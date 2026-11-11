# Refusal revamp: refuse only on physical failure, measured on the cell

Written 10 November 2026 against 5.49.0. Grounded in audit/CONSULT-over-refusal.md and its response.

## The rule

A claim is withheld only when that cell's own data shows the claim cannot be measured: the DNA is
not there, the calls are at chance, or a second genome is present where the claim needs one. Every
other weakness is reported as a confidence on a claim that is made, never as a silent refusal.

Every withheld claim carries exactly one physical cause from a closed list, and the number that
showed it. A check enforces that. This is what ends "why was this rejected": the answer is always
printed next to the refusal.

## The closed list of causes

| cause | physical meaning | measured by |
|---|---|---|
| `no-template` | the DNA is not there, on the array or on this chromosome | foreign-allele rate at same-homozygote markers; call collapse |
| `no-genotype-signal` | calls are at chance: hybridisation or clustering failed | information ratio R; parent-free drop-in |
| `second-genome` | another person's DNA is in the reaction | foreign-allele enrichment E (once validated) |
| `no-intensity-reference` | nothing in this genome can anchor intensity | allele-disomic chromosome count |
| `parents-uninformative` | the family offers too few informative markers here | offered count, before the sample is read |
| `parent-not-loaded` | the claim needs a parent array that was not supplied | input |

Anything else that currently refuses is a proxy, and becomes either a confidence on the claim or a
statistical test whose result is shown.

## What gates claims today (inventory)

**A. Cell-level gates that erase every claim.**
`stage.ts` `inferStage` returns `failed` on heterozygosity over `MAX_DIPLOID_HET` (0.25, scaled) or
call rate under `QC_CALL_FLOOR` (0.40). `scoreSample` then clears `aneuploidy`, `lostParent`,
`aneuploidyParent`, `mixed`, `alleleGain` on every chromosome, empties `findings`, nulls relatedness
and withholds the copy-number scan. `integrity.ts` adds three `uninterpretable` levels: stage failed,
`COMPLEX_CALL_RATE` (0.70), `COMPLEX_DEVIANT_FRACTION` (0.30).
Heterozygosity is a PROXY: the consult measured that the refused cells it was blamed for are single
genomes with drop-in of 0.220 to 0.307, whose allele calls agree with published truth at 133 of 138.

**B. Channel gates.** Copy-number scan withheld on a chaotic genome (`HET_SPREAD_LOCALISED`).
Intensity really does break there (AUC 0.67), so the gate is right in kind; its trigger is wrong.

**C. Per-chromosome floors.** `obligateHet` `MIN_INFORMATIVE` 200, the undecided band between
`HAPLOID_MAX` 0.20 and `DIPLOID_MIN` 0.45, `SURVIVOR_MIN`, `MIXED_SHARE`; `gainOrigin`
`MIN_INFORMATIVE_*`, `SHARE_MARGIN`; `parentage` `CALL_RATE_FLOOR` 0.60, `CALL_COLLAPSE`, the
`HET_BAND_*` family; `dosageOrigin` `MAX_HET_BAF_SD` 0.11, which consult C showed refuses the wrong
group.

**D. Statistical tests.** `SEGMENT_LRT`, the `Z_*` constants, posterior bands. These are evidence,
not floors. They stay, and the value is shown when they decide.

## Phases

### 0. Corrections (now)
- "Two cells" was wrong in five places: it is 28 calls over five cells published chaotic before the
  gain arm, 31 over seven after. Fix CHANGELOG 5.49.0, the result file, the consult (as a correction
  note, since it was sent), and Joy's packet.
- Joy's packet gets the consult's finding on the refused cells.

### 1. Refusal census: measure before changing anything
- `result.refusals`: one record per withheld claim, `{ claim, chrom, gate, cause, measured,
  withheld }`, pushed at the point each gate clears a value, so the withheld value is kept and the
  census can score it. The UI renders cause and number only.
- Census harness built on the REAL `scoreSample` (pattern of audit/zuccaro-table-s3.ts), never a
  re-implementation, because a harness that diverged from the product produced every gains figure
  this project recorded until 10 November.
- Scored against Zuccaro Table S3 with the consult's metric: by claim type and direction; a false
  refusal is a withheld claim that published truth confirms AND that the claim's own channel made;
  wrong parent is its own column; an adjudicate bin for claims that two independent channels support
  against the hand reading.
- Output: per gate, claims withheld, confirmed by truth, contradicted by truth, unscorable.
- TESTS: a source-level pin that every assignment clearing a claim in the refusal blocks sits beside
  a `refusals.push`, and that `cause` is a member of the closed list. A unit check that a constructed
  failed cell produces one refusal per cleared claim with the withheld value attached.

### 2. Answer keys: the physical measurements, per cell
New module, each statistic with a known expectation:
- drop-in: heterozygous calls at markers where both parents are the same homozygote (trio)
- parent-free drop-in: the same at markers homozygous-identical across a founder panel. Founders:
  laboratory bulk arrays S23-00 and S20-18 on tesla plus unrelated GSE148488 parental bulk
  arrays, substituting for the consult's two missing founders. The substitution is VALIDATED, not
  assumed: correlation with trio drop-in at 0.99 or better over the trio cells, as the consult got
  0.998 over 76.
- allele contrast k: 1 - 2 x foreign BAF at same-homozygote markers
- template-absent signature per chromosome: foreign rate against the cell's own median
- allele-disomic fraction: the replacement for "normal majority", which flat noise cannot pass
- corrected heterozygosity: het minus drop-in, the ploidy-relevant quantity
- E (mixture) is NOT shipped until validated by construction with a floor under 5% admixture
- TESTS: each on constructed arrays with a known planted value, recovered within tolerance; and the
  consult's measured values as known-answer targets (5.42.1 cells drop-in 0.220 to 0.307; corrected
  het 0.015 to 0.053; the five noise cells drop-in over 0.50).

### 3. Per-claim gating
Replace the erase-all block with each claim checking only what it consumes:

| claim | withheld only on |
|---|---|
| which parent, whole-chromosome loss | `no-genotype-signal`, `second-genome`, `no-template`, `parents-uninformative`, `parent-not-loaded` |
| which parent, whole-chromosome gain | the above, plus no copy evidence |
| whole-chromosome copy number | `no-intensity-reference` |
| segments | `no-intensity-reference`, `no-template` |
| stage and ploidy | `no-genotype-signal` |

OPEN CONFLICT, to settle by measurement and not by argument: the consult's allele calls on chaotic
cells are referenced to 0.5 with k, scoring 190 of 196 against 183 of 196 own-median. gainOrigin's
docstring says referencing to 0.5 against a RECONSTRUCTED parent is biased maternal by up to 0.077.
Both are true in their own setting. Plan: 0.5 with k when both parents are measured arrays; own
median when either is reconstructed; each path tested.
- TESTS: independence checks. A cell at R near 1 keeps intensity claims where anchored. A chaotic
  cell keeps allele claims. The five noise cells keep refusing every genotype claim. A constructed
  high-drop-in single genome gets its allele claims back with the right parent.

### 4. Proxy floors become confidence
Group C floors that are proxies stop hiding calls. A call under `MIN_INFORMATIVE` is made with its
uncertainty from the markers it has; the undecided band stays undecided but shows the posterior. A
floor survives only if the census shows it withholds claims truth contradicts.
- TESTS: each converted floor gets a boundary check that a call either side carries a confidence
  and never silently disappears.

### 5. Pre-registered validation
The consult's P1 to P6, plus the census bars below, written in their own commit before any is run.

### 6. Ship, and one page for Joy
Release; and a single page listing the six causes, what each means physically, and the number
printed beside it, so she can read any refusal herself.

## Bars, fixed now
- wrong parent: 0 on the Zuccaro arm
- false segments: 1 or fewer on the 72 laboratory arrays with no hand call
- contradicted claims on the 477 published-normal autosomes: 34 or fewer, with the adjudicate list
  published
- false refusals in the census: strictly fewer than before
- true refusals kept: the five cells at R 1.01 to 1.12 refused for every genotype claim
- every withheld claim carries one cause from the closed list and its measurement, check-enforced

## Order
0 and 1 now. 2 next, since every later phase consumes it. 3 to 6 each pre-registered in its own
commit before it is scored. Nothing in 3 or 4 ships before the census exists to measure it against.

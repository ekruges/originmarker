<!-- Laboratory sample identifiers redacted: this file is tracked in a public repository.
     Cells of the laboratory series are pseudonymised per audit/local/barcode-map.tsv (S12, S13).
     GSM accessions are public GEO identifiers and are given as they are. -->

# Consult: the tool refuses cells it could say things about, and every fix has missed why

Written 10 November 2026 against version 5.49.0, commit c71b64a.

## The request

Our collaborator in the Egli laboratory has pushed, repeatedly and across several months, for the
tool to say more about more cells. She is right, and the measurements below show it. What I am
asking for is not a lower threshold. It is a structure in which a true refusal stays a refusal and
everything a cell CAN support still gets said.

The developer's own framing, which I think is the correct one: "not forcing true refusals through
but helping say more stuff about more cells."

## 1. The architecture that produces the over-refusal

One stage verdict gates every claim. `inferStage` (web/src/stage.ts) returns `failed` on either of:

- heterozygosity over `MAX_DIPLOID_HET` = 0.25 (scaled by material), or
- call rate under `QC_CALL_FLOOR` = 0.40.

When it does, `scoreSample` (web/src/scoreSample.ts, the `if (uninterpretable)` block) clears, for
every chromosome: `aneuploidy`, `lostParent`, `aneuploidyParent`, `mixed`, `alleleGain`. It also
withholds the whole copy-number scan (`const copy = (uninterpretable || chaoticGenome) ? [] : ...`).
So a cell that fails the heterozygosity test loses its intensity-based claims, its allele-based
claims and its segment claims together, although those rest on different evidence and fail for
different reasons.

## 2. The history, which is the reason for this consult

Every round asked whether a threshold was right and answered that correctly. None asked whether a
single upstream gate should erase every downstream claim.

| release | date | what it did | effect on refusal |
|---|---|---|---|
| 4.8.4 | 14 Aug | introduced `MAX_DIPLOID_HET = 0.25` because a 52%-heterozygous mixture was reading as bulk DNA | fixed an UNDER-refusal, and made "too heterozygous" mean "withhold everything" |
| 4.8.1 | 14 Aug | `QC_CALL_FLOOR`, so a failed amplification is not called haploid | more refusal, correctly |
| 5.21.0 | 31 Aug | scaled the ceiling and dropout anchor by material | threshold tuned |
| 5.42.1 | 26 Oct | **tested the ceiling as an over-refusal** (see below) | defended; no change |
| 5.46.0 | 2 Nov | every withheld call states why; a diploid cell says diploid | explains, does not reduce |
| 5.47.0 | 4 Nov | says whether the array failed or the chromosome had nothing to read | explains |
| 5.48.0 | 9 Nov | information ratio R: says when calls carry no information | explains, and makes a NEW separation possible that nothing yet uses |
| 5.49.0 | 10 Nov | gain nomination on intensity; found the benchmark was scanning refused arrays | exposed how much the gate costs |

**5.42.1 is the one to read first** (CHANGELOG.md, entry 5.42.1). Five cells of GSE148488 refused
at heterozygosity 0.318 to 0.376, whose published karyotypes read two of them "normal for other
chromosomes" and one "46,XX", carrying three hand-called `del(6q)(pat)` deletions lost to the
refusal. The hypothesis tested was that drop-in on amplified material is larger than on bulk, which
would license raising the ceiling. Measured on GSE60909 (twenty amplified single cells of two
siblings' lymphoblastoid lines against unamplified controls of the same genomes), it was not, so
the ceiling stood and the five cells stayed refused. That conclusion is correct for the question it
asked. It is the wrong question.

## 3. What is measured today

### 3a. The 21 laboratory cells carrying hand calls

All 8 refusals fire on the heterozygosity branch; none on the call-rate floor. Scored against
5.48.0's information ratio R, which does not depend on heterozygosity (R near 1 is a sample whose
calls carry no information about any genome; a two-person mixture reads BELOW a single genome):

| cell | stage today | call | het | R | hand gains | hand losses |
|---|---|---|---|---|---|---|
| S12-11 | failed | 0.508 | 0.636 | 1.010 | 7 | 0 |
| S12-22 | failed | 0.608 | 0.619 | 1.107 | 1 | 0 |
| S13-39 | failed | 0.622 | 0.616 | 1.082 | 2 | 0 |
| S12-20 | failed | 0.605 | 0.615 | 1.122 | 1 | 0 |
| S12-19 | failed | 0.612 | 0.596 | 1.115 | 1 | 0 |
| **S12-60** | failed | 0.569 | 0.458 | **0.429** | 1 | 0 |
| **S13-47** | failed | 0.547 | 0.353 | **0.410** | 9 | 14 |
| **S12-05** | failed | 0.562 | 0.262 | **0.323** | 1 | 9 |
| S12-67 | trophectoderm | 0.587 | 0.202 | 0.356 | 0 | 8 |
| S13-44 | trophectoderm | 0.616 | 0.198 | 0.363 | 3 | 11 |
| S13-06 | trophectoderm | 0.717 | 0.150 | 0.305 | 0 | 1 |
| S13-08 | blastomere | 0.767 | 0.108 | 0.330 | 0 | 5 |

(Nine further cells scan today at R 0.220 to 0.334 and hold 0 gains and 18 losses between them.)

The first five are correct refusals by an independent statistic. **The bolded three are not**: S12-05
is refused at R 0.323 while S13-06 and S13-08 pass at 0.305 and 0.330. Those three hold 11 of the 26
hand-called gains and 23 of the 67 losses.

### 3b. The replacement I tried, and why it failed

The copy-number scan normalises against the sample's own per-chromosome median, so its real
requirement is a normal majority. I measured that: the share of autosomes whose median log R lies
within 0.15 of the genome median, with the majority threshold at 0.5. Script:
audit/joy-normal-majority.ts.

Against Zuccaro Table S3 (31 cells, published per-chromosome karyotypes) it is nearly perfect:
12 of 12 cells published chaotic are refused, 18 of 19 others admitted, the one error a cell
published 45,XY,-6 sitting at exactly 0.500.

On the laboratory cells it fails in both directions. It ADMITS the four cells at R 1.08 to 1.12,
whose calls are noise and whose intensity happens to be flat (concordance 0.955), and REFUSES eight
cells that scan today, including one at R 0.220 and heterozygosity 0.136 (concordance 0.273). On
heavily amplified single cells, per-chromosome intensity swings without copy number. Intensity
concordance and genotype information are orthogonal, so neither substitutes for the other.

### 3c. The published arm, and an under-refusal

Zuccaro 2020 Table S3 has 477 autosomes published normal. The tool raises 30 segments on them
before 5.49.0's gain arm and 34 after. 28 and 31 of those sit on two cells the paper hand-called
chaotic aneuploid, which the tool does NOT refuse. So the same gate is too strict on the
laboratory cells and too loose on these two.

## 4. What I am asking

**Q1. The structure.** Is per-claim gating the right replacement for one stage gate, and what does
each claim actually require? The tool's claims, roughly:

- whole-chromosome copy number from intensity
- which parent a whole-chromosome loss came from, from alleles at opposite-homozygote markers
- which parent supplied a whole-chromosome gain, from paternal share
- segmental copy number from the intensity scan
- stage and ploidy

For each: which of these break it, and which do not: (a) calls carrying no genotype information,
(b) a two-person mixture, (c) a triploid, (d) a genome with no normal majority, (e) heavy
amplification noise in intensity. My guess is that (a) breaks every allele claim and no intensity
claim, (d) breaks only claims normalised against the sample's own median, and (b) and (c) break
"which parent" in ways that depend on the mixing fraction. I would like that checked rather than
assumed.

**Q2. The over-ceiling, informative class.** Cells over the heterozygosity ceiling whose calls DO
carry information (R 0.32 to 0.43 here; 0.318 to 0.376 heterozygous in 5.42.1's five GSE148488
cells) are candidates for a two-person mixture, a triploid, or something else. These should be
separable POSITIVELY from BAF band structure rather than refused: a triploid shows bands near 1/3
and 2/3, a mixture shows extra bands whose positions give the mixing fraction. Is that sound on
amplified single cells at these call rates, and what would the measurement be? Test material: the
three laboratory cells above, the five 5.42.1 cells (published normal), and any known mixture or
triploid you can find in public data on this platform.

**Q3. Segmental calls without a normal majority.** Consult A (CONSULT-RESPONSE-segmental-gain-
statistic.md, section 2) showed a per-chromosome offset built from one series' clean arrays is real
and improves 30 of 36 arrays within a laboratory, and does not transfer between laboratories. The
laboratory has 72 arrays with no hand call on the same platform. Would an external per-series
reference remove the scan's need for the sample's own majority, and does it survive amplified
single cells? Validation I can run: leave-one-out on the 72, and the Zuccaro arm.

**Q4. The two chaotic cells the tool admits.** What should it say about them?

**Q5. Making "says more" measurable.** I want a metric to optimise that cannot be gamed by
forcing refusals through. My proposal: for each cell, the claims made, each scored against
published truth where truth exists; over-refusal is then the number of claims withheld that
published truth would have confirmed, and over-calling the number made that it contradicts. Is
that the right pair, and what is the right bar between them?

## 5. Constraints

- The 72 laboratory arrays with no hand call are the specificity arm and are never traded for
  recall. One false segment across them is the standing bar.
- The published arm (477 normal autosomes of Zuccaro Table S3) must not get worse than today.
- Every change is pre-registered with its pass criteria fixed before it is scored, and a failed
  prediction is reported as a failure.
- The laboratory region table's provenance is unsettled (all 93 rows reproduce from this
  project's own output by filtering), so recall against it is reported but specificity rests on
  published data.
- Bulk reference arrays are not available in the laboratory corpus, which is why 5.48.0's
  "informative" band does not ship as a positive claim.

## 6. Files

- web/src/stage.ts (`inferStage`, the `failed` branches), web/src/scoreSample.ts (the
  `uninterpretable` block and the copy-number withholding), web/src/informationRatio.ts (R)
- CHANGELOG.md entries 4.8.4, 5.21.0, 5.42.1, 5.46.0, 5.48.0, 5.49.0
- audit/gain-nomination-result-2026-11-10.txt and its addendum
- audit/joy-normal-majority.ts, audit/joy-public-specificity.ts, audit/zuccaro-table-s3.tsv
- audit/noise-information-ratio-2026-11-08.txt

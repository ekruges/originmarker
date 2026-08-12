# Accuracy audits

Two, one per feature that makes a call from array data. Each runs against public arrays whose
answer was established at the bench rather than by any statistic the tool computes.

| | correct | refused | incorrect |
|---|---|---|---|
| [Syngamy](AUDIT.txt) | 22 | 4 | **0** |
| [Progenitor](PROGENITOR-AUDIT.txt) | 35 | 6 | **0** |
| [Breakpoints and gain origin](BREAKPOINTS-AUDIT.txt) | 3 of 4 parts | 1 blocked | **0** |

The third audit is newer and narrower. Three of its four parts run on the public arrays in this
repo: breakpoint accuracy against constructed events with known truth (12.2x better than the
window edge), a positive control that recovers both directions of a gain on a real array, and
specificity. The split-half holdout has nothing to run on, because no segment clears the shipped
floor in a subset this sparse. It also carries an independent corroboration of a real chr1 event
from a statistic unrelated to the one that first found it. What is missing is named in the file
rather than left as an absence.

---

# Syngamy

Syngamy, run against 27 public arrays whose answer was established at the bench rather than by
any statistic this tool computes.

**[AUDIT.txt](AUDIT.txt)** is the record: every case, what the experiment says, what the tool
said, the numbers behind it, and a remark. **[reports/](reports)** holds the report PDF for each
of the 26 scored cases, produced by the same builder as the app's Report (PDF) button.

## Result

| | |
|---|---|
| correct | 22 |
| refused | 4 (declined to call, asserted nothing) |
| **incorrect** | **0** |

22 of 23 correct on samples that pass every one of the tool's own quality gates. The
egg-donor identification matrix resolved 7 of 7 maternal pronuclei to exactly one of four
donors, with the runner-up 10 to 18 times further away.

## Why this data

Zuccaro et al. separated the two pronuclei of a fertilised human zygote by micromanipulation and
arrayed each one alone. A "paternal nucleus" sample contains one parental genome and a "maternal
nucleus" the other, physically, before any statistic is applied. That is the answer to the
question Syngamy asks, established independently of it.

The series also carries four arrays of one sperm donor and eight of four egg donors, which give
identity and unrelated controls on the same platform from the same laboratory.

Source: GEO [GSE148488](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148488).
Zuccaro MV et al. *Cell* 2020;183(6):1650-1664. doi:10.1016/j.cell.2020.10.025

## What it changed

The first run returned **3 incorrect calls**, not 0. Three isolated paternal pronuclei were
reported as biparental: they are haploid and therefore homozygous by construction, but each
showed an 18-27% heterozygous BAF band, which pushed zygosity to "diploid". All three sat below
the 60% call rate at which erroneous heterozygous calls become common, and the tool had already
excluded all three on that gate, and its het-to-hom gate had already printed *"SUSPENDED: below
60% call rate erroneous heterozygous calls are common"* - and the classifier used the
heterozygous band anyway.

Zygosity is now withheld below that floor, in both `origin.py` and `web/src/parentage.ts`, and
the class falls to unclear with the reason stated. Absence is unaffected and still called,
because it is Mendelian. Three wrong answers became three refusals; every other case is
unchanged.

## Reproducing it

The harness drives `web/src/ingest.ts`, `web/src/parentage.ts` and `web/src/syngamyPdf.ts`, which
is what the browser runs. Only the file plumbing differs: a gzipped file off disk rather than a
dropped one, and `favicon.svg` read rather than fetched.

```sh
mkdir -p raw && cd raw
# the 27 GSM accessions are listed under SOURCES in AUDIT.txt
curl -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4472nnn/GSM4472397/suppl/GSM4472397_sperm_DNA_71.CEL.txt.gz
cd ..
node --experimental-strip-types --max-old-space-size=8192 audit/cases.ts raw out
node --experimental-strip-types audit/report.ts out
```

## What it does not show

One laboratory, one platform, one sperm donor. Single pronuclei are harder material than a
trophectoderm biopsy, so this is a stress test rather than a representative one. 26 scored cases
cannot support an accuracy percentage: 22 correct with 0 incorrect puts a 95% lower bound near
87% on the proportion that are not incorrect, and that is as far as it goes.

Research use only. Not a clinical diagnostic.

---

# Progenitor

A man's genotype reconstructed from eight of his sons' paternal pronuclei, then scored against
his own array, which the reconstruction never saw.

**[PROGENITOR-AUDIT.txt](PROGENITOR-AUDIT.txt)** is the record.

## Result

| | |
|---|---|
| correct | 35 |
| refused | 6 (three products set aside on call rate, three sub-floor cases recorded rather than scored) |
| **incorrect** | **0** |

He reads present at 0.12x and 0.18x on his two replicate arrays. All nine genomes that are not
his read absent, at 4.71x to 9.11x. Membership put the five usable products in one group and
split them cleanly from seven maternal pronuclei, with zero misread pairs out of 45.

## What makes it checkable

The series carries the sperm donor's own bulk DNA alongside the pronuclei. That turns three
modelled quantities into measured ones: ascertainment against his real 16.66% heterozygosity,
contamination by asking his array which asserted homozygotes are really heterozygous, and
absence on nine genomes known not to be his. His two replicate arrays bound what the ground
truth itself can be wrong about, at 4.15%.

## What it found

Predicted contamination sits below observed at every depth, by 0.6 to 1.1 percentage points.
That gap is inside the replicate floor, so this data cannot separate a model that understates it
from a truth that overstates it. It is reported rather than resolved, and the direction is
stated: if the model is the one that is wrong, it is optimistic.

The density sweep rebuilds the reference on disjoint slices of the array, down to every 16th
marker: 31 references, 310 classifications, all clean. One array density is one observation, and
a reconstruction that survives a sixteenth of this one is not leaning on density.

## Reproducing it

```sh
node --experimental-strip-types --max-old-space-size=8192 audit/progenitor.ts raw out
node --experimental-strip-types audit/progenitor-report.ts out
```

## What it does not show

One laboratory, one platform, one man, and five usable products is the floor rather than a
comfortable margin. The five-product floor was set on true offspring inverting below it; this
series has no offspring of these products, so the sub-floor cases score the parent instead and
do not reproduce that measurement. Nothing here separates this man from a close relative, since
none is in the series.

Research use only. Not a clinical diagnostic.

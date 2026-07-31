# Accuracy audit

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

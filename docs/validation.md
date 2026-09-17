# Validation record

What this tool has been measured against, what each measurement establishes, and where each one
stops. Every answer below was established independently of this tool: by a karyotype, by an
experimental design, by a dissection performed at a bench, or by a submitting laboratory's own
metadata. Nothing here is scored against a number this tool produced.

Supersedes `validation-2.0.md`, which covers the 2.0 inference path only.

> **Research use only. Not a clinical diagnostic.**

## The corpus

| series | n | material | where the answer comes from |
|---|---|---|---|
| GSE148488 | 135 | zygotes, pronuclei, embryos, donor gDNA | pronuclei separated by micromanipulation; EYS targeted |
| GSE186407 | 166 | single blastomeres | chromosome 16 targeted, stated per file |
| GSE290961 | 158 | single blastomeres | preimplantation embryo, no per-sample answer |
| GSE19247 | 622 | sperm, lymphoblasts, blood, cleavage embryos | karyotype-confirmed trisomy 21; euploid reference pedigrees |
| GSE117672 | 13 | complete hydatidiform moles | androgenetic, dispermic; placental sex; 3 carry a trisomy |
| GSE60909 | 164 | PGD families, blood and blastomeres | mother, father, sibling labelled by the submitters |
| GSE20975 | 82 | single cleavage-stage blastomeres | 24-chromosome karyotype per blastomere, from the submitters' array method |
| GSE18932 | 206 | blastocyst biopsies, Coriell cell lines, two-person mixtures | karyotype per sample; Coriell isodisomy karyotypes; stated mixing proportions |
| GSE12713 | 100 | complete hydatidiform moles | monospermic, stated by the series for its whole collection |
| GSE207887 | 189 of 564 | miscarriage chorionic villus and fetal skin, BULK | one stated abnormality per sample, trisomies, 45,X, triploidy, or free of variation |
| GSE163799 | 23 | fetal amniotic and villus, BULK | a stated microduplication or microdeletion per sample, with its band |
| GSE21732 | 51 | embryos of a balanced t(2;20) carrier, plus both parents and the newborn | one ISCN karyotype per sample, agreed by real-time PCR, SNP array and FISH |
| laboratory corpus | 884 | mixed | targeted loci, from experimental design |

Roughly 2,800 arrays on six array platforms (GPL28377, GPL6985, GPL8855, GPL13829, GPL3718,
GPL18637), spanning about twenty years of instrument generations, from four laboratories on three
continents.

## What is established

### Detection does not invent findings

| set | observations | false calls |
|---|---|---|
| adult gamete donors, whole-chromosome | 242 autosomes | **0** |
| euploid reference lines and donors, both gates | 1,584 autosomes | **1** (0.63 per 1,000) |
| structural validity across the corpus | ~2,500 arrays | **0** malformed events, **0** crashes |
| one-complement genomes, ploidy call withheld | 144 arrays | **0** heterozygosity-loss findings |
| monospermic complete moles, whole-chromosome | 100 arrays | **0** |
| normal-karyotype blastocysts the tool accepts | 2 arrays | **0** |

A living gamete donor carries no autosomal whole-chromosome loss, so every call on those arrays is a
false one. Structural validity is checked separately from correctness: a start past its own end,
coordinates off the end of the chromosome they claim, a non-finite number in a reported field, a
call with no markers behind it, a duplicated locus. None occurred on any chemistry.

**Heterozygosity-loss findings stand on the array, not on the ploidy call.** Copy-neutral LOH and
uniparental disomy exclude a genome the tool calls uniparental. A genome whose ploidy was refused
passes that exclusion, so `audit/upd-premise.ts` drives both detectors with the ploidy call withheld,
on genomes whose answer biology fixes:

| set | arrays | findings, 5.29.0 | findings, 5.30.0 |
|---|---|---|---|
| GSE12713 monospermic moles | 100 | 2,200 | **0** |
| GSE19247 single sperm | 29 | 384 | **0** |
| GSE148488 pronuclei | 15 | 163 | **0** |
| GSE18932 stated isodisomy, stated chromosome | 8 | 8 | 7 |

Every finding on a one-complement genome is false. Through 5.29.0 the run detector measured a run
against the run-finder's fixed 5% tolerance rather than against the array, so an array near that
heterozygosity read as runs end to end. It now carries the three bounds the copy-neutral detector
already had. The isodisomy array that drops expects 380 heterozygotes on the stated chromosome at
its own background, under the measured bound of 400, and the stage inference refuses it at a 38.1%
call rate before anything is reported.

### Detection finds events where an experiment put them

| series | targeted | rank among 22 autosomes | enrichment |
|---|---|---|---|
| GSE148488, EYS | chromosome 6 | **1st** of 22, over 108 embryos | 2.9x |
| GSE186407 | chromosome 16 | **1st** of 22, over 128 blastomeres | 7.9x |

The tool is never told where the experiment aimed. The other 21 autosomes are a specificity set on
the same arrays at the same time. The chromosome 6 figure is taken over every embryo in the series;
the first 40 alone read 7.3x, which overstated the effect.

### Attribution names the right parent

| measurement | result |
|---|---|
| segmental breakages, trio-resolved | **8 of 8 correct, 0 wrong** |
| across all parental calls made to date | **0 wrong in 263** |
| one-complement genomes (pronuclei) | **0 parents named on 15** |
| mirrored constructed losses | 12/12 paternal, 10/12 maternal, 0/12 on the untouched control |

Trio truth is Mendelian and computed outside the tool: where one parent is homozygous and the other
homozygous for the opposite allele, a normal child must be heterozygous, and which allele survives
says which copy went. Margins were 479:1, 1115:8 and 235:1.

That derivation independently reproduces the source study's own finding, that the chromosome 6
losses are paternal, without being told it.

**Both directions, on twelve families.** GSE60909 labels mother, father and embryos per array.
Where one parent is homozygous and the embryo heterozygous, the embryo's other allele must be in the
other parent's array:

| direction | informative markers | carried by the other parent |
|---|---|---|
| paternal | 145,591 | 78.52% |
| maternal | 146,599 | 78.67% |

The two agree within every family as well as in aggregate: 47.06% against 46.94% in the family with
the most dropout, 94.70% against 94.61% in the cleanest. The spread between families is embryo
quality, not a parental asymmetry. Pairs drawn from different families called relatives: 0 of 12.

### The panel builder keeps the markers that hold up

`informativity.ts` labels which markers can testify about an embryo's paternal contribution. The
label is computed from the father's first array and checked against his second, so his genotype is
not on both sides of the comparison. A contradiction is the second array lacking the allele the
mother could not have supplied:

| label | kept, contradicted | discarded, contradicted | ratio |
|---|---|---|---|
| provesPaternalPresence | 1,314,496 markers, 1.70% | 528,654 markers, 93.14% | 0.018 |
| keySnp | 725,378 markers, 2.91% | 1,117,772 markers, 44.16% | 0.066 |

The father's two arrays disagree about heterozygosity at 4.152% of jointly called markers, which
bounds how clean any kept set can read.

### A parent can be reconstructed from haploid products

Built from 8 paternal pronuclei at a marker depth chosen without ever seeing him, then scored
against his own arrays, which were never given to the reconstruction:

| | |
|---|---|
| verdicts | **35 correct, 6 refused, 0 incorrect** |
| his own two arrays | present, 0.0277 and 0.0431 |
| nine genomes that are not his | absent, 0.1379 to 0.1482 |
| genotype concordance at the chosen depth | 96.637%, opposite homozygotes 0.1140% |

The record states its own ceiling: his two arrays disagree with each other about heterozygosity at
4.152% of jointly called markers, and no comparison against him can resolve below that.

### The same DNA gives the same answer

| | |
|---|---|
| zygosity, integrity, whole-chromosome set | **14/14** |
| segment set | 13/14 |
| named parents, by chromosome | 13/14 |
| one locus given two different parents | **0** |
| informative-marker calls stable across replicate arrays | 99.6981% |

### Isodisomy is found on a stated chromosome

Coriell GM11496 carries isodisomy of chromosome 7 and GM15603 of chromosome 8, four arrays each in
GSE18932, genotype calls only:

| | |
|---|---|
| stated chromosome least heterozygous of 22, in the data | 8/8 |
| arrays refused by the call-rate floor, at 36.4% to 39.4% | 4/8 |
| stated chromosome flagged, on the arrays the tool accepts | **4/4** |
| false findings on those four arrays, 5.29.0 then 5.30.0 | 17, then **0** |

A refused array reports nothing and says so. The floor is not moved to reach the other four, because
moving a quality gate to fit known answers is fitting the test. Without intensity a single copy reads
the same as isodisomy, and the finding says it is either. The 17 false findings were 13 segmental
uniparental disomies on other chromosomes and 4 complex-genome calls made from call rate alone; a
low call rate is now reported as a damaged array rather than as a disturbed genome.

### Karyotyped single cells sit under the quality floor

GSE20975 and GSE18932 publish a karyotype per blastomere or biopsy, called by the submitters' own
array method. Their genotype calls on GPL3718 were lifted from NCBI Build 35 to GRCh37, with 0
illegal placements of 262,208 markers.

| series | listed monosomies | refused | normal arrays | refused | findings on accepted normals |
|---|---|---|---|---|---|
| GSE20975 | 15 | 15 | 60 | 60 | 0 of 0 |
| GSE18932 | 32 | 32 | 131 | 129 | **0 of 2** |

The monosomic chromosome is plainly in the data, least heterozygous of 22 on 44 of 47, but every one
sits on an array under the 40% call-rate floor. Sensitivity cannot be scored on this material under
the tool's own gate, and no claim is made from it. The two accepted normal arrays each carried a
complex-genome call through 5.29.0. Sex is not callable on this platform, which carries too few
chromosome Y markers, and the tool refuses rather than guessing.

### Mixed samples are not detected

Four known-proportion mixtures of GM00323 (male, Cellosaurus CVCL_7280) and GM00321 (female,
CVCL_7279), with each line pure as its own control:

| | arrays | stage refused | sex call ambiguous | either |
|---|---|---|---|---|
| pure | 2 | 0 | 0 | 0 |
| mixed | 4 | 0 | 3 | 3 |

One signal responds to a mixture without identifying one, and the mixture that is 75% female trips
nothing. Two of these mixtures were refused by the stage inference until 5.31.0, and that refusal is
not evidence of anything: the ceiling refusing them also refuses undamaged bulk arrays on this
panel. No mixture gate ships, because none can be validated here: checks for excess heterozygosity
detect only contamination above 5 to 10%, and the array method reads B-allele frequencies against
population allele frequencies (Jun et al. 2012, Am J Hum Genet 91:839). Every sample reports the
mixture question as not tested.

### Sex is called correctly

| set | expected | observed |
|---|---|---|
| hydatidiform moles | 7 male, 6 female | **7, 6** |
| adult donors | 4 male, 8 female | **4, 8** |
| maternal pronuclei | zero Y-bearing | **0 of 7** |

An egg carries an X and no Y, so a Y call on a maternal pronucleus is an error and the count is the
error count.

One deviation is recorded rather than explained away: 23 sperm cells on the converted GSE19247
platform read 17 Y-bearing against the 1:1 that meiosis fixes, exact two-sided p = 0.0347.

### What each class emits, on 600 real arrays

`audit/coverage.ts` drives 600 arrays from five sets through the whole pipeline and counts what each
class of the taxonomy emits. These are counts, not accuracy: most emitting arrays carry no
per-sample answer for these classes. Between 5.29.0 and 5.30.0 only the three classes this release
changed moved:

| class | arrays, 5.29.0 | arrays, 5.30.0 | instances, 5.29.0 | instances, 5.30.0 |
|---|---|---|---|---|
| segmental uniparental disomy | 282 | 36 | 2,976 | 178 |
| whole-chromosome isodisomy | 61 | 15 | 244 | 35 |
| complex or chaotic genome | 170 | 12 | 170 | 12 |

Arrays of refused ploidy emitting a disomy fell from 165 to 3, and GSE19247, the series of single
sperm, lymphoblasts and cleavage embryos, went from 181 disomy and 158 complex-genome arrays to none.
Every other class is identical in both releases: monosomy on 72 arrays, trisomy on 5, segmental
deletion on 162, copy-neutral LOH on 12, haploidy on 26. The 12 adult donors carry no class a
healthy adult cannot, in either release. Three classes still fire on no real array: segmental
duplication, triploidy, and a segmental change introduced by the gamete.

The full regression battery, 19 harness runs on the same inputs, reads identically across the two
releases in 17 of 19 logs. The other two differ only by the new mixture gate and by a smaller
serialised result where the removed findings used to be.

### Two constants were platform properties

Both were measured on one chemistry and applied to every other, and each refused real material:

| panel | bulk heterozygosity |
|---|---|
| Axiom GPL28377 | 0.141 to 0.170 |
| Affymetrix 250K Nsp GPL3718 | 0.201 to 0.308 |
| Illumina CytoSNP-12 GPL13829 | 0.301 to 0.316 |

The stage inference refused any array over 25% heterozygous as "not a genome", a figure derived from
the Axiom rate plus drop-in. On the other two panels an undamaged bulk array sits above it: a
newborn's own gDNA, 96.9% called and 30.8% heterozygous, was refused. The ceiling and the dropout
anchor now follow the panel, taken from the parental array, which is on the same panel as the
sample. They scale up only, because a parent's own dropout can only depress its heterozygosity, and
they are capped at twice the anchor.

The stage RUNGS are deliberately left alone. Scaling those as well moved one bulk array of the
stress corpus onto an amplified rung, where the relationship verdict is withheld, and the tool
stopped warning that an array had been placed in its own parental slot on one run in eight. That
regression was caught by the battery and is pinned by a check.

### A published method's own calls, on its own data

GSE21732 is a PGD case report: 48 embryos of a mother carrying a balanced t(2;20), both parents
arrayed, one ISCN karyotype per sample agreed by real-time PCR, SNP array and FISH. It is the only
set here whose answers were produced by another method rather than by a bench dissection.

| | |
|---|---|
| arrays the tool refuses, at 19% to 43% call rate | 45 of 51 |
| listed whole-chromosome events called, on accepted arrays | 2 of 6 |
| normal and balanced-carrier arrays, no call anywhere | 4 of 4 |
| derivative chromosomes, which are maternal by the parents' karyotypes | 0 of 4 reported |

The comparison is narrow because the material is not: these are 2010 arrays of trophectoderm
biopsies and arrested whole embryos, and the tool refuses most of them on call rate. What it does
establish is that on the arrays it accepts it invents nothing, and that the derivative segments,
which are 20 to 90 Mb of imbalance, are not reported: this platform's segment floor is 26.1 Mb and
those events sit near it.

### Bulk material, with the answer stated per sample

Every sensitivity figure in this record used to come from amplified single cells. GSE207887 is
miscarriage material on Affymetrix CytoScan 750K, unamplified, with one stated abnormality per
sample:

| stated | called |
|---|---|
| trisomy 13 | 20 of 21 |
| trisomy 16 | 18 of 20 |
| trisomy 18 | 7 of 7 |
| trisomy 21 | 18 of 18 |
| trisomy 22 | 20 of 20 |
| **whole-chromosome gains, all** | **83 of 86, 96.5% [90.2%, 98.8%]** |
| triploidy | 20 of 23, 87.0% [67.9%, 95.5%] |
| stated free of clinical variation | **58 of 58 clean**, 2 refused |
| 45,X | 17 of 20, 85.0% [64.0%, 94.8%] |

Triploidy had never been emitted on a real array before this series, and it is called here with no
false triploidy on any of the 58 negatives.

The 45,X row read 0 of 20 through 5.30.1, for two reasons that were both fixable. The converter
mapped the file's numeric chromosome codes with a fixed table and wrote every chrX marker out as
chrY; the file states its own names and the converter now reads them. And the sex call read chrX
heterozygosity alone, which is the same on a 46,XY male and a 45,X female, so one X always read
male. chrX is now measured for copy number and chrY is consulted, against the array's own chrX
rather than an absolute level, because both sex chromosomes are single-copy in a male and an
absolute cut taken from another panel called four normal males a monosomy.

| measured on the 189 stated-karyotype arrays | log2 |
|---|---|
| chrX, one copy (20 stated 45,X and 26 normal males) | -0.49 to -0.63 |
| chrX, two copies (33 normal females) | -0.05 to +0.01 |
| chrY minus chrX, normal male | +0.09 |
| chrY minus chrX, stated 45,X | -1.78 |
| chrY minus chrX, normal female | -2.25 |

No array of the other 144 is called a monosomy. The three stated 45,X not called read two X copies
with female heterozygosity, 0.79 to 0.85 of their own autosomal rate, so those three arrays are not
reading one X whatever the pregnancy was karyotyped as.

### 43 arrays stated to carry a pathogenic CNV

The same series states "Chromosomal deletions/duplications" on 43 further arrays and says neither
which chromosome nor how big, so the only scoreable question is whether the tool reports an
imbalance at all. It does on 8 of 43, 18.6% [9.7%, 32.6%], and only ONE of those eight is a
segment. The other seven are five whole-chromosome gains and two sex chromosome constitutions,
which may be the stated event or may be a different event on the same array; the label cannot say.

That single segment is the honest measure of this class on bulk material. Clinical CNVs of this
kind run 0.6 to 2.6 Mb where a companion series states their size, and the smallest segment this
panel can report is 26 Mb. The arm measures the floor, not the detector.

Through 5.30.1 the gains row read 1 of 86. The magnitude gate demanded 0.40 log2 before calling a
gain, a figure measured from four arrays of one biopsy disagreeing with each other by 0.33, which is
an amplification artefact. On this chemistry a stated trisomy shifts by 0.299 to 0.357, the spread
between one array's own chromosomes is 0.0155, and the largest drift on any autosome of the 60
negatives is 0.070. The floor now follows the material.

### A segment has to be big enough to be seen

The segment scanner reports nothing under 2,400 markers, which is a different number of megabases on
every panel:

| panel | marker spacing | smallest reportable segment |
|---|---|---|
| Axiom GPL28377 | 3.5 kb | 8.3 Mb |
| Illumina CytoSNP-12 GPL13829 | 9.9 kb | 23.9 Mb |
| Affymetrix 250K Nsp GPL3718 | 10.9 kb | 26.1 Mb |
| CytoScan 750K GPL18637, genotyped markers | 14.3 kb | 34.3 Mb |

GSE163799 states a microduplication per sample, 0.6 to 2.6 Mb, which is three to forty times under
that floor: **0 of 15 are reported**, and none could be. Segmental duplication therefore remains a
class this tool has never emitted on real material. The one absence-of-heterozygosity case in that
series, stated on a whole chromosome, is reported.

## The detection limit

Sensitivity was never quantified before. It is now, on the only material with a per-sample answer
for both classes: chromosome 21 on karyotype-confirmed trisomy 21 arrays as positives, every
autosome of a euploid reference line or adult donor as negatives.

Both gates swept together, because a call must clear both and quoting either alone hides the other.
Each cell is sensitivity, then false calls per 1,000 known-negative chromosomes.

| z \ floor | 0.10 | 0.20 | 0.30 | **0.40** | 0.50 |
|---|---|---|---|---|---|
| 2 | 42% / 135.1 | 42% / 99.7 | 42% / 80.8 | 40% / 68.8 | 25% / 61.9 |
| 3 | 3% / 48.0 | 3% / 44.2 | 3% / 35.4 | 3% / 27.8 | 2% / 25.3 |
| 4 | 0% / 22.7 | 0% / 22.1 | 0% / 18.9 | 0% / 15.8 | 0% / 14.5 |
| **5.49** | 0% / 2.5 | 0% / 2.5 | 0% / 1.9 | **0% / 0.63** | 0% / 0.6 |

**At the shipped operating point on this material: sensitivity 0 of 60, false calls 0.63 per
1,000.**

That is the weaker of the series' two chemistries. Over every karyotype-confirmed trisomy 21 array
rather than the 60 the curve draws on:

| platform | chr21 called | false calls on euploid material |
|---|---|---|
| GPL6985 | 16 of 82, 19.5% | 0 of 198 |
| GPL8855 | 11 of 232, 4.7% | 5 of 2,134 |

The curve has to be drawn per platform; pooled, it hides a four-fold difference.

The signal is present and the statistic cannot reach it. Missed positives carry a median shift of
0.518, against the 0.585 a trisomy produces in theory, at a median z of 1.9 against a threshold of
5.49. **58 of the 60 missed positives carry a real shift**; only 2 have no signal at all. Lowering z
to 2 buys 40% sensitivity at a 6.9% false-call rate, which is not an operating point anyone should
use.

This is a property of the material and the design together. Within one array, a chromosome is
measured against that array's own 21 others, and on amplified single cells that spread is wide
enough to swallow a real trisomy. Between arrays the same data separates almost perfectly. A
per-array instrument cannot reach a cohort statistic, and closing this needs a reference cohort
rather than a different constant.

**The curve covers amplified single cells only.** The corpus holds no known-positive bulk arrays:
the 12 bulk donors are negatives, and the mole trisomies are 3 of 13 without the submitters saying
which three. No claim is made here about bulk sensitivity.

## Uncertainty

Wilson 95% intervals on the rates above. At 0 or 1 the interval is what bounds the claim. Where
observations share an array they are not independent draws, so both ends of what the interval could
be are given: the independent end assumes no sharing, the by-array end assumes an array moves as
one, and the truth lies between them.

| claim | k/n | independent | by array |
|---|---|---|---|
| trio-resolved attributions correct | 8/8 | 67.6% to 100% | 4/4 arrays, 51.0% to 100% |
| all parental calls correct | 263/263 | 98.6% to 100% | cluster count not recorded |
| false positives, adult donors | 0/242 | 0% to 1.6% | 0/11 arrays, 0% to 25.9% |
| false calls, euploid lines and donors | 1/1584 | 0% to 0.4% | 1/72 arrays, 0.2% to 7.5% |
| heterozygosity-loss findings, one-complement genomes | 0/144 | 0% to 2.6% | one per array |
| sensitivity, single cells, shipped point | 0/60 | 0% to 6.0% | one per array |
| reconstruction verdicts correct | 35/41 | 71.6% to 93.1% | one per array |
| whole-chromosome gains, bulk | 83/86 | 90.2% to 98.8% | one per array |
| triploidy, bulk | 20/23 | 67.9% to 95.5% | one per array |
| bulk negatives clean | 58/58 | 93.8% to 100% | one per array |

Two arms that were marked clustered are not: the mirrored-loss arms are twelve arrays carrying one
constructed event each, so their intervals stand as printed. Narrowing the genuinely clustered rows
needs a cluster bootstrap over arrays, and that needs each array's own numerator and denominator.
The harnesses emit the pooled figure today, which is the open item rather than a correction factor.

## Stated limits

- **The panel builder's selection is validated on one father with two arrays.** A second
  independent father is needed before the ratio generalises beyond him.
- **Eight trio-resolved breakages are consistent with a true accuracy as low as 67.6%.** The claim
  rests on the twelve-family arm, not on those eight.
- **The zygosity boundary is validated on one chemistry.** `HET_BAND_EXCESS` subtracts the mid-band
  mass at homozygous calls, which exists only where genotype and allele frequency were derived by
  separate algorithms. Cluster-file and GenomeStudio exports are structurally not that, confirmed on
  three file formats.
- **The relationship verdict is withheld on amplified material.** It read `unrelated` for a
  confirmed father 9 times out of 9 and for a stranger 9 times out of 9. Identical words for
  opposite facts are worse than silence. Two independent measurements disagree about whether a
  per-material threshold exists, so none is set.
- **Direction between two arrays is not recoverable.** Under Hardy-Weinberg the parent-child joint
  distribution is exchangeable. A third array does not rescue it: measured over 8 complete trios,
  both the opposite-homozygote rate and the window floor put the inverted configuration inside the
  correct range.
- **Sensitivity for whole-chromosome gains on single cells is low and platform-dependent**: 19.5%
  on GPL6985 and 4.7% on GPL8855 across every confirmed trisomy 21 array, at a false-call rate that
  stays near zero.
- **Mixed samples are not detected.** Two signals respond on 3 of 4 real two-person mixtures, and
  neither identifies a mixture.
- **A file without intensity shows no gains,** and a loss there cannot be told from isodisomy. The
  run states both.
- **Whole-chromosome sensitivity cannot be scored on the karyotyped GPL3718 single cells.** All 47
  listed monosomies sit under the call-rate floor.
- **A sex chromosome constitution is reported only on material expected to be diploid**, and only
  where chrY can be measured. Where a panel carries no chrY the tool reports one X and states that
  a male karyotype and a monosomy X are not separable on that file. 47,XXX and 47,XXY are emitted by
  the same rule and have no stated case in any corpus here, so neither has been measured.
- **A segment under 2,400 markers is not reported**, which is 8.3 Mb on the densest panel here and
  34.3 Mb on the sparsest. Clinical microduplications are 0.6 to 2.6 Mb and cannot be seen at all.
- **Segmental duplication has still never been emitted on real material**, and neither has a
  segmental change introduced by the gamete.

## How the evidence is produced

Every number above comes from a harness in `audit/`, each run against files on disk rather than
against fixtures. The harnesses that carry a scored answer are `detection-truth`,
`attribution-truth`, `detection-limit`, `handoff-truth`, `mole-truth`, `blastomere-truth`,
`replicates`, `zygote-origin`, `parental-power`, `progenitor`, `meiotic-segregation`,
`duplicate-truth`, `pgd-families`, `karyotype-truth`, `cellline-truth`, `mole-genotype`,
`upd-premise`, `coverage`, `cytoscan-truth`, `gse21732-truth`, and the two on the second public
series.

Fifteen harness defects were found and corrected during this work, most of them the same mistake:
a negative set built from material that was not actually negative. Two of them produced clean-looking
results that were reported before being caught. Of the last two, one scored a refusal as a miss,
which read 47 refused monosomies as 47 missed ones, and the other took a cell line's sex from a
series field that contradicts itself rather than from Cellosaurus. Two more were in a format
converter rather than a harness. The reader for Affymetrix CytoScan files emitted empty marker
names, so every row inherited one position; and it mapped the numeric chromosome codes with a fixed
table, writing every chrX marker out as chrY on 214 arrays. The first was caught on the first file,
because a converter is checked against the file it came from before anything is scored. The second
was not, and it stood until chromosome X was measured for its own sake: a mislabelled chromosome
reads as a plausible result on every array that does not depend on it, and no autosomal number
moved. Both are why the converter now takes the chromosome names from the file's own header rather
than from an assumption about the format. The rule that now governs every arm is that a
dataset enters the corpus only if its answer is retrievable per sample, from the source, before any
analysis runs. A 5,062-sample series was rejected under it, because its answers live in a paper
rather than in the data.

## What no public dataset here supplies

Four things the tool would be measured better against, none of which exist in the corpora above.

- **Embryo arrays with both parents.** Of the three public embryo series, only GSE148488 genotyped
  the parents, so parent of origin is answerable on 69 arrays of 355. Every other embryo array here
  can be read for copy number and not for whose copy it is.
- **A sample sheet for the lab's own 884 arrays**, stating stage, embryo and parents. Without it
  those arrays can be scored one at a time and cannot be grouped: mechanism needs two units known
  to be one embryo, and the array cannot say which embryo it came from.
- **A stated 47,XXX or 47,XXY.** Both are emitted by the sex chromosome rule and neither has a
  stated case in any corpus here, so neither has been measured.
- **Bulk arrays with stated segmental imbalances and coordinates, at a size this panel can see.**
  The stated CNVs available are 0.6 to 2.6 Mb against a 26 Mb floor, and the 43 arrays whose CNV is
  stated without a size measure presence only.
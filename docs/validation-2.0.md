# 2.0 validation against public data

The 2.0 inference path decides whether an embryo carries the paternal allele, and if not,
whether the absence is copy-reducing. This file records what that path has been tested
against, using published data with an answer established independently of this tool.

Everything here is reproducible from public files. No patient data is involved and none is
required to repeat it.

> **Research use only.** Not a clinical diagnostic.

## Why these controls

Until this pass the positive path had only ever run against markers written by hand. Every
run on real data had ended in a refusal, correctly, so nothing confirmed the tool could
detect a real event rather than merely decline to.

The controls need a documented answer that does not come from this tool. Two sources give
one:

- **Biology.** A male child's X chromosome is entirely maternal. His father contributes no X
  to him at all, so across the whole non-pseudoautosomal X the paternal allele is absent by
  construction, at one copy. No paper is needed to know the answer.
- **An orthogonal channel.** The run statistic reads genotypes only and never touches
  intensity. So a deletion found by the run statistic can be confirmed, or refuted, by the
  Log R Ratio at the same markers, which the statistic did not see.

## Data

The PennCNV example trio, `WGLab/PennCNV/example/{father,mother,offspring}.txt`. Illumina
GenomeStudio export with GType, Log R Ratio and B Allele Freq; bulk DNA; chromosomes 3, 11,
20 and X; 93,129 markers per sample.

Confirmed from the data before use: father male, mother female, **offspring male**
(chrX/autosomal heterozygosity 0.000, 0.964, 0.000).

```bash
for f in father.txt mother.txt offspring.txt; do curl -sLO "https://raw.githubusercontent.com/WGLab/PennCNV/master/example/$f"; done
```

## Results

Fitted dropout across the trio was 0.014 from 44,131 father-heterozygous by mother-homozygous
markers, so `r_min` is 3 on these chromosomes.

| Chromosome | Truth | Longest run / informative | p | State called |
|---|---|---|---|---|
| X | Paternal X absent, one copy (male child) | **4,324 / 4,324** | 0 | withheld, see below |
| 3 | De novo hemizygous deletion, 71 kb | 13 | 1.2e-21 | `PAT0_MAT1`, 24.8 sigma |
| 11 | Deletion inherited from the father, 8 kb | 3 | 4.0e-03 | `PAT0_MAT1`, 16.2 sigma |
| 20 | Nothing | 2 | 0.16 | none, correctly quiet |

Every informative marker on the X shows paternal absence: the run is not merely long, it is
the entire chromosome, which is the correct answer and the strongest form it could take.

The two autosomal deletions were not known in advance. They were found by the run statistic
and then confirmed independently by intensity:

| | chr3 span, 34 markers | chr11 span, 6 markers |
|---|---|---|
| Offspring mean LRR | **-0.643**, 0% heterozygous | **-0.826**, 0% heterozygous |
| Father | +0.018, 29.4% het (normal) | **-0.564**, 0% het (carries it) |
| Mother | -0.006, 52.9% het (normal) | -0.033, 66.7% het (normal) |

So chr3 is de novo in the child and chr11 was transmitted by the father, and in both cases
the tool's claim that the *paternal* copy is the missing one is correct. Neither parent's
data alone distinguishes these, and no single-sample CNV caller can: parental origin is the
part that needs the trio.

## Compression, measured two ways

The chr3 deletion is a one-copy autosomal segment, so it measures LRR compression directly.
Two independent estimates:

- offspring, chr3, n = 34: **c = 0.594**
- father, chr11, n = 6: **c = 0.578**

They agree to 3%. As a null, the father's LRR across the chr3 span gives c = -0.009, which is
the correct answer for someone who does not carry that deletion and is a check that the
method is measuring what it claims.

## Defects this found

Five, all fixed in this pass. Listed because the point of a control is what it catches.

1. **`score_paternal` ignored the mother.** Presence was scored whenever the embryo carried
   the father's allele. But presence is an identity claim: it only holds if the mother could
   not have supplied that allele. At a paternal deletion the embryo is hemizygous for the
   maternal allele, so wherever a heterozygous mother transmitted the same allele the father
   carries, the call still contains his allele and read as "present". This broke real runs at
   roughly half of all mother-heterozygous markers. The X chromosome came back as a longest
   run of 21 instead of 4,324.

2. **Markers that cannot say were counted as presence.** `None` was folded into the flag list
   as `False`, which both broke runs and inflated `n`. Excluded now.

3. **The reported segment was the longest anywhere in the window,** not the one at the
   variant. On a chromosome with more than one event that quotes a real sigma for an
   unrelated locus. Worse, it counted every marker carrying that label anywhere rather than a
   contiguous stretch, and separation grows as sqrt(n): a normal chromosome 3 reported 7,499
   markers at 359 sigma.

4. **The CLI defaulted `--product` to a named array.** The chrX baseline offset is a measured
   per-array constant and the calibrator correctly refuses an array it has no measurement for.
   Defaulting to a specific product defeated that refusal and applied the Axiom constant
   (0.1908) to Illumina data, putting c at 0.456 against a measured 0.59.

5. **The chrX probe offset was applied when deriving c but not when interpreting chrX.** It
   is most of the chrX-to-autosome shift, so without it a male's single X reads as two
   copies, and it did: `PAT0_MAT2`. Copy-number states on chrX are now refused unless the
   array's offset has been measured. Genotype statistics are untouched by this, which is why
   the X result above still stands.

## Where the run statistic fails, and the test that catches it

A second control was added to settle this: the HapMap CEU trio (NA12891 father, NA12892
mother, NA12878 offspring), bulk DNA, 3,908,761 markers across all 22 autosomes, from the
2010-08 phase II+III forward-strand release. Fitted dropout 0.0076. At opposite-homozygote
parents the child is non-heterozygous at 0.737% of markers, which lands independently beside
the Kothiyal floor of 0.63% the tool already uses.

Genome-wide it produces 35 runs reaching `r_min`, and a whole-genome-amplified blastomere
(GSE148488 A8, dropout 0.331) produces 824. Neither number means anything on its own: a normal
human carries hundreds of deletions, so dozens of genuine paternal-origin losses are expected,
and a per-locus p-value never claimed otherwise.

What separates them is a test needing no external data. Across a real paternal deletion the
child is hemizygous for the maternal allele, so every informative marker shows the paternal
allele absent and **none** shows the maternal allele absent. Amplification dropout has no such
asymmetry: it drops alleles at random, independently of which parent supplied them, so it
scatters violations in both directions through the same region. The maternal test is the same
scoring function with the parents exchanged, so it costs nothing to run.

| | runs at `r_min` | one-directional | both directions |
|---|---|---|---|
| HapMap CEU trio, bulk DNA | 35 | **34** | 1 |
| A8 blastomere, amplified | 824 | 1 | **823** |

Inside A8's runs the maternal-violation rate per informative marker is about 0.35, against a
fitted dropout of 0.331. Those regions are not deletions; they are dropout, measured as
dropout.

Power runs the right way round. On A8 there is a median of 16 maternal-informative markers
inside every run and none has zero, so the test is fully powered exactly where the artefacts
are. On the clean trio's small events the median is 1 to 2, so "34 of 35 clean" is consistent
with real deletions rather than proof of them, and a small event can leave the test with
nothing to say. It should therefore report insufficient evidence rather than a pass when too
few maternal-informative markers fall inside a span, and it is unavailable in the no-mother
degraded mode.

## A separate defect: runs ignore physical distance

Run length counted markers that are adjacent in the informative subsequence, with no constraint
on how far apart they are on the chromosome. Three violations spanning 5,580 kb on chromosome 7
counted as a run of 3 exactly like three violations inside 1 kb. A 5.6 Mb hemizygous deletion in
a healthy adult would be a major clinical finding, not a quiet data point.

Against the genome-wide informative density of one marker per 18.6 kb, 3 of the 35 clean-trio
runs span 30 to 150 times what their length implies. Runs are now split where consecutive
absent markers lie more than ten times the local informative spacing apart, on the grounds that
the intervening sequence was never interrogated.

## Both gates, measured

| | runs at `r_min`, ungated | surviving as paternal loss |
|---|---|---|
| A8 blastomere, amplified | 824 | **2** |
| HapMap chr3 / 7 / 11 / 15 / 19, bulk | 18 | 12 |

On A8 the chromosome 6 run of 41 markers at p = 6.9e-17, which sits in the MHC, is refuted with
27 of 33 maternal-informative markers in the same span also missing the maternal allele, a
likelihood ratio of 10^45 favouring dropout. The two survivors are at chr6:28.3 Mb and
chr10:6.7 Mb; neither is at the edited locus.

On the PennCNV trio all four controls behave as before: the chr3 and chr11 deletions survive
with `PAT0_MAT1` intact, chr20 stays quiet, and chrX still reports the whole paternal
chromosome absent.

That last one needed a fix. The compactness rule splits at centromeres and heterochromatin, so
a male's X fragments into pieces even though it is genuinely absent end to end, and the longest
run fell from 4,324 to 116. The whole-chromosome test was keyed to the longest run; it is now
keyed to whether every informative marker on the chromosome shows absence. Fragmentation
describes where the array put its markers, not what is missing from the sample.

Only a positive refutation blocks a result: "insufficient" leaves the run standing with the
caveat attached, because on clean DNA the test is powerless by construction and withholding
every clean-DNA result would discard the cases that work best.

## Parent of origin from sperm and sample alone

The question a lab asks first, on the pairing a lab actually produces: one sperm array and the
sample that came with it. Three measurements, no mother, no variant, no phase.

| axis | what it answers | paternal-only | unrelated | biparental |
|---|---|---|---|---|
| paternal alleles absent | is the paternal genome here | 0.16-0.22% | 5.1-5.6% | low |
| alleles the father lacks | is another genome here | 3.4-5.2% | 14.2% | 34.7% |
| heterozygous BAF band | is it diploid at all | 1.3-3.4% | - | 15-16% |

The second axis is the one that matters most, because paternal absence alone cannot separate a
paternal-only genome from a biparental one: the father's alleles are present in both. Alleles he
cannot supply had to come from somewhere else, so they measure the other parent directly.

The third is an independent check through intensity rather than genotypes, and it is not fooled
by dropout: a biparental embryo at 33% dropout still reads 22.2% in the heterozygous band,
because dropout smears that cluster rather than emptying it.

On a five-sample Egli-protocol case (one sperm, four samples):

| sample | absent | father lacks | BAF band | call | chrY |
|---|---|---|---|---|---|
| 02 | 0.16% | 3.39% | 1.27% | androgenetic, Y-bearing | present |
| 03 | 6.81% | 11.69% | 3.43% | gynogenetic | absent |
| 04 | 0.22% | 5.15% | 2.02% | androgenetic, X-bearing | absent |
| 05 | 0.17% | 4.07% | 1.57% | androgenetic, X-bearing | absent |

Three of the four have no chrY, so a chrY test sees 03, 04 and 05 as one indistinguishable
group. The SNP rates separate them by thirty to forty fold. Sperm type is inferred from the chrX
SNP rate, so chrY is a cross-check and never the method, and where the two disagree the sperm
type is reported as unsettled rather than guessed.

## Why the decision boundary is not a population constant

The first version compared the absence rate against a measured "unrelated genome" rate of 5.5%.
Measuring properly across 50+ pairs on this array showed that number cannot exist. It ranged
from **6.8% to 50.3%**, driven by two things: the ancestry the two people share (same-study pairs
sit near 7%, cross-study pairs at 8 to 15%) and how noisy the sample used as the father is (a
sample with a 33% no-call rate, used as the father, drives it to 48 to 50%).

Worse, the two populations overlap. Zuccaro's A8 embryo against **its own father** shows 9.69%
absence, which is higher than genuinely unrelated pairs at 5.1 to 5.6%. Any fixed threshold
between those numbers calls a real father-offspring pair unrelated.

The cause is mechanistic. Dropout only manufactures paternal absence by turning a heterozygous
call homozygous and discarding the paternal allele, so the inflation scales with heterozygosity
TIMES dropout, and neither alone. A homozygous genome is immune however much drops out, which is
why androgenotes hold at 0.16% despite a 13.6% no-call rate while a diploid embryo at 46.9%
reaches 9.69%.

So the boundary is now the sample's own noise, `no-call rate x heterozygous fraction`, which
needs no population reference and no ancestry assumption. It bounds every related pair measured,
and tightly:

| pair | no-call | het | bound | observed absence |
|---|---|---|---|---|
| PennCNV offspring | 0.4% | 31.8% | 0.13% | 0.05% |
| 52461 androgenote 02 | 13.6% | 1.3% | 0.18% | 0.16% |
| 52461 androgenote 04 | 13.7% | 2.0% | 0.27% | 0.22% |
| Zuccaro A8 | 46.9% | 22.2% | 10.4% | 9.69% |

Against nine pairs of known relationship, including the adversarial A8 case, the rule makes
**zero wrong calls**: three androgenotes, one gynogenote, three unrelated pairs and two bulk
trios all resolve correctly, and a rate landing between the bound and three times it is left
uncalled rather than guessed.

The 5.5% figure is retained for reporting context only and is explicitly not a boundary.

## The second axis is derived too, and its margin is narrow

A second parent makes the sample carry an allele the father lacks only where he is homozygous,
at the sum over markers of `p^2 q + q^2 p = pq`. His own heterozygosity is the sum of `2pq` over
those same markers. So the expected second-parent contribution is exactly **half the father's
heterozygosity**, measurable from his own array with no reference panel and no ancestry
assumption. The hardcoded 7% is gone.

| pair | composition | half the father's het | observed |
|---|---|---|---|
| 52461 androgenotes 02 / 04 / 05 | paternal-only | 8.48% | 3.39 / 5.15 / 4.07% |
| Zuccaro A8 | biparental | 8.33% | 34.66% |
| PennCNV offspring | biparental | 16.10% | 16.35% |
| HapMap NA12878 | biparental | 9.72% | 10.65% |

The prediction holds, but note what the biparental rows show: a clean biparental sample lands
almost exactly ON the prediction, because the prediction is its expectation. So this axis
separates by roughly 1.6x where paternal absence separates by 30x, and noise pushes a
paternal-only genome upward toward the boundary.

The classifier is therefore restructured so this axis carries the decision only where it must.
A homozygous genome has one allele per locus, and if those alleles are the father's there is no
room for a maternal complement, so zygosity settles it and the axis-2 residual is error by
construction. The axis is consulted only for a DIPLOID sample carrying the father's genome,
which is the genuinely ambiguous case: biparental, or paternal-only from two sperm. There the
number, its derived boundary and the margin are all printed, and a call within 40% of the
boundary is declared non-decisive.

Where a file carries no B-allele frequencies at all, zygosity falls back to genotype
heterozygosity against the same derived reference, flagged as the weaker measure.

Against nine pairs of known composition, three paternal-only, three biparental and three
unrelated, the classifier now makes **zero wrong calls with no hardcoded population constant
anywhere in it**.

## Segmental parent of origin

The genome-wide call answers whether a paternal genome arrived. This answers which parts of one
did, which is the question for a mosaic or partially-lost sample.

Resolution is derived rather than chosen. One informative marker per 5 kb on the Axiom array
makes 200 markers about 1 Mb, and 200 markers distinguish 0.16% from 6.8% at odds near 10^17. A
two-state chain runs at marker level with transitions scaled by physical distance, so boundaries
are not quantised to a window grid.

**The failure it was built through.** The first version estimated the two emission rates as the
10th and 90th percentiles of the observed window rates. That finds two populations even in a
uniform genome, because those are just the tails of one noise distribution: a sample whose
paternal genome is present everywhere at 0.16% came back as 28 "absent" segments running at 0.3
to 1.5%, and a uniformly gynogenetic sample was carved into 720. The rates are now anchored to
the calibrated noise bound, and segmentation is skipped entirely when no window exceeds it or
when nearly all do.

**The independent check.** On sample 02, an androgenote from a Y-bearing sperm, the one confident
segment is chrX:2,700,151-155,233,115. It begins just past 2,699,520, the PAR1 boundary. The tool
recovered the pseudoautosomal boundary from the data alone: PAR1 is present because the father
transmits it, and the non-PAR X is absent because this sample carries his Y instead. Nothing in
the code knows where PAR1 is.

**Two filters, because small spans are where correlated probe failure shows.** Every segment
carries its margin above the sample's own noise floor: the real chrX finding sits at 18x while
every spurious subtelomeric stretch sat at 1 to 2.5x. And across a plate, a stretch recurring in
sample after independent sample is a property of the probes there rather than of any genome. On
this plate every marginal call recurred across all three paternal samples while chrX appeared in
exactly the one carrying the Y.

Known limitation: small subtelomeric and pericentromeric spans remain the false-positive mode.
They are reported with their margin and recurrence rather than suppressed by a region blacklist,
which would need validating in its own right.

## Input normalisation

One reader, validated against real public files rather than against a specification.

| format | source | samples/file | coding | BAF | LRR | build |
|---|---|---|---|---|---|---|
| Axiom `.CEL.probes`, tab | lab 52461 | 1 | numeric | yes | yes | GRCh37 |
| Axiom `.CEL.txt`, comma | GSE148488 | 1 | numeric | yes | yes | GRCh37 |
| FinalReport wide | PennCNV example | 1 | AB | yes | yes | undetermined |
| FinalReport long, `[Header]` block | GSE16912 | many | nucleotide | yes | yes | undetermined |
| HapMap raw, space-delimited | HapMap r28 | **174** | nucleotide | no | no | undetermined |
| WGS VCF, GT/AD/DP | GIAB HG001 | 1 | numeric | from AD | from DP | GRCh38 |

Four things vary independently and each had to be handled rather than assumed.

**Where the header is.** GenomeStudio writes a `[Header]` metadata block before `[Data]`, putting
the columns on line ten. Rather than special-casing that marker, the reader takes the first line
that RESOLVES to the fields it needs, which covers that case, a comment preamble and a plain
header without knowing which it is looking at.

**How samples are laid out.** One column per sample suffixed `<sample>.GType`; one column per
sample named for nothing but the sample, as HapMap writes it; or one ROW per sample and SNP with
a `Sample ID` column. All three appear in the files above. The middle one yields 174 samples from
a single chr22 file.

**How the call is spelled.** A single column, or two allele columns (`Allele1 - Forward` and its
pair), coded numerically, in AB space, or as nucleotides.

**Nucleotide harmonisation, which cannot be done per file.** A marker whose alleles are A and G
shows only "AA" in a sample homozygous for the first, so that file alone cannot know G exists,
and two files mapping independently would disagree about which allele is A. The alleles seen at
each marker are pooled across every sample being compared and one mapping is fixed for all of
them. Which allele becomes A never matters, only that it is the same everywhere, and that makes
the step a no-op on data already in AB space.

Validated end to end on the CEU trio read straight out of the raw HapMap release with no
converter: father to child 0.009% paternal absence against 0.78% explainable, father to mother
4.74%.

**A bug this found.** A wide export with more than one `.GType` column read only the first sample
and dropped the rest silently. The test written for the harmonisation no-op is what exposed it.

## PLINK

Binary `.bed`/`.bim`/`.fam` and text `.ped`/`.map`. The binary layout is the documented one: two
magic bytes, a mode byte that must be SNP-major, then two bits per sample per marker, least
significant pair first, with **01 as the missing code rather than a genotype**. Reading that pair
as a call would invent a genotype for every sample that failed, so it is tested explicitly.

`.bim` names A1 and A2 for every marker, so the AB mapping is exact here rather than pooled
across samples. `.ped` has no such column and goes through the same harmonisation as any other
nucleotide input. Individual-major files and bad magic bytes are refused rather than misread.

Caveat worth stating: the encoding was implemented from the specification and round-tripped, not
checked against a file written by PLINK itself, because no public `.bed` was reachable during
development. The bit-level decode is the part that would benefit from one real file.

## Both parents

With only the father, whether a second parent contributed has to be read off the rate of alleles
he lacks, an axis that separates by about 1.6x and needs his own heterozygosity as its reference.
With the mother present, each parent is tested the same way on the axis that separates by
thirty-fold, and the class falls out of the pair rather than being inferred.

| paternal | maternal | class |
|---|---|---|
| present | present | biparental |
| present | absent | androgenetic |
| absent | present | gynogenetic |
| absent | absent | neither parent: check the pairing before the biology |

Validated on every real trio available. PennCNV and HapMap CEU both classify as biparental;
substituting an unrelated individual as the father correctly yields gynogenetic. Zuccaro A8
returns **unclear**, which is right: at a 46.9% no-call rate its maternal side lands inside the
ambiguous band, and the tool declines rather than guessing.

Two guards. A maternal X reported absent is a real loss rather than sex determination, since a
mother transmits an X to a child of either sex, so the paternal-side exemption is not applied.
And two declared parents whose genotypes agree at over 99% are flagged as a duplicate or
relabelled file, because both tests then pass against the same person and produce a confident
"biparental" that means nothing.

## Read-level structural evidence

Two mechanisms leave no trace on an array. A large deletion at the edited site can fall between
markers, and inside it a hemizygous call is indistinguishable from a homozygous one. An insertion
or balanced rearrangement changes no copy number at all, so intensity is flat and genotypes are
unremarkable.

Both come from a **structural-variant VCF** rather than from a BAM, and the distinction is
deliberate: an SV VCF is the output of a caller that already did the read analysis and carries
the split-read and paired-read counts behind each call. Parsing BAM would mean reimplementing
that caller with no BAM on hand to check the result against.

Validated on the GIAB HG002 Tier 1 benchmark: 3,192 calls parsed, 567 passing and carried, 218
deletions with a median span of 167 bp and a largest of 32,197 bp.

Supplying it moves `H3b` and `H3f` off "not tested" in the mechanism ledger. Three guards: a call
present in the file but genotyped `0/0` is not carried by the sample; a call that failed the
caller's filters does not count; and **absence of a call is never reported as absence of an
event**, because a caller run with a size threshold above the event, or without the type enabled,
produces exactly the same empty result as a genome with nothing there.

## Assembly matching

The bundled deCODE genetic maps are GRCh38, and the scaffold's confidence comes from
centimorgan distances read off them. Coordinates from another assembly index the map at the
wrong genomic locations and return a plausible number from the wrong locus, which is the worst
kind of wrong: nothing looks broken.

The build is determined from marker POSITIONS alone. A marker cannot lie inside an assembly
N-gap or past the end of a chromosome, so counting illegal placements under each candidate
separates them, and the correct build gives zero. No rsIDs, no manifest, no chain file.

| dataset | verdict | illegal under GRCh37 | under GRCh38 |
|---|---|---|---|
| lab 52461 Axiom, 825,656 markers | **GRCh37** | 0 | 6,315 |
| GIAB HG001 WGS VCF, 89,974 | **GRCh38** | 761 | 0 |
| PennCNV Illumina, 93,129 | undetermined | 405 | 280 |
| HapMap chr11, 207,012 | undetermined | 620 | 146 |

The two undetermined ones are hg18, which these tables do not carry. They are reported as
undetermined rather than forced to whichever is closer, because the closer one would be a silent
error.

**This mattered.** The lab's own Axiom exports are GRCh37 while the maps are GRCh38, so the
scaffold had been reading the map at wrong coordinates for exactly the data it was built for. On
a build mismatch the map is now not consulted at all: distances fall back to a uniform 1 cM/Mb,
which ignores hotspots, and the refusal says so. The homologue calls did not change; the stated
confidence did, from 10^8 to 10^6, which is the honest direction when you know less than you
thought.

There is no liftOver here by design. UCSC gap and chromosome-length tables are free for any use;
UCSC CHAIN files are free for non-commercial use only, and a field-of-use restriction cannot be
sublicensed under Apache 2.0 however little the file costs. So the tool can tell you the build
and cannot convert between them.

The Python and TypeScript tables are transcoded from one generated source, and a test asserts
they are identical so a regeneration of one cannot silently diverge from the other.

## Sequencing input

A VCF is read directly, and where the file carries `AD` and `DP` it supplies two things an array
cannot.

`AD` gives the reads behind each allele, so the alternate fraction is a measurement rather than a
three-way call. A genotype rounds 8% alternate reads to homozygous reference and the evidence of a
minority lineage is gone; the fraction keeps it.

`DP` is linear in copy number, so `log2(DP / median DP)` is an LRR with a compression of exactly
one. The entire per-array calibration that intensity requires, the chrX probe offset, the
compression constant, the credible bounds and the noise ceiling, simply does not arise. That is
the single biggest advantage of sequencing input here.

Checked against the GIAB HG001 benchmark VCF: 89,974 biallelic records, alternate fraction 0.495
at heterozygous calls where 0.5 is expected, and a depth-derived LRR at sigma 0.232, tighter than
the 0.30 to 0.43 the array gives on bulk DNA.

**One trap, found in that same file.** It contains **zero** homozygous-reference records out of
89,974, because a single-sample variant-only VCF lists only sites where the sample differs from
the reference. Those are exactly the markers a paternal-absence test needs. The analysis would
still run, silently, on a marker set conditioned on which sample happened to vary. Any VCF whose
homozygous-reference fraction is under 1% is now flagged, with the fix stated: a joint call across
the family, or a gVCF.

## The per-locus deletion test is the one thing that needs the mother

Worth separating clearly, because the scope is narrow. Parent of origin, sperm type, segmental
loss and the homologue scaffold are all decided from RATES across hundreds of thousands of
markers, or from co-inheritance between siblings, and none of them needs maternal genotypes.

What needs her is the per-locus test below, because attributing a single allele requires someone
who could not have supplied it. So the oocyte donor's array is a strong recommendation rather
than a requirement: it measures dropout directly, restores that test, and raises confidence
across every borderline call elsewhere.

## Father plus embryo cannot support the per-locus deletion statistic

The section above and this one are about different questions, and the difference is the whole
point. A RATE across half a million markers has its own null and needs no mother. A claim about
one locus is an attribution, and attribution needs someone who could not have supplied the
allele.

For the per-locus question, testing the degraded mode found it structurally broken rather than
merely weaker, and it is now refused rather than answered.

Three failures compound. Paternal PRESENCE cannot be established without someone who could not
have supplied the allele, so with no mother the informative set contains **absences only**:
nothing can ever break a run, and run length becomes a count of Mendelian violations wearing a
p-value. Dropout cannot be fitted either, since the estimator needs father-heterozygous by
mother-homozygous markers, so `q` falls back to the 0.63% error floor on a sample whose real
dropout is 33%, and `r_min` collapses with it. And the maternal cross-check above, the one
thing that catches dropout, is unavailable by definition.

Measured on the same embryo:

| A8 blastomere, chromosome 6 | with mother | without mother |
|---|---|---|
| fitted dropout | 0.331 | unmeasurable, falls back to 0.0063 |
| `r_min` | 12 | 3 |
| longest run | 41 | **610** |
| verdict | refuted as dropout, 10^45 | run of 610, p = 0 |

The run gets *longer* without the mother, because nothing can break it. On the PennCNV trio's
normal chromosome 20 the same mode produced a run of 3 spanning 35 Mb at p = 2.5e-07.

A mother-free dropout proxy was tried and does not work: comparing embryo to father
heterozygosity gives -0.917 on A8 where the truth is 0.331.

So no significance is claimed without maternal genotypes. Violation counts are still reported,
because a cluster of them is worth a human look, but father plus embryo can raise a suspicion
and cannot settle one. On single-cell input a deletion and amplification dropout are the same
picture.

## Which paternal homologue: the sibling scaffold

The other half of the question, "was the disease allele inherited and corrected or was the
wild-type inherited", is about WHICH of the father's two chromosomes came through. Nothing in a
genotype file says which one carries the mutation.

Several embryos from the same father break the deadlock. Each received exactly one paternal
homologue, so embryos that received the same one agree wherever the father is heterozygous and
embryos that received different ones disagree there. That splits them into two groups without
naming either. One embryo whose variant site was typed outside the array then names both groups
at once, and the label propagates to every sibling.

Note this reads the OPPOSITE marker class to the presence test. Presence needs the father
homozygous; homologue identity needs him heterozygous. A panel informative for one can be
useless for the other.

Measured on the real CEU trio, chromosome 11, 207,012 shared markers:

| | |
|---|---|
| father heterozygous | 39,807 (19.2%) |
| paternal allele deducible | 31,587 (79.4% of those) |
| deducible without the mother | 20,766 (52.2%) |
| usable within 1 cM of a locus | **559 markers** |

So marker supply is not the constraint. Recombination is.

### The flaw this found

The first implementation summed per-marker log-odds and reported roughly 10^50 confidence.
That is marker count masquerading as evidence. Markers on the same side of a crossover all
report one fact, not many, so within 1 cM there are one or two independent observations rather
than 559, and no amount of density reduces the real error source, which is whether a crossover
fell between the markers and the variant.

It is now a two-state chain along the chromosome, with transitions from the deCODE map and
forward-backward giving the posterior at the variant, which is inserted as an unobserved node.
Confidence saturates instead of accumulating: at a 0.01 error rate a 0.1 cM window and a 1 cM
window both report 10^11.3 despite ten times the markers.

### Crossover behaviour

Tested with real markers, real parental genotypes and the real map, transmission simulated so
the true homologue and the true crossover positions are known. Real siblings could not give
that, since nobody knows where a real meiosis recombined.

| embryos | deduction error | correct |
|---|---|---|
| 2, 4, 8 | 0.01, 0.10, 0.33 | **100%** |

The summed-odds version got one call wrong at 8 embryos and 0.33 error, and reported no
uncertainty on any call in any condition. The chain version gets that case right.

Caveat on this test: HapMap r28 is on b36 and the bundled map is GRCh38, so the cM readings sit
at the wrong coordinates. It cancels here because the same map drives both the simulation and
the inference, keeping it a valid test of the algorithm. Real use needs build-matched positions.

### What it requires

Two or more embryo lines from the same father, and one of them typed at the variant site
outside the array. Without the anchor the partition is still produced and still useful, but the
groups stay unnamed and H1 versus H2 remains unseparated.

## Reproducing

```bash
python -m originmarker.origin --father father.txt --mother mother.txt --embryo offspring.txt --chrom X --pos 75000000 --compression 0.594
```

Swap `--chrom X` for `3`, `11` or `20`, centring `--pos` on 4,035,000 and 81,185,000 for the
two deletions. Omit `--compression` to see the calibrator refuse an array whose chrX offset
has never been measured.

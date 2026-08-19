# Consult: parental balance, and whether any of this fits a lab

Two questions, and the second matters as much as the first. Sample identifiers are replaced with
letters throughout; every number below is measured, not illustrative.

## What the tool is

OriginMarker reads SNP-array exports from human embryo material and answers, per chromosomal
abnormality: where it is, and which parent's genome it came from. It is research-use-only. The
material is mostly amplified single cells (single blastomere, trophectoderm biopsy, single
embryonic stem cell), plus bulk DNA where available. Typical array: 825,656 markers, GRCh37.

Runs load one parent's array (usually the sperm donor) plus N sample arrays. A second parent is
sometimes available and is worth roughly threefold in detectable mosaic fraction.

The driving question, from the lab: **"if breakpoints and losses affect paternal or maternal
genomes differentially, or equally."**

### Channels that can name a parent

1. **Dosage.** Self-referenced allele-dosage shift at markers where the loaded parent is
   homozygous, compared against the rest of that same array's genome. Class-marginal posterior:
   copy-number class and mosaic fraction are both marginalised, because a gain inverts the sign map
   that loss and copy-neutral LOH share. Calibrated by isotonic regression on our own injections,
   leave-one-array-out. Four bands: A >= 0.985, B >= 0.90, C >= 0.75, D < 0.75. Measured accuracy
   across bulk / ES-single / trophectoderm / blastomere: band A 0.9952 to 0.9980, band D 0.6038 to
   0.6360.
2. **Obligate heterozygosity.** Mendelian, not self-referenced.
3. **Untransmitted.** Markers where the loaded parent is heterozygous and the sample reads
   homozygous.
4. **Zygosity (new, and question 1 below).** Described under "The inference we just shipped".

All reported confidences are capped by a systematic error bound of 0.206, derived from 13
independent validation units with zero observed events.

### The detection floors that caused the problem

Minimum mosaic fraction detectable at 1% false positives, one genotyped parent:

| Material | loss (whole chr / segment) | CNN-LOH (whole / segment) | gain (whole / segment) |
|---|---|---|---|
| bulk | 0.040 / 0.044 | 0.040 / 0.040 | 0.044 / 0.044 |
| ES-single | 0.348 / **none** | 0.399 / **none** | **none** / **none** |
| trophectoderm | 0.186 / 0.327 | 0.135 / 0.186 | 0.511 / **none** |
| blastomere | 0.628 / **none** | 0.399 / **none** | 0.620 / **none** |

"None" means no floor exists at any fraction up to 0.70. With two parents ES-single segments are
still none. **Sub-chromosomal intervals on amplified material are unreachable by dosage, at any
quality, with either parent count.**

Related: the fraction of detected events whose ORIGIN resolves while the CLASS does not is 0.296 on
bulk, 0.759 on ES-single, and 1.000 on both trophectoderm and blastomere.

## The run that prompted this

18 samples against one sperm-donor array. Genome-level calls: 8 gynogenetic (maternal genome only),
7 androgenetic (paternal only), 3 unclear. The taxonomy detected **1,272 regions**. Parental origins
called: **zero**. Every row read "not evaluable", correctly, because every region was
sub-chromosomal on ES-single material.

Sample A, gynogenetic: absence 9.04% against a 0.70% ceiling that already includes dropout and
error, i.e. **12.9x**, on 580,133 informative markers. Sample B, unclear: 11.63% against a 22.43%
ceiling, i.e. 0.5x, on 282,055 informative markers. Samples with detected change: one gynogenetic
at 5.0x (chr13 loss, chr19 loss), one gynogenetic at 10.1x (chr4 loss, 42 Mb), one unclear at 0.5x
(five gains).

## The inference we just shipped, and question 1

A gynogenetic embryo carries maternal DNA and no paternal complement. So a chromosome lost from it
was maternal, because there was no other chromosome there to lose. Androgenetic is the mirror. This
needs no detection floor. We now emit a parental origin for every change in a confidently
uniparental sample, from the genome-level call rather than from a per-event measurement.

Guards: consulted only where dosage returned nothing, so a measured answer always wins; silent on
biparental samples, unclear ones, unresolved zygosity, and genome-level calls that did not clear
their own ceiling.

Confidence: `reportedConfidence(log(foldOverCeiling), bound = 0.206, halfEvidence = log(3.0))`,
floor 0.5, cap 1 - 0.206 = 0.794. The half-evidence is anchored to `ABSENCE_MARGIN`, the threshold
at which a genome is declared uniparental at all, so a sample sitting exactly there lands halfway
from chance to the cap. 1.5x over ceiling gives 0.579, 3x gives 0.647, 12.9x gives 0.706. The
channel cannot reach band B by construction, because it has no injection series of its own.

**Q1a.** Is the inference sound as stated? Specifically: is there a mechanism by which a change
detected in a confidently uniparental genome could belong to the absent parent? Residual paternal
fragments below the absence ceiling, chimaerism, and contamination are the ones we thought of.

**Q1b. The band scale makes a differently-validated channel look uniformly weak, and we think this
is the sharpest question here.** Bands are A >= 0.985, B >= 0.90, C >= 0.75, D >= 0.55. They were
set for the DOSAGE channel, which has an injection series behind it and a measured accuracy of
0.9972 in band A. Any channel capped at 1 - 0.206 = 0.794 by this project's 13-unit validation can
therefore only ever occupy band D, with the top 5.5% of its range in C, and bands A and B
mathematically unreachable however decisive its own evidence is.

Concretely: a genome called gynogenetic with 9.04% absence against a 0.70% ceiling on 580,133
informative markers, which is 12.9x and not a marginal call by any reading, reports 0.706 and prints
as "weak, not for reporting". That is the number the lab sees on every row.

So: is the cap the right instrument, applied to the right channel? The 0.206 bound is a zero-event
bound on the ONE-PARENT MENDELIAN ORIGIN call, from 251 of 251 correct across 13 independent units.
The uniparental inference's failure modes are different ones: the genome-level call being wrong,
contamination, chimaerism. None of those is what those 13 units measured. We reached for the
nearest available bound rather than a bound on this claim.

Three ways out that we can see, and we do not know which is right. Bound this channel by something
that actually describes it, if such a bound can be derived. Report it on its own scale rather than
forcing it onto bands calibrated for a channel with thousands of validation points. Or accept that
it belongs in band D and that the lab's answer really is "weak" until the channel has a validation
series of its own.

Related: is a bounded, monotone, explicitly-uncalibrated number the right thing to publish at all,
or should this channel emit a verdict with no number? The lab's requirement was that every
abnormality carry a confidence-scored parental call, which pushed us toward a number.

We have already fixed one clear error found while asking this: the confidence floor was 1/3, the
chance level for the three-hypothesis channel the helper was written for, applied to a two-way
maternal-or-paternal question where chance is 0.5. It reported ignorance as worse than a coin flip
and dragged every number above it down. Corrected, 12.9x moved from 0.665 to 0.706. Still band D.

**Q1c.** Every change in one sample gets the same answer and the same confidence, because they all
rest on one genome-level call. We say so in the row text. Is that sufficient to stop a reader
treating 1,272 rows as 1,272 independent confirmations?

## The aggregation, and question 2

`parentalBalance.ts` compares uniparental groups directly: maternal-only genomes against
paternal-only genomes. No per-event assignment, so no per-event detection bias. Counts normalised
per 100,000 informative markers. Null by shuffling which genome is whose, 10,000 permutations,
keeping each sample's event and marker counts together. Two-sided. Per class as well as overall,
Bonferroni over classes tested. Seeded, so a report regenerates identically.

Refusals built in:
- Median informative markers between groups may not differ by more than **1.5x**, else
  `not-comparable`.
- At least **3** genomes per group, else `underpowered`.
- Where group sizes admit no p below the corrected threshold, the class is marked `underpowered`
  rather than reading as no difference. 4 vs 4 genomes admit no p below 0.042, so a planted 7.5x
  effect cannot be flagged at all. 8 vs 7 resolves it.

**Q2a. The unit problem, which we think is the biggest threat.** One wholly homozygous genome
produces hundreds of "regions" from the sliding-window LOH detector: one sample yielded 183
copy-neutral findings. Those are one biological fact sliced by a window, not 183 events. Our rate
is therefore partly a measure of how finely the detector cut each genome. The permutation shuffles
whole samples, so the test's type-I error should be protected, but the effect size is not
interpretable. What is the right unit? Candidates we see: independent breakpoints rather than
regions; total affected megabases; a per-genome binary "any change"; segment count after merging to
maximal intervals. Or is the whole rate framing wrong here?

**Q2b. The denominator.** Informative markers, or assessable megabases, or called markers, or
something that accounts for call-rate collapse regions being both a detection target and a hole in
coverage?

**Q2c. The test.** Label permutation on 8 vs 7 genomes with several classes. Is that the best
available, or should this be a mixed-effects or negative-binomial model with sample as a random
effect? The counts are heavily overdispersed and within-genome correlated.

**Q2d. The thresholds.** `MAX_POWER_SKEW = 1.5` and `MIN_PER_GROUP = 3` are ours, not measured.
They decide whether a lab sees a result at all. What should they be, and on what basis?

**Q2e.** Is comparing gynogenetic against androgenetic embryos even the right comparison for the
lab's question? These are abnormal conceptuses by definition. A difference between them may not
generalise to breakage in biparental embryos, which is what the question is presumably about. Is
there a defensible way to use the biparental samples, given per-event calls on them are band C/D at
best on this material?

## Question 3: does any of this fit a lab

This is the part we have the least basis for and would value most.

The tool is a browser page. Files are read locally, nothing is transmitted or retained. Output is a
12-page PDF plus a genome viewer, plus an optional feature-comparison addon that scores called
regions against common fragile sites, genes over 500 kb, centromeres and telomeres, and ENCODE
Repli-seq late-replication valleys.

**Q3a.** Who actually runs something like this, and at what point? Embryologist, genetic counsellor,
lab director, bioinformatician? Our assumption has been a research group re-analysing arrays after
the fact, not anyone in a clinical decision path. Is that assumption right, and does it change what
the output should be?

**Q3b.** What decision, if any, does a parental-origin call feed? We can name a parent for a
chromosomal change. We do not know what a lab does with that.

**Q3c.** Bands C and D name a parent at 0.75 and below, with measured accuracy around 0.60 for D.
We render them dimmed and label them "not for reporting". Is showing a number at all defensible
there, or does any number get quoted downstream regardless of its label?

**Q3d.** The run log was recently found to be printing 319 identical warnings and re-rendering per
line. Fixed. But it raises the general question: what does an operator actually need to see while a
run proceeds, versus what belongs only in the report?

**Q3e.** How does this sit beside existing PGT-A and PGT-M workflows? Is parent-of-origin something
those already deliver, and if so by what method, and is this redundant?

**Q3f.** What would make a lab distrust the output on sight? We would rather hear it now.

## What we would find most useful back

1. A ruling on Q2a, since it decides whether the headline number means anything.
2. Numbers for Q2d, or a principled way to derive them.
3. Blunt answers on Q3. If the honest answer to Q3b is "nothing, this is a research curiosity",
   that is worth knowing before we build the reporting layer around it.

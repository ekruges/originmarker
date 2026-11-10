<!-- Laboratory sample identifiers redacted: this file is tracked in a public repository.
     See scripts/release-check.sh, which refused the first commit of an earlier consult. -->

# Consult: saying "this array is too noisy" rather than "this chromosome had few markers"

Written 8 November 2026 against version 5.47.0, commit 47e48d4. Smallest of three consults sent
together. The others cover the segmental gain statistic and the sister-cell prior.

## The request

Our collaborator in the Egli laboratory, 7 November, in her words: when they run SNP arrays some
come out with much higher background than others, and in those cases the program should say it is
not making a call because the array is too noisy or poor quality, rather than because there were
not enough informative markers. She asked for this as a distinct reason, not a severity label.

## 1. What already ships

5.47.0 splits a withheld call two ways, and neither is hers.

**Marker-count branch.** `HetTally.offered` is the number of markers the two parents make
informative on a chromosome, counted before the sample's own calls are examined, so it is an upper
bound on what any array could return there. Under the floor, no array could do better and the
refusal says the parents offered too little. This branch has never fired on any autosome of 459
arrays: 200 informative markers is low against the ten thousand or so a chromosome offers, so a
sample reaches the other branch first.

**Call-rate branch.** A chromosome genotyped at 17 percent on an array that reached 95 percent is a
failure confined to that chromosome. The same 17 percent on an array that reached 22 percent is the
array. The comparison is against the array's own rate rather than a fixed threshold, because array
chemistry varies between runs.

Measured over 459 arrays in seven corpora: 3,065 chromosomes withheld, every one with a stated
reason. 2,469 are the collapsed-call-rate refusal, of which 2,398 name the array and 71 name the
chromosome alone. 596 are sex chromosomes, which the opposite-homozygote channel does not cover.
42 had previously produced no row at all.

So the existing split answers "did the array fail here, or did this chromosome have nothing to
read". It does not answer "did the array produce calls that should not be trusted", which is what
she asked. A high-background array can hold a respectable call rate and still be unreliable, and
nothing in the shipped refusal path looks at dispersion.

## 2. The statistic already exists and nothing consumes it

`SampleProfile.hetBafSd` in `web/src/ingest.ts`: the spread of B-allele frequency at heterozygous
calls, computed in the same streaming pass as everything else. Its docstring records the
calibration, which is the part that makes this consult cheap:

- unamplified genomic DNA reads about 0.088;
- whole-genome amplified material reads 0.21 to 0.30;
- NaN where the file carries no B-allele frequencies.

Its sibling `hetBand`, the fraction of autosomal markers inside the heterozygous BAF band, IS
consumed: 15 to 16 percent for a diploid genome, 1.0 to 7.8 percent across confirmed haploid
meiotic products, above 30 percent for a failed array, and it drives the ploidy and zygosity calls.

`hetBafSd` is reported nowhere and gates nothing.

So we are not asking you to invent a metric. We are asking whether this one is the right one, and
what to do with it.

## 3. The three things we cannot settle

**(a) Does dispersion mean anything on this material.** Every sample in the corpus that matters is
whole-genome amplified single-cell DNA, where dropout has already inflated the BAF spread by a
factor of two to three over unamplified DNA. The published array-QC literature we know of (LRR
standard deviation, BAF drift, waviness) was built on bulk DNA where a threshold of about 0.3 on
LRR SD marks a poor sample. We do not know whether any of those thresholds transfer, or whether on
WGA single cells dispersion measures amplification rather than array quality, in which case
refusing on it would refuse the material rather than the run.

**Question 1: on WGA single-cell SNP arrays, which dispersion statistic separates a bad RUN from
heavily amplified but usable material? If the honest answer is that none does and the number should
be reported without a threshold, say so. We would rather print a number with no gate than invent a
gate.**

**(b) The level is wrong for the refusal we need.** `hetBafSd` is per SAMPLE. The refusal it would
feed is per CHROMOSOME. `ChromStats` carries markers, called, het, nocall, unreadable and
nucleotide counts, and no BAF quantity at all, so a per-chromosome dispersion statistic does not
exist and would have to be added to the streaming pass.

**Question 2: should the noise reason be a per-sample statement that applies to every withheld
chromosome on that array, or a per-chromosome one? A per-sample statement is honest about what we
measure and useless for localising a problem. A per-chromosome one is what she would want and
invites reading a 200-marker dispersion as a measurement when it is not.**

**(c) The threshold, and what it must not disturb.** 3,065 chromosomes currently carry a reason. If
a noise branch is added, some of them will change reason. Our view is that a reclassification is
acceptable only when the new branch is a better description of the same array, that every such
change must be counted and listed, and that no chromosome may move from refused to called or back.

**Question 3: how should the threshold be calibrated given that we have no labelled set of noisy
arrays? Options we see: the corpus's own distribution with a quantile fixed in advance; the
unamplified 0.088 reference as a floor with a multiple on it; the published bulk thresholds
transferred and stated as transferred; or a per-platform calibration, since the corpus spans
Affymetrix-style probe exports and Illumina HumanCytoSNP-12. We have no principled way to choose
among these and would take a recommendation.**

## 4. What we think the answer looks like, so you can contradict it

A third refusal branch reading: the array's BAF spread at heterozygous calls is X, against Y for
unamplified DNA and a range of 0.21 to 0.30 for amplified material, so calls from this array carry
more noise than the channel can use. Threshold set from the corpus distribution at a quantile fixed
before measuring. Reported per sample, applied to every withheld chromosome of that array, with the
existing per-chromosome reasons kept underneath rather than replaced.

We think this is right and small. We also think there is a decent chance the correct answer is that
`hetBafSd` on single cells measures amplification and nothing else, in which case her ask is
answered by reporting it beside every refusal and refusing nothing on it.

## 5. One thing we are not asking

Whether to tell her. We will, either way, including if the answer is that the distinction she asked
for cannot be drawn from these files. She has been explicit that she wants to know how a claim was
validated before she uses it.

# Consult: an instrument that cannot see, and a report card that was never earned

Sample identifiers are replaced with letters throughout. Every number below is measured, and each
is marked with how far we trust it. Public GEO accessions are given as-is.

## What the tool is

OriginMarker Syngamy reads SNP-array exports from human embryo material and answers, per
chromosomal abnormality: where it is, and which parent's genome it came from. Research use only.
The material is mostly amplified single cells (single blastomere, trophectoderm biopsy, single
embryonic stem cell) plus bulk DNA where available. Typical array 825,656 markers, GRCh37. Runs
load one parent's array plus N samples; a second parent is sometimes available.

The driving question, from the lab: **"if breakpoints and losses affect paternal or maternal
genomes differentially, or equally."**

A sweep over the whole codebase has just found that most of the machinery has never run on a real
event, and that the numbers it prints its confidence in were measured by an experiment that does not
measure what it claims. This consult is about the small number of decisions we cannot make
ourselves without repeating the same class of mistake.

## What we have already settled, so you need not spend effort on it

These were open questions three hours ago. They are now closed, by measurement, and they are
engineering rather than science. They are listed so you can see the shape of the failure and skip
them.

- **No gain is ever detected. 0 of 235 constructed gains, any material, any fraction, any width.**
  Cause: a chromosome is classified as a gain when its intensity shift clears `LRR_SHIFT = 1.0`, a
  doubling. A trisomy is 2 to 3 copies, log2(3/2) = 0.585, so it never clears it and is filed as a
  loss. It is also excluded earlier: a chromosome is only considered at all once its genotype call
  rate collapses below `CALL_COLLAPSE = 0.60`, and three copies genotype perfectly well. Verified
  directly.
- **0 of 271 constructed whole-chromosome events detected, 0 of 753 segment-scanner runs fired.**
  Same entry gate: detection keys on genotype call-rate collapse, never on intensity. A chromosome
  that has genuinely lost a copy but still genotypes is invisible. Verified directly by constructing
  a full copy loss with intensity dropped 1.0 and the call rate intact: `aneuploid = []`.
- **Every detected loss is reported under a class that asserts the opposite fact.** 62 of 62
  constructed losses surfaced as `isodisomy` or `segmental-upd`, whose glosses are "both copies from
  one parent" and "segment from one parent only". The event has one copy, not two.
- **An array the tool has itself rejected still emits origin calls.** The guard empties the taxonomy
  findings only; 30 aneuploidy calls and 86 segment calls survived it across 186 lab arrays.
- **`om link` cannot tell a parent from a child.** The statistic is symmetric to four decimals in
  both directions, and 58 of 60 constructed full siblings pass its 0.020 gate. A separator does
  exist and we will implement it: presence of IBD0 tracts gave 0 percent error over 100
  pedigree and constructed pairs.
- **Loading a second parent silenced the tool completely.** Both the dosage block and the Mendelian
  block were scoped to one-parent runs, so a two-parent run emitted nothing at all. Fixed and
  verified on a public trio.

## What is sound

Worth stating, because it is the part worth protecting.

**The Mendelian channel works.** At a marker where the loaded parent is homozygous, losing that
parent's copy leaves the sample carrying an allele the parent does not have; dropout removes alleles
and cannot invent one. Under a construction that does not leak the answer, it scores **1.0000 and
0.9773** with a true parent (44 calls, 22 chromosome clusters) and **collapses to exactly 0.5000**
with an unrelated one, returning a constant verdict in that arm. That is what a working instrument
looks like: the wrong-parent arm falls to chance.

**Position is sound.** It comes from intensity and needs no parent.

## The evidence base, measured

This is the context for questions 3 and 4, and it is worse than we expected.

- `audit/blastomere-origin.ts`, which produces the shipped `OBVIOUS_EVENT_ACCURACY = 0.920`, **does
  not measure parent of origin.** On a genetically confirmed public trio it scores 44/44 with the
  true father, 44/44 with the true mother, and 44/44 with two arrays that are certainly not parents.
  Cluster-robust CI 0.0000 over 22 clusters in every arm. Mechanism, confirmed arithmetically: its
  construction builds the child's post-loss genotype out of the loaded parent's own alleles, so the
  signal exists whoever is loaded. The arm difference equals the child's heterozygous-call rate to
  three decimals (0.1475 against 0.1472).
- Worse: the parent the harness auto-picks on the real corpus is not a parent of any array it
  scores (opposite-homozygote 0.0406 or higher against a true parent-child 0.0032 to 0.0056), so the
  committed 0.920 **was already a wrong-parent measurement.**
- `audit/bands_measured.csv`, the sole source of the four `BAND_ACCURACY` rows the tool prints most
  often, **has no producer in any of 183 commits.** The two named producer scripts do not exist in
  the tree or in any commit in its history.
- Where band A could be measured against genetic truth it is **0.8933, not the shipped 0.9952**,
  n = 506 over 22 chromosome clusters, 95 percent CI [0.8375, 0.9574], which excludes the shipped
  point estimate and its whole interval.
- **Five of six committed calibration scripts inject the tool's own forward model** as the observed
  quantity and then measure whether the model agrees with it. The sixth leaks the true class.
- **Loading a stranger as the parent does not make the tool refuse.** It named a parent on 452 of
  1,320 rows (0.342), reached band A on 224 of them at confidence up to 1.0000, and was right
  0.3527 of the time [0.2319, 0.5641], against 0.8933 with the true parent.

## The corpus

- 884 array files, 726 unique after removing 146 name-and-size duplicates and 13 byte-identical
  copies under different names.
- Of 324 later-stage arrays (bulk, trophectoderm, blastomere), 10 get a confirmed parent from the
  linkage test, and after an independent Mendelian check **exactly one survives**. Zero samples have
  both parents. Two of the four study groups contain no unamplified array at all.
- Public data: one genetically confirmed trio with a trophectoderm child (GSM4472415 egg donor,
  GSM4472397 sperm donor, GSM4472424 embryo4_TE), plus haploid meiotic products GSM4774680 to
  GSM4774685. GSM4472397 and GSM4472398 are technical replicates of one man.

---

# The questions

## Q1. The intensity error model, which is the one that names the wrong parent

The tool derives a standard error for the intensity shift of a region against the rest of the same
array. **That error is 0.69 to 6.8 times narrower than the spread the same array actually shows
between its own chromosomes** (median 4.3x, n = 6 lab arrays). The detection-stage error is worse:
5.1 to 50.4 times too narrow, median 31.9x.

The consequence is measured, not hypothesised. With the error that narrow, a positive intensity
excursion of ordinary array drift resolves the copy-number class as a gain; a gain inverts the sign
map that loss and copy-neutral loss of heterozygosity share; and the call then **names the opposite
parent at confidence 1.0000, band A.** The flip threshold is 1.18 to 1.90 drift standard deviations
at a shift of 0.03 to 0.05.

**Q1a.** What is the right error model for the mean intensity of a region measured against the rest
of the same array, on whole-genome-amplified single-cell material? We currently treat markers as
independent, which is plainly wrong: amplification produces long correlated runs, so the effective
number of independent observations is far below the marker count.

**Q1b.** Is the right correction an effective-sample-size deflation, an explicit variance component
for between-chromosome drift, or a per-array empirical null built from the chromosomes themselves?
We can measure a between-chromosome spread per array; what we do not know is whether using the
array's own chromosomes as the null is legitimate when some of those chromosomes may carry the
events we are looking for.

**Q1c. This is the part we are most likely to get wrong by ourselves.** Widening the error is not
obviously safe. Too narrow and the tool inverts parents at maximum confidence, which is what it does
today. Too wide and every call becomes a refusal, which on this material is most of them already.
Is there a principled place to put it, or is the honest answer that the intensity channel cannot
support a parental call on amplified material at all and should be used for detection only?

## Q2. A threshold derived from the wrong material, deciding half the corpus

`HET_BAND_DIPLOID = 0.08` separates a genome carrying two parental contributions from one carrying
a single contribution, using the fraction of B-allele frequencies falling in the heterozygous band.
Its stated basis in the source is "diploid genomes run 15 to 16 percent in the heterozygous BAF
band, uniparental ones 1.3 to 3.4 percent", measured on unamplified DNA.

Applied to amplified arrays, **diploid embryo material on this platform reads a median band
fraction of 0.0730**, which is below the threshold. Measured consequence across 726 unique arrays:

- 316 of 726 (43.5 percent) are read as carrying one parental contribution.
- Among arrays the tool itself stages as diploid material, 62.5 percent are read that way.
- By stage: 157 of 181 single-cell (86.7 percent), 74 of 81 blastomere (91.4 percent), 85 of 225
  trophectoderm (37.8 percent).
- 122 of 726 (16.8 percent) sit within 10 percent of the threshold, which is the densest part of the
  distribution.

That single read switches off the copy-neutral detector, the uniparental-disomy detector and the
entire Mendelian origin channel. On the public single-cell array, forcing the gate to diploid makes
the same detectors find the same constructed events on the same array.

**Q2a.** Given an allele dropout rate d, what fraction of B-allele frequencies should fall within a
band of half-width w around 0.5 for genuinely diploid material? We want to derive the boundary from
a dropout model rather than refit it, because refitting a threshold to the corpus it will then judge
is exactly how the current situation arose. Is that derivation tractable, and what does it need
beyond d, w and the marker allele-frequency spectrum?

**Q2b.** Should this be one threshold at all, or should the quantity be dropout-adjusted per array
before comparison to a single boundary? The tool already estimates a dropout per array, though see
Q5 and the note below on its reliability.

**Q2c.** Is the heterozygous-band fraction even the right statistic here, against the alternatives:
the genotype heterozygosity, the count of markers where the sample carries an allele neither the
band nor the genotype can explain, or a mixture model fitted to the BAF distribution directly?

## Q3. A fix that looks too good, in a project where that has meant circular

The dosage channel resolves the copy-number class, uses it for the detection floor, for the
direction rule and for the implied mosaic fraction, and then calls the posterior **without passing
it**, so the posterior marginalises over a uniform 1/3 prior across loss, gain and copy-neutral. The
posterior accepts a `classPrior` argument that has no caller anywhere in the codebase.

Passing the class the caller already holds:

```
uniform prior:   9/48 named,  8/9 correct
class prior:    48/48 named, 48/48 correct, 40 at band A
```

It also removes an inversion in which a gain names the opposite parent.

**Q3. We do not trust this and we would like you to tell us whether we are right not to.** The class
was resolved from the same data the posterior then scores. Five of six calibration scripts in this
project failed in exactly this shape, injecting the model's own prediction and then measuring
agreement. 48 of 48 at band A has the same smell.

Specifically: is conditioning the parental posterior on a class estimated from the same intensity
and allele data a legitimate hierarchical step, or is it double-counting the evidence? If it is
legitimate in principle, what has to be true of the class estimate for it to be safe, and how would
we demonstrate that rather than assert it? If the class is uncertain, is the correct treatment a
soft prior from the class posterior rather than a hard conditioning, and does that recover any of
the 39 answers?

## Q4. What may honestly be printed

The lab's requirement was that every abnormality carry a confidence-scored parental call. Given the
section above, we cannot currently justify the numbers we print, on any material.

What we can and cannot measure:

- One material can be measured against genetic truth at all, trophectoderm, using the single public
  trio. One band within it, band A, has enough calls. That is one array.
- The lab corpus yields one array with a genuine confirmed parent, so it cannot supply a
  calibration series.
- The Mendelian channel has a sound construction and a wrong-parent control that collapses to
  0.5000, so its accuracy is measurable. The dosage channel's is not, on this data.

**Q4a.** With one array of genetic truth, is any calibrated confidence defensible, or is the honest
output a verdict with no number? The lab's stated requirement pushes toward a number, and we think
that pressure is how the current unfounded ones came to be printed.

**Q4b.** If a number should not be printed, what should? Candidates: a verdict plus the evidence
count; a bound rather than a point estimate; a two-level report where the channel is named and the
reader is told what validation exists behind it.

**Q4c.** Is there a construction we have not thought of that would generate genuine calibration data
from material we actually have? The held-out-parent construction is sound for the Mendelian channel
because removing a copy is a genetic operation with a known answer. We cannot see the equivalent for
the intensity channel, because constructing an intensity shift means assuming the very noise model
Q1 asks about.

**Q4d.** The tool caps every reported confidence at 1 - 0.206 using a rule-of-three bound from 13
zero-event validation units. Given that the underlying accuracies are now known to be unfounded, is
that cap doing anything useful, or is it lending a false precision to numbers that have no basis?

## Q5. An allele frequency used where it cannot be estimated

`DEFAULT_Q = 0.30` is the assumed frequency of the allele the known parent lacks, used when it
cannot be estimated per marker. Measured empirically on this panel it is **0.1472 and 0.1407**,
about half.

The error is one-sided. On a euploid chromosome with the true parent loaded, the shipped value makes
the Mendelian channel call "the other parent's copy is lost" on **8 to 17 of 22 chromosomes**
depending on the assumed dropout, against **0 of 22 when q is estimated from the data**.

**Q5a.** Is a single genome-wide q defensible at all, or must it be per-marker from the panel's own
allele frequencies? We have the panel and can compute per-marker frequencies.

**Q5b.** If a single value is used, should it be the mean allele frequency, the mean of the
frequency of the specific allele the parent lacks conditional on the parent being homozygous, or
something else? Those differ, and the conditioning is what we are least sure of.

**Q5c.** How should ancestry enter? The panel's expected heterozygosity is documented in this project
as 0.99, 0.93 and 0.81 of a European anchor in South Asian, African and East Asian samples. A
genome-wide q fixed for one population is a systematic error for the others, in a tool whose output
names a parent.

## Q6. Whether the driving question is answerable at all

The lab wants to know whether breakpoints and losses affect paternal or maternal genomes
differentially. The current implementation compares gynogenetic against androgenetic genomes.

**Q6a.** Those are abnormal conceptuses, arguably two different classes of fertilisation failure
with different replication histories. Is a difference between them evidence about breakage in normal
biparental embryos, or is it confounded with the mechanism that produced the class?

**Q6b.** The alternative is a within-embryo comparison on biparental material, which removes
between-sample confounding entirely. It needs per-event parental calls on that material. Given that
the Mendelian channel is the only sound one, that it needs a genotyped parent, and that exactly one
array in this corpus has a confirmed parent, is that design reachable? What is the minimum data that
would make it reachable, stated as a number of embryos and what must be arrayed alongside each?

**Q6c.** What is the right countable unit? One genome carrying a single parental contribution
produced 183 copy-neutral regions from a sliding-window detector; that is one biological fact sliced
by a window, not 183 events. Candidates: independent breakpoints, total affected megabases, a
per-genome binary, segments after merging to maximal intervals.

**Q6d.** Is there a design that answers the lab's question without needing per-event parental calls
at all? We would rather be told the question is unanswerable with this data than build a second
apparatus on the first one's foundations.

---

## What we would find most useful back

1. **Q1c and Q3.** Those two decide whether the intensity channel is repairable or should be demoted
   to detection only, and whether the one large recovery we found is real or is the same circularity
   in a new place. Everything else we can proceed on.
2. **Q2a**, because a derivation would let us set that boundary from theory instead of refitting it
   to the data it judges, and 43.5 percent of the corpus turns on it.
3. **A blunt answer to Q4a.** If the honest position is that this tool cannot print a calibrated
   confidence on amplified material, we would rather hear it now and change what we report than keep
   a number that four separate measurements have failed to support.
4. Anything in the "already settled" section that you think we have settled wrongly. We would rather
   be corrected there than proceed confidently on a bad reading, and one prior consult on this
   project overturned a conclusion we had reached on our own that had cost us a great deal of time.

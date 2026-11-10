<!-- Laboratory sample identifiers redacted: this file is tracked in a public repository.
     See scripts/release-check.sh, which refused the first commit of an earlier consult. -->

# Consult: the segmental gain statistic, its gate, and its floor

Written 8 November 2026 against version 5.47.0, commit 47e48d4. Every number below was produced
today or read out of the shipped source. Say which of our claims are wrong; three days of work rest on two of
them and we would rather lose them now than later.

## The request behind this

Our collaborator in the Egli laboratory (Columbia) asked on 7 November why the tool matches only
11 of 26 hand-called segmental gains, and offered two ideas of her own. This consult is about the
first: use marker classes beyond the opposite-homozygote one, parent-heterozygous markers
included, so a small event has enough informative sites. Her second idea, a sister cell's
breakpoint as a prior, is a separate consult.

## 1. What the detector is, exactly

`web/src/segments.ts`, `scanCopyNumber`. One function, three gates, and the indicator is the
part that matters:

1. Scan the NO-CALL indicator along the chromosome against an external null, the median
   per-chromosome absence rate over the sample's OTHER autosomes. Threshold `SEGMENT_LRT = 250`,
   empirical, taken from genomes known to carry nothing, because absence artefact on amplified
   single cells is spatially clustered and a closed-form tail produces 4 to 30 false segments per
   clean genome.
2. Require at least `MIN_SEGMENT_MARKERS = 2_400` called markers. Titrated by splicing a block of
   a real maternal pronucleus into a clean paternal pronucleus of the same series: 1,200 markers
   scored 152 to 197 against a null maximum of 139 and was not detected; 2,400 scored 431 to 556
   and was.
3. Require `|median log2R inside - genome median log2R| >= SEGMENT_LRR_SHIFT = 1.0`. The sign of
   that deviation is what labels the segment a gain or a loss.

So the detector finds a gain only where the array STOPPED GENOTYPING, and then checks the
intensity agrees. For a loss that indicator is the event: absent DNA produces no calls. For a
gain the DNA is more present, not less, and we can see no reason the no-call rate should rise at
all.

The constructed gain controls in that docstring are four of four correct at 12.6x to 100x the
threshold. We now think they establish less than they appear to. The spliced donor block was
taken from a chromosome the tool calls a whole-chromosome gain at 0.34x to 0.41x CALL RATE and
+1.60 to +1.95 log2. Those controls therefore carry a no-call excess (83% against a 12%
background) and an intensity excursion far above a third copy. They are amplification-distorted
material, not a clean trisomy. The identical construction from a euploid array returns nothing,
so the controls do measure something real; our claim is that what they measure is not the
biological event we need to detect.

**Claim 1, please contradict if wrong: the no-call indicator is the wrong indicator for a gain,
and the gain controls passed because the donor material was dropout-rich rather than because it
carried three copies.**

## 2. The gate is above the biology

A third copy is log2(3/2) = 0.585. The gate is 1.0, applied symmetrically to both directions. A
hemizygous loss is log2(1/2) = -1.0, so the gate is exactly calibrated for losses and the
measured real losses sit at -0.97 to -1.94.

Measured today over the 93 arrays of that laboratory series against 93 hand-called regions
(26 gains, 67 losses, 21 cells; the other 72 arrays carry no hand call and are the specificity
arm). Reciprocal overlap of 0.5 required for a match. Full output in
`audit/joy-dieter-baseline-2026-11-08.txt`.

```
recall            gain 11 of 26        loss 44 of 67
false-call arm    1 segment on 1 of the 72 arrays he called nothing on

the 15 gain misses, by cause
  under the 1.0 shift gate            4     deviations 0.815, 0.954, 0.994, 0.998
  boundary disagreement               6     best overlap 0.33 0.35 0.44 0.47 0.48 0.49
  under the 2,400-marker floor        2     both 3 Mb, 1,569 and 1,326 markers
  no gate explains it                 3     chr9 26 Mb, chr7 16 Mb, chr17 7 Mb

                     GAIN found      GAIN missed     LOSS found      LOSS missed
size, Mb             median 69       median 22       median 77       median 11
log2R deviation      median 1.90     median 1.14     median 1.55     median 1.16
markers in region    median 22,063   median 6,440    median 24,677   median 6,076
sample call rate     median 0.547    median 0.547    median 0.616    median 0.616
```

One of those four is a 145 Mb gain carrying 34,297 markers, refused because its deviation is
0.998. Call rate is identical between found and missed to three decimals in both directions, so
array quality separates nothing.

**Claim 2: the gate should not be symmetric. Our instinct is a direction-specific gate derived
from the copy-number prediction and the chemistry compression the code already grids (0.55 to
1.00 of the ideal log ratio), rather than a single absolute number. Is that right, and what is
the defensible form? If the honest answer is that single-cell log R cannot support a gain gate at
all without an allele channel, say so, because that reorders everything below.**

## 3. The model we already have, and where it does not reach

`web/src/gainOrigin.ts`, `callPair`, is wired and shipped. It weighs the allele share against
intensity over six states:

```
states        disomy, gain-paternal, gain-maternal, loss-paternal, loss-maternal, unexplained
share pred    gain  d = +/- f / (2 * (2 + f))        loss  d = +/- f / (2 * (2 - f))
lrr pred      gain  k * log2((2 + f) / 2)            loss  k * log2((2 - f) / 2)
grids         f in 0.1 .. 1.0 by 0.1        k in {0.55, 0.70, 0.85, 1.00}
'unexplained' lets the share float free and draws intensity from the cell's own spread
output        state, posterior, runner-up, best-fit f and k, plus a TILT reported separately
```

Validated under a pre-registration committed before scoring
(`audit/gain-origin-split-preregistration.txt`): split by CELL not chromosome, even accessions to
tuning and odd to held-out, real hand-called events excluded from both halves and read once at the
end. Over 300 constructed builds it names a state on 267 and the tilt is right on 262, including
94 of 98 at a cell fraction of one third where the shipped thresholds detect nothing at all. On
held-out cells the direction was wrong on 1 of 121 named rows against a 5 percent bar fixed in
advance; the tilt was wrong on 4 of 287 against a 1 percent bar, which it misses, so the parent is
not taken from the model. Every one of the 22 wrong-direction calls sat on a chromosome whose own
baseline intensity was between +0.18 and +0.49, against -0.05 to +0.10 on every arm that never
erred.

**`callPair` runs per chromosome. It has never been given a segment.** The segment path has no
allele channel whatsoever. That looks to us like the whole gap, and it is a reuse rather than a
new model.

The obstacle is marker supply inside an interval. The allele channel's informative class is
markers where the two parents are OPPOSITE HOMOZYGOTES, which is 3.4 percent of calls, about
1,000 per chromosome, roughly 8 per Mb at this density. The floor is
`ALLELE_GAIN_MIN_INFORMATIVE = 200`, with `ALLELE_GAIN_Z = 6`, `ALLELE_GAIN_SHIFT = 0.08`,
`ALLELE_GAIN_MIN_REFERENCE = 6` reference chromosomes, and an intensity gate of
`ALLELE_GAIN_LRR_Z = 2` of the cell's own chromosome-intensity spreads. At 8 informative markers
per Mb a segment needs about 25 Mb before it can reach 200. Half the hand gains are smaller than
that. That is our collaborator's point restated as arithmetic, and it is why her suggestion is
not merely a convenience.

## 4. What we already know about the other marker classes

`hetCallOneParent` is implemented, tested with eight assertions, and deliberately NOT wired into
any verdict (5.44.3 and 5.44.4). On 1,021 chromosome-observations against a published karyotype it
scores better than the shipped class, 0.9667 against 0.9068, and better on the hard arm, 0.777
against 0.490 on one-contribution cases.

It stays out for a measured reason. A sample carrying two haplotypes of ONE parent is
heterozygous wherever a crossover separated them, so it reads high on the parent-heterozygous
class while no second parent contributed anything. On real arrays of a series with a genotyped egg
donor and no father, the two samples carrying her genome alone read 0.3408 and 0.3734 against
0.3974 for the biparental sample in the same series. Under 0.03 of separation, so no boundary
divides them, and the substitution calls both biparental: it asserts a paternal contribution that
does not exist.

**Claim 3: that failure is about asserting a second parent's PRESENCE, and it does not transfer to
a DOSAGE use inside an interval already known to carry an event, where the question is how many
copies rather than from how many parents. If that reading is wrong, the whole of our collaborator's
first idea is closed and we would rather hear it now.**

Our intended construction, stated so it can be corrected. At a marker where one parent is AB and
the other is homozygous, the hard genotype call carries little, but the embryo's CONTINUOUS BAF
carries dosage: with the mother AA and the father AB, a normal cell reads near 0 or near 0.5 and a
cell with an extra paternal copy reads near 1/3 or 2/3. So the evidence is in the BAF, not in the
AA/AB/BB call, and the class is roughly an order of magnitude larger than the opposite-homozygote
one.

Three things we do not know how to settle:
- what dropout does to those band positions, and whether the shift survives it;
- whether `hetBafSd`, the measured spread of BAF at heterozygous calls (0.088 for unamplified DNA,
  0.21 to 0.30 for whole-genome amplified material), already swamps a 1/6 band separation;
- whether to model per-marker likelihoods or a windowed share, given that this project's one
  recorded disaster came from centring a share statistic on the median of a distribution whose
  median IS a boundary of its support.

## 5. The floor, which we believe should be derived rather than tabulated

An earlier consult on 25 October established two things we are relying on:

- marker floors should be derived as N = ln(R) / KL for required posterior odds R, not tabulated;
- the Bernoulli KL is asymmetric, the project's recorded rule wrote KL(h0 || h1), and the quantity
  that governs whether a true event clears the threshold is KL(h1 || h0), because the scan's score
  has expectation w * KL(h1 || h0). Applied to the opposite-homozygote set that reached window 300
  at LRT 250, valid for material up to roughly 43 percent dropout.

**Request: derive the floor PER MARKER CLASS for the gain hypothesis.** A gain is a weaker
contrast than a loss in both channels, and this number decides whether the parent-heterozygous
class can carry a 14 Mb event or whether nothing can. If the derivation says the smallest
detectable gain is 40 Mb whatever we do, that is a publishable limit and we will write it down
instead of building a feature.

## 6. Measured limits we are not asking you to overturn, only to work within

- Mosaic gains at cell fraction 0.3 and 0.5 are 0 of 120 detected. The share moves +0.012 rather
  than +0.10 because dropout has already spread the markers. 300 constructed builds.
- The ratio lrr / |share shift| is +2.0 for gains at every mosaic fraction tested (2.04, 1.96, 2.01
  at f = 0.3, 0.5, 1) and -0.75 to -0.97 for losses. This is the fact that makes a pair likelihood
  work at all.
- A second platform transfers the statistic but not the tool: on Illumina HumanCytoSNP-12 data
  (Zamani Esteki 2015, PMID 25983246), of 8 published parent-named events 5 move by |z| >= 2 and
  all 5 name the right parent, including one cell carrying extra maternal chr22 and extra paternal
  chr6 in opposite directions. The shipped pipeline scores 0 chromosomes on that chip because the
  stage call reads haploid or failed.

## 7. The questions, in priority order

1. Is Claim 1 right, and if so what indicator should a segmental gain scan use?
2. What is the defensible form of the gain intensity gate, given compression of 0.55 to 1.00 and
   per-chromosome baseline offsets of +0.18 to +0.49 that a single cell cannot measure?
3. Is Claim 3 right: does the 5.44.4 blindness argument leave a dosage use of parent-het markers
   untouched?
4. Derive N = ln(R) / KL per marker class for the gain hypothesis, orientation resolved, and say
   what the smallest honestly detectable segmental gain is in each class.
5. Per-marker likelihood or windowed share, and if a share, how to centre it without repeating the
   median-is-a-boundary error.
6. Anything in sections 1 to 6 that is wrong.

## 8. What we will do with each answer

A floor that admits a 14 Mb event in the parent-het class: wire that class as a dosage statistic,
pre-register the bar (recall rises, the 1 false segment across the 72 quiet arrays does not), and
give `callPair` a segment. A floor that does not: write the limit into `docs/validation.md`, report
the counts per region so our collaborator can see for herself which events could never be called,
and spend the time on the sister-cell prior instead.

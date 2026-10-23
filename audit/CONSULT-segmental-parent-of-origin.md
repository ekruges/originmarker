# Consult: naming the parent on a segment, and a statistic that failed its bar by one chromosome

Public data throughout. No laboratory sample is named. Every number is measured on real arrays
with a published by-eye answer key, and each is marked with what it rests on.

## The instrument, in one paragraph

OriginMarker reads SNP-array exports from human embryo material and answers, per chromosomal
abnormality: where it is, and whose copy it was. Affymetrix Axiom, 825,656 markers, GRCh37.
Material is mostly amplified single cells, so dropout is the dominant noise. A run loads one
parent's array plus N samples, and sometimes a second parent.

The core statistic for a whole-chromosome absence is the opposite-homozygote rate against a
parental array: a haploid product of that parent reads 0.17 to 0.75 percent, a non-product 6.5 to
16.6. Roughly 24,000 of the 825,656 markers are informative, a median 1,034 per chromosome.

## The data this consult rests on

```
arrays    GSE148488, Zuccaro et al., Cell 2020;183(6):1650-1664
truth     Table S3, the authors' own by-eye karyotype of every arrayed cell,
          committed as audit/zuccaro-table-s3.tsv: 682 chromosome calls over 31 cells
parents   both arrayed in bulk for all 31, donor assignment read from GEO's own
          !Sample_source_name_ch1, rebuilt by audit/geo_trio_manifest.py
no simulation anywhere. Every genotype is measured; every truth label is the published hand
reading of that same array.
```

## Question 1. Two channels name the parent on a segment. Which one wins?

A `parental-absence` segment is found by scanning, at markers where ONE parent is homozygous,
whether that parent's allele is missing across a window. Until this week that scan ran against the
first parent only, so a deletion of the second parent's copy read present at every marker and was
never looked for: 11 of 17 paternal segmental events found and **0 of 10 maternal**. The scan now
runs a second time against the other parent's own homozygous markers, which is a different marker
set of the same size, so the window minimum and the threshold carry over unchanged. Maternal
detection goes to 5 of 10.

That leaves two statements about whose copy it was, on the same window:

- **Which run found it.** The marker set is defined by one parent's homozygous genotypes, so the
  run that produced the segment names the parent by construction. It carries no confidence of its
  own: if the segment is real, the parent is that parent.
- **The allele share inside the window.** A separate statistic over the same region, with a
  calibrated confidence and band, which reads a loss the same way round (the under-represented
  parent is the one that lost the copy).

**What we did, and want challenged.** The run now wins, the share supplies the confidence where it
agrees, and where the share names the OTHER parent the row names no parent and says the two
readings disagree. The reasoning is that the run is definitional rather than inferential, and that
printing the share's confidence beside a parent half the evidence contradicts would be a
confidence calibrated on runs that had no conflict.

Three ways that could be wrong, and we cannot tell which:

1. A segment can be real while its *extent* is wrong. If the window overruns the event, the run
   still names its parent while the share, averaging over the true event plus flanking normal
   sequence, goes quiet or tilts the other way. Then the conflict is about the edges, not the
   parent, and declining is an overreaction.
2. Dropout is parent-blind, but the two marker sets are not the same markers. If one parent's
   homozygous set is systematically denser in a region of poor amplification, the two runs are not
   equally powered there, and preferring the run silently prefers whichever set is denser.
3. The share has a measured accuracy; the run does not have one at all, because there is nothing
   to measure. "Correct by construction" is an argument, not a number, and this codebase has been
   wrong before in exactly that register.

**What the 31 cells say about it.** The ladder's conflict branch has never fired on real data.

```
parental-absence segments over 31 cells          69
carrying which run found them                    69   (43 paternal run, 26 maternal run)
with BOTH a run and a window allele share        23
  the two channels agree                         23 of 23
```

So the field is always present where both parents are loaded, and the two channels have never once
named opposite parents. The conflict rule is therefore asserted by a unit test and by nothing else.
That is not a reason to remove it, but it is a reason not to claim it has been validated.

**What we are asking.** Is the ladder right? If not, what is the measurement that would decide
it, on 31 cells with a published hand karyotype and nothing else?

## Question 1b. 35 of 69 segments are windows that BOTH runs called

This is the finding we most want looked at, because we have withheld a parent on the basis of an
argument rather than a measurement.

Both runs calling one region cannot be a single-copy event: a lost copy takes one parent's alleles
with it. Both parents' alleles depleted over one window is either a region missing both copies,
which the intensity channel calls as a copy loss on its own evidence, or dropout wide enough to
clear the window minimum of 2,400 markers. Measured:

```
one window called by BOTH runs                      18 pairs
segments involved                                   35 of 69
cells carrying any                                  6 of 31
those cells carrying a hand-called segmental event  4 of 6
hand-called events that lose their parent to this   0
detection, before and after                         11 of 17 paternal, 5 of 10 maternal
```

Half the parental-absence segments in the corpus are in one of these pairs. The flag costs nothing
measurable: no hand-called event loses its parent to it, and detection is identical before and
after, because the four affected cells that do carry a hand-called segmental event carry no
DETECTED one. That is a statement about what the flag costs, not a defence of the cells.

The spans are the tell. They are nearly chromosome-length: chr3 30.7 to 176.4 Mb against 39.9 to
176.4 Mb on one cell, chr8 0.2 to 144.8 against 0.2 to 130.5 on the same one, and so on across
chr11, chr12, chr13, chr14, chr15, chr16, chr17 and chr20. One cell carries four such pairs, and
one segment is in two of them. That reads as a globally depleted genome rather than as located
events.

**What we did.** The region is kept, the parent is withheld, and a limit line says both runs called
it and names the two candidate explanations. Nothing is filtered: the laboratory has asked
explicitly that chaotic cells stay in rather than being screened out, so suppressing the cell was
not an option we considered.

**What we cannot settle.** Three different things are consistent with the measurement.

1. **Withhold the parent, keep the region.** What we shipped. Honest, and it costs no hand-called
   event, but it leaves a near-chromosome-length region in the segment list whose only content is
   that it could not be read.
2. **Withhold the segment calls on such a cell entirely.** If a cell produces four of these, the
   scan's null is wrong for that cell and every segment it reports on it is suspect, including the
   34 that are not contested. We have no measurement that separates those two populations, and the
   fact that contested segments are half the corpus rather than a handful is the strongest argument
   for this option.
3. **Resolve it by intensity.** If the log R over the window says DNA is missing, it is a copy
   loss and belongs on that channel; if it does not, it is dropout and belongs nowhere. The
   intensity channel already exists and is already run. We did not route the contested windows
   through it because choosing that route after seeing which cells were affected is the same
   mistake this directory keeps recording.

**What we are asking.** Which of the three, and what fixed-in-advance rule decides it? If the
answer is the intensity route, what is the statement of it that can be committed before the
contested windows are scored, given that we already know which seven cells carry them?

## Question 2. A pre-registered statistic met one bar and missed the other by one chromosome

The absence scan reads markers where ONE parent is homozygous, where a lost copy still reads
present about half the time because the other parent may be heterozygous. The alternative reads
only markers where BOTH parents are homozygous for OPPOSITE alleles, where a lost copy shows at
every marker, on a set roughly fifteen times smaller: about 1,034 markers per chromosome against
about 16,000.

Whether the undiluted signal pays for the marker count is not settled by arithmetic. A power
calculation written in advance said the diluted set should win. A first look said the opposite,
and that first look chose its operating point on the cells it then reported. So the split, the
rule for picking the operating point, and the bar were committed before the held-out half was
scored: `audit/segment-statistic-preregistration.txt`. 31 cells, 16 tuning and 15 held out, split
by cell and sorted by accession, no reassignment.

```
held out, at the point the rule picked   deletions found   segmental only   on hand-normal
shipped, window 2400 LRT 250                   76 of 102        12 of 14       12 of 224
opposite-homozygote, window 150 LRT 80         85 of 102        13 of 14       13 of 224
```

The bar required at least as many deletions AND no more chromosomes the hand reads as normal. It
met the first and missed the second by one chromosome, so the shipped statistic stays.

Three things went wrong with the design, all recorded rather than patched:

- **The halves came out badly uneven in what matters**: 50 hand-called deletions on the tuning half
  against 102 on the held-out half, and 253 hand-normal chromosomes against 224. Splitting by cell
  cannot control for how many events each cell carries. The secondary endpoint that asked whether
  the halves differ is what caught it.
- **The rule walked to the corner.** It picked window 150, LRT 80, the most aggressive point on the
  grid, because the cap it was given was the shipped statistic's own rate on the loose tuning half.
  "Detect the most subject to that cap" then has nothing holding it back. That is a defect in the
  rule, not in the data.
- **A point that dominates is sitting in the held-out results and cannot be adopted.** Window 300,
  LRT 250 finds 75 of 102 against the shipped 76, and calls 5 of 224 hand-normal chromosomes
  against the shipped 12. One fewer event for less than half the false calls. It is read off the
  held-out half, so a point chosen that way carries none of the protection the pre-registration
  exists to provide, and it is recorded and NOT adopted.

**What we are asking.**

1. Is there any legitimate route to window 300 / LRT 250 on this corpus, or is a second corpus the
   only answer? Our own reading is that it is the only answer, and we would like that contradicted
   if it is wrong.
2. If a second corpus it is: GSE60909 is Illumina HumanCytoSNP-12, both parents, a published
   per-cell answer key, a second platform and a second laboratory. What should the
   pre-registration say so that the rule does not walk to a corner again? Specifically, what caps
   a specificity constraint when the only available reference rate is the incumbent's own, measured
   on a half that is loose?
3. A call on a chromosome the hand karyotype reads as normal is not proven wrong: the authors
   tabulated the losses they were looking for, so a real untabulated event lands there. We use it
   as a relative measure between two statistics on the same chromosomes of the same cells, never
   as an error rate. Is that use sound, and is there a better one available on public data?

## Question 3. What limits detection

16 of 27 hand-called segmental events are detected, 11 of 17 paternal and 5 of 10 maternal. The
shipped whole-chromosome channel separately reaches most whole-chromosome losses: on the same 31
cells it agrees on 475 intact and 78 parent calls, misses 29, and names the wrong parent once. Of
the 29 misses, a majority are on cells whose stage reads `failed`.

We have not established what the 11 missed segmental events have in common, and we are not going
to guess at it in a changelog. If the answer is "those events are below the window minimum", that
is a different fix from "those cells are too noisy", and both are different from "the hand
karyotype is reading something the array cannot carry".

**What we are asking.** Which of those three is testable on this corpus without choosing the test
after seeing the answer?

## What has already shipped from this work, so you need not spend effort on it

- The second absence run, against both parents. Defect fix, measured separately: maternal
  segmental detection 0 of 10 to 5 of 10, with the window allele share naming the depleted parent
  correctly on every event the scan reaches. Shipped 5.38.0.
- Three display surfaces were labelling every parental-absence segment with the loaded parent's
  role, one of them hardcoded to `paternal`. Correct while one run existed, wrong on every segment
  the second run finds. All three now read the segment's own field, and a seam test asserts that
  any file carrying that label reads it, so a fourth surface fails instead of shipping a wrong
  label.
- The trio manifest every public audit depends on existed only as a hand-written file beside the
  arrays and could not be rebuilt by a reviewer. It is now generated from GEO's own
  `!Sample_source_name_ch1` by `audit/geo_trio_manifest.py`, which refuses a row it cannot source.

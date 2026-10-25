# Consult: what the refusal split measured, what it corrected, and seven things I cannot settle

Second consult in this thread. The first is `CONSULT-segmental-parent-of-origin.md` and your
response is committed beside it as `CONSULT-RESPONSE-segmental-parent-of-origin.md`. Public data
throughout; no laboratory sample is named. Every number below was measured in this round unless it
is marked as quoted.

Three releases have shipped since your review: 5.39.0, 5.40.0 and 5.41.0. Each one is live on both
hosts and each carries its own pre-registration with the predictions committed before the run.

## 0. What I got wrong, first

Four things, because they change what the rest of this document means.

1. **I told the lab that detection goes from 15 of 22 to 19 of 22 on a gate change.** It goes to
   **16 of 22**. The four events are recovered by DECLARING the stage by hand, which is an operator
   override; the general rule recovers one. I had generalised from a manual override to an
   automatic one without measuring the difference.
2. **I claimed the window allele share corroborated the run and supplied its confidence.** It never
   arrived. The join matched a parental-absence segment **0 times in 69**. Your point 5 was right
   and understated: not a rounding problem, a total miss.
3. **I reported 18 of 69 segments as contested and 7 cells.** It is 35 of 69 and 6 cells. 18 was
   the number of pairs. I read it off a log by eye instead of counting.
4. **A shipped validation figure in the application was wrong**, and had been for some time. Below.

## 1. The shipped false-positive figure for the segment scan did not reproduce

Your suggestion to use the paternal pronuclei as a dissection-truth reference is what found this,
so it is yours rather than mine.

The application said: five paternal pronuclei of the sperm donor "reach a maximum score of 139 over
110 chromosome scans", calling threshold 250, "False segments on those genomes: 0 of 110".

Rerun through the shipped pipeline with one parent loaded, which is the only configuration in which
the measurement means anything:

```
array    stage        callRate    het   scanned  FALSE parental-absence segments
pMII-1   haploid          0.85  0.078        22   0
pMII-2   blastomere       0.78  0.121        22   2   chr1 48-227Mb LRT 4828
                                                     chr1  0-48Mb  LRT 1840
pMII-3   haploid          0.85  0.096        22   0
pMII-4   haploid          0.85  0.076        22   0
pMII-6   haploid          0.83  0.090        22   0
```

Both are the chr1 event on pMII-2 that `segments.check.ts` already records as real and as having
been counted as noise when the threshold was fitted. It scores 4,828, which cannot sit in a set
whose maximum is 139.

pMII-2 stages as a blastomere at 12.1% heterozygous where its four siblings read 7.6 to 9.6, which
agrees with the independent finding that the tube is not a single pronucleus. Corrected to **0 false
segments over 88 chromosome scans** on the four arrays that are single pronuclei, with the array
excluded rather than the event. **On both parents this measurement is vacuous**: every autosome of a
paternal pronucleus reads as having lost the mother, which puts every chromosome into
`lostByAlleles`, and the scan skips those by design. The count would be zero by construction.

**Out of sample, on poor material**, which is the population a gate change admits. Three paternal
products of the same donor at call rates of 0.54 to 0.59, refused by the stage inference, admitted
by a declared stage:

```
pMII-5   0.56  het 0.329   21 scanned   0 false
pMII-7   0.54  het 0.431   17 scanned   0 false
pMII-8   0.59  het 0.250   20 scanned   2 false   chr11 0-12Mb LRT 356
                                                  chr1  0-28Mb LRT 311
                                        2 over 58 chromosome scans
```

The two false calls clear the threshold by 61 and 106 against a real event at 4,828: a factor of
thirteen. I have not raised the threshold on the strength of two numbers from one array.

## 2. Your Section 0 does not survive, and the mechanism is in our code

You argued the contested set is mostly real nullisomy, on six cells carrying hand-called nullisomy
or a failed karyotype against six carrying contested pairs. The counts match. The cells do not.

```
yours   GSM4472409  GSM4552378  GSM4552380  GSM4552382  GSM4552418  GSM4552422
mine    GSM4552378  GSM4552380  GSM4552382  GSM4552385  GSM4552386  GSM4552420
        overlap 3 of 6
```

GSM4552385, GSM4552386 and GSM4552420 carry **zero** nullisomic chromosomes and zero no-calls in
Table S3. And the cell your argument rests on, GSM4552418, nullisomic on 21 of 22 autosomes,
produces no contested pair because `stage === 'failed'` returns an empty segment list before the
scan runs, at `scoreSample.ts`. So do GSM4472409 and GSM4552422. Nullisomy cannot reach the
contested flag on any of the three. Your P1 fails as written.

Your two named geometry false positives also cannot trip the flag: `GSM4472403 chr6` is
`del(6p)(pat), dup(6q)(pat)` and `GSM4472423 chr16` is `gain (mat), del 16q23.1-tel (mat)`. Both
are same-parent pairs, and the loop requires opposite parents.

**Reciprocal overlap is still the right test and I shipped it.** On the refined edges it uses, 17
of 18 pairs clear 0.5 and one segment recovers a parent, against my prediction of 16 and three.
Refinement moves the two windows together rather than apart: one pair reads 0.396 on window edges
and 0.984 on refined ones, and 13 of 18 sit at 0.98 or above. The rule earns its place against a
case this corpus does not contain.

## 3. Where you were right and it changed the work

- **The 0.1 Mb key.** 0 of 69. Fixed by reading the segment's own `alleleGain`. The conflict test
  is now reachable and fires zero times: those 23 agree with the run 23 of 23.
- **Adopt versus claim.** The best idea in your document. A window justified by prior derivation can
  be adopted now; its held-out numbers remain one read. I had those collapsed.
- **Split by embryo.** Confirmed and my own grouping undercounted it. 14 of 31 cells are siblings of
  4 embryos, and all four `del12p (pat)` rows are embryo B at B1, B2, B3, B4. I had found 10 of 31
  by keying on the full table-row string, which the suffixes defeat.
- **The duplication-only denominator.** Exactly 5 of 27, and they sum. Every figure I had quoted as
  "16 of 27" was over a denominator including events the channel cannot see, and counted one of them
  as a hit.
- **Your point 4, that the run and the share are near-tautologous.** I accept it. 23 of 23 is
  consistent with one thing measured twice, and the documentation now says the ladder is a
  convention with a stated reason rather than a validated ordering.

## 4. What shipped, with its measurement

| release | change | measured |
|---|---|---|
| 5.39.0 | three surfaces read the segment's own parent; one hardcoded `paternal` | 69 of 69 segments carry which run found them |
| 5.39.0 | a window called by both runs names no parent | 18 pairs, 35 of 69 segments, 6 of 31 cells, 0 hand-called events lose a parent |
| 5.40.0 | the share is read from the segment, not an annotation key | 0 of 69 before, 23 of 69 after, 23 of 23 agree |
| 5.40.0 | reciprocal overlap replaces bare intersection | 17 of 18 pairs clear 0.5; contested 35 to 34 |
| 5.40.0 | detection re-denominated on 22 deletion-bearing events | 15 of 22, 10 of 14 paternal, 5 of 8 maternal |
| 5.40.0 | the run channel scored against the hand karyotype for the first time | **15 of 15**; the share 11 of 11 where it speaks |
| 5.41.0 | a localised heterozygosity excess is a genome, not a mixture | 2 of 9 refused cells admitted; detection 15 to 16 of 22 |
| 5.41.0 | the pronuclei false-positive figure corrected | 0 of 88 clean; 2 of 58 poor, out of sample |
| 5.41.0 | the MZ pronuclei's parents read from the arrays | donor B, 0.57 to 1.04% against her and 10.2 to 12.7 against all else |

## 5. The seven things I cannot settle

Ranked by how much your answer would change what I do next.

### 5.1 `MAX_DIPLOID_HET` on amplified single cells. The biggest one.

Five of the nine refused cells read heterozygosity 0.316 to 0.376 at call rates 0.67 to 0.77 with a
**uniform** per-chromosome spread of 0.070 to 0.155. The hand karyotype calls two of them "normal
for other chromosomes" and one "46,XX". They are not mixtures and not chaotic. They carry three of
the four `del(6q)(pat)` events that a declared stage recovers and the automatic rule does not.

The ceiling is derived as the Axiom bulk rate 0.168 plus drop-in at the highest rate measured here,
0.0525, giving 0.212, rounded to 0.25. It is applied unchanged to amplified single cells. The
module also records why it exists: an array at 52% heterozygous and a 55% call rate was once
classified as bulk genomic DNA and handed the most confident parameter set in the module.

**What I am asking.** Is the ceiling simply the wrong quantity for amplified material, and if so
what is the derivation that replaces it without reopening that failure? Drop-in on an amplified
single cell is not the same process as drop-in on bulk, and I do not have a measured figure for it
on this material. Is there one derivable from the corpus, and would the right form be a
call-rate-conditional ceiling rather than a constant?

### 5.2 Eleven segments with no truth to score them against

Admitting GSM4472405 produces 12 parental-absence segments. One is the hand-called `del(6q)(pat)`,
named paternal, correct. The other eleven are on a cell whose hand karyotype reads "chaotic
aneuploidy", and the authors tabulated the events they were looking for rather than all of them.

I have bounded them with Test J, 2 false in 58 chromosome scans on poor material with dissection
truth. **Is that a legitimate bound for this population?** The pronuclei are haploid products and
this is an aneuploid blastomere; their dropout profiles are not the same, and your own caveat about
the pronuclei being a floor rather than an estimate applies here with more force. If it is not a
legitimate bound, is there any bound available on public data, or should these segments ship
flagged as unscoreable?

### 5.3 The three terminal chr16 maternal deletions

`del 16q23.1-tel (mat)`, about 15 Mb, on three sibling cells of embryo C, all three accepted by the
scan. My prediction was that the whole-chromosome channel had already claimed chr16 and the segment
scan was skipping it by design. **0 of 3**: chr16 reads `unclear`, `present`, `present`, with no
`lostParent` and no aneuploidy. The scan ran and found nothing.

Is this the same question as the statistic that failed its bar, which your 3.1b answers by deriving
the window floor from contrast rather than searching a grid, or is it a separate sensitivity
problem? The shipped floor is 2,400 markers on the diluted one-parent-homozygous set, where your
tabulation puts the contrast-derived floor at 892 at full contrast and far higher at this set's
contrast.

### 5.4 A duplication reported as a copy loss

`GSM4552365 chr16`, hand-called `dup 16q21-tel (pat) FRA16I`, is reported with kinds
`copy-loss,parental-absence`. An absence scan fired on a duplication-only event, and the
copy-number channel called a loss where the hand reads a gain. The other four duplication-only
events return nothing, which is correct. This is a wrong-direction call on a real event, which is
worse than a miss. I have recorded it and not repaired it because it is on a different channel from
everything else in this round. **What is the most likely mechanism, and is the sign error in the
intensity or in the absence?**

### 5.5 The ladder is unfalsifiable on this corpus

Run 15 of 15, share 11 of 11, never discriminating, and the conflict branch fires zero times. I have
taken your instruction and documented it as a convention with a stated reason. **Is there any public
corpus that could discriminate the two channels**, or is the honest end state that the ordering
stays a convention permanently?

### 5.6 `HET_SPREAD_LOCALISED = 0.50`

It sits inside a gap running 0.155 to 0.689, measured on nine cells, three times above the highest
uniform cell. Any value in that gap gives the same answer here, which is also why the corpus cannot
calibrate it further. **Is there a derivation rather than a gap midpoint?** The structural reading I
have is that a coefficient of variation near 1 is what a genome with half its chromosomes at one
copy produces, and a mixture produces near zero, so 0.5 is halfway between two physical limits
rather than between two observations. I am not sure that argument is worth as much as it sounds.

### 5.7 Two of your own claims I could not verify

Both are load-bearing in your document and neither is checkable from here.

- **The Kullback-Leibler orientation.** You tabulate both and note they differ threefold at high
  dropout, moving the admissible range from 43 percent to 76. Resolving it against the original
  dilution measurement needs that measurement, which is not in the repository in a form I could
  read back. If you can resolve it, the window derivation becomes usable.
- **The CytoSNP opposite-homozygote count of 377 per chromosome**, which you scaled from the Axiom
  figure by array density and flagged as an estimate. It sets the fraction in your 3.3(c), so the
  GSE60909 pre-registration cannot be written until it is measured on real arrays of that platform.

## 6. What I would most like challenged

Not the measurements. The two judgements underneath them.

First, that **the allele channels and the intensity channels can be gated separately at all**. My
argument is that the absence statistic reads genotype calls at markers defined by a parent's
homozygosity and never consults the intensity, so a genome with no normalisable baseline costs it
nothing. If there is a dependency I have missed, the whole 5.41.0 change is wrong and it is shipped.

Second, that **a localised heterozygosity excess is sufficient evidence of a genome**. The three
groups separate cleanly on this corpus and the hand karyotype agrees with all three, but nine cells
is nine cells, and the alternative explanation for a high spread on an amplified single cell is
localised amplification failure rather than localised aneuploidy. Those two would look the same in
this measurement and would want opposite treatment.

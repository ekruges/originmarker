<!-- Laboratory sample identifiers redacted: this file is tracked in a public repository.
     See scripts/release-check.sh, which refused the first commit of it. -->

# Response: segmental parent-of-origin, the ladder, and the statistic that missed by one

What this rests on. The repo at /Users/ezrakruger/claudecodegeneralworkspace/originmarker,
commit state of 2026-10-24. audit/zuccaro-table-s3.tsv recounted from scratch (682 chromosome
calls, 31 cells, 27 del/dup rows all parent-named, 17 pat / 10 mat - all four figures in the
consult reproduce). audit/segment-statistic-preregistration.txt read in full. The scan itself
read in web/src/segments.ts, the contested flag in web/src/scoreSample.ts, the ladder in
web/src/defects.ts.

What I could NOT do: the GSE148488 arrays are not on this machine (the local array directory holds the
laboratory series, not GSMs) and node is not installed, so segment-origin.ts and
segment-sensitivity.ts could not be re-run. Every detection number below is quoted from your
documents, not reproduced. Everything I assert as new is derived from the truth table or from
the shipped source, both of which I did read.

Zuccaro et al., Cell 2020;183(6):1650-1664, DOI 10.1016/j.cell.2020.10.025 - checked for
retraction/correction: none found. It carries a same-issue commentary (PMID 33306952), which is
commentary, not a correction.

Headline: two of your three questions have a defect upstream of the decision you are asking me
to make, and in both cases the defect is in the shipped code rather than in the argument.
Question 1b is the one you flagged as most wanting review and it is the one where I think you
are closest to right - but for a reason you did not give, and the rule I would fix in advance is
not any of your three.

And on Question 2 you asked to be contradicted if "a second corpus is the only route" was wrong.
It is wrong, and the contradiction comes from your own prior records: the project already
established that the segment floor is contrast-dependent and that floors should be derived as
N = ln(R)/KL rather than tabulated. Applied to the opposite-homozygote set that rule reaches
window 300 at LRT 250 from measurements predating this pre-registration - though only for
material up to roughly 43 percent dropout once the KL orientation is taken conservatively, which
is a scope limit this corpus's `failed` cells make real. The deeper issue is
that your document treats "adopt this operating point" and "claim this operating point's
numbers" as one act; they separate, and only the second needs GSE60909. See 3.1b - I reached the
wrong answer first and have left both in.

---

## 0. The finding that reorders everything: your contested set is mostly real biology

You describe the 18 both-run pairs as "either a region missing both copies ... or dropout wide
enough to clear the window minimum", and you read the near-chromosome-length spans as "a
globally depleted genome rather than as located events."

There is a third explanation you did not list, and the truth table says it is the dominant one.

A NULLISOMIC chromosome has both parental copies absent. Both absence runs SHOULD call it. That
is not an artefact, not dropout, and not a failure of the scan's null - it is the correct
behaviour of two correct runs on a chromosome that genuinely lost both copies.

Counted in audit/zuccaro-table-s3.tsv:

```
hand-called nullisomic chromosome-cells                42
cells carrying at least one                             5
cells that are entirely 'no call' (hand: "failed")      1
cells matching either description                       6 of 31
your contested-pair cells                               6 of 31
```

Six and six. And the per-cell detail lines up with the spans you quote:

```
cell         nullisomic  no-call  monosomy  normal   hand label
GSM4552418       21         0         0        0     "chr6p fragment only"
GSM4472409       10         0        12        0     chaotic aneuploid
GSM4552378        5         0         7        9     chaotic aneuploid
GSM4552380        5         0         8        7     chaotic aneuploid
GSM4552382        1         0         9       10     chaotic aneuploid
GSM4552422        0        22         0        0     "failed"
```

GSM4552418 is nullisomic on 21 of 22 autosomes and the authors' own karyotype column calls it
"chr6p fragment only." A cell that retains one arm and nothing else SHOULD produce both-run
calls on chr3, chr8, chr11, chr12, chr13, chr14, chr15, chr16, chr17 and chr20 - which is very
nearly the chromosome list you report. "One cell carries four such pairs" is the expected
reading of a cell whose hand karyotype is a single surviving arm, not evidence that the null is
broken.

So the contested flag is not firing on noise. It is firing, mostly, on the one event class where
both runs are simultaneously right, and it is converting a correct and agreeing pair of readings
into "no parent named." On these cells the honest answer is not "whose copy was lost" but
"both were" - which is a statement your data structure can already carry (`verdict:
'contradiction'` in bothParentsOrigin.ts, explicitly documented as nullisomy) and which the
segment channel throws away.

One caveat I want to state rather than bury: the chromosome lists overlap but I cannot prove
pair-by-pair identity without running the scan. The 6-and-6 cell count is an exact match and the
span pattern matches GSM4552418's karyotype specifically. Treat it as a strong and testable
prediction, not a closed result. The test is one run of segment-origin.ts with the hand label
joined on - see the rule in section 1 below.

### The second defect, in the same block

web/src/scoreSample.ts, the contested loop:

```js
if (a.startBp > b.endBp || b.startBp > a.endBp) continue
a.contested = true
```

Bare interval intersection. No reciprocal-overlap fraction, no requirement that the two windows
describe the same event. At near-chromosome spans, any two events on one chromosome overlap
trivially - so a genuine compound event trips the flag. The truth table has exactly those:

```
GSM4472403 chr6   del(6p)(pat),dup(6q)(pat)              two events, one chromosome
GSM4472423 chr16  Chr16 gain (mat), del 16q23.1-tel (mat)  two events, one chromosome
```

GSM4472403 chr6 is a paternal p-arm deletion plus a paternal q-arm duplication. Two real,
oppositely-signed, non-overlapping-in-biology events that a bare-intersection test on scanned
window coordinates will mark contested and strip of a parent. Note this one is NOT explained by
nullisomy - it is a straightforward geometry bug, and it is independent of whether section 0 is
right.

Both defects point the same way: the flag as written cannot distinguish "two runs describing one
region" from "two runs describing two different regions that happen to abut."

---

## 1. Question 1b: which of the three, and the rule fixed in advance

None of the three as stated. Your option 2 (withhold segment calls on such a cell entirely) is
the one I would most strongly argue against, and it is the one you say is best supported - so
that is the disagreement worth being explicit about.

Option 2 reasons: "if a cell produces four of these, the scan's null is wrong for that cell." On
GSM4552418 the scan's null is not wrong. The cell really is nullisomic almost everywhere; the
scan is reporting that correctly in both runs. Suppressing every segment on that cell would
discard correct findings to protect against a failure mode that cell does not exhibit. It also
contradicts the standing instruction that chaotic cells stay in - not just in letter (you would
keep the cell, dropping only its segments) but in substance, because the segments are the only
per-region content such a cell has.

Option 3 (resolve by intensity) is closer to right but cannot be adopted as stated, for a reason
that is in the code rather than in the epistemics. scanCopyNumber in segments.ts does not run
an independent measurement - it calls scanChromosome on `absent: !m.called`, i.e. on the same
no-call pattern, then filters on a median log2R shift of at least SEGMENT_LRR_SHIFT = 1.0. So:

- it is not independent of dropout. The window-finding step IS a dropout scan; only the filter
  is intensity.
- a shift floor of 1.0 in log2 is a severe requirement on amplified single cells. For a
  one-copy-of-two loss the expectation is about -1.0 exactly; for a both-copy loss the DNA is
  gone and there is no reliable intensity at all.
- and on these cells the question is not "is DNA missing" but "which parent" - intensity
  explicitly informs the state and never the origin (your own dosageOrigin.ts says so in those
  words).

So routing contested windows through intensity answers a different question from the one the
contested flag suppresses.

What I would do instead. Split the flag by what the two runs are actually saying, which is
available at the moment the flag is set and requires no new channel:

```
RULE, fixable in advance, no held-out data needed.

A pair of overlapping opposite-parent parental-absence windows on one chromosome is classified
BEFORE any parent is named, by two tests on the pair itself:

  (i) RECIPROCAL OVERLAP. Let r = (length of intersection) / (length of union), on the refined
      marker-resolution edges, not the scanning window.
        r >= 0.5  -> the two runs describe ONE region.   -> go to (ii)
        r <  0.5  -> the two runs describe TWO regions.  -> NOT contested. Each window keeps
                     its own run's parent. This is a compound event (del one arm, dup/del the
                     other) and is reported as two segments.

  (ii) For a single region called against both parents, the report states BOTH PARENTAL COPIES
       ABSENT OVER THIS REGION, with no parent named, using the existing nullisomy language from
       bothParentsOrigin.ts rather than the current "could not be read" limit line. The region,
       its interval, and both runs' scores are reported.

Threshold justification, so it is not a free parameter chosen later: r >= 0.5 is the standard
reciprocal-overlap criterion for calling two CNV intervals the same event, and 0.5 is its
conventional value. It is NOT tuned here and must not be tuned on these cells.
```

Why this is the honest resolution rather than a relabelling. The current flag conflates two
situations that have opposite correct outputs - and in the more common one, both runs agree and
are both right. "Both copies absent" is a positive finding about the region; "could not be read"
is a confession about the instrument. Printing the second where the first is true understates
what the tool knows, and does so on exactly the chaotic cells the laboratory asked to keep.

What it costs, and the pre-commitment you asked for. You already know which six cells are
affected, so a measurement chosen now is not blind - which is precisely why the rule above is
specified in terms of the hand karyotype's own categories rather than in terms of a result.
Commit this before running:

```
PREDICTION, committed before segment-origin.ts is re-run with the hand label joined on:

  P1. Every contested pair on GSM4552418, GSM4472409, GSM4552378, GSM4552380, GSM4552382 lands
      on a chromosome the hand karyotype calls nullisomic, or on GSM4552422 which the hand calls
      failed throughout.
  P2. Under test (i), the chr6 pair on GSM4472403 and the chr16 pair on GSM4472423 - if the scan
      produces them - separate (r < 0.5) and recover a parent each.
  P3. Detection of hand-called segmental events is unchanged: 11 of 17 paternal, 5 of 10
      maternal. The rule only redistributes parents among already-called regions.

If P1 fails on any pair - a contested pair on a chromosome the hand reads as normal, on a cell
with no nullisomy - then the dropout explanation is live for that pair and section 0 is wrong
about it. Record the failures rather than widening the rule.
```

That is a rule that can be falsified by the run, states what it expects, and does not depend on
knowing which cells are affected - because its content is "classify by reciprocal overlap and
report nullisomy as nullisomy," not "treat these six cells specially."

---

## 2. Question 1: the ladder, and the three failure modes you missed

Your three are all real. 1 (extent error) and 2 (unequal marker density) are the two I would
rank highest. Here are four more.

(4) THE SHARE AND THE RUN ARE NOT INDEPENDENT, so 23 of 23 is weaker evidence than it looks.
Both read the same sample's genotype calls over the same window. The run asks "are this parent's
alleles absent"; the share asks "what fraction of the window's informative markers carry each
parent's allele." A window where the parent's alleles are absent is a window where that parent's
share is low - nearly by construction. Agreement between them is close to a tautology on
clean data, and it is the dropout-driven cases where they could come apart. So 23 of 23 does not
license "these channels cross-validate." It is consistent with "we measured one thing twice."
This is the same in-build/out-of-build distinction the leave-one-out work on a held-out
laboratory sample already
caught in a different part of this codebase: agreement with yourself is not corroboration.

(5) THE CONFLICT BRANCH IS PARTLY UNREACHABLE BY CONSTRUCTION, which is a better explanation of
0 of 23 than "no conflicts exist." defects.ts computes the conflict against `ann`, keyed by
`byWhere.get(where)` where `where` is built from `segmentCoords(sg)` rounded to 0.1 Mb. The
share annotation has to land on the same rounded string as the segment. Any segment whose
annotation key does not match gets `origin = 'unclear'` and never reaches the conflict test. 23
of 69 carrying a share is itself the symptom: two thirds of segments have no annotation to
conflict with. Before concluding the branch never fires, check how many of the 46 share-less
segments are share-less because the key missed rather than because the share was not computed.
That is a one-line instrumentation change, not an experiment.

(6) THE RUN IS NOT ACTUALLY DEFINITIONAL WHEN THE WINDOW IS WRONG - which is your point 1, but
it has a sharper form. "Correct by construction" is a statement about the MARKER SET, not about
the WINDOW. The marker set guarantees: if these markers show absence, it is this parent's
alleles that are absent. It does not guarantee that the window the peak-picker drew is the
region where that absence lives. scanChromosome picks windows at scales w, 2w, 4w... and then
refines edges; refinement can move the boundary but cannot change which scale won. So the run's
parent attaches to a region that may be the wrong region - and on the near-chromosome-length
spans this corpus produces, "wrong region" can mean "mostly normal sequence." The construction
argument is sound for the alleles and silent about the coordinates, and the ladder uses it as
though it covered both.

(7) THE LADDER IS UNTESTED IN THE DIRECTION THAT MATTERS. Your unit test (defects.check.ts:378)
constructs a contested segment and asserts the origin goes unclear. That tests the branch's
plumbing. It does not test the branch's PREMISE - that when the two channels disagree, the run
is the one more likely to be right. Nothing in the corpus speaks to that, because the branch has
never fired on real data. So the ladder's ordering is not merely unvalidated: it is
unfalsifiable on this corpus as currently instrumented.

What would decide it, on 31 cells and nothing else. Not the conflict rate - that is zero and will
stay zero. Decide the PREMISE instead, on the 16 detected hand-called segmental events where the
hand names a parent:

```
For each detected segmental event with a hand-named parent, record both channels independently:
  - the run that found it (parent by construction)
  - the sign of the window allele share (parent by dosage)
and score BOTH against the hand label.

The run's accuracy is currently assumed to be 1.0 and has never been printed. The share's is
quoted as correct on every event the scan reaches. If both read 16 of 16, the ladder's ordering
is untestable here and should be documented as a convention rather than a finding. If the run is
ever wrong where the share is right, the ladder is backwards and that single event is worth more
than the 23 agreements.
```

This is cheap - it is a column in a report you already generate - and it converts "correct by
construction, which is an argument not a number" into a number, which is the gap your own point 3
identifies. Keep the ladder as it is in the meantime; your reasoning for not printing a
contradicted confidence is sound and I would not change the user-facing behaviour on present
evidence. Change the claim, not the code: the ladder is a convention adopted for a stated reason,
not a validated ordering.

---

## 3. Question 2: you are right that a second corpus is needed, and the pre-registration needs
       to fix something other than what you think

### 3.1 Is there a legitimate route to 300/250 on this corpus? No - and I tried to find one

I went looking for a geometry-only argument that reaches window 300 without touching the
held-out results, on the theory that the corner is a multiple-testing artefact. The reasoning
was: the window minimum is not just a floor, it is the base of the scan's scale ladder
(`for (let w = minMarkers; w <= n; w *= 2)`), so lowering it adds scales and multiplies the
number of windows tested. Measured on the real geometry:

```
opposite-hom set, 1,034 markers/chrom     scales   windows tested
  minWindow 2400                            0        0   (scan returns [] : n < minWindow)
  minWindow 1200                            0        0   (same)
  minWindow  600                            1        3
  minWindow  300                            2       13
  minWindow  150                            3       37
```

Two things fall out of that table and both are worth keeping.

First, a real finding: on the opposite-homozygote set, the SHIPPED setting is unreachable. 1,034
< 2,400, and scanChromosome opens with `if (n < minMarkers) return []`. The grid rows at 2400 and
1200 do not represent a conservative version of the new statistic - they represent no scan at
all. Any row of your grid at those settings is structurally empty, and the comparison "shipped at
2400/250 vs opposite at 150/80" is not two settings of one statistic, it is the only reachable
region of the new statistic against the only reachable region of the old one. That is worth a
line in the record independent of everything else.

Second, the honest negative: the multiple-testing argument does NOT get you to 300. I computed
the exact binomial tail for the scan's own likelihood-ratio statistic, summed over the real
window counts above, at null absence rates of 1, 2 and 5 percent. Expected false windows per
chromosome comes out around 1e-110 at 250 and 1e-35 at 80. Both are so far below one that the
2.85x difference in windows tested between 150 and 300 (37 against 13) is irrelevant. The grid's
corner is not a multiple-testing failure, and a family-wise-error cap would not have held it back.

That negative is informative: it means the 12 of 224 and 13 of 224 calls on hand-normal
chromosomes are not the binomial null leaking. At 1e-110 per window the binomial null leaks
nothing. Those calls come from the null being MISSPECIFIED - real absence structure in the
sample that the per-chromosome external null does not model (localised dropout, amplification
bias, mis-clustering). Which means no threshold chosen by arithmetic on the null can calibrate
this statistic. Your instinct that the operating point has to come from data is not conservatism;
it is forced by the fact that the only tractable null is the wrong one.

### 3.1b CORRECTION TO MY OWN ANSWER: there IS a legitimate route, and it is in your own records

I wrote "no legitimate route" above and then found one while checking this project's prior
measurements. I am leaving the wrong answer in place above rather than deleting it, because the
route only makes sense against the failed attempt.

You asked to be contradicted if you were wrong. You are, but not about the thing you expected -
and the correction splits your question into two claims that you have been treating as one.

Two results from this project's own prior work, both predating this pre-registration and both
measured on DIFFERENT material (dilution series, 32 replicates/cell, 12 clean haploid genomes,
0 false segments):

  - the segment floor is CONTRAST-dependent, not a single marker count
  - marker floors should be derived as N = ln(R) / KL for required posterior odds R, not
    tabulated

The scan's LRT threshold is already in log units, so ln(R) = T exactly.

FIRST, AN ORIENTATION PROBLEM THAT CHANGES THE NUMBERS. Bernoulli KL is asymmetric, and the
project's recorded rule is written KL(h0 || h1) while the quantity that governs whether a true
event clears the threshold is KL(h1 || h0). The scan's score is
logLike(k, w, observed) - logLike(k, w, null); under a true event of in-window rate h1 its
expectation is w * KL(h1 || h0). So KL(h1 || h0) is the correct orientation FOR THIS PURPOSE and
the recorded rule's order is either a different intended use or a transcription slip. I cannot
tell which from the record, so both are tabulated below and the weaker one governs the
conclusion. h0 = 0.003 throughout:

```
                              KL(h1||h0)              KL(h0||h1)  [as recorded]
dropout    h1       KL     N at T=250      KL     N at T=250
    0%  0.999    5.795           43      6.867           36
   30%  0.700    3.456           72      1.181          212
   50%  0.500    2.213          113      0.673          372
   70%  0.300    1.134          221      0.339          738

one-parent-hom, full contrast (prior work): in-segment 0.106 vs 0.003 outside
                        KL(h1||h0) = 0.280 -> N at T=250 = 892
```

Window 300 clears the floor under BOTH orientations at low-to-moderate dropout, and the two
diverge where it matters:

```
                       300 markers suffice while dropout stays below
  KL(h1||h0)                              76%
  KL(h0||h1)  [as recorded]               43%
```

So the claim has to be scoped rather than stated flat, and I had it too strong on first writing:
window 300 at LRT 250 is derivable in advance **for material up to roughly 43 percent dropout**,
which is the conservative reading and the one I would use. Beyond that the recorded orientation
requires 372 markers at 50 percent and 738 at 70 percent, and window 300 does NOT clear it. Since
this corpus explicitly includes cells whose stage reads `failed` and pronuclei documented at 53.8
to 59.1 percent call rate, the high-dropout regime is populated here, so the scope limit is real
and not hypothetical.

What survives regardless of orientation: the opposite-homozygote set needs FEWER markers than the
diluted set at equal dropout, because full contrast is its defining property - under either order
its floor is far below the 892 the diluted set needs. And window 150 fails the floor under both
orientations at 70 percent dropout (221 and 738 against 150), which is a second reason the corner
was wrong that owes nothing to the held-out scores.

So window 300 at LRT 250 is derivable, in advance, from measurements that predate this
pre-registration and were taken on different material. That is an external justification, not a
read off the held-out half. The grid was the wrong instrument: a floor that your own prior work
says must be derived from contrast was instead chosen by searching a tabulated grid, which is
exactly the practice that work concluded against.

WHAT THIS DOES AND DOES NOT LICENSE - and this is the distinction your question collapses:

  - It licenses ADOPTING window 300 / LRT 250 as the pre-specified operating point, derived
    rather than searched. No second corpus needed for that.
  - It does NOT license CLAIMING the (75 of 102, 5 of 224) performance. Those two numbers are
    held-out observations and remain held-out observations whatever justifies the window.

Your document treats "adopt the point" and "claim its numbers" as the same act, which is why the
only route you could see was a second corpus. They are separable. You may adopt a point on prior
derivation and then report its performance as what it is: a single read on a half that has now
been seen, quoted as the reason for the change but not as a validated operating characteristic.
The second corpus is still needed - for the CLAIM, not for the CHOICE.

Three honest limits on this, since it is now load-bearing:

  1. The 0.003 null is the prior work's measured outside-segment rate, not the scan's actual
     externalNull, which is this cell's own per-chromosome median and runs higher on chaotic
     cells. Recompute the KL floor with the real per-chromosome nulls before fixing the number.
  2. The KL floor is a DETECTION-POWER floor: it answers "can a true event clear the threshold",
     not "will false windows stay out". The specificity gap (5 vs 12) is a null-misspecification
     problem and the KL argument says nothing about it. So this route justifies the window on
     power grounds and leaves the specificity claim exactly where it was.
  3. Deriving a floor (221 or 738 at 70 percent dropout, depending on orientation) and then
     choosing 300 because 300 is the grid point nearby is still a grid artefact. Derive the
     number and use the derived number, stated as a fraction of informative markers per
     chromosome, rather than snapping to the grid.

One route I should name and reject explicitly, since it looks legitimate: re-splitting the 31
cells (different parity, randomisation, k-fold) to get a better-balanced tuning half. It is not
legitimate. You have now seen the held-out results; any re-split is chosen in their light, and
the 102-vs-50 imbalance that motivates it is itself a held-out observation. Record the imbalance
as a design defect - which you did - and leave the split alone.

### 3.2 The honest intermediate, before a second corpus arrives

There is one thing the 31 cells can still support, and it is not an operating point. The
dominating point's advantage is concentrated in specificity: 5 of 224 against 12 of 224. If the
7-chromosome difference is mostly the 'failed'/chaotic cells, then a CELL-LEVEL GATE - not a
threshold change - captures most of the gain, and a gate is defensible on grounds that predate
the held-out scoring (you already refuse arrays on call-rate and BAF-extremity grounds
elsewhere; origin.py documents the mis-clustering gate). Ask of the shipped statistic's 12:

```
  how many of the 12 hand-normal calls fall on cells whose stage reads `failed`, or whose
  BAF-extremity or call rate would already refuse the array on existing gates?
```

If most do, the specificity gap is a cell-admission question rather than a statistic question,
and tightening admission is a change you can justify from the existing gate rationale rather
than from the held-out curve. If few do, that is also worth knowing: it says the two statistics
differ on genuinely clean material and the second corpus is the only arbiter. Either answer is
useful and neither requires adopting a point read off held-out data.

### 3.3 What the GSE60909 pre-registration must say

Your question - what caps a specificity constraint when the only reference rate is the
incumbent's own, measured on a loose half - has a direct answer: nothing, and that is why the cap
must not be a rate at all.

The defect in the first rule is not that the cap was too loose. It is that the rule was
LEXICOGRAPHIC: maximise sensitivity, subject to a cap. A lexicographic rule always sits on its
constraint boundary, and when the constraint is slack it walks to the grid corner. Re-writing it
with a tighter cap just moves the corner. The fix is to stop optimising one axis subject to the
other.

```
WHAT THE GSE60909 PRE-REGISTRATION SHOULD FIX

(a) THE SELECTION RULE IS NOT LEXICOGRAPHIC. Score every grid point by a single objective that
    prices both axes, and name the exchange rate in advance:

      maximise   (hand-called deletions detected) - k * (calls on hand-normal chromosomes)

    with k committed before scoring. k is the number of false region-calls you would trade for
    one true detection, and it is a clinical judgement rather than a statistic: for a tool whose
    output is "whose copy was lost" on material heading toward embryo work, a wrong region named
    with a parent is more costly than a missed region, because the missed one is silent and the
    wrong one is actionable. So k should be well above 1. Worked on this corpus's own numbers:

        k      shipped (76,12)   150/80 (85,13)   300/250 (75,5)   winner
        0.50        70.00            78.50           72.50         150/80
        1.00        64.00            72.00           70.00         150/80
        1.25        61.00            68.75           68.75         tie
        1.50        58.00            65.50           67.50         300/250
        2.00        52.00            59.00           65.00         300/250
        3.00        40.00            46.00           60.00         300/250

    Note what this does NOT say, correcting an arithmetic error in my first draft: k = 1 does not
    reach the dominating point. At k = 1 the corner still wins, 72 against 70. The breakeven
    between 300/250 and 150/80 is exactly k = 1.25, and 300/250 only wins outright above it. So
    k must be committed at 1.5 or higher for this objective to select the point I am arguing for,
    and a k chosen at 1 would have reproduced the corner under a new name. Commit k >= 1.5 on the
    stated clinical asymmetry, or acknowledge that the objective does not by itself avoid the
    corner. (Against the shipped statistic alone the breakeven is k = 0.14, so almost any k
    prefers 300/250 to shipped; the contest is entirely between 300/250 and the corner.)

(b) TIES AND NEAR-TIES RESOLVE TOWARD THE LARGER WINDOW, which is already your tie rule - but
    extend it: among points within one detection of the optimum, take the largest window. This
    is the explicit anti-corner clause the first rule was missing.

(c) THE WINDOW MINIMUM IS DECLARED AS A FRACTION OF A CHROMOSOME'S INFORMATIVE MARKERS, NOT AS A
    MARKER COUNT. This is the cross-platform defect, and it is the one that would silently wreck
    the transfer. HumanCytoSNP-12 carries 301,232 markers against Axiom's 825,656 (2.74x fewer),
    so a window fixed in markers is a 2.74x larger REGION on the second platform:

        window (markers)   span on Axiom   span on HumanCytoSNP-12
                    2400          9.0 Mb                  24.7 Mb
                     600          2.3 Mb                   6.2 Mb
                     300          1.1 Mb                   3.1 Mb
                     150          0.6 Mb                   1.5 Mb

    Worse, the opposite-homozygote set scales with array density too - roughly 377 markers per
    chromosome on CytoSNP against 1,034 on Axiom - so the grid's reachability changes: at
    minWindow 600 the scan returns [] on the second platform, and at 300 it tests 2 windows per
    chromosome against 13. The same nominal setting is not the same test. Declare the window as
    "29 percent of a chromosome's informative markers" (which is what 300 means on Axiom) and
    let the marker count follow from the platform, or the comparison measures array density
    rather than statistic quality.

(d) THE SPECIFICITY REFERENCE IS THE INCUMBENT MEASURED ON THE SAME CELLS, scored in the same
    run, never a rate carried over from GSE148488. Platform, laboratory, amplification and
    answer-key conventions all differ; a cross-corpus rate comparison would confound all four.
    The relative use you already defend is the only sound one - see (e).

(e) THE SPLIT IS BY EMBRYO, NOT BY CELL. This corpus taught that lesson and it should be carried
    forward explicitly. 14 of the 31 cells are sibling blastomeres of 4 embryos, and siblings
    share both their zygote and their amplification lineage:

        5 cells  MII injection single 'wt' blastomere embryo   (chr6 event in 4 of them)
        4 cells  MII injection blastomere embryo B             (del12p (pat) in ALL FOUR)
        3 cells  2-cell injection, blastomere embryo one
        2 cells  2-cell injection, blastomere embryo

    All four `del12p (pat)*` rows are one embryo. Splitting by cell puts correlated copies of one
    event on both sides of the split, which inflates apparent agreement and is part of why the
    halves came out 50-vs-102. Split by embryo, and state the expected imbalance up front rather
    than discovering it in a secondary endpoint.

(f) THE PRIMARY ENDPOINT IS PRE-DECLARED AS THE OBJECTIVE IN (a), SINGLE NUMBER, ONE READ. No
    bar-with-two-conditions, because a two-condition bar is what produced "met one, missed the
    other by one chromosome" - an outcome with no decision attached to it.

(g) THE WINDOW MINIMUM IS DERIVED, NOT SEARCHED. Per 3.1b and this project's own prior finding
    that marker floors should come from N = ln(R)/KL rather than a table - fixing the KL
    orientation explicitly, since the two orders differ threefold at high dropout: compute the
    floor from the contrast the opposite-homozygote set actually delivers on CytoSNP material,
    at the LRT you intend to ship, at the worst dropout you will admit. Put that single derived
    number in the pre-registration and do not grid-search the window at all. A grid over a
    quantity your own method says is derivable is how the corner happened; removing the axis
    removes the corner more reliably than any cap on it. Grid only the LRT, if anything.
```

### 3.4 Question 2.3: is the hand-normal measure sound?

Yes, with one qualification you have not stated, and there is a better measure available.

Sound: as a paired relative measure between two statistics scored on the same chromosomes of the
same cells, it is valid, and your refusal to call it an error rate is correct. The authors
tabulated the losses they were looking for, so untabulated real events land in that column.

The qualification: it is only valid if the two statistics are scored on the SAME SET of
chromosomes. They are not, necessarily. The shipped statistic at 2400 cannot scan an
opposite-homozygote chromosome at all (1,034 < 2,400), and both statistics skip chromosomes
excluded by `measured`, `lostByAlleles` and the uninterpretable-array gate. If the two arms
scored different denominators, 12 vs 13 is not a paired comparison. Report the scored
denominator per arm, not just "of 224."

The better measure, which exists on public data and in your own repo: the five paternal pronuclei
of GSE148488 (GSM4774680/681/682/683/685, already bundled as germlineExamples.ts). Each is a
single meiotic product of the sperm donor, so every autosome carries exactly one parental
complement BY DISSECTION rather than by hand reading - and critically, these arrays are not among
the 31 hand-karyotyped cells (I checked: zero overlap). On that material, "the other parent's
alleles are absent genome-wide" is ground truth from the dissection, and any SEGMENTAL call
against the dissection-known present parent is a false positive with no "untabulated event"
escape hatch. That converts your relative measure into an absolute false-positive rate on
material whose truth does not come from the same by-eye reading being tested. Your own
citations.ts already makes this argument for the absence-detection step; it extends to the
segment scan unchanged. It is also a cleaner specificity reference than the incumbent's own rate,
which answers the cap question in 3.3 from a different direction: measure the floor on material
where the truth is structural.

Caveat to keep attached: a dissected pronucleus has been through neither syngamy nor a mitosis,
so its dropout profile is not a blastomere's, and a false-positive rate measured there is a floor
rather than an estimate. Your citations.ts already says exactly this about sensitivity; the same
discount applies here.

---

## 4. Question 3, briefly, since you did not rank it

All three of your candidate explanations are testable on this corpus without choosing the test
after seeing the answer, because the truth table stratifies them in advance. Specifically: 5 of
the 27 hand-called segmental events are DUPLICATION-ONLY, and an absence scan structurally cannot
see a copy that is present twice - your own pre-registration says so. So the denominator for
"what the absence scan missed" is 22, not 27, and those five belong in a different channel's
accounting:

```
GSM4552380 chr16  dup 16q23.1-tel (mat)*
GSM4552382 chr6   dup(6q)(pat)
GSM4472422 chr6   dup (6q)(pat)
GSM4552365 chr16  dup 16q21-tel (pat) FRA16I
GSM4472403 chr18  dup18p (mat) FRA18E
```

Pre-specifiable strata, all readable off the table's own text before any scan is run: deletion-
bearing (22) vs duplication-only (5); terminal ...-tel (10) vs whole-arm (15); FRA fragile-site
labelled (3); asterisked (10). Cross the 11 misses against those four strata and the three
explanations separate: concentration in duplication-only or terminal events means "below the
window minimum / wrong detector"; concentration in the chaotic and failed cells means "those
cells are too noisy"; concentration in the FRA and asterisked rows means the hand karyotype is
reading cytogenetic context the array cannot carry. None of those strata were chosen by looking
at which events missed.

---

## 5. What I would change, in order

1. Fix the contested flag's overlap test: reciprocal overlap r >= 0.5, not bare intersection.
   This is a bug with a known false-positive case in the corpus (GSM4472403 chr6) and is
   independent of everything else here.
2. Report both-copies-absent AS nullisomy on single regions called by both runs, using the
   language bothParentsOrigin.ts already has, instead of "could not be read." Commit the P1-P3
   predictions first, then run segment-origin.ts with the hand label joined on.
3. Record in the pre-registration that the opposite-homozygote statistic is UNREACHABLE at
   windows 2400 and 1200 (1,034 < minWindow, scan returns []), so those grid rows were empty
   rather than conservative.
4. Print the run channel's own accuracy against the hand label beside the share's, on the 16
   detected events. Until then, describe the ladder as a convention with a stated rationale, not
   a validated ordering.
5. Instrument how many of the 46 share-less segments are share-less because the 0.1 Mb
   annotation key missed, before concluding the conflict branch has never had the chance to fire.
6. Re-denominate segmental detection on 22 deletion-bearing events rather than 27.
7. Re-derive the window minimum from contrast, per 3.1b, using this project's own N = ln(R)/KL
   rule with the KL orientation resolved, the scan's real per-chromosome nulls, and the worst
   dropout you admit. If it lands near 300 on the opposite-homozygote set, you may adopt that
   point now on
   prior-derivation grounds - while reporting its held-out numbers as one read rather than as a
   validated operating characteristic. The second corpus is for the claim, not the choice.
8. Leave the shipped statistic's PERFORMANCE CLAIM and the withheld-parent user-facing behaviour
   alone. What does not survive this review is the claim that 1b was a choice between three
   options, the belief that the grid corner was a tuning accident rather than a lexicographic
   rule walking to its boundary, and - on my own second look - the conclusion that a second
   corpus is the only route to the dominating point.

## 6. Where I am least certain

- Section 0's pair-by-pair identity. The cell counts match exactly (6 and 6) and GSM4552418's
  karyotype explains the span pattern, but I could not run the scan. P1 is the falsifier.
- Whether the share and the run can come apart at all (point 4). I argue they are near-tautologous
  on clean data; I did not measure their correlation, which would need the arrays.
- 3.1b's KL floor, on two counts. First, it uses 0.003 as the outside-segment null, which is the
  prior dilution work's measured figure on different material, not the scan's live externalNull;
  on chaotic cells the real null is higher and the floor rises. Second, the project's recorded
  rule writes KL(h0||h1) whereas the scan's own score has expectation w*KL(h1||h0), and the two
  disagree by a factor of 3 at high dropout (221 vs 738 markers at 70 percent). I argue the
  second orientation is the right one for a detection floor and have taken the first as the
  conservative bound, but I cannot tell from the record which the prior work intended, and that
  single choice sets whether window 300 clears at 43 percent dropout or at 76 percent. Resolve
  the orientation against the original dilution measurement before this becomes load-bearing.
  The direction of the argument is robust under both (full contrast needs far fewer markers than
  the diluted set); the specific thresholds are not.
- k in 3.3(a). The exchange rate is a judgement about clinical cost, not something I can measure.
  I have argued for small k and named why; the number is yours.
- The CytoSNP opposite-homozygote count (377/chrom) is scaled from the Axiom figure by array
  density. It is an estimate for planning the pre-registration, not a measurement; compute it
  from the actual arrays before fixing the fraction in (c).

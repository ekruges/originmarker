# Response 3: the reporting layer, and the defect underneath the questions it is asking

What this rests on. Repo read 2026-10-26: `web/src/segments.ts`, `parentage.ts`, `obligateHet.ts`,
`scoreSample.ts`, `stage.ts`, `ingest.ts`. Node is present (v26.8.1); GSE148488 is still not on
this machine, so nothing was re-run. Every detection figure is quoted from your documents.
Everything asserted as new is derived from the shipped source or from algebra, and the algebra is
shown so it can be checked rather than trusted.

Headline, and it reframes the design request. **Your section 4 is close to right, and the reason
it is right is not the reason you give.** More importantly, four of your seven questions - 5.2,
5.6, 5.7, and the disqualification half of 5.1 - are the same defect seen from four sides, and it
is not a reporting problem. It is in `externalNull`. The reporting layer should be designed after
that is fixed, because the fix changes what there is to report.

I will take 5.1 first since you asked for it as a derivation, then show the defect, then answer the
rest against it.

---

## 1. Question 5.1: the partition, derived

The classifying question is: **what is in the denominator?** Four reference classes, and only the
first two survive a genome with no normal majority.

```
CLASS                       survives?   why
this chromosome's own markers   YES     a ratio within one chromosome needs no outside baseline
the PARENTS' genotypes          YES     the parents are bulk; their arrays are not the sick one
this sample's other chromosomes NO      a chaotic genome has no normal majority to normalise to
this sample as a whole          PARTLY  valid as an ARRAY-level statement, never per chromosome
```

Applied to everything the pipeline computes:

```
quantity                                     reference class          reportable
informative markers (per chrom, per parent)   PARENT only                 yes
absence count and rate (per chrom, per parent) PARENT only                yes
call rate (per chrom)                         THIS CHROM only             yes
no-call fraction (per window)                 THIS CHROM only             yes
het share on the opposite-hom set             THIS CHROM + PARENTS        yes
onlyA / onlyB  (survivor share)               THIS CHROM + PARENTS        yes
contribution (uniparental / biparental)       THIS CHROM + PARENTS        yes
bafExtreme (share outside 0.15-0.85)          THIS CHROM only             yes
nonParental rate                              PARENT only                 yes
aneuploidy from CALL_COLLAPSE                 THIS CHROM only             yes
hetRate, callRate, hetBafSd, homInBand        WHOLE SAMPLE                array-level only
------------------------------------------------------------------------------------------
hetSpread (CV across chromosomes)             CHROM vs OTHER CHROMS       no
mosaicZ                                       CHROM vs OTHER CHROMS       no
externalNull                                  CHROM vs OTHER CHROMS       no
segment LRT score                             CHROM vs OTHER CHROMS       no
log2R median, genomeLrr                       CHROM vs GENOME             no
copy-number segment (LRR shift)               CHROM vs GENOME             no
```

**What the partition puts on the reportable side that your section 4 omits.** Your table proposes
four columns: informative markers, absence vs father, absence vs mother, call rate. The partition
adds five more, all already computed, none requiring a new measurement to *print*:

1. **`onlyA` / `onlyB`, and the survivor share derived from them.** This is the answer to 5.2. See
   section 3.
2. **Het share on the opposite-homozygote set.** On that set a biparental diploid is obligate
   heterozygous, so this reads one-versus-two contributions directly, per chromosome, with the
   parents as the only reference.
3. **`bafExtreme`.** Print it because it is the gate's own input. A chromosome suppressed as
   `not_measured` should show the number that suppressed it; otherwise the reader cannot see which
   chromosomes the tool declined to measure versus measured and found nothing.
4. **Aneuploidy from `CALL_COLLAPSE`.** A within-chromosome call-rate collapse references nothing
   outside that chromosome. It is the one aneuploidy channel that survives.
5. **`nonParental` rate.** An allele neither parent carries. Parent-referenced, and it is the
   measure that separates a contaminated reaction from a chaotic genome - which bears on 5.3.

**What the partition disqualifies that you assumed safe.** Nothing in your four columns, but the
boundary is in a different place than you drew it. You wrote that the per-chromosome absence rate
qualifies "because it is a ratio within one chromosome" - that is true of the *rate*, and the
reason is that the markers are defined by the parent, not that the ratio is internal. But the
**segment** built from that rate does not qualify, because `scanChromosome` scores it against
`externalNull`, a median over the sample's other chromosomes (`segments.ts:233-244`).

So your section-4 boundary - numbers ship, segments do not - is correct. The reason you gave is
not the reason it holds, and the correct reason has a consequence you have not drawn: **the
segments on the ten chaotic cells you already analyse rest on that same disqualified reference.**
That is 5.6, and it is not a judgement call.

A caution on your phrasing. You say "she is right that the allele channel needs no baseline". On
the partition she is right about the *rate* and about `obligateHet`, and wrong about the
*segments* - and the `clustered` dependency I traced last round is a third, separate route by
which intensity enters. Three different things have been called "the allele channel" in this
thread; the partition is worth keeping precisely because it stops that conflation.

---

## 2. The defect: `externalNull` breaks down on exactly these cells

`segments.ts:233` computes the null as the **median per-chromosome absence rate over the cell's
other autosomes**. Its own comment says why a median: *"Taking the median rather than the mean is
what stops several events from doing it collectively."*

A median tolerates a minority of contaminated values. Its breakdown point is 50 percent. Your
chaotic cells carry **12 to 15 abnormal autosomes of 22** - 55 to 68 percent. On exactly these
cells the median *is* an affected chromosome, so the null is set by the rate real events produce.

And the second half, which I had wrong last round and which makes this decisive: **a real
single-parent loss does not drive absence to 1.0.** On markers where the mother is homozygous and
her copy is lost, the child carries only the paternal allele - which is the *opposite* allele only
at the fraction of markers where the father transmitted one. Call that h:

```
h = P(child truly heterozygous | this parent homozygous)
  computed over an allele-frequency spectrum calibrated to the panel's measured 0.168 bulk het:
h = 0.100
```

So in-event absence is about 0.100, not 1.0. Last round I tabulated the LRT "even at 10% in-event
absence" and called it the pessimistic corner. **It is the expectation.** That was an error and it
hid this.

Now put the two together. The scan's first test is `if (rate <= nullRate) continue`
(`segments.ts:274`) - before the LRT is computed at all.

```
null (median of a chaotic cell)   event rate   rate > null?   LRT at w=2400
                          0.003        0.100        yes             620.5
                          0.020        0.100        yes             202.3
                          0.050        0.100        yes              49.6
                          0.080        0.100        yes               6.1
                          0.100        0.100        NO      rejected before scoring
                          0.150        0.100        NO      rejected before scoring
```

On a cell whose median absence rate reaches h, **every window on every chromosome is rejected
before it is scored.** The marker count is irrelevant - which is why chr16 with 14,322 to 14,991
maternal-homozygous markers, six times the floor, finds nothing. That is 5.7, and it is not a
sensitivity problem, not a thin marker set, and not the window floor. My previous first-choice
explanation was wrong and so was my second.

It also explains the contested windows with no extra assumption. On a cell whose bar is already
at event level, the windows that *do* clear are the ones where dropout pushes absence above it -
and dropout is symmetric between parents, so it clears **both runs at once**. Observed signature:
35 of 69 contested, all on chaotic cells, zero on clean cells. A clean cell's median is ~0.003, so
its bar is low and real asymmetric events clear it easily.

One figure for both halves:

![The segment null breaks down on chaotic cells; the discriminator already exists](art_e805ead0-2476-413d-9c29-c74194e50ca6)

Panel a is a schematic - a two-point idealisation where affected chromosomes all sit at h and the
rest at a clean 0.003. Real per-chromosome rates are a distribution and the transition will be a
ramp rather than a step. The breakdown point and its direction are the claim; the sharpness is not.

**Predictions, falsifiable, one instrumented run, no new data.**

```
P1  On the three chr16 cells, externalNull for the MATERNAL run is >= ~0.10, i.e. at or
    above the rate a true maternal loss produces.
P2  On those cells the 16q23.1-tel window's maternal absence rate is ~h and fails
    `rate <= nullRate`, so it is rejected before the LRT.
P3  On the six contested cells, both runs' nullRate values are elevated, and the contested
    windows' two absence rates are approximately EQUAL (survivor share near 0.5), not both
    near h.
FALSIFIER: if P3 fails and the contested windows are strongly asymmetric, the dropout
    explanation is wrong and so is the 5.2 answer below.
```

---

## 3. Question 5.2: option 3, and the discriminator is already in your codebase

You say you cannot find a signature in the allele data alone. There is one, it is exact, it needs
no intensity and no new free parameter, and `obligateHet.ts` already computes it.

`addTwoParent` works on the **opposite-homozygote** set - father AA and mother BB - where a
biparental diploid is **obligate heterozygous**. It tallies `het`, `onlyA`, `onlyB` per chromosome.
With per-allele dropout at rate d:

```
                    probability        conditioned on being called
called AB           (1-d)^2            (1-d)/(1+d)
called AA (only A)  d(1-d)             d/(1+d)        father's allele alone
called BB (only B)  d(1-d)             d/(1+d)        mother's allele alone
no call             d^2
```

`onlyA` and `onlyB` are **equal in expectation at every d**. They have to be: the array has no way
to know which allele came from which parent, so dropout cannot be parent-specific. Therefore

```
survivor share = max(onlyA, onlyB) / (onlyA + onlyB)

   symmetric dropout, ANY d        -> 0.5
   real single-parent loss         -> 1.0
```

And the gap is **already measured on this corpus**: `SURVIVOR_MIN = 0.70`, from 595 autosomes of
GSE148488, where 105 hand-called lost chromosomes never read below 0.735 and 490 intact ones never
reached 0.627. `obligateHet.ts` even states the mechanism in words: *"Dropout cannot cross it,
because dropout removes a heterozygote's alleles at random rather than one parent's allele at
every marker of a chromosome."*

So the answer to 5.2 is your option 3, and the third quantity exists, is self-normalising within
one chromosome, and has a threshold sitting in a measured gap. Why it is not already doing this
work: the segment scan uses the **one-parent-homozygous** set with a cross-chromosome null, while
`obligateHet` uses the **opposite-homozygote** set with a within-chromosome reference. The second
survives a chaotic genome; the first does not.

**The limit, and it is a real one.** `SURVIVOR_MIN` is calibrated *per chromosome*, with
`MIN_INFORMATIVE = 200`. The opposite-homozygote set is roughly 1/16th the size of the
one-parent-homozygous set, so chr16's ~14,650 one-parent markers correspond to only about 947
opposite-homozygote markers, and the 16q23.1-tel window to about **171** - below the 200 floor.

```
chr16 whole                  ~947 opposite-hom markers   -> test available
chr16 16q23.1-tel (18.0%)    ~171                        -> BELOW MIN_INFORMATIVE
```

So the survivor share can answer "did chr16 lose a maternal copy anywhere" but **not** "is the
terminal 16.3 Mb the part that went". That is a genuine resolution limit of the allele channel,
it should be stated in the layer rather than worked around, and re-deriving the floor at window
scale is a measurement to pre-register, not an assumption to ship. (The 1/16 ratio is scaled from
the project's own median figures, so the 171 is an estimate; count it directly before relying on
it.)

**How to present it.** Not your option 1, and not quite your option 2:

```
Per chromosome, print: informative (both sets), call rate, absence vs father, absence vs
mother, onlyA, onlyB, survivor share, bafExtreme.

Attach ONE of three named states, from the survivor share and the call rate alone:

  survivor share >= SURVIVOR_MIN, call rate not collapsed
      -> "one parent's allele alone over this chromosome" and NAME the parent.
         This is a parent call and it is defensible on a chaotic genome, because
         nothing in it references another chromosome.

  survivor share ~ 0.5, both absence rates elevated
      -> "both parents' alleles depleted; consistent with dropout at this call rate"
         No parent, no region. This is the 35-of-69 case and it is correctly a
         non-finding rather than a finding about a region.

  call rate collapsed, both denominators thin
      -> "too few called markers to measure" and print how few.
```

The middle state is your option 2 minus a claim I would not make. You proposed naming it "both
parents depleted over this window", which you call "a positive statement about the region". It is
not positive - symmetric depletion is what dropout *predicts*, so the honest statement names the
artefact, not the region. Naming it as a property of the region is the one thing in your design I
would strike.

---

## 4. Question 5.6: yes, and it is not a judgement call

You least wanted to ask it. The answer is that the ten chaotic cells should stop emitting segment
verdicts, and the reason is section 1's partition, not a contested-fraction threshold.

Those segments are scored against `externalNull`, which the partition disqualifies, and section 2
shows the disqualification is not theoretical on these cells specifically - it is a median past
its breakdown point. A verdict computed against a reference the sample cannot supply is not a
verdict with a wide error bar; it is a verdict about nothing.

So: not "drop to measurements-only above some contested fraction" - that would be a new free
parameter in a gap you have not measured, and your own constraint forbids it. The rule follows
from the partition:

```
Segment verdicts require a null the sample can supply. Emit them only where a majority of
autosomes are unaffected. Where they are not, emit the per-chromosome measurements and the
survivor-share state instead.
```

"A majority of autosomes unaffected" is not a tuned threshold - it is the median's breakdown
point, which is 50 percent by construction. That is a derivation, not a gap midpoint, and it
satisfies your constraint 2 without a new measurement.

Note what this does **not** do: it does not reduce what the lab sees on those cells. It replaces
69 segments, a third self-contradictory, with per-chromosome measurements plus parent calls
wherever the survivor share clears 0.70 - which is the quantity the lab actually described in the
earlier call ("if it's the mother's SNPs that are missing, then it's the mother's chromosome that
was lost"). She described `obligateHet`, not the segment scan.

---

## 5. Question 5.4: a real risk, and the mitigation is already in your code

Yes, it is real, and your framing of why is right: the intended reader is precisely the person who
will make the call the tool declined.

But the hazard is not "numbers without a verdict". It is **numbers without the reason the verdict
was withheld**. A reader who sees two elevated absence rates and no verdict will infer the tool was
being cautious. A reader who sees two elevated absence rates, a survivor share of 0.52, and the
line "symmetric depletion is what dropout produces at this call rate" has been told the actual
finding, which is that there is no finding.

So the mitigation is: **every suppressed verdict prints the quantity that suppressed it.** Not a
caveat - the number. `bafExtreme` for a `not_measured` chromosome; the survivor share for a
both-parents case; the call rate for a thin one. The reader can then overrule the tool *on the
evidence*, which is what the lab asked for, rather than in the absence of it.

`scoreSample.ts:552` already does exactly this for the mixed-signal case: it logs the counts, says
it was not called, and names no parent. That is the precedent, it is shipped, and the taxonomy row
you quoted - *"Report the ORIGIN and withhold the CLASS as separate fields"* - is the same
principle. The design should generalise that line rather than invent a pattern.

One addition. Order the table so the suppressing quantity sits **beside** the suppressed field,
not in a footnote. A reader scanning for "which parent" should hit "survivor share 0.52" in the
adjacent column, not a legend.

---

## 6. Question 5.5: three tests that score something, none requiring the laboratory

"The numbers are correct" is indeed a unit test. But a measurements layer makes three claims that
*are* scoreable, and you have truth for two of them.

**Test 1, the layer's states against the hand karyotype.** The survivor-share state is a call -
parent-named, both-depleted, or unmeasurable - so it scores exactly like everything else in the
project. Pre-register: on the 595 autosomes with hand karyotypes, what fraction of hand-called
single-parent losses get a parent named, what fraction of hand-normal chromosomes get one wrongly,
and what fraction land in each non-call state. This is the layer's primary validation and it needs
no new data.

**Test 2, negative controls that must read "no finding".** The dissected pronuclei. A haploid
product is uniparental by dissection, so its survivor share is 1.0 everywhere **by construction** -
which makes it a positive control for the parent-naming state. And the 5 bulk diploid arrays are
the negative control: every autosome should read biparental, survivor share near 0.5, with no
parent named. You already use both sets this way for other statistics; nothing new is needed.

**Test 3, the one that tests "not misleading" rather than "correct".** Pre-register the prediction
that **no state transition is possible without the printed number moving across its stated
threshold.** That is a property of the display, scoreable mechanically over the corpus: for every
chromosome, recompute the state from the printed columns alone and check it reproduces the state
shown. If a reader cannot re-derive the verdict from the numbers beside it, the layer is
misleading by construction, regardless of whether the numbers are right. This is the test your
question asks for and it does not require showing anything to anyone.

What none of these test is whether the lab finds it *useful*. That genuinely requires asking her,
and I would not pretend otherwise - but "useful" was never the thing at risk. "Misleading" was,
and Test 3 scores it.

---

## 7. Question 5.3: no sound discriminator for three of the four, and what to print

Direct answer: you have one real discriminator and cannot get the others from public data.

What separates what:

```
chaotic genome        call-rate-conditioned spread. Validated: 0 of 5 bulk diploid, 0 of 15
                      pronuclei read localised. This one works.
contaminated reaction the nonParental rate - alleles NEITHER parent carries. A contaminant
                      contributes alleles from a third individual; a chaotic genome cannot
                      invent an allele. This is parent-referenced, so it survives the
                      partition, and it is already computed (`tallyRow`, nonParentalDen).
                      I do not think you have scored it, and it is the one addition to 5.3
                      that is free.
two individuals mixed BAF modality. A 50:50 mixture is trimodal at 0.25/0.50/0.75. Your
                      46-percent-outside-the-bands figure is consistent with it but does not
                      separate it from a destroyed array, because both smear.
destroyed array       no positive discriminator available.
```

So: two of four are separable with what you have (one validated, one free and unscored), and mixed
versus destroyed is not separable on public data.

**What to print for a cause that cannot be narrowed.** Not a list of three possibilities with the
measurement that would narrow it - that is better than today's lumped sentence, but it asks the
reader to run an experiment the repository cannot specify. Print instead:

```
what was measured      the value and the threshold it crossed
what it rules OUT      "this is not a genome reading this way" - the strong, certain claim
what remains           the named causes still live, as a set
what would narrow it   only where a named, available measurement exists. Where none does,
                       say "no measurement available on this platform separates these."
```

The last line is the honest form and it is better than a remedy column that lists something the
lab cannot do. Your taxonomy has eleven rows wired to nothing; wire the ones whose remedy is real
and mark the rest as unseparable rather than inventing remedies to fill the column.

---

## 8. Revised design

Your section 4 structure survives. Three changes, one addition, one deletion.

```
A. NAME THE REFUSAL.  Keep, with the four-field form in section 7, and mark unseparable
   causes as unseparable rather than supplying a remedy.

B. PER-CHROMOSOME MEASUREMENTS TABLE.  Keep, and run it regardless of verdict. Columns:
     informative (one-parent-hom, per parent)
     informative (opposite-hom)
     call rate
     absence vs father, absence vs mother
     onlyA, onlyB, survivor share          <- ADDED: the 5.2 discriminator
     bafExtreme                            <- ADDED: the gate's own input
     state: parent-named / both-depleted / unmeasurable, from survivor share + call rate
   Every suppressed field sits beside the number that suppressed it.

C. SEGMENTS ON CHAOTIC CELLS.  CHANGED, and this is new work rather than reporting:
   stop emitting segment verdicts where a majority of autosomes are affected, because
   externalNull is past its breakdown point there. Per section 4.

D. DELETED from your proposal: naming "both parents depleted over this window" as a
   positive statement about a region. Symmetric depletion is what dropout predicts; name
   the artefact, not the region.
```

**And a prior step, before any of it.** Run P1 to P3. If the null is what is suppressing chr16 and
manufacturing the contested windows, then fixing `externalNull` - a trimmed or parent-referenced
null rather than a within-sample median - may recover the three missed deletions and dissolve a
third of the contested segments outright. That changes what the reporting layer has to report, so
it comes first. It is one instrumented run on data you already have.

I would not guess at the replacement null here. The obvious candidates - a trimmed mean, the
minimum over autosomes, a null derived from h and the window's own call rate - have different
failure modes and the choice should be pre-registered against the 595-autosome truth, not argued.

---

## 9. Where I am least certain

- **h = 0.100** comes from an allele-frequency spectrum (Beta) fitted to reproduce the panel's
  measured 0.168 bulk heterozygosity, not from the panel's actual spectrum. The *form* h/2 and the
  symmetry argument are exact; the number is an estimate. Compute h directly from the two parental
  arrays, where it is a sum over known genotype pairs with no spectrum assumption, before it is
  load-bearing.
- **The 171-marker figure** for 16q23.1-tel on the opposite-hom set is scaled by the project's own
  median set-size ratio, not counted. It decides whether the survivor share works at window scale,
  so count it directly.
- **Panel a is a schematic.** Real per-chromosome absence rates are a distribution and the
  breakdown will be a ramp; the direction and the 50-percent breakdown point are the claims.
- **P3 is the load-bearing prediction.** If the contested windows turn out strongly asymmetric, my
  whole 5.2 answer is wrong and the dropout reading with it.
- **The `nonParental` suggestion for 5.3 is unscored.** I am confident it is parent-referenced and
  therefore survives the partition; I have not seen it separate a contaminated reaction from a
  chaotic genome on real data, and it may be confounded by genotyping error at the rates these
  arrays show.
- I have not re-derived whether a **trimmed** null fixes the problem or merely moves its breakdown
  point. On a cell with 15 of 22 affected, no robust location estimator of the sample's own
  chromosomes can work, so the replacement may have to be parent-referenced rather than
  within-sample. That is the open design question this response leaves.

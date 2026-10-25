# Response 2: the refusal split, the two judgements, and seven open questions

What this rests on. Repo at `/Users/ezrakruger/claudecodegeneralworkspace/originmarker`, read
2026-10-25: `web/src/stage.ts`, `scoreSample.ts`, `parentage.ts`, `segments.ts`, `ingest.ts`, and
`audit/zuccaro-table-s3.tsv` recounted independently. chr16 cytobands verified against UCSC hg19
(`api.genome.ucsc.edu`, cytoBand track): 16q21 starts 57,400,000; 16q23.1 starts 74,100,000;
chr16 ends 90,354,753. So `16q23.1-tel` is 16.3 Mb = 18.0% of chr16, and `16q21-tel` is 33.0 Mb
= 36.5%.

What I could not do: GSE148488 is still not on this machine, so no audit was re-run. Node is
present (v26.8.1) but there is nothing to run it against. Every detection figure below is quoted
from your document. Everything I assert as new comes from the source or from the truth table.

Headline. Your Section 0 refutation is correct and I have confirmed it independently. On
Section 6, the first judgement is **wrong as stated** - there is a dependency, I traced it, and
it is load-bearing for 5.41.0 - though the conclusion survives for a reason different from your
argument. The second judgement is also wrong, but in the opposite direction from the one you
feared. And 5.1 has an answer that is not a new constant.

---

## 0. My Section 0 was wrong. The mechanism, not just the fact.

Confirmed from `zuccaro-table-s3.tsv` directly. Your six contested cells:

```
cell         nullisomic  no-call  karyotype            in my predicted list?
GSM4552378        5         0     chaotic aneuploid    yes
GSM4552380        5         0     chaotic aneuploid    yes
GSM4552382        1         0     chaotic aneuploid    yes
GSM4552385        0         0     chaotic aneuploid    NO
GSM4552386        0         0     chaotic aneuploid    NO
GSM4552420        0         0     chaotic aneuploidy   NO
```

Three of my six carry zero nullisomy and zero no-calls. P1 fails as written.

The three I named that you do not have - GSM4472409 (10 nullisomic), GSM4552418 (21), GSM4552422
(22 no-call, hand `failed`) - are the three most degraded cells in the corpus, and that is
exactly why they are absent: `stage==='failed'` empties the chromosome list before any scan runs.

**My error, named precisely.** I predicted nullisomy would drive the contested flag. Nullisomy
and stage-failure are the same cells, and the failed-stage gate fires first. I reasoned about the
biology and never checked the control flow - the same error I had just warned about in my own
point 6, that "correct by construction" covers the marker set and not the window. I checked the
overlap test's arithmetic in `scoreSample.ts` and never asked what reached it.

What replaces it, and it supports your original reading rather than mine. Your six cells are
**mid-burden**: 12 to 15 abnormal chromosomes of 22, 7 to 12 monosomies each, none nullisomic
enough to be refused.

```
GSM4552378  abnormal 13/22  monosomy  7
GSM4552380  abnormal 15/22  monosomy  8
GSM4552382  abnormal 12/22  monosomy  9
GSM4552385  abnormal 13/22  monosomy 12
GSM4552386  abnormal 13/22  monosomy 12
GSM4552420  abnormal 12/22  monosomy 12
```

A both-run overlap selects cells degraded enough for wide dropout and not degraded enough to be
refused. That is the dropout explanation you gave originally. The nullisomy hypothesis was not
merely unproven, it was testing a population the gate had already removed.

---

## 1. Section 6, first judgement: there IS a dependency, and it is in the path you rely on

Your argument: the absence statistic reads genotype calls at parent-homozygous markers and never
consults intensity, so a genome with no normalisable baseline costs it nothing.

The statistic does not consult intensity. **The gate in front of it does.** Traced in shipped
code:

```
scoreSample.ts:781   const measured = new Set(result.chroms
                       .filter((c) => c.verdict !== 'not_measured').map((c) => c.chrom))
   ... and the absence scan runs only on chromosomes in `measured`.

parentage.ts:1244    const [ext, nBaf] = t.bafByChrom.get(c) ?? ...
parentage.ts:1245    const extreme = nBaf ? ext / nBaf : NaN
parentage.ts:1253    const collapsed = (callFrac.get(c) ?? NaN) < CALL_COLLAPSE
parentage.ts:1254    const clustered = !Number.isFinite(extreme)
                       || extreme >= BAF_EXTREME_FLOOR || collapsed
   ... and `verdict: !clustered ? 'not_measured' : ...`
```

`bafByChrom` is accumulated from `row.baf` (`parentage.ts:575-579`), which is a B-allele
frequency - an intensity-derived quantity, not a genotype call. So the chain is:

**absence scan → `measured` → `verdict` → `clustered` → `bafExtreme` → BAF → intensity.**

The module says so in its own words at `parentage.ts:241-245`: mis-clustering is something "no
allelic statistic can see... Only the intensity-derived allelic ratio shows the cause." That is
the gate's entire purpose. It cannot do its job without intensity, and the absence scan cannot
run on a chromosome it fails.

Three consequences, in increasing order of how much they matter:

1. The comment at `scoreSample.ts:606-611` is wrong where it says the absence channels "need
   informative markers and a call rate and never the intensity baseline." They need the BAF
   column too, upstream.
2. The dependency is on **per-chromosome BAF distribution**, not on the **per-chromosome
   intensity median**. Those are different quantities. Your argument is right that a genome with
   no normalisable baseline costs the absence scan nothing, because `bafExtreme` is a
   within-chromosome share with no cross-chromosome normalisation - unlike `scanCopyNumber`,
   which normalises against the sample's own per-chromosome median and genuinely does need a
   normal majority.
3. **So 5.41.0 survives, but not for the reason you gave.** The defensible statement is narrower:
   *the absence channel needs BAF to be readable per chromosome, but not comparable across
   chromosomes; the intensity channel needs both.* That distinction is what licenses gating them
   separately. "Never consults the intensity" does not, and it is the sentence currently in the
   code.

One thing to check before relying on this. `bafByChrom` is accumulated **before** the no-call
return (`parentage.ts:575`, deliberately, per the comment at 561-564), so **no-calls are in the
denominator**. On a chromosome with heavy dropout the extreme share is diluted by markers that
have intensity but no genotype, which pushes `extreme` toward the 0.40 floor from above. That
makes `not_measured` more likely on exactly the degraded cells 5.41.0 is designed to admit. It is
worth measuring what fraction of the newly-admitted cells' chromosomes then fail `bafExtreme`,
because the gate you removed at the genome level may still be removing them one chromosome at a
time.

### And a second, sharper problem with that gate

On well-clustered data, homozygous calls sit at BAF ~0 or ~1 (extreme) and heterozygous calls at
~0.5 (not extreme), so `extreme ≈ 1 - het`:

```
apparent het   expected extreme share   vs BAF_EXTREME_FLOOR = 0.40
       0.168                   0.832    pass
       0.250                   0.750    pass
       0.380                   0.620    pass
       0.600                   0.400    pass (at the boundary)
       0.690                   0.310    not_measured
       0.870                   0.130    not_measured
```

The documented mis-clustered case, GSM4774681 chr1 at 0.130, corresponds to a het-equivalent of
0.87 - no genome sits there, which is why it reads as broken. But the mapping means
`BAF_EXTREME_FLOOR = 0.40` is approximately **a second heterozygosity ceiling at het ≤ 0.60,
applied per chromosome**, while `MAX_DIPLOID_HET = 0.25` is the first, applied genome-wide. They
are the same axis at two scales with two different numbers and two separate justifications.
Worth reconciling explicitly: a chromosome can pass the per-chromosome test at het 0.55 on a cell
the genome-wide ceiling refuses at 0.32.

---

## 2. Section 6, second judgement: wrong, but the opposite way round

You fear that a localised heterozygosity **excess** and localised amplification **failure** look
identical in your measurement. They do - but not because the statistic cannot tell them apart.
Because **only one of them can reach the gate at all.**

`hetSpread` is a coefficient of variation, sd/mean. It is symmetric: it measures dispersion with
no sign. Computing what each mechanism can produce, at the panel's 0.168 diploid rate across 22
autosomes:

```
mechanism                                   reachable hetSpread
EXCESS, trisomy-realistic (1.5x on k chrom)   ceiling 0.20
EXCESS, generous (2.0x on k chrom)            ceiling 0.35
DEFICIT (k chrom driven to ~0)                0.22 at k=1, 0.54 at k=5, 1.00 at k=11
```

`HET_SPREAD_LOCALISED = 0.50` is **unreachable by elevation** and easily reached by depression -
it clears 0.50 at k=5 of 22. And the two cells the rule classifies as chaotic read 0.69 and 0.78,
far above anything elevation can produce.

This is arithmetic, not a measurement, and it decides the question: **whatever the rule is firing
on, it is not a heterozygosity excess.** It is a deficit on a subset of chromosomes. A rule named
"a localised heterozygosity excess is a genome" is detecting localised heterozygosity *collapse*.

Your own justification for the constant says this without quite saying it: "a coefficient of
variation near 1 is what a genome with half its chromosomes at one copy produces." Half the
chromosomes at one copy is the **deficit** branch - het goes to zero there, not up. The
constant's derivation and the rule's name describe opposite signs.

Why the rule still works on this corpus: a chaotic aneuploid genome has many monosomies (7 to 12
per cell in your six), each collapsing het on its chromosome, which is a real deficit from a real
genome. A uniform mixture elevates every chromosome equally and reads CV ≈ 0. So the three groups
do separate - but on *deficit versus uniformity*, not *excess versus mixture*.

**Which means your feared confound is real and the rule has no defence against it.** Localised
amplification failure also collapses het locally, also reads high CV, and is not a genome. The
corpus cannot distinguish them because every high-CV cell in it is hand-called chaotic aneuploid.
The hand karyotype agreeing three times is not evidence the statistic separates the two causes;
it is evidence that in nine cells only one cause occurred.

### What would separate them, using quantities already computed

The discriminator is **where the deficit's markers went**, and both destinations are already
measured per chromosome:

- **monosomy**: the chromosome is present at one copy. Markers still genotype. Call rate stays
  near normal; het collapses because there is one haplotype. → `callByChrom` normal, het low.
- **amplification failure**: the locus did not amplify. Markers do not genotype at all. → call
  rate collapses with het.

So: **condition the spread on call rate per chromosome.** A chromosome with low het and *normal*
call rate is hemizygous - a genome. A chromosome with low het and *collapsed* call rate is a
dropout region - not evidence of a genome. `callByChrom` is already accumulated
(`parentage.ts:543-546`) and `CALL_COLLAPSE` already exists as a threshold.

A one-sided, call-rate-conditioned spread would be strictly better than the CV on both counts:
it would stop being sign-blind, and it would stop counting dropout as evidence of ploidy. I would
not ship another corpus-free constant before measuring that, but the statistic is computable
today from columns already in the tally.

---

## 3. Question 5.1: the ceiling is the wrong *shape*, and the fix is not a new constant

Direct answers first.

**Is the ceiling the wrong quantity for amplified material?** Partly - but not in the direction
you expect, and on your five cells it is arguably too *generous*, not too strict.

Work the arithmetic the module's own model implies. Dropout destroys heterozygotes (AB → AA or
BB): it is monotone **down**. Only drop-in creates them. At call rates 0.67-0.77 these cells have
substantial dropout, so their *observed* het is already suppressed relative to their true het. To
observe 0.32-0.38 there:

```
call rate   implied dropout   het after dropout   drop-in needed for 0.32   for 0.38
     0.67              0.33              0.113                     0.234      0.301
     0.72              0.28              0.121                     0.226      0.295
     0.77              0.23              0.129                     0.219      0.288
```

The module's measured drop-in ceiling is 0.0525, from 113 same-genome pairs. The required
drop-in is 0.22-0.30, five to six times that. **On the module's own numbers, a single diploid
plus measured drop-in cannot read 0.32-0.38 at these call rates.**

So if these cells are diploid, the measured drop-in figure does not apply to amplified single
cells - which is your own suspicion, and it is well founded, because drop-in on WGA material is a
different process (chimeric amplicons, template switching) from drop-in on bulk. You do not have
that figure. I cannot derive it from the corpus either; it needs same-genome replicate pairs **on
amplified material**, which the dilution series would have if those replicates were WGA rather
than diluted bulk.

**But there is a hypothesis the module does not list, and it is a genome.** The ceiling's caveat
names four causes: two individuals mixed, contaminated reaction, failed array, chaotic genome.
Missing: **triploidy**. Calibrating an allele-frequency spectrum to the panel's measured diploid
rate of 0.168:

```
state                                  expected het    + measured drop-in
diploid (2 copies)                            0.166                 0.210
TRIPLOID (3 copies)                           0.249                 0.288
tetraploid / two unrelated genomes            0.302                 0.339
```

A triploid sits essentially **on** the 0.25 ceiling and is refused by construction. This material
produces triploids routinely - MII injection with retained second polar body (digynic) or
dispermy (diandric). And a by-eye SNP-array karyotype **cannot see uniform triploidy**: it reads
relative copy number across chromosomes, and a triploid is flat. "46,XX" and "normal for other
chromosomes" are exactly what a hand reader would write for a triploid.

I want to be careful here: at 0.32-0.38 your five cells read above the triploid expectation and
closer to the tetraploid/two-genome figure. So I am **not** claiming they are triploid. I am
claiming the ceiling cannot distinguish a triploid from a mixture, that triploidy is absent from
its stated list of causes, and that this is a genome the tool will refuse whenever it meets one.

**The form of the fix, and it is not a call-rate-conditional constant.** You asked whether a
call-rate-conditional ceiling is the right shape. I do not think so, for a specific reason: the
call rate enters with the wrong sign. Dropout lowers het, so conditioning the ceiling on call
rate would *raise* the ceiling on low-call-rate arrays, where the observed het is already
suppressed - making the gate looser exactly where the original failure occurred (an array at 52%
het and 55% call rate). That reopens the failure you do not want reopened.

The right discriminator is **BAF band structure at heterozygous calls**, and it is already
computed:

```
state                      BAF bands at het calls      0.5 band occupied?
diploid                    0.50                        yes
triploid                   0.33, 0.67                  NO - empty
tetraploid / 50:50 mix     0.25, 0.50, 0.75            yes
two genomes, uneven ratio  continuum                   smeared
```

A triploid has **no 0.5 band at all**; a mixture retains one. These are not subtle differences,
and the occupancy of [0.45, 0.55] separates them in one pass over a column `ingest.ts` already
reads (`accumulateBaf`, band = share in [0.35, 0.65]).

Crucially, **the precedent exists in the same function.** `inferStage` already lets the band
override heterozygosity in the haploid direction (`BAND_DIPLOID_CERTAIN`, line 400: "Intensity
outranks the genotype calls in this one direction"). The proposal is to let it override in the
other direction too. That is not a new mechanism, a new corpus, or a new constant - it is the
existing override applied above the ceiling as well as below it.

Concretely:

```
Replace the unconditional refusal at `h > maxDiploid` with a two-way test:
  h > maxDiploid AND the 0.5 band is OCCUPIED      -> failed (mixture/contamination). Unchanged.
  h > maxDiploid AND the 0.5 band is EMPTY,
    with density at 1/3 and 2/3                    -> triploid. A genome. Report it as one.
  h > maxDiploid AND BAF is smeared with no modes  -> failed (array destroyed). Unchanged.
```

The 52%-het-at-55%-call-rate array that prompted the ceiling would still be refused: at 52% het
with two individuals mixed, the 0.5 band is occupied. The failure does not reopen.

Three caveats I will not bury. The band thresholds need measuring on this platform before they
are constants - I have given the expected positions, not measured ones. Amplification noise
smears all three patterns and the amount of smearing is exactly what `BULK_MAX_BAF_SD` says
varies 0.21-0.30 on this material, which may be enough to fill a triploid's 0.5 gap. And this
does not explain your five cells on its own, since they read above the triploid expectation. It
makes the gate answer the right question; it does not promise those five pass.

---

## 4. Question 5.3: not a sensitivity problem, and not the same question as 3.1b

Direct answer: **separate problem, and the window floor is not the binding constraint.**

Using the verified UCSC coordinates, `16q23.1-tel` is 16.3 Mb = 18.0% of chr16. Scaling the
one-parent-homozygous informative count by chromosome length (median autosome 134.4 Mb, chr16
90.4 Mb, so ~10,750 informative on chr16), the event contains roughly 1,935 informative markers
against a floor of 2,400. So the event is below the floor, and the smallest admissible window
(2,400 markers ≈ 20.2 Mb) is 1.24x the event.

That is the natural reading and it is wrong. The question is not whether the window fits the
event; it is whether a window **overlapping** the event scores high enough. Computing the LRT for
a 2,400-marker window of which 1,935 markers are inside a maternal deletion, against the scan's
own null of 0.003:

```
in-event absence (after dropout)       k     rate       LRT    clears 250?
                            1.00    1936    0.807   10069.4    yes
                            0.60    1162    0.484    5091.6    yes
                            0.30     581    0.242    2052.3    yes
                            0.10     194    0.081     459.7    yes
```

Even at 10% in-event absence the diluted window clears the threshold by a factor of 1.8. **The
floor is not what is stopping these three.** Dilution costs power here, and there is far more
power than the task needs.

So this is not the same question as 3.1b, and raising the window floor's derivation will not
recover these events. Three candidates that the code does support, all checkable in one
instrumented run:

1. **Peak-picking discard.** `scanChromosome` keeps the highest-scoring window and drops
   everything overlapping it. If a larger or stronger window on chr16 also clears, the terminal
   event is absorbed and never reported. Your note that chr16 shows no `lostParent` rules out a
   whole-chromosome call absorbing it, but not a larger segment window.
2. **The maternal run's marker set.** These are maternal deletions, scanned against
   `matAbsenceByChrom`. If chr16's maternal-homozygous informative count is thin, the floor may
   be unreachable on *that set specifically* even though it is comfortable on the paternal one.
   This is the one I would check first, because it is the only candidate that explains why three
   sibling cells all fail the same way.
3. **The `measured` gate.** Confirm `bafExtreme ≥ 0.40` on chr16 for all three. One reads
   `unclear`, which is a different field, but the gate is upstream of the scan and worth ruling out.

The length-scaling is an estimate and I will flag it: Axiom spacing is severely non-uniform, so
the 1,935 figure could be off by ~1.5x either way. It does not change the conclusion, because the
LRT clears by far more than that margin. Counting informative markers in 16q23.1-tel directly
would settle it.

---

## 5. Question 5.4: two defects on one chromosome, not one sign error

Direct answer: **the sign error is in neither channel. The absence call is correct and
mis-worded; the copy-loss is a separate false positive.**

**The absence call.** On `dup 16q21-tel (pat)`, the region is 2 paternal + 1 maternal copies.

- *Paternal run* (markers where the father is homozygous): the father's allele is present on two
  copies. Absence goes down. The paternal run cannot fire. Correct.
- *Maternal run* (markers where the mother is homozygous): at a marker where the mother is BB and
  the father AA, the sample is A,A,B - allelic ratio 2:1, so BAF sits at 1/3 rather than 1/2. The
  caller's AB cluster is centred on 0.5 and the AA/AB boundary sits near 0.25; a 1/3 point is
  only ~0.33 sd from that boundary at this material's het BAF spread of 0.21-0.30. A large
  fraction of those sites miscall to AA, and **the mother's allele reads absent.**

So the absence scan fired correctly on a real allelic imbalance. "The mother's alleles are absent
here" is *true*. It means "the father's copy is duplicated", and the tool has no vocabulary for
that. This is a **semantics** defect, not a sign error: a one-sided statistic reporting a
two-sided phenomenon.

Falsifiable prediction: that segment should carry `absentParent: 'maternal'` - the run against
the **non-duplicated** parent - and its span should coincide with 16q21-tel (33.0 Mb, 36.5% of
chr16, per UCSC). If it reads paternal, this mechanism is wrong.

**The copy-loss is unrelated.** `scanCopyNumber` scans no-call runs and keeps a window only if
median log2R - genomeLrr ≤ -1.0. Copy states:

```
one copy   log2R -1.000   passes a <= -1.0 filter
two        log2R  0.000   filtered out
three      log2R +0.585   filtered out - cannot pass in EITHER sign convention
```

A 3-copy region sits at +0.585 and cannot clear a -1.0 threshold however the sign is read. The
copy-loss is therefore **not** the duplication inverted. GSM4552365's hand karyotype is otherwise
normal, so `genomeLrr ≈ 0` and there is no baseline shift to flip it. Your own code notes that
rejected failures span -0.72 to +0.38 while kept ones span -0.97 to -1.94, so a severe
amplification failure does reach the kept range: this is a dropout false positive elsewhere on
chr16 that the 1.0 shift filter was meant to catch and did not. Note 64% of chr16 lies outside
the duplication, so there is ample room for it.

Two defects, two repairs, and the fact that they landed on one chromosome is coincidence.

---

## 6. Question 5.2: Test J is not a valid bound for this population

Direct answer: **no, and I would ship the eleven flagged as unscoreable.**

My own pronuclei caveat applies with more force than you state, and for a reason beyond dropout
profile. A pronucleus is **uniparental by dissection**: every marker where the other parent is
homozygous-different is informative, and the absence signal is genome-wide and maximal. An
aneuploid blastomere is biparental with localised losses. The false-positive mechanisms are not
the same population, and 2 in 58 on material whose truth comes from dissection does not transfer
to a cell where the question is whether a *local* excursion is real.

There is also a selection problem: Test J's 58 scans are on "poor" material selected by a quality
criterion, and GSM4472405 is selected by being hand-called chaotic. Those are different poorness.

Is any bound available on public data? Not a sound one, for the same reason the specificity
question has been hard all along: you need material that is (a) aneuploid, (b) amplified, (c)
exhaustively truthed rather than truthed for the events someone was looking for. The corpus's own
answer key is explicitly not exhaustive, which is what creates the problem.

So: ship them flagged. "Eleven segments on a cell whose hand karyotype is not exhaustive; not
scoreable against this answer key" is an honest statement and a useful one. The alternative -
quoting a bound measured on dissected haploid products - would be the same category of error as
the shipped 0-of-110 figure you just corrected.

---

## 7. Questions 5.5, 5.6, 5.6b, 5.7, briefly

**5.5 - can any public corpus discriminate the ladder?** Probably not, and the reason is
structural rather than empirical. The two channels are near-tautologous by construction: a window
where a parent's alleles are absent is a window where that parent's allele share is low. They
come apart only where the share is low for a reason *other* than absence - allelic imbalance
without loss, which is precisely the duplication case in 5.4. So the discriminating material is
**duplication-bearing cells with parental arrays**, where the share drops and the absence is not
a loss. GSE60909 may contain some; this corpus has five duplication-only events and you have
already found the channels disagreeing on one of them in substance if not in label. That is the
test, and it is worth saying that 5.4 is the first real evidence the ordering matters.

**5.6 - is there a derivation for 0.50?** Not the one you have, and your instinct that the
argument is worth less than it sounds is correct - the CV=1 limit comes from the *deficit*
branch, as section 2 above shows, while the rule is named for an excess. Any derivation should
come from the one-sided, call-rate-conditioned statistic instead, where the threshold has a
meaning: how many chromosomes must show hemizygous-level het at normal call rate before the
genome is chaotic rather than noisy. That is a count with a binomial null, not a gap midpoint.

**5.6b - does coverage change the exchange rate?** It does something worse: it makes the raw
counts incomparable across grid points, and it explains the corner better than my cap argument
did. Coverage runs 1.84 → 36.62 → 83.75 → 97.47 → 99.49 percent as the window falls - monotone in
exactly the direction the rule walked. So the objective was monotone in an **unrecorded nuisance
axis**, and "the rule walked to the corner because the cap was slack" is incomplete. It walked
partly because smaller windows let it look at more chromosomes, which is not the same thing as
being more sensitive.

For the exchange rate specifically: between 300 (97.47%) and 150 (99.49%) coverage is nearly
matched, so the k = 1.25 breakeven I gave stands for that pair. Against 600 or 1200 it does not,
and those points cannot be priced without a per-scanned-chromosome denominator. Add coverage as a
**third pre-registered endpoint** and define detections and false calls **per scanned
chromosome**; then the grid is comparable and the coverage gradient stops masquerading as
sensitivity.

**5.7a - the KL orientation. Resolved, by derivation, and the dilution file is not needed.** The
scan's per-marker score is log[P1/P0]. Its expectation under the true event distribution P1 is
E_P1[log(P1/P0)] = KL(P1 || P0) - that is the definition of KL, not a modelling choice. So a
**detection** floor uses KL(h1 || h0). Monte Carlo over 200,000 markers confirms it: at h1 = 0.30
the empirical mean score per marker is 1.1391 against KL(h1||h0) = 1.1340 and KL(h0||h1) =
0.3388. KL(h0||h1) is the expected score under the *null*, which governs false positives - a real
quantity, answering a different question. **So the admissible range is 76% dropout, not 43%, and
my conservative hedge in 3.1b can be withdrawn.** Window 300 clears the floor across the whole
usable range of this material.

**5.7b - my CytoSNP figure of 377 per chromosome is wrong by 3.4x, and in the wrong direction.**
I scaled by array density alone. The opposite-homozygote yield depends on the **allele-frequency
spectrum**: for two unrelated parents, P(opposite homozygotes) = 2p²q² = H²/2, so per-marker
yield scales as H², not as density. Your own `stage.ts` already records the panel heterozygosity
for both platforms - Axiom 0.141-0.170, CytoSNP-12 0.301-0.316:

```
H ratio (0.308/0.168)              1.833   ->  per-marker yield ratio H^2 = 3.361
density ratio (301,232/825,656)    0.365
net per-chromosome ratio                       1.226

Axiom measured opposite-hom per chromosome     1,034
CytoSNP-12, corrected estimate                 1,268     (I said 377)
```

**CytoSNP has MORE opposite-homozygote markers per chromosome than Axiom**, because its markers
are far more polymorphic, which more than offsets its 2.7x lower density. This reverses the
implication of my 3.3(c): the platform transfer is more favourable than I said, not less. At
~1,268 markers per chromosome, minWindow 300 gives 3 scales and 19 windows per chromosome;
minWindow 2400 remains unreachable.

Still an estimate - it uses mean H, and E[H²] ≥ (E[H])², so both panels' yields are understated
and the ratio less so. But it is now anchored on a measured panel property rather than on density
alone, and it is computable exactly from any one CytoSNP trio.

---

## 8. What I would change, in order

1. Fix the comment at `scoreSample.ts:606-611`. The absence channel does depend on intensity,
   through `measured` → `bafExtreme` → BAF. Replace the claim with the narrower true one: it
   needs BAF readable **per chromosome**, not comparable **across** chromosomes. 5.41.0 stands on
   that; it does not stand on "never the intensity baseline."
2. Measure how many of the newly-admitted cells' chromosomes then fail `bafExtreme`. No-calls are
   in that denominator, so the genome-level gate you removed may be re-applying per chromosome.
3. Rename and re-found `HET_SPREAD_LOCALISED`. It fires on a heterozygosity **deficit**;
   elevation cannot reach 0.50. Replace the CV with a one-sided, call-rate-conditioned count -
   chromosomes at hemizygous-level het **and** normal call rate - which distinguishes a genome
   from localised amplification failure. Both inputs are already in the tally.
4. Make the ceiling two-way using the BAF 0.5 band, per section 3. Triploidy is a genome this
   material produces and the current gate cannot see it. The 52%/55% failure stays refused.
5. For 5.3, check the maternal-homozygous informative count on chr16 before anything else. The
   window floor is not the constraint; the LRT clears by 1.8x at 10% in-event absence.
6. For 5.4, confirm `absentParent: 'maternal'` on that segment, then treat the copy-loss as a
   separate dropout false positive. Two defects, not one sign error.
7. Ship the eleven segments flagged unscoreable. Do not bound them with Test J.
8. Add coverage as a third pre-registered endpoint and denominate detections and false calls per
   scanned chromosome.
9. Withdraw my 43% hedge: the KL orientation is settled by derivation at KL(h1||h0).
10. Replace my CytoSNP figure of 377 with ~1,268 pending direct measurement on a real array of
    that platform.

## 9. Where I am least certain

- The triploid hypothesis for 5.1. Your five cells read 0.316-0.376, above the triploid
  expectation of ~0.25 and nearer the two-genome figure. I am confident the ceiling *cannot see*
  triploidy and that triploidy is missing from its list of causes; I am **not** claiming these
  five are triploid. The BAF band test decides it and costs one pass.
- The BAF band positions in section 3 are expected values, not measured ones, and amplification
  smearing of 0.21-0.30 sd may partly fill a triploid's 0.5 gap. Measure before making them
  constants.
- The chr16 marker count (1,935 in the event) is length-scaled on a non-uniform array; good to
  about 1.5x. The conclusion survives because the LRT margin is much larger than that.
- The `extreme ≈ 1 - het` mapping assumes well-clustered data. It is exact enough to show the two
  ceilings share an axis, but the per-chromosome boundary of 0.60 het-equivalent is approximate.
- The allele-frequency spectrum behind the triploid and CytoSNP arithmetic is a Beta fitted to
  reproduce the panel's measured diploid rate, not the panel's actual spectrum. Ratios are robust
  to that; absolute figures less so.

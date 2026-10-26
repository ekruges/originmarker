# Consult: designing what a refused or chaotic array may report

Third consult in this thread. The first two are `CONSULT-segmental-parent-of-origin.md` and
`CONSULT-refusal-split-and-attribution.md`, with your responses committed beside them. This one is
a **design** request rather than an audit: I want the scientific design of a reporting layer
decided before I write it, because the last two things I built in this area were both reverted and
both for the same reason, which is that I designed from an argument and measured afterwards.

Public data throughout. No laboratory sample is named.

## 0. The requirement, in the laboratory's own words

From the 21 October call, verbatim, after I said I was filtering chaotic cells out:

> **"But aren't those cells the ones that we're interested in, the ones that are aneuploid?"**
>
> **"I think the ones that are chaotic, so there's a lot of aneuploidies, those are still important
> for biology... even in some early embryo samples, I can see there's so much aneuploidy, so many
> chromosome gains and losses, but knowing which one is which is still going to be important."**

And from the earlier call, which is the same request put as a method:

> "if all the B alleles were lost, then it would be the mother's... Then you would just call which
> SNPs is missing. If it's the mother's SNPs that are missing, then it's the mother's chromosome
> that was lost."
>
> **"Can you print it out... are we able to check humanly these it's making the right call?"**
> "I think for some of these things, it's hard to really figure out what the program might not be
> doing correctly because we can't see what it's going, what it's doing. If we can check, then we
> can see: okay, maybe this thing is too stringent."

So the ask is not "call more". It is "show me the readings so I can check you", on exactly the
cells the tool is least able to call. She is right that the allele channel needs no baseline, and
she identified the dependency before I did, in the same breath as asking what log R normalises to.

## 1. What I had wrong until an hour ago, because it changes the problem

I have been saying the tool filters out chaotic cells. It does not.

```
cells the published karyotype calls chaotic      12
of those, refused by the tool                     2
of those, analysed and emitting calls            10
```

And the nine refused cells are mostly not chaotic:

| cell | abnormal autosomes | normal | karyotype |
|---|---|---|---|
| GSM4472409 | 22 | 0 | chaotic aneuploid |
| GSM4552418 | 21 | 0 | chr6p fragment only |
| GSM4472405 | 9 | 13 | chaotic aneuploidy |
| GSM4472422 | 4 | 18 | 45,XY,-5,-13,+15 |
| GSM4552389 | 2 | 20 | monosomy 5 |
| GSM4552387 | 1 | 21 | normal for other chromosomes |
| GSM4552388 | 1 | 21 | normal for other chromosomes |
| GSM4552422 | 0 | 0 | failed |
| GSM4472421 | 0 | 22 | 46,XX |

So the live problem is not the refusal. It is that on the ten chaotic cells we already analyse,
**a third of what we emit contradicts itself**: 35 of 69 parental-absence segments are windows the
scan called against BOTH parents, and every one of those 35 is on a chaotic cell. Zero occur on a
clean cell. A single lost copy takes one parent's alleles with it, so both cannot be the missing
ones over one window.

And the three missed deletions nobody has been able to explain, `del 16q23.1-tel (mat)` on three
sibling cells of one embryo, are **also** all on chaotic cells we analyse, with 14,322 to 14,991
maternal-homozygous markers available on chr16 against a window floor of 2,400. Six times the
floor, and the scan finds nothing.

## 2. What a refused array reports today

It reports a diagnosis and no analysis.

The diagnosis is real. `StageCallout` renders for a refused array carrying the stage, heterozygous
rate, BAF spread at heterozygous calls, whether it reads amplified, dropout, the marker floor, and
a `why` that is specific: *"34.5% heterozygous exceeds the 25% ceiling a single diploid can reach
on this platform, which is its 16.8% rate plus drop-in at the highest rate measured here. No
genome reads this way, so this is not a stage."*

The analysis is nothing. The refusal clears every per-chromosome field:

```js
for (const c of result.chroms) {
  c.aneuploidy = undefined; c.lostParent = undefined; c.aneuploidyParent = undefined
  c.mixed = undefined; c.alleleGain = undefined
}
result.segments = []
```

Two weaknesses in the diagnosis half. The `caveat` lumps four causes into one sentence, *"two
individuals mixed, a contaminated reaction and a failed array are not distinguishable from each
other here, only from a genome"*, so it says "not a genome" without saying which, and offers no
remedy. And `audit/refusal_taxonomy.csv` already holds the richer form, eleven rows of trigger,
measurement and remedy, including one row that is exactly the governing principle here, **"Report
the ORIGIN and withhold the CLASS as separate fields"**. It is wired to nothing.

## 3. What is settled, so you need not spend effort on it

All measured this round, all recorded with pre-registrations.

- **The diploid ceiling is right and stays.** Tested against 20 single cells from lymphoblastoid
  lines of two siblings against multi-cell controls of the same genomes, 279,998 autosomal
  markers. Amplification **lowers** heterozygosity: MDA reads 0.048 to 0.121 and PicoPlex 0.147 to
  0.264 against their own controls at 0.314 to 0.318, and median amplified drop-in is 0.033 against
  the 0.0525 bulk figure the ceiling is built on. A cell at twice its platform's bulk rate is not
  an amplified diploid.
- **The heterozygosity spread statistic is conditioned on per-chromosome call rate** and validated
  against four classes: 0 of 5 bulk diploid arrays and 0 of 15 dissected pronuclei read localised,
  where the unconditioned version read 0.774 to 0.982 on three pronuclei and would have called
  dissected haploid products chaotic genomes. It admits 0 of 9 refused cells.
- **The absence channel does depend on intensity**, through `measured` to `verdict` to `clustered`
  to the per-chromosome extreme-BAF share. Your trace was right. The defensible statement is that
  it needs BAF readable **per chromosome** and not comparable **across** chromosomes.
- **Detection**, measured positionally against UCSC hg19 cytobands rather than by chromosome name:
  deletions 15 of 22 with any overlap and 14 of 22 at reciprocal overlap 0.5, paternal 10 of 14
  against maternal 4 of 8; duplications 0 of 6.
- **Attribution**, where an event is detected: 15 of 15 against the published karyotype.
- **False positives** on dissection truth: 0 over 88 chromosome scans on clean pronuclei, 2 over
  58 on poor ones at call rates 0.54 to 0.59.
- **Triploidy cannot be tested on public data.** Zero of the 682 truth rows mention it; GEO carries
  no human triploid embryo genotype array with parents; the one published study with triploid trios
  states that GDPR prevents release. Your BAF half-band discriminator has no positive class here
  and I have not built it.

## 4. The design I am proposing, which I want replaced if it is wrong

Two pieces, neither of which loosens a gate.

**A. Name the refusal.** Replace the one lumped caveat with the class actually fired, from the
taxonomy: what was measured, which of the causes it is or which remain possible, and what the
laboratory does about it.

**B. A per-chromosome measurements table that runs regardless of the verdict.** Per chromosome:
informative markers, absence rate against the father, absence rate against the mother, call rate.
**Numbers with no call attached.** A chromosome where both rates are elevated then shows as both
rates elevated, which is the truth, rather than being resolved into a parent or vanishing behind a
refusal.

My reasoning for B over a per-chromosome parent call: a third of what we would call on these cells
is self-contradictory, so a call would be wrong a third of the time, while the numbers are never
wrong. And it is what she asked for twice: print it so I can check you.

## 5. The questions I want answered before I write any of it

Ranked by how much the answer changes the build.

### 5.1 What is the complete set of quantities computable on a genome with no normal majority?

I have asserted that the per-chromosome opposite-homozygote rate against each parent qualifies,
because it is a ratio within one chromosome, and that copy number does not, because log R
normalises to the sample's own per-chromosome median and a chaotic genome has no such median.

**Derive the set rather than checking my two.** Every quantity the pipeline computes should be
classifiable by what it references: this chromosome only, this chromosome against the genome, this
sample against a parent, this sample against itself. I expect that partition to put things on the
reportable side that I have not thought of, and possibly to disqualify something I assumed was
safe. Note that `clustered` already shows my intuition is not reliable here: I called the absence
scan baseline-free and it is not.

### 5.2 How should a both-parents-absent chromosome be presented?

This is 35 of 69 segments, so it is the common case on exactly the cells in question, not an edge.
The possibilities I can see, and I do not know which is right:

1. Both numbers, no interpretation. Honest, and leaves the reader to notice that both being high
   is impossible for one lost copy.
2. Both numbers plus a named class, "both parents depleted over this window", which is a positive
   statement about the region even though it names no parent.
3. A third quantity that discriminates. Dropout wide enough to clear the window against both
   parents should have a signature distinguishable from a region genuinely missing both copies,
   and the intensity would be the obvious place, except that the intensity is the channel that
   does not work on these cells. If there is a signature in the allele data alone I cannot find it.

### 5.3 What measurably separates the four refusal causes, given no positive class for one of them?

The taxonomy wants to name which cause fired. Two individuals mixed, a contaminated reaction, a
destroyed array, and a chaotic genome. I have one discriminator that works, the per-chromosome
spread, which separates localised from uniform and is validated on 20 negative controls. I have
nothing separating mixed from contaminated from destroyed, and the measurements I do have on those
five cells are that they carry roughly 46 percent of their heterozygous mass outside the half,
third and two-third BAF bands, against 30 percent on euploid controls.

Is there a sound discriminator available, and if there is not, **what is the right thing to print**
for a cause that cannot be narrowed? A list of three possibilities with the one measurement that
would narrow it, or something better.

### 5.4 Does a measurement without a verdict become a verdict in the reader's hands?

The hazard I can see in my own proposal. A table of absence rates per chromosome, printed on a
cell the tool has declined to call, invites the reader to make the call the tool refused, with none
of the gating that refusal represented. The laboratory's position is that their manual call is the
benchmark, so a reader doing exactly this is the intended use, and I am not sure whether that makes
the hazard acceptable or worse.

Is this a real risk, and what in the design mitigates it without defeating the purpose?

### 5.5 How do you validate a display that deliberately makes no calls?

Everything in this project is validated by scoring calls against a truth. A measurements table
makes no calls, so there is nothing to score, and "the numbers are correct" is a unit test rather
than a validation. What is the test that this layer is *useful and not misleading*, and is there
one that does not require showing it to the laboratory and asking?

### 5.6 Should the ten already-analysed chaotic cells keep emitting verdicts?

The question I least want to ask, because it is about shipped behaviour rather than new work. On
those ten we currently emit parent calls, and a third of the segments contradict each other. The
options are to keep calling and flag the contradictions, to drop to measurements-only on any cell
whose contested fraction passes some level, or to leave it and say so loudly.

There is no truth to decide it with: Zuccaro tabulated the events they were looking for, so of the
twelve segments on one of those cells, one matches a hand call and eleven are unscoreable. Nothing
public is simultaneously aneuploid, amplified and exhaustively truthed, which you established last
round and I have not been able to get past.

### 5.7 And the three missed deletions, since they are in this population

`del 16q23.1-tel (mat)`, 16.3 Mb terminal, three sibling cells of embryo C, all three chaotic, all
three analysed, 14,000 to 15,000 maternal-homozygous markers on chr16, window floor 2,400, chr16
verdicts `unclear`, `present`, `present`. Your previous first-choice explanation was a thin
maternal marker set and it is not thin. Your arithmetic also showed the window floor is not the
constraint, the LRT clearing 250 by 1.8x even at 10 percent in-event absence.

Is this the same problem as the contested windows, a chaotic cell's dropout swamping a real
terminal event, or a separate one? If it is the same problem, it is evidence for whatever answer
5.2 gets.

## 6. Constraints that are not negotiable

- Nothing ships without a measurement, and the measurement is pre-registered before the run. The
  two things I built here without that were both reverted.
- No new free parameter unless its value is derived or sits in a measured gap, and either way the
  derivation goes in the file beside it.
- No laboratory sample identifier in the repository; a gate enforces it and has caught me once.
- A refused array must not become a called array by the back door. Whatever is printed is printed
  as a measurement, and the refusal stands.

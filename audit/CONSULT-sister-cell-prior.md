<!-- Laboratory sample identifiers redacted: this file is tracked in a public repository.
     See scripts/release-check.sh, which refused the first commit of an earlier consult. -->

# Consult: a sister cell's breakpoint as a prior on its sibling's gain

Written 8 November 2026 against version 5.47.0, commit 47e48d4. Companion to
`audit/CONSULT-segmental-gain-statistic.md`, which covers the detector itself. Read that one first
if only one is read, because this construction inherits its floors.

## The request

Our collaborator in the Egli laboratory, 7 November, in her words: cleavage-stage cells usually
carry reciprocal gains and losses, so if one cell lost a piece the sister gained it. Where a loss
breakpoint can be called definitively in one cell, use it to pinpoint the gain in the sister.

She is pointing at the right half of our error table. Of 15 missed gains, 6 are boundary
disagreements rather than absences, and 4 of those 6 sit within 0.06 of the 0.5 reciprocal-overlap
bar. On losses the pattern is stronger: 13 of 23 misses are boundary disagreements. The detector
is seeing these events and drawing them differently.

## 1. How often the prior would be available

Counted today from the laboratory's own manual region table, 93 hand-called regions over 21 cells:

**19 of the 26 hand-called gains have an overlapping hand-called LOSS in another cell of the same
embryo group, on the same chromosome.**

Two caveats we cannot resolve ourselves and which bound how much that number is worth:

- The grouping is the barcode prefix, which is the same rule `audit/joy-common.ts` applies. That
  harness carries its own warning that a prefix group is not always an embryo: in one public series
  the two members of a group are a paternal and a maternal pronucleus from different oocytes, so a
  reciprocal pair between them would mean nothing. We are confirming the groups against the
  laboratory sheet before quoting the figure to anyone.
- Both sides of the 19 are hand calls. In the tool's own output the loss has to be detected and
  refined first. Losses are our stronger channel but they are 44 of 67 at the 0.5 bar, so the
  availability under tool-detected losses will be lower than 19 and we do not yet know by how much.

## 2. What exists to build on

**`web/src/refine.ts`.** Breakpoint refinement by likelihood drop. `LR_DROP = 12`,
`SEARCH_RADIUS = 3_000`, `REFINE_CALL_FLOOR = 0.70`, `MIN_INTERVAL_MARKERS = 100`. The drop is
empirical, calibrated to 96.9 percent coverage on 470 edges. A variance-inflation correction for
clustered dropout was tried and reaches only 92.0 percent, so it is recorded as not an improvement.
`scanChromosome` refines its edges and `scanCopyNumber` spreads the result, so copy-number segment
edges are already refined.

**`web/src/siblingOrigin.ts`.** The construction that removes the parent requirement: multiple
blastomeres of one embryo share both parents exactly, so where an event is present in some cells
and absent in others, the unaffected siblings define that embryo's heterozygous sites and the
affected cell is read against them. Its docstring states four traps, three of which bear directly
on this consult:

1. Heterozygosity must come from the SIBLINGS, never from the cell being scored. Selecting on the
   scored cell's own heterozygous calls selects the markers where both alleles survived, which
   forces the share to 0.5 at every retained marker: simulated mean 0.5000 with zero markers at
   either tail. Sibling-established heterozygosity keeps both tails at 15.3 percent each.
2. There is no centre, deliberately. Each marker contributes a likelihood under three hypotheses
   and the region's evidence is their product. The project's one-directional artefacts came from
   centring a share on the median of a distribution whose median is a boundary of its own support.
3. It does not phase and cannot name a parent alone. Observations must already be oriented to a
   consistent haplotype. Given unphased markers, counting AA against BB cancels to noise, so the
   honest decomposition is that it answers deletion versus dropout without orientation and
   answers which side only when phase is supplied.
4. False heterozygosity is the only one-directional term: a marker truly homozygous but called
   heterozygous by the siblings always retains the reference allele, so it always reads as loss of
   the other copy. Its rate phi travels with every call, and above `PHI_WALL = 0.15` the call is
   refused outright.

**`audit/joy-sisters.ts`.** Already groups arrays by embryo, runs the shipped copy-number scan on
each cell, and sorts overlaps into reciprocal and shared. It does not yet join the gain side to the
loss side, which is the work this consult is meant to specify.

## 3. What we intend to build, stated so it can be corrected

Take a loss segment in cell A whose refined edges clear the existing refinement criteria. Treat
its interval as FIXED and evaluate the gain hypothesis in sister cell B at that interval, with no
search over boundaries.

The claim we are relying on: removing the search removes a multiple-comparison penalty, so the
same evidence clears a lower bar and the detectable event size falls. The scan's threshold
`SEGMENT_LRT = 250` is empirical precisely because absence artefact is spatially clustered and a
closed-form tail produced 4 to 30 false segments per clean genome. We do not know how much of that
250 is paying for the search and how much for the clustering, and therefore we cannot say what a
fixed interval buys.

**Question 1: how much detectable size does a fixed interval buy here, and what is the correct way
to state it? If the answer is that an empirically calibrated threshold cannot be decomposed that
way, we need to know before we claim any improvement.**

## 4. What we think is double-counted, and are not sure

Two cells of one embryo share their parents and their genome, so the markers that are informative
in one are informative in the other. They do not share amplification: each cell was amplified
separately, and dropout is the dominant error in this material.

So our reading is that using A's breakpoint on B costs nothing in independence terms for the
DETECTION in B, because the two dropout realisations are independent, but that the two events are
not independent as biological observations, since reciprocity is exactly the hypothesis that makes
them one mis-segregation rather than two.

**Question 2: what is actually double-counted, and does it affect the detection claim, the
reciprocity claim, or both? And if we later report a per-event parental label across a cohort, is a
reciprocal pair one event or two?** The deliverable this project owes Dr Egli is per-event parental
labels with a defensible error rate, so the unit matters.

## 5. Whether the premise transfers

Reciprocity in cleavage-stage embryos comes to us from Vanneste 2009, which reported segmental
deletions and duplications reciprocal between sister blastomeres. Our material is different in
three ways we can name: it is largely pronuclei and blastomeres from gene-editing experiments
rather than untouched IVF embryos; the three public series are Cas9 or base-editor studies, so no
count from them is a natural aneuploidy rate; and the laboratory series carries engineered cut
sites at known loci.

**Question 3: does reciprocity transfer to this material, and what fraction of events should be
expected reciprocal rather than shared? A shared event in the same direction points to something
before the split or in the gamete, which is a different biological claim and not a failure of the
method. If the expected reciprocal fraction is low, the prior is available far less often than our
19 of 26 suggests.**

Also worth your view: a CRISPR-induced double-strand break at a known site should produce
reciprocal products more often than spontaneous mis-segregation does, which would make this
corpus unusually favourable to the construction and therefore a poor guide to what it would do
elsewhere. Say if that is backwards.

## 6. Reporting

Our instinct is that a prior-assisted call is a different kind of claim from an unaided one and
must be reported as its own evidence class, with unaided recall quoted unchanged beside it, so
neither we nor a reader can mistake one for the other. The row would carry: the sister cell it
came from, the refined interval and its coverage, the evidence in this cell at that interval, and
the fact that the interval was not independently found here.

**Question 4: is that the right boundary, and what else belongs in the row? In particular, should a
prior-assisted call be allowed to stand alone in a published count, or only to corroborate?**

## 7. The specificity arm, which we will not trade

Fixed before any of this is built: across the 72 arrays the laboratory called nothing on, the
detector raises 1 segment, on 1 array. Any version of this construction that raises that number
fails, whatever it does for recall. A prior supplied by a sister cell is a new way to find a
segment and therefore a new way to find a false one, and the obvious failure mode is that a
dropout-rich interval in cell A nominates the same dropout-rich interval in cell B, since marker
density and probe quality are properties of the array design shared by both.

**Question 5: how should the null for a prior-assisted test be built so that shared marker-level
artefact cannot manufacture a reciprocal pair? This is the question we are least confident about
and the one most likely to invalidate the whole construction.**

## 8. What we will do with each answer

A real gain in detectable size plus a clean null: build it, pre-registered, as a separate evidence
class, and report availability under tool-detected losses rather than hand ones. No real gain, or
no clean null: write it down as attempted and measured, tell our collaborator why her idea does not
reach in this material, and put the time into the marker classes in the companion consult.

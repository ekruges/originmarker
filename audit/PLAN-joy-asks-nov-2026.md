# Plan: Joy's asks from the 7 November call, due Wednesday 11 November

Written 8 November, revised the same day after three consult responses came back. Every repo fact
below was read out of the source or the data, not recalled. Section 4 records what the consults
overturned, including two findings that change what the gains number means.

## 1. What she asked

From the call, in her order:

1. The 0.33 / 0.67 slide. She asked twice what the number IS and how it is computed, and said it
   reads like a B allele frequency and is not one. "You would have to look into the program and see
   how it's done."
2. The dot figure. A maternal loss and a paternal loss are drawn identical except for red dots.
   She said the colours confuse it.
3. The one-page schematic. She asked for it two calls running. The presenter's own verdict was
   that it is not simple enough.
4. A noisy array must say so. Her words: the program should say "I am not making a call not
   because there aren't enough informative markers, but rather the SNP array is too noisy."
5. Gains. 11 of 26 hand-called gains matched. Her two suggestions:
   a. use marker classes beyond the opposite-homozygote one, parent-heterozygous markers included,
      to get enough sites inside a small event;
   b. use a confidently called loss breakpoint in one cell as a prior for the gain in its sister
      cell, since cleavage-stage events are reciprocal.
6. GUARDIAN. Data security, settled before a half-finished product reaches her supervisors.

## 2. What the code does today

Read in `web/src`, commit 47e48d4, version 5.47.0.

- `scanCopyNumber` (segments.ts) finds segments from the NO-CALL run against an external null, then
  keeps one only if `|median log2R deviation| >= SEGMENT_LRR_SHIFT = 1.0`, and the sign labels gain
  or loss. `MIN_SEGMENT_MARKERS = 2_400` called markers.
- A real trisomy is log2(3/2) = 0.585. The 1.0 gate therefore sits above the biology.
- `gainOrigin.callPair` exists and is wired at gainOrigin.ts:699. It weighs the allele share against
  intensity over six states, grids cell fraction and chemistry compression, and has a pre-registered
  held-out split. It runs PER CHROMOSOME and has never been given a segment.
- `paternalShare` (gainOrigin.ts:118) is the quantity on slide 5. It is BAF flipped per marker by the
  parents' genotypes, at opposite-homozygote markers only: father AA / mother BB returns `1 - baf`,
  father BB / mother AA returns `baf`. The region statistic is the MEDIAN of that over informative
  markers, judged against the sample's own genome-wide median, not against 0.50.
- The allele channel uses opposite-homozygote markers only, with a floor of 200 informative markers.
- 5.47.0 splits a refusal into "the array failed" and "this chromosome offered nothing", by call
  rate, in `obligateHet.whyNotInformative`. No array-noise statistic exists in the shipped path.

## 3. The numbers that set the targets

From `audit/local/features/regions-dieter.tsv` and the 8 November baseline run:

- 26 hand-called gains, 67 losses, over 21 cells of 93. The other 72 arrays carry no hand call and
  are the false-call arm.
- Gain spans, Mb: min 3, Q1 14, median 29, Q3 69, max 245.
- Baseline at 0.5 reciprocal overlap: gains 11 of 26, losses 44 of 67, 1 false segment on 1 of the
  72 clean arrays. Recorded in `audit/joy-dieter-baseline-2026-11-08.txt`.
- The 15 gain misses by cause: 4 under the 1.0 shift gate, 6 boundary disagreement, 2 under the
  2,400-marker floor, 3 unexplained.
- Opposite-homozygote marker density is 11.1 per Mb with both parents called, not the 8 per Mb the
  first plan assumed. Corrected by consult A against GSE148488.

## 4. What the consults answered, and what they overturned

Three consults went out on 8 November and all three came back the same day. Responses are in
`audit/CONSULT-RESPONSE-segmental-gain-statistic.md`, `-sister-cell-prior.md` and
`-array-noise-refusal.md`. What follows is only what changes the work.

### Two findings that block any gains number reaching her

- **The reference table may not be hand calls.** All 93 rows of `regions-dieter.tsv` are an exact
  subset of `regions-dieter-sensitive.tsv`, `log2R_dev` identical to 4 dp, 100% on a 1 Mb grid, and
  `STATE.txt` lists both as offline script output. Filtering the sensitive table at |dev| >= 0.815
  and >= 3 Mb reproduces it with 0 missed and 6 extra rows. Either the hand calls were transcribed
  onto the script's segments, or "recall against his calls" is the tool scored against itself. Not
  resolvable from the repo. Ask Dieter.
- **The 11 of 26 may be agreement on distorted material.** 10 of the 11 found gains sit in cells
  with zero hand-called losses, at deviation +1.61 to +2.83, which is the amplification-distorted
  range rather than a third copy. 12 of the 13 gains in cells that do carry losses are missed.
  Testable with one number per hand gain: homozygous-pole fraction inside the region against that
  array's median. Prediction: found gains below 0.5, missed gains near 0.83.

### What is dead

- **Lowering the shift gate.** Real gains measure +0.49, +0.27 and +0.33 on a measured trio. The
  clean-window null at 10 Mb has p99.9 = 0.30 and a maximum of 0.53. Single-cell log2R cannot carry
  a gain gate at anything near 1-false-in-72. Allele channel first, intensity second.
- **The sister prior as specified.** Fixing the sister's interval rescues 0 of the 15 misses and
  buys 4 to 17 nats against a bar of 250. The 19 of 26 measures chromosome co-occurrence, not
  sisterhood: non-sister cells overlap more often at r >= 0.5 (12 of 26) than same-prefix ones
  (8 of 26). The barcode prefix is 38 and 32 arrays, not an embryo.
- **hetBafSd as a noise gate.** The 29 failed arrays read 0.159 to 0.170, below most usable
  amplified ones at 0.102 to 0.299. Any upper cut refuses good cells and passes every failure.

### What replaces them

- **Ask 4 is solved.** The statistic is the information ratio R: the sample's opposite-homozygote
  rate against a reference, over the rate expected if its calls carried no information. Parent-free,
  carries its own null. Informative arrays read <= 0.26, failed arrays >= 0.72, with an empty gap
  from 0.424 to 0.908 across 135 arrays. Cuts 0.30 and 0.60. It changes reason text only and never
  the predicate, so no chromosome moves between called and refused.
- **Joy's 5a is right in information terms and not yet usable.** Parent-heterozygous markers do
  carry dosage signal (KL 0.08 to 0.14 per marker), but their null is misspecified about 3x beyond
  what a variance-inflation correction captures, so adding them lowers the derived floor and raises
  the empirical one. It is a null-calibration problem, not a marker-supply problem. Both parents must
  be genotyped, and no parental array for these cells exists locally.
- **Joy's 5b survives in a different form.** A hemizygous sister supplies PHASE over exactly the
  interval in question. Orienting the sibling's BAF by the sister's retained allele raises
  per-marker information 7 to 12x and needs no parent. It also comes with a free null: score the
  sibling for a gain of the sister's RETAINED haplotype, which artefact cannot tell apart from the
  real orientation but a real reciprocal gain can.
- **Reciprocity is rarer than assumed.** Turocy 2026, this lab's own Cas9 series, finds 22 reciprocal
  segmental gains against 96 losses, 22.9%. Expect a sister gain for about 1 loss in 4 or 5, not 19
  of 26. Breakpoint transfer is side-specific: the acentric edge is conserved to kb, the centric edge
  is moved by attrition.

### Code defects the consults found

- `callPair`'s compression grid runs 0.55 to 1.00 and excludes both corpora. Measured losses read
  k = 0.475 and 0.535; the hand losses here read k about 1.43.
- `MAX_HET_BAF_SD = 0.11` is documented as adopted from MoChA. It is in neither MoChA README, whose
  filters are call_rate < 0.97 and baf_auto > 0.03. Provenance unverified.
- The shipped refusal says "a rerun is what fixes it". On the only data that test it, replicate
  re-hybridisations of one amplification product fail every time, and the group we send off to rerun
  reproduced 8 of 8 collapsed chromosomes. The wording is false.
- `ALLELE_GAIN_MIN_INFORMATIVE = 200` is not the binding constraint; the derived floors are 56 to
  175 markers.

### Defects found in the deck while reading the code for ask 1

- Slide 5 names the quantity "B allele frequency, averaged over about a thousand markers". It is the
  paternal share, it is a median, and it is judged against the sample's own centre.
- Slide 5's two gain labels are swapped. `gainOrigin.check.ts` asserts share 0.6664 is **paternal**
  and 0.3331 is **maternal**; the slide has 0.33 as paternal and 0.67 as maternal.
- The dot figure's own caption says "The last two plots are identical. The red is the only thing
  that tells them apart." On a raw allele-frequency axis that is true, which is the defect: the
  figure never draws the axis the program actually uses. It also reds only one of the two
  informative classes, leaving the other drawn identically in black.

## 5. Tracks

Each track states its hypothesis, its fixed bar and its stop rule. Bars are fixed before scoring.

### Track A. Unblock the gains question (today)

A1. Draft the provenance question to Dieter: is `regions-dieter.tsv` his reading, or script output
he reviewed? Nothing in ask 5 is quotable until this is answered.
A2. Pole-fraction harness: one number per hand gain, homozygous-pole fraction inside the region
against that array's median. Runs on tesla. Bar fixed in advance: if found gains read below 0.5 and
missed gains near 0.83, the honest baseline is recall on pole-intact gains and the headline number
changes. Stop rule: report the split either way, including if the prediction fails.

### Track B. The three figures (today and Monday)

B1. Slide 5 rebuilt: name the quantity, show the per-marker flip, correct the swapped labels, say
median and say centre. This is ask 1.
B2. The dot figure rebuilt on the paternal-share axis, where a maternal loss piles at 1.0 and a
paternal loss at 0.0 and no colour is needed to tell them apart. Keep one raw-BAF panel as the
"before" so the point of the flip is visible. This is ask 2, and it is the same fix as ask 1.
B3. The schematic cut from seven boxes to three: pick the markers, read two channels, combine or
refuse. Every threshold leaves the figure. This is ask 3.

### Track C. Ask 4, the noise refusal (Monday)

C1. `audit/noise-information-ratio.ts`: per-array and per-autosome R beside call rate, het rate and
hetBafSd, using the shipped parser, naming no sample.
C2. Run the 459 on tesla against the pre-registered predictions: P1 every array at R >= 0.60 is
already fully withheld, P2 arrays in 0.30 to 0.60 are <= 2% and all under 0.50 genome call rate,
P3 their hetBafSd sits below the informative median, P4 any informative HumanCytoSNP-12 array over
0.30 gets the number reported and no gate.
C3. Wire the reason text, per sample, with the per-chromosome reasons kept underneath. Count and
list every chromosome whose reason text changes; zero may change verdict, and a check pins that.
C4. Correct the shipped "a rerun is what fixes it" wording at the same time. It is wrong today and
it is wrong independently of the new branch.

### Track D. Code defects (Tuesday)

D1. Extend `callPair`'s compression grid to cover 0.4 to 1.6, or fit k per cell. Re-run the
pre-registered held-out split; the bar is that the held-out direction error does not get worse.
D2. Trace `MAX_HET_BAF_SD = 0.11` or restate its docstring to say it is a project constant with an
unverified attribution.

### Track E. Rescoring, honestly labelled (Tuesday)

E1. chr17 in two S13 cells is split by the reference into three pieces against one 0 to 73 Mb loss
in a third. Merged at a 2 Mb gap, one of them scores r = 0.95 against that loss. Rescore gains on a
merged event unit against the detector's own coordinates. Bar: report the change as a RESCORING in
its own row, never folded into recall.

### Track F. GUARDIAN (Monday)

F1. One page: what data the tool holds, where it goes, what never leaves the machine, and the
barcode pseudonymisation rule. Her condition is that this is settled before her supervisors see
anything.

### Track G. Her packet (Tuesday into Wednesday)

G1. The corrected scoreboard, with the losses correction: 44 of 67 at the 0.5 bar, and the 78 of 85
she heard is monosomy parent-naming on a different series and a different question.
G2. The three rebuilt figures.
G3. A one-page note on ask 5 saying what was tried, what was measured and what was not won.

## 6. Order of work

- Sunday: A1 draft, A2 launched on tesla, B1, B2, B3.
- Monday: C1, C2 launched, F1, C4.
- Tuesday: C3, D1, D2, E1, G1, G2, G3.
- Wednesday: deliver.

## 7. What is certain by Wednesday, and what is not

Certain: asks 1, 2, 3, 4 and 6. All five are either a figure rebuild, a writing task, or a reason
branch whose statistic is already derived and whose cuts are already fixed.

Not certain: ask 5 as a recall improvement. No honest version of it lands by Wednesday. What lands
is the diagnosis, the corrected baseline, the merged-event rescoring labelled as a rescoring, and a
written statement of what the two suggestions are worth. Joy gets told that plainly rather than
shown a number that moved for the wrong reason.

## 8. Rules that hold throughout

- No bare laboratory barcodes in tracked files. Per-array tables stay under `audit/local/`.
- Do not re-derive analysis outside the shipped modules; harnesses import `web/src` directly.
- The 72 clean arrays are the specificity arm and are never traded for recall.
- Bars are fixed and written down before scoring, and a failed prediction is reported as a failure.
- Every number that reaches Joy traces to a run, a file or a paper.

# Changelog

Releases are named after the biology of crossing over. All dates 2026.

The entries below are unusually detailed about bugs. That is deliberate: this tool's output
informs which parental allele an embryo inherited, so a plausible wrong answer is worse than
an error. An error gets investigated; a plausible wrong answer gets used. Anyone deciding
whether to trust a panel from an older build deserves to know exactly what it got wrong.

---

## 5.2.1 "Interference"

The taxonomy's findings are now scored for parental origin through the same path as everything
else, which is what 5.2.0 claimed and did not do.

A copy-neutral event and an isodisomy carry a parental origin exactly as a deletion does, but they
arrived in the defect list with the origin unresolved because nothing had scored them. The dosage
scorer was chromosome-scoped and inline; it now takes the interval as a predicate and is used
twice, so a finding is measured over its OWN interval rather than borrowing its chromosome's
answer, and gets the same posterior, the same four bands and the same class-inversion veto.

Classes whose origin is blocked by the class itself are skipped rather than scored and discarded,
and a stray call arriving on one is ignored rather than trusted: a triploidy must not acquire a
parent from a number that means nothing for it.

Two log lines compared a dosage verdict against 'refused', which belongs to the genotype channel
and never appears there, so both reported success for every outcome including a withheld parent.
The trisomy mechanism gate is now also conditioned on the interval being a whole chromosome, since
copy number three is a whole-chromosome property and asking it of a sub-chromosomal interval
returns the by-exclusion answer every time.

---

## 5.2.0 "Interference"

The full taxonomy of chromosomal abnormality, including the classes that cannot be answered, with
every new class wired into the same machinery the old ones use.

**What was wrong before.** This tool found three things: a parental absence, a copy loss and a copy
gain. Everything else an embryologist looks at was not looked for, and silence read as absence. A
reader had no way to separate "this array has no uniparental disomy" from "this tool has never
checked for uniparental disomy". Fourteen classes are now enumerated, each saying what it rests on
and what a parental origin for it rests on.

**Four classes are newly detected, and need no parental array at all.** Copy-neutral loss of
heterozygosity, from heterozygosity depletion against the array's own background with copy number
unchanged. Uniparental isodisomy and segmental uniparental disomy, from runs of homozygosity past
the 13.5 Mb reporting length. Triploidy, from allele fractions occupying the one-third and
two-thirds bands while vacating the half band, which a diploid genome cannot do. And a genome too
disturbed to self-reference, which is reported as a finding in its own right rather than as a
refusal: the array may be perfect and the genome is what is disturbed.

**The copy-neutral detector will not fire on a deletion.** Depleted heterozygosity is equally the
signature of a lost copy, and calling that copy-neutral would name the wrong event and then invert
its origin, since loss and gain read opposite signs. A window whose intensity has moved belongs to
the deletion detector and is left to it. Where no intensity is supplied the class is reported as
loss-or-copy-neutral rather than assumed. The calling threshold also sits clear of the false-positive
floor: dropout that clusters rather than scattering already puts 1% of event-free windows at 40%
depletion.

**No run found is not no disomy, and the negative says so.** Runs find isodisomy, where both copies
are the same homologue. Heterodisomy leaves copy number at two and heterozygosity normal and is
invisible from an embryo alone, and 32% of confirmed uniparental disomy carries no significant run
at all. Both facts travel with every result rather than sitting in documentation. Runs across four
or more chromosomes are the pattern of shared parental ancestry, and those calls are flagged rather
than withdrawn: the runs are real, what they imply is what changes.

**Three limits are named so nobody chases them.** Uniparental heterodisomy needs both parents
alongside the embryo. Tandem versus inserted duplication carries no positional information in
either channel. Reverse segregation is a blind spot rather than a power problem: in 20 of 26
observed cases the copy number was normal, and it was the most frequent non-canonical pattern
across 23 complete meioses, so no marker density fixes it. Two further walls are conditional and
lift when the user supplies more, and the report lists which still stand for the run actually done.

**Everything goes through one list, not a second panel.** The new classes enter the same defect
display, the same headline with its confidence band, the same PDF, with a readable phrase for each
so a clinician need not know the code's vocabulary. Where a class carries no parental origin at all,
that is stated as a property of the class rather than left as an empty field, because a structural
impossibility and a failed attempt must not look alike. `om taxonomy` prints the whole table.

---

## 5.1.0 "Disjunction"

Every origin channel now emits a confidence and a band, not just the dosage one. Measuring what
those numbers are worth turned up a defect in the two channels that were already shipping.

**The gap this closes.** 5.0.0 gave the dosage channel a calibrated posterior and left the other
two alone, so a two-parent genotype call, which is the strongest evidence this tool has, printed a
bare parent name while a weak dosage call printed 0.62 beside it. The formatting told a reader the
opposite of the truth. All three channels now carry a number, a band and a reason, through the
defect list, the callout and all three PDF tables.

**The defect that found.** Both genotype channels return a posterior that does not vary. Measured
on the shipped obligate-het channel: 1.0000 at 0, 2, 8, 40, 120 and 300 exclusive markers out of
400. It flips its verdict between 8 and 40 and never moves its number. The two-parent share model
behaves the same way once given a realistic error, reading 0.9974 at the calling margin and 1.0000
everywhere above it. A number that does not vary is not a confidence.

This is not a fault in either likelihood. It is what a ratio over hundreds of near-independent
Mendelian markers does when taken at face value. The real error rate on those channels is set by
things no likelihood here represents: contamination, a mis-specified dropout rate, a sample that is
not the one on the label. So on those two channels the BAND is the output and the digits are not,
both are capped below the top band, and both report themselves uncalibrated. Only the dosage
posterior, which marginalises its nuisance parameters rather than plugging them in and floors its
error with drift that does not average down, earned its bands from measurement.

**The standard error on the two-parent channel is floored.** Sampling error alone is about 0.002
over the marker counts this channel sees, which would make any deviation clearing the 0.056 margin
read as certainty. The floor is the margin itself, the project's own long-standing estimate of
where this statistic stops being trustworthy, with sampling added in quadrature beneath it. This is
the same lesson the dosage channel learned from drift.

**What is still outstanding.** Giving the genotype channels a confidence that discriminates needs
the treatment the dosage channel got: marginalise the nuisance parameters, and validate against a
parent-child truth set. No such truth set exists in the corpus to hand. The families its filenames
suggest are not families under test, one nominal parent reading 55% heterozygous against a
diploid's 16.8%, so no pair from it can serve.

Detection is unchanged in this release. Nothing new is found; what is found is now reported with a
confidence on every channel.

---

## 5.0.1 "Anaphase"

The PDF said the dosage table carried a confidence column and it did not. The column is now there,
printing the calibrated number and its band, grey for bands C and D so a weak number cannot be read
as a strong one. A parent printed without its confidence is the exact failure 5.0.0 set out to
remove, and the report was still doing it.

---

## 5.0.0 "Anaphase"

Parental origin is no longer read from the sign of the allele-dosage shift. It is a calibrated
posterior that marginalises the copy-number class, and it comes with a confidence number and a
band on every call.

**The bug this fixes, stated plainly.** The dosage channel named a parent by taking the sign of the
self-referenced centroid shift, inverting it through the loss formula whenever the copy-number
class was unresolved. That default was not conservative. Gain inverts the map that loss and
copy-neutral LOH share: at a mosaic fraction of 0.10 the loaded parent's copy reads +0.0263 under a
loss and -0.0238 under a gain, so the same positive shift names opposite parents depending on a
class the tool usually could not resolve. On trophectoderm and blastomere material 89 to 100% of
detected events resolve an origin without a class, so the assumption was untested in the majority
regime rather than at the edges. Measured wrong-parent rate for a true low-fraction gain scored
that way: 0.551 to 0.580 at f = 0.05, across four material classes. Worse than chance, and
systematically so, because a sign inversion is not noise. Anyone who read a parental origin for a
gain or a possible gain from build 4.x should re-run it.

**What replaces it.** The class and the mosaic fraction are both marginalised into a posterior
probability that the affected copy came from the un-genotyped parent. The intensity channel feeds
the class and never the origin: at a fixed class and fraction its term is identical under both
parental hypotheses and cancels exactly. Where intensity resolves the class the posterior sharpens;
where it does not, it hedges across the inverting map and lands near 0.5, which is the honest
answer instead of a confident wrong one.

**Four bands, and every one of them carries its number.** A very confident, B confident, C weak and
direction only, D weak and not for reporting. Band D is not a coin flip: measured 0.604 to 0.636
with array-clustered intervals that all exclude 0.50, and calibrated within itself to within 1.2
points, so suppressing its number would hide a calibration that can be demonstrated. What changes
across bands is the words beside the number, not whether a number appears. Bands C and D render
dimmed so a weak number cannot be mistaken for a strong one at a glance.

**Two thresholds retired.** The 0.30 sign-security bound is gone: it was calibrated against
wrong-parent rates that were themselves produced by assuming the class, and on blastomere material
it refused 24.1% of events that sit in a band measured at 0.9971. A single threshold on fraction
cannot express this in any case, since fraction explains only 51.7% of the variance in achievable
confidence. The detection floor also no longer gates assignment. By the time a parent is named a
detection has already been made at z >= 2.576, and re-applying a power threshold refused events the
array did in fact see. The floor stays where it belongs, on the question of whether any array of
this kind at this width could answer, and as reported information on the call.

**One cell still withholds a parent.** Amplified material, an implied gain fraction under 0.15, and
an intensity channel that does not resolve the direction at 99%. There a small loss and a small
gain both explain the observation and name opposite parents, and the measured accuracy of a stated
parent is 0 of 34 on trophectoderm and 0 of 18 on blastomere. The event, its interval and its
evidence are all still reported. Only the parent is withheld, with the remedy stated: a second
parental array resolves the class categorically rather than by threshold.

**Calibration, measured here rather than asserted.** The recalibration maps this design calls for
were not available, and could not be honestly refitted, because the corpus carries no material
labels and the staging function separates ploidy rather than amplification. So the question that
decides whether the posterior is safe without them was measured directly: 140,000 injections
carrying the real per-chromosome noise of 35 arrays. Expected calibration error 0.0039, erring
under-confident in the middle bands by 2 to 3 points, which is the conservative direction. Band A
reads 0.99955 over 64,522 rows with no intensity channel at all. Every posterior the tool emits
until the maps arrive says in its own reason text that it is not recalibrated. Method, numbers and
limits in `audit/calibration/FINDINGS.txt`.

That measurement also reproduced the one defect independently. With no intensity supplied, true
gains reaching the top bands are correct 0 of 60. The algebraic gate was then tested and does not
work: a gain cannot displace more than 0.1296 at the fraction ceiling, but noise carries true gains
above it, so gating there demotes zero rows and catches zero errors. It is an identifiability limit
rather than a tuning failure, and it is now confirmed rather than assumed.

**Breaking.** The dosage channel's verdicts are `loaded-parent` and `other-parent`, replacing
`known-parent-lost` and `other-parent-lost`. The old names presumed a loss, which is exactly what
that channel cannot presume. The obligate-het channel keeps its own vocabulary, where a missing
haplotype really is missing and no gain can invert it. A dosage verdict in the retired vocabulary
is no longer honoured rather than being silently read as an answer. The `--sign-secure-f` dial is
removed with the threshold it controlled.

**Trisomy mechanism.** `SPH` is renamed `not-BPH`. Meiosis II non-disjunction with a crossover
produces both-homologue tracts distal to the centromere, so the both-versus-single homologue
distinction is not meiotic versus mitotic, and the second category pools non-recombinant meiosis II
with post-zygotic duplication. The call now states what it pools instead of inviting a reader to
hear a mechanism the data does not choose.

**Release names.** The ladder in `build_info.py` had been spent, and 4.14.0 and 4.15.0 took
"Kinetochore" and "Recombinase" back from 2.0.0 and 3.0.0. Six new names are added and
`release-check.sh` now fails a build that reuses one.

---

## 4.15.0 "Recombinase"

Two channels built on markers the tool has always discarded.

### Added

- **The untransmitted-haplotype channel.** Every origin call this project has made used markers
  where the loaded parent is HOMOZYGOUS. `untransmitted.ts` uses the complement, a disjoint set of
  comparable size that was simply unused: 27,302 to 58,571 markers per array. Where that parent is
  heterozygous and the sample reads homozygous, the transmission is determined, so the parent gave
  the allele shown and withheld the other.

  The advantage is not a cleverer statistic. The marker set is CLEAN BY CONSTRUCTION: in the
  parent-homozygous window only 32 to 90% of markers are truly heterozygous in the child, because
  the other parent often transmits the same allele, and those carry no information while still
  counting toward the denominator. Here every marker carries information. Measured 1.40 to 1.98x
  signal-to-noise above the obligate-het channel, on every amplified material at every fraction.

  **It is the only channel that gives a single blastomere a floor at all.** Against one parent that
  case has none at any fraction to 0.70, and still none at eight cells. Through this one it has
  0.628. That is above the callable bound, so this converts a material impossibility into a
  measurable quantity rather than into an answer, which is a different sentence to a reader.

  It also yields a per-array error rate needing no external truth: readings whose dosage sits at
  the extreme opposite to their own genotype call, measured 1.7 to 13.8% on this platform.

- **Trisomy mechanism from band occupancy.** Both parental homologues puts the other parent's share
  at 1/3 in one band; a duplicated single homologue puts it at 1/3 or 1.0 in two, neither at 2/3.
  The diagnostic is OCCUPANCY rather than a mean, which is what makes it work here: a per-marker
  three-band assignment needs BAF spread under 0.053 and misassigns 13 to 22% on every amplified
  class, while a fractional occupancy over hundreds of markers is insensitive to per-marker
  dispersion. AUC 1.000 from 400 markers.

  Reported positively for both-homologues only. The single-homologue bands are also populated by an
  ordinary euploid genome, so occupancy there has power 0.134 against euploid and means nothing
  alone.

### Fixed

- **The mechanism call is gated on copy number three being established**, and the gate was added
  because the ungated version got it wrong on real data: a euploid chromosome of a confirmed
  parent-child pair came back "SPH", confidently, by exclusion. It answers that for every normal
  chromosome, because a euploid genome populates both of the bands it reasons from. Asking the
  question at all requires the copy number first.

---

## 4.14.0 "Kinetochore"

Fewer refusals, by asking the right questions separately instead of one question that fails as a
whole. Twelve of twenty-four material-state-parent combinations are now callable where the module
previously recognised four.

### Changed

- **Origin and class are separate verdicts.** This is the largest single source of refusals removed.
  Power to detect an allelic imbalance exceeds power to resolve which state produced it, and on the
  material this tool targets the gap is total: at 400 informative markers 100% of detected events on
  trophectoderm and blastomere resolve an ORIGIN and none resolves a CLASS. A caller emitting one
  verdict had to refuse every one of them, discarding an origin it could support. `classVerdict` and
  `classWhy` are now their own fields, and an unresolvable class no longer takes the origin with it.
  Lifting the class needs a three to fivefold narrower window log2R spread, measured at 0.17-0.22
  against the 0.029-0.081 required, and four times the markers buys 1.2-1.4x. It is not reachable by
  collecting more of the same, so it is reported rather than chased.

- **Floors are state-aware and parent-count aware.** The module shipped with loss-only floors and so
  refused the class it handles BEST. Copy-neutral LOH is the largest-signal state at every material,
  because copy number stays at 2 and the genotype caller never degrades: median call rate 0.867 at
  copy-neutral log2R against 0.287 below -1.5. On a trophectoderm biopsy with ONE parent its origin
  floor is 0.186, comfortably callable, where a loss on the same array is 0.625.

- **A second parent is now modelled as an input, not an aspiration.** On the same single file it
  moves a TE whole chromosome from 0.625 to 0.186, a factor of 3.36, where one cell to five moves
  0.186 to 0.135, a factor of 1.38. The app passes the count automatically; the CLI takes
  `--parents 2`.

- **A trophectoderm loss reports an unassigned imbalance rather than nothing.** It has a floor of
  0.625, too high for the sign to be secure but not absent, so the event is reported without a
  parent instead of being refused as unevaluable. Only one combination has no floor at all and it
  stays refused: a single blastomere loss against one parent, at any fraction to 0.70 and still none
  at eight cells. That is a material limit and no file the user can load moves it.

---

## 4.13.3 "Cohesin"

### Fixed

- **The mosaic-fraction inversion was the loss formula applied to every state.** A methods review
  derived the full algebra and verified it by brute-force copy counting to 1.1e-16. The three
  inversions differ: f = 4d/(1+2d) for a loss, 4d/(1-2d) for a gain, 2d for copy-neutral LOH. At
  d = 0.04 those are 0.148, 0.174 and 0.080. Sharing the loss form was correct on losses and wrong
  in the OPPOSITE direction on gains, which is the worse way to be wrong: the error inflates the
  fraction rather than shrinking it. `fractionFromShift` now takes the state.

### Corrected

- **Copy-neutral LOH is the easiest class to assign a parent to, not the hardest.** Rank order by
  deviation is CNN-LOH > loss > gain at every mosaic fraction, because copy number stays at 2 so
  the genotype caller does not degrade: median call rate 0.867 at copy-neutral log2R against 0.287
  below -1.5. The intuition that a class with no copy-number signal must be hardest is backwards.
  Detection is the hard half there, 1.75 to 17.1x harder than assignment.

- **"A gain shifts about a third as much as a loss" is the f = 1 endpoint.** The pooled ratio is
  (2-f)/(2+f), which is 0.905 at f = 0.1 and 0.818 at f = 0.2, so the mosaic-range penalty is 1.1
  to 1.7x rather than 3x.

- **"The signal that finds the event destroys the evidence that assigns it" is half wrong.** At a
  vendor no-call the discrete genotype carries 0.0000 bits about origin but the continuous BAF
  carries 0.0471, and BAF is present at 88-91% of them; the cause is low intensity, not an
  off-cluster reading. Discretisation throws the information away, not the event. The refusal on
  whole-chromosome losses survives for a different reason: in this dataset a loss large enough to
  see in one cell co-occurs with an array too damaged to self-reference, and the two are perfectly
  confounded. All 70 whole-chromosome losses sit on arrays with 40-100% of autosomes deviant.

- The review is committed at `audit/ORIGIN-CONSULT.md` with `boundary_of_certainty.csv` and
  `refusal_taxonomy.csv`, so the 22 verdicts and 14 refusal conditions can be read rather than
  paraphrased.

---

## 4.13.2 "Cohesin"

### Added

- **Continuous integration.** Dependabot has been opening grouped pull requests against nothing, so
  the merge-if-green design in its config was unavailable and every bump was checked by hand. Three
  jobs, so a failure names itself: build and lint, the self-checks plus the bundled examples, and
  release metadata.

- **`scripts/checks.sh`**, which runs all 31 self-checks each from the directory it needs. A
  `for f in src/*.check.ts` from the repo root reports PanelTable and PrimerOptions as failing,
  because they resolve their component through Vite at `/src/...` against the working directory.
  The failure is a false alarm that costs a re-run to establish, and a runner that cries wolf
  teaches you to skim its output.

- **`scripts/release-check.sh`**, asserting the version the build reports, the newest changelog
  entry, the codename and the citation table all agree. It found five releases with no changelog
  entry at all on its first run: 4.12.1 through 4.13.1 were committed and deployed while the
  edits meant to document them silently matched nothing. Those entries are now written.

- **`scripts/deploy.sh`**, replacing six steps run by hand every release, with the health check as
  a gate that exits non-zero rather than output someone has to read. A release once shipped
  reporting the previous version and was caught only because the output happened to be read.

- **`cli/examples.check.ts`**, the only check here that guards an answer rather than a build. Each
  bundled example carries an `expect` string making specific numeric claims, and nothing verified
  them, so any threshold moved in this project could falsify one silently. It runs the shipped
  modules over the files a user actually downloads and asserts the opposite-homozygote rates,
  second parental contributions and inferred dropouts still match the prose.

---

## 4.13.1 "Cohesin"

### Fixed

- **The per-locus deletion test did nothing when only an oocyte donor was loaded.** Its form was
  gated on the oocyte alone, which was correct until 4.8.3 made oocyte-only runs possible. After
  that a run with only an oocyte rendered the form and enabled the Run button while the handler
  returned immediately for want of the sperm donor: no result, no error, no spinner. The gate asks
  for both parents now and names whichever is missing. A control that is live in a state where it
  cannot act is worse than one that is absent. The test itself was never broken; it needs a marker
  where one parent could not have supplied the allele, so it needs two people.

---

## 4.13.0 "Cohesin"

### Added

- **A second family in the bundled examples**, found by searching the public deposits rather than
  reusing the series already in hand. A donor and two of its children from GSE290961, on the same
  Axiom platform. The pair is the point: the same relationship read from bulk DNA and from a single
  amplified cell. One links at 1.11% opposite-homozygote with a 9.0% second parental contribution
  and reads dropout 0.008; the other at 1.82% with 8.7% and reads 0.193. Nothing is declared for
  either. No example showed that contrast before.

- **The subset's cost is stated where the examples are defined.** Linkage, stage inference and
  whole-chromosome events survive every-eighth-marker sampling because each uses a whole genome or
  a whole chromosome. A segmental change does not: the resolution floor is 2,400 informative
  markers and a 12 Mb region holds about 390 at that density.

### Searched and not found

No sample carrying a chromosomal change attributable to a parent. GSE186407 is 166 arrays with no
parental DNA. GSE290961 has 158 arrays, nine bulk-quality candidates, six with confirmed children,
and 22 arrays carrying a whole-chromosome event. The two sets do not intersect, and not by
accident: every array with an event has a call rate between 0.587 and 0.868 and a BAF spread
between 0.219 and 0.300, against 0.061 to 0.085 for the bulk arrays. In this material a detected
whole-chromosome loss travels with a collapsed call rate.

Provenance corrected: the examples header claimed every file came from one series, which three no
longer do. The second is cited by accession, no publication having been confirmed for it.

---

## 4.12.3 "Tetrad"

### Added

- **A Syngamy example that passes.** Every example was a refusal or a negative, so the set read as
  though the tool only ever says no. GSM4472424 is a trophectoderm biopsy of the donor's own child,
  confirmed at a 0.66% opposite-homozygote rate with a second contribution at 14.7% of 74,399
  markers and no chromosomal change. It also carries the intensity column, which none of the
  others do: without it the copy-number channel is dead on every bundled example.

---

## 4.12.2 "Tetrad"

### Added

- **Both origin channels in the report.** A genotype table carrying the exclusive-marker count,
  which is the one number checkable without trusting the model, and a dosage table for whole
  chromosomes with the shift, its z, the implied fraction and the central window. Separate tables,
  because the two answer the same question from different measurements and a reader needs to know
  which one answered.

- **A documentation section on the dosage channel**: why self-referencing is the method rather than
  a refinement, why drift and not marker count sets the uncertainty, the order the four questions
  are asked in, and why no parent is named below a mosaic fraction of 0.30.

### Fixed

- **The documentation nav scrolls inside itself.** At 35 sections the list outgrew the viewport,
  and a sticky element taller than the screen hides its own tail, so the last sections could not be
  reached without scrolling the article. Mobile opts out, being a full-width header there.

---

## 4.12.1 "Tetrad"

### Added

- **The measured correlation between the dosage and intensity channels, and a joint term using it.**
  The review fitted this on technical-replicate differences (-0.055 to +0.036) and means (up to
  -0.46); neither applies to the shipped statistic, which is self-referenced, and that subtraction
  removes exactly the artefact separating those two figures. Measured on the shipped quantity over
  81 arrays: bulk -0.058, blastomere +0.486, esc-single +0.508, trophectoderm +0.633. Bulk is
  independent and every amplified material is strongly correlated, the opposite of what quadrature
  assumes, which would overstate the joint z by 1.27x on trophectoderm.

  Intensity enters DETECTION and never DIRECTION. On haploid products of a known parent log2R
  cannot distinguish maternal from paternal (p = 0.54) while oriented dosage separates them at
  1.3e-15, so the joint term decides whether an event is present and the dosage sign alone decides
  whose it is.

---

## 4.12.0 "Tetrad"

The dosage channel rebuilt to the methods review's specification. 4.11.0 stopped it naming a
parent; this is the statistic that earns the right back, and mostly declines to use it.

### Changed

- **The statistic is a self-referenced centroid shift, not a band-occupancy likelihood.** Orient so
  the loaded parent's own allele is low, keep markers where that parent is homozygous, take the
  central window 0.20 to 0.80, and subtract the same quantity computed on the array's OWN other
  chromosomes. Self-referencing is what removes the directional null: the raw one-parent null sits
  at -0.031 on trophectoderm and -0.023 on blastomeres under no event, pointing at the parent that
  was not genotyped, and on TE that offset is the shift a real mosaic fraction of 0.117 would
  produce. Subtracting the array's own genome brings the null median to -0.001..+0.006.

- **The standard error is floored by within-array drift, not by sampling.** Drift does not average
  down with more markers: against a sampling term of 0.0019-0.0102 it runs 1.44-1.90x larger on
  amplified material. Omitting it is how a z of 15 arrives on a chromosome independently verified
  diploid. Variance inflation from spatial correlation is applied per material (bulk 0.94,
  blastomere 3.63 per chromosome), because effective independent markers saturate near 250 per
  chromosome out of a nominal 900-1,100 and adding markers past that buys nothing.

- **The mosaic fraction is inverted correctly**, f = 4d/(1+2d) from E[BAF] = 1/(2-f), not f = 2d
  from the per-cell form, which understates it about 1.8x.

- **The guard's two jobs are split.** Array quality is a property of the array, asked once, using
  MoChA's gate of BAF spread above 0.11 at heterozygous sites. Evidence is the centroid shift. The
  old single guard conflated them and keyed on a proxy for call rate.

- **Evaluability is decided before the data is read**, from measured detection floors: a 12 Mb
  interval is out of reach on every amplified material, and a blastomere against one parent has no
  floor at any fraction to 0.70. Those return `not-evaluable` rather than a refusal, and the
  distinction is the point. A refusal says this array was bad; not-evaluable says no array of this
  kind at this width could have answered, which sends a reader to change the design rather than
  re-run a sample. That question is asked FIRST for the same reason: nearly every amplified array
  fails the MoChA gate, so asking quality first would report a QC failure for a study-design limit.

- **A parent is named only above an implied fraction of 0.30.** Below it, between half and all
  detections name the wrong parent. Between detection and that bound the verdict is
  `imbalance-unassigned`: the event is reported, the class is withheld.

### Verified against the arrays that produced the wrong answers

- The pronucleus chr1 that once returned both-copies-present at posterior 1.0000 is now
  `array-excluded`, on BAF spread of 0.256 at heterozygous sites against the 0.11 gate, which is a
  property of the array rather than a proxy for it.
- The two chr16 trophectoderm biopsies now read `not-evaluable`, naming the interval width and the
  material, which is the output the review specified for them.

---

## 4.11.0 "Diakinesis"

The dosage channel stops naming a parent. A methods review measured its null and the result is
disqualifying as written.

### Fixed

- **The one-parent null is off-centre, and it points at the parent that was NOT genotyped.** Raw
  centroid null medians are -0.031 on trophectoderm and -0.023 on blastomeres under NO event at
  all. On TE that offset is the shift a genuine mosaic fraction of 0.117 would produce. The
  mechanism is the orientation step itself: the loaded parent's allele is the one always
  identifiable, so dropout of the other parent's allele moves readings toward it more often than
  the reverse. This module's check asserts symmetry between the two DIRECTIONS of the statistic and
  that assertion holds; it does not protect against a null that is off-centre to begin with, which
  is a distinction worth keeping in mind about symmetry tests generally.

- **And the error is confident.** On the honest null the wrong-sign tail reaches z = -15.9 at the
  0.1st percentile, with 0.6% of null units beyond |z| = 10. Below a mosaic fraction of 0.3, the
  proportion of detections naming the WRONG parent runs 0.50 to 1.00. Posterior 1.0000 on noise is
  not a bug that was fixed in 4.10.0; it is reachable whenever the statistic is not self-referenced
  against the array's own clean genome, which this one is not.

  The channel now reports the EVENT and withholds the CLASS. That is the published response to this
  exact dilemma: MoChA left 29% of its events unassigned because power to detect an imbalance
  exceeded power to resolve which imbalance it was.

### Corrected

- **The mosaic expectation was wrong by a factor of (2-f).** An earlier note reasoned that a mosaic
  fraction f displaces the expected dosage to 0.5 + 0.5f. That is the mean of PER-CELL frequencies
  and assumes a monosomic cell contributes as much DNA as a disomic one. An array reads POOLED
  dosage and the loss removes template from the denominator too, so the correct expectation is
  1/(2-f). The old form overstates displacement up to 1.95x at low f, puts the 0.65 band edge at
  f = 0.30 where it truly sits at 0.46, and understates an inverted fraction about 1.8x. In logit
  space the shift is additive, -log(1-f), which is the form to build a likelihood in. Loss WITH
  reduplication is the exception, where total DNA is unchanged and (1+f)/2 is right, so the two
  mechanisms differ about twofold in implied f and must be separated by log2R before any fraction
  is quoted.

- **Mosaicism does not live in the unresolved readings**, which was the premise the guard added in
  4.10.0 rests on. Measured across stages, the unresolved zone's share of band redistribution never
  reaches a majority at any f. Below f = 0.3 in a blastomere it is 0.02 to 0.20 while the extremes
  take 0.57 to 0.92: at low f the shift is a change in the allele-dropout BALANCE, readings tipping
  from one extreme to the other, not intermediate readings accumulating.

- **Two biopsies of one embryo agreeing is not evidence.** Under the null with no event, this
  statistic correlates 0.73 between replicate arrays and the unresolved fraction correlates 0.91,
  because the drift is a reproducible property of the DNA and its amplification. The agreement
  between the two chr16 biopsies was quoted here as an internal consistency check; it is not one.

- **The blastomere in the 4.10.0 note is a paternal pronucleus** (GSM4774681), not a blastomere,
  and its chr1 is genuinely absent rather than noise: only 10.0% of its chr1 markers reach
  log2R > -0.5 against 77.2% elsewhere. A region can be a complete loss and still produce a
  middle-band-rich profile that a band-occupancy likelihood reads as "both copies present". The 3x
  guard refused it for the wrong reason, keying on a proxy for call rate.

The review is committed at `audit/MOSAIC-AUDIT.txt` with its own limits section.

---

## 4.10.3 "Diplotene"

### Reverted

- **4.10.2's concurrent Ensembl prefetch made resolution slower and is withdrawn.** Measured on the
  deployed build: 69s against a 35s baseline, and one request reached 125s and was cut off by the
  tunnel with a 524. The reasoning behind it was wrong in a way worth writing down. Ensembl's
  throttle is not a simple per-second cap that concurrency can hide under: concurrent requests from
  one IP queue behind each other or are penalised, so three at once cost more than three in a row.
  Two of the three were speculative as well, since the gene lookup is only a fallback and the
  assembly call is not always reached, so it was adding load to a throttled endpoint to warm a
  cache that often went unread.

  The helper is kept in place, disabled, carrying its numbers, because the next person to read that
  page will have the same idea and should not have to re-run the experiment.

  The measurement that stands: NCBI answers in 0.09s, gnomAD moves 928 KB in 0.14s, Ensembl takes
  12.8s for a 2 KB record and 18.9s for a 426-byte one, and Ensembl's own /info/ping answers in
  0.30s from the same container. The route is fine and the payloads are tiny. The cost is per-IP
  throttling at Ensembl and it is not something this end can optimise around by rearranging calls.

---

## 4.10.2 "Diplotene"

### Fixed

- **Variant resolution was spending its whole wall time waiting on Ensembl, one call at a time.**
  Measured from the deployed container: NCBI E-utilities answers in 0.09s and gnomAD moves 928 KB
  in 0.14s, while rest.ensembl.org takes 12.8s for a 2 KB variation record and 18.9s for a
  426-byte gene lookup. Its own `/info/ping` answers in 0.30s from that same container, so the
  route is fine and the payloads are tiny: this is per-IP throttling, which charges for round
  trips rather than bytes. Those Ensembl calls depend only on the record just resolved and not on
  each other, so they were costing the sum of their latencies for no reason.

  `prefetch_for` now issues them together as soon as the record exists, so the sequential code
  below finds the cache warm and pays the slowest instead of the sum.

  **Nothing about any result changes.** The prefetch computes nothing, decides nothing and returns
  nothing. It issues the same GETs the existing call sites are about to issue, and every one of
  those keeps its own retry count, its own error handling and its own wording, including the
  distinction between "Ensembl says no such variant" and "Ensembl did not answer". A prefetch
  failure is discarded and the real call then behaves exactly as it did before, same outcome and
  same error text. The only observable difference is latency.

  Not the diagnosis that was tried first. An NCBI API key was added on the theory that the 3
  requests per second cap was the constraint; NCBI was never the bottleneck and the key changes
  nothing here, though it is harmless and remains set.

---

## 4.10.1 "Diplotene"

### Added

- **The dosage channel is wired into the app, not only the command line.** 4.10.0 shipped the
  module and exposed it through `om origin --channel`, which left the browser tool unable to answer
  the events the channel exists for. Syngamy now retains allele dosage for every marker rather than
  every called marker, since a no-call still has an intensity reading and on a collapsed chromosome
  those are the only readings left, and calls whole-chromosome origin from it against a background
  taken from the rest of the genome.

- **The defect callout gains `channel` and `excluded dosage` chips.** The two channels answer the
  same question from different measurements and are not interchangeable, so a reader deciding how
  much to trust a call needs to know which one produced it.

- **Dosage fills what genotypes cannot reach and never overrides them.** Where the genotype channel
  answered, its answer stands; dosage only speaks where a whole chromosome left it with no
  evidence. `defects.check.ts` asserts both directions, and that a refused dosage call names nobody.

---

## 4.10.0 "Diplotene"

Parental origin from allele dosage, for the events genotypes structurally cannot answer.

### Added

- **`dosageOrigin.ts` calls which parent's copy is missing from the B-allele frequency** rather
  than from genotype calls. 4.8.5 recorded why this is needed: a whole-chromosome loss is DETECTED
  by the collapse of its genotype call rate and ASSIGNED from genotypes, so the signal that finds
  the event destroys the evidence that would name its parent. All four segmental losses in
  GSE148488 scored and all three whole-chromosome losses refused, which is not a coincidence of
  three. Dosage is read whether or not a genotype is emitted, so it survives the collapse.

  The Mendelian fact is the same one the genotype channel rests on. At a marker where the loaded
  parent is homozygous, oriented so its own allele is low: a HIGH dosage requires that parent's
  copy to be absent, since a present copy would contribute its allele and hold the dosage at or
  below the middle. A MIDDLE dosage requires both copies, since one copy cannot be heterozygous.
  Each is impossible under the alternatives, which is where the power is.

  **It is not a share and it is not centred.** Two earlier attempts formed a paternal allele share
  and centred it on the sample's own median, and for a one-parent informative set the median IS the
  ceiling, so upward headroom was zero by construction and six one-directional results followed.
  What is counted here is the rate of an event impossible under the alternative, in each direction
  separately, by the same code with the orientation flipped. The check asserts the two loss
  directions are equally reachable from the same weight of evidence, which is the property those
  six results violated.

- **`om origin --channel auto|genotype|dosage|both`**, defaulting to dosage for whole chromosomes
  and genotypes elsewhere, and naming which channel produced each answer.

### Fixed

- **The dosage channel returned both-copies-present at posterior 1.0000 from unusable intensity**,
  which is the same failure as 4.8.5 in the other channel and was caught before shipping only by
  comparing a chromosome against the rest of its own genome. A middle-band reading is impossible
  with one copy, so it is enormously informative for two copies being present, and a region whose
  intensity is noise scatters readings into that band. On a real blastomere:

      chr1              call 30.7%   own  7.7%   middle 49.7%   excluded 5.4%   between 37.2%
      its other chroms  call 83.7%   own 88.1%   middle  4.0%   excluded 0.9%   between  7.1%

  The tell is the BETWEEN band, which no hypothesis predicts: a genome resolves its alleles into
  clusters and a region that does not is not being measured. The region is now compared against the
  sample's own genome, because amplification sets that baseline per array.

  **The guard cannot separate unmeasurable from mosaic.** In a single cell an unresolved dosage is
  noise, since one genome cannot be partly anything. In a multi-cell biopsy a mosaic loss genuinely
  puts the average between the bands. Both are refused, which is right for a blastomere and
  conservative for a biopsy. Two trophectoderm biopsies of one embryo showed 29.9% and 23.7%
  unresolved over a chr16 loss and were refused here while the genotype channel called them.

### Not yet demonstrated

The channel is validated on synthetic regions and guarded against the failure that reached
posterior 1.000, but every real region tested so far has been refused by the unresolved-dosage
guard. It has not yet produced a positive call on laboratory material, and should not be quoted as
though it has.

---

## 4.9.0 "Bivalent"

Everything the browser tool does, available on the command line, with the thresholds exposed.

### Added

- **`cli/om.ts`.** The web tool is a mindless drag and drop by design: it infers what it can,
  refuses what it cannot, and offers no knob whose correct setting the user would have to know.
  This is the opposite end. Eight commands covering stage inference, linkage, one-parent origin,
  whole-directory cohort runs, the haploid census, parental reconstruction, positional enrichment,
  and a listing of every tunable with its provenance. `--json` on any of them.

  It imports the modules the app runs, so it cannot drift from the app or skip a guard it applies.

- **The thresholds are now parameters rather than constants read from module scope.** `inferStage`
  takes an optional boundary set, `callOneParentOrigin` takes an optional threshold set, and both
  default to the measured values, so the web tool passing nothing is exactly the shipped
  configuration. `cli/om.check.ts` asserts that equivalence on every branch of both functions,
  which is what stops the two front ends from quietly disagreeing.

  It also asserts that each dial actually moves the answer. The first version of the CLI
  advertised `--max-region-het` and silently ignored it, because that threshold lived inside the
  module and was not a parameter: the flag parsed, the run completed, and the number came out
  unchanged. A dial that does nothing is worse than no dial.

  Loosening a threshold can admit a region the defaults refuse. It cannot manufacture a call from
  nothing, and that is asserted too: with no informative markers the refusal survives every dial
  being opened at once.

- **Any command run off its defaults says so**, in the text output and as a `tuning` key in the
  JSON, naming each constant and both values. A number produced under a changed constant is not
  comparable with one produced under the defaults, and the validation figures do not apply to it.

---

## 4.8.5 "Holliday"

### Fixed

- **The one-parent caller issued confident verdicts on genotypes that were not measuring the
  region.** Found on a real array rather than by inspection: a blastomere with a call-rate collapse
  on chromosome 1, 69% no-call, read 59.7% heterozygous across the markers where the father is
  homozygous. Where the loaded parent is homozygous a biparental sample can only be heterozygous
  when the OTHER parent transmitted the allele this one lacks, which is the panel's allele
  frequency, so 59.7% is not a genome. The likelihood took the heterozygote count at face value,
  because a heterozygote is near-impossible under either deletion hypothesis, and returned
  both-copies-present at **posterior 1.000** from an input carrying no information.

  `MAX_REGION_HET = 0.40` now refuses the region instead. The ceiling clears `q` plus drop-in at
  the top of the measured range, so a real region is never refused, and the check asserts that
  directly. This is the region-level form of the array-level ceiling added in 4.8.4, and it exists
  for the same reason: an impossible rate is evidence about the reaction, not about the genome.
  A confident verdict from meaningless input is worse than a refusal, because only the refusal is
  visible to the reader.

---

## 4.8.4 "Holliday"

### Fixed

- **An array reading 52% heterozygous was classified as bulk genomic DNA.** The stage ladder was
  closed at the bottom, where a failed amplification drives heterozygosity toward zero, and open at
  the top, where `BOUNDS.find` matched bulk for any value at or above 0.158. So a mixture of two
  individuals, or a contaminated reaction, was called bulk DNA and handed a dropout of 0.008, the
  most confident parameter set in the module, and every downstream likelihood then treated it as
  pristine material. This is the same failure the review found at the haploid end, at the other end
  and with worse consequences.

  Found on real arrays rather than by inspection: four in one experiment read 0.53-0.56 heterozygous
  at 0.45-0.59 call rate, which is over twice what any single genome can reach.

  The ceiling is `MAX_DIPLOID_HET = 0.25`. A diploid on this panel reads 0.168, the European anchor
  is already the highest of the ancestries measured against it, and drop-in at the top of the range
  measured here adds (1 - 0.168) x 0.0525, giving 0.212 as the most a real diploid can show. Above
  the ceiling the array is reported as failed, and the caveat says what it cannot distinguish:
  mixture, contamination and outright failure look alike there, and only differ from a genome.

---

## 4.8.3 "Holliday"

Everything the detail box could say, and a run that starts from either parent.

### Fixed

- **A run with only an oocyte donor did nothing at all.** The sample loop was gated on a sperm
  donor being present, so an oocyte-only run left the Run button disabled behind a notice saying
  one file must be labelled sperm. Nothing below that gate was paternal except the label it was
  given: the tally, the segments and the Mendelian exclusion all ask whether the LOADED parent's
  copy is present, and which parent that is only decides what the answer is called. Either parent
  alone now runs.

- **A single-parent verdict was read as paternal whichever parent was loaded.** `defectsFrom` has
  always taken the loaded parent as an argument and has always used it correctly; the caller passed
  the literal `paternal`, because the result did not carry its own role and there was nothing
  truthful to pass. On an oocyte-only run that would have printed PATERNAL over a maternal loss,
  which is the one failure this display's own documentation calls the most damaging thing it can
  do. `ParentageResult` now carries `role` and every label in the result view reads it: the
  aneuploidy callout, the segment section heading and its explanation, and the kind label that
  previously read "paternal alleles absent" as a fixed string.

- **The exclusive-marker count was computed and discarded.** It is the one number on a
  single-parent run a reader can check without trusting the model: a parent who is AA has no B to
  give, and dropout removes alleles without inventing one, so markers carrying an allele the loaded
  parent lacks are Mendelian exclusion rather than a statistic. It now has a chip. Measured
  5,130-5,172 when the loaded parent's copy was removed on a real trio and exactly 0 when the other
  parent's was.

- **The `false-het` chip could never render.** `phi` was declared on the defect and never populated
  by anything, so the chip was unreachable. Removed rather than left as a field that looks
  available.

### Added

- **Material and dropout are chips on the defect box**, so the stage the run inferred travels with
  the finding rather than living only in the log and the report. A failed array contributes no
  dropout chip rather than one reading NaN, which a reader would take for a measurement that came
  out strange.

- **`defects.check.ts`**, which the module's own header has promised since it was split out of the
  component for exactly this reason. It pins the direction of a single-parent verdict in both
  roles, that two-parent dosage is never overwritten by a single-parent call, that a run with no
  parent still reports its coordinates, and that a failed stage arrives as absent rather than NaN.

---

## 4.8.2 "Holliday"

### Fixed

- **A heavily dropped-out diploid can no longer be called haploid.** The call-rate gate added in
  4.8.1 catches an array that failed outright, but not a blastomere that amplified acceptably and
  simply lost most of its heterozygotes: at 0.40 dropout it sits at 0.101 heterozygous, below the
  0.105 boundary, and was called one genome. The BAF band now vetoes that. Intensity sees a
  heterozygous locus whether or not the genotype caller kept it, and intermediate signal at a large
  fraction of markers cannot come from a single template, so a band at the diploid line decides
  against haploid however few heterozygous calls survived.

  **The veto is one-directional, deliberately.** The review asked for band position as a required
  second signal in both directions; this project has the measurement that it does not work that
  way here. Across 120 non-haploid arrays spanning blastomeres, trophectoderm biopsies and ESC
  lines, 91 fall inside the haploid band range and 22 blastomeres sit below the diploid-exclusion
  floor, because amplification widens the heterozygote band and fills the homozygote band until the
  two are one distribution. A low band is therefore treated as no evidence of one genome, and the
  recommendation was adopted only in the direction it survives contact with this material.

---

## 4.8.1 "Holliday"

The stage mechanism asserted in 4.8.0 was put to an external review. Three of its six claims did not
survive, and one of them was a failure mode rather than a wording problem: an array whose
amplification failed outright was being classified as haploid and handed a confident parameter set.
This release retracts what was wrong, fixes what was broken, and states the remaining uncertainty in
the output rather than in a footnote.

### Fixed

- **A failed amplification can no longer be classified as haploid.** Near-total dropout drives
  heterozygosity toward zero, which is exactly where a polar body sits, so on heterozygosity alone
  the two are indistinguishable. 4.8.0 read one as the other and attached a confident dropout of
  0.02 to what was noise. Call rate is now gated first, at `QC_CALL_FLOOR = 0.40`, and an array
  below it is reported as `failed` with no dropout figure at all.

- **A `failed` array no longer voids the one-parent caller silently.** Its dropout is not a number,
  and passing that into the likelihood would have produced `NaN` posteriors that compare false
  against every threshold, so a failed array would have refused every call for an invisible reason.
  It now takes the most conservative dropout measured on any stage.

### Changed

- **The ladder is documented as calibration constants, not as a copy-number model.** If each
  template failed independently with probability f, calibrating f on the single-cell rung predicts
  0.0000055 dropout for a 15-template biopsy against a measured 0.050, four orders of magnitude out,
  and exactly 0 for bulk against a measured 0.013. The bulk and biopsy rungs are therefore not
  template loss: 0.013 sits inside this platform's own replicate error, measured at 3.31% marker
  disagreement between technical replicates of one bulk sample, 98.4% of it heterozygote to
  homozygote. The numbers are unchanged and still measured; the mechanism claimed for them is
  withdrawn.

- **The chromatin explanation for blastomere against single ES cell is withdrawn.** A published
  experiment varying reaction conditions across more than 3,000 single-cell amplifications found
  amplicon size, DNA degradation, freeze-thaw and cell number all affected dropout while cell type
  had little or no effect (Piyamongkol 2003), which is the mechanism 4.8.0 asserted and the one that
  study looked for and did not find. A same-laboratory comparison also reports first polar bodies,
  carrying one genome, with the lowest dropout of three cell types (Rechitsky 1998), so template
  count does not order these data either. The difference is unresolved and is now stated as such.

- **First polar bodies are no longer treated as homozygous.** A PB1 carries a dyad, and distal to
  each crossover its two sister chromatids carry different haplotypes: roughly 44% of a PB1 genome
  is genuinely heterozygous, expected heterozygosity about 0.074 with no dropout whatsoever. 4.8.0
  called every heterozygous site in a haploid sample an error, which mis-parameterises exactly this
  material. A sample in that band is now labelled `1 homologue, 2 chromatids` and says so.

- **Every heterozygosity-shortfall estimate now carries its confounds in its own output.**
  Consanguinity, copy-neutral LOH, uniparental disomy and ancestry differing from the anchor all
  depress heterozygosity and are absorbed into the figure: first-cousin parents bias it +0.050, 20%
  of the genome in LOH +0.160, and an East Asian sample against this European-derived anchor +0.155,
  which is enough to classify bulk DNA as a single cell with nothing in the output indicating a
  problem. `StageCall` gains `basis` and `caveat`, both printed in the report.

### Added

- **`dropoutFromReplicates` estimates dropout with no population anchor.** Among markers called
  heterozygous in at least one of two amplifications of one genome, the discordant fraction is
  2d/(1+d), so d = phi/(2-phi). It is immune to consanguinity, LOH, UPD and ancestry, because those
  change which markers are heterozygous and not the chance a heterozygous one survives twice. Two
  limits are documented with it: correlated failure biases it low in proportion, to 0.90x at 10%
  shared failures and 0.50x at 50%, so it is a lower bound; and it has a floor set by the platform's
  own genotyping error, returning about 0.10 on bulk replicates where true dropout is near zero.

- The review and the model behind it are committed at `audit/STAGE-AUDIT.txt` and
  `audit/stage_audit_core.py`, so the arithmetic above can be re-run rather than taken on trust.

---

## 4.8.0 "Holliday"

Stage inferred from physics rather than fitted, carried in every output, and regions written the
way a genome browser writes them.

### Added

- **`stage.ts` infers the developmental stage of each array, grounded in template copy number.**
  Allele dropout is a sampling failure during amplification and its rate is set by how many
  template molecules the reaction started from. A locus in bulk DNA is present in millions of
  copies, so losing every copy of one allele is impossible; the same locus in a single cell is
  present in exactly two, and a heterozygote survives only if both amplify. Nothing about the
  genome changed, the reaction lost it.

      bulk genomic DNA      >10^6 cells    ~2x10^6 templates/locus    dropout 0.013   floor 100
      trophectoderm         5-10 cells     10-20 templates            0.050           floor 100
      single ES cell, WGA   1 cell         2 templates                0.199           floor 200
      cleavage blastomere   1 cell         2 templates                0.308           floor 200

  The last two both start from two molecules yet differ nearly two-fold, which is chromatin rather
  than copy number: a blastomere is in a rapid cell cycle with a decondensed, replication-active
  genome that primes less reliably.

  **RETRACTED in 4.8.1.** The dropout figures are measured and stand; the copy-number mechanism and
  the chromatin explanation do not. See 4.8.1.

- **Haploid material is separated before any diploid stage.** A polar body, pronucleus or sperm
  carries one genome and is homozygous by construction, so on heterozygosity alone it resembles a
  catastrophically dropped-out diploid; reading it as one would attach a nonsensical dropout to it.
  The boundary sits in the empty gap between measured populations, 0.002-0.10 for haploid products
  against 0.116 for the lowest diploid.

  **PARTLY RETRACTED in 4.8.1.** A FIRST polar body is not homozygous, and a failed amplification
  reached this branch and was called haploid. Both corrected in 4.8.1.

- **Stage is bundled into the result** and printed in the report, with its template count, the
  dropout it implies and the marker floor it sets. A failed array is called unknown and given the
  most conservative dropout rather than the most flattering.

- **The marker floor now varies by stage**, which the previous release recorded as outstanding:
  `MIN_INFORMATIVE_TROPHECTODERM` and `MIN_INFORMATIVE_BULK` existed and were never used.

### Changed

- **Regions are written `chr6:39,302-294,904`**, the conventional genome-browser form, in the
  callout chips and in the report's segment table. Thousands separators included because these are
  read by people, and an unseparated nine-digit coordinate is where a reader loses an order of
  magnitude.

- Documentation gains two sections, on parental origin from one parent and on material,
  amplification and stage, both stating the mechanism rather than the setting.

---

## 4.7.1 "Pachytene"

Stage handled without asking, and breakpoint intervals readable in the unit they live in.

### Fixed

- **Refined intervals print in kilobases below a megabase.** A refined edge is a median 151
  markers, which is tens of kilobases, and rendering that as `+/-0.06 Mb` threw away the precision
  the refinement bought. Now `+/-60 kb`, and where the two edges agree it prints one figure rather
  than the same number twice separated by a slash, which read as more precision than was claimed.

- **Allele dropout is inferred from the sample instead of proxied by its no-call rate.** Those are
  different quantities: a no-call is a marker with no genotype, dropout is a heterozygote read as a
  homozygote. Heterozygosity is the correct readout, since dropout removes one allele of a
  heterozygote and depresses the observed rate below the bulk figure for the platform.

  This is what makes stage handling automatic. A trophectoderm biopsy and a blastomere differ
  six-fold in dropout, 0.050 against 0.308, and the difference is visible in the arrays, so the
  user is never asked which they dropped. The check pins the inference across four stages, that it
  is monotone, that degenerate input falls back to the worst stage rather than the most optimistic,
  and that it stays bounded.

### Known

`MIN_INFORMATIVE_TROPHECTODERM` and `MIN_INFORMATIVE_BULK` are exported and unused; the marker
floor is still the strictest value at every stage. Dropout now adapts, the floor does not.

---

## 4.7.0 "Pachytene"

**Parental origin from one parent.** A run with only the sperm donor now names which parent's copy
is missing, per region.

### The evidence

Validated on a real HapMap CEPH trio, chr21, NA12891 x NA12892 -> NA12878, with the trio confirmed
from the genotypes rather than a pedigree file. The child is euploid, so each loss is constructed
from the family's own alleles, which is what makes the truth exact: both parents are genotyped, so
the surviving allele under each loss is known rather than assumed. **The mother is then hidden and
only the father is given to the caller.**

    scenario                  verdict              posterior   exclusive markers
    ADO 0.050  keep mother    known-parent-lost        1.000              5172
    ADO 0.050  keep father    other-parent-lost        1.000                 0
    ADO 0.050  keep both      both-present             1.000               147
    ADO 0.450  keep mother    known-parent-lost        1.000              5160
    ADO 0.450  keep father    other-parent-lost        1.000                 0
    ADO 0.450  keep both      both-present             1.000              1163

    12 of 12 correct, across dropout from 0.050 to 0.450, with only the father loaded.

The `exclusive` column is the mechanism. It counts markers where the loaded parent is homozygous
and the sample carries the allele that parent does not have. When his copy is gone the count is
5,130 to 5,172; when the OTHER parent's copy is gone it is exactly 0, at every dropout rate.

### Why it works where two attempts failed

**A parent who is AA has no B to give.** So a child reading BB at that marker cannot be carrying
his copy, and dropout cannot manufacture the observation, because dropout removes an allele and
never invents one. The estimator is that fact counted properly, with the two ways it can be faked
priced in.

There is no allele share and no centre, which is what killed the earlier attempts: a share centred
on its own median could only move one way, because for every q below 0.5 the median IS the ceiling.
Each marker here contributes a likelihood and the region evidence is their product, so neither
direction is privileged. The check pins symmetry, that the caller is not a constant function, and
that raising dropout to 0.6 cannot turn an intact region into a loss.

### Added

- `oneParentOrigin.ts`, and its verdict feeds the chromosomal-change callout, so a sperm-only run
  shows a named parent with `basis: one parent`, the informative marker count and the posterior.
- The display mapping is pinned in both parent roles: with the sperm donor loaded
  `known-parent-lost` reads paternal, with a maternal reference it reads maternal, and two-parent
  dosage still outranks the single-parent call where both exist.

### Still true

The sibling-referenced caller from 4.6.0 remains inert and unused. This release does not repair it;
it makes it unnecessary for the single-parent case.

---

## 4.6.1 "Synapsis"

A red callout that leads with WHOSE copy, and an honest report that 4.6.0's sibling caller does not
work.

### Added

- **Chromosomal change now has its own callout, and origin is the headline.** A reader wants whose
  copy first and the coordinates second, so the parent is the large text and the evidence sits
  beneath it as chips: interval, span, refined edge, event type, informative markers, posterior,
  false-heterozygous fraction, and which channel the call rests on. Where origin cannot be
  determined the headline says so in the same place and the same size, because a defect with no
  origin is a result rather than a missing field.

- **`defects.ts` carries the mapping, separately from the markup.** Which origin is shown, and what
  is said when there is none, is the layer where a mistake does real damage, so it is testable on
  its own. The check pins the case that matters: a region with no annotation must not borrow a
  neighbour's origin, must not claim a basis, and must keep its coordinates, since position needs
  no parental genotype.

### Fixed

- **A homologue call is no longer displayed as a named parent.** The uniparental path writes
  "this parent, other homologue (meiotic)", which names no parent; a prefix match would have shown
  it as one.

### Does not work

- **The sibling-referenced caller shipped in 4.6.0 is inert.** Measured against the intensity
  channel on real arrays, 8 embryos with 3 or more arrayed cells: of 74 regions the intensity
  channel called a copy-loss, the sibling caller reported no deletion for all 74, and reported no
  deletion for all 74 matched negatives as well. Sensitivity 0.000, specificity 1.000. It is a
  constant function and carries no information. Two candidate causes are open, an allele-dropout
  parameter fed the array's no-call rate instead, and sibling-established heterozygous markers that
  are not reliably heterozygous, and they have not been distinguished. Stated here rather than
  quietly patched.

---

## 4.6.0 "Synapsis"

Parental origin without a parental array. The blocker was never the method.

### Added

- **`siblingOrigin.ts` calls which parental copy is missing using the embryo's OWN event-free
  cells as the reference.** Multiple blastomeres of one embryo share both parents exactly, so where
  an event is present in some cells and absent in others, the unaffected siblings establish that
  embryo's heterozygous sites and the affected cell is read against them. No parental array is
  involved at any point.

  Three hypotheses per region: the reference-side copy lost, the other copy lost, and no deletion.
  Posterior over the three, called only above 0.95, refused otherwise.

- **There is no centre in it, deliberately.** Each marker contributes a likelihood and the region's
  evidence is their product. The quantity that was previously mis-centred does not exist here.

- **Drop-in measured on this platform**, on 113 same-genome array pairs across four experiments:
  medians 0.0390 to 0.0525. False heterozygosity is the only one-directional term in this
  construction, because a falsely heterozygous marker always retains the reference allele, so it is
  reported with every call and bounds it. Above 0.15 the call is refused rather than downgraded.

- **The sibling agreement rule scales with the panel**, `>=ceil(k/2)` with a floor of 2. Under a
  fixed rule of two the false-het fraction RISES with panel size, 0.0019 at two siblings to 0.0675
  at ten, because more cells give a spurious call more chances to appear. Adding arrays would have
  made the bias worse while appearing to make the call more reliable. Under the majority rule it
  falls to 0.0000 at ten.

### Fixed

- **`recentre` is documented as unusable for a one-parent marker set, and why.** The one-parent
  asymmetry recorded in earlier audits was not a geometric limit. The two loss hypotheses are
  symmetric about the disomic expectation; the asymmetry came from centring on a median that IS the
  boundary of that distribution's support, since for every q below 0.5 the majority of markers sit
  exactly at 1.0. A statistic whose centre is a boundary point can only move one way, which is what
  produced this project's one-directional results. Dropout was never the cause: expected share is
  invariant in it, 0.8754 to 0.8748 across ADO 0.000 to 0.600.

- **`paternalShareOneParent` is marked do-not-use.** It selects on the scored cell's own
  heterozygosity, which selects the markers where both alleles survived and forces the share to 0.5
  at every one of them. No information survives the conditioning. Retained rather than deleted so
  that callers fail loudly.

### Not done

The module is not wired into Syngamy, which streams one sample at a time and would need a two-pass
restructure to hold a sibling panel. The audit harness is the route that reads many arrays at once.

---

## 4.5.0 "Chiasma"

4.0 said where a chromosomal change is and which parent an extra copy came from. This release
answers the question actually asked of the data: when a piece of a genome is MISSING, whose was it,
and where does it sit.

### Added

- **Losses now name which parent's copy is missing.** Only gains carried a parental annotation
  before, and in this material a loss is the more common event. The direction INVERTS between the
  two: an extra paternal copy raises the paternal allele share and a missing one removes it, so the
  same deviation names opposite parents depending on the event. Passing a loss through the gain
  call named the wrong parent every time. `callLossOrigin` writes that down once, and a test now
  fails if the two calls ever agree.

- **Each chromosome reports one parental contribution or two**, from heterozygosity at markers
  where a parent is homozygous. Reported, never used as a gate: the boundaries hold for a bulk
  reference parent, and against a single-cell reference they do not separate at all.

- **Called regions are scored for WHERE they sit**, against common fragile sites, genes over
  500 kb, centromeres and telomeres, and ENCODE Repli-seq late-replication valleys. Needs only
  breakpoint position and no parental genotype, so it runs on every sample.

  The null is the substance of it. A region can only be called where the array carries markers and
  where amplification produced calls, and marker density tracks gene density, so intervals drawn
  uniformly along the genome report an enrichment for almost any feature. Each null interval is
  drawn on the same chromosome carrying the same number of informative markers as the observed
  region. The self-check constructs a genome where marker density and feature position are
  deliberately confounded and requires no enrichment to be reported, because that failure is
  invisible in the output: it produces a confident p from a biased comparison.

- **The report plots the null**, not just a p value. Observed against the middle 95% of the null
  per feature, and for anything clearing 0.05 the null distribution itself with the observation
  marked. On this material a ratio of 0.97 reached p 0.002 because the null was narrow, and a
  reader shown only the p would take that for a finding. A methods section is included that can be
  pasted rather than paraphrased.

### Fixed

- **The one-parent ploidy boundaries were wrong, and quietly so.** They were 0.12 and 0.30, scaled
  from the two-parent figures and marked provisional in the source rather than measured. Against a
  bulk parent across nine of his children, his one-genome products read 0.0428 to 0.0587 and his
  biparental children 0.0894 upward. BOTH CLASSES SAT UNDER THE OLD 0.12 BOUND, so every biparental
  child of a single genotyped parent was called uniparental. Now 0.065 and 0.085.

  Anyone who ran a single-parent reference on a 4.0.x build should re-run: the affected calls read
  "one parent's genome and nothing else" for cells that had two.

### Changed

- **The informative-marker floor halves, 400 to 200**, which halves the smallest region that can
  carry a parental label from a median 47.6 Mb to 23.8 Mb. Re-measured on the material it governs
  by the criterion it was originally set by, and both halves were required, since specificity alone
  cannot justify a floor: a caller that never calls has no error.

      specificity   0.0030 two-way error at 200 markers, against a 1% bar, on five biparental
                    real WGA arrays scored against a bulk-genotyped father, where every window is
                    euploid so any directional call is an error
      sensitivity   1.00 at 200 markers with half the alleles dropping, 0.87 at 50, on a real CEPH
                    trio with a trisomy built from the parents' own alleles and dropout drawn from
                    a Markov model fitted to real arrays (a marker beside a dropped one is 2.18x
                    more likely to drop)

### Validation

- **Gain detection on 330 karyotype-confirmed trisomy 21 cells**, 100 re-karyotyped blind by an
  independent laboratory. AUC 0.9779 against 114 confirmed euploid on one platform, 0.9991 against
  27 on the other. Inside the trisomy cells the affected chromosome ranks first of 22 autosomes on
  both, which is what shows the statistic reads dosage and not amplification.
- **Parental origin on a real CEPH trio**, its structure established from the genotypes rather than
  from a pedigree file: 0.0002 opposite homozygotes for parent-child against 0.0677 for the couple.
  Nine of nine calls correct, including refusing the euploid child.
- **1,189 copy-number regions across four experiments**, 318 gains and 871 losses, each scored
  against every feature set.
- **Per-allele calling bias measured against meiosis**, where segregation fixes the true ratio at
  exactly one half: 0.4934 over 127,826 markers, a bias of -0.0066.

### Not claimed

No parental split is admissible from this material. Six one-directional results were produced while
building this release and every one was an artefact; each is recorded in the audit with the
measurement that exposed it, including that a self-derived null turns platform failure into a
parental result, that a single-parent reference cannot see a paternal loss at all, and that
reconstructing a father from the cells being scored is circular.

The parental origin of a biologically real trisomy remains unvalidated: no public series carries
one with both parents genotyped. Detection is validated on 330 real cases; naming the parent of one
is not.

---

## 4.0.1

A positive control for the parental direction of a gain. 4.0.0 shipped that annotation shown only
to refuse correctly and to be orientation-safe, never to fire correctly, because no confirmed
trisomy with a known parent exists in either body of data. That gap is now partly closed.

### Added

- **The direction is recovered on a real array, both ways.** Constructed by the same method an
  external review used and labelled the same way: re-weighting real per-marker allele intensities
  rather than simulating genotypes. The array is the sperm donor's own zygote, 12,185 informative
  markers, one parent genotyped and the mother absent, which is the one-parent path the tool marks
  provisional.

      euploid, untouched    deviation +0.0000  ->  unclear     correct
      paternal gain 2:1     deviation +0.1502  ->  paternal    correct
      maternal gain 1:2     deviation -0.1714  ->  maternal    correct

  Runs as `node --experimental-strip-types audit/gain-positive-control.ts`.

- **And it proves the recentring on real data rather than by argument.** That array's own centre
  is 0.5993, not the theoretical 0.5000. Judged against the theoretical centre its UNTOUCHED
  euploid genome sits +0.0993 away, over the margin, and would be reported as a paternal gain that
  is not there. Recentring on the sample's own median puts it at exactly +0.0000. Both numbers are
  now pinned in the module's check, so the one design decision that would invert calls cannot be
  undone silently.

Still not checked against a real trisomy, and the tool still says so everywhere it makes the call.

## 4.0.0 "Crossover"

Chromosomal change is specific now: every breakpoint is measured rather than read off a scanning
window, and every extra copy is asked where it came from.

Breakpoints are measured, not read off a scanning window. This is the precondition for
everything the lab wants to do next with them.

### Added

- **Marker-resolution breakpoint refinement.** The scan slides windows of 2,400 markers and up,
  stepping a quarter of a window, and used to report the winning window's EDGES as the event's
  coordinates. Those edges are an artefact of the scan: the step alone quantises them to a median
  1.96 Mb. Measured against 848 events spliced from real arrays at known marker positions, the
  window edge lands a median 373.5 markers from the truth, 2.51 Mb, with a p95 of 15.46 Mb.

  Once the window has found roughly where the event is, coordinate ascent on the two edges at
  marker resolution finds where it actually starts and stops, evaluating the same
  segment-versus-background likelihood ratio the scan already uses. Off a cumulative sum each
  candidate is O(1). Median error 9 markers, 0.063 Mb: **40x better**, and it localises more
  events than the window edge did, 629 of 848 against 412. Circular binary segmentation and binary
  segmentation were measured on the same constructs; both beat the window edge, both lost to this,
  and both cost more.

- **A confidence interval on every breakpoint, and no bare coordinates.** A point estimate would
  claim 151x more precision than the calibration supports. The interval is a profile
  likelihood-ratio drop, and the drop is 12 rather than the nominal chi-squared 1.92: the nominal
  value gives 75.5% coverage of a nominal 95%, and deflating by the measured variance-inflation
  factor only reaches 92.0%. Twelve was calibrated empirically to 96.9% on 470 edges. Every
  bootstrap fails outright, 64.9% to 66.7%, and widening the block length does not rescue them,
  because the failure is model misspecification rather than the resampling scheme.

  On the lab's own chr4 event the intervals come out at 0.23 Mb and 0.40 Mb.

### Fixed

- **One event was being reported as two.** Peak-picking ran on window coordinates, so an event
  straddling a window boundary survived as two adjacent non-overlapping hits. On a real array a
  chr4 loss was reported as 0.06-5.84 Mb plus 5.85-38.29 Mb; refinement resolves both to one
  region running to 42 Mb, which is also what the laboratory record says, singular. The pick is
  now re-run on the sharpened coordinates, so an event that straddles a boundary is counted once.
  This would have double-counted events in exactly the enrichment analysis the precision was built
  for.

### Added: where an extra copy came from

- **Every gain is annotated with its origin, or with the reason there is none.** A gain is two
  different questions and the tool answers whichever one the material supports.

  In a cell carrying ONE parent's genome, which is what a pronucleus is, the parent is not in
  question. What is asked instead is whether the extra copy is that parent's OTHER homologue,
  which is meiotic, or a duplicate of the same one, which is mitotic. A haploid genome cannot be
  heterozygous, so the other homologue shows as heterozygosity far above the array's error rate:
  AUC 1.000 from 50 markers, holding even on degraded arrays at 24-45% spurious heterozygosity. A
  duplicate of the SAME homologue is bit-identical to a single copy: AUC 0.037 to 0.119, below
  chance, and intensity does not rescue it. So the meiotic case is reported and the identical case
  is reported as indistinguishable from a normal single copy, which is NOT the same as reporting
  no gain, and the wording says so.

  In a cell carrying both parents, allele dosage names the parent, at 400 informative markers by
  default. Two details are load-bearing and both are checked. Each marker is converted to a
  paternal allele share before anything is averaged, because in raw B-allele frequency an extra
  paternal copy pushes the value down at father-AA markers and up at father-BB markers, so a
  statistic that mixes them cancels the signal to exactly nothing. And the centre is the sample's
  own median, never the theoretical 0.5, which against a reconstructed parent is biased toward
  maternal by up to 0.077, or 46-72% of the whole band separation, and inverts calls.

  Medians, never means. On a single amplified cell the mean band separation compresses 36% to
  0.106, which would put the bands inside the decision margin and refuse every single-cell gain;
  the medians do not compress. That is asserted as a test rather than left as a comment.

### Fixed while building it

- **The gain background was self-referential.** The first implementation judged a chromosome's
  heterozygosity against the genome-wide rate, so where several chromosomes carry an event the bar
  is lifted by the very events under test and each one reads as ordinary. This is the same failure
  the segment scan already solves with an external per-chromosome null, made again. The background
  is now the median of the OTHER chromosomes, floored, with thin chromosomes excluded.

### Refused, and stated in the code

- **Intensity plays no part.** log2R localises 6% of these events against 74% for the absence
  channel, because half of what this tool looks for is an origin switch, a haplotype substitution
  with no dosage step to find at all. Combining the two is worse than absence alone, 53% against
  74%. The same conclusion the copy-number scan reached independently.
- **Refinement below a 0.70 call rate is not reported.** Positional error rises from 6-9 markers
  to 160-222 and coverage falls to 0.81.
- **An event that does not localise says so** rather than being given a wide interval.
  Unconditional coverage never exceeded 0.743 for any method measured; the shortfall is a
  localisation failure and no interval width fixes it.
- **Naming a parent for a gain needs a biparental cell and both parents loaded.** On a
  uniparental product the question does not arise and the tool says so rather than reporting a
  parent it inferred from nothing.
- **Focal gains cannot be annotated on a single amplified cell.** 400 informative markers span a
  median 47.6 Mb; a 5 Mb region carries a median 51. Whole chromosomes and large arms can be
  annotated, focal events cannot, and the refusal names the marker count.
- **NOT VALIDATED ON A TRUE POSITIVE.** The bands and marker requirements are measured, but on
  constructed contrasts: no confirmed gain with an independently known parent of origin exists in
  the material behind them, where every dosage event is a loss. This has been shown to refuse
  correctly and has never been shown to fire correctly. The one array in the set carrying
  whole-chromosome gains is a 42% call-rate QC failure and the tool refuses all five of its gains,
  which is the right answer and not a positive control. Stated on the page, in the report and in
  the documentation, not only here.
- **The 2,400-marker floor is unchanged.** An external review measured a full-contrast floor of
  800, which conflicts with this project's standing 2,428-marker figure, and the conflict cannot
  be resolved from the record. Changing it would trade a measured floor for an unreconciled one.

## 3.7.2

Hotfix. An external review of the chromosomal-change work found two things wrong with what this
repo says about itself. Neither changes a call the tool makes today; both would have if the
statements had been trusted one step further.

### Fixed

- **A real event was sitting in the clean-genome calibration set.** `segments.check.ts` described
  the five paternal pronuclei the segment null was calibrated on as "present on every autosome".
  They are not: GSM4774681 (pMII-2) carries a chr1 event, measured at absence 0.190 against 0.0105
  genome-wide with log2R -1.79 and a multiscale LRT of 5,411 against a leave-one-out threshold of
  134. Checked here against the subset this repo ships, scoring that product against a reference
  built from the other four: chr1 is its worst chromosome at 9.78% absence against 1.04%
  genome-wide, 9.4x, which corroborates the direction on 460 informative markers. The same check
  flags pMII-6 chr20 at 12x, so the number of affected arrays is not established either.

  A real event counted as noise makes every false-positive rate derived from that set optimistic
  by an unknown margin. The thresholds are NOT changed here: swapping them without recalibrating
  on a verified-clean set trades a known bias for an unknown one. The claim is corrected, the
  consequence is stated, and the recalibration is an open item.

- **The ploidy gate is documented for what it actually covers.** The BAF-band gate separates a
  haploid meiotic product from a bulk diploid adult, which is the only distinction it was measured
  on, and it does that correctly on all 46 arrays this tool has run. It does not generalise: on 120
  post-zygotic arrays, 91 fall inside the haploid band range and 22 blastomeres sit below the
  diploid-exclusion floor, because whole-genome amplification widens the heterozygote band and
  fills the homozygote band until they are one distribution. It was never wrong here; it would
  have become wrong the moment it was reused as a general ploidy branch.

### Added

- **`obligateHet.ts`**, the statistic that does separate there, against the day the tool handles
  cleavage-stage material. Heterozygous fraction at markers where a parent is homozygous: haploid
  products measure 0.047 to 0.100 and the lowest post-zygotic array 0.4213, a gap four times the
  width of the entire haploid range, where the BAF band has no gap at all. Boundaries at 0.20 and
  0.45 sit inside that gap rather than on either edge, which means the lowest observed diploid is
  refused rather than called, and that is the intended trade.

  It carries its own limits. With one parent instead of two the diploid signal is diluted to the
  markers where the unseen parent differs, so the boundaries narrow, every call is marked
  provisional, and the same fraction can be a call with two parents and a refusal with one. Under
  200 informative markers it refuses outright. Nothing is wired to it yet.

## 3.7.1

Progenitor's tagline read "builds a parent\u2019s SNP array" on the live page, with the escape
printed rather than the apostrophe. A JSX attribute in quotes is raw text and does not process
string escapes, so the six characters went straight to the screen. Found while screenshotting the
page for a presentation, which is the only reason anyone looked at that line closely.

## 3.7.0

The two halves of a run are joined up. Progenitor builds the array; Syngamy is what you point it
at; between them sat a file picker.

### Added

- **Download and open in Syngamy**, under the inferred-array export. Saves the file, opens
  Syngamy in a new tab, and hands the array across already loaded as the donor. The saved copy is
  still written, because a run you cannot reproduce from disk is not a run; the handoff is the
  convenience on top of it.

  The transfer is a `postMessage` between the two tabs rather than storage. `localStorage` caps
  out around 5MB and this file is three times that, `sessionStorage` is per tab and its
  inheritance by an opened tab is not something to rely on, and IndexedDB would work but leaves a
  copy of an identifiable person's reconstructed genotype sitting in the browser after the run,
  which is the one thing this tool has always refused to do. The new tab asks and the opener
  answers, in that order, because the opener cannot know when the receiver's listener is mounted
  and a message sent before then is simply lost.

  The new tab says what it is waiting for while it waits, and says what to do if nothing arrives.

### Fixed

- **A blocked popup no longer destroys the run it came from.** The obvious way to write this is
  `window.open(url, '_blank')`, and it has a bad failure mode that showed up on the first test:
  the target differs from the current page only by its hash, so a browser blocking the popup can
  fall back to navigating THAT tab to it instead. Measured: `window.open` returned null, meaning
  the code believed it had failed, and the current tab navigated to Syngamy anyway, taking the
  finished Progenitor run and the array it was handing over with it.

  The blank tab is now opened first, with no URL at all, and only pointed at Syngamy once it is
  known to exist. An empty target has nothing to fall back to, so a blocked open stays blocked:
  the page stays put, the file is still saved, and the button says so.

## 3.6.1

One button at a time. 3.6.0 put all three steps in the toolbar at once, where two of them sat
disabled until their turn came: three buttons together, two of them dead, is a worse instruction
than one.

Run is the only control in the toolbar now. The build button appears at the foot of the
reconstruction section once there is a group worth building, one per buildable group since with
more than one that is the actual decision, and it names the parent it is about to build. The
call button appears below the reference once there is an array to score against, with a line
saying it re-reads every file and why. Each disappears when its step is done.

## 3.6.0

Progenitor is a SNP array builder that also calls parental origin, rather than the other way
round, and it no longer freezes the tab while it works.

### Fixed

- **The page stopped locking up.** Two passes had no yield in them and nothing to show, so a run
  spent tens of seconds looking hung. Both now report progress and hand the main thread back
  between chunks.

  The arithmetic underneath was also doing several times the work it needed to. Choosing the
  agreement depth called the full build once per candidate depth, and a build materialises a
  600,000-entry map, so picking a threshold cost more than the build it was picking for; there is
  now a path returning the one scalar that choice needs. Contamination was a SECOND full pass
  over every probe, with a map lookup each, to recover a number the first pass already had. And
  every build rebuilt an 825,657-element array of probe pairs to iterate over. Measured on an
  18-array experiment, the blocking arithmetic went from roughly ten seconds to 2.4, and what
  remains is chunked.

  What is NOT optimised away is the re-read: scoring an array needs its B-allele frequencies, and
  the allele codes held in memory carry only genotypes, so every file streams a second time. That
  is most of the wall clock. It now sits behind its own button and reports per file.

### Changed

- **Three buttons instead of one.** Compare the products, build the reference, then call the
  arrays. Each states what it is about to do, and the expensive step is no longer reached by
  surprise: a run that turns out to have three parent groups can stop after the first.

- **Every step logs what it is doing.** Each pair with its rate, shared marker count and verdict;
  each group with its members, its Y-bearing products and whether it clears the floor; each rung
  of the depth ladder with the ascertainment it retains; each array as it is scored, and for
  reference members the depth and marker count of the leave-one-out build behind it.

- **The groups are stated, not only plotted.** Membership led with a concordance matrix, which is
  the evidence and is unreadable past about eight products. A roster above it now names the
  members of each group, which one is the father and on what basis, which group the reference
  came from, and which groups sit under the floor. The report carries the same, with Y-bearing
  counts per group.

- **Framed as what it is.** The page builds a parent's SNP array out of the cells that parent
  produced. Parental origin is what scoring arrays against that array gives you, not the other
  way round.

### Added

- **A methods paragraph written to be pasted into a paper.** Prose rather than documentation:
  what was gated and at what threshold, how products were grouped and why chained links are not
  accepted, how chromosome Y named the parent and when that naming is withheld, how the depth was
  chosen, what contamination means, and how a call is read. Bracketed where the text depends on
  the samples rather than on the software.

- **The Egli Lab and Columbia University Irving Medical Center are named** in the README and at
  the foot of every documentation page. The Columbia Stem Cell Initiative mark was already on the
  home page; nothing said it in words.

## 3.5.3

3.5.0 documented a safety property it had not implemented, and left the exports and pages
describing the feature as it was before it. Both fixed.

### Fixed

- **The inferred-array mark now does something.** The export banner said, and the docs repeated,
  that a reconstructed genotype carries a machine-readable mark "which the tools that ingest
  arrays look for". Nothing looked for it. Two things do now.

  Syngamy reads the mark on drop, before anything is profiled. It warns in the run log and states
  it in the first section of the report, naming the file and saying that absence measured against
  a reconstruction is not the same quantity as absence against a measured array. Nothing is
  blocked: using one as the donor is the intended workflow and the reason the file is written.

  Progenitor REFUSES such a file as a product. This is the direction that could have been quietly
  wrong. The ploidy gate does reject the file as currently written, but only because the export
  carries no BAF column and an undefined band lands on "borderline" - a column happening to be
  absent, not a reconstruction being recognised. On its genotypes the file is homozygous at every
  marker, which is 0% heterozygosity and the cleanest haploid product ever submitted. Add a BAF
  column or convert the file and it would pass, fold a reference into itself, agree with itself
  everywhere, and report a contamination of nothing.

### Changed

- **The run manifest records how the parent was named.** A `naming` block: which side was
  reconstructed, which group is the father's, whether naming was withheld and for which of the
  three reasons, and the Y verdict for every product. The per-sample rows already carried the
  origin call; there was no record of what the naming rested on.

- **Documentation caught up with 3.5.0.** The Progenitor exports table lists six artefacts rather
  than five and describes the new one. The samples CSV and run manifest entries say they carry the
  origin call. "Using the page" describes what Run actually does now. Syngamy's docs gained a
  section on using a reconstructed donor: how it is detected, what is stated where, and the two
  ways its numbers are not comparable with a measured array's. The README and the home page card
  described Progenitor as reconstructing a genotype, which was its 3.4 scope.

## 3.5.2

Laboratory sample identifiers had reached three source files in 3.5.0 and 3.5.1, in comments
citing the two arrays that set the chromosome Y thresholds and the one that turned out not to be
a parent array. The measurements are what those comments are for and they are unchanged; the
identifiers are gone. Nothing else changed.

## 3.5.1

The validation record for 3.5.0, written down. No behaviour changed.

### Validation

3.5.0 shipped with one experiment's answer checked against a held-back sperm array. Three further
tests were run afterwards, on everything else available, and are now on the docs page.

- **A second held-back parent array, on an independent public dataset.** Reconstructed from five
  pronuclei, the reference reads the father's own bulk array at 0.03x and a SECOND array of the
  same man at 0.09x, both present, and his own zygote at 0.75x. Two unrelated oocyte donors read
  3.06x, absent, and 2.44x, refused. No wrong answer.

- **Specificity, on 26 arrays of other people.** One experiment's paternal reference against every
  array of the two unrelated experiments: not one reads as carrying him. Nineteen decisive at
  1.75x to 10.84x, seven uncalled, five of those already excluded by the gates.

- **Robustness to depth.** Rebuilding from five of the eight products moves no call: 17 of 18
  identical, one degrading to a refusal, no inversion. Contamination rises from 0.48% to 2.24%
  across that drop, which is where the refusal comes from.

- **The Y as an independent witness.** Y-bearing is read from chromosome Y call rate and intensity
  and shares nothing with the autosomal agreement that produces the call. Across the twelve
  Y-bearing arrays in three experiments the two never disagree: seven agree, five are refused by
  the autosomal side, all of them low-quality arrays.

### Measured and deliberately not acted on

- **The BestProbeset column stays ignored.** These arrays flag 7,180 of 825,656 probes as not
  best, and 95% of those sit at a locus that already has a best probeset, so honouring the flag
  looked like an obvious improvement. It is not one. Filtering to best probes only moves no call
  of 18, shifts every ratio by at most 0.07x, raises contamination from 0.483% to 0.492% rather
  than lowering it, and costs 5,795 markers. The flag is parsed and ignored on purpose, and this
  is the measurement saying so.

- **The one array with no row in the laboratory record is not a parent array.** It sits in the
  same folder as an experiment that has no sperm array, and would have anchored it. It does not:
  every one of that experiment's 17 arrays reads no parental contribution against it at 2.17x to
  10.22x. It is an unrelated female diploid.

## 3.5.0

Progenitor calls parental origin. Drop an experiment's arrays in and it works out which came from
which parent, with no array of either parent anywhere in the run.

### Added

- **Parental origin for every array, from products alone.** The reconstruction was already the
  hard part and it stopped one step short: it established which arrays share a parent and how
  contaminated a reference built from them would be, then left the question the whole thing exists
  for unanswered. Now every array that goes in is scored against the reconstructed genotype, not
  only the ones it was built from. An array carrying that parent's genome came from that parent;
  one decisively lacking it came from the other.

  Which parent gets reconstructed is not chosen and does not need to be. The largest group that
  clears the five-product floor is built, whichever parent it belongs to, because the other side
  reads off it. An experiment whose paternal products are too few but whose maternal ones are not
  is answered by reconstructing the mother.

  Arrays the gates excluded are called too, and their rows say so. A product has to be a clean
  haploid cell to go INTO a reference; that is not a reason to withhold an answer ABOUT it, and a
  fused or half-failed zygote is the case this gets asked about most.

- **Naming the reconstructed parent, from chromosome Y.** Grouping products by shared parent says
  nothing about which parent: a set of siblings looks the same either way. A maternal cell cannot
  carry a Y, so a group with one Y-bearing product is the father's.

  TWO measurements are required and neither is sufficient. Each rule alone inverts a real
  experiment. Call rate alone is wrong on an array that genotypes 86.2% of its Y probes while its
  Y intensity sits a full log2 below its own autosomes, exactly where arrays with no Y sit: an
  absent chromosome still produces calls, and taking that at its word names a maternal group
  paternal and inverts every call in it. Intensity alone is wrong in the other direction, on an
  array reading -0.10 log2 while calling not one Y probe.

  Requiring both separates cleanly over 46 arrays: Y-bearing arrays call 93.7% to 97.3% of their Y
  probes at +0.16 to +0.43 log2, and every other array either calls 0.0% or sits at -0.81 to
  -1.25. Naming is withheld rather than guessed where no group carries a Y, where more than one
  does, and where the file's Y panel is too thin to ask. The split between the two sides holds in
  every case; only the naming is withheld.

- **The reconstructed genotype is exported as an array file.** It was withheld for several
  releases on the reasoning that a file of homozygous calls belonging to an identifiable person
  could be re-imported as though it had been measured. That reasoning was right about the hazard
  and wrong about the gain: reconstructing a parent exists in order to have something to call
  parental origin against, and an experiment with no array of that parent has nothing else to use.
  It is written now, and the hazard is handled where it lives. Every such file opens with a banner
  saying it is not a measured array and carries a machine-readable mark saying the same. The
  banner lines are comments the header sniffer skips, so it loads anywhere a real export does, and
  it can be dropped into Syngamy as the donor to call arrays the Progenitor run never saw.

### Changed

- **Run does the whole run.** It used to compare pairs and stop, then wait again for a choice of
  which group to reconstruct, which is not a decision the person holding the folder can make
  better than the arrays can. One press now groups, names, reconstructs and calls. Sample chips
  appear as each file finishes rather than when the whole batch does, they carry the origin call
  once there is one, and the answer sits at the top of the page rather than four sections down.

### Validation

Three experiments, 46 arrays, no array of either parent used anywhere in the pipeline. On the one
experiment where the answer is known independently, the sperm donor's own array was held back and
used only to mark the result: 16 of 18 arrays were callable and all 16 agree with it, including
two the laboratory record has the wrong way round. The other two are below the call-rate floor and
read unclear against both. No inversions.

The remaining two experiments have no such array. Both are internally consistent: the group named
paternal by its Y is the group the arrays labelled paternal fall in, and in one of them the
buildable group is the mother's, which the run reports as such and reads the paternal side off.

## 3.4.2

The segmental gain direction is no longer unvalidated, and the report no longer tells the reader
that gains are not called.

### Added

- **A positive class for segmental GAINS, constructed from real arrays.** 3.4.0 shipped the gain
  direction with its status stated plainly: no segmental gain occurs anywhere in the 46 arrays, so
  it had never been shown to fire on a true positive. It now has been, by the method already used
  to set the marker floor rather than by simulation. A block of one real array is spliced into
  another of the same series, on the same platform from the same lab, so it carries genuine
  trisomic intensity and genuine amplification artefact. The block is taken from a chromosome this
  tool already calls as a whole-chromosome gain (0.34x to 0.41x call rate at +1.60 to +1.95 log2)
  and spliced into the same chromosome of an array carrying no event.

      spliced markers   span      no-call       score    verdict
           2,400         6.9 Mb   83% vs 12%    3,162    copy-gain
           4,800        14.0 Mb   82% vs 12%    6,256    copy-gain
           9,600        30.6 Mb   82% vs 12%   12,312    copy-gain
          19,200        62.8 Mb   82% vs 12%   25,007    copy-gain

  Four of four, the right direction every time, at 12.6x to 100x the threshold. The IDENTICAL
  construction with the block taken from a euploid array returns nothing at all four sizes, so
  what is being measured is the gain and not the act of splicing, and the recipient unspliced
  returns nothing on every autosome. Localisation was checked at three offsets: the reported
  interval lands within 3.5 Mb of the spliced one on a 40 to 47 Mb block. Below the marker floor
  the scan refuses rather than answering weakly, at 600 and at 1,200 markers.

  What this does NOT establish is a false-positive rate for gains in the wild, because a spliced
  boundary is sharper than most biological ones. That is stated wherever the result is quoted.

- **`origin.py` gained the copy-number scan.** `scan_copy_number` and `SEGMENT_LRR_SHIFT` mirror
  the TypeScript, and `Segment` carries the `kind` that distinguishes the two channels. The
  self-check pins both languages to the same fixture: 9,600 markers, 8,009 no-calls, 47,995,000 bp
  span, score 14,298.1603, kind copy-loss, identical to sixteen significant figures. A divergence
  here would be invisible to a user, who only ever runs one of the two.

### Fixed

- **The report told the reader that gains are not called on this platform.** That sentence was
  correct when it was written and became false in 3.4.0, when the copy-number channel shipped. The
  per-sample callout now names which of the two findings it is looking at, and says what each one
  means, since `parental alleles absent` and `copy loss` are different claims about a region and
  were being introduced by the same paragraph.

- **The Methods section described only the allelic channel.** A reader reconstructing the analysis
  from the report would not have learnt that a second scan runs on a different indicator, nor that
  the intensity requirement is what separates a real deletion from a region the array merely failed
  on. Both are now stated, along with the constant and its measurement in the constants table.

## 3.4.1

The report's first section now says whether a chromosome is missing, rather than leaving it to a
per-sample block further down.

### Changed

- **Chromosomal change is flagged in the Result summary.** A reader who takes only the first
  section of the PDF was previously told the origin call and the absence rates, and had to reach
  the per-sample detail to learn that a chromosome was gone. Whole-chromosome and segmental
  changes now appear twice in that first section: as a banded line above the table naming the
  affected samples and what was found in each, and as a per-sample column in the table itself.
  The coordinates and the evidence stay where they were, below.

  Samples with nothing found read "none detected" rather than being left blank, and the note
  under the table says what that does and does not mean: a region under the size floor, or one on
  a chromosome whose calls are not measuring it, is not reported either way.

- The summary table was over-wide and ran past the right margin when an oocyte donor was loaded.
  Every column is narrower; the table now ends exactly at the margin in both layouts.

## 3.4.0

Segmental copy number, on a second indicator. The aneuploidy work is now complete on both scales
and in both directions.

### Added

- **Segmental gains and losses, with their coordinates and length.** A second segment scan, on a
  different indicator from the first and deliberately so. The existing one asks where a PARENT's
  alleles are missing, which a region can do with its DNA entirely present, if it came from the
  other parent alone. This asks where the DNA itself is not there, read from the array no longer
  calling: a region that is gone stops producing genotypes. Both events are real, both are
  reported, and they are labelled rather than merged.

  A region is called only where the no-call rate rises AND the intensity agrees by at least 1.0
  log2. That requirement is load-bearing rather than cautious. Across 46 arrays the scan without
  it returns 31 regions, one of them a 44.9 Mb loss on the sperm donor's own bulk DNA, who is a
  diploid adult and cannot have one. With it, 17 survive and NONE fall on any of the six bulk
  diploid arrays, whose truth is zero. Rejected regions carry shifts of -0.72 to +0.38; kept ones
  carry -0.97 to -1.94.

  The largest call is a contiguous 38 Mb loss at the start of one chromosome, reaching 85% no-call
  against a 14% background at -1.94 log2.

  EXTRA COPIES use the same machinery with the shift reversed, and the status of that direction is
  stated plainly rather than implied: no segmental gain occurs anywhere in those 46 arrays, so it
  has never been shown to fire on a true positive, only never to fire on a true negative.
  Whole-chromosome gains ARE validated, where the effect is an order of magnitude larger.

## 3.3.0

Whole-chromosome aneuploidy, called from the call rate. And a correction: three chromosomes this
tool was suppressing as artefact are real losses.

### Added

- **Aneuploidy calling.** A chromosome that has been lost yields no DNA, so the array reads
  nothing at its probes and cannot genotype it. The call rate on that chromosome collapses while
  every allelic statistic there merely looks noisy. That is a genotype-level signal needing no
  intensity, which matters because fine copy-number work on amplified material is refused
  elsewhere in this tool for measured reasons.

  Over 1,012 chromosome observations from 46 arrays across three experiments, eleven chromosomes
  sat at 0.20x to 0.41x their genome's median call rate and the other 1,001 at 0.78x to 1.16x. The
  gap is empty by a factor of 1.9. All eleven also carried an intensity shift of 1.59 to 2.04 log2
  where the rest never left -0.79 to +0.42, so two channels that fail in unrelated ways agree on
  every event.

  The sign of that shift says which way it went, used only after the call rate has established
  that something is wrong: six sat at -1.59 to -2.04 and five at +1.60 to +1.95. So gains ARE
  called at whole-chromosome scale, where the effect is an order of magnitude larger than the
  noise that forced the refusal of finer copy-number work.

  Whose copy went is read from what survives: if the declared parent's alleles are still present
  on the remainder, the copy that went was the other parent's; if absent, it was theirs. Where
  nothing remains, or the sample carries none of that parent's genome anywhere, no parent is
  attached and the report says so.

### Fixed

- **Three real chromosome losses were being suppressed as array mis-clustering.** The 3.0.4 gate
  withheld any chromosome whose allelic ratio was depressed, on the reasoning that a paternal
  pronucleus is its father's on every autosome by construction. That reasoning was wrong: the
  series this tool is validated against exists precisely because chromosomes ARE lost in these
  embryos, so a pronucleus can genuinely lack one.

  Both a loss and a mis-clustering depress the allelic ratio. The CALL RATE separates them,
  because a mis-clustered chromosome still calls and an absent one cannot. The gate now applies
  only where the chromosome is still being genotyped, and a collapsed one is reported as the loss
  it is. Every chromosome the gate had suppressed carries a call rate at 0.28x to 0.36x and an
  intensity shift near -1.8, which is a lost chromosome and not a clustering failure.

## 3.2.1

Documentation catching up with 3.2.0. The mosaic contrast shipped after the report's Methods and
constants table were written, so it reached the PDF only as a limit and was not described or
sourced there.

### Added

- A Methods paragraph for the mosaic scan in the report, and the mosaic contrast in the report's
  constants table with the measurement behind it: four bulk diploid arrays with no mosaic ran
  -1.65 to 5.18 over 88 chromosome observations, against 9.7 to 22.7 at a third of cells.

### Changed

- The README said Syngamy reports which parental genome is present and on which chromosomes. It
  now reports which REGIONS within them too, and has since 3.1.0.

## 3.2.0

The last two items in the review. A mosaic is now detected from the allelic ratio, and the
segment scan exists in both implementations rather than one.

### Added

- **Mosaic detection, from the continuous allelic ratio rather than the genotype.** A genotype
  call is a threshold on that ratio, the AB-against-BB boundary sits near 0.917, and a chromosome
  present in a fraction f of cells has a true ratio of 1/(2-f), which does not cross that boundary
  until f = 0.909. Below it the call never changes, so no number of markers helps. The ratio
  itself is not thresholded and moves from the first percent.

  Each chromosome's mean deviation from 0.5 over its heterozygous sites is contrasted against the
  sample's OWN other chromosomes. That being internal is load-bearing: a degraded array has a
  globally shifted ratio, so an absolute deviation measure calls the whole genome mosaic, and a
  contrast cannot, because the shift sits in both terms.

  Measured on four bulk diploid arrays with no mosaic anywhere, 88 chromosome observations, the
  contrast runs -1.65 to 5.18. Titrated by shifting the observed ratio on one chromosome and
  keeping each array's own noise:

      fraction of cells   contrast, four arrays        detected
          0.15             3.4   1.1   7.9   4.0        none
          0.20             5.2   2.1  12.2   6.3        1 of 4
          0.30             9.7   4.6  22.7  12.5        3 of 4
          0.50            21.3  11.7  48.5  28.6        all 4

  So the floor is half the cells for a reliable call and a third for a partial one, not the 0.20 an
  earlier estimate suggested. Still far below the genotype route's structural floor of nine tenths.

  NO FRACTION IS REPORTED, only that a mixture is present: inverting the statistic is biased low by
  roughly half, because the truncation that hides a mosaic from the genotype also removes the
  most-shifted markers from the mean. The axis is withheld entirely on a uniparental genome, which
  is homozygous by construction and has no heterozygous sites for a mixture to shift.

- **The segment scan now exists in Python as well**, pinned to the TypeScript by a self-check on
  the same planted event: 9,600 markers, 48.0 Mb, rate 0.108, score 2244 in both. The two
  implementations had diverged when the scan shipped in the browser only, and a divergence there is
  invisible to a user, who runs one of them and never the other.

## 3.1.3

The last three items from the independent review. One deletion and two refusals that had no
measurement written beside them.

### Removed

- **The insertion state is gone from the copy-number state space**, along with the hypothesis only
  it could reach. It was there so the mechanism ledger could report an insertion as not tested
  rather than excluded, but the ledger never read the state space: that entry is built from
  whether reads and a construct sequence were supplied. So the state bought nothing, and it cost
  something real. An insertion adds no probeset and alters no genotype at any of the 825,656 fixed
  markers, so its emissions were bit-identical to the normal diploid state in both channels, and
  two states with identical emissions split the posterior between them. The normal state was being
  reported at roughly half its true probability wherever the prior allowed the other.

  Nine states, not ten. A test now pins that nothing in the space shares emissions with the normal
  state, which is the property that made this worth removing rather than merely tidy.

### Changed

- **Mosaic detection from genotypes is refused, with the reason.** No such feature existed, but
  nothing said why one should not be added. A genotype call is a threshold on the allelic ratio,
  and the AB-against-BB boundary sits at 0.9168 over 115,199 heterozygous and 372,656 homozygous
  calls. A chromosome at mosaic fraction f has a true ratio of 1/(2-f), which does not cross that
  boundary until f = 0.909, so below that the genotype does not change and the marker carries no
  information at all. Power is 0.000 to 0.250 at every fraction from 0.02 to 0.80. Adding markers
  buys nothing against a structural floor.

  Recorded alongside it: the continuous B-allele frequency reaches f = 0.20 as a per-chromosome
  contrast, four to six times better, using data already in the files. That is a different
  evidence channel with its own validation burden and is noted as available rather than claimed.

- **Meiosis I versus meiosis II is refused near the centromere, with the marker counts.** Within
  2 Mb of the centromere the informative markers, father heterozygous and mother homozygous, number
  ZERO on chromosome 9, one on chromosome 16 and three on chromosome 15, against a median of 30.
  At 1 Mb the power is zero on twelve of twenty-two autosomes. That is an absence of data rather
  than a noise problem and no threshold recovers it; the regions are probe-poor because they are
  repeat-rich, which is a property of the sequence rather than the vendor.

  Also on record: at a real trisomy an AAB locus is called AB by the array, so a meiosis II
  isodisomy is not reliably called homozygous by the genotype channel to begin with.

## 3.1.2

The other half of the copy-number work was built, measured, and refused. This entry is the
measurement, so nobody attempts it again from the same evidence.

### Changed

- **Attributing an extra chromosome to a parent is refused, on measurement rather than caution.**
  The question was: a sample has two copies of a chromosome where it should have one, so whose is
  the extra? Three mechanisms produce that, and they were constructed from real arrays rather than
  simulated: a haploid paternal pronucleus as background, with one chromosome merged either with a
  second pronucleus of the SAME father, which is two independent draws of his genome, or with a
  maternal pronucleus, which is a different person.

  The first mechanism, his one copy replicated, is invisible by construction: duplicating a
  chromosome changes no genotype at any marker.

  The other two were measured, self-normalised against the sample's own other chromosomes:

      construction                    het at father-het    het at father-hom
      nothing added, 88 observations  median 1.01 max 2.03  median 1.01 max 1.91
      both his homologues             1.46 to 7.25          1.19 to 2.17
      a chromosome from someone else  3.74 to 4.99          2.18 to 3.17

  The null runs through the signal. A chromosome carrying nothing extra reaches 2.03 while the
  weakest real second-genome construction reaches 1.46. And the two mechanisms that are both real
  events separate by 2.17 against 2.18, half a percent. Neither call is available.

  An earlier pass looked far more promising because the control range was taken from four
  observations of one chromosome, 0.79 to 1.27. Measured properly across 88 chromosome
  observations the null maximum is 2.03. The gap was in the sample size, not in the genome.

- **Recorded alongside it**, from the independent review: even with clean separation the statistic
  detects material that is NOT PATERNAL, never material that is MATERNAL. The true mother and an
  unrelated adult are indistinguishable on it at matched array quality. Dispermy, a second sperm,
  lab contamination, a swapped sample and a genuinely retained maternal chromosome all produce the
  same signature, so the attribution would need the oocyte donor's own array to be a measurement
  rather than an inference from an assumed pedigree.

- The refusal appears in the segments documentation with its table, and in scope and limits.

## 3.1.1

Chromosomal change is stated where a reader will see it, rather than left to be inferred from a
table of chromosomes.

### Changed

- **A segmental loss now announces itself.** It was reported only in a table below the fold, which
  is the wrong place for it: the whole-chromosome verdict for a partly lost chromosome is
  "unclear", and that reads as an absence of information rather than a located event, so a reader
  skimming the result has no reason to look further.

  It now appears three times over. As a filled badge in the sample headline beside the sperm type
  and the zygosity, so it is visible without opening anything. As a bordered callout at the top of
  the detail, naming the chromosomes, the total span and each region's coordinates and rate. And
  as a banded block in the report PDF directly under the sample's class, before the axis numbers.

  Each says the same three things: the paternal genome is missing over these regions, the
  whole-chromosome verdict cannot show it, and it is a LOSS rather than a statement about physical
  copy number or a gain, since gains are not called on this platform in either channel.

### Added

- **A Methods paragraph for the scan** in the report, and the segment floor, the segment threshold
  and the allelic-ratio floor in the report's constants table, each with the measurement behind it.

- **A documentation section** covering why the null is external, why the threshold is empirical
  rather than a closed-form tail, the titration that set the 2,400-marker floor, and why a marker
  count is not a resolution on an array whose spacing runs from 1 bp to 21 kb.

## 3.1.0

Syngamy now reports WHERE along a chromosome the paternal genome is missing, not only whether a
whole chromosome is. First half of the copy-number work; the half an independent review found
recoverable.

### Added

- **Segmental loss detection.** A chromosome that is partly lost reads as neither present nor
  absent: the rate lands between the two references, the chromosome comes back unclear, and the
  part that IS missing goes unreported. A multiscale scan now finds those regions and reports each
  with its span, its marker count, its local absence rate, the rate it was scored against, and its
  score.

  Three things had to be right, each measured rather than chosen.

  The null is EXTERNAL: the median per-chromosome rate over the sample's OTHER chromosomes. Scoring
  a window against the sample's own genome-wide rate misses the largest events, because a big event
  inflates the rate it is being tested against.

  The threshold is EMPIRICAL, never a closed-form tail. A Bonferroni exact binomial and an
  Erdos-Renyi run-length scan both fail here for the same reason the run-length p-value did:
  absence artefact on amplified single cells is spatially clustered, so a tail computed under
  independence calls artefact significant. Measured on five genomes carrying no event, 110
  chromosome scans, the maximum score reached was 139. The threshold sits at 250.

  The floor is 2,400 called informative markers, titrated on real material: a block of a maternal
  pronucleus spliced into a clean paternal pronucleus of the same series, so the segment is a
  genuine alternative genome carrying real amplification artefact. Twelve constructions per size.
  At 1,200 markers the weakest scores 152 against a null maximum of 139, which is 1.09x and not a
  detection; at 2,400 the weakest scores 431, which is 3.1x. The old 200-marker floor was not a
  weak test, it was no test at all.

  False segments on the five genomes known to carry none: **0 of 110 chromosome scans.** That
  threshold is fitted on the same null it is evaluated against, so it awaits an out-of-sample
  clean cohort, and the score is printed beside every segment so a marginal call reads marginal.

  A chromosome withheld by the mis-clustering gate is not scanned either: the same broken calls
  would produce a confident segment inside it.

### Fixed

- **Localisation, caught during titration.** Merging every overlapping window across scales
  reported a 12 Mb loss as spanning 231 Mb, because a whole-chromosome window carrying a diluted
  version of the same event overlaps the tight one. Replaced with peak-picking: the highest-scoring
  window wins and everything overlapping it is discarded as another view of the same event.

### Changed

- Segments appear in the result page, in the report PDF, and in the run log. The span is stated as
  the resolution rather than the event, because a marker count is not a resolution on an array
  whose spacing runs from 1 bp to 21 kb.

## 3.0.4

A chromosome whose genotype calls are not measuring it now gets no verdict. Found while scoping
the copy-number feature, in shipped code rather than in the new work.

### Fixed

- **Array mis-clustering was being reported as a chromosome-scale paternal loss.** When a
  chromosome's probes cluster badly its genotype calls come back systematically wrong, absence
  rises to something indistinguishable from a real loss, and every genotype-derived measure agrees
  with it, because they are all reading the same broken calls. Nothing in the tool caught it.

  On the bundled series this was live: **GSM4774681 chromosome 1 was reported ABSENT at 18.36%**
  on a sample whose genome-wide verdict is "parent genome present". That sample is one of the
  sperm donor's own paternal pronuclei, so every autosome of it is his by construction. It is also
  one of the five arrays this project certifies as clean and validates against.

  Only the allelic ratio shows the cause, so a per-chromosome gate now runs before any verdict:
  at least 40% of a chromosome's B-allele frequencies must sit outside 0.15 to 0.85. Measured
  across the five paternal pronuclei, 110 chromosomes: correctly clustered ones run **0.752 to
  0.976**, and the false event ran **0.130**. The gap is empty by a factor of 5.8. With the gate
  in place exactly one chromosome of 110 is withheld and it is the one that was wrong; nothing
  real is suppressed.

  A withheld chromosome reports its measured extremeness and the reason in place of a verdict,
  and its ratio against the ceiling is blanked in both the table and the report, because a ratio
  is a finding and a chromosome that was never measured has none. A genuinely trisomic chromosome
  fails this gate too, since its allelic ratios cluster at a third and two thirds rather than at
  the extremes. That is correct rather than a cost: gains are not called on this platform in
  either channel, so the alternative to withholding is a wrong answer.

## 3.0.3

An independent adversarial review of the per-locus run-length statistic found its null false on
the material this tool targets. The p-value is withheld until a sample earns it.

### Fixed

- **The run-length p-value assumed something that is not true of amplified single cells.** The
  statistic's whole principle, stated in its own module, is that "artefacts are independent
  across markers, real structural events are contiguous", and the p-value is the Erdos-Renyi tail
  under that assumption. It is false here. Measured on the public series, maximal runs of two or
  more against the independence prediction at the same fitted rate: **6.1x to 10.3x** on haploid
  pronuclei, and **2.0x to 2.1x** even on bulk gDNA from unrelated adults. Every material class
  exceeds the model, single cells by roughly an order of magnitude, so a run of ordinary artefact
  reaches a length the model calls impossible and the tail reports significance for it.

  Independence is now measured rather than assumed, on the sample's own absence indicator across
  every chromosome except the one under test, because a real event is contiguous by definition
  and would otherwise inflate the null it is about to be judged against. A sample that does not
  demonstrate independence gets no p-value at all: the run length is still reported as the
  observation it is, and the new `independence_not_demonstrated` verdict says why, with the
  measured excess in the note. A sample that supplies no off-locus markers gets no p-value
  either, because nothing has shown the null holds.

  Both implementations are fixed together. The shared fixture pins the arithmetic of the tail,
  which is unchanged; the gate is a separate property with its own checks.

  This is the second time a statistic in this project was validated on the wrong material class.
  The first was the same assumption, checked on bulk gDNA, where it does hold.

### Changed

- **`UNIFORM_CV` may not be extended to mosaicism**, and now says so with the measurement.
  Against uniform artefact its separation ratio is 0.58 at f=0.10, 0.83 at 0.30 and 1.00 at 0.80.
  A ratio at or below 1 is no separation at any fraction. The shipped use, telling a blend of two
  whole genomes from a partial loss, is a different question and is unaffected. Also recorded:
  the 1.106 figure in its calibration is a same-individual replicate, the easiest control in the
  set, and carries no weight.

- **A number in `emissions.py` was overstated.** Its mosaic-BAF correction was described as
  "nearly three standard deviations at a typical sigma_BAF of 0.03". Heterozygous BAF sd on the
  arrays this tool reads is 0.088, so the error is 0.95 sd, not 2.8. The correction itself stands.

## 3.0.2

One line of throat-clearing removed from the References section of both documentation pages.
No code changes.

---

## 3.0.1

An accuracy audit for Progenitor, and the report sections the other two tools already had.

### Added

- **A second accuracy audit, for Progenitor**, on 19 public arrays from the same GEO series the
  Syngamy audit uses. A man's genotype is reconstructed from eight of his sons' paternal
  pronuclei and graded against his own array, which the reconstruction never sees: 35 correct,
  6 refused, 0 incorrect. He reads present at 0.12x and 0.18x on his two replicate arrays, and
  the nine genomes that are not his read absent at 4.71x to 9.11x. Membership put the five usable
  products in one group and split them from seven maternal pronuclei with zero misread pairs
  out of 45.

  The series carries his own bulk DNA beside the pronuclei, which turns three modelled quantities
  into measured ones: ascertainment against his real 16.66% heterozygosity, contamination by
  asking his array which asserted homozygotes are really heterozygous, and absence on nine
  genomes known not to be his. His two replicate arrays bound what the ground truth itself can be
  wrong about, at 4.15%.

  A density sweep rebuilds the reference on disjoint slices of the array, down to every sixteenth
  marker: 31 references, 310 classifications, all clean. One array density is one observation.

  Two results are recorded rather than resolved. Predicted contamination sits 0.6 to 1.1
  percentage points below observed at every depth, which is inside the replicate floor, so this
  data cannot separate a model that understates it from a truth that overstates it; the direction
  is stated, since if the model is the wrong one it is optimistic. And the sub-floor cases score
  the parent rather than an offspring, because this series contains no offspring of these
  products, so they do not reproduce the measurement the five-product floor was set from.

- **Methods, Files read and Citation in the Progenitor report**, and **Methods and Citation in
  the panel report**, both written from the actual run rather than from a template, and both
  naming the refusals that fired. Parity with the Syngamy report, which has carried all three
  since 2.0.

### Changed

- The ploidy gate moves out of `ProductPanel.tsx` into `productTriage.ts`. It is the gate the
  whole feature stands on, and living in a `.tsx` put it out of reach of every check that runs
  under node, including the audit, which then could not exercise the shipped rule at all.

- Syngamy's run log sits at the top of the page, where Progenitor's already was.

## 3.0.0 "Recombinase"

The enzyme class that catalyses strand invasion and exchange.

Progenitor: a parent's genotype reconstructed from the haploid cells that parent produced, with
no array of the parent. The first release in which a reference the tool built itself can carry a
verdict, which is why the release is major and why most of the work below is the refusals.

### Added

- **Progenitor.** Drop several haploid meiotic products of one parent, a pronucleus, a polar
  body, a single sperm, and the parent's genotype is reconstructed from them. Where the parent is
  homozygous every product carries that allele and the site is recovered exactly. Where the parent
  is heterozygous each product carries one at random, so the site is recognised only by observing
  both alleles, and a site where every product happened to agree enters the reference as
  homozygous. That single error is the whole method's uncertainty, it is quantified per run, and
  it is reported next to every number derived from it.

  Depth is chosen by measurement rather than by a rule of thumb. A marker enters when at least
  `m` products called it and agreed; the tool takes the deepest `m` whose retained marker set
  still holds 90% of the heterozygosity the products imply. An earlier rule fixing `m = n - 1` is
  right at five products and wrong at eight, where it demanded seven agreeing calls and left the
  marker set at roughly 62% of the genome's heterozygosity. Against a father whose real array
  gives 16.66% heterozygosity,
  the m>=2 estimate lands at 16.93% with no array of him in the calculation.

  The agreement probability is measured, not assumed. Products are not independent at a
  heterozygous site, because a probe with allele-specific dropout calls the same homozygote in
  every product. Over 133,631 heterozygous markers the excess above 2^(1-m) runs 1.24x at m=2 to
  1.42x at m=4 and vanishes by m=6.

- **Membership before reconstruction, from concordance alone.** Every pair of products is
  compared before anything is built, because a reference built from two people is worse than no
  reference. Measured across two experiments: 4.68% to 9.70% opposite homozygotes between products
  of one father over 46 pairs, and 9.88% to 16.10% between products of different fathers over 45.
  The separation is real and the margin is 0.18 of a percentage point, and one genuine
  cross-father pair sits at 9.88%, under the same-parent cut. So a per-pair label is not evidence
  of membership: grouping requires every pair inside a group to agree, as an exact maximum clique.
  Single linkage merges two men through that one edge. A greedy all-pairs pass returns the wrong
  partition outright and was replaced after a synthetic case caught it, which the real products,
  all of one father, could not have.

  This is not hypothetical. A series of fourteen products split into groups of nine and five on
  concordance alone, with no reference and no prompting: it holds two sperm donors.

- **A refusal below five products.** At three and fewer, every true offspring tested inverted to a
  decisive wrong answer rather than to a refusal, 24 of 24 across two experiments. The tool
  declines to build rather than returning a weak answer.

- **Six things a reconstructed reference will not report**, each measured coming out wrong rather
  than guessed at, and each travelling with the results into every export. Chief among them:
  androgenetic against biparental, because that axis derives from the parent's own heterozygosity
  and a reconstruction is homozygous everywhere by construction.

- **Five more public example arrays load with one click**, from the same GEO series as
  Syngamy's, GSE148488 (Zuccaro et al. 2020). Every sixteenth marker, 51,604 each, served from
  this instance and subsampled without altering any value. They are the five paternal pronuclei
  of one man, chosen because they are a positive case: they pass every gate, group as one,
  reconstruct and verify. Someone seeing the feature for the first time should see it work
  before they see it refuse. The terms of use name them alongside Syngamy's.

- **Documentation for Progenitor**, and five exports: a Letter PDF in the same format as the
  Syngamy report, pairwise concordance and the concordance matrix as CSV, per-sample results
  including every excluded file with its reason, and the whole run as JSON. Every artefact repeats
  on every page that the reference was inferred rather than measured, so a page lifted out of
  context still says so. The CSV provenance sits in a leading `#` block, which pandas and R skip,
  so the file still loads in one call.

  The reconstructed genotype itself is not exported. The build is deterministic, so the product
  files and their checksums are the reference. A file of homozygous calls belonging to an
  identifiable person, re-importable as though it were a measured array, is a hazard with no
  matching gain.

### Fixed

- **The noise ceiling ignored the reference's own false absence, so true products read unclear.**
  At a contaminated marker the reference asserts a homozygote where the parent is really
  heterozygous, and a haploid product carries the other allele half the time. That is absence
  through no fault of the sample, it equals half the contamination, and the tool already computed
  and displayed it while leaving it out of the ceiling that decides present against absent. Every
  true product was therefore judged against a ceiling roughly its own size: on the five public
  pronuclei, three of five read unclear at 1.10x to 1.32x. Admitting it puts all five at 0.12x to
  0.50x. The error direction is the dangerous one, calling a true relative unrelated.

  The term is passed only where it applies. A diploid sample needs the minor allele from its other
  parent as well, so its added absence is far smaller, and the diploid path is calibrated without
  one: an unrelated adult reads 3.37x against the reconstruction and 3.35x against the father's
  real array. Nothing about a run against a measured parental array changes.

- **Two products of one parent could be told apart by an unstripped filename.** Membership keys on
  the product name, and the name was taken with its extension in one place and without it in
  another, so a group's members failed to match the loaded files.

### Changed

- **Validated against a second father held out entirely.** 18 arrays, 16 agreeing with the answer
  his real array gives, with 0 inversions and 0 false relationships. His own replicate array is
  recognised at 0.14x by a reference built from nothing but his sons' pronuclei.

- **The browser and Python implementations are pinned to each other.** The arithmetic exists twice
  on purpose, and the storage differs because a tab cannot afford what a workstation can: Python
  holds per-marker observation lists, the browser holds one `Uint8Array` of allele codes per
  product, so eight products of 825,657 probes cost about 6.6 MB rather than tens of millions of
  boxed tuples. A cross-implementation check compares every quantity a reader is shown, at every
  depth, plus the grouping and the synthetic case that broke the greedy implementation. It reports
  0 divergences.

- The three documentation pages now share one shell, so their numbering, deep links and furniture
  cannot drift. Each links to the other two as buttons rather than as a sentence of prose.

- The landing page names all three tools rather than two.

---

## 2.4.2

README trimmed to what a repo README is for. No code changes.

---

## 2.4.1 "Kinetochore"

Two ways the noise ceiling could be trusted when it was not measuring anything.

### Fixed

- **A failed array can no longer manufacture a relationship from its own noise ceiling.** The
  ceiling is the no-call rate times the heterozygous band, so an array whose band is artefact
  gets a ceiling wide enough to swallow an unrelated genome's absence entirely and report it as
  present. Three arrays with bands of 36.9% to 43.0% produced ceilings near 19% and were called
  as carrying a parental genome they are unrelated to, at 6.4% to 11.9% absence. A band above
  30% is now refused outright: a diploid genome reads 15% to 16%, and dropout inflates that only
  into the high twenties.

  The bound is on the band rather than the call rate, and the public series settles why. Zuccaro
  A8 sits at a 53.0% call rate with a 22.2% band, and the tool reads it correctly in both
  directions: present against its own father at 0.91x, refused against two unrelated egg donors
  at 1.16x and 1.18x. A call-rate gate would have discarded that. The failing arrays sit within
  three points of A8 on call rate and fourteen points away on the band.

- **A parent with no heterozygosity no longer resolves the second-parent axis.** That axis is
  derived from the parent's own heterozygosity, so a parent showing none drives the expectation
  to 0.00% and the comparison then reads biparental for any nonzero rate, on no evidence. The
  expectation is now withheld below `PANEL_HET_FLOOR` and the call goes to unclear. A genotype
  reconstructed from haploid meiotic products is homozygous everywhere by construction, so this
  was reachable from ordinary use rather than only from a malformed file.

### Added

- The Columbia Stem Cell Initiative mark on the landing page, linking to the lab this was built
  for. Landing page only, and above the research-use notice rather than inside it.

### Changed

- The test fixtures built parents who were homozygous at every marker, which is what let the
  second defect survive: the second-parent axis was asserted on in four tests and exercised in
  none of them. Parents now carry the 17% heterozygosity a real person does, assigned
  independently of the sample's genotype so the informative denominator is not shrunk. Two
  synthetic checks also set every B-allele frequency to 0.5, a 100% band that no genome
  produces; both now use 16%.

---

## 2.4.0 "Kinetochore"

Three quantities the tool measured and did not report. None of these changes a verdict; they
change what a reader can see behind one.

### Added

- **The pseudoautosomal region is scored apart from the rest of chrX.** PAR1 and PAR2 sit on
  both the X and the Y and recombine between them, so a Y-bearing sperm delivers paternal PAR
  alleles while the rest of the X is legitimately absent. Pooled into one chrX bucket that
  positive control was invisible: the row now appears as `X:PAR` and reads present on exactly
  the samples where `X` reads expected_absent. Boundaries from the GRC assembly region reports,
  and the split is skipped entirely when the assembly could not be called, since PAR2 is at a
  different address in the two builds and guessing would put ordinary chrX markers into a
  control. On a full Axiom array roughly 976 informative markers fall in the PAR; on the
  shipped one-in-eight example subsets 122 do, below the 200-marker floor, so the row appears
  on real arrays rather than on the demonstrations.
- **The noise ceiling now reports the two factors it is the product of**, the no-call rate and
  the heterozygous fraction, next to the ratio in the browser and in the report. A ratio that
  fell between two runs can mean the signal shrank or the ceiling grew, and those have opposite
  implications. The factors say which.
- **Dispersion between chromosomes is reported**, as the coefficient of variation of the
  per-chromosome absence rate together with the cleanest chromosome. Measured on the shipped
  arrays: 1.11 for two arrays of one man, 0.76 for a degraded but true parent-offspring pair,
  0.11 for an unrelated adult, and 0.10 for a 50:50 blend of two unrelated genomes. Where a
  sample already sits in the uncalled band and the absence is uniform across every chromosome
  with none of them reaching the ceiling, the report now says the shape is a blend rather than
  a partial loss, and that no reanalysis separates the two genomes.

  Worth stating plainly, because the reverse is easy to assume: uniformity does not identify a
  mixture on its own. An unrelated adult and a blend measure 0.112 against 0.104, which is no
  separation at all. What the statistic distinguishes is a genome-wide difference from one
  confined to part of the genome, and it is consulted only where the rate is already ambiguous.

### Fixed

- `parental_origin` accepts a known assembly rather than always inferring one, for callers
  holding a manifest.

---

## 2.3.1 "Kinetochore"

The rest of the audit's list, and a test of its central proposal that falsifies it.

### Added

- **The two implementations are now pinned against each other**, which is the guarantee a
  docstring had claimed for months without a test behind it and which the audit found violated on
  every pair it tried. `tools/gen_parentage_fixture.py` writes what `origin.py` computes on the
  shipped examples; `parentage.check.ts` asserts the browser reproduces the class, zygosity,
  sperm type, verdict and six rates to within 1e-6. Verified live rather than assumed: perturbing
  one constant by 0.0001 fails it with the sample and quantity named. After the 2.2.x and 2.3.0
  fixes the two agree exactly on all four cases.
- The CLI's verdicts are role-neutral, `parent_genome_present` and `no_parental_contribution`,
  matching the browser. The old names asserted "paternal" even when the call was made against the
  oocyte donor.

### Changed

- The shipped description of example A8 said biparental. Since 2.2.1 withheld zygosity below a
  60% call rate it returns unclear, which is the correct answer on a 53% call rate and is now what
  the example and the documentation say.

### Tested and not adopted

- **The audit's answer to open question A does not survive contact with the data.** Its proposal
  was to separate structural missingness from dropout by the heterozygous fraction: dropout
  converts heterozygous calls to homozygous so the fraction falls, while markers that were never
  assayed leave it unchanged. The audit derived this by simulation and said testing it against the
  real files was the obvious next step.

  Measured. On the HapMap individuals typed on a different panel, which is the structural case,
  the heterozygous fraction of called markers is 31.08% against 19.34% for the well-called
  individuals on the same platform, a ratio of 1.61. On the amplified Zuccaro material, which is
  the dropout case, it is 31.98% against 16.43%, a ratio of 1.95. Both rise, and the dropout case
  rises further, so the statistic does not separate them and points the wrong way if used.

  The reason is already recorded in this codebase: below a 60% call rate the dominant artefact on
  amplified material is spurious heterozygous calls, not clean allele dropout. The textbook
  het-to-hom conversion is real but it is not what moves the measured fraction. Open question A
  stays open, with one proposed answer now closed off.

---

## 2.3.0 "Kinetochore"

Responds to an independent scientific audit of the 2.x inference path. Every finding below was
re-derived from the shipped data before being acted on; one recommended fix was implemented and
removed, with the reason recorded in the code.

### Fixed

- **The chrY denominator collapsed to zero against a female parent.** Counting chrY probes only
  where the *parent* was called meant an oocyte donor, who calls no chrY, gave a denominator of
  zero, so a Y-bearing sample read as not Y-bearing. Whether a sample carries a Y is a property of
  that sample, so the denominator is now chrY probes on its own array. Introduced in 2.2.1 by the
  fix that made the chrX exemption require a chrY measurement.
- **The heterozygosity gate was genome-wide when the failure is per chromosome.** A real
  blastomere carried ten of twenty-two autosomes with no homozygous B-allele population at all,
  68.9% heterozygous, while twelve intact chromosomes diluted the genome-wide figure to 32% and
  nothing fired. Evaluated per chromosome it fires immediately. Re-running the pronucleus audit,
  it also flags one broken chromosome each on three haploid pronuclei the audit had passed
  (chr1 60%, chr13 61%, chr21 55%) and does not fire on the clean bulk-DNA donor.
- **The CLI inferred a Y-bearing sperm from a missing paternal X without consulting chrY**, which
  is circular: X loss and chrX assay failure produce the same observation. It now requires the
  chrY measurement and reports `unknown` when the two do not corroborate, which is what the
  browser has always done.
- **The duplicate-parent check used 0.99 where the ingestion layer already used 0.90.** Real
  replicate arrays of one person concord at 95.8% and a parent-offspring pair at 54.9%, so the
  threshold sat above the thing it was meant to catch and missed every genuine duplicate.
- **The pseudoautosomal region was stated as 1.99% of chrX.** Measured on the shipped array it is
  6.80% (178 of 2,617 markers, GRCh37 PAR1 and PAR2), understated 3.4-fold. The sex bands rest on
  the measured chrX-to-autosome ratio, not on that figure, so no call was affected.
- **A docstring claimed the two implementations are pinned by a shared fixture.** They are not.
  The per-locus layer is; the parentage layer is not, and the audit found the two sides differing
  on every pair it tried. The claim is removed rather than restated.

### Investigated and not taken

- **A goodness-of-fit gate on the per-chromosome likelihood ratio.** The audit is right that the
  ratio ranks two hypotheses without asking whether the winner fits. Implemented, it refused a
  genuine whole-chromosome loss: at 31.5% absence the data sits 51 sigma from the 5.5% unrelated
  rate, because a missing chromosome is a third state more extreme than unrelated, not a poor fit
  to it. The ratio only reaches its threshold when the observation is at or beyond one of the two
  rates, so its decisive calls are directionally sound. The audit's concrete case, fifteen
  chromosomes called on a mis-clustered array, is caught by the per-chromosome heterozygosity gate
  above, which is where that failure actually lives. Reasoning recorded at the call site.

---

## 2.2.6 "Kinetochore"

### Fixed

- **The run log stopped following itself.** It scrolled to the bottom on each new line and then
  read its own geometry back in the `onScroll` that assignment triggers. Landing mid-append, that
  read often computed a false "the reader has scrolled up", after which every further line
  arrived below the fold in silence. The component now marks its own scrolls and ignores them,
  scrolls in a layout effect so the move lands in the same frame the line paints, and treats
  anything within 24px of the bottom as still at the bottom rather than 4px, which sub-pixel
  positions alone could exceed. Scrolling up to read still holds position, and returning to the
  bottom resumes following.

---

## 2.2.5 "Kinetochore"

### Changed

- **Served at `originmarker.app` as well as `ezrakruger.cc/originmarker`, from one container.**
  No second image and no per-host build was needed: `app/main.py` reads `ROOT_PATH` into
  `FastAPI(root_path=...)`, which strips the prefix only when the request carries it, so one
  process answers both `/api/health` and `/originmarker/api/health`; and the frontend is built
  with vite `base: './'`, so its assets resolve at either depth. Tunnel ingress rules added for
  the apex and `www`.
- The terms say the service runs at both hostnames and that they are one instance sharing a
  cache, rate limits and downtime, so choosing one over the other buys nothing.
- README, the deploy runbook, `app/main.py` and `vite.config.ts` name both.

---

## 2.2.4 "Kinetochore"

A second audit suite, covering the readers the pronucleus audit never touched: 163 cases across
the HapMap and VCF paths, against published pedigrees. 132 correct, 6 refused, 25 incorrect, and
the 25 are one failure mode that is documented rather than patched.

### Added

- **`audit/FORMATS.txt` and `audit/formats.py`.** HapMap CEU chr22 (International HapMap Project,
  174 individuals, 52 complete trios, nucleotide genotypes, build 36) and the 1000 Genomes phase 3
  chr22 related-samples VCF, with father and mother taken from the published 1000G pedigree.
  Every one of the 52 documented father-child pairs and 52 mother-child pairs is called correctly,
  at 0.00x to 0.06x of its ceiling.
- A limit fires when the parent's heterozygosity says the marker set is not a polymorphic panel.
  Every rate here is calibrated on common-SNP arrays, where a parent runs 15-19% heterozygous; the
  1000G VCF runs 3.2%, because it is a whole-genome variant callset whose sites are mostly rare
  and cannot show absence in anyone. Those sites dilute the denominator, and the unrelated pair
  reading 4.1% absence on the array reads 0.69% in the callset. Annotated, never adjusted for.

### Known limitation, found here and deliberately not patched

- **A presence call is only as strong as the ceiling it cleared.** The test is one-sided: present
  means the observed absence is no larger than this sample's own noise could manufacture. Where
  the call rate is low that bound reaches 5-22%, and an unrelated adult showing 7-10% absence sits
  under it. All 25 incorrect calls are that, on HapMap individuals typed on one panel only and so
  carrying 35% call rates.

  A guard refusing presence below the call-rate floor was written and reverted: it breaks the case
  the method exists for, Zuccaro embryo A8, a true father-offspring pair at a 53% call rate whose
  9.69% absence is higher than genuinely unrelated pairs at 5.1-5.6%. The second axis does not
  separate them either, since a real biparental embryo sits on the same side as an unrelated
  adult. The margin is printed beside every call for this reason; `audit/FORMATS.txt` sets out the
  evidence in full.

---

## 2.2.3 "Kinetochore"

### Fixed

- **The browser could not read AB-space genotypes at all.** The parser recognised only Axiom's
  numeric codes 0/1/2/-1, so a file spelling its calls `AA`, `AB`, `BB` (which is what most
  exports write) came back as 100% no-call, was excluded on its call rate, and said nothing about
  why. The file was fine; only the dialect was unknown. AB space, `BA`, and the usual no-call
  tokens (`NC`, `--`, `---`, `?`, `.`) are now read.
- **Nucleotide genotypes were silently swallowed rather than refused.** A file writing `AG` also
  read as no-call with no explanation. Which nucleotide is allele A is a per-marker convention
  that has to be resolved by pooling across every sample in a run, which a page reading one file
  at a time cannot do, so these are now refused explicitly on a new **genotype format** gate that
  names the dialect and points at the command line, which does pool.
- The documentation claimed the browser accepted nucleotide pairs. It did not, and now says so.

  Missed by the accuracy audit because all 27 of its arrays are Axiom numeric exports: the audit
  covered the statistics thoroughly and the ingestion layer through exactly one path.

---

## 2.2.2 "Kinetochore"

### Fixed

- **The terms did not contain the promise the rest of the project cites them for.** Five modules
  and the Syngamy documentation tell a reader that "the terms promise nothing about a family is
  submitted or retained". The terms said no such thing, and section 4 said the opposite: "the
  tool takes a variant identifier and nothing else". That has been false since 1.5.0, when the
  carrier genotype loader arrived, and 2.0 widened it from one carrier's VCF to a family's
  arrays. There is now a section 5 that makes the promise in full: files are read in the browser,
  no endpoint accepts a genotype, the report is typeset locally for the same reason, and the
  reader is told how to confirm it in their own network panel. Two of its sentences are pinned by
  `tests/test_terms.py` so a rewrite cannot drop them again.
- **The terms said the tool "cannot determine which parental allele an embryo inherited."** That
  was written for 1.x and Syngamy does exactly that. Section 1 now describes both halves and
  draws the line the docs draw: which parent contributed is reported, which of that parent's two
  chromosomes came through is not.

### Changed

- The one-line description everywhere it appears (meta description, landing page, README) now
  names both halves: marker panels before an experiment, parent of origin after one.
- Terms section 6 names the bundled example arrays as public GEO data served under NCBI's terms.

---

## 2.2.1 "Kinetochore"

A self-audit of the 2.0 surface, going after what had been deferred, asserted without a check, or
left matching a backend that was itself wrong. Five defects, each with a check that fails if it
returns. The accuracy audit re-ran unchanged: 22 correct, 4 refused, 0 incorrect, no case moved.

### Fixed

- **The chrX exemption was inferred from the absent X itself, not measured.** A male offspring has
  no paternal X, so that absence is exempted rather than reported as a loss. The browser decided
  "he sent a Y" from the chrX rate alone, which is circular: on a sample calling no chrY at all it
  exempted a genuine paternal X deletion as ordinary sex determination and reported
  `Y_bearing sperm` on the strength of a chromosome the file does not show. A silent false
  negative on exactly the event this tool exists to detect. Now gated on a chrY measurement, as
  `origin.py` already did; where the X is absent and no Y is called, the loss is reported and no
  sperm type is claimed.
- **The per-locus test used a lower bound as if it were the rate.** The Kothiyal floor is
  explicitly a lower bound on the spurious-violation rate, and the browser left `q` there for
  every sample. A lower `q` shrinks the required run, so the test over-called: at a 5% spurious
  rate it demanded 2 markers where 4 are needed, and at 30% it demanded 2 where 8 are. It now
  fits `q` from the sample's own genome-wide absence rate, which measures exactly that quantity,
  and falls back to the floor only where the paternal genome is absent and the rate would not
  mean that.
- **`r_min` could be 1, making a single marker a significant run.** Below n = 3 the threshold came
  back as 1, so one absent marker cleared it at p = 0.0063 on the strength of being the only
  marker in the window. Contiguity is the entire discriminator and a run of one has none, so the
  threshold is now never below 2, in both implementations.
- **The locus box accepted `chr0` and `chr23` through `chr99`**, which parsed and then matched no
  marker. Restricted to 1-22, X and Y.
- **An unresolved pair discarded the parent it had resolved.** With one parent confirmed and the
  other in the uncalled band, the class must stay unclear, but the note said only "at least one
  parent is unresolved" and dropped a contribution that had been measured with margin. It now
  names which parent is settled and which way.

### Changed

- A self-check fixture asserted the old chrX behaviour: it built a sample with an absent paternal
  X and no chrY and expected the exemption to fire. The fixture now carries the chrY it was
  implicitly claiming.

---

## 2.2.0 "Kinetochore"

### Added

- **An accuracy audit against 27 public arrays**, in `audit/`: the harness, the record as a text
  file, and a report PDF for each of the 26 scored cases. The data is GEO GSE148488, where
  Zuccaro et al. separated the two pronuclei of a fertilised zygote by micromanipulation and
  arrayed each alone, so the answer exists physically before any statistic is applied. Result:
  22 correct, 4 refused, 0 incorrect, and 22 of 23 correct on samples passing every quality
  gate. The egg-donor identification matrix resolved 7 of 7 maternal pronuclei to exactly one of
  four donors, with the runner-up 10 to 18 times further away. The harness drives the same
  modules the browser runs.

### Fixed

- **Zygosity was asserted from heterozygous calls the tool had already declared untrustworthy.**
  Below a 60% call rate erroneous heterozygous calls are common, which is why the ingestion
  gates exclude there and the het-to-hom gate prints "SUSPENDED"; the classifier then read the
  heterozygous BAF band anyway. On three isolated paternal pronuclei at 53.8%, 55.6% and 59.1%
  call rate, each haploid and therefore homozygous by construction, the spurious band of 18-27%
  pushed zygosity to "diploid" and each was reported **biparental**. Zygosity is now withheld
  below that floor in both `origin.py` and `web/src/parentage.ts`, and the class falls to unclear
  with the reason stated. Absence is unaffected and still called, since it is Mendelian: on the
  paternal axis all three were already right, at 0.44x, 0.49x and 0.45x of their ceilings. Found
  by the audit above, which is what it was for.

---

## 2.1.1 "Kinetochore"

### Added

- The report carries the per-locus deletion test when one has been run: the locus, the declared
  event size, and per sample a row per window with the L3 marker count, how many of those the
  sample was called at, total absences, longest run, r_min, p, the maternal mirror run and the
  window's resolution floor, then each window's own note. Methods says the test ran and over
  which windows, since it is written from the run rather than from a template.

---

## 2.1.0 "Kinetochore"

### Added

- **The per-locus deletion test runs in the browser**, as a box under the results taking a
  variant position. Each sample is re-read, one chromosome of it, and scored marker by marker;
  the statistic is the longest run of consecutive paternal absences against the longest run
  independent genotyping error would produce at the same marker count. Three windows, 25 kb,
  10 Mb and the whole chromosome. An optional event size turns an absent run into an explicit
  "below resolution" rather than letting it read as no event.

  It requires the oocyte donor and refuses without her, for the reason the command line refuses:
  absence is Mendelian and needs nothing from her, but presence is an identity claim and needs
  her homozygous for the other allele. Without that nothing scores as present, nothing breaks a
  run, and the statistic has no null.

  The sample is re-read rather than retained, because its genotypes were streamed and discarded
  to keep memory flat.
- Its own documentation section, and the marks linking to the source and the author now appear
  on the Syngamy documentation as they do on the panel documentation.

### Fixed

- **`scoreMarker` scored paternal presence without consulting the mother**, the same defect
  `origin.score_paternal` was fixed for and documents at length. At a paternal deletion the
  embryo is hemizygous for the maternal allele, so wherever a heterozygous mother transmitted
  the allele the father also carries, the call still contains his allele and read as "present".
  That breaks the run at roughly half of all mother-heterozygous markers. The browser statistic
  was unreachable from the UI until now, so nothing shipped a wrong answer, but it would have
  the moment the test was exposed.
- The role control and the stat line shared one row on a 318px chip, so "female" ran off the
  edge and "ambiguous" further still. They stack.
- The two links under the search box had no separator.

---

## 2.0.1 "Kinetochore"

### Fixed

- **The landing page had no route to Syngamy.** The header carries the only link to it and
  the landing page hides the header, so the front page could not reach the other half of the
  tool at all. The hero now links to both halves.
- **The Syngamy documentation route rendered with no site header**, because it fell through
  the `atHome` test and was treated as the landing page.
- **Syngamy was unusable on a phone.** The chip is sized for two-up on a desktop and at 375px
  was wider than the drop zone holding it, so chips spilled past the dashed border and clipped
  their own right-hand column. They now fill the width they have, and the role control and the
  stat line stack rather than fighting for one row.
- **Every table in the Syngamy detail view overflowed its card**, since none had a scroll box.
  All 21 of them now scroll inside their own width, the way the panel documentation's tables
  already did.
- **The footer squeezed the disclaimer into a one-word column.** It and the version and links
  sat at either end of one nowrap row, and the links cannot give way, so at 375px the
  disclaimer was left about 100px. They stack below 700px.
- **The header wrapped to three lines at 375px**, leaving the wordmark, the nav and the status
  badge on separate rows. "Documentation" shortens to "Docs" below 700px, which fits the
  wordmark and both links on one line.

---

## 2.0.0 "Kinetochore"

The structure that couples a chromosome to the spindle and pulls it to a pole.

1.x builds a marker panel before an experiment. 2.0 reads the experiment afterwards: given a
sperm donor's array and a sample's, it reports which parental genome is present and on which
chromosomes. Cas9 cleavage in a human embryo frequently removes the cut chromosome rather than
repairing it, and a wild-type read at the variant site cannot tell a correction from a loss, so
the question is answered away from the cut site instead.

Research use only. Not a clinical diagnostic.

### Added

- **Syngamy**, at `#/syngamy`. Drop array exports, label one the sperm donor, run. Each sample
  is classified androgenetic, gynogenetic, biparental or unclear, per chromosome as well as
  genome-wide. Files are read in the browser and there is no endpoint to send them to.
- **Three measurements, each against a reference derived from the data in front of it.**
  Paternal absence against the rate that sample's own noise can produce, taken as its no-call
  rate times its heterozygous fraction; alleles the donor lacks against half his heterozygosity,
  which is what a second parent contributes; zygosity from the fraction of B-allele frequencies
  in the heterozygous band. No population constant appears anywhere in the calculation.
- **Sperm type from the chrX rate, not chrY.** An X-bearing sperm delivering a complete paternal
  genome leaves no Y, so a chrY test reads it identically to no paternal contribution at all.
  On one four-sample set, three of four had no chrY and the SNP measurement separated them by 15
  to 40-fold. The chrX exemption for a Y-bearing sample is paternal only: a mother transmits an
  X to a child of either sex, so an absent maternal X is a real finding.
- **An oocyte donor array is accepted alongside the sperm donor.** Each sample is tallied
  against both in one pass and classified from the pair. That measures the maternal side instead
  of inferring it, which separates a sample lacking a paternal contribution from one belonging
  to `neither_parent`. The two declared parents are also compared against each other, since two
  arrays of one person pass both tests and yield a confident biparental call.
- **An uncalled band between the noise bound and three times it.** A noisy diploid sample can
  reach the unrelated range against its own parent, so that region is left uncalled with the
  reason stated rather than guessed. The other parent's array measures dropout directly and
  settles it.
- **Six quality gates on every file**, reported and never applied silently: call rate against
  the only published threshold measured on amplified material, suspension of the het-to-hom
  asymmetry below that threshold, genome-wide LOH, an upper bound on heterozygosity that a
  diploid genome cannot exceed, detection of inverted numeric genotype coding, and a sex call
  from the chrX to autosome heterozygosity ratio.
- **Assembly detection from marker positions** against the UCSC chromInfo and gap tables. No
  liftOver is performed and no chain file is shipped: those carry a non-commercial field-of-use
  restriction Apache 2.0 cannot sublicense. A mismatch is reported, and the steps that depend on
  coordinates refuse rather than proceed.
- **A citable report**, typeset in the browser so producing it needs no upload. Per sample: the
  class, every measured rate with the reference it was compared against, a full per-chromosome
  table per parent, the quality read, all six gates, findings and limits. Then Methods written
  from the run, why no confidence percentage is reported, every constant with its provenance,
  and each input file with its SHA-256.
- **Five public example arrays load with one click**, from GEO GSE148488. Every eighth marker,
  no value altered, so the calls are the ones the full 825,657-marker files give.
- **Syngamy documentation** at `#/syngamy-docs`, cross-linked with the panel documentation in
  both directions. Twenty-five sections from the biology through the method, the validation
  record, scope and limits.
- **Backend beyond the browser**: a per-locus deletion test with a run-length statistic, a
  bidirectional cross-check that a real paternal deletion cannot produce maternal-absence calls,
  a sibling haplotype scaffold whose confidence saturates on the genetic map rather than on
  marker count, structural-variant VCF corroboration, and readers for VCF, Illumina FinalReport,
  HapMap raw and PLINK.

### Notes on what it refuses to assume

- **No confidence percentage.** The per-sample likelihood ratio runs past 10^10000, which is
  not a probability any reader should be handed. The honest figure is the method's empirical
  accuracy, 9 of 9 correct on pairs of known relationship, whose 95% lower bound is roughly 72%.
  About 300 consecutive correct calls would be needed to claim 99%. The margin is reported
  instead: the observed rate, the reference, and the ratio between them.
- **The per-locus test refuses to run without the mother.** Paternal presence cannot be
  established at any single marker without someone who could not have supplied the allele. The
  informative set holds absences only, nothing can break a run, and the statistic has no null to
  be significant against. Dropout cannot be fitted either. In that mode a normal chromosome 20
  produced a run of 3 across 35 Mb at p = 2.5e-07, so it is refused rather than reported. Parent
  of origin, sperm type and segmental loss need no mother and are unaffected.
- **An "unrelated" mismatch rate is not a constant.** Measured across more than fifty pairs it
  ranged from 6.8% to 50.3%, and the populations overlap: one embryo scores 9.69% against its
  own father while genuinely unrelated pairs score 5.1% to 5.6%. No threshold separates them,
  which is why the reference is derived per sample.
- **A run is split where markers are too far apart.** Consecutive absent markers separated by
  more than ten times the local informative spacing are not one run: adjacency in the
  informative subsequence is not adjacency on the chromosome. Whole-chromosome loss is keyed to
  every informative marker being absent, not to the longest run, since marker deserts fragment
  runs.
- **Segment rates are anchored to the calibrated bound**, not estimated from the distribution of
  observed window rates. Percentile estimation finds two populations in a uniform genome: it
  carved a sample present everywhere at 0.16% into 28 "absent" segments, and a gynogenetic
  sample into 720.
- **Nucleotide-to-AB harmonisation pools across every file in a run**, never per file, or two
  files can be assigned opposite conventions at the same marker.
- **A no-call is excluded, not counted as evidence** in either direction. The heterozygous BAF
  band is counted before that exclusion, because a marker whose genotype failed still has an
  intensity reading and a dropped heterozygote sits mid-band.

### Fixed

- **Presence was scored without consulting the mother.** At a paternal deletion the embryo is
  hemizygous-maternal, so wherever a heterozygous mother transmitted the same allele the marker
  read as present. This broke runs at roughly half of maternal heterozygous markers and reported
  a chrX loss of 4,324 markers as a run of 21.
- **Markers that could not speak were folded in as presence** rather than excluded, which broke
  runs and inflated the marker count at the same time.
- **The segment reported was the longest anywhere on the chromosome**, not the one at the
  variant, and counted scattered markers. A normal chromosome 3 was reported at 7,499 markers
  and 359 standard deviations.
- **A runner re-derived significance downstream of the gates**, bypassing the no-mother refusal
  and reporting paternal loss on three samples whose paternal genome was present at 0.2%.
- **The chrX probe offset was applied when deriving the compression factor but not when
  interpreting chrX**, so a male's single X read as two maternal copies.
- **`--product` defaulted to a named array**, applying one platform's chrX offset to another's
  data and giving 0.456 against a measured 0.59. It defaults to unspecified.
- **A wide multi-sample genotype file was read as its first sample only**, dropping the rest
  silently.
- **The column alias resolver did not strip hyphens**, so `Allele1 - Forward` missed its alias,
  the tab split failed, and parsing fell through to whitespace.
- **Findings were listed as limitations** in the report's Methods. Limits are now held
  separately from notes.
- **Assembly detection ran only inside the scaffold branch**, so a parentage-only run reported
  "undetermined" for data proven GRCh37.
- **The browser and Python disagreed on the BAF band.** The browser returned early on a no-call
  before counting it, giving 0.92% against Python's 1.27% on the same file and shifting the
  noise ceiling with it.
- **The Kothiyal 2019 reference named the wrong journal.** It is Journal of Computational
  Biology 26:405-419, not BMC Genomics.

---

## 1.5.1 "Cohesin"

### Changed

- The carrier's genotype file can be chosen in Manual input, alongside the query, instead of
  only after the panel is built. Only the file reference is taken there: nothing is read
  until the build lands and there is a window to filter it against, and it is applied to the
  results without a second action.

---

## 1.5.0 "Cohesin"

The ring complex holding sister chromatids together, without which crossovers could not be
resolved correctly.

### Added

- **Load the carrier's own genotypes.** Every figure this panel ranks on is a population
  prior, and the lab protocol's first two steps are to genotype the carrier and drop the
  markers where that person turns out to be homozygous. Given their VCF or SNP-array export,
  the panel now does that step: each marker is reported as heterozygous, homozygous, no-call
  or not in the file, the shortlist can be filtered to the markers that are genuinely
  informative, and the both-sides coverage count is recomputed on a real person instead of an
  average.

  **The file never leaves the browser.** There is no endpoint that receives it. A carrier's
  sequencing file is the most identifying artefact in this workflow, and the terms promise
  nothing about a family is submitted or retained, so it is parsed in the page and forgotten
  when the panel is cleared or rebuilt. A whole-genome VCF is streamed and discarded down to
  the panel's window as it reads, so it never has to fit in memory.

- **A documentation chapter on choosing between SNP array, WGS and PCR**, and on what the
  carrier's data does and does not settle. The short version: sequencing the carrier tells
  you which markers are informative and finds rare ones no chip carries, but it does not give
  you phase. Short reads span hundreds of bases and these markers sit tens of kilobases from
  the variant, so which marker allele rides on the pathogenic chromosome still comes from an
  informative relative, from reads long enough to span the distance, or from typing gametes.

### Notes on what it refuses to assume

- **Not in the file is kept separate from homozygous reference.** A variants-only VCF omits
  every site where the sample matches the reference, so absence usually means homozygous and
  in an uncalled region means nothing was measured. Only an all-sites or gVCF file
  distinguishes them, so absent markers are left uncounted rather than assumed either way.
- **A file on another assembly is detected rather than silently mis-joined.** Markers are
  matched by rsID first, because an rsID survives an assembly difference and a coordinate
  does not. If the rsID-matched sites systematically disagree about position, the page says
  the file looks like a different build and reports the offset.
- **Heterozygosity is strand-invariant**, which is what makes array exports safe to read
  here: a file reported on the opposite strand turns A/G into T/C and both are still one
  allele of each. The A/T and C/G ambiguity that plagues array merging cannot turn a
  heterozygote into a homozygote.
- Loaded genotypes are cleared on a rebuild. They were filtered to the previous window as
  they streamed, so against a wider one every new marker would read as "not in the file" when
  the file may well carry it.

---

## 1.4.1 "Tetrad"

### Changed

- The genome check moved to the top of the primer dialog, above the settings grid: it is the
  one control in there anyone acts on, and the knobs below it are mostly provenance.
- Starting a check now scrolls the page back to the build log and opens it, from either
  button. The log sits above the table, so a reader who started the run from a marker row was
  being shown progress somewhere they were not looking.

---

## 1.4.0 "Tetrad"

The four chromatids of a paired chromosome, the unit crossover acts on.

Responds to an independent scientific audit of the methodology. The audit found no errors in
the statistics or genetics that is implemented; every change here closes a gap between what
the method assumes and the inputs it silently accepted, or improves a displayed figure.

### Added

- **mtDNA variants are refused at resolve time.** Mitochondrial DNA is maternally inherited,
  does not recombine, and is often heteroplasmic, so "which parental chromosome an embryo
  inherited" is undefined and flanking-SNP linkage does not apply. A chrM variant used to
  resolve like any point variant and hand back a thin panel that looked valid; it now returns
  a clear statement that the method does not apply.
- **A "cases this tool does not handle, or cannot detect" section** in Scope and limits:
  de novo variants, uniparental disomy, mosaicism, consanguinity, mtDNA and repeat expansions,
  each with the consequence for a panel. The tool never meets the family, so most of these it
  cannot detect and names as the reader's to rule out.

### Fixed

- **gnomAD exome frequencies were fetched and then dropped.** The single-variant query asked
  for both the genome and exome callsets but read only the 76k genomes. A coding pathogenic
  variant, which is most of this tool's input, carries its frequency in the 730k exomes and
  can be sparse or absent from the genomes, so the rarity card could read blank while a good
  exome answer sat unused. The verdict and the headline frequency now take whichever callset
  observed more chromosomes at the site, and both callsets are shown on the card and in every
  export. The LD-usable verdict is unchanged: it remains gated on the 1000 Genomes count.

### Changed

- The primer field reference notes that the 69 C Tm default is deliberately stringent, higher
  than routine genotyping, and can be lowered toward 60 C for a standard single-anneal assay.
- The star's 1 Mb clause is explained as inert at the default window (every candidate is
  already well inside 1 Mb) and binding only if the window is widened; ESHRE's 2 Mb
  "acceptable but not advisable" allowance is noted for that case.

---

## 1.3.4 "Diakinesis"

### Added

- A dismissable nudge at the top of the candidate list recommending a genome check on the
  primer pairs. Its button reopens the build log and streams the verification into it.
- A download button on the build log that exports it as a .txt, with a footer line carrying
  the release, job, instance URL, timestamp and data versions.

---

## 1.3.3 "Diakinesis"

### Fixed

- **Free-text descriptions never resolved.** The intent prompt forbade the model to recall
  an rsID, so "the sickle cell mutation, in Africans" returned no variant at all rather than
  rs334. It now names the standard identifier for a described variant, while still never
  emitting a coordinate; the identifier is confirmed by live lookup as before.

---

## 1.3.2 "Diakinesis"

Terms of use, brought up to what the app actually does. No behaviour changed.

### Changed

- **The terms disclose the two things that leave the server carrying your input**, which they
  did not. Section 5 was called "Third-party data" and described only what is retrieved:
  - **Free text with no identifier in it is sent verbatim to Anthropic's API**, to be read by
    a small model. Text containing an rsID, HGVS or a ClinVar accession is read here by a
    regular expression and reaches no model at all. The distinction is the whole privacy
    story of that box, the terms already asked readers to keep identifying information out of
    it, and they never said why.
  - **Primer sequences are sent to UCSC** when the check is asked for. They come from the
    reference genome rather than any sample, and nothing else from the query goes with them.
- **Primers are in the terms at all now.** Section 1 described a tool that proposes markers,
  which stopped being the whole of it in 1.3.0. A primer pair is a candidate in the same
  sense a marker is: nothing here has run a PCR, a reference is not a patient's genome, and a
  private variant under a primer site causes exactly the dropout the design cannot see.
- **Section 6 no longer lets "no tracking" imply "nothing leaves".** They are different
  claims and only the first was true. It says which is which, and that browser history stays
  in the browser.
- **Section 7 states the one dependency that is not Apache 2.0.** primer3-py is GPLv2, which
  is why it is optional and absent from the default image. Redistributing an image built with
  it switched on makes the combined work GPLv2. Previously the terms said "redistribute,
  including commercially" with no qualification.
- **Section 4 says whose UCSC quota the check spends.** It is the operator's key answering for
  UCSC's published limit, which is why it is rate limited per client.

### Added

- A check that fails when the app reaches a host the terms do not name, and when the terms
  stop saying a load-bearing thing ("not a clinical diagnostic", "candidate", "wet-lab",
  "GPLv2"). Prose drifts from code silently and in one direction: someone adds an outbound
  call, every test passes because the call works, and the terms keep saying it does not
  happen. It reads the rendered prose with the hrefs stripped, because the first version of
  it passed a gutted disclosure on the strength of a leftover URL.

---

## 1.3.1 "Diakinesis"

A patch, so it keeps 1.3.0's name.

### Added

- **The primer design is now reachable from Manual input**, which is where it was always
  meant to be. Pick which markers get a pair, and set every constraint they are designed
  under: the melting temperatures, the lengths, the composition, the product window, the
  reaction conditions the Tm is only meaningful beside, and the mask. The form seeds from the
  server's own numbers, so what is on screen is what the build is asked for.
- **A checkbox to check every pair against the genome as part of the build**, for anyone who
  would rather wait once than build and then press a second button. It states the cost beside
  itself: UCSC allows one query every 15 seconds, so it adds about 15 s per designed pair
  while the panel alone takes 20 to 60 s. The build log names each verdict as it lands, so a
  long run is not a blank wait, and it is off unless ticked, every time.

### Fixed

- **The Manual input primer form has never rendered for anyone.** It draws only against a
  server that states its defaults, which is right: the numbers must come from the engine that
  will use them rather than a copy in the browser. `/api/health` never sent them. So the
  section was gated on a field that did not exist, the flag beside it said primers were
  enabled, and the form was simply absent with nothing anywhere reporting a problem.
- **The build log claimed no pair had been checked, in builds that then checked them.** It
  said "none has been checked against the genome by this build", which was true when the
  design emitted it and false three lines later once the bundled check ran in the same job,
  about the same pairs, into the same console. It says "not yet" now, which is true on both
  paths.
- **A bundled check held a build slot for the whole of its wait**, found by watching the live
  site refuse an ordinary build while two ticked boxes sat waiting on UCSC. The slot bounds
  BUILDS, as its `MAX_CONCURRENT_BUILDS` name says, and the default is 2: two people ticking
  the box therefore blocked panel builds for everyone, for minutes, over work that had
  finished. The slot goes back when the build ends now. Verification stays bounded by the
  per-client budget and by the process-wide gate, which is what was keeping it polite anyway.
- **`npx tsc --noEmit` type-checks nothing, and the README told you to run it.** `tsconfig.json`
  is a solution file with `"files": []` and two references, so bare tsc has no inputs and
  exits 0 over a codebase that does not compile: only `tsc -b` follows the references. Proven
  rather than reasoned: a deliberate type error exits 0 under `--noEmit` and 2 under `-b`.
  Nothing shipped broken, because `npm run build` runs `tsc -b` and was in the same gate, but
  the line above it was pure reassurance. The README now documents the command that checks.

### Security

- **The bundled check spends the verification budget, not the build budget.** Both routes to
  UCSC reach one published daily quota, and the budgets differ: 20 builds a client per window
  against 4 verification runs. Charged only to the build budget, the checkbox would have been
  a five-fold rate-limit bypass on someone else's server. It is charged to the same key the
  button spends, and a test fails if that stops being true.

---

## 1.3.0 "Diakinesis"

The final condensation of prophase I, chiasmata still holding the homologues.

### Added

- **Primer design for the markers that meet the flanking criteria.** Each one gets a
  candidate FWD/REV pair for genotyping it by PCR, designed by primer3 against a reference
  template fetched around the marker, and folded into the PDF, CSV, XLSX and JSON. Defaults
  are 20 to 35 bases, GC 40 to 60%, Tm 69 C, product under 600 bp, GRCh38. Every field is
  settable, from the panel's own primer box or through the API, and the settings a panel was
  built under travel with it as provenance.
- **Common variants are masked out of both primer sites.** This is the part that matters. A
  primer sitting on a common SNP fails to bind in exactly the carriers who have it, their
  allele goes unamplified, and a heterozygote is read as a homozygote: allele dropout, which
  is silent and yields a genotype that looks clean and is wrong. The pool of markers is also
  the pool of hazards, so every gnomAD variant at or above the mask floor is excluded from
  under both primers, and the marker itself sits in the product under neither of them.
- **Optional verification against the whole genome**, through UCSC In-Silico PCR, behind a
  button and never part of a build: UCSC publishes one request every 15 seconds, so a public
  URL that verified on its own would spend the owner's quota on visitors who never looked at
  the result. It needs a UCSC API key (`deploy/README-deploy.md` says where to get one) and
  reports itself unavailable without one rather than offering a button that cannot run.
- **A primers chapter in the documentation**: every field with its bounds and what it
  constrains, what each warning means, what a clean in-silico result is and is not worth, and
  the steps to obtain a key. The field table is generated from the form's own field list, so
  a knob that reaches the form and not the docs fails the typecheck.

### Fixed

- **A verified pair kept the warning saying it had never been verified.** The verdict was
  written onto the pair as a bare state code, so a pair UCSC had called dangerous still
  carried "NOT CHECKED AGAINST THE GENOME" on the same row, and the reader was being asked to
  trust exactly one of two statements the same document made. The finding itself, naming the
  loci and their positions, reached only the build log: never the table, never the PDF. The
  verdict's own words are welded onto the pair now, and the caveat with them.
- **The check asked UCSC a narrower question than UCSC asks itself, and reported the answer
  as genome-wide.** Max product size went out at 1000 bp against hgPcr's own default of 4000.
  Measured against the live endpoint with a pair known to give one 549 bp product: at 400 the
  product is not reported at all. So the field bounds the search rather than filtering its
  result, and a pair whose second locus amplified between 1001 and 4000 bp came back holding
  one product, classified clean, and printed VERIFIED CLEAN (in silico) on a filed PDF. That
  is the multi-locus pass this lane exists to prevent, reached without a single component
  behaving incorrectly. It cut the other way too, since a design may ask for a 3000 bp
  product: that product could not be reported, and its absence classified as "found no
  product, do not order", which was our own request accusing a good pair. The question is
  UCSC's default now, the design has a server-side ceiling below it, and a check fails if
  either moves.
- **A page listing two loci could come back holding one.** The parser decided a line was a
  FASTA header before stripping the HTML around it, so a header whose ">" UCSC had wrapped or
  escaped was skipped as quietly as a sequence line, without being marked unreadable. The
  remaining product then classified as a clean single band. Tags come off and entities decode
  before that test now, so the test and the regex read the same text.
- **A UCSC timeout was rendered as a UCSC finding.** The page had no name for the state that
  means "asked, and the answer could not be read", so a timeout or a spent quota drew a red
  DANGER badge and a banner reading "In-silico PCR contradicts N primer pairs in this panel",
  over pairs UCSC had never answered about. The PDF from the same job id said NOT VERIFIED
  for all of them: one job, two documents, and two different instructions. Not verified is
  now neither dangerous nor clean, which are separate questions and were fixed together, as
  answering only the first renders a quota stop green.
- **The daily verdict on alt scaffolds overclaimed.** A hit on chr6 and one on a chr6 alt
  haplotype are usually one locus reported twice, and the note asserted flatly that the pair
  amplified more than one locus and must be redesigned. It says products rather than loci
  now, and names the ambiguity where alt or fix scaffolds are among the hits. The state stays
  DANGER: hgPcr cannot separate a redundant alt copy from a real second locus on that
  haplotype, and guessing toward clean is this tool's worst direction to guess in.
- **A UCSC key in `.env` never reached the container.** Compose reads that file for variable
  substitution and passes nothing it was not asked to, and the key was not named in the
  environment block. The key could be correctly generated, correctly stored, and the feature
  would still report itself unavailable, with no error anywhere: every layer behaving exactly
  as written. A test now fails if any key-shaped variable the app reads is not forwarded.
- **The JSON export shipped the shortlist twice in two different shapes.** The file writes it
  once inside `candidates` and again whole under `recommended`, and only the first was given
  its primers and its ancestry column, so a reader who took the shortlist rather than
  filtering the candidates got the lesser copy of the same marker.

### Changed

- **Every repeated block on the panel is one sentence and a link now.** Each warning carries
  a short form and a long one, written beside each other in the engine so they cannot drift.
  The table takes the short one; the exports keep the full wording, because a filed page
  cannot follow a link. The not-checked warning goes from 242 characters to 67, the mask note
  from 244 to 71, and the primer form's six paragraphs of constraint move to the docs.
- **The primer box is collapsed behind a thin line.** Four lines of detail under every marker
  is a table nobody reads. What collapsing never does is hide a finding: a dangerous pair and
  a failed design open themselves, the verdict sits on the summary line either way, and the
  panel-level alert lists every dangerous pair regardless of what is open.
- **The star's key is a key again.** It said what the criteria were, next to every star,
  beside a hover that already said it; it names the claim now and links to the chapter.
- Documentation section numbers derive from the section list, the way citation numbers always
  have. Inserting a chapter used to renumber eight headings and strand every cross-reference
  that named one.
- primer3-py stays an optional dependency and is not in `requirements.txt`: it is GPLv2 and
  this repo is Apache 2.0. Absent, panels build exactly as before and simply carry no
  primers. Nothing from kent/isPcr is vendored; verification calls UCSC's hosted service.

### Known limitations

- Verification is bounded at 4000 bp, which is UCSC's own bound and not an exhaustive search:
  a second locus amplifying wider than that is still invisible to it.
- A clean in-silico result is not a wet-lab validation. It does not model cycling conditions,
  and it cannot see a carrier's private variants under a primer site, which cause dropout
  exactly where the reference cannot show it.

---

## 1.2.0 "Zygotene"

The substage where homologues find each other and begin to pair.

### Added

- **Local query history.** Click the search box and your previous queries drop down, each
  with its candidate count and an x to forget it, plus Clear all. It lives in your browser
  and goes nowhere: no account, no server, nothing about anyone leaves the page. Storage
  that is full, disabled, or holding another version's garbage degrades to an empty list
  rather than taking the search box down.
- **A star on markers that meet ESHRE's flanking criteria**, in the table and in all four
  exports. What it means, exactly: within 1 Mb of the variant, no recombination hotspot in
  between on an assessed map, and no position disagreement between gnomAD and Ensembl. It
  is a structural check, not a ranking, and not a claim about any carrier's genotype.

### Why the star is a predicate and not a top-three

**There is no convention for a "top 3".** ESHRE's PGT-M recommendations (doi:
10.1093/hropen/hoaa018) do give real numbers, but a different shape: *at least three SNPs
proximal and three SNPs distal*, within 1 Mb of the variant, avoiding known hotspots. That
is six, per side, as a minimum count, not a ranking. ESHRE's own informativity rank is a
function of the couple's actual genotypes, which this tool structurally cannot have: it
proposes candidates and has no genotypes (R3), and 2pq is a population prior, not a
genotype (R4). So "strongest" in the only sense the convention defines is not computable
here, and claiming it would be an invented recommendation printed on a filed PDF.

**Starring the top 3 by the existing rank would have been actively harmful.** On the
reference ABCC8 panel, the three top-ranked markers are all on the same side and all three
have a recombination hotspot between them and the variant. Starring them would endorse
exactly the markers the tool's own coverage flag warns about, and violate the both-sides
rule on the same page that states it.

**So the star is a predicate, and it is not capped.** Every marker meeting the criteria
gets one; the count per side is reported against ESHRE's minimum of three. A cap would
force an ordering, and no source gives an exchange rate between heterozygosity and
distance: on the reference panel, capping at three would star a marker 24 kb away over one
125 bp away on a 0.01 difference in a population prior. The five unstarred qualifiers are
not worse than the three starred ones, and a star saying so would be a fabrication.

The label says "meets ESHRE's structural criteria". It never says strongest, best, or
recommended.

### Fixed

- **The history dropdown could submit a variant you never named.** Highlight a row, type a
  character that filters it out, press Enter, and the panel built for whichever row had
  slid into that index: a plausible wrong answer, produced silently. The selection index is
  a reference into a list that re-filters as you type, and resetting it in an effect runs
  after paint, which is too late for the keypress that follows. It resets in the same
  handler that changes the list now.
- **The arrow keys did nothing.** The input's own key handler and the combobox's collided,
  and the result was neither: navigation silently never ran. All keys are handled in one
  place now.
- **The export contradicted itself on 349 rows.** The criteria are evaluated over the
  shortlist, and the column is written for every candidate, so ~1,200 markers nobody judged
  shipped a False verdict. Several were nearer the variant than the starred ones, and the
  printed note routes a reader looking for a fourth marker straight into them. The field is
  None where it was not assessed, and the column prints empty: False is an assertion,
  absence is not.
- Four documentation cross-references named the wrong section number after the previous
  release inserted a section. A check now derives them from the section list.

---

## 1.1.0 "Leptotene"

The beginning of prophase I, where chromosomes first condense into threads.

### Added

- **A build log you can watch.** A dropdown beside the progress bar streams what the engine
  is actually doing, one tagged line per event: `[FETCH]` a request going out, `[CACHE]` one
  answered from disk and how old that answer is, `[INFO]` a count, `[SKIP]` sites dropped
  and why, `[WARN]` something worth reading, `[DONE]` the summary. It survives the build
  finishing, which is when it is most useful, and a late subscriber replays the whole log
  rather than joining midway.
- **The monogram on the PDF and XLSX masthead**, rendered from `web/public/favicon.svg` at
  export time. Page 1 only.

### Fixed

- **A provider URL, with its query string, reached the browser.** A failed call raised an
  error built from the raw URL, and the job's error text is both shown to the user and
  appended to the build log. The NCBI api_key travels in that query string. It did not in
  fact leak, but only because the URL was truncated at 80 characters and the key happens to
  be appended last: two couplings nothing asserted, one reordered dict away from publishing
  the key. The error is now built from the same scrubbed label the log uses, and a check
  drives a real failing call and asserts the key is absent.
- **`NM_000518.5(HBB):c.20A>T, MAF at least 0.1` did not resolve.** The search box decided
  client-side whether text was an identifier, and its HGVS test was a substring search, so
  any text containing `:c.` was posted whole, modifiers included, as the variant. The
  predicate is deleted: everything goes to the one parser, which is the server's. This also
  fixes free text being switched off refusing `rs334 in Europeans`, which needs no model.

### Changed

- The monogram's geometry is checked across its three copies. `favicon.svg` feeds the tab
  icon, the PDF and the spreadsheet; `Mark.tsx` holds a copy because its colours are CSS
  variables and the favicon's must be literal. Changing one used to leave the other alone,
  silently, with every gate green. The colours are still allowed to differ; the geometry is
  not, and a check says so.
- The log tag set is likewise pinned across Python and TypeScript.

### Performance

**No speedup shipped, and the measurements say why.** A cold build is 6.8 s, about 70% of it
waiting on the network. The candidates worth trying were tried and rejected on evidence:

- **The fan-out is CPU-bound, not network-bound.** Against a warm cache, one worker takes
  1,027 ms and eight take 887 ms: eight times the workers buys 1.16x. `json.loads` is 40 ms
  per 3.7 MB chunk and does not release the GIL, so ~890 ms is a floor no worker count
  removes. Sixteen workers is slower than eight.
- **Bigger requests are not fewer costs.** 40 kb chunks (7 requests instead of 13) measured
  *slower* than 20 kb (1.35 s vs 1.07 s): the cost is bytes and parsing, not round trips.
  10 kb was far worse at 5.97 s.
- **Nothing unread is worth removing.** The region query reads every field it asks for.
  Populations are 68% of the payload and gnomAD cannot filter them server-side.
- `annotate` is linear, and gnomAD's own run-to-run variance on an identical request is 8x
  (1.07 s vs 8.67 s), which is the noise floor any future claim here has to clear.

---

## 1.0.0 "Synaptonemal"

The protein scaffold that zips paired homologues together along their length.

**First public release.** Open source under Apache 2.0. Versioning restarts at semver here:
the 1.x/2.x numbers below were pre-release build counters, not compatibility promises.

### Added

- **Free-text input.** Two paths, and the difference is not cosmetic. Name an identifier
  (rsID, HGVS, ClinVar accession) and a regex reads it: no model runs, and it costs
  nothing. Describe the variant in words and a model is asked what you meant, answering
  from its own knowledge.
- **The model cannot touch a coordinate.** The typed query it fills has no field for a
  chromosome, position, strand or allele, so the coordinate on every panel came from a live
  lookup regardless of how the variant was named. This is a property of the code rather
  than a promise, and it was attacked specifically: extra coordinate keys are ignored,
  coordinate-shaped identifiers are refused by an allow-list.
- **A gene cross-check.** Any gene symbol in your own text is compared against the gene of
  the record that actually resolved. Ask for "the ABCC8 splice mutation" and get a variant
  in HBB, and it refuses rather than warns. It reads the resolved record, never the model's
  claim about itself: a model wrong about the variant can be wrong about its gene in the
  same breath.
- **Model provenance in every export.** A panel whose variant a model chose says so, in the
  CSV, JSON, XLSX and on the PDF variant card, and it says the two things separately: the
  model chose which variant, and it did not supply the coordinate. Previously a panel built
  from prose was byte-identical to one typed by hand.
- **An intent cache**, so identical prose is not billed twice.
- Documentation of all of the above, including what it cannot do.
- The monogram: a geometric sans ring around the wordmark's own serif M, outlined from the
  Merriweather glyph rather than redrawn. It is the site icon and the corner mark.
- Rotating input examples cycling the accepted forms.

### Fixed

- **The search box accepted one character at a time.** A component declared inside the
  render body took a fresh identity every render, so React remounted the whole subtree on
  each keystroke and destroyed the input's focus.
- **The gene extractor could not read `C9orf72`.** HGNC's Cxorfy convention carries
  lowercase, which a capitals-only pattern reads straight past. The commonest genetic cause
  of ALS was invisible to the safeguard.
- `MB` was excluded as a unit, but it is also myoglobin. Units are now stripped before gene
  extraction, so `500MB` is a window and `MB` is a gene.
- Ordinary report jargon (NIPT, CVS, WES, MLPA and others) read as gene symbols and refused
  correct answers, after billing for them.
- Legitimate aliases (`SUR1` for ABCC8, `ND1` for MT-ND1) were refused with no way past it.

### Known limitations

A gene symbol you lowercase reads as prose, and the cross-check then stays silent exactly
as if you had named no gene. Capitalisation is what separates a symbol from an ordinary
word here, and relaxing that makes every English word a symbol. Documented rather than
papered over.

---

## 2.4 "Bivalent"

A pair of synapsed homologues, held together by the chiasmata between them.

### Fixed

- **Default ranking used a best-case-over-ancestries statistic.** With no ancestry selected
  (the default), the primary sort key was the maximum expected heterozygosity across the
  eight gnomAD populations. That is not a prior for any carrier: it is the order statistic
  E[max], upward-biased by construction (+0.09 absolute, +26% relative on a representative
  locus), and each marker's figure came from a different population, so the ranking assumed
  a different ancestry for every row. One marker held a core slot advertising 2pq = 0.4884,
  the African value; for a Northern European family the same marker is 0.1418, and for an
  East Asian family 0.0241. Ranking is now on the global 2pq. On the reference panel this
  changed 15 of 20 shortlisted markers, took the count of markers with 2pq < 0.10 from 2 to
  0, and improved the worst-case marker in all eight populations.
- **Every export declared a ranking basis the engine did not use.** Exports stated "global
  2pq prior" while the sort keyed on the max across populations. The engine now names the
  quantity that produced the order, and exports render that name verbatim rather than
  restating it, so the two cannot drift apart again.
- **The PDF still labelled flanks `tel`/`cen`.** The engine moved to `lower coord` /
  `higher coord` in 2.3, but the PDF recomputed the label from the sign of the distance.
  This is the artifact that gets printed and filed, and it was anatomically inverted for
  every gene on a q arm.
- **`rs334` (sickle cell) resolved to a haplotype record.** ClinVar's relevance ranking
  returns eight haplotypes ahead of the HbS allele, each pairing it with a second HBB
  variant. A haplotype record resolves cleanly by returning one constituent's position with
  the haplotype's own classification, so the canonical sickle cell rsID produced a panel
  labelled "other" instead of "Pathogenic". Records describing combinations of variants
  (haplotype, compound heterozygote, diplotype) are now recognised by object type and
  refused with an explanation.
- **The client-side ranking mirror diverged from the engine.** `rank.ts` still fell back to
  the max-across-populations statistic, so the table headed "Engine rank" was ordered on it
  for all candidates while the PDF reported a different basis. Screen and export now agree
  for every row.
- **Concurrent builds shared one fetch ledger.** The ledger records the age of the data and
  supplies `queried_utc`, which is the data date printed on a filed export. It was a module
  global: thread-safe, but not build-safe. Two builds in one process shared it, so an export
  could carry the other panel's data date. Now scoped per build.
- **A failed ClinVar efetch was diagnosed as a structural variant.** The failure was
  swallowed, and the empty alleles fell through to the copy-number guard, refusing a routine
  SNV with a message that contradicted itself: "is single nucleotide variant (...), which
  has no single reference and alternate allele". An efetch failure now reports as a failed
  lookup.
- **A position disagreement between sources produced no panel-level warning.** The
  cross-check runs after selection, so a marker gnomAD and Ensembl place differently was
  already shortlisted when the disagreement surfaced, and it appeared only in that marker's
  own cell. It now raises a coverage flag, which renders in the UI and in all four exports.
- **Panels reported the live server's Ensembl release, not their own.** One frozen panel
  could render two different release numbers as the server moved on beneath it. The release
  is captured at build time; when it is unknown the panel says so rather than borrowing the
  current value.
- **PDF printed an allele frequency at full float precision** (`0.012718820176763365`) four
  lines above its own banner rendering the same quantity as `1.27e-2`.

### Changed

- Module self-checks now run under `pytest`. They existed and were never executed: with the
  ranking defect reintroduced, the suite still reported 54 passed. It now fails.
- Both test factories (Python and TypeScript) defaulted the global and max-across-population
  heterozygosities to the same number, making the two quantities indistinguishable by
  construction. They now default apart.
- Tests: 49 to 69.

---

## 2.3 "Holliday"

The four-armed branched junction intermediate that resolves into a crossover.

### Fixed

- **Genetic distances beyond the end of a chromosome map were clamped and reported as
  measurements.** The bundled maps stop short of the telomeres (chr11 by ~11 kb), and
  positions past the last map point collapsed onto a single value, reporting cM = 0 and
  recombination fraction = 0 sourced as deCODE. A fabricated theta of zero does not merely
  mislead: it ranks as a perfect marker. Positions off the map are now extrapolated at the
  documented 1 cM/Mb fallback and flagged as approximate.
- **`chrX` was labelled sex-averaged and is not.** The bundled chrX reads 1.297 cM/Mb
  against 1.175 (chr7) and 1.160 (chr8). A genuinely sex-averaged X must sit well below the
  autosomes, since paternal meioses are a third of X transmissions and contribute no
  crossovers outside the pseudoautosomal regions. It is the female map: correct for a female
  carrier, and an overstatement for a male one. It now says which.
- **Hotspot detection was off by one on both bounds**, so a hotspot outside the queried span
  could be reported for it.
- **Coverage was blind to recombination.** Eight of twenty shortlisted markers on the
  reference panel have a recombination hotspot between them and the variant, and coverage
  reported no flags. Sides are now judged on markers clear of an intervening hotspot, with
  "not assessed" kept distinct from "clear": an unmapped chromosome can never reach the
  hotspot threshold, so a negative there means unmeasured, not safe.
- **`tel`/`cen` flank labels were anatomically inverted for every gene on a q arm** (roughly
  60% of the genome, including CFTR, BRCA2, SMN1, MECP2, FMR1, and every acrocentric, whose
  p arm carries no genes). No centromere table exists anywhere in the tool, so it was
  asserting an anatomical fact it never looked up. Sides are now named by coordinate.
- **gnomAD quality filters were never requested**, so 33 QC-failed sites in the reference
  window were offered as candidates. Artifacts have inflated heterozygosity, so the ranking
  actively preferred them.
- **A failed 1000 Genomes join was reported as a count of zero**, i.e. a lie about the data
  rather than an admission that the lookup failed. Now reported as unavailable.
- **ClinVar significance was read from the first RCV rather than the variation-level
  aggregate**, and the significance and review status were scraped by independent passes, so
  they could describe different submissions.
- **`queried_utc` was wall-clock "now" even for a panel assembled entirely from cache.** The
  panel now reports both when it was built and how old the data in it is.
- **The natural-language rate limiter billed for refusals.** The global cap was tested after
  the model call, so a rejected request still cost money.
- **The free-text parser took the first identifier and discarded the rest**, so "not
  rs1801133, I mean rs151344623" built the wrong panel. It now refuses and names both.
- **The locus track's y-axis named an ancestry while plotting the global statistic.**
- **Rate limiting could be bypassed** by supplying an `X-Forwarded-For` header, since
  proxies append to it and the first entry is client-controlled.
- Error paths no longer leak module paths or raw provider errors to the browser.

### Added

- Mobile support.
- Recombination hotspot annotation and per-marker genetic distance from the bundled deCODE
  2019 map.

---

## 2.2 "Pachytene"

The substage of prophase I in which crossing over actually occurs.

### Fixed

- **Free text resolved to an arbitrary real variant.** ClinVar's esearch is a full-text
  search, and the code took the first hit on trust, so the phrase "a pathogenic variant"
  returned a real gene at real coordinates with no indication anything was wrong. Only
  identifiers are accepted now, and every candidate record must reconcile against the query
  before it is used.
- **The first ClinVar hit was not always the right one.** `rs1801133` (MTHFR, chr1) returned
  four records, the first of which is a CPS1 variant on chr2. All candidates are now fetched
  in one call and the one that reconciles is used.
- A copy-number variant resolved to empty alleles and a nonsense gnomAD identifier rather
  than being refused.
- The variant could appear as a marker for itself.
- Strand was reported as plus when it was unknown.
- A failed gnomAD chunk returned an empty list instead of raising.

### Added

- Genome-wide verification sweep against live data across all 23 chromosomes, both strands,
  and SNV/insertion/deletion/duplication classes, cross-checking every coordinate against a
  second source.

---

## 2.1 "Synapsis"

The lengthwise pairing of homologous chromosomes that has to happen before they can exchange
anything.

### Changed

- Two-tone wordmark.
- Removed the decorative helix.

---

## 1.0 "Chiasma"

The X-shaped point where two chromatids have crossed over: the visible evidence of
recombination, and the event every flanking marker is chosen to bracket.

First complete release.

### Added

- Panel builder: resolve a variant, enumerate common flanking SNPs from gnomAD, rank by
  expected heterozygosity and proximity, select a balanced panel covering both sides.
- Coordinates from live APIs only (ClinVar, Ensembl, gnomAD). No hardcoded positions.
- CSV, JSON, XLSX and PDF exports, each self-describing: build, both variant forms, source
  versions, timestamps, the wet-lab protocol, and the disclaimer.
- Progress streaming, docs with verified citations, optional LDlink annotation.
- Landing page, documentation.

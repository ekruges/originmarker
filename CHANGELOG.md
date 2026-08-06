# Changelog

Releases are named after the biology of crossing over. All dates 2026.

The entries below are unusually detailed about bugs. That is deliberate: this tool's output
informs which parental allele an embryo inherited, so a plausible wrong answer is worse than
an error. An error gets investigated; a plausible wrong answer gets used. Anyone deciding
whether to trust a panel from an older build deserves to know exactly what it got wrong.

---

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

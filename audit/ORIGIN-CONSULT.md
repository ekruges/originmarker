# Parental origin as a general feature: what it would take, and where the boundary is

Consult answering the eight questions. Substrate: GSE148488 (Zuccaro 2020), 135 Affymetrix UK Biobank
Axiom arrays x 825,656 markers, GRCh37, 1 sperm donor / 4 egg donors / 31 embryos / 66 units / 15 haploid
pronuclei of known parent. Five independent measurement tracks, each with its own memo and CSVs; this
document is the synthesis and states where the tracks disagreed.

Build verified empirically, not assumed: all 22 autosomal maximum marker positions fit hg19 lengths
(margin 14-248 kb) and 14 of 22 exceed hg38 lengths by 0.1-2.7 Mb. Genotype encoding verified by
cross-tabulation against BAF on arrays spanning the material range, independently in four tracks.

## The short answer

Your instinct to rebuild was right, but not in the direction you expected, and the single most useful
result is a negative one that closes off the most attractive-looking rebuild.

1. Mendelian phase from one parent is free, switch-error-free, and does NOT break the variance floor.
   Three tracks measured this independently and all three agree. The floor is set by noise that lives in
   the parent-allele frame, which is the frame a signed statistic works in, so signing cannot cancel it.
   Measured common-mode variance fraction 0.079-0.208; SE scaling exponents -0.107 to -0.223 on amplified
   material against the -0.500 the hypothesis needs; effective markers per chromosome rise 53 -> 55
   (blastomere) and 102 -> 113 (TE), a factor 1.04-1.11, not the 60-200x hoped for. Your doubt that "more
   markers alone" is the answer was correct, and it is now explained rather than just observed.
2. Signing is still worth doing, for a different reason: it fixes the DIRECTION. The blastomere
   wrong-parent rate falls from 1.00 to 0.25 at f=0.10 and to 0.042 at f=0.30. That is the fix for the
   defect you flagged as most dangerous, and it is not a power gain.
3. The largest available gain remains the second parent, and it is a marker-count gain, not a noise gain.
   Floor ratio 2.8x (f50) to 3.4x (f80) on TE chromosomes, reproducing your ~3x by a different route.
   Design the tool around asking for the second parent.
4. Abandon the mean-of-deviations statistic, not the channel. What replaces it is a per-marker likelihood
   over (BAF, log2R) with self-anchored nuisance parameters and an absolute goodness-of-fit gate. Its
   value is not raw accuracy (0.9455 ungated vs 0.9303 for the discretised rule) but calibrated refusal:
   gated it is 1.0000 on 251/330 known-parent tests where the discrete rule on those same rows is 0.9203.
5. Two things you are not using are worth more than most of what you are. The untransmitted parental
   haplotype (1.40-1.98x SNR, and the only channel that gives a single blastomere a defined floor at all,
   f80 0.628). And CNN-LOH/UPD, which you handle hardly at all, is the EASIEST class to assign a parent
   to, not the hardest.
6. The two-parent version of this problem was solved in clinic during 2023-2025 and your prior survey
   missed it. PH-trace 2025 reaches AUC 0.999 on WGA TE biopsies with an explicit refusal band; AJHG 2023
   reports 99-100% truth-model concordance on 2,277 blastocysts. Adopt for the two-parent case. The
   one-parent case is genuinely unoccupied, and that is where the novel work is.

## 1. What is the best achievable approach

Keep the channels, keep the Mendelian sign, replace the estimator, and add two channels you do not have.

The statistic. At a marker where the loaded parent is homozygous, define the parent-allele share
  PS = BAF if the parent is BB, 1 - BAF if AA
so PS is the fraction of signal attributable to the loaded parent's transmitted allele. The sign comes
from the PARENT's genotype, which is known without reference to the sample, so there is no switch error
and no cohort phasing. Expectations, in the corrected pooled form:

  event                                  E[PS at a truly het marker]   deviation from 0.5
  loss of loaded parent's copy           (1-f)/(2-f)                   -f/(4-2f)
  loss of other parent's copy            1/(2-f)                       +f/(4-2f)
  CNN-LOH, loaded parent's copy lost     (1-f)/2                       -f/2
  CNN-LOH, other parent's copy lost      (1+f)/2                       +f/2
  gain of loaded parent's copy           (1+f)/(2+f)                   +f/(4+2f)
  gain of other parent's copy            1/(2+f)                       -f/(4+2f)

Verified by brute-force copy counting at n=200,000 cells, max |analytic - counted| 1.11e-16. Three
consequences that correct the framing you inherited:

- Rank order by deviation is CNN-LOH > loss > gain at every f. CNN-LOH is the LARGEST-deviation state,
  which inverts the intuition that a class with no copy-number signal must be the hardest.
- The "gains shift about a third as much as a loss" statement is the f=1 endpoint. The pooled ratio is
  (2-f)/(2+f): 0.905 at f=0.1, 0.818 at f=0.2, 0.333 only at f=1. In the mosaic range the gain penalty
  is 1.1-1.7x, not 3x.
- Inversions differ by state and must not be shared: f_loss = 4d/(1+2d), f_gain = 4d/(1-2d),
  f_cnnloh = 2d. At d=0.04 these give 0.148 / 0.174 / 0.080 against 0.240 for the per-cell f=6d form. A
  tool that fixed the loss case by inspection is still wrong on gains, and in the opposite direction.

The estimator. Replace mean-of-deviations with a per-marker likelihood over the joint observation
(BAF, log2R) given (copy-number state x parent of origin). States: CN1_pat, CN1_mat, CN2_norm,
CN2_updpat, CN2_updmat, CN3_gainpat, CN3_gainmat, plus a null-template state (below). Each state implies
q centres with weights obtained by marginalising which alleles the unknown parent transmitted and which
survived amplification. Density Gaussian, REFLECTED at 0 and 1 so no mass is lost off the support, plus a
uniform slab absorbing drop-in and off-cluster probes; separate dispersions at boundary and interior
centres. This dispatches your four known defects structurally: there is no mosaic-fraction inversion
anywhere (copy number is a state, not a displacement to invert), no recentring on a median that sits at
the ceiling, a self-anchored dropout estimator, and an absolute fit test per call.

The nuisance parameters, self-anchored. Among called parent-homozygous markers under CN2, the three band
masses give delta = r/(1+r) with r = a0/ah, and beta = (a0+ah)(1+delta). dropout is identified by a RATIO
of two observed masses on the array itself, with no external allele frequency and no population anchor.
Validated against the mother's own held-back array: median |beta_self - beta_mother| 0.0044 at BAF
dispersion below 0.06, 0.0169 at 0.06-0.12, breaking to 0.0996-0.1121 above 0.20. Negative control on the
father's own genotype, where the answer must be exactly zero, returns 0.00000. This replaces the
d = 1 - h/0.168 estimator, which was EUR-anchored and returned negative dropout at the measured drop-in.
Boundary of validity: sigma = 0.12.

The absolute goodness-of-fit gate, and its honest limit. A parametric chi-square against the mixture
rejects 100% of units including verifiably clean bulk ESC (median chi2/df 84.4): it measures the mixture
model's misfit, not the hypothesis's, and is useless as a guard. Reported so it is not tried again. What
works on the haploid truth set is a bootstrap on the binned multinomial deviance at FIXED subsample size
2,000, which fixes power by design and makes p-values comparable across segments of very different marker
counts. Calibrated: gate at deviance/marker <= 1.20 with dLL >= 20 retains 251/330 tests at 1.0000
accuracy and catches all 18 errors the ungated likelihood makes.
  Limit, stated plainly: a nonparametric shape test with the array as its own reference has real power
only on bulk material (0.793 detection of a misspecified extent at 1% FPR) and is powerless on TE and
blastomere (0.010-0.021), because the null spread on amplified material swamps the shape change. On
single ESC WGA it is worse than powerless, rejecting the CORRECT model 47.7% of the time. On WGA few-cell
material the honest substitute for a fit test is a refusal, not a weaker test. I could not construct one.

Add a null-template state. pat:PNzyg2 chr1 (72.2% no-call, 60.4% het among called, self-referenced log2R
-1.852) is fitted better by a symmetric Beta(1.498, 1.498) than by the best allelic state, by 41,135
log-likelihood units, and 34.4% of its mass sits in zones no allelic state can reach. Your remark that
this "is not a genome" is correct and now has a number. It is what a caller does with intensity and no
template: both allele intensities are background, so BAF lands near 0.5 and the caller either no-calls or
emits a spurious het. The correct output is neither a copy number nor an origin. Fit it as a symmetric
Beta with a near 1.5 and route any segment it wins to refusal. Threshold not calibrated (one chromosome
plus qualitative recurrence in the highest-no-call blastomeres).

Two channels to add:
- The untransmitted parental haplotype. At markers where the loaded parent is HETEROZYGOUS and the other
  homozygous, the child's genotype identifies which parental allele was transmitted, so the untransmitted
  one follows by subtraction. Disjoint from the parent-homozygous set. Yield 27,302-58,571 markers per
  array, 86.2-99.8% of informative markers, with 1.7-13.8% impossible readings giving a direct per-array
  error rate. SNR 1.40-1.98x above the obligate-het channel on every amplified material at every f
  tested. Mechanism is a cleaner marker set (100% truly heterozygous by construction against 32-90% for
  the parent-homozygous window), not drift cancellation.
- BPH vs SPH from parent-HETEROZYGOUS / other-homozygous markers. Budget 42,881-65,044 per genome,
  695-5,426 per chromosome, which is LARGER than the obligate-het origin class. Under BPH the other
  parent's allele share is always 2/3 (one band); under SPH it is 1.0 or 1/3 (two bands, neither at 2/3).
  So the diagnostic is band OCCUPANCY, not a mean. AUC 1.000 from k=400 with 0.0000 SPH contamination.
  Note this reverses the framing in your question: it is not a 3-band resolution problem defeated by
  dispersion 0.21-0.30. Per-marker 3-band assignment does need BAF sd < 0.053 and fails on every WGA
  class (misassignment 0.133-0.218), but a fractional-occupancy statistic over hundreds of markers is
  insensitive to per-marker dispersion. Report the occupancy, never a per-marker classification.

Browser compute budget. Nothing above needs a panel, a cohort, or a server: the sign comes from a
parental genotype, the nuisance parameters from ratios of band masses on the array itself, the thresholds
from the array's own other chromosomes. It is a per-marker pass plus per-window sums over 825k markers, a
few hundred MB of typed arrays. Precedent for the delivery: kana runs single-cell analysis fully in
browser via WebAssembly, ViralWasm reports 2-3x slowdown versus native. Note that both treat the browser
as a delivery convenience; any local execution gives the identical privacy property.

## 2. The boundary of near-certainty

boundary_of_certainty.csv, 22 rows, three verdicts only. Criterion for callable: f80 <= 0.35 AND
wrong-parent < 0.05 at that f AND an absolute goodness-of-fit gate with measured power exists.

  input      material               event class                        verdict            f80
  1 parent   bulk gDNA / ESC        segmental loss >=12 Mb             callable           0.050
  1 parent   bulk gDNA / ESC        whole-chromosome loss              callable           0.050
  1 parent   bulk gDNA / ESC        CNN-LOH / UPD, origin              callable           0.040
  1 parent   bulk gDNA / ESC        CNN-LOH / UPD, detection           withheld           0.679
  1 parent   bulk gDNA / ESC        single-copy gain, origin           callable           0.044
  1 parent   single ESC, WGA        whole-chromosome loss              withheld           0.348
  1 parent   single ESC, WGA        segmental loss >=12 Mb             withheld           none
  1 parent   TE biopsy 5-10         whole-chromosome loss              withheld           0.625
  1 parent   TE biopsy 5-10         CNN-LOH / UPD, origin              callable           0.186
  1 parent   TE biopsy 5-10         single-copy gain, origin           withheld           0.511
  1 parent   single blastomere      whole-chromosome loss              never on this data none
  1 parent   single blastomere      segmental loss >=12 Mb             never on this data none
  1 parent   single blastomere      CNN-LOH / UPD, origin              withheld           0.399
  1 parent   haploid product        which parent contributed           callable           n/a (0.9911)
  2 parents  single ESC, WGA        whole-chromosome loss              callable           0.232
  2 parents  TE biopsy 5-10         whole-chromosome loss              callable           0.186
  2 parents  TE biopsy 5-10         segmental loss >=12 Mb             withheld           0.327
  2 parents  single blastomere      whole-chromosome loss              withheld           0.628
  2 parents  any                    BPH vs SPH given CN3               callable           n/a (AUC 1.000)
  0 parents  any                    any origin call                    NEVER              n/a
  any        >30% deviant autosomes or call rate <0.70                 NEVER              n/a
  any        TE / blastomere, k<=800  copy-number CLASS                NEVER              n/a

What governs it, in order of size:

1. Whether a parent is genotyped at all. This is not a power axis, it is an identifiability axis. Origin
   is defined by reference to a parent, and the nuisance parameter (the probability that the unknown
   parent transmitted the non-parental allele) is confounded with the estimand. Supplied from the other
   parent's own array, the false origin-call rate on euploid chromosomes is 0.0790 (n=937). Supplied from
   a population median of four egg donors, it is 1.0000 (n=475) - every euploid chromosome falsely
   called. A factor of 12.7 on the same statistic, same arrays, same markers. Self-anchoring beta from
   CN2 territory sits between them at 0.2917 (n=24). Note: the whole-chromosome track's memo reports
   0.0024 / n=825 for the mother-array arm; that value is NOT reproducible from the saved
   false_positive_rate.csv, which gives 0.0790 / n=937 under every gate combination I tried. The figures
   here are the ones the artifact supports. The direction and the design conclusion are unchanged and do
   not depend on which is right - a population fallback fails completely.
2. Whether the array retains a self-reference. All 70 whole-chromosome losses in this series sit on
   arrays with 40-100% of autosomes deviant. Forced through, 28 score and 27 of 28 return the same parent
   (27 maternal / 1 paternal) at deviance 1.431-23.640, every one of the 28 above the 1.2 gate - and the
   primary paper predicts the opposite parent (paternal, EYS 6q12). Near-unanimity toward one parent while
   every call's absolute fit test rejects is the signature of a biased statistic, not of biology.
3. Amplification chemistry, which sets a floor offset spanning 23x across materials on the SAME platform
   (tau 0.00037 bulk to 0.00853 blastomere).
4. Event class, through the deviation algebra: CNN-LOH 2.1-3.0x a gain at the same f.
5. Segment length, through the informative-marker count.

## 3. Minimum input

Design the tool around asking for the second parent. It is worth more than four extra cells of the same
embryo (TE f80: 1 -> 5 cells moves 0.186 to 0.135, a factor 1.38; one parent -> two parents on the SAME
single file moves 0.625 to 0.186, a factor 3.36).

  event class                material            minimum input that yields a call     f80
  whole-chrom loss or gain   bulk gDNA / ESC     1 parent, 1 file                     0.050
  whole-chrom loss or gain   single ESC, WGA     2 parents, 1 file                    0.232
  whole-chrom loss or gain   TE biopsy 5-10      2 parents, 1 file                    0.186
  whole-chrom loss or gain   single blastomere   not reachable (2 parents -> 0.628)   none
  segmental >=12 Mb          bulk gDNA / ESC     1 parent, 1 file                     0.050
  segmental >=12 Mb          TE biopsy 5-10      2 parents, >=4 independent cells      <0.20
  segmental >=12 Mb          single ESC / blastomere  not reachable                   none
  CNN-LOH / UPD origin       TE biopsy 5-10      1 parent, 1 file                     0.186
  which parent, haploid      polar body / PN     1 parent, 1 file                     n/a

Mosaic-sensitive calling (f80 <= 0.10) is reachable ONLY on bulk gDNA / ESC material, at both scales,
with one parent. A single blastomere never becomes callable for loss at any k from 1 to 8 (k=8 limit
f80 0.501) on any channel. That is the sharp boundary you asked for: it is a material boundary, not an
input boundary, and no file the user can load moves it.

Second-parent material may be amplified. Reconstructed from >=5 of her own haploid products, an egg donor
is indistinguishable from her bulk gDNA array on this statistic (SNR 13.81 vs 13.94, null sd 0.0105 vs
0.0105). A second parent inferred from diploid SIBLINGS is not acceptable - section 4.

## 4. What the tool can earn from files the user already has

Measured per file structure, as an effective floor reduction. multifile_gain.csv.

  file structure                        what it buys                       ceiling
  technical replicates of one DNA       1.02-1.06x at k=4                  84-95% of null variance shared
                                        (1.03-1.09x at k=inf)
  independent cells of one embryo       1.20x per doubling                 1.60-1.87x
  multiple embryos of one couple        1.00-1.20x on nuisance terms       ICC of drift 0.004-0.300
  untransmitted haplotype (1 parent)    1.40-1.98x SNR                     the only blastomere channel
  haploid products of missing parent    equals a bulk gDNA array at n>=5   -
  diploid siblings as a parent proxy    NOT SHIPPABLE on amplified material -

Four findings worth acting on:

Technical replicate arrays of one WGA product are not additional evidence about copy number. 84-95% of the
null variance is a property of the DNA and its amplification and reproduces on every array run from it. Four
replicate arrays buy a factor 1.020-1.065 in null width against sqrt(4) = 2.0 if they were independent, and
the ceiling at unlimited replicates is only 1.027-1.089.
The dropout already happened before the array. Corollary the tool must not get wrong: agreement between
replicate arrays is NOT evidence that a signal is real.

Independent cells do average down, but the model matters more than the averaging. Fitting shared-event
when the truth is cell-specific inflates implied f by a factor tracking k to within 20% across the grid
(1.9x at k=2, 4.8x at k=5). An f reported without stating which multi-cell model was fitted is wrong by
up to 6x. The two models are also not freely separable at low f, since cell-specific at f and shared at
f/k predict the same statistic.

Within-sibship parent reconstruction is not shippable. Leave-one-out on the same couple's other embryos:
on amplified material the reconstructed parent flips the origin sign on 29-38% of chromosomes and five
siblings reduce that only from 0.376 to 0.292. It does not converge. Raw bias -0.027 BAF against a true
displacement of 0.0072 at f=0.20, i.e. 3.8x the signal it is meant to measure, in the wrong direction.
Diagnosis: at one sibling on the harder donor, 51% of false positives come from markers where the mother
is genuinely HETEROZYGOUS, and more siblings barely touch that class (51% -> 38%), because a sibling
carrying allele B tells you the mother has B, never that she lacks A. It is safe only where the tool
cannot verify the precondition (bulk-quality material both sides, >=3 siblings, dropout below 3%,
assignment margin above 1.5). If offered at all, offer it as a marker-set expansion with its own error
rate displayed, never silently in place of a parent.
  Note this is the mechanism behind the asymmetry you already found: a haploid product of a parent reports
HER alleles directly, a diploid sibling reports a mixture of both. That is why >=5 haploid products give a
usable reference at the same nominal precision that gives a 29-38% sign error from siblings.

Crossover-based phase extension is a bulk-material capability. On clean ESC lines sibling comparison
decisively phases 93% of a chromosome in blocks to 52 Mb with 4 apparent switches, close to the real
meiotic count. On amplified material only 11% is phased, longest block 2.9 Mb, and 22 apparent switches
per chromosome, far above any plausible crossover count - those are dropout artifacts. And 2.9 Mb is
already the BAF autocorrelation length, so it adds no independent blocks.

Sample-parent consistency is a free feature you are not exposing. With one parent genotyped, embryo arrays
read 0.003-0.02 against that parent and 0.04-0.10 against a non-parent, so 23 of 31 embryos assign to an
egg donor with a minimum margin of 67.1 binomial SE. That is a sample-swap and mis-assignment detector
that needs no phasing and no panel.

## 5. Platforms and densities

Two quantities govern the floors, not one. Chemistry sets a floor OFFSET that density cannot cross, and it
also sets how steeply density still buys: the fitted exponent itself ranges -0.46 to -0.18 across materials,
so 'chemistry sets the offset, density sets the slope' is a shorthand and both terms are chemistry-dependent.
Measured by
subsampling 825,656 markers down to 25,000 with chemistry held fixed, then refitting.

  material          power-law exponent   chemistry floor tau   density where tau dominates
  bulk gDNA / ESC        -0.46               0.00037                3.7e6 markers
  single ESC, WGA        -0.34               0.00445                3.6e5
  TE biopsy              -0.30               0.00300                2.6e5
  single blastomere      -0.18               0.00853                8.8e4

This corrects "more markers is not the answer" into a material-specific statement. It holds for TE and
blastomere at your platform's density - both are already past their own saturation points, TE by 3.2x -
and it is wrong for bulk material, still on the sqrt slope at 825,656 markers with a fitted saturation at
3.7 million. Random versus spatially uniform thinning agree to 0.81-1.16x, so the law is a function of
count, not layout.

What transfers: the pooled-dosage algebra; the orientation and self-referencing requirement; intensity
being parent-blind; the sign-reliability threshold in f; the variance decomposition as a fraction.
What must be re-measured per platform: null width in BAF units; f80/f50; the per-array 1% FPR threshold;
effective independent marker count.

The runtime estimator, which is what makes degradation honest. From the user's own file with one parent
loaded, no reference data:
  log sd_null = 1.171 log sd(BAF at het calls) - 0.261 log n_window - 0.354      R2 = 0.930 over 40 cells
sd(BAF at het calls) carries the chemistry axis and is nearly invariant to density (0.0481 vs 0.0479 for
bulk at 825k vs 100k); the window count carries the density axis. Leave-one-chemistry-out, the honest
test for an unseen platform: median ratio 0.75-1.90, worst 2.30, and it errs in BOTH directions -
conservative on bulk, OPTIMISTIC on TE and blastomere, which is the dangerous direction. So: predict the
null width, then multiply the threshold by 2.3 unless that chemistry has been calibrated, and state the
widened floor. The 2.3x rests on four chemistries from one laboratory on one array family and is the
number I would most want widened by external data.

Practical note: HumanKaryomap-12 at ~300,000 markers is what PGT-M laboratories actually run. Measured TE
one-parent f80 at 300k is 0.676 against 0.625 at 826k (1.08x worse) and two-parent 0.191 against 0.186
(1.03x). A Karyomap-density file is within 10% of your dev platform for TE work.

## 6. The classes you handle badly

CNN-LOH and UPD invert the usual difficulty. Copy number stays 2, so the genotype caller does not
degrade: median call rate 0.8674 at copy-neutral log2R against 0.2872 below -1.5 (2,376 progeny plus 330
pronucleus chromosome-arrays). Under maternal UPD
every paternal-homozygous marker yields Mendelian exclusion at full strength. Origin accuracy given
detection is 1.0000 across 2,060 windows at k=5 and every larger window, for every UPD and loss
construct - there is no wrong-parent regime once the deviation clears threshold, only a power constraint.
99% power at k=12 (0.94 Mb) for heterodisomy, k=25 (2.21 Mb) for isodisomy.
  But DETECTING the LOH is 1.75-17.10x harder than assigning it at k=800 (1.28-17.10x across k=100-800),
and the gap is worst where the data is best (17.10x on bulk ESC). Detection cannot use dosage by definition, so it falls to het depletion, which
competes with WGA dropout that depletes heterozygosity identically. Het depletion is a THRESHOLD
statistic: at f=0.5 the allele share is only 0.25 or 0.75, still inside the het band, so depletion is
essentially zero up to f=0.5 even in bulk. Nuisance floor 0.40-0.45 in every material, i.e. one window in
a hundred of a clean array already reads 40% het-depleted from clustered dropout.
  Design consequence: with a parent available, do NOT detect CNN-LOH by het depletion. The allelic
deviation detects the same event at f=0.040-0.469 and should be the trigger; het depletion is a
consequence to report. With NO parent, het depletion is the only channel, the floor is f=0.65-0.90, and
that is a refusal rather than a floor.
  Isodisomy vs heterodisomy: both assignable at full strength; discriminating them needs the
other-parent-het / loaded-parent-hom class (62,798-138,633 markers per donor) rather than the obligate-het
class, which does not separate them. Not measured at power here - the construct sums two independent
amplifications and overstates its own het recovery. MI vs MII origin is refused: the pericentromeric
budget is 2.5-121 markers at +/-1 to +/-10 Mb, chr9 has 0-1 out to +/-5 Mb, and a window wide enough to
hold the markers is wide enough to contain a crossover.

Single-copy gains need 400 informative markers (47.5 Mb median) for 95% power on a complete event, and
f80 0.511-0.691 on TE, 0.620 on blastomere. Mosaic gain origin is out of reach on your target material.
BPH is callable positively at AUC 1.000 from k=400; SPH only by exclusion ("CN3 established by dosage,
and not BPH"), because SPH's bands at 1/3 and 1.0 are both populated by a euploid genome and its power at
1% FPR against euploid is 0.134.

State resolution fails on your material, and this is the most important design consequence in this
section. Power to detect an allelic imbalance exceeds power to resolve which state produced it:

  material          k=400 origin-but-no-class   k=1600   fully classed at k=1600
  bulk gDNA / ESC          0.296                 0.224        0.776
  single ESC, WGA          0.759                 0.643        0.357
  TE biopsy                1.000                 0.899        0.101
  single blastomere        1.000                 0.746        0.254

On bulk material 22-30%, strikingly close to MoChA's 29% unclassified. On your target material 75-100%.
Emit ORIGIN and STATE as separate fields with separate confidence. The log2R channel does not average
down (4x the markers buys 1.2-1.4x, window sd 0.15-0.22 on every WGA class against a needed 0.029-0.081),
so this is not fixable with more markers.
  The real chr9q event in the ESC line demonstrates it: deviation 0.36 on 627-654 markers across four
concordant replicate arrays, maternal, origin beyond any threshold - and every single-state hypothesis
rejected, pure loss at z=124 and pure CNN-LOH at z=10.3, with a pure gain excluded by algebra (f=5.65).
A tool that must emit a class would emit a wrong one here. Note this corrects the prior record, which
labelled it "loss of the paternal chr9": the intensity observation that analysis reported (+0.037, here
+0.072 +/- 0.007 against a matched control set) is itself what rejects a pure loss. Same measurement,
different label.

Breakpoints: origin uncertainty dominates for every class except heterodisomic UPD, where the two are
comparable (0.94 Mb origin vs 0.91 Mb breakpoint interval). For gains it dominates by 50-100x. Breakpoint
error acts on the origin call only through the dilution fraction 2e/k, and power holds at 1.000 until the
mis-assigned fraction reaches 0.10-0.12, so it matters for small segments and is harmless for large ones.

## 7. What the tool should refuse, and how to say it

refusal_taxonomy.csv, 14 conditions, each with the triggering quantity, the measurement behind it, and
the remedy. A refusal that names no remedy is not acceptable output. The load-bearing ones:

  no parent genotyped -> false origin-call rate 0.0790 (n=937) with the second parent's own array vs
    1.0000 (n=475) from a population median. Remedy: load one parental genotype.
  array has no self-reference (>30% deviant autosomes or call rate <0.70) -> remedy: re-amplify and
    re-array, or load a sibling unit of the same embryo, or supply a matched euploid array of the same WGA
    chemistry and accept a stated loss of self-referencing.
  no template (symmetric Beta wins) -> the DNA is absent, not aneuploid. No analytic remedy; re-biopsy.
  goodness of fit rejected (deviance/marker >1.20 or dLL <20) -> outside the model's state space. Report
    unresolved, never as an origin.
  class unresolvable -> report the origin, withhold the class. Needs a 3-5x reduction in window log2R sd,
    which more markers does not deliver.
  implied f >0.5 with own-phase -> switch to an external phase source; above this gate the sample cannot
    phase itself (wrong-parent rises to 0.681 at f=0.70).
  implied f <0.30 on blastomere -> report unassigned; load the second parent or more independent cells.
  unseen chemistry -> multiply the threshold by 2.3 and state the widened floor.

Two reporting disciplines that follow from the measurements rather than from taste:

Report implied f with its bias, or not at all. Median f_hat / f_true is 0.97-0.98 on bulk but 0.34-0.70
on blastomere and 0.55-0.70 on single ESC, worsening with f (blastomere 0.96 at f=0.05 down to 0.34 at
f=0.50). The correct inversion fixes the (2-f) factor; it does not fix dilution-driven under-recovery,
which is calibrated here and not solved.

Frame the operating point as a refusal rate, not as sensitivity. At f=0.10 the phased blastomere statistic
names the wrong parent on 25% of detections instead of 100%, but its power to name the correct parent is
0.0129. The honest description is "it detects almost nothing, and what it detects it gets right." That is
exactly the operating point you asked for, and it must be stated that way.

## 8. Prior art

Full survey in PRIOR_ART.md (27 new records, two-channel retraction check, no retractions, one Author
Correction on the mLOX 2024 record). Three things you should act on.

The two-parent problem is solved in clinic and your prior survey missed it. PH-trace 2025: 1,221
consecutive PGT results, mBAF ratio statistic, ROC critical values >2.408 and <0.380 both at AUC 0.999,
explicit balanced/refuse band, exclusion below 30 informative SNPs per aberration, 52.6% of embryos called
mosaic reclassified euploid, 94.1% of predicted false positives absent on re-biopsy. AJHG 2023: 99-100%
truth-model concordance on 2,277 blastocysts. Both need both parents (PH-trace's statistic is a
maternal:paternal ratio, so with one parent the denominator does not exist). Adopt their thresholds and
refusal discipline for your two-parent path rather than inventing.
  - PH-trace - https://doi.org/10.1093/hropen/hoaf075
  - AJHG 2023 - https://doi.org/10.1016/j.ajhg.2023.03.003

The one-parent architecture is genuinely novel, and both its ingredients are published separately.
hapLOH-AI 2012 is the signed per-marker construction, scored as agreement over a local window against a
known haplotype, with the null measured empirically in known-negative material (0.5005 in a pure normal
sample) rather than assumed. Delaneau 2022 uses duo phase at parent-homozygous sites as GROUND TRUTH
against which panel inference is scored, error assumed ~0, Mendel-inconsistent sites excluded. Nobody has
joined them. Duo phase has no switch-error process, so the mirrored-state transitions MoChA and Numbat
need become unnecessary - a genuine simplification.
  - hapLOH-AI - https://doi.org/10.1101/gr.141374.112
  - duo phase as ground truth - https://doi.org/10.1038/s41467-022-34383-6
  Implication for validation burden: no external benchmark exists for the one-parent case, so every
threshold must be derived and shown in PH-trace's style rather than asserted, and your own 135 arrays are
the only validation available.

The literature predicted the negative result, from its own numbers. hapLOH-AI reports concordance 0.65 at
10% aberrant cells on 103,556 het sites and extrapolates to power >50% at 2-3 aberrant cells per 1,000 on
1.25M signed sites. If those signed observations were independent, that operating point would sit at
z=8.4, not at 50% power; matching the paper's own claim requires SE 5.1x above the independent value, i.e.
~26x variance inflation and ~2,200 effective markers per chromosome. On BULK DNA with panel phase at het
BAF spread 0.03-0.10. Treat as an order-of-magnitude expectation, not a result: the linearity in f is the
paper's extrapolation, "power >50%" is read as exactly 50%, and the 5M figures are themselves
extrapolations. But the direction is the one measured here.

Also worth noting, from the same track: two components are liftable regardless. Numbat's Beta-Binomial
overdispersion is the principled replacement for an assumed binomial when allelic noise is inflated, and
its per-interval phase confidence is the right shape for a gated phase source. SCAN-SNV estimates
amplification imbalance from NEIGHBOURING het markers, which is the closest published instance of the
local-versus-per-probe question, and for MDA sequencing its answer is local - at a 1-10 kb amplicon scale,
a different length scale from your 1.0-2.2 Mb BAF autocorrelation. That mismatch is worth a direct test.
  - Numbat - https://doi.org/10.1038/s41587-022-01468-y
  - SCAN-SNV - https://doi.org/10.1038/s41467-019-11857-8

Gaps that remain in the literature: no duo-based UPD detection method paper (so a one-parent UPD call has
no external benchmark); no published error model for BAF on WGA few-cell material measured on an ARRAY;
no application of a phase-based mosaic caller to embryo or polar-body material, confirming your prior
survey. GENType's single-parent claim resolves against you - its abstract states parents-only haplotyping
requires prior diagnosis of at least one reference embryo by independent technology, so it substitutes a
diagnosed sibling for the missing parent rather than phasing from one parent alone.

## Where the tracks disagreed, and how it resolves

Two tracks reported the phase-selection circularity with opposite sign, and both are right about
different quantities. Comparing the SAME estimator with the marker set frozen before injection versus
re-derived after it, the deviation SHRINKS (blastomere: to 0.320 of truth for UPD, 0.160 for loss, 0.013
for gain), because conditioning on "the sample called this heterozygous" conditions on "the allele share
landed near 0.5". Comparing ACROSS set definitions - parent-defined versus own-het-defined - the apparent
displacement is INFLATED 1.02-1.58x, because a het call selects markers with good template where the
odds shift moves furthest before saturating, while discarding 35% of the window. Both are measured; they
are different contrasts. The operational rule is the same either way and is not sensitive to which
dominates: define the informative set from PARENTAL genotypes; if the sample's own calls must be used,
freeze them outside the candidate interval and never re-derive inside it. Freezing removes the bias
(0.84-1.05 across all materials).

Smaller disagreements, resolved: two tracks quote different one-parent TE floors (0.306 vs 0.625) because
one applies the pooled shift to all parent-homozygous markers and the other only to markers that are
truly heterozygous in the child. The second is correct - 68% of the TE central window is truly
homozygous and carries no dosage signal at all - so the honest TE one-parent floor is the worse number.
Two tracks measured the signed SE exponent on TE as -0.144 and -0.219 with different window
constructions; both are far from -0.500 and the conclusion is unaffected.

## Corrections to the record you gave me

- "Gains shift about a third as much as losses" is the f=1 endpoint. In the mosaic range it is 1.1-1.7x.
- Intensity p=0.54 on haploid products does not reproduce: 0.232 (Mann-Whitney) or 0.107 (Welch) on median
  log2R, and 0.029 on the MEAN, the mean being parent-associated only because it inherits three degraded
  arrays that all happen to be paternal. Conclusion unchanged - use the median, gate on call rate, and the
  channel is parent-blind. Oriented dosage separates at 2.8e-19 here against your 1.3e-15.
- The reported INVERSION below 5 haploid products does not reproduce: 0 inversions in 20-70 leave-one-out
  tests at every n from 1 to 5, with the residual bias POSITIVE (over-reporting the loaded parent's
  presence). The n>=5 gate still holds, but on margin grounds (1.22x at n=1 to 3.30x at n=5), not
  inversion grounds. If the inversion is real in your hands the likely cause is a majority rather than
  unanimity consensus rule, or retention of markers where fewer than n products spoke. That is a code
  question, not a data question, and worth checking before the gate is defended on those grounds.
- "The signal that finds the event destroys the evidence that assigns it" is half wrong as a mechanism.
  At vendor no-call markers the discrete genotype carries 0.0000 bits about origin and continuous BAF
  carries 0.0471 bits, with BAF present at 88-91% of them; the cause is low intensity, not off-cluster.
  Discretisation throws away half the information at called markers and all of it at no-call markers. The
  refusal on whole-chromosome losses still SURVIVES, but because a loss large enough to detect in a single
  cell co-occurs with the array being too damaged to self-reference. In this dataset those two conditions
  are perfectly confounded.
- Effective markers per chromosome: your ~250 and the ~55/113 measured here are different quantities
  (segment-mean scatter on obligate-het markers versus a drift floor resolved from the full density
  ladder). The transferable number is the RATIO between signed and unsigned, 1.04-1.11x.

## Gaps, and what would settle each

- No confirmed UPD or CNN-LOH exists in GSE148488. Verified over 2,376 progeny chromosome-arrays (108 arrays x 22 autosomes): maximum
  paternal absence 0.1365 against an informativity ceiling of 0.104-0.163, and every high-absence
  chromosome is a complete loss. Every UPD number here is from constructions built by recombining real
  allele intensities. Settled by: one array from a karyotype-confirmed UPD case with both parents.
- Whole-chromosome loss and array damage are perfectly confounded here (all 70 losses on arrays with >=40%
  deviant chromosomes). Settled by: a single blastomere carrying ONE whole-chromosome loss on an otherwise
  clean array, e.g. titrated WGA of a known monosomic line at single-cell input. This is the single most
  valuable additional experiment, because it decides whether your 3/3 refusal is intrinsic to the event
  class or incidental to array quality.
- Losses are injected, not observed at known f. The injection is validated in prior project work against
  real trisomies to within 0.007-0.013 in median BAF; the loss direction rests on the pooled algebra.
  Settled by: a real mosaic monosomy of independently known f on this platform.
- One father only. Every loaded-parent result on the paternal side is one individual, so the sign balance,
  consensus call rate and drop-in floor are single-donor quantities.
- ESC_single_WGA rests on 4 arrays from 2 lines and its GOF threshold is demonstrably miscalibrated
  (rejects the correct model 47.7%). Treat every single-ESC number as indicative.
- n is arrays, not embryos, throughout. 16 ESC arrays are ~5 lines x replicates; 4 blastomeres are 2
  embryos. Per-material sd is understated.
- The null-template state rests on one chromosome. Threshold not calibrated.
- The per-marker likelihood is WORSE than the discrete count in the 20-35% no-call band (3/8 correct on 8
  tests). Too few tests to know if that is real. Check before the likelihood replaces the discrete rule as
  a default rather than supplementing it.
- The 2.3x unseen-platform safety factor rests on four chemistries from one laboratory on one array family.
- Zuccaro 2020 full text was not machine-retrievable (bronze OA), so per-embryo truth was not read from
  its supplementary tables. All accuracy numbers rest on the 15 known-parent pronuclei plus internal
  concordance. Reading that table would let every embryo event be scored against external truth.
- X and Y excluded throughout (803,906 of 825,656 markers are autosomal).
- No absolute goodness-of-fit test with adequate power exists for WGA few-cell material. I could not
  construct one; the honest substitute is refusal.

## Ground truth and retraction check

Dataset: GSE148488, platform GPL28377, verified from the GEO record (GSM4472397 series id).
Primary paper: Zuccaro et al. 2020, "Allele-Specific Chromosome Removal after Cas9 Cleavage in Human
Embryos", Cell 183(6). Target locus verified from the GEO sample record, not from memory: source name
reads "semen from subject with homozygous frame shift mutation at rs758109813", which resolves in dbSNP to
chr6, GRCh38 6:63999115, GRCh37 6:64709008, gene EYS, pathogenic frameshift,
NM_001292009.2:c.6794del / p.Pro2265fs. EYS span GRCh37 chr6:64,429,876-66,417,118 (Ensembl REST),
cytoband 6q12 (UCSC hg19), distal to the chr6 centromere at 58.7-63.3 Mb. The affected parental allele is
PATERNAL. Internal observation consistent with the paper's arm-level prediction: chr6 self-referenced
log2R in embryo C blastomeres reads 6p -0.071 / 6q -1.757 and 6p -0.046 / 6q -1.725, i.e. 6q loss sparing
6p.
  Retraction check, method stated: NCBI E-utilities elink and efetch on PMID 33125898, reading every
CommentsCorrections element, plus a targeted ESearch on "retracted publication[pt]", plus Crossref
update-to and updated-by relations on the DOI, plus the OpenAlex is_retracted flag. Result: no retraction,
no erratum, no expression of concern, no corrected-and-republished relation. One linked item of RefType
CommentIn: a peer commentary at Cell 183(6):1464-1466, PMID 33306952. The field has published general
methodological critiques of embryo-editing assessment (Nat Commun 2023, s41467-023-36820-6) which are not
indexed as corrections to this paper. A very recent notice not yet indexed would not appear.
  - GSE148488 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148488
  - Zuccaro 2020 - https://doi.org/10.1016/j.cell.2020.10.025
  - commentary - https://doi.org/10.1016/j.cell.2020.11.022
  - rs758109813 - https://www.ncbi.nlm.nih.gov/snp/rs758109813

## Files

Synthesis: CONSULT_parental_origin_general.md (this document), boundary_of_certainty.csv,
refusal_taxonomy.csv, fig_hypothesis_verdict.png, fig_boundary_and_classes.png.
Track memos and their CSVs: PHASED_FRAME.md, WHOLE_CHROM.md, NEGLECTED_CLASSES.md, MULTIFILE.md,
PRIOR_ART.md.

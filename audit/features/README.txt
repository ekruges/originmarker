BREAKPOINT POSITION AGAINST GENOMIC FEATURES - hg19
====================================================
Asks where called regions fall, not which parent they came from. It therefore needs only
breakpoint POSITION, which the dosage channel produces without any parental genotype, and is the
one part of this work not blocked by the parental-reference problem.

  python3 features.py <regions.tsv>

BUILD. The laboratory arrays are hg19/GRCh37, established rather than assumed: the highest marker
position on chr1 is 249,222,527, which exceeds hg38's chr1 length of 248,956,422 and fits hg19's
249,250,621. All tracks below are hg19.

SOURCES, all public
  cytoBand    UCSC hg19. Used to give coordinates to the common fragile sites, which the
              literature defines by band rather than by base pair.
  refGene     UCSC hg19, 77,614 entries. Used three ways: any gene, genes over 500 kb, and gene
              density. The long-gene class is the one named in this project's correspondence -
              WWOX, DPP10 and the long neuronal genes.
  gap         UCSC hg19. Centromere and telomere gaps, reported because proximity to either is
              both a real breakage correlate and a known artefact source, and must be visible so
              the two are not confused.
  Repli-seq   ENCODE wgEncodeUwRepliSeq wavelet signal, downloaded. NOT yet wired in: it is
              bigWig and needs a parser. Replication timing is the strongest remaining track.

  24 common fragile sites are mapped, the classical aphidicolin-induced set.

THE NULL IS THE POINT
---------------------
A region can only be called where the array carries markers AND where amplification produced
calls. Drawing random intervals uniformly along the genome compares callable genome against
uncallable and reports an enrichment for anything correlated with marker density, which includes
gene density and therefore fragile sites. Every null interval here is drawn on the SAME chromosome
and matched on NUMBER OF INFORMATIVE MARKERS rather than on span, so an interval of 2,000 markers
is compared against other stretches of 2,000 markers.

FIRST RUN IS A NEGATIVE CONTROL, AND IT PASSES
-----------------------------------------------
Run on the 29 regions this audit establishes are amplification dropout rather than copy loss
(see ../regions/FINDINGS.txt, confound 5):

  feature                observed   null    ratio      p
  common fragile site       0.069   0.093    0.74   0.486
  any gene                  1.000   0.999    1.00   0.981
  gene over 500 kb          0.448   0.416    1.08   0.426
  centromeric gap           0.034   0.040    0.87   0.682
  telomeric gap             0.000   0.000    0.00   1.000
  gene density           39.4/Mb  36.6/Mb    1.08

Nothing is enriched, which is the correct answer for regions that are artefact. A feature tool
that reported fragile-site enrichment on those would be broken, so this is the check that the
matched null works rather than a finding about biology.

WHAT IT NEEDS NEXT
------------------
Real breakpoint calls. The dosage channel produces them on all 884 arrays without parents, and
that run has not been done. Until it is, this pipeline has nothing biological to score. Wiring the
Repli-seq bigWig is the other open item, and replication timing is the track most likely to carry
a real signal, with the caveat that Xu 2024 reports the late-replicating fragile compartment is
established symmetrically on both parental genomes, so it is unlikely to explain a PARENTAL
asymmetry even if it explains position.


FIRST RUN ON REAL COPY-NUMBER REGIONS
--------------------------------------
Regions from dosage.py on one experiment, 93 arrays, segmented from log2R in 1-Mb windows against
a panel-of-normals null with each array held out of its own. The deposited copy_number column is
NOT usable for this: 55% of markers read 3.0 and 8% read 2.0 on a normal array, so it is a
per-probe estimate rather than a segmented call.

  93 regions across 21 of the 93 arrays: 26 gains, 67 losses.

Both directions appear, which is the first result in this audit that is not one-directional.

  feature                     all    losses   gains        verdict
  common fragile site        1.09      1.09    1.08     not significant
  gene over 500 kb           1.02      1.04    0.97     not significant
  centromere / telomere      0.75***   0.74***  0.80    DEPLETED, p 0.0005
  late-repl valley (ES)      1.00      1.00    1.00     SATURATED, see below
  late-repl valley (const)   0.92      0.87*   1.11     depleted in losses, p 0.027

NO FRAGILE-SITE ENRICHMENT. The classical aphidicolin-induced sites are touched by 51.6% of
regions against 47.5% for marker-matched intervals, a ratio of 1.09 at p 0.135. On this experiment
the answer is negative, and it is a real answer rather than an absence of one: the same machinery
detects a 0.75 depletion at p 0.0005 in the same run, so it is not simply underpowered.

TWO METHODOLOGICAL FAULTS, BOTH NOW FIXED AND MEASURED.

  1. The ES valley track was SATURATED: 4,843 valleys against regions of 3 Mb and up meant every
     region and every null interval contained one, so the 1.00 carried no information. Rescored as
     valley DENSITY per Mb rather than binary overlap:

       all       1.70/Mb against a null of 1.74   ratio 0.97   p 0.0020
       losses    1.69                    1.74     ratio 0.97   p 0.0075
       gains     1.71                    1.74     ratio 0.98   p 0.0985

     WORTH READING TWICE. That is significant at p 0.002 and it is a three per cent difference.
     With 2,000 permutations and a tight null the p value detects effects far below the size worth
     reporting, so effect size governs here and this one is not a finding.

  2. The centromere depletion is PARTLY THE CALLER. Segmentation required three consecutive
     windows, which is hardest to satisfy near a centromere, and the marker-matched null controls
     for marker density but not for that run-length requirement. Re-run with the requirement
     dropped to one window, which yields 324 regions instead of 93:

       centromere / telomere    0.75 at p 0.0005   ->   0.86 at p 0.031

     The depletion weakens substantially and should be quoted at the weaker figure, if at all.

  THE FRAGILE-SITE NEGATIVE SURVIVES BOTH SETTINGS, which is what makes it a result:

       MIN_WIN 3,  93 regions    1.08   p 0.156
       MIN_WIN 1, 324 regions    0.95   p 0.330

  Different region counts, different caller sensitivity, same answer. On this experiment the
  called copy-number changes are not enriched at common fragile sites.


FULL COHORT: 1,189 REGIONS ACROSS ALL FOUR EXPERIMENTS
-------------------------------------------------------
  DIETER   93 regions across  21 of  93 arrays      26 gains,  67 losses
  JENNA   164 regions across  16 of 178 arrays      14 gains, 150 losses
  ROBLES  417 regions across  65 of 264 arrays     123 gains, 294 losses
  TREFF   515 regions across  73 of 349 arrays     155 gains, 360 losses
  TOTAL 1,189 regions, 175 of 884 arrays carry one, 318 gains and 871 losses

  feature                      all      losses     gains
  common fragile site         0.96*      0.96       0.95
  gene over 500 kb            1.04***    1.04***    1.03
  centromere / telomere       0.81***    0.78***    0.89*
  late-replication valley     0.77***    0.76***    0.79***

THE ANSWER ON FRAGILE SITES IS STILL NO, AND NOW AT n=1,189. The ratio is 0.96, which is a slight
DEPLETION rather than an enrichment, and it does not survive as a finding either: at this sample
size the p is 0.03 for a four per cent difference.

READ THE EFFECT SIZES, NOT THE STARS. Every ratio here lies between 0.77 and 1.04. With 1,189
regions the permutation null is narrow enough that a four per cent difference reaches p 0.0005,
and three of the four features are flagged at that level while none of them moves the rate by more
than a quarter. The largest effect in the table is the late-replication depletion at 0.77, and
"regions avoid the latest-replicating positions" is the opposite of the hypothesis this analysis
was built to test.

The centromere depletion is partly the caller, as recorded above, and the same caution applies to
the late-replication and long-gene figures: none has been re-run at MIN_WIN 1.

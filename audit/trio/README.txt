PARENTAL ORIGIN ON A REAL TRIO WITH EXACT TRUTH - HapMap chr21, CEPH
=====================================================================
No public series carries a confirmed trisomy WITH both parents genotyped, so the origin call had
never been checked against a real answer. This supplies one, on public data.

  ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/forward/
  genotypes_chr21_CEU_r28_nr.b36_fwd.txt.gz     174 samples, 51,010 chr21 markers, CEPH-60-trios

THE TRIO IS CONFIRMED FROM THE GENOTYPES, not taken from a pedigree file. A parent and child share
one allele at every locus, so opposite homozygotes are near zero; two unrelated people are not.

  NA12891 vs NA12878     0.0002 over 33,021 markers      parent-child
  NA12892 vs NA12878     0.0002 over 32,996 markers      parent-child
  NA12891 vs NA12892     0.0677 over 31,797 markers      the couple, unrelated to each other

A 338-fold separation between parent-child and the couple. The pedigree is therefore established
by measurement here rather than asserted.

  2,121 obligate-het markers on chr21, both parents homozygous and opposite.
  The child is heterozygous at 99.53% of them, where Mendel requires 100%. The missing 0.47% is
  genotyping error, and it is the error floor for this panel.

THE TRISOMY IS BUILT FROM THE PARENTS' REAL ALLELES, so which parent contributed the extra
chromosome is exact rather than inferred. An extra copy from one parent makes the allele ratio 2:1
in that parent's favour at every obligate-het marker.

  state                ADO     paternal share     deviation    call        truth
  euploid             0.00           0.5000        +0.0000     unclear     correct
  paternal trisomy    0.00           0.6667        +0.1667     paternal    correct
  maternal trisomy    0.00           0.3333        -0.1667     maternal    correct
  ... repeated at ADO 0.20 and 0.35, nine of nine correct

WHAT THIS DOES AND DOES NOT ESTABLISH. It establishes that the orientation and the decision rule
are correct end to end on real trio genotypes: the euploid child is refused rather than called,
and each direction names the right parent, at the 0.056 margin the module ships.

It does NOT establish noise robustness, and the table shows why: the median over 2,121 markers is
unmoved by symmetric dropout, reading 0.6667 identically at 0%, 20% and 35%. That is a real
property worth knowing - the statistic is robust to symmetric allele dropout by construction - but
it means this test does not stress the caller. The cases that do stress it are elsewhere in this
audit: the marker-count floor, and the constructed control on a real WGA array at 53% call rate.

The remaining untested case is a trisomy that is BIOLOGICALLY real rather than constructed, with
both parents genotyped. That material is not public; see AUDIT-4.1.txt.


STRESS TEST: WHERE THE CALLER BREAKS
-------------------------------------
The table above does not stress anything, so the caller was swept across marker count, allele
dropout, and BIASED dropout, where one allele survives more often than the other, which is the
real platform effect. 200 replicates per cell, calls scored against the known truth.

  markers   ADO  bias   euploid   paternal   maternal   overall
       50  0.30  0.00      1.00       0.98       0.98      0.99
       50  0.30  0.20      1.00       0.96       0.95      0.97
       50  0.50  0.20      0.83       0.80       0.80      0.81
      100  0.30  0.20      1.00       0.99       0.99      0.99
      100  0.50  0.20      0.96       0.92       0.92      0.93
      200  0.50  0.20      1.00       0.99       0.98      0.99
      400  0.50  0.20      1.00       0.98       0.99      0.99
      800  0.50  0.20      1.00       1.00       1.00      1.00

IT BREAKS AT 50 MARKERS UNDER HALF DROPOUT, and not before. At 200 markers it is at 0.99 even with
half the alleles dropping and a 20% asymmetry in which allele survives. The failure at 50 is
symmetric across all three states, so it degrades into refusal rather than into a wrong parent,
which is the safe direction.

  THE SHIPPED FLOOR IS CONSERVATIVE. MIN_INFORMATIVE_DEFAULT is 400 and 200 markers reach 0.99
  under conditions harsher than this material presents. Halving the floor would halve the smallest
  region that can carry a parental label, which is the binding limit on resolution.

  DO NOT ACT ON THAT YET, for a stated reason. These are clean array genotypes with dropout
  simulated as independent per allele. Whole-genome amplification drops markers in spatially
  clustered runs, so neighbouring markers fail together and the effective marker count in a region
  is lower than its nominal one. The floor should be re-derived on WGA material with clustered
  dropout before it is moved, and the number to beat is the 0.99 above.


RE-DERIVING THE FLOOR WITH MEASURED CLUSTERED DROPOUT
------------------------------------------------------
The sweep above simulates dropout independently per allele, which understates the difficulty:
amplification drops markers in runs. Measured on six real WGA arrays, chr1, in marker order:

  no-call rate 0.082 to 0.668, and P(no-call | previous marker no-call) 0.162 to 0.699
  mean LIFT 2.18x, mean run length 1.19 to 3.33 markers

A marker beside a dropped one is 2.18 times more likely to drop. Independent dropout assumes 1.00.
Re-running the sweep with a Markov dropout model at that lift, holding the marginal rate fixed:

  markers   ADO 0.30   ADO 0.50
       50       1.00       0.87
      100       1.00       0.99
      200       1.00       1.00
      400       1.00       1.00

So under a MODEL of clustered dropout the caller needs about 200 informative markers, against the
400 the tool ships.

THE FLOOR HAS NOT BEEN CHANGED, and the reason is a difference in what was measured, not caution
for its own sake. The shipped 400 was derived on REAL WGA blastomeres as the smallest count at
which the worst array in that stage is under 1% two-way error. This is clean array genotypes with
dropout SIMULATED from a two-parameter Markov model fitted to a lift statistic. Those are not the
same measurement, and the simulated one is the weaker evidence: it captures the run-length
structure of dropout but not genotyping error, not mis-clustering, and not whatever else makes a
real blastomere harder than its no-call rate suggests.

  WHAT WOULD JUSTIFY HALVING IT. The same two-way error criterion, computed on real WGA arrays at
  200 markers, against a genotyped parent. The material for that now exists in this project: a
  confirmed maternal reference with biparental children on real arrays. Until that is run the
  floor stays at 400, and the cost of it stays as stated in the module: 400 informative markers
  span a median 47.6 Mb, so a focal gain cannot be annotated and a chromosome arm can.

#!/bin/bash
# Per-chromosome null shift: the real dispersion the self-referenced statistic faces, measured on
# an array with NO event injected. Markers are those the sample called heterozygous with BAF in the
# central window, which is the set the shipped statistic runs on (parent homozygous, child het).
# Orientation is arbitrary here because dispersion, not direction, is what this measures.
f="$1"
gzcat "$f" | awk -v name="$f" '
NR==1 { FS=(index($0,"\t")>0)?"\t":","; $0=$0
  for(i=1;i<=NF;i++){ if($i=="baf")bi=i; else if($i=="genotype")gi=i; else if($i=="chr")ci=i }
  next }
$gi==1 && $bi!="" && $bi>=0.20 && $bi<=0.80 { s[$ci]+=$bi; n[$ci]++; gs+=$bi; gn++ }
END{
  if(gn<1000) exit
  gm=gs/gn
  for(c in s) if(n[c]>=200 && c+0>=1 && c+0<=22) printf "%s,%s,%d,%.6f\n", name, c, n[c], s[c]/n[c]-gm
}'

#!/bin/bash
# Per-array census. The corpus carries TWO layouts: a comma-separated reference export
# (probeset_id,copy_number,chr,position,probe_classification,baf,genotype,normalized_intensity)
# and tab-separated cohort exports (probeset_id,chr,position,log2R,baf,copy_number,genotype).
# Column order differs between them, so columns are located BY NAME from the header rather than
# by position. Assuming one layout silently produced a census of a single array.
f="$1"
gzcat "$f" | awk -v name="$f" '
NR==1 {
  FS = (index($0, "\t") > 0) ? "\t" : ","
  $0 = $0
  for (i=1; i<=NF; i++) {
    if ($i=="baf") bi=i
    else if ($i=="genotype") gi=i
    else if ($i=="log2R" || $i=="normalized_intensity") li=i
  }
  next
}
{
  g[$gi]++
  if ($gi==1 && $bi!="") { s+=$bi; ss+=$bi*$bi; h++ }
  if (li && $li!="") { ls+=$li; lss+=$li*$li; ln++ }
}
END{
  n=g[0]+g[1]+g[2]; nc=g["-1"]+0
  if (n+nc==0) { printf "%s,PARSE_FAILED\n", name; exit }
  sd  = h>1  ? sqrt((ss  - s*s/h)/(h-1))    : -1
  lsd = ln>1 ? sqrt((lss - ls*ls/ln)/(ln-1)): -1
  printf "%s,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f\n", name, g[0],g[1],g[2],nc, n/(n+nc), g[1]/n, sd, lsd
}'

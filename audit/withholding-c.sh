#!/usr/bin/env bash
# Part C of audit/withholding-preregistration.txt: each S23 cell withheld, the father rebuilt by the
# shipped `om reconstruct` from the other seventeen, and the withheld cell scored against that
# rebuild and against the measured sperm donor.
#
# Run from the repository root: audit/withholding-c.sh <dir of S23 .probes> <out dir> [workers]
# with OM_S23_SERIES set to the series number (audit/local/barcode-map.tsv; it stays out of the
# repository). Cells are <series>-01 to -18; <series>-00 is the sperm donor's bulk array.
# One rebuild from seventeen full arrays holds about 6 GB; one worker is the default for that reason.
# A cell whose two scores are already written is skipped, so a stopped run resumes.
set -uo pipefail
S23=$1; OUT=$2; JOBS=${3:-1}
SER=${OM_S23_SERIES:?set OM_S23_SERIES to the series number}
NODE=${NODE:-node}
mkdir -p "$OUT"
one() {
  k=$1
  [ "$(wc -l < "$OUT/score-$k.jsonl" 2>/dev/null || echo 0)" -ge 2 ] && return
  rm -f "$OUT/score-$k.jsonl" "$OUT/score-$k.err"
  others=$(ls "$S23"/"$SER"-*.probes | grep -v -e "$SER-00" -e "$SER-$k")
  nice -n 19 $NODE --experimental-strip-types --max-old-space-size=8192 cli/om.ts reconstruct $others \
    --json --out "$OUT/recon-$k.probes" > "$OUT/recon-$k.json" 2> "$OUT/recon-$k.err"
  for p in "$OUT/recon-$k.probes" "$S23/$SER-00.probes"; do
    nice -n 19 $NODE --experimental-strip-types audit/withholding-score.ts "$S23/$SER-$k.probes" "$p" \
      >> "$OUT/score-$k.jsonl" 2>> "$OUT/score-$k.err"
  done
  rm -f "$OUT/recon-$k.probes"
}
export -f one; export S23 OUT NODE SER
printf '%s\n' 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 | xargs -P "$JOBS" -I{} bash -c 'one {}'
touch "$OUT/C.DONE"

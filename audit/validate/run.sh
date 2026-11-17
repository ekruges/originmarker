#!/usr/bin/env bash
# The validation suite, one command. Run from the repository root on the machine holding the arrays.
#
#   audit/validate/run.sh <out dir> [previous scorecard.json]
#
# Environment:
#   OM_TRIOS     GSE148488 arrays and trio_manifest_full.csv            (required)
#   OM_FATHER    the father reconstructed from his eight pronuclei      (required: the configuration
#                the laboratory runs, so every one-parent check is read against it)
#   OM_S23, OM_S23_SERIES   the laboratory series and its number, for the hold-out rebuild (optional)
#   OM_SURVEY_LIST          a file of every array path, for the no-parent survey (optional)
#   OM_UNITS     the truth units list: with it the constructed gains are also scored with the run's other
#                arrays given, each checked open there first (optional)
#   NODE         node 22 or later                                       (default: node)
#
# Parts run one after another, never together: each holds whole arrays in memory, and the machine
# the corpus lives on hosts other guests.
set -uo pipefail
OUT=$1; PREV=${2:-}
NODE=${NODE:-node}
N="nice -n 19 $NODE --experimental-strip-types"
mkdir -p "$OUT"
shards() { # shards <n> <label> <command...>: run n shards of a harness and wait for all of them
  local n=$1 label=$2; shift 2
  for s in $(seq 0 $((n - 1))); do "$@" "$s" "$n" "$OUT/$label-$s.out" > "$OUT/$label-$s.log" 2>&1 & done
  wait
}
echo "== census: both parents, the measured father alone, the rebuilt father alone"
P=web/public/panel/founders-gse186441.tsv.gz
for arm in two meas1 recon; do
  for s in 0 1 2; do
    case $arm in
      two) F= ;; meas1) F=$OM_TRIOS/GSM4472397.probes.gz ;; recon) F=$OM_FATHER ;;
    esac
    OM_FATHER=$F OM_PANEL_FILE=$P OM_TRIOS=$OM_TRIOS OM_SHARDS=3 OM_SHARD=$s OM_TSV=$OUT/census-$arm-$s.tsv \
      $N audit/refusal-census.ts > "$OUT/census-$arm-$s.log" 2>&1 &
  done
done
wait
echo "== invariance against the rebuilt father"
OM_FATHER=$OM_FATHER OM_TRIOS=$OM_TRIOS shards 3 invariance $N audit/validate/invariance.ts
echo "== constructed gains, the original six chromosomes and the held-out sixteen"
if [ -n "${OM_UNITS:-}" ]; then export OM_RUN=1 OM_UNITS; fi
OM_FATHER=$OM_FATHER OM_TRIOS=$OM_TRIOS shards 3 gains $N audit/validate/constructed-gains.ts
OM_CHROMS=1,3,4,6,7,8,10,11,12,14,15,16,17,19,20,22 OM_FATHER=$OM_FATHER OM_TRIOS=$OM_TRIOS \
  shards 3 gains-heldout $N audit/validate/constructed-gains.ts
if [ -n "${OM_S23:-}" ] && [ -n "${OM_S23_SERIES:-}" ]; then
  echo "== the laboratory series, each cell withheld and its father rebuilt from the rest"
  OM_S23_SERIES=$OM_S23_SERIES NODE=$NODE audit/withholding-c.sh "$OM_S23" "$OUT/holdout" 1
fi
if [ -n "${OM_SURVEY_LIST:-}" ]; then
  echo "== every array, read without a parent"
  for s in 0 1 2; do $N audit/nocall-survey.ts "$OM_SURVEY_LIST" $s 3 "$OUT/survey-$s.tsv" > "$OUT/survey-$s.log" 2>&1 & done
  wait
fi
echo "== scorecard"
python3 audit/validate/scorecard.py \
  --census "$OUT/census-two" "$OUT/census-meas1" "$OUT/census-recon" \
  --invariance "$OUT"/invariance-*.out \
  --constructed "$OUT"/gains-[0-9].out --heldout "$OUT"/gains-heldout-*.out \
  ${PREV:+--previous "$PREV"} --out "$OUT/scorecard.json" | tee "$OUT/report.txt"

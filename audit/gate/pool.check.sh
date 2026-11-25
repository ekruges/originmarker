#!/usr/bin/env bash
# The pool starts every job exactly once: with two schedulers alive at once, and with one killed
# mid-queue and started again. Lanes are served in turn.
#
# Run: bash audit/gate/pool.check.sh
set -uo pipefail
POOL=$(cd "$(dirname "$0")" && pwd)/pool.sh
T=$(mktemp -d)
trap 'kill $(jobs -p) 2>/dev/null; rm -rf "$T"' EXIT
export POOL_DISK_GB=0 POOL_DIR=$T/pool POOL_JOBS=3 POOL_GAP=0 POOL_TICK=1 POOL_FREE_KB=99999999 POOL_RESERVE_GB=1 POOL_JOB_GB=1
fail() { echo "FAIL: $*"; exit 1; }

for lane in a b; do
  : > "$T/q-$lane"
  for i in $(seq 1 6); do printf '%s%s\t%s\techo %s%s >> %s/ran; sleep 1\n' "$lane" "$i" "$T" "$lane" "$i" "$T" >> "$T/q-$lane"; done
  bash "$POOL" submit "$lane" "$T/q-$lane" > /dev/null
done

# Two schedulers; the first is killed after its first jobs start, and a third starts after that.
bash "$POOL" run & s1=$!
sleep 2
kill -9 "$s1"; wait "$s1" 2>/dev/null
POOL_DRAIN=1 bash "$POOL" run & s2=$!
POOL_DRAIN=1 bash "$POOL" run & s3=$!
for _ in $(seq 1 60); do [ -z "$(find "$POOL_DIR/queue" -name '*.job')" ] && [ -z "$(find "$POOL_DIR/running" -name '*.job')" ] && break; sleep 1; done
sleep 2
kill "$s2" "$s3" 2>/dev/null

[ "$(wc -l < "$T/ran")" -eq 12 ] || fail "12 jobs, ran $(wc -l < "$T/ran") lines"
[ "$(sort -u "$T/ran" | wc -l)" -eq 12 ] || fail "a job ran twice: $(sort "$T/ran" | uniq -d | tr '\n' ' ')"
[ "$(find "$POOL_DIR/done" -name '*.exit' | wc -l)" -eq 12 ] || fail "not every job is in done/"
grep -h . "$POOL_DIR"/done/*.exit | grep -qv '^0$' && fail "a job did not exit 0"
# Lanes in turn: the first four starts alternate between the lanes.
order=$(awk '{print substr($4,1,1)}' "$POOL_DIR/pool.log" | head -4 | tr -d '\n')
[ "$order" = abab ] || [ "$order" = baba ] || fail "lanes not in turn: $order"

# A job claimed by a scheduler that died before launching it is given up, never started.
mkdir -p "$POOL_DIR/running"
printf 'x\t%s\techo lost >> %s/ran\n' "$T" "$T" > "$POOL_DIR/running/a__lost.job"
echo $(( $(date +%s) - 100 )) > "$POOL_DIR/running/a__lost.job.claimed"
POOL_DRAIN=1 bash "$POOL" run
grep -q lost "$T/ran" && fail "a claimed job was started again"
[ "$(cat "$POOL_DIR/done/a__lost.job.exit")" = LOST ] || fail "a claimed job not marked LOST"
# Weights: a lane of weight 2 gets two starts for each of a lane of weight 1.
export POOL_DIR=$T/pool2
for lane in w1 w2; do
  : > "$T/q-$lane"; for i in $(seq 1 6); do printf '%s%s\t%s\ttrue\n' "$lane" "$i" "$T" >> "$T/q-$lane"; done
  bash "$POOL" submit "$lane" "$T/q-$lane" > /dev/null
done
bash "$POOL" priority w2 2 > /dev/null
POOL_JOBS=1 POOL_DRAIN=1 bash "$POOL" run
order=$(awk '{print $4}' "$POOL_DIR/pool.log" | head -6 | tr -d '\n')
[ "$(echo "$order" | grep -o w2 | wc -l)" -eq 4 ] || fail "weights not kept: $order"

# Memory: a class whose measured peak does not fit is held; one that fits starts.
export POOL_DIR=$T/pool3
printf 'big\t%s\tsleep 0\nsmall\t%s\tsleep 0\n' "$T" "$T" > "$T/q-m"
bash "$POOL" submit m "$T/q-m" > /dev/null
printf 'm/big\t5000000\n' >> "$POOL_DIR/peaks.tsv"
POOL_FREE_KB=4194304 POOL_RESERVE_GB=1 POOL_JOB_GB=1 timeout 8 bash "$POOL" run
grep -q 'start m small' "$POOL_DIR/pool.log" || fail "a job that fits was not started"
grep -q 'start m big' "$POOL_DIR/pool.log" && fail "a job larger than the memory left was started"
# wait returns once a lane is finished, with each job's exit code, and fails on a failed job.
export POOL_DIR=$T/pool4
printf 'good\t%s\ttrue\nbad\t%s\tfalse\n' "$T" "$T" > "$T/q-x"
bash "$POOL" submit x "$T/q-x" > /dev/null
POOL_DRAIN=1 bash "$POOL" run
out=$(bash "$POOL" wait x); rc=$?
[ "$rc" -eq 1 ] && echo "$out" | grep -q 'good exit 0' && echo "$out" | grep -q 'bad exit 1' || fail "wait: $rc $out"
# A new lane name starts level with the busiest lane; a lane named first goes before the rest.
export POOL_DIR=$T/pool5
printf 'o1\t%s\ttrue\no2\t%s\ttrue\no3\t%s\ttrue\n' "$T" "$T" "$T" > "$T/q-old"
bash "$POOL" submit old "$T/q-old" > /dev/null
POOL_JOBS=1 POOL_DRAIN=1 bash "$POOL" run
printf 'z1\t%s\ttrue\n' "$T" > "$T/q-new"; printf 'p1\t%s\ttrue\np2\t%s\ttrue\n' "$T" "$T" > "$T/q-old2"
bash "$POOL" submit old "$T/q-old2" > /dev/null; bash "$POOL" submit zz "$T/q-new" > /dev/null
POOL_JOBS=1 POOL_DRAIN=1 bash "$POOL" run
order=$(awk '{print $5}' "$POOL_DIR/pool.log" | tail -3 | tr '\n' ' ')
[ "$order" = "p1 z1 p2 " ] || fail "a new lane jumped the queue: $order"
printf 'f1\t%s\ttrue\nf2\t%s\ttrue\n' "$T" "$T" > "$T/q-f"; printf 'q1\t%s\ttrue\n' "$T" > "$T/q-o3"
bash "$POOL" submit old "$T/q-o3" > /dev/null; bash "$POOL" submit late "$T/q-f" > /dev/null
echo 99 > "$POOL_DIR/started/late"; bash "$POOL" first late > /dev/null
POOL_JOBS=1 POOL_DRAIN=1 bash "$POOL" run
[ "$(awk '{print $5}' "$POOL_DIR/pool.log" | tail -3 | tr '\n' ' ')" = "f1 f2 q1 " ] || fail "first lane not served first"
echo "pool checks pass"

#!/usr/bin/env bash
# One job pool for every lane on the scoring machine.
#
#   pool.sh submit <lane> <queue file>   queue lines "<label>\t<dir>\t<command>", the old sched.sh form
#   pool.sh run                          the scheduler (one is enough; a second waits on the lock)
#   pool.sh status                       queued, running and finished jobs, by lane
#
# A job starts only by moving its file out of queue/ into running/, and mv is atomic, so no job
# starts twice, whatever number of schedulers are alive or however often one is killed and started
# again. A job in running/ is never started again: one whose process is gone is moved to done/ with
# its exit code, or as LOST when the scheduler died before its process began.
#
# A job starts while fewer than POOL_JOBS run, while the machine keeps POOL_RESERVE_GB available
# beyond the POOL_JOB_GB the job may take, and no sooner than POOL_GAP seconds after the last start,
# so the reading of memory includes the job before it. The next job is taken from the lanes in
# turn. Jobs run detached (setsid, where there is one) under nice 19 and outlive the scheduler.
set -uo pipefail
P=${POOL_DIR:-/root/om-pool}
JOBS=${POOL_JOBS:-10}
RESERVE_KB=$(( ${POOL_RESERVE_GB:-6} * 1048576 ))
JOB_KB=$(( ${POOL_JOB_GB:-3} * 1048576 ))
GAP=${POOL_GAP:-30}
TICK=${POOL_TICK:-10}
mkdir -p "$P/queue" "$P/running" "$P/done" "$P/logs"
DETACH=$(command -v setsid || true)

avail_kb() { if [ -n "${POOL_FREE_KB:-}" ]; then echo "$POOL_FREE_KB"; else awk '/^MemAvailable:/{print $2}' /proc/meminfo; fi; }
alive() { [ -n "$1" ] && kill -0 "$1" 2>/dev/null; }

# Move finished and lost jobs to done/; print how many still run.
reap() {
  local n=0 j pid
  for j in "$P"/running/*.job; do
    [ -e "$j" ] || continue
    pid=$(cat "$j.pid" 2>/dev/null)
    if [ -z "$pid" ]; then
      # Claimed, never launched: a scheduler died between the two steps. Wait a tick in case it is
      # launching now, then give it up; it is never started from here.
      if [ -e "$j.claimed" ] && [ $(( $(date +%s) - $(cat "$j.claimed") )) -gt $(( TICK * 3 )) ]; then
        echo LOST > "$j.exit"; mv "$j" "$j.exit" "$j.claimed" "$P/done/" 2>/dev/null
      else n=$((n + 1)); fi
    elif alive "$pid"; then n=$((n + 1))
    else
      [ -e "$j.exit" ] || echo LOST > "$j.exit"
      mv "$j" "$j.pid" "$j.exit" "$P/done/" 2>/dev/null; rm -f "$j.claimed"
    fi
  done
  echo "$n"
}

# The next queued job, lanes in turn after the lane that started last.
next_job() {
  local last l f k=0 o i=0 L=()
  last=$(cat "$P/last-lane" 2>/dev/null)
  while read -r l; do L+=("$l"); [ "$l" = "$last" ] && k=$((i + 1)); i=$((i + 1)); done < <(ls "$P/queue")
  for ((o = 0; o < ${#L[@]}; o++)); do
    l=${L[$(((k + o) % ${#L[@]}))]}
    f=$(ls "$P/queue/$l" 2>/dev/null | grep '\.job$' | sort | head -1)
    [ -n "$f" ] && { echo "$l/$f"; return; }
  done
}

start() { # start <lane/file>
  local lane=${1%%/*} name="${1%%/*}__${1#*/}" j
  j="$P/running/$name"
  mv "$P/queue/$1" "$j" 2>/dev/null || return 1   # another scheduler took it
  date +%s > "$j.claimed"
  echo "$lane" > "$P/last-lane"
  local dir cmd label
  IFS=$'\t' read -r label dir cmd < "$j"
  $DETACH bash -c 'cd "$1" && nice -n 19 bash -c "$2" > "$3" 2>&1 < /dev/null; echo $? > "$4"' _ \
    "$dir" "$cmd" "$P/logs/$lane-$label.log" "$j.exit" &
  echo $! > "$j.pid"
  echo "$(date '+%F %T') start $lane $label pid $!" >> "$P/pool.log"
}

case "${1:-}" in
  submit)
    lane=$2; q=$3; mkdir -p "$P/queue/$lane"
    seq=$(date +%s)-$$
    i=0
    while IFS=$'\t' read -r label dir cmd; do
      [ -z "$label" ] && continue
      i=$((i + 1)); f=$(printf '%s-%04d-%s.job' "$seq" "$i" "$label")
      printf '%s\t%s\t%s\n' "$label" "$dir" "$cmd" > "$P/queue/$lane/.$f.tmp" && mv "$P/queue/$lane/.$f.tmp" "$P/queue/$lane/$f"
    done < "$q"
    echo "$i jobs queued for $lane"
    ;;
  run)
    # One scheduler at a time by a lock; the lock is a courtesy, the mv in start() is the guard.
    while ! mkdir "$P/lock" 2>/dev/null; do
      holder=$(cat "$P/lock/pid" 2>/dev/null)
      if [ -n "$holder" ] && ! alive "$holder"; then rm -rf "$P/lock"; else sleep "$TICK"; fi
      [ -n "${POOL_ONCE:-}" ] && exit 0
    done
    echo $$ > "$P/lock/pid"
    trap 'rm -rf "$P/lock"' EXIT
    lastStart=0
    while :; do
      n=$(reap)
      now=$(date +%s)
      if [ "$n" -lt "$JOBS" ] && [ $((now - lastStart)) -ge "$GAP" ] && [ "$(avail_kb)" -ge $((RESERVE_KB + JOB_KB)) ]; then
        j=$(next_job)
        if [ -n "$j" ] && start "$j"; then lastStart=$now; continue; fi
      fi
      if [ -n "${POOL_DRAIN:-}" ] && [ "$n" -eq 0 ] && [ -z "$(next_job)" ]; then exit 0; fi
      sleep "$TICK"
    done
    ;;
  status)
    for d in queue running done; do
      echo "$d: $(find "$P/$d" -name '*.job' | wc -l)"
    done
    for l in $(ls "$P/queue" 2>/dev/null); do echo "  queued $l: $(ls "$P/queue/$l" | grep -c '\.job$')"; done
    for j in "$P"/running/*.job; do [ -e "$j" ] && echo "  running $(basename "$j") pid $(cat "$j.pid" 2>/dev/null)"; done
    for e in "$P"/done/*.exit; do [ -e "$e" ] && echo "  done $(basename "$e" .exit) exit $(cat "$e")"; done | tail -20
    ;;
  *) echo "usage: pool.sh submit <lane> <queue> | run | status" >&2; exit 2 ;;
esac

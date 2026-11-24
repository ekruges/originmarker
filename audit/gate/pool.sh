#!/usr/bin/env bash
# One job pool for every lane on the scoring machine.
#
#   pool.sh submit <lane> <queue file>   queue lines "<label>\t<dir>\t<command>", the old sched.sh form
#   pool.sh priority <lane> <weight>     a lane's share of starts (default 1; 2 starts twice as often)
#   pool.sh run                          the scheduler (one is enough; a second waits on the lock)
#   pool.sh status                       queued, running and finished jobs, by lane
#   pool.sh wait <lane>                  returns when the lane has nothing queued or running; prints
#                                        each finished job's exit code; exits 1 if any is not 0
#
# A job starts only by moving its file out of queue/ into running/, and mv is atomic, so no job
# starts twice, whatever number of schedulers are alive or however often one is killed and started
# again. A job in running/ is never started again: one whose process is gone is moved to done/ with
# its exit code, or as LOST when the scheduler died before its process began.
#
# Cores: a job starts while fewer than POOL_JOBS of the pool's jobs run and the one-minute load is
# under POOL_JOBS, so work started outside the pool counts too. POOL_JOBS defaults to the cores
# less one, left to the host and the other guests.
#
# Memory, measured: the scheduler reads each running job's resident memory (its whole process
# group) every tick and keeps its peak; a finished job's peak is kept for its class (lane and label
# without the trailing shard number). A job may start when available memory, less POOL_RESERVE_GB,
# covers its class's largest peak seen (POOL_JOB_GB where none is), with a tenth more, plus the
# growth still to come of every running job, each counted up to its own class's peak.
#
# Disk: no job starts while the pool's disk has under POOL_DISK_GB free (default 3).
#
# Lanes: the next job is taken from the lane with the fewest starts for its weight. Jobs run
# detached (setsid, where there is one) under nice 19 and outlive the scheduler.
set -uo pipefail
P=${POOL_DIR:-/root/om-pool}
CORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)
JOBS=${POOL_JOBS:-$((CORES - 1))}
RESERVE_KB=$(( ${POOL_RESERVE_GB:-6} * 1048576 ))
JOB_KB=$(( ${POOL_JOB_GB:-3} * 1048576 ))
GAP=${POOL_GAP:-30}
TICK=${POOL_TICK:-10}
mkdir -p "$P/queue" "$P/running" "$P/done" "$P/logs" "$P/prio" "$P/started"
touch "$P/peaks.tsv"
DETACH=$(command -v setsid || true)

avail_kb() { if [ -n "${POOL_FREE_KB:-}" ]; then echo "$POOL_FREE_KB"; else awk '/^MemAvailable:/{print $2}' /proc/meminfo; fi; }
disk_gb() { if [ -n "${POOL_DISK_FREE_GB:-}" ]; then echo "$POOL_DISK_FREE_GB"; else df -Pk "$P" | awk 'NR==2{print int($4 / 1048576)}'; fi; }
load1() { if [ -n "${POOL_LOAD:-}" ]; then echo "$POOL_LOAD"; elif [ -r /proc/loadavg ]; then cut -d' ' -f1 /proc/loadavg; else echo 0; fi; }
alive() { [ -n "$1" ] && kill -0 "$1" 2>/dev/null; }
class_of() { echo "${1##*/}" | sed -E 's/^(.*)__[0-9]+-[0-9]+-[0-9]{4}-(.*)\.job$/\1\/\2/; s/-[0-9]+$//'; }
# Largest peak seen for a class, in KB, with a tenth more; the default where none is.
need_kb() { awk -F'\t' -v c="$1" -v d="$JOB_KB" '$1==c && $2>m{m=$2} END{print (m ? int(m * 1.1) : d)}' "$P/peaks.tsv"; }

# Resident KB of every process group, read once a tick.
rss_by_group() { ps -e -o pgid=,rss= 2>/dev/null | awk '{s[$1]+=$2} END{for (g in s) print g, s[g]}'; }

# Move finished and lost jobs to done/, update peaks; print running count and memory still to come.
reap() {
  local n=0 grow=0 j pid rss peak need groups
  groups=$(rss_by_group)
  for j in "$P"/running/*.job; do
    [ -e "$j" ] || continue
    pid=$(cat "$j.pid" 2>/dev/null)
    if [ -z "$pid" ]; then
      # Claimed, never launched: a scheduler died between the two steps. Wait in case it is
      # launching now, then give it up; it is never started from here.
      if [ -e "$j.claimed" ] && [ $(( $(date +%s) - $(cat "$j.claimed") )) -gt $(( TICK * 3 )) ]; then
        echo LOST > "$j.exit"; mv "$j" "$j.exit" "$j.claimed" "$P/done/" 2>/dev/null
      else n=$((n + 1)); grow=$((grow + JOB_KB)); fi
    elif alive "$pid"; then
      n=$((n + 1))
      rss=$(echo "$groups" | awk -v g="$pid" '$1==g{print $2}'); rss=${rss:-0}
      peak=$(cat "$j.peak" 2>/dev/null || echo 0)
      [ "$rss" -gt "$peak" ] && { peak=$rss; echo "$peak" > "$j.peak"; }
      need=$(need_kb "$(class_of "$j")")
      [ "$need" -gt "$rss" ] && grow=$((grow + need - rss))
    else
      [ -e "$j.exit" ] || echo LOST > "$j.exit"
      [ -e "$j.peak" ] && printf '%s\t%s\n' "$(class_of "$j")" "$(cat "$j.peak")" >> "$P/peaks.tsv"
      mv "$j" "$j.pid" "$j.exit" "$P/done/" 2>/dev/null; rm -f "$j.claimed" "$j.peak"
    fi
  done
  echo "$n $grow"
}

# The next queued job: the lane with the fewest starts for its weight, ties by name, and in it the
# first job, in queue order, whose class fits the memory given (no limit where none is given).
next_job() {
  local l f g best='' bestScore='' s w budget=${1:-}
  for l in $(ls "$P/queue"); do
    f=''
    for g in $(ls "$P/queue/$l" 2>/dev/null | grep '\.job$' | sort); do
      if [ -z "$budget" ] || [ "$(need_kb "$(class_of "${l}__$g")")" -le "$budget" ]; then f=$g; break; fi
    done
    [ -n "$f" ] || continue
    w=$(cat "$P/prio/$l" 2>/dev/null || echo 1)
    s=$(awk -v n="$(cat "$P/started/$l" 2>/dev/null || echo 0)" -v w="$w" 'BEGIN{printf "%012.4f", n / w}')
    if [ -z "$best" ] || [[ "$s" < "$bestScore" ]]; then best="$l/$f"; bestScore=$s; fi
  done
  echo "$best"
}

start() { # start <lane/file>
  local lane=${1%%/*} name="${1%%/*}__${1#*/}" j
  j="$P/running/$name"
  mv "$P/queue/$1" "$j" 2>/dev/null || return 1   # another scheduler took it
  date +%s > "$j.claimed"
  echo $(( $(cat "$P/started/$lane" 2>/dev/null || echo 0) + 1 )) > "$P/started/$lane"
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
  priority)
    echo "$3" > "$P/prio/$2"; echo "$2 weight $3"
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
      read -r n grow < <(reap)
      now=$(date +%s)
      j=$(next_job $(( $(avail_kb) - RESERVE_KB - grow )))
      if [ -n "$j" ] && [ "$n" -lt "$JOBS" ] && [ $((now - lastStart)) -ge "$GAP" ] \
        && awk -v l="$(load1)" -v c="$JOBS" 'BEGIN{exit !(l < c)}' && [ "$(disk_gb)" -ge "${POOL_DISK_GB:-3}" ]; then
        if start "$j"; then lastStart=$now; continue; fi
      fi
      if [ -n "${POOL_DRAIN:-}" ] && [ "$n" -eq 0 ] && [ -z "$(next_job)" ]; then exit 0; fi
      sleep "$TICK"
    done
    ;;
  wait)
    lane=$2
    while [ -n "$(ls "$P/queue/$lane" 2>/dev/null | grep '\.job$')" ] || ls "$P/running/${lane}__"*.job > /dev/null 2>&1; do sleep "$TICK"; done
    bad=0
    for e in "$P"/done/"${lane}"__*.exit; do
      [ -e "$e" ] || continue
      c=$(cat "$e"); echo "$(basename "$e" .job.exit | sed "s/^${lane}__[0-9]*-[0-9]*-[0-9]*-//") exit $c"
      [ "$c" = 0 ] || bad=1
    done
    exit $bad
    ;;
  status)
    for d in queue running done; do
      echo "$d: $(find "$P/$d" -name '*.job' | wc -l)"
    done
    for l in $(ls "$P/queue" 2>/dev/null); do echo "  queued $l: $(ls "$P/queue/$l" | grep -c '\.job$') weight $(cat "$P/prio/$l" 2>/dev/null || echo 1)"; done
    for j in "$P"/running/*.job; do [ -e "$j" ] && echo "  running $(basename "$j") pid $(cat "$j.pid" 2>/dev/null) peak $(cat "$j.peak" 2>/dev/null || echo ?) KB"; done
    for e in "$P"/done/*.exit; do [ -e "$e" ] && echo "  done $(basename "$e" .exit) exit $(cat "$e")"; done | tail -20
    echo "available $(( $(avail_kb) / 1048576 )) GB, load $(load1), at most $JOBS jobs"
    ;;
  *) echo "usage: pool.sh submit <lane> <queue> | priority <lane> <weight> | run | status | wait <lane>" >&2; exit 2 ;;
esac

#!/usr/bin/env bash
# Ship to the NAS and prove it landed. Replaces six manual steps run by hand every release.
#
# THE HEALTH CHECK IS A GATE, NOT A PRINTOUT. A release once went out reporting the previous
# version and it was noticed only because someone read the output. Anything this script verifies,
# it verifies by exiting non-zero.
#
#   usage: scripts/deploy.sh [--skip-checks]
set -euo pipefail
cd "$(dirname "$0")/.."

NAS_DIR=/volume1/docker/originmarker
# The same host and port this script deploys to, reused by the verify step's fallback so an
# environment that cannot resolve the public URL can still confirm what is actually being served.
NAS_SSH=nas
APP_PORT=8091
HOSTS=(https://originmarker.app https://ezrakruger.cc/originmarker)

want=$(grep -oE 'VERSION = "[^"]+"' originmarker/build_info.py | head -1 | sed 's/.*"\(.*\)"/\1/')
echo "deploying $want"

if [ "${1:-}" != "--skip-checks" ]; then
  echo "== release metadata"; ./scripts/release-check.sh
  echo "== self-checks";      ./scripts/checks.sh | tail -2
  echo "== build";            (cd web && npm run build >/dev/null 2>&1) && echo "  ok   vite build"
fi

# A DIRTY TREE IS A STOP, NOT A WARNING. The checks above run against the working tree, so an
# uncommitted file can make them pass while the same checks fail on what was actually committed.
# That is not hypothetical: comparison.ts shipped referring to a field only an uncommitted edit
# added, the suite passed locally on every release for five releases, and CI was red the whole
# time because it checks out the commit. A warning was printed each time and read past.
if [ -n "$(git status --porcelain)" ]; then
  if [ "${ALLOW_DIRTY:-}" = "1" ]; then
    echo "WARNING: working tree is dirty, deploying uncommitted state (ALLOW_DIRTY=1)"
    git status --porcelain | sed 's/^/         /'
  else
    echo "working tree is dirty, so these checks did not test what is committed:"
    git status --porcelain | sed 's/^/  /'
    echo
    echo "commit or stash first, or re-run with ALLOW_DIRTY=1 to deploy this state anyway."
    exit 1
  fi
fi

echo "== copy"
# COPYFILE_DISABLE stops macOS tar smuggling ._* AppleDouble files into the build context, which
# once doubled an apparent file count. data/maps is included: the app reads it from disk.
COPYFILE_DISABLE=1 tar czf - \
  --exclude=.git --exclude=.venv --exclude=node_modules --exclude=tests/fixtures \
  --exclude=.panelbuilder_cache --exclude=__pycache__ --exclude=.pytest_cache \
  --exclude='*.tgz' --exclude=.DS_Store . 2>/dev/null \
  | ssh nas "cat > $NAS_DIR.tgz"

echo "== extract and start"
ssh nas "tar xzf $NAS_DIR.tgz -C $NAS_DIR 2>/dev/null; rm -f $NAS_DIR.tgz; \
  cd $NAS_DIR && docker compose up -d --build" 2>&1 | tail -1

echo "== verify"
fails=0
for host in "${HOSTS[@]}"; do
  got=""
  # Up to thirty tries at 4s. The container is recreated, so the first request can beat it to the
  # port, and the tunnel-routed host takes materially longer to answer than the direct one: 5.0.0
  # verified green on originmarker.app and FAILED on ezrakruger.cc inside a 40s window, then
  # answered correctly by hand moments later. A gate that cries wolf is a gate people stop reading.
  for _ in $(seq 1 30); do
    got=$(curl -s --max-time 20 "$host/api/health" \
      | python3 -c 'import json,sys; print(json.load(sys.stdin).get("version",""))' 2>/dev/null || true)
    [ -n "$got" ] && break
    sleep 4
  done
  if [ "$got" = "$want" ]; then
    printf '  ok   %s reports %s\n' "$host" "$got"
    continue
  fi
  # AN UNREACHABLE HOST IS NOT A WRONG HOST, and conflating them made this gate lie. Outbound DNS
  # is blocked from some environments this script runs in, so the tunnel-routed URL returns nothing
  # while the site is serving perfectly. That produced a FAIL on 5.4.0, 5.5.0 and 5.5.1, all three
  # of which were live and correct, which is how a gate teaches people to ignore it.
  #
  # So an EMPTY answer falls back to asking the host itself over the connection this script already
  # used to deploy. A WRONG answer never falls back: that is a real failure and stays one.
  if [ -z "$got" ]; then
    onhost=$(timeout 30 ssh -o BatchMode=yes -o ConnectTimeout=10 "$NAS_SSH" \
      "curl -s --max-time 8 http://127.0.0.1:${APP_PORT}/api/health" 2>/dev/null \
      | python3 -c 'import json,sys; print(json.load(sys.stdin).get("version",""))' 2>/dev/null || true)
    if [ "$onhost" = "$want" ]; then
      printf '  ok   %s reports %s, verified ON THE HOST (the public URL was unreachable from here)\n' \
        "$host" "$onhost"
      continue
    fi
    printf '  FAIL %s unreachable from here, and the host itself reports "%s"\n' "$host" "$onhost"
    fails=$((fails + 1))
    continue
  fi
  printf '  FAIL %s reports "%s", expected %s\n' "$host" "$got" "$want"
  fails=$((fails + 1))
done

echo
[ "$fails" -eq 0 ] && echo "$want is live on ${#HOSTS[@]} host(s)" || echo "$fails host(s) wrong"
exit "$fails"

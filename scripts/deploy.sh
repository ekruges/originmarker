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
HOSTS=(https://originmarker.app https://ezrakruger.cc/originmarker)

want=$(grep -oE 'VERSION = "[^"]+"' originmarker/build_info.py | head -1 | sed 's/.*"\(.*\)"/\1/')
echo "deploying $want"

if [ "${1:-}" != "--skip-checks" ]; then
  echo "== release metadata"; ./scripts/release-check.sh
  echo "== self-checks";      ./scripts/checks.sh | tail -2
  echo "== build";            (cd web && npm run build >/dev/null 2>&1) && echo "  ok   vite build"
fi

if [ -n "$(git status --porcelain)" ]; then
  echo "WARNING: working tree is dirty, deploying uncommitted state"
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
  else
    printf '  FAIL %s reports "%s", expected %s\n' "$host" "$got" "$want"
    fails=$((fails + 1))
  fi
done

echo
[ "$fails" -eq 0 ] && echo "$want is live on ${#HOSTS[@]} host(s)" || echo "$fails host(s) wrong"
exit "$fails"

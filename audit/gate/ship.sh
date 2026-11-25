#!/usr/bin/env bash
# Copy a commit's tree to the scoring machine with its GATE-KEYS.json, the closures the store keys
# on (audit/gate/plan.ts --keys). The bundled example arrays are left out.
#
#   audit/gate/ship.sh <commit> <dir on tesla>
set -euo pipefail
c=$1; dest=$2
root=$(cd "$(dirname "$0")/../.." && pwd)
keys=$(mktemp)
trap 'rm -f "$keys"' EXIT
node --experimental-strip-types "$root/audit/gate/plan.ts" --keys "$c" 2>/dev/null > "$keys"
ssh -o ControlPath=none tesla "test ! -e '$dest/GATE-KEYS.json'" || { echo "$dest already holds a tree" >&2; exit 1; }
git -C "$root" archive "$c" audit cli originmarker scripts web ':(exclude)web/public/examples' \
  | ssh -o ControlPath=none tesla "mkdir -p '$dest' && tar xf - -C '$dest'"
scp -q -o ControlPath=none "$keys" "tesla:$dest/GATE-KEYS.json"
echo "$(git -C "$root" rev-parse --short "$c") at tesla:$dest"

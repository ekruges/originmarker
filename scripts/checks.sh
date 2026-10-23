#!/usr/bin/env bash
# Every self-check in the repo, each from the directory it needs.
#
# WHY THIS EXISTS RATHER THAN A ONE-LINE LOOP. A `for f in src/*.check.ts` from the repo root
# reports PanelTable and PrimerOptions as FAILING, because those load their component through
# Vite's SSR resolver at `/src/...`, which resolves against the working directory. The failure is
# a false alarm and it costs a re-run to establish that every time. Worse, a runner that cries
# wolf trains you to skim its output, which is the opposite of what a check suite is for.
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
fails=0
run() { # run <label> <dir> <file>
  if (cd "$2" && node --experimental-strip-types "$3" >/dev/null 2>&1); then
    printf '  ok   %s\n' "$1"
  else
    printf '  FAIL %s\n' "$1"
    (cd "$2" && node --experimental-strip-types "$3" 2>&1 | tail -20 | sed 's/^/       /')
    fails=$((fails + 1))
  fi
}

echo "web/src self-checks"
for f in "$ROOT"/web/src/*.check.ts; do
  run "$(basename "$f")" "$ROOT/web" "src/$(basename "$f")"
done

echo "cli self-checks"
for f in "$ROOT"/cli/*.check.ts; do
  run "$(basename "$f")" "$ROOT" "cli/$(basename "$f")"
done

echo "requirements, against the shipped code"
# THE REQUIREMENTS AUDIT, run against the shipped code rather than read off it. Every other check
# here tests a unit; this one exercises the paths a user's run takes and reports what each of the
# original requirements actually produces. It exists because this project has repeatedly shipped
# features that were complete, documented and never wired, each passing its own unit tests because
# each test exercised a fixture the application never sees.
run 'requirements audit' . audit/requirements.check.ts

# TYPESCRIPT IS TYPE-CHECKED HERE, NOT ONLY AT DEPLOY TIME. The self-checks above run through
# tsx, which strips types without checking them, so this script reported "all checks passed" on a
# tree that does not compile. 5.38.0 was tagged and pushed in that state. compileall below has
# covered the Python half since it was written; this is the other half.
echo "typescript type-checks"
if out=$(cd web && npx tsc -b 2>&1); then
  printf '  ok   tsc -b\n'
else
  printf '%s\n' "$out" | sed 's/^/    /'
  printf '  FAIL tsc -b\n'; fails=$((fails + 1))
fi

echo "python modules compile"
if python3 -m compileall -q originmarker app tools >/dev/null 2>&1; then
  printf '  ok   compileall\n'
else
  printf '  FAIL compileall\n'; fails=$((fails + 1))
fi

echo
if [ "$fails" -eq 0 ]; then
  echo "all checks passed"
else
  echo "$fails check(s) failed"
fi
exit "$fails"

#!/usr/bin/env bash
# The version the build reports, the version the changelog documents, and the tag, must agree.
#
# WHY. A release once shipped whose /api/health reported the PREVIOUS version, because the code
# was committed and build_info.py was not. It took a second commit to correct, and it was caught
# only because someone happened to read the deploy output. The version string is what a reader
# uses to work out which build produced a report, so a build that misidentifies itself is worse
# than one that fails to start.
set -uo pipefail
cd "$(dirname "$0")/.."
fails=0
say() { printf '  %s\n' "$1"; }

build=$(grep -oE 'VERSION = "[^"]+"' originmarker/build_info.py | head -1 | sed 's/.*"\(.*\)"/\1/')
codename=$(grep -oE 'CODENAME = "[^"]+"' originmarker/build_info.py | head -1 | sed 's/.*"\(.*\)"/\1/')
# The newest changelog heading, which is the first "## <version>" after the preamble.
changelog=$(grep -m1 -oE '^## [0-9]+\.[0-9]+\.[0-9]+' CHANGELOG.md | sed 's/^## //')

say "build_info.py  $build \"$codename\""
say "CHANGELOG.md   $changelog"

# A CODENAME MUST NEVER BE REUSED. It was, twice: 4.14.0 took "Kinetochore" back from 2.0.0 and
# 4.15.0 took "Recombinase" back from 3.0.0, because the ladder in build_info.py had been spent and
# nothing checked. A codename that names two releases cannot identify either of them.
prior=$(grep -oE '^## [0-9]+\.[0-9]+\.[0-9]+ "[^"]+"' CHANGELOG.md | sed 's/.*"\(.*\)"/\1/' | tail -n +2)
if printf '%s\n' "$prior" | grep -qxF "$codename"; then
  say "REUSED CODENAME: \"$codename\" already names an earlier release. Take the next unused name"
  say "from RELEASES in originmarker/build_info.py."
  fail=1
fi

if [ "$build" != "$changelog" ]; then
  say "MISMATCH: the build reports $build and the newest changelog entry is $changelog"
  fails=$((fails + 1))
fi

# The codename must appear on that changelog heading if the heading carries one at all.
head_line=$(grep -m1 -E '^## [0-9]+\.[0-9]+\.[0-9]+' CHANGELOG.md)
case "$head_line" in
  *'"'*)
    if ! printf '%s' "$head_line" | grep -q "\"$codename\""; then
      say "MISMATCH: build codename \"$codename\" is not the one on $head_line"
      fails=$((fails + 1))
    fi ;;
esac

# Every citation the app ships must be one the verifier actually resolves against Crossref.
missing=0
while read -r id; do
  grep -q "\"$id\"" deploy/verify-citations.py || { say "citation $id is not in verify-citations.py"; missing=$((missing + 1)); }
done < <(grep -oE '^  [a-z0-9_]+: \{' web/src/citations.ts | sed 's/[ :{]//g')
[ "$missing" -gt 0 ] && fails=$((fails + 1))
[ "$missing" -eq 0 ] && say "citations: every shipped id is in the verifier"

if [ "${1:-}" = "--tag" ]; then
  if git rev-parse "v$build" >/dev/null 2>&1; then
    say "tag v$build exists"
  else
    say "NO TAG: v$build has not been created"
    fails=$((fails + 1))
  fi
fi

echo
[ "$fails" -eq 0 ] && echo "release metadata is consistent" || echo "$fails release problem(s)"
exit "$fails"

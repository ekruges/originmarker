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

# THE CODENAME MUST EXIST ON THE LADDER, which the app enforces at RUNTIME and this did not check.
# build_info.current() raises when CODENAME is absent from RELEASES, and /api/health calls it on
# every request, so an unknown name is not a cosmetic slip: it returns 500 for the whole service.
# Two releases shipped that way, 5.8.0 and 5.9.0, and both passed this gate. A check that validates
# how a name is SPELLED while ignoring whether it EXISTS is the wrong half of the problem.
if ! python3 -c "
import sys
sys.path.insert(0, 'originmarker')
import build_info as b
sys.exit(0 if b.CODENAME in {r.name for r in b.RELEASES} else 1)
" 2>/dev/null; then
  say "CODENAME \"$codename\" is not in RELEASES in originmarker/build_info.py."
  say "build_info.current() raises on it, and /api/health calls that on every request, so this"
  say "would ship a service that returns 500. Add the name to the ladder or fix the typo."
  fails=$((fails + 1))
fi

# A CODENAME MUST NOT BE REUSED ACROSS RELEASE LINES. It was, twice: 4.14.0 took "Kinetochore" back
# from 2.0.0 and 4.15.0 took "Recombinase" back from 3.0.0, because the ladder in build_info.py had
# been spent and nothing checked. A codename that names two unrelated releases identifies neither.
#
# A PATCH SHARING ITS MINOR LINE'S NAME IS NOT REUSE. 5.0.1 is part of 5.0.0 "Anaphase", so the
# comparison is against OTHER major.minor lines only. The first version of this guard failed a
# legitimate patch, and set the wrong variable while doing it, so it printed a warning that could
# never fail the build. Both are why this reads the way it does.
line=$(printf '%s' "$build" | cut -d. -f1,2)
clash=$(grep -oE '^## [0-9]+\.[0-9]+\.[0-9]+ "[^"]+"' CHANGELOG.md \
  | awk -v n="$codename" -v l="$line" '{
      v=$2; q=index($0,"\""); name=substr($0,q+1,length($0)-q-1)
      split(v,a,"."); ml=a[1] "." a[2]
      if (name == n && ml != l) print v
    }')
if [ -n "$clash" ]; then
  say "REUSED CODENAME: \"$codename\" already names $(printf '%s' "$clash" | tr '\n' ' ')"
  say "Take the next unused name from RELEASES in originmarker/build_info.py."
  fails=$((fails + 1))
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

# NO LAB DATA IN THE REPOSITORY, AND NO SAMPLE IDENTIFIERS.
#
# The arrays this tool reads are a family's. They get copied into the tree during testing, because
# the browser can only fetch what the dev server serves, and a copy left behind is committed by the
# next `git add -A`. The identifiers are the same problem in a smaller package: they appear in run
# logs, reports and screenshots, and they name real samples.
#
# This is a gate rather than a note because the failure is unrecoverable: a public repository does
# not forget, and rewriting history does not un-publish.
arrays=$(git ls-files | grep -iE '\.(probes|cel)(\.|$)' || true)
if [ -n "$arrays" ]; then
  say "ARRAY DATA IS TRACKED. These are read from disk at run time and never belong in the tree:"
  printf '%s\n' "$arrays" | sed 's/^/    /'
  fails=$((fails + 1))
fi
# Five digits, a dash, two digits, an underscore, two digits, which is the shape these carry.
ids=$(git grep -lE '[0-9]{5}-[0-9]{2}_[0-9]{2}' -- . ':!CHANGELOG.md' 2>/dev/null || true)
if [ -n "$ids" ]; then
  say "A LAB SAMPLE IDENTIFIER APPEARS IN TRACKED FILES:"
  printf '%s\n' "$ids" | sed 's/^/    /'
  fails=$((fails + 1))
fi
[ -z "$arrays" ] && [ -z "$ids" ] && say "no array data and no sample identifiers are tracked"

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

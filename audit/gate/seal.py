"""The open list: which array files of a gate's lists the store and the binary copies may read.

An array is named as audit/truth/build.py names it. One that units.tsv lists may be read only if it
is open; sealed and exposed ones are refused. One that units.tsv does not list may be read only if it
is a laboratory array: a barcode in one of the laboratory folders (build.FOLDER_RANK's first three)
or in GSE288845, the laboratory's own deposit, which is not public. Every published array this
project holds a fact about is in units.tsv; any other file is refused rather than guessed, except a
file named with --derived: one the tool built from open arrays, such as a rebuilt parent.

Writes the allowed paths, one a line, to --out; prints counts only, never a path.

Run: python3 audit/gate/seal.py --units <units.tsv> --out <open.txt> [--derived <file>...] <list.tsv|jobs.jsonl>...
"""
import argparse
import json
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'truth'))
import build  # noqa: E402

LAB = build.FOLDER_RANK[:3] + ['GSE288845']


def paths_of(f):
    """Every array path a run list (tsv) or a known-answer job list (jsonl) names."""
    out = []
    for line in open(f):
        line = line.rstrip('\n')
        if not line:
            continue
        if line.startswith('{'):
            j = json.loads(line)
            out += [j[k] for k in ('path', 'parent', 'other', 'buildParent', 'buildOther') if j.get(k)]
        else:
            c = line.split('\t')
            out += c[1:3] + c[4:5] + (c[6].split(';') if len(c) > 6 else [])
    return [p for p in out if p and p != '-' and build.ARRAY_EXT.search(p)]


def classify(paths, split, derived=()):
    """{path: (array id, 'open' | 'refused')}."""
    with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as t:
        t.write(''.join(f'{p}\t0\n' for p in sorted(set(paths))))
    inv = build.load_inventory(t.name, {})
    os.unlink(t.name)
    out = {}
    for aid, a in inv.items():
        for _, p in a['paths']:
            if aid in split:
                ok = split[aid] == 'open'
            else:
                ok = bool(re.fullmatch(r'\d{5}-\d{2}', aid)) and any(f'/{f}/' in p for f in LAB)
            out[p] = (aid, 'open' if ok or p in derived else 'refused')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--units', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--derived', nargs='*', default=[])
    ap.add_argument('lists', nargs='+')
    a = ap.parse_args()
    split = {}
    for line in open(a.units):
        c = line.rstrip('\n').split('\t')
        if len(c) >= 4 and c[0] != 'array':
            split[c[0]] = c[3]
    paths = [p for f in a.lists for p in paths_of(f)]
    got = classify(paths, split, set(a.derived))
    allowed = sorted(p for p, (_, s) in got.items() if s == 'open')
    with open(a.out, 'w') as f:
        f.write(''.join(f'{p}\n' for p in allowed))
    refused = len(got) - len(allowed)
    print(f'{len(got)} array files named, {len(allowed)} open, {refused} refused')
    return 1 if refused else 0


if __name__ == '__main__':
    sys.exit(main())

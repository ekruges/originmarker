"""The functions a stage's runs called, from node's coverage files (NODE_V8_COVERAGE), as
`path<TAB>offset` lines: each called function's start, in the files of one tree, relative to it.

audit/gate/plan.ts --ran maps the offsets to declarations of the same tree. Prints counts only.

Run: python3 audit/gate/coverage.py <tree root> <out.tsv> <coverage dir>...
"""
import glob
import json
import os
import sys


def main():
    root, out, dirs = os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3:]
    seen = set()
    files = 0
    for d in dirs:
        for f in glob.glob(os.path.join(d, '*.json')):
            files += 1
            for s in json.load(open(f)).get('result', []):
                u = s.get('url', '')
                if not u.startswith(f'file://{root}/'):
                    continue
                rel = u[len(f'file://{root}/'):]
                for fn in s.get('functions', []):
                    r = fn['ranges'][0]
                    if r['count'] > 0:
                        seen.add((rel, r['startOffset'], r['endOffset']))
    with open(out, 'w') as o:
        for rel, a, b in sorted(seen):
            o.write(f'{rel}\t{a}\t{b}\n')
    print(f'{files} coverage files, {len(seen)} called functions')


if __name__ == '__main__':
    main()

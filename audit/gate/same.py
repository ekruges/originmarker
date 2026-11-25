"""Two record sets compared field for field: how many records match, differ, or are missing.

Records are JSON lines keyed by `id` and `tag` (run.ts) or `key`, `kind` and `id` (known-answer-run.ts);
lines that are not JSON are ignored. --skip names top-level fields set aside, such as run times.
Prints counts and, for differing records, the field names that differ, never a sample id.

Run: python3 audit/gate/same.py [--skip f,g] <a.jsonl|dir> <b.jsonl|dir>
"""
import argparse
import collections
import glob
import json
import os
import sys


def records(path):
    files = sorted(glob.glob(os.path.join(path, '*.jsonl'))) if os.path.isdir(path) else [path]
    out = {}
    for f in files:
        for line in open(f):
            if line.startswith('{'):
                r = json.loads(line)
                k = (r.get('id'), r.get('tag'), r.get('key'), r.get('kind'))
                out[k] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--skip', default='')
    ap.add_argument('a')
    ap.add_argument('b')
    a = ap.parse_args()
    skip = set(filter(None, a.skip.split(',')))
    A, B = records(a.a), records(a.b)
    same = differ = 0
    fields = collections.Counter()
    errors = sum(1 for r in B.values() if 'error' in r)
    for k, r in A.items():
        if k not in B:
            continue
        x = {f: v for f, v in r.items() if f not in skip}
        y = {f: v for f, v in B[k].items() if f not in skip}
        if x == y:
            same += 1
        else:
            differ += 1
            for f in set(x) | set(y):
                if x.get(f) != y.get(f):
                    fields[f] += 1
    print(f'{len(A)} and {len(B)} records: {same} identical, {differ} differ, '
          f'{len(set(A) - set(B))} missing from the second, {len(set(B) - set(A))} only in the second, {errors} errors in the second')
    for f, n in fields.most_common():
        print(f'  {f}: {n}')
    return 0 if differ == 0 and len(A) == len(B) == same else 1


if __name__ == '__main__':
    sys.exit(main())

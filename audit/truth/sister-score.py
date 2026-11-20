"""
Score the sister-cell gain reading of plan item 3.1 against the published and sheet tables.

Reads the readings audit/truth/sisters.ts wrote, decides which pairs read as one genome by the
rule web/src/sisterGain.ts ships (parentage.agreement off the chromosomes either cell reports a
change on, above SISTER_AGREE_MIN), applies the call rule, and counts every call against
truth.tsv.

A reading is only taken on a chromosome the cell itself reports no change on, which is what the
module does, so every call counted here is a gain the tool would not otherwise print.

Run: python3 audit/truth/sister-score.py --readings <jsonl> --pairs <pair-chromosome tsv>
       --truth audit/local/truth --runs <jsonl> --tag <tag> [--agree 0.86] [--shift 0.25] [--z 3]
"""
import argparse
import collections
import csv
import json
import os

GAIN = ('gain:whole', 'gain:segmental', 'gain:both',
        'gain+loss:whole', 'gain+loss:segmental', 'gain+loss:both')
LOSS = ('loss:whole', 'loss:segmental', 'loss:both')


def changed_sets(runs, tag):
    """Each cell's changed autosomes in one arm, which is what the comparison leaves out."""
    out = collections.defaultdict(set)
    for line in open(runs):
        if not line.startswith('{'):
            continue
        r = json.loads(line)
        if r.get('tag') != tag or not r.get('id'):
            continue
        for c in r.get('chroms') or []:
            if c.get('aneuploidy') and str(c.get('chrom')).isdigit():
                out[r['id']].add(str(c['chrom']))
        for s in r.get('segments') or []:
            if str(s.get('chrom', '')).isdigit():
                out[r['id']].add(str(s['chrom']))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--readings', required=True)
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--runs', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--agree', type=float, default=0.86)
    ap.add_argument('--shift', type=float, default=0.25)
    ap.add_argument('--z', type=float, default=3.0)
    a = ap.parse_args()

    stated = {}
    for r in csv.DictReader(open(os.path.join(a.truth, 'truth.tsv')), delimiter='\t'):
        if r['kind'] == 'chrom_state':
            stated[(r['array'], r['chrom'])] = r['value']
    units = {r['array']: r['split'] for r in
             csv.DictReader(open(os.path.join(a.truth, 'units.tsv')), delimiter='\t')}
    changed = changed_sets(a.runs, a.tag)

    # Per pair per chromosome, so the aggregate can leave chromosomes out.
    pc = collections.defaultdict(dict)
    for r in csv.DictReader(open(a.pairs), delimiter='\t'):
        pc[tuple(sorted((r['a'], r['b'])))][r['chrom']] = (int(r['shared']), int(r['same']))

    def one_genome(x, y):
        rows = pc.get(tuple(sorted((x, y))))
        if not rows:
            return None
        skip = changed[x] | changed[y]
        shared = same = 0
        for c, (s, m) in rows.items():
            if c in skip:
                continue
            shared += s
            same += m
        if shared <= 1000:          # the floor parentage.agreement refuses below
            return False
        return same / shared > a.agree

    t = collections.Counter()
    gains = set()
    wrong = []
    skipped = collections.Counter()
    for line in open(a.readings):
        if not line.startswith('{'):
            continue
        d = json.loads(line)
        if d.get('changed'):
            continue
        cell, src = d['cell'], d['from']
        if units.get(cell) != 'open' or units.get(src) != 'open':
            skipped['not open'] += 1
            continue
        ok = one_genome(cell, src)
        if ok is None:
            skipped['no pair measured'] += 1
            continue
        if not ok:
            skipped['not one genome'] += 1
            continue
        w = stated.get((cell, d['chrom']))
        truth = ('gain' if w in GAIN else 'loss' if w in LOSS
                 else 'normal' if w == 'normal' else 'unstated')
        kind = 'whole' if d['whole'] else 'piece'
        t['tests'] += 1
        t[f't:{truth}'] += 1
        t[f't:{kind}'] += 1
        if d['shift'] >= a.shift and d['z'] >= a.z:
            t['calls'] += 1
            t[f'c:{truth}'] += 1
            t[f'c:{kind}'] += 1
            if truth == 'gain':
                gains.add((cell, d['chrom'], kind))
            elif truth in ('normal', 'loss'):
                wrong.append((cell, d['chrom'], truth, round(d['shift'], 3), round(d['z'], 2)))

    print(f'arm {a.tag}: readings skipped {dict(skipped)}')
    print(f"  sister readings {t['tests']} (whole {t['t:whole']}, piece {t['t:piece']}); "
          f"stated gain {t['t:gain']}, normal {t['t:normal']}, loss {t['t:loss']}, "
          f"unstated {t['t:unstated']}")
    print(f"  calls {t['calls']}: stated gain {t['c:gain']}, normal {t['c:normal']}, "
          f"loss {t['c:loss']}, unstated {t['c:unstated']}")
    print(f"  distinct stated gains found: {len(gains)} "
          f"({sum(1 for g in gains if g[2] == 'whole')} whole, "
          f"{sum(1 for g in gains if g[2] == 'piece')} piece)")
    if wrong:
        print(f'  WRONG calls: {wrong}')
    return t


if __name__ == '__main__':
    main()

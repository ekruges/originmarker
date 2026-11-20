"""
Score the sister-cell gain reading of plan item 3.1 against the published and sheet tables.

Reads the readings audit/truth/sisters.ts wrote, re-decides which pairs are sisters from a pair
table (audit/truth/sister-pairs.ts), applies the call rule, and counts every call against
truth.tsv. Prints both candidate same-genome measures so the shipped one can be seen beside the
other, and the label the sheets carry so the size of the fault being fixed is visible.

A reading is only taken on a chromosome the cell itself reports no change on, which is what the
module does, so every call counted here is a gain the tool would not otherwise print.

Run: python3 audit/truth/sister-score.py --readings <jsonl> --pairs <tsv> --truth audit/local/truth
       [--embryos <tsv>] [--agree 0.90] [--opp 0.02] [--shift 0.25] [--z 3]
"""
import argparse
import collections
import csv
import json
import os

GAIN = ('gain:whole', 'gain:segmental', 'gain:both',
        'gain+loss:whole', 'gain+loss:segmental', 'gain+loss:both')
LOSS = ('loss:whole', 'loss:segmental', 'loss:both')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--readings', required=True)
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--embryos')
    ap.add_argument('--agree', type=float, default=0.90)
    ap.add_argument('--opp', type=float, default=0.02)
    ap.add_argument('--shift', type=float, default=0.25)
    ap.add_argument('--z', type=float, default=3.0)
    a = ap.parse_args()

    # The stated state of every chromosome of every array, and the seal state.
    stated = {}
    for r in csv.DictReader(open(os.path.join(a.truth, 'truth.tsv')), delimiter='\t'):
        if r['kind'] == 'chrom_state':
            stated[(r['array'], r['chrom'])] = r['value']
    units = {r['array']: r['split'] for r in
             csv.DictReader(open(os.path.join(a.truth, 'units.tsv')), delimiter='\t')}

    pair = {}
    for r in csv.DictReader(open(a.pairs), delimiter='\t'):
        key = tuple(sorted((r['a'], r['b'])))
        pair[key] = (float(r['agree']), float(r['opp_rate']))

    label = {}
    if a.embryos:
        for line in open(a.embryos):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2 and f[0]:
                label[f[0]] = f[1]

    # Every reading, with which rules call its pair a sister.
    tally = collections.defaultdict(collections.Counter)
    missing = collections.Counter()
    calls = collections.defaultdict(list)
    for line in open(a.readings):
        if not line.startswith('{'):
            continue
        d = json.loads(line)
        if d.get('changed'):
            continue
        cell, src = d['cell'], d['from']
        if units.get(cell) != 'open' or units.get(src) != 'open':
            missing['not open'] += 1
            continue
        m = pair.get(tuple(sorted((cell, src))))
        if m is None:
            missing['no pair measured'] += 1
            continue
        agree, opp = m
        called = d['shift'] >= a.shift and d['z'] >= a.z
        want = stated.get((cell, d['chrom']))
        truth = ('gain' if want in GAIN else 'loss' if want in LOSS
                 else 'normal' if want == 'normal' else 'unstated')
        for rule, ok in (('agree', agree > a.agree), ('opp', opp <= a.opp),
                         ('label', label.get(cell) is not None and label.get(cell) == label.get(src))):
            if not ok:
                continue
            tally[rule]['tests'] += 1
            tally[rule][f'test:{truth}'] += 1
            if called:
                tally[rule]['calls'] += 1
                tally[rule][f'call:{truth}'] += 1
                calls[rule].append((cell, d['chrom'], truth, d['shift'], d['z'], d['whole']))

    print(f'readings skipped: {dict(missing)}')
    for rule in ('agree', 'opp', 'label'):
        t = tally[rule]
        print(f"\n-- sister by {rule} --")
        print(f"  pairs tested: {t['tests']} readings "
              f"(stated gain {t['test:gain']}, normal {t['test:normal']}, "
              f"loss {t['test:loss']}, unstated {t['test:unstated']})")
        print(f"  calls: {t['calls']} = gain {t['call:gain']}, normal {t['call:normal']}, "
              f"loss {t['call:loss']}, unstated {t['call:unstated']}")
    return tally, calls


if __name__ == '__main__':
    main()

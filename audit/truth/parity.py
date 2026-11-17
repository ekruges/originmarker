"""One parent against two: the same arrays read with both parents, one measured parent and one
reconstructed parent, each line of the answer key side by side.

Rows are run.ts records tagged `two`, `father` and `recon`; only arrays present under all three
tags are counted, so every column is over the same arrays. Printed calls only, what a reader sees,
but for the one line that says so: a claim computed and withheld is in the record with its reason.

  normal    stated-normal chromosomes with no loss or gain printed on them
  whole     stated whole-chromosome losses printed as a loss
  segment   stated segment losses printed as a loss, by the parent the table names
  parent    stated losses with a stated parent, where the tool names a parent for that chromosome:
            right or wrong, over every stated loss and over the whole-chromosome ones alone

With --scorecard, the validation suite's scorecard.json adds the constructed gains: an extra copy
spliced into a real array, and whether the tool names its parent. The suite reads those with both
parents and with the rebuilt father; it has no measured-father arm.

Run: python3 audit/truth/parity.py --truth audit/local/truth --runs <jsonl>... [--scorecard scorecard.json]
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calls import chrom_calls  # noqa: E402

ARMS = (('two', 'both parents'), ('father', 'measured father'), ('recon', 'rebuilt father'))


def named_parent(rec, chrom, loaded):
    """The parent the record names as missing on a chromosome, from the row or a segment on it.

    A parental-absence segment that carries no parent is the loaded parent's: with one parent
    loaded that is what the page prints."""
    for c in rec['chroms']:
        if c['chrom'] == chrom and c['lostParent']:
            return c['lostParent']
    named = {s.get('absentParent') or loaded for s in rec['segments']
             if s['chrom'] == chrom and s['kind'] == 'parental-absence'} - {None}
    return named.pop() if len(named) == 1 else None


def table(truth_dir, runs):
    states = collections.defaultdict(dict)
    parent = {}
    for r in csv.DictReader(open(os.path.join(truth_dir, 'truth.tsv')), delimiter='\t'):
        if r['kind'] == 'chrom_state':
            states[r['array']][r['chrom']] = r['value']
        elif r['kind'] == 'parent':
            parent[(r['array'], r['chrom'])] = r['value']
    units = {r['array']: r for r in csv.DictReader(open(os.path.join(truth_dir, 'units.tsv')), delimiter='\t')}
    recs = collections.defaultdict(dict)
    for path in runs:
        for line in open(path):
            if line.startswith('{'):
                r = json.loads(line)
                if 'error' not in r and r['tag'] in dict(ARMS):
                    recs[r['tag']][r['id']] = r
    shared = set.intersection(*(set(recs[a]) for a, _ in ARMS)) & set(states)
    if any(units.get(i, {}).get('split') == 'sealed' for i in shared):
        sys.exit('sealed arrays in these runs. Refused.')
    out = {}
    for arm, _ in ARMS:
        c = collections.Counter()
        for i in sorted(shared):
            rec = recs[arm][i]
            calls = chrom_calls(rec, alleles=True, printed_only=True)
            computed = chrom_calls(rec, alleles=True, printed_only=False)
            for chrom, value in states[i].items():
                cls, _, extent = value.partition(':')
                got = calls.get(chrom, frozenset())
                who = parent.get((i, chrom))
                if cls == 'normal':
                    c['normal'] += 1
                    c['normal.left'] += not got
                elif cls == 'loss':
                    key = 'whole' if extent == 'whole' else 'segment' if extent.startswith('seg') else None
                    if not key:
                        continue
                    c[key] += 1
                    c[key + '.hit'] += 'loss' in got
                    if key == 'segment' and who in ('paternal', 'maternal'):
                        c[f'segment.{who}'] += 1
                        c[f'segment.{who}.hit'] += 'loss' in got
                        c[f'segment.{who}.computed'] += 'loss' in computed.get(chrom, frozenset())
                    if who in ('paternal', 'maternal'):
                        said = named_parent(rec, chrom, None if arm == 'two' else 'paternal')
                        c['parent.right'] += said == who
                        c['parent.wrong'] += said is not None and said != who
                        if key == 'whole':
                            c['parent.whole.right'] += said == who
                            c['parent.whole.wrong'] += said is not None and said != who
        out[arm] = c
    return len(shared), out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True)
    ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--scorecard')
    a = ap.parse_args()
    n, t = table(a.truth, a.runs)
    print(f'{n} arrays read in all three arms')
    print(f"{'':44s}" + ''.join(f'{name:>18s}' for _, name in ARMS))
    for label, hit, tot in (
        ('Stated-normal chromosomes left alone', 'normal.left', 'normal'),
        ('Whole-chromosome losses printed', 'whole.hit', 'whole'),
        ('Losses on a segment printed', 'segment.hit', 'segment'),
        ("  of those, the father's", 'segment.paternal.hit', 'segment.paternal'),
        ("  of those, the mother's", 'segment.maternal.hit', 'segment.maternal'),
        ("    the mother's, printed or withheld", 'segment.maternal.computed', 'segment.maternal'),
    ):
        print(f'{label:44s}' + ''.join(f"{t[arm][hit]:>11d} of {t[arm][tot]:<3d}" for arm, _ in ARMS))
    for label, k in (('Parent of a stated loss, right / wrong', 'parent'),
                     ('  whole-chromosome losses only', 'parent.whole')):
        print(f'{label:44s}'
              + ''.join(f"{t[arm][k + '.right']:>12d} / {t[arm][k + '.wrong']:<3d}" for arm, _ in ARMS))
    if a.scorecard:
        card = json.load(open(a.scorecard))
        for label, key in (('Parent of a constructed gain named', 'constructed'),
                           ('  the same, held-out chromosomes', 'constructed_heldout')):
            c = card[key]
            meas = (f"{c['meas1.named']:>11d} of {c['meas1.constructed']:<3d}" if 'meas1.named' in c
                    else f"{'not in the suite':>18s}")
            print(f'{label:44s}' + f"{c['two.named']:>11d} of {c['two.constructed']:<3d}" + meas
                  + f"{c['one.named']:>11d} of {c['one.constructed']:<3d}")
            print(f"{'    wrong':44s}" + f"{c['two.wrong']:>18d}"
                  + (f"{c['meas1.wrong']:>18d}" if 'meas1.wrong' in c else f"{'':>18s}") + f"{c['one.wrong']:>18d}")
            if 'one.lost_parent_on_gain' in c:
                print(f"{'    built gain printed as a lost parent':44s}"
                      + ''.join(f"{c.get(k + '.lost_parent_on_gain', 0):>18d}" for k in ('two', 'meas1', 'one')))

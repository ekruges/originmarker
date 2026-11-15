"""Score run.ts records against the stated facts in audit/local/truth.

Four questions, each against what a published table states about that array:

  material   the stage read with no parent, against the stated material
  ploidy     any copy-number call on a cell, against euploid or aneuploid
  chromosome a loss or a gain on a chromosome, against its per-chromosome karyotype
  break      the nearest copy-segment edge, against a published break interval

A row with no parent, or an unrelated adult standing in for one, is asked what the call rate and the
intensity compute and the parent-free stage; --alleles adds what a true parent's alleles name, for
rows that have one. By default a claim counts whether the tool printed it or withheld it; --printed
counts what a reader sees. `readable` rows are those the parent-free stage does not refuse. Every count is given for
the arrays this tool was developed on (`seen`) and the ones it was not (`fresh`) apart. Sealed arrays
are refused unless --sealed is passed, which is for the one scoring a major version gets.

Run: python3 audit/truth/score.py --truth audit/local/truth --runs <jsonl>... [--tag standin]
       [--alleles] [--printed] [--out scorecard.json]
"""
import argparse
import collections
import csv
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calls import chrom_calls as _chrom_calls, copy_segments as _copy_segments  # noqa: E402

# What each stated material may be read as. A single cell is 'blastomere' or 'single-cell' to this
# tool by its dropout; several cells are 'trophectoderm' or 'bulk'.
CONSISTENT = {
    'blastomere': {'blastomere', 'single-cell'},
    'zygote': {'blastomere', 'single-cell'},
    'pronucleus': {'haploid'},
    'trophectoderm-biopsy': {'trophectoderm', 'bulk'},
    'trophectoderm-outgrowth': {'trophectoderm', 'bulk'},
    'esc-cells': {'trophectoderm', 'bulk'},
    'multi-cell': {'trophectoderm', 'bulk'},
    'bulk-gdna': {'bulk'},
}
NEAR_BP = (1_000_000, 5_000_000)


def load(a):
    truth = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(a.truth, 'truth.tsv')), delimiter='\t'):
        truth[r['array']].append(r)
    units = {r['array']: r for r in csv.DictReader(open(os.path.join(a.truth, 'units.tsv')), delimiter='\t')}
    recs = {}
    for path in a.runs:
        for line in open(path):
            if line.startswith('{'):
                r = json.loads(line)
                if r.get('tag') == a.tag and 'error' not in r:
                    recs[r['id']] = r
    sealed = [i for i in recs if units.get(i, {}).get('split') == 'sealed']
    if sealed and not a.sealed:
        sys.exit(f'{len(sealed)} sealed arrays in these runs. Refused: pass --sealed for the one scoring.')
    burned = set(a.burned.split(',')) | {'GSE148488'}
    arm = {i: ('seen' if units[i]['series'] in burned else 'fresh') for i in recs if i in units}
    return truth, units, recs, arm


def tally():
    return collections.defaultdict(collections.Counter)


def score(a):
    truth, units, recs, arm = load(a)
    chrom_calls = lambda r: _chrom_calls(r, alleles=a.alleles, printed_only=a.printed)  # noqa: E731
    copy_segments = lambda r: _copy_segments(r, printed_only=a.printed)  # noqa: E731
    card = {}
    out = []
    say = out.append
    say(f"{len(recs)} arrays scored: {sum(1 for v in arm.values() if v == 'fresh')} fresh, "
        f"{sum(1 for v in arm.values() if v == 'seen')} seen")
    readable = {i for i, r in recs.items() if r['stageFree']['stage'] != 'failed'}
    for which in ('fresh', 'seen'):
        n = [i for i in recs if arm.get(i) == which]
        say(f"  {which}: {len(n)} arrays, {sum(1 for i in n if i in readable)} readable, "
            f"{sum(1 for i in n if i not in readable)} refused by the stage")
        card[f'{which}.arrays'] = len(n)
        card[f'{which}.readable'] = sum(1 for i in n if i in readable)

    # ---- material
    say('\nMATERIAL: the stage read with no parent, against the stated material')
    table = tally()
    for i, r in recs.items():
        for t in truth[i]:
            if t['kind'] == 'material' and t['source'] in ('jerabek2026', 'turocy2026'):
                m = t['value']
                m = 'blastomere' if m.startswith('blastomere') else m
                table[(arm[i], m)][r['stageFree']['stage']] += 1
    for which in ('fresh', 'seen'):
        ok = tot = fail = 0
        for (w, m), c in sorted(table.items()):
            if w != which:
                continue
            n = sum(c.values())
            good = sum(v for k, v in c.items() if k in CONSISTENT.get(m, set()))
            if m in CONSISTENT:
                ok += good; tot += n - c['failed']; fail += c['failed']
            say(f"  {which:5s} {m:26s} n={n:3d}  " + ', '.join(f'{k} {v}' for k, v in c.most_common()))
        say(f"  {which}: consistent with the stated material {ok} of {tot} read, {fail} refused")
        card[f'{which}.material'] = f'{ok}/{tot}'
        card[f'{which}.material_refused'] = fail

    # ---- sample ploidy
    say('\nPLOIDY: any copy-number call on a cell, against the stated euploid or aneuploid')
    for which in ('fresh', 'seen'):
        for scope, keep in (('all', lambda i: True), ('readable', lambda i: i in readable)):
            c = collections.Counter()
            for i, r in recs.items():
                if arm.get(i) != which or not keep(i):
                    continue
                stated = next((t['value'] for t in truth[i] if t['kind'] == 'ploidy'), None)
                if stated in ('euploid', 'aneuploid'):
                    c[(stated, bool(chrom_calls(r)))] += 1
            an, eu = c[('aneuploid', True)] + c[('aneuploid', False)], c[('euploid', True)] + c[('euploid', False)]
            say(f"  {which:5s} {scope:8s} aneuploid cells with a call {c[('aneuploid', True)]} of {an}; "
                f"euploid cells with no call {c[('euploid', False)]} of {eu}")
            card[f'{which}.{scope}.aneuploid_called'] = f"{c[('aneuploid', True)]}/{an}"
            card[f'{which}.{scope}.euploid_clean'] = f"{c[('euploid', False)]}/{eu}"

    # ---- per chromosome
    say('\nCHROMOSOME: a loss or a gain on a chromosome, against its stated karyotype')
    for which in ('fresh', 'seen'):
        for scope, keep in (('all', lambda i: True), ('readable', lambda i: i in readable)):
            c = collections.Counter()
            for i, r in recs.items():
                if arm.get(i) != which or not keep(i):
                    continue
                calls = chrom_calls(r)
                for t in truth[i]:
                    if t['kind'] != 'chrom_state':
                        continue
                    cls, _, extent = t['value'].partition(':')
                    got = calls.get(t['chrom'], frozenset())
                    if cls == 'normal':
                        c['normal'] += 1
                        c['normal.called'] += bool(got)
                        c['normal.called_gain'] += 'gain' in got
                    elif cls in ('loss', 'gain'):
                        key = f"{cls}.{extent or 'unspecified'}"
                        c[key] += 1
                        c[key + '.hit'] += cls in got
                        c[key + '.opposite'] += (cls not in got) and bool(got)
            if not c:
                continue
            say(f"  {which:5s} {scope:8s} stated normal: {c['normal'] - c['normal.called']} of {c['normal']} left alone "
                f"({c['normal.called']} called, {c['normal.called_gain']} of them a gain)")
            for key in ('loss.whole', 'loss.segmental', 'gain.whole', 'gain.segmental'):
                if c[key]:
                    say(f"  {'':5s} {'':8s} stated {key.replace('.', ', '):16s} found {c[key + '.hit']} of {c[key]}"
                        f"  (called the other way {c[key + '.opposite']})")
                    card[f'{which}.{scope}.{key}'] = f"{c[key + '.hit']}/{c[key]}"
            card[f'{which}.{scope}.normal_left_alone'] = f"{c['normal'] - c['normal.called']}/{c['normal']}"

    # ---- breaks
    say('\nBREAK: the nearest copy-segment edge on that chromosome, against a published break interval')
    for which in ('fresh', 'seen'):
        for source in ('jerabek2026', 'turocy2026'):
            d, none = [], 0
            for i, r in recs.items():
                if arm.get(i) != which or i not in readable:
                    continue
                segs = copy_segments(r)
                for t in truth[i]:
                    if t['kind'] != 'break' or t['source'] != source or not t['start']:
                        continue
                    lo, hi = int(t['start']), int(t['end'])
                    edges = [e for ch, _, s, e2 in segs if ch == t['chrom'] for e in (s, e2)]
                    if not edges:
                        none += 1
                        continue
                    d.append(min(0 if lo <= e <= hi else min(abs(e - lo), abs(e - hi)) for e in edges))
            n = len(d) + none
            if not n:
                continue
            near = [sum(1 for x in d if x <= k) for k in NEAR_BP]
            med = f'{statistics.median(d) / 1e6:.2f} Mb' if d else 'none'
            say(f"  {which:5s} {source:12s} {n} breaks on readable arrays: an edge within 1 Mb {near[0]}, "
                f"within 5 Mb {near[1]}, no segment on that chromosome {none}; median distance where one exists {med}")
            card[f'{which}.{source}.break_within_1Mb'] = f'{near[0]}/{n}'
            card[f'{which}.{source}.break_within_5Mb'] = f'{near[1]}/{n}'
            card[f'{which}.{source}.break_no_segment'] = none
    return '\n'.join(out), card


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True)
    ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--tag', default='standin')
    ap.add_argument('--burned', default='S12,S13,S20,S21,S23,S24')
    ap.add_argument('--sealed', action='store_true')
    ap.add_argument('--alleles', action='store_true')
    ap.add_argument('--printed', action='store_true')
    ap.add_argument('--out')
    args = ap.parse_args()
    text, card = score(args)
    print(text)
    if args.out:
        json.dump(card, open(args.out, 'w'), indent=1)

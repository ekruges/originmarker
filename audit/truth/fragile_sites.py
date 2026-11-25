"""Printed breaks against fragile-site catalogues, with a redraw null matched on chromosome and marker count.

The rule is written in joy-tables/dev/FRAGILE-SITES-plan-item.md (sections 2 to 4) and follows the
owner's decision D9. In short:
  break      a printed start or end of a non-whole gain or loss piece on the result card, not within
             0.5 Mb of the chip's first or last marker on its chromosome, and not the centromere-side
             edge of an arm call (no chip marker between it and the hg19 acen band)
  at a site  start <= printed coordinate <= end
  guide      a break within 1 Mb of a guide site is left out and counted apart
  embryo     breaks of one embryo on one chromosome within 1 Mb of each other are one break; the
             first member in list order (np, lab, h; line; position) stands for it
  null       each break redrawn as the same edge of a stretch of k consecutive chip markers (k = its
             piece's markers) on its chromosome, uniform over stretches whose edge is a break;
             10,000 rounds; one-sided share (1 + rounds >= observed) / (rounds + 1)
Prints counts and positions, never a sample name.

Run: python3 -I audit/truth/fragile_sites.py <records dir> <np list> <lab list> <h list> <chip probes>
       <truth.tsv> <embryo-lab.tsv> <guides-s1.tsv> <cytoband.tsv> <catalogue.tsv>...
With no arguments it runs its own known-answer checks.
"""
import bisect
import collections
import glob
import gzip
import json
import os
import sys

import numpy as np

ROUNDS = 10_000
SEED = 20261125
END_BP = 500_000
GUIDE_BP = 1_000_000
LINK_BP = 1_000_000
CHROMS = [str(i) for i in range(1, 23)] + ['X']
PIECE_KINDS = {'copy-loss', 'copy-gain', 'parental-absence'}


def catalogue(path):
    sites = collections.defaultdict(list)
    for line in open(path):
        if line.startswith(('#', 'name\t')) or not line.strip():
            continue
        name, c, s, e = line.rstrip('\n').split('\t')[:4]
        sites[c].append((int(s), int(e), name))
    return {c: sorted(v) for c, v in sites.items()}


def site_at(sites, chrom, x):
    """The name of the site holding x (ends included), or None."""
    for s, e, name in sites.get(chrom, ()):
        if s <= x <= e:
            return name
    return None


def hit_table(pos, sites, chrom):
    """For each marker index, 1 where the marker lies in a site of this chromosome."""
    h = np.zeros(len(pos), dtype=np.int8)
    for s, e, _ in sites.get(chrom, ()):
        h[np.searchsorted(pos, s, 'left'):np.searchsorted(pos, e, 'right')] = 1
    return h


def allowed_edges(pos, k, side):
    """Marker indices the redrawn edge may take: the same edge of a stretch of k markers that fits on
    the chromosome, the edge not within END_BP of the first or last marker."""
    n = len(pos)
    if k < 1 or k > n:
        return np.empty(0, dtype=np.int64)
    idx = np.arange(0, n - k + 1) if side == 'start' else np.arange(k - 1, n)
    keep = (pos[idx] - pos[0] > END_BP) & (pos[-1] - pos[idx] > END_BP)
    return idx[keep]


def redraw(breaks, pos_of, hits_of, rounds, rng):
    """Total breaks at a site in each round. breaks: (chrom, k, side)."""
    tot = np.zeros(rounds, dtype=np.int64)
    p = []
    for c, k, side in breaks:
        a = allowed_edges(pos_of[c], k, side)
        draw = a[rng.integers(0, len(a), size=rounds)]
        got = hits_of[c][draw]
        tot += got
        p.append(got.mean())
    return tot, p


def by_length(breaks, pos_of, sites):
    """Expected breaks at a site when each is redrawn as the same edge of a stretch of the same span
    in bases, its edge uniform over the stretch places on the chromosome; computed exactly."""
    e = 0.0
    for b in breaks:
        pos = pos_of[b['chrom']]
        lo, hi = pos[0] + END_BP, pos[-1] - END_BP
        if b['side'] == 'start':
            hi = min(hi, pos[-1] - b['span'])
        else:
            lo = max(lo, pos[0] + b['span'])
        if hi <= lo:
            continue
        e += sum(max(0, min(hi, s) - max(lo, a)) for a, s, _ in sites.get(b['chrom'], ())) / (hi - lo)
    return e


def share(tot, observed):
    return (1 + int((tot >= observed).sum())) / (len(tot) + 1)


def chip_markers(path):
    op = gzip.open if path.endswith('.gz') else open
    pos = collections.defaultdict(set)
    with op(path, 'rt') as f:
        head = f.readline().rstrip('\n').split('\t')
        ic, ip = head.index('chr'), head.index('position')
        for line in f:
            r = line.split('\t')
            if r[ic] in CHROMS:
                pos[r[ic]].add(int(r[ip]))
    return {c: np.array(sorted(v), dtype=np.int64) for c, v in pos.items()}


def acen(path):
    out = collections.defaultdict(list)
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        c, s, e, _, stain = line.rstrip('\n').split('\t')[:5]
        if stain == 'acen':
            out[c.replace('chr', '')] += [int(s), int(e)]
    return {c: (min(v), max(v)) for c, v in out.items()}


def is_break(pos, cen, x):
    if x - pos[0] <= END_BP or pos[-1] - x <= END_BP:
        return False, 'chromosome end'
    if cen:
        lo, hi = cen
        if lo <= x <= hi:
            return False, 'arm edge at the centromere'
        # next to the centromere band with no chip marker between: an arm call's unplaced edge
        i = bisect.bisect_right(pos, x)
        if x < lo and (i >= len(pos) or pos[i] >= lo):
            return False, 'arm edge at the centromere'
        j = bisect.bisect_left(pos, x) - 1
        if x > hi and (j < 0 or pos[j] <= hi):
            return False, 'arm edge at the centromere'
    return True, ''


def list_lines(path):
    out = {}
    for n, line in enumerate(open(path), start=1):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4:
            out[(p[0], p[3])] = n
    return out


def main(rec_dir, np_list, lab_list, h_list, chip, truth, emb_lab, guides_s1, cyto, *cats):
    pos_of = chip_markers(chip)
    cen = acen(cyto)
    order = {'np': list_lines(np_list), 'lab': list_lines(lab_list), 'h': list_lines(h_list)}
    rank = {'np': 0, 'lab': 1, 'h': 2}

    embryo, guides = {}, collections.defaultdict(list)
    for line in open(truth):
        p = line.rstrip('\n').split('\t')
        if len(p) < 8:
            continue
        if p[3] == 'embryo':
            embryo[p[0]] = p[7]
        elif p[3] == 'break' and p[1] == 'jerabek2026' and p[5] and p[7] in (
                'on-target', 'proximal to predicted off-target site'):
            guides[p[4]].append((int(p[5]), int(p[6])))
    for line in open(emb_lab):
        a, e = line.rstrip('\n').split('\t')
        embryo[a] = e
    n_s1 = 0
    for line in open(guides_s1):
        c, s, e = line.split('\t')[:3]
        guides[c].append((int(s), int(e)))
        n_s1 += 1
    n_jer = sum(len(v) for v in guides.values()) - n_s1

    # every printed break
    recs = {}
    for f in sorted(glob.glob(os.path.join(rec_dir, '*.jsonl'))):
        if 'crashed' in f:
            continue
        setup = os.path.basename(f).split('-')[0]
        for line in open(f):
            if line.startswith('{'):
                r = json.loads(line)
                recs[(r['id'], r['tag'])] = (setup, r)
    seen = collections.Counter()
    raw = []
    for (rid, tag), (setup, r) in recs.items():
        if 'error' in r or 'chip' in r:
            seen['rows not scored'] += 1
            continue
        line = order[setup][(rid, tag)]
        for g in r['card']['counted']:
            if g['whole']:
                continue
            for p in g['pieces']:
                c = g['chrom']
                if p['kind'] not in PIECE_KINDS or c not in pos_of:
                    seen[f'pieces of kind {p["kind"]} on chr{c} left out'] += 1
                    continue
                pos = pos_of[c]
                k = int(np.searchsorted(pos, p['end'], 'right') - np.searchsorted(pos, p['start'], 'left'))
                for side, x in (('start', p['start']), ('end', p['end'])):
                    ok, why = is_break(pos, cen.get(c), x)
                    if not ok:
                        seen[f'edges not a break: {why}'] += 1
                        continue
                    raw.append(dict(emb=embryo.get(rid, f'array:{rid}'), chrom=c, x=x, side=side, k=k, span=p['end'] - p['start'],
                                    key=(rank[setup], line, x), setup=setup))
    seen['printed edges that are breaks (all cells and setups)'] = len(raw)

    # one break per embryo: link breaks on one chromosome within LINK_BP, first member stands for it
    units = []
    by = collections.defaultdict(list)
    for b in raw:
        by[(b['emb'], b['chrom'])].append(b)
    for v in by.values():
        v.sort(key=lambda b: b['x'])
        group = [v[0]]
        for b in v[1:]:
            if b['x'] - group[-1]['x'] <= LINK_BP:
                group.append(b)
            else:
                units.append(min(group, key=lambda b: b['key']))
                group = [b]
        units.append(min(group, key=lambda b: b['key']))
    units.sort(key=lambda b: (CHROMS.index(b['chrom']), b['x']))

    def near_guide(b):
        return any(s - GUIDE_BP <= b['x'] <= e + GUIDE_BP for s, e in guides.get(b['chrom'], ()))

    at_guide = [b for b in units if near_guide(b)]
    pool = [b for b in units if not near_guide(b)]
    usable = [b for b in pool if len(allowed_edges(pos_of[b['chrom']], b['k'], b['side']))]
    dropped = len(pool) - len(usable)

    P = print
    P('FRAGILE SITES: printed breaks on the open rows against two catalogues (plan 3.4b, decision D9)')
    P('Rule: joy-tables/dev/FRAGILE-SITES-plan-item.md, sections 2 to 4, written before this run.')
    P('Records: 6.27.1, tesla:/root/om-fastgate/S17/rec (np 330, lab 271, h 252 rows). Run once.')
    P(f'Chip markers: {sum(len(v) for v in pos_of.values())} on chr1-22 and X. Seed {SEED}, {ROUNDS} rounds.')
    P(f'Guide sites: {n_s1} of the lab\'s Table S1 (hg19) and {n_jer} published on-target or off-target '
      f'break intervals of Jerabek 2026; a break within {GUIDE_BP / 1e6:g} Mb is left out.\n')
    P('Counting')
    for k, v in sorted(seen.items()):
        P(f'  {k}: {v}')
    P(f'  breaks, one per embryo: {len(units)} (embryos with a break: {len({b["emb"] for b in units})})')
    P(f'  left out at a guide site: {len(at_guide)}')
    P(f'  dropped, no stretch to redraw on: {dropped}')
    P(f'  breaks in the test: {len(usable)} '
      f'({collections.Counter(b["setup"] for b in usable)["np"]} first seen on np rows, '
      f'{collections.Counter(b["setup"] for b in usable)["lab"]} lab, '
      f'{collections.Counter(b["setup"] for b in usable)["h"]} h)')
    P(f'  by side: {dict(collections.Counter(b["side"] for b in usable))}; '
      f'marker count per piece median {int(np.median([b["k"] for b in usable])) if usable else 0}')

    rng = np.random.default_rng(SEED)
    for path in cats:
        sites = catalogue(path)
        n_sites = sum(len(v) for v in sites.values())
        hits_of = {c: hit_table(pos_of[c], sites, c) for c in pos_of}
        obs = [site_at(sites, b['chrom'], b['x']) for b in usable]
        o = sum(1 for s in obs if s)
        tot, p = redraw([(b['chrom'], b['k'], b['side']) for b in usable], pos_of, hits_of, ROUNDS, rng)
        e = float(tot.mean())
        P(f'\n== {os.path.basename(path)}: {n_sites} sites ==')
        P(f'  observed {o} of {len(usable)} breaks at a site; expected {e:.2f}; '
          f'observed over expected {o / e if e else float("nan"):.2f}; one-sided share {share(tot, o):.4f}')
        P(f'  chance rounds: 5th to 95th percentile {np.percentile(tot, 5):.0f} to {np.percentile(tot, 95):.0f}')
        per = collections.defaultdict(lambda: [0, 0.0])
        for s in obs:
            if s:
                per[s][0] += 1
        # expected per site: the redraw chance of each break landing in that site
        for b in usable:
            c = b['chrom']
            a = allowed_edges(pos_of[c], b['k'], b['side'])
            for s0, e0, name in sites.get(c, ()):
                inside = ((pos_of[c][a] >= s0) & (pos_of[c][a] <= e0)).mean()
                if inside:
                    per[name][1] += float(inside)
        P('  per site, observed and expected (sites with either above zero):')
        for name in sorted(per, key=lambda n: (-per[n][0], -per[n][1])):
            if per[name][0] or per[name][1] >= 0.005:
                P(f'    {name:<10} observed {per[name][0]:>3}  expected {per[name][1]:6.2f}')
        P('  breaks at a site (chromosome, printed coordinate, site):')
        for b, s in zip(usable, obs):
            if s:
                P(f'    chr{b["chrom"]} {b["x"]:>11,}  {s}')
    P('\n== A comparison, not the rule: the expectation by length ==')
    P('D9 redraws by marker count, and the result above is the result. Beside it, each break redrawn as')
    P('the same edge of a stretch of the same length in bases, uniform over every place on its chromosome')
    P(f'where that edge is not within {END_BP / 1e6:g} Mb of the first or last marker:')
    for path in cats:
        sites = catalogue(path)
        P(f'  {os.path.basename(path)}: expected by length {by_length(usable, pos_of, sites):.2f}')
    c6, fr = pos_of['6'], [x for x in catalogue(cats[0]).get('6', ()) if x[2] == 'FRA6H']
    if fr:
        s0, e0, _ = fr[0]
        m = ((c6 >= s0) & (c6 <= e0)).mean()
        L = (e0 - s0) / (c6[-1] - c6[0])
        P(f'  FRA6H holds {m:.1%} of chromosome 6\'s chip markers in {L:.1%} of its length '
          f'({m / L:.1f} times the density of the chromosome as a whole).')
    per_emb = collections.Counter(b['emb'] for b in usable)
    P(f'  Breaks cluster by embryo: the five embryos with most breaks carry '
      f'{sum(n for _, n in per_emb.most_common(5))} of {len(usable)}.')
    P('\nBreaks left out at a guide site (chromosome, printed coordinate):')
    for b in at_guide:
        P(f'    chr{b["chrom"]} {b["x"]:>11,}')


def demo():
    """Known answers on made-up data."""
    sites = {'1': [(100, 200, 'A')]}
    assert site_at(sites, '1', 150) == 'A'
    assert site_at(sites, '1', 100) == 'A' and site_at(sites, '1', 200) == 'A'  # ends included
    assert site_at(sites, '1', 99) is None and site_at(sites, '1', 201) is None
    assert site_at(sites, '2', 150) is None

    # A chromosome of 101 markers, one per Mb from 0 to 100 Mb; a site from 40 to 49 Mb (10 markers).
    # Start edges of 10-marker stretches: index 0..91; not within 0.5 Mb of the ends: 1..91 (91 places).
    # 10 of them (40..49) are in the site, so the chance is 10/91 per break.
    pos = np.arange(0, 101, dtype=np.int64) * 1_000_000
    hits = {'1': hit_table(pos, {'1': [(40_000_000, 49_000_000, 'S')]}, '1')}
    a = allowed_edges(pos, 10, 'start')
    assert a[0] == 1 and a[-1] == 91 and len(a) == 91
    b = allowed_edges(pos, 10, 'end')
    assert b[0] == 9 and b[-1] == 99 and len(b) == 91
    rng = np.random.default_rng(1)
    tot, p = redraw([('1', 10, 'start')] * 20, {'1': pos}, hits, 20000, rng)
    want = 20 * 10 / 91
    assert abs(tot.mean() - want) < 0.05, (tot.mean(), want)
    # one-sided share: with every round at 0, observing 1 gives 1/(n+1); observing 0 gives 1
    z = np.zeros(99, dtype=np.int64)
    assert share(z, 1) == 1 / 100 and share(z, 0) == 1.0
    # chromosome ends and the arm edge next to the centromere band
    cen = (55_500_000, 56_500_000)
    assert is_break(pos, cen, 400_000)[0] is False
    assert is_break(pos, cen, 100_000_000)[0] is False
    assert is_break(pos, cen, 55_000_000) == (False, 'arm edge at the centromere')
    assert is_break(pos, cen, 57_000_000) == (False, 'arm edge at the centromere')
    assert is_break(pos, cen, 54_000_000) == (True, '')
    assert is_break(pos, cen, 56_000_000) == (False, 'arm edge at the centromere')  # inside the band
    # by length: one start edge of a 10 Mb stretch on a 0 to 100 Mb chromosome, edge uniform over
    # 0.5 to 90 Mb; a site of 40 to 49 Mb holds 9 of 89.5 Mb.
    want = 9 / 89.5
    got = by_length([dict(chrom='1', side='start', span=10_000_000)], {'1': pos}, {'1': [(40_000_000, 49_000_000, 'S')]})
    assert abs(got - want) < 1e-9, (got, want)
    print('fragile_sites demo ok')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        demo()
    else:
        main(*sys.argv[1:])

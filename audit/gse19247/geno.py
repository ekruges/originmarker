"""Genotypes and B-allele frequency from the cluster positions, plus the checks that say whether
they are trustworthy.

The clusters come from cluster.py. Everything below is downstream of them, so if the checks at the
bottom fail there is no point running anything else on this data.
"""
import json, os, sys
import numpy as np
from cluster import read_idat_np, sample_theta, MIN_R_FRAC

# A call is made only when one cluster is a clear winner. d1/d2 is the distance to the nearest
# cluster over the distance to the next: at 1.0 the point sits on the midpoint between two clouds
# and means nothing, so the bar is well inside that.
CALL_RATIO = 0.55

AA, AB, BB, NC = 0, 1, 2, -1


def load():
    z = np.load('clusters_gpl6985.npz', allow_pickle=True)
    return z['addr'], z['chrom'], z['pos'], z['mu']


def call(grn, red, addr, mu):
    """genotype (AA/AB/BB/NC) and BAF for one sample."""
    th, tot = sample_theta(grn, red, addr)
    med = np.median(tot[tot > 0]) if (tot > 0).any() else 1.0
    live = tot >= MIN_R_FRAC * med

    d = np.abs(th[:, None] - mu)                    # markers x 3
    order = np.argsort(d, axis=1)
    near = order[:, 0]
    d1 = np.take_along_axis(d, order[:, :1], 1)[:, 0]
    d2 = np.take_along_axis(d, order[:, 1:2], 1)[:, 0]
    gt = np.where(live & (d1 <= CALL_RATIO * np.maximum(d2, 1e-9)), near, NC).astype(np.int8)

    # BAF: piecewise linear through the three cluster positions, which is how a cluster file
    # defines it. Outside the homozygote clouds it saturates rather than extrapolating.
    lo, mid, hi = mu[:, 0], mu[:, 1], mu[:, 2]
    baf = np.where(
        th <= lo, 0.0,
        np.where(th >= hi, 1.0,
                 np.where(th < mid,
                          0.5 * (th - lo) / np.maximum(mid - lo, 1e-6),
                          0.5 + 0.5 * (th - mid) / np.maximum(hi - mid, 1e-6))))
    return gt, baf.astype(np.float32), live


def paths(s):
    g = ['idat/' + f for f in s['idats'] if 'Grn' in f]
    r = ['idat/' + f for f in s['idats'] if 'Red' in f]
    return (g[0], r[0]) if g and r and os.path.exists(g[0]) and os.path.exists(r[0]) else None


def main():
    addr, chrom, pos, mu = load()
    S = json.load(open('gpl6985_samples.json'))
    by = lambda pred: [s for s in S if pred(s)]
    fam = lambda s, n: f'family {n}' in s['src']

    print('=== CHECK 1: bulk DNA heterozygosity ===')
    print('a real bulk sample runs near 0.30; the global-threshold attempt gave 0.415 and 0.528')
    bulk = by(lambda s: 'buccal' in s['src'])
    gts = {}
    for s in bulk:
        p = paths(s)
        gt, baf, live = call(*p, addr, mu)
        gts[s['gsm']] = gt
        called = gt >= 0
        het = (gt == AB).sum() / max(called.sum(), 1)
        print(f"  {s['gsm']}  {s['src'][:34]:36} call {called.mean():.3f}  het {het:.4f}")

    print()
    print('=== CHECK 2: which buccal is the father, by his own sperm ===')
    print('the father cannot be an opposite homozygote to his own sperm; the mother can and will')
    for n in ('15', '28'):
        sp = [s for s in by(lambda s: 'sperm' in s['src']) if fam(s, n)][:8]
        spg = []
        for s in sp:
            p = paths(s)
            if p:
                spg.append(call(*p, addr, mu)[0])
        par = [s for s in bulk if fam(s, n)]
        out = []
        for s in par:
            gt = gts[s['gsm']]
            opp = tot = 0
            for sg in spg:
                both = (gt != NC) & (sg != NC) & (gt != AB) & (sg != AB)
                opp += int((both & (gt != sg)).sum())
                tot += int(both.sum())
            out.append((s['gsm'], opp / max(tot, 1), tot))
        out.sort(key=lambda x: x[1])
        sep = out[1][1] / max(out[0][1], 1e-9)
        print(f'  family {n}: ' + '  '.join(f'{g}={v:.4f}' for g, v, _ in out)
              + f'   -> father {out[0][0]}, separation {sep:.1f}x')

    print()
    print('=== CHECK 3: trisomy 21 vs euploid, on chr21 intensity ===')
    print('83 confirmed trisomy 21 single cells and the euploid cells, same chip, same pass')


if __name__ == '__main__':
    main()

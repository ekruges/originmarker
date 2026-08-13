"""Segment-level parental regions, with an empirical null that already contains WGA clustering.

Two failures preceded this. Testing 5-Mb bins independently against a MAD null gave 8 detections in
a cell that is uniparental on every bin. Replacing that with a binomial likelihood ratio gave 419
regions across nearly every cell, because the binomial assumes markers within a bin drop out
independently and amplification drops them in clustered runs.

Neither problem is the threshold. A real event spans many bins, so the evidence has to be POOLED
across them, and the null has to be built from the data rather than from a distribution that
assumes independence.

  score      per bin, the standardised shortfall of the heterozygous fraction against the same bin
             in the OTHER cells, so bin difficulty and cell quality are already divided out
  segment    the maximum-sum contiguous run of (score - PENALTY) along a chromosome, which is the
             classical maximum-segment-sum and finds the best interval without pre-set edges
  null       the same statistic on the SAME scores with bin order shuffled within the cell. That
             preserves each cell's marginal noise exactly and destroys only the spatial structure,
             which is the thing an event has and clustered dropout has less of.

Truth for sensitivity is material type from the submission sheet: polar bodies, MII oocyte biopsies
and isolated pronuclei carry one parental genome, so every chromosome in them is a true positive.
"""
import numpy as np

CACHE = "cells_14801.npz"
AUTOSOMES = [str(i) for i in range(1, 23)]
NC, AA, AB, BB = -1, 0, 1, 2
MIN_COV, HET_FRAC = 10, 0.30
BIN_MB, MIN_CALLED, MIN_COHORT = 5.0, 60, 8
PENALTY = 1.0          # per bin, so a segment must average better than this to grow
NPERM = 400
DIR_MARGIN = 0.60

# One parental genome by construction, from the GEO submission sheet.
UNIPARENTAL = {'14801-01', '14801-02', '14801-03', '14801-04', '14801-06',
               '14801-12', '14801-13', '14801-14', '14801-16'}


def max_segment(x):
    """Largest contiguous sum, with its bounds. Kadane, returning (score, i, j) inclusive."""
    best, bi, bj = 0.0, -1, -1
    cur, ci = 0.0, 0
    for i, v in enumerate(x):
        if cur <= 0:
            cur, ci = v, i
        else:
            cur += v
        if cur > best:
            best, bi, bj = cur, ci, i
    return best, bi, bj


def main():
    z = np.load(CACHE, allow_pickle=True)
    names, G, chrom, pos, mgt = z['names'], z['G'], z['chrom'], z['pos'], z['mgt']
    n = len(names)
    cov, het = (G != NC).sum(0), (G == AB).sum(0)
    pgt = np.where(mgt == AA, BB, AA)

    bounds, step = [], int(BIN_MB * 1e6)
    for c in AUTOSOMES:
        w = np.where(chrom == c)[0]
        if not len(w):
            continue
        p = pos[w]
        for lo in range(0, int(p.max()) + step, step):
            sel = w[(p >= lo) & (p < lo + step)]
            if len(sel) >= MIN_CALLED:
                bounds.append((c, sel[0], sel[-1]))
    F = np.full((n, len(bounds)), np.nan)
    P = np.full((n, len(bounds)), np.nan)
    for i in range(n):
        c = cov - (G[i] != NC)
        h = het - (G[i] == AB)
        obl = (c >= MIN_COV) & (h / np.maximum(c, 1) >= HET_FRAC)
        called = obl & (G[i] != NC)
        isab, ispat, ismat = called & (G[i] == AB), called & (G[i] == pgt), called & (G[i] == mgt)
        for j, (_, a, b) in enumerate(bounds):
            nn = int(called[a:b + 1].sum())
            if nn < MIN_CALLED:
                continue
            F[i, j] = isab[a:b + 1].sum() / nn
            nh = int(ispat[a:b + 1].sum() + ismat[a:b + 1].sum())
            if nh >= 20:
                P[i, j] = ispat[a:b + 1].sum() / nh

    S = np.zeros_like(F)
    for j in range(len(bounds)):
        ok = np.isfinite(F[:, j])
        if ok.sum() < MIN_COHORT:
            continue
        for i in np.where(ok)[0]:
            o = F[ok & (np.arange(n) != i), j]
            m = np.median(o)
            mad = max(np.median(np.abs(o - m)) * 1.4826, 0.02)
            S[i, j] = (m - F[i, j]) / mad          # positive = heterozygosity lost

    rng = np.random.default_rng(11)
    rows, seen = [], {}
    for i in range(n):
        allsc = S[i][np.isfinite(F[i])]
        for c in AUTOSOMES:
            idx = [j for j, b in enumerate(bounds) if b[0] == c and np.isfinite(F[i, j])]
            if len(idx) < 3:
                continue
            x = S[i, idx] - PENALTY
            sc, a, b = max_segment(x)
            if sc <= 0:
                continue
            null = np.array([max_segment(rng.permutation(allsc)[:len(x)] - PENALTY)[0]
                             for _ in range(NPERM)])
            p = (1 + (null >= sc).sum()) / (NPERM + 1)
            seen.setdefault(names[i], []).append(p)
            if p <= 0.01:
                jj = idx[a:b + 1]
                pv = float(np.nanmedian(P[i, jj]))
                call = ('MATERNAL copy missing' if pv >= DIR_MARGIN else
                        'PATERNAL copy missing' if pv <= 1 - DIR_MARGIN else
                        'uniparental, parent unresolved')
                s0, e0 = bounds[jj[0]][1], bounds[jj[-1]][2]
                rows.append((names[i], c, int(pos[s0]), int(pos[e0]), sc, p, pv, call,
                             float(np.nanmedian(F[i, jj]))))

    print(f'{n} cells, {len(bounds)} bins, {sum(len(v) for v in seen.values())} chromosome scans')
    uni = [p for k, v in seen.items() if k in UNIPARENTAL for p in v]
    bip = [p for k, v in seen.items() if k not in UNIPARENTAL for p in v]
    print(f'\nchromosomes in one-genome material: {len(uni)}   in the rest: {len(bip)}')
    for cut in (0.05, 0.01, 0.005):
        print(f'  p <= {cut:<6} sensitivity {np.mean(np.array(uni) <= cut):.3f}   '
              f'false-positive {np.mean(np.array(bip) <= cut):.4f}')

    print(f'\n{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"span":>7} {"het":>6} '
          f'{"score":>7} {"p":>7} {"his%":>6}  call')
    for r in sorted(rows, key=lambda x: x[5]):
        print(f'{r[0]:<11} {r[1]:>3} {r[2]/1e6:>9.2f} {r[3]/1e6:>9.2f} {(r[3]-r[2])/1e6:>6.1f}M '
              f'{r[8]:>6.3f} {r[4]:>7.1f} {r[5]:>7.4f} {r[6]:>6.2f}  {r[7]}')
    mat = sum(1 for r in rows if r[7].startswith('MATERNAL'))
    pat = sum(1 for r in rows if r[7].startswith('PATERNAL'))
    print(f'\n{len(rows)} regions: {mat} maternal-missing, {pat} paternal-missing, '
          f'{len(rows)-mat-pat} unresolved')
    with open('segment_regions.tsv', 'w') as f:
        f.write('cell\tchr\tstart_bp\tend_bp\tspan_Mb\thet_frac\tscore\tp\this_allele_frac\tcall\n')
        for r in rows:
            f.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{(r[3]-r[2])/1e6:.3f}\t{r[8]:.4f}\t'
                    f'{r[4]:.1f}\t{r[5]:.4f}\t{r[6]:.3f}\t{r[7]}\n')
    print('wrote segment_regions.tsv')


if __name__ == '__main__':
    main()

"""Parent-labelled regions by binomial likelihood ratio, and the sensitivity that buys.

The bin-fraction z-score discards the count. With about 110 obligate-het markers in a 5-Mb bin the
sampling deviation of the fraction is large, a MAD across cells absorbs it into the null, and the
result was 8 flagged bins of about 560 in a cell that is paternal-absent EVERYWHERE. Direction was
right and sensitivity was not.

A binomial likelihood ratio uses the counts, and two nuisance factors are divided out first rather
than left in the null:

  cell quality   a poorly amplified cell reads fewer heterozygotes on every bin
  bin difficulty some bins genotype badly in every cell

so the expected heterozygous count for cell i in bin j is that cell's own genome-wide rate scaled
by how this bin behaves in the OTHER cells. Departure from that is the statistic.

  expected_ij = n_ij * r_i * (q_j / q)      r_i cell rate, q_j cohort rate here, q cohort overall

Truth for measuring sensitivity is 14801-02 and 14801-16, which carry no paternal genome anywhere,
so every scoreable bin in them is a true positive and every bin in a biparental cell is a true
negative for the paternal direction.
"""
import numpy as np

CACHE = "cells_14801.npz"
AUTOSOMES = [str(i) for i in range(1, 23)]
NC, AA, AB, BB = -1, 0, 1, 2
MIN_COV, HET_FRAC = 10, 0.30
BIN_MB, MIN_CALLED = 5.0, 60
MIN_COHORT = 8
LRT = 12.0                # about p 5e-4 on 1 df, before multiplicity
DIR_MARGIN = 0.60
MATERNAL_ONLY = {'14801-02', '14801-16'}     # established independently, see FINDINGS


def lrt(k, n, e):
    """2 log likelihood ratio for k of n against expectation e, one-sided on a DEFICIT."""
    if n <= 0 or k >= e:
        return 0.0
    p, p0 = max(k / n, 1e-9), min(max(e / n, 1e-9), 1 - 1e-9)
    t = k * np.log(p / p0)
    if k < n:
        t += (n - k) * np.log((1 - p) / (1 - p0))
    return max(2.0 * t, 0.0)


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

    K = np.zeros((n, len(bounds)))      # heterozygous count
    N = np.zeros((n, len(bounds)))      # called obligate-het count
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
            N[i, j], K[i, j] = nn, int(isab[a:b + 1].sum())
            nh = int(ispat[a:b + 1].sum() + ismat[a:b + 1].sum())
            if nh >= 20:
                P[i, j] = ispat[a:b + 1].sum() / nh

    ok = N >= MIN_CALLED
    r = np.where(N.sum(1) > 0, K.sum(1) / np.maximum(N.sum(1), 1), np.nan)   # cell rate
    q_all = K[ok].sum() / N[ok].sum()
    print(f'{n} cells, {len(bounds)} bins, cohort heterozygous rate {q_all:.4f}')

    S = np.zeros_like(K)
    for j in range(len(bounds)):
        m = ok[:, j]
        if m.sum() < MIN_COHORT:
            continue
        for i in np.where(m)[0]:
            o = m & (np.arange(n) != i)
            qj = K[o, j].sum() / max(N[o, j].sum(), 1)          # bin difficulty, cell held out
            e = N[i, j] * r[i] * (qj / q_all)
            S[i, j] = lrt(K[i, j], N[i, j], e)

    # sensitivity against the two cells that carry no paternal genome anywhere
    tp = np.array([nm in MATERNAL_ONLY for nm in names])
    pos_bins = ok & tp[:, None]
    neg_bins = ok & ~tp[:, None]
    print(f'\nscoreable bins: {int(pos_bins.sum())} in maternal-only cells (all true positives), '
          f'{int(neg_bins.sum())} in biparental cells')
    print(f'{"LRT cut":>8} {"sensitivity":>12} {"false-positive rate":>21}')
    for cut in (6, 9, 12, 16, 20, 30):
        se = (S[pos_bins] >= cut).mean()
        fp = (S[neg_bins] >= cut).mean()
        print(f'{cut:>8} {se:>12.3f} {fp:>21.4f}')

    print(f'\nregions at LRT >= {LRT}')
    print(f'{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"het":>6} {"exp":>6} '
          f'{"LRT":>7} {"his%":>6}  call')
    rows = []
    for i in range(n):
        j = 0
        while j < len(bounds):
            if S[i, j] < LRT:
                j += 1; continue
            k = j
            while k + 1 < len(bounds) and S[i, k + 1] >= LRT and bounds[k + 1][0] == bounds[j][0]:
                k += 1
            kk, nn = K[i, j:k + 1].sum(), N[i, j:k + 1].sum()
            pv = float(np.nanmedian(P[i, j:k + 1]))
            call = ('MATERNAL copy missing' if pv >= DIR_MARGIN else
                    'PATERNAL copy missing' if pv <= 1 - DIR_MARGIN else
                    'uniparental, parent unresolved')
            a, b = bounds[j][1], bounds[k][2]
            print(f'{names[i]:<11} {bounds[j][0]:>3} {pos[a]/1e6:>9.2f} {pos[b]/1e6:>9.2f} '
                  f'{kk/max(nn,1):>6.3f} {r[i]:>6.3f} {S[i, j:k+1].max():>7.1f} {pv:>6.2f}  {call}')
            rows.append((names[i], bounds[j][0], int(pos[a]), int(pos[b]),
                         kk / max(nn, 1), float(S[i, j:k + 1].max()), pv, call))
            j = k + 1
    mat = sum(1 for x in rows if x[7].startswith('MATERNAL'))
    pat = sum(1 for x in rows if x[7].startswith('PATERNAL'))
    print(f'\n{len(rows)} regions: {mat} maternal-missing, {pat} paternal-missing, '
          f'{len(rows)-mat-pat} unresolved')
    with open('lrt_regions.tsv', 'w') as f:
        f.write('cell\tchr\tstart_bp\tend_bp\thet_frac\tLRT\this_allele_frac\tcall\n')
        for x in rows:
            f.write(f'{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]:.4f}\t{x[5]:.1f}\t{x[6]:.3f}\t{x[7]}\n')
    print('wrote lrt_regions.tsv')


if __name__ == '__main__':
    main()

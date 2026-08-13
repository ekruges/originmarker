"""Parent-labelled regions from GENOTYPE calls, which is the channel that survives WGA.

A median B-allele frequency on this material reports the dropout rate: only 21-53% of obligate-het
markers are called heterozygous and the rest sit at zero, so no threshold recovers dosage from it.
Genotype calls are insensitive to where BAF lands inside a heterozygous call, and they carry both
pieces of information needed:

  IS there an event    the fraction of obligate-het markers CALLED heterozygous falls, because a
                       cell with one parental copy cannot be heterozygous there
  WHICH parent         among the markers that are NOT heterozygous, whether they carry the
                       mother's allele or the father's. A missing MATERNAL copy leaves his allele.

Both are read against the same window in the other cells, each cell held out of its own null, so a
region has to depart from where the cohort sits rather than from the cell's own average. That is
what stops recurrent platform failure being reported as biology.
"""
import numpy as np

CACHE = "cells_14801.npz"
AUTOSOMES = [str(i) for i in range(1, 23)]
NC, AA, AB, BB = -1, 0, 1, 2

MIN_COV, HET_FRAC = 10, 0.30
# 0.30 admits some heterozygous fathers, whose markers a normal cell reads AB only half the time.
# That lowers the baseline heterozygous fraction UNIFORMLY across cells, which the cohort null
# absorbs, and it buys the marker density a genomic bin needs: 78k markers rather than 14k.
# Bins are genomic, not marker-index, because a marker-index window over all maternal-homozygous
# markers contains only about 2% obligate-het ones and never reaches a usable count. Fixed Mb bins
# are also what a region has to be reported as.
# 5 Mb rather than 2: at 2 Mb a bin holds about 43 obligate-het markers and the sampling standard
# deviation of the heterozygous fraction is 0.076, which is the same size as a real event's signal.
# 5 Mb roughly doubles the markers and halves that noise.
BIN_MB = 5.0
MIN_CALLED = 60          # called obligate-het markers in a bin before it is scored
MIN_Z_AB = 4.0           # departure of the heterozygous fraction from the cohort
MIN_DROP = 0.15          # and an absolute fall, so a tight cohort cannot manufacture one
# Single bins are allowed: 5 Mb is already a reportable region, and requiring two adjacent
# ones discarded every call at z >= 4 on this cohort.
MIN_WIN, MIN_COHORT = 1, 8
DIR_MARGIN = 0.60        # of non-het calls carrying one parent's allele before naming that parent


def main():
    z = np.load(CACHE, allow_pickle=True)
    names, G, chrom, pos, mgt = z['names'], z['G'], z['chrom'], z['pos'], z['mgt']
    n = len(names)
    cov = (G != NC).sum(0)
    het = (G == AB).sum(0)
    pgt = np.where(mgt == AA, BB, AA)          # father's allele at an obligate-het marker

    bounds = []
    step = int(BIN_MB * 1e6)
    for c in AUTOSOMES:
        w = np.where(chrom == c)[0]
        if not len(w):
            continue
        p = pos[w]
        for lo in range(0, int(p.max()) + step, step):
            sel = w[(p >= lo) & (p < lo + step)]
            if len(sel) >= MIN_CALLED:
                bounds.append((c, sel[0], sel[-1]))
    print(f'{n} cells, {len(bounds)} bins of {BIN_MB} Mb')

    HET = np.full((n, len(bounds)), np.nan)    # heterozygous fraction
    PAT = np.full((n, len(bounds)), np.nan)    # of non-het calls, fraction carrying HIS allele
    for i in range(n):
        c = cov - (G[i] != NC)
        h = het - (G[i] == AB)
        obl = (c >= MIN_COV) & (h / np.maximum(c, 1) >= HET_FRAC)
        called = obl & (G[i] != NC)
        isab = called & (G[i] == AB)
        ispat = called & (G[i] == pgt)
        ismat = called & (G[i] == mgt)
        for j, (_, a, b) in enumerate(bounds):
            nc = int(called[a:b + 1].sum())
            if nc < MIN_CALLED:
                continue
            HET[i, j] = isab[a:b + 1].sum() / nc
            nh = int(ispat[a:b + 1].sum() + ismat[a:b + 1].sum())
            if nh >= 20:
                PAT[i, j] = ispat[a:b + 1].sum() / nh

    rows = []
    flags = np.zeros_like(HET, dtype=bool)
    for j in range(HET.shape[1]):
        col = HET[:, j]; ok = np.isfinite(col)
        if ok.sum() < MIN_COHORT:
            continue
        for i in np.where(ok)[0]:
            other = col[ok & (np.arange(n) != i)]
            med = np.median(other)
            mad = max(np.median(np.abs(other - med)) * 1.4826, 0.02)
            # one-sided: an event REMOVES heterozygosity, it never adds it
            flags[i, j] = (med - col[i]) / mad >= MIN_Z_AB and (med - col[i]) >= MIN_DROP

    print(f'\n{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"span":>7} '
          f'{"het":>6} {"cohort":>7} {"z":>6} {"his%":>6}  parental call')
    print('-' * 104)
    for i in range(n):
        j = 0
        while j < HET.shape[1]:
            if not flags[i, j]:
                j += 1; continue
            k = j
            while k + 1 < HET.shape[1] and flags[i, k + 1] and bounds[k + 1][0] == bounds[j][0]:
                k += 1
            if k - j + 1 >= MIN_WIN:
                hv = float(np.nanmedian(HET[i, j:k + 1]))
                coh = float(np.nanmedian([np.nanmedian(np.delete(HET[:, q], i))
                                          for q in range(j, k + 1)]))
                mad = max(np.median(np.abs(HET[:, j][np.isfinite(HET[:, j])]
                                           - np.nanmedian(HET[:, j]))) * 1.4826, 0.02)
                zz = (coh - hv) / mad
                pv = float(np.nanmedian(PAT[i, j:k + 1]))
                if not np.isfinite(pv):
                    call = 'parent not determined, too few homozygous calls'
                elif pv >= DIR_MARGIN:
                    call = 'MATERNAL copy missing'
                elif pv <= 1 - DIR_MARGIN:
                    call = 'PATERNAL copy missing'
                else:
                    call = 'uniparental, parent unresolved'
                a, b = bounds[j][1], bounds[k][2]
                print(f'{names[i]:<11} {bounds[j][0]:>3} {pos[a]/1e6:>9.2f} {pos[b]/1e6:>9.2f} '
                      f'{(pos[b]-pos[a])/1e6:>6.1f}M {hv:>6.3f} {coh:>7.3f} {zz:>6.1f} '
                      f'{pv:>6.2f}  {call}')
                rows.append((names[i], bounds[j][0], int(pos[a]), int(pos[b]), hv, coh, zz, pv, call))
            j = k + 1
    print('-' * 104)
    mat = sum(1 for r in rows if r[8].startswith('MATERNAL'))
    pat = sum(1 for r in rows if r[8].startswith('PATERNAL'))
    print(f'{len(rows)} regions: {mat} maternal-missing, {pat} paternal-missing, '
          f'{len(rows)-mat-pat} unresolved')
    with open('geno_regions.tsv', 'w') as f:
        f.write('cell\tchr\tstart_bp\tend_bp\tspan_Mb\thet_frac\tcohort_het\tz\this_allele_frac\tcall\n')
        for r in rows:
            f.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{(r[3]-r[2])/1e6:.3f}\t{r[4]:.4f}\t'
                    f'{r[5]:.4f}\t{r[6]:.1f}\t{r[7]:.3f}\t{r[8]}\n')
    print('wrote geno_regions.tsv')


if __name__ == '__main__':
    main()

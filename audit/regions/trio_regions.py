"""Both-direction parental regions on experiment 14801, with the three confounds closed.

  1  external per-window null, so recurrent platform failure is not read as a parental result
  2  obligate-het markers, so a paternal loss has the same headroom as a maternal one
  3  the father called LEAVE-ONE-OUT and from HETEROZYGOUS children only, so no cell contributes
     to selecting the markers it is then scored on

Only AB in a child shows the father supplying an allele the mother could not. Counting a child's
opposite homozygote as evidence about the father instead selects markers where the MATERNAL copy
is absent, which is the circularity this replaces: it put paternal share at 0.96-0.98 in 27 of 33
cells where a biparental cell reads about 0.5.

Pass 1 caches genotypes and BAF aligned to one marker index. Everything after is arithmetic.
"""
import glob, os, sys
import numpy as np
from regions import read_probes

USB = "/Volumes/SANDISK USB/SNP array data"
CACHE = "cells_14801.npz"
AUTOSOMES = [str(i) for i in range(1, 23)]

MAX_OPP = 0.15
MIN_COV = 10        # other cells covering a marker before the father is called there
# Calibrated on this data, not assumed. At a TRUE obligate-het marker a cell reads AB only when
# both alleles survive amplification, so heavy dropout puts the observed fraction near 0.30 rather
# than near 1. Sweeping the threshold, 0.30 recovers 77,993 markers against the ~80,000 the two
# parents' genotype frequencies predict, while 0.60 keeps 4,401 and throws the class away.
HET_FRAC = 0.30
WIN, STEP = 200, 100
MIN_Z, MIN_ABS, MIN_WIN = 5.0, 0.12, 2
MIN_COHORT = 8

NC, AA, AB, BB = -1, 0, 1, 2


def build():
    mumf = glob.glob(f'{USB}/ROBLES/**/14801-05_*.probes', recursive=True)[0]
    mum = read_probes(mumf)
    keys = [k for k, v in mum.items() if v[3] in ('0', '2') and v[0] in AUTOSOMES]
    keys.sort(key=lambda k: (AUTOSOMES.index(mum[k][0]), mum[k][1]))
    idx = {k: i for i, k in enumerate(keys)}
    chrom = np.array([mum[k][0] for k in keys])
    pos = np.array([mum[k][1] for k in keys], dtype=np.int64)
    mgt = np.array([0 if mum[k][3] == '0' else 2 for k in keys], dtype=np.int8)
    print(f'mother: {len(keys)} autosomal homozygous markers', flush=True)

    names, G, B = [], [], []
    for cf in sorted(f for f in glob.glob(f'{USB}/ROBLES/**/14801-*.probes', recursive=True)
                     if '14801-05_' not in f):
        name = cf.split('/')[-1].split('_')[0]
        s = read_probes(cf)
        g = np.full(len(keys), NC, dtype=np.int8)
        b = np.full(len(keys), np.nan, dtype=np.float32)
        for k, i in idx.items():
            w = s.get(k)
            if w and w[3] != '-1':
                g[i] = int(w[3])
                if np.isfinite(w[2]):
                    b[i] = w[2]
        called = g != NC
        opp = np.sum(called & (g != AB) & (g != mgt))
        tot = np.sum(called & (g != AB))
        if not tot or opp / tot > MAX_OPP:
            print(f'  {name} excluded ({opp/max(tot,1):.3f})', flush=True)
            continue
        names.append(name); G.append(g); B.append(b)
        print(f'  {name} in ({called.sum()} called)', flush=True)
    np.savez_compressed(CACHE, names=np.array(names), G=np.array(G), B=np.array(B),
                        chrom=chrom, pos=pos, mgt=mgt)
    print(f'cached {len(names)} cells', flush=True)


def main():
    if not os.path.exists(CACHE) or '--rebuild' in sys.argv:
        build()
    z = np.load(CACHE, allow_pickle=True)
    names, G, B, chrom, pos, mgt = z['names'], z['G'], z['B'], z['chrom'], z['pos'], z['mgt']
    n = len(names)
    print(f'\n{n} cells, {G.shape[1]} maternal-homozygous markers')

    cov = (G != NC).sum(0)
    het = (G == AB).sum(0)

    # paternal share per cell at that cell's own leave-one-out obligate-het markers
    share = np.full(G.shape, np.nan, dtype=np.float32)
    obl_n = []
    for i in range(n):
        c = cov - (G[i] != NC)
        h = het - (G[i] == AB)
        ok = (c >= MIN_COV) & (h / np.maximum(c, 1) >= HET_FRAC)
        obl_n.append(int(ok.sum()))
        # Orient to the FATHER. At an obligate-het marker he carries the allele the mother does
        # not, so where she is AA his allele is B and his share is the B-allele frequency itself;
        # where she is BB his allele is A and his share is its complement. Taking the maternal
        # orientation here and calling it paternal inverts every call that follows.
        s = np.where(mgt == AA, B[i], 1.0 - B[i])
        share[i] = np.where(ok & (G[i] != NC) & np.isfinite(B[i]), s, np.nan)
    print(f'obligate-het markers per cell (leave-one-out): '
          f'median {int(np.median(obl_n))}, range {min(obl_n)}-{max(obl_n)}')

    print(f'\n{"cell":<11} {"informative":>11} {"centre":>8} {"p1":>8} {"p99":>8}  symmetry')
    centres = np.full(n, np.nan)
    for i in range(n):
        v = share[i][np.isfinite(share[i])]
        if len(v) < 2000:
            continue
        centres[i] = np.median(v)
        d = v - centres[i]
        lo, hi = np.percentile(d, 1), np.percentile(d, 99)
        print(f'{names[i]:<11} {len(v):>11} {centres[i]:>8.4f} {lo:>+8.4f} {hi:>+8.4f}  '
              f'{"SYMMETRIC" if abs(hi/min(lo,-1e-9)) > 0.4 else "skewed"}')

    # window medians, one fixed grid so the cohort null compares like with like
    bounds = []
    for c in AUTOSOMES:
        w = np.where(chrom == c)[0]
        for a in range(0, max(len(w) - WIN + 1, 0), STEP):
            bounds.append((c, w[a], w[min(a + WIN, len(w) - 1)]))
    M = np.full((n, len(bounds)), np.nan)
    for i in range(n):
        if not np.isfinite(centres[i]):
            continue
        for j, (_, a, b) in enumerate(bounds):
            seg = share[i][a:b]
            seg = seg[np.isfinite(seg)]
            if len(seg) >= WIN // 4:
                M[i, j] = np.median(seg) - centres[i]

    print(f'\n{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"span":>7} {"dev":>7} '
          f'{"z":>7}  parental call')
    print('-' * 92)
    out = []
    flags = np.zeros_like(M, dtype=bool)
    for j in range(M.shape[1]):
        col = M[:, j]; ok = np.isfinite(col)
        if ok.sum() < MIN_COHORT:
            continue
        for i in np.where(ok)[0]:
            other = col[ok & (np.arange(n) != i)]
            med = np.median(other)
            mad = max(np.median(np.abs(other - med)) * 1.4826, 0.02)
            flags[i, j] = abs((col[i] - med) / mad) >= MIN_Z and abs(col[i] - med) >= MIN_ABS
    for i in range(n):
        j = 0
        while j < M.shape[1]:
            if not flags[i, j]:
                j += 1; continue
            k = j
            while k + 1 < M.shape[1] and flags[i, k + 1] and bounds[k + 1][0] == bounds[j][0]:
                k += 1
            if k - j + 1 >= MIN_WIN:
                dev = float(np.nanmedian(M[i, j:k + 1]))
                coh = float(np.nanmedian([np.nanmedian(np.delete(M[:, q], i))
                                          for q in range(j, k + 1)]))
                mad = max(np.median(np.abs(M[:, j][np.isfinite(M[:, j])]
                                           - np.nanmedian(M[:, j]))) * 1.4826, 0.02)
                zz = (dev - coh) / mad
                call = 'MATERNAL copy missing' if dev > coh else 'PATERNAL copy missing'
                a, b = bounds[j][1], bounds[k][2]
                print(f'{names[i]:<11} {bounds[j][0]:>3} {pos[a]/1e6:>9.2f} {pos[b]/1e6:>9.2f} '
                      f'{(pos[b]-pos[a])/1e6:>6.1f}M {dev:>+7.3f} {zz:>+7.1f}  {call}')
                out.append((names[i], bounds[j][0], pos[a], pos[b], dev, zz, call))
            j = k + 1
    print('-' * 92)
    mat = sum(1 for o in out if o[6].startswith('MATERNAL'))
    print(f'{len(out)} regions: {mat} maternal-missing, {len(out)-mat} paternal-missing')
    with open('trio_regions.tsv', 'w') as f:
        f.write('cell\tchr\tstart_bp\tend_bp\tspan_Mb\tdev\tz\tcall\n')
        for o in out:
            f.write(f'{o[0]}\t{o[1]}\t{o[2]}\t{o[3]}\t{(o[3]-o[2])/1e6:.3f}\t'
                    f'{o[4]:+.4f}\t{o[5]:+.1f}\t{o[6]}\n')
    print('wrote trio_regions.tsv')


if __name__ == '__main__':
    main()

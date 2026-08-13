"""Parent-labelled regions with a PER-REGION EXTERNAL NULL.

Scoring each cell against its own genome-wide centre produces 154 regions on this experiment, all
154 of them 'maternal copy missing' and none paternal, with chr16:0-1Mb appearing in 18 of 32
cells. A biological process does not run one way in every cell and land on the same megabase in
half of them. Those are places where genotyping fails on this platform - GC-rich subtelomeres, the
HLA region at chr6:32Mb, the imprinted cluster at chr11:2Mb - and the failure drops maternal-
matching signal, so a self-derived null converts it into a parental result.

The null here is therefore the SAME WINDOW IN THE OTHER CELLS, not the rest of the same cell. A
region survives only if it departs from where that window sits across the cohort, which is what
removes an artefact common to all of them. Cells are held out of their own null so a real event
cannot suppress its own evidence.

Two passes: window medians are cached once, then the null is applied. The cache makes the second
half re-runnable without re-reading 27 million marker rows.
"""
import glob, os, sys
import numpy as np
from regions import read_probes

USB = "/Volumes/SANDISK USB/SNP array data"
CACHE = "windows_14801.npz"
AUTOSOMES = [str(i) for i in range(1, 23)]

WIN, STEP = 300, 150
MAX_OPP = 0.15
MIN_MARKERS = 20_000
# A window must clear BOTH: this many robust SDs from the cohort at that position, and this much
# raw share. The absolute floor stops a very tight cohort window turning noise into significance.
MIN_Z = 5.0
MIN_ABS = 0.10
MIN_WIN = 2
MIN_COHORT = 8      # cells needed at a window before it has a usable null


def build_cache():
    mumf = glob.glob(f'{USB}/ROBLES/**/14801-05_*.probes', recursive=True)[0]
    mum = read_probes(mumf)
    hom = {k: v for k, v in mum.items() if v[3] in ('0', '2') and v[0] in AUTOSOMES}
    print(f'maternal reference: {len(hom)} autosomal homozygous markers', flush=True)

    # A single fixed window grid, shared by every cell, or the cohort null compares
    # different stretches of genome to each other.
    grid = {}
    for c in AUTOSOMES:
        ps = sorted(p for (ch, p, _, _) in hom.values() if ch == c)
        grid[c] = np.array(ps, dtype=np.int64)

    cells, mats = [], []
    for cf in sorted(f for f in glob.glob(f'{USB}/ROBLES/**/14801-*.probes', recursive=True)
                     if '14801-05_' not in f):
        name = cf.split('/')[-1].split('_')[0]
        s = read_probes(cf)
        opp = tot = 0
        for k, v in hom.items():
            w = s.get(k)
            if w and w[3] in ('0', '2'):
                tot += 1
                opp += (w[3] != v[3])
        if not tot or opp / tot > MAX_OPP:
            print(f'  {name}: excluded, {opp/max(tot,1):.3f} opposite-homozygote', flush=True)
            continue
        per = {}
        n = 0
        for k, (c, pos, _, gm) in hom.items():
            w = s.get(k)
            if not w or w[3] == '-1' or not np.isfinite(w[2]):
                continue
            per.setdefault(c, {})[pos] = (1.0 - w[2]) if gm == '0' else w[2]
            n += 1
        if n < MIN_MARKERS:
            continue
        allv = np.array([v for d in per.values() for v in d.values()])
        centre = float(np.median(allv))
        row = []
        for c in AUTOSOMES:
            g = grid[c]
            sh = np.array([per.get(c, {}).get(int(p), np.nan) for p in g])
            for i in range(0, max(len(g) - WIN + 1, 0), STEP):
                seg = sh[i:i + WIN]
                seg = seg[np.isfinite(seg)]
                row.append(np.median(seg) - centre if len(seg) >= WIN // 3 else np.nan)
        cells.append(name); mats.append(row)
        print(f'  {name}: {n} markers, centre {centre:.4f}', flush=True)

    meta = []
    for c in AUTOSOMES:
        g = grid[c]
        for i in range(0, max(len(g) - WIN + 1, 0), STEP):
            meta.append((c, int(g[i]), int(g[min(i + WIN, len(g) - 1)])))
    np.savez_compressed(CACHE, cells=np.array(cells), mat=np.array(mats, dtype=np.float64),
                        chrom=np.array([m[0] for m in meta]),
                        start=np.array([m[1] for m in meta]),
                        end=np.array([m[2] for m in meta]))
    print(f'cached {len(cells)} cells x {len(meta)} windows', flush=True)


def report():
    z = np.load(CACHE, allow_pickle=True)
    cells, mat, chrom, start, end = z['cells'], z['mat'], z['chrom'], z['start'], z['end']
    print(f'\n{len(cells)} cells, {mat.shape[1]} windows\n')

    flags = np.zeros_like(mat, dtype=bool)
    zs = np.full_like(mat, np.nan)
    for j in range(mat.shape[1]):
        col = mat[:, j]
        ok = np.isfinite(col)
        if ok.sum() < MIN_COHORT:
            continue
        for i in np.where(ok)[0]:
            other = col[ok & (np.arange(len(col)) != i)]   # hold this cell out of its own null
            med = np.median(other)
            mad = np.median(np.abs(other - med)) * 1.4826
            if mad < 0.02:
                mad = 0.02          # floor, or a tight window turns noise into significance
            zz = (col[i] - med) / mad
            zs[i, j] = zz
            flags[i, j] = abs(zz) >= MIN_Z and abs(col[i] - med) >= MIN_ABS

    print(f'{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"span":>7} {"dev":>7} '
          f'{"cohort":>7} {"z":>7}  parental call')
    print('-' * 96)
    total = 0
    for i, cell in enumerate(cells):
        j = 0
        while j < mat.shape[1]:
            if not flags[i, j]:
                j += 1; continue
            k = j
            while k + 1 < mat.shape[1] and flags[i, k + 1] and chrom[k + 1] == chrom[j]:
                k += 1
            if k - j + 1 >= MIN_WIN:
                dev = float(np.nanmedian(mat[i, j:k + 1]))
                coh = float(np.nanmedian([np.nanmedian(np.delete(mat[:, q], i))
                                          for q in range(j, k + 1)]))
                zz = float(np.nanmedian(zs[i, j:k + 1]))
                call = 'PATERNAL copy missing' if dev > coh else 'MATERNAL copy missing'
                print(f'{cell:<11} {chrom[j]:>3} {start[j]/1e6:>9.2f} {end[k]/1e6:>9.2f} '
                      f'{(end[k]-start[j])/1e6:>6.1f}M {dev:>+7.3f} {coh:>+7.3f} {zz:>+7.1f}  {call}')
                total += 1
            j = k + 1
    print('-' * 96)
    print(f'{total} regions survive a per-window external null '
          f'(>= {MIN_Z} robust SD and >= {MIN_ABS} share from the cohort at that position)')


if __name__ == '__main__':
    if not os.path.exists(CACHE) or '--rebuild' in sys.argv:
        build_cache()
    report()

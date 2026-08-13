"""Copy-number regions from log2R, with a panel-of-normals null. No parental genotype needed.

The deposited copy_number column is not a segmented call - 55% of markers read 3.0 and 8% read 2.0
on a normal array - so segmentation is done here from log2R.

TWO LESSONS FROM THIS AUDIT ARE BUILT IN.
  Chromosomes carry their own intensity baseline on this chemistry: GC-rich chr19 and chr22 read
  -0.57 and -0.52 with no dosage change at all, so a genome-wide baseline calls a loss in every
  cell. Every window is therefore scored against the SAME WINDOW ACROSS THE OTHER ARRAYS, each
  array held out of its own null.
  Recurrent positions are artefact, not biology. The cohort median removes them by construction.
"""
import glob, sys
import numpy as np

USB = "/Volumes/SANDISK USB/SNP array data"
AUT = [str(i) for i in range(1, 23)]
BIN = 1_000_000
MIN_MARK = 100
MIN_WIN = 3            # consecutive windows, so one noisy bin is not a region
Z = 5.0
ABS = 0.25             # log2 units, on top of the z, so a tight cohort cannot make noise loud


def read(path, step=4):
    ch, po, lr = [], [], []
    with open(path) as f:
        f.readline()
        for i, l in enumerate(f):
            if i % step:
                continue
            p = l.rstrip('\n').split('\t')
            if len(p) < 5 or p[1] not in AUT:
                continue
            try:
                v = float(p[3])
            except ValueError:
                continue
            if np.isfinite(v):
                ch.append(p[1]); po.append(int(p[2])); lr.append(v)
    return np.array(ch), np.array(po, np.int64), np.array(lr, np.float32)


def main(pattern, out):
    files = sorted(glob.glob(pattern, recursive=True))
    print(f'{len(files)} arrays', flush=True)
    ch0, po0, _ = read(files[0])
    grid = []
    for c in AUT:
        w = np.where(ch0 == c)[0]
        if not len(w):
            continue
        p = po0[w]
        for lo in range(0, int(p.max()) + BIN, BIN):
            s = w[(p >= lo) & (p < lo + BIN)]
            if len(s) >= MIN_MARK:
                grid.append((c, lo, lo + BIN))
    print(f'{len(grid)} windows of {BIN/1e6:.0f} Mb', flush=True)

    names, M = [], []
    for i, f in enumerate(files):
        ch, po, lr = read(f)
        if len(lr) < 20000:
            continue
        base = float(np.median(lr))
        row = []
        for c, a, b in grid:
            m = (ch == c) & (po >= a) & (po < b)
            row.append(float(np.median(lr[m])) - base if m.sum() >= MIN_MARK // 2 else np.nan)
        names.append(f.split('/')[-1].split('_')[0]); M.append(row)
        if i % 20 == 0:
            print(f'  {i}/{len(files)}', flush=True)
    M = np.array(M)
    print(f'{len(names)} arrays scored', flush=True)

    rows = []
    for j in range(M.shape[1]):
        col = M[:, j]
        ok = np.isfinite(col)
        if ok.sum() < 8:
            continue
    for i, nm in enumerate(names):
        flags = np.zeros(M.shape[1], bool)
        devs = np.full(M.shape[1], np.nan)
        for j in range(M.shape[1]):
            col = M[:, j]
            ok = np.isfinite(col) & (np.arange(len(col)) != i)
            if ok.sum() < 8 or not np.isfinite(M[i, j]):
                continue
            med = np.median(col[ok])
            mad = max(np.median(np.abs(col[ok] - med)) * 1.4826, 0.05)
            d = M[i, j] - med
            devs[j] = d
            flags[j] = abs(d / mad) >= Z and abs(d) >= ABS
        j = 0
        while j < len(flags):
            if not flags[j]:
                j += 1; continue
            k = j
            while k + 1 < len(flags) and flags[k + 1] and grid[k + 1][0] == grid[j][0]:
                k += 1
            if k - j + 1 >= MIN_WIN:
                d = float(np.nanmedian(devs[j:k + 1]))
                rows.append((nm, grid[j][0], grid[j][1], grid[k][2], d,
                             'gain' if d > 0 else 'loss'))
            j = k + 1
    with open(out, 'w') as f:
        f.write('cell\tchr\tstart_bp\tend_bp\tlog2R_dev\tkind\n')
        for r in rows:
            f.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]}\n')
    g = sum(1 for r in rows if r[5] == 'gain')
    print(f'\n{len(rows)} regions: {g} gains, {len(rows)-g} losses, '
          f'across {len({r[0] for r in rows})} of {len(names)} arrays')
    print(f'wrote {out}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

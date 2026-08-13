"""Resolution against a BULK-genotyped father, which is what the method was specified for.

Every resolution figure so far was measured against a father reconstructed from siblings, and the
contaminated marker class that produced was shown to be the limit. The EYS sperm donor is bulk
gDNA on this platform at full density and five arrays of the 32477 series are his children, so the
obligate-het class can be DEFINED here rather than inferred.

With one parent the informative set is his homozygous markers: a cell carrying only his genome
cannot be heterozygous there, and a biparental cell reads heterozygous wherever the mother differs.
That is exactly obligateHet.ts's one-parent path, whose stated precondition is a bulk reference,
and section E2 measured that a single-cell reference cannot substitute.

Reported per array, per chromosome and per 5-Mb bin, so the resolution ceiling measured earlier
against a reconstructed father can be compared like for like.
"""
import gzip, glob
import numpy as np

USB = "/Volumes/SANDISK USB/SNP array data"
AUTOSOMES = [str(i) for i in range(1, 23)]
BIN_MB, MIN_BIN = 5.0, 40


def read(p, skip=1):
    op = gzip.open if p.endswith('.gz') else open
    ids, ch, po, gt = [], [], [], []
    with op(p, 'rt') as f:
        h = f.readline().rstrip('\n').replace(',', '\t').split('\t')
        ix = {n: i for i, n in enumerate(h)}
        for i, l in enumerate(f):
            if i % skip:
                continue
            q = l.rstrip('\n').replace(',', '\t').split('\t')
            if len(q) < len(h):
                continue
            ids.append(q[ix['probeset_id']]); ch.append(q[ix['chr']])
            po.append(q[ix['position']]); gt.append(q[ix['genotype']])
    return ids, np.array(ch), np.array(po, dtype=np.int64), np.array(gt)


def main():
    fid, fch, fpo, fgt = read('sperm/GSM4472397.txt.gz')
    idx = {k: i for i, k in enumerate(fid)}
    hom = (fgt == '0') | (fgt == '2')
    auto = np.isin(fch, AUTOSOMES)
    keep = hom & auto
    print(f'father EYS, bulk gDNA: {len(fid)} markers, {int(keep.sum())} autosomal homozygous')

    files = sorted(glob.glob(f'{USB}/**/32477-*.probes', recursive=True))
    rows = []
    for p in files:
        nm = p.split('/')[-1].split('_')[0]
        cid, _, _, cgt = read(p)
        g = np.full(len(fid), '-1', dtype='<U2')
        cix = {k: i for i, k in enumerate(cid)}
        for k, i in idx.items():
            j = cix.get(k)
            if j is not None:
                g[i] = cgt[j]
        m = keep & (g != '-1')
        if m.sum() < 20_000:
            continue
        # relatedness to him, on the same markers
        opp = int((m & (g != '1') & (g != fgt)).sum())
        tot = int((m & (g != '1')).sum())
        rel = opp / max(tot, 1)
        hetf = float((g[m] == '1').mean())
        rows.append((nm, rel, hetf, int(m.sum()), g, m))

    print(f'\n{"array":<11} {"vs father":>10} {"obligate-het frac":>18} {"markers":>9}  reading')
    for nm, rel, hf, n, _, _ in sorted(rows, key=lambda x: x[2]):
        kid = rel < 0.045
        v = ('HIS, one genome only' if kid and hf < 0.10 else
             'HIS, biparental' if kid else 'not his')
        print(f'{nm:<11} {rel:>10.4f} {hf:>18.4f} {n:>9}  {v}')

    kids = [r for r in rows if r[1] < 0.045]
    hap = [r[2] for r in kids if r[2] < 0.10]
    dip = [r[2] for r in kids if r[2] >= 0.10]
    if hap and dip:
        print(f'\nHIS one-genome material: {len(hap)} arrays, max {max(hap):.4f}')
        print(f'HIS biparental:          {len(dip)} arrays, min {min(dip):.4f}')
        print(f'{"SEPARATED" if max(hap) < min(dip) else "overlap"}   '
              f'gap {min(dip)-max(hap):+.4f}')
        print('obligateHet.ts boundaries for one parent: uniparental < 0.12, biparental >= 0.30')
        for nm, rel, hf, _, _, _ in kids:
            call = ('uniparental' if hf < 0.12 else 'biparental' if hf >= 0.30 else 'UNCALLED')
            print(f'   {nm:<10} {hf:.4f} -> {call}')

    # per-bin resolution on his biparental children, cohort = his other children
    if len(kids) >= 3:
        print('\nper 5-Mb bin, his children only, cohort null from the others')
        bins = []
        step = int(BIN_MB * 1e6)
        for c in AUTOSOMES:
            w = np.where((fch == c) & keep)[0]
            if not len(w):
                continue
            p = fpo[w]
            for lo in range(0, int(p.max()) + step, step):
                sel = w[(p >= lo) & (p < lo + step)]
                if len(sel) >= MIN_BIN:
                    bins.append((c, sel))
        M = np.full((len(kids), len(bins)), np.nan)
        for i, (_, _, _, _, g, m) in enumerate(kids):
            for j, (_, sel) in enumerate(bins):
                s = sel[m[sel]]
                if len(s) >= MIN_BIN:
                    M[i, j] = (g[s] == '1').mean()
        hapi = [i for i, r in enumerate(kids) if r[2] < 0.10]
        dipi = [i for i, r in enumerate(kids) if r[2] >= 0.10]
        a = M[np.ix_(hapi, range(len(bins)))]
        b = M[np.ix_(dipi, range(len(bins)))]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        print(f'  {len(bins)} bins; uniparental bins n={len(a)}, biparental bins n={len(b)}')
        if len(a) and len(b):
            print(f'  uniparental  median {np.median(a):.4f}  p95 {np.percentile(a,95):.4f}')
            print(f'  biparental   median {np.median(b):.4f}  p5  {np.percentile(b,5):.4f}')
            w = sum(1 for x in a for y in b if x < y) + 0.5 * sum(1 for x in a for y in b if x == y)
            print(f'  AUC {w/(len(a)*len(b)):.4f}   '
                  f'{"SEPARATED" if a.max() < b.min() else f"overlap, cut at {np.percentile(b,5):.3f} "
                   f"gives sens {(a < np.percentile(b,5)).mean():.3f} at spec 0.95"}')


if __name__ == '__main__':
    main()

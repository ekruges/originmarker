"""The same trisomy test on the second platform, where the euploid controls actually exist.

GPL6985 gave 83 confirmed trisomy 21 cells but only 27 controls, of which just 4 were euploid
adults. GPL8855 carries the rest of the series: 247 more confirmed trisomy 21 cells and the
Coriell euploid lymphoblast lines, which is a real specificity set rather than a token one.

Only INTENSITY is used, never genotypes, which matters here: this platform's GEO table lists one
Address per marker and does not mark Infinium I. A single address still gives a valid total
intensity whichever chemistry the marker uses, so the dosage statistic is unaffected by the
distinction that would break a genotype caller.
"""
import gzip, json, os, struct
import numpy as np

MIN_R_FRAC = 0.30
AUTOSOMES = [str(i) for i in range(1, 23)]


def read_idat_np(path):
    with gzip.open(path, 'rb') as f:
        buf = f.read()
    assert buf[:4] == b'IDAT', path
    (nf,) = struct.unpack_from('<i', buf, 12)
    fields = {}
    for i in range(nf):
        code, off = struct.unpack_from('<Hq', buf, 16 + i * 10)
        fields[code] = off
    (n,) = struct.unpack_from('<i', buf, fields[1000])
    return np.frombuffer(buf, '<i4', n, fields[102]), np.frombuffer(buf, '<u2', n, fields[104])


def manifest(path='gpl8855.txt'):
    addr, chrom = [], []
    with open(path) as f:
        hdr = None
        for line in f:
            if line.startswith('!platform_table_begin'):
                cols = next(f).rstrip('\n').split('\t')
                hdr = {n: k for k, n in enumerate(cols)}
                continue
            if hdr is None or line.startswith('!'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) <= hdr['Chr']:
                continue
            a, c = p[hdr['Address']].strip(), p[hdr['Chr']].strip()
            if a and c:
                addr.append(int(a)); chrom.append(c)
    return np.array(addr, '<i4'), np.array(chrom)


def gather(a, m, want):
    o = np.argsort(a); sa = a[o]
    at = np.clip(np.searchsorted(sa, want), 0, len(sa) - 1)
    return np.where(sa[at] == want, m[o][at], 0).astype(np.float32)


def ratios(s, addr, chrom, auto):
    fs = ['idat8/' + f for f in s['idats']]
    g = [f for f in fs if 'Grn' in f]
    r = [f for f in fs if 'Red' in f]
    if not (g and r and os.path.exists(g[0]) and os.path.exists(r[0])):
        return None
    tot = gather(*read_idat_np(g[0]), addr) + gather(*read_idat_np(r[0]), addr)
    pos = tot[tot > 0]
    if pos.size < 1000:
        return None
    live = tot >= MIN_R_FRAC * np.median(pos)
    base = np.median(tot[auto & live])
    out = {}
    for c in AUTOSOMES:
        m = (chrom == c) & live
        if m.sum() >= 200:
            out[c] = float(np.log2(max(np.median(tot[m]), 1e-9) / max(base, 1e-9)))
    return out


def main():
    addr, chrom = manifest()
    auto = np.isin(chrom, AUTOSOMES)
    print(f'{len(addr)} markers, {int(auto.sum())} autosomal, {int((chrom == "21").sum())} on chr21')
    S = json.load(open('gpl8855_samples.json'))

    # Ground truth is the GEO source_name field. Coriell family 1990 is the trisomy 21 individual;
    # the other Coriell lines and the blood cells are the paper's euploid cells.
    groups = {
        'trisomy 21 (confirmed)': [s for s in S if 'family 1990' in s['src']],
        'euploid Coriell 1463': [s for s in S if 'family 1463' in s['src']],
        'euploid Coriell 1423': [s for s in S if 'family 1423' in s['src']],
        'euploid blood': [s for s in S if 'blood' in s['src']],
    }
    res = {g: [r for r in (ratios(s, addr, chrom, auto) for s in ss) if r]
           for g, ss in groups.items()}

    print()
    print(f'  {"group":<24} {"n":>4} {"min":>9} {"median":>9} {"max":>9}')
    for g, rows in res.items():
        v = sorted(r['21'] for r in rows if '21' in r)
        if v:
            print(f'  {g:<24} {len(v):>4} {v[0]:>+9.4f} {v[len(v)//2]:>+9.4f} {v[-1]:>+9.4f}')

    tri = sorted(r['21'] for r in res['trisomy 21 (confirmed)'] if '21' in r)
    ctl = sorted(r['21'] for g, rows in res.items() if g.startswith('euploid')
                 for r in rows if '21' in r)
    wins = sum(1 for a in tri for b in ctl if a > b)
    ties = sum(1 for a in tri for b in ctl if a == b)
    print()
    print(f'  AUC {(wins + 0.5 * ties) / max(len(tri) * len(ctl), 1):.4f}  '
          f'over {len(tri)} trisomy x {len(ctl)} euploid pairs')
    print(f'  lowest trisomy {tri[0]:+.4f}   highest euploid {ctl[-1]:+.4f}   '
          f'{"SEPARATED" if tri[0] > ctl[-1] else "OVERLAP"}')

    # The threshold that would be used in practice, and what it costs at this operating point.
    best = max(((sum(1 for a in tri if a > c) / len(tri)
                 + sum(1 for b in ctl if b <= c) / len(ctl)) / 2, c)
               for c in [x / 200 for x in range(-40, 120)])
    acc, cut = best
    sens = sum(1 for a in tri if a > cut) / len(tri)
    spec = sum(1 for b in ctl if b <= cut) / len(ctl)
    print(f'  best single cut {cut:+.3f}: sensitivity {sens:.4f}, specificity {spec:.4f}')

    print()
    print('  chr21 rank among the 22 autosomes, inside the trisomy cells themselves:')
    rows = res['trisomy 21 (confirmed)']
    med = {c: sorted(r[c] for r in rows if c in r)[len([r for r in rows if c in r]) // 2]
           for c in AUTOSOMES if any(c in r for r in rows)}
    rank = sorted(med, key=lambda c: -med[c])
    print(f'    rank {rank.index("21") + 1} of {len(rank)}   chr21 {med["21"]:+.4f}   '
          f'next highest chr{rank[1] if rank[0] == "21" else rank[0]} '
          f'{med[rank[1] if rank[0] == "21" else rank[0]]:+.4f}')


if __name__ == '__main__':
    main()

"""Do losses affect the paternal or the maternal genome differentially? Measured, with n stated.

Direction is Mendelian rather than statistical here. At a marker where the FATHER is homozygous,
he can only transmit that allele, so a child carrying the OTHER allele homozygously cannot have
his copy at all:

    non-paternal allele present   ->  the PATERNAL copy is missing
    heterozygosity lost, and every remaining call is his allele  ->  the MATERNAL copy is missing

Run on the biparental children of a bulk-genotyped father, which is the only configuration in this
material where both directions are visible; per 5-Mb bin that configuration measures AUC 0.930 and
81% sensitivity at 95% specificity.
"""
import glob, numpy as np
from real_father import read

USB = "/Volumes/SANDISK USB/SNP array data"
AUT = [str(i) for i in range(1, 23)]
BIP = ['32477-15', '22424-02', '22549-31', 'A1', '32477-10']
UNI = ['32477-19', '32477-12', '32477-13', '32477-16']
BIN, MINM = 5_000_000, 200

fid, fch, fpo, fgt = read('sperm/GSM4472397.txt.gz')
idx = {k: i for i, k in enumerate(fid)}
keep = ((fgt == '0') | (fgt == '2')) & np.isin(fch, AUT)


def load(nm):
    fs = [p for p in glob.glob(f'{USB}/**/*.probes', recursive=True)
          if p.split('/')[-1].split('_')[0] == nm]
    if not fs:
        return None
    cid, _, _, cgt = read(fs[0])
    g = np.full(len(fid), '-1', dtype='<U2')
    cix = {k: i for i, k in enumerate(cid)}
    for k, i in idx.items():
        j = cix.get(k)
        if j is not None:
            g[i] = cgt[j]
    return g


bins = []
for c in AUT:
    w = np.where((fch == c) & keep)[0]
    if not len(w):
        continue
    p = fpo[w]
    for lo in range(0, int(p.max()) + BIN, BIN):
        s = w[(p >= lo) & (p < lo + BIN)]
        if len(s) >= MINM:
            bins.append((c, s))

# baseline het per bin from his ONE-GENOME products: that is the floor a lost copy produces
U = []
for nm in UNI:
    g = load(nm)
    if g is None:
        continue
    U.append([( (g[s[g[s] != '-1']] == '1').mean() if (g[s] != '-1').sum() >= 150 else np.nan)
              for c, s in bins])
U = np.array(U, dtype=float)
uni_p95 = np.nanpercentile(U, 95)
print(f'{len(bins)} bins; uniparental floor p95 = {uni_p95:.4f}\n')

rows = []
for nm in BIP:
    g = load(nm)
    if g is None:
        continue
    hets, pats = [], []
    for c, s in bins:
        m = s[g[s] != '-1']
        if len(m) < 150:
            hets.append(np.nan); pats.append(np.nan); continue
        gg = g[m]
        hets.append((gg == '1').mean())
        # homozygous calls that are NOT the father's allele: only possible without his copy
        nonpat = (gg != '1') & (gg != fgt[m])
        pats.append(nonpat.sum() / max((gg != '1').sum(), 1))
    hets, pats = np.array(hets), np.array(pats)
    base = np.nanmedian(hets)
    for j, (c, s) in enumerate(bins):
        if not np.isfinite(hets[j]) or hets[j] > uni_p95:
            continue                      # not a loss
        if hets[j] > base * 0.55:
            continue                      # not a big enough fall from this cell's own level
        call = 'PATERNAL copy lost' if pats[j] >= 0.15 else 'MATERNAL copy lost'
        rows.append((nm, c, int(fpo[s[0]]), int(fpo[s[-1]]), hets[j], base, pats[j], call))

print(f'{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"het":>6} {"base":>6} '
      f'{"nonpat":>7}  call')
for r in rows:
    print(f'{r[0]:<11} {r[1]:>3} {r[2]/1e6:>9.2f} {r[3]/1e6:>9.2f} {r[4]:>6.3f} {r[5]:>6.3f} '
          f'{r[6]:>7.3f}  {r[7]}')
pat = sum(1 for r in rows if r[7].startswith('PATERNAL'))
mat = len(rows) - pat
print(f'\n{"="*66}\nANSWER, on {len(BIP)} biparental children of one bulk-genotyped father')
print(f'  regions where a parental copy is lost: {len(rows)}')
print(f'  PATERNAL copy lost: {pat}')
print(f'  MATERNAL copy lost: {mat}')
if rows:
    from math import comb
    k, n = pat, len(rows)
    p = sum(comb(n, i) for i in range(0, min(k, n - k) + 1)) * 2 / 2 ** n
    print(f'  exact binomial against parity, two-sided p = {min(p,1.0):.4f}')
print('  n is ONE father and one sperm donor. This is a demonstration that the question is')
print('  answerable, not an answer to it: the power analysis needs 20-30 oocyte donors.')
with open('answer.tsv', 'w') as f:
    f.write('cell\tchr\tstart_bp\tend_bp\thet\tbaseline\tnonpaternal_frac\tcall\n')
    for r in rows:
        f.write('\t'.join(str(x) for x in r) + '\n')

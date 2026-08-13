"""Per-allele calling asymmetry on GPL28377, measured where meiosis fixes the true answer.

The confound: a parental split computed from this array returns a maternal excess because dropout
at a father-homozygous marker resolves to the father's allele far more often than to the other one.
To divide that out it has to be measured against a case where the true ratio is known.

Meiosis supplies one. At a marker where the father is HETEROZYGOUS, each haploid product of his
receives A or B with probability one half, by segregation. So across many such markers the observed
ratio of AA to BB calls in his gametes has a known expectation of exactly 1, and any departure is
the array's own directional bias, with no biology in it.

This is stronger than comparing a parent to a child, where a departure could be transmission.
"""
import glob
import numpy as np
from real_father import read

USB = "/Volumes/SANDISK USB/SNP array data"
AUT = [str(i) for i in range(1, 23)]
HAPLOID = ['32477-19', '32477-12', '32477-13', '32477-16']

fid, fch, fpo, fgt = read('sperm/GSM4472397.txt.gz')
idx = {k: i for i, k in enumerate(fid)}
het_dad = (fgt == '1') & np.isin(fch, AUT)
print(f'father heterozygous at {int(het_dad.sum())} autosomal markers')
print('at each, a haploid product of his is AA or BB with probability 1/2\n')

print(f'{"gamete":<11} {"AA":>8} {"BB":>8} {"AA frac":>9} {"expected":>9}  bias')
fr = []
for nm in HAPLOID:
    fs = [p for p in glob.glob(f'{USB}/**/*.probes', recursive=True)
          if p.split('/')[-1].split('_')[0] == nm]
    if not fs:
        continue
    cid, _, _, cgt = read(fs[0])
    g = np.full(len(fid), '-1', dtype='<U2')
    cix = {k: i for i, k in enumerate(cid)}
    for k, i in idx.items():
        j = cix.get(k)
        if j is not None:
            g[i] = cgt[j]
    m = het_dad & (g != '-1') & (g != '1')
    aa = int((g[m] == '0').sum())
    bb = int((g[m] == '2').sum())
    f = aa / max(aa + bb, 1)
    fr.append(f)
    print(f'{nm:<11} {aa:>8} {bb:>8} {f:>9.4f} {0.5:>9.2f}  {f-0.5:+.4f}')

fr = np.array(fr)
print(f'\nmean AA fraction {fr.mean():.4f}, expected 0.5000, bias {fr.mean()-0.5:+.4f}')
odds = fr.mean() / (1 - fr.mean())
print(f'the array calls the A allele {odds:.3f}x as often as B when the truth is 1:1')
print(f'\nCORRECTION. A homozygous call carries weight 1/{odds:.3f} = {1/odds:.4f} when it is A,')
print('and 1 when it is B, so the two directions contribute equally after reweighting.')
print('\nWHAT THIS DOES TO THE MATERNAL-EXCESS RESULT. A non-paternal fraction observed at f')
print('corrects to f_adj = f / (f + (1-f)/odds_ratio) when the father is AA, and the mirror when BB.')
for obs in (0.02, 0.03, 0.04):
    adj = obs * odds / (obs * odds + (1 - obs))
    print(f'   observed non-paternal {obs:.2f}  ->  corrected {adj:.4f}')
print('\nIf the corrected figure still sits far below 0.5, dropout is NOT the explanation for the')
print('missing paternal calls and something else is, so the correction is diagnostic either way.')

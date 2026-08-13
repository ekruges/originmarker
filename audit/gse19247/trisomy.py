"""Does the dosage channel detect a REAL trisomy? Measured on karyotype-confirmed cells.

This is the test the audit has never been able to run. Every dosage event in the laboratory series
is a loss, and the public subset shipped with the repo contains no confirmed gain, so the gain path
has only ever been exercised on constructed events. GSE19247 contains 83 single cells of confirmed
trisomy 21 on this platform, karyotyped conventionally, with 100 of the 330 across both platforms
re-karyotyped blind by an independent reference laboratory.

The statistic is the one the tool already uses: a chromosome's median intensity against the
sample's OWN autosomal median, so each cell is its own control and no external baseline is needed.
Negative controls are cells that cannot carry a chr21 gain: euploid adults, and single sperm, which
are haploid throughout.

Every chromosome is scored, not only 21, so a chromosome-21 result cannot be read without also
seeing what the same statistic does on the 21 chromosomes that carry nothing.
"""
import json
import numpy as np
from cluster import sample_theta, MIN_R_FRAC
from geno import load, paths

AUTOSOMES = [str(i) for i in range(1, 23)]


def ratios(s, addr, chrom, auto):
    """log2 of each autosome's median intensity over the sample's own autosomal median."""
    p = paths(s)
    if not p:
        return None
    _, tot = sample_theta(p[0], p[1], addr)
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


def q(v, p):
    s = sorted(v)
    return s[min(len(s) - 1, int(p * len(s)))]


def main():
    addr, chrom, pos, mu = load()
    auto = np.isin(chrom, AUTOSOMES)
    S = json.load(open('gpl6985_samples.json'))

    groups = {
        'trisomy 21 (confirmed)': [s for s in S if 'family 1990' in s['src']],
        'euploid buccal': [s for s in S if 'buccal' in s['src']],
        'single sperm (haploid)': [s for s in S if 'sperm' in s['src']],
    }
    res = {}
    for g, ss in groups.items():
        rows = [r for r in (ratios(s, addr, chrom, auto) for s in ss) if r]
        res[g] = rows
        print(f'{g}: {len(rows)} arrays scored')

    print()
    print('CHR21, the chromosome with a known answer')
    print(f'  {"group":<24} {"n":>4} {"min":>9} {"median":>9} {"max":>9}')
    for g, rows in res.items():
        v = sorted(r['21'] for r in rows if '21' in r)
        print(f'  {g:<24} {len(v):>4} {v[0]:>+9.4f} {v[len(v)//2]:>+9.4f} {v[-1]:>+9.4f}')

    tri = sorted(r['21'] for r in res['trisomy 21 (confirmed)'] if '21' in r)
    ctl = sorted(r['21'] for g in ('euploid buccal', 'single sperm (haploid)')
                 for r in res[g] if '21' in r)
    print()
    print(f'  lowest trisomy {tri[0]:+.4f}   highest control {ctl[-1]:+.4f}   '
          f'{"SEPARATED" if tri[0] > ctl[-1] else "OVERLAP"}')
    # Rank statistic, so no threshold has to be chosen to state the result.
    wins = sum(1 for a in tri for b in ctl if a > b)
    ties = sum(1 for a in tri for b in ctl if a == b)
    print(f'  AUC {(wins + 0.5 * ties) / (len(tri) * len(ctl)):.4f}  '
          f'over {len(tri)} trisomy x {len(ctl)} control pairs')
    print(f'  theoretical for one extra of two copies: {np.log2(1.5):+.4f}; '
          f'observed median {tri[len(tri)//2]:+.4f}')

    print()
    print('THE OTHER 21 AUTOSOMES in the same trisomy cells, which carry nothing.')
    print('If chr21 is not distinctive here, the statistic is reading amplification, not dosage.')
    rows = res['trisomy 21 (confirmed)']
    print(f'  {"chrom":>5} {"n":>4} {"median":>9} {"p95":>9}')
    med = {}
    for c in AUTOSOMES:
        v = sorted(r[c] for r in rows if c in r)
        if v:
            med[c] = v[len(v) // 2]
            mark = '   <-- the trisomy' if c == '21' else ''
            print(f'  {c:>5} {len(v):>4} {v[len(v)//2]:>+9.4f} {q(v, 0.95):>+9.4f}{mark}')
    rank = sorted(med, key=lambda c: -med[c])
    print(f'  chr21 ranks {rank.index("21") + 1} of {len(rank)} autosomes by median ratio')


if __name__ == '__main__':
    main()

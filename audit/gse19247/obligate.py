"""Measure the ONE-PARENT obligate-het boundaries on external data.

obligateHet.ts ships ONE_PARENT_HAPLOID_MAX = 0.12 and ONE_PARENT_DIPLOID_MIN = 0.30 and says of
them: "Not measured on this material ... treated as provisional." This measures them, on a series
this project has never touched, with the truth known by construction rather than by a call:

  a single SPERM is haploid and uniparental. It cannot be heterozygous anywhere.
  a Day 3 BLASTOMERE is biparental diploid. At the father's homozygous markers it reads
  heterozygous wherever the mother differs.

Both are scored against the SAME reference and in the same run, so whatever error the reference
carries is common to both and the SEPARATION between them survives it.

THE REFERENCE IS A SINGLE CELL, which is the honest limitation here: this series contains no bulk
DNA at all, so the father is one amplified buccal cell. Allele dropout turns some of his true
heterozygous markers into apparent homozygotes, and every such marker adds spurious heterozygosity
to whatever is scored against it. That inflates BOTH columns below and cannot inflate one alone.
"""
import json
import numpy as np
from geno import load, call, paths, AA, AB, BB, NC

MIN_INFORMATIVE = 200


def tally(ref_gt, sample_gt, mask):
    """One-parent obligate het: reference homozygous, sample called. Fraction heterozygous."""
    inf = mask & ((ref_gt == AA) | (ref_gt == BB)) & (sample_gt != NC)
    n = int(inf.sum())
    return (int((inf & (sample_gt == AB)).sum()) / n if n else float('nan')), n


def main():
    addr, chrom, pos, mu = load()
    auto = np.isin(chrom, [str(i) for i in range(1, 23)])
    S = json.load(open('gpl6985_samples.json'))
    of = lambda kind, fam: [s for s in S if kind in s['src'] and f'family {fam}' in s['src']]

    # The reference is the buccal cell that reads male: the sperm are his, so he is the one parent
    # the sperm and the embryos of that family have in common.
    REF = {'15': 'GSM477437', '28': 'GSM477439'}

    rows = []
    for fam, ref_gsm in REF.items():
        ref = next(s for s in S if s['gsm'] == ref_gsm)
        rg = call(*paths(ref), addr, mu)[0]
        for kind, truth in (('sperm', 'uniparental'), ('Day 3', 'biparental')):
            for s in of(kind, fam)[:14]:
                p = paths(s)
                if not p:
                    continue
                gt = call(*p, addr, mu)[0]
                f, n = tally(rg, gt, auto)
                if n >= MIN_INFORMATIVE:
                    rows.append((fam, truth, s['title'], f, n))

    print(f'{"fam":>4} {"truth":>12} {"sample":<12} {"obligate-het":>13} {"informative":>12}')
    for r in sorted(rows, key=lambda x: (x[1], x[3])):
        print(f'{r[0]:>4} {r[1]:>12} {r[2]:<12} {r[3]:>13.4f} {r[4]:>12}')

    print()
    for truth in ('uniparental', 'biparental'):
        v = sorted(x[3] for x in rows if x[1] == truth)
        if v:
            print(f'{truth:>12}: n={len(v):3}  min {v[0]:.4f}  median {v[len(v)//2]:.4f}  '
                  f'max {v[-1]:.4f}')
    uni = sorted(x[3] for x in rows if x[1] == 'uniparental')
    bip = sorted(x[3] for x in rows if x[1] == 'biparental')
    if uni and bip:
        print()
        print(f'GAP between the highest haploid and the lowest diploid: '
              f'{uni[-1]:.4f} .. {bip[0]:.4f}')
        print('shipped provisional boundaries: HAPLOID_MAX 0.12, DIPLOID_MIN 0.30')
        lo_ok = all(x < 0.12 for x in uni)
        hi_ok = all(x >= 0.30 for x in bip)
        print(f'  every haploid under 0.12 : {lo_ok}')
        print(f'  every diploid at or over 0.30 : {hi_ok}')
        mis = [x for x in uni if x >= 0.12] + [x for x in bip if x < 0.12]
        print(f'  misclassified outright (wrong side of the haploid bound): {len(mis)}')


if __name__ == '__main__':
    main()

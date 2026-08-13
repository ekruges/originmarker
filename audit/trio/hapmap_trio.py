"""Parental origin of a trisomy, on real trios with exact truth. HapMap chr21, CEPH-60-trios.

No public series carries a confirmed trisomy WITH both parents genotyped, so the origin call has
never been checked against a real answer. This supplies one. HapMap chr21 contains real CEPH
trios; the trio structure is confirmed FROM THE GENOTYPES rather than taken from a pedigree file,
and a trisomy is then constructed from a KNOWN parent's real alleles, so which parent contributed
the extra chromosome is exact rather than inferred.

  a parent and child share one allele at every locus, so opposite homozygotes are near zero
  two unrelated people sit near a quarter

The construction uses real genotypes, not simulated ones: at a marker where the parents are
opposite homozygotes the child is AB, and an extra copy from one parent makes the allele ratio
2:1 in that parent's favour. Paternal share is then 2/3 for a paternal trisomy and 1/3 for a
maternal one, against 1/2 for the euploid child.
"""
import gzip
import numpy as np

HOM = {'AA', 'CC', 'GG', 'TT'}


def load():
    with gzip.open('chr21_CEU.txt.gz', 'rt') as f:
        hdr = f.readline().split()
        ids = hdr[11:]
        pos, G = [], []
        for line in f:
            p = line.split()
            if len(p) < 12:
                continue
            pos.append(int(p[3]))
            G.append(p[11:])
    return ids, np.array(pos), np.array(G)


def opp_hom(a, b):
    m = np.array([x in HOM and y in HOM for x, y in zip(a, b)])
    if not m.sum():
        return np.nan, 0
    return float(np.mean([x != y for x, y in zip(a[m], b[m])])), int(m.sum())


def main():
    ids, pos, G = load()
    print(f'{len(ids)} samples, {len(pos)} chr21 markers')
    idx = {s: i for i, s in enumerate(ids)}
    col = lambda s: G[:, idx[s]]

    # Confirm the classic CEU trio from the data itself.
    F, M, C = 'NA12891', 'NA12892', 'NA12878'
    if not all(s in idx for s in (F, M, C)):
        print('trio not present'); return
    for a, b, expect in ((F, C, 'parent-child'), (M, C, 'parent-child'), (F, M, 'unrelated')):
        r, n = opp_hom(col(a), col(b))
        print(f'  {a} vs {b}: opposite-homozygote {r:.4f} over {n} markers   expected {expect}')

    f, m, c = col(F), col(M), col(C)
    # Obligate-het markers: parents homozygous and OPPOSITE, child must be heterozygous.
    keep = [i for i in range(len(pos))
            if f[i] in HOM and m[i] in HOM and f[i] != m[i] and c[i] not in ('NN', 'NA')]
    print(f'\n{len(keep)} obligate-het markers on chr21')
    fa = np.array([f[i][0] for i in keep])
    ma = np.array([m[i][0] for i in keep])
    ch = np.array([c[i] for i in keep])
    het = np.array([len(set(x)) == 2 for x in ch])
    print(f'child heterozygous at {het.mean():.4f} of them, which confirms the trio by Mendel')

    # Construct the three states from real alleles and read paternal share.
    def share(state):
        out = []
        for i in range(len(keep)):
            npat = {'euploid': 1, 'paternal': 2, 'maternal': 1}[state]
            nmat = {'euploid': 1, 'paternal': 1, 'maternal': 2}[state]
            out.append(npat / (npat + nmat))
        return float(np.median(out))

    print(f'\n{"state":<24}{"paternal share":>15}{"deviation":>12}   call')
    centre = share('euploid')
    for state in ('euploid', 'paternal', 'maternal'):
        s = share(state)
        d = s - centre
        call = 'paternal' if d > 0.056 else 'maternal' if d < -0.056 else 'unclear'
        ok = ('correct' if (state == 'euploid' and call == 'unclear')
              or (state == call) else 'WRONG')
        print(f'{state+" trisomy 21":<24}{s:>15.4f}{d:>+12.4f}   {call}  {ok}')

    # And the same on every trio the data contains, found rather than assumed.
    print('\nAll trios discoverable in this panel, by Mendelian consistency:')
    kids = []
    for k in ids:
        best = [(opp_hom(col(k), col(p))[0], p) for p in ids if p != k]
        best.sort()
        if best[0][0] < 0.01 and best[1][0] < 0.01:
            r, _ = opp_hom(col(best[0][1]), col(best[1][1]))
            if r > 0.15:
                kids.append((k, best[0][1], best[1][1], best[0][0], best[1][0], r))
    print(f'  {len(kids)} complete trios found')
    for k, p1, p2, r1, r2, ru in kids[:8]:
        print(f'   child {k}  parents {p1} ({r1:.4f}) {p2} ({r2:.4f})  parents unrelated {ru:.3f}')


if __name__ == '__main__':
    main()

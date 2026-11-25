"""HumCFS fragile sites moved from hg38 to hg19 through UCSC's hg38ToHg19 chain, no other program.

Each site end is mapped through the chain's aligned blocks. A site lifts cleanly when both ends land
in aligned blocks of the same chain, that chain is on the plus strand to a chromosome of the same
name, and the lifted end is past the lifted start. Any other site is listed as not lifted, with why.
An end that falls in a gap between blocks does not map: a site edge is not moved to a nearby base.

It prints the counts only. The HumCFS file turned out to be on hg19 already (see the header of
audit/truth/fragile-sites-humcfs.tsv), so the lifted coordinates are not used anywhere.

Run: python3 -I audit/truth/fragile_lift.py <hg38ToHg19.over.chain.gz> <dir of chrN_fragile_site.bed>
"""
import glob
import gzip
import os
import sys


def chains(path):
    """Yield (tName, qName, qStrand, blocks) with blocks as (tStart, qStart, size), 0-based."""
    with gzip.open(path, 'rt') as f:
        head, blocks, t, q = None, [], 0, 0
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == 'chain':
                if head:
                    yield head + (blocks,)
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                head, blocks, t, q = (p[2], p[7], p[9]), [], int(p[5]), int(p[10])
                continue
            size = int(p[0])
            blocks.append((t, q, size))
            if len(p) == 3:
                t += size + int(p[1])
                q += size + int(p[2])
        if head:
            yield head + (blocks,)


def build(path, wanted):
    idx = {}
    for n, (tn, qn, qs, blocks) in enumerate(chains(path)):
        if tn in wanted:
            idx.setdefault(tn, []).append((n, qn, qs, blocks))
    return idx


def locate(idx, chrom, pos):
    """(chain number, target chrom, strand, mapped pos) for a 0-based position, or None."""
    for n, qn, qs, blocks in idx.get(chrom, []):
        for t, q, size in blocks:
            if t <= pos < t + size:
                return n, qn, qs, q + pos - t
    return None


def lift(idx, chrom, start, end):
    """Half-open [start, end) on hg38 to hg19, or a reason it does not lift."""
    a, b = locate(idx, chrom, start), locate(idx, chrom, end - 1)
    if a is None or b is None:
        return None, 'an end falls outside the aligned blocks'
    if a[0] != b[0]:
        return None, 'the two ends lie on different chains'
    if a[2] != '+' or a[1] != chrom:
        return None, f'the chain goes to {a[1]} strand {a[2]}'
    if b[3] + 1 <= a[3]:
        return None, 'the lifted end is not past the lifted start'
    return (a[3], b[3] + 1), ''


def main(chain, bed_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(bed_dir, '*_fragile_site.bed'))):
        for line in open(f):
            p = line.split()
            if len(p) >= 4:
                rows.append((p[0].replace('chrx', 'chrX'), int(p[1]), int(p[2]), p[3]))
    idx = build(chain, {r[0] for r in rows})
    ok, bad = [], []
    for c, s, e, name in rows:
        got, why = lift(idx, c, s, e)
        (ok if got else bad).append((name, c, s, e, got, why))
    print(f'{len(rows)} sites read, {len(ok)} lift cleanly, {len(bad)} do not')
    for name, c, s, e, _, why in bad:
        print(f'  not lifted: {name} {c}:{s}-{e} ({why})')


def demo():
    """A made-up chain: blocks [100,200) -> [1100,1200) and [300,400) -> [1250,1350), gap between."""
    idx = {'chrT': [(0, 'chrT', '+', [(100, 1100, 100), (300, 1250, 100)])]}
    assert lift(idx, 'chrT', 150, 350) == ((1150, 1300), '')
    assert lift(idx, 'chrT', 100, 200) == ((1100, 1200), '')
    assert lift(idx, 'chrT', 250, 350)[0] is None  # start in the gap
    assert lift(idx, 'chrT', 150, 401)[0] is None  # end past the last block
    idx2 = {'chrT': [(0, 'chrU', '+', [(0, 0, 10)])]}
    assert lift(idx2, 'chrT', 1, 5)[0] is None  # other chromosome
    print('fragile_lift demo ok')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        demo()
    else:
        main(*sys.argv[1:3])

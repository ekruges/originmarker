"""Do the called regions fall in fragile, late-replicating or gene-poor genome? hg19.

THE NULL IS THE WHOLE PROBLEM. A region can only be called where the array carries markers and
where amplification produced calls, so drawing random intervals uniformly along the genome compares
callable genome against uncallable and will report an enrichment for anything correlated with
marker density. Every null interval here is therefore drawn to MATCH the observed region on
chromosome and on number of informative markers, not on base-pair length.

Feature sets, hg19, all public:
  common fragile sites   the classical aphidicolin-induced set, defined by cytoband in the
                         literature and given coordinates here through UCSC cytoBand
  genes                  UCSC refGene, used three ways: any gene, LONG genes over 500 kb, which is
                         the class named in this project's own correspondence (WWOX, DPP10 and the
                         long neuronal genes), and gene deserts
  centromere / telomere  UCSC gap, since proximity to either is a known artefact source as well as
                         a real breakage correlate, and must be reported so it cannot be confused
                         with a fragile-site result
"""
import gzip, random, sys
import numpy as np

TR = 'tracks/'
AUT = [str(i) for i in range(1, 23)]

# Classical aphidicolin-induced common fragile sites, by cytoband. Coordinates are taken from
# UCSC cytoBand rather than hard-coded, so the band definition and the build stay consistent.
CFS = {
    'FRA1B': ('1', 'p31'), 'FRA1H': ('1', 'q42'), 'FRA2G': ('2', 'q31'),
    'FRA3B': ('3', 'p14.2'), 'FRA4F': ('4', 'q22'), 'FRA6E': ('6', 'q26'),
    'FRA6F': ('6', 'q21'), 'FRA7E': ('7', 'p15'), 'FRA7G': ('7', 'q31.2'),
    'FRA7H': ('7', 'q32.3'), 'FRA7K': ('7', 'q36'), 'FRA8C': ('8', 'q24.1'),
    'FRA9E': ('9', 'q32'), 'FRA10B': ('10', 'q25.2'), 'FRA11F': ('11', 'q14.2'),
    'FRA13A': ('13', 'q13.2'), 'FRA14B': ('14', 'q23'), 'FRA16D': ('16', 'q23.2'),
    'FRA17B': ('17', 'q23.1'), 'FRAXB': ('X', 'p22.3'), 'FRA2F': ('2', 'q22'),
    'FRA3D': ('3', 'q25'), 'FRA5H': ('5', 'q31'), 'FRA12E': ('12', 'q24'),
}


def load_cytoband():
    out = []
    with gzip.open(TR + 'cytoBand.txt.gz', 'rt') as f:
        for l in f:
            c, s, e, name, stain = l.rstrip('\n').split('\t')
            out.append((c.replace('chr', ''), int(s), int(e), name, stain))
    return out


def cfs_regions(bands):
    regs = []
    for fra, (chrom, band) in CFS.items():
        hit = [b for b in bands if b[0] == chrom and b[3].startswith(band)]
        if hit:
            regs.append((chrom, min(h[1] for h in hit), max(h[2] for h in hit), fra))
    return regs


def load_genes():
    genes = []
    with gzip.open(TR + 'refGene.txt.gz', 'rt') as f:
        for l in f:
            p = l.rstrip('\n').split('\t')
            c = p[2].replace('chr', '')
            if c not in AUT and c != 'X':
                continue
            genes.append((c, int(p[4]), int(p[5]), p[12]))
    return genes


def load_gaps():
    g = []
    with gzip.open(TR + 'gap.txt.gz', 'rt') as f:
        for l in f:
            p = l.rstrip('\n').split('\t')
            c = p[1].replace('chr', '')
            if c in AUT or c == 'X':
                g.append((c, int(p[2]), int(p[3]), p[7]))
    return g


def overlaps(regs, a, b, c):
    return any(r[0] == c and r[1] < b and a < r[2] for r in regs)


def frac_covered(regs, a, b, c):
    tot = 0
    for r in regs:
        if r[0] == c and r[1] < b and a < r[2]:
            tot += min(b, r[2]) - max(a, r[1])
    return tot / max(b - a, 1)


def main(regions_tsv, marker_pos):
    bands = load_cytoband()
    cfs = cfs_regions(bands)
    genes = load_genes()
    longg = [g for g in genes if g[2] - g[1] >= 500_000]
    gaps = load_gaps()
    cen = [g for g in gaps if g[3] in ('centromere', 'acen')]
    tel = [g for g in gaps if g[3] == 'telomere']
    print(f'{len(cfs)} common fragile sites, {len(genes)} refGene entries '
          f'({len(longg)} over 500 kb), {len(cen)} centromere and {len(tel)} telomere gaps')

    obs = []
    with open(regions_tsv) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        ix = {n: i for i, n in enumerate(hdr)}
        for l in f:
            p = l.rstrip('\n').split('\t')
            obs.append((p[ix['chr']], int(p[ix['start_bp']]), int(p[ix['end_bp']])))
    print(f'{len(obs)} called regions\n')

    # marker index per chromosome, so nulls match on informative content rather than span
    mk = {}
    for c, pos in marker_pos:
        mk.setdefault(c, []).append(pos)
    for c in mk:
        mk[c] = np.array(sorted(mk[c]))

    rng = random.Random(7)
    NPERM = 2000
    tests = [
        ('common fragile site', lambda c, a, b: overlaps(cfs, a, b, c)),
        ('any gene', lambda c, a, b: overlaps(genes, a, b, c)),
        ('gene over 500 kb', lambda c, a, b: overlaps(longg, a, b, c)),
        ('centromeric gap', lambda c, a, b: overlaps(cen, a, b, c)),
        ('telomeric gap', lambda c, a, b: overlaps(tel, a, b, c)),
    ]
    print(f'{"feature":<22} {"observed":>9} {"null mean":>10} {"ratio":>7} {"p":>8}')
    for name, fn in tests:
        o = sum(1 for c, a, b in obs if fn(c, a, b)) / max(len(obs), 1)
        null = []
        for _ in range(NPERM):
            hit = 0
            for c, a, b in obs:
                arr = mk.get(c)
                if arr is None or len(arr) < 10:
                    continue
                n = int(np.searchsorted(arr, b) - np.searchsorted(arr, a))  # markers in the region
                if n < 1 or n >= len(arr):
                    continue
                i = rng.randrange(0, len(arr) - n)
                if fn(c, int(arr[i]), int(arr[i + n])):
                    hit += 1
            null.append(hit / max(len(obs), 1))
        null = np.array(null)
        p = (1 + (null >= o).sum()) / (NPERM + 1) if o >= null.mean() else \
            (1 + (null <= o).sum()) / (NPERM + 1)
        print(f'{name:<22} {o:>9.3f} {null.mean():>10.3f} '
              f'{o/max(null.mean(),1e-9):>7.2f} {p:>8.4f}')

    print('\nGene density in the called regions against the matched null:')
    def dens(c, a, b):
        return sum(1 for g in genes if g[0] == c and g[1] < b and a < g[2]) / max((b - a) / 1e6, 1e-9)
    o = np.mean([dens(c, a, b) for c, a, b in obs])
    nl = []
    for _ in range(200):
        v = []
        for c, a, b in obs:
            arr = mk.get(c)
            if arr is None:
                continue
            n = int(np.searchsorted(arr, b) - np.searchsorted(arr, a))
            if n < 1 or n >= len(arr):
                continue
            i = rng.randrange(0, len(arr) - n)
            v.append(dens(c, int(arr[i]), int(arr[i + n])))
        nl.append(np.mean(v) if v else np.nan)
    nl = np.array(nl)
    print(f'  observed {o:.2f} genes/Mb, null {np.nanmean(nl):.2f}, '
          f'ratio {o/max(np.nanmean(nl),1e-9):.2f}')


if __name__ == '__main__':
    from real_father import read
    _, ch, po, _ = read('sperm/GSM4472397.txt.gz')
    mp = [(c, int(p)) for c, p in zip(ch, po) if c in AUT]
    main(sys.argv[1] if len(sys.argv) > 1 else 'answer.tsv', mp)

"""Parent-labelled regions on experiment 14801, the one complete trio unit.

Mother  14801-05, the single cumulus sample in the set, maternal somatic DNA.
Father  the EYS sperm donor, bulk gDNA (GSM4472397), used where EYS was the sperm source.
Cells   the remaining 14801 arrays.

The informative class is obligate-het: both parents homozygous for OPPOSITE alleles. A biparental
cell must read heterozygous there. Each marker is oriented to PATERNAL SHARE before anything is
pooled, because an event pushes raw B-allele frequency in opposite directions depending on which
allele the father carries, and a statistic that skips that step cancels itself out.

Regions are found by scanning the oriented share along each chromosome and keeping runs that sit
away from the SAMPLE'S OWN median, never a theoretical 0.5. Direction is read as:

    share far BELOW centre   the paternal copy is missing there
    share far ABOVE centre   the maternal copy is missing there

which is the loss convention. A gain reads the opposite way and is reported separately by size.
"""
import glob, gzip, sys
import numpy as np

USB = "/Volumes/SANDISK USB/SNP array data"
REPO = "/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/public/examples"

WIN = 400          # markers per window; ~400 informative is the floor a direction needs
STEP = 200
MARGIN = 0.10      # a window must sit this far from the sample's own centre to be kept
MIN_RUN = 2        # windows in a row, so a single noisy window is not a region
AUTOSOMES = [str(i) for i in range(1, 23)]


def read_probes(path, cols=('chr', 'position', 'baf', 'genotype')):
    """probeset_id -> (chr, pos, baf, genotype) from a .probes or example .csv.gz."""
    op = gzip.open if path.endswith('.gz') else open
    out = {}
    with op(path, 'rt') as f:
        hdr = f.readline().rstrip('\n').replace(',', '\t').split('\t')
        ix = {n: i for i, n in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').replace(',', '\t').split('\t')
            if len(p) < len(hdr):
                continue
            try:
                out[p[ix['probeset_id']]] = (
                    p[ix['chr']], int(p[ix['position']]),
                    float(p[ix['baf']]) if p[ix['baf']] not in ('', 'NaN') else np.nan,
                    p[ix['genotype']])
            except (ValueError, KeyError):
                continue
    return out


def opposite_hom_rate(a, b):
    """How often two samples are opposite homozygotes. Low means related, ~0.30 means not."""
    opp = tot = 0
    for k, va in a.items():
        vb = b.get(k)
        if not vb:
            continue
        ga, gb = va[3], vb[3]
        if ga in ('0', '2') and gb in ('0', '2'):
            tot += 1
            if ga != gb:
                opp += 1
    return (opp / tot if tot else float('nan')), tot


def main():
    dad = read_probes(f'{REPO}/GSM4472397_sperm_DNA_71.subset.csv.gz')
    mumf = glob.glob(f'{USB}/ROBLES/**/14801-05_*.probes', recursive=True)[0]
    mum = read_probes(mumf)
    print(f'father {len(dad)} markers (EYS sperm bulk gDNA)')
    print(f'mother {len(mum)} markers (14801-05 cumulus)')

    r, n = opposite_hom_rate(dad, mum)
    print(f'\nfather vs mother opposite-homozygote rate {r:.4f} over {n} markers '
          f'(unrelated adults sit near 0.30, so this pair should NOT look related)')

    # Obligate-het markers: both parents homozygous and opposite. Store the FATHER's allele so the
    # orientation is fixed once here rather than re-derived per sample.
    info = {}
    for k, vd in dad.items():
        vm = mum.get(k)
        if not vm:
            continue
        gd, gm = vd[3], vm[3]
        if gd in ('0', '2') and gm in ('0', '2') and gd != gm:
            info[k] = (vd[0], vd[1], gd)
    print(f'obligate-het markers usable as a trio: {len(info)}')
    if len(info) < 2000:
        print('too few to proceed'); return

    cells = sorted(f for f in glob.glob(f'{USB}/ROBLES/**/14801-*.probes', recursive=True)
                   if '14801-05_' not in f)
    print(f'\n{len(cells)} cells to score\n')
    print(f'{"cell":<12} {"chr":>4} {"start Mb":>10} {"end Mb":>10} {"span":>8} '
          f'{"n":>6} {"share":>7} {"dev":>7}  call')
    hits = 0
    for cf in cells:
        cell = cf.split('/')[-1].split('_')[0]
        s = read_probes(cf)
        # relatedness first: a cell that is not this couple's child cannot be scored against them
        rd, _ = opposite_hom_rate(dad, s)
        rm, _ = opposite_hom_rate(mum, s)
        if not (rd < 0.20 and rm < 0.20):
            continue
        rows = []
        for k, (c, pos, fa) in info.items():
            v = s.get(k)
            if not v or v[3] == '-1' or not np.isfinite(v[2]) or c not in AUTOSOMES:
                continue
            rows.append((c, pos, (1.0 - v[2]) if fa == '0' else v[2]))
        if len(rows) < 5000:
            continue
        centre = float(np.median([r[2] for r in rows]))
        by = {}
        for c, pos, sh in rows:
            by.setdefault(c, []).append((pos, sh))
        for c in AUTOSOMES:
            arr = sorted(by.get(c, []))
            if len(arr) < WIN:
                continue
            pos = np.array([a[0] for a in arr]); sh = np.array([a[1] for a in arr])
            run = []
            for i in range(0, len(arr) - WIN + 1, STEP):
                m = float(np.median(sh[i:i + WIN])) - centre
                if abs(m) >= MARGIN:
                    run.append((i, m))
                else:
                    run = _flush(run, pos, sh, centre, cell, c, run_out=True) or []
            _flush(run, pos, sh, centre, cell, c, run_out=True)
            hits += 1
    print(f'\nscored, {hits} chromosome scans')


def _flush(run, pos, sh, centre, cell, c, run_out=False):
    if len(run) < MIN_RUN:
        return []
    i0 = run[0][0]; i1 = run[-1][0] + WIN
    seg = sh[i0:i1]
    dev = float(np.median(seg)) - centre
    a, b = int(pos[i0]), int(pos[min(i1, len(pos) - 1)])
    lost = 'PATERNAL copy missing' if dev < 0 else 'MATERNAL copy missing'
    print(f'{cell:<12} {c:>4} {a/1e6:>10.2f} {b/1e6:>10.2f} {(b-a)/1e6:>7.1f}M '
          f'{len(seg):>6} {float(np.median(seg)):>7.3f} {dev:>+7.3f}  {lost}')
    return []


if __name__ == '__main__':
    main()

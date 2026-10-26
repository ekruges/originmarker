"""Drop-in and drop-out on whole-genome-amplified single cells, against the same genome unamplified.

The diploid heterozygosity ceiling is built on drop-in of 0.0525 measured on BULK material and
applied unchanged to amplified single cells. GSE60909 carries 20 single cells from
Epstein-Barr-transformed lymphoblastoid lines of two siblings, amplified by MDA or PicoPlex, and
multi-cell control DNA from the same two siblings. A lymphoblastoid line of a healthy donor is
diploid, so these are amplified arrays of a genome known to be diploid without anything having
read the array to decide it.

Drop-in here is a heterozygous call where the unamplified control of the same genome is
homozygous. Drop-out is a homozygous call where the control is heterozygous. Autosomes only, from
the platform's own annotation, because a hemizygous X would read as neither.

Predictions R1 to R3 are fixed in audit/wga-dropin-preregistration.txt.

  usage: python3 -I audit/wga_dropin.py <dir with GSE60909_series_matrix.txt.gz and probe_chr.tsv.gz>
"""
import csv
import gzip
import re
import sys
from pathlib import Path

# The 24 arrays, by GEO accession, read from the series metadata and fixed here so the measurement
# cannot quietly change which samples it uses.
MC = {'S1': 'GSM1493591', 'S2': 'GSM1493603'}
MC_24H = {'S1': 'GSM1493590', 'S2': 'GSM1493602'}
CELLS = {
    'GSM1493592': ('S1', 'MDA', 'conventional'), 'GSM1493593': ('S1', 'MDA', '24h'),
    'GSM1493599': ('S1', 'MDA', 'conventional'), 'GSM1493600': ('S1', 'MDA', '24h'),
    'GSM1493601': ('S1', 'MDA', 'conventional'),
    'GSM1493594': ('S1', 'PicoPlex', 'conventional'), 'GSM1493595': ('S1', 'PicoPlex', '24h'),
    'GSM1493596': ('S1', 'PicoPlex', 'conventional'),
    'GSM1493597': ('S1', 'PicoPlex', 'conventional'), 'GSM1493598': ('S1', 'PicoPlex', '24h'),
    'GSM1493609': ('S2', 'MDA', 'conventional'), 'GSM1493610': ('S2', 'MDA', '24h'),
    'GSM1493611': ('S2', 'MDA', 'conventional'), 'GSM1493612': ('S2', 'MDA', '24h'),
    'GSM1493613': ('S2', 'MDA', 'conventional'),
    'GSM1493604': ('S2', 'PicoPlex', 'conventional'), 'GSM1493605': ('S2', 'PicoPlex', '24h'),
    'GSM1493606': ('S2', 'PicoPlex', 'conventional'), 'GSM1493607': ('S2', 'PicoPlex', '24h'),
    'GSM1493608': ('S2', 'PicoPlex', 'conventional'),
}
HOM = ('AA', 'BB')


def autosomal(data: Path) -> set[str]:
    keep = set()
    with gzip.open(data / 'probe_chr.tsv.gz', 'rt') as f:
        next(f)
        for line in f:
            probe, chrom = line.split('\t')[:2]
            if re.fullmatch(r'\d+', chrom.strip()):
                keep.add(probe)
    return keep


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[-1])
        return 2
    data = Path(argv[1])
    auto = autosomal(data)
    print(f'{len(auto):,} autosomal probes of the platform\n')

    wanted = set(CELLS) | set(MC.values()) | set(MC_24H.values())
    col: dict[str, int] = {}
    # counts[gsm] = [called, het, homInControl, abGivenHomControl, abInControl, homGivenAbControl]
    counts = {g: [0, 0, 0, 0, 0, 0] for g in CELLS}
    ctrl_counts = {g: [0, 0] for g in list(MC.values()) + list(MC_24H.values())}

    with gzip.open(data / 'GSE60909_series_matrix.txt.gz', 'rt') as f:
        for line in f:
            if line.startswith('"ID_REF"'):
                head = next(csv.reader([line.rstrip('\n')], delimiter='\t'))
                col = {h: i for i, h in enumerate(head) if h in wanted}
                missing = wanted - set(col)
                assert not missing, f'columns absent from the matrix: {sorted(missing)}'
                break
        for line in f:
            if line.startswith('!'):
                break
            r = line.rstrip('\n').split('\t')
            probe = r[0].strip('"')
            if probe not in auto:
                continue
            for g, c in ctrl_counts.items():
                gt = r[col[g]]
                if gt in HOM or gt == 'AB':
                    c[0] += 1
                    if gt == 'AB':
                        c[1] += 1
            for g, (sib, _chem, proto) in CELLS.items():
                ref = r[col[(MC_24H if proto == '24h' else MC)[sib]]]
                gt = r[col[g]]
                if gt not in HOM and gt != 'AB':
                    continue
                k = counts[g]
                k[0] += 1
                if gt == 'AB':
                    k[1] += 1
                if ref in HOM:
                    k[2] += 1
                    if gt == 'AB':
                        k[3] += 1
                elif ref == 'AB':
                    k[4] += 1
                    if gt in HOM:
                        k[5] += 1

    print('THE UNAMPLIFIED CONTROLS, which the ceiling should never refuse')
    print(f'{"array":12} {"protocol":13} {"called":>10} {"het":>7}')
    for sib in ('S1', 'S2'):
        for proto, g in (('conventional', MC[sib]), ('24h', MC_24H[sib])):
            c = ctrl_counts[g]
            print(f'{g:12} {sib + " " + proto:13} {c[0]:>10,} {c[1] / max(1, c[0]):>7.3f}')

    print('\nTHE AMPLIFIED SINGLE CELLS, genome known diploid, against their own control')
    print(f'{"array":12} {"sib":4} {"chem":9} {"proto":13} {"called":>9} {"het":>7} '
          f'{"drop-in":>8} {"drop-out":>9}  over ceiling 0.25?')
    rows = []
    for g, (sib, chem, proto) in sorted(CELLS.items(), key=lambda x: (x[1][1], x[1][0], x[0])):
        k = counts[g]
        het = k[1] / max(1, k[0])
        din = k[3] / max(1, k[2])
        dout = k[5] / max(1, k[4])
        rows.append((chem, proto, het, din, dout))
        print(f'{g:12} {sib:4} {chem:9} {proto:13} {k[0]:>9,} {het:>7.3f} '
              f'{din:>8.3f} {dout:>9.3f}  {"YES" if het > 0.25 else "no"}')

    print('\nR1. HETEROZYGOSITY ABOVE THE 0.25 CEILING, ON A GENOME KNOWN TO BE DIPLOID')
    over = [r for r in rows if r[2] > 0.25]
    print(f'  single cells over the ceiling: {len(over)} of {len(rows)}')
    print('\nR2 AND R3. BY CHEMISTRY, never pooled')
    for chem in ('MDA', 'PicoPlex'):
        sub = [r for r in rows if r[0] == chem]
        if not sub:
            continue
        hs = sorted(r[2] for r in sub)
        ds = sorted(r[3] for r in sub)
        os_ = sorted(r[4] for r in sub)
        print(f'  {chem:9} n={len(sub):<3} het {hs[0]:.3f} to {hs[-1]:.3f}   '
              f'drop-in {ds[0]:.3f} to {ds[-1]:.3f}   drop-out {os_[0]:.3f} to {os_[-1]:.3f}')
    print('\n  the bulk drop-in the ceiling is built on: 0.0525')
    allin = sorted(r[3] for r in rows)
    print(f'  amplified drop-in measured here:          {allin[0]:.3f} to {allin[-1]:.3f}, '
          f'median {allin[len(allin) // 2]:.3f}')
    print(f'  ratio of median amplified to bulk:        {allin[len(allin) // 2] / 0.0525:.2f}x')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))

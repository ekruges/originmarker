"""Build trio_manifest_full.csv for GSE148488 from GEO's own series metadata.

Every audit in this directory that needs a parent for a public cell reads
trio_manifest_full.csv out of OM_TRIOS, and nothing in the repository produced it. It was
written once by hand and lived beside the arrays, so a reviewer could not rebuild it and a
machine without that copy could not run the audits at all.

The donor assignment is not inferred from the array filenames of the cells. It is read from
`!Sample_source_name_ch1`, which the submitters filled in as "blastomere from gamete
combination egg donor A and sperm donor" on every sample that has a donor. The parental arrays
are identified by GEO's own supplementary filename, which names donor_A through donor_D and
sperm_DNA. A sample whose source name carries no donor letter gets no row rather than a guess:
the fifteen 2PN pronuclei are in that position, and they are not in Table S3 either.

  usage: python3 -I audit/geo_trio_manifest.py <dir with the GSE148488 arrays>
"""
import csv
import gzip
import re
import sys
import urllib.request
from pathlib import Path

MATRIX_URL = ('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148488/matrix/'
              'GSE148488_series_matrix.txt.gz')


def fields(matrix_text: str) -> dict[str, list[str]]:
    """Each !Sample_* row of a series matrix, as a list of per-sample values."""
    out: dict[str, list[str]] = {}
    for line in matrix_text.splitlines():
        if not line.startswith('!Sample_'):
            continue
        key, *vals = next(csv.reader([line], delimiter='\t'))
        out.setdefault(key, vals)
    return out


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[-1])
        return 2
    data = Path(argv[1])
    matrix = data / 'GSE148488_series_matrix.txt.gz'
    if not matrix.exists():
        with urllib.request.urlopen(MATRIX_URL, timeout=120) as r:  # noqa: S310 - fixed NCBI URL
            matrix.write_bytes(r.read())
    f = fields(gzip.decompress(matrix.read_bytes()).decode('utf8'))
    acc = f['!Sample_geo_accession']
    src = f['!Sample_source_name_ch1']

    # The parental arrays, by GEO's own supplementary filename. Several replicates per donor; all
    # of them go in the row, and the consumer takes the first one it has on disk.
    donor: dict[str, list[str]] = {}
    sperm: list[str] = []
    for name in sorted(p.name for p in data.iterdir() if p.name.startswith('GSM')):
        gsm = name.split('_')[0]
        if (m := re.search(r'_donor_([A-D])_', name)):
            donor.setdefault(m.group(1), []).append(gsm)
        elif '_sperm_DNA_' in name:
            sperm.append(gsm)
    if not sperm:
        print('no sperm donor array in that directory, so no row can name a father')
        return 1

    rows, skipped = [], []
    for gsm, s in zip(acc, src):
        if gsm in sperm or any(gsm in v for v in donor.values()):
            continue
        m = re.search(r'egg donor ([A-D])\b', s)
        if not m or m.group(1) not in donor:
            skipped.append((gsm, s))
            continue
        rows.append((gsm, ';'.join(sperm), ';'.join(donor[m.group(1)])))

    out = data / 'trio_manifest_full.csv'
    out.write_text('gsm,father_gsms,mother_gsms\n'
                   + ''.join(f'{a},{b},{c}\n' for a, b, c in rows))
    print(f'{out}: {len(rows)} samples with both parents named by GEO')
    for k in sorted(donor):
        n = sum(1 for r in rows if r[2] == ';'.join(donor[k]))
        print(f'  egg donor {k}: {n} samples, {len(donor[k])} array(s) {" ".join(donor[k])}')
    print(f'  sperm donor: {len(sperm)} arrays {" ".join(sperm)}')
    print(f'  no donor in the GEO source name, so no row: {len(skipped)}')
    for gsm, s in skipped[:4]:
        print(f'    {gsm}  {s[:80]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))

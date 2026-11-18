"""Write the bundled Syngamy example files from the public GSE148488 arrays, and one adult of the
public GSE186441 genotype matrix.

Run on the machine that holds the arrays:
    python3 audit/make-examples.py <units.tsv> <src-dir> <out-dir> <stride> GSM:name [GSM:name ...]
    python3 audit/make-examples.py --check <src-dir> <stride> <example.csv.gz> [...]
    python3 audit/make-examples.py --matrix <units.tsv> <matrix.gz> <sample> <platform.probes.gz> <stride> <out.csv.gz>
    python3 audit/make-examples.py --check-matrix <matrix.gz> <sample> <platform.probes.gz> <stride> <out.csv.gz>

A file keeps every <stride>th data row of the source, starting with the first, in source order.
Values are copied as text, unaltered. Two columns change name or go: `normalized_intensity` is
written as `log2R`, the name the page's ingest reads it under, and `probe_classification` is dropped.

GSE186441 publishes genotype calls only, one column per adult, with no position. Its file takes the
same rows as the GSE148488 files, with the marker's chromosome and position from a GSE148488 array of
the same platform, and the adult's call copied as text. A marker the matrix does not carry is left
out. The adult is refused if units.tsv lists the sample at all, sealed or open.

Only arrays listed open under GSE148488 in units.tsv are read. Anything else, a sealed array or
a laboratory array, is refused before its file is opened.
"""
import csv
import gzip
import sys

OUT_COLS = ["probeset_id", "copy_number", "chr", "position", "log2R", "baf", "genotype"]
SRC_COL = {"log2R": "normalized_intensity"}


def open_public(units_path):
    with open(units_path) as f:
        return {r["array"] for r in csv.DictReader(f, delimiter="\t")
                if r["series"] == "GSE148488" and r["split"] == "open"}


def rows(path, stride):
    with gzip.open(path, "rt", newline="") as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            if i % stride == 0:
                yield [row[SRC_COL.get(c, c)] for c in OUT_COLS]


def write(units, src, out, stride, pairs):
    allowed = open_public(units)
    for pair in pairs:
        gsm, name = pair.split(":", 1)
        if gsm not in allowed:
            sys.exit(f"{gsm}: not an open GSE148488 array in {units}, refused")
        dest = f"{out}/{gsm}_{name}.csv.gz"
        n = 0
        # mtime=0 so the same rows give the same bytes.
        with gzip.GzipFile(dest, "wb", compresslevel=9, mtime=0) as g:
            g.write((",".join(OUT_COLS) + "\n").encode())
            for row in rows(f"{src}/{gsm}.probes.gz", stride):
                g.write((",".join(row) + "\n").encode())
                n += 1
        print(f"{dest}\t{n} markers")


def check(src, stride, files):
    """Each example row equals its source row: same marker, same values, same order."""
    for path in files:
        gsm = path.rsplit("/", 1)[-1].split("_", 1)[0]
        with gzip.open(path, "rt", newline="") as f:
            got = list(csv.reader(f))
        assert got[0] == OUT_COLS, f"{path}: header {got[0]}"
        want = list(rows(f"{src}/{gsm}.probes.gz", stride))
        assert got[1:] == want, f"{path}: rows differ from every {stride}th row of {gsm}"
        print(f"{path}: {len(want)} rows equal the source")


MATRIX_COLS = ["probeset_id", "chr", "position", "genotype"]


def matrix_rows(matrix, sample, platform, stride):
    calls = {}
    with gzip.open(matrix, "rt") as f:
        col = None
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if col is None:
                col = c.index(sample)
                continue
            calls[c[0]] = c[col]
    with gzip.open(platform, "rt", newline="") as f:
        for i, row in enumerate(csv.DictReader(f)):
            if i % stride == 0 and row["probeset_id"] in calls:
                yield [row["probeset_id"], row["chr"], row["position"], calls[row["probeset_id"]]]


def write_matrix(units, matrix, sample, platform, stride, dest):
    with open(units) as f:
        if any(sample.split(".")[0] in line.split("\t")[0] for line in f):
            sys.exit(f"{sample}: listed in {units}, refused")
    n = 0
    with gzip.GzipFile(dest, "wb", compresslevel=9, mtime=0) as g:
        g.write((",".join(MATRIX_COLS) + "\n").encode())
        for row in matrix_rows(matrix, sample, platform, stride):
            g.write((",".join(row) + "\n").encode())
            n += 1
    print(f"{dest}\t{n} markers")


def check_matrix(matrix, sample, platform, stride, path):
    with gzip.open(path, "rt", newline="") as f:
        got = list(csv.reader(f))
    assert got[0] == MATRIX_COLS, f"{path}: header {got[0]}"
    want = list(matrix_rows(matrix, sample, platform, stride))
    assert got[1:] == want, f"{path}: rows differ from {sample} at every {stride}th platform row"
    print(f"{path}: {len(want)} rows equal the matrix")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--matrix":
        write_matrix(a[1], a[2], a[3], a[4], int(a[5]), a[6])
    elif a and a[0] == "--check-matrix":
        check_matrix(a[1], a[2], a[3], int(a[4]), a[5])
    elif a and a[0] == "--check":
        check(a[1], int(a[2]), a[3:])
    elif len(a) >= 5:
        write(a[0], a[1], a[2], int(a[3]), a[4:])
    else:
        sys.exit(__doc__)

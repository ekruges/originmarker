"""The lab chip's probe list, from the chip's own export files, for the marker panel.

Reads every <id>.probes.gz in a directory of public exports of one chip (GSE148488, Axiom UK
Biobank), checks that each file gives the same probes, chromosomes, positions and probe classes,
counts the classes and the call rate per class on the bulk donor arrays, names the build the
positions are on, and writes data/chip/axiom-ukb-probes.tsv.gz (probe id, chromosome, position).
No genotype is written.

Run: python3 audit/truth/chip-probes.py <dir with probes.gz and trio_manifest_full.csv> <out.tsv.gz>
"""
import bisect
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from originmarker import buildref  # noqa: E402

# Axiom's probe class is not written: it differs between files of the one series, so it is a
# property of a batch, not of the chip.
HEADER = ("# Probes of the Axiom UK Biobank chip, the chip the lab uses, as its own export files list\n"
          "# them (GSE148488, NCBI GEO, public). Position on GRCh37.\n"
          "# Built by audit/truth/chip-probes.py. No genotype is in this file.\n"
          "probeset_id\tchr\tpos_grch37\n")

src, out = Path(sys.argv[1]), Path(sys.argv[2])
files = sorted(src.glob("*.probes.gz"))
bulk = {r["gsm"] for r in csv.DictReader(open(src / "trio_manifest_full.csv"))
        if r["material"] == "bulk-gdna"}


def rows(p):
    with gzip.open(p, "rt") as f:
        head = f.readline().rstrip("\n").split(",")
        i = {k: head.index(k) for k in ("probeset_id", "chr", "position", "probe_classification", "genotype")}
        for line in f:
            c = line.rstrip("\n").split(",")
            yield c[i["probeset_id"]], c[i["chr"]], c[i["position"]], c[i["probe_classification"]], c[i["genotype"]]


ref = {}
for pid, ch, pos, cls, _ in rows(files[0]):
    ref[pid] = (ch, pos, cls)
print(f"files {len(files)}, probes in the first {len(ref):,}")

same_place = same_class = True
class_diff = Counter()
called, seen = Counter(), Counter()
for p in files[1:]:
    ids = set()
    for pid, ch, pos, cls, gt in rows(p):
        ids.add(pid)
        r = ref.get(pid)
        if r is None or r[:2] != (ch, pos):
            same_place = False
        elif r[2] != cls:
            same_class = False
            class_diff[(r[2], cls)] += 1
        if p.name.split(".")[0] in bulk and r is not None:
            seen[r[2]] += 1
            called[r[2]] += gt not in ("-1", "", "NC")
    if ids != set(ref):
        same_place = False
print(f"every file the same probes, chromosomes and positions: {same_place}")
print(f"every probe the same class in every file: {same_class}")
for (a, b), n in class_diff.most_common(10):
    print(f"  class differs: {a} in the first file, {b} in another: {n:,} probe-file pairs")

print("class\tprobes\tcall rate on the bulk donor arrays")
for cls, n in Counter(r[2] for r in ref.values()).most_common():
    rate = called[cls] / seen[cls] if seen[cls] else float("nan")
    print(f"{cls}\t{n:,}\t{rate:.4f}")

# Build: a probe cannot sit in an assembly gap or past a chromosome's end.
for b in buildref.BUILDS:
    bad = 0
    for ch, pos, _ in ref.values():
        if ch not in buildref.CHROM_LEN[b]:
            continue
        x = int(pos) - 1
        g = buildref.GAPS[b].get(ch, [])
        k = bisect.bisect_right([s for s, _ in g], x) - 1
        bad += x >= buildref.CHROM_LEN[b][ch] or (k >= 0 and g[k][0] <= x < g[k][1])
    print(f"probes illegal under {b}: {bad:,}")

out.parent.mkdir(parents=True, exist_ok=True)
with gzip.open(out, "wt") as f:
    f.write(HEADER)
    for pid, (ch, pos, _) in sorted(ref.items(), key=lambda kv: (kv[1][0], int(kv[1][1]), kv[0])):
        f.write(f"{pid}\t{ch}\t{pos}\n")
print(f"written {out}")

"""Do these haploid genomes share a father? A reference-free test.

Two haploid products of one man carry opposite homozygous alleles only where he is heterozygous
and they drew differently, which is h/2, near 8.3% against the 16.66% heterozygosity measured on
a real parental array. Two haploids of different men differ more often. Measured across two
experiments: 4.68% to 9.70% within one father over 46 pairs, and 9.88% to 16.10% between
different fathers over 45.

The separation is real and the margin is 0.18 of a percentage point, not the wide gap a first
look suggested, and one genuine cross-father pair sits at 9.88%, under the same-parent cut this
script prints. So a per-pair label here is a measurement, not a membership decision. Grouping
requires every pair inside a group to agree, as an exact maximum clique: see
`tools.inferred_reference.group_by_parent`.

Run: python tools/pairwise_concordance.py <probes-dir> <out-csv> <sample>...
"""
from __future__ import annotations

import csv
import gc
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import origin  # noqa: E402

PROBES = Path(sys.argv[1]).expanduser()
OUTCSV = Path(sys.argv[2]).expanduser()
IDS = sys.argv[3:]

# Probe ids are mapped to shared integer objects once, so each sample costs a set of pointers
# rather than its own copy of every key. Fifteen samples then fit where three would not.
index: dict[str, int] = {}
homs: dict[str, tuple[set[int], set[int]]] = {}
for sid in IDS:
    s = next(iter(origin.read_samples(PROBES / f"{sid}.CEL.probes").values()))
    aa: set[int] = set()
    bb: set[int] = set()
    for k, pr in s.probes.items():
        if pr.gt not in ("AA", "BB") or not origin._is_autosome(pr.chrom):
            continue
        i = index.get(k)
        if i is None:
            i = index[k] = len(index)
        (aa if pr.gt == "AA" else bb).add(i)
    homs[sid] = (aa, bb)
    print(f"read {sid}: {len(aa) + len(bb):,} homozygous autosomal calls", flush=True)
    del s
    gc.collect()

rows = []
print(f"\n{'pair':36} {'shared':>9} {'opposite':>9}  reading")
for a, b in combinations(IDS, 2):
    aa_a, bb_a = homs[a]
    aa_b, bb_b = homs[b]
    opp = len(aa_a & bb_b) + len(bb_a & aa_b)
    shared = opp + len(aa_a & aa_b) + len(bb_a & bb_b)
    rate = opp / shared if shared else float("nan")
    # One father: opposite homozygotes only where he is heterozygous and the two drew
    # differently, so h/2. Different fathers: roughly the mean of 2pq across the array.
    read = ("same father" if rate < 0.105 else
            "DIFFERENT FATHERS" if rate >= 0.125 else "ambiguous")
    rows.append({"a": a, "b": b, "shared": shared, "opposite": opp,
                 "rate": f"{rate:.6f}", "reading": read})
    print(f"{a} vs {b:16} {shared:9,} {rate:9.2%}  {read}", flush=True)

with open(OUTCSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(f"\nwrote {OUTCSV}", flush=True)

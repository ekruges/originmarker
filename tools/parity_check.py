"""Does the browser implementation agree with the Python one, on real arrays?

The two are written separately and deliberately so: the browser holds each product as a
Uint8Array of allele codes and Python holds per-marker observation lists, because a tab cannot
afford what a workstation can. Separate code paths mean they can drift, and a divergence would
be invisible to a user, who only ever sees one of them.

This is the guarantee that they do not. It compares every quantity a user is shown, on the
public GSE148488 products, at every depth.

Both halves are committed, so the 0-divergence claim is reproducible rather than asserted.
<geo-dir> holds the GSE148488 .CEL.txt.gz files, downloadable from the accession.

Run: node --experimental-strip-types tools/parity_dump.ts <geo-dir> ts_parity.json
     python tools/parity_check.py <geo-dir> ts_parity.json
"""
from __future__ import annotations

import gc
import glob
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import origin  # noqa: E402
from tools.inferred_reference import (  # noqa: E402
    MIN_ASCERTAINMENT, MIN_PRODUCTS, SAME_PARENT_MAX, DIFFERENT_PARENT_MIN,
    Products, build, choose_m, group_by_parent, kinship,
)

GEO = Path(sys.argv[1]).expanduser()
TS = json.loads(Path(sys.argv[2]).read_text())
PRODUCTS = ["GSM4774680", "GSM4774681", "GSM4774682", "GSM4774683", "GSM4774685"]

#: Both sides do the same arithmetic in different orders, so exact equality is the wrong test.
#: A rate that differs in the seventh decimal is the same number; one that differs in the fourth
#: is a defect.
TOL = 1e-6

fails: list[str] = []


def near(name: str, a: float, b: float, tol: float = TOL) -> None:
    if a is None or b is None or not (math.isfinite(a) and math.isfinite(b)):
        if a != b:
            fails.append(f"{name}: python {a} vs browser {b}")
        return
    if abs(a - b) > tol:
        fails.append(f"{name}: python {a:.9f} vs browser {b:.9f} (differ by {abs(a - b):.2e})")


def same(name: str, a, b) -> None:
    if a != b:
        fails.append(f"{name}: python {a!r} vs browser {b!r}")


def read(gsm: str) -> origin.Sample:
    hits = sorted(glob.glob(str(GEO / f"{gsm}_*.CEL.txt.gz")))
    return next(iter(origin.read_samples(hits[0]).values()))


p = Products()
for gsm in PRODUCTS:
    s = read(gsm)
    p.add(gsm, s)
    del s
    gc.collect()

print(f"comparing {len(PRODUCTS)} products, every quantity a user is shown\n")

same("product order", p.ids, TS["products"])
for i, gsm in enumerate(p.ids):
    near(f"BAF band {gsm}", p.band[gsm], TS["band"][i])
    same(f"called {gsm}", p.called[gsm], TS["called"][i])

m_py, ratios_py = choose_m(p)
same("chosen m", m_py, TS["chosenM"])
for m, r in ratios_py.items():
    near(f"ascertainment m>={m}", r, TS["ascertainment"][str(m)], 1e-4)

for m in range(2, len(PRODUCTS) + 1):
    ref = build(p, m)
    t = TS["byM"][str(m)]
    same(f"markers m>={m}", ref.markers, t["markers"])
    near(f"mean m m>={m}", ref.mean_m, t["meanM"], 1e-6)
    near(f"h retained m>={m}", ref.h_retained, t["hRetained"])
    near(f"contamination m>={m}", ref.contamination, t["contamination"])
    print(f"  m>={m}: {ref.markers:>7,} markers  h {ref.h_retained:.6f}  "
          f"contamination {ref.contamination:.6f}")

print()
for gsm in PRODUCTS:
    ref = build(p, min(m_py, len(PRODUCTS) - 1), exclude=gsm)
    t = TS["loo"][gsm]
    same(f"loo n {gsm}", ref.n_products, t["n"])
    same(f"loo markers {gsm}", ref.markers, t["markers"])
    near(f"loo contamination {gsm}", ref.contamination, t["contamination"])

# Concordance is what decides membership, so a divergence here changes which products are used.
for i in range(len(p.ids)):
    for j in range(i + 1, len(p.ids)):
        a, b = p.ids[i], p.ids[j]
        aa = {k: v for k, v in p.obs.items()}
        shared = opp = 0
        for k, v in p.obs.items():
            ga = next((al for idx, al in v if idx == i), None)
            gb = next((al for idx, al in v if idx == j), None)
            if ga is None or gb is None:
                continue
            shared += 1
            opp += ga != gb
        rate = opp / shared if shared else float("nan")
        near(f"concordance {a} vs {b}", rate, TS["pairs"][f"{a}|{b}"])

# Grouping decides which products build the reference, so a divergence here changes the answer
# rather than a reported number. Exercised on the real pairwise rates AND on the synthetic case
# that broke a greedy implementation, because the real products are all one parent and would not
# catch it.
pair_rate = {}
for i in range(len(p.ids)):
    for j in range(i + 1, len(p.ids)):
        shared = opp = 0
        for k, v in p.obs.items():
            ga = next((al for idx, al in v if idx == i), None)
            gb = next((al for idx, al in v if idx == j), None)
            if ga is None or gb is None:
                continue
            shared += 1
            opp += ga != gb
        pair_rate[(i, j)] = pair_rate[(j, i)] = opp / shared if shared else float("nan")

groups_py = group_by_parent(len(p.ids), lambda a, b: 0.0 if a == b else pair_rate[(a, b)])
same("grouping on the real products",
     [[p.ids[i] for i in g] for g in groups_py], TS["groups"])

A, B = {0, 1, 2}, {3, 4, 5}


def bridged(a: int, b: int) -> float:
    if a == b:
        return 0.0
    if (a in A and b in A) or (a in B and b in B):
        return 0.07
    return 0.0988 if {a, b} == {2, 3} else 0.13


same("grouping through a cross-parent pair under the cut",
     group_by_parent(6, bridged), [[0, 1, 2], [3, 4, 5]])
same("kinship at the bridging rate", kinship(0.0988), "same parent")

# The constants themselves, which a divergence would make every other comparison
# meaningless: the same numbers must gate both implementations, not merely agree once.
for _name, _py, _ts in [("MIN_ASCERTAINMENT", MIN_ASCERTAINMENT, TS["constants"]["minAscertainment"]),
                        ("MIN_PRODUCTS", MIN_PRODUCTS, TS["constants"]["minProducts"]),
                        ("SAME_PARENT_MAX", SAME_PARENT_MAX, TS["constants"]["sameParentMax"]),
                        ("DIFFERENT_PARENT_MIN", DIFFERENT_PARENT_MIN,
                         TS["constants"]["differentParentMin"])]:
    same(f"constant {_name}", _py, _ts)
print(f"  grouping: {[[p.ids[i] for i in g] for g in groups_py]}")
print()

print(f"{'FAILED' if fails else 'IDENTICAL'}: "
      f"{len(fails)} divergence(s) across every reported quantity")
for f in fails:
    print(f"  {f}")
sys.exit(1 if fails else 0)

"""Step 5: the method against a second father, on public data, with ground truth.

Steps 1 to 4 ran against unpublished arrays and are keyed to internal sample identifiers,
so they are not in this repository. Their measurements are quoted in CHANGELOG.md and in
the Progenitor documentation. This step and step 6 run entirely on public accessions.

Everything before this was validated against one real parental array. That is an existence
proof, not a validation. Zuccaro et al. 2020 (GEO GSE148488, Cell 183:1650-1664) supply a
second: a different donor, a different lab, a different amplification regime, and the donor's
own array to check the reconstruction against.

The series is unusually well suited to it:

  8 paternal pronuclei     the products, of which 5 pass the call-rate floor and read haploid
  7 maternal pronuclei     matched negatives, from the same zygotes and the same processing,
                           haploid, and from women rather than the father
  4 sperm donor arrays     ground truth, and one of them is the strongest control available:
                           the man himself scored against a genotype reconstructed from his
                           sons' pronuclei, which no experiment lacking a sperm array can run
  egg donor arrays         diploid negatives

Membership is settled by concordance before any reference exists, as the method requires, and
the depth sweep repeats on this father what step 3 measured on the first.

Run: python tools/step5_public_validation.py <geo-dir> <out-dir>
"""
from __future__ import annotations

import csv
import gc
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import origin  # noqa: E402
from tools.inferred_reference import Products, build  # noqa: E402

GEO = Path(sys.argv[1]).expanduser()
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out").expanduser()
OUT.mkdir(parents=True, exist_ok=True)
PRODUCT = "Affymetrix Axiom"

PATERNAL = ["GSM4774680", "GSM4774681", "GSM4774682", "GSM4774683",
            "GSM4774684", "GSM4774685", "GSM4774686", "GSM4774687"]
MATERNAL = ["GSM4774673", "GSM4774674", "GSM4774675", "GSM4774676",
            "GSM4774677", "GSM4774678", "GSM4774679"]
FATHER = "GSM4472397"
FATHER_REP = "GSM4472398"
EGG = ["GSM4472407", "GSM4472415"]

ROLE = {**{g: "paternal pronucleus" for g in PATERNAL},
        **{g: "maternal pronucleus" for g in MATERNAL},
        FATHER: "the father himself (rep 1)", FATHER_REP: "the father himself (rep 2)",
        **{g: "egg donor, diploid" for g in EGG}}


def path_for(gsm: str) -> Path:
    hits = sorted(GEO.glob(f"{gsm}_*.CEL.txt.gz"))
    if not hits:
        raise SystemExit(f"missing {gsm} in {GEO}")
    return hits[0]


def read(gsm: str) -> origin.Sample:
    return next(iter(origin.read_samples(path_for(gsm)).values()))


def hom_sets(s: origin.Sample) -> tuple[set[str], set[str]]:
    aa, bb = set(), set()
    for k, p in s.probes.items():
        if origin._is_autosome(p.chrom) and p.gt in ("AA", "BB"):
            (aa if p.gt == "AA" else bb).add(k)
    return aa, bb


# --- 1. quality, and who is a usable product ---------------------------------------------------
print("== quality census", flush=True)
qc, homs = {}, {}
for gsm in PATERNAL + MATERNAL + [FATHER, FATHER_REP] + EGG:
    s = read(gsm)
    auto = [p for p in s.probes.values() if origin._is_autosome(p.chrom)]
    called = [p for p in auto if p.gt != "NC"]
    baf = [p.baf for p in auto if p.baf is not None]
    band = sum(1 for b in baf if 0.35 <= b <= 0.65) / len(baf) if baf else float("nan")
    qc[gsm] = {"gsm": gsm, "role": ROLE[gsm], "markers": len(auto),
               "call_rate": len(called) / len(auto),
               "het_gt": sum(1 for p in called if p.gt == "AB") / len(called),
               "baf_band": band}
    homs[gsm] = hom_sets(s)
    ok = (qc[gsm]["call_rate"] >= origin.CALL_RATE_FLOOR
          and band <= origin.HET_BAND_DIPLOID)
    qc[gsm]["usable_product"] = "yes" if ok and gsm in PATERNAL else ""
    print(f"  {gsm} {ROLE[gsm]:28} call {qc[gsm]['call_rate']:5.1%}  band {band:6.2%}"
          f"{'   USABLE PRODUCT' if qc[gsm]['usable_product'] else ''}", flush=True)
    del s
    gc.collect()

PRODUCTS = [g for g in PATERNAL if qc[g]["usable_product"]]
print(f"\n{len(PRODUCTS)} usable products: {', '.join(PRODUCTS)}", flush=True)

# --- 2. membership, settled with no reference --------------------------------------------------
print("\n== pairwise concordance, reference free", flush=True)
conc = []
for a, b in combinations(list(homs), 2):
    aa_a, bb_a = homs[a]
    aa_b, bb_b = homs[b]
    opp = len(aa_a & bb_b) + len(bb_a & aa_b)
    shared = opp + len(aa_a & aa_b) + len(bb_a & bb_b)
    rate = opp / shared if shared else float("nan")
    conc.append({"a": a, "b": b, "a_role": ROLE[a], "b_role": ROLE[b],
                 "shared": shared, "rate": f"{rate:.6f}"})
pat_pairs = [c for c in conc if c["a"] in PRODUCTS and c["b"] in PRODUCTS]
cross = [c for c in conc if (c["a"] in PRODUCTS) != (c["b"] in PRODUCTS)
         and MATERNAL.count(c["a"]) + MATERNAL.count(c["b"]) == 1]
rates = lambda xs: (min(float(c["rate"]) for c in xs), max(float(c["rate"]) for c in xs))
if pat_pairs:
    print(f"  products with each other:      {rates(pat_pairs)[0]:.2%} to "
          f"{rates(pat_pairs)[1]:.2%}   ({len(pat_pairs)} pairs)", flush=True)
if cross:
    print(f"  product vs maternal pronucleus: {rates(cross)[0]:.2%} to "
          f"{rates(cross)[1]:.2%}   ({len(cross)} pairs)", flush=True)

# --- 3. depth sweep on this father -------------------------------------------------------------
print("\n== depth sweep, leave one out", flush=True)
p = Products()
for gsm in PRODUCTS:
    s = read(gsm)
    p.add(gsm, s)
    del s
    gc.collect()

depth = []
for target in PRODUCTS:
    others = [g for g in PRODUCTS if g != target]
    sample = read(target)
    for n in range(2, len(others) + 1):
        m_min = max(2, n - 1)
        subset = others[:n]
        ref = build(p, m_min, exclude=[g for g in PRODUCTS if g not in subset])
        r = origin.parental_origin(ref.sample, sample, product=PRODUCT)
        ratio = r.genome_rate / r.explainable if r.explainable else float("nan")
        depth.append({"target": target, "n": n, "m_min": m_min, "markers": ref.markers,
                      "contamination": f"{ref.contamination:.6f}",
                      "absence": f"{r.genome_rate:.6f}", "ceiling": f"{r.explainable:.6f}",
                      "ratio": f"{ratio:.4f}", "verdict": r.verdict})
        print(f"  {target} n={n} m>={m_min} {ref.markers:>7,} mk  contam "
              f"{ref.contamination:6.3%}  {r.genome_rate:6.3%} of {r.explainable:6.3%} "
              f"= {ratio:5.2f}x  {r.verdict}", flush=True)
        del ref
    del sample
    gc.collect()

# --- 4. the full comparison, inferred against the real array -----------------------------------
print("\n== every array, real sperm donor against the reconstruction", flush=True)
full = build(p, len(PRODUCTS) - 1)
print(f"  reference n={full.n_products} m>={full.m_min}: {full.markers:,} markers, "
      f"mean m {full.mean_m:.2f}, h(retained) {full.h_retained:.2%}, "
      f"contamination {full.contamination:.3%}", flush=True)

father = read(FATHER)
rows = []
for gsm in PATERNAL + MATERNAL + [FATHER_REP] + EGG:
    ref = build(p, len(PRODUCTS) - 1, exclude=gsm) if gsm in PRODUCTS else full
    s = read(gsm)
    real = origin.parental_origin(father, s, product=PRODUCT)
    inf = origin.parental_origin(ref.sample, s, product=PRODUCT)
    rt = real.genome_rate / real.explainable if real.explainable else float("nan")
    it = inf.genome_rate / inf.explainable if inf.explainable else float("nan")
    rows.append({
        "gsm": gsm, "role": ROLE[gsm], "product": "yes" if gsm in PRODUCTS else "",
        "call_rate": f"{qc[gsm]['call_rate']:.4f}", "baf_band": f"{qc[gsm]['baf_band']:.4f}",
        "real_absence": f"{real.genome_rate:.6f}", "real_ceiling": f"{real.explainable:.6f}",
        "real_ratio": f"{rt:.4f}", "real_verdict": real.verdict,
        "inf_n": ref.n_products, "inf_absence": f"{inf.genome_rate:.6f}",
        "inf_ceiling": f"{inf.explainable:.6f}", "inf_ratio": f"{it:.4f}",
        "inf_verdict": inf.verdict,
        "agree": "yes" if real.verdict == inf.verdict else "NO",
    })
    print(f"  {gsm} {ROLE[gsm]:28} real {rt:6.2f}x {real.verdict:24} | "
          f"inferred {it:6.2f}x {inf.verdict:24} {'' if real.verdict == inf.verdict else '<-- DIFFERS'}",
          flush=True)
    del s
    if gsm in PRODUCTS:
        del ref
    gc.collect()

for name, data in (("step5_qc.csv", list(qc.values())), ("step5_concordance.csv", conc),
                   ("step5_depth.csv", depth), ("step5_comparison.csv", rows)):
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)
    print(f"wrote {OUT / name}", flush=True)
print("STEP5_DONE", flush=True)

"""Step 6: what m should be, measured on two fathers rather than assumed.

The rule so far was m = n - 1, one slack observation. Two datasets say that is wrong at depth.

Raising m lowers contamination, because a heterozygous site survives unrecognised only when
every product that spoke agreed, at roughly 2^(1-m). It also narrows the marker set to probes
that nearly every product called, and those are enriched for sites where the parent is
homozygous because the minor allele is rare. Unrelated people carry the common allele there
too, so the reference loses its grip on them. The two effects pull in opposite directions and
m = n - 1 puts a deep reference on the wrong side of the trade:

  n=6  m>=5   retained heterozygosity 93% of true   unrelated diploid held  4.15% -> 4.06%
  n=5  m>=4   retained heterozygosity 93% of true   unrelated diploids improved
  n=8  m>=7   retained heterozygosity ~62% of true  unrelated diploid fell  4.15% -> 2.49%

This sweeps m across its whole range on both fathers whose true heterozygosity is known from a
real array, and scores the three things that decide a rule: do true products come back, do
negatives stay out, and how far has the marker set drifted from the genome.

Run: python tools/step6_ascertainment.py <geo-dir> <out-dir>
"""
from __future__ import annotations

import csv
import gc
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from originmarker import origin  # noqa: E402
from tools.inferred_reference import Products, build  # noqa: E402

GEO = Path(sys.argv[1]).expanduser()
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out").expanduser()
OUT.mkdir(parents=True, exist_ok=True)
PRODUCT = "Affymetrix Axiom"

PRODUCTS = ["GSM4774680", "GSM4774681", "GSM4774682", "GSM4774683", "GSM4774685"]
FATHER = "GSM4472397"
#: Everything that is not a product, labelled by what it should read. The father's second array
#: is the strongest positive control available: a genotype reconstructed from his sons must
#: recognise the man himself.
PANEL = [("GSM4472398", "father rep 2", "present"),
         ("GSM4774673", "maternal pronucleus", "absent"),
         ("GSM4774676", "maternal pronucleus", "absent"),
         ("GSM4774678", "maternal pronucleus", "absent"),
         ("GSM4472407", "egg donor, diploid", "absent"),
         ("GSM4472415", "egg donor, diploid", "absent")]


def path_for(gsm: str) -> Path:
    hits = sorted(GEO.glob(f"{gsm}_*.CEL.txt.gz"))
    if not hits:
        raise SystemExit(f"missing {gsm}")
    return hits[0]


def read(gsm: str) -> origin.Sample:
    return next(iter(origin.read_samples(path_for(gsm)).values()))


father = read(FATHER)
auto = [p for p in father.probes.values() if origin._is_autosome(p.chrom) and p.gt != "NC"]
H_TRUE = sum(1 for p in auto if p.gt == "AB") / len(auto)
print(f"father's true autosomal heterozygosity {H_TRUE:.4%}\n", flush=True)
del auto

p = Products()
for gsm in PRODUCTS:
    s = read(gsm)
    p.add(gsm, s)
    del s
    gc.collect()

# Read the panel once. Each sample is scored against every m, so re-reading per m would be
# the whole cost of the run.
panel = {gsm: (read(gsm), what, expect) for gsm, what, expect in PANEL}

rows = []
for m_min in range(2, len(PRODUCTS) + 1):
    ref = build(p, m_min)
    ascertainment = ref.h_retained / H_TRUE if H_TRUE else float("nan")
    print(f"== m>={m_min}: {ref.markers:,} markers, mean m {ref.mean_m:.2f}, "
          f"h(retained) {ref.h_retained:.2%} = {ascertainment:.0%} of true, "
          f"contamination {ref.contamination:.3%}", flush=True)

    # Products leave one out, so a member is never judged by a reference containing itself.
    prod_ok = prod_wrong = 0
    for gsm in PRODUCTS:
        loo = build(p, min(m_min, len(PRODUCTS) - 1), exclude=gsm)
        s = read(gsm)
        r = origin.parental_origin(loo.sample, s, product=PRODUCT)
        ratio = r.genome_rate / r.explainable if r.explainable else float("nan")
        prod_ok += r.verdict == "parent_genome_present"
        prod_wrong += r.verdict == "no_parental_contribution"
        rows.append({"m_min": m_min, "gsm": gsm, "what": "product (leave one out)",
                     "expect": "present", "n": loo.n_products,
                     "markers": loo.markers, "h_retained": f"{loo.h_retained:.6f}",
                     "ascertainment": f"{loo.h_retained / H_TRUE:.4f}",
                     "contamination": f"{loo.contamination:.6f}",
                     "absence": f"{r.genome_rate:.6f}", "ceiling": f"{r.explainable:.6f}",
                     "ratio": f"{ratio:.4f}", "verdict": r.verdict})
        del s, loo
        gc.collect()

    neg_min = float("inf")
    neg_bad = 0
    for gsm, (s, what, expect) in panel.items():
        r = origin.parental_origin(ref.sample, s, product=PRODUCT)
        ratio = r.genome_rate / r.explainable if r.explainable else float("nan")
        got = ("present" if r.verdict == "parent_genome_present"
               else "absent" if r.verdict == "no_parental_contribution" else "unclear")
        if expect == "absent":
            neg_min = min(neg_min, ratio)
            neg_bad += got != "absent"
        rows.append({"m_min": m_min, "gsm": gsm, "what": what, "expect": expect,
                     "n": ref.n_products, "markers": ref.markers,
                     "h_retained": f"{ref.h_retained:.6f}",
                     "ascertainment": f"{ascertainment:.4f}",
                     "contamination": f"{ref.contamination:.6f}",
                     "absence": f"{r.genome_rate:.6f}", "ceiling": f"{r.explainable:.6f}",
                     "ratio": f"{ratio:.4f}", "verdict": r.verdict})
        print(f"   {gsm} {what:22} expect {expect:8} got {got:8} {ratio:7.2f}x", flush=True)
    print(f"   -> products present {prod_ok}/{len(PRODUCTS)}, inverted {prod_wrong}, "
          f"negatives wrong {neg_bad}, closest negative {neg_min:.2f}x\n", flush=True)
    del ref
    gc.collect()

with open(OUT / "step6_ascertainment.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(f"wrote {OUT / 'step6_ascertainment.csv'}", flush=True)
print("STEP6_DONE", flush=True)

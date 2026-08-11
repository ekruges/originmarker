"""Generate the cross-implementation fixture from origin.py.

Run: python3 tools/gen_parentage_fixture.py

The two implementations of the parentage arithmetic are meant to agree exactly. Nothing compared
them until an independent audit found them differing on every pair it tried, because
`parentage.check.ts` pinned the browser against measured values, `origin.py` had its own suite,
and neither looked at the other.

This writes what the CLI computes on the shipped example arrays; `parentage.check.ts` asserts the
browser reproduces it. A divergence now fails a test rather than shipping.
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from originmarker import origin  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "web" / "public" / "examples"
DONOR = "GSM4472397_sperm_DNA_71.subset.csv.gz"
SAMPLES = ["GSM4472409_A8_45", "GSM4472407_donor_A_47", "GSM4472415_donor_C_70",
           "GSM4472398_sperm_DNA_79"]


def load(name: str) -> origin.Sample:
    with tempfile.NamedTemporaryFile("wb", suffix=".csv", delete=False) as out:
        with gzip.open(EX / name, "rb") as fh:
            shutil.copyfileobj(fh, out)
        tmp = out.name
    try:
        return origin.read_sample(tmp)
    finally:
        Path(tmp).unlink(missing_ok=True)


def main() -> None:
    father = load(DONOR)
    auto = [p for p in father.probes.values() if origin._is_autosome(p.chrom) and p.gt != "NC"]
    het = sum(1 for p in auto if p.gt == "AB") / len(auto)
    cases = []
    for stem in SAMPLES:
        s = load(next(p.name for p in EX.iterdir() if p.name.startswith(stem)))
        r = origin.parental_origin(father, s)
        cases.append({
            "sample": stem,
            "verdict": r.verdict,
            "origin_class": r.origin_class,
            "zygosity": r.zygosity,
            "sperm_type": r.sperm_type,
            "genome_rate": round(r.genome_rate, 9),
            "explainable": round(r.explainable, 9),
            "nonpaternal_rate": round(r.nonpaternal_rate, 9),
            "second_parent_expected": round(r.second_parent_expected, 9),
            "het_band": round(r.het_band, 9),
            "dispersion": round(r.dispersion, 9),
            "min_chrom_rate": round(r.min_chrom_rate, 9),
            "no_call_rate": round(r.no_call_rate, 9),
        })
    out = {"generated_by": "tools/gen_parentage_fixture.py", "donor": DONOR,
           "donor_heterozygosity": round(het, 9), "cases": cases}
    dst = ROOT / "tests" / "fixtures" / "parentage_cross.json"
    dst.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {dst}")
    for c in cases:
        print(f"  {c['sample']:26} {c['origin_class']:13} {c['genome_rate']:.6f} "
              f"of {c['explainable']:.6f}  {c['sperm_type']}")


if __name__ == "__main__":
    main()

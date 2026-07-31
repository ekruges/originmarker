"""Format-coverage audit: the same question, asked through every reader.

Run: python3 audit/formats.py <raw-dir> <out-dir>

The pronucleus audit in `cases.ts` tests the statistics hard and the ingestion layer through
exactly one path: 27 Affymetrix Axiom exports, all writing numeric genotype codes. That gap hid a
defect where an ordinary AB-space file read as 100% no-call, so this suite drives the other
readers against public files with documented relationships.

Truth comes from published pedigrees, not from this tool:

  HapMap CEU  International HapMap Project phase II+III, chr22, forward strand, build 36. 174
              individuals, of which 52 are children with both parents in the same file. Genotypes
              are nucleotide pairs, so this also exercises the pooled A/B assignment. Pedigree
              from the 1000 Genomes sample-info file, which carries the HapMap families.
  1000G VCF   Phase 3 related-samples release, chr22. A whole-genome variant callset rather than
              a polymorphic panel, which is a case worth having on the record for what it does to
              the rates.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import origin  # noqa: E402

RAW = Path(sys.argv[1] if len(sys.argv) > 1 else "raw")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out")
OUT.mkdir(parents=True, exist_ok=True)

HAPMAP = RAW / "hapmap_chr22_CEU.txt.gz"
VCF = RAW / "g1k_chr22_related.vcf.gz"
PED = RAW / "g1k.ped"

rows: list[dict] = []


def pedigree() -> dict[str, tuple[str, str, str]]:
    """child -> (father, mother, population), from the published sample-info file."""
    out = {}
    for line in PED.read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 7 and f[2] != "0" and f[3] != "0":
            out[f[1]] = (f[2], f[3], f[6])
    return out


def record(**kw) -> None:
    rows.append(kw)
    verdict = kw["got"]
    ok = {True: "ok  ", False: "FAIL", None: "...."}[kw["pass"]]
    print(f"  {ok} {kw['id']:38} {verdict:26} {kw['absence']:7.3%} of {kw['ceiling']:7.3%}"
          f"  {kw['margin']:6.2f}x")


def score(parent, sample, role="paternal"):
    r = origin.parental_origin(parent, sample, role=role)
    return r, (r.genome_rate / r.explainable if r.explainable else float("nan"))


# --- HapMap CEU: 52 documented trios, nucleotide space, bare sample columns --------------------
print(f"reading {HAPMAP.name}")
d = origin.read_samples(HAPMAP)
print(f"  {len(d)} samples, {len(next(iter(d.values())).probes)} markers on chr22")
ped = pedigree()
trios = [(c, p, m, pop) for c, (p, m, pop) in ped.items()
         if c in d and p in d and m in d]
trios.sort()
print(f"  {len(trios)} complete trios present\n")

print("== HapMap CEU: child against its own father, mother, and an unrelated adult")
# Harmonised together, which is the point: the A/B assignment at each marker is decided once
# across every sample rather than per file.
members = sorted({x for t in trios for x in t[:3]})
harm = {s.sample_id: s for s in origin.harmonise([d[k] for k in members])}

for child, father, mother, pop in trios:
    r, x = score(harm[father], harm[child])
    record(id=f"hapmap-{child}-vs-father", suite="HapMap CEU trios", fmt="HapMap raw (nucleotide)",
           source=f"child {child}, father {father} ({pop})",
           expect="paternal genome present", got=r.verdict,
           absence=r.genome_rate, ceiling=r.explainable, margin=x,
           informative=r.informative_markers if hasattr(r, "informative_markers") else 0,
           note="", **{"pass": r.verdict == "paternal_genome_present"})

    r, x = score(harm[mother], harm[child], role="maternal")
    record(id=f"hapmap-{child}-vs-mother", suite="HapMap CEU trios", fmt="HapMap raw (nucleotide)",
           source=f"child {child}, mother {mother} ({pop})",
           expect="maternal genome present", got=r.verdict,
           absence=r.genome_rate, ceiling=r.explainable, margin=x, informative=0,
           note="Role-neutral in the browser; the Python still words this field paternally.",
           **{"pass": r.verdict == "paternal_genome_present"})

    r, x = score(harm[father], harm[mother])
    record(id=f"hapmap-{child}-parents-unrelated", suite="HapMap CEU trios",
           fmt="HapMap raw (nucleotide)",
           source=f"the two parents of {child}: {father} and {mother} ({pop})",
           expect="no contribution: spouses are unrelated", got=r.verdict,
           absence=r.genome_rate, ceiling=r.explainable, margin=x, informative=0, note="",
           **{"pass": r.verdict == "no_paternal_contribution"})

# --- 1000 Genomes VCF -------------------------------------------------------------------------
print("\n== 1000 Genomes phase 3 VCF, chr22 related samples")
with gzip.open(VCF, "rt") as fh:
    names = next(l.split() for l in fh if l.startswith("#CHROM"))[9:]
pairs = [(c, p, "father") for c, (p, m, _) in ped.items() if c in names and p in names]
pairs += [(c, m, "mother") for c, (p, m, _) in ped.items() if c in names and m in names]

cache = {}


def vcf(name):
    if name not in cache:
        cache[name] = origin.read_vcf(VCF, name)
    return cache[name]


for child, parent, which in pairs:
    r, x = score(vcf(parent), vcf(child), role="paternal" if which == "father" else "maternal")
    record(id=f"vcf-{child}-vs-{which}", suite="1000 Genomes VCF", fmt="VCF (GT)",
           source=f"child {child}, {which} {parent}",
           expect=f"{which} genome present", got=r.verdict,
           absence=r.genome_rate, ceiling=r.explainable, margin=x, informative=0,
           note=" ".join(r.limits)[:400],
           **{"pass": r.verdict == "paternal_genome_present"})

others = [n for n in names if n not in {p for _, p, _ in pairs} | {c for c, _, _ in pairs}][:6]
base = pairs[0][1]
for o in others:
    r, x = score(vcf(base), vcf(o))
    record(id=f"vcf-{base}-vs-{o}-unrelated", suite="1000 Genomes VCF", fmt="VCF (GT)",
           source=f"{base} against {o}, no documented relationship",
           expect="no contribution", got=r.verdict,
           absence=r.genome_rate, ceiling=r.explainable, margin=x, informative=0,
           note=" ".join(r.limits)[:400],
           **{"pass": r.verdict == "no_paternal_contribution"})

scored = [r for r in rows if r["pass"] is not None]
ok = sum(1 for r in scored if r["pass"])
refused = sum(1 for r in scored if not r["pass"] and r["got"] == "unclear")
print(f"\n===== {ok} correct, {refused} refused, {len(scored) - ok - refused} incorrect, "
      f"of {len(scored)}")
(OUT / "formats.json").write_text(json.dumps({"rows": rows}, indent=1))
print(f"wrote {OUT / 'formats.json'}")

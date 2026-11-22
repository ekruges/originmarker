"""Recount the run's overview table from the per-sample records, sharing no code with it.

Reads the table as web/src/overview.ts prints it (audit/truth/overview.ts) and the records it was
built from, recounts every line from the records with the rule as written in
audit/truth/overview-counts-2026-11-22.txt, and exits 1 on any mismatch, on any line it cannot
account for, or on a class that does not sum to its total.

Run: python3 audit/truth/overview_check.py <table.tsv> <records.jsonl> [more.jsonl ...]
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
MARGIN = 5_000_000
GAIN = {"copy-gain", "trisomy", "segmental-duplication"}
LOSS = {"copy-loss", "parental-absence", "monosomy", "segmental-deletion"}
WORDS = {
    "meiotic": "from the gamete",
    "consistent-with-gamete": "consistent with the gamete (two cells read)",
    "mitotic-reciprocal": "after fertilisation",
    "mitotic-not-all": "after fertilisation",
}
STAGE_WORD = {"single-cell": "single cell", "": "no stage called"}


def centromeres():
    cen, end = {}, {}
    for chrom, start, stop, kind in json.loads((ROOT / "web/public/hg19_features.json").read_text())["gaps"]:
        end[chrom] = max(end.get(chrom, 0), stop)
        if kind == "centromere":
            cen[chrom] = (start, stop)
    return cen, end


def recount(records):
    cen, end = centromeres()

    def extent(row):
        if row["whole"]:
            return "whole"
        if row["chrom"] not in cen:
            return "piece"
        lo, hi = cen[row["chrom"]][0] - MARGIN, cen[row["chrom"]][1] + MARGIN
        near = lambda x: lo <= x <= hi
        if row["start"] <= MARGIN and near(row["end"]):
            return "arm"
        if near(row["start"]) and row["end"] >= end[row["chrom"]] - MARGIN:
            return "arm"
        return "piece"

    def direction(row):
        if row["gain"]:
            return "gain"
        return "loss" if any(p["kind"] in LOSS for p in row["pieces"]) else None

    recs = [r for r in records if r.get("card")]
    gates = lambda r: {x["gate"] for x in r.get("refusals", []) if x["claim"] == "sample"}
    cls = lambda r: ("not reportable" if "integrity.uninterpretable" in gates(r)
                     else "limited" if "integrity.disturbed" in gates(r) else "reportable")
    want = {"samples": len(recs)}
    stage = lambda r: (r.get("stage") or {}).get("stage") or ""
    for st in ["single-cell", "blastomere", "trophectoderm", "haploid", "bulk", "failed"] + sorted({stage(r) for r in recs}):
        want["stage: " + STAGE_WORD.get(st, st)] = sum(stage(r) == st for r in recs)
    want["reportable"] = sum(cls(r) == "reportable" for r in recs)
    want["limited (integrity disturbed)"] = sum(cls(r) == "limited" for r in recs)
    want["not reportable (integrity uninterpretable)"] = sum(cls(r) == "not reportable" for r in recs)
    read = [r for r in recs if cls(r) != "not reportable"]
    rows = [g for r in read for g in r["card"]["counted"]]
    changes = [g for g in rows if direction(g)]
    want["changes (as the result cards count them)"] = len(changes)
    for d, label in (("loss", "losses"), ("gain", "gains")):
        of = [g for g in changes if direction(g) == d]
        want[label] = len(of)
        for e in ("whole", "arm", "piece"):
            want[f"{d}, {e}"] = sum(extent(g) == e for g in of)
        for p in ("maternal", "paternal"):
            want[f"{d}, {p}"] = sum(g["origin"] == p for g in of)
        want[f"{d}, parent not named"] = sum(g["origin"] not in ("maternal", "paternal") for g in of)
    want["findings (copy-neutral, uneven shares, complex genome)"] = (
        sum(1 for g in rows if not direction(g)) + sum(1 for r in read if r["card"]["disturbed"]))
    want["stretches seen once, not counted"] = sum(r["card"]["once"] for r in read)
    sizes = {}
    for r in recs:
        if r.get("embryo") is not None:
            sizes[r["embryo"]] = sizes.get(r["embryo"], 0) + 1
    want["embryos of two or more units"] = len(sizes)
    want["samples in them"] = sum(sizes.values())
    want["samples in no embryo group"] = sum(r.get("embryo") is None for r in recs)
    for k in set(sizes.values()):
        want[f"embryos of {k} units"] = sum(v == k for v in sizes.values())
    said = lambda g: {WORDS[p["mechanism"]] for p in g["pieces"] if p.get("mechanism") in WORDS}
    for w in ("from the gamete", "consistent with the gamete (two cells read)", "after fertilisation"):
        want["origin: " + w] = sum(said(g) == {w} for g in changes)
    want["origin: pieces read differently"] = sum(len(said(g)) > 1 for g in changes)
    want["origin: no origin word"] = sum(not said(g) for g in changes)

    sums = [
        ("stage lines", sum(v for k, v in want.items() if k.startswith("stage: ")), want["samples"]),
        ("reportable classes", want["reportable"] + want["limited (integrity disturbed)"]
         + want["not reportable (integrity uninterpretable)"], want["samples"]),
        ("losses and gains", want["losses"] + want["gains"], want["changes (as the result cards count them)"]),
        ("embryo lines", want["samples in them"] + want["samples in no embryo group"], want["samples"]),
        ("origin lines", sum(v for k, v in want.items() if k.startswith("origin: ")),
         want["changes (as the result cards count them)"]),
    ]
    for d, label in (("loss", "losses"), ("gain", "gains")):
        sums.append((f"{label} by extent", sum(want[f"{d}, {e}"] for e in ("whole", "arm", "piece")), want[label]))
        sums.append((f"{label} by parent", want[f"{d}, maternal"] + want[f"{d}, paternal"]
                     + want[f"{d}, parent not named"], want[label]))
    unreportable_changes = sum(1 for r in recs if cls(r) == "not reportable" for g in r["card"]["counted"])
    return want, sums, unreportable_changes


def main(argv):
    table = {}
    for line in pathlib.Path(argv[1]).read_text().splitlines()[1:]:
        label, count = line.rsplit("\t", 1)
        table[label] = int(count)
    records = [json.loads(l) for f in argv[2:] for l in pathlib.Path(f).read_text().splitlines() if l.startswith("{")]
    want, sums, held = recount(records)
    bad = [f"{k}: table {table.get(k)} recount {v}" for k, v in want.items() if table.get(k) != v]
    bad += [f"{k}: in the table, not in the rule" for k in table if k not in want]
    bad += [f"{name} sum to {a}, total {b}" for name, a, b in sums if a != b]
    for b in bad:
        print("MISMATCH", b)
    print(f"overview check: {len(table)} lines, {len(bad)} mismatches; "
          f"{held} counted rows on samples not reportable, none in the table")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

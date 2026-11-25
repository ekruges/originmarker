"""Gate G2 of the shipped GRCh38 column: the old live-map method and the shipped column give the
same chip set ({GRCh38 position: probe}) on every recorded assembly-map answer and on live windows.

The old method: Ensembl's map of the GRCh38 window to GRCh37 blocks, each probe's GRCh37 position
carried through its block. The new: the shipped pos_grch38 column, read inside the window.

Run: ENSEMBL_ADDR=<address> python3 audit/truth/chip-grch38-check.py tests/fixtures <chrom:lo-hi> ...
"""
import bisect
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
table: dict = {}
with gzip.open(ROOT / "data" / "chip" / "axiom-ukb-probes.tsv.gz", "rt") as f:
    for line in f:
        if line.startswith(("#", "probeset_id")):
            continue
        pid, ch, p37, p38 = line.rstrip("\n").split("\t")
        table.setdefault(ch, []).append((int(p37), pid, int(p38) if p38 else None))
for v in table.values():
    v.sort()


def old(chrom: str, maps: list) -> dict:
    rows = table.get(chrom, []); p37 = [r[0] for r in rows]; out = {}
    for mp in maps:
        o, m = mp["original"], mp["mapped"]
        if (m.get("seq_region_name") != chrom or m.get("assembly") != "GRCh37"
                or m["end"] - m["start"] != o["end"] - o["start"]):
            continue
        for i in range(bisect.bisect_left(p37, m["start"]), bisect.bisect_right(p37, m["end"])):
            off = p37[i] - m["start"]
            out[o["start"] + off if m["strand"] == o["strand"] else o["end"] - off] = rows[i][1]
    return out


def new(chrom: str, lo: int, hi: int) -> dict:
    return {p38: pid for _, pid, p38 in table.get(chrom, []) if p38 is not None and lo <= p38 <= hi}


def compare(label: str, chrom: str, lo: int, hi: int, maps: list) -> bool:
    a, b = old(chrom, maps), new(chrom, lo, hi)
    same = a == b
    print(f"{label}\tchr{chrom}:{lo:,}-{hi:,}\told {len(a)}\tnew {len(b)}\t{'identical' if same else 'DIFFER'}")
    if not same:
        for p in sorted(set(a) ^ set(b) | {p for p in set(a) & set(b) if a[p] != b[p]}):
            print(f"  {p:,}\told {a.get(p)}\tnew {b.get(p)}")
    return same


ok = True
for fx in sorted(Path(sys.argv[1]).glob("*.json")):
    try:
        d = json.loads(fx.read_text())
    except ValueError:
        continue                      # an XML answer, not an assembly map
    maps = d.get("mappings") if isinstance(d, dict) else None
    if not (maps and isinstance(maps[0], dict) and isinstance(maps[0].get("original"), dict)
            and maps[0]["original"].get("assembly") == "GRCh38"):
        continue
    chrom = maps[0]["original"]["seq_region_name"]
    lo = min(m["original"]["start"] for m in maps); hi = max(m["original"]["end"] for m in maps)
    ok &= compare("recorded", chrom, lo, hi, maps)
pin = ["--resolve", f"rest.ensembl.org:443:{os.environ['ENSEMBL_ADDR']}"] if os.environ.get("ENSEMBL_ADDR") else []
for w in sys.argv[2:]:
    chrom, span = w.split(":"); lo, hi = (int(x) for x in span.split("-"))
    url = f"https://rest.ensembl.org/map/human/GRCh38/{chrom}:{lo}..{hi}:1/GRCh37?content-type=application/json"
    r = subprocess.run(["curl", "-s", "-m", "60", *pin, url], capture_output=True, text=True, check=True)
    ok &= compare("live", chrom, lo, hi, json.loads(r.stdout)["mappings"])
print("all identical" if ok else "NOT all identical")
sys.exit(0 if ok else 1)

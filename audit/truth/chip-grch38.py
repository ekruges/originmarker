"""Each chip probe's GRCh38 position, from Ensembl's assembly map, taken once. A one-off.

Reads data/chip/axiom-ukb-probes.tsv.gz (probe id, chromosome, GRCh37 position), asks Ensembl's
/map/human/GRCh37/<chrom>:1..<length>:1/GRCh38 once per chromosome, and writes the same file with a
fourth column, pos_grch38. A probe gets a GRCh38 position only where it lies inside exactly one
aligned block that maps to its own chromosome; a probe in no block, in a block mapped to another
chromosome, in a block whose two sides differ in length, or in more than one block, is written with
the column empty and is never offered. Prints the Ensembl release and the counts by reason.

The host name may be pinned to an address (ENSEMBL_ADDR), as curl --resolve does: the name and the
certificate check stay. One call per chromosome, a pause between, and a stop on repeated 429 or 5xx.

Run: ENSEMBL_ADDR=<address> python3 audit/truth/chip-grch38.py data/chip/axiom-ukb-probes.tsv.gz
"""
import bisect
import gzip
import json
import os
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from originmarker import buildref  # noqa: E402

HOST = "rest.ensembl.org"
ADDR = os.environ.get("ENSEMBL_ADDR")
LEN37 = dict(buildref.CHROM_LEN["GRCh37"], MT=16569)


def get(path: str) -> dict:
    """One GET by curl, retried with a growing pause; stops after three 429 or 5xx answers."""
    url = f"https://{HOST}{path}{'&' if '?' in path else '?'}content-type=application/json"
    pin = ["--resolve", f"{HOST}:443:{ADDR}"] if ADDR else []
    for attempt in range(3):
        r = subprocess.run(["curl", "-s", "-m", "120", "-w", "\n%{http_code}", *pin, url],
                           capture_output=True, text=True)
        body, _, code = r.stdout.rpartition("\n")
        if code == "200":
            return json.loads(body)
        print(f"  {path}: answered {code or 'nothing'}, attempt {attempt + 1}", file=sys.stderr)
        time.sleep(10 * (attempt + 1))
    sys.exit(f"stopped: {path} did not answer 200 three times")


src = Path(sys.argv[1])
with gzip.open(src, "rt") as f:
    lines = f.read().splitlines()
head = [l for l in lines if l.startswith("#")]
rows = [l.split("\t") for l in lines if l and not l.startswith(("#", "probeset_id"))]
release = get("/info/data")["releases"]
print(f"Ensembl release {release}")

by_chrom: dict = {}
for i, r in enumerate(rows):
    by_chrom.setdefault(r[1], []).append(i)

pos38 = [""] * len(rows)
why = Counter()
for chrom, idx in sorted(by_chrom.items()):
    idx.sort(key=lambda i: int(rows[i][2]))
    p37 = [int(rows[i][2]) for i in idx]
    maps = get(f"/map/human/GRCh37/{chrom}:1..{LEN37[chrom]}:1/GRCh38")["mappings"]
    time.sleep(1)
    hits: dict = {}
    for b, mp in enumerate(maps):
        o = mp["original"]
        for k in range(bisect.bisect_left(p37, o["start"]), bisect.bisect_right(p37, o["end"])):
            hits.setdefault(k, []).append(b)
    for k, i in enumerate(idx):
        h = hits.get(k, [])
        if not h:
            why["no mapping"] += 1
            continue
        if len(h) > 1:
            why["mapped to more than one place"] += 1
            continue
        o, m = maps[h[0]]["original"], maps[h[0]]["mapped"]
        if m.get("seq_region_name") != chrom:
            why["mapped to another chromosome"] += 1
            continue
        if m["end"] - m["start"] != o["end"] - o["start"]:
            why["in a block whose two sides differ in length"] += 1
            continue
        off = p37[k] - o["start"]
        pos38[i] = str(m["start"] + off if m["strand"] == o["strand"] else m["end"] - off)
        why["placed"] += 1
    print(f"chr{chrom}: {len(idx):,} probes, {len(maps):,} blocks")

with gzip.open(src, "wt", compresslevel=9) as f:
    for h in head:
        if not h.startswith("# Built by"):
            f.write(h + "\n")
    f.write(f"# pos_grch38 from Ensembl's assembly map, release {release[0]}, by audit/truth/chip-grch38.py;\n"
            "# empty where the probe has no single GRCh38 position, and such a probe is never offered.\n")
    f.write("# Built by audit/truth/chip-probes.py and audit/truth/chip-grch38.py. No genotype is in this file.\n")
    f.write("probeset_id\tchr\tpos_grch37\tpos_grch38\n")
    for r, p in zip(rows, pos38):
        f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{p}\n")
print(f"probes {len(rows):,}")
for k, v in why.most_common():
    print(f"  {k}\t{v:,}")

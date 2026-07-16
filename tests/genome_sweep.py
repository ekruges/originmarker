"""
Genome-wide sweep against live data. Not part of the offline suite; run it by hand.

    .venv/bin/python tests/genome_sweep.py            # resolve sweep, all chromosomes
    .venv/bin/python tests/genome_sweep.py --builds 3 # + N full panel builds

The offline suite covers one variant. This asks whether the tool holds across the genome,
against real records, without fooling itself:

1. Variants are sampled live from ClinVar per chromosome, never written from memory, so
   the sweep cannot encode a wrong rsID and then "confirm" it.
2. Coordinates are cross-checked against Ensembl, a different database on a different
   endpoint. The app agreeing with itself is not corroboration.
3. A refusal is not a failure. The failures worth finding are silent: a coordinate that
   disagrees with Ensembl, a strand asserted without evidence, a transcript-sense reading
   containing a fabricated base.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
os.environ.setdefault("PANELBUILDER_CACHE", "/tmp/om_genome_sweep")

import panelbuilder as pb  # noqa: E402

CHROMS = [str(c) for c in range(1, 23)] + ["X"]

# NCBI asks for <=3 req/s unauthenticated; exceeding it gets the caller's IP blocked.
THROTTLE = 0.5


def sample_clinvar(chrom: str, extra: str, n: int = 1) -> list[str]:
    """Live-sample ClinVar accessions for a chromosome. Never invents an identifier."""
    term = f'{extra} AND {chrom}[Chromosome] AND "clinsig pathogenic"[Properties]'
    try:
        r = json.loads(pb._get(pb.EUTILS, "/esearch.fcgi", pb._eutils_params(
            {"db": "clinvar", "term": term, "retmode": "json", "retmax": str(n)})))
        return r["esearchresult"]["idlist"][:n]
    except Exception:
        return []


def accession_of(uid: str) -> str | None:
    try:
        s = json.loads(pb._get(pb.EUTILS, "/esummary.fcgi", pb._eutils_params(
            {"db": "clinvar", "id": uid, "retmode": "json"})))["result"][uid]
        return s.get("accession")
    except Exception:
        return None


def ensembl_pos(rsid: str) -> int | None:
    """Independent coordinate for an rsID, from Ensembl rather than ClinVar."""
    try:
        rec = json.loads(pb._get(pb.ENSEMBL, f"/variation/homo_sapiens/{rsid}",
                                 {"content-type": "application/json"}, tries=2))
        m = next((m for m in rec.get("mappings", [])
                  if m.get("assembly_name") == "GRCh38"), None)
        return int(m["start"]) if m else None
    except Exception:
        return None


def ensembl_strand(gene: str) -> int | None:
    try:
        g = json.loads(pb._get(pb.ENSEMBL, f"/lookup/symbol/homo_sapiens/{gene}",
                               {"content-type": "application/json"}, tries=2))
        return int(g.get("strand"))
    except Exception:
        return None


def check(query: str, label: str) -> dict:
    """Resolve one variant and interrogate the answer. Returns a row for the report."""
    row: dict = {"query": query, "label": label, "status": "?", "notes": []}
    try:
        v = pb.resolve_variant(query)
    except pb.ApiError as e:
        row["status"] = "refused"
        row["detail"] = str(e)[:90]
        return row
    except Exception as e:  # noqa: BLE001
        row["status"] = "CRASH"
        row["detail"] = f"{type(e).__name__}: {e}"[:90]
        return row

    row.update(status="resolved", rsid=v.rsid, gene=v.gene, chrom=v.chrom,
               pos=v.pos_grch38, ref=v.vcf_ref, alt=v.vcf_alt, strand=v.strand,
               sig=v.clinical_significance, ts=v.transcript_sense_change())

    if v.rsid:
        e = ensembl_pos(v.rsid)
        if e is None:
            row["notes"].append("ensembl: no GRCh38 mapping to compare")
        elif e != v.pos_grch38:
            row["notes"].append(f"COORD MISMATCH clinvar={v.pos_grch38} ensembl={e}")
            row["status"] = "MISMATCH"
        else:
            row["notes"].append("coord confirmed by ensembl")

    if v.gene and v.strand is not None:
        es = ensembl_strand(v.gene)
        if es is not None and es != v.strand:
            row["notes"].append(f"STRAND MISMATCH clinvar={v.strand} ensembl={es}")
            row["status"] = "MISMATCH"
        elif es is not None:
            row["notes"].append("strand confirmed by ensembl")

    if "?" in (row["ts"] or ""):
        row["notes"].append("FABRICATED BASE in transcript sense")
        row["status"] = "MISMATCH"
    if v.strand is None and "unknown" not in (row["ts"] or ""):
        row["notes"].append("strand unknown but not declared unknown")
        row["status"] = "MISMATCH"

    if not v.vcf_ref or not v.vcf_alt:
        row["notes"].append(f"EMPTY ALLELE ref={v.vcf_ref!r} alt={v.vcf_alt!r}")
        row["status"] = "MISMATCH"
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--builds", type=int, default=0, help="also run N full panel builds")
    args = ap.parse_args()

    print("Sampling ClinVar live, one pathogenic variant per chromosome.")
    print("Nothing here is written from memory.\n")

    targets: list[tuple[str, str]] = []
    for c in CHROMS:
        uids = sample_clinvar(c, '"single nucleotide variant"[Type]')
        time.sleep(THROTTLE)
        if not uids:
            print(f"  chr{c:<3} no ClinVar sample returned")
            continue
        acc = accession_of(uids[0])
        time.sleep(THROTTLE)
        if acc:
            targets.append((acc, f"chr{c} SNV"))
            print(f"  chr{c:<3} {acc}")

    # Variant classes the golden case never exercises.
    for typ, lbl in [('"Deletion"[Type]', "deletion"),
                     ('"Duplication"[Type]', "duplication"),
                     ('"Insertion"[Type]', "insertion"),
                     ('"Indel"[Type]', "indel")]:
        uids = sample_clinvar("1", typ) or sample_clinvar("2", typ)
        time.sleep(THROTTLE)
        if uids:
            acc = accession_of(uids[0])
            time.sleep(THROTTLE)
            if acc:
                targets.append((acc, lbl))
                print(f"  {lbl:<12} {acc}")

    print(f"\nResolving {len(targets)} real records\n")
    rows = []
    for q, lbl in targets:
        r = check(q, lbl)
        rows.append(r)
        if r["status"] == "resolved":
            print(f"  ok       {lbl:<14} {r.get('rsid') or r['query']:<14} "
                  f"{r.get('gene') or '?':<12} chr{r['chrom']}:{r['pos']} "
                  f"{r['ref']}>{r['alt']} strand={r['strand']}")
        elif r["status"] == "refused":
            print(f"  refused  {lbl:<14} {r['query']:<14} {r['detail'][:52]}")
        else:
            print(f"  {r['status']:<8} {lbl:<14} {r['query']:<14} {r.get('detail','')}")
        for n in r["notes"]:
            if n.isupper() or "MISMATCH" in n or "FABRICATED" in n or "EMPTY" in n:
                print(f"           ^^ {n}")
        time.sleep(THROTTLE)

    ok = [r for r in rows if r["status"] == "resolved"]
    ref = [r for r in rows if r["status"] == "refused"]
    bad = [r for r in rows if r["status"] in ("MISMATCH", "CRASH")]
    corroborated = [r for r in ok if any("confirmed by ensembl" in n for n in r["notes"])]

    print("\n" + "=" * 72)
    print(f"resolved {len(ok)}/{len(rows)}   refused {len(ref)}   "
          f"WRONG {len(bad)}")
    print(f"coordinates independently corroborated by Ensembl: {len(corroborated)}")
    print(f"chromosomes covered: "
          f"{sorted({r['chrom'] for r in ok}, key=lambda c: (c.isdigit() is False, c))}")
    strands = {r["strand"] for r in ok}
    print(f"strands seen: {sorted(s for s in strands if s is not None)}"
          f"{' + unknown' if None in strands else ''}")
    print(f"variant classes: {sorted({r['label'].split()[-1] for r in ok})}")

    if bad:
        print("\nFAILURES THAT ARE ACTUALLY WRONG (not refusals):")
        for r in bad:
            print(f"  {r['query']}: {r.get('detail','')} {r['notes']}")
        return 1
    print("\nNo silent wrongness found: every resolved coordinate that Ensembl could")
    print("corroborate matched, no fabricated bases, no undeclared strands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

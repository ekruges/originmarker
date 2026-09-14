"""Genotype calls from a GPL3718 series matrix into the tab-separated shape OriginMarker reads.

The Affymetrix 250K Nsp series in GEO ship their genotype calls inside the series matrix, one
column per sample, and the platform table places markers on NCBI Build 35. OriginMarker works in
GRCh37 or GRCh38, so every position is lifted to GRCh37 with the UCSC hg17ToHg19 chain before it is
written. A marker that does not lift, or lifts onto a different chromosome, is dropped rather than
placed.

These matrices carry calls only. The log2R and BAF columns are written empty: the intensity
channels have nothing to read and must see nothing rather than a stand-in.

Per-sample metadata goes to truth.json exactly as the series states it, every characteristics field
included, so a harness reads the submitters' answer and never this file's interpretation of it.

Run inside a Python with pyliftover:
  python convert.py <series_matrix.txt.gz> <annot_b35.tsv> <hg17ToHg19.over.chain.gz> <outdir>
"""
import gzip
import json
import os
import sys

from pyliftover import LiftOver

HEAD = "probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset"
CALLS = {"AA": "AA", "AB": "AB", "BB": "BB"}


def load_positions(annot, chain):
    """probe -> (chrom, GRCh37 position), for markers that lift onto the same chromosome."""
    lo = LiftOver(chain)
    out = {}
    dropped = 0
    with open(annot) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or not p[2].isdigit():
                continue
            pid, c, pos = p[0], p[1].strip(), int(p[2])
            if c in ("", "---", "0"):
                continue
            src = "chr" + c
            hits = lo.convert_coordinate(src, pos - 1)
            if not hits or hits[0][0] != src:
                dropped += 1
                continue
            out[pid] = (c, hits[0][1] + 1)
    return out, dropped


def main():
    matrix, annot, chain, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    pos, dropped = load_positions(annot, chain)
    print("markers placed on GRCh37: %d, dropped by liftover: %d" % (len(pos), dropped), flush=True)

    gsms, titles, chars = [], [], []
    handles = None
    written = 0
    with gzip.open(matrix, "rt") as f:
        in_table = False
        for line in f:
            if not in_table:
                cells = [x.strip('"') for x in line.rstrip("\n").split("\t")]
                if line.startswith("!Sample_geo_accession"):
                    gsms = cells[1:]
                elif line.startswith("!Sample_title"):
                    titles = cells[1:]
                elif line.startswith("!Sample_characteristics_ch1"):
                    chars.append(cells[1:])
                elif line.startswith("!series_matrix_table_begin"):
                    in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            p = line.rstrip("\n").split("\t")
            pid = p[0].strip('"')
            if pid == "ID_REF":
                handles = [open(os.path.join(outdir, g + ".probes"), "w") for g in gsms]
                for h in handles:
                    h.write(HEAD + "\n")
                continue
            hit = pos.get(pid)
            if hit is None or handles is None:
                continue
            pre = "%s\t%s\t%d\t\t\t\t" % (pid, hit[0], hit[1])
            for i, v in enumerate(p[1:len(handles) + 1]):
                handles[i].write(pre + CALLS.get(v.strip('"'), "NC") + "\t1\n")
            written += 1
    for h in handles or []:
        h.close()

    truth = []
    for i, g in enumerate(gsms):
        row = {"gsm": g, "title": titles[i] if i < len(titles) else ""}
        for c in chars:
            if i < len(c) and ":" in c[i]:
                k, v = c[i].split(":", 1)
                row[k.strip()] = v.strip()
        truth.append(row)
    with open(os.path.join(outdir, "truth.json"), "w") as w:
        json.dump(truth, w, indent=1)
    print("samples: %d, markers written per sample: %d" % (len(gsms), written), flush=True)


if __name__ == "__main__":
    main()

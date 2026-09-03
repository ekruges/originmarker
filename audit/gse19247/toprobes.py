"""Emit `.probes` files from GSE19247, so the actual tool can be run on a second public study.

WHY THIS EXISTS. The scripts beside it validate the DOSAGE CHANNEL on this series by
reimplementing the arithmetic in numpy. That answered one question and left the rest untouched:
nothing here has ever gone through `scoreSample`, so the stage inference, the zygosity call, the
segment scan, the parental channels and every guard have only ever seen one study and one platform.
This converter closes that. It writes the same tab-separated shape the tool already reads, so every
surface can be attacked on material from a different laboratory, a different array chemistry and a
different decade.

WHAT THE TRUTH IS, from the GEO `source_name` field and nothing else:

  lymphoblast from Coriell family 1990   karyotype-confirmed TRISOMY 21
  sperm cells from Gene Security ...     HAPLOID by biology, so one parental complement
  blood cell from Gene Security 231      diploid somatic, and the same family as that family's
                                         sperm, so a genotype reconstructed from the sperm can be
                                         checked against a real array of the man himself
  Day 3 cleavage stage embryo NNNN       several single cells of ONE embryo, so one genome read
                                         many times

INTENSITY. The tool reads a log2 ratio, which this platform does not publish. It is built here the
standard way: each marker's total intensity against the median of that marker across every array in
the series, then the whole sample recentred on its own median so the array's own loading drops out.
That is a real LRR and it is the quantity the intensity channel was designed for. It is NOT the
publisher's own normalisation, and any number that comes out of it should be read as this
conversion's as much as the tool's.

CALLS COME FROM cluster.py AND geno.py UNCHANGED. Nothing about genotyping is re-decided here.

Run, after cluster.py has written its clusters:
  ./.venv/bin/python toprobes.py --platform 6985 --out ~/gse19247/probes
"""
import argparse
import json
import os
import sys

import numpy as np

from geno import call as geno_call
from cluster import sample_theta, MIN_R_FRAC

CALL_NAME = {0: 'AA', 1: 'AB', 2: 'BB', -1: 'NC'}
HEADER = 'probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset'


def idat_paths(rec, root):
    """The Grn and Red file for one sample, in that order."""
    grn = [f for f in rec['idats'] if 'Grn' in f]
    red = [f for f in rec['idats'] if 'Red' in f]
    if not grn or not red:
        return None
    g, r = os.path.join(root, grn[0]), os.path.join(root, red[0])
    return (g, r) if os.path.exists(g) and os.path.exists(r) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--platform', default='6985', choices=['6985', '8855'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--limit', type=int, default=0, help='stop after N samples, for a smoke run')
    a = ap.parse_args()

    idat_root = 'idat' if a.platform == '6985' else 'idat8'
    index = f'gpl{a.platform}_samples.json'
    recs = json.load(open(index))
    recs = [r for r in recs if idat_paths(r, idat_root)]
    if a.limit:
        recs = recs[:a.limit]
    if not recs:
        sys.exit(f'no readable IDAT pairs under {idat_root}')

    # Cluster positions and the address-to-marker map. Both from the scripts already here.
    z = np.load(f'clusters_gpl{a.platform}.npz', allow_pickle=True)
    addr, chrom, pos, mu = z['addr'], z['chrom'], z['pos'], z['mu']
    print(f'{len(recs)} samples, {len(addr)} markers, platform GPL{a.platform}', flush=True)

    # ---- PASS ONE: the per-marker reference intensity.
    #
    # The median across every array in the series. A per-marker reference is what makes the ratio a
    # RATIO: probe affinity varies by orders of magnitude between markers and would otherwise sit in
    # the signal, swamping the copy-number differences this is for.
    tot_sum = None
    n_seen = 0
    mats = []
    for i, r in enumerate(recs):
        g, red = idat_paths(r, idat_root)
        # sample_theta reads the IDATs itself: it takes PATHS, not arrays.
        _, tot = sample_theta(g, red, addr)
        mats.append(tot.astype(np.float32))
        n_seen += 1
        if (i + 1) % 25 == 0:
            print(f'  intensity pass {i + 1}/{len(recs)}', flush=True)
    stack = np.vstack(mats)
    ref = np.median(stack, axis=0)
    ref[ref <= 0] = np.nan
    print(f'  reference built from {n_seen} arrays', flush=True)

    os.makedirs(a.out, exist_ok=True)
    written = 0
    for i, r in enumerate(recs):
        g, red = idat_paths(r, idat_root)
        # Genotype and BAF from the shared caller, untouched. It takes paths too.
        # geno.call returns a third value, the live mask, which this converter does not need.
        gt, baf, _live = geno_call(g, red, addr, mu)
        tot = mats[i]
        with np.errstate(divide='ignore', invalid='ignore'):
            lrr = np.log2(tot / ref)
        # Recentre on the array's own median, so how much DNA went on the chip drops out and what
        # is left is the difference between this marker and the rest of this array.
        finite = np.isfinite(lrr)
        if finite.any():
            lrr = lrr - np.median(lrr[finite])

        out = os.path.join(a.out, f"{r['gsm']}.probes")
        with open(out, 'w') as fh:
            fh.write(HEADER + '\n')
            for k in range(len(addr)):
                c = str(chrom[k])
                if not c or c in ('0', 'nan'):
                    continue
                l = lrr[k]
                b = baf[k]
                fh.write(
                    f'{addr[k]}\t{c}\t{int(pos[k])}\t'
                    f'{"" if not np.isfinite(l) else f"{l:.4f}"}\t'
                    f'{"" if not np.isfinite(b) else f"{b:.4f}"}\t'
                    f'\t{CALL_NAME[int(gt[k])]}\t1\n')
        written += 1
        if written % 25 == 0:
            print(f'  written {written}/{len(recs)}', flush=True)

    # The truth table, so a harness never has to parse a filename.
    truth = [{'gsm': r['gsm'], 'title': r['title'], 'src': r['src']} for r in recs]
    with open(os.path.join(a.out, 'truth.json'), 'w') as fh:
        json.dump(truth, fh, indent=1)
    print(f'wrote {written} probes files and truth.json to {a.out}', flush=True)


if __name__ == '__main__':
    main()

"""Fill the log2R column of GPL3718 probes files from the series' own copy-number exports.

The 250K Nsp series in GEO carry genotype calls in the series matrix, and some of them also ship a
per-sample CNAT copy-number export with one Log2Ratio per marker. convert.py writes an empty
intensity column because a matrix carries none; where the exports exist this joins them by
probeset, so the intensity channels have something to read. Positions stay as convert.py lifted
them and only the intensity column changes.

Run: python cn_merge.py <probes_dir> <cn_dir>
"""
import gzip
import os
import re
import sys


def log2_by_probe(path):
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > 3 and p[0].startswith("SNP_A-"):
                out[p[0]] = p[3]
    return out


def main():
    probes_dir, cn_dir = sys.argv[1:3]
    cn = {}
    for f in os.listdir(cn_dir):
        m = re.match(r"(GSM\d+)_.*\.cn\.txt\.gz$", f)
        if m and "_1Mb" not in f:
            cn[m.group(1)] = os.path.join(cn_dir, f)
    done, missing = 0, []
    for f in sorted(os.listdir(probes_dir)):
        if not f.endswith(".probes"):
            continue
        gsm = f[: -len(".probes")]
        src = cn.get(gsm)
        if not src:
            missing.append(gsm)
            continue
        vals = log2_by_probe(src)
        path = os.path.join(probes_dir, f)
        tmp = path + ".tmp"
        filled = rows = 0
        with open(path) as r, open(tmp, "w") as w:
            w.write(r.readline())
            for line in r:
                p = line.rstrip("\n").split("\t")
                rows += 1
                v = vals.get(p[0])
                if v:
                    p[3] = v
                    filled += 1
                w.write("\t".join(p) + "\n")
        os.replace(tmp, path)
        done += 1
        print("%s: %d of %d markers carry intensity" % (gsm, filled, rows), flush=True)
    print("filled %d samples; no copy-number export for: %s"
          % (done, ", ".join(missing) or "none"), flush=True)


if __name__ == "__main__":
    main()

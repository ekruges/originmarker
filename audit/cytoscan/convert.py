"""Affymetrix CytoScan CYCHP files into the tab-separated shape OriginMarker reads.

A CYCHP is a Command Console Generic Data File: a binary container of typed data sets. The three
this needs are ProbeSets/CopyNumber (Log2Ratio per marker), Genotyping/Calls (the genotype and the
raw A and B signals) and ProbeSets/AllelicData. Call codes are the Affymetrix ones, 6 AA, 7 BB,
8 AB, 11 no call, confirmed on this data by mean contrast: -1.96, +1.93 and -0.02.

WHY B-ALLELE FREQUENCY IS COMPUTED RATHER THAN READ. The BAF column of these files is empty on
every marker. The A and B signals are not, so the allele fraction b = B / (A + B) is recoverable,
but a raw fraction carries each probe's own bias and puts homozygotes near 0.2 and 0.8 instead of 0
and 1. GenomeStudio removes that bias by clustering each marker across a cohort, and this does the
same: per marker, the median raw fraction of the samples called AA, AB and BB anchors 0, 0.5 and 1,
and every sample is mapped onto those anchors. A marker whose cohort never shows a class is left
without the corresponding anchor and interpolates from the ones it has.

Three stages, because the cluster anchors need every sample before any sample can be written:
    stage1   <cychp> <outdir>            one compact binary per sample, and the marker list once
    clusters <outdir>                    per-marker anchors from every binary present
    probes   <outdir> <gsm> <out.probes> one sample, written against the anchors

Run under any Python 3, no dependencies.
"""
import array
import os
import struct
import sys

HEAD = "probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset"
CALL = {6: "AA", 7: "BB", 8: "AB"}
# Fallback only. The numeric chromosome codes are NOT a fixed convention: a CytoScan 750K file
# numbers X as 24 and Y as 25, with no 23 at all, so a hardcoded 23=X map writes every X marker
# out as chrY. Every file states its own names in Chromosomes/Summary and that is what stage1
# records; this map is used only for a staged directory written before it did.
CHROM = {i: str(i) for i in range(1, 23)}
CHROM.update({23: "X", 24: "Y", 25: "MT"})


def chrom_names(f, ds):
    """Code to display name, as the file itself states them."""
    out = {}
    key = ("Chromosomes", "Summary")
    if key in ds:
        for row in rows(f, ds[key]):
            if len(row) >= 2 and str(row[1]).strip():
                out[int(row[0])] = str(row[1]).strip()
    return out


class Reader:
    def __init__(s, f):
        s.f = f

    def u8(s):
        return struct.unpack(">B", s.f.read(1))[0]

    def i32(s):
        return struct.unpack(">i", s.f.read(4))[0]

    def u32(s):
        return struct.unpack(">I", s.f.read(4))[0]

    def astr(s):
        n = s.i32()
        return s.f.read(n).decode("ascii", "replace") if n > 0 else ""

    def wstr(s):
        n = s.i32()
        return s.f.read(n * 2).decode("utf-16-be", "replace") if n > 0 else ""

    def blob(s):
        return s.f.read(s.i32())


def _header(r):
    r.astr(), r.astr(), r.wstr(), r.wstr()
    for _ in range(r.i32()):
        r.wstr(), r.blob(), r.wstr()
    for _ in range(r.i32()):
        _header(r)


def datasets(path):
    """{(group, set): descriptor} for one file, with the file handle to read rows from."""
    f = open(path, "rb")
    r = Reader(f)
    r.u8(), r.u8()
    ngroups, pos = r.i32(), r.u32()
    _header(r)
    out = {}
    for _ in range(ngroups):
        f.seek(pos)
        nxt, first_set, nsets = r.u32(), r.u32(), r.i32()
        gname = r.wstr()
        sp = first_set
        for _ in range(nsets):
            f.seek(sp)
            first_el, next_set = r.u32(), r.u32()
            sname = r.wstr()
            for _ in range(r.i32()):
                r.wstr(), r.blob(), r.wstr()
            cols = []
            for _ in range(r.u32()):
                cols.append((r.wstr(), r.u8(), r.i32()))
            out[(gname, sname)] = {"cols": cols, "rows": r.u32(), "at": first_el}
            sp = next_set
        pos = nxt
    return f, out


FMT = {1: "B", 2: "h", 3: "H", 4: "i", 5: "I", 6: "f", 0: "b"}


def rows(f, ds):
    """Every row of a data set, as tuples. Columns are fixed width, strings included.

    A string column is its own little record: a four-byte length, then that many characters, then
    padding out to the column's declared width. Reading the field as text and trimming at the first
    NUL loses the name entirely, because the length prefix starts with NULs.
    """
    fmt = ">"
    types = []
    for _, ct, cs in ds["cols"]:
        fmt += FMT[ct] if ct in FMT else "%ds" % cs
        types.append(ct)
    st = struct.Struct(fmt)
    f.seek(ds["at"])
    buf = f.read(st.size * ds["rows"])
    for i in range(ds["rows"]):
        vals = st.unpack_from(buf, i * st.size)
        out = []
        for ct, v in zip(types, vals):
            if not isinstance(v, bytes):
                out.append(v)
                continue
            ln = struct.unpack(">i", v[:4])[0]
            wide = ct == 8
            out.append(v[4:4 + ln * (2 if wide else 1)]
                       .decode("utf-16-be" if wide else "ascii", "replace"))
        yield out


def stage1(cychp, outdir):
    """One sample into a compact binary: raw allele fraction, call code and log2 ratio per marker."""
    os.makedirs(outdir, exist_ok=True)
    gsm = os.path.basename(cychp).split("_")[0].split(".")[0]
    f, ds = datasets(cychp)
    names_path = os.path.join(outdir, "chroms.tsv")
    if not os.path.exists(names_path):
        with open(names_path, "w") as w:
            for code, name in sorted(chrom_names(f, ds).items()):
                w.write("%d\t%s\n" % (code, name))
    log2 = {}
    # Every copy number probe, not only the genotyped ones: a chromosome the genotyping algorithm
    # never calls is invisible in the join below, and chrY is exactly that chromosome on this
    # platform. It is what separates one X in a male from one X in a monosomy, so its probes are
    # carried as intensity with no genotype rather than dropped.
    cn = {}
    for name, chrom, pos, l2, _w, _s, _n in rows(f, ds[("ProbeSets", "CopyNumber")]):
        cn[name] = (chrom, pos, l2)
        if name.startswith("S"):
            log2[name] = l2
    names, frac, calls, l2s, chrs, poss = [], array.array("f"), array.array("B"), array.array("f"), array.array("B"), array.array("I")
    allelic = {}
    for name, chrom, pos, _ad, _baf in rows(f, ds[("ProbeSets", "AllelicData")]):
        allelic[name] = (chrom, pos)
    for _i, name, call, _c, _fc, a, b, _ss, _ct in rows(f, ds[("Genotyping", "Calls")]):
        hit = allelic.get(name)
        if hit is None:
            continue
        tot = (a or 0) + (b or 0)
        names.append(name)
        frac.append(b / tot if tot > 0 else float("nan"))
        calls.append(call if call in CALL else 11)
        l2s.append(log2.get(name, float("nan")))
        chrs.append(hit[0])
        poss.append(hit[1])
    genotyped = set(chrs)
    for name, (chrom, pos, l2) in cn.items():
        if chrom in genotyped or name in allelic:
            continue
        names.append(name)
        frac.append(float("nan"))
        calls.append(11)
        l2s.append(l2)
        chrs.append(chrom)
        poss.append(pos)
    f.close()
    markers = os.path.join(outdir, "markers.tsv")
    if not os.path.exists(markers):
        with open(markers, "w") as w:
            for i, n in enumerate(names):
                w.write("%s\t%d\t%d\n" % (n, chrs[i], poss[i]))
    with open(os.path.join(outdir, gsm + ".bin"), "wb") as w:
        w.write(struct.pack(">I", len(names)))
        w.write(frac.tobytes()); w.write(calls.tobytes()); w.write(l2s.tobytes())
    print("%s: %d markers" % (gsm, len(names)), flush=True)


def load_bin(path):
    with open(path, "rb") as f:
        n = struct.unpack(">I", f.read(4))[0]
        frac = array.array("f"); frac.frombytes(f.read(4 * n))
        calls = array.array("B"); calls.frombytes(f.read(n))
        l2 = array.array("f"); l2.frombytes(f.read(4 * n))
    return n, frac, calls, l2


def clusters(outdir):
    """Per-marker anchors: the cohort's mean raw fraction for each genotype class.

    A mean rather than a median, because a median needs every sample's value for every marker held
    at once, which is 38 million floats on this cohort. Running sums are three numbers per marker
    per class, and a wrong genotype call moves an anchor by at most one sample in the count.
    """
    files = sorted(x for x in os.listdir(outdir) if x.endswith(".bin"))
    n = load_bin(os.path.join(outdir, files[0]))[0]
    sums = array.array("d", bytes(8 * 3 * n))
    counts = array.array("I", bytes(4 * 3 * n))
    used = 0
    for fn in files:
        m, frac, calls, _l2 = load_bin(os.path.join(outdir, fn))
        if m != n:
            print("skipped %s: %d markers against %d" % (fn, m, n), flush=True)
            continue
        used += 1
        for i in range(n):
            c = calls[i]
            v = frac[i]
            if v == v and c in (6, 7, 8):
                k = 3 * i + (0 if c == 6 else (2 if c == 7 else 1))
                sums[k] += v
                counts[k] += 1
    out = array.array("f")
    for k in range(3 * n):
        out.append(sums[k] / counts[k] if counts[k] else float("nan"))
    files = files[:used]
    with open(os.path.join(outdir, "clusters.bin"), "wb") as w:
        w.write(struct.pack(">I", n))
        w.write(out.tobytes())
    print("anchors from %d samples over %d markers" % (len(files), n), flush=True)


def probes(outdir, gsm, out_path):
    """One sample written against the cohort anchors."""
    with open(os.path.join(outdir, "clusters.bin"), "rb") as f:
        n = struct.unpack(">I", f.read(4))[0]
        anchors = array.array("f"); anchors.frombytes(f.read(12 * n))
    named = dict(CHROM)
    names_path = os.path.join(outdir, "chroms.tsv")
    if os.path.exists(names_path):
        with open(names_path) as f:
            named = {}
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) == 2:
                    named[int(p[0])] = p[1]
    marks = []
    with open(os.path.join(outdir, "markers.tsv")) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            marks.append((p[0], named.get(int(p[1]), ""), p[2]))
    m, frac, calls, l2 = load_bin(os.path.join(outdir, gsm + ".bin"))
    written = 0
    with open(out_path, "w") as w:
        w.write(HEAD + "\n")
        for i in range(min(m, n, len(marks))):
            name, chrom, pos = marks[i]
            if not chrom:
                continue
            aa, ab, bb = anchors[3 * i], anchors[3 * i + 1], anchors[3 * i + 2]
            v = frac[i]
            baf = ""
            if v == v and ab == ab:
                # Piecewise linear between the anchors this marker has, clamped to the unit range.
                if v <= ab and aa == aa and ab > aa:
                    baf = 0.5 * (v - aa) / (ab - aa)
                elif v > ab and bb == bb and bb > ab:
                    baf = 0.5 + 0.5 * (v - ab) / (bb - ab)
                if baf != "":
                    baf = "%.4f" % min(1.0, max(0.0, baf))
            g = CALL.get(calls[i], "NC")
            v2 = l2[i]
            w.write("%s\t%s\t%s\t%s\t%s\t\t%s\t1\n"
                    % (name, chrom, pos, ("%.4f" % v2) if v2 == v2 else "", baf, g))
            written += 1
    print("%s: %d markers written" % (gsm, written), flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "stage1":
        stage1(sys.argv[2], sys.argv[3])
    elif cmd == "clusters":
        clusters(sys.argv[2])
    elif cmd == "probes":
        probes(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        raise SystemExit("stage1 | clusters | probes")

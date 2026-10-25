"""Turn Table S3's karyotype strings into GRCh37 coordinates, from the UCSC cytoBand track.

The audits were joining a called segment to a hand-called event by CHROMOSOME alone, so a
segment anywhere on a chromosome counted as detecting an event anywhere else on it. One such
case is in this corpus and it was published as a detection.

Every region here is derived from `ucsc-hg19-cytoband.tsv`. Nothing is estimated from chromosome
length, and a string this module cannot parse returns None rather than a guess.

  usage: python3 -I audit/karyotype_regions.py          # self-check over Table S3
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def bands() -> dict[str, list[tuple[int, int, str]]]:
    """Per chromosome, (start, end, band) for every cytoband, in order."""
    out: dict[str, list[tuple[int, int, str]]] = {}
    for line in (HERE / 'ucsc-hg19-cytoband.tsv').read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        chrom, start, end, band, _stain = line.split('\t')
        c = chrom.removeprefix('chr')
        out.setdefault(c, []).append((int(start), int(end), band))
    return out


BANDS = bands()


def arm_span(chrom: str, arm: str) -> tuple[int, int] | None:
    """The whole p or q arm."""
    hits = [(s, e) for s, e, b in BANDS.get(chrom, []) if b.startswith(arm)]
    return (min(s for s, _ in hits), max(e for _, e in hits)) if hits else None


def band_start(chrom: str, band: str) -> int | None:
    """Start of a named band. '16q23.1' -> the start of q23.1 on chr16."""
    for s, _e, b in BANDS.get(chrom, []):
        if b == band:
            return s
    # A parent band names the span of its sub-bands: 6q21 covers 6q21 alone here, but a table may
    # name a band that UCSC only carries split. Fall back to the first sub-band of that prefix.
    subs = [s for s, _e, b in BANDS.get(chrom, []) if b.startswith(band)]
    return min(subs) if subs else None


def chrom_end(chrom: str) -> int:
    return max(e for _s, e, _b in BANDS[chrom])


def region(chrom: str, text: str) -> tuple[int, int] | None:
    """The coordinates a karyotype clause names, or None if it names no located region.

    Handles the forms Table S3 actually uses: `del 16q23.1-tel`, `del(6p)`, `del12p`,
    `del (6q21-tel)`, `dup 16q21-tel`, `del8p23.1`, `del7q`. A whole-chromosome call
    (`monosomy 5`, `nullisomy`, `trisomy 15`, `normal`) returns None: it is not a segment.
    """
    t = text.replace(' ', '')
    # band-to-telomere, e.g. 16q23.1-tel or 6q21-tel
    m = re.search(r'(\d+)([pq][\d.]+)-tel', t, re.I)
    if m and m.group(1) == chrom:
        s = band_start(chrom, m.group(2).lower())
        return (s, chrom_end(chrom)) if s is not None else None
    # a single named band, e.g. 8p23.1
    m = re.search(r'(?:del|dup)\(?(\d+)([pq][\d.]+)\)?', t, re.I)
    if m and m.group(1) == chrom:
        b = m.group(2).lower()
        for s, e, nm in BANDS.get(chrom, []):
            if nm == b:
                return (s, e)
        subs = [(s, e) for s, e, nm in BANDS.get(chrom, []) if nm.startswith(b)]
        if subs:
            return (min(s for s, _ in subs), max(e for _, e in subs))
        return None
    # a whole arm, e.g. del(6q), del12p, del7q
    m = re.search(r'(?:del|dup)\(?(\d+)([pq])\)?(?![\d.])', t, re.I)
    if m and m.group(1) == chrom:
        return arm_span(chrom, m.group(2).lower())
    # an arm with the chromosome implied by the row, e.g. del(6q) already covered; bare (6q21-tel)
    m = re.search(r'\((\d+)([pq][\d.]*)-tel\)', t, re.I)
    if m and m.group(1) == chrom:
        s = band_start(chrom, m.group(2).lower())
        return (s, chrom_end(chrom)) if s is not None else None
    return None


def clauses(text: str) -> list[str]:
    """A karyotype cell can name two events: 'del(6p)(pat),dup(6q)(pat)'."""
    return [c for c in re.split(r',(?![^()]*\))', text) if re.search(r'del|dup', c, re.I)]


def main() -> int:
    rows = [l.rstrip('\n').split('\t')
            for l in (HERE / 'zuccaro-table-s3.tsv').read_text().splitlines()
            if l.strip() and not l.startswith('#')]
    hdr, rows = rows[0], rows[1:]
    gi, ci, mi = (hdr.index(x) for x in ('gsm', 'chrom', 'manual'))
    parsed = unparsed = 0
    print(f"{'cell':12} {'chr':>3}  {'region (Mb)':>18}  {'%chrom':>6}  clause")
    for r in rows:
        if not re.search(r'del|dup', r[mi], re.I):
            continue
        for cl in clauses(r[mi]):
            reg = region(r[ci], cl)
            if reg is None:
                unparsed += 1
                print(f"{r[gi]:12} {r[ci]:>3}  {'UNPARSED':>18}  {'':>6}  {cl.strip()}")
                continue
            parsed += 1
            span = (reg[1] - reg[0]) / chrom_end(r[ci])
            print(f"{r[gi]:12} {r[ci]:>3}  "
                  f"{reg[0] / 1e6:7.1f}-{reg[1] / 1e6:<8.1f}  {span:6.1%}  {cl.strip()}")
    print(f"\nparsed {parsed} clauses, unparsed {unparsed}")
    assert unparsed == 0, 'every segmental clause in Table S3 must resolve to coordinates'
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

"""structural - the evidence a genotyping array structurally cannot carry.

An array reports at marker positions. Two of the mechanisms this tool has to weigh leave no
trace there at all:

    H3b  a large deletion at the edited site. The array's nearest markers may sit outside it, and
         inside it a hemizygous call is indistinguishable from a homozygous one.
    H3f  an insertion, or a balanced rearrangement. Neither changes copy number, so intensity is
         flat and genotypes are unremarkable, and the event is invisible by construction.

Reads see both. This module reads them from a STRUCTURAL VARIANT VCF rather than from a BAM, and
the distinction is worth stating: an SV VCF is the output of a caller that already did the read
analysis, and it carries the split-read and paired-read counts behind each call. Parsing BAM here
would mean reimplementing that caller, with no BAM on hand to check the result against.

What this module will not do is turn absence of an SV call into evidence of absence. A caller
that was never run, or was run with a size threshold above the event, produces exactly the same
empty result as a genome with no structural variant, and those are reported differently.

Self-check:  python -m structural
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional, Sequence

#: Types this tool reasons about. Anything else is read and reported but not interpreted.
DELETION = ("DEL", "CNV_DEL", "CON")
INSERTION = ("INS", "DUP", "INV", "BND", "CPX")


@dataclass
class StructuralVariant:
    chrom: str
    start: int
    end: int
    svtype: str
    svlen: Optional[int]
    variant_id: str
    filter: str
    genotype: str = "NC"
    #: Split reads and paired reads behind the call, where the caller recorded them.
    split_reads: Optional[int] = None
    paired_reads: Optional[int] = None

    @property
    def span(self) -> int:
        return max(0, self.end - self.start)

    @property
    def passes(self) -> bool:
        return self.filter in ("PASS", ".", "")

    @property
    def carried(self) -> bool:
        """Whether the sample actually has it. A call in the file is not a call in the sample."""
        return self.genotype in ("AB", "BB")


def _info(field_text: str) -> dict[str, str]:
    out = {}
    for part in field_text.split(";"):
        k, _, v = part.partition("=")
        out[k] = v
    return out


def read_sv_vcf(path: str | Path, sample_name: Optional[str] = None) -> list[StructuralVariant]:
    """Structural calls for one sample, from a VCF written by an SV caller."""
    path = Path(path)
    opener = gzip.open if path.name.endswith(".gz") else open
    out: list[StructuralVariant] = []
    col: Optional[int] = None
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = f[9:] if len(f) > 9 else []
                if names:
                    col = 9 + (names.index(sample_name) if sample_name in names else 0)
                continue
            if len(f) < 8:
                continue
            info = _info(f[7])
            svtype = info.get("SVTYPE", "")
            if not svtype:
                continue
            try:
                start = int(f[1])
            except ValueError:
                continue
            try:
                svlen = int(info["SVLEN"].split(",")[0])
            except (KeyError, ValueError):
                svlen = None
            try:
                end = int(info["END"])
            except (KeyError, ValueError):
                # Without END, a deletion's extent is its length and an insertion has none.
                end = start + abs(svlen) if (svlen and svtype in DELETION) else start

            gt, sr, pr = "NC", None, None
            if col is not None and len(f) > col:
                rec = dict(zip(f[8].split(":"), f[col].split(":")))
                g = rec.get("GT", "./.").replace("|", "/")
                a, _, b = g.partition("/")
                if a in ("0", "1") and b in ("0", "1"):
                    gt = {"00": "AA", "01": "AB", "10": "AB", "11": "BB"}[a + b]
                for key, target in (("SR", "sr"), ("PR", "pr")):
                    v = rec.get(key)
                    if v and "," in v:
                        try:
                            val = int(v.split(",")[1])
                        except ValueError:
                            val = None
                        if target == "sr":
                            sr = val
                        else:
                            pr = val
            out.append(StructuralVariant(
                chrom=f[0].removeprefix("chr"), start=start, end=end, svtype=svtype,
                svlen=svlen, variant_id=f[2], filter=f[6], genotype=gt,
                split_reads=sr, paired_reads=pr))
    return out


@dataclass
class StructuralEvidence:
    """What the reads say about one locus, and how far that goes."""

    supplied: bool
    chrom: str
    pos: int
    deletions: list[StructuralVariant] = field(default_factory=list)
    insertions: list[StructuralVariant] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def spans_locus(self) -> bool:
        return bool(self.deletions or self.insertions)

    def summary(self) -> str:
        if not self.supplied:
            return ("No structural calls supplied. A large deletion at the locus and an insertion "
                    "or balanced rearrangement are NOT TESTED: an array cannot see either.")
        if not self.spans_locus:
            return (f"Structural calls supplied and none covers chr{self.chrom}:{self.pos:,}. "
                    "That excludes an event the caller could detect, and not one below its size "
                    "threshold or outside its type set.")
        bits = []
        for v in self.deletions:
            bits.append(f"{v.svtype} {v.span:,} bp spanning the locus, {v.genotype}"
                        + (f", {v.split_reads} split reads" if v.split_reads else ""))
        for v in self.insertions:
            bits.append(f"{v.svtype} at chr{v.chrom}:{v.start:,}, {v.genotype}"
                        + (f", {v.split_reads} split reads" if v.split_reads else ""))
        return "; ".join(bits)


def at_locus(
    svs: Sequence[StructuralVariant], chrom: str, pos: int, *, flank: int = 0,
    supplied: bool = True, carried_only: bool = True,
) -> StructuralEvidence:
    """Structural calls covering a position, split into the two mechanisms they bear on."""
    want = str(chrom).removeprefix("chr")
    ev = StructuralEvidence(supplied=supplied, chrom=want, pos=pos)
    if not supplied:
        return ev
    for v in svs:
        if v.chrom != want or not v.passes:
            continue
        if carried_only and not v.carried:
            continue
        if v.svtype in DELETION and v.start - flank <= pos <= v.end + flank:
            ev.deletions.append(v)
        elif v.svtype in INSERTION and abs(v.start - pos) <= max(flank, 1_000):
            ev.insertions.append(v)
    if svs and not ev.spans_locus:
        ev.notes.append(
            "Absence of a call is not absence of an event. A caller run with a minimum size "
            "above this event, or without the relevant type enabled, gives the same empty "
            "result as a genome with nothing there."
        )
    return ev


if __name__ == "__main__":
    import tempfile

    vcf = "\n".join([
        "##fileformat=VCFv4.2",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "6\t65399000\tdel1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=65402000;SVLEN=-3000\tGT:SR\t0/1:12,7",
        "6\t65400500\tins1\tN\t<INS>\t50\tPASS\tSVTYPE=INS;SVLEN=340\tGT:SR\t0/1:9,4",
        "6\t70000000\tdel2\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=70010000;SVLEN=-10000\tGT\t0/0",
        "6\t65401000\tlow\tN\t<DEL>\t50\tLowQual\tSVTYPE=DEL;END=65402000\tGT\t1/1",
    ])
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "sv.vcf"
        p.write_text(vcf + "\n")
        svs = read_sv_vcf(p)

    assert len(svs) == 4
    d1 = svs[0]
    assert (d1.svtype, d1.span, d1.genotype, d1.split_reads) == ("DEL", 3000, "AB", 7)

    ev = at_locus(svs, "6", 65_400_000)
    assert len(ev.deletions) == 1 and ev.deletions[0].variant_id == "del1"
    assert len(ev.insertions) == 1 and ev.insertions[0].variant_id == "ins1"
    assert "spanning the locus" in ev.summary()

    # A call the sample does not carry, and one that failed filters, are both excluded.
    assert all(v.variant_id != "del2" for v in ev.deletions), "0/0 is not carried"
    assert all(v.variant_id != "low" for v in ev.deletions), "LowQual does not pass"

    clear = at_locus(svs, "6", 90_000_000)
    assert not clear.spans_locus
    assert any("not absence of an event" in n for n in clear.notes)

    none = at_locus([], "6", 65_400_000, supplied=False)
    assert not none.supplied and "NOT TESTED" in none.summary()

    print("structural self-check OK")

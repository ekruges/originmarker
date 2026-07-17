"""primers - candidate FWD/REV pairs for genotyping a marker, with common variants masked.

The novel part is the mask. Markers are SNPs, and the caller already holds every variant
gnomAD reports in the window, so the pool of candidates is also the pool of hazards. A
primer sitting on a common variant fails to bind in exactly the carriers who have it, and a
heterozygote read as a homozygote is allele dropout: the worst error this app can make. So
every variant at or above `mask_maf` is excluded from both primer footprints, and the target
sits inside the product under neither primer.

This module never fetches. Sequence, coordinates and variants arrive as arguments, which
keeps the network in one place (R1) and keeps the design testable offline.

primer3-py is GPLv2 and this repo is Apache 2.0, so it is an OPTIONAL dependency and is not
in requirements.txt: see requirements-primers.txt. Absent, design() returns a result whose
`error` says so and whose pair is None. Never raises for a missing dependency.

The pair is a CANDIDATE (R3). It has not been run against the genome unless a verification
lane has set `insilico_pcr`, and `warnings` says so on every result.

Self-check:  python primers.py
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, NamedTuple, Optional, Sequence

try:
    import primer3
except ImportError:                                  # optional by licence, see module head
    primer3 = None

# Where the primer explanation lives. One route, so a note gains a link by existing.
PRIMER_DOCS = "#/docs/primers"


@dataclass(frozen=True)
class Note:
    """One thing wrong with a pair, at two lengths, from one place.

    A warning has two readers and they cannot share a wording. The table shows a note beside
    every pair, so the same paragraph repeated down the page is a paragraph nobody reads by
    the third row; the PDF is filed and read once, away from any link, so it needs the whole
    thing. Written twice they drift, and the short one is the one on screen: the wrong half
    to have go stale. So both are written here, beside each other, and `docs` is where the
    reader who wants more than either goes.
    """
    code: str
    short: str
    long: str
    docs: str = PRIMER_DOCS

    def __str__(self) -> str:
        # Print is the long form: a filed page cannot follow a link, so the sentence that
        # defers to the docs would defer to nothing.
        return self.long


UNAVAILABLE = Note(
    code="primer3_missing",
    short="Primer design is switched off on this server: primer3 is not installed.",
    long="Primer design is unavailable: primer3 is not installed on this server. It is an "
         "optional GPLv2 dependency and is deliberately absent from the default image "
         "(see requirements-primers.txt). The panel is unaffected; it simply carries no "
         "primers.",
)

# Carried by every pair this module returns, which is every pair before a verification lane
# has looked at it. It states only what this module can know: the pair was designed against
# one template and checked for products in that template and nowhere else. Whether a
# verifier is configured is app/ispcr.py's fact to report, not this module's to guess.
NOT_CHECKED_WARNING = Note(
    code="not_checked",
    short="Not checked against the genome: it may amplify more than one locus.",
    long="NOT CHECKED AGAINST THE GENOME. This pair was designed against a single reference "
         "template and has not been tested for products anywhere else in GRCh38. It may "
         "amplify more than one locus, which cannot be genotyped. Verify it before ordering.",
)

# primer3's built-in cap is 36, but its manual scopes the melting temperature model to 35.
# A Tm past that is computed, plausible and out of model, which is the one kind of answer
# this app refuses.
SIZE_CAP = 35

# The widest product a design may ask for. Not primer3's limit: it is the width of the
# question app/ispcr.py can put to UCSC (ispcr.REPORT_MAX_BP, 4000, hgPcr's own default for
# Max Product Size). A product past that cannot be REPORTED by UCSC, and hgPcr reports
# nothing rather than reporting it truncated, so the pair's own on-target product goes
# missing and reads as "found no product, do not order": our own request accusing a good
# pair. The number lives here rather than there because this module must not import a module
# that fetches; ispcr's self-check asserts its own width still clears this one.
PRODUCT_CAP = 3000

# A product needs room to sit inside the template AND room for primer3 to search either
# side of it. Pinning max product to the template length leaves no search space.
FLANK_MARGIN = 100

# libprimer3's PR_MAX_INTERVAL_ARRAY. Compile-time, so it cannot be raised from here.
MAX_MASK_INTERVALS = 200


def available() -> bool:
    """True if primer3 is importable. Design is optional; the app works without it."""
    return primer3 is not None


class MaskSite(NamedTuple):
    """One variant to keep out of the primer footprints, in GRCh38 plus-strand VCF terms.

    The caller decides what goes in this list, and its reach is the mask's reach: a list
    built from SNPs alone cannot mask an indel, and a common indel under a primer drops an
    allele just as surely as a SNP does.
    """
    pos: int            # 1-based GRCh38
    ref: str
    alt: str
    maf: float

    @property
    def span(self) -> int:
        """Reference bases this variant occupies. len(ref), so an insertion covers 1."""
        return len(self.ref)


@dataclass(frozen=True)
class PrimerSettings:
    """Every knob the design exposes, with its defaults. Overrides arrive as a replace().

    A Tm means nothing without the reaction it was computed for, so the salt, dNTP and DNA
    concentrations are stated here rather than inherited: primer3 ships two default sets
    that disagree on them, and an inherited default is not a reproducible one.
    """
    # Tm. Length floats to reach it: no 20-mer at 40-60% GC comes near 69 C.
    opt_tm: float = 69.0
    min_tm: float = 67.0
    max_tm: float = 71.0
    # primer3's own default is 100.0, which is no constraint at all. A pair several degrees
    # apart has no annealing temperature that suits both.
    max_pair_diff_tm: float = 5.0
    # Size.
    min_size: int = 20
    opt_size: int = 26
    max_size: int = SIZE_CAP
    # GC.
    min_gc: float = 40.0
    max_gc: float = 60.0
    gc_clamp: int = 1
    # Reaction conditions. These four equal primer3's defaults and are stated anyway.
    salt_monovalent: float = 50.0
    salt_divalent: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0
    # Product. Coupled to the template the caller must supply: see flank_needed().
    min_product: int = 250
    max_product: int = 600
    # Structure.
    max_poly_x: int = 4
    max_self_any_th: float = 47.0
    max_self_end_th: float = 47.0
    max_hairpin_th: float = 47.0
    # An N is an unknown base and a primer over one is a fabricated primer (R1).
    max_ns_accepted: int = 0
    tm_formula: int = 1                 # SantaLucia
    salt_corrections: int = 1           # SantaLucia
    # Lower than the marker floor on purpose. A marker needs heterozygosity to be
    # informative; dropout does not care, and a 1% variant under a primer still drops the
    # allele in the carriers who have it.
    mask_maf: float = 0.01
    # Keeps the target off the primers by a margin rather than by luck, and clear of the
    # unreadable trace immediately downstream of a sequencing primer.
    target_pad: int = 50

    def __post_init__(self):
        if not 0.0 < self.min_tm <= self.opt_tm <= self.max_tm:
            raise ValueError(f"Tm window out of order: {self.min_tm}/{self.opt_tm}/{self.max_tm}")
        if not 1 <= self.min_size <= self.opt_size <= self.max_size:
            raise ValueError(f"size window out of order: "
                             f"{self.min_size}/{self.opt_size}/{self.max_size}")
        if self.max_size > SIZE_CAP:
            raise ValueError(f"max_size {self.max_size} exceeds {SIZE_CAP}, the longest oligo "
                             f"primer3's Tm model covers. Past it the Tm is still computed "
                             f"and no longer valid.")
        if not 0.0 <= self.min_gc <= self.max_gc <= 100.0:
            raise ValueError(f"GC window out of order: {self.min_gc}/{self.max_gc}")
        if not 0 < self.min_product <= self.max_product:
            raise ValueError(f"product window out of order: {self.min_product}/{self.max_product}")
        if self.max_product > PRODUCT_CAP:
            raise ValueError(f"max_product {self.max_product} exceeds {PRODUCT_CAP}, the "
                             f"widest product the verification lane can ask UCSC about. A "
                             f"pair designed past it cannot have its own product reported, "
                             f"and that silence classifies as DANGEROUS.")
        if not 0.0 <= self.mask_maf < 0.5:
            raise ValueError(f"mask_maf out of range: {self.mask_maf}")
        if self.target_pad < 0:
            raise ValueError(f"target_pad out of range: {self.target_pad}")
        # A product must be able to hold the padded target and a primer at each end (R5).
        floor = 2 * self.min_size + 2 * self.target_pad + 1
        if self.max_product < floor:
            raise ValueError(f"max_product {self.max_product} cannot hold a {self.target_pad}bp"
                             f"-padded target between two {self.min_size}bp primers "
                             f"(needs >= {floor})")


DEFAULTS = PrimerSettings()


def flank_needed(s: PrimerSettings = DEFAULTS) -> int:
    """Sequence the caller must fetch each side of the marker, so an inclusive fetch of
    marker_pos +/- this is the template design() wants. Raising max_product without
    widening the fetch starves the search; design() enforces that rather than trust it.
    """
    return s.max_product // 2 + FLANK_MARGIN


@dataclass(frozen=True)
class Primer:
    """One oligo. Only ever reached through a PrimerResult, which carries the warnings."""
    seq: str
    pos: int            # 1-based GRCh38 of its leftmost template base, plus strand
    idx: int            # 0-based offset into the template
    length: int
    tm: float
    gc: float

    @property
    def end_idx(self) -> int:
        """0-based offset of its rightmost template base, inclusive."""
        return self.idx + self.length - 1


@dataclass(frozen=True)
class PrimerResult:
    """A pair with its warnings welded on, or a failure that says why.

    `fwd`/`rev` are None whenever `error` is set, and `warnings` is never empty for a pair:
    a caller cannot render a primer from this without also holding what is wrong with it.
    """
    fwd: Optional[Primer] = None
    rev: Optional[Primer] = None
    product_size: Optional[int] = None
    product_start: Optional[int] = None          # 1-based GRCh38, FWD's 5' base
    masked: tuple = ()                           # tuple[MaskSite], what was kept clear
    mask_note: Optional[Note] = None
    warnings: tuple = ()
    # Multi-valued, like ensembl_pos_check and hotspot_between: "not_checked" is not
    # "passed". A two-state flag renders an unverified pair as a clean one. This module
    # only ever sets the default; app/ispcr.py owns the vocabulary and the other states,
    # and importing it here would drag HTTP into a module that must not fetch.
    insilico_pcr: str = "not_checked"
    error: Optional[str] = None

    @property
    def ok(self) -> bool:
        return self.error is None

    @property
    def product_end(self) -> Optional[int]:
        """1-based GRCh38 of REV's 5' base, inclusive. The amplicon is start..end."""
        if self.product_start is None or self.product_size is None:
            return None
        return self.product_start + self.product_size - 1


def _merge(spans: Iterable[tuple[int, int]]) -> list[list[int]]:
    """Overlapping and abutting [start, length] intervals combined, sorted by start.

    Merging is not a nicety: adjacent variants are common and the interval array is capped,
    so unmerged input hits the cap on windows that merge to a handful of ranges.
    """
    out: list[list[int]] = []
    for start, length in sorted(spans):
        if out and start <= out[-1][0] + out[-1][1]:
            end = max(out[-1][0] + out[-1][1], start + length)
            out[-1][1] = end - out[-1][0]
        else:
            out.append([start, length])
    return out


def design(template: str, region_start: int, marker_pos: int, marker_ref: str,
           snps: Sequence[MaskSite] = (), settings: PrimerSettings = DEFAULTS) -> PrimerResult:
    """Design a FWD/REV pair around `marker_pos`, clear of every masked variant.

    `template` is GRCh38 plus-strand sequence whose first base is `region_start` (1-based);
    `snps` is the caller's variant list, and its reach is the mask's reach. Failure returns
    a result with `error` set: nothing here relaxes a constraint to force a pair, because a
    pair outside the stated Tm or over a common variant is worse than no pair at all.
    """
    if not available():
        # `error` is prose, not a Note: it is the one field a caller renders whole, and every
        # other thing that sets it is a sentence primer3 handed back about this one marker.
        return PrimerResult(error=UNAVAILABLE.long, warnings=(UNAVAILABLE,))

    seq = template.upper()
    bad = set(seq) - set("ACGTN")
    if bad:
        raise ValueError(f"template is not DNA: unexpected {sorted(bad)!r} (R7)")
    idx = marker_pos - region_start
    if not 0 <= idx < len(seq):
        raise ValueError(f"marker at {marker_pos:,} is outside the template "
                         f"{region_start:,}-{region_start + len(seq) - 1:,}")

    # The whole safety story of this module. A mismatch means the template is not what the
    # marker describes: an off-by-one, the wrong build, a minus-strand fetch or a stale
    # cache entry all land here, and any of them would place primers on the wrong sequence.
    if seq[idx] != marker_ref.upper():
        raise ValueError(f"template base at {marker_pos:,} is {seq[idx]!r}, but the marker's "
                         f"GRCh38 ref is {marker_ref.upper()!r}. The template does not match "
                         f"the marker; no primer is designed on it (R6, R7).")

    # R5: room for a primer on BOTH sides, or the pair cannot flank the marker. Enforced
    # here rather than left to primer3, whose failure would read as a hard window.
    need = flank_needed(settings)
    if idx < need or len(seq) - idx - 1 < need:
        raise ValueError(f"template gives {idx} bp before and {len(seq) - idx - 1} bp after the "
                         f"marker; a {settings.max_product} bp max product needs {need} bp each "
                         f"side. Widen the fetch or lower max_product.")

    lo = max(idx - settings.target_pad, 0)
    hi = min(idx + settings.target_pad, len(seq) - 1)
    target = [lo, hi - lo + 1]

    masked, spans = [], []
    for s in snps:
        if s.maf < settings.mask_maf or s.pos == marker_pos:
            continue
        start = max(s.pos - region_start, 0)                       # clamp to the template:
        end = min(s.pos - region_start + s.span, len(seq))         # an indel may overhang
        if start >= end:
            continue                                              # falls outside entirely
        masked.append(s)
        # A site inside the padded target is already off both primers, and the interval
        # array is capped: spending budget on it can cost a pair elsewhere in the window.
        if not (lo <= start and end - 1 <= hi):
            spans.append((start, end - start))

    excluded = _merge(spans)
    n_indel = sum(1 for s in masked if len(s.ref) != 1 or len(s.alt) != 1)
    # Two lengths, for the same reason the warnings have two: the counts differ per marker
    # and belong under every pair, while the sentences after them are identical on every row
    # of the table and are the docs' job to carry.
    note = Note(
        code="mask",
        short=f"Primer sites clear of {len(masked)} variants "
              f"({len(masked) - n_indel} SNPs, {n_indel} indels) at MAF >= "
              f"{settings.mask_maf:.3%}.",
        long=f"Primer-binding regions are clear of {len(masked)} variants "
             f"({len(masked) - n_indel} SNPs, {n_indel} indels) at MAF >= "
             f"{settings.mask_maf:.3%} across this {len(seq):,} bp template. The marker sits "
             f"inside the product with {settings.target_pad} bp of padding, under neither "
             f"primer. Variants the caller did not supply are not masked.",
    )

    if len(excluded) > MAX_MASK_INTERVALS:
        return PrimerResult(
            masked=tuple(masked), mask_note=note,
            error=f"{len(masked)} variants at MAF >= {settings.mask_maf:.3%} merge to "
                  f"{len(excluded)} regions in this {len(seq):,} bp window, over primer3's "
                  f"built-in limit of {MAX_MASK_INTERVALS}. The mask is NOT truncated to fit: "
                  f"a truncated mask designs primers over unmasked common variants while "
                  f"still reporting a mask. Raise mask_maf or narrow the window.")

    globals_ = {
        "PRIMER_OPT_TM": settings.opt_tm,
        "PRIMER_MIN_TM": settings.min_tm,
        "PRIMER_MAX_TM": settings.max_tm,
        "PRIMER_PAIR_MAX_DIFF_TM": settings.max_pair_diff_tm,
        "PRIMER_MIN_SIZE": settings.min_size,
        "PRIMER_OPT_SIZE": settings.opt_size,
        "PRIMER_MAX_SIZE": settings.max_size,
        "PRIMER_MIN_GC": settings.min_gc,
        "PRIMER_MAX_GC": settings.max_gc,
        "PRIMER_GC_CLAMP": settings.gc_clamp,
        "PRIMER_SALT_MONOVALENT": settings.salt_monovalent,
        "PRIMER_SALT_DIVALENT": settings.salt_divalent,
        "PRIMER_DNTP_CONC": settings.dntp_conc,
        "PRIMER_DNA_CONC": settings.dna_conc,
        "PRIMER_PRODUCT_SIZE_RANGE": [[settings.min_product, settings.max_product]],
        "PRIMER_MAX_POLY_X": settings.max_poly_x,
        "PRIMER_MAX_SELF_ANY_TH": settings.max_self_any_th,
        "PRIMER_MAX_SELF_END_TH": settings.max_self_end_th,
        "PRIMER_MAX_HAIRPIN_TH": settings.max_hairpin_th,
        "PRIMER_MAX_NS_ACCEPTED": settings.max_ns_accepted,
        "PRIMER_TM_FORMULA": settings.tm_formula,
        "PRIMER_SALT_CORRECTIONS": settings.salt_corrections,
        "PRIMER_NUM_RETURN": 1,
    }
    seq_args = {"SEQUENCE_ID": f"chr_{marker_pos}", "SEQUENCE_TEMPLATE": seq,
                "SEQUENCE_TARGET": target}
    if excluded:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = excluded

    try:
        out = primer3.design_primers(seq_args, globals_)
    except (OSError, ValueError) as e:
        # primer3's own errors quote the input back, template included. Truncate: this
        # string is browser-bound and an 800 bp dump is not a message.
        raise ValueError(f"primer3 rejected the design: {str(e)[:200]}") from e

    if not out.get("PRIMER_PAIR_NUM_RETURNED"):
        gc = 100.0 * sum(c in "GC" for c in seq) / len(seq)
        return PrimerResult(
            masked=tuple(masked), mask_note=note,
            error=f"No pair found. No {settings.min_size}-{settings.max_size} bp primer reaches "
                  f"{settings.min_tm:.0f}-{settings.max_tm:.0f} C at {settings.min_gc:.0f}-"
                  f"{settings.max_gc:.0f}% GC while clearing {len(masked)} masked variants in "
                  f"this {len(seq):,} bp window ({gc:.0f}% GC). Nothing was relaxed to force a "
                  f"pair: a pair outside the stated Tm or over a common variant is worse than "
                  f"none. Adjust the Tm window, the GC window or mask_maf. "
                  f"primer3: left={out.get('PRIMER_LEFT_EXPLAIN')} "
                  f"right={out.get('PRIMER_RIGHT_EXPLAIN')} "
                  f"pair={out.get('PRIMER_PAIR_EXPLAIN')}")

    l_idx, l_len = out["PRIMER_LEFT_0"]
    r_idx, r_len = out["PRIMER_RIGHT_0"]        # r_idx is the RIGHTMOST base of the oligo
    fwd = Primer(seq=out["PRIMER_LEFT_0_SEQUENCE"], pos=region_start + l_idx, idx=l_idx,
                 length=l_len, tm=round(out["PRIMER_LEFT_0_TM"], 2),
                 gc=round(out["PRIMER_LEFT_0_GC_PERCENT"], 1))
    rev = Primer(seq=out["PRIMER_RIGHT_0_SEQUENCE"], pos=region_start + r_idx - r_len + 1,
                 idx=r_idx - r_len + 1, length=r_len, tm=round(out["PRIMER_RIGHT_0_TM"], 2),
                 gc=round(out["PRIMER_RIGHT_0_GC_PERCENT"], 1))
    return PrimerResult(
        fwd=fwd, rev=rev,
        product_size=out["PRIMER_PAIR_0_PRODUCT_SIZE"], product_start=region_start + l_idx,
        masked=tuple(masked), mask_note=note, warnings=(NOT_CHECKED_WARNING,),
    )


if __name__ == "__main__":
    import random

    if not available():
        print(f"primers self-check SKIPPED: {UNAVAILABLE}")
        raise SystemExit(0)

    # Synthetic, and named so: R1 forbids sequence from a table, and a real genomic string
    # committed here would be exactly that. The logic under test is the mask geometry and
    # the guards, none of which care whether the bases are real.
    rng = random.Random(11)
    # The geometry a caller fetching +/- flank_needed() around the marker actually gets:
    # an inclusive region is 2*flank + 1 bases, and the marker sits on the middle one.
    POS, REGION_START = 17_396_930, 17_396_930 - flank_needed()
    tpl = "".join(rng.choice("ACGT") for _ in range(2 * flank_needed() + 1))
    tpl = tpl[:POS - REGION_START] + "C" + tpl[POS - REGION_START + 1:]

    assert flank_needed() == 400, flank_needed()
    assert len(tpl) == 801 and tpl[POS - REGION_START] == "C"

    # A pair is found, and it carries its warning.
    r = design(tpl, REGION_START, POS, "C")
    assert r.ok, r.error
    assert r.fwd and r.rev, r
    assert NOT_CHECKED_WARNING in r.warnings, r.warnings
    assert r.insilico_pcr == "not_checked", r.insilico_pcr
    assert DEFAULTS.min_product <= r.product_size <= DEFAULTS.max_product, r.product_size
    for p in (r.fwd, r.rev):
        assert DEFAULTS.min_tm <= p.tm <= DEFAULTS.max_tm, p
        assert DEFAULTS.min_gc <= p.gc <= DEFAULTS.max_gc, p
        assert DEFAULTS.min_size <= p.length <= DEFAULTS.max_size, p
        assert len(p.seq) == p.length, p
    assert abs(r.fwd.tm - r.rev.tm) <= DEFAULTS.max_pair_diff_tm, (r.fwd.tm, r.rev.tm)
    # Primer3 measures the product from the left primer's 5' end to the right primer's.
    assert r.product_size == r.rev.end_idx - r.fwd.idx + 1, r
    assert r.fwd.pos == REGION_START + r.fwd.idx and r.product_start == r.fwd.pos

    # The target is inside the product and under neither primer, padded.
    tgt = POS - REGION_START
    assert r.fwd.end_idx < tgt < r.rev.idx, (r.fwd.end_idx, tgt, r.rev.idx)
    assert tgt - r.fwd.end_idx > DEFAULTS.target_pad, "target too close to FWD"
    assert r.rev.idx - tgt > DEFAULTS.target_pad, "target too close to REV"

    # Masked sites are avoided. Seed them right where the unmasked pair sits.
    hits = [Primer(seq="", pos=0, idx=i, length=1, tm=0, gc=0)
            for i in (r.fwd.idx + 2, r.rev.end_idx - 2)]
    snps = [MaskSite(REGION_START + p.idx, "A", "G", 0.20) for p in hits]
    r2 = design(tpl, REGION_START, POS, "C", snps)
    assert r2.ok, r2.error
    assert len(r2.masked) == 2 and r2.masked[0].maf == 0.20, r2.masked
    for p in (r2.fwd, r2.rev):
        for s in r2.masked:
            i = s.pos - REGION_START
            assert not (p.idx <= i <= p.end_idx), f"primer {p} sits on masked site {s}"
    assert (r2.fwd.idx, r2.rev.idx) != (r.fwd.idx, r.rev.idx), "mask changed nothing"
    # The counts belong in BOTH lengths: they are the per-marker fact, and the short form
    # is the only one the table shows.
    for words in (r2.mask_note.short, r2.mask_note.long):
        assert "2 variants (2 SNPs, 0 indels)" in words, words
        assert "MAF >= 1.000%" in words, words

    # An indel masks its whole reference footprint, and is counted as one.
    r3 = design(tpl, REGION_START, POS, "C",
                [MaskSite(REGION_START + r.fwd.idx, "ACGTACGTAC", "A", 0.10)])
    assert r3.ok and "(0 SNPs, 1 indels)" in r3.mask_note.short, r3.mask_note
    for p in (r3.fwd, r3.rev):
        assert not (p.idx <= r.fwd.idx <= p.end_idx), f"{p} sits on the masked indel"

    # An indel overhanging the template start is clamped, not dropped: its tail still sits
    # under a primer. One starting past the end is dropped, having no footprint here.
    over_l = design(tpl, REGION_START, POS, "C",
                    [MaskSite(REGION_START - 4, "A" * (r.fwd.idx + 8), "A", 0.3)])
    assert len(over_l.masked) == 1, over_l.masked
    assert over_l.fwd.idx >= r.fwd.idx + 4, "clamped indel did not push FWD clear"
    assert design(tpl, REGION_START, POS, "C",
                  [MaskSite(REGION_START + len(tpl), "A", "G", 0.3)]).masked == ()

    # Below mask_maf is not masked; the same site above it is.
    rare = [MaskSite(REGION_START + r.fwd.idx + 2, "A", "G", 0.001)]
    assert design(tpl, REGION_START, POS, "C", rare).masked == (), "0.1% masked at a 1% floor"
    assert len(design(tpl, REGION_START, POS, "C", rare,
                      replace(DEFAULTS, mask_maf=0.0005)).masked) == 1, "floor not honoured"

    # The marker itself is never masked: it is the thing being amplified.
    assert design(tpl, REGION_START, POS, "C",
                  [MaskSite(POS, "C", "T", 0.4)]).masked == (), "the target must not be masked"

    # A fully-masked template fails loudly rather than returning a primer over a variant.
    # Spaced under min_size, so no legal primer fits between two masks anywhere.
    gap = DEFAULTS.min_size - 5
    wall = [MaskSite(REGION_START + i, "A", "G", 0.3) for i in range(0, len(tpl), gap)]
    dead = design(tpl, REGION_START, POS, "C", wall)
    assert not dead.ok and dead.fwd is None and dead.rev is None, dead
    assert "No pair found" in dead.error and "worse than none" in dead.error, dead.error
    assert len(dead.masked) == len(wall) and dead.mask_note, "a failure still reports its mask"

    # Past the interval cap it fails on the cap, and says the mask was not truncated.
    flood = [MaskSite(REGION_START + i, "A", "G", 0.3) for i in range(0, 800, 2)]
    over = design(tpl, REGION_START, POS, "C", flood)
    assert not over.ok and "NOT truncated" in over.error, over.error
    assert str(MAX_MASK_INTERVALS) in over.error, over.error
    # ...and merging is what keeps an ordinary window under it.
    assert _merge([(10, 5), (12, 1), (15, 3), (40, 1)]) == [[10, 8], [40, 1]]
    assert _merge([(5, 1), (5, 1)]) == [[5, 1]]

    # Guards. Each of these would otherwise put primers on the wrong sequence.
    def raises(fn, needle):
        try:
            fn()
        except ValueError as e:
            assert needle in str(e), f"wrong error: {e}"
            return
        raise AssertionError(f"expected ValueError containing {needle!r}")

    raises(lambda: design(tpl, REGION_START, POS, "A"), "does not match the marker")
    raises(lambda: design(tpl.replace("A", "U", 1), REGION_START, POS, "C"), "not DNA")
    raises(lambda: design(tpl, REGION_START, REGION_START + 900, "C"), "outside the template")
    raises(lambda: design(tpl[:600], REGION_START, POS, "C"), "each side")
    # R5 again, from the other end: a marker 100 bp into a long template has no room behind.
    raises(lambda: design(tpl + tpl, REGION_START, REGION_START + 100, tpl[100]), "each side")
    # A minus-strand template is caught by the same ref check, which is the point of it.
    comp = str.maketrans("ACGT", "TGCA")
    raises(lambda: design(tpl.translate(comp)[::-1], REGION_START, POS, "C"),
           "does not match the marker")

    # Settings are validated, not trusted: an out-of-model Tm is refused at the door.
    raises(lambda: replace(DEFAULTS, max_size=36), "Tm model covers")
    raises(lambda: replace(DEFAULTS, min_tm=72.0), "Tm window out of order")
    raises(lambda: replace(DEFAULTS, min_size=30, opt_size=26), "size window out of order")
    raises(lambda: replace(DEFAULTS, min_product=250, max_product=100), "product window")
    raises(lambda: replace(DEFAULTS, min_product=100, max_product=140), "padded target")
    raises(lambda: replace(DEFAULTS, mask_maf=0.7), "mask_maf out of range")
    # The coupling is enforced, not advisory: a wider product needs a wider template.
    wide = replace(DEFAULTS, max_product=900)
    assert flank_needed(wide) == 550, flank_needed(wide)
    raises(lambda: design(tpl, REGION_START, POS, "C", settings=wide), "Widen the fetch")

    # Overrides reach primer3 rather than being decoration.
    lo_tm = replace(DEFAULTS, opt_tm=60.0, min_tm=58.0, max_tm=62.0)
    r4 = design(tpl, REGION_START, POS, "C", settings=lo_tm)
    assert r4.ok and 58.0 <= r4.fwd.tm <= 62.0, r4

    # The honesty invariants the three consuming lanes rest on.
    assert not any(x.ok and not x.warnings for x in (r, r2, r3, r4)), "a pair without warnings"
    assert all(x.fwd is None and x.rev is None for x in (dead, over)), "a primer with an error"
    assert "NOT CHECKED" in NOT_CHECKED_WARNING.long and "GRCh38" in NOT_CHECKED_WARNING.long

    # Both lengths, on every note, and actually two lengths. A `short` copied from `long` is
    # the failure this shape exists to prevent: it renders down a page nobody then reads.
    for n in (NOT_CHECKED_WARNING, UNAVAILABLE, r.mask_note, *r.warnings):
        assert n.code and n.docs and n.short and n.long, n
        assert len(n.short) <= 120, (n.code, len(n.short))
        assert len(n.long) > len(n.short), n.code
        assert str(n) == n.long, "print must be the long form"

    print(f"primers self-check OK  |  primer3 {primer3.__version__}  "
          f"FWD {r.fwd.length}bp Tm {r.fwd.tm} GC {r.fwd.gc}%  "
          f"REV {r.rev.length}bp Tm {r.rev.tm} GC {r.rev.gc}%  "
          f"product {r.product_size}bp  |  masked pair shifted to "
          f"{r2.fwd.idx}/{r2.rev.idx} from {r.fwd.idx}/{r.rev.idx}")

"""origin - the end-to-end answer for one embryo: which paternal homologue came through.

This is the integrator. Every layer below it is tested in isolation; this is the path from
parental and embryo genotypes to a single verdict, with the refusals attached.

The chain, and what each link can and cannot settle:

    normalize   the intensity channel onto this sample's own scale, via a male chrX reference
    score       per marker, is the paternal allele represented, absent, or unable to say
    run         do the absences RUN, which separates a structural event from independent dropout
    states      what copy configuration explains the segment, which is the ONLY thing that
                separates a gene-conversion tract from a paternal deletion
    verdict     H2-vs-H3 always; H1-vs-H2 only with a declared phase source, never otherwise

The last line is the honest boundary of the whole tool. "Was the disease allele inherited and
corrected, or was the wild-type inherited" is H1 versus H2, and those two differ ONLY in which of
the father's homologues was transmitted. Nothing in any genotype file says which of his two
chromosomes carried the mutation, so absent an external phase source the answer is refused rather
than guessed - and refused loudly, because a confident wrong answer here is exactly the error
that produced the Ma/Egli dispute.

The run-length statistic is reimplemented here rather than imported from the TypeScript layer.
That duplication is deliberate and bounded: the browser needs phases 1 and 2 without loading a
20 MB runtime, and the CLI needs them without a JavaScript one. Both sides are pinned by the same
two external artefacts - `informativity_table.csv` and the 12-case golden fixture - so a
divergence fails a test rather than reaching a result.

Self-check:  python -m origin          (no arguments; the CLI needs --father/--embryo/...)
Full tests:  python -m pytest tests/test_origin.py
"""

from __future__ import annotations

import gzip
import math
import statistics as st
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional, Sequence

import buildref
import emissions as em
import hmm
import normalize as nz
import structural

#: Kothiyal 2019's per-trio Mendelian-inconsistency floor, 0.63% of variant sites across 1,314
#: nuclear families. A LOWER BOUND on the spurious-violation rate, never the value itself.
KOTHIYAL_FLOOR = 0.0063


# --- the run-length statistic (mirrors runlength.ts, pinned by the same fixture) -----------


def run_length_p(r: int, n: int, q: float) -> float:
    """P(longest run >= r) among n independent Bernoulli(q) trials, Erdos-Renyi.

    Written as -expm1(m * log1p(-x)) because the direct form `1 - (1-q^r)^m` cancels
    catastrophically: at q = 0.0063 and r = 10 the truth is 9.8e-23 but `1 - x` rounds to exactly
    1.0 and the subtraction returns 0. Reporting p = 0 for a real event would be unfalsifiable.
    """
    if r <= 0:
        return 1.0
    if r > n or q <= 0.0:
        return 0.0
    if q >= 1.0:
        return 1.0
    x = q ** r
    if x >= 1.0:
        return 1.0
    return -math.expm1((n - r + 1) * math.log1p(-x))


def r_min(n: int, q: float, alpha: float = 0.05, windows: int = 3) -> Optional[int]:
    """Smallest run length clearing the window-wise error budget. Never hardcoded."""
    budget = alpha / max(1, windows)
    for r in range(1, n + 1):
        if run_length_p(r, n, q) <= budget:
            # Never below 2: contiguity is the whole discriminator and a run of one has none.
            # Below n = 3 the threshold comes back as 1, so a single absent marker would read
            # as a significant run on the strength of being the only marker in the window.
            return max(2, r)
    return None


def q_hat(empirical: Optional[float]) -> tuple[float, str]:
    """max(empirical, Kothiyal floor). The max, not the sum: Kothiyal already contains error
    classes the empirical estimator also captures, so adding them double-counts."""
    if empirical is None or not math.isfinite(empirical) or empirical <= KOTHIYAL_FLOOR:
        return KOTHIYAL_FLOOR, "kothiyal_floor"
    return empirical, "empirical"


def score_paternal(father: str, mother: Optional[str], embryo: str) -> Optional[int]:
    """1 if the paternal allele is provably represented, 0 if provably absent, None if silent.

    Needs a HOMOZYGOUS father either way: a heterozygous one is compatible with either observed
    allele and so testifies to nothing.

    The two directions then need different things from the mother, which is the whole reason she
    appears here. ABSENCE is Mendelian: a homozygous father must transmit his allele, so an embryo
    lacking it is missing the paternal contribution whatever the mother's genotype. That is why
    the no-mother mode keeps this half. PRESENCE is an identity claim: the embryo carrying his
    allele only shows the allele came from HIM if the mother could not have supplied it, so it
    needs her homozygous for the other allele.

    Scoring presence without that check is not a conservative error. At a paternal deletion the
    embryo is hemizygous for the maternal allele, and wherever a heterozygous mother transmitted
    the same allele the father carries, the call still contains his allele and reads as "present".
    That breaks the run at roughly half of all mother-heterozygous markers, which is how a
    whole-chromosome absence came back as a run of 21.
    """
    if father not in ("AA", "BB") or embryo in ("NC", ""):
        return None
    allele = father[0]
    if allele not in embryo:
        return 0
    other = "B" if allele == "A" else "A"
    return 1 if mother == other * 2 else None


#: A run is split where consecutive absent markers lie further apart than this multiple of the
#: local informative spacing. Generous on purpose: inside a real deletion the absent markers sit
#: at ordinary array spacing, so a true event never approaches this.
GAP_FACTOR = 10.0


@dataclass
class AbsenceRun:
    """A stretch of consecutive paternal-absence markers, and the ground it actually covers."""

    start: int
    length: int
    lo: int
    hi: int


def absence_runs(
    scored: Sequence[tuple["Marker", Optional[int]]], max_gap_bp: Optional[float] = None
) -> list[AbsenceRun]:
    """Maximal runs of paternal absence, split where the chromosome was never interrogated.

    Adjacency in the informative subsequence is not adjacency on the chromosome. On a real trio
    three violations spread across 5.6 Mb of chromosome 7 counted as a run of 3, exactly like
    three violations inside 1 kb, because nothing informative happened to lie between them. A
    5.6 Mb hemizygous deletion in a healthy adult is a major clinical finding, not a quiet data
    point, and claiming one on three observations asserts that the intervening megabases are
    absent on no evidence at all.
    """
    runs: list[AbsenceRun] = []
    start: Optional[int] = None
    prev: Optional[int] = None
    for i, (m, s) in enumerate(scored):
        if s != 0:
            if start is not None:
                runs.append(AbsenceRun(start, i - start, scored[start][0].pos, prev))
            start = prev = None
            continue
        if start is not None and max_gap_bp is not None and m.pos - prev > max_gap_bp:
            runs.append(AbsenceRun(start, i - start, scored[start][0].pos, prev))
            start = i
        elif start is None:
            start = i
        prev = m.pos
    if start is not None:
        runs.append(AbsenceRun(start, len(scored) - start, scored[start][0].pos, prev))
    return runs


BidirectionalVerdict = Literal["one_directional", "artefact", "insufficient", "unavailable"]


def _log_binom_pmf(k: int, n: int, p: float) -> float:
    if not 0.0 < p < 1.0:
        return 0.0 if (p <= 0.0 and k == 0) or (p >= 1.0 and k == n) else -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def bidirectional_test(
    informative: int, violations: int, q: float,
    floor: float = KOTHIYAL_FLOOR, decisive: float = 1.0,
) -> tuple[BidirectionalVerdict, float]:
    """Real paternal loss, or amplification dropout wearing its shape?

    Across a real paternal deletion the embryo is hemizygous for the MATERNAL allele, so nothing
    in that span can show the maternal allele missing: the maternal copy is physically present
    and called correctly. Dropout has no such asymmetry. It removes alleles without regard to
    which parent supplied them, so it scatters violations in BOTH directions through the same
    region. The maternal-violation count inside the span therefore separates the two.

    Under dropout that count is Binomial(k, q); under a real loss it is Binomial(k, floor), the
    residual genotyping-error rate. Returns the verdict and log10 of the likelihood ratio,
    positive favouring artefact.

    The test loses power exactly as it stops being needed. When q falls to the error floor the
    two hypotheses coincide, the ratio goes to one, and the answer is "insufficient" - which is
    correct, not a shortcoming, because a sample with no excess dropout cannot be faked by
    dropout. On a blastomere at q = 0.33 it is decisive; on bulk DNA at q = 0.008 it says
    nothing, and needs to say nothing.
    """
    if informative <= 0:
        return "insufficient", 0.0
    lr = (_log_binom_pmf(violations, informative, q)
          - _log_binom_pmf(violations, informative, floor)) / math.log(10.0)
    if not math.isfinite(lr):
        return "insufficient", 0.0
    if lr >= decisive:
        return "artefact", lr
    if lr <= -decisive:
        return "one_directional", lr
    return "insufficient", lr


# --- inputs -------------------------------------------------------------------------------


@dataclass
class Marker:
    rsid: str
    chrom: str
    pos: int
    father: str
    mother: Optional[str]
    embryo: str
    baf: Optional[float] = None
    lrr: Optional[float] = None


@dataclass
class PhaseSource:
    """Which of the father's homologues carries the pathogenic allele, and how that was learnt.

    Required for any H1-vs-H2 claim. There is no default and no inference: both of the father's
    homologues are simply "his alleles" until something external distinguishes them.
    """

    mutant_allele: Literal["A", "B"]
    route: str


# --- the verdict ---------------------------------------------------------------------------


@dataclass
class WindowVerdict:
    window: str
    n_l3: int
    z_sum: int
    longest_run: int
    run_p: float
    r_min: Optional[int]
    q: float
    q_source: str
    #: Whether the run survived every gate: length, the maternal cross-check, the no-mother
    #: refusal. Read this rather than recomputing from longest_run and r_min, which bypasses
    #: the refusals.
    significant: bool = False
    #: Physical extent of the longest run, and the gap cap that was allowed to break it.
    run_lo: Optional[int] = None
    run_hi: Optional[int] = None
    max_gap_bp: Optional[float] = None
    #: The bidirectional test over that span: does the maternal side also show violations?
    bidirectional: BidirectionalVerdict = "unavailable"
    bidirectional_log10_lr: float = 0.0
    maternal_informative: int = 0
    maternal_violations: int = 0
    #: The state the HMM settled on across the affected segment, if any.
    segment_state: Optional[str] = None
    segment_markers: int = 0
    #: Separation between copy-loss and copy-neutral over THAT segment length, not per marker.
    separation_sigma: float = math.nan
    hypotheses: dict[str, float] = field(default_factory=dict)
    ledger: list[hmm.LedgerEntry] = field(default_factory=list)
    refusals: list[str] = field(default_factory=list)


@dataclass
class EmbryoVerdict:
    embryo_id: str
    calibration: nz.Calibration
    dropout: Optional[float]
    windows: list[WindowVerdict] = field(default_factory=list)
    refusals: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """One paragraph, in the language section 14 requires: bounded exclusion, not a verdict."""
        lines = [f"{self.embryo_id}"]
        lines.append(
            f"  calibration: c={self.calibration.lrr_compression:.3f} "
            f"sigma={self.calibration.sigma_lrr:.3f} via {self.calibration.route}"
        )
        if self.dropout is not None:
            lines.append(f"  fitted dropout: {self.dropout:.3f}")
        for w in self.windows:
            if w.hypotheses:
                top = max(w.hypotheses.items(), key=lambda kv: kv[1])
                tail = f"-> {top[0]} {top[1]:.2f}"
            else:
                # Never print a "none" that could be read as "nothing found" beside a run this
                # extreme. No hypotheses means the state layer did not run, not that it ran empty.
                tail = "(no hypothesis: state layer did not run)"
            lines.append(
                f"  [{w.window}] {w.n_l3} L3 markers, longest paternal-absence run "
                f"{w.longest_run} (r_min {w.r_min}, p={w.run_p:.2e}) {tail}"
            )
            if w.longest_run and w.run_lo is not None:
                kb = (w.run_hi - w.run_lo) / 1000.0
                lines.append(f"      run spans {w.run_lo:,}-{w.run_hi:,} ({kb:,.1f} kb)")
            if w.bidirectional != "unavailable" and w.maternal_informative:
                verdict = {
                    "artefact": "DROPOUT, not a paternal loss",
                    "one_directional": "one-directional, consistent with real paternal loss",
                    "insufficient": "cannot say: too little maternal evidence in this span",
                }[w.bidirectional]
                lines.append(
                    f"      maternal check: {w.maternal_violations}/{w.maternal_informative} "
                    f"markers also missing the MATERNAL allele -> {verdict} "
                    f"(10^{w.bidirectional_log10_lr:+.1f})"
                )
            if w.segment_state:
                lines.append(
                    f"      segment: {w.segment_state} over {w.segment_markers} markers, "
                    f"copy-loss vs copy-neutral separated by {w.separation_sigma:.1f} sigma"
                )
            for e in w.ledger:
                if e.verdict != "excluded":
                    lines.append(f"      {e.mechanism}: {e.verdict.upper()}")
        for r in self.refusals:
            lines.append(f"  REFUSED: {r}")
        return "\n".join(lines)


def analyse_embryo(
    embryo_id: str,
    markers: Sequence[Marker],
    variant_chrom: str,
    variant_pos: int,
    *,
    calibration: nz.Calibration,
    dropout: Optional[float] = None,
    phase: Optional[PhaseSource] = None,
    sigma_baf: float = 0.05,
    eps: float = 0.02,
    structural_evidence: Optional["structural.StructuralEvidence"] = None,
) -> EmbryoVerdict:
    """Run the whole chain for one embryo.

    Markers may arrive unsorted; they are ordered here, because runs are physical and an unsorted
    input would silently fragment them and understate every event.
    """
    v = EmbryoVerdict(embryo_id=embryo_id, calibration=calibration, dropout=dropout)

    if phase is None:
        v.refusals.append(
            "No phase source declared. H1 (paternal wild-type inherited) and H2 (paternal mutant "
            "inherited and repaired) differ ONLY in which paternal homologue was transmitted, and "
            "nothing in these files says which of the father's chromosomes carried the mutation. "
            "That half of the question is REFUSED. H2-vs-H3 - whether the paternal allele is there "
            "at all, and whether its absence is copy-reducing - is reported below."
        )
    if structural_evidence is not None:
        v.refusals.append("structural evidence: " + structural_evidence.summary())
        for n in structural_evidence.notes:
            v.refusals.append(n)
    if not math.isfinite(calibration.lrr_compression):
        # Quote the calibrator's own reason rather than restating a guess at it. A male WGA'd
        # embryo fails the noise ceiling with a perfectly good chrX reference available, and
        # saying "no male chrX reference" there would send the operator after the wrong fix.
        why = calibration.warnings[-1] if calibration.warnings else "no reference available"
        v.refusals.append(
            f"Intensity channel is centred but UNSCALED ({calibration.route}). {why} "
            "Copy-number states cannot be distinguished, so copy-neutral LOH is not tested."
        )

    q, q_src = q_hat(dropout)
    params = em.EmissionParams(
        lrr_compression=calibration.lrr_compression
        if math.isfinite(calibration.lrr_compression) else 0.0,
        sigma_baf=sigma_baf,
        sigma_lrr=calibration.sigma_lrr if math.isfinite(calibration.sigma_lrr) else 1.0,
        eps=eps,
    )

    want = variant_chrom.replace("chr", "").lower()
    # chrX probes carry a per-array intensity offset unrelated to copy number, and without it a
    # male's single X reads as two copies. Where it has not been measured the copy-number layer
    # is refused on chrX. The genotype statistics below are untouched.
    on_sex_chrom = want in ("x", "23")
    offset_missing = on_sex_chrom and calibration.chrx_offset is None
    if offset_missing:
        v.refusals.append(
            "chrX copy-number states NOT TESTED: no measured chrX probe offset for this array. "
            "The offset is most of the chrX-to-autosome shift, so without it a single X reads as "
            "two copies. Paternal presence and absence below are genotype-only and unaffected."
        )
    on_chrom = sorted(
        (m for m in markers if m.chrom.replace("chr", "").lower() == want and m.pos != variant_pos),
        key=lambda m: m.pos,
    )
    # Without the mother presence is unobservable, so the informative set holds absences only,
    # nothing breaks a run, and the statistic has no null: a normal chromosome 20 gave a run of 3
    # across 35 Mb at p = 2.5e-07. Dropout is unfittable too, so q collapses to the error floor.
    has_mother = any(m.mother is not None for m in on_chrom)
    if not has_mother:
        v.refusals.append(
            "No maternal genotypes, so THIS per-locus test is not run. Paternal presence cannot "
            "be established at any single marker without someone who could not have supplied the "
            "allele, the informative set holds absences only, and the run-length statistic has "
            "no null to be significant against. Dropout cannot be fitted either. Violation counts "
            "are still reported as counts.\n"
            "      This limit is specific to the per-locus deletion test. Parent of origin, "
            "sperm type and segmental loss are all decided from RATES across hundreds of "
            "thousands of markers and need no mother at all. Adding the oocyte donor is a strong "
            "recommendation rather than a requirement: it measures dropout directly, restores "
            "this test, and raises confidence throughout."
        )

    for name, half in (("local", 25_000), ("segmental", 10_000_000),
                       ("whole_chromosome", math.inf)):
        in_win = [m for m in on_chrom if abs(m.pos - variant_pos) <= half]
        # A marker that cannot say must be EXCLUDED, not counted as evidence of presence.
        # Folding None into the flag list breaks runs and inflates n at the same time.
        l3 = [(m, s) for m in in_win
              if (s := score_paternal(m.father, m.mother, m.embryo)) is not None]
        flags = [s == 0 for _, s in l3]
        n = len(flags)
        z = sum(flags)
        spacing = _median_spacing([m.pos for m, _ in l3])
        max_gap = GAP_FACTOR * spacing if spacing else None
        runs = absence_runs(l3, max_gap)
        top = max(runs, key=lambda r: r.length, default=None)
        best = top.length if top else 0
        thr = r_min(n, q) if n else None

        wv = WindowVerdict(window=name, n_l3=n, z_sum=z, longest_run=best,
                           r_min=thr, q=q if has_mother else math.nan,
                           q_source=q_src if has_mother else "uncalibrated_no_mother",
                           run_p=run_length_p(best, n, q) if (n and has_mother) else math.nan,
                           run_lo=top.lo if top else None, run_hi=top.hi if top else None,
                           max_gap_bp=max_gap)
        significant = thr is not None and best >= thr and has_mother

        # Does the maternal side violate over the same ground? Only dropout can do that, so this
        # is what separates a real paternal loss from amplification noise shaped like one.
        if significant and top is not None:
            if any(m.mother is not None for m in in_win):
                mat = [score_paternal(m.mother, m.father, m.embryo)
                       for m in in_win if top.lo <= m.pos <= top.hi]
                wv.maternal_informative = sum(1 for s in mat if s is not None)
                wv.maternal_violations = sum(1 for s in mat if s == 0)
                wv.bidirectional, wv.bidirectional_log10_lr = bidirectional_test(
                    wv.maternal_informative, wv.maternal_violations, q)
            else:
                wv.bidirectional = "unavailable"

        # Only a positive refutation blocks. "Insufficient" leaves the run standing with the
        # caveat attached, because on clean DNA the test is powerless by construction and
        # withholding every clean-DNA result would discard the cases that work best.
        refuted = wv.bidirectional == "artefact"
        if refuted:
            wv.refusals.append(
                f"Paternal-absence run REFUTED as amplification dropout: {wv.maternal_violations} "
                f"of {wv.maternal_informative} maternal-informative markers inside "
                f"{wv.run_lo:,}-{wv.run_hi:,} also show the MATERNAL allele missing "
                f"(likelihood ratio 10^{wv.bidirectional_log10_lr:.1f} favouring dropout over a "
                "real loss). A deletion removes one parent's copy, never both, so this span is "
                "not evidence that the paternal allele is absent."
            )
        significant = significant and not refuted
        wv.significant = significant

        # The signal layer runs only on the window that actually shows something: a state call
        # over a window with no candidate event is noise fitting.
        if significant and math.isfinite(calibration.lrr_compression) and not offset_missing:
            shift = calibration.center + (calibration.chrx_offset or 0.0 if on_sex_chrom else 0.0)
            obs = [
                hmm.MarkerObs(m.rsid, m.pos, m.father, m.mother, baf=m.baf,
                              lrr=None if m.lrr is None else m.lrr - shift,
                              genotype=m.embryo)
                for m in in_win
            ]
            tp = hmm.TransitionParams(event_length_bp=max(20_000.0, half / 50)
                                      if math.isfinite(half) else 10_000_000.0)
            inf = hmm.infer(obs, variant_chrom, params, tp)
            path = [hmm.STATES[j].name for j in inf.viterbi]
            at_variant = min(range(len(obs)), key=lambda i: abs(obs[i].pos - variant_pos)) \
                if obs else None
            # The stretch containing the variant, not the longest anywhere in the window, and
            # contiguous rather than scattered: counting scattered markers read a normal
            # chromosome 3 as 7,499 markers at 359 sigma.
            if at_variant is not None and path[at_variant] != "PAT1_MAT1":
                state = path[at_variant]
                lo = hi = at_variant
                while lo > 0 and path[lo - 1] == state:
                    lo -= 1
                while hi + 1 < len(path) and path[hi + 1] == state:
                    hi += 1
                wv.segment_state = state
                wv.segment_markers = hi - lo + 1
                wv.separation_sigma = calibration.separation_sigma(wv.segment_markers)
            if at_variant is not None:
                hp = hmm.hypotheses_at(
                    inf.gamma[at_variant],
                    # z == n, not the longest run: compactness splits runs at marker deserts, so
                    # a genuinely absent chromosome fragments into pieces.
                    whole_chromosome_absent=(name == "whole_chromosome" and z == n and n > 0),
                    informative_markers=n,
                )
                wv.hypotheses = hp.posterior
                wv.refusals += hp.refusals

        wv.ledger = hmm.mechanism_ledger(
            run_significant=significant,
            copy_neutral_separation_sigma=wv.separation_sigma
            if math.isfinite(wv.separation_sigma) else 0.0,
            marker_spacing_bp=_median_spacing([m.pos for m, _ in l3]),
            r_min=thr,
            has_lrr=math.isfinite(calibration.lrr_compression) and not offset_missing
            and any(m.lrr is not None for m in in_win),
            # An array cannot see a large on-target deletion or an insertion at all. Reads can,
            # and structural calls are how they get here.
            has_reads=bool(structural_evidence and structural_evidence.supplied),
            construct_sequence_supplied=bool(
                structural_evidence and structural_evidence.spans_locus),
        )
        v.windows.append(wv)
    return v


def _median_spacing(positions: Sequence[int]) -> Optional[float]:
    if len(positions) < 2:
        return None
    gaps = [b - a for a, b in zip(positions, positions[1:]) if b > a]
    return st.median(gaps) if gaps else None


def estimate_dropout(
    father: dict[str, str], mother: dict[str, str], embryo: dict[str, str]
) -> Optional[tuple[float, float, int]]:
    """d = 1 - 2h at father-het x mother-homozygous markers. Returns (d, se, n).

    Frequency-free: Mendel gives exactly half AA and half AB there regardless of allele
    frequency, so no reference panel and no assumed population heterozygosity are needed.

    Only valid when `embryo` really is the offspring of the other two. Handing it a parent
    returns a clamped zero, which is not a measurement - a mistake worth naming because I made it.
    """
    n = het = 0
    for rsid, fa in father.items():
        if fa != "AB":
            continue
        mo = mother.get(rsid)
        if mo not in ("AA", "BB"):
            continue
        e = embryo.get(rsid)
        if e is None or e == "NC":
            continue
        n += 1
        if e == "AB":
            het += 1
    if n == 0:
        return None
    h = het / n
    return min(1.0, max(0.0, 1.0 - 2.0 * h)), 2.0 * math.sqrt(h * (1 - h) / n), n


# --- reading files --------------------------------------------------------------------------
#
# Only the CLI needs this; the browser parses its own drop box. Pinned against the real
# published files, the same way the run statistic is.

#: Genotype coding. Axiom writes 0/1/2, Illumina writes the AB-space letters directly. Anything
#: else is a no-call, never guessed: a nucleotide-space report degrades to nearly all no-calls
#: and trips the quality gate rather than being silently misread as AB space.
_GT = {"0": "AA", "1": "AB", "2": "BB",
       "AA": "AA", "AB": "AB", "BA": "AB", "BB": "BB"}

#: Column aliases in normalised form, mirroring `ingest.ts`.
_ALIASES: dict[str, tuple[str, ...]] = {
    "rsid": ("probesetid", "snpname", "name", "rsid", "id", "markername", "rs#", "snp"),
    "chrom": ("chr", "chrom", "chromosome", "chrid"),
    "pos": ("position", "pos", "physicalposition", "mapinfo", "coordinate"),
    "gt": ("genotype", "call", "gtype", "geno"),
    "baf": ("baf", "ballelefreq", "ballelefrequency"),
    "lrr": ("log2r", "normalizedintensity", "lrr", "logrratio"),
    "cn": ("copynumber", "cn"),
    # GenomeStudio splits the call across two columns; the strand suffix is stripped by _norm.
    "a1": ("allele1", "allele1forward", "allele1top", "allele1ab", "allele1plus"),
    "a2": ("allele2", "allele2forward", "allele2top", "allele2ab", "allele2plus"),
    # Long-format exports carry every sample in one file, one row per sample and SNP.
    "sample": ("sampleid", "sample", "sampleName", "samplename", "iid"),
    #: HapMap-style files name both alleles for the marker, which makes AB mapping exact.
    "alleles": ("alleles",),
}

#: Delimiters tried in order of specificity. Space last: a space-delimited file usually has no
#: tabs or commas, but a tab-delimited one often has spaces inside field values.
_DELIMS = ("\t", ",", None)


def _norm(col: str) -> str:
    """Normalise a column name for alias matching.

    GenomeStudio prefixes every per-sample column with the sample id, so `Log R Ratio` arrives as
    `99HI0698C.Log R Ratio`. Dropping everything before the last dot and squashing separators
    makes the Illumina and Axiom headers resolve through one table.
    """
    out = col.strip().rsplit(".", 1)[-1].lower()
    for ch in (" ", "_", "-"):
        out = out.replace(ch, "")
    return out


@dataclass
class Probe:
    chrom: str
    pos: int
    gt: str
    baf: Optional[float] = None
    lrr: Optional[float] = None
    cn: Optional[float] = None


@dataclass
class Sample:
    sample_id: str
    probes: dict[str, Probe]
    columns: list[str]
    notes: list[str] = field(default_factory=list)

    def autosomal(self, attr: str) -> list[float]:
        return [v for p in self.probes.values()
                if _is_autosome(p.chrom) and (v := getattr(p, attr)) is not None]

    def chrx(self, attr: str) -> list[float]:
        return [v for p in self.probes.values()
                if p.chrom in ("X", "23") and (v := getattr(p, attr)) is not None]


def _is_autosome(chrom: str) -> bool:
    return chrom.isdigit() and 1 <= int(chrom) <= 22


def _resolve(cols: list[str]) -> dict[str, int]:
    idx: dict[str, int] = {}
    for field_name, names in _ALIASES.items():
        for n in names:
            if n in cols:
                idx[field_name] = cols.index(n)
                break
    return idx


def _header_line(lines: Sequence[str]) -> tuple[int, str, list[str], dict[str, int]]:
    """Find the row that is actually the column header, and how the file is delimited.

    GenomeStudio writes a `[Header]` block of metadata before `[Data]`, so the columns are on
    line ten rather than line one. Rather than special-casing that marker, take the first line
    that RESOLVES to the fields needed. It handles the marker, files with a comment preamble, and
    files with neither, without knowing which it is looking at.
    """
    best: tuple[int, str, list[str], dict[str, int]] | None = None
    for i, raw in enumerate(lines):
        if not raw.strip():
            continue
        for d in _DELIMS:
            cols = [_norm(c) for c in raw.rstrip("\n").split(d)]
            if len(cols) < 3:
                continue
            idx = _resolve(cols)
            if {"rsid", "chrom", "pos"} <= set(idx):
                # A genotype column is preferred but not required here: a wide export can name
                # its sample columns by nothing but the sample, and that is only visible once the
                # layout is known. `read_samples` refuses if no route to genotypes exists.
                if "gt" in idx or ("a1" in idx and "a2" in idx):
                    return i, d, cols, idx
                if best is None:
                    best = (i, d, cols, idx)
            elif best is None and len(idx) >= 3:
                best = (i, d, cols, idx)
    if best is not None:
        return best
    raise ValueError("no row resolves to a recognisable column header")


def read_plink(prefix: str | Path) -> dict[str, Sample]:
    """PLINK, binary or text. Give the prefix, or any one of the files.

    `.bim` names A1 and A2 for every marker, so the AB mapping is exact here rather than pooled
    across samples: A1 becomes A and A2 becomes B, per marker, from the file itself.

    The binary layout is the documented one: two magic bytes, a mode byte that must be 1 for
    SNP-major, then ceil(n/4) bytes per marker holding two bits per sample, least significant
    pair first. The codes are 00 homozygous A1, 01 MISSING, 10 heterozygous, 11 homozygous A2.
    Note that 01 is the missing code and not a genotype: reading it as one silently invents a
    call for every sample that failed.
    """
    prefix = Path(str(prefix))
    for ext in (".bed", ".bim", ".fam", ".ped", ".map"):
        if prefix.name.endswith(ext):
            prefix = prefix.with_suffix("")
            break
    bim, fam, bed = (prefix.with_suffix(e) for e in (".bim", ".fam", ".bed"))
    ped, mp = prefix.with_suffix(".ped"), prefix.with_suffix(".map")

    if bim.exists() and fam.exists() and bed.exists():
        markers = []
        for line in bim.read_text().splitlines():
            f = line.split()
            if len(f) >= 6:
                markers.append((f[1], f[0].removeprefix("chr"), int(f[3]), f[4], f[5]))
        ids = [ln.split()[1] for ln in fam.read_text().splitlines() if ln.split()]
        raw = bed.read_bytes()
        if raw[:2] != b"\x6c\x1b":
            raise ValueError(f"{bed.name}: not a PLINK .bed (magic bytes wrong)")
        if raw[2] != 1:
            raise ValueError(f"{bed.name}: individual-major .bed is not supported; "
                             "re-export with a current PLINK, which writes SNP-major")
        stride = (len(ids) + 3) // 4
        out: dict[str, dict[str, Probe]] = {i: {} for i in ids}
        for m, (rsid, chrom, pos, a1, a2) in enumerate(markers):
            block = raw[3 + m * stride: 3 + (m + 1) * stride]
            if len(block) < stride:
                break
            for j, sid in enumerate(ids):
                code = (block[j >> 2] >> ((j & 3) * 2)) & 3
                gt = {0: "AA", 1: "NC", 2: "AB", 3: "BB"}[code]
                out[sid][rsid] = Probe(chrom, pos, gt)
        return {k: Sample(k, v, ["plink-bed"]) for k, v in out.items()}

    if ped.exists() and mp.exists():
        markers = []
        for line in mp.read_text().splitlines():
            f = line.split()
            if len(f) >= 4:
                markers.append((f[1], f[0].removeprefix("chr"), int(f[3])))
        out = {}
        for line in ped.read_text().splitlines():
            f = line.split()
            if len(f) < 6:
                continue
            sid, calls = f[1], f[6:]
            probes = {}
            for m, (rsid, chrom, pos) in enumerate(markers):
                a, b = calls[2 * m: 2 * m + 2] or ("0", "0")
                g = "NC" if a in ("0", "N") or b in ("0", "N") else "".join(sorted(a + b))
                probes[rsid] = Probe(chrom, pos, g)
            out[sid] = Sample(sid, probes, ["plink-ped"])
        # .ped carries nucleotides with no A1/A2 column, so the mapping is pooled as elsewhere.
        return {s.sample_id: s for s in harmonise(list(out.values()))}

    raise ValueError(f"{prefix.name}: no .bed/.bim/.fam or .ped/.map alongside this prefix")


def read_samples(path: str | Path) -> dict[str, Sample]:
    """Every sample in a file, keyed by id. Single-sample files return one entry.

    Handles the format families a genotyping lab actually produces: wide exports with one column
    per sample, long exports with a `Sample ID` column and one row per sample and SNP, calls in a
    single column or split across two allele columns, and genotypes coded numerically, in AB
    space or as nucleotides.

    Nucleotide calls are NOT converted here. Mapping them to AB needs the two alleles seen at a
    marker, and a single homozygous sample shows only one of them, so the mapping has to be
    agreed across every sample being compared. `harmonise` does that; doing it per file would let
    two files disagree about which allele is A.
    """
    path = Path(path)
    if any(path.name.endswith(x) for x in (".vcf", ".vcf.gz")):
        s = read_vcf(path)
        return {s.sample_id: s}
    if any(path.name.endswith(x) for x in (".bed", ".bim", ".fam", ".ped", ".map")):
        return read_plink(path)
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        lines = fh.read().splitlines()
    if not any(l.strip() for l in lines):
        raise ValueError(f"{path.name}: empty file")
    try:
        hdr_i, delim, cols, idx = _header_line(lines[:200])
    except ValueError as e:
        raise ValueError(f"{path.name}: {e}. Found: {', '.join(lines[0].split()[:12])}") from None
    missing = [r for r in ("rsid", "chrom", "pos") if r not in idx]
    if missing:
        raise ValueError(f"{path.name}: no column for {', '.join(missing)}. "
                         f"Found: {', '.join(cols[:20])}")

    raw_cols = lines[hdr_i].rstrip("\n").split(delim)
    width = len(cols)
    long_mode = "sample" in idx
    # A wide export names each per-sample column `<sample>.GType`, so the sample ids are the
    # prefixes; a long export names them in a column.
    wide_ids: list[str] = []
    wide_gt_col: dict[str, int] = {}
    # Per-sample column map for a wide export. Without this only the first sample was read and
    # the rest were silently dropped, which a two-sample PennCNV-style file makes immediate.
    wide_cols: dict[str, dict[str, int]] = {}
    if not long_mode:
        for j, c in enumerate(raw_cols):
            if "." not in c:
                continue
            field = _norm(c)
            if field not in ("gtype", "genotype", "logrratio", "ballelefreq", "copynumber"):
                continue
            sid = c.rsplit(".", 1)[0].strip()
            if sid not in wide_ids:
                wide_ids.append(sid)
            key = {"gtype": "gt", "genotype": "gt", "logrratio": "lrr",
                   "ballelefreq": "baf", "copynumber": "cn"}[field]
            wide_cols.setdefault(sid, {})[key] = j
        if not wide_ids and "gt" not in idx and not ("a1" in idx and "a2" in idx):
            # A third layout: one column per sample named simply by the sample, as HapMap and
            # several public releases write it. Any column that is not a recognised field and
            # whose values look like two-character calls is taken to be a sample.
            claimed = set(idx.values())
            probe_row = next((l.split(delim) for l in lines[hdr_i + 1:]
                              if len(l.split(delim)) >= width), None)
            if probe_row:
                for j, c in enumerate(raw_cols):
                    if j in claimed or j >= len(probe_row):
                        continue
                    v = probe_row[j].strip()
                    if len(v) == 2 and (set(v) <= set("ACGTN") or v in ("AA", "AB", "BB", "--")):
                        wide_ids.append(c.strip())
                        wide_gt_col[c.strip()] = j

    out: dict[str, dict[str, Probe]] = {}
    for line in lines[hdr_i + 1:]:
        f = line.rstrip("\n").split(delim)
        if len(f) < width:
            continue

        def num(key: str) -> Optional[float]:
            if key not in idx:
                return None
            try:
                v = float(f[idx[key]])
            except (ValueError, IndexError):
                return None
            return v if math.isfinite(v) else None

        try:
            pos = int(float(f[idx["pos"]]))
        except (ValueError, IndexError):
            continue
        chrom = f[idx["chrom"]].strip().removeprefix("chr")
        rsid = f[idx["rsid"]]
        if wide_gt_col:
            for sid, j in wide_gt_col.items():
                g = f[j].strip()
                out.setdefault(sid, {})[rsid] = Probe(
                    chrom=chrom, pos=pos,
                    gt=_GT.get(g, g if _is_nucleotide(g) else "NC"))
            continue
        if len(wide_ids) > 1 and wide_cols:
            for sid, cmap in wide_cols.items():
                if "gt" not in cmap:
                    continue
                g = f[cmap["gt"]].strip()

                def wnum(key: str) -> Optional[float]:
                    if key not in cmap:
                        return None
                    try:
                        v = float(f[cmap[key]])
                    except (ValueError, IndexError):
                        return None
                    return v if math.isfinite(v) else None

                out.setdefault(sid, {})[rsid] = Probe(
                    chrom=chrom, pos=pos,
                    gt=_GT.get(g, g if _is_nucleotide(g) else "NC"),
                    baf=wnum("baf"), lrr=wnum("lrr"), cn=wnum("cn"))
            continue
        if "gt" in idx:
            raw_gt = f[idx["gt"]].strip()
        elif "a1" in idx and "a2" in idx:
            raw_gt = (f[idx["a1"]].strip() + f[idx["a2"]].strip()).replace("-", "")
        else:
            continue
        sid = f[idx["sample"]].strip() if long_mode else (wide_ids[0] if wide_ids
                                                          else path.name.split(".")[0])
        out.setdefault(sid, {})[rsid] = Probe(
            chrom=chrom, pos=pos,
            gt=_GT.get(raw_gt, raw_gt if _is_nucleotide(raw_gt) else "NC"),
            baf=num("baf"), lrr=num("lrr"), cn=num("cn"),
        )

    if not out or not any(out.values()):
        raise ValueError(
            f"{path.name}: found marker columns but no genotypes. Expected a genotype column, "
            f"a pair of allele columns, per-sample columns named <sample>.GType, or columns "
            f"named for the samples themselves. Header: {', '.join(cols[:20])}")
    return {k: Sample(k, v, cols) for k, v in out.items()}


def harmonise(samples: Sequence[Sample]) -> list[Sample]:
    """Put every sample into one consistent AB space, agreed across all of them.

    Nucleotide calls cannot be mapped file by file. A marker whose alleles are A and G shows only
    "AA" in a sample homozygous for the first, so that file alone cannot know G exists, and two
    files mapping independently would disagree about which allele is A. Pooling the alleles seen
    at each marker across every sample being compared fixes one mapping for all of them.

    Which allele becomes A never matters, only that it is the same everywhere, because every test
    downstream is invariant under a consistent relabel. That also makes this a no-op on data
    already in AB space: the alleles there are A and B, sorted order maps A to A and B to B.
    """
    alleles: dict[str, set[str]] = {}
    for s in samples:
        for rsid, p in s.probes.items():
            if p.gt in ("NC", "", "--", "00", "NN"):
                continue
            alleles.setdefault(rsid, set()).update(p.gt)

    out = []
    for s in samples:
        probes = {}
        for rsid, p in s.probes.items():
            seen = sorted(alleles.get(rsid, ()))
            if len(seen) > 2 or p.gt in ("NC", "", "--", "00", "NN") or not seen:
                gt = "NC"
            else:
                m = {seen[0]: "A", (seen[1] if len(seen) > 1 else seen[0]): "B"}
                if len(seen) == 1:
                    m = {seen[0]: "A"}
                letters = sorted(m.get(c, "?") for c in p.gt)
                gt = "".join(letters) if "?" not in letters else "NC"
                gt = {"AA": "AA", "AB": "AB", "BB": "BB"}.get(gt, "NC")
            probes[rsid] = Probe(p.chrom, p.pos, gt, baf=p.baf, lrr=p.lrr, cn=p.cn)
        out.append(Sample(s.sample_id, probes, s.columns, list(s.notes)))
    return out


def _is_nucleotide(g: str) -> bool:
    return len(g) == 2 and all(c in "ACGT" for c in g)


def read_sample(path: str | Path, sample_id: Optional[str] = None) -> Sample:
    """One sample out of a file. See `read_samples` for the format families handled."""
    path = Path(path)
    if any(path.name.endswith(x) for x in (".vcf", ".vcf.gz")):
        return read_vcf(path, sample_id)
    got = read_samples(path)
    if sample_id and sample_id in got:
        return got[sample_id]
    if sample_id and len(got) > 1:
        raise ValueError(f"{Path(path).name}: no sample {sample_id!r}. "
                         f"Present: {', '.join(sorted(got))}")
    s = next(iter(got.values()))
    return Sample(sample_id or s.sample_id, s.probes, s.columns, s.notes)


def read_vcf(path: str | Path, sample_name: Optional[str] = None) -> Sample:
    """Read one sample out of a VCF, with read-level evidence where the file carries it.

    Sequencing supplies two things an array cannot.

    AD gives the reads supporting each allele, so the alternate fraction is a MEASUREMENT rather
    than a three-way call. A genotype rounds 8% alternate reads to homozygous reference and the
    evidence of a minority lineage is gone; the fraction keeps it.

    DP gives depth, and depth is LINEAR in copy number. `log2(DP / median DP)` is an LRR with a
    compression of exactly one, so the whole per-array calibration that intensity needs - the
    chrX probe offset, the compression constant, the credible bounds - simply does not arise.

    REF becomes A and ALT becomes B. Which is which never matters, only that it is consistent
    within a marker, because every test downstream is invariant under that relabel.
    """
    path = Path(path)
    opener = gzip.open if path.name.endswith(".gz") else open
    gts: dict[str, tuple[str, str, int, Optional[float], Optional[float]]] = {}
    col: Optional[int] = None
    sid = sample_name or path.name.split(".")[0]
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                if len(f) < 10:
                    raise ValueError(f"{path.name}: no genotype columns in this VCF")
                names = f[9:]
                col = 9 + (names.index(sample_name) if sample_name in names else 0)
                sid = sample_name if sample_name in names else names[0]
                continue
            if col is None or len(f) <= col or "," in f[4]:
                continue          # header not seen yet, or multi-allelic: this tool is biallelic
            keys = f[8].split(":")
            vals = f[col].split(":")
            rec = dict(zip(keys, vals))
            g = rec.get("GT", "./.").replace("|", "/")
            a, _, b = g.partition("/")
            if a in ("0", "1") and b in ("0", "1"):
                gt = {"00": "AA", "01": "AB", "10": "AB", "11": "BB"}[a + b]
            else:
                gt = "NC"
            frac = depth = None
            ad = rec.get("AD")
            if ad and "," in ad:
                try:
                    ref, alt = (int(x) for x in ad.split(",")[:2])
                except ValueError:
                    ref = alt = 0
                if ref + alt > 0:
                    frac = alt / (ref + alt)
            try:
                depth = float(rec["DP"])
            except (KeyError, ValueError):
                depth = None
            rsid = f[2] if f[2] not in (".", "") else f"{f[0]}:{f[1]}"
            gts[rsid] = (f[0].removeprefix("chr"), gt, int(f[1]), frac, depth)

    if not gts:
        raise ValueError(f"{path.name}: no usable biallelic records")
    depths = [d for *_, d in gts.values() if d and d > 0]
    med = st.median(depths) if depths else None
    probes = {}
    for rsid, (chrom, gt, pos, frac, depth) in gts.items():
        lrr = math.log2(depth / med) if (depth and med and depth > 0) else None
        probes[rsid] = Probe(chrom, pos, gt, baf=frac, lrr=lrr)

    notes = []
    hom_ref = sum(1 for p in probes.values() if p.gt == "AA")
    if hom_ref < 0.01 * len(probes):
        notes.append(
            f"Only {hom_ref:,} of {len(probes):,} records are homozygous reference, so this looks "
            "like a variant-only VCF rather than an all-sites or jointly-called one. Parent-of-"
            "origin needs the sites where the father is homozygous, including homozygous "
            "REFERENCE, and those are missing here. Use a joint call across the family, or a "
            "gVCF, or the marker set is conditioned on which sample happened to vary."
        )
    return Sample(sid, probes, ["vcf"], notes)


@dataclass
class BuildCall:
    build: Optional[str]
    tested: int
    illegal: dict[str, int]
    note: str


def detect_build(sample: Sample, *, min_markers: int = 1_000) -> BuildCall:
    """Which assembly these coordinates belong to, from the positions alone.

    A marker cannot lie inside an assembly N-gap or past the end of a chromosome, so counting
    illegal placements under each candidate build separates them: the right build gives zero.
    No rsIDs, no manifest, no chain file.
    """
    illegal = {b: 0 for b in buildref.BUILDS}
    tested = 0
    for p in sample.probes.values():
        c = p.chrom.removeprefix("chr")
        if c not in buildref.CHROM_LEN["GRCh37"]:
            continue
        tested += 1
        for b in buildref.BUILDS:
            if p.pos > buildref.CHROM_LEN[b][c]:
                illegal[b] += 1
                continue
            for lo, hi in buildref.GAPS[b].get(c, ()):
                if lo < p.pos <= hi:
                    illegal[b] += 1
                    break
    if not tested:
        return BuildCall(None, 0, illegal, "no markers on a recognised primary chromosome")
    clean = [b for b in buildref.BUILDS if illegal[b] == 0]
    if len(clean) == 1:
        if tested < min_markers:
            return BuildCall(None, tested, illegal,
                             f"only {tested} markers tested: too few to call")
        return BuildCall(clean[0], tested, illegal,
                         f"{clean[0]}: every one of {tested:,} markers is legal there, "
                         f"{illegal[[b for b in buildref.BUILDS if b != clean[0]][0]]:,} are not "
                         "under the alternative")
    if len(clean) > 1:
        return BuildCall(None, tested, illegal,
                         "markers are legal under both builds, so positions cannot separate them")
    best = min(buildref.BUILDS, key=lambda b: illegal[b])
    other = [b for b in buildref.BUILDS if b != best][0]
    return BuildCall(None, tested, illegal,
                     f"no build is clean ({illegal[best]:,} illegal under {best}, "
                     f"{illegal[other]:,} under {other}): coordinates may be a build this table "
                     "does not carry, or a different reference altogether")


def call_sex(s: Sample) -> tuple[str, float]:
    """chrX heterozygosity as a fraction of autosomal. Needs no reference panel.

    The middle band returns "unknown" rather than picking a side, which matters because the chrX
    calibration route is only valid in a male and a wrong sex call silently doubles the scale.
    """
    def het(values: Sequence[str]) -> float:
        called = [g for g in values if g != "NC"]
        return sum(g == "AB" for g in called) / len(called) if called else math.nan

    x = het([p.gt for p in s.probes.values() if p.chrom in ("X", "23")])
    a = het([p.gt for p in s.probes.values() if _is_autosome(p.chrom)])
    ratio = x / a if math.isfinite(x) and math.isfinite(a) and a > 0 else math.nan
    if not math.isfinite(ratio):
        return "unknown", ratio
    return ("male" if ratio < 0.15 else "female" if ratio > 0.55 else "unknown"), ratio


def calibrate_sample(
    s: Sample, *, product: str, declared_compression: Optional[float] = None
) -> tuple[nz.Calibration, str, float]:
    cn_obs = [(int(p.cn), p.lrr) for p in s.probes.values()
              if p.cn is not None and p.lrr is not None and p.cn == int(p.cn)]
    sex, ratio = call_sex(s)
    cal = nz.calibrate(
        s.autosomal("lrr"), s.chrx("lrr"),
        sex=sex, product=product, copy_number_obs=cn_obs or None,
        declared_compression=declared_compression,
    )
    return cal, sex, ratio


def build_markers(
    father: Sample,
    embryo: Sample,
    mother: Optional[Sample] = None,
    *,
    chrom: Optional[str] = None,
) -> list[Marker]:
    """Inner join on probe id, sorted by position, which is what the run statistic requires."""
    out = []
    for rsid, e in embryo.probes.items():
        if chrom is not None and e.chrom != chrom:
            continue
        fa = father.probes.get(rsid)
        if fa is None:
            continue
        mo = mother.probes.get(rsid) if mother else None
        out.append(Marker(rsid, e.chrom, e.pos, fa.gt,
                          mo.gt if mo else None, e.gt, e.baf, e.lrr))
    out.sort(key=lambda m: (len(m.chrom), m.chrom, m.pos))
    return out


# --- parental origin from sperm and sample alone ---------------------------------------------
#
# Did this sample get a paternal genome, and on which chromosomes. A RATE over hundreds of
# thousands of markers, so it needs no mother, no variant position and no phase. It also covers
# the case chrY cannot reach: an X-bearing sperm leaves no Y, so chrY alone cannot separate a full
# paternal genome from none.

#: For REPORTING only, never a decision boundary: measured across 50+ pairs it ranges from 6.8%
#: to 50%, since it depends on the ancestry the two people share and on the father file's noise.
UNRELATED_ABSENCE = {"UK Biobank Axiom Array": 0.055}

#: Above this a chromosome is treated as having lost the paternal contribution outright. Sits far
#: from both populations seen in practice: 0.1-0.6% when present, 5-8% when absent.
DEFAULT_UNRELATED_ABSENCE = 0.055

def second_parent_signal(father_heterozygosity: float) -> float:
    """Rate at which a SECOND parent makes the sample carry an allele the father lacks.

    Derived, not fitted. The father lacks an allele only where he is homozygous, and a second
    parent supplies the other one at the population frequency, so the rate is the sum over
    markers of p^2 q + q^2 p = pq. The father's own heterozygosity is the sum of 2pq over the
    same markers, so the quantity wanted is exactly HALF HIS HETEROZYGOSITY, measurable from his
    own array with no reference panel and no ancestry assumption.

    Checked against six pairs of known composition: 8.48% predicted against 3.4-5.2% observed for
    three paternal-only genomes, and 8.33/16.10/9.72% predicted against 34.66/16.35/10.65%
    observed for three biparental ones.

    The margin here is far narrower than on the absence axis, roughly 1.6x rather than 30x,
    because a clean biparental sample lands almost exactly on the prediction while noise pushes a
    paternal-only sample up toward it. So this axis is used only where it has to be, which is a
    DIPLOID sample: a homozygous genome carrying the father's alleles has one allele per locus and
    that allele is his, leaving no room for a maternal complement whatever the residual says.
    """
    if not math.isfinite(father_heterozygosity):
        return math.nan
    return father_heterozygosity / 2.0

#: Fraction of autosomal BAF in the heterozygous band, which says whether the genome is diploid.
#: 15-16% for diploids, 1.3-3.4% for uniparental. Dropout smears that band rather than emptying
#: it: a biparental embryo at 33% dropout still reads 22.2%.
HET_BAND_DIPLOID = 0.08

#: Call rate below which heterozygous calls stop being trustworthy, so nothing derived from them
#: may be asserted. Natesan 2014, the only published threshold measured on amplified material,
#: and the same figure the ingestion gates already exclude on.
#:
#: Absence is unaffected and still called: it is Mendelian, and a homozygous father must transmit
#: his allele whatever the call rate. What is withheld is ZYGOSITY, read from the heterozygous
#: band, and with it the androgenetic/biparental split that rests on it. Measured on three
#: isolated paternal pronuclei (GSM4774684/686/687) at 53.8%, 55.6% and 59.1% call rate: each is
#: haploid and therefore homozygous by construction, each showed an 18-27% heterozygous band that
#: a haploid genome cannot produce, and each was called biparental on the strength of it.
CALL_RATE_FLOOR = 0.60

#: NOTE on the per-chromosome likelihood ratio below, from an independent audit: it ranks two
#: hypotheses and never asks whether the winner fits, so in principle it can report a decisive
#: verdict on data neither describes. A goodness-of-fit gate was written for it and removed. The
#: ratio only reaches +/-3 when the observation is at or BEYOND one of the two rates, and beyond
#: "unrelated" in the absent direction is a missing chromosome, which is the finding this whole
#: layer exists to make: the gate refused a real whole-chromosome loss at 31.5% absence for not
#: resembling the 5.5% unrelated rate closely enough. The audit's concrete case, fifteen
#: chromosomes called on a mis-clustered array, is caught instead by the per-chromosome
#: heterozygosity gate in ingest.ts, which is where that failure actually lives.

#: Parent heterozygosity below which the marker set is not a polymorphic panel.
#:
#: Every rate here is calibrated on common-SNP arrays, where a parent runs 15-19% heterozygous. A
#: whole-genome variant callset runs far lower because most of its sites are rare: measured at
#: 3.2% on a 1000 Genomes chr22 VCF, where the parent is homozygous reference at 95% of markers
#: against 75% on an array of the same chromosome. Those sites cannot show absence in anyone, so
#: they dilute the denominator and pull an unrelated pair toward a related one: the pair reading
#: 4.1% absence on the array reads 0.69% in the callset. Annotated, never adjusted for.
PANEL_HET_FLOOR = 0.10

#: Residual paternal-absence rate on clean data, from genotyping error alone. Measured at 0.03%
#: and 0.05% on two bulk trios.
ABSENCE_ERROR_FLOOR = 0.005

#: How far above the explainable level the observed rate must sit before absence is called.
ABSENCE_MARGIN = 3.0


def absence_explainable(no_call_rate: float, het_fraction: float) -> float:
    """The highest paternal-absence rate a TRUE father-offspring pair can produce by noise alone.

    Dropout turns a heterozygous call homozygous, and at a marker where the father is homozygous
    the embryo's heterozygous genotype is exactly paternal-allele plus maternal-allele. Half the
    time the allele dropout discards is his, and the call then reads as his allele being absent.
    So the inflation scales with heterozygosity TIMES dropout, and neither alone: a homozygous
    genome is immune however much drops out, which is why androgenotes hold at 0.16% despite a
    13.6% no-call rate, while a diploid embryo at 46.9% no-call reaches 9.69% against its own
    father and overtakes genuinely unrelated pairs at 5.56%.

    Without a mother the dropout rate cannot be measured, so the no-call rate stands in for it.
    The product bounds every related pair measured, and tightly: predicted 0.13% against observed
    0.05%, 0.18% against 0.16%, 0.27% against 0.22%, and 10.4% against 9.69%.

    This replaces comparing against an "unrelated" rate, which cannot be a constant: it depends on
    the ancestry shared by the two people and ranged from 6.8% to 50% across the pairs measured.
    Asking whether the absence exceeds what this sample's own noise can produce needs no
    population reference at all.
    """
    if not (math.isfinite(no_call_rate) and math.isfinite(het_fraction)):
        return math.nan
    return max(0.0, no_call_rate) * max(0.0, het_fraction)


@dataclass
class ChromOrigin:
    chrom: str
    informative: int
    absent: int
    rate: float
    verdict: Literal["paternal_present", "paternal_absent", "expected_absent", "unclear"]
    log10_lr: float
    note: str = ""


@dataclass
class ParentageReport:
    sample_id: str
    father_id: str
    genome_rate: float
    informative: int
    baseline: float
    unrelated_rate: float
    y_called: int
    y_total: int
    y_bearing_sperm: bool
    verdict: Literal["paternal_genome_present", "no_paternal_contribution", "unclear"]
    #: Second axis: alleles the father cannot supply, so a non-paternal genome is present.
    #: Things that CONSTRAIN the answer, as opposed to notes, which explain it. Kept apart so a
    #: Methods section does not list a finding as a limitation.
    limits: list[str] = field(default_factory=list)
    nonpaternal_rate: float = math.nan
    #: Third axis: is the genome diploid at all, from the BAF heterozygous band.
    het_band: float = math.nan
    zygosity: Literal["diploid", "uniparental_homozygous", "unknown"] = "unknown"
    #: Highest absence rate this sample's own noise could produce in a TRUE father-offspring
    #: pair, plus the error floor. The decision boundary, and it needs no population reference.
    explainable: float = math.nan
    no_call_rate: float = math.nan
    #: Half the father's heterozygosity: what a second parent would contribute on axis 2.
    second_parent_expected: float = math.nan
    origin_class: Literal["androgenetic", "gynogenetic", "biparental", "unclear"] = "unclear"
    #: Inferred from the chrX SNP rate, not from chrY.
    sperm_type: Literal["X_bearing", "Y_bearing", "unknown"] = "unknown"
    chroms: list[ChromOrigin] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"{self.sample_id} vs {self.father_id}: {self.origin_class.upper()}",
                 f"  paternal alleles absent      {self.genome_rate:>7.2%}   "
                 f"(this sample's noise can explain up to {self.explainable:.2%}; "
                 f"an unrelated genome runs {self.unrelated_rate:.0%}+)",
                 f"  alleles the father lacks     {self.nonpaternal_rate:>7.2%}   "
                 f"(a second parent would contribute {self.second_parent_expected:.2%}, "
                 f"half the father's heterozygosity)",
                 f"  heterozygous BAF band        {self.het_band:>7.2%}   "
                 f"-> {self.zygosity}",
                 f"  based on {self.informative:,} autosomal markers where the father is "
                 f"homozygous"]
        if self.sperm_type != "unknown":
            lines.append(f"  sperm type {self.sperm_type}, from the chrX SNP rate. "
                         f"chrY {self.y_called:,}/{self.y_total:,} called, as a cross-check only.")
        for c in self.chroms:
            if c.verdict == "paternal_present":
                continue
            lines.append(f"  chr{c.chrom}: {c.rate:.1%} absent over {c.informative:,} markers "
                         f"-> {c.verdict.upper()}" + (f"  ({c.note})" if c.note else ""))
        for n in self.notes:
            lines.append(f"  {n}")
        return "\n".join(lines)


def parental_origin(
    father: Sample, sample: Sample, *, product: str = "unspecified",
    decisive_log10: float = 3.0, role: Literal["paternal", "maternal"] = "paternal",
) -> ParentageReport:
    """Did this sample inherit a paternal genome, and on which chromosomes.

    At every marker where the father is homozygous he must transmit that allele, so a sample
    lacking it is missing his contribution there. One such marker means nothing; the rate across
    a chromosome means a great deal, and it needs neither the mother nor a variant position.
    """
    unrelated = UNRELATED_ABSENCE.get(product, DEFAULT_UNRELATED_ABSENCE)
    fa = {k: p for k, p in father.probes.items() if p.gt in ("AA", "BB")}

    per: dict[str, list[int]] = {}
    nonpat = nonpat_den = 0
    for k, p in sample.probes.items():
        fp = father.probes.get(k)
        if fp is None or fp.gt == "NC" or p.gt == "NC":
            continue
        if _is_autosome(fp.chrom):
            nonpat_den += 1
            # An allele the father does not carry at all had to come from somewhere else. This
            # is the only measure that separates a paternal-only genome from a biparental one:
            # both have his alleles present, so paternal absence cannot tell them apart.
            nonpat += bool(set(p.gt) - set(fp.gt))
        f = fa.get(k)
        if f is None:
            continue
        d = per.setdefault(f.chrom, [0, 0])
        d[0] += 1
        d[1] += f.gt[0] not in p.gt

    band = [p.baf for p in sample.probes.values()
            if p.baf is not None and _is_autosome(p.chrom)]
    het_band = (sum(1 for b in band if 0.35 <= b <= 0.65) / len(band)) if band else math.nan

    y_total = sum(1 for p in father.probes.values() if p.chrom == "Y")
    y_called = sum(1 for k, p in sample.probes.items() if p.chrom == "Y" and p.gt != "NC")
    y_bearing = y_total > 0 and y_called / y_total > 0.5

    autosomes = {c: v for c, v in per.items() if _is_autosome(c) and v[0] >= 200}
    n_tot = sum(v[0] for v in autosomes.values())
    a_tot = sum(v[1] for v in autosomes.values())
    genome_rate = a_tot / n_tot if n_tot else math.nan
    # The present-rate baseline comes from the sample's own best chromosomes, so it absorbs this
    # sample's genotyping quality instead of assuming a fixed error rate.
    rates = sorted(v[1] / v[0] for v in autosomes.values())
    baseline = max(rates[len(rates) // 4] if rates else 0.002, 5e-4)

    rep = ParentageReport(
        sample_id=sample.sample_id, father_id=father.sample_id, genome_rate=genome_rate,
        informative=n_tot, baseline=baseline, unrelated_rate=unrelated,
        y_called=y_called, y_total=y_total, y_bearing_sperm=y_bearing,
        verdict="unclear",
        nonpaternal_rate=(nonpat / nonpat_den) if nonpat_den else math.nan,
        het_band=het_band,
    )
    if not n_tot:
        rep.limits.append("No markers where the father is homozygous and the sample is called.")
        return rep

    # The baseline comes from the sample's own quietest chromosomes, which has one total failure
    # mode: a sample with no paternal genome has no quiet chromosome to learn from. Detect that
    # rather than quoting a "when present" rate above the unrelated one. The decision itself is
    # whether absence exceeds what this sample's own noise can manufacture, which needs no
    # population reference.
    called = [p.gt for p in sample.probes.values() if p.gt != "NC"]
    nc_rate = 1.0 - (len(called) / len(sample.probes)) if sample.probes else math.nan
    het_frac = (het_band if math.isfinite(het_band)
                else (sum(1 for g in called if g == "AB") / len(called) if called else math.nan))
    explainable = absence_explainable(nc_rate, het_frac) + ABSENCE_ERROR_FLOOR
    rep.explainable = explainable
    rep.no_call_rate = nc_rate
    # Used by both the presence guard and the zygosity guard below.
    call_rate = 1.0 - nc_rate if math.isfinite(nc_rate) else math.nan
    rep.limits.extend(sample.notes)
    rep.limits.extend(f"father: {n}" for n in father.notes)

    if not math.isfinite(explainable):
        rep.verdict = "unclear"
        rep.limits.append("Sample quality could not be measured, so there is no bound on "
                         "noise-driven absence and no call is made.")
    elif genome_rate <= explainable:
        rep.verdict = "paternal_genome_present"
    elif genome_rate >= explainable * ABSENCE_MARGIN:
        rep.verdict = "no_paternal_contribution"
        rep.baseline = math.nan
    else:
        rep.verdict = "unclear"
        rep.limits.append(
            f"Absence of {genome_rate:.2%} sits between what this sample's own noise can "
            f"explain ({explainable:.2%}) and {ABSENCE_MARGIN:g}x that. A noisy diploid sample "
            "can reach the unrelated range against its own father, so this is left uncalled "
            "rather than guessed. Strongly recommended rather than required: the oocyte donor's "
            "array measures dropout directly and would settle this case, and every other "
            "borderline one, without changing what the tool can already call without her."
        )
    no_clean_chromosome = rep.verdict != "paternal_genome_present"

    for c in (() if no_clean_chromosome else sorted(per, key=lambda x: (len(x), x))):
        n, a = per[c]
        if n < 200:
            continue
        rate = a / n
        # A male offspring has no paternal X at all: his father sent a Y instead. Flagging that
        # as a loss would report ordinary sex determination as an anomaly.
        expected = (c in ("X", "23") and y_bearing
                    and rep.verdict == "paternal_genome_present")
        clr = (_log_binom_pmf(a, n, baseline) - _log_binom_pmf(a, n, unrelated)) / math.log(10.0)
        if expected and rate > baseline * 5:
            verdict, note = "expected_absent", "no paternal X: this sample carries the father's Y"
        elif clr <= -decisive_log10:
            verdict, note = "paternal_absent", ""
        elif clr >= decisive_log10:
            verdict, note = "paternal_present", ""
        else:
            verdict, note = "unclear", "rate sits between the two references"
        rep.chroms.append(ChromOrigin(c, n, a, rate, verdict, clr, note))

    # --- the three axes, resolved into one parent-of-origin call ----------------------------
    #
    # Absence cannot separate a paternal-only genome from a biparental one, since his alleles are
    # in both. The second axis does that, and the third says whether the genome is diploid.
    pat_present = rep.verdict == "paternal_genome_present"
    # The second-parent signal is DERIVED from the father's own heterozygosity, not fitted.
    father_het = (sum(1 for p in father.probes.values()
                      if p.gt == "AB" and _is_autosome(p.chrom))
                  / max(1, sum(1 for p in father.probes.values()
                               if p.gt != "NC" and _is_autosome(p.chrom))))
    rep.second_parent_expected = second_parent_signal(father_het)

    if math.isfinite(father_het) and father_het < PANEL_HET_FLOOR:
        rep.limits.append(
            f"The parent is heterozygous at {father_het:.1%} of called autosomal markers. Every "
            "rate here is calibrated on common-SNP arrays, where that figure runs 15-19%. A marker "
            "set this monomorphic is a whole-genome variant callset rather than a polymorphic "
            "panel: most of its sites are rare, they cannot show absence in anyone, and they "
            "dilute the denominator. Related and unrelated pairs move closer together, so a call "
            "here is weaker than the same call on an array and an unrelated pair may read as "
            "unclear. Restrict such a file to common variants before relying on it."
        )

    # The gate the ingestion layer already applies, applied to the axis that depends on it:
    # below this call rate the heterozygous band is artefact, so zygosity is withheld.
    if math.isfinite(call_rate) and call_rate < CALL_RATE_FLOOR:
        rep.limits.append(
            f"Call rate {call_rate:.1%} is below {CALL_RATE_FLOOR:.0%}, where erroneous "
            "heterozygous calls become common. Zygosity is read from the heterozygous band, so "
            "it is not reported here and the genome cannot be separated into uniparental or "
            "biparental. The presence or absence of this parent's contribution is unaffected: "
            "that is Mendelian and does not rest on heterozygous calls."
        )
    else:
        if math.isfinite(het_band):
            rep.zygosity = "diploid" if het_band > HET_BAND_DIPLOID else "uniparental_homozygous"
        elif called and math.isfinite(father_het) and father_het > 0:
            # No intensity channel at all, so fall back to genotype heterozygosity against the same
            # derived reference: a biparental sample resembles the father, a uniparental one tends to
            # zero. Weaker than the BAF band, which is why it is flagged.
            gt_het = sum(1 for g in called if g == "AB") / len(called)
            rep.zygosity = ("diploid" if gt_het > father_het / 2.0
                            else "uniparental_homozygous")
            rep.limits.append(
                f"No B-allele frequencies in this file, so zygosity comes from genotype "
                f"heterozygosity ({gt_het:.1%} against the father's {father_het:.1%}) rather than "
                "the BAF band. That is the weaker of the two measures."
            )

    if rep.zygosity == "uniparental_homozygous":
        # One allele per locus. If those alleles are the father's there is no room for a maternal
        # complement, so this axis is not consulted and its residual is error by construction.
        rep.origin_class = ("androgenetic" if pat_present else
                            "gynogenetic" if rep.verdict == "no_paternal_contribution"
                            else "unclear")
    elif rep.zygosity == "diploid" and pat_present:
        # The only case where the second axis has to carry the decision, and its margin is
        # narrow, so the number and its reference are both reported.
        if not math.isfinite(rep.second_parent_expected):
            rep.origin_class = "unclear"
        elif rep.nonpaternal_rate > rep.second_parent_expected:
            rep.origin_class = "biparental"
        else:
            rep.origin_class = "androgenetic"
            rep.notes.append(
                f"Diploid but carrying {rep.nonpaternal_rate:.2%} alleles the father lacks, "
                f"below the {rep.second_parent_expected:.2%} a second parent would contribute "
                "(half the father's heterozygosity). Consistent with a paternal-only genome that "
                "is heterozygous rather than duplicated, meaning two sperm. The margin on this "
                "axis is about 1.6x, far narrower than the 30x on paternal absence, so treat a "
                "near-boundary call as provisional."
            )
        margin = (rep.nonpaternal_rate / rep.second_parent_expected
                  if rep.second_parent_expected else math.nan)
        if math.isfinite(margin) and 0.7 < margin < 1.4:
            rep.limits.append(
                f"Second-parent signal sits within 40% of its boundary ({rep.nonpaternal_rate:.2%} "
                f"against {rep.second_parent_expected:.2%}). A clean biparental sample lands "
                "almost exactly on this prediction, so this call is not decisive on its own."
            )
    elif rep.verdict == "no_paternal_contribution":
        rep.origin_class = "gynogenetic"
    else:
        rep.origin_class = "unclear"

    # Sperm type from the chrX SNP rate, so chrY is confirmation rather than the method. An
    # X-bearing sperm delivers the father's X and leaves no Y, which is exactly the case chrY
    # cannot call.
    xr = next((c for c in rep.chroms if c.chrom in ("X", "23")), None)
    # Sperm type is a paternal question only. A mother transmits an X to a child of either sex,
    # so an absent maternal X is a real loss rather than ordinary sex determination.
    if role == "maternal" and xr is not None and xr.verdict == "expected_absent":
        xr.verdict, xr.note = "paternal_absent", "maternal X absent: not explained by sex"
    if role == "paternal" and pat_present and xr is not None:
        # An absent paternal X has at least three causes: a Y-bearing sperm, X loss, and local
        # assay failure. Concluding the first from the absence alone is circular, so the chrY
        # measurement decides it and "unknown" is a legal answer. The browser has always required
        # this; the CLI used to infer Y_bearing from the X rate and never consult chrY.
        x_gone = xr.rate > (baseline + unrelated) / 2
        rep.sperm_type = ("X_bearing" if not x_gone
                          else "Y_bearing" if y_bearing else "unknown")
        if rep.sperm_type == "unknown":
            rep.notes.append(
                "The paternal X is absent but no chrY was called, so this is not ordinary sex "
                "determination and no sperm type is claimed. A Y-bearing sperm, loss of the X, "
                "and assay failure on chrX all produce this and the array cannot separate them."
            )
        if rep.sperm_type == "X_bearing" and not y_bearing:
            rep.notes.append(
                "No chrY, yet the father's alleles are present on chrX as well as the autosomes: "
                "an X-bearing sperm. chrY alone cannot distinguish this from no paternal "
                "contribution at all, which is the gap the SNP rate closes."
            )
        if rep.sperm_type == "Y_bearing" and not y_bearing:
            rep.notes.append(
                "chrX shows no paternal contribution but chrY was not called either. Those "
                "disagree, so the sperm type is not settled."
            )
    if rep.origin_class == "gynogenetic":
        rep.notes.append(
            "No paternal genome, and the genome is "
            + ("homozygous, consistent with a duplicated maternal complement."
               if rep.zygosity == "uniparental_homozygous" else "diploid.")
            + " chrY is silent here too, so chrY alone would not have separated this from an "
              "X-bearing sperm carrying a full paternal genome."
        )
    if rep.origin_class == "androgenetic":
        rep.notes.append(
            "Every allele traces to the father and the genome is homozygous: a paternal-only "
            "complement, duplicated. Nothing here required chrY."
        )
    return rep


# --- segmental parent of origin ----------------------------------------------------------------
#
# Which PARTS of a paternal genome arrived, for a mosaic or partially-lost sample. The same
# measurement, resolved along the chromosome instead of summed over it.


@dataclass
class OriginSegment:
    chrom: str
    start: int
    end: int
    markers: int
    absent: int
    rate: float
    state: Literal["paternal_present", "paternal_absent"]
    #: How far the rate sits above the level this sample's own noise can reach. A real loss runs
    #: far above it; repeat-rich subtelomeric and pericentromeric regions sit just over the line,
    #: which is where correlated probe failure shows up, so the two must not be reported alike.
    margin: float = math.nan

    @property
    def span_mb(self) -> float:
        return (self.end - self.start) / 1e6

    @property
    def confident(self) -> bool:
        return math.isfinite(self.margin) and self.margin >= 3.0


def _viterbi_two_state(
    obs: Sequence[tuple[int, int]], p_lo: float, p_hi: float, length_bp: float
) -> list[int]:
    """States: 0 = rate p_lo (paternal present), 1 = rate p_hi (paternal absent)."""
    logs = []
    for p in (p_lo, p_hi):
        logs.append((math.log(max(p, 1e-12)), math.log(max(1.0 - p, 1e-12))))
    best = [math.log(0.5), math.log(0.5)]
    back: list[tuple[int, int]] = []
    prev_pos = obs[0][0]
    for i, (pos, a) in enumerate(obs):
        gap = max(0, pos - prev_pos)
        prev_pos = pos
        # Longer physical gaps make a state change more likely, which is what keeps boundaries
        # tied to the chromosome rather than to marker density.
        sw = -math.expm1(-gap / length_bp) if i else 0.5
        sw = min(max(sw, 1e-9), 0.5)
        stay, move = math.log(1.0 - sw), math.log(sw)
        cur, ptr = [0.0, 0.0], [0, 0]
        for s in (0, 1):
            same, other = best[s] + stay, best[1 - s] + move
            ptr[s] = s if same >= other else 1 - s
            cur[s] = max(same, other) + (logs[s][0] if a else logs[s][1])
        best, _ = cur, back.append((ptr[0], ptr[1]))
    path = [0 if best[0] >= best[1] else 1]
    for ptr in reversed(back[1:]):
        path.append(ptr[path[-1]])
    return list(reversed(path))


def segmental_origin(
    father: Sample, sample: Sample, *, expected_segment_mb: float = 20.0,
    min_markers: int = 100, window: int = 200, explainable: Optional[float] = None,
) -> tuple[list[OriginSegment], dict]:
    """Where along the genome the paternal contribution is present, and where it is not.

    The two emission rates are estimated from this sample rather than assumed: the quiet windows
    give what "present" looks like at this sample's noise level, and the loud ones give what
    "absent" looks like for this pair. That keeps the ancestry and quality dependence out, exactly
    as on the genome-wide call.

    Resolution is set by how many markers a rate needs to separate the two, not chosen. On the
    Axiom array one informative marker per 5 kb makes 200 markers about 1 Mb, which distinguishes
    0.16% from 6.8% at odds near 10^17.
    """
    per: dict[str, list[tuple[int, int]]] = {}
    for k, p in sample.probes.items():
        f = father.probes.get(k)
        if f is None or f.gt not in ("AA", "BB") or p.gt == "NC":
            continue
        per.setdefault(f.chrom, []).append((f.pos, int(f.gt[0] not in p.gt)))
    for v in per.values():
        v.sort()

    # Anchored to the calibrated boundary, not to percentiles of this sample. Percentiles find
    # two populations in a uniform genome, since those are the tails of one noise distribution: it
    # carved a sample present everywhere at 0.16% into 28 "absent" segments.
    if explainable is None:
        explainable = parental_origin(father, sample).explainable
    if not math.isfinite(explainable):
        return [], {"reason": "no noise bound for this sample, so nothing can be called absent"}
    threshold = explainable * ABSENCE_MARGIN

    rates = []
    windows: list[tuple[str, int, int, float]] = []
    for chrom, v in per.items():
        for i in range(0, len(v) - window, window):
            chunk = v[i:i + window]
            r = sum(a for _, a in chunk) / len(chunk)
            rates.append(r)
            windows.append((chrom, chunk[0][0], chunk[-1][0], r))
    if len(rates) < 4:
        return [], {"reason": "too few markers to segment", "track": windows}
    loud = [r for r in rates if r > threshold]
    info = {"explainable": explainable, "threshold": threshold, "windows": len(rates),
            "windows_above": len(loud), "window_markers": window, "track": windows}
    if not loud:
        info["reason"] = (f"no window exceeds {threshold:.2%}, the level this sample's own noise "
                          "can reach, so the paternal genome is present throughout")
        return [], info
    if len(loud) >= 0.9 * len(rates):
        info["reason"] = (f"{len(loud)} of {len(rates)} windows exceed {threshold:.2%}, so the "
                          "paternal genome is absent throughout and there is no boundary to "
                          "find. Segmental loss is only meaningful against a genome that has "
                          "the rest of it.")
        return [], info

    p_lo = max(explainable / 2.0, 5e-4)
    p_hi = max(st.median(loud), p_lo * 5.0)
    info["p_present"], info["p_absent"] = p_lo, p_hi

    out: list[OriginSegment] = []
    for chrom, v in sorted(per.items(), key=lambda kv: (len(kv[0]), kv[0])):
        if len(v) < min_markers:
            continue
        path = _viterbi_two_state(v, p_lo, p_hi, expected_segment_mb * 1e6)
        i = 0
        while i < len(path):
            j = i
            while j + 1 < len(path) and path[j + 1] == path[i]:
                j += 1
            seg = v[i:j + 1]
            if len(seg) >= min_markers:
                a = sum(x for _, x in seg)
                rate = a / len(seg)
                out.append(OriginSegment(
                    chrom=chrom, start=seg[0][0], end=seg[-1][0], markers=len(seg),
                    absent=a, rate=rate,
                    state="paternal_absent" if path[i] else "paternal_present",
                    margin=rate / threshold if threshold > 0 else math.nan))
            i = j + 1
    return out, info


# --- one experiment: a sperm array and the plate of samples that came with it ------------------


@dataclass
class BiparentalOrigin:
    """What a sample inherited, when BOTH parents are on hand.

    Two direct measurements beat one measurement plus an inference. With only the father, whether
    a second parent contributed has to be read off the rate of alleles he lacks, an axis that
    separates by about 1.6x and needs the father's own heterozygosity as its reference. With the
    mother present each parent is tested the same way, on the axis that separates by thirty-fold,
    and the class falls out of the pair.
    """

    sample_id: str
    paternal: ParentageReport
    maternal: ParentageReport
    origin_class: Literal["biparental", "androgenetic", "gynogenetic",
                          "neither_parent", "unclear"]
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"{self.sample_id}: {self.origin_class.upper()}",
                 f"  paternal alleles absent {self.paternal.genome_rate:>7.2%}  "
                 f"(noise ceiling {self.paternal.explainable:.2%}) -> {self.paternal.verdict}",
                 f"  maternal alleles absent {self.maternal.genome_rate:>7.2%}  "
                 f"(noise ceiling {self.maternal.explainable:.2%}) -> {self.maternal.verdict}"]
        lines += [f"  {n}" for n in self.notes]
        return "\n".join(lines)


def both_parents(
    father: Sample, mother: Sample, sample: Sample, *, product: str = "unspecified",
) -> BiparentalOrigin:
    """Test each parent the same way and read the class off the pair."""
    pat = parental_origin(father, sample, product=product, role="paternal")
    mat = parental_origin(mother, sample, product=product, role="maternal")
    P = pat.verdict == "paternal_genome_present"
    M = mat.verdict == "paternal_genome_present"
    unresolved = "unclear" in (pat.verdict, mat.verdict)

    if unresolved:
        cls = "unclear"
    elif P and M:
        cls = "biparental"
    elif P:
        cls = "androgenetic"
    elif M:
        cls = "gynogenetic"
    else:
        cls = "neither_parent"

    out = BiparentalOrigin(sample.sample_id, pat, mat, cls)
    # Two declared parents that are the same genome is a labelling accident, not a pedigree, and
    # it produces a confident "biparental" because both tests pass against the same person.
    shared = [(father.probes[k].gt, mother.probes[k].gt) for k in father.probes
              if k in mother.probes and father.probes[k].gt != "NC"
              and mother.probes[k].gt != "NC"]
    if len(shared) > 1_000:
        same = sum(1 for a, b in shared if a == b) / len(shared)
        # 0.90, not 0.99, matching the ingestion layer's own measured figure: real replicate
        # arrays of one person concord at 95.8% and a parent-offspring pair at 54.9%, so 0.99 sat
        # above the thing it was meant to catch.
        if same > 0.90:
            out.notes.append(
                f"The two declared parents agree at {same:.1%} of shared markers, which is a "
                "duplicate or a relabelled file rather than two people. Every conclusion below "
                "rests on them being different, so resolve this before reading it."
            )
    if cls == "neither_parent":
        out.notes.append(
            "Neither declared parent accounts for this genome. Before reading that as biology, "
            "check the pairing: a mislabelled sample or the wrong donor pair produces exactly "
            "this, and so does contamination."
        )
    if cls == "androgenetic":
        out.notes.append("The maternal complement is absent and the mother's own array confirms "
                         "it directly, rather than it being inferred from allele sharing.")
    if cls == "unclear":
        out.notes.append("At least one parent is unresolved, so the pair cannot be classified. "
                         "The per-parent rates above say which one and by how much.")
    return out


@dataclass
class SampleVerdict:
    sample_id: str
    parentage: ParentageReport
    variant: Optional[EmbryoVerdict] = None
    #: Stretches lacking the paternal contribution, and those unique to this sample.
    segments: list["OriginSegment"] = field(default_factory=list)
    private_segments: list["OriginSegment"] = field(default_factory=list)
    #: Which paternal homologue, from co-inheritance across the siblings.
    relationship: Optional[str] = None
    homologue_odds: float = math.nan
    inherited: Optional[str] = None
    notes: list[str] = field(default_factory=list)


@dataclass
class ExperimentReport:
    father_id: str
    variant_chrom: Optional[str] = None
    variant_pos: Optional[int] = None
    build: Optional[str] = None
    assembly_note: str = ""
    samples: list[SampleVerdict] = field(default_factory=list)
    scaffold: Optional["scaffold.Scaffold"] = None
    refusals: list[str] = field(default_factory=list)

    def summary(self) -> str:
        head = f"experiment: sperm {self.father_id}, {len(self.samples)} samples"
        if self.variant_chrom:
            head += f", variant chr{self.variant_chrom}:{self.variant_pos:,}"
        lines = [head, ""]
        lines.append(f"  {'sample':<20s}{'origin':<15s}{'sperm':<11s}"
                     f"{'paternal locus':<18s}{'homologue':<14s}{'inherited'}")
        for s in self.samples:
            p = s.parentage
            locus = "-"
            if s.variant is not None:
                assessed = any(w.q_source != "uncalibrated_no_mother" for w in s.variant.windows)
                locus = ("ABSENT" if any(w.significant for w in s.variant.windows)
                         else "present" if assessed else "not tested")
            rel = s.relationship or "-"
            if math.isfinite(s.homologue_odds):
                rel += f" 10^{s.homologue_odds:+.0f}"
            lines.append(f"  {s.sample_id:<20s}{p.origin_class:<15s}{p.sperm_type:<11s}"
                         f"{locus:<18s}{rel:<18s}{s.inherited or '-'}")
        for s in self.samples:
            for n in s.notes:
                lines.append(f"  {s.sample_id}: {n}")
        if self.scaffold is not None:
            for r in self.scaffold.refusals:
                lines.append(f"  SCAFFOLD REFUSED: {r}")
        for r in self.refusals:
            lines.append(f"  REFUSED: {r}")
        return "\n".join(lines)


def run_experiment(
    father: Sample,
    samples: Sequence[Sample],
    *,
    mother: Optional[Sample] = None,
    variant_chrom: Optional[str] = None,
    variant_pos: Optional[int] = None,
    product: str = "unspecified",
    compression: Optional[float] = None,
    anchor: Optional[str] = None,
    anchor_carries: Optional[Literal["mutant", "wildtype"]] = None,
) -> ExperimentReport:
    """The whole plate at once: parent of origin, then which homologue, then the locus.

    The order is a gate, not a convenience. Which paternal homologue an embryo received is only
    a question if it received one, so a sample with no paternal genome is answered and set aside.
    And a paternal segment that is ABSENT at the locus makes the homologue call void there: you
    cannot have inherited a chromosome that is not present, and reporting a homologue anyway is
    exactly how a deletion gets read as a correction.
    """
    import scaffold as _scaffold

    rep = ExperimentReport(father_id=father.sample_id, variant_chrom=variant_chrom,
                           variant_pos=variant_pos)

    for s in samples:
        p = parental_origin(father, s, product=product)
        sv = SampleVerdict(sample_id=s.sample_id, parentage=p)
        if p.origin_class == "gynogenetic":
            sv.notes.append("no paternal genome, so no paternal homologue to identify")
        rep.samples.append(sv)

    paternal = [sv for sv in rep.samples
                if sv.parentage.origin_class in ("androgenetic", "biparental")]

    # The assembly is a property of the input, not of whether a variant was asked about. It was
    # previously detected only inside the scaffold branch, so a parentage-only run reported the
    # build as undetermined even where the coordinates say plainly which one it is.
    detected = detect_build(father)
    rep.build = detected.build
    rep.assembly_note = detected.note

    # Segmental loss, then the plate-level filter: a stretch recurring across independent samples
    # is a property of the probes there, not of any one genome.
    by_id_all = {s.sample_id: s for s in samples}
    for sv in paternal:
        segs, _info = segmental_origin(father, by_id_all[sv.sample_id],
                                       explainable=sv.parentage.explainable)
        sv.segments = [x for x in segs if x.state == "paternal_absent"]
    for sv in rep.samples:
        for seg in sv.segments:
            others = sum(1 for o in rep.samples if o is not sv and any(
                x.chrom == seg.chrom and x.start <= seg.end and seg.start <= x.end
                for x in o.segments))
            if others == 0:
                sv.private_segments.append(seg)
        strong = [x for x in sv.private_segments if x.confident]
        weak = len(sv.private_segments) - len(strong)
        for seg in strong:
            sv.notes.append(
                f"paternal contribution absent chr{seg.chrom}:{seg.start:,}-{seg.end:,} "
                f"({seg.span_mb:.1f} Mb, {seg.rate:.1%} absent, {seg.margin:.0f}x this sample's "
                f"noise floor), and in no other sample on the plate"
            )
        if weak:
            sv.notes.append(
                f"{weak} further stretch(es) sit within 3x the noise floor and are not called. "
                "Small subtelomeric and pericentromeric spans land there through correlated probe "
                "failure, which is why margin and plate-recurrence are both reported."
            )
    if variant_chrom is None:
        return rep

    by_id = {s.sample_id: s for s in samples}
    dropout = None
    if mother is not None:
        est = estimate_dropout({k: p.gt for k, p in father.probes.items()},
                               {k: p.gt for k, p in mother.probes.items()},
                               {k: p.gt for k, p in by_id[paternal[0].sample_id].probes.items()}
                               ) if paternal else None
        dropout = est[0] if est else None

    # --- the locus: is the paternal segment actually there -----------------------------------
    cal, _sex, _r = calibrate_sample(
        by_id[paternal[0].sample_id], product=product, declared_compression=compression
    ) if paternal else (nz.Calibration(0.0, math.nan, math.nan, "none"), None, None)
    for sv in paternal:
        s = by_id[sv.sample_id]
        c, _s2, _r2 = calibrate_sample(s, product=product, declared_compression=compression)
        ms = build_markers(father, s, mother, chrom=variant_chrom)
        if ms:
            sv.variant = analyse_embryo(sv.sample_id, ms, variant_chrom, variant_pos,
                                        calibration=c, dropout=dropout, phase=None)

    # --- which homologue, from the siblings ---------------------------------------------------
    if len(paternal) < 2:
        rep.refusals.append(
            "Fewer than two samples carry a paternal genome, so there is no co-inheritance to "
            "read. The paternal homologue stays unnamed."
        )
        return rep

    ids = [sv.sample_id for sv in paternal]
    fp = father.probes
    marks = []
    for k, f in fp.items():
        if f.gt != "AB" or f.chrom.replace("chr", "").lower() != variant_chrom.lower():
            continue
        gts = {i: by_id[i].probes[k].gt for i in ids if k in by_id[i].probes}
        if len(gts) >= 2:
            marks.append(_scaffold.Marker(k, f.chrom, f.pos, f.gt,
                                          mother.probes[k].gt if mother and k in mother.probes
                                          else None, gts))
    # Detect once, from the father, and hand it down: the scaffold is the only consumer whose
    # answer changes with the assembly, because it is the only one reading a genetic map.
    rep.scaffold = _scaffold.build_scaffold(
        marks, variant_chrom, variant_pos, anchor=anchor, anchor_carries=anchor_carries,
        dropout=dropout, build=detected.build)

    calls = {c.embryo_id: c for c in rep.scaffold.calls}
    for sv in paternal:
        c = calls.get(sv.sample_id)
        if c is None:
            continue
        sv.relationship = c.relationship
        sv.homologue_odds = c.log10_odds
        absent_here = sv.variant is not None and any(w.significant for w in sv.variant.windows)
        if absent_here:
            sv.relationship = "void"
            sv.homologue_odds = math.nan
            sv.inherited = None
            sv.notes.append(
                "the paternal segment is ABSENT at this locus, so which homologue was inherited "
                "is undefined here. An absent segment reads as a correction if the homologue is "
                "reported anyway, which is the error this ordering exists to prevent."
            )
        else:
            sv.inherited = c.inherited
    return rep


def _self_check() -> int:
    """Prove the whole chain runs standalone, with no data files and no network."""
    cal = nz.Calibration(center=0.0, lrr_compression=0.594, sigma_lrr=0.11, route="declared")
    absent = range(20, 41)
    ms = [Marker(f"m{i}", "11", 1_000_000 + i * 2_000, "AA", "BB",
                 "BB" if i in absent else "AB",
                 baf=1.0 if i in absent else 0.5,
                 lrr=-0.594 if i in absent else 0.0) for i in range(60)]
    v = analyse_embryo("selfcheck", ms, "11", 1_060_000, calibration=cal, dropout=0.01)
    w = v.windows[-1]
    assert w.n_l3 == 59, w.n_l3
    assert w.longest_run == 20, w.longest_run
    assert w.segment_state == "PAT0_MAT1", w.segment_state
    assert w.run_p < 1e-20, w.run_p
    assert any("phase source" in r for r in v.refusals), v.refusals
    assert run_length_p(10, 10, KOTHIYAL_FLOOR) > 0, "the run form must not underflow to zero"
    assert score_paternal("AA", "AB", "AB") is None, "presence needs a mother who cannot supply"
    assert estimate_dropout({"r": "AB"}, {"r": "AA"}, {"r": "AB"})[0] == 0.0
    print("origin self-check OK")
    return 0


def _main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    if argv is None and len(sys.argv) == 1:
        return _self_check()

    ap = argparse.ArgumentParser(
        prog="origin",
        description="Which paternal allele did this embryo inherit, and is a correction real. "
                    "Research use only; not a clinical diagnostic.",
    )
    ap.add_argument("--father", required=True, help="paternal array export")
    ap.add_argument("--embryo", help="one embryo or biopsy array export")
    ap.add_argument("--samples", nargs="+", metavar="FILE",
                    help="the whole plate at once: every sample that came with this sperm array")
    ap.add_argument("--anchor", help="sample id whose variant-site genotype is known externally")
    ap.add_argument("--anchor-carries", choices=("mutant", "wildtype"),
                    help="what that sample was typed as. Names both homologue groups at once.")
    ap.add_argument("--mother", help="maternal array export. Omit to run the degraded mode.")
    # Optional: with neither, the tool answers the question a lab asks first, which is whether
    # this sample carries a paternal genome at all. That needs no variant and no mother.
    ap.add_argument("--chrom", help="variant chromosome, e.g. 6. Omit to scan parental origin.")
    ap.add_argument("--pos", type=int, help="variant position, in the array's build")
    # No platform default. The chrX baseline offset is a MEASURED per-array constant, and
    # defaulting to a named array silently asserts one array's constant for another's data. On
    # the PennCNV Illumina trio that mistake put c at 0.456 against a measured 0.59.
    ap.add_argument("--product", default="unspecified",
                    help="array product, which selects the measured chrX baseline offset. "
                         "Without one the chrX route is refused rather than approximated.")
    ap.add_argument("--compression", type=float,
                    help="externally measured LRR compression, e.g. from a known one-copy "
                         "segment. Used only if no internal reference is available, and flagged "
                         "as unverified in the output.")
    ap.add_argument("--mutant-allele", choices=("A", "B"),
                    help="paternal allele carrying the variant. Without it H1 vs H2 is refused.")
    ap.add_argument("--phase-route", help="how --mutant-allele was established, recorded verbatim")
    args = ap.parse_args(argv)

    if bool(args.mutant_allele) != bool(args.phase_route):
        ap.error("--mutant-allele and --phase-route must be given together: an unattributed "
                 "phase claim is the failure mode this tool exists to prevent")

    if bool(args.chrom) != bool(args.pos is not None):
        ap.error("--chrom and --pos go together; give both to analyse a variant, or neither to "
                 "scan parental origin genome-wide")
    if bool(args.embryo) == bool(args.samples):
        ap.error("give either --embryo for one sample or --samples for the whole plate")
    if bool(args.anchor) != bool(args.anchor_carries):
        ap.error("--anchor and --anchor-carries go together: an anchor without its typed status "
                 "names nothing, and a status without a sample belongs to nobody")

    if args.samples:
        father = read_sample(args.father)
        mother = read_sample(args.mother) if args.mother else None
        plate = [read_sample(p) for p in args.samples]
        rep = run_experiment(
            father, plate, mother=mother, variant_chrom=args.chrom, variant_pos=args.pos,
            product=args.product, compression=args.compression,
            anchor=args.anchor, anchor_carries=args.anchor_carries)
        print(rep.summary())
        return 0

    father = read_sample(args.father)
    embryo = read_sample(args.embryo)
    mother = read_sample(args.mother) if args.mother else None

    if args.chrom is None:
        print(parental_origin(father, embryo, product=args.product).summary())
        return 0

    cal, sex, ratio = calibrate_sample(embryo, product=args.product,
                                       declared_compression=args.compression)
    print(f"embryo {embryo.sample_id}: chrX/autosomal heterozygosity {ratio:.3f} -> sex {sex}")
    print(f"calibration: c={cal.lrr_compression:.3f} sigma={cal.sigma_lrr:.3f} via {cal.route}")
    for w in cal.warnings:
        print(f"  warning: {w}")

    dropout = None
    if mother:
        est = estimate_dropout(
            {k: p.gt for k, p in father.probes.items()},
            {k: p.gt for k, p in mother.probes.items()},
            {k: p.gt for k, p in embryo.probes.items()},
        )
        if est:
            dropout, se, n = est
            print(f"fitted dropout: {dropout:.3f} +/- {se:.4f} from n={n:,}")
        else:
            print("dropout: no father-het x mother-hom markers, cannot fit")
    else:
        print("dropout: no mother supplied, falling back to the Kothiyal floor")

    markers = build_markers(father, embryo, mother, chrom=args.chrom)
    print(f"chr{args.chrom} markers: {len(markers):,}")
    if not markers:
        print("REFUSED: no shared markers on that chromosome")
        return 1

    phase = (PhaseSource(mutant_allele=args.mutant_allele, route=args.phase_route)
             if args.mutant_allele else None)
    verdict = analyse_embryo(embryo.sample_id, markers, args.chrom, args.pos,
                             calibration=cal, dropout=dropout, phase=phase)
    print()
    print(verdict.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())

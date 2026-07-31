"""emissions - what an array marker looks like under each physical state of the locus.

Phase 3 of OriginMarker 2.0 (build spec 6, 7.1-7.2). The layer that finally separates
copy-neutral from copy-loss LOH, which is the distinction Ma's reply could not make: gene
conversion and paternal segmental deletion predict identical loss of paternal SNPs and differ
ONLY in copy number.

This module is emissions and per-sample fitting. The HMM, segmentation and posterior sit on top
of it and are separate.

Three corrections to spec 7.2 are implemented here rather than the text as written. Each is
checked against `array_signatures.csv`, the spec's own companion table, which has them right.

1. TWO channels, not three. The spec multiplies P(GType) * P(BAF) * P(LRR) and calls them
   "distinct measurements". They are not. An Infinium array measures two numbers per marker, a
   normalised allelic ratio theta and a total intensity R. BAF is an interpolation in theta; LRR
   is log2(R_obs / R_exp) with R_exp interpolated from the same theta; GType is the cluster call
   from that same pair. So (BAF, LRR) is a reparameterisation of the two real measurements and
   factorising it is an ordinary conditional-independence approximation, while multiplying GType
   in as a third factor counts theta twice. GType enters ONLY where BAF is missing, as a coarser
   readout of the same allelic quantity, never alongside it.

   Which channel is the fallback is settled by the data, not by preference: in real Axiom output
   BAF is present at 95.6% of markers and at 87% of NO-CALL markers, against a genotype present
   at 90.7%. BAF is the more complete channel, so it leads.

2. MOSAIC BAF is a ratio of expectations, not a linear interpolation of BAF. The spec writes
   `f * mu_b(S) + (1-f) * mu_b(normal)`. At mosaic loss with f = 0.5 that gives 0.75 where the
   truth is 0.667, an error of 0.083 - nearly three standard deviations at a typical
   sigma_BAF of 0.03, applied systematically to every mosaic state.

3. MOSAIC LRR likewise. The spec writes `f * mu_LRR(S) + (1-f) * 0`, linear in LRR. At CN1 and
   f = 0.5 that gives -0.500 against a true -0.415.

Both of those are the same mistake: the observation is a nonlinear function of copy number, so
the mixture has to be taken in copy number and pushed through the nonlinearity afterwards, not
the other way round.

Self-check:  python -m pytest tests/test_emissions.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, NamedTuple, Optional, Sequence

# Kothiyal's floor lives in the run-length layer; nothing here needs it.

# --- the state space (spec 7.1) ------------------------------------------------------------


@dataclass(frozen=True)
class CopyState:
    """A local copy/structural configuration, counted per parent.

    Counting per PARENT rather than as a total is what makes the two isodisomies distinct
    states. PAT2_MAT0 and PAT0_MAT2 are both CN2 with flat LRR and both show loss of
    heterozygosity, and they are opposite findings: one is paternal isodisomy, the other
    maternal. A total-copy-number state space cannot tell them apart.
    """

    name: str
    n_pat: int
    n_mat: int
    #: An insertion at the cut site leaves copy number nominal, so it is flagged rather than
    #: encoded in the counts. Array data cannot see it at all; this state exists so the
    #: mechanism ledger can say "not tested" instead of "excluded".
    insertion: bool = False

    @property
    def cn(self) -> int:
        return self.n_pat + self.n_mat


STATES: tuple[CopyState, ...] = (
    CopyState("PAT1_MAT1", 1, 1),
    CopyState("PAT0_MAT1", 0, 1),
    CopyState("PAT1_MAT0", 1, 0),
    CopyState("PAT0_MAT2", 0, 2),          # maternal isodisomy, CN2, copy-neutral
    CopyState("PAT2_MAT0", 2, 0),          # PATERNAL isodisomy, CN2, copy-neutral
    CopyState("PAT2_MAT1", 2, 1),
    CopyState("PAT1_MAT2", 1, 2),
    CopyState("PAT2_MAT2", 2, 2),
    CopyState("PAT0_MAT0", 0, 0),
    CopyState("PAT1_MAT1_INS", 1, 1, insertion=True),
)

BY_NAME = {s.name: s for s in STATES}

NORMAL = BY_NAME["PAT1_MAT1"]


# --- allelic composition ------------------------------------------------------------------

#: A genotype in AB space. NC is a no-call and is never treated as a homozygote.
Genotype = str  # 'AA' | 'AB' | 'BB' | 'NC'


class Band(NamedTuple):
    """One BAF band the state predicts at this marker, with its prior weight."""

    baf: float
    lrr: float
    weight: float
    #: B-allele copies and total copies BEFORE mosaic mixing, kept for auditing.
    n_b: int
    cn: int


def _transmissions(
    father: Genotype, mother: Optional[Genotype], paternal_allele: Optional[str]
) -> list[tuple[str, str, float]]:
    """Every (paternal allele, maternal allele, weight) this marker could have transmitted.

    A homozygous parent transmits a known allele. A heterozygous one is marginalised over, at
    half weight each, because which allele he or she passed on is simply not observed. That
    marginalisation is the honest treatment and it is also why a both-parents-heterozygous
    marker carries so little information: every state predicts the same mixture there.

    `paternal_allele` is the allele the declared phase says came through. Supplying it collapses
    the paternal marginalisation; withholding it keeps both branches, which is what happens when
    no phase source was declared.
    """
    if father == "NC":
        return []
    pat: list[tuple[str, float]]
    if father in ("AA", "BB"):
        pat = [(father[0], 1.0)]
    elif paternal_allele in ("A", "B"):
        pat = [(paternal_allele, 1.0)]
    else:
        pat = [("A", 0.5), ("B", 0.5)]

    mat: list[tuple[str, float]]
    if mother is None or mother == "NC":
        # No maternal genotype: her contribution is unconstrained, so both alleles stay live.
        mat = [("A", 0.5), ("B", 0.5)]
    elif mother in ("AA", "BB"):
        mat = [(mother[0], 1.0)]
    else:
        mat = [("A", 0.5), ("B", 0.5)]

    return [(p, m, wp * wm) for p, wp in pat for m, wm in mat]


def bands(
    state: CopyState,
    father: Genotype,
    mother: Optional[Genotype],
    *,
    paternal_allele: Optional[str] = None,
    mosaic_fraction: float = 1.0,
    lrr_compression: float = 1.0,
) -> list[Band]:
    """The BAF bands and LRR mean this state predicts at this marker.

    `mosaic_fraction` is the fraction of cells carrying the state; the rest are normal diploid.
    f = 1.0 is a clonal event.

    The mosaic mixing is a ratio of EXPECTATIONS, per `array_signatures.csv`:

        BAF = E[B copies] / E[total copies]
        LRR = c * log2( E[total copies] / 2 )

    and NOT a linear interpolation of BAF or of LRR, which is what spec 7.2 writes. See the
    module docstring for the size of that error.
    """
    if not 0.0 <= mosaic_fraction <= 1.0:
        raise ValueError(f"mosaic fraction out of range: {mosaic_fraction}")

    f = mosaic_fraction
    e_cn = f * state.cn + (1.0 - f) * NORMAL.cn
    if e_cn <= 0:
        # A homozygous deletion at f = 1 has no signal to interpolate theta from, so BAF is
        # noise rather than a band. Reported as no bands; the caller must fall back to the
        # intensity channel alone.
        return []

    lrr = lrr_compression * math.log2(e_cn / NORMAL.cn)

    out: list[Band] = []
    for p, m, w in _transmissions(father, mother, paternal_allele):
        n_b = state.n_pat * (p == "B") + state.n_mat * (m == "B")
        n_b_normal = (p == "B") + (m == "B")
        e_b = f * n_b + (1.0 - f) * n_b_normal
        out.append(Band(baf=e_b / e_cn, lrr=lrr, weight=w, n_b=n_b, cn=state.cn))
    return out


# --- the two channels ---------------------------------------------------------------------

_SQRT2 = math.sqrt(2.0)
_LOG_SQRT_2PI = 0.5 * math.log(2.0 * math.pi)


def _norm_logpdf(x: float, mu: float, sigma: float) -> float:
    z = (x - mu) / sigma
    return -0.5 * z * z - math.log(sigma) - _LOG_SQRT_2PI


def _norm_cdf(x: float, mu: float, sigma: float) -> float:
    return 0.5 * (1.0 + math.erf((x - mu) / (sigma * _SQRT2)))


def _truncated_norm_logpdf(x: float, mu: float, sigma: float, lo: float, hi: float) -> float:
    """Normal density truncated to [lo, hi].

    Truncation is not cosmetic here and omitting it biases the whole model. BAF is bounded to
    [0, 1] and the homozygous bands sit ON the boundaries, so an untruncated normal spends about
    half its mass outside the support at band 0 and band 1 while the heterozygous band at 0.5
    spends almost none. The result is a systematic factor-of-two penalty against every
    homozygous state, which is to say against exactly the states that indicate allele loss.
    """
    if not lo <= x <= hi:
        return -math.inf
    mass = _norm_cdf(hi, mu, sigma) - _norm_cdf(lo, mu, sigma)
    if mass <= 0.0:
        return -math.inf
    return _norm_logpdf(x, mu, sigma) - math.log(mass)


def _genotype_from_baf(baf: float) -> Genotype:
    """The discrete call an array would make from a band. A coarsening, by design."""
    if baf < 0.25:
        return "AA"
    if baf > 0.75:
        return "BB"
    return "AB"


def genotype_logp(observed: Genotype, expected: Genotype, eps: float) -> float:
    """Discrete allelic channel, used ONLY when BAF is missing.

    The dropout kernel is one-directional: a heterozygote can be read as either homozygote, and
    a homozygote cannot be read as a heterozygote. That asymmetry is what makes a het call
    robust to every H3 mechanism and a hom call robust to none.

    Note the asymmetry is an approximation rather than a theorem on amplified material -
    erroneous heterozygous calls become common below 60% call rate - so the caller must gate
    this on call rate rather than rely on it universally.
    """
    if observed == "NC" or expected == "NC":
        return 0.0  # no information, contributes nothing
    if observed == expected:
        return math.log1p(-eps) if eps < 1.0 else -math.inf
    if expected == "AB" and observed in ("AA", "BB"):
        return math.log(eps / 2.0) if eps > 0.0 else -math.inf
    # A homozygote read as the OTHER homozygote, or as a heterozygote, needs two independent
    # errors. Kept at eps^2 rather than zero so a single bad marker cannot veto a state outright.
    return 2.0 * math.log(eps) if eps > 0.0 else -math.inf


@dataclass(frozen=True)
class EmissionParams:
    """Per-sample fitted parameters. None of these is a literature constant (spec 7.4)."""

    #: LRR compression. Fitted, per sample: real amplified material runs 0.45-0.58 against the
    #: 0.70-0.85 spec 7.2 expected from Illumina vendor figures, and PennCNV's own shipped
    #: defaults sit at 0.57-0.68. There is no single right value to hardcode.
    lrr_compression: float
    sigma_baf: float
    sigma_lrr: float
    #: Per-marker genotype error/dropout rate, used only by the discrete fallback.
    eps: float


def emission_logp(
    state: CopyState,
    father: Genotype,
    mother: Optional[Genotype],
    *,
    baf: Optional[float],
    lrr: Optional[float],
    genotype: Genotype = "NC",
    params: EmissionParams,
    paternal_allele: Optional[str] = None,
    mosaic_fraction: float = 1.0,
) -> float:
    """log P(observation | state) for one marker, over TWO channels.

    Allelic channel: BAF when it is present, the discrete genotype only when it is not. Never
    both, because they are two readings of the same theta and multiplying them counts it twice.

    Intensity channel: LRR, which carries R and so is a genuinely separate measurement.

    A channel that is absent contributes nothing rather than contributing a penalty, so a marker
    with no LRR is simply less informative rather than evidence against every state.
    """
    predicted = bands(
        state,
        father,
        mother,
        paternal_allele=paternal_allele,
        mosaic_fraction=mosaic_fraction,
        lrr_compression=params.lrr_compression,
    )
    if not predicted:
        # No signal to predict: a homozygous deletion. Only the intensity channel can speak.
        if lrr is None:
            return 0.0
        mu = params.lrr_compression * math.log2(max(1e-3, 0.0) / NORMAL.cn)
        return _norm_logpdf(lrr, mu, params.sigma_lrr)

    total = 0.0

    # --- allelic ---
    if baf is not None:
        terms = [
            math.log(b.weight) + _truncated_norm_logpdf(baf, b.baf, params.sigma_baf, 0.0, 1.0)
            for b in predicted
            if b.weight > 0.0
        ]
        total += _logsumexp(terms)
    elif genotype != "NC":
        terms = [
            math.log(b.weight) + genotype_logp(genotype, _genotype_from_baf(b.baf), params.eps)
            for b in predicted
            if b.weight > 0.0
        ]
        total += _logsumexp(terms)

    # --- intensity ---
    if lrr is not None:
        total += _norm_logpdf(lrr, predicted[0].lrr, params.sigma_lrr)

    return total


def _logsumexp(terms: Sequence[float]) -> float:
    finite = [t for t in terms if t > -math.inf]
    if not finite:
        return -math.inf
    hi = max(finite)
    return hi + math.log(sum(math.exp(t - hi) for t in finite))


# --- per-sample fitting -------------------------------------------------------------------


class CompressionFit(NamedTuple):
    """Fitted LRR compression, pooled and per copy-number state."""

    pooled: float
    per_cn: dict[int, float]
    n: dict[int, int]
    note: str


def fit_lrr_compression(observations: Iterable[tuple[int, float]]) -> CompressionFit:
    """Fit c from (copy_number, lrr) pairs.

    c is the ratio of observed LRR to its theoretical value log2(CN/2). It is fitted rather than
    assumed because the observed value is not a property of the method: it depends on the array,
    the chemistry and the input mass, and on real amplified material it runs roughly half the
    theoretical figure.

    Reported per copy-number state as well as pooled, because the compression is NOT constant
    across states - PennCNV's shipped defaults already encode per-state means rather than one
    scale factor. A large spread between states is a finding about the sample, not noise to
    average away.

    CN2 contributes nothing: log2(2/2) is zero, so it constrains no scale. CN0 is excluded
    because log2(0) diverges.
    """
    sums: dict[int, list[float]] = {}
    for cn, lrr in observations:
        if cn <= 0 or cn == NORMAL.cn or not math.isfinite(lrr):
            continue
        sums.setdefault(cn, []).append(lrr)

    per_cn: dict[int, float] = {}
    counts: dict[int, int] = {}
    num = den = 0.0
    for cn, vals in sums.items():
        x = math.log2(cn / NORMAL.cn)
        mean = sum(vals) / len(vals)
        per_cn[cn] = mean / x
        counts[cn] = len(vals)
        # Least squares through the origin, weighted by how many markers each state contributes.
        num += len(vals) * x * mean
        den += len(vals) * x * x

    if den == 0.0:
        return CompressionFit(
            pooled=math.nan,
            per_cn={},
            n={},
            note="no markers outside CN2, so no scale is constrained and c cannot be fitted "
            "from this input. It must come from the HMM fit or be declared.",
        )

    pooled = num / den
    spread = (max(per_cn.values()) - min(per_cn.values())) if len(per_cn) > 1 else 0.0
    note = (
        f"pooled c = {pooled:.3f} from {sum(counts.values())} markers across CN "
        f"{sorted(per_cn)}; per-state spread {spread:.3f}"
    )
    if spread > 0.10:
        note += (
            ". That spread is large enough that a single pooled c misfits some states: prefer "
            "the per-state values"
        )
    return CompressionFit(pooled=pooled, per_cn=per_cn, n=counts, note=note)


def fit_sigma(residuals: Iterable[float]) -> float:
    """Robust scale from residuals, via the median absolute deviation.

    MAD rather than the standard deviation because the residuals are contaminated by exactly the
    events being looked for: a real segmental loss puts a block of markers far from the normal
    mean, and a plain standard deviation would absorb them into the noise estimate and then be
    unable to see them.
    """
    vals = [r for r in residuals if math.isfinite(r)]
    if len(vals) < 2:
        return math.nan
    vals.sort()
    mid = len(vals) // 2
    med = vals[mid] if len(vals) % 2 else 0.5 * (vals[mid - 1] + vals[mid])
    dev = sorted(abs(v - med) for v in vals)
    mad = dev[mid] if len(dev) % 2 else 0.5 * (dev[mid - 1] + dev[mid])
    return 1.4826 * mad  # scaled so it estimates sigma for a normal


# --- the distinction the whole layer exists for -------------------------------------------


def separates_copy_neutral(
    father: Genotype,
    mother: Optional[Genotype],
    *,
    lrr_compression: float,
    sigma_lrr: float,
    segment_markers: int = 1,
) -> dict[str, float]:
    """How far apart, in LRR standard deviations, are the copy-loss and copy-neutral states.

    This is the number that decides whether the layer can do its job on a given sample.
    PAT0_MAT1 (paternal loss, CN1) and PAT0_MAT2 (maternal isodisomy, CN2) predict the SAME
    allelic pattern - paternal alleles gone - and differ only in intensity. Gene conversion and
    paternal segmental deletion are the biological versions of that pair, and Ma's reply rested
    on LOH without copy number, which is why it could not choose between them.

    `segment_markers` is load-bearing and defaulting it to 1 is a trap I fell into: the states
    differ by c in EXPECTED LRR, and the uncertainty on a segment MEAN falls as sigma/sqrt(n), so
    the per-marker figure understates the available power by sqrt(n). On real data c = 0.596
    against sigma = 0.417 is 1.4 sigma per marker - which I wrongly reported as "this sample
    cannot make the call" - and 14 sigma over a hundred markers. Always pass the length of the
    segment actually being called.

    Returns the separation in sigma units. Below about 2 the segment cannot support the call, and
    the mechanism ledger must report copy-neutral LOH as NOT TESTED rather than excluded.
    """
    loss = bands(BY_NAME["PAT0_MAT1"], father, mother, lrr_compression=lrr_compression)
    neutral = bands(BY_NAME["PAT0_MAT2"], father, mother, lrr_compression=lrr_compression)
    if not loss or not neutral:
        return {"delta_lrr": math.nan, "sigmas": math.nan}
    delta = abs(loss[0].lrr - neutral[0].lrr)
    n = max(1, int(segment_markers))
    effective = sigma_lrr / math.sqrt(n)
    return {
        "delta_lrr": delta,
        "sigmas": delta / effective if sigma_lrr > 0 else math.nan,
        "sigmas_per_marker": delta / sigma_lrr if sigma_lrr > 0 else math.nan,
        "segment_markers": float(n),
    }

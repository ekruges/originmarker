"""hmm - one joint inference over the latent physical state of the locus.

Phase 3 of OriginMarker 2.0 (build spec 7.5-7.6, 6.5). Sits on `emissions`.

The structural point, and the reason this replaced an earlier design: the latent state is a PAIR,
(T, C). T is which paternal homologue was transmitted and changes only by meiotic crossover. C is
the local copy/structural configuration and changes only by a somatic breakpoint. Keeping them
separate is what lets the model tell a crossover (T flips, C stays diploid) from a deletion
(T unchanged, C drops a paternal copy) - the confusion a smoothing step cannot resolve, because
on the observations alone both look like "the markers changed".

Everything is one joint likelihood over markers. An earlier version multiplied per-layer
likelihoods aggregated over many markers, which counted a single segmental loss three times as
independent evidence, in the direction of overconfidence.

Deliberately stdlib-only, like `panelbuilder` and `genetic_map`. numpy is present in the dev
environment but is NOT in requirements.txt, so importing it here would build a container that
dies on startup. A window holds thousands of markers, not millions, and 20 states over that is
milliseconds in plain Python.

Self-check:  python -m pytest tests/test_hmm.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Literal, NamedTuple, Optional, Sequence

import emissions as em
import genetic_map

# --- the state space -----------------------------------------------------------------------

#: The state space is the copy states, and nothing else. Which paternal homologue was
#: transmitted is `scaffold`'s question, answered from co-inheritance across siblings; a
#: per-marker A/B letter convention cannot carry it, because which allele is labelled A is
#: arbitrary until the father is phased.
STATES: tuple[em.CopyState, ...] = em.STATES


# --- transitions (spec 7.5) ---------------------------------------------------------------


@dataclass(frozen=True)
class TransitionParams:
    """Per-window transition structure, with ASYMMETRIC entry and exit rates.

    Spec 7.5 gives a single hazard, `P(C changes) = 1 - exp(-L / L_event)`, applied in both
    directions. That cannot be right: with one rate the expected event length equals the expected
    normal length, so the prior asserts that about half the genome is aberrant. Real structural
    events are rare and bounded, which needs two rates and is what CNV HMMs use.

      entry, per bp:  events_per_mb / 1e6      how often an event STARTS
      exit,  per bp:  1 / event_length_bp      how long one LASTS once started

    Both are user-set rather than fitted. No calibration exists for edited human embryos, and the
    spec's three windows genuinely expect different things - kilobases locally, megabases
    segmentally - so a single fitted value would be wrong for two of the three. [N]
    """

    #: Expected extent of an event once it has started. Sets the exit rate.
    event_length_bp: float
    #: Expected number of independent events per megabase. Sets the entry rate. The default is
    #: deliberately low: an event is an excursion, not the background.
    events_per_mb: float = 0.05


class MapUse(NamedTuple):
    """Which genetic map actually supplied the crossover probability, and whether it is the
    one the spec asks for."""

    cm: float
    theta: float
    source: str
    approx: bool
    #: True whenever the map is not male-specific. Male recombination is strongly
    #: telomere-biased, and no licence-clear male map was found, so this fallback is labelled
    #: and travels with every affected output.
    sex_averaged_fallback: bool


def crossover(chrom: str, pos_a: int, pos_b: int) -> MapUse:
    """P(T flips) between two adjacent markers, from the bundled map via Haldane.

    Reuses `genetic_map` rather than reimplementing interpolation. That map is SEX-AVERAGED for
    the autosomes, and 2.0 wants the male map, so the result is flagged. The direction of the
    error is worth knowing: male recombination is telomere-biased, so a sex-averaged rate
    OVERSTATES crossover risk near the centromere and UNDERSTATES it near the telomeres. For a
    near-centromeric locus the fallback is therefore conservative - it makes haplotype calls less
    confident, not more - and for a subtelomeric one it is anti-conservative.
    """
    gm = genetic_map.load(chrom)
    cm = abs(gm.cm_between(min(pos_a, pos_b), max(pos_a, pos_b)))
    covered = gm.covers(pos_a) and gm.covers(pos_b)
    return MapUse(
        cm=cm,
        theta=genetic_map.haldane_theta(cm),
        source=gm.source,
        approx=gm.approx or not covered,
        sex_averaged_fallback=True,
    )


class Hazards(NamedTuple):
    """Entry and exit probabilities across one inter-marker gap."""

    enter: float
    exit: float


def breakpoint_hazard(gap_bp: int, params: TransitionParams) -> Hazards:
    """P(an event starts) and P(a running event ends) across an inter-marker gap.

    Two rates, not one: see TransitionParams. Both saturate at 1 across a gap much larger than
    their characteristic length, which is correct for a Poisson process and is why a very wide
    gap carries essentially no information about whether the state persisted across it.
    """
    if gap_bp <= 0:
        return Hazards(0.0, 0.0)
    enter = (-math.expm1(-gap_bp * params.events_per_mb / 1e6)
             if params.events_per_mb > 0 else 0.0)
    exit_ = (-math.expm1(-gap_bp / params.event_length_bp)
             if params.event_length_bp > 0 else 0.0)
    return Hazards(enter=enter, exit=exit_)


NORMAL_C = em.NORMAL


def transition_logp(
    frm: em.CopyState, to: em.CopyState, *, hazards: Hazards, n_aberrant: int
) -> float:
    """log P(to | frm), factorised as P(T'|T) * P(C'|C).

    The factorisation is an assumption and a defensible one: a meiotic crossover in the father's
    germline and a somatic breakpoint in the embryo are physically independent processes. It is
    NOT the discarded layer-independence assumption, which multiplied three readings of the same
    markers.

    Structurally, events are excursions FROM normal and back. An aberrant state transitions only
    to normal, never straight to a different aberrant state. Two adjacent distinct events are
    therefore represented as event-normal-event, which costs one marker of normal between them -
    a real limitation, and the alternative is inventing a transition prior over event pairs that
    nothing in the literature constrains.
    """
    from_normal = frm == NORMAL_C
    to_normal = to == NORMAL_C

    if from_normal:
        p_c = (1.0 - hazards.enter) if to_normal else hazards.enter / max(1, n_aberrant)
    elif to == frm:
        p_c = 1.0 - hazards.exit
    elif to_normal:
        p_c = hazards.exit
    else:
        p_c = 0.0                      # aberrant to a different aberrant: via normal only

    if p_c <= 0.0:
        return -math.inf
    return math.log(p_c)


# --- observations --------------------------------------------------------------------------


@dataclass(frozen=True)
class MarkerObs:
    """One marker's parental context and observation."""

    rsid: str
    pos: int
    father: str
    mother: Optional[str]
    baf: Optional[float] = None
    lrr: Optional[float] = None
    genotype: str = "NC"


def emission_row(
    obs: MarkerObs, params: em.EmissionParams, *, mosaic_fraction: float = 1.0
) -> list[float]:
    """log P(obs | S) for every state, in STATES order.

    A homozygous father transmits a known allele and it is used. A heterozygous one is
    marginalised over: which of his homologues came through is not a question this layer can
    answer, and the convention that used to answer it here tracked an arbitrary A/B labelling
    rather than a chromosome.
    """
    out: list[float] = []
    for st in STATES:
        # A homozygous father transmits a known allele; a heterozygous one is marginalised over,
        # since which homologue came through is `scaffold`'s question.
        pat_allele = obs.father[0] if obs.father in ("AA", "BB") else None
        out.append(
            em.emission_logp(
                st,
                obs.father,
                obs.mother,
                baf=obs.baf,
                lrr=obs.lrr,
                genotype=obs.genotype,
                params=params,
                paternal_allele=pat_allele,
                mosaic_fraction=mosaic_fraction,
            )
        )
    return out


# --- forward-backward and Viterbi ---------------------------------------------------------


def _logsumexp(xs: Iterable[float]) -> float:
    vals = [x for x in xs if x > -math.inf]
    if not vals:
        return -math.inf
    hi = max(vals)
    return hi + math.log(sum(math.exp(v - hi) for v in vals))


class Inference(NamedTuple):
    #: Per-marker posterior over states, each row summing to 1.
    gamma: list[list[float]]
    #: Most likely single path, as indices into STATES.
    viterbi: list[int]
    log_likelihood: float
    #: The map fallback used at each step, so the label reaches the output.
    map_uses: list[MapUse]


def infer(
    markers: Sequence[MarkerObs],
    chrom: str,
    params: em.EmissionParams,
    tp: TransitionParams,
    *,
    mosaic_fraction: float = 1.0,
) -> Inference:
    """Forward-backward and Viterbi over the joint state, in log space.

    Log space throughout is not defensive: with a few thousand markers the forward variable
    underflows float64 within a few hundred steps, and a scaled-linear implementation would
    silently return zeros for exactly the long windows that matter most.
    """
    n = len(markers)
    k = len(STATES)
    if n == 0:
        return Inference(gamma=[], viterbi=[], log_likelihood=0.0, map_uses=[])

    n_aberrant = len(em.STATES) - 1
    logB = [emission_row(m, params, mosaic_fraction=mosaic_fraction) for m in markers]

    # Flat prior: no state is privileged before the data speaks. Anchoring on normal diploid
    # would bias the first marker's posterior toward "nothing happened".
    log_pi = -math.log(k)

    uses: list[MapUse] = []
    trans: list[list[list[float]]] = []
    for i in range(1, n):
        mu = crossover(chrom, markers[i - 1].pos, markers[i].pos)
        hz = breakpoint_hazard(markers[i].pos - markers[i - 1].pos, tp)
        uses.append(mu)
        trans.append([
            [transition_logp(a, b, hazards=hz, n_aberrant=n_aberrant)
             for b in STATES]
            for a in STATES
        ])

    # forward
    alpha = [[log_pi + logB[0][j] for j in range(k)]]
    for i in range(1, n):
        prev, A = alpha[-1], trans[i - 1]
        alpha.append([
            _logsumexp([prev[a] + A[a][j] for a in range(k)]) + logB[i][j] for j in range(k)
        ])
    loglik = _logsumexp(alpha[-1])

    # backward
    beta = [[0.0] * k for _ in range(n)]
    for i in range(n - 2, -1, -1):
        A, nxt, nb = trans[i], logB[i + 1], beta[i + 1]
        beta[i] = [
            _logsumexp([A[a][j] + nxt[j] + nb[j] for j in range(k)]) for a in range(k)
        ]

    gamma: list[list[float]] = []
    for i in range(n):
        row = [alpha[i][j] + beta[i][j] for j in range(k)]
        z = _logsumexp(row)
        gamma.append([math.exp(v - z) if z > -math.inf else 1.0 / k for v in row])

    # Viterbi
    delta = [log_pi + logB[0][j] for j in range(k)]
    back: list[list[int]] = []
    for i in range(1, n):
        A = trans[i - 1]
        nd, bp = [0.0] * k, [0] * k
        for j in range(k):
            best, arg = -math.inf, 0
            for a in range(k):
                v = delta[a] + A[a][j]
                if v > best:
                    best, arg = v, a
            nd[j], bp[j] = best + logB[i][j], arg
        delta, _ = nd, None
        back.append(bp)
    path = [max(range(k), key=lambda j: delta[j])]
    for bp in reversed(back):
        path.append(bp[path[-1]])
    path.reverse()

    return Inference(gamma=gamma, viterbi=path, log_likelihood=loglik, map_uses=uses)


# --- segmentation, the continuous second opinion (spec 7.6 step 4) -------------------------


def binary_segment(
    values: Sequence[float], sigma: float, *, min_size: int = 3, penalty: float = 3.0
) -> list[int]:
    """Boundary indices from recursive binary segmentation of a real-valued series.

    Stands in for the spec's PCF-at-gamma-10 and CBS. Binary segmentation is chosen deliberately
    over PCF: PCF at that parameter is the one remaining haplarithmisis-derived constant in the
    design, and the haplarithmisis BAF machinery is patented (WO/2015/028576, licensed to
    Agilent). Binary segmentation is textbook, unencumbered, and does the same job in the
    agreement rule - which only needs SOME continuous method to corroborate a boundary.

    Only ONE continuous method is implemented. The spec asks for two so their agreement is
    informative rather than tautological; the second is not built, and the agreement rule below
    reports how many corroborated rather than pretending two did.
    """
    bounds: list[int] = []

    def split(lo: int, hi: int) -> None:
        n = hi - lo
        if n < 2 * min_size:
            return
        total = sum(values[lo:hi])
        best_stat, best_at = 0.0, -1
        left_sum = 0.0
        for t in range(lo + min_size, hi - min_size + 1):
            left_sum = sum(values[lo:t]) if t == lo + min_size else left_sum + values[t - 1]
            nl, nr = t - lo, hi - t
            ml, mr = left_sum / nl, (total - left_sum) / nr
            stat = abs(ml - mr) * math.sqrt(nl * nr / n)
            if stat > best_stat:
                best_stat, best_at = stat, t
        # Accept only if the shift is large against the noise scale.
        if best_at < 0 or sigma <= 0 or best_stat < penalty * sigma:
            return
        bounds.append(best_at)
        split(lo, best_at)
        split(best_at, hi)

    split(0, len(values))
    return sorted(bounds)


# --- the agreement rule (spec 7.6 step 5) -------------------------------------------------


class Boundary(NamedTuple):
    index: int
    #: How many of the three conditions held.
    viterbi_change: bool
    continuous_support: bool
    run_significant: bool

    @property
    def accepted(self) -> bool:
        return self.viterbi_change and self.continuous_support and self.run_significant


def viterbi_boundaries(path: Sequence[int]) -> list[int]:
    """Indices where the Viterbi path changes COPY state.

    Deliberately ignores T-only changes: a crossover is not a structural boundary and folding
    the two together is what made the earlier design unable to tell them apart.
    """
    out = []
    for i in range(1, len(path)):
        if STATES[path[i]] != STATES[path[i - 1]]:
            out.append(i)
    return out


def reconcile(
    viterbi_bounds: Sequence[int],
    continuous_bounds: Sequence[int],
    run_significant: bool,
    *,
    tolerance: int = 2,
) -> list[Boundary]:
    """Accept a state change only where the discrete and continuous paths agree AND the run test
    is significant. Anything else is DISCORDANT and reported as a candidate, never as a call.

    The rule is conservative by construction: sensitivity is the product of three methods'
    sensitivities. That is the intended trade - a missed boundary is reported as discordant and
    stays visible, whereas a boundary accepted on one method's word alone would be invisible.
    """
    out: list[Boundary] = []
    for b in viterbi_bounds:
        near = any(abs(b - c) <= tolerance for c in continuous_bounds)
        out.append(Boundary(b, True, near, run_significant))
    # Continuous boundaries with no Viterbi partner are candidates too, and are reported rather
    # than dropped: the discrete path missing one is itself information.
    for c in continuous_bounds:
        if not any(abs(b - c) <= tolerance for b in viterbi_bounds):
            out.append(Boundary(c, False, True, run_significant))
    return sorted(out, key=lambda x: x.index)


# --- the reportable hypotheses (spec 7.1) -------------------------------------------------

Hypothesis = Literal[
    "H1_paternal_wt_inherited",
    "H2_repaired",
    "H3_paternal_segment_absent",
    "H3_paternal_chromosome_absent",
    "insufficient_evidence",
]

_ABSENT = ("PAT0_MAT1", "PAT0_MAT2", "PAT0_MAT0")


@dataclass
class HypothesisPosterior:
    posterior: dict[str, float] = field(default_factory=dict)
    refusals: list[str] = field(default_factory=list)
    #: Marker support that produced it, so a number never travels without its basis.
    informative_markers: int = 0


def hypotheses_at(
    gamma_row: Sequence[float],
    *,
    whole_chromosome_absent: bool,
    informative_markers: int,
    min_markers_per_side: int = 2,
) -> HypothesisPosterior:
    """Map the state posterior at the variant onto the reportable hypotheses.

    Two refusals are structural rather than thresholds:

    H1 versus H2 differ ONLY in T, so without a declared phase source the two cannot be
    separated at all and their mass is reported pooled. That is the same refusal the
    informativity and emission layers make, arriving here from the signal side.

    And a posterior resting on too few informative markers is reported as insufficient whatever
    its numeric value, because a confident number from two markers is still a number from two
    markers.
    """
    out = HypothesisPosterior(informative_markers=informative_markers)
    mass: dict[str, float] = {}
    for s, p in zip(STATES, gamma_row):
        if p <= 0.0:
            continue
        if s.name in _ABSENT:
            key = ("H3_paternal_chromosome_absent" if whole_chromosome_absent
                   else "H3_paternal_segment_absent")
        elif s.name == "PAT1_MAT1":
            # Present at one copy. Whether that is H1 or H2 turns on which homologue came
            # through, which this layer cannot see; `scaffold` answers it.
            key = "H1_or_H2_paternal_allele_present"
        else:
            # Gains, maternal loss and tetrasomy are real states but none of them is one of the
            # reportable hypotheses. Pooling them into H2 would be the convenient answer.
            key = "insufficient_evidence"
        mass[key] = mass.get(key, 0.0) + p

    if mass.get("H1_or_H2_paternal_allele_present"):
        out.refusals.append(
            "The paternal allele is present, but H1 and H2 differ only in WHICH paternal "
            "homologue was transmitted, which genotypes at this locus cannot show. That is "
            "answered by co-inheritance across siblings with one externally typed anchor, in "
            "`scaffold`, not here."
        )

    if informative_markers < 2 * min_markers_per_side:
        out.refusals.append(
            f"{informative_markers} informative markers is below the floor of "
            f"{min_markers_per_side} per side; reported as insufficient evidence regardless of "
            "the numeric posterior."
        )
        mass = {"insufficient_evidence": 1.0}

    total = sum(mass.values()) or 1.0
    out.posterior = {k: v / total for k, v in sorted(mass.items(), key=lambda kv: -kv[1])}
    return out


# --- the mechanism ledger (spec 6.5, 14) --------------------------------------------------

LedgerVerdict = Literal["excluded", "not_excluded", "not_tested"]


class LedgerEntry(NamedTuple):
    mechanism: str
    verdict: LedgerVerdict
    detail: str


def mechanism_ledger(
    *,
    run_significant: bool,
    copy_neutral_separation_sigma: float,
    marker_spacing_bp: Optional[float],
    r_min: Optional[int],
    has_lrr: bool,
    has_reads: bool,
    construct_sequence_supplied: bool,
) -> list[LedgerEntry]:
    """Per-mechanism verdicts, with the resolution each was tested at.

    The rule this enforces is the one the whole design turns on: a mechanism that could not be
    tested is reported NOT TESTED, never excluded. An array cannot see an insertion, so on array
    input H3f is never excluded no matter how clean the data is.
    """
    floor = (r_min * marker_spacing_bp) if (r_min and marker_spacing_bp) else None
    out: list[LedgerEntry] = []

    out.append(LedgerEntry(
        "H3a stochastic dropout",
        "excluded" if run_significant else "not_excluded",
        "a contiguous run separates a structural event from independent dropout"
        if run_significant else
        "no significant run, so an isolated absence is indistinguishable from dropout",
    ))

    # The array has a floor: a deletion shorter than r_min markers cannot make a significant run,
    # so below it the array is blind however clean the data. Reads have no such floor, which is
    # the whole reason structural calls are worth supplying.
    if has_reads:
        b_verdict = "excluded" if (run_significant or construct_sequence_supplied) \
            else "not_excluded"
        b_detail = ("structural calls cover this locus, so an event below the array's resolution "
                    "floor is visible here where the array alone could not exclude it")
    elif floor is None:
        b_verdict, b_detail = "not_tested", ("no resolution floor computable without r_min and "
                                             "marker spacing")
    else:
        b_verdict = "excluded" if run_significant else "not_excluded"
        b_detail = (f"array resolution floor is r_min x marker spacing = {floor:.0f} bp; events "
                    f"below that cannot produce a significant run and are NOT excluded")
    out.append(LedgerEntry("H3b large on-target deletion", b_verdict, b_detail))

    for name in ("H3c segmental paternal loss", "H3d whole-chromosome loss"):
        out.append(LedgerEntry(
            name,
            "not_tested" if not has_lrr else ("excluded" if run_significant else "not_excluded"),
            "needs the intensity channel to distinguish loss from copy-neutral change"
            if not has_lrr else "tested on BAF and LRR jointly",
        ))

    # The entry the real data turned into a live constraint.
    if not has_lrr:
        detail = "no LRR channel, so copy-neutral LOH cannot be distinguished from copy loss"
        verdict: LedgerVerdict = "not_tested"
    elif copy_neutral_separation_sigma < 2.0:
        detail = (
            f"copy-loss and copy-neutral states are only "
            f"{copy_neutral_separation_sigma:.2f} sigma apart on this sample, below the 2 sigma "
            "needed to choose between them. Gene conversion and paternal deletion both remove "
            "paternal alleles contiguously and differ ONLY here, so this is the distinction the "
            "sample cannot make"
        )
        verdict = "not_tested"
    else:
        detail = (f"separated by {copy_neutral_separation_sigma:.2f} sigma in LRR")
        verdict = "excluded" if run_significant else "not_excluded"
    out.append(LedgerEntry("H3e copy-neutral LOH", verdict, detail))

    if not has_reads:
        f_detail = "array input cannot see an insertion or a balanced rearrangement at all"
    elif not construct_sequence_supplied:
        f_detail = "reads available but no delivery-construct sequence, so a vector insertion " \
                   "cannot be identified"
    else:
        f_detail = "tested on split reads and discordant pairs"
    out.append(LedgerEntry(
        "H3f insertion or balanced rearrangement",
        "not_tested" if not (has_reads and construct_sequence_supplied) else "not_excluded",
        f_detail,
    ))

    return out

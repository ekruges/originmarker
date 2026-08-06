"""Phase-3 inference: transitions, forward-backward, Viterbi, segmentation, agreement, ledger.

The properties worth pinning here are the ones a passing-but-wrong implementation would still
satisfy: that the joint state really does separate a crossover from a deletion, that T goes flat
rather than arbitrary where it is unidentifiable, that the agreement rule refuses rather than
resolves, and that a mechanism which could not be tested is never reported as excluded.
"""

from __future__ import annotations

import math

import pytest

import emissions as em
import hmm

FA, MO = "AA", "BB"
PARAMS = em.EmissionParams(lrr_compression=0.78, sigma_baf=0.03, sigma_lrr=0.12, eps=0.01)
TP = hmm.TransitionParams(event_length_bp=200_000.0)


def markers(specs, chrom_pos_step=10_000, start=1_000_000):
    """Build a marker run from (baf, lrr) pairs at fixed spacing."""
    return [
        hmm.MarkerObs(rsid=f"rs{i}", pos=start + i * chrom_pos_step, father=FA, mother=MO,
                      baf=b, lrr=l)
        for i, (b, l) in enumerate(specs)
    ]


# --- state space ---------------------------------------------------------------------------


def test_entry_and_exit_hazards_are_asymmetric():
    """A single hazard would make expected event length equal expected normal length, i.e. a
    prior asserting half the genome is aberrant. Entry must be far rarer than exit."""
    h = hmm.breakpoint_hazard(10_000, TP)
    # The rate ratio is exit/enter = (1/event_length) / (events_per_mb/1e6), which is exactly 100
    # for these defaults; at a finite gap the exponential curvature makes it marginally less.
    assert h.exit / h.enter > 90, f"entry {h.enter} must be far rarer than exit {h.exit}"
    assert hmm.breakpoint_hazard(0, TP) == (0.0, 0.0)
    # Exit saturates at the characteristic length, which is what makes an event bounded.
    assert hmm.breakpoint_hazard(int(TP.event_length_bp), TP).exit == pytest.approx(1 - math.exp(-1))
    # Both rise monotonically with the gap and saturate at 1.
    gaps = (1_000, 100_000, 10_000_000)
    enters = [hmm.breakpoint_hazard(g, TP).enter for g in gaps]
    assert enters[0] < enters[1] < enters[2] <= 1.0


def test_crossover_uses_the_bundled_map_and_flags_the_sex_averaged_fallback():
    """Spec 9.4 wants the MALE map and no redistributable licence-clear one was found, so every
    affected output has to carry the label rather than look quietly correct."""
    mu = hmm.crossover("11", 17_147_055, 17_397_055)
    assert mu.cm > 0.0
    assert 0.0 < mu.theta < 0.5
    assert mu.sex_averaged_fallback is True
    assert "deCODE" in mu.source or mu.approx
    # A zero-length step cannot recombine.
    assert hmm.crossover("11", 17_147_055, 17_147_055).theta == pytest.approx(0.0)


def test_transition_rows_are_proper_distributions():
    """Every row must sum to one, from normal AND from an aberrant state, or the forward pass
    silently leaks probability mass."""
    n_ab = len(em.STATES) - 1
    if True:
        for hz in (hmm.Hazards(0.001, 0.05), hmm.Hazards(0.3, 0.9)):
            for a in hmm.STATES:
                total = sum(
                    math.exp(hmm.transition_logp(a, b, hazards=hz, n_aberrant=n_ab))
                    for b in hmm.STATES
                )
                assert total == pytest.approx(1.0, abs=1e-9), f"row {a.name}"


def test_an_aberrant_state_returns_to_normal_rather_than_hopping_sideways():
    """Events are excursions from normal and back. Two adjacent distinct events are represented
    as event-normal-event, which is a stated limitation rather than an invented prior over
    event pairs."""
    kw = dict(hazards=hmm.Hazards(0.001, 0.05), n_aberrant=len(em.STATES) - 1)
    loss = em.BY_NAME["PAT0_MAT1"]
    gain = em.BY_NAME["PAT2_MAT1"]
    assert hmm.transition_logp(loss, gain, **kw) == -math.inf
    assert hmm.transition_logp(loss, em.NORMAL, **kw) > -math.inf


# --- inference -----------------------------------------------------------------------------


def test_forward_backward_posteriors_are_normalised():
    ms = markers([(0.5, 0.0)] * 12)
    inf = hmm.infer(ms, "11", PARAMS, TP)
    assert len(inf.gamma) == 12
    for row in inf.gamma:
        assert sum(row) == pytest.approx(1.0, abs=1e-9)
    # A log DENSITY over continuous BAF and LRR, not a log probability, so it is routinely
    # positive: at sigma_BAF 0.03 the peak density alone is about 13.
    assert math.isfinite(inf.log_likelihood)
    assert len(inf.viterbi) == 12
    assert len(inf.map_uses) == 11


def test_a_clean_diploid_run_is_called_normal_throughout():
    ms = markers([(0.5, 0.0)] * 20)
    inf = hmm.infer(ms, "11", PARAMS, TP)
    assert all(hmm.STATES[j].name == "PAT1_MAT1" for j in inf.viterbi)


def test_a_paternal_loss_segment_is_localised():
    """Markers 8-13 lose the paternal copy: BAF to 1.0 and LRR to the CN1 mean."""
    specs = [(0.5, 0.0)] * 8 + [(1.0, -0.78)] * 6 + [(0.5, 0.0)] * 8
    inf = hmm.infer(markers(specs), "11", PARAMS, TP)
    called = [hmm.STATES[j].name for j in inf.viterbi]
    assert called[:8] == ["PAT1_MAT1"] * 8
    assert all(n == "PAT0_MAT1" for n in called[8:14]), called
    assert called[14:] == ["PAT1_MAT1"] * 8


def test_copy_neutral_and_copy_loss_are_told_apart_by_lrr_alone():
    """The distinction the layer exists for, end to end. Identical BAF, different LRR."""
    base = [(0.5, 0.0)] * 6
    loss = hmm.infer(markers(base + [(1.0, -0.78)] * 8 + base), "11", PARAMS, TP)
    neutral = hmm.infer(markers(base + [(1.0, 0.0)] * 8 + base), "11", PARAMS, TP)
    assert {hmm.STATES[j].name for j in loss.viterbi[6:14]} == {"PAT0_MAT1"}
    assert {hmm.STATES[j].name for j in neutral.viterbi[6:14]} == {"PAT0_MAT2"}


def test_infer_handles_the_empty_window():
    inf = hmm.infer([], "11", PARAMS, TP)
    assert inf.gamma == [] and inf.viterbi == [] and inf.log_likelihood == 0.0


# --- segmentation and the agreement rule ---------------------------------------------------


def test_binary_segmentation_finds_a_real_shift_and_ignores_noise():
    flat = [0.0, 0.01, -0.01, 0.005] * 6
    assert hmm.binary_segment(flat, sigma=0.12) == []
    stepped = [0.0] * 12 + [-0.78] * 12
    bounds = hmm.binary_segment(stepped, sigma=0.12)
    assert bounds and min(abs(b - 12) for b in bounds) <= 1


def test_agreement_requires_all_three_conditions():
    both = hmm.reconcile([10], [10], True)
    assert len(both) == 1 and both[0].accepted

    # Continuous method disagrees -> candidate, not a call.
    disc = hmm.reconcile([10], [40], True)
    assert not any(b.accepted for b in disc)
    assert len(disc) == 2, "both the unmatched discrete and continuous boundaries are reported"

    # The run test failing is enough on its own to withhold acceptance.
    assert not hmm.reconcile([10], [10], False)[0].accepted
    # Within tolerance is agreement; outside it is not.
    assert hmm.reconcile([10], [12], True)[0].accepted
    assert not hmm.reconcile([10], [13], True)[0].accepted


def test_a_continuous_only_boundary_is_reported_not_dropped():
    """The discrete path missing a boundary the continuous one found is itself information."""
    out = hmm.reconcile([], [25], True)
    assert len(out) == 1
    assert out[0].index == 25 and not out[0].viterbi_change and out[0].continuous_support
    assert not out[0].accepted


# --- reportable hypotheses ------------------------------------------------------------------


def gamma_on(state_name: str) -> list[float]:
    return [1.0 if s == em.BY_NAME[state_name] else 0.0 for s in hmm.STATES]


def test_paternal_absence_maps_to_segment_or_whole_chromosome():
    seg = hmm.hypotheses_at(gamma_on("PAT0_MAT1"), whole_chromosome_absent=False, informative_markers=20)
    assert seg.posterior["H3_paternal_segment_absent"] == pytest.approx(1.0)
    whole = hmm.hypotheses_at(gamma_on("PAT0_MAT1"), whole_chromosome_absent=True, informative_markers=20)
    assert whole.posterior["H3_paternal_chromosome_absent"] == pytest.approx(1.0)
    # Maternal isodisomy also removes the paternal allele, so it belongs to the same hypothesis.
    assert hmm.hypotheses_at(gamma_on("PAT0_MAT2"), whole_chromosome_absent=False,
                             informative_markers=20).posterior["H3_paternal_segment_absent"] == 1.0


def test_gains_are_not_pooled_into_h2():
    # There is no insertion state and no H2_with_local_artefact hypothesis. Both were removed in
    # 3.1.3: an array cannot see an insertion at any level, and a state whose emissions match
    # PAT1_MAT1 exactly only takes posterior mass away from it.
    assert "H2_with_local_artefact" not in hmm.Hypothesis.__args__
    # A gain is a real state but not a reportable hypothesis. Folding it into H2 would be the
    # convenient answer rather than the true one.
    gain = hmm.hypotheses_at(gamma_on("PAT2_MAT1"), whole_chromosome_absent=False, informative_markers=20)
    assert gain.posterior["insufficient_evidence"] == pytest.approx(1.0)
    assert hmm.hypotheses_at(gamma_on("PAT1_MAT0"), whole_chromosome_absent=False,
                             informative_markers=20).posterior["insufficient_evidence"] == 1.0


def test_too_few_markers_is_insufficient_whatever_the_number_says():
    thin = hmm.hypotheses_at(gamma_on("PAT0_MAT1"), whole_chromosome_absent=False, informative_markers=3)
    assert thin.posterior == {"insufficient_evidence": 1.0}
    assert any("below the floor" in r for r in thin.refusals)
    assert thin.informative_markers == 3


# --- the mechanism ledger -------------------------------------------------------------------


def ledger(**kw):
    base = dict(run_significant=True, copy_neutral_separation_sigma=6.5,
                marker_spacing_bp=20_000.0, r_min=2, has_lrr=True, has_reads=False,
                construct_sequence_supplied=False)
    base.update(kw)
    return {e.mechanism: e for e in hmm.mechanism_ledger(**base)}


def test_an_untestable_mechanism_is_never_reported_as_excluded():
    """The rule the whole design turns on."""
    led = ledger()
    assert led["H3f insertion or balanced rearrangement"].verdict == "not_tested"
    assert "cannot see an insertion" in led["H3f insertion or balanced rearrangement"].detail
    # Even with perfect data, array input never excludes H3f.
    assert ledger(copy_neutral_separation_sigma=99.0)[
        "H3f insertion or balanced rearrangement"].verdict == "not_tested"


def test_copy_neutral_loh_is_not_tested_when_the_sample_cannot_separate_it():
    """The real sample's fitted parameters put this at 1.32 sigma, below the 2 sigma needed.

    This is the constraint that matters: gene conversion and paternal deletion both remove
    paternal alleles contiguously and differ ONLY in copy number, so a sample that cannot resolve
    that difference cannot answer the question Ma's reply stalled on.
    """
    poor = ledger(copy_neutral_separation_sigma=1.32)
    assert poor["H3e copy-neutral LOH"].verdict == "not_tested"
    assert "1.32 sigma" in poor["H3e copy-neutral LOH"].detail
    good = ledger(copy_neutral_separation_sigma=6.5)
    assert good["H3e copy-neutral LOH"].verdict == "excluded"
    assert ledger(has_lrr=False)["H3e copy-neutral LOH"].verdict == "not_tested"


def test_dropout_is_excluded_only_by_contiguity():
    assert ledger(run_significant=True)["H3a stochastic dropout"].verdict == "excluded"
    assert ledger(run_significant=False)["H3a stochastic dropout"].verdict == "not_excluded"


def test_the_deletion_entry_carries_its_resolution_floor():
    led = ledger(r_min=2, marker_spacing_bp=20_000.0)
    e = led["H3b large on-target deletion"]
    assert "40000 bp" in e.detail
    # A noisier sample needs a longer run, so its floor is worse and the entry must say so.
    assert "200000 bp" in ledger(r_min=10, marker_spacing_bp=20_000.0)[
        "H3b large on-target deletion"].detail
    assert ledger(r_min=None)["H3b large on-target deletion"].verdict == "not_tested"


def test_segmental_and_whole_chromosome_need_the_intensity_channel():
    no_lrr = ledger(has_lrr=False)
    for m in ("H3c segmental paternal loss", "H3d whole-chromosome loss"):
        assert no_lrr[m].verdict == "not_tested"
        assert "intensity channel" in no_lrr[m].detail


def test_every_mechanism_appears_exactly_once():
    entries = hmm.mechanism_ledger(
        run_significant=False, copy_neutral_separation_sigma=1.0, marker_spacing_bp=None,
        r_min=None, has_lrr=False, has_reads=False, construct_sequence_supplied=False)
    names = [e.mechanism for e in entries]
    assert len(names) == len(set(names)) == 6
    assert all(e.verdict in ("excluded", "not_excluded", "not_tested") for e in entries)
    assert all(e.detail for e in entries), "a verdict without a reason is not reportable"


def test_the_homologue_axis_is_gone_and_not_merely_hidden():
    """It used to derive the transmitted allele from a state T by the convention "the mutant
    homologue carries the A allele". Which of a father's alleles is labelled A at each marker is
    arbitrary until he is phased, so that axis tracked a relabelling rather than a chromosome.
    `scaffold` answers the question properly, from co-inheritance across siblings."""
    assert not hasattr(hmm, "HOMOLOGUES") and not hasattr(hmm, "State")
    assert not hasattr(hmm, "t_identifiable")
    assert hmm.STATES == em.STATES, "the states ARE the copy states"
    # Nine, not ten: the insertion state was removed in 3.1.3 as array-invisible.
    assert len(hmm.STATES) == 9


def test_h1_and_h2_are_never_split_here_whatever_is_passed():
    """The old code split them on T and pooled them only when no phase was declared, so a caller
    that declared one got an answer resting on a letter convention. There is now no way to ask."""
    got = hmm.hypotheses_at(gamma_on("PAT1_MAT1"), whole_chromosome_absent=False,
                            informative_markers=50)
    assert "H1_or_H2_paternal_allele_present" in got.posterior
    assert not any(k in got.posterior for k in ("H1_paternal_wt_inherited", "H2_repaired"))
    assert any("scaffold" in r for r in got.refusals)

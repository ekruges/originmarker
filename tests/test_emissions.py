"""Phase-3 emission model.

The externally-anchored tests are the mosaic ones. `array_signatures.csv` tabulates the band
positions and LRR shifts independently of this code, and it disagrees with spec 7.2's own
formulas, so reproducing the TABLE rather than the text is the whole point. The magnitudes of
those disagreements are asserted too, so nobody quietly reverts to the spec's version.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

import emissions as em

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"

# A marker where both parents are homozygous for opposite alleles, so the transmitted pair is
# fully determined and every band position is unambiguous. This is the fully-informative class.
FA, MO = "AA", "BB"


def band_of(state_name: str, f: float = 1.0, c: float = 1.0) -> em.Band:
    got = em.bands(em.BY_NAME[state_name], FA, MO, mosaic_fraction=f, lrr_compression=c)
    assert len(got) == 1, "both parents homozygous means exactly one band"
    return got[0]


# --- clonal states, against array_signatures.csv -------------------------------------------


@pytest.mark.parametrize(
    "state,expected_baf,expected_lrr",
    [
        ("PAT1_MAT1", 0.5, 0.0),          # normal diploid, het band
        ("PAT0_MAT1", 1.0, -1.0),         # paternal loss: only the maternal B survives
        ("PAT1_MAT0", 0.0, -1.0),         # the mirror
        ("PAT2_MAT0", 0.0, 0.0),          # PATERNAL isodisomy: LOH with FLAT LRR
        ("PAT0_MAT2", 1.0, 0.0),          # maternal isodisomy: same, opposite allele
        ("PAT2_MAT1", 1 / 3, 0.5849625),  # CN3, log2(3/2)
        ("PAT1_MAT2", 2 / 3, 0.5849625),
        ("PAT2_MAT2", 0.5, 1.0),          # CN4, log2(4/2)
    ],
)
def test_clonal_bands_match_the_signature_table(state, expected_baf, expected_lrr):
    b = band_of(state)
    assert b.baf == pytest.approx(expected_baf, abs=1e-9)
    assert b.lrr == pytest.approx(expected_lrr, abs=1e-6)


def test_copy_neutral_loh_is_flat_in_lrr_and_that_is_the_whole_point():
    """The distinction the layer exists for.

    Paternal loss and maternal isodisomy predict the SAME allelic pattern - paternal alleles
    gone - and differ ONLY in intensity. Gene conversion and paternal segmental deletion are the
    biological versions of that pair, and Ma's reply rested on LOH without copy number, which is
    exactly why it could not choose between them.
    """
    loss = band_of("PAT0_MAT1")
    neutral = band_of("PAT0_MAT2")
    assert loss.baf == neutral.baf == 1.0, "identical allelic reading"
    assert loss.lrr == pytest.approx(-1.0)
    assert neutral.lrr == pytest.approx(0.0)
    assert abs(loss.lrr - neutral.lrr) == pytest.approx(1.0)


def test_homozygous_deletion_has_no_bands_to_predict():
    """CN0 leaves no signal to interpolate theta from, so BAF is noise rather than a band."""
    assert em.bands(em.BY_NAME["PAT0_MAT0"], FA, MO) == []


# --- mosaic: the ratio of expectations, NOT a linear interpolation --------------------------

# Verbatim from array_signatures.csv, mosaic_loss_CN1_fraction_f.
MOSAIC_LOSS = {0.05: (0.487, 0.513, -0.037), 0.10: (0.474, 0.526, -0.074),
               0.20: (0.444, 0.556, -0.152), 0.30: (0.412, 0.588, -0.235),
               0.50: (0.333, 0.667, -0.415), 1.00: (0.000, 1.000, -1.000)}

# mosaic_gain_CN3_fraction_f.
MOSAIC_GAIN = {0.05: (0.488, 0.512, 0.036), 0.10: (0.476, 0.524, 0.070),
               0.20: (0.455, 0.545, 0.138), 0.30: (0.435, 0.565, 0.202),
               0.50: (0.400, 0.600, 0.322), 1.00: (0.333, 0.667, 0.585)}

# mosaic_CN_neutral_LOH_fraction_f: LRR is exactly zero at every f.
MOSAIC_CNLOH = {0.05: (0.475, 0.525), 0.10: (0.450, 0.550), 0.20: (0.400, 0.600),
                0.30: (0.350, 0.650), 0.50: (0.250, 0.750), 1.00: (0.000, 1.000)}


@pytest.mark.parametrize("f,expected", sorted(MOSAIC_LOSS.items()))
def test_mosaic_loss_reproduces_the_table(f, expected):
    lo, hi, lrr = expected
    # (1-f)/(2-f) is losing the maternal copy, 1/(2-f) is losing the paternal one.
    assert band_of("PAT1_MAT0", f=f).baf == pytest.approx(lo, abs=5e-4)
    assert band_of("PAT0_MAT1", f=f).baf == pytest.approx(hi, abs=5e-4)

    # Assert the closed form exactly, and the table only at the precision it is quoted to.
    # The table's f = 0.30 row reads -0.235 where log2(0.85) is -0.234465, which rounds to
    # -0.234: a rounding slip in that one row. The other five round correctly, and the gain
    # side at f = 0.30 does too, so it is an isolated typo rather than a different formula.
    exact = math.log2((f * 1 + (1 - f) * 2) / 2)
    assert band_of("PAT0_MAT1", f=f).lrr == pytest.approx(exact, abs=1e-12)
    assert band_of("PAT0_MAT1", f=f).lrr == pytest.approx(lrr, abs=1e-3)


@pytest.mark.parametrize("f,expected", sorted(MOSAIC_GAIN.items()))
def test_mosaic_gain_reproduces_the_table(f, expected):
    lo, hi, lrr = expected
    assert band_of("PAT2_MAT1", f=f).baf == pytest.approx(lo, abs=5e-4)
    assert band_of("PAT1_MAT2", f=f).baf == pytest.approx(hi, abs=5e-4)
    assert band_of("PAT2_MAT1", f=f).lrr == pytest.approx(lrr, abs=5e-4)


@pytest.mark.parametrize("f,expected", sorted(MOSAIC_CNLOH.items()))
def test_mosaic_copy_neutral_loh_reproduces_the_table(f, expected):
    lo, hi = expected
    assert band_of("PAT2_MAT0", f=f).baf == pytest.approx(lo, abs=5e-4)
    assert band_of("PAT0_MAT2", f=f).baf == pytest.approx(hi, abs=5e-4)
    # Zero at EVERY mosaic fraction. This is why a depth-based caller can never see this state.
    assert band_of("PAT2_MAT0", f=f).lrr == pytest.approx(0.0, abs=1e-12)
    assert band_of("PAT0_MAT2", f=f).lrr == pytest.approx(0.0, abs=1e-12)


def test_the_specs_linear_mosaic_forms_are_wrong_and_by_how_much():
    """Spec 7.2 interpolates linearly in the OBSERVED quantity. Both observations are nonlinear
    functions of copy number, so the mixture has to be taken in copy number instead.

    Asserted with magnitudes because the BAF error is large enough to matter: at a typical
    sigma_BAF of 0.03 it is nearly three standard deviations, applied to every mosaic state.
    """
    f = 0.50
    loss = em.BY_NAME["PAT0_MAT1"]

    spec_baf = f * 1.0 + (1.0 - f) * 0.5           # linear in BAF, as written
    correct_baf = band_of("PAT0_MAT1", f=f).baf
    assert spec_baf == pytest.approx(0.750)
    assert correct_baf == pytest.approx(0.667, abs=5e-4)
    assert abs(spec_baf - correct_baf) == pytest.approx(0.0833, abs=5e-4)
    assert abs(spec_baf - correct_baf) / 0.03 > 2.7, "over 2.7 sigma at sigma_BAF = 0.03"

    spec_lrr = f * math.log2(loss.cn / 2) + (1.0 - f) * 0.0   # linear in LRR, as written
    correct_lrr = band_of("PAT0_MAT1", f=f).lrr
    assert spec_lrr == pytest.approx(-0.500)
    assert correct_lrr == pytest.approx(-0.415, abs=5e-4)
    assert abs(spec_lrr - correct_lrr) == pytest.approx(0.0849, abs=5e-4)


def test_mosaic_fraction_is_validated():
    with pytest.raises(ValueError):
        em.bands(em.NORMAL, FA, MO, mosaic_fraction=1.5)
    with pytest.raises(ValueError):
        em.bands(em.NORMAL, FA, MO, mosaic_fraction=-0.1)


# --- compression is fitted, never assumed --------------------------------------------------


def test_lrr_compression_fits_the_real_samples_values():
    """The real `.CEL.probes` sample gives CN1 at mean log2R -0.454 and CN3 at +0.339.

    Those are the numbers the fitter has to recover, and they are roughly half the 0.70-0.85
    spec 7.2 says to expect from Illumina vendor figures. c cannot be a constant.
    """
    obs = [(1, -0.454)] * 200 + [(3, 0.339)] * 400 + [(2, -0.045)] * 1000
    fit = em.fit_lrr_compression(obs)
    assert fit.per_cn[1] == pytest.approx(0.454, abs=1e-3)
    assert fit.per_cn[3] == pytest.approx(0.339 / math.log2(1.5), abs=1e-3)
    assert 2 not in fit.per_cn, "CN2 constrains no scale: log2(2/2) is zero"
    # The spread between states is large on this sample, and the fitter must say so rather than
    # hand back one averaged number that misfits both.
    assert "prefer the per-state values" in fit.note


def test_lrr_compression_refuses_when_nothing_constrains_the_scale():
    fit = em.fit_lrr_compression([(2, 0.01)] * 500)
    assert math.isnan(fit.pooled)
    assert "cannot be fitted" in fit.note
    assert em.fit_lrr_compression([]).per_cn == {}


def test_sigma_is_robust_to_the_events_being_looked_for():
    """A real segmental loss puts a block of markers far from the mean. A plain standard
    deviation absorbs them into the noise estimate and is then blind to them; the MAD does not.
    """
    clean = [0.01, -0.02, 0.00, 0.015, -0.01] * 40
    contaminated = clean + [-1.0] * 20            # a real event, 10% of the markers
    s_clean = em.fit_sigma(clean)
    s_dirty = em.fit_sigma(contaminated)
    assert s_dirty == pytest.approx(s_clean, rel=0.35), "MAD barely moves"
    # A standard deviation would roughly triple over the same contamination.
    def sd(xs):
        m = sum(xs) / len(xs)
        return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    assert sd(contaminated) > 3 * sd(clean)
    assert math.isnan(em.fit_sigma([0.1]))


# --- two channels, not three ----------------------------------------------------------------

PARAMS = em.EmissionParams(lrr_compression=0.78, sigma_baf=0.03, sigma_lrr=0.12, eps=0.01)


def logp(state: str, **kw) -> float:
    return em.emission_logp(em.BY_NAME[state], FA, MO, params=PARAMS, **kw)


def test_genotype_is_ignored_when_baf_is_present():
    """The allelic channel is BAF or genotype, never both: they are two readings of one theta,
    and multiplying them counts the same measurement twice. That was the error in spec 7.2."""
    with_gt = logp("PAT1_MAT1", baf=0.5, lrr=0.0, genotype="AB")
    without = logp("PAT1_MAT1", baf=0.5, lrr=0.0, genotype="NC")
    assert with_gt == pytest.approx(without), "genotype must not contribute alongside BAF"
    # Even a CONTRADICTORY genotype cannot move the answer while BAF is present.
    contradictory = logp("PAT1_MAT1", baf=0.5, lrr=0.0, genotype="AA")
    assert contradictory == pytest.approx(without)


def test_genotype_carries_the_allelic_channel_when_baf_is_missing():
    """Real Axiom files leave BAF empty at about 4% of markers, so the fallback is not academic."""
    het = logp("PAT1_MAT1", baf=None, lrr=0.0, genotype="AB")
    hom = logp("PAT1_MAT1", baf=None, lrr=0.0, genotype="AA")
    assert het > hom, "a normal diploid predicts the heterozygote"
    assert logp("PAT0_MAT1", baf=None, lrr=-0.78, genotype="BB") > \
           logp("PAT0_MAT1", baf=None, lrr=-0.78, genotype="AB")


def test_a_missing_channel_costs_nothing_rather_than_penalising_every_state():
    """A marker with no LRR is less informative, not evidence against everything."""
    both = logp("PAT1_MAT1", baf=0.5, lrr=0.0)
    baf_only = logp("PAT1_MAT1", baf=0.5, lrr=None)
    assert baf_only > -math.inf
    assert both != baf_only
    # With no channels at all the marker contributes exactly zero to every state, so it cannot
    # tilt the posterior.
    for s in ("PAT1_MAT1", "PAT0_MAT1", "PAT2_MAT0"):
        assert logp(s, baf=None, lrr=None, genotype="NC") == 0.0


def test_truncation_stops_a_systematic_bias_against_allele_loss():
    """BAF is bounded to [0,1] and the homozygous bands sit ON the boundary, so an untruncated
    normal loses about half its mass off the support there while the het band at 0.5 loses none.
    Without truncation every allele-loss state is penalised by roughly a factor of two - which
    is to say, exactly the states the tool exists to detect."""
    # Observed exactly on the band, so only the normalisation differs between the two states.
    at_hom = logp("PAT0_MAT1", baf=1.0, lrr=-0.78)
    at_het = logp("PAT1_MAT1", baf=0.5, lrr=0.0)
    # Truncation roughly doubles the density at a boundary band: log 2 = 0.693.
    assert at_hom - at_het == pytest.approx(math.log(2.0), abs=0.02)
    assert logp("PAT1_MAT1", baf=1.5, lrr=0.0) == -math.inf, "outside [0,1] is impossible"


# --- state discrimination on real parameters -------------------------------------------------


def test_the_correct_state_wins_on_its_own_observation():
    cases = {
        "PAT1_MAT1": (0.50, 0.00),
        "PAT0_MAT1": (1.00, -0.78),
        "PAT1_MAT0": (0.00, -0.78),
        "PAT0_MAT2": (1.00, 0.00),
        "PAT2_MAT0": (0.00, 0.00),
        "PAT2_MAT1": (1 / 3, 0.456),
        "PAT1_MAT2": (2 / 3, 0.456),
    }
    for truth, (baf, lrr) in cases.items():
        scores = {s.name: logp(s.name, baf=baf, lrr=lrr) for s in em.STATES if s.cn > 0}
        best = max(scores, key=scores.get)
        assert best == truth, f"{truth} observation was best explained by {best}: {scores}"


def test_copy_neutral_and_copy_loss_are_separated_by_lrr_alone():
    """The C3-versus-C5 pair from the golden fixture, at the model level: identical BAF, and only
    the intensity channel can choose. Drop LRR and the two become indistinguishable, which is
    what the mechanism ledger must then report as NOT TESTED rather than excluded."""
    loss = logp("PAT0_MAT1", baf=1.0, lrr=-0.78)
    neutral = logp("PAT0_MAT2", baf=1.0, lrr=-0.78)
    assert loss > neutral + 5, "with LRR, copy-loss wins decisively"

    loss_nolrr = logp("PAT0_MAT1", baf=1.0, lrr=None)
    neutral_nolrr = logp("PAT0_MAT2", baf=1.0, lrr=None)
    assert loss_nolrr == pytest.approx(neutral_nolrr), "without LRR they are the same hypothesis"


def test_separation_is_reported_in_sigma_units_so_a_sample_can_be_judged():
    """Whether the layer can work at all on a sample is a number, and it should be surfaced."""
    good = em.separates_copy_neutral(FA, MO, lrr_compression=0.78, sigma_lrr=0.12)
    assert good["delta_lrr"] == pytest.approx(0.78, abs=1e-9)
    assert good["sigmas_per_marker"] == pytest.approx(6.5, abs=0.1)

    # The metric that matters is per SEGMENT, and defaulting to one marker is the trap: on the
    # real chrX-calibrated values (c 0.596, sigma 0.417) a single marker gives 1.4 sigma, which I
    # wrongly reported as "cannot make the call", while 25 markers give 7.
    real = dict(lrr_compression=0.596, sigma_lrr=0.417)
    assert em.separates_copy_neutral(FA, MO, **real)["sigmas"] < 2.0
    assert em.separates_copy_neutral(FA, MO, **real, segment_markers=25)["sigmas"] > 5.0
    # sqrt scaling, not linear.
    a = em.separates_copy_neutral(FA, MO, **real, segment_markers=100)["sigmas"]
    b = em.separates_copy_neutral(FA, MO, **real, segment_markers=400)["sigmas"]
    assert b == pytest.approx(2 * a, rel=1e-9)


# --- marginalisation over unobserved transmissions ------------------------------------------


def test_a_heterozygous_mother_is_marginalised_not_guessed():
    """Which allele a heterozygous parent transmitted is simply not observed, so both branches
    stay live at half weight. Picking one would invent information."""
    got = em.bands(em.NORMAL, "AA", "AB")
    assert len(got) == 2
    assert sorted(round(b.baf, 6) for b in got) == [0.0, 0.5]
    assert all(b.weight == 0.5 for b in got)


def test_a_heterozygous_father_collapses_once_phase_is_declared():
    """Without a phase source both paternal branches are live; with one, the declared allele
    picks the branch. This is the same refusal the informativity layer makes, at the signal
    level."""
    undeclared = em.bands(em.NORMAL, "AB", "AA")
    assert len(undeclared) == 2 and all(b.weight == 0.5 for b in undeclared)
    declared = em.bands(em.NORMAL, "AB", "AA", paternal_allele="B")
    assert len(declared) == 1 and declared[0].weight == 1.0
    assert declared[0].baf == pytest.approx(0.5)


def test_no_maternal_genotype_leaves_her_contribution_unconstrained():
    got = em.bands(em.NORMAL, "AA", None)
    assert len(got) == 2, "no mother: both maternal alleles stay possible"
    assert em.bands(em.NORMAL, "NC", "BB") == [], "father no-call: nothing to predict"


def test_both_parents_heterozygous_is_nearly_uninformative():
    """Every state predicts almost the same mixture there, which is the signal-level statement of
    why that marker class carries so little information."""
    normal = em.bands(em.NORMAL, "AB", "AB")
    assert len(normal) == 4
    assert sum(b.weight for b in normal) == pytest.approx(1.0)


# --- the golden fixture, at the level this layer can be held to ------------------------------


def test_golden_fixture_c3_and_c5_differ_only_in_lrr():
    """The fixture's design lesson, checked against the file rather than restated: identical
    genotypes and identical run statistics, separated only by mean LRR inside the run."""
    path = FIXTURES / "golden_test_vectors.csv"
    rows = list(csv.DictReader(path.open()))
    c3 = [r for r in rows if r["case"] == "C3_H3c_segmental_loss"]
    c5 = [r for r in rows if r["case"] == "C5_H3e_copyneutral_LOH"]
    assert [r["embryo_gt"] for r in c3] == [r["embryo_gt"] for r in c5]
    assert [r["LRR"] for r in c3] != [r["LRR"] for r in c5]

    def mean_affected(rs):
        vals = [float(r["LRR"]) for r in rs if r["true_CN"] and r["true_state"] != "PAT1_MAT1"]
        return sum(vals) / len(vals)

    # Copy-loss sits near the CN1 mean; copy-neutral sits at zero. Fitted c = 0.78 in the fixture.
    assert mean_affected(c3) == pytest.approx(-0.78, abs=0.05)
    assert mean_affected(c5) == pytest.approx(0.0, abs=0.05)

    # And the model prefers the right state for each, from the fixture's own numbers.
    assert logp("PAT0_MAT1", baf=1.0, lrr=mean_affected(c3)) > \
           logp("PAT0_MAT2", baf=1.0, lrr=mean_affected(c3))
    assert logp("PAT0_MAT2", baf=1.0, lrr=mean_affected(c5)) > \
           logp("PAT0_MAT1", baf=1.0, lrr=mean_affected(c5))


def test_state_space_covers_the_signature_table():
    """Every copy-number state array_signatures.csv enumerates must exist here, or the model
    would have to force a real configuration into the nearest wrong one."""
    cns = {s.cn for s in em.STATES}
    assert {0, 1, 2, 3, 4} <= cns
    # Both isodisomies, which a total-copy-number state space would collapse into one.
    assert em.BY_NAME["PAT2_MAT0"].cn == em.BY_NAME["PAT0_MAT2"].cn == 2
    assert em.BY_NAME["PAT2_MAT0"].n_pat != em.BY_NAME["PAT0_MAT2"].n_pat
    # The insertion state is copy-number-nominal and array-invisible by design.
    ins = em.BY_NAME["PAT1_MAT1_INS"]
    assert ins.cn == 2 and ins.insertion
    assert band_of("PAT1_MAT1_INS").baf == pytest.approx(band_of("PAT1_MAT1").baf)
    assert band_of("PAT1_MAT1_INS").lrr == pytest.approx(band_of("PAT1_MAT1").lrr)

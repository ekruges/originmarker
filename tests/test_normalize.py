"""In-flight intensity normalisation.

Every number here was measured on real public data: Zuccaro 2020 GSE148488 (two known-male sperm
donor replicates, two known-female egg donors) and Turocy 2026 GSE186407, both UK Biobank Axiom,
825,656 markers. They are pinned as measurements rather than as fixtures because the files are
40 MB each and patient-derived even though public.

The load-bearing tests are the two that caught real errors: that the female chrX baseline must be
subtracted, and that separation is a per-SEGMENT quantity rather than a per-marker one.
"""

from __future__ import annotations

import math

import pytest

import normalize as nz

UKB = "UK Biobank Axiom Array"

# Measured medians. A male's chrX sits 0.4075 below his autosomes; a known FEMALE's sits 0.1908
# ABOVE hers despite carrying two copies, which is the array's own offset.
MALE_RAW_SHIFT = -0.4075
FEMALE_RAW_SHIFT = +0.1908


def series(median: float, spread: float, n: int) -> list[float]:
    """A symmetric series with a known median and MAD-scale, so calibration is exactly checkable."""
    half = n // 2
    step = spread / max(1, half)
    # Symmetric about `median`, with the centre element only at odd n. Including it at even n
    # leaves the list one long, and trimming from one end shifts the median by half a step.
    out = [median] if n % 2 else []
    for i in range(1, half + 1):
        out += [median - i * step, median + i * step]
    return out[:n]


def test_the_female_chrx_baseline_is_a_real_offset_not_biology():
    """A female carries two X copies, so her chrX-to-autosome shift should be zero. It is +0.191.

    That is probe design and reference-panel composition. Ignoring it is not a rounding error: it
    biases the male-derived compression by a third, from 0.598 down to 0.407.
    """
    assert nz.CHRX_CN2_OFFSET[UKB] == pytest.approx(0.1908, abs=1e-4)
    uncorrected = -MALE_RAW_SHIFT
    corrected = -(MALE_RAW_SHIFT - nz.CHRX_CN2_OFFSET[UKB])
    assert uncorrected == pytest.approx(0.4075, abs=1e-4)
    assert corrected == pytest.approx(0.5983, abs=1e-3)
    assert abs(corrected - uncorrected) / corrected > 0.30, "a third of the value"


def test_chrx_calibration_recovers_the_published_compression_range():
    """c from a male's chrX lands in PennCNV's shipped range for CN probes, 0.572-0.627.

    An independent measurement agreeing with published defaults is the strongest validation
    available here, because nothing in this derivation used PennCNV's numbers.
    """
    auto = series(0.044, 0.4, 20_000)
    x = series(0.044 + MALE_RAW_SHIFT, 0.4, 2_000)
    c, raw, warns = nz.compression_from_chrx(auto, x, sex="male", product=UKB)
    assert c is not None
    assert 0.572 <= c <= 0.627, f"c was {c}"
    assert raw == pytest.approx(MALE_RAW_SHIFT, abs=1e-3)
    assert warns == []


def test_a_female_cannot_calibrate_and_saying_so_is_the_point():
    auto = series(0.06, 0.4, 20_000)
    x = series(0.06 + FEMALE_RAW_SHIFT, 0.4, 2_000)
    c, raw, warns = nz.compression_from_chrx(auto, x, sex="female", product=UKB)
    assert c is None
    assert any("two copies" in w for w in warns)
    # Ambiguous sex is refused too: guessing male would invent a one-copy reference.
    assert nz.compression_from_chrx(auto, x, sex="ambiguous", product=UKB)[0] is None


def test_an_unmeasured_array_refuses_rather_than_borrowing_an_offset():
    """The offset is a per-array constant. Borrowing another array's would bias c by a third,
    which is exactly the error this route exists to avoid."""
    auto = series(0.0, 0.4, 20_000)
    x = series(MALE_RAW_SHIFT, 0.4, 2_000)
    c, _, warns = nz.compression_from_chrx(auto, x, sex="male", product="Some Other Axiom")
    assert c is None
    assert any("per-array constant" in w for w in warns)


def test_too_few_markers_is_refused_not_averaged():
    c, _, warns = nz.compression_from_chrx([0.0] * 100, [0.0] * 10, sex="male", product=UKB)
    assert c is None
    assert any("too few markers" in w for w in warns)


# --- the fallback, and why it is a fallback -------------------------------------------------


def test_the_copy_number_route_reproduces_its_own_real_world_failure():
    """On the real sample this route returned 0.19 where chrX returned 0.60.

    It is not independent of the channel it calibrates: the delivered calls are derived from the
    same LRR, so a systematically wrong caller yields a self-consistent and wrong scale.
    """
    obs = [(1, -0.192)] * 2_000 + [(3, 0.192 * math.log2(1.5))] * 4_000
    c, warns = nz.compression_from_copy_number(obs)
    assert c == pytest.approx(0.192, abs=0.01)
    assert any("NOT independent" in w for w in warns)
    assert any("credible floor" in w for w in warns), "0.19 is below the floor and must say so"


def test_route_preference_is_by_independence():
    auto = series(0.044, 0.4, 20_000)
    x = series(0.044 + MALE_RAW_SHIFT, 0.4, 2_000)
    # A CREDIBLE copy-number column, for the fallback case.
    cn_ok = [(1, -0.60)] * 2_000
    # The real one, which returned 0.192 - below the credible floor.
    cn_bad = [(1, -0.192)] * 2_000

    # A usable male chrX wins over a copy-number column, even when both are present.
    cal = nz.calibrate(auto, x, sex="male", product=UKB, copy_number_obs=cn_ok)
    assert cal.route == "chrx_male"
    assert 0.572 <= cal.lrr_compression <= 0.627

    # Without a male, a credible copy-number column is used AND flagged.
    fem_x = series(0.044 + FEMALE_RAW_SHIFT, 0.4, 2_000)
    cal2 = nz.calibrate(auto, fem_x, sex="female", product=UKB, copy_number_obs=cn_ok)
    assert cal2.route == "copy_number_column"
    assert any("NOT independent" in w for w in cal2.warnings)

    # But an INCREDIBLE one is rejected rather than used. On the real sample this route returned
    # 0.192 where chrX returned 0.598, so accepting it would have scaled every state call by a
    # third of the truth. Falling through to unscaled is the honest outcome.
    cal_bad = nz.calibrate(auto, fem_x, sex="female", product=UKB, copy_number_obs=cn_bad)
    assert cal_bad.route == "none"
    assert math.isnan(cal_bad.lrr_compression)
    assert any("also rejected" in w for w in cal_bad.warnings)

    # With neither, a declared value is used and flagged as unverified.
    cal3 = nz.calibrate(auto, fem_x, sex="female", product=UKB, declared_compression=0.60)
    assert cal3.route == "declared"
    assert any("unverified" in w for w in cal3.warnings)

    # With nothing at all the channel is centred but UNSCALED, and says so.
    cal4 = nz.calibrate(auto, fem_x, sex="female", product=UKB)
    assert cal4.route == "none"
    assert math.isnan(cal4.lrr_compression)
    assert any("must be reported as not tested" in w for w in cal4.warnings)


def test_centring_uses_the_samples_own_autosomal_median():
    """Zero must mean two copies IN THIS SAMPLE. The delivered channel arrives centred against
    somebody else's cluster file, and the two real papers' samples sit at medians from -0.009 to
    +0.060 - small, systematic, and enough to bias a state call."""
    auto = series(0.0597, 0.4, 20_000)     # egg donor A's measured autosomal median
    cal = nz.calibrate(auto, [], sex="female", product=UKB, declared_compression=0.6)
    assert cal.center == pytest.approx(0.0597, abs=1e-3)
    centred = nz.apply_center([0.0597, 0.0597 - 0.6, None], cal)
    assert centred[0] == pytest.approx(0.0, abs=1e-9)
    assert centred[1] == pytest.approx(-0.6, abs=1e-9)
    # A missing intensity stays missing: centring it to zero would turn absent data into a
    # confident claim of normal copy number.
    assert centred[2] is None


# --- the metric I got wrong ------------------------------------------------------------------


def test_separation_is_per_segment_not_per_marker():
    """The error worth pinning. c = 0.596 against sigma = 0.417 is 1.4 sigma per MARKER, and I
    reported that as 'this sample cannot make the call'. The HMM aggregates over a segment, where
    uncertainty on the mean falls as sigma/sqrt(n), so the same data gives 14 sigma over a hundred
    markers. Quoting the per-marker figure understates the available power by sqrt(n)."""
    cal = nz.Calibration(center=0.0, lrr_compression=0.596, sigma_lrr=0.417, route="chrx_male")
    assert cal.separation_sigma(1) == pytest.approx(1.43, abs=0.02)
    assert cal.separation_sigma(1) < 2.0, "per marker it genuinely does not separate"
    assert cal.separation_sigma(2) >= 2.0, "and two markers already do"
    assert cal.separation_sigma(100) == pytest.approx(14.3, abs=0.2)
    # sqrt scaling, not linear.
    assert cal.separation_sigma(400) == pytest.approx(2 * cal.separation_sigma(100), rel=1e-9)


def test_markers_needed_is_the_honest_resolution_statement():
    cal = nz.Calibration(center=0.0, lrr_compression=0.596, sigma_lrr=0.417, route="chrx_male")
    assert cal.markers_needed(2.0) == 2
    assert cal.markers_needed(3.0) == 5
    assert cal.markers_needed(5.0) == 13
    # A noisier sample needs a longer segment, monotonically.
    noisy = nz.Calibration(center=0.0, lrr_compression=0.596, sigma_lrr=1.2, route="chrx_male")
    assert noisy.markers_needed(3.0) > cal.markers_needed(3.0)
    # An uncalibrated sample cannot state a resolution at all.
    none = nz.Calibration(center=0.0, lrr_compression=math.nan, sigma_lrr=0.4, route="none")
    assert math.isnan(none.separation_sigma(100))


def test_sigma_comes_from_the_mad_and_survives_a_real_event():
    """A real segmental loss puts a block of markers far from the mode. A standard deviation
    absorbs them into the noise estimate and is then blind to them."""
    clean = series(0.0, 0.4, 2_000)
    contaminated = clean + [-0.6] * 200
    assert nz.robust_sigma(contaminated) == pytest.approx(nz.robust_sigma(clean), rel=0.30)
    assert math.isnan(nz.robust_sigma([0.1]))


def test_the_real_samples_measured_values_round_trip():
    """End to end on the two known males, from their measured medians."""
    for auto_med, x_med, spread in ((0.0440, -0.3635, 0.417), (-0.0093, -0.4125, 0.428)):
        cal = nz.calibrate(
            series(auto_med, spread, 20_000),
            series(x_med, spread, 2_000),
            sex="male", product=UKB,
        )
        assert cal.route == "chrx_male"
        assert 0.57 <= cal.lrr_compression <= 0.63, f"c was {cal.lrr_compression}"
        assert cal.center == pytest.approx(auto_med, abs=1e-3)
        # And with that calibration the distinction the whole layer exists for IS available.
        assert cal.separation_sigma(25) > 5.0
        assert cal.markers_needed(3.0) <= 8


# --- the two ceilings a real WGA'd embryo exposed -------------------------------------------


def test_compression_above_one_is_impossible_and_is_rejected():
    """A real WGA'd embryo (Zuccaro A8) gave c = 1.486 from its chrX. Compression is observed
    shift over THEORETICAL shift, so above 1 the array would have to amplify the copy-number
    signal rather than compress it. Every published value sits below 0.7. The ceiling was missing
    until end-to-end integration on real data produced the impossible number."""
    auto = series(0.0, 0.4, 20_000)
    # A chrX median far below the autosomes: noise, not one copy.
    x = series(-1.7, 0.4, 2_000)
    c, _, warns = nz.compression_from_chrx(auto, x, sex="male", product=UKB)
    assert c is not None and c > nz.MAX_CREDIBLE_COMPRESSION
    assert any("impossible" in w for w in warns)
    # And calibrate() must refuse it rather than pass it on.
    cal = nz.calibrate(auto, x, sex="male", product=UKB)
    assert cal.route == "none"
    assert math.isnan(cal.lrr_compression)


def test_an_unusably_noisy_channel_is_declared_unusable_not_fitted():
    """sigma 1.36 on that same embryo, against 0.42 on clean bulk DNA. A channel that noisy
    carries no copy-number information whatever scale is fitted to it, so copy-number states must
    be reported NOT TESTED rather than called from a number that happens to exist."""
    noisy = series(0.0, 3.0, 20_000)          # MAD well above the usable ceiling
    x = series(-0.6, 3.0, 2_000)
    cal = nz.calibrate(noisy, x, sex="male", product=UKB)
    assert cal.sigma_lrr > nz.MAX_USABLE_SIGMA_LRR
    assert cal.route == "none"
    assert any("NOT TESTED" in w for w in cal.warnings)
    # A clean channel at the same compression is still accepted.
    clean = nz.calibrate(series(0.0, 0.4, 20_000), series(-0.4075, 0.4, 2_000),
                         sex="male", product=UKB)
    assert clean.route == "chrx_male"
    assert clean.sigma_lrr < nz.MAX_USABLE_SIGMA_LRR

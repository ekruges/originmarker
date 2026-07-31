"""normalize - bring a sample's intensity channel onto a standard this tool can reason with.

The delivered `log2R` or `normalized_intensity` column was normalised by somebody else's
pipeline, against a cluster file built from samples that may bear no relation to this one. Track
4 of the research said so plainly: the .egt cluster file is load-bearing, and a sample processed
under conditions that differ from the ones used to build the clusters produces LRR that does not
mean what it says. So the channel arrives on an unknown scale with an unknown centre, and using
it as though it were calibrated is how a real distinction gets declared untestable.

This module re-derives the two things that matter, from the sample itself:

  centre  the autosomal median, subtracted, so 0 means two copies IN THIS SAMPLE
  scale   the LRR compression c, from an internal one-copy reference

The internal reference is chrX in a male. A male has one X, so the shift between his chrX and his
autosomes IS the one-copy shift, measured with no assumptions about copy number anywhere else.

That needs one correction, and it is not optional. Known FEMALES on the UK Biobank Axiom array
show chrX sitting +0.191 above their autosomes despite carrying two copies of each - probe design
and reference-panel composition, not biology. Ignoring it makes the male shift look like -0.407
and yields c = 0.407; correcting for it gives -0.598 and c = 0.598, which is where PennCNV's own
shipped defaults for CN probes sit (0.572-0.627). The uncorrected figure is wrong by a third.

Measured, not assumed: two known males (Zuccaro sperm donor, replicates) and two known females
(egg donors A and C), GSE148488.

Self-check:  python -m pytest tests/test_normalize.py
"""

from __future__ import annotations

import math
import statistics as st
from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, Sequence

#: chrX LRR offset at two copies, on the UK Biobank Axiom array, from known females.
#: A per-array constant. Any other array needs its own, and until one is measured the chrX
#: calibration route must be refused rather than run with this number.
CHRX_CN2_OFFSET = {"UK Biobank Axiom Array": 0.1908}

#: Below this the compression is not credible: it would mean a one-copy state shifts LRR by less
#: than a fifth of a log2 unit, which no published array does. Seeing it means the reference used
#: to derive it was wrong - a delivered copy_number column, usually.
MIN_CREDIBLE_COMPRESSION = 0.25

#: Above this it is not credible either, and this ceiling was missing until a real WGA'd embryo
#: produced c = 1.486. Compression is observed shift over THEORETICAL shift, so a value above 1
#: means the array amplified the copy-number signal rather than compressing it. No hybridisation
#: array does that; every published value sits below 0.7 (PennCNV 0.57-0.68, ASCAT 0.55). A value
#: over 1 means the reference was not a one-copy state - on that embryo, chrX was simply noisy.
MAX_CREDIBLE_COMPRESSION = 1.0

#: Above this per-marker noise the intensity channel carries no usable copy-number information,
#: whatever c says. Measured for scale: clean bulk DNA on this array runs 0.42, and the WGA'd
#: embryo that exposed the ceiling above ran 1.36 - three times the state separation it would need
#: to resolve. A channel that noisy has to be declared unusable rather than fitted.
MAX_USABLE_SIGMA_LRR = 0.80

CalibrationRoute = Literal["chrx_male", "copy_number_column", "declared", "none"]


def robust_sigma(values: Sequence[float]) -> float:
    """Scale from the median absolute deviation.

    MAD rather than the standard deviation because the residuals are contaminated by exactly the
    events being looked for: a real segmental loss puts a block of markers far from the mode, and
    a standard deviation would absorb them into the noise estimate and then be unable to see them.
    """
    vals = [v for v in values if math.isfinite(v)]
    if len(vals) < 2:
        return math.nan
    med = st.median(vals)
    return 1.4826 * st.median([abs(v - med) for v in vals])


@dataclass
class Calibration:
    """What the intensity channel means, for this sample."""

    #: Subtract this so that zero means two copies in this sample.
    center: float
    #: LRR compression: observed shift divided by the theoretical log2(CN/2).
    lrr_compression: float
    #: Per-marker noise, after centring.
    sigma_lrr: float
    route: CalibrationRoute
    warnings: list[str] = field(default_factory=list)
    #: Raw chrX-minus-autosome shift, before the female-baseline correction, for auditing.
    chrx_shift_raw: Optional[float] = None
    #: The measured chrX probe offset for this array, or None if this array has none on file.
    #: It is needed twice and the second use is easy to miss: once to DERIVE c from a male chrX,
    #: and again to INTERPRET any chrX intensity afterwards. Skipping the second use makes a
    #: male's single X read as two copies, because the offset is most of the shift being measured.
    chrx_offset: Optional[float] = None

    def separation_sigma(self, segment_markers: int) -> float:
        """How far apart copy-loss and copy-neutral are, over a segment of this many markers.

        The metric that matters, and not the per-marker one. Copy-loss (CN1) and copy-neutral
        (CN2 with LOH) differ by exactly c in expected LRR, and the uncertainty on a segment MEAN
        falls as sigma / sqrt(n). Quoting the per-marker figure understates the available power by
        a factor of sqrt(n) and will declare a perfectly callable distinction untestable: on real
        data c = 0.60 against sigma = 0.42 is 1.4 sigma per marker and 14 sigma over a hundred.
        """
        if segment_markers < 1 or not math.isfinite(self.sigma_lrr) or self.sigma_lrr <= 0:
            return math.nan
        return (self.lrr_compression / self.sigma_lrr) * math.sqrt(segment_markers)

    def markers_needed(self, target_sigma: float = 3.0) -> Optional[int]:
        """Smallest segment that reaches `target_sigma`. The honest resolution statement."""
        if not math.isfinite(self.sigma_lrr) or self.sigma_lrr <= 0:
            return None
        per = self.lrr_compression / self.sigma_lrr
        if per <= 0:
            return None
        return max(1, math.ceil((target_sigma / per) ** 2))


def compression_from_chrx(
    autosomal_lrr: Sequence[float],
    chrx_lrr: Sequence[float],
    *,
    sex: str,
    product: str,
) -> tuple[Optional[float], Optional[float], list[str]]:
    """Derive c from chrX in a male. Returns (c, raw_shift, warnings).

    Only a male gives a one-copy reference. A female's chrX is two copies, so her shift measures
    the array's own offset and nothing about compression - which is precisely why she is the
    control that makes this route trustworthy rather than a second guess.
    """
    warns: list[str] = []
    if sex != "male":
        return None, None, [
            f"chrX calibration needs a male sample; this one is {sex}. A female's chrX carries "
            "two copies, so her chrX-to-autosome shift measures the array's probe offset rather "
            "than the one-copy compression."
        ]
    if product not in CHRX_CN2_OFFSET:
        return None, None, [
            f"no measured chrX baseline offset for {product!r}. The offset is a per-array "
            "constant and using another array's would bias c by roughly a third, so the chrX "
            "route is refused rather than approximated."
        ]
    auto = [v for v in autosomal_lrr if math.isfinite(v)]
    x = [v for v in chrx_lrr if math.isfinite(v)]
    if len(auto) < 1_000 or len(x) < 200:
        return None, None, [
            f"too few markers for a stable median: {len(auto)} autosomal, {len(x)} chrX"
        ]

    raw = st.median(x) - st.median(auto)
    corrected = raw - CHRX_CN2_OFFSET[product]
    c = -corrected
    if c < MIN_CREDIBLE_COMPRESSION:
        warns.append(
            f"chrX-derived compression {c:.3f} is below the credible floor "
            f"{MIN_CREDIBLE_COMPRESSION}: a one-copy state cannot shift LRR that little. Either "
            "the sex call is wrong or chrX is itself aberrant in this sample."
        )
    elif c > MAX_CREDIBLE_COMPRESSION:
        warns.append(
            f"chrX-derived compression {c:.3f} exceeds 1.0, which is impossible: the observed "
            "shift cannot be larger than the theoretical one. chrX is not behaving as a one-copy "
            "reference in this sample - too noisy, itself aberrant, or the sex call is wrong."
        )
    return c, raw, warns


def compression_from_copy_number(
    observations: Iterable[tuple[int, float]],
) -> tuple[Optional[float], list[str]]:
    """Fallback: derive c from a delivered copy_number column.

    Weaker than it looks, and the warning is the point. On real Axiom output this route returned
    c = 0.19 where chrX calibration on the same sample returned 0.60, because the delivered calls
    placed half the genome at CN3 and only 15% at CN2 - self-consistent with its own LRR and still
    not a genome. A caller derived from the same channel it is calibrating cannot be independent
    of it.
    """
    by_cn: dict[int, list[float]] = {}
    for cn, lrr in observations:
        if cn <= 0 or cn == 2 or not math.isfinite(lrr):
            continue
        by_cn.setdefault(cn, []).append(lrr)
    if not by_cn:
        return None, ["no markers outside CN2, so nothing constrains the scale"]

    num = den = 0.0
    for cn, vals in by_cn.items():
        xk = math.log2(cn / 2)
        num += len(vals) * xk * st.fmean(vals)
        den += len(vals) * xk * xk
    if den == 0.0:
        return None, ["no scale constrained"]
    c = num / den
    warns = [
        "compression derived from the delivered copy_number column, which is NOT independent of "
        "the channel being calibrated. On real data this route disagreed with chrX calibration by "
        "a factor of three. Prefer a male chrX reference where one exists."
    ]
    if c < MIN_CREDIBLE_COMPRESSION:
        warns.append(
            f"and the result {c:.3f} is below the credible floor {MIN_CREDIBLE_COMPRESSION}, "
            "which is the signature of unreliable copy-number calls rather than a flat array"
        )
    return c, warns


def calibrate(
    autosomal_lrr: Sequence[float],
    chrx_lrr: Sequence[float],
    *,
    sex: str,
    product: str,
    copy_number_obs: Optional[Iterable[tuple[int, float]]] = None,
    declared_compression: Optional[float] = None,
) -> Calibration:
    """Normalise one sample's intensity channel, preferring the most independent reference.

    Route order is by independence, not convenience:

      1. chrX in a male, corrected by the female-derived baseline. Internal to the sample and
         independent of any copy-number call.
      2. the delivered copy_number column, flagged, because it is derived from the very channel
         being calibrated.
      3. a declared value, flagged as unverified.

    A sample with none of the three is returned uncalibrated with the scale left at NaN, so
    downstream code cannot mistake an unscaled channel for a scaled one.
    """
    auto = [v for v in autosomal_lrr if math.isfinite(v)]
    center = st.median(auto) if auto else 0.0
    sigma = robust_sigma([v - center for v in auto])

    c, raw, warns = compression_from_chrx(auto, chrx_lrr, sex=sex, product=product)
    route: CalibrationRoute = "chrx_male"

    def credible(x: Optional[float]) -> bool:
        return (x is not None and math.isfinite(x)
                and MIN_CREDIBLE_COMPRESSION <= x <= MAX_CREDIBLE_COMPRESSION)

    if not credible(c):
        if c is not None and math.isfinite(c):
            warns.append(f"chrX route rejected (c = {c:.3f}); falling back")
        if copy_number_obs is not None:
            c2, w2 = compression_from_copy_number(copy_number_obs)
            if credible(c2):
                c, route, warns = c2, "copy_number_column", warns + w2
            elif c2 is not None:
                warns += w2 + [f"copy_number route also rejected (c = {c2:.3f})"]
                c = None
            else:
                warns += w2
                c = None
        else:
            c = None
        if not credible(c) and declared_compression is not None:
            c, route = declared_compression, "declared"
            warns.append(
                f"using the declared compression {declared_compression:.3f}; it is unverified "
                "against this sample and no internal reference confirmed it"
            )
    # A channel this noisy cannot resolve a state whatever the scale says, so it is declared
    # unusable rather than fitted. Discovered on a real WGA'd embryo at sigma 1.36.
    too_noisy = math.isfinite(sigma) and sigma > MAX_USABLE_SIGMA_LRR
    if too_noisy:
        c, route = math.nan, "none"
        warns.append(
            f"per-marker intensity noise {sigma:.3f} exceeds the usable ceiling "
            f"{MAX_USABLE_SIGMA_LRR}: this channel carries no reliable copy-number information, "
            "so copy-number states must be reported as NOT TESTED regardless of any fitted scale."
        )
    if not credible(c) and not too_noisy:
        # Name which routes were unavailable versus tried and rejected. The distinction decides
        # what the operator can do about it: supply a male reference, or distrust the array's own
        # copy-number calls.
        tried = [r for r, avail in (("male chrX", sex == "male"),
                                    ("copy_number column", copy_number_obs is not None),
                                    ("declared value", declared_compression is not None)) if avail]
        c, route = math.nan, "none"
        warns.append(
            ("every available reference was rejected (" + ", ".join(tried) + ")"
             if tried else "no calibration reference available at all")
            + ". The intensity channel is CENTRED but UNSCALED, so copy-number states cannot be "
            "distinguished and must be reported as not tested."
        )

    return Calibration(center=center, lrr_compression=c, sigma_lrr=sigma, route=route,
                       warnings=warns, chrx_shift_raw=raw,
                       chrx_offset=CHRX_CN2_OFFSET.get(product))


def apply_center(values: Iterable[Optional[float]], cal: Calibration) -> list[Optional[float]]:
    """Re-centre a channel so that zero means two copies in this sample.

    None passes through as None: a missing intensity is missing, and centring it to zero would
    turn absent data into a confident statement of normal copy number.
    """
    return [None if v is None or not math.isfinite(v) else v - cal.center for v in values]

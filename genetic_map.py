"""
genetic_map - recombination context for candidate markers.

Interpolates cM locally from the bundled deCODE-derived GRCh38 map (provenance in
data/maps/README.md), so there is no network call at request time. Annotates each marker's
genetic distance from the pathogenic variant and flags hotspots between them.

Two cases are approximations rather than map readings, and both raise map_approx: a
chromosome with no bundled map file (chrY and chrM among them), and a position beyond
either end of the mapped span. Both are modelled at a uniform 1 cM/Mb.

The autosomal maps are sex-averaged and the chrX map is female-specific; neither describes
a particular carrier's meiosis. The MAP_SOURCE strings carry that caveat to the reader.
"""
from __future__ import annotations

import bisect
import gzip
import math
from array import array
from functools import lru_cache
from pathlib import Path
from typing import Optional

MAPS_DIR = Path(__file__).parent / "data" / "maps"

# Ships to the reader: stamped into export provenance and shown in the UI beside the cM
# and theta columns.
MAP_SOURCE = (
    "deCODE 2019 sex-averaged (Beagle GRCh38 liftover, plink format). Averaged over male "
    "and female meioses, so it is the map of neither. PGT-M linkage runs through one "
    "carrier of known sex: female recombination exceeds male genome-wide and the ratio "
    "varies by region, so cM and theta here can err in either direction for the carrier "
    "being tested."
)

# chrX is the FEMALE map, not a sex-averaged one: the bundled file runs ~1.30 cM/Mb
# against ~1.16-1.18 cM/Mb for the sex-averaged autosomes, and a sex-averaged X would have
# to sit below them (paternal meioses are a third of X transmissions and contribute no
# crossovers outside the PARs).
X_MAP_SOURCE = (
    "deCODE 2019 chrX FEMALE map (Beagle GRCh38 liftover, plink format). Female-specific, "
    "not sex-averaged like the autosomes, because the X does not recombine outside the "
    "pseudoautosomal regions in male meiosis. It applies to a female carrier. For a male "
    "carrier it overstates recombination: outside the PARs his X passes to each daughter "
    "as a single haplotype, so theta here does not describe his transmission."
)

APPROX_SOURCE = "uniform 1 cM/Mb APPROXIMATION (no map bundled for this chromosome)"

# Every chromosome's map stops short of the telomere, so a variant near an end lands here.
OFF_MAP_SOURCE = (
    "position beyond the end of the mapped span: cM EXTRAPOLATED at uniform 1 cM/Mb from "
    "the last map point, not a map reading"
)

# Rate used wherever the map cannot answer. Always labelled, never presented as a reading.
APPROX_CM_PER_MB = 1.0

HOTSPOT_CM_PER_MB = 10.0


class GeneticMap:
    """Interpolated cM lookup for one chromosome."""

    def __init__(self, chrom: str, pos: array, cm: array, approx: bool):
        self.chrom = chrom
        self._pos = pos
        self._cm = cm
        self.approx = approx

    @property
    def source(self) -> str:
        if self.approx:
            return APPROX_SOURCE
        return X_MAP_SOURCE if self.chrom == "X" else MAP_SOURCE

    def covers(self, position: int) -> bool:
        """True if position is inside the mapped span, i.e. cm_at() is a real reading."""
        return bool(len(self._pos)) and self._pos[0] <= position <= self._pos[-1]

    def cm_at(self, position: int) -> float:
        """Genetic position (cM) at a physical position, linearly interpolated.

        Beyond either end of the map the value is extrapolated at APPROX_CM_PER_MB rather
        than clamped to cm[0]/cm[-1]: clamping would report two distinct off-map positions
        as 0 cM apart. Callers must label any position where covers() is False.
        """
        if self.approx:
            return position / 1e6 * APPROX_CM_PER_MB

        pos, cm = self._pos, self._cm
        if position <= pos[0]:
            return cm[0] - (pos[0] - position) / 1e6 * APPROX_CM_PER_MB
        if position >= pos[-1]:
            return cm[-1] + (position - pos[-1]) / 1e6 * APPROX_CM_PER_MB

        # strictly inside the span now, so bisect lands on 1..len-1 and i-1 exists
        i = bisect.bisect_left(pos, position)
        if pos[i] == position:
            return cm[i]
        p0, p1 = pos[i - 1], pos[i]
        c0, c1 = cm[i - 1], cm[i]
        if p1 == p0:
            return c0
        return c0 + (c1 - c0) * (position - p0) / (p1 - p0)

    def cm_between(self, a: int, b: int) -> float:
        """Absolute genetic distance (cM) between two physical positions."""
        return abs(self.cm_at(b) - self.cm_at(a))

    def max_rate_between(self, a: int, b: int) -> float:
        """Peak local rate (cM/Mb) over the interval - used for hotspot detection.

        Scans the actual map points inside [a, b] rather than averaging the endpoints, so
        a narrow hotspot between two distant markers is not smoothed away.
        """
        lo, hi = (a, b) if a <= b else (b, a)
        if self.approx:
            return APPROX_CM_PER_MB

        pos, cm = self._pos, self._cm
        # cm_at models off-map stretches at 1 cM/Mb, so that is their floor rate here.
        peak = 0.0 if self.covers(lo) and self.covers(hi) else APPROX_CM_PER_MB
        # bisect_right on the low end, bisect_left on the high end: admits exactly the
        # intervals overlapping [lo, hi], never one past either side.
        i = max(bisect.bisect_right(pos, lo) - 1, 0)
        j = min(bisect.bisect_left(pos, hi), len(pos) - 1)
        for k in range(i, j):
            dp = pos[k + 1] - pos[k]
            if dp <= 0:
                continue
            rate = (cm[k + 1] - cm[k]) / (dp / 1e6)
            peak = max(peak, rate)
        return peak


def haldane_theta(cm: float) -> float:
    """Recombination fraction from genetic distance via Haldane's map function.

    theta = 0.5 * (1 - e^(-2d)) with d in Morgans. Haldane assumes no interference and so
    yields a slightly higher theta than Kosambi at the same cM: the conservative direction
    here, since it does not understate marker risk.
    """
    d = cm / 100.0
    return 0.5 * (1.0 - math.exp(-2.0 * d))


@lru_cache(maxsize=32)
def load(chrom: str) -> GeneticMap:
    """Load and cache one chromosome's map. Falls back to a labelled approximation."""
    c = str(chrom).replace("chr", "")
    f = MAPS_DIR / f"chr{c}.pos_cm.gz"
    if not f.exists():
        return GeneticMap(c, array("l"), array("d"), approx=True)

    pos, cm = array("l"), array("d")
    with gzip.open(f, "rt") as fh:
        for line in fh:
            a, _, b = line.partition(" ")
            if not b:
                continue
            pos.append(int(a))
            cm.append(float(b))
    if not pos:
        return GeneticMap(c, array("l"), array("d"), approx=True)
    return GeneticMap(c, pos, cm, approx=False)


def annotate_distance(chrom: str, variant_pos: int, marker_pos: int) -> dict:
    """cM distance, recombination fraction and hotspot flag between variant and marker."""
    gm = load(chrom)
    cm = gm.cm_between(variant_pos, marker_pos)
    peak = gm.max_rate_between(variant_pos, marker_pos)
    off_map = not (gm.covers(variant_pos) and gm.covers(marker_pos))
    return {
        "cm": round(cm, 5),
        # Significant figures, not decimal places: round() would flatten a small but real
        # theta to 0.0, and zero recombination is a claim this map cannot make.
        "recomb_fraction": float(f"{haldane_theta(cm):.3g}"),
        "hotspot_between": bool(peak >= HOTSPOT_CM_PER_MB),
        "peak_cm_per_mb": round(peak, 2),
        "map_approx": gm.approx or off_map,
        "map_source": OFF_MAP_SOURCE if off_map and not gm.approx else gm.source,
    }


if __name__ == "__main__":
    gm = load("11")
    assert not gm.approx, "chr11 map should be bundled"

    a = gm.cm_at(17_147_055)
    b = gm.cm_at(17_397_055)
    c = gm.cm_at(17_647_055)
    assert a < b < c, f"cM not monotonic: {a} {b} {c}"

    # ABCC8 locus: rs757110 sits 125 bp from the variant, so ~0 cM and ~0 theta.
    near = annotate_distance("11", 17_397_055, 17_396_930)
    assert near["cm"] < 0.001, near
    assert near["recomb_fraction"] < 0.001, near
    assert near["map_approx"] is False

    far = annotate_distance("11", 17_397_055, 17_606_055)  # ~209 kb out
    assert far["cm"] > near["cm"], (near, far)

    # Haldane: 0 cM -> 0; 50 cM -> ~0.316; asymptote at 0.5.
    assert haldane_theta(0) == 0.0
    assert abs(haldane_theta(50) - 0.31606) < 1e-4
    assert haldane_theta(10_000) > 0.4999

    fake = load("99")
    assert fake.approx and "APPROXIMATION" in fake.source
    assert abs(fake.cm_at(2_000_000) - 2.0) < 1e-9  # 1 cM/Mb

    # Past the last map point: must not collapse to a single cM value.
    end = gm._pos[-1]
    assert not gm.covers(end + 1), "covers() must exclude positions past the map"
    off = annotate_distance("11", end + 300_000, end + 500_000)
    assert off["cm"] > 0.0019, off                    # ~200 kb at 1 cM/Mb
    assert off["recomb_fraction"] > 0.0, off
    assert off["map_approx"] is True, off
    assert "EXTRAPOLATED" in off["map_source"], off
    assert off["peak_cm_per_mb"] > 0.0, off
    # Same below the first map point, and a straddling span is approx too.
    assert annotate_distance("11", 1_000, 50_000)["map_approx"] is True
    assert annotate_distance("11", gm._pos[0] - 1_000, 5_000_000)["map_approx"] is True
    assert gm.cm_at(1) < gm.cm_at(gm._pos[0]) < gm.cm_at(50_000_000)  # monotonic off-map too

    # A marker fully inside the map is still a real reading.
    assert near["map_approx"] is False and "deCODE" in near["map_source"]

    # max_rate_between must scan only the intervals overlapping the query.
    synth = GeneticMap("t", array("l", [100, 200, 300, 400]),
                       array("d", [0.0, 0.0, 0.0, 5.0]), approx=False)
    assert synth.max_rate_between(100, 200) == 0.0, "hotspot in [300,400] is not in [100,200]"
    assert synth.max_rate_between(300, 400) == 50_000.0, "hotspot in [300,400] must be found"
    below = GeneticMap("t", array("l", [100, 200, 300, 400]),
                       array("d", [0.0, 5.0, 5.0, 5.0]), approx=False)
    assert below.max_rate_between(200, 300) == 0.0, "hotspot in [100,200] is not in [200,300]"
    assert below.max_rate_between(150, 300) == 50_000.0, "a straddled hotspot must be found"

    # The sex caveat has to reach the reader: it ships in the provenance string.
    assert "carrier" in MAP_SOURCE and "sex-averaged" in MAP_SOURCE

    # chrX must outpace a sex-averaged autosome, or the bundled X has been scaled and the
    # FEMALE label on X_MAP_SOURCE is a lie.
    x, a7 = load("X"), load("7")
    rate = lambda m: m._cm[-1] / ((m._pos[-1] - m._pos[0]) / 1e6)
    assert rate(x) > rate(a7), f"chrX {rate(x):.3f} <= chr7 {rate(a7):.3f}: X may be scaled"
    assert x.source is X_MAP_SOURCE and "FEMALE" in x.source
    assert "male carrier" in x.source, "a male carrier's X does not recombine outside the PARs"
    assert load("11").source is MAP_SOURCE, "autosomes must not carry the chrX label"

    print(f"genetic_map self-check OK  |  chr11 rows={len(gm._pos)}  "
          f"cM@variant={b:.4f}  rs757110 theta={near['recomb_fraction']:.2e}  "
          f"chrX {rate(x):.3f} cM/Mb vs chr7 {rate(a7):.3f} cM/Mb (female, unscaled)")

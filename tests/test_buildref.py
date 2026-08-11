"""Assembly detection from marker positions, and the guard against the two tables drifting."""

from __future__ import annotations

import json
import re
from pathlib import Path

from originmarker import buildref
from originmarker import origin
ROOT = Path(__file__).resolve().parent.parent


def test_the_python_and_typescript_tables_are_the_same_tables():
    """`buildref.py` is transcoded from `buildref.ts`. Both are generated from UCSC and neither
    is hand-edited, so any divergence means one was regenerated without the other."""
    ts = (ROOT / "web" / "src" / "buildref.ts").read_text()

    def literal(name):
        i = ts.index(f"export const {name}")
        i = ts.index("{", ts.index("=", i))
        depth, j = 0, i
        while True:
            depth += (ts[j] == "{") - (ts[j] == "}")
            if depth == 0:
                break
            j += 1
        body = re.sub(r"\b(GRCh3[78])\s*:", r'"\1":', ts[i:j + 1]).replace("'", '"')
        return json.loads(re.sub(r",(\s*[}\]])", r"\1", body))

    assert literal("CHROM_LEN") == buildref.CHROM_LEN
    assert {b: {c: [list(x) for x in v] for c, v in d.items()}
            for b, d in literal("GAPS").items()} == buildref.GAPS
    assert sum(len(v) for v in buildref.GAPS["GRCh37"].values()) == 357
    assert sum(len(v) for v in buildref.GAPS["GRCh38"].values()) == 603


def sample_at(positions, chrom="1"):
    return origin.Sample("s", {f"m{i}": origin.Probe(chrom, p, "AA")
                               for i, p in enumerate(positions)}, [])


def test_a_position_past_the_end_of_a_chromosome_rules_that_build_out():
    """chr1 is 249,250,621 in GRCh37 and 248,956,422 in GRCh38, so a marker between the two can
    only be GRCh37. No rsIDs, no manifest, no chain file: positions alone."""
    between = 249_100_000  # past GRCh38 chr1, and in a legal GRCh37 stretch
    call = origin.detect_build(sample_at([between] * 1_500))
    assert call.build == "GRCh37"
    assert call.illegal["GRCh38"] == 1_500 and call.illegal["GRCh37"] == 0


def test_a_position_inside_an_assembly_gap_rules_that_build_out():
    """A marker cannot sit inside a run of Ns. GRCh37 chr1 has a gap at 121,535,434-124,535,434
    that GRCh38 does not place there."""
    call = origin.detect_build(sample_at([123_000_000] * 1_500))
    assert call.illegal["GRCh37"] == 1_500
    assert call.build == "GRCh38", call.note


def test_too_few_markers_is_undetermined_rather_than_a_coin_flip():
    """A single clean build is only a call once the alternative had a real chance to look clean
    too. Below the floor the honest answer is that positions have not separated them."""
    call = origin.detect_build(sample_at([249_100_000] * 50))
    assert call.build is None and "too few to call" in call.note
    assert call.tested == 50


def test_an_unlisted_assembly_is_reported_as_such_not_forced_into_one():
    """Real hg18 data: neither table is clean. Measured on the HapMap b36 export, 620 illegal
    under GRCh37 and 146 under GRCh38, and forcing the closer one would be a silent error."""
    call = origin.detect_build(sample_at([121_000_000, 123_000_000] * 750 + [249_100_000] * 100))
    assert call.build is None
    assert "no build is clean" in call.note


def test_markers_off_the_primary_chromosomes_are_not_counted():
    call = origin.detect_build(sample_at([1_000] * 100, chrom="GL000191.1"))
    assert call.tested == 0 and call.build is None
    assert "recognised primary chromosome" in call.note

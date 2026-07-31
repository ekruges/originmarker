"""The citable artifact: headline, figure, provenance, and a Methods paragraph from the run."""

from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import origin
import report


def plate(tmp_path, *, absent=0.002, het_band=0.02, nonpaternal=0.04, n=3_000):
    fp, sp = {}, {}
    for i in range(n):
        k, pos = f"1:{i}", 1_000 + i * 1_000
        fp[k] = origin.Probe("1", pos, "AB" if i % 5 == 0 else "AA", baf=0.5 if i % 5 == 0 else 0.0)
        scatter = (i * 7919) % 1000
        gt = "BB" if scatter < absent * 1000 else ("AB" if scatter >= 1000 - nonpaternal * 1000
                                                   else "AA")
        sp[k] = origin.Probe("1", pos, gt, baf=0.5 if scatter < het_band * 1000 else 0.0)
    f = tmp_path / "sperm.txt"
    f.write_text("Name,Chr,Position,GType\nrs1,1,100,AA\n")
    return origin.Sample("sperm", fp, []), origin.Sample("kid", sp, []), f


def build(tmp_path, **kw):
    fa, kid, f = plate(tmp_path, **kw)
    exp = origin.run_experiment(fa, [kid])
    pr = origin.parental_origin(fa, kid)
    tracks = {"kid": origin.segmental_origin(fa, kid, explainable=pr.explainable)}
    return report.build_report(exp, f, {"kid": f}, tracks=tracks)


def test_the_headline_states_the_verdict_and_the_number_that_produced_it(tmp_path):
    rep = build(tmp_path)
    h = rep.samples[0].headline
    assert h.state == "result"
    assert h.verdict in ("ANDROGENETIC", "BIPARENTAL", "GYNOGENETIC")
    assert any("paternal alleles absent" in f and "own noise could produce" in f for f in h.facts)
    assert any("informative markers" in f for f in h.facts)


def test_no_confidence_percentage_is_attached_to_the_verdict(tmp_path):
    """Deliberate. The statistical evidence runs to 10^11,660 because it multiplies across half a
    million markers, while the tool's own error rate is bounded by the handful of samples it has
    been checked against. A percentage would report the first as though it were the second."""
    rep = build(tmp_path)
    h = rep.samples[0].headline
    blob = " ".join([h.verdict, h.gloss] + h.facts).lower()
    for phrase in ("confidence", "% confident", "probability", "certainty"):
        assert phrase not in blob, phrase


def test_an_inconclusive_result_is_as_prominent_as_a_positive_one(tmp_path):
    """A report that whispers what it could not determine is the dangerous kind."""
    rep = build(tmp_path, absent=0.012, het_band=0.30, nonpaternal=0.30)
    h = rep.samples[0].headline
    assert h.state == "inconclusive" and h.verdict == "INCONCLUSIVE", h
    assert any("oocyte donor" in f for f in h.facts), "and what would settle it"


def test_provenance_carries_checksums_and_every_constant_with_its_source(tmp_path):
    rep = build(tmp_path)
    p = rep.provenance
    assert p.version and p.codename and p.tool == "OriginMarker"
    assert all(len(i["sha256"]) == 64 for i in p.inputs), "inputs are pinned by content"
    assert {i["role"] for i in p.inputs} == {"sperm", "sample"}
    names = {c["name"] for c in p.constants}
    assert "Mendelian inconsistency floor" in names
    assert "second-parent signal" in names
    assert all(len(c["provenance"]) > 30 for c in p.constants), "a value without a source is not"
    assert any("Kothiyal" in c["provenance"] for c in p.constants), "literature is cited"
    assert any("field-of-use" in c["provenance"] for c in p.constants), "licence limits stated"


def test_the_methods_paragraph_is_written_from_the_run_not_boilerplate(tmp_path):
    rep = build(tmp_path)
    m = rep.methods
    assert rep.provenance.version in m and rep.provenance.codename in m
    assert "no-call rate and heterozygous fraction" in m, "the derived bound is named"
    assert "half his heterozygosity" in m, "and so is the second axis"
    assert "chrX SNP rate rather than from chrY" in m
    assert "no liftOver was performed" in m
    # A clean run has no limits, so the section is absent rather than padded. One with a limit
    # names it: findings live in notes and must never appear here.
    assert "reported rather than resolved" not in m
    noisy = build(tmp_path, absent=0.012, het_band=0.30, nonpaternal=0.30)
    assert "reported rather than resolved" in noisy.methods
    assert "sits between what this sample" in noisy.methods
    assert "Every allele traces to" not in noisy.methods, "a finding is not a limitation"


def test_the_figure_is_valid_vector_with_the_thresholds_drawn(tmp_path):
    """Drawn rather than described, so the figure stands alone in a paper. Nothing is encoded in
    colour alone, so it survives greyscale printing."""
    rep = build(tmp_path, absent=0.002)
    sr = rep.samples[0]
    sr.track = report.Track(
        points=[report.TrackPoint("1", i * 1_000, i * 1_000 + 999, 0.002 if i < 40 else 0.30)
                for i in range(60)],
        noise_ceiling=0.0067, call_threshold=0.0201,
        called=[{"chrom": "1", "start": 40_000, "end": 59_999, "span_mb": 0.02,
                 "rate": 0.30, "margin": 15.0, "confident": True}])
    svg = report.track_svg(sr)
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert "noise ceiling" in svg and "calling threshold" in svg
    assert svg.count("<polyline") == 1, "one track line"
    assert 'fill-opacity="0.07"' in svg, "the called segment is shaded"
    assert "stroke-dasharray" in svg, "thresholds are dashed, not colour-coded"
    assert sr.sample_id in svg and sr.compared_with in svg


def test_a_uniform_genome_yields_a_figure_that_says_so_rather_than_a_blank(tmp_path):
    rep = build(tmp_path)
    sr = rep.samples[0]
    sr.track = None
    svg = report.track_svg(sr)
    assert ET.fromstring(svg).tag.endswith("svg")
    assert "no track" in svg


def test_every_number_is_in_the_json_so_a_reader_can_replot(tmp_path):
    rep = build(tmp_path)
    d = json.loads(rep.to_json())
    assert d["plate_summary"] and d["methods"]
    s = d["samples"][0]
    assert {"headline", "evidence", "track", "refusals"} <= set(s)
    assert len(s["evidence"]) == 3
    assert all({"axis", "observed", "reference", "basis"} == set(r) for r in s["evidence"])
    if s["track"]:
        assert s["track"]["points"] and "noise_ceiling" in s["track"]

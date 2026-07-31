"""report - a run turned into something that can be cited.

Everything a reader needs to check the result, and nothing that would let them mistake a refusal
for an answer. Three deliverables from one object:

    headline    the verdict, and equally prominent when there ISN'T one. A report that whispers
                what it could not determine is the dangerous kind.
    figure      one genome-wide track carrying the whole finding, as SVG: vector for print,
                native in a browser, and no plotting dependency to generate.
    provenance  version, input checksums, detected assembly, and every constant with where it
                came from, plus a Methods paragraph written from the actual run.

No percentage is attached to the verdict, deliberately. The statistical evidence for parent of
origin runs to 10^11,660 because it multiplies across half a million markers, while the tool's
own error rate is bounded by the samples it has been checked against, which is a far smaller
number. Printing a confidence would be reporting the first as though it were the second.

Self-check:  python -m report
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional, Sequence

import build_info
import origin

HeadlineState = Literal["result", "inconclusive", "not_determined"]

#: Every constant the parent-of-origin path can use, with where its value comes from. Printed in
#: the appendix so a reader can audit the ones that mattered rather than taking them on trust.
CONSTANTS: tuple[tuple[str, str, str], ...] = (
    ("Mendelian inconsistency floor", f"{origin.KOTHIYAL_FLOOR:.4f}",
     "Kothiyal 2019, 0.63% of variant sites across 1,314 nuclear families. A lower bound on the "
     "spurious-violation rate, never the value itself."),
    ("paternal-absence noise bound", "no-call rate x heterozygous fraction",
     "Derived per sample, not fitted. Dropout only manufactures paternal absence by turning a "
     "heterozygous call homozygous and discarding the paternal allele, so it scales with the "
     "product and neither term alone. Bounds every related pair measured, within 2x."),
    ("absence calling margin", f"{origin.ABSENCE_MARGIN:g}x",
     "How far above the noise bound an absence must sit before it is called."),
    ("second-parent signal", "half the father's heterozygosity",
     "Derived. A second parent supplies an allele the father lacks only where he is homozygous, "
     "at sum(pq); his heterozygosity is sum(2pq) over the same markers."),
    ("diploid heterozygous BAF band", f"{origin.HET_BAND_DIPLOID:.2f}",
     "Measured at 15-16% for diploid genomes and 1.3-3.4% for uniparental ones. Dropout does not "
     "mimic this: a biparental embryo at 33% dropout still reads 22.2%."),
    ("segment gap factor", f"{origin.GAP_FACTOR:g}x local spacing",
     "A run is split where consecutive markers lie further apart than this, because adjacency in "
     "the informative subsequence is not adjacency on the chromosome."),
    ("assembly tables", "UCSC chromInfo and gap, hg19 and hg38",
     "Freely redistributable. No liftOver chain file is used or shipped: those carry a "
     "non-commercial field-of-use restriction Apache 2.0 cannot sublicense."),
)


@dataclass
class Headline:
    state: HeadlineState
    verdict: str
    gloss: str
    facts: list[str] = field(default_factory=list)


@dataclass
class EvidenceRow:
    axis: str
    observed: str
    reference: str
    basis: str


@dataclass
class TrackPoint:
    chrom: str
    start: int
    end: int
    rate: float


@dataclass
class Track:
    points: list[TrackPoint]
    noise_ceiling: float
    call_threshold: float
    called: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class SampleReport:
    sample_id: str
    compared_with: str
    headline: Headline
    evidence: list[EvidenceRow]
    track: Optional[Track]
    refusals: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


@dataclass
class Provenance:
    tool: str
    version: str
    codename: str
    assembly: str
    inputs: list[dict[str, str]] = field(default_factory=list)
    constants: list[dict[str, str]] = field(default_factory=list)


@dataclass
class Report:
    plate_summary: str
    samples: list[SampleReport]
    provenance: Provenance
    methods: str

    def to_json(self, **kw: Any) -> str:
        """Every number in the report, so a reader can re-plot rather than trust a raster."""
        return json.dumps(asdict(self), indent=2, **kw)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _headline(sv: origin.SampleVerdict) -> Headline:
    p = sv.parentage
    GLOSS = {
        "androgenetic": "Every allele traces to the sperm donor. No maternal complement present.",
        "gynogenetic": "No paternal contribution anywhere. The genome is maternal in origin.",
        "biparental": "Both parents contributed. The paternal genome is present.",
    }
    if p.verdict == "unclear" or p.origin_class == "unclear":
        # An inconclusive result is stated as loudly as a positive one, with what would settle it.
        return Headline(
            "inconclusive", "INCONCLUSIVE",
            "The evidence does not separate the possibilities for this sample.",
            [f"paternal alleles absent {p.genome_rate:.2%}, against {p.explainable:.2%} this "
             f"sample's own noise could produce and {origin.ABSENCE_MARGIN:g}x that to call it",
             "the oocyte donor's array would measure dropout directly and settle it"],
        )
    facts = [
        f"paternal alleles absent {p.genome_rate:.2%}, against {p.explainable:.2%} this sample's "
        f"own noise could produce",
        f"{p.informative:,} informative markers where the sperm donor is homozygous",
    ]
    if p.sperm_type != "unknown":
        facts.append(f"{p.sperm_type.replace('_', '-')} sperm, from the chrX SNP rate rather "
                     f"than chrY")
    return Headline("result", p.origin_class.upper(),
                    GLOSS.get(p.origin_class, ""), facts)


def _evidence(sv: origin.SampleVerdict) -> list[EvidenceRow]:
    p = sv.parentage
    return [
        EvidenceRow("paternal alleles absent", f"{p.genome_rate:.2%}",
                    f"{p.explainable:.2%}",
                    "this sample's noise ceiling: no-call rate x heterozygous fraction"),
        EvidenceRow("alleles the father lacks", f"{p.nonpaternal_rate:.2%}",
                    f"{p.second_parent_expected:.2%}",
                    "half the father's heterozygosity, what a second parent would contribute"),
        EvidenceRow("heterozygous BAF band", f"{p.het_band:.2%}",
                    f"{origin.HET_BAND_DIPLOID:.0%}",
                    "diploid genomes run 15-16%, uniparental 1.3-3.4%"),
    ]


def build_report(
    experiment: origin.ExperimentReport,
    father_path: str | Path,
    sample_paths: dict[str, str | Path],
    *,
    tracks: Optional[dict[str, tuple[list, dict]]] = None,
) -> Report:
    """Turn a run into the citable object. `tracks` maps sample id to segmental_origin output."""
    tracks = tracks or {}
    samples: list[SampleReport] = []
    for sv in experiment.samples:
        tr = None
        segs, info = tracks.get(sv.sample_id, ([], {}))
        if info.get("track"):
            tr = Track(
                points=[TrackPoint(c, a, b, r) for c, a, b, r in info["track"]],
                noise_ceiling=info.get("explainable", math.nan),
                call_threshold=info.get("threshold", math.nan),
                called=[{"chrom": x.chrom, "start": x.start, "end": x.end,
                         "span_mb": round(x.span_mb, 3), "rate": x.rate,
                         "margin": x.margin, "confident": x.confident}
                        for x in segs if x.state == "paternal_absent"],
            )
        samples.append(SampleReport(
            sample_id=sv.sample_id, compared_with=experiment.father_id,
            headline=_headline(sv), evidence=_evidence(sv), track=tr,
            refusals=list(sv.parentage.limits), notes=list(sv.parentage.notes) + list(sv.notes),
        ))

    counts: dict[str, int] = {}
    for sr in samples:
        counts[sr.headline.verdict] = counts.get(sr.headline.verdict, 0) + 1
    summary = ", ".join(f"{v} {k.lower()}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1]))

    inputs = [{"role": "sperm", "file": Path(father_path).name, "sha256": _sha256(Path(father_path))}]
    for sid, p in sample_paths.items():
        inputs.append({"role": "sample", "id": sid, "file": Path(p).name,
                       "sha256": _sha256(Path(p))})

    prov = Provenance(
        tool="OriginMarker", version=build_info.VERSION, codename=build_info.CODENAME,
        assembly=experiment.build or "undetermined",
        inputs=inputs,
        constants=[{"name": n, "value": v, "provenance": s} for n, v, s in CONSTANTS],
    )
    return Report(summary, samples, prov, methods_paragraph(experiment, prov))


def methods_paragraph(experiment: origin.ExperimentReport, prov: Provenance) -> str:
    """A Methods section written from the run, including the refusals that actually fired."""
    n = len(experiment.samples)
    parts = [
        f"Parent of origin was determined from SNP array genotypes using {prov.tool} "
        f"{prov.version} ({prov.codename}). For each of {n} sample(s), the rate at which the "
        "sperm donor's obligate allele was absent was computed across autosomal markers where he "
        "is homozygous, and compared against the rate that sample's own genotyping noise can "
        "produce, taken as the product of its no-call rate and heterozygous fraction; a sample "
        f"was called as lacking a paternal contribution only where the observed rate exceeded "
        f"that bound by {origin.ABSENCE_MARGIN:g}-fold. Whether a second parent contributed was "
        "assessed from the rate at which the sample carries alleles the sperm donor does not "
        "possess, against half his heterozygosity, which is the contribution a second parent "
        "makes. Zygosity was taken from the fraction of B-allele frequencies in the heterozygous "
        "band. Sperm type was inferred from the chrX SNP rate rather than from chrY, so that "
        "X-bearing sperm are distinguishable from absent paternal contribution.",
        (f"Assembly was determined from marker positions against UCSC chromInfo and gap tables "
         f"and found to be {prov.assembly}; no liftOver was performed."
         if prov.assembly != "undetermined" else
         "Assembly could not be determined from marker positions against UCSC chromInfo and gap "
         "tables, which carry GRCh37 and GRCh38 only; no liftOver was performed and no genetic "
         "map was read at these coordinates."),
    ]
    fired = [r for sv in experiment.samples for r in sv.parentage.limits] + experiment.refusals
    if fired:
        parts.append("The following limits applied to this run and are reported rather than "
                     "resolved: " + " ".join(f"({i + 1}) {r.rstrip('.')}."
                                             for i, r in enumerate(dict.fromkeys(fired))))
    return "\n\n".join(parts)


# --- the figure ------------------------------------------------------------------------------


def track_svg(sr: SampleReport, *, width: int = 1000, height: int = 260) -> str:
    """Genome-wide paternal absence, as vector.

    The idiom a cytogeneticist already reads: position along the genome against a rate, with the
    thresholds drawn rather than described. Linear and auto-scaled, so a sample with one event
    shows a floor and one excursion, and a sample with no paternal genome shows everything above
    the line. Nothing is encoded in colour alone, so it survives greyscale printing.
    """
    t = sr.track
    if t is None or not t.points:
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="60" '
                f'role="img"><text x="8" y="34" font-family="sans-serif" font-size="13">'
                f'{sr.sample_id}: no track (genome uniform, or too few markers)</text></svg>')

    order = sorted({p.chrom for p in t.points}, key=lambda c: (len(c), c))
    spans = {c: max(p.end for p in t.points if p.chrom == c) for c in order}
    total = sum(spans.values()) or 1
    offset, acc = {}, 0
    for c in order:
        offset[c] = acc
        acc += spans[c]

    pad_l, pad_r, pad_t, pad_b = 52, 12, 14, 26
    plot_w, plot_h = width - pad_l - pad_r, height - pad_t - pad_b
    top = max([p.rate for p in t.points] + [t.call_threshold * 1.4, 0.01]) * 1.08

    def x(chrom: str, pos: int) -> float:
        return pad_l + (offset[chrom] + pos) / total * plot_w

    def y(rate: float) -> float:
        return pad_t + plot_h * (1.0 - min(rate, top) / top)

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="paternal allele absence across the genome for {sr.sample_id}">',
           '<style>text{font-family:ui-sans-serif,system-ui,sans-serif}</style>',
           f'<rect x="{pad_l}" y="{pad_t}" width="{plot_w}" height="{plot_h}" fill="none" '
           'stroke="#999" stroke-width="1"/>']

    for c in order:
        xs = x(c, spans[c])
        out.append(f'<line x1="{xs:.1f}" y1="{pad_t}" x2="{xs:.1f}" y2="{pad_t + plot_h}" '
                   'stroke="#e5e5e5" stroke-width="1"/>')
        mid = x(c, spans[c] // 2)
        out.append(f'<text x="{mid:.1f}" y="{height - 9}" font-size="9" fill="#444" '
                   f'text-anchor="middle">{c}</text>')

    for seg in t.called:
        x0, x1 = x(seg["chrom"], seg["start"]), x(seg["chrom"], seg["end"])
        out.append(f'<rect x="{x0:.1f}" y="{pad_t}" width="{max(1.0, x1 - x0):.1f}" '
                   f'height="{plot_h}" fill="#000" fill-opacity="0.07"/>')

    for rate, dash, label in ((t.noise_ceiling, "2,3", "noise ceiling"),
                              (t.call_threshold, "6,3", "calling threshold")):
        if not (isinstance(rate, float) and math.isfinite(rate)):
            continue
        yy = y(rate)
        out.append(f'<line x1="{pad_l}" y1="{yy:.1f}" x2="{pad_l + plot_w}" y2="{yy:.1f}" '
                   f'stroke="#000" stroke-width="1" stroke-dasharray="{dash}"/>')
        out.append(f'<text x="{pad_l + plot_w - 4:.1f}" y="{yy - 3:.1f}" font-size="9" '
                   f'fill="#000" text-anchor="end">{label} {rate:.2%}</text>')

    pts = " ".join(f'{x(p.chrom, (p.start + p.end) // 2):.1f},{y(p.rate):.1f}'
                   for p in sorted(t.points, key=lambda p: (offset[p.chrom], p.start)))
    out.append(f'<polyline points="{pts}" fill="none" stroke="#000" stroke-width="0.7" '
               'stroke-opacity="0.75"/>')

    for frac in (0.0, 0.5, 1.0):
        yy = pad_t + plot_h * (1 - frac)
        out.append(f'<text x="{pad_l - 6}" y="{yy + 3:.1f}" font-size="9" fill="#444" '
                   f'text-anchor="end">{top * frac:.1%}</text>')
    out.append(f'<text x="6" y="{pad_t + plot_h / 2:.0f}" font-size="9" fill="#444" '
               f'transform="rotate(-90 10 {pad_t + plot_h / 2:.0f})">paternal allele absent</text>')
    out.append(f'<text x="{pad_l}" y="{pad_t - 4}" font-size="10" fill="#000">'
               f'{sr.sample_id} vs {sr.compared_with}</text>')
    out.append("</svg>")
    return "\n".join(out)


if __name__ == "__main__":
    fa = origin.Sample("sperm", {f"1:{i}": origin.Probe("1", i * 1_000, "AA", baf=0.0)
                                 for i in range(3_000)}, [])
    kid = origin.Sample("kid", {f"1:{i}": origin.Probe("1", i * 1_000,
                                                       "BB" if 1_000 <= i < 1_400 else "AA",
                                                       baf=0.0)
                                for i in range(3_000)}, [])
    exp = origin.run_experiment(fa, [kid])
    segs, info = origin.segmental_origin(fa, kid, explainable=0.01)
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "sperm.txt"
        f.write_text("Name,Chr,Position,GType\nrs1,1,100,AA\n")
        rep = build_report(exp, f, {"kid": f}, tracks={"kid": (segs, info)})

    assert rep.samples and rep.samples[0].headline.verdict
    assert rep.provenance.inputs[0]["sha256"] and len(rep.provenance.inputs[0]["sha256"]) == 64
    assert rep.provenance.constants and all(c["provenance"] for c in rep.provenance.constants)
    assert "no liftOver was performed" in rep.methods
    json.loads(rep.to_json())
    svg = track_svg(rep.samples[0])
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "noise ceiling" in svg or "no track" in svg
    print("report self-check OK")

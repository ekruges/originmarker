"""Exports for a pb.PanelResult: CSV, JSON, XLSX, PDF.

Every export is self-describing: build, both variant forms (genomic-VCF and transcript
sense), source versions, timestamps, the Layer-B protocol and pb.DISCLAIMER verbatim.
Heterozygosity columns are named "..._prior": 2pq is a population prior, never a claim
about this carrier.

    to_csv/to_json/to_xlsx/to_pdf(result) -> bytes
    FILENAME(result, ext) -> 'originmarker_rs151344623_GRCh38_2026-07-15.csv'
"""
from __future__ import annotations

import csv
import io
import json
import re

import panelbuilder as pb

# gnomAD ancestry codes, in the engine's own order (AFR AMR ASJ EAS FIN NFE SAS MID).
POPS = list(pb.GNOMAD_POPS.values())

def _columns(anc=None) -> list:
    """Column names, parallel to _row(): keep the two in step, in order.

    Every quantity that decides in_recommended_panel must be a column, so the
    ancestry-matched 2pq joins the table whenever one was selected.
    """
    return (["rsid", "chrom", "pos_grch38", "ref", "alt", "signed_dist_bp", "side", "tier",
             "maf", "af", "gnomad_an", "het_2pq_prior_global", "het_2pq_prior_max_pop"]
            + ([f"het_2pq_prior_{anc}"] if anc else [])
            + [f"maf_{p}" for p in POPS]
            + [f"an_{p}" for p in POPS]
            + ["cm_to_variant", "recomb_fraction", "hotspot_between", "map_approx",
               "ensembl_pos_check", "in_recommended_panel"])


CSV_COLUMNS = tuple(_columns())            # the shape when no ancestry was selected


# --------------------------------------------------------------------------- #
# Shared facts. One source of truth for the header block / provenance sheet /
# PDF variant card, so the four formats cannot drift apart.
# --------------------------------------------------------------------------- #

def _ranking_key(result) -> str:
    """The engine's own name for the quantity that produced the sort order.

    Only panelbuilder knows the sort key, so exports render the name it gives and never
    restate the basis in their own words.
    """
    return result.provenance.get("ranking_key") or "not reported by this build"


def _rank_pop(result):
    """The ancestry the panel was ranked for, or None."""
    return result.provenance.get("ancestry_rank")


def _het(maf):
    """2pq from a minor allele frequency. The one place this formula is written."""
    return None if maf is None else round(2 * maf * (1 - maf), 4)


def _rank_het(m, anc):
    """The ancestry-matched 2pq prior for this marker, or None where that population has
    no gnomAD frequency here (the engine then falls back to the global 2pq)."""
    return _het(m.per_pop_maf.get(anc)) if anc else None


def _decider(m, anc):
    """The number that decided this marker's rank.

    Mirrors panelbuilder._rank_key across a seam: the engine exposes no per-marker key.
    The self-check pins the two together against the engine's real key.
    """
    h = _rank_het(m, anc)
    return m.het if h is None else h


def _map_approx(result) -> bool:
    """True if any cM value came from the 1 cM/Mb fallback rather than the map."""
    return any(m.map_approx for m in result.candidates)


def _cm_note(result) -> str:
    src = result.provenance["sources"]["genetic_map"]
    if _map_approx(result):
        return (f"cM/recomb_fraction are APPROXIMATE: some markers fall outside the "
                f"genetic map and use a 1 cM/Mb fallback (see map_approx column). "
                f"Map: {src}")
    return f"cM/recomb_fraction from genetic map: {src}"


def _genomic(v) -> str:
    return f"chr{v.chrom}:{v.pos_grch38} {v.vcf_ref}>{v.vcf_alt} ({v.build})"


def _facts(result) -> list[tuple[str, str]]:
    """Ordered (label, value) pairs describing identity, evidence and provenance."""
    v, rar, prov, cov = result.variant, result.rarity, result.provenance, result.coverage
    src = prov["sources"]
    strand = {1: "plus (+1)", -1: "minus (-1)"}.get(v.strand, "unknown")
    f = [
        ("Reference build", v.build),
        ("Query", v.query),
        ("rsID", v.rsid or "n/a"),
        ("Gene", v.gene or "n/a"),
        # R7: the two forms are always shown side by side and labelled.
        ("Genomic (VCF)", _genomic(v)),
        ("Transcript sense (HGVS c.)", v.transcript_sense_change()),
        ("Gene strand", strand),
        ("ClinVar significance", v.clinical_significance or "n/a"),
        ("ClinVar review status", v.review_status or "n/a"),
        ("ClinVar accession", v.clinvar_accession or "n/a"),
    ]
    if v.pos_grch37:
        # R6: carried for display only; every computation above is on pos_grch38.
        f.append(("GRCh37 position (display only)", f"chr{v.chrom}:{v.pos_grch37}"))
    if v.build_note:
        f.append(("Build note", v.build_note))
    f += [
        ("gnomAD genome AF", _num(rar.gnomad_af_genome)),
        ("gnomAD genome AC / AN", f"{_num(rar.gnomad_ac_genome)} / {_num(rar.gnomad_an_genome)}"),
        ("1000 Genomes AC", _num(rar.thousand_genomes_ac)),
        ("Population LD usable", str(rar.population_LD_usable)),
        ("LD verdict reason", rar.reason),
        ("Ranking key", _ranking_key(result)),
        ("Ranking exclusion (R2)",
         "LD with the pathogenic variant is never a ranking key."),
        ("Requested build", prov["requested_build"]),
        ("Window (bp, each side)", str(prov["window_bp"])),
        ("Common MAF floor", str(prov["common_maf"])),
        # Which ancestry was asked for, and nothing more: never name a 2pq here, the sort
        # key is _ranking_key's to report.
        ("Ancestry re-rank", prov["ancestry_rank"] or "none"),
        ("Candidates (common pool)", str(prov["candidate_n"])),
        ("Recommended panel size", str(len(result.recommended))),
        ("Coverage lower / higher coord", f"{cov['lower_count']} / {cov['higher_count']} markers; "
                               f"core-near {cov['lower_core_near']} / {cov['higher_core_near']}"),
        ("Coverage flags (R5)", "; ".join(cov["flags"]) or "none"),
        ("Source: ClinVar", src["clinvar"]),
        ("Source: Ensembl", src["ensembl"]),
        ("Source: gnomAD", src["gnomad"]),
        ("Source: genetic map", src["genetic_map"]),
        ("Genetic map note", _cm_note(result)),
        ("Source data as of (UTC)", prov["queried_utc"]),
        ("Panel built (UTC)", prov.get("built_utc") or prov["queried_utc"]),
        ("Source responses", f"{prov.get('source_responses_from_network', 0)} fetched live, "
                             f"{prov.get('source_responses_from_cache', 0)} from local cache"),
        ("Build elapsed (s)", str(prov["elapsed_s"])),
    ]
    return f


def _num(x) -> str:
    return "" if x is None else repr(x) if isinstance(x, float) else str(x)


def _row(m, recommended_ids: set, anc=None) -> list:
    """Values in _columns(anc) order. Keep the two in step."""
    return ([m.rsid, m.chrom, m.pos, m.ref, m.alt, m.dist, m.side, m.tier,
             m.maf, m.af, m.an, m.het, m.het_max_pop]
            + ([_rank_het(m, anc)] if anc else [])
            + [m.per_pop_maf.get(p) for p in POPS]
            + [m.per_pop_an.get(p) for p in POPS]
            + [m.cm, m.recomb_fraction, m.hotspot_between, m.map_approx,
               m.ensembl_pos_check, m.variant_id in recommended_ids])


def FILENAME(result, ext: str) -> str:
    v = result.variant
    stem = v.rsid or re.sub(r"\W+", "_", v.query).strip("_")[:40] or "variant"
    date = result.provenance["queried_utc"][:10]
    return f"originmarker_{stem}_{v.build}_{date}.{ext}"


# --------------------------------------------------------------------------- #
# CSV
# --------------------------------------------------------------------------- #

def to_csv(result) -> bytes:
    buf = io.StringIO(newline="")
    w = csv.writer(buf, lineterminator="\n")
    c = lambda s="": buf.write(f"# {s}\n".rstrip(" ") if s else "#\n")

    c("OriginMarker candidate linkage-marker panel - RESEARCH USE ONLY")
    c(pb.DISCLAIMER)                                   # R8, verbatim
    c()
    for k, val in _facts(result):
        c(f"{k}: {val}")
    c()
    c("CANDIDATE MARKERS ONLY (R3). These are not usable as-is:")
    for i, step in enumerate(pb.LAYER_B_STEPS, 1):     # R3, verbatim
        c(f"  {i}. {step}")
    c()
    c("het_2pq_prior_* are POPULATION PRIORS (R4): 2pq is the expected heterozygosity")
    c("in the population, NOT evidence that this carrier is heterozygous.")
    if _rank_pop(result):
        anc_ = _rank_pop(result)
        c(f"het_2pq_prior_{anc_} is the number the ranking keyed on. It is empty where gnomAD")
        c(f"reports no {anc_} frequency at that site (see an_{anc_}), and the ranking then")
        c("fell back to het_2pq_prior_global for that marker.")
    c("signed_dist_bp is relative to the pathogenic variant; negative = lower GRCh38 coordinate.")
    c("side is the genomic axis (lower / higher GRCh38 coordinate). It is NOT an arm label:")
    c("nothing here knows where the centromere is, so side never names a chromosome arm.")
    c(f"gnomad_an / an_* are allele NUMBERS: what each frequency rests on. Sites below")
    c(f"AN {pb.CALL_RATE_AN_FLOOR} are excluded entirely; a population MAF is reported only")
    c(f"above AN {pb.MIN_POP_AN} (~{pb.MIN_POP_AN // 2} people). gnomAD QC-failed sites")
    c("(AC0 / AS_VQSR) are excluded and never appear here.")
    c(f"All positions are {result.variant.build} (R6).")
    c(_cm_note(result))
    c()

    anc = _rank_pop(result)
    w.writerow(_columns(anc))
    rec = {m.variant_id for m in result.recommended}
    for m in result.candidates:
        w.writerow(_row(m, rec, anc))
    return buf.getvalue().encode("utf-8")


# --------------------------------------------------------------------------- #
# JSON
# --------------------------------------------------------------------------- #

def to_json(result) -> bytes:
    d = result.to_dict()
    rec = {m.variant_id for m in result.recommended}
    anc = _rank_pop(result)
    # zip relies on to_dict() being asdict(), which preserves list order; the assert checks it.
    for m, md in zip(result.candidates, d["candidates"]):
        assert md["variant_id"] == m.variant_id, "to_dict() reordered candidates"
        md["in_recommended_panel"] = md["variant_id"] in rec
        # Via _rank_het, the same expression the other three formats use: only _rank_het is
        # pinned to the engine's ranking key by the self-check.
        if anc:
            md[f"het_2pq_prior_{anc}"] = _rank_het(m, anc)
    prov = dict(d["provenance"], genetic_map_note=_cm_note(result),
                ranking_key=_ranking_key(result),
                ranking_excludes="LD with the pathogenic variant is never a ranking key (R2)",
                het_2pq_semantics="population prior, not a per-carrier genotype claim (R4)",
                transcript_sense=result.variant.transcript_sense_change())
    payload = {
        "provenance": prov,
        "variant": d["variant"],
        "rarity": d["rarity"],
        "coverage": d["coverage"],
        "candidates": d["candidates"],
        "recommended": d["recommended"],
        "layer_b_steps": list(pb.LAYER_B_STEPS),
        "disclaimer": pb.DISCLAIMER,
    }
    return json.dumps(payload, indent=2, ensure_ascii=False).encode("utf-8")


# --------------------------------------------------------------------------- #
# XLSX
# --------------------------------------------------------------------------- #

def to_xlsx(result) -> bytes:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font

    bold = Font(bold=True)
    wrap = Alignment(wrap_text=True, vertical="top")
    rec = {m.variant_id for m in result.recommended}
    anc = _rank_pop(result)

    wb = Workbook()
    wb.remove(wb.active)

    def sheet(title):
        """New sheet with pb.DISCLAIMER pinned in row 1 (R8, on EVERY sheet)."""
        ws = wb.create_sheet(title)
        ws["A1"] = pb.DISCLAIMER
        ws["A1"].font = Font(bold=True, italic=True)
        return ws

    def table(title, markers):
        ws = sheet(title)
        headers = _columns(anc)
        if _map_approx(result):
            headers[headers.index("cm_to_variant")] = "cm_to_variant (APPROX - see map_approx)"
        ws.append([])                       # row 2 spacer keeps A1 readable
        ws.append(headers)                  # row 3
        for c in ws[3]:
            c.font = bold
        for m in markers:
            ws.append(_row(m, rec, anc))
        ws.freeze_panes = "A4"              # rows 1-3: disclaimer + header stay visible
        _autosize(ws)
        return ws

    ws = table("Recommended panel", result.recommended)
    ws["A2"] = (f"Balanced subset, both sides (R5). CANDIDATES ONLY - genotype the carrier, "
                f"keep heterozygous markers, then phase (R3). het_2pq_prior_* are population "
                f"priors (R4). Positions {result.variant.build} (R6).")
    table("All candidates", result.candidates)["A2"] = (
        f"Full common pool (MAF >= {result.provenance['common_maf']}), ranked on "
        f"{_ranking_key(result)} - never on LD with the pathogenic variant (R2). "
        f"side is the genomic axis, not a chromosome arm: nothing here knows where the "
        f"centromere is. Positions {result.variant.build} (R6).")

    ws = sheet("Variant + provenance")
    ws.append([])
    ws.append(["Field", "Value"])
    for c in ws[3]:
        c.font = bold
    for k, v in _facts(result):
        ws.append([k, v])
    for row in ws.iter_rows(min_row=4, min_col=2, max_col=2):
        row[0].alignment = wrap
    ws.freeze_panes = "A4"
    _autosize(ws, cap=95)

    # Sheet name avoids "Layer B": that is this repo's word for the wet-lab half, not the
    # field's, and a reader opening the workbook has never seen the spec.
    ws = sheet("Using these markers")
    ws.append([])
    ws.append(["Step", "Action"])
    for c in ws[3]:
        c.font = bold
    for i, step in enumerate(pb.LAYER_B_STEPS, 1):
        ws.append([i, step])
    ws.append([])
    ws.append(["", "This app builds CANDIDATE panels only. It cannot determine phase, "
                   "and it makes no claim that any given carrier is heterozygous at any "
                   "marker - 2pq is a population prior (R3/R4)."])
    for row in ws.iter_rows(min_row=4, min_col=2, max_col=2):
        row[0].alignment = wrap
    ws.freeze_panes = "A4"
    _autosize(ws, cap=110)

    out = io.BytesIO()
    wb.save(out)
    return out.getvalue()


def _autosize(ws, cap=28):
    """Width from the longest cell, ignoring the pinned disclaimer/note rows."""
    for col in ws.columns:
        letter = col[0].column_letter
        widest = max((len(str(c.value)) for c in col[2:] if c.value is not None),
                     default=10)
        ws.column_dimensions[letter].width = min(max(widest + 2, 9), cap)


# --------------------------------------------------------------------------- #
# PDF
# --------------------------------------------------------------------------- #

def to_pdf(result) -> bytes:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfbase.pdfmetrics import stringWidth
    from reportlab.pdfgen import canvas as rl_canvas

    W, H = letter
    L, R, TOP, BOT = 40, W - 40, H - 42, 78      # margins; BOT clears the footer
    v = result.variant
    out = io.BytesIO()
    c = rl_canvas.Canvas(out, pagesize=letter)
    c.setTitle(FILENAME(result, "pdf"))
    state = {"y": TOP, "page": 0}

    def footer():
        """R8: pb.DISCLAIMER verbatim in the footer of EVERY page."""
        c.setFont("Helvetica-Oblique", 7)
        c.setFillColor(colors.black)
        yy = 44
        for line in _wrap(pb.DISCLAIMER, "Helvetica-Oblique", 7, R - L, stringWidth):
            c.drawString(L, yy, line)
            yy -= 8.5
        c.setFont("Helvetica", 6.5)
        c.setFillColor(colors.grey)
        c.drawString(L, 26, f"OriginMarker | {v.build} | queried {result.provenance['queried_utc']} "
                            f"| {result.provenance['sources']['gnomad']} | page {state['page']}")
        c.setFillColor(colors.black)

    def new_page():
        if state["page"]:
            footer()
            c.showPage()
        state["page"] += 1
        state["y"] = TOP

    def need(h):
        if state["y"] - h < BOT:
            new_page()

    def text(s, size=8, font="Helvetica", gap=2, colour=colors.black):
        need(size + gap)
        c.setFont(font, size)
        c.setFillColor(colour)
        for line in _wrap(s, font, size, R - L, stringWidth):
            need(size + gap)
            c.drawString(L, state["y"] - size, line)
            state["y"] -= size + gap
        c.setFillColor(colors.black)

    new_page()

    # --- variant card -------------------------------------------------------
    text("OriginMarker - candidate linkage markers for PGT-M", 14, "Helvetica-Bold", 4)
    text(f"{v.rsid or v.query}   {v.gene or ''}   {v.clinical_significance or ''}"
         f"   {v.clinvar_accession or ''}", 9, "Helvetica-Bold", 5)
    state["y"] -= 2
    card = [
        # R7: both forms, always labelled.
        ("Genomic (VCF):", _genomic(v)),
        ("Transcript sense (HGVS c.):", v.transcript_sense_change()),
        ("Query:", v.query),
        ("ClinVar review:", v.review_status or "n/a"),
        # fmt_af, not _num: this is prose, and the exact value sits beside it as AC/AN.
        # Machine-readable fields (the CSV data columns) must never be formatted this way.
        ("gnomAD genome AF:", f"{pb.fmt_af(result.rarity.gnomad_af_genome)}  "
                              f"(AC {_num(result.rarity.gnomad_ac_genome)} / "
                              f"AN {_num(result.rarity.gnomad_an_genome)})"),
        ("1000G AC:", _num(result.rarity.thousand_genomes_ac)),
        ("Population LD usable:", f"{result.rarity.population_LD_usable} - {result.rarity.reason}"),
        ("Window / MAF floor:", f"+/-{result.provenance['window_bp']} bp | MAF >= "
                                f"{result.provenance['common_maf']} | ancestry "
                                f"{_rank_pop(result) or 'none'}"),
        # tel_/cen_ are coverage-dict key names only: the label a reader sees names the
        # genomic axis, never a chromosome arm.
        ("Coverage (R5):", f"lower coord {result.coverage['lower_count']} / higher coord "
                           f"{result.coverage['higher_count']}"
                           f" | flags: {'; '.join(result.coverage['flags']) or 'none'}"),
    ]
    kw = 118
    for k, val in card:
        need(11)
        c.setFont("Helvetica-Bold", 8)
        c.drawString(L, state["y"] - 8, k)
        c.setFont("Helvetica", 8)
        for i, line in enumerate(_wrap(str(val), "Helvetica", 8, R - L - kw, stringWidth)):
            if i:
                need(10)
            c.drawString(L + kw, state["y"] - 8, line)
            state["y"] -= 10
    text(f"Ranked on {_ranking_key(result)}. LD with the pathogenic variant is never a "
         f"ranking key (R2).", 7, "Helvetica-Oblique", 3, colors.grey)

    # --- locus figure -------------------------------------------------------
    state["y"] -= 6
    _figure(c, result, L, R, state, need, colors)

    # --- recommended table --------------------------------------------------
    state["y"] -= 10
    text(f"Recommended panel ({len(result.recommended)} candidate markers, both sides)",
         10, "Helvetica-Bold", 4)
    approx = " *approx" if _map_approx(result) else ""
    anc = _rank_pop(result)
    # The head must name the population the deciding 2pq is for; a bare "2pq prior" on a
    # page stamped with an ancestry reads as the global figure. The global stays alongside,
    # since it is the fallback.
    dec_col = f"2pq {anc}" if anc else "2pq global"
    cols = ([("rsID", 58, "l"), ("pos " + v.build, 58, "r"), ("dist bp", 46, "r"),
             ("side", 42, "l"), ("tier", 52, "l"), ("2pq global", 42, "r")]
            + ([(dec_col, 46, "r")] if anc else [])
            + [("MAF", 32, "r"), ("cM" + approx, 36, "r"),
               ("hotspot", 34, "l"), ("pos chk", 72, "l")])
    assert sum(w for _, w, _ in cols) <= R - L, "PDF table wider than the page"

    def header_row():
        need(14)
        c.setFont("Helvetica-Bold", 6.8)
        x = L
        for name, w_, al in cols:
            c.drawString(x, state["y"] - 7, name)
            x += w_
        state["y"] -= 9
        c.setLineWidth(0.4)
        c.line(L, state["y"], R, state["y"])
        state["y"] -= 2

    header_row()
    for m in sorted(result.recommended, key=lambda m: m.dist):
        if state["y"] - 9 < BOT:
            new_page()
            text(f"Recommended panel (cont.) - {v.build}", 10, "Helvetica-Bold", 4)
            header_row()
        # Render m.side; never derive an arm from the sign of m.dist. Nothing here knows
        # where the centromere is, and below a q-arm locus is centromeric, not telomeric.
        rh = _rank_het(m, anc)
        vals = ([m.rsid, f"{m.pos:,}", f"{m.dist:+,}", m.side, m.tier, f"{m.het:.3f}"]
                + ([("-" if rh is None else f"{rh:.3f}")] if anc else [])
                + [f"{m.maf:.3f}",
                   "n/a" if m.cm is None else f"{m.cm:.4f}",
                   "yes" if m.hotspot_between else "no", m.ensembl_pos_check or "-"])
        c.setFont("Helvetica", 6.8)
        c.setFillColor(colors.red if m.ensembl_pos_check and
                       m.ensembl_pos_check.startswith("MISMATCH") else colors.black)
        x = L
        for (name, w_, al), val in zip(cols, vals):
            if al == "r":
                c.drawRightString(x + w_ - 6, state["y"] - 7, str(val))
            else:
                c.drawString(x, state["y"] - 7, str(val))
            x += w_
        c.setFillColor(colors.black)
        state["y"] -= 9
    note = (f"2pq = expected heterozygosity in the population (R4), NOT evidence that this "
            f"carrier is heterozygous. {dec_col} is the figure the ranking keyed on. "
            f"dist is signed from the variant; - = lower coordinate. side is the genomic "
            f"axis, not a chromosome arm: nothing here knows where the centromere is.")
    if anc:
        note += (f" A dash means gnomAD reports no {anc} frequency at that site, and the "
                 f"ranking fell back to the global 2pq for that marker.")
    if _map_approx(result):
        note += " *cM approximate: 1 cM/Mb fallback outside the map."
    text(note, 6.5, "Helvetica-Oblique", 2, colors.grey)

    # --- the wet-lab hand-off -----------------------------------------------
    state["y"] -= 8
    need(30)
    text("Using these markers: lab steps this app cannot do for you",
         10, "Helvetica-Bold", 4)
    for i, step in enumerate(pb.LAYER_B_STEPS, 1):     # R3, verbatim
        need(11)
        c.setFont("Helvetica-Bold", 7.6)
        c.drawString(L, state["y"] - 8, f"{i}.")
        c.setFont("Helvetica", 7.6)
        for j, line in enumerate(_wrap(step, "Helvetica", 7.6, R - L - 14, stringWidth)):
            if j:
                need(9.6)
            c.drawString(L + 14, state["y"] - 8, line)
            state["y"] -= 9.6
    text(f"Sources: ClinVar {result.provenance['sources']['clinvar']} | Ensembl "
         f"{result.provenance['sources']['ensembl']} | gnomAD "
         f"{result.provenance['sources']['gnomad']} | map "
         f"{result.provenance['sources']['genetic_map']} | queried "
         f"{result.provenance['queried_utc']} | build {v.build} (R6)",
         6.3, "Helvetica", 2, colors.grey)

    footer()
    c.showPage()
    c.save()
    return out.getvalue()


def _figure(c, result, L, R, state, need, colors):
    """Locus map: axis, dashed line at the variant, lollipops for recommended markers.

    Stem height comes from _decider, so the picture matches the table's deciding 2pq
    column. Colour = side.
    """
    # Bands, top-down: title, legend, stems, axis, ticks. axis_y is low enough that a
    # full-height stem (2pq=0.5) stops short of the legend row.
    FH = 104                                 # figure block height
    need(FH)
    top = state["y"]
    legend_y = top - 20
    axis_y = top - FH + 28
    stem_h = 42
    win = result.provenance["window_bp"]
    span = max([win] + [abs(m.dist) for m in result.recommended]) or 1
    mid = (L + R) / 2
    half = (R - L) / 2 - 14
    x_of = lambda d: mid + (d / span) * half

    anc = _rank_pop(result)
    dec_col = f"2pq {anc}" if anc else "2pq global"
    c.setFont("Helvetica-Bold", 8)
    c.setFillColor(colors.black)
    c.drawString(L, top - 8, f"Locus map ({result.variant.build}) - lollipop height = "
                             f"{dec_col} population prior (R4), colour = side")

    lo, hi = colors.HexColor("#2c6fbb"), colors.HexColor("#d1741f")
    # Keyed on m.side, not the sign of dist. An unknown side goes grey rather than
    # silently taking the other side's colour.
    side_col = {"lower coord": lo, "higher coord": hi}
    for m in sorted(result.recommended, key=lambda m: _decider(m, anc)):
        x = x_of(m.dist)
        h = (_decider(m, anc) / 0.5) * stem_h    # 2pq maxes out at 0.5
        col = side_col.get(m.side, colors.grey)
        c.setStrokeColor(col)
        c.setLineWidth(0.7)
        c.line(x, axis_y, x, axis_y + h)
        c.setFillColor(col)
        c.circle(x, axis_y + h, 2.1, stroke=0, fill=1)

    c.setStrokeColor(colors.black)
    c.setLineWidth(0.8)
    c.line(L, axis_y, R, axis_y)
    for d in (-span, -span / 2, 0, span / 2, span):
        x = x_of(d)
        c.setLineWidth(0.6)
        c.line(x, axis_y, x, axis_y - 3)
        c.setFont("Helvetica", 6)
        c.setFillColor(colors.grey)
        c.drawCentredString(x, axis_y - 11, f"{d/1000:+.0f} kb" if d else "0")

    c.setStrokeColor(colors.HexColor("#b12222"))
    c.setLineWidth(1.0)
    c.setDash(2, 2)
    c.line(x_of(0), axis_y - 5, x_of(0), axis_y + stem_h + 6)
    c.setDash()
    c.setFillColor(colors.HexColor("#b12222"))
    c.setFont("Helvetica-Bold", 6.5)
    c.drawCentredString(x_of(0), legend_y, result.variant.rsid or "pathogenic variant")

    c.setFont("Helvetica", 6)
    c.setFillColor(lo)
    c.drawString(L, legend_y, "lower GRCh38 coordinate")
    c.setFillColor(hi)
    c.drawRightString(R, legend_y, "higher GRCh38 coordinate")
    c.setFillColor(colors.grey)
    c.drawString(L, axis_y - 21, f"{dec_col} scale: 0 to 0.50 (max for a biallelic SNP)")
    c.setFillColor(colors.black)
    state["y"] = axis_y - 26


def _wrap(text: str, font: str, size: float, width: float, stringWidth) -> list[str]:
    lines, cur = [], ""
    for word in text.split():
        trial = f"{cur} {word}".strip()
        if stringWidth(trial, font, size) <= width or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = word
    return lines + [cur] if cur else lines or [""]


# --------------------------------------------------------------------------- #
# Self-check: PYTHONPATH=<repo root> PANELBUILDER_CACHE=tests/fixtures python app/exports.py
# --------------------------------------------------------------------------- #

def _pdf_text(data: bytes) -> str:
    """The strings a reader actually sees, pulled back out of the PDF.

    Not a general PDF parser: it reads the ASCII85+Flate streams and text-showing
    operators reportlab happens to emit, which is what drawString put there.
    """
    import base64
    import zlib
    out = []
    for s in re.findall(rb"stream\r?\n(.*?)~>\s*endstream", data, re.S):
        d = zlib.decompress(base64.a85decode(s))
        out += [m.decode("latin-1") for m in re.findall(rb"\((.*?)\) Tj", d)]
    return "\n".join(out)


def _pdf_stems(data: bytes) -> list[float]:
    """Heights of the locus figure's lollipop stems, read back off the page.

    Stems are the only vertical segments rising from the axis (ticks drop below it, the
    variant's dashed line starts 5pt lower), so they are the up-segments sharing the
    commonest bottom edge. Never assumes where axis_y is, so re-layout cannot break it.
    """
    import base64
    import collections
    import zlib
    segs = []
    for s in re.findall(rb"stream\r?\n(.*?)~>\s*endstream", data, re.S):
        d = zlib.decompress(base64.a85decode(s)).decode("latin-1")
        for x1, y1, x2, y2 in re.findall(r"([\d.]+) ([\d.]+) m\s+([\d.]+) ([\d.]+) l", d):
            if abs(float(x1) - float(x2)) < 0.01 and float(y2) > float(y1):
                segs.append((round(float(y1), 2), float(y2) - float(y1)))
    if not segs:
        return []
    axis_y = collections.Counter(y for y, _ in segs).most_common(1)[0][0]
    return [h for y, h in segs if y == axis_y]


if __name__ == "__main__":
    import os
    import pathlib
    import sys

    os.environ.setdefault("PANELBUILDER_CACHE",
                          str(pathlib.Path(__file__).resolve().parent.parent / "tests" / "fixtures"))
    r = pb.build("NM_000352.6(ABCC8):c.3989-9G>A")
    assert r.variant.rsid == "rs151344623", r.variant.rsid

    blobs = {ext: fn(r) for ext, fn in
             [("csv", to_csv), ("json", to_json), ("xlsx", to_xlsx), ("pdf", to_pdf)]}
    for ext, data in blobs.items():
        name = FILENAME(r, ext)
        assert data, f"{ext} is empty"
        pathlib.Path("/tmp", name).write_bytes(data)
        print(f"  /tmp/{name:52} {len(data):>9,} bytes")

    # R8: the disclaimer must survive verbatim into the text formats.
    for ext in ("csv", "json"):
        assert pb.DISCLAIMER.encode() in blobs[ext], f"disclaimer missing from {ext}"
    assert blobs["xlsx"][:2] == b"PK", "xlsx magic bytes"
    assert blobs["pdf"][:4] == b"%PDF", "pdf magic bytes"

    # R8 in the binaries too: xlsx zip entries are deflated, so read it back.
    from openpyxl import load_workbook
    wb = load_workbook(io.BytesIO(blobs["xlsx"]))
    assert wb.sheetnames == ["Recommended panel", "All candidates",
                             "Variant + provenance", "Using these markers"], wb.sheetnames
    for ws in wb:
        assert ws["A1"].value == pb.DISCLAIMER, f"disclaimer missing from sheet {ws.title!r}"
    assert wb["All candidates"].max_row == len(r.candidates) + 3, "candidate row count"

    # R4/R6/R7 spot checks on the CSV header + table.
    head, table = blobs["csv"].decode().split("rsid,chrom,pos_grch38", 1)
    for must in ("2pq", "prior", "GRCh38", "transcript sense", "minus strand",
                 pb.LAYER_B_STEPS[0], "C>T"):
        assert must in head, f"CSV header missing {must!r}"
    assert "het_2pq_prior_global" in ",".join(CSV_COLUMNS)
    assert table.count("\n") >= len(r.candidates)

    js = json.loads(blobs["json"])
    assert js["disclaimer"] == pb.DISCLAIMER
    assert js["layer_b_steps"] == pb.LAYER_B_STEPS
    assert js["provenance"]["build"] == "GRCh38"
    # Pinned to the result, not a literal: the candidate count moves legitimately whenever
    # upstream QC filters change.
    assert len(js["candidates"]) == len(r.candidates) > 1000
    assert sum(m["in_recommended_panel"] for m in js["candidates"]) == len(r.recommended)
    assert js["variant"]["pos_grch38"] == 17397055
    assert js["rarity"]["population_LD_usable"] is False

    # --- branches the golden case does not hit -------------------------------
    import copy

    # 1. Approximate genetic map => every format must say so next to the cM column.
    approx = copy.deepcopy(r)
    approx.candidates[0].map_approx = True
    assert "APPROXIMATE" in _cm_note(approx)
    assert b"APPROXIMATE" in to_csv(approx)
    assert "APPROX" in json.loads(to_json(approx))["provenance"]["genetic_map_note"].upper()
    wb2 = load_workbook(io.BytesIO(to_xlsx(approx)))
    assert any("APPROX" in str(c.value) for c in wb2["All candidates"][3]), "xlsx cM note"
    assert to_pdf(approx)[:4] == b"%PDF"

    # 2. No rsID (bare-HGVS resolve) => filename falls back to the sanitised query.
    anon = copy.deepcopy(r)
    anon.variant.rsid = None
    # Pin the stem, not the calendar: the date comes from the result's own provenance.
    fn = FILENAME(anon, "csv")
    assert fn.startswith("originmarker_NM_000352_6_ABCC8_c_3989_9G_A_GRCh38_"), fn
    assert fn.endswith(".csv") and fn[-14:-4] == r.provenance["queried_utc"][:10], fn
    assert to_pdf(anon)[:4] == b"%PDF"      # figure labels the variant without an rsID

    # 3. A panel too tall for one page must paginate, not overflow.
    big = copy.deepcopy(r)
    big.recommended = (big.recommended * 6)[:110]
    pdf = to_pdf(big)
    assert pdf.count(b"/Type /Page\n") >= 2 or pdf.count(b"/Type/Page") >= 2, "expected 2+ pages"
    assert len(pdf) > len(blobs["pdf"])

    # 4. The PDF says what the MODEL says about which side a marker is on, and never
    # names a chromosome arm.
    pdf_text = _pdf_text(blobs["pdf"])
    assert "rs757110" in pdf_text and "Recommended panel" in pdf_text, "pdf text extraction"
    for m in r.recommended:
        assert m.side in pdf_text, f"PDF does not render m.side {m.side!r} for {m.rsid}"
    for fmt in ("pdf", "csv"):
        body = pdf_text if fmt == "pdf" else blobs["csv"].decode()
        arm = re.findall(r"\b(?:tel|cen)\b", body)
        assert not arm, f"{fmt} still claims a chromosome arm: {arm[:4]}"

    # ...and renders that field rather than deriving it. The engine sets side from the sign
    # of dist, so only a marker whose side CONTRADICTS its sign can tell the two apart.
    q_arm = copy.deepcopy(r)
    q_arm.recommended[0].side = "SIDE_FROM_MODEL"
    assert q_arm.recommended[0].dist < 0
    assert "SIDE_FROM_MODEL" in _pdf_text(to_pdf(q_arm)), "PDF recomputes side, not renders it"
    assert "SIDE_FROM_MODEL" in to_csv(q_arm).decode()

    # 5. Every quantity that decides in_recommended_panel is a column in every export.
    # With an ancestry the sort keys on that population's 2pq, so it must be readable
    # alongside the global one.
    eas = pb.build(pb.StructuredQuery(variant="NM_000352.6(ABCC8):c.3989-9G>A", ancestry="EAS"))
    assert _rank_pop(eas) == "EAS"
    m0 = max(eas.recommended, key=lambda m: abs((_rank_het(m, "EAS") or 0) - m.het))
    assert _rank_het(m0, "EAS") is not None and abs(_rank_het(m0, "EAS") - m0.het) > 0.05, \
        "need a marker whose EAS 2pq differs from global, or this check proves nothing"
    eas_pdf = _pdf_text(to_pdf(eas))
    assert "2pq EAS" in eas_pdf, "PDF must head the deciding column with the ancestry"
    assert f"{_rank_het(m0, 'EAS'):.3f}" in eas_pdf, "PDF must print the EAS 2pq it ranked on"
    # Read the deciding column back BY NAME and check the VALUE under it, in each format:
    # presence-only and length-only asserts both pass while the artifact is wrong, since
    # _columns()/_row() are parallel lists that can drift out of order at equal length.
    eas_csv = to_csv(eas).decode().split("\n")
    hdr_i = next(i for i, l in enumerate(eas_csv) if l.startswith("rsid,"))   # find it, don't count back
    eas_cols = eas_csv[hdr_i].split(",")
    eas_xl = load_workbook(io.BytesIO(to_xlsx(eas)))["All candidates"]
    xl_hdr = [c.value for c in eas_xl[3]]
    # Keyed by variant_id, NOT rsid: multi-allelic sites share an rsID, so an rsid-keyed
    # dict drops every alt but the last and compares the wrong row.
    eas_json = {c["variant_id"]: c for c in json.loads(to_json(eas))["candidates"]}
    for fmt, cols in (("csv", eas_cols), ("xlsx", xl_hdr)):
        assert "het_2pq_prior_EAS" in cols, f"{fmt} has no het_2pq_prior_EAS column"
    for m, csv_line, xl_row in zip(eas.candidates, eas_csv[hdr_i + 1:],
                                   eas_xl.iter_rows(min_row=4, values_only=True)):
        want = _rank_het(m, "EAS")
        by_name = dict(zip(_columns("EAS"), _row(m, set(), "EAS")))
        assert by_name["het_2pq_prior_EAS"] == want, "_columns/_row are out of order"
        assert by_name["rsid"] == m.rsid, "_columns/_row are out of order"
        got_csv = csv_line.split(",")[eas_cols.index("het_2pq_prior_EAS")]
        assert got_csv == ("" if want is None else str(want)), \
            f"csv: het_2pq_prior_EAS is {got_csv!r} under its own header, want {want!r}"
        assert xl_row[xl_hdr.index("het_2pq_prior_EAS")] == want, "xlsx: wrong value under header"
        assert eas_json[m.variant_id]["het_2pq_prior_EAS"] == want, \
            f"json: het_2pq_prior_EAS is {eas_json[m.variant_id]['het_2pq_prior_EAS']!r}, want {want!r}"
    assert len(_columns("EAS")) == len(_row(m0, set(), "EAS")) == len(CSV_COLUMNS) + 1
    # No ancestry: the decider is the global 2pq, which is on the page under its own name.
    assert "2pq global" in pdf_text and "2pq prior" not in pdf_text

    # _decider mirrors panelbuilder._rank_key across a seam: pin it to the engine's actual
    # key so a change of ranking quantity fails here, not in a shipped export.
    for res, a in ((r, None), (eas, "EAS")):
        for m in res.recommended:
            engine_het = -pb._rank_key(a)(m)[0]
            # Tolerance covers _het's round-to-4dp and nothing more. NOT 5e-5: the engine
            # keeps 2p(1-p) unrounded, and ordinary MAFs land on that exact boundary just
            # over a strict < 5e-5. A change of ranking quantity is orders bigger.
            assert abs(engine_het - _decider(m, a)) < 1e-4, \
                f"exports and the engine disagree on what ranked {m.rsid}: " \
                f"{_decider(m, a)} vs {engine_het}"

    # 6. Seam contract: exports render provenance["ranking_key"] verbatim and never
    # restate the ranking basis in their own words.
    keyed = copy.deepcopy(r)
    keyed.provenance["ranking_key"] = "SENTINEL-KEY 2pq in XYZ"
    assert _ranking_key(keyed) == "SENTINEL-KEY 2pq in XYZ"
    assert "SENTINEL-KEY 2pq in XYZ" in _pdf_text(to_pdf(keyed))
    assert "SENTINEL-KEY 2pq in XYZ" in to_csv(keyed).decode()
    assert json.loads(to_json(keyed))["provenance"]["ranking_key"] == "SENTINEL-KEY 2pq in XYZ"
    for body in (pdf_text, blobs["csv"].decode()):
        assert "none (global 2pq prior)" not in body, "exports restating the ranking basis"
    if not r.provenance.get("ranking_key"):
        print("  NOTE: panelbuilder sets no provenance['ranking_key'] yet; exports render "
              "the marked fallback.", file=sys.stderr)

    # 7. The FIGURE draws the number it ranked on: every check above reads text, none the
    # geometry. Scale-free on purpose, asserting heights are PROPORTIONAL to _decider
    # rather than equal to a layout constant, so re-sizing the figure does not fail this.
    dec = sorted(_decider(m, "EAS") for m in eas.recommended)
    stems = sorted(_pdf_stems(to_pdf(eas)))
    assert len(stems) == len(dec), f"{len(stems)} lollipops for {len(dec)} markers"
    ratios = [h / d for h, d in zip(stems, dec)]
    assert max(ratios) - min(ratios) < 1.0, \
        f"lollipop heights are not proportional to the 2pq that ranked them: {ratios[:3]}"

    # 8. The no-frequency fallback: "-" in the column, global 2pq in the ranking. The
    # marker is synthesised because no fixture in this repo reaches this branch: ABCC8 has
    # a gnomAD frequency for every recommended marker in all 8 populations.
    gap = copy.deepcopy(eas)
    m_gap = gap.recommended[0]
    del m_gap.per_pop_maf["EAS"]
    assert _rank_het(m_gap, "EAS") is None
    assert _decider(m_gap, "EAS") == m_gap.het, "fallback must be the global 2pq, as the engine does"
    assert -pb._rank_key("EAS")(m_gap)[0] == m_gap.het, "engine disagrees about the fallback"
    gap_row = [l for l in _pdf_text(to_pdf(gap)).split("\n")]
    i = gap_row.index(m_gap.rsid)
    assert gap_row[i + 6] == "-", f"no-frequency marker must print '-', printed {gap_row[i + 6]!r}"
    assert "no EAS frequency" in _pdf_text(to_pdf(gap)), "PDF must footnote what the dash means"
    assert _row(m_gap, set(), "EAS")[13] is None, "csv/xlsx cell must be empty, never a number"

    # An allele frequency is prose on the PDF and a data cell in the CSV: the PDF renders
    # it via fmt_af, the CSV column keeps the exact float.
    pdf_txt = _pdf_text(to_pdf(r))
    assert repr(r.rarity.gnomad_af_genome) not in pdf_txt, "float repr leaked into the PDF"
    assert pb.fmt_af(r.rarity.gnomad_af_genome) in pdf_txt, "PDF should render the AF via fmt_af"
    af_col = _columns(None).index("af")
    assert _row(r.recommended[0], set())[af_col] == r.recommended[0].af, \
        "the CSV af column must stay an exact float, never a formatted string"

    print(f"\nself-check OK: {len(r.candidates)} candidates, {len(r.recommended)} recommended, "
          f"4/4 formats, disclaimer verbatim in all; map-approx/no-rsid/multipage branches OK; "
          f"PDF renders m.side and the deciding 2pq, and draws it at the right height.",
          file=sys.stderr)

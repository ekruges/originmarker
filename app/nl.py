"""Natural-language INTENT parser for OriginMarker (R1).

MAY DECIDE: which variant the user named, and the search knobs they asked for (window,
MAF floor, ancestry).

MAY NEVER DECIDE: a coordinate. No chromosome, position, strand or ref/alt allele, from
regex, LLM or memory. `variant` leaves this module as an opaque identifier that
pb.resolve_variant() looks up live; every coordinate in the app is produced there.

The boundary is enforced twice: pb.StructuredQuery has no coordinate field, so this
layer has nowhere to put one, and _reject_coordinates() screens whatever the model
returns. The system prompt asks the model not to invent coordinates, but the prompt is
not the control: the regex is.

The regex fast path handles any text containing an rsID or HGVS for zero tokens.
Modifiers are always parsed locally, even on the LLM path. The LLM is a last resort for
prose that names no identifier at all.

Self-check:  PYTHONPATH=<repo root> python app/nl.py
"""

from __future__ import annotations

import json
import logging
import os
import re
from typing import Optional, Tuple

import panelbuilder as pb

MODEL = "claude-haiku-4-5-20251001"
MAX_TOKENS = 200

_log = logging.getLogger(__name__)

# Every dead end on this path has the same way out, so it is written once.
_NAME_IT = ("Name the variant directly, for example rs151344623 or "
            "NM_000352.6(ABCC8):c.3989-9G>A.")


# --------------------------------------------------------------------------- #
# Shapes: what an identifier looks like, and what a coordinate looks like.
# --------------------------------------------------------------------------- #

_RSID = re.compile(r"\brs\d+\b", re.I)

_HGVS = re.compile(
    r"\b(?:N[MCGRP]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?)"   # accession
    r"\s*(?:\(([^()]{1,32})\))?"                            # optional (GENE)
    r"\s*:\s*[cgmnpr]\.\S+",                                # :c.3989-9G>A
    re.I,
)

_CHROM = r"(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT?)"

# Only for a friendlier error message, NOT the control. The control is _ID_SHAPES.
_LOOKS_COORDINATE = (
    re.compile(rf"{_CHROM}[-:_]\d+[-:_][ACGTN]+[-:_][ACGTN]+", re.I),   # 11-17397055-C-T
    re.compile(rf"{_CHROM}\s*[:\-]\s*[\d,]{{4,}}", re.I),               # chr11:17,397,055
    re.compile(r"^[\d,]{4,}$"),                                         # bare 17397055
    re.compile(r"(?:GRCh3[78]|hg\d+)", re.I),                           # GRCh38:11:17397055
)

# An ALLOW-LIST of identifier shapes, and the R1 tripwire. It must stay an allow-list: a
# blacklist of coordinate formats can never be complete, whereas anything that is not a
# recognised identifier is refused here by construction.
_ID_SHAPES = (
    _RSID,                                                     # rs151344623
    _HGVS,                                                     # NM_000352.6(ABCC8):c.3989-9G>A
    re.compile(r"(?:VCV|RCV)\d+(?:\.\d+)?", re.I),             # ClinVar accession
    re.compile(r"[A-Z][A-Z0-9-]{1,14}\s*[: ]\s*[cp]\.\S+", re.I),  # ABCC8:c.3989-9G>A / ABCC8 p.R1215Q
)


def _reject_coordinates(s: str) -> str:
    """Return `s` only if it is a recognised variant IDENTIFIER; raise otherwise (R1).

    The enforcement point for R1 on LLM output: a coordinate must never reach
    StructuredQuery, where it would look like a resolved locus no API confirmed. Raises
    rather than dropping, since a coordinate here means the intent layer misread its job.

    Accession-anchored HGVS ('NC_000011.10:g.17397055C>T') is an identifier, not a
    coordinate: resolve_variant() still looks it up, so the number is never trusted.
    """
    t = (s or "").strip().strip("`'\"")
    if any(shape.fullmatch(t) for shape in _ID_SHAPES):
        return t
    # These messages reach the user through /api/nl's 400 detail: say what to type
    # instead, and never cite internal rule numbering.
    if any(shape.search(t) for shape in _LOOKS_COORDINATE):
        raise ValueError(
            f"{s!r} looks like a genomic position. Enter an rsID or an HGVS expression "
            "instead, and the position will be looked up."
        )
    raise ValueError(
        f"{s!r} is not a recognised variant identifier. Expected an rsID, an HGVS "
        "expression, or a ClinVar accession."
    )


def _require_identifier(s: str) -> str:
    v = _reject_coordinates(s)
    if not (_RSID.fullmatch(v) or _HGVS.fullmatch(v)):
        raise ValueError(
            f"{s!r} is not an rsID or HGVS identifier. Name the variant the way it "
            "appears in the report (e.g. rs151344623 or NM_000352.6(ABCC8):c.3989-9G>A)."
        )
    return v


# --------------------------------------------------------------------------- #
# Local (free) extraction
# --------------------------------------------------------------------------- #

_UNITS = {"bp": 1, "kb": 1_000, "mb": 1_000_000}
_WINDOW = re.compile(r"(\d+(?:\.\d+)?)\s*(kb|mb|bp)\b", re.I)
_MAF = re.compile(
    r"\bmaf\b\s*(?:>=|<=|>|<|=|of|:|at\s+least|above|over)?\s*(\d+(?:\.\d+)?)\s*(%?)",
    re.I,
)

# Bare gnomAD codes, case-SENSITIVE and sourced from the engine so the two can't
# drift. Case matters: it keeps \bMID\b off "middle" and \bFIN\b off "finding".
_CODES = re.compile(r"\b(" + "|".join(sorted(pb.GNOMAD_POPS.values())) + r")\b")

# Ordered: more specific phrases first, so "non-Finnish European" beats "Finnish".
# No trailing \b - these must survive the plural ("Europeans", "South Asians").
_ANCESTRY_WORDS = (
    (re.compile(r"\bnon[-\s]?finnish\s+european", re.I), "NFE"),
    (re.compile(r"\bashkenazi", re.I), "ASJ"),
    (re.compile(r"\bfinnish", re.I), "FIN"),
    (re.compile(r"\beuropean", re.I), "NFE"),
    (re.compile(r"\bafrican", re.I), "AFR"),
    (re.compile(r"\beast[-\s]asian", re.I), "EAS"),
    (re.compile(r"\bsouth[-\s]asian", re.I), "SAS"),
    (re.compile(r"\bmiddle[-\s]eastern", re.I), "MID"),
    (re.compile(r"\b(?:latino|hispanic|admixed\s+american)", re.I), "AMR"),
)


def _local_variants(text: str) -> list[Tuple[str, Optional[str]]]:
    """Every DISTINCT identifier in the text, in the order typed. No LLM.

    Returns all of them, never just the first: the caller has to see every identifier to
    know it cannot choose between them.
    """
    hgvs = [(m.span(), m.group(0).rstrip(".,;"), (m.group(1) or "").strip() or None)
            for m in _HGVS.finditer(text)]
    hits = [(span[0], v, gene) for span, v, gene in hgvs]
    for m in _RSID.finditer(text):
        # HGVS's trailing \S+ is greedy, so an rsID falling inside an HGVS match is the
        # same variant restated, not a second one.
        if any(a <= m.start() < b for (a, b), _, _ in hgvs):
            continue
        hits.append((m.start(), m.group(0), None))

    out: list[Tuple[str, Optional[str]]] = []
    seen: set[str] = set()
    for _, v, gene in sorted(hits, key=lambda h: h[0]):
        if v.upper() not in seen:          # 'rs1 ... RS1' is one variant, not a conflict
            seen.add(v.upper())
            out.append((v, gene))
    return out


def _local_modifiers(text: str) -> dict:
    """Window / MAF / ancestry straight off the text. No LLM, no network.

    Runs on both paths; on the LLM path the model's values only fill the gaps.
    """
    out: dict = {}

    m = _WINDOW.search(text)
    if m:
        out["window_bp"] = int(round(float(m.group(1)) * _UNITS[m.group(2).lower()]))

    m = _MAF.search(text)
    if m:
        v = float(m.group(1))
        if m.group(2):                      # "MAF >= 10%"
            v /= 100.0
        # A bare value > 1 is deliberately not reinterpreted as a percentage:
        # StructuredQuery rejects it and the user restates it.
        out["common_maf"] = v

    m = _CODES.search(text)
    if m:
        out["ancestry"] = m.group(1)
    else:
        for pat, code in _ANCESTRY_WORDS:
            if pat.search(text):
                out["ancestry"] = code
                break

    return out


def _build(variant: str, gene: Optional[str], mods: dict) -> pb.StructuredQuery:
    """Hand off to the engine. StructuredQuery.__post_init__ owns validation and raises
    ValueError on a bad window/MAF/ancestry. Do not re-validate here: two validators
    drift."""
    return pb.StructuredQuery(variant=variant, gene=gene or None, **mods)


def _describe(mods: dict) -> str:
    bits = []
    if "window_bp" in mods:
        bits.append(f"window +/-{mods['window_bp']:,} bp")
    if "common_maf" in mods:
        bits.append(f"MAF floor {mods['common_maf']}")
    if "ancestry" in mods:
        bits.append(f"ancestry {mods['ancestry']}")
    return "; ".join(bits)


# --------------------------------------------------------------------------- #
# LLM fallback: intent only, and only when no identifier is present
# --------------------------------------------------------------------------- #

# Prompt caching needs a >=4096-token prefix on Haiku 4.5, so the cache_control marker
# below is silently a no-op at this size. It is set for the day the prompt outgrows the
# minimum; the regex fast path, not caching, is what makes this cheap.
_SYSTEM = (
    "You extract search INTENT from a geneticist's free-text request for a PGT-M "
    "linkage marker panel. Reply with ONE JSON object and nothing else:\n"
    '{"variant": str|null, "gene": str|null, "window_bp": int|null, '
    '"ancestry": str|null, "common_maf": float|null}\n\n'
    "HARD RULES:\n"
    "- `variant` must be an identifier the user ACTUALLY TYPED: an rsID (rs...) or an "
    "HGVS expression (NM_/NC_/NG_/NR_... with :c. :g. :n. or :p.). Copy it verbatim.\n"
    "- NEVER invent, recall, complete, or infer a genomic coordinate, chromosome, "
    "position, rsID, or ref/alt allele. You do not know them and must not guess. "
    "Never output anything shaped like 11-17397055-C-T, chr11:17397055, GRCh38:11:..., "
    "or a bare position number. Such output is rejected by the application.\n"
    "- If the user did not type an rsID or HGVS identifier, set variant to null. "
    "Returning null is CORRECT and expected: the app then fails loudly and asks the "
    "user to name the variant. A guess is a serious error; null is not.\n"
    "- gene: HGNC symbol if one is named, else null. It is a hint only, never used to "
    "derive coordinates.\n"
    "- ancestry: one of AFR AMR ASJ EAS FIN NFE SAS MID, else null.\n"
    "- window_bp: flank size in base pairs. common_maf: minor allele frequency floor "
    "as a fraction (0-0.5). Use null for anything the user did not state."
)


def _llm_intent(text: str) -> dict:
    """One tiny Haiku call. Returns the raw parsed JSON dict; the caller validates it
    (_require_identifier), because the model is never trusted. Raises ValueError."""
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise ValueError(f"No rsID or HGVS expression was found in your text, and "
                         f"free-text parsing is switched off here. {_NAME_IT}")
    try:
        import anthropic
    except ImportError as e:
        _log.warning("intent parser unavailable: %s", e)
        raise ValueError(f"Free-text parsing is unavailable. {_NAME_IT}") from e

    client = anthropic.Anthropic()
    try:
        resp = client.messages.create(
            model=MODEL,
            max_tokens=MAX_TOKENS,
            temperature=0,
            system=[{
                "type": "text",
                "text": _SYSTEM,
                "cache_control": {"type": "ephemeral"},
            }],
            messages=[{"role": "user", "content": text}],
        )
    except anthropic.APIError as e:
        # Every failure here is the same to the caller: intent could not be parsed, so
        # 400, never 500. The provider's message goes to the log, not to the visitor.
        _log.warning("intent parser unavailable: %s", e)
        raise ValueError(f"Free-text parsing is unavailable. {_NAME_IT}") from e

    raw = "".join(b.text for b in resp.content if b.type == "text").strip()
    try:
        return json.loads(raw)
    except ValueError:
        m = re.search(r"\{.*\}", raw, re.S)          # tolerate ``` fencing
        if not m:
            _log.warning("intent parser returned unusable output: %r", raw[:200])
            raise ValueError(f"Free-text parsing did not return a usable answer. {_NAME_IT}")
        return json.loads(m.group(0))


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def needs_llm(text: str) -> bool:
    """Would this text require a model call? Regex only: no network, no tokens.

    Lets the caller price a request before serving it. Answers off the same
    _local_variants() that parse() branches on, so the two cannot disagree about what
    costs money. Ambiguous text (several identifiers) is False: parse() raises on it
    without reaching the model.
    """
    return not _local_variants((text or "").strip())


def parse(text: str) -> Tuple[pb.StructuredQuery, bool, str]:
    """Free text -> (StructuredQuery, used_llm, note). Raises ValueError -> 400.

    The returned query carries an opaque `variant` identifier only. Nothing in this
    module knows or asserts where that variant sits in the genome (R1).
    """
    text = (text or "").strip()
    if not text:
        raise ValueError("empty query")

    mods = _local_modifiers(text)
    found = _local_variants(text)

    if len(found) > 1:
        # Which one the user meant is not recoverable from the text: "not rs1, use rs2"
        # and "rs1, or maybe rs2" are the same to a regex. Ask, never guess.
        names = ", ".join(v for v, _ in found)
        raise ValueError(
            f"Your text names more than one variant ({names}), so it is not clear which "
            "panel to build. Send just the one you want."
        )

    if found:
        variant, gene = found[0]
        q = _build(_require_identifier(variant), gene, mods)
        note = "Matched a variant identifier in your text."
        if _describe(mods):
            note += f" Read: {_describe(mods)}."
        return q, False, note

    data = _llm_intent(text)
    variant = _require_identifier(data.get("variant") or "")
    gene = (data.get("gene") or None)

    # Local reading wins wherever we have one; the model only fills gaps.
    for key in ("window_bp", "ancestry", "common_maf"):
        if key not in mods and data.get(key) is not None:
            mods[key] = data[key]
    if mods.get("window_bp") is not None:
        mods["window_bp"] = int(mods["window_bp"])

    q = _build(variant, gene, mods)
    # This path runs only because no identifier could be read from the text, so every
    # character of `variant` came from the model. The note must disclose that.
    note = (
        f"No rsID or HGVS expression could be read from your text, so {MODEL} was asked "
        f"what you meant and it supplied {variant!r}. Check that this is the variant you "
        "intend before you use the panel: it was not taken from anything you typed. "
        "Position, strand and alleles are not taken from the model; they come from a live "
        "lookup of that identifier."
    )
    if _describe(mods):
        note += f" Read: {_describe(mods)}."
    return q, True, note


# --------------------------------------------------------------------------- #
# Self-check
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    GOLDEN = "NM_000352.6(ABCC8):c.3989-9G>A"

    # 1. The golden HGVS parses with no LLM, verbatim, with the gene picked up free.
    q, used_llm, note = parse(GOLDEN)
    assert q.variant == GOLDEN, q.variant
    assert q.gene == "ABCC8", q.gene
    assert used_llm is False, note
    assert q.build == "GRCh38" and q.window_bp == 250_000

    # 2. Bare rsID, also free. Also inside a sentence.
    q, used_llm, _ = parse("rs151344623")
    assert (q.variant, used_llm) == ("rs151344623", False)
    q, used_llm, _ = parse("please build a panel around rs151344623 for this family")
    assert (q.variant, used_llm) == ("rs151344623", False)

    # 2b. TWO identifiers is a question, not a coin flip.
    for text in ["not rs1801133, I mean rs151344623",
                 "rs1 was wrong, use rs2",
                 f"{GOLDEN} or rs151344623?",
                 "rs151344623 vs NM_000352.6:c.3989-9G>A"]:
        try:
            parse(text)
        except ValueError as e:
            assert "more than one" in str(e), e
        else:
            raise AssertionError(f"parse() guessed between two identifiers: {text!r}")
    # ...and ambiguous text is never charged: parse() refuses before it would call out.
    assert needs_llm("not rs1801133, I mean rs151344623") is False

    # 2c. The same identifier twice is one variant, not an ambiguity.
    assert parse("rs151344623 (aka RS151344623)")[0].variant == "rs151344623"
    assert len(_local_variants(f"{GOLDEN} confirmed, {GOLDEN} again")) == 1

    # 3. Prose with no identifier: modifiers still come out locally, no token spent.
    mods = _local_modifiers(
        "markers near the ABCC8 splice mutation with MAF>0.1 in Europeans"
    )
    assert mods["ancestry"] == "NFE", mods
    assert mods["common_maf"] == 0.1, mods

    # 4. Window phrases.
    assert _local_modifiers("+/-100kb")["window_bp"] == 100_000
    assert _local_modifiers("250 kb")["window_bp"] == 250_000
    assert _local_modifiers("1Mb")["window_bp"] == 1_000_000
    assert _local_modifiers("500000 bp")["window_bp"] == 500_000

    # 5. Modifiers ride along the fast path end-to-end, still with no LLM.
    q, used_llm, _ = parse("rs151344623 +/-100kb with MAF>=0.1 in Europeans")
    assert (q.window_bp, q.common_maf, q.ancestry, used_llm) == (
        100_000, 0.1, "NFE", False)

    # 6. Ancestry synonyms, incl. the ordering trap and the plural.
    for text, code in [
        ("in Europeans", "NFE"), ("non-Finnish European", "NFE"),
        ("a Finnish family", "FIN"), ("Ashkenazi Jewish", "ASJ"),
        ("East Asian", "EAS"), ("South Asians", "SAS"), ("Middle Eastern", "MID"),
        ("African American", "AFR"), ("admixed American", "AMR"), ("Latino", "AMR"),
        ("in AFR", "AFR"), ("pop MID", "MID"),
    ]:
        assert _local_modifiers(text).get("ancestry") == code, (text, code)
    # ...and words that merely contain a code must not trigger one.
    assert "ancestry" not in _local_modifiers("markers in the middle of the gene")

    # 7. THE R1 BOUNDARY. Coordinate shapes are refused, whatever they look like.
    for bad in ["11-17397055-C-T", "chr11:17397055", "17397055",
                "GRCh38:11:17397055", "11:17397055", "chr11:17397055-17400000",
                "X:12345", "11_17397055_C_T",
                "chr11:17,397,055", "chr11:17,397,055 C>T",
                "GRCh38 11 17397055 C T", "11 17397055 C T",
                "the variant at chr11 position 17397055",
                '{"chrom":"11","pos":17397055}']:
        try:
            _reject_coordinates(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"_reject_coordinates let a coordinate through: {bad}")

    # ...while real identifiers pass: the allow-list must not be so tight that it refuses
    # input a researcher would legitimately type.
    for ok in [GOLDEN, "rs151344623", "NC_000011.10:g.17397055C>T",
               "NM_000352.6:c.3989-9G>A", "VCV000009088", "ABCC8:c.3989-9G>A",
               "ABCC8 p.R1215Q", "ENST00000389817.8:c.100A>G"]:
        assert _reject_coordinates(ok) == ok, ok

    # 8. A coordinate from the model is a raise, never a pass-through.
    for bad in ["11-17397055-C-T", "chr11:17397055", "ABCC8 splice mutation", ""]:
        try:
            _require_identifier(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"_require_identifier accepted {bad!r}")

    # 9. Engine validation is not duplicated here - it fires through us.
    for bad in ["rs1 +/- 1Mb and 900kb", "rs1 with MAF 0.9"]:
        try:
            parse("rs151344623 with MAF 0.9")
        except ValueError:
            break
    else:
        raise AssertionError("out-of-range MAF should have raised")

    # 10. The LLM-path note must disclose that the model supplied an identifier the user
    #     never typed.
    _real, _llm_intent = _llm_intent, lambda text: {"variant": "rs151344623"}
    try:
        q, used_llm, note = parse("markers near the ABCC8 splice mutation in Europeans")
    finally:
        _llm_intent = _real
    assert used_llm is True and q.variant == "rs151344623"
    assert q.ancestry == "NFE", "local modifiers still win on the LLM path"
    assert "echoed" not in note, note
    assert "supplied" in note and "rs151344623" in note, note
    assert "R1" not in note, "internal rule IDs must not reach the user"

    print("nl.py self-check OK (no LLM calls, no network)")

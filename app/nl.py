"""Natural-language INTENT parser for OriginMarker (R1).

MAY DECIDE: which variant the user means -- named outright, or described in words ("the
sickle cell mutation") and resolved to its standard rsID/HGVS identifier from the model's
own knowledge -- and the search knobs they asked for (window, MAF floor, ancestry).

MAY NEVER DECIDE: a coordinate. No chromosome, position, strand or ref/alt allele, from
regex, LLM or memory. `variant` leaves this module as an opaque identifier that
pb.resolve_variant() looks up live; every coordinate in the app is produced there. An rsID
the model recalled is still only an identifier: the app confirms it and derives its
coordinate by live lookup, and the user is told the model supplied it.

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
from functools import lru_cache
from typing import NamedTuple, Optional, Tuple

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

_VCV = re.compile(r"\b(?:VCV|RCV)\d+(?:\.\d+)?", re.I)

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
    _VCV,                                                      # ClinVar accession
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
    """Narrow _reject_coordinates to the shapes resolve_variant() can actually look up."""
    v = _reject_coordinates(s)
    if not (_RSID.fullmatch(v) or _HGVS.fullmatch(v) or _VCV.fullmatch(v)):
        raise ValueError(
            f"{s!r} is not an rsID, HGVS or ClinVar identifier. Name the variant the way "
            "it appears in the report (e.g. rs151344623 or "
            "NM_000352.6(ABCC8):c.3989-9G>A)."
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
    for pat in (_RSID, _VCV):
        for m in pat.finditer(text):
            # HGVS's trailing \S+ is greedy, so an identifier falling inside an HGVS match
            # is the same variant restated, not a second one.
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


# HGNC symbols are uppercase alphanumeric, 2-10 chars, occasionally hyphenated (MT-ND1).
# The lookbehind keeps HGVS suffixes out ('p.R1215Q' is not a gene called R1215Q); the
# 10-char cap keeps accessions out (VCV000009088, ENST00000389817).
# Two shapes: an all-caps symbol (ABCC8, MT-ND1), and HGNC's Cxorfy convention, whose
# lowercase "orf" a caps-only pattern reads straight past (C9orf72 is the commonest genetic
# cause of ALS, so missing it is not academic).
_GENE = re.compile(r"(?<![a-z]\.)\b(?:[A-Z]\d+orf\d+|[A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]{1,4})?)\b")

# A hand-written exclusion set, not an HGNC download: fetch the real list only if this
# becomes a losing game. Everything here is jargon a geneticist writes in caps that is not
# a gene. Note what is absent: MB is myoglobin, a real symbol, so units are stripped before
# extraction rather than blacklisted here.
_NOT_GENES = frozenset(pb.GNOMAD_POPS.values()) | {
    "DNA", "RNA", "SNP", "SNV", "CNV", "INDEL", "PGT", "PGD", "PGS", "IVF",
    "HGVS", "VCF", "MAF", "LD", "GRCH37", "GRCH38", "HG19", "HG38", "OMIM",
    "ACMG", "ID", "OK",
    # Assay and report jargon. Every one of these read as a symbol and produced a
    # confident complaint that the user had named a gene they had not.
    "NIPT", "CVS", "PDF", "REI", "WES", "WGS", "NGS", "MLPA", "FISH", "QPCR",
    "PCR", "IUI", "ART", "ICSI", "TE", "ADO", "SNV", "VUS", "LOH",
}


def _named_genes(text: str) -> list[str]:
    """Gene symbols the USER typed, uppercase, in the order typed. Regex only, no LLM.

    Deliberately read from the user's own words and never from the model's `gene` field:
    this is the string the model gets checked against, so it may not come from the model.
    """
    # Units first: "500MB" is a window, and MB is also myoglobin. Blacklisting MB would
    # lose the gene; letting it through would read every window as one. _WINDOW already
    # knows which is which, so remove what it claims and read genes from the rest.
    text = _WINDOW.sub(" ", text)
    out: list[str] = []
    for m in _GENE.finditer(text):
        g = m.group(0)
        # The stem too, so PGT-M and PGT-A fall out with PGT while MT-ND1 stays.
        if g.upper() in _NOT_GENES or g.split("-")[0].upper() in _NOT_GENES or g in out:
            continue
        out.append(g)
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
    "You read a geneticist's free-text request for a PGT-M linkage marker panel and "
    "identify which pathogenic variant they mean. Reply with ONE JSON object and nothing "
    'else: {"variant": str|null, "gene": str|null, "window_bp": int|null, '
    '"ancestry": str|null, "common_maf": float|null}\n\n'
    "variant - the standard identifier for the variant the user names or describes:\n"
    "- If they typed an rsID (rs...) or an HGVS expression (NM_/NC_/NG_/NR_/ENST... with "
    ":c. :g. :n. or :p.), copy it verbatim.\n"
    "- If they described it in words ('the sickle cell mutation', 'factor V Leiden', 'the "
    "common CF deletion'), give the well-known rsID or HGVS identifier for it from your "
    "own knowledge. This is the reason you are being asked: the identifier is not in their "
    "text, and naming it from what they described is your job.\n"
    "- If you cannot identify ONE specific variant with confidence, set variant to null. "
    "Null is a correct answer: the app then asks the user to name it. A confident wrong "
    "guess is far worse than null.\n\n"
    "NEVER output a genomic COORDINATE. Give the IDENTIFIER, never the position it sits "
    "at: no chromosome, position, strand, or ref/alt allele, and nothing shaped like "
    "11-5227002-T-A, chr11:5227002, GRCh38:11:5227002, or a bare position number. The app "
    "looks the coordinate up from the identifier and rejects a coordinate here. An rsID or "
    "an accession-anchored HGVS expression is an identifier, not a coordinate.\n\n"
    "gene - HGNC symbol if the variant's gene is clear, else null (a hint only, never used "
    "to derive coordinates).\n"
    "ancestry - one of AFR AMR ASJ EAS FIN NFE SAS MID if the user names a population, "
    "else null.\n"
    "window_bp - flank size in base pairs, ONLY if the user stated one, else null.\n"
    "common_maf - minor allele frequency floor as a 0-0.5 fraction, ONLY if the user "
    "stated one, else null.\n"
    "Use null for any modifier the user did not state; do not fill modifiers from the "
    "variant you identified."
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

    # max_retries, because a 529 "Overloaded" is weather, not a verdict. Observed in
    # production: three consecutive parses died on it and told the user free text was
    # unavailable, while the key was fine. Every other network call in this app retries.
    client = anthropic.Anthropic(max_retries=4)
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
        # 400, never 500: the provider's message goes to the log, not to the visitor. But
        # a busy model and a broken one are different facts, and "unavailable" reads as the
        # feature being off rather than as something to retry.
        _log.warning("intent parser failed: %s", e)
        status = getattr(e, "status_code", None)
        if status in (429, 529) or isinstance(e, anthropic.APIStatusError) and status and status >= 500:
            raise ValueError(
                f"The model is busy right now, and retrying did not clear it. This is a "
                f"transient failure and says nothing about your text. Try again in a "
                f"moment, or skip the wait: {_NAME_IT}"
            ) from e
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


def _norm(text: str) -> str:
    """Cache key: the same prose asked twice must bill once."""
    return " ".join(text.split()).casefold()


@lru_cache(maxsize=256)
def _cached_intent(key: str) -> dict:
    """_llm_intent keyed on normalised prose. The model sees the NORMALISED text.

    Looks _llm_intent up on the module so a monkeypatched stub is honoured; call
    _cached_intent.cache_clear() when swapping the stub, or the swap is invisible.
    lru_cache does not memoise exceptions, so a failed call is retried, not pinned.
    """
    return _llm_intent(key)


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

class NLResponse(NamedTuple):
    """What parse() returns.

    named_genes are the symbols the USER typed. parse() cannot resolve, so it cannot
    check them against the variant's real gene: the caller must, once it has resolved.
    """
    query: pb.StructuredQuery
    used_llm: bool
    note: str
    named_genes: list[str]

def needs_llm(text: str) -> bool:
    """Would this text require a model call? Regex only: no network, no tokens.

    Lets the caller price a request before serving it. Answers off the same
    _local_variants() that parse() branches on, so the two cannot disagree about what
    costs money. Ambiguous text (several identifiers) is False: parse() raises on it
    without reaching the model.
    """
    return not _local_variants((text or "").strip())


def parse(text: str) -> NLResponse:
    """Free text -> NLResponse. Raises ValueError -> 400.

    The returned query carries an opaque `variant` identifier only. Nothing in this
    module knows or asserts where that variant sits in the genome (R1).
    """
    text = (text or "").strip()
    if not text:
        raise ValueError("empty query")

    mods = _local_modifiers(text)
    named = _named_genes(text)
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
        return NLResponse(q, False, note, named)

    data = dict(_cached_intent(_norm(text)))   # copy: the cached dict is shared
    variant = _require_identifier(data.get("variant") or "")
    gene = (data.get("gene") or None)

    # The gene cross-check is NOT done here, against the model's own `gene` field. That
    # asks the model to validate itself: a model wrong about the variant can be wrong about
    # its gene in the same breath and pass. It is also unable to distinguish a real
    # disagreement from an alias the user typed (SUR1 is ABCC8, ND1 is MT-ND1), so it
    # refused correct answers with no way past it, after billing for them. The check that
    # means something compares named_genes against the RESOLVED record's gene, which is
    # authoritative and arrives later. It lives at the resolve step.

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
    return NLResponse(q, True, note, named)


# --------------------------------------------------------------------------- #
# Self-check
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    GOLDEN = "NM_000352.6(ABCC8):c.3989-9G>A"

    # 1. The golden HGVS parses with no LLM, verbatim, with the gene picked up free.
    q, used_llm, note, named = parse(GOLDEN)
    assert q.variant == GOLDEN, q.variant
    assert q.gene == "ABCC8", q.gene
    assert used_llm is False, note
    assert q.build == "GRCh38" and q.window_bp == 250_000
    assert named == ["ABCC8"], named

    # 2. Bare rsID, also free. Also inside a sentence.
    q, used_llm, _, _ = parse("rs151344623")
    assert (q.variant, used_llm) == ("rs151344623", False)
    q, used_llm, _, _ = parse("please build a panel around rs151344623 for this family")
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
    q, used_llm, _, _ = parse("rs151344623 +/-100kb with MAF>=0.1 in Europeans")
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
    calls: list[str] = []

    def _stub(payload):
        def f(text):
            calls.append(text)
            return payload
        return f

    def _swap(fn):
        """Install a stub AND drop the cache: a warm entry hides the swap."""
        global _llm_intent
        _llm_intent = fn
        _cached_intent.cache_clear()

    _real = _llm_intent
    try:
        _swap(_stub({"variant": "rs151344623"}))
        q, used_llm, note, named = parse(
            "markers near the ABCC8 splice mutation in Europeans")
        assert used_llm is True and q.variant == "rs151344623"
        assert q.ancestry == "NFE", "local modifiers still win on the LLM path"
        assert "echoed" not in note, note
        assert "supplied" in note and "rs151344623" in note, note
        assert "R1" not in note, "internal rule IDs must not reach the user"
        assert named == ["ABCC8"], named

        # 11. THE GENE CROSS-CHECK IS NOT MADE HERE. The user said ABCC8 and rs334 is HBB,
        #     but the only gene this layer could compare against is the model's own claim,
        #     and a model wrong about the variant can be wrong about its gene in the same
        #     breath. parse() surfaces what the user named and lets the resolver, which
        #     has the authoritative gene, make the call.
        _swap(_stub({"variant": "rs334", "gene": "HBB"}))
        r = parse("the ABCC8 splice mutation we discussed")
        assert r.query.variant == "rs334", r.query
        assert r.named_genes == ["ABCC8"], r.named_genes

        # ...an agreeing gene passes, and so does a model that names no gene at all.
        _swap(_stub({"variant": "rs151344623", "gene": "ABCC8"}))
        assert parse("the ABCC8 splice mutation").query.variant == "rs151344623"
        _swap(_stub({"variant": "rs151344623"}))
        assert parse("that splice mutation from the report").used_llm is True

        # 12. COST. Identical prose bills once; whitespace and case must not defeat it.
        _swap(_stub({"variant": "rs151344623"}))
        calls.clear()
        for text in ["the  splice   mutation", "The Splice Mutation",
                     "  the splice mutation  ", "the splice mutation"]:
            assert parse(text).used_llm is True
        assert len(calls) == 1, f"intent cache billed {len(calls)} times for one query"

        # ...and the fast path never reaches the model at all, however chatty the text.
        calls.clear()
        for text in [GOLDEN, "rs151344623", f"panel around {GOLDEN} please",
                     "rs151344623 +/-100kb MAF>=0.1 in Europeans",
                     "VCV000009088 for the ABCC8 family",
                     "NC_000011.10:g.17397055C>T",
                     "rs151344623 (aka RS151344623)"]:
            assert needs_llm(text) is False, text
            assert parse(text).used_llm is False, text
        assert calls == [], f"regex fast path reached the model: {calls}"
    finally:
        _swap(_real)

    # 13. Gene extraction reads the USER, so it must not invent genes out of HGVS,
    #     accessions, ancestry codes or jargon.
    assert _named_genes("markers near the ABCC8 splice mutation in Europeans") == ["ABCC8"]
    assert _named_genes("MT-ND1 and BRCA1 and HBB") == ["MT-ND1", "BRCA1", "HBB"]
    assert _named_genes("ABCC8 p.R1215Q") == ["ABCC8"], "p.R1215Q is not a gene"
    assert _named_genes(GOLDEN) == ["ABCC8"]
    for quiet in ["rs151344623", "VCV000009088", "ENST00000389817.8:c.100A>G",
                  "NC_000011.10:g.17397055C>T", "SNP panel, GRCh38, MAF 0.1, 250 KB, AFR",
                  "PGT-M DNA HGVS", "a splice mutation in the middle of the gene"]:
        assert _named_genes(quiet) == [], (quiet, _named_genes(quiet))

    # Gene extraction, every case an adversarial pass found. Each was a live miss or a
    # live false alarm, not a hypothetical.
    assert _named_genes("the ABCC8 splice mutation") == ["ABCC8"]
    assert _named_genes("C9orf72 repeat expansion carrier") == ["C9orf72"]   # HGNC Cxorfy
    assert _named_genes("the MB variant") == ["MB"]                          # myoglobin
    assert _named_genes("rs334 with a 500MB window") == []                   # MB the unit
    assert _named_genes("markers near HBB, 250kb window") == ["HBB"]
    assert _named_genes("the splice mutation from the NIPT report") == []
    assert _named_genes("that mutation from the CVS sample") == []
    assert _named_genes("the MT-ND1 variant") == ["MT-ND1"]
    assert _named_genes("PGT-M panel for ABCC8") == ["ABCC8"]
    # A symbol the user lowercases reads as prose. Known and documented: the caps rule IS
    # the detector, and relaxing it makes every English word a symbol.
    assert _named_genes("the abcc8 splice mutation") == []

    # The model's own gene claim is never the thing checked against: a model wrong about
    # the variant can be wrong about its gene too, and an alias (SUR1 for ABCC8) is not a
    # disagreement. parse() must not refuse on it.
    _orig = globals().get("_llm_intent")
    globals()["_llm_intent"] = lambda t: {"variant": "rs151344623", "gene": "ABCC8"}
    _cached_intent.cache_clear()
    try:
        r = parse("the SUR1 splice mutation we discussed")     # SUR1 is ABCC8's alias
        assert r.query.variant == "rs151344623", r.query
        assert r.named_genes == ["SUR1"], r.named_genes
    finally:
        globals()["_llm_intent"] = _orig
        _cached_intent.cache_clear()

    # 14. The prompt must ask the model to do the one thing this path exists for: name an
    #     identifier for a variant the user only DESCRIBED. An earlier prompt forbade the
    #     model to "recall ... an rsID", which is self-defeating -- it made every advertised
    #     free-text example ("the sickle cell mutation") return null, since the identifier is
    #     never in the user's text. The regex is the R1 control, not this prompt, but a
    #     prompt that re-locks the feature fails silently and only on the live path, so pin
    #     it here. Coordinates must still be forbidden in the same breath.
    _p = _SYSTEM.lower()
    assert "from your own knowledge" in _p, "the prompt no longer lets the model name a described variant"
    assert "never output a genomic coordinate" in _p, "the prompt no longer forbids coordinates"
    assert "actually typed" not in _p, "the prompt re-locked the model to only echo typed identifiers"

    print("nl.py self-check OK (no LLM calls, no network)")

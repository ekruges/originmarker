"""
panelbuilder - candidate linked-SNP marker panels for parental-origin determination.

Given a pathogenic variant, build a ranked, reusable menu of common flanking SNPs that
can - after the carrier is genotyped and the markers are phased in a lab - report which
parental chromosome an embryo inherited. This is PGT-M linkage / karyomapping.

Data sources (free, GRCh38): ClinVar (E-utilities), Ensembl REST, gnomAD v4 (GraphQL).

Domain rules. Comments across the codebase cite these by number:
  R1  coordinates/rsIDs/strand/ref-alt come ONLY from live API calls, never from memory.
  R2  markers ranked by heterozygosity + proximity, NEVER by LD with the pathogenic variant.
  R3  output is CANDIDATE markers; per-family genotyping + phasing still required (caller's job).
  R4  expected heterozygosity (2pq) is a population prior, labelled as such.
  R5  both genomic sides must be covered; under-coverage is flagged.
  R6  reference build is explicit on every record (GRCh38).
  R7  strand handled: transcript-sense vs genomic ref/alt reconciled for minus-strand genes.
  R8  every page footer and every export carries pb.DISCLAIMER verbatim, never paraphrased.

This module is pure/importable and has no web framework or LLM dependency.
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import os
import re
import contextvars
import gzip
import logging
import threading
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Union

import genetic_map

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL = "https://rest.ensembl.org"
GNOMAD = "https://gnomad.broadinstitute.org/api"

# Ancestry population ids in gnomAD v4 (top-level only; ignore ids containing "_").
GNOMAD_POPS = {"afr": "AFR", "amr": "AMR", "asj": "ASJ", "eas": "EAS",
               "fin": "FIN", "nfe": "NFE", "sas": "SAS", "mid": "MID"}

DEFAULT_WINDOW = 250_000
CALL_RATE_AN_FLOOR = 10_000     # require gnomAD genome AN >= this (decent call rate)
COMMON_MAF = 0.05               # candidate pool floor
# Per-population AN floor: the global floor says the SITE is well called, not that a given
# population was sampled there. AN 200 is 100 people, +/-7pp at p=0.5, low enough to admit
# gnomAD's smaller groups (MID, ASJ).
MIN_POP_AN = 200
REGION_CHUNK = 20_000           # gnomAD region slice size (bp)
MAX_WORKERS = 8                 # concurrent API fetches
CACHE_DIR = Path(os.environ.get("PANELBUILDER_CACHE", ".panelbuilder_cache"))
# ClinVar publishes weekly and gnomAD/Ensembl quarterly, so a week cannot serve a
# classification that changed two releases ago. PANELBUILDER_CACHE_TTL=0 disables expiry:
# the recorded-fixture test runs need that, since their cache IS the source of truth and
# expiring it would turn an offline suite into a live one.
CACHE_TTL_S = float(os.environ.get("PANELBUILDER_CACHE_TTL", 7 * 86_400))
CONTACT_EMAIL = os.environ.get("NCBI_CONTACT_EMAIL")   # optional etiquette
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")          # optional, raises rate limit


# --------------------------------------------------------------------------- #
# HTTP layer: retrying, cached, thread-safe-enough
# --------------------------------------------------------------------------- #

class ApiError(RuntimeError):
    pass


log = logging.getLogger(__name__)

# The Ensembl release that answered THIS build. app/main warms it off the request path and
# writes it here; build() stamps it into provenance. It must not be re-read at render time:
# a frozen panel would then re-render with whatever release the live server has moved on
# to. None renders as "unknown", which is an honest gap.
_ENSEMBL_RELEASE: dict = {"value": None}


def ensembl_release() -> Optional[int]:
    """The warmed Ensembl release, or None if it has not landed yet."""
    return _ENSEMBL_RELEASE["value"]


def set_ensembl_release(v: Optional[int]) -> None:
    _ENSEMBL_RELEASE["value"] = v

# What THIS build fetched versus what it read off disk; build() resets it and reads it back
# to stamp provenance with the data's age. A ContextVar over a MUTABLE dict, not a module
# global (concurrent builds in one process would share one ledger) and not a thread-local
# (_note_fetch is called from inside enumerate_candidates' worker pool, which a
# thread-local does not survive). The lock guards those workers against each other.
_fetch_lock = threading.Lock()
_fetch_ledger: contextvars.ContextVar[dict] = contextvars.ContextVar("fetch_ledger")


def _new_fetch_log() -> dict:
    return {"oldest": None, "from_cache": 0, "from_network": 0}


def _fetch_log() -> dict:
    """This build's ledger, or a fresh private one if we are not inside a build, so a bare
    _http call cannot contaminate a build's data date."""
    try:
        return _fetch_ledger.get()
    except LookupError:
        d = _new_fetch_log()
        _fetch_ledger.set(d)
        return d


def _utc(epoch: Optional[float]) -> Optional[str]:
    """Epoch seconds to the ISO-Z form used throughout provenance."""
    return None if epoch is None else time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(epoch))


def _reset_fetch_log() -> None:
    """Start a fresh ledger for this build. Rebinds rather than mutating in place: a
    concurrent build may be holding the old dict."""
    _fetch_ledger.set(_new_fetch_log())


def _note_fetch(epoch: float, cached: bool) -> None:
    flog = _fetch_log()
    with _fetch_lock:
        if flog["oldest"] is None or epoch < flog["oldest"]:
            flog["oldest"] = epoch
        flog["from_cache" if cached else "from_network"] += 1


# The span enumerate_candidates actually queried, per side, which is NOT always the window
# asked for: a variant near a telomere gets a truncated flank, and R5 requires that be
# disclosed rather than quietly shortened. Thread-local, not a module global, because
# concurrent builds share this process; sound only because _note_span and build()'s
# read-back both run on the build's own thread, before the fan-out to the worker pool.
# _fetch_log cannot use this trick, see _fetch_ledger above.
_span_state = threading.local()


def _span() -> dict:
    """This build's queried span. Empty before enumerate_candidates has run."""
    return getattr(_span_state, "span", {})


def _note_span(pos: int, lo: int, hi: int, window: int, chrom: str,
               contig_len: Optional[int]) -> None:
    flags = []
    if pos - lo < window:
        flags.append(f"Left flank truncated to {pos - lo:,} bp of the {window:,} bp "
                     f"requested: chr{chrom} starts at position 1.")
    if contig_len is None:
        flags.append(f"Length of chr{chrom} was unavailable from Ensembl, so the right "
                     f"flank was not checked against the end of the chromosome.")
    elif hi - pos < window:
        flags.append(f"Right flank truncated to {hi - pos:,} bp of the {window:,} bp "
                     f"requested: chr{chrom} ends at {contig_len:,}.")
    # No lock: nothing outside this thread can see or touch this build's span.
    _span_state.span = dict(queried_start=lo, queried_stop=hi, window_bp=window,
                            left_bp=pos - lo, right_bp=hi - pos,
                            contig_length=contig_len, flags=flags)


def _cache_key(method: str, url: str, body: Optional[bytes]) -> Path:
    h = hashlib.sha256()
    h.update(method.encode()); h.update(url.encode())
    if body:
        h.update(body)
    return CACHE_DIR / f"{h.hexdigest()}.json"


def _body_is_expected(url: str, txt: str) -> bool:
    """Is this response the content type this endpoint is supposed to return?

    efetch is the only XML endpoint here; esearch, esummary, Ensembl and gnomAD answer
    JSON. Checked on READ as well as on write: a poisoned entry outlives the process that
    fetched it. Must stay a check on the RAW TEXT, not a parse: a 502 page can be
    well-formed XML, so ET.fromstring succeeding proves nothing about who answered.
    """
    s = txt.lstrip()
    return s.startswith("<?xml") if "efetch" in url else s.startswith(("{", "["))


def _cache_read(ck: Path) -> Optional[tuple[str, float, Path]]:
    """(text, mtime, path) for a cache entry, or None if there is no hit.

    A .gz beside the key reads too. The committed test fixtures are stored gzipped
    because recorded gnomAD region responses are ~200 MB of JSON at rest and compress
    about 9:1. Live writes stay plain, so a working cache directory needs no tooling.
    """
    gz = ck.with_name(ck.name + ".gz")
    if ck.exists():
        return ck.read_text(), ck.stat().st_mtime, ck
    if gz.exists():
        with gzip.open(gz, "rt") as f:
            return f.read(), gz.stat().st_mtime, gz
    return None


def _http(method: str, url: str, *, body: Optional[bytes] = None,
          headers: Optional[dict] = None, timeout: int = 60,
          tries: int = 4, use_cache: bool = True) -> str:
    """One retrying HTTP call with optional on-disk caching. Returns response text."""
    ck = _cache_key(method, url, body) if use_cache else None
    hit = _cache_read(ck) if ck else None
    if hit:
        # mtime, not now: a cache hit is data of some age and provenance reports it as such.
        txt, mtime, path = hit
        if not _body_is_expected(url, txt):
            path.unlink(missing_ok=True)        # poison: drop it, never serve it
        elif CACHE_TTL_S and time.time() - mtime > CACHE_TTL_S:
            path.unlink(missing_ok=True)        # stale: refetch below
        else:
            _note_fetch(mtime, cached=True)
            return txt

    hdrs = {"User-Agent": "panelbuilder", "Accept": "application/json"}
    if headers:
        hdrs.update(headers)
    last: Optional[Exception] = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers=hdrs, method=method)
            with urllib.request.urlopen(req, timeout=timeout) as f:
                txt = f.read().decode()
            # An empty or HTML body under load is retryable, and must never be cached.
            if not _body_is_expected(url, txt):
                raise ApiError(f"expected {'XML' if 'efetch' in url else 'JSON'}, got "
                               f"{txt.lstrip()[:40]!r}")
            if ck:
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                # Write beside the key, then rename over it: a reader sees the old entry or
                # the new one, never a half-written one.
                tmp = ck.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")
                tmp.write_text(txt)
                os.replace(tmp, ck)
            _note_fetch(time.time(), cached=False)
            return txt
        except Exception as e:  # noqa: BLE001 - urllib raises a zoo of types
            last = e
            time.sleep(1.5 * (i + 1))
    raise ApiError(f"{method} {url[:80]} failed after {tries} tries: {last}")


def _get(base: str, path: str, params: Optional[dict] = None, **kw) -> str:
    url = base + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    return _http("GET", url, **kw)


def _graphql(query: str, variables: dict, **kw) -> dict:
    body = json.dumps({"query": query, "variables": variables}).encode()
    txt = _http("POST", GNOMAD, body=body,
                headers={"Content-Type": "application/json"}, **kw)
    return json.loads(txt)


def _eutils_params(extra: dict) -> dict:
    p = dict(extra)
    if CONTACT_EMAIL:
        p["email"] = CONTACT_EMAIL
        p["tool"] = "panelbuilder"
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _revcomp(allele: str) -> Optional[str]:
    """Reverse complement of a VCF allele, or None if it is not plain A/C/G/T.

    None rather than a placeholder: a symbolic or empty allele has no transcript-sense
    reading, and the caller says so in words instead of printing a "?" base (R7).
    Orients alleles only; it does not re-left-align an indel in transcript space.
    """
    out = [_COMP.get(b) for b in (allele or "").upper()]
    return "".join(reversed(out)) if out and all(out) else None


@dataclass
class VariantRecord:
    """Canonical, API-resolved identity of the pathogenic variant (R1, R6, R7)."""
    query: str
    rsid: Optional[str]
    gene: Optional[str]
    strand: Optional[int]
    chrom: str
    pos_grch38: int
    vcf_ref: str
    vcf_alt: str
    clinical_significance: Optional[str]
    review_status: Optional[str]
    clinvar_accession: Optional[str]
    build: str = "GRCh38"
    # Display and conversion labelling only: every coordinate the pipeline computes on is
    # pos_grch38, and the two are never mixed (R6).
    pos_grch37: Optional[int] = None
    build_note: Optional[str] = None

    def transcript_sense_change(self) -> str:
        """Genomic ref>alt mapped to transcript sense (reverse complement if minus) (R7).

        Unknown strand is its own answer and is said in words: it must never fall through
        to "plus strand", which is exactly backwards for a minus-strand gene. Alleles are
        reverse complemented, not complemented base-by-base: reversal is a no-op for an
        SNV but not for an indel.
        """
        if self.strand == -1:
            ref, alt = _revcomp(self.vcf_ref), _revcomp(self.vcf_alt)
            if ref and alt:
                return f"{ref}>{alt} (transcript sense; minus strand)"
            return (f"{self.vcf_ref}>{self.vcf_alt} (genomic; alleles are not plain "
                    f"A/C/G/T, so this is NOT converted to transcript sense)")
        if self.strand == 1:
            return f"{self.vcf_ref}>{self.vcf_alt} (plus strand)"
        return (f"{self.vcf_ref}>{self.vcf_alt} (genomic; strand unknown, so this is NOT "
                f"converted to transcript sense)")


def fmt_af(af: Optional[float]) -> str:
    """Render an allele frequency the way the UI does: 1.77e-4, not Python's zero-padded
    1.77e-04. Shared so the server, the exports and the frontend agree."""
    if af is None:
        return "unknown"
    return f"{af:.2e}".replace("e-0", "e-").replace("e+0", "e+")


@dataclass
class Rarity:
    gnomad_af_genome: Optional[float]
    gnomad_ac_genome: Optional[int]
    gnomad_an_genome: Optional[int]
    thousand_genomes_ac: Optional[int]
    population_LD_usable: bool
    reason: str
    # "defined" | "undefined" | "unknown". Three states, because population_LD_usable is a
    # two-way flag and LD nobody could look up is neither "exists" nor "cannot". Kept here
    # so the UI colours its badge from one string rather than re-deriving the verdict and
    # disagreeing with `reason` beside it.
    ld_status: str


@dataclass
class Marker:
    rsid: str
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float
    maf: float
    het: float                         # global expected het 2pq (R4: a prior)
    het_max_pop: float
    dist: int                          # signed bp from variant
    side: str                          # "lower coord" | "higher coord" (see annotate)
    tier: str
    per_pop_maf: dict = field(default_factory=dict)
    an: Optional[int] = None                 # gnomAD genome AN at this site
    per_pop_an: dict = field(default_factory=dict)
    ensembl_pos_check: Optional[str] = None  # "ok" | "MISMATCH:<pos>" | None
    # Recombination context from the bundled genetic map (see genetic_map.py).
    cm: Optional[float] = None               # genetic distance to variant
    recomb_fraction: Optional[float] = None  # Haldane theta for that distance
    hotspot_between: Optional[bool] = None   # a recombination hotspot sits in between
    # True => cm is a 1 cM/Mb approximation, not a map reading. Must be labelled as such.
    map_approx: Optional[bool] = None


@dataclass
class PanelResult:
    variant: VariantRecord
    rarity: Rarity
    candidates: list            # list[Marker], full common pool, ranked
    recommended: list           # list[Marker], balanced subset
    coverage: dict              # per-side/tier counts + under-coverage flags
    params: dict
    provenance: dict

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


@dataclass
class StructuredQuery:
    """The typed front door to build(), and the ONLY shape a natural-language layer may
    produce.

    Note what is absent: no chrom, no pos, no ref/alt, no strand. The variant's identity is
    an opaque HGVS/rsID string that resolve_variant() must look up live, so an NL layer
    cannot inject a coordinate even if it hallucinates one. R1 by construction. The two
    nl_* fields below are prose and a model id, so they carry no coordinate either.
    """
    variant: str                        # HGVS or rsID, verbatim from the user
    gene: Optional[str] = None          # hint only; never used to derive coordinates
    window_bp: int = DEFAULT_WINDOW
    build: str = "GRCh38"
    ancestry: Optional[str] = None      # re-rank by ancestry-matched 2pq
    common_maf: float = COMMON_MAF
    cross_check: bool = True
    # Provenance, never input: build() must not read these to decide anything. They record
    # that a model, not the user, chose `variant`, so the exports can say so.
    nl_text: Optional[str] = None       # the user's original prose, verbatim
    nl_model: Optional[str] = None      # the model id that proposed the variant

    def __post_init__(self):
        if not (self.variant or "").strip():
            raise ValueError("StructuredQuery.variant is required")
        if self.build not in ("GRCh38", "GRCh37"):
            raise ValueError(f"unsupported build {self.build!r} (R6)")
        if not 1_000 <= self.window_bp <= 2_000_000:
            raise ValueError(f"window_bp out of range: {self.window_bp}")
        if not 0.0 <= self.common_maf < 0.5:
            raise ValueError(f"common_maf out of range: {self.common_maf}")
        if self.ancestry:
            self.ancestry = self.ancestry.upper()
            if self.ancestry not in GNOMAD_POPS.values():
                raise ValueError(f"unknown ancestry {self.ancestry!r}; "
                                 f"expected one of {sorted(GNOMAD_POPS.values())}")


# --------------------------------------------------------------------------- #
# Step 1: resolve variant (ClinVar -> fallback Ensembl for bare rsID)  (R1)
# --------------------------------------------------------------------------- #

_RSID_RE = re.compile(r"^rs\d+$", re.I)
_HGVS_RE = re.compile(
    r"^(?:N[MCGRP]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?)"   # accession
    r"\s*(?:\([^()]{1,32}\))?"                             # optional (GENE)
    r"\s*:\s*[cgmnpr]\.\S+$",
    re.I,
)
_VCV_RE = re.compile(r"^(?:VCV|RCV)\d+(?:\.\d+)?$", re.I)


def _looks_like_identifier(q: str) -> bool:
    q = q.strip()
    return bool(_RSID_RE.match(q) or _HGVS_RE.match(q) or _VCV_RE.match(q))


def _norm_hgvs(s: str) -> str:
    """Comparable form: drop the (GENE) parenthetical, the accession version, space, case.

    ClinVar's title carries the gene and whatever accession version is current; a user's
    paste often carries neither. The base accession and the change still have to match
    exactly, so NM_000352 vs NM_000353 and c.3989-9G>A vs c.3989-9G>T are still caught.
    """
    s = re.sub(r"\([^()]*\)", "", s)                 # (GENE)
    s = re.sub(r"^([A-Za-z_]+\d+)\.\d+", r"\1", s)   # accession version
    return s.replace(" ", "").lower()


# ClinVar obj_type values describing a COMBINATION of variants rather than one. Matched on
# a substring because ClinVar's casing and spacing vary across records.
_COMBINATION_OBJ = ("haplotype", "compound", "diplotype", "phase", "distinct chromosome")


def _is_combination(summ: dict) -> bool:
    """True when the record describes several variants together, not a single allele."""
    obj = (summ.get("obj_type") or "").lower()
    return any(k in obj for k in _COMBINATION_OBJ)


def _assert_record_matches(query: str, summ: dict, vset: dict) -> None:
    """Confirm the ClinVar record we got back is the one that was asked for.

    ClinVar's esearch is a FULL-TEXT search, not a lookup: it fuzzy-matches prose, so a hit
    can be a real record for an unrelated gene. Every downstream R1 check still passes on
    it, because the coordinate did come from a live API. It is just the wrong variant, and
    nothing about the result looks wrong, so the hit is reconciled before it is trusted.
    """
    q = query.strip()
    title = summ.get("title") or vset.get("variation_name") or ""
    accession = summ.get("accession") or ""

    if _RSID_RE.match(q):
        want = q.lower().removeprefix("rs")
        got = {str(x.get("db_id")) for x in (vset.get("variation_xrefs") or [])
               if x.get("db_source") == "dbSNP"}
        if want not in got:
            raise ApiError(
                f"ClinVar returned {accession or 'a record'} ({title!r}) for {query!r}, but that "
                f"record's dbSNP ids are {sorted(got) or 'none'}. Refusing to report a variant "
                f"that was not asked for.")
        return

    if _VCV_RE.match(q):
        if q.lower().split(".")[0] not in accession.lower():
            raise ApiError(f"ClinVar returned accession {accession!r} for {query!r}. "
                           f"Refusing to report a variant that was not asked for.")
        return

    # HGVS: the record's own name must be the same variant.
    nq, nt = _norm_hgvs(q), _norm_hgvs(title)
    if nq and nt and (nq == nt or nq in nt or nt in nq):
        return
    raise ApiError(
        f"ClinVar's best match for {query!r} is {accession or 'a record'} ({title!r}), which is "
        f"not the same variant. ClinVar's search is full-text and will return a near miss; "
        f"refusing to build a panel on it. Check the identifier, or paste the rsID.")


def resolve_variant(query: str, build: str = "GRCh38") -> VariantRecord:
    """Resolve an HGVS/rsID to a canonical GRCh38 record via live API only (R1).

    `build` describes the assembly the *user* is thinking in. The pipeline always works
    in GRCh38; if the caller said GRCh37 we still return GRCh38 coordinates and attach an
    explicit conversion note plus the GRCh37 position, so the two are labelled rather
    than silently mixed (R6).
    """
    # Gate before searching: esearch is full-text, so prose does not fail, it matches
    # something. Only an identifier may be looked up.
    if not _looks_like_identifier(query):
        raise ApiError(
            f"{query!r} is not a variant identifier. Enter an rsID (rs151344623), an HGVS "
            f"expression (NM_000352.6(ABCC8):c.3989-9G>A) or a ClinVar accession "
            f"(VCV000009088).")

    # retmax is explicit: the default is 20, and relevance order is not usefulness order
    # (see the reconcile loop below). This widens what is CONSIDERED, never what is
    # ACCEPTED, since every candidate still has to clear _assert_record_matches.
    found = json.loads(_get(EUTILS, "/esearch.fcgi",
                            _eutils_params({"db": "clinvar", "term": query,
                                            "retmode": "json",
                                            "retmax": "50"})))["esearchresult"]
    ids = found["idlist"]
    total = int(found.get("count") or len(ids))
    if not ids:
        return _resolve_via_ensembl(query, build=build)

    # Take the hit that RECONCILES, not the first one, and consider EVERY id esearch
    # returned: relevance rank is not a ranking of "is this the variant you named", so
    # neither idlist[0] nor a prefix of the list is a shortlist. All candidates come back
    # in ONE esummary call however many there are, so checking all of them is one request.
    summaries = json.loads(_get(EUTILS, "/esummary.fcgi",
                                _eutils_params({"db": "clinvar", "id": ",".join(ids),
                                                "retmode": "json"})))["result"]
    # Relevance decides ORDER; order must not decide what KIND of object answers the query.
    # For a well-studied rsID ClinVar ranks combination records first, and a haplotype
    # resolves silently: it hands back one constituent variant's position carrying the
    # haplotype's own classification. This pipeline anchors on ONE point variant.
    ordered = sorted(ids, key=lambda i: _is_combination(summaries.get(i) or {}))
    summ = vset = None
    cid = None
    combos_only = []
    for candidate_id in ordered:
        cand = summaries.get(candidate_id)
        if not cand or not cand.get("variation_set"):
            continue
        cvset = cand["variation_set"][0]
        try:
            _assert_record_matches(query, cand, cvset)
        except ApiError:
            continue
        if _is_combination(cand):
            # Reconciles but is not a point variant: remembered so the refusal below can
            # name what was found rather than claim nothing matched.
            combos_only.append(cand)
            continue
        summ, vset, cid = cand, cvset, candidate_id
        break

    if summ is None and combos_only:
        titles = ", ".join((c.get("title") or c.get("accession") or "?")
                           for c in combos_only[:3])
        raise ApiError(
            f"{query!r} matched only combination records in ClinVar ({len(combos_only)}: "
            f"{titles}), which describe several variants inherited together rather than a "
            f"single one. OriginMarker builds a panel around one point variant, so there "
            f"is nothing here to anchor on. Use the HGVS expression or the VCV accession "
            f"for the individual allele.")

    if summ is None:
        # Nothing ClinVar offered is the variant asked for. An rsID may still be real and
        # simply have no ClinVar submission, so let Ensembl answer before giving up.
        if _RSID_RE.match(query.strip()):
            return _resolve_via_ensembl(query, build=build)
        # Report the SEARCH as having missed, not the user as having mistyped: a full-text
        # search failing to surface a record is not evidence the identifier is wrong. Every
        # id was examined unless esearch had more than retmax, and if it did, say so rather
        # than let "none of them" imply an exhaustive search that did not happen.
        seen = ("" if total <= len(ids)
                else f", of which only the {len(ids)} most relevant were examined")
        raise ApiError(
            f"ClinVar's full-text search returned {total} record(s) for {query!r}{seen}, "
            f"and none of them is this variant. That search is a relevance ranking rather "
            f"than a lookup, so this does not mean the expression is malformed: it means "
            f"ClinVar did not surface a matching record. If you have the rsID or the VCV "
            f"accession, either one resolves directly and is not subject to this.")
    loc38 = next((l for l in vset.get("variation_loc", [])
                  if l.get("assembly_name") == "GRCh38" and l.get("chr")), None)
    if loc38 is None:
        raise ApiError(f"no GRCh38 mapping for {query}")
    loc37 = next((l for l in vset.get("variation_loc", [])
                  if l.get("assembly_name") == "GRCh37" and l.get("chr")), None)
    pos37 = int(loc37["start"]) if loc37 and loc37.get("start") else None
    note = None
    if build == "GRCh37":
        note = (f"Input interpreted as GRCh37; converted to GRCh38 via ClinVar's own "
                f"assembly mapping. GRCh37 chr{loc37['chr']}:{pos37} -> "
                f"GRCh38 chr{loc38['chr']}:{loc38['start']}. All output is GRCh38 (R6)."
                if loc37 else
                "Input flagged GRCh37 but ClinVar has no GRCh37 mapping; "
                "output is GRCh38 (R6).")
    gene, strand = _gene_and_strand(summ, vset)

    # esummary's variation_loc often omits ref/alt; the efetch SequenceLocation (GRCh38,
    # forDisplay) carries positionVCF/referenceAlleleVCF/alternateAlleleVCF.
    sig = rev = acc = None
    vcf_ref = loc38.get("ref") or ""
    vcf_alt = loc38.get("alt") or ""
    pos = int(loc38["start"])
    # Two contracts share this call. The CLASSIFICATION is best-effort: a panel builds
    # fine without ClinVar's verdict. The VCF ALLELES are load-bearing: every downstream
    # step keys on them. So the failure is remembered rather than swallowed, and only the
    # best-effort half may come back empty.
    efetch_err: Optional[str] = None
    try:
        vcv = _get(EUTILS, "/efetch.fcgi",
                   _eutils_params({"db": "clinvar", "rettype": "vcv",
                                   "id": cid, "is_variationid": "true"}))
        root = ET.fromstring(vcv)
        acc = next((e.get("Accession") for e in root.iter("VariationArchive")), None)
        sig, rev = _aggregate_classification(root)
        sl = next((s for s in root.iter("SequenceLocation")
                   if s.get("Assembly") == "GRCh38" and s.get("referenceAlleleVCF")), None)
        if sl is not None:
            vcf_ref = sl.get("referenceAlleleVCF") or vcf_ref
            vcf_alt = sl.get("alternateAlleleVCF") or vcf_alt
            pos = int(sl.get("positionVCF") or pos)
    except Exception as e:  # noqa: BLE001 - urllib and ElementTree raise a zoo of types
        efetch_err = f"{type(e).__name__}: {e}"
        log.warning("clinvar efetch failed for %s: %s", cid, efetch_err)

    # An efetch failure must never be diagnosed as a property of the variant: without this
    # branch, empty alleles fall into the structural-variant test below and a routine SNV
    # is refused for having no alleles.
    if efetch_err and (not vcf_ref or not vcf_alt):
        raise ApiError(
            f"ClinVar's efetch could not be read for VariationID {cid} ({efetch_err}), so "
            f"the reference and alternate alleles for {query!r} are unknown. This is a "
            f"failed lookup, not a statement about the variant. Retry.")

    # A structural variant has no VCF alleles, and everything downstream needs them: the
    # gnomAD lookup keys on <chrom>-<pos>-<ref>-<alt>, the rarity verdict needs an allele
    # to count, and transcript sense needs bases to complement.
    obj = (summ.get("obj_type") or "").lower()
    if not vcf_ref or not vcf_alt or "copy number" in obj or "structural" in obj:
        raise ApiError(
            f"{query!r} is {summ.get('obj_type') or 'a variant'} "
            f"({summ.get('title') or acc or 'no title'}), which has no single reference and "
            f"alternate allele. OriginMarker builds panels around a point variant, so it "
            f"has nothing to anchor on here. Use a variant with explicit VCF alleles.")

    return VariantRecord(
        query=query, rsid=_first_rsid(vset), gene=gene, strand=strand,
        chrom=str(loc38["chr"]), pos_grch38=pos,
        vcf_ref=vcf_ref, vcf_alt=vcf_alt,
        clinical_significance=sig, review_status=rev, clinvar_accession=acc,
        pos_grch37=pos37, build_note=note)


def _aggregate_classification(root) -> tuple[Optional[str], Optional[str]]:
    """ClinVar's VARIATION-level classification and its review status, as a pair.

    A VCV carries one GermlineClassification per RCV, i.e. per variant+condition pairing,
    and those disagree with each other by design, so no document-order scan may pick one.
    The aggregate ClinVar publishes sits at ClassifiedRecord/Classifications, outside
    RCVList. Both values are read off the one node, so the pair describes one record.
    """
    for cr in root.iter("ClassifiedRecord"):
        for cls in cr.findall("Classifications"):
            for kind in ("GermlineClassification", "SomaticClinicalImpact",
                         "OncogenicityClassification"):
                node = cls.find(kind)
                if node is None:
                    continue
                d = node.findtext("Description")
                if d and "removed" not in d.lower():
                    return d, node.findtext("ReviewStatus")
    # No aggregate block: say nothing rather than reach into RCVList and present one
    # condition's opinion as the variant's classification.
    return None, None


def _first_rsid(vset: dict) -> Optional[str]:
    for x in vset.get("variation_xrefs", []) or []:
        if x.get("db_source") == "dbSNP" and x.get("db_id"):
            return "rs" + str(x["db_id"])
    return None


def _resolve_via_ensembl(rsid: str, build: str = "GRCh38") -> VariantRecord:
    """Fallback for an identifier ClinVar does not carry (an rsID with no submission).

    Only reachable for identifier-shaped input, since resolve_variant() gates first.
    "Ensembl says no such variant" and "Ensembl did not answer" are different facts and
    must not collapse into one sentence: only Ensembl's own "not found" may be reported as
    absence. Transient 5xx from rest.ensembl.org is ordinary weather.
    """
    tries = 2
    try:
        rec = json.loads(_get(ENSEMBL, f"/variation/homo_sapiens/{rsid}",
                              {"content-type": "application/json"}, tries=tries))
    except ApiError as e:
        # Ensembl reports an unknown id as 400 with {"error": "<id> not found for
        # homo_sapiens"}, NOT as 404, so 404 alone is the wrong discriminator. Anything
        # that is not Ensembl answering "no" (503, timeout, DNS, a reset) is a question
        # that never got asked, and is reported as one. The status is parsed back out of
        # _http's message because _http raises a bare ApiError.
        m = re.search(r"HTTP Error (\d{3})", str(e))
        status = int(m.group(1)) if m else None
        if status not in (400, 404):
            raise ApiError(
                f"Ensembl could not be reached ({status or 'no response'} after {tries} "
                f"tries), so it was never asked about {rsid!r}. This is NOT a statement "
                f"about whether {rsid!r} exists: the lookup did not happen. Retry in a "
                f"moment.") from None
        raise ApiError(
            f"{rsid!r} was not found in ClinVar, and Ensembl answered that it does not "
            f"know it either. Check the identifier: a transcript version or a base change "
            f"that does not exist will not resolve, and nothing is guessed from a near "
            f"match.") from None
    m = next((m for m in rec.get("mappings", [])
              if m.get("assembly_name") == "GRCh38"), None)
    if m is None:
        raise ApiError(f"{rsid!r} exists but has no GRCh38 mapping, so there is no "
                       f"coordinate to build a panel around.")
    ref, alt = (m["allele_string"].split("/") + ["", ""])[:2]
    note = ("Input flagged GRCh37; resolved via Ensembl to GRCh38 coordinates. "
            "All output is GRCh38 (R6)." if build == "GRCh37" else None)
    return VariantRecord(
        query=rsid, rsid=rec.get("name"), gene=None, strand=None,
        chrom=str(m["seq_region_name"]), pos_grch38=int(m["start"]),
        vcf_ref=ref, vcf_alt=alt, clinical_significance=None,
        review_status=None, clinvar_accession=None, build_note=note)


def _gene_strand(symbol: str) -> Optional[int]:
    try:
        g = json.loads(_get(ENSEMBL, f"/lookup/symbol/homo_sapiens/{symbol}",
                            {"content-type": "application/json"}))
        return int(g.get("strand"))
    except Exception:  # noqa: BLE001
        return None


_TITLE_GENE_RE = re.compile(r"\(([A-Za-z0-9_.\-]{1,32})\)\s*:\s*[cgmnpr]\.", re.I)


def _gene_and_strand(summ: dict, vset: dict) -> tuple[Optional[str], Optional[int]]:
    """The gene a ClinVar record is curated against, and that gene's strand.

    genes[] is ordered ALPHABETICALLY, not by relevance, so genes[0] is not the disease
    gene and may sit on the opposite strand, which would render transcript sense exactly
    backwards while looking well-formed; gene_sort is the same alphabetical key. The HGVS
    title does name the curated gene, so that selects the entry. The strand then comes from
    ClinVar's own field on it (R1: the esummary the caller already fetched), with the
    Ensembl symbol lookup only as a fallback, since a second round-trip is a second thing
    that can fail.
    """
    genes = summ.get("genes") or []
    m = _TITLE_GENE_RE.search(summ.get("title") or vset.get("variation_name") or "")
    want = m.group(1).lower() if m else None
    pick = (next((g for g in genes if (g.get("symbol") or "").lower() == want), None)
            or (genes[0] if genes else None))
    sym = (pick or {}).get("symbol") or summ.get("gene_sort")
    if not sym:
        return None, None
    strand = {"+": 1, "-": -1}.get((pick or {}).get("strand"))
    return sym, (strand if strand is not None else _gene_strand(sym))


# --------------------------------------------------------------------------- #
# Step 2: rarity / LD-usability  (R2 driver)
# --------------------------------------------------------------------------- #

_VARIANT_Q = """query V($id:String!,$ds:DatasetId!){
  variant(variantId:$id,dataset:$ds){ genome{ac an af} exome{ac an af} }}"""

# What "common enough for population LD to be defined" means, stated once per source.
# 1000G phase 3 ALL is 5008 haplotypes, so >5 copies and AF > ~1e-3 are the same claim;
# deriving the second from the first stops the two sources disagreeing about one allele.
# Raise _KG_MIN_AC to be stricter, and the gnomAD side follows automatically.
_KG_HAPLOTYPES = 5008
_KG_MIN_AC = 5
_GNOMAD_MIN_AF = _KG_MIN_AC / _KG_HAPLOTYPES


def assess_rarity(v: VariantRecord) -> Rarity:
    vid = f"{v.chrom}-{v.pos_grch38}-{v.vcf_ref}-{v.vcf_alt}"
    g_af = g_ac = g_an = None
    try:
        j = _graphql(_VARIANT_Q, {"id": vid, "ds": "gnomad_r4"})
        gv = (j.get("data") or {}).get("variant") or {}
        gen = gv.get("genome") or {}
        g_af, g_ac, g_an = gen.get("af"), gen.get("ac"), gen.get("an")
    except Exception:  # noqa: BLE001
        pass

    kg_ac = None
    if v.rsid:
        try:
            rec = json.loads(_get(ENSEMBL, f"/variation/homo_sapiens/{v.rsid}",
                                  {"pops": "1", "content-type": "application/json"}))
            kg_all = [p for p in rec.get("populations", [])
                      if p.get("population") == "1000GENOMES:phase_3:ALL"]
            alt_entries = [p for p in kg_all if p.get("allele") == v.vcf_alt]
            if alt_entries:
                kg_ac = min(p.get("allele_count") for p in alt_entries)
            elif kg_all:
                # The panel HAS this variant, but none of its alleles string-equals our VCF
                # alt: a failed join, not a count of zero. Systematic for indels, where we
                # say C>CA and Ensembl says '-' and 'A'. None means unknown.
                kg_ac = None
            else:
                # Ensembl knows the rsID but carries no 1000G ALL frequencies for it:
                # genuinely absent from the panel.
                kg_ac = 0
        except Exception:  # noqa: BLE001
            pass

    # Displayed and exported, so format rather than letting a float repr leak in beside a
    # cell rendering the same number as 1.77e-4. Neither figure may be interpolated raw:
    # both are Optional, and a None rendered into prose reads as a measurement.
    af_txt = fmt_af(g_af)
    ev = (f"1000G allele count {kg_ac if kg_ac is not None else 'unavailable'}, "
          f"gnomAD genome AF {af_txt}")

    # Three states, not two: a source that FAILED TO ANSWER is not a source that answered
    # "rare". Every branch below must be reached from evidence that EXISTS, never from a
    # gap, and neither may a missing source alone force "unknown", which would drop a true
    # rarity finding on an allele only gnomAD could answer for. Note the 1000G join drops
    # out deterministically on indel notation, so this is not only an outage path.
    kg_common = kg_ac is not None and kg_ac > _KG_MIN_AC
    g_common = g_af is not None and g_af > _GNOMAD_MIN_AF
    said_rare = ((kg_ac is not None and not kg_common)
                 or (g_af is not None and not g_common))

    # 1000G-only and conservative: gnomAD calling an allele common is good reason to
    # believe LD estimates exist but is not confirmation, and app/ldlink.py hard-gates on
    # this flag. Unconfirmed must read as "unknown", not "yes".
    usable = kg_common
    if kg_common:
        ld_status = "defined"
        reason = (f"Common enough in reference panels ({ev}) for population-LD estimates "
                  f"to exist. They are not used here: markers are ranked by expected "
                  f"heterozygosity and proximity, and origin still comes from per-family "
                  f"phasing.")
    elif g_common:
        ld_status = "unknown"
        reason = (f"gnomAD reports this allele as common ({ev}), so population-LD "
                  f"estimates for it most likely exist, but the 1000 Genomes count did "
                  f"not confirm that, so it is not asserted here. This is not a rarity "
                  f"finding. LD is not used for ranking or origin either way, and origin "
                  f"comes from per-family phasing.")
    elif said_rare:
        ld_status = "undefined"
        reason = (f"Too rare in reference panels for linkage disequilibrium to be defined "
                  f"({ev}). Per-family phasing required.")
    else:
        ld_status = "unknown"
        reason = ("Frequency data could not be retrieved for this allele, so its rarity is "
                  "unknown and population LD cannot be evaluated either way. Treated as "
                  "unusable for LD, and per-family phasing is required regardless. This is "
                  "a missing-data verdict, not a rarity finding: retry, or check the "
                  "variant in gnomAD directly.")
    return Rarity(g_af, g_ac, g_an, kg_ac, usable, reason, ld_status)


# --------------------------------------------------------------------------- #
# Step 3: enumerate + annotate candidates (concurrent region pull)
# --------------------------------------------------------------------------- #

_REGION_Q = """query R($chrom:String!,$start:Int!,$stop:Int!,$ds:DatasetId!){
  region(chrom:$chrom,start:$start,stop:$stop,reference_genome:GRCh38){
    variants(dataset:$ds){ variant_id rsid pos ref alt genome{af an filters populations{id ac an}} }
  }
}"""


def _fetch_region_chunk(chrom: str, start: int, stop: int) -> list[dict]:
    """One gnomAD region slice. Raises on failure rather than returning nothing.

    Returning [] is indistinguishable from "this slice genuinely has no variants", so an
    outage would produce a smaller panel instead of an error: markers missing, coverage
    computed over a hole. A panel built on a partial pull is unsound, so the build fails.
    _http has already retried with backoff before we get here.
    """
    j = _graphql(_REGION_Q, {"chrom": chrom, "start": start, "stop": stop, "ds": "gnomad_r4"})
    if j.get("errors"):
        msg = "; ".join(str(e.get("message", e)) for e in j["errors"])[:200]
        raise ApiError(f"gnomAD rejected the region chr{chrom}:{start}-{stop}: {msg}")
    region = (j.get("data") or {}).get("region")
    if region is None:
        raise ApiError(f"gnomAD returned no region object for chr{chrom}:{start}-{stop}; "
                       f"the pull is incomplete, so no panel is built on it.")
    return region.get("variants") or []


def _contig_length(chrom: str) -> Optional[int]:
    """GRCh38 length of one chromosome, live from Ensembl (R1: never from memory).

    Advisory, deliberately not fatal: gnomAD accepts a stop past the telomere and answers
    with an empty region rather than an error, so failing to clamp hi costs a footnote, not
    a coordinate. Returns None on any failure, which leaves the span marked unverified
    rather than failing a build that is otherwise sound.
    """
    try:
        j = json.loads(_get(ENSEMBL, f"/info/assembly/homo_sapiens/{chrom}",
                            {"content-type": "application/json"}))
        # Only GRCh38's length may clamp a GRCh38 coordinate (R6): a wrong-assembly length
        # would not error, it would restate the telomere, and this reply both moves hi and
        # gets printed as fact. Mismatch takes the advisory path, same as an outage.
        return int(j["length"]) if j.get("assembly_name") == "GRCh38" else None
    except (ApiError, KeyError, ValueError, TypeError):
        return None


def enumerate_candidates(v: VariantRecord, window: int = DEFAULT_WINDOW,
                         on_progress: Optional[Callable[[int, int], None]] = None
                         ) -> dict[str, dict]:
    # A flank must not run off the end of a chromosome: gnomAD answers a negative start
    # with an HTTP 500. Clamping alone is not enough, though, since a clamped flank is
    # SHORTER than the window asked for, and R5 exists to stop a short flank shipping
    # labelled as a full one. So the span actually queried per side is recorded, and any
    # side that got less than it asked for says so.
    clen = _contig_length(v.chrom)
    lo = max(v.pos_grch38 - window, 1)
    hi = min(v.pos_grch38 + window, clen) if clen else v.pos_grch38 + window
    _note_span(v.pos_grch38, lo, hi, window, v.chrom, clen)
    chunks = [(s, min(s + REGION_CHUNK, hi)) for s in range(lo, hi, REGION_CHUNK)]
    out: dict[str, dict] = {}
    done = 0
    with cf.ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        # copy_context() must be evaluated on the submitting thread: pool workers start with
        # an empty context and would tally into a ledger nobody reads.
        futs = {ex.submit(contextvars.copy_context().run, _fetch_region_chunk, v.chrom, s, e):
                (s, e) for s, e in chunks}
        for fut in cf.as_completed(futs):
            for var in fut.result():
                out[var["variant_id"]] = var
            done += 1
            if on_progress:
                on_progress(done, len(chunks))
    return out


def _is_snp(var: dict) -> bool:
    return len(var["ref"]) == 1 and len(var["alt"]) == 1


def _per_pop_maf(var: dict) -> tuple[dict, dict]:
    """(per-population MAF, per-population AN), both keyed by ancestry code.

    The AN floor applies per population, not just globally: CALL_RATE_AN_FLOOR says the
    SITE is well called, not that a given population was sampled there, and the ancestry
    ranking can select on a MAF drawn from four people. AN is returned alongside so a
    reader can see what each figure rests on.
    """
    g = var.get("genome") or {}
    maf, ans = {}, {}
    for p in (g.get("populations") or []):
        pid = p["id"]
        if "_" in pid or pid not in GNOMAD_POPS:   # skip sex/subset splits
            continue
        an, ac = p.get("an") or 0, p.get("ac") or 0
        if an < MIN_POP_AN:
            continue
        a = ac / an
        maf[GNOMAD_POPS[pid]] = round(min(a, 1 - a), 4)
        ans[GNOMAD_POPS[pid]] = an
    return maf, ans


def _tier(dist: int) -> str:
    a = abs(dist)
    if a < 2_000:
        return "A_core(<2kb)"
    if a < 30_000:
        return "B_near(2-30kb)"
    return "C_flank(30kb+)"


def annotate(variants: dict[str, dict], v: VariantRecord,
             common_maf: float = COMMON_MAF) -> list[Marker]:
    markers: list[Marker] = []
    for var in variants.values():
        if not _is_snp(var):
            continue
        g = var.get("genome") or {}
        af = g.get("af")
        if af is None or (g.get("an") or 0) < CALL_RATE_AN_FLOOR:
            continue
        # gnomAD's own QC verdict: a non-empty filters list means gnomAD failed the site
        # (AC0, or an AS_VQSR artifact). These cannot merely be tolerated: a mismapped site
        # collects spurious heterozygotes, which inflates 2pq, and 2pq is the ranking key,
        # so the sort is drawn toward exactly the sites gnomAD already rejected.
        if g.get("filters"):
            continue
        maf = min(af, 1 - af)
        if maf < common_maf:
            continue
        dist = var["pos"] - v.pos_grch38
        # A variant does not flank itself: a COMMON pathogenic allele clears the MAF floor
        # and is otherwise recruited into its own panel at +0 bp, where it has no side and
        # is counted on neither, breaking the R5 arithmetic.
        if dist == 0:
            continue

        pop, pop_an = _per_pop_maf(var)
        het = 2 * af * (1 - af)
        het_pop = max((2 * a * (1 - a) for a in pop.values()), default=het)
        gm = genetic_map.annotate_distance(v.chrom, v.pos_grch38, var["pos"])
        markers.append(Marker(
            rsid=var.get("rsid") or var["variant_id"], variant_id=var["variant_id"],
            chrom=v.chrom, pos=var["pos"], ref=var["ref"], alt=var["alt"],
            af=round(af, 4), maf=round(maf, 4), het=round(het, 4),
            het_max_pop=round(het_pop, 4), dist=dist,
            an=g.get("an"), per_pop_an=pop_an,
            # Named for the axis, not the arm: a chromosome runs telomere, p, CENTROMERE,
            # q, telomere, so a lower coordinate is toward the telomere on a p-arm and
            # toward the centromere on a q-arm. Nothing here knows where the centromere
            # is, so nothing here says.
            side="higher coord" if dist > 0 else "lower coord",
            tier=_tier(dist), per_pop_maf=pop,
            cm=gm["cm"], recomb_fraction=gm["recomb_fraction"],
            hotspot_between=gm["hotspot_between"], map_approx=gm["map_approx"]))
    return markers


# --------------------------------------------------------------------------- #
# Step 4: rank + select balanced panel (R2, R5) + cross-check (R1)
# --------------------------------------------------------------------------- #

def _rank_key(ancestry: Optional[str]):
    """Sort key: 2pq for the ancestry asked about, else the global 2pq. Never het_max_pop.

    het_max_pop is the MAXIMUM 2pq across the eight gnomAD populations: an order statistic,
    biased upward, with each marker's figure drawn from a different population. Sorting on
    it assumes the family belongs to whichever population flatters each row. The global 2pq
    is computed on the whole cohort, and is also the fallback when a selected ancestry has
    no frequency here, so an EAS query is never ordered on some other population's figure.
    """
    def key(m: Marker):
        p = m.per_pop_maf.get(ancestry) if ancestry else None
        het = 2 * p * (1 - p) if p is not None else m.het
        # R2: het then proximity, never LD. pos and variant_id are not genetics; they make
        # the order total, since het saturates at 0.5 for many sites and two markers can be
        # equidistant on opposite sides. Without a final tiebreak, a stable sort settles an
        # exact tie by gnomAD's response order, and two labs get different panels.
        return (-het, -m.het, abs(m.dist), m.pos, m.variant_id)
    return key


def _ranking_key_label(ancestry: Optional[str]) -> str:
    """Names the quantity that ACTUALLY produced the sort order. Exports render it
    verbatim and must never restate the ranking basis in their own words."""
    if not ancestry:
        return "global 2pq (het)"
    return f"2pq in {ancestry}, global 2pq where {ancestry} has no frequency"


def select_panel(markers: list[Marker], ancestry: Optional[str] = None,
                 per_bin: int = 2) -> tuple[list[Marker], list[Marker], dict]:
    # 2pq = 0 cannot say which parental haplotype an embryo inherited, so it is not a
    # marker. common_maf=0.0 is a legal floor, and would otherwise let monomorphic sites
    # fill a band and satisfy the coverage rule below carrying no information.
    usable = [m for m in markers if max(m.het, m.het_max_pop) > 0]
    ranked = sorted(usable, key=_rank_key(ancestry))

    # One marker per site: a genotyping assay reads the site once, and multi-allelic sites
    # share an rsID, so this must key on pos rather than rsid. ranked is best-first, so
    # setdefault keeps the better alt.
    by_site: dict[int, Marker] = {}
    for m in ranked:
        by_site.setdefault(m.pos, m)
    pool = list(by_site.values())

    bands = [(0, 2_000), (2_000, 10_000), (10_000, 30_000),
             (30_000, 80_000), (80_000, 10**12)]
    recommended: list[Marker] = []
    for want_higher in (False, True):
        for lo, hi in bands:
            recommended += [m for m in pool if (m.dist > 0) == want_higher
                            and lo <= abs(m.dist) < hi][:per_bin]
    recommended.sort(key=lambda m: m.dist)   # total: one marker per site, so dist is unique

    lower = [m for m in recommended if m.dist < 0]
    higher = [m for m in recommended if m.dist > 0]
    coverage = {
        "lower_count": len(lower), "higher_count": len(higher),
        "lower_core_near": sum(1 for m in lower if abs(m.dist) < 30_000),
        "higher_core_near": sum(1 for m in higher if abs(m.dist) < 30_000),
        "flags": [],
    }
    # The 1 cM/Mb fallback can never reach the hotspot threshold, so its False means "not
    # assessed", which is not "no hotspot" and must never be counted as one.
    def assessed(m: Marker) -> bool:
        return m.hotspot_between is not None and m.map_approx is False

    for side, label in ((lower, "lower-coordinate"), (higher, "higher-coordinate")):
        near = [m for m in side if abs(m.dist) < 30_000]
        judged = [m for m in near if assessed(m)]
        clear = [m for m in judged if not m.hotspot_between]
        if not side:
            coverage["flags"].append(
                f"No markers at all on the {label} side: a recombination between the "
                f"variant and the panel cannot be detected on that side.")
        elif len(near) < 2:
            coverage["flags"].append(
                f"Fewer than 2 markers within 30 kb on the {label} side.")
        elif judged and len(clear) < 2:
            coverage["flags"].append(
                f"Fewer than 2 markers within 30 kb on the {label} side are clear of an "
                f"intervening recombination hotspot ({len(clear)} of {len(judged)}).")

    hot = [m for m in recommended if assessed(m) and m.hotspot_between]
    if hot:
        coverage["flags"].append(
            f"{len(hot)} of the {len(recommended)} shortlisted markers have a "
            f"recombination hotspot between them and the variant, so they are the most "
            f"likely to have lost phase with it.")
    unjudged = [m for m in recommended if not assessed(m)]
    if unjudged:
        coverage["flags"].append(
            f"Recombination hotspots were not assessed for {len(unjudged)} of the "
            f"{len(recommended)} shortlisted markers.")
    return ranked, recommended, coverage


def cross_check_ensembl(markers: list[Marker], top_n: int = 8) -> None:
    """Verify top markers' GRCh38 positions against Ensembl (independent source) (R1)."""
    subset = sorted(markers, key=lambda m: abs(m.dist))[:top_n]
    ids = [m.rsid for m in subset if m.rsid.startswith("rs")]
    if not ids:
        return
    body = json.dumps({"ids": ids}).encode()
    try:
        resp = json.loads(_http("POST", ENSEMBL + "/variation/homo_sapiens",
                                body=body,
                                headers={"Content-Type": "application/json"}))
    except Exception:  # noqa: BLE001
        return
    epos = {}
    for rid, rec in resp.items():
        m = next((mm for mm in rec.get("mappings", [])
                  if mm.get("assembly_name") == "GRCh38"), None)
        if m:
            epos[rid] = int(m["start"])
    for m in subset:
        if m.rsid in epos:
            m.ensembl_pos_check = "ok" if epos[m.rsid] == m.pos else f"MISMATCH:{epos[m.rsid]}"


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

# Required wording (R8), reproduced verbatim in the footer and in every export.
DISCLAIMER = ("Research use only. Candidate markers require validation and per-family "
              "phasing in a qualified genetics laboratory. Not a clinical diagnostic.")

# Fixed Layer-B protocol (R3): wet-lab work the app cannot do. Surfaced in the UI drawer
# and stamped into every export.
LAYER_B_STEPS = [
    "Genotype the carrier parent at these candidate markers.",
    "Keep only the markers where that carrier is actually heterozygous. Expected "
    "heterozygosity is a population average, not this individual's genotype.",
    "Phase the retained markers against an informative relative (an affected or "
    "unaffected child, a proband, or a grandparent), or by read-based sequencing. "
    "OriginMarker cannot determine phase.",
    "Genotype the embryo biopsy at the phased markers.",
    "Require at least two concordant markers per side before calling parental origin.",
    "Use markers on both sides: one side alone cannot reveal a recombination between "
    "the marker and the locus, and per-side redundancy guards against allele dropout.",
]


def build(query: Union[str, "StructuredQuery"], window: int = DEFAULT_WINDOW,
          ancestry: Optional[str] = None, common_maf: float = COMMON_MAF,
          cross_check: bool = True,
          on_progress: Optional[Callable[[str, float], None]] = None) -> PanelResult:
    """Full Layer-A build. Returns a PanelResult (R1-R7).

    `query` is either an HGVS/rsID string or a StructuredQuery, which supplies every other
    parameter and takes precedence over the keyword args. `on_progress(stage, fraction)` is
    called as the build advances, for a web layer to stream progress from.
    """
    if isinstance(query, StructuredQuery):
        q = query
    else:
        q = StructuredQuery(variant=query, window_bp=window, ancestry=ancestry,
                            common_maf=common_maf, cross_check=cross_check)

    def progress(stage: str, frac: float) -> None:
        if on_progress:
            on_progress(stage, frac)

    t0 = time.time()
    _reset_fetch_log()          # provenance reports the age of THIS build's data
    progress("resolving variant", 0.05)
    v = resolve_variant(q.variant, build=q.build)     # R1, R6

    progress("confirming rarity", 0.15)
    rarity = assess_rarity(v)                         # R2 driver

    progress("pulling population variants", 0.25)
    raw = enumerate_candidates(v, q.window_bp,
                               on_progress=lambda done, total: progress(
                                   f"pulling variants ({done}/{total} regions)",
                                   0.25 + 0.45 * (done / max(total, 1))))

    progress(f"annotating {len(raw)} variants", 0.75)
    markers = annotate(raw, v, q.common_maf)          # R4

    progress("ranking + selecting panel", 0.85)
    ranked, recommended, coverage = select_panel(markers, q.ancestry)  # R2, R5

    if q.cross_check:
        progress("cross-checking positions against Ensembl", 0.92)
        cross_check_ensembl(recommended)              # R1

        # The cross-check runs AFTER selection, since it costs a request per marker and is
        # spent on the shortlist rather than on the whole pool. So a marker whose sources
        # disagree about where it IS has already been shortlisted by the time we find out.
        # That is not a cell value, it is a reason to distrust the row: it goes to
        # coverage.flags, which renders in the UI alert and in every export.
        bad = [m for m in recommended
               if (m.ensembl_pos_check or "").startswith("MISMATCH")]
        if bad:
            coverage.setdefault("flags", []).append(
                f"{len(bad)} of the {len(recommended)} shortlisted markers are placed "
                f"differently by gnomAD and Ensembl "
                f"({', '.join(f'{m.rsid} at {m.pos:,} vs {m.ensembl_pos_check.split(':')[1]}' for m in bad[:3])}"
                f"{', and others' if len(bad) > 3 else ''}). A marker its own sources "
                f"cannot place is not safe to genotype against: confirm the position "
                f"before ordering it, or drop it.")

    # R5: a clamped flank is a COVERAGE fact, so it belongs on the coverage card beside the
    # per-side counts, which is where a reader checks whether both sides are covered.
    # enumerate_candidates detects and words the truncation; this is the pickup.
    span = dict(_span())
    coverage.setdefault("flags", []).extend(span.get("flags", []))
    # One read of the ledger, so the data date and the two counts below describe the same
    # snapshot rather than three reads a concurrent fetch could slide between.
    _flog = dict(_fetch_log())
    prov = {
        # What was actually asked of gnomAD; window_bp is what was requested. They differ
        # exactly when a flank is clamped.
        "queried_span": span,
        "ensembl_release": ensembl_release(),
        "sources": {"clinvar": "NCBI E-utilities", "ensembl": ENSEMBL,
                    "gnomad": "v4 (gnomad_r4) GraphQL",
                    "genetic_map": genetic_map.load(v.chrom).source},
        "build": "GRCh38", "window_bp": q.window_bp, "common_maf": q.common_maf,
        "ancestry_rank": q.ancestry, "candidate_n": len(markers),
        "requested_build": q.build,
        "ranking_key": _ranking_key_label(q.ancestry),
        "queried_utc": _utc(_flog["oldest"]) or _utc(time.time()),
        "built_utc": _utc(time.time()),
        "source_responses_from_cache": _flog["from_cache"],
        "source_responses_from_network": _flog["from_network"],
        "elapsed_s": round(time.time() - t0, 1),
        "disclaimer": DISCLAIMER,
        "layer_b_steps": LAYER_B_STEPS,
        # Nothing above was decided by these. Both keys are always present so a consumer
        # can read one shape; nl_model None means the user named the variant themselves,
        # and every render must then stay silent rather than say "none".
        "nl_text": q.nl_text,
        "nl_model": q.nl_model,
    }
    progress("done", 1.0)
    return PanelResult(v, rarity, ranked, recommended, coverage,
                       {"window": q.window_bp, "ancestry": q.ancestry,
                        "common_maf": q.common_maf, "build": q.build},
                       prov)


if __name__ == "__main__":
    import sys

    # Ranking key, offline. Asserts the key VALUE rather than an emergent ordering: which
    # quantity came out is the thing that has to be pinned.
    _flattered = Marker(                      # low globally, high in exactly one group
        rsid="rs_flattered", variant_id="1-100-A-G", chrom="1", pos=100, ref="A", alt="G",
        af=0.18, maf=0.18, het=0.2952, het_max_pop=0.5000, dist=100, side="higher coord",
        tier=_tier(100), per_pop_maf={"AFR": 0.5, "EAS": 0.01})
    _honest = Marker(                         # high globally, no population flatters it
        rsid="rs_honest", variant_id="1-200-A-G", chrom="1", pos=200, ref="A", alt="G",
        af=0.43, maf=0.43, het=0.4902, het_max_pop=0.4950, dist=200, side="higher coord",
        tier=_tier(200), per_pop_maf={"AFR": 0.45, "NFE": 0.44})

    _k = _rank_key(None)
    assert _k(_flattered)[0] == -0.2952, _k(_flattered)   # global 2pq, NOT het_max_pop 0.5
    assert _k(_honest)[0] == -0.4902, _k(_honest)
    assert sorted([_flattered, _honest], key=_k)[0].rsid == "rs_honest"

    _k = _rank_key("EAS")
    assert abs(_k(_flattered)[0] + 2 * 0.01 * 0.99) < 1e-12, _k(_flattered)  # EAS's number
    # _honest has no EAS frequency, so the fallback is the global 2pq, not het_max_pop.
    assert _k(_honest)[0] == -0.4902, _k(_honest)

    assert _ranking_key_label(None) == "global 2pq (het)"
    assert _ranking_key_label("EAS").startswith("2pq in EAS")

    # Cache poisoning, offline. Points at a dead port so the refetch is guaranteed to fail:
    # the assertion is that the poison is DROPPED rather than served, not that a retry
    # succeeds.
    import tempfile
    _saved_dir = CACHE_DIR
    CACHE_DIR = Path(tempfile.mkdtemp())
    _dead = "http://127.0.0.1:9/efetch.fcgi?db=clinvar&id=9088"
    _ck = _cache_key("GET", _dead, None)
    _ck.parent.mkdir(parents=True, exist_ok=True)
    _poison = "<html><head><title>502 Bad Gateway</title></head><body>nginx</body></html>"
    _ck.write_text(_poison)
    try:
        _served = _http("GET", _dead, tries=1, timeout=1)
    except ApiError:
        _served = None                      # the refetch failed, which is the correct end
    assert _served != _poison, "a cached error page was served as if it were ClinVar"
    assert not _ck.exists(), "poisoned cache entry survived a read"
    # The 502 page above parses cleanly under ET.fromstring, which is why the check is on
    # the raw text.
    assert not _body_is_expected(_dead, _poison)
    assert _body_is_expected(_dead, "<?xml version='1.0'?><VariationArchive/>")
    assert _body_is_expected(GNOMAD, '{"data":{}}')
    assert not _body_is_expected(GNOMAD, "<html>502</html>")
    CACHE_DIR = _saved_dir

    # Telomere clamp, offline via stubs. The stubbed coordinates below are the only ones in
    # this file that did not come from an API this run: they are here so the ARITHMETIC can
    # be checked without a network, and nothing downstream ever sees them (R1).
    _asked: list = []
    _saved_chunk, _saved_len = _fetch_region_chunk, _contig_length
    _fetch_region_chunk = lambda chrom, s, e: (_asked.append((s, e)), [])[1]  # noqa: E731
    _contig_length = lambda chrom: 90_338_345          # chr16, as Ensembl reports it
    _hba2 = VariantRecord(query="HBA2-synthetic", rsid=None, gene="HBA2", strand=1,
                          chrom="16", pos_grch38=172_021, vcf_ref="G", vcf_alt="A",
                          clinical_significance=None, review_status=None,
                          clinvar_accession=None)
    enumerate_candidates(_hba2, window=250_000)
    assert min(s for s, _ in _asked) >= 1, f"negative coordinate sent to gnomAD: {_asked}"
    assert _span()["left_bp"] == 172_020, _span()
    # ...and the short flank is DISCLOSED, not merely survived.
    assert any("Left flank truncated" in f for f in _span()["flags"]), _span()
    assert not any("Right flank" in f for f in _span()["flags"]), _span()

    _asked.clear()
    _qter = VariantRecord(query="qter-synthetic", rsid=None, gene=None, strand=1,
                          chrom="16", pos_grch38=90_300_000, vcf_ref="G", vcf_alt="A",
                          clinical_significance=None, review_status=None,
                          clinvar_accession=None)
    enumerate_candidates(_qter, window=250_000)
    assert max(e for _, e in _asked) <= 90_338_345, f"queried past the telomere: {_asked}"
    assert _span()["right_bp"] == 38_345, _span()
    assert any("Right flank truncated" in f for f in _span()["flags"]), _span()
    _fetch_region_chunk, _contig_length = _saved_chunk, _saved_len

    # ...and the length doing the clamping has to be GRCh38's (R6). The stub lengths are
    # synthetic on purpose: what is asserted is that a wrong-assembly reply is REFUSED, not
    # that any particular number is right, so no remembered coordinate goes in.
    _saved_get = _get
    _get = lambda *a, **kw: json.dumps({"length": 1, "assembly_name": "GRCh37"})  # noqa: E731
    assert _contig_length("16") is None, "a non-GRCh38 length was accepted as a flank clamp"
    _get = lambda *a, **kw: json.dumps({"length": 12_345, "assembly_name": "GRCh38"})  # noqa: E731
    assert _contig_length("16") == 12_345
    _get = _saved_get

    # ...and the truncation reaches a READER, not just _span(): a flag nothing picks up
    # discloses nothing, and asserting on _span() alone passes in exactly that world. So
    # this asserts on what build() hands out, and specifically on the R5 coverage card,
    # which is where "is each side covered" is answered.
    _mid = VariantRecord(query="mid-synthetic", rsid=None, gene="MID", strand=1,
                         chrom="16", pos_grch38=45_000_000, vcf_ref="G", vcf_alt="A",
                         clinical_significance=None, review_status=None,
                         clinvar_accession=None)
    _recs = {"HBA2-synthetic": _hba2, "mid-synthetic": _mid}
    _saved_res, _saved_rar = resolve_variant, assess_rarity
    _saved_chunk2, _saved_len2 = _fetch_region_chunk, _contig_length
    resolve_variant = lambda q, **kw: _recs[q]                        # noqa: E731
    assess_rarity = lambda v: Rarity(None, None, None, None, False, "stub", "unknown")  # noqa: E731
    _fetch_region_chunk = lambda chrom, s, e: (time.sleep(0.02), [])[1]   # noqa: E731
    _contig_length = lambda chrom: 90_338_345                         # noqa: E731
    try:
        _clamped = build("HBA2-synthetic", window=250_000, cross_check=False)
        _cov_flags = _clamped.coverage.get("flags", [])
        assert any("Left flank truncated" in f for f in _cov_flags), _cov_flags
        assert _clamped.provenance["queried_span"]["left_bp"] == 172_020, \
            _clamped.provenance["queried_span"]
        # The disclosure has to CONTRADICT window_bp, or it discloses nothing.
        assert _clamped.provenance["window_bp"] == 250_000
        assert _clamped.provenance["queried_span"]["left_bp"] < 250_000

        # ...and it survives a CONCURRENT build: app/jobs.py runs MAX_CONCURRENT builds in
        # one process, and a serial assertion cannot see one build clobbering another's
        # span, which is why this one runs two at once.
        _out: dict = {}
        def _go(name: str) -> None:
            r = build(name, window=250_000, cross_check=False)
            _out[name] = r.provenance["queried_span"]["left_bp"]
        _threads = [threading.Thread(target=_go, args=(n,)) for n in _recs]
        for t in _threads:
            t.start()
        for t in _threads:
            t.join()
        assert _out["HBA2-synthetic"] == 172_020, f"span raced: {_out}"   # not 250,000
        assert _out["mid-synthetic"] == 250_000, f"span raced: {_out}"    # no phantom clamp
    finally:
        resolve_variant, assess_rarity = _saved_res, _saved_rar
        _fetch_region_chunk, _contig_length = _saved_chunk2, _saved_len2

    q = sys.argv[1] if len(sys.argv) > 1 else "NM_000352.6(ABCC8):c.3989-9G>A"
    r = build(q)
    assert r.provenance["ranking_key"] == _ranking_key_label(None), r.provenance["ranking_key"]

    # The nl_* fields are provenance and nothing else: they must reach prov and leave the
    # resolved variant and the selected panel identical. A build that read them diverges here.
    assert r.provenance["nl_text"] is None and r.provenance["nl_model"] is None
    _nl = build(StructuredQuery(variant=q, nl_text="prose a user typed",
                                nl_model="claude-test-model-1"))
    assert _nl.provenance["nl_text"] == "prose a user typed", _nl.provenance["nl_text"]
    assert _nl.provenance["nl_model"] == "claude-test-model-1", _nl.provenance["nl_model"]
    assert ([m.variant_id for m in _nl.recommended]
            == [m.variant_id for m in r.recommended]), "build() branched on nl_* provenance"
    assert _nl.variant.pos_grch38 == r.variant.pos_grch38, "build() branched on nl_* provenance"
    print(f"{r.variant.rsid}  {r.variant.gene}  chr{r.variant.chrom}:{r.variant.pos_grch38}"
          f"  {r.variant.vcf_ref}>{r.variant.vcf_alt}  {r.variant.clinical_significance}")
    print("LD usable:", r.rarity.population_LD_usable, "|", r.rarity.reason)
    print("candidates:", len(r.candidates), "recommended:", len(r.recommended),
          "| coverage flags:", r.coverage["flags"] or "none")
    print("nearest recommended:")
    for m in sorted(r.recommended, key=lambda m: abs(m.dist))[:6]:
        print(f"  {m.rsid:12} {m.pos:>10} {m.dist:>+7} bp  2pq={m.het:.2f}"
              f"  {m.side:20} chk={m.ensembl_pos_check}")

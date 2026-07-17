"""Optional in-silico PCR verification of a primer pair against UCSC (hg38, R6 explicit).

Primer3 designs against an 800bp template and cannot see the other 3.1Gb. A pair that is
perfect on that template may still prime a second locus: the trace is then a mixture of
loci, and a heterozygote at the marker is indistinguishable from two loci differing at
that base. This module is the only thing here that looks at the rest of the genome.

Three things it is never allowed to be:

  - a claim of validation. One product in silico against a reference is not a wet-lab
    result, and PASS says that much and no more.
  - a reason to hide a primer. Every finding keeps the primer and warns. A hidden primer
    is a decision made for the user; a warned primer is information.
  - a pass by default. No key, an unreadable page, a CAPTCHA, a timeout and a quota stop
    are all NOT VERIFIED, never "clean". Only a parsed single product, on the expected
    chromosome, at the expected size, is a pass.

UCSC publishes one hit per 15 seconds and 5000 per day for program-driven use. The gate
below is module-level and process-wide because that limit is UCSC's, not any one caller's:
11 pairs therefore take about 165 seconds. That is the module working, not hanging.

Never bypass the CAPTCHA. The API key is UCSC's sanctioned programmatic route; without one
hgPcr answers with a Turnstile challenge, which classifies as NOT VERIFIED and must stay
that way.

SERVER SIDE ONLY: UCSC_API_KEY never reaches the browser.

Self-check:  PYTHONPATH=<repo root> python app/ispcr.py
"""

from __future__ import annotations

import html
import os
import re
import threading
import time
import urllib.parse
import urllib.request
from typing import Optional

import panelbuilder as pb

UCSC = "https://genome.ucsc.edu"

# R6: hg38 stated on every request. UCSC's default db is not ours to inherit.
DB = "hg38"

# UCSC's published program-driven limits. Not tuning knobs: raising either one is a
# decision about someone else's server, so they are not read from the environment.
MIN_INTERVAL_S = 15.0
MAX_PER_DAY = 5000

# UCSC's form minimum for both match thresholds.
PERFECT = 15

# hgPcr's own default for Max Product Size, and the width of the question this module asks.
#
# It is NOT a filter on the answer. Measured against the live endpoint with a pair known to
# give exactly one 549bp product on chr11: wp_size=400 returns ZERO products, 1000 and 4000
# both return the one. So gfPcr takes this as a bound on the search, and a product wider
# than it is never generated to be reported.
#
# That makes a low value the quietest way to break this module. At the 1000 this once
# defaulted to, a pair whose second locus amplifies at 1001-4000bp came back holding one
# product, classified ONE_PRODUCT, and printed VERIFIED CLEAN (in silico) on a filed PDF:
# the multi-locus pass the module exists to prevent, reached by asking UCSC a narrower
# question than UCSC asks itself and reporting the answer as genome-wide. It cut the other
# way too, since the design's own max_product goes to 3000: a 2500bp product could not be
# reported at all, and "found no product, do not order" was our cap talking, not UCSC.
#
# 4000 does not make the check exhaustive: it moves the blind window to 4001+bp. It makes
# the window UCSC's own rather than one this file invented.
REPORT_MAX_BP = 4000

# The four states. Only ONE_PRODUCT is a pass; the other three are all "not verified",
# and they are kept apart because their remedies differ: configure a key, investigate the
# pair, or retry later.
NOT_CHECKED = "not_checked"     # never asked
ONE_PRODUCT = "one_product"     # asked, and the answer was clean
DANGER = "danger"               # asked, and the answer was bad
UNKNOWN = "unknown"             # asked, and the answer could not be read

CAVEAT = (
    "IN-SILICO ONLY, NOT A WET-LAB VALIDATION. UCSC In-Silico PCR aligns the two primers "
    "against the GRCh38 reference and reports where they would amplify it. It is not a "
    "PCR. It does not see this carrier's genome: a private variant under a primer site "
    "will not appear here and can still cause allele dropout. It does not model the "
    "cycling conditions, so a single reference product is not a promise of a single band "
    "on a gel. A clean result here means the pair is not obviously multi-locus against "
    "one reference sequence, and nothing more. Confirm every pair in the lab before use."
)

# The same caveat in one sentence. It rides every verdict, a pass included, so on the page
# it is the line under every checked pair: at CAVEAT's length nobody reads the fifth one.
CAVEAT_SHORT = (
    "In silico against the reference only, not a wet-lab validation, and blind to this "
    "carrier's own variants."
)

NOT_CHECKED_NOTE = (
    "Not checked against the genome. This pair was designed against an 800bp Ensembl "
    "template and has not been tested for products anywhere else in GRCh38. It may "
    "amplify more than one locus. UCSC In-Silico PCR verification is optional and is not "
    "configured on this server."
)

NOT_CHECKED_SHORT = (
    "Not checked against the genome: verification is not configured on this server."
)


def available() -> bool:
    """True if a UCSC API key is configured. Verification is optional; unset is NOT VERIFIED."""
    return bool(os.environ.get("UCSC_API_KEY"))


_KEY_RE = re.compile(r"(apiKey=)[^&\s]*", re.I)


def _scrub(msg: str) -> str:
    """Strip the API key out of anything browser-bound.

    urllib embeds the request URL in some of its errors, and the key rides in the query
    string. `note` reaches the browser, so scrub explicitly rather than trusting a caller.
    """
    out = _KEY_RE.sub(r"\1***", str(msg))
    key = os.environ.get("UCSC_API_KEY")
    if key:
        out = out.replace(key, "***")       # belt and braces: any other embedding
    return out


def _state(state: str, reason: str, products: Optional[list] = None,
           short: str = "") -> dict:
    # Every return path funnels through here, so scrubbing once covers them all, and
    # `available` cannot drift out of step with `state`: only a real answer from UCSC,
    # good or bad, counts as having asked.
    #
    # `short` is the same finding in one sentence, for the table; `note` is the whole of it,
    # for the page that gets filed. Both are scrubbed: either can carry a URL back from a
    # failure, and a key leaks the same either way.
    return {"state": state, "note": _scrub(reason), "short": _scrub(short or reason),
            "products": products or [], "db": DB, "caveat": CAVEAT,
            "available": state in (ONE_PRODUCT, DANGER)}


# --------------------------------------------------------------------------- #
# Rate gate - process-wide, because the limit belongs to UCSC and not to a caller
# --------------------------------------------------------------------------- #

_gate = threading.Lock()
_last = 0.0             # monotonic: a wall clock that steps backwards must not open the gate
_day_start = 0.0
_day_count = 0


class QuotaExhausted(RuntimeError):
    """Raised by _gate_pass when UCSC's daily budget is spent. Caller degrades to UNKNOWN."""


def _gate_pass() -> None:
    """Block until UCSC may be called again. Serialises every caller in this process.

    Holds the lock across the sleep on purpose: concurrency above one would breach the
    published interval no matter how any caller's pool is sized.
    """
    global _last, _day_start, _day_count
    with _gate:
        now = time.time()
        if now - _day_start >= 86_400:
            _day_start, _day_count = now, 0
        if _day_count >= MAX_PER_DAY:
            raise QuotaExhausted(
                f"UCSC's published daily limit of {MAX_PER_DAY:,} requests is spent for "
                f"today. Verification stops here rather than exceed it."
            )
        wait = MIN_INTERVAL_S - (time.monotonic() - _last)
        if wait > 0:
            pb._emit(pb.Tag.INFO, f"UCSC In-Silico PCR: waiting {wait:.0f}s for the "
                                  f"published rate limit of one request per "
                                  f"{MIN_INTERVAL_S:.0f}s")
            time.sleep(wait)
        # ponytail: counter is process-local and resets on restart, so a redeploy forgives
        # the day's tally. Move to a shared store only if this ever runs multi-worker.
        _last = time.monotonic()
        _day_count += 1


# --------------------------------------------------------------------------- #
# Fetch and parse
# --------------------------------------------------------------------------- #

# pb._http cannot carry this: its _body_is_expected requires JSON or XML for every
# endpoint, and hgPcr answers HTML, so every reply would be dropped as poison and retried.
# Hence a private fetch, with no retry: a retry spends another 15s slot and another of the
# 5000 to re-ask a question that already failed. One try, then UNKNOWN.
def _fetch(fwd: str, rev: str, max_product: int, key: str, timeout: int) -> str:
    """One hgPcr request. Rate-gated. Returns the raw page."""
    url = UCSC + "/cgi-bin/hgPcr?" + urllib.parse.urlencode({
        "db": DB, "wp_target": "genome", "wp_f": fwd, "wp_r": rev,
        "wp_size": max_product, "wp_perfect": PERFECT, "wp_good": PERFECT,
        "boolshad.wp_flipReverse": 0, "apiKey": key,
    })
    _gate_pass()
    # Never the url: it carries the key, and _emit's line is shown to a user.
    pb._emit(pb.Tag.FETCH, f"UCSC In-Silico PCR ({DB})")
    req = urllib.request.Request(url, headers={"User-Agent": "panelbuilder",
                                               "Accept": "text/html"})
    with urllib.request.urlopen(req, timeout=timeout) as f:
        return f.read().decode(errors="replace")


_PRE_RE = re.compile(r"<PRE>(.*?)</PRE>", re.I | re.S)
_TAG_RE = re.compile(r"<[^>]+>")
# hgPcr wraps the position in an anchor INSIDE the FASTA header, so the tags come off
# before the header is read. A naive parser takes '<A HREF=...' for the sequence name.
_HEAD_RE = re.compile(r"^>\s*(chr[\w]+):(\d+)([+-])(\d+)\s+(\d+)bp", re.I)
_NOMATCH_RE = re.compile(r"No matches to", re.I)


def _parse(txt: str) -> Optional[list]:
    """Amplicons from an hgPcr page. None means the shape was not recognised.

    None is not zero amplicons: a CAPTCHA page, an error page and a changed layout all
    contain no products, and reporting them as "no product found" would claim UCSC
    answered a question it never ran.
    """
    products, unreadable = [], False
    for block in _PRE_RE.findall(txt):
        for raw in block.splitlines():
            # Tags off and entities decoded BEFORE deciding this is a header, so the test
            # and the regex read the same text. Ordered the other way, a header whose '>'
            # UCSC had moved inside the anchor, or escaped to &gt;, would not look like a
            # header at all: it would be skipped as quietly as a sequence line, `unreadable`
            # would stay False, and a page listing two loci would come back holding one.
            # That is this module's one forbidden outcome, a multi-locus pair called clean,
            # reached by markup drift alone. hgPcr grew a Turnstile CAPTCHA recently, so
            # drift is this endpoint's observed behaviour rather than a hypothetical.
            line = html.unescape(_TAG_RE.sub("", raw)).strip()
            if not line.startswith(">"):
                continue
            m = _HEAD_RE.match(line)
            if not m:
                unreadable = True       # a header we cannot read is not zero amplicons
                continue
            products.append({"chrom": m.group(1), "start": int(m.group(2)),
                             "strand": m.group(3), "end": int(m.group(4)),
                             "size": int(m.group(5))})
    if unreadable:
        return None
    if products:
        return products
    # Only UCSC's own no-match sentence proves the search ran and found nothing.
    return [] if _NOMATCH_RE.search(txt) else None


def _norm(chrom: str) -> str:
    """'chr11', 'Chr11' and '11' are one chromosome. Compare through here, never raw."""
    return str(chrom).strip().casefold().removeprefix("chr")


# The two suffixes that mark a REDUNDANT representation of sequence already in the primary
# assembly: an alternate haplotype, or a patch correcting it. No genome carries a locus and
# its own alt copy, so a hit on each is UCSC reporting the same place twice.
#
# _random and chrUn_* are deliberately absent. They look alike and are not: they are real
# primary-assembly sequence, present in every genome, and a hit there is a genuine second
# locus. A suffix rule that swept them in would talk a real multi-locus pair down.
_REDUNDANT_SUFFIXES = ("_alt", "_fix")


def _alt_caveat(products: list) -> str:
    """One clause where alt/fix scaffolds are among the hits, otherwise nothing.

    Wording, never state: this reports what UCSC cannot tell us, and stays DANGER.
    A pair whose extra hits are all on alt scaffolds is PROBABLY one locus reported twice,
    and probably is not a verdict. hgPcr's output cannot separate a redundant alt copy from
    a real second locus on that haplotype, so promoting it to a pass would guess toward
    clean, which is this app's worst direction to guess in. The reader gets the scaffold
    names, which are already in the note, plus the one fact that recontextualises them.
    """
    if not any(str(p.get("chrom", "")).endswith(_REDUNDANT_SUFFIXES) for p in products):
        return ""
    return (" Some of these are on alternate-haplotype or fix scaffolds, which may be "
            "redundant representations of the same locus rather than separate loci; UCSC "
            "cannot distinguish those here.")


def _classify(products: Optional[list], chrom: str, expect_bp: int, tol_bp: int) -> dict:
    """Products to a state. Everything that is not exactly one clean product warns."""
    if products is None:
        return _state(UNKNOWN,
                      "UCSC In-Silico PCR returned a page this app could not read, so "
                      "this pair is NOT VERIFIED. It has not been tested for products "
                      "elsewhere in GRCh38 and may amplify more than one locus. An "
                      "unreadable answer is not a clean answer.",
                      short="UCSC sent back a page this app could not read, so the pair "
                            "is still unverified.")

    if not products:
        return _state(DANGER,
                      f"DANGEROUS: UCSC In-Silico PCR found no product for this pair in "
                      f"{DB}. Primer3 designed it against an Ensembl GRCh38 template, so "
                      f"a null result means the two disagree. Do not order without "
                      f"investigating.",
                      short=f"UCSC finds no product at all in {DB}, where this pair was "
                            f"designed to give one. Do not order it.")

    if len(products) > 1:
        where = ", ".join(f"{p['chrom']}:{p['start']} ({p['size']}bp)" for p in products[:5])
        return _state(DANGER,
                      f"DANGEROUS: UCSC In-Silico PCR found {len(products)} products for "
                      f"this pair in {DB}, not one ({where}).{_alt_caveat(products)} A pair "
                      f"that amplifies more than one locus cannot be genotyped: the trace "
                      f"is a mixture of loci, and a heterozygote at the marker is "
                      f"indistinguishable from two loci differing at that base. Do not "
                      f"order this pair without redesigning.", products,
                      short=f"UCSC finds {len(products)} products in {DB}, not one. Check "
                            f"the placements before ordering.")

    p = products[0]
    if _norm(p["chrom"]) != _norm(chrom):
        return _state(DANGER,
                      f"DANGEROUS: UCSC In-Silico PCR placed a product on {p['chrom']}, "
                      f"not the marker's {chrom}. This pair does not amplify the marker. "
                      f"Do not order it.", products,
                      short=f"UCSC puts the product on {p['chrom']}, not the marker's "
                            f"{chrom}: this pair does not amplify the marker.")

    if abs(p["size"] - expect_bp) > tol_bp:
        return _state(DANGER,
                      f"DANGEROUS: UCSC In-Silico PCR reports a product of {p['size']}bp; "
                      f"Primer3 designed for {expect_bp}bp (tolerance {tol_bp}bp). The "
                      f"pair is not binding where it was designed to bind. Confirm before "
                      f"ordering.", products,
                      short=f"UCSC reports {p['size']}bp where {expect_bp}bp was designed "
                            f"for: the pair is not binding where it should.")

    return _state(ONE_PRODUCT,
                  f"UCSC In-Silico PCR: one product, {p['chrom']}, {p['size']}bp, as "
                  f"designed.", products,
                  short=f"UCSC finds one product, {p['chrom']}, {p['size']}bp, as designed.")


def verify(fwd: str, rev: str, chrom: str, expect_bp: int,
           max_product: int = REPORT_MAX_BP, size_tol_bp: int = 10,
           timeout: int = 30) -> dict:
    """Check one primer pair against the whole of hg38. -> {state, note, products, caveat}.

    Takes at least MIN_INTERVAL_S per call by design, so 11 pairs take about 165 seconds.
    Never raises for operational reasons: no key, quota spent, service down and an
    unreadable page all degrade to a NOT VERIFIED state with a readable `note`.

    `max_product` is the question's width, not a filter on its answer: see REPORT_MAX_BP.
    Lower it and this function answers a narrower question than the one its verdict claims.
    """
    if not fwd or not rev:
        raise ValueError("verify() needs both primers")

    key = os.environ.get("UCSC_API_KEY")
    if not key:
        return _state(NOT_CHECKED, NOT_CHECKED_NOTE, short=NOT_CHECKED_SHORT)

    try:
        txt = _fetch(fwd, rev, max_product, key, timeout)
    except QuotaExhausted as e:
        return _state(UNKNOWN, f"{e} This pair is NOT VERIFIED: it has not been tested "
                               f"for products elsewhere in GRCh38.",
                      short="Still unverified: this server has spent its daily UCSC quota.")
    except Exception as e:  # noqa: BLE001 - urllib raises a zoo of types
        return _state(UNKNOWN,
                      f"UCSC In-Silico PCR unavailable ({type(e).__name__}: {e}), so this "
                      f"pair is NOT VERIFIED. It has not been tested for products "
                      f"elsewhere in GRCh38 and may amplify more than one locus.",
                      short=f"Still unverified: UCSC could not be reached "
                            f"({type(e).__name__}).")

    return _classify(_parse(txt), chrom, expect_bp, size_tol_bp)


# --------------------------------------------------------------------------- #
# Self-check - must pass with no key and no network
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    # Recorded 2026-07-16 from a keyless GET to hgPcr, trimmed. This is what an unkeyed
    # request actually gets: it must never parse as a result.
    CAPTCHA = (
        "<html><head>\n<script>\nfunction showWidget() { \nturnstile.render('#myWidget', "
        "{\nsitekey: '0x4AAAAAABgxfUUSzYakb-Pd',\ntheme: 'light',\n});\n}\n</script>\n"
        "</head><body>\n<h1>Verifying you are human</h1>\n</body></html>"
    )
    # SYNTHETIC, not recorded: built from hgPcr.c's <TT><PRE> wrapper and gfPcrLib.c's
    # outputFa() label, because no key was available to observe a real one. The anchor
    # inside the FASTA header is the detail most likely to be wrong. Re-record against a
    # real keyed response before trusting the parser.
    ONE = ("<HTML><BODY><TT><PRE>\n"
           '><A HREF="../cgi-bin/hgTracks?db=hg38&position=chr22">chr22:31000551+31001000'
           "</A> 449bp GAGGGGGTCTCACTGTGTTG TTGACCTTGGCCTTGAACTC\n"
           "acgtacgtacgtacgtacgt\nACGTACGTACGTACGTACGT\n"
           "</PRE></TT></BODY></HTML>")
    TWO = ("<HTML><BODY><TT><PRE>\n"
           '><A HREF="x">chr22:31000551+31001000</A> 449bp AAAA TTTT\nacgt\n'
           '><A HREF="x">chr7:5530601+5531050</A> 449bp AAAA TTTT\nacgt\n'
           "</PRE></TT></BODY></HTML>")
    NOMATCH = ("<HTML><BODY><p>No matches to GAGGGGGTCTCACTGTGTTG TTGACCTTGGCCTTGAACTC "
               "in hg38 UCSC Genes.</p></BODY></HTML>")

    key = os.environ.pop("UCSC_API_KEY", None)      # make the check deterministic
    try:
        # 1. No key -> NOT VERIFIED, and it says so instead of raising or passing.
        assert available() is False
        out = verify("GAGGGGGTCTCACTGTGTTG", "TTGACCTTGGCCTTGAACTC", "chr22", 449)
        assert out["state"] == NOT_CHECKED, out
        assert out["available"] is False
        assert out["state"] != ONE_PRODUCT
        assert "not configured" in out["note"], out["note"]
        assert set(out) >= {"state", "note", "short", "products", "caveat", "available"}

        # 2. The caveat says what a pass is worth, on every return path.
        for word in ["in-silico", "not a", "reference", "dropout", "lab"]:
            assert word.lower() in out["caveat"].lower(), word

        # 3. An unreadable body is UNKNOWN, never a pass. The CAPTCHA is the live case.
        for body in [CAPTCHA, "", "<html><body>502 Bad Gateway</body></html>",
                     '{"error":"nope"}', "<TT><PRE>\n>garbage header\nacgt\n</PRE></TT>"]:
            assert _parse(body) is None, body[:40]
            got = _classify(_parse(body), "chr22", 449, 10)
            assert got["state"] == UNKNOWN, (body[:40], got)
            assert got["state"] != ONE_PRODUCT
            assert got["available"] is False
            assert "NOT VERIFIED" in got["note"], got["note"]

        # 4. Parsing the shapes hgPcr emits.
        assert _parse(NOMATCH) == []                    # searched, found nothing
        one = _parse(ONE)
        assert one == [{"chrom": "chr22", "start": 31000551, "strand": "+",
                        "end": 31001000, "size": 449}], one
        assert len(_parse(TWO)) == 2

        # 4b. A header this parser cannot read must never come back as ONE FEWER PRODUCT.
        # Losing a locus silently is the one failure with no floor under it: the remaining
        # product classifies as a clean single band, so a pair that amplifies two loci is
        # handed over as verified. Every shape below hides the '>' from a naive line test,
        # which is why the tags come off and the entities decode BEFORE that test runs.
        for hidden in [
            '<A HREF="x">&gt;chr7:5530601+5531050</A> 449bp AAAA TTTT',   # '>' inside the tag
            '&gt;chr7:5530601+5531050 449bp AAAA TTTT',                   # '>' escaped
            '<SPAN>></SPAN>chr7:5530601+5531050 449bp AAAA TTTT',         # '>' wrapped alone
        ]:
            page = ("<HTML><BODY><TT><PRE>\n"
                    '><A HREF="x">chr22:31000551+31001000</A> 449bp AAAA TTTT\nacgt\n'
                    f"{hidden}\nacgt\n</PRE></TT></BODY></HTML>")
            got = _parse(page)
            assert got is None or len(got) == 2, (
                f"a page listing 2 loci parsed as {len(got)}: {hidden[:40]}")
            assert _classify(got, "chr22", 449, 10)["state"] != ONE_PRODUCT, \
                f"a two-locus page passed as one clean product: {hidden[:40]}"

        # 5. Classification. Exactly one clean product passes; everything else warns.
        assert _classify(_parse(ONE), "chr22", 449, 10)["state"] == ONE_PRODUCT
        assert _classify(_parse(ONE), "22", 449, 10)["state"] == ONE_PRODUCT   # chr-prefix
        assert _classify(_parse(ONE), "chr22", 445, 10)["state"] == ONE_PRODUCT  # in tol
        assert _classify(_parse(TWO), "chr22", 449, 10)["state"] == DANGER      # 2 bands
        assert _classify(_parse(NOMATCH), "chr22", 449, 10)["state"] == DANGER  # no band
        assert _classify(_parse(ONE), "chr7", 449, 10)["state"] == DANGER       # elsewhere
        assert _classify(_parse(ONE), "chr22", 900, 10)["state"] == DANGER      # wrong size

        # Every danger keeps the primer and says DANGEROUS out loud.
        for bad in [_classify(_parse(TWO), "chr22", 449, 10),
                    _classify(_parse(NOMATCH), "chr22", 449, 10),
                    _classify(_parse(ONE), "chr7", 449, 10),
                    _classify(_parse(ONE), "chr22", 900, 10)]:
            assert bad["state"] == DANGER
            assert "DANGEROUS" in bad["note"], bad["note"]

        # 5c. Alt/fix scaffolds are reported, never reasoned away. The clause appears, the
        #     state stays DANGER, and the pair keeps its sequences. _random and chrUn are
        #     real primary sequence: a hit there is a real second locus and must NOT be
        #     softened, which is the whole risk in a suffix rule.
        def _hits(*chroms):
            return [{"chrom": c, "start": 1000 + i, "strand": "+", "end": 1450, "size": 449}
                    for i, c in enumerate(chroms)]

        alt = _classify(_hits("chr6", "chr6_GL000251v2_alt"), "chr6", 449, 10)
        assert alt["state"] == DANGER, "an alt hit is not reasoned into a pass"
        assert "redundant representations" in alt["note"], alt["note"]
        fix = _classify(_hits("chr1", "chr1_KN196472v1_fix"), "chr1", 449, 10)
        assert "redundant representations" in fix["note"], fix["note"]
        for real in (_hits("chr1", "chr1_KI270711v1_random"),
                     _hits("chr1", "chrUn_GL000195v1"),
                     _hits("chr1", "chr7")):
            got = _classify(real, "chr1", 449, 10)
            assert got["state"] == DANGER, got
            assert "redundant" not in got["note"], \
                f"{real[1]['chrom']} is real sequence in every genome: a hit there is a " \
                f"second locus and must not be talked down"
        # The note claims products, and only the clause speaks about loci being distinct.
        assert "not one" in alt["short"] and "loci" not in alt["short"], alt["short"]

        # Every branch carries both lengths, and they are two lengths. The short one is what
        # the table shows beside the pair, so a branch that forgot it falls back to `note`
        # and silently prints a paragraph per row: green everywhere, unreadable on screen.
        for out in [_classify(_parse(ONE), "chr22", 449, 10),
                    _classify(_parse(TWO), "chr22", 449, 10),
                    _classify(_parse(NOMATCH), "chr22", 449, 10),
                    _classify(_parse(ONE), "chr7", 449, 10),
                    _classify(_parse(ONE), "chr22", 900, 10),
                    _classify(None, "chr22", 449, 10)]:
            assert out["short"], out["state"]
            assert len(out["short"]) <= 140, (out["state"], len(out["short"]))
            assert out["short"] != out["note"], f"{out['state']} has no short form"
        assert len(CAVEAT_SHORT) <= 140 and CAVEAT_SHORT != CAVEAT

        # 5b. The question is at least as wide as UCSC's own, and wider than any product
        #     the design can ask for. Narrowing this is invisible: every page still parses,
        #     every verdict still reads clean, and a second locus past the cap is simply
        #     never on the page to be found. Pinned against primers.DEFAULTS rather than a
        #     number, so raising the design's ceiling without raising this fails here.
        import inspect

        import primers
        _default_max = inspect.signature(verify).parameters["max_product"].default
        assert _default_max == REPORT_MAX_BP >= 4000, _default_max
        assert REPORT_MAX_BP > primers.PRODUCT_CAP, (
            f"the design may ask for a {primers.PRODUCT_CAP}bp product and UCSC is only "
            f"asked about {REPORT_MAX_BP}bp: the on-target product could not be reported, "
            f"and its absence classifies as DANGEROUS")

        # 6. The gate honours UCSC's published interval, and the constant IS the published
        #    one. Asserted apart so neither a fast mechanism nor a slow constant hides.
        assert MIN_INTERVAL_S >= 15.0 and MAX_PER_DAY <= 5000
        _saved, MIN_INTERVAL_S = MIN_INTERVAL_S, 0.05   # keep the check quick, not the limit
        try:
            _gate_pass()
            t0 = time.monotonic()
            _gate_pass()
            _gate_pass()
            spacing = time.monotonic() - t0
            assert spacing >= 2 * 0.05, f"gate let calls through {spacing:.3f}s apart"
        finally:
            MIN_INTERVAL_S = _saved

        # 7. The daily cap stops rather than overrun, and degrades to UNKNOWN not a pass.
        _saved_count, _day_count = _day_count, MAX_PER_DAY
        os.environ["UCSC_API_KEY"] = "selfcheck-not-a-real-key"
        try:
            out = verify("ACGT", "TGCA", "chr22", 449)
            assert out["state"] == UNKNOWN, out
            assert out["available"] is False
            assert "daily limit" in out["note"], out["note"]
        finally:
            _day_count = _saved_count
            os.environ.pop("UCSC_API_KEY", None)

        # 8. Key present, service unreachable -> degrades, never raises, and the key is not
        #    in the browser-bound dict. Pins the guarantee, not urllib's message format.
        SECRET = "selfcheck-not-a-real-key"
        os.environ["UCSC_API_KEY"] = SECRET
        try:
            _saved_ucsc, UCSC = UCSC, "http://127.0.0.1:1"     # nothing listening
            _saved, MIN_INTERVAL_S = MIN_INTERVAL_S, 0.0       # do not wait out the gate
            out = verify("ACGT", "TGCA", "chr22", 449)
            assert out["state"] == UNKNOWN, out
            assert out["available"] is False
            assert "NOT VERIFIED" in out["note"], out["note"]
            import json
            assert SECRET not in json.dumps(out), "KEY LEAKED into browser-bound dict"
        finally:
            UCSC, MIN_INTERVAL_S = _saved_ucsc, _saved
            os.environ.pop("UCSC_API_KEY", None)

        # ...and _scrub redacts independently of urllib.
        os.environ["UCSC_API_KEY"] = SECRET
        try:
            assert _scrub("GET https://x/hgPcr?db=hg38&apiKey=abc123&z=1") == (
                "GET https://x/hgPcr?db=hg38&apiKey=***&z=1")
            assert SECRET not in _scrub(f"boom: key {SECRET} rejected")
        finally:
            os.environ.pop("UCSC_API_KEY", None)

        print("ispcr.py self-check OK (no key, no network)")
    finally:
        if key is not None:
            os.environ["UCSC_API_KEY"] = key

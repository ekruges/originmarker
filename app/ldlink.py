"""Optional LD annotation between two COMMON SNPs. Labelled prior, nothing more (R2).

LD looks like it answers the question the user actually has ("which allele came from
which parent?"), and it does not. Two things it is never allowed to be:

  - a ranking key. Markers are ranked by heterozygosity plus proximity, full stop.
    Ranking by LD prefers markers that correlate in a reference panel over markers
    that are informative in this family.
  - an origin call. r2/D' are population averages over unrelated reference samples; a
    carrier's haplotype is not the population's, and phase needs the lab (R3).

Against a rare pathogenic variant LD is undefined, not merely weak: at AC<=a few the
2x2 haplotype table is nearly empty and any value returned is noise. `allow_rare`
exists to be refused, not enabled: see ld_between().

So this module returns a number with a caveat welded to it, or it returns nothing.
Failures degrade to an unavailable signal rather than an exception, so the UI survives
LDlink being down.

SERVER SIDE ONLY: LDLINK_TOKEN never reaches the browser.

Self-check:  PYTHONPATH=<repo root> python app/ldlink.py
"""

from __future__ import annotations

import json
import os
import re
from typing import Optional

import panelbuilder as pb

LDLINK = "https://ldlink.nih.gov/LDlinkRest"

# R6: LDlink defaults to GRCh37. Everything else in this app is GRCh38, and a
# silently mixed build is exactly the failure R6 exists to prevent. Always explicit.
GENOME_BUILD = "grch38"

CAVEAT = (
    "POPULATION PRIOR, NOT A MEASUREMENT IN THIS FAMILY. r2/D' here are computed over "
    "unrelated 1000 Genomes reference samples and describe that population, not this "
    "carrier. Valid only between two COMMON SNPs: against a rare pathogenic variant "
    "(allele count <= a few) LD is undefined - the haplotype table is empty and any "
    "value returned would be noise. This is never a basis for a parental-origin call, "
    "and it is never used to rank or select markers: ranking is expected "
    "heterozygosity plus proximity only. Phase still requires genotyping the carrier "
    "and phasing against an informative relative; this app cannot determine phase."
)


def available() -> bool:
    """True if an LDlink token is configured. LD is optional; the app works without it."""
    return bool(os.environ.get("LDLINK_TOKEN"))


_TOKEN_RE = re.compile(r"(token=)[^&\s]*", re.I)


def _scrub(msg: str) -> str:
    """Strip the token out of anything browser-bound.

    pb's ApiError embeds the request URL, which carries &token=<secret>, and that message
    lands in `note`. pb's URL truncation is not a guarantee: scrub explicitly.
    """
    out = _TOKEN_RE.sub(r"\1***", str(msg))
    tok = os.environ.get("LDLINK_TOKEN")
    if tok:
        out = out.replace(tok, "***")       # belt and braces: any other embedding
    return out


def _unavailable(reason: str, pop: str) -> dict:
    # Every degraded path funnels through here, so scrubbing once covers them all.
    return {"r2": None, "dprime": None, "pop": pop, "note": _scrub(reason),
            "caveat": CAVEAT, "available": False}


def _num(v) -> Optional[float]:
    """LDlink returns stats as strings, sometimes as 'NA'."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _parse(txt: str) -> tuple[Optional[float], Optional[float]]:
    """(r2, dprime) from an LDpair response. Missing stats come back as None.

    Handles both shapes the endpoint serves, JSON and plain text, and raises pb.ApiError
    on its 'error' payload.
    """
    try:
        data = json.loads(txt)
    except ValueError:
        r2 = re.search(r"R2:\s*([\d.]+)", txt)
        dp = re.search(r"D'\s*:\s*([\d.]+)", txt)
        return (_num(r2.group(1)) if r2 else None, _num(dp.group(1)) if dp else None)

    if not isinstance(data, dict):
        return None, None
    if data.get("error"):
        raise pb.ApiError(str(data["error"]))
    stats = data.get("statistics") or data
    return (_num(stats.get("r2")),
            _num(stats.get("d_prime", stats.get("dprime"))))


def ld_between(a: str, b: str, pop: str = "CEU", allow_rare: bool = False) -> dict:
    """LD between two COMMON SNPs, as a labelled prior. -> {r2, dprime, pop, note, caveat}.

    `allow_rare` is a tripwire, not a feature flag: there is no correct value other than
    False. A caller who knows a variant is rare (pb.assess_rarity ->
    population_LD_usable=False) says so and gets a ValueError instead of a
    plausible-looking number. Refusing is the feature (R2/R3).

    Never raises for operational reasons: no token, service down, timeout and bad input
    all degrade to {..., available: False} with a readable `note`. The R2 guard above is
    the only raise.
    """
    if allow_rare:
        raise ValueError(
            "refusing LD against a variant flagged rare: LD is undefined at AC<=a few "
            "and is never a basis for an origin call or for ranking (R2). Rank markers "
            "by heterozygosity and proximity; establish phase in the lab (R3)."
        )

    token = os.environ.get("LDLINK_TOKEN")
    if not token:
        return _unavailable(
            "LD annotation is disabled: no LDLINK_TOKEN configured on the server. "
            "This is optional - it never affects marker ranking or selection (R2).",
            pop,
        )

    try:
        txt = pb._get(LDLINK, "/ldpair", {
            "var1": a, "var2": b, "pop": pop,
            "genome_build": GENOME_BUILD, "token": token,
        }, tries=2, timeout=10)
        r2, dprime = _parse(txt)
    except Exception as e:  # noqa: BLE001 - every failure here means "no LD today"
        return _unavailable(f"LDlink unavailable: {e}", pop)

    if r2 is None and dprime is None:
        return _unavailable(f"LDlink returned no LD statistics for {a}/{b} in {pop}.", pop)

    return {
        "r2": r2,
        "dprime": dprime,
        "pop": pop,
        "note": (f"LDlink LDpair: {a} vs {b}, 1000 Genomes {pop}, {GENOME_BUILD.upper()} "
                 f"(build stated explicitly, R6). Annotation only - not used in ranking."),
        "caveat": CAVEAT,
        "available": True,
    }


# --------------------------------------------------------------------------- #
# Self-check - must pass with no token and no network
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    token = os.environ.pop("LDLINK_TOKEN", None)   # make the check deterministic
    try:
        # 1. No token -> disabled, and it says so instead of raising.
        assert available() is False
        out = ld_between("rs757110", "rs151344623")
        assert out["available"] is False
        assert out["r2"] is None and out["dprime"] is None
        assert "LDLINK_TOKEN" in out["note"], out["note"]
        assert out["pop"] == "CEU"
        assert set(out) >= {"r2", "dprime", "pop", "note", "caveat"}

        # 2. The caveat carries the whole R2/R3 story, on every single return path.
        for key in ["prior", "COMMON", "undefined", "origin", "rank", "phase"]:
            assert key.lower() in out["caveat"].lower(), key

        # 3. The R2 guard is the one thing that raises.
        try:
            ld_between("rs757110", "rs151344623", allow_rare=True)
        except ValueError as e:
            assert "undefined" in str(e)
        else:
            raise AssertionError("allow_rare=True must be refused (R2)")

        # 4. Token present but service unreachable -> still degrades, never raises, and
        #    the token is not in the browser-bound dict. The assert pins the guarantee,
        #    not pb's URL truncation, so it holds however pb formats the message.
        SECRET = "selfcheck-not-a-real-token"
        os.environ["LDLINK_TOKEN"] = SECRET
        try:
            saved, LDLINK = LDLINK, "http://127.0.0.1:1"   # nothing listening
            out = ld_between("rs151344623", "rs757110")
            assert out["available"] is False, out
            assert "unavailable" in out["note"].lower(), out["note"]
            blob = json.dumps(out)
            assert SECRET not in blob, f"TOKEN LEAKED into browser-bound dict: {blob}"
        finally:
            LDLINK = saved
            os.environ.pop("LDLINK_TOKEN", None)

        # ...and _scrub itself redacts, independent of pb.
        os.environ["LDLINK_TOKEN"] = SECRET
        try:
            assert _scrub("GET https://x/ldpair?var1=rs1&token=abc123&z=1") == (
                "GET https://x/ldpair?var1=rs1&token=***&z=1")
            assert SECRET not in _scrub(f"boom: Authorization {SECRET} rejected")
        finally:
            os.environ.pop("LDLINK_TOKEN", None)

        # 5. Response parsing: JSON, legacy text, and the error shape.
        assert _parse('{"statistics":{"r2":"0.83","d_prime":"0.95"}}') == (0.83, 0.95)
        assert _parse('{"r2":"0.5","dprime":"0.7"}') == (0.5, 0.7)
        assert _parse("R2: 0.42\nD': 0.88\n") == (0.42, 0.88)
        assert _parse('{"statistics":{"r2":"NA","d_prime":"NA"}}') == (None, None)
        try:
            _parse('{"error":"Invalid or expired API token"}')
        except pb.ApiError:
            pass
        else:
            raise AssertionError("LDlink error payload should surface as ApiError")

        print("ldlink.py self-check OK (no token, no network)")
    finally:
        if token is not None:
            os.environ["LDLINK_TOKEN"] = token

"""The terms make factual claims about where data goes. Code moves; the page does not.

Every claim below is a promise to a reader who cannot see the source, about a tool used for
PGT-M work. The failure mode is silent and one-directional: someone adds an outbound call,
every test passes because the call works, and the terms keep saying it does not happen. So
the recipients are asserted against the code that reaches them, not against a list.

This does not read the prose. It answers one question: does the app talk to a host the terms
do not name?
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TERMS = ROOT / "web" / "src" / "TermsPage.tsx"

# Hosts the app reaches, and the section 5 wording that must account for each. A new one is a
# new disclosure, which is the point: this list is short because the app's reach is small.
EXPECTED = {
    "api.anthropic.com": "anthropic",       # free text with no identifier, sent verbatim
    "genome.ucsc.edu": "ucsc",              # the two primer sequences, on request
    "eutils.ncbi.nlm.nih.gov": "clinvar",   # the variant identifier
    "rest.ensembl.org": "ensembl",
    "gnomad.broadinstitute.org": "gnomad",
}

# Where a URL constant lives. app/nl.py reaches Anthropic through their SDK rather than a
# literal, so it is matched on the import instead.
_SOURCES = [ROOT / "panelbuilder.py", *(ROOT / "app").glob("*.py")]

# Loopback is this process talking to itself (the container healthcheck), never a third party.
_LOOPBACK = re.compile(r"^(127\.|localhost$|0\.0\.0\.0$|\[?::1)")


def _hosts_the_app_reaches() -> set:
    """Hosts in the SHIPPING code. Self-check fixtures are excluded on purpose.

    A recorded page under `if __name__ == "__main__"` carries whatever hrefs the real
    response carried, and those are evidence of nothing: the block does not run in the app.
    Counting them would make this fail on a fixture, which teaches the next person to widen
    the ignore list, which is how a check like this dies.
    """
    found = set()
    for path in _SOURCES:
        src = path.read_text()
        src = src.split('if __name__ == "__main__":')[0]
        for m in re.finditer(r"https?://([a-z0-9.\-]+)", src):
            if not _LOOPBACK.match(m.group(1)):
                found.add(m.group(1))
        # The SDK carries the host, so the import is the evidence of the egress.
        if re.search(r"^\s*import anthropic", src, re.M):
            found.add("api.anthropic.com")
    return found


def _terms_prose() -> str:
    """The terms with every href stripped: what a reader actually reads.

    Searching the raw file instead would let a URL stand in for a disclosure. It did: gutting
    the word Anthropic from the prose left `https://www.anthropic.com` in an href and this
    check went green over terms that no longer named the recipient.
    """
    return re.sub(r'href="[^"]*"', "", TERMS.read_text()).lower()


def test_the_terms_name_every_host_the_app_sends_to() -> None:
    terms = _terms_prose()
    reached = _hosts_the_app_reaches()
    unnamed = []
    for host, must_appear in EXPECTED.items():
        if host in reached and must_appear not in terms:
            unnamed.append(f"{host} (terms never say {must_appear!r})")
    assert not unnamed, (
        "the app sends to hosts the terms do not disclose: " + "; ".join(unnamed)
    )


def test_a_new_outbound_host_forces_a_terms_decision() -> None:
    """An egress this file has never heard of. Add it to EXPECTED and to the terms, or, if it
    is not an egress, add it to the ignore list with a reason.
    """
    # Documentation links, licence URLs and provider homepages that are printed rather than
    # called. Each is a place the app POINTS at, never one it sends to.
    NOT_EGRESS = {
        "www.ncbi.nlm.nih.gov",         # dbSNP/ClinVar links printed into exports
        "www.apache.org", "creativecommons.org", "doi.org", "www.ebi.ac.uk",
        "ldlink.nih.gov", "ldlink.nci.nih.gov",   # optional, and already in section 5's scope
        "www.internationalgenome.org", "github.com", "ezrakruger.cc",
        "hgdownload.soe.ucsc.edu",      # the bundled map's provenance, not fetched at runtime
        "schema.org", "www.w3.org",
    }
    unknown = _hosts_the_app_reaches() - set(EXPECTED) - NOT_EGRESS
    assert not unknown, (
        f"new hosts in the code: {sorted(unknown)}. If the app SENDS to one, name it in "
        f"TermsPage section 5 and add it to EXPECTED. If it is only a printed link, add it "
        f"to NOT_EGRESS."
    )


def test_the_terms_still_say_the_load_bearing_things() -> None:
    """The claims a reader acts on, and that a rewrite could quietly drop."""
    t = _terms_prose()
    for claim in [
        "not a clinical diagnostic",     # the whole point
        "candidate",                     # markers and primers alike
        "wet-lab",                       # what the in-silico check is not
        "verbatim",                      # what happens to free text
        "gplv2",                         # the one non-Apache dependency
        "apache license 2.0",
    ]:
        assert claim in t, f"the terms no longer say {claim!r}"

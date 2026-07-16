"""The monogram has one geometry, and three consumers that each hold their own copy of it.

web/public/favicon.svg is the tab icon AND the source app/exports.py renders into the PDF
and XLSX masthead at runtime. web/src/Mark.tsx is the corner mark on the site, and carries
a hardcoded copy: it cannot read the svg, because its colours are CSS custom properties so
the mark tracks the page theme, while the favicon must stand alone with literal hex.

So the colours differ on purpose and the GEOMETRY must not. Without this, changing the
favicon turned the icon, the PDF and the spreadsheet red while the site header stayed blue,
and every gate passed.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FAVICON = ROOT / "web" / "public" / "favicon.svg"
MARK_TSX = ROOT / "web" / "src" / "Mark.tsx"


def _grab(text: str, pattern: str, what: str) -> str:
    m = re.search(pattern, text)
    assert m, f"could not read {what} out of the source; the check itself is stale"
    return m.group(1).strip()


def _num(s: str) -> float:
    return float(s)


@pytest.fixture(scope="module")
def sources() -> tuple[str, str]:
    for p in (FAVICON, MARK_TSX):
        if not p.exists():
            pytest.skip(f"{p.relative_to(ROOT)} not present")
    return FAVICON.read_text(), MARK_TSX.read_text()


def test_ring_matches(sources) -> None:
    svg, tsx = sources
    assert _num(_grab(svg, r'<circle[^>]*\br="([\d.]+)"', "favicon ring radius")) == \
           _num(_grab(tsx, r'<circle[^>]*\br="([\d.]+)"', "Mark.tsx ring radius"))
    assert _num(_grab(svg, r'stroke-width="([\d.]+)"', "favicon stroke width")) == \
           _num(_grab(tsx, r'strokeWidth="([\d.]+)"', "Mark.tsx stroke width"))
    for axis in ("cx", "cy"):
        assert _num(_grab(svg, rf'<circle[^>]*\b{axis}="([\d.]+)"', f"favicon {axis}")) == \
               _num(_grab(tsx, rf'<circle[^>]*\b{axis}="([\d.]+)"', f"Mark.tsx {axis}"))


def test_m_glyph_matches(sources) -> None:
    """The M is an outlined Merriweather glyph: 246 characters nobody can proofread, which
    is exactly why it needs a machine to compare rather than a promise in a comment."""
    svg, tsx = sources
    assert _grab(svg, r'<path d="([^"]+)"', "favicon M path") == \
           _grab(tsx, r"const D = '([^']+)'", "Mark.tsx M path")
    assert _grab(svg, r'<path[^>]*transform="([^"]+)"', "favicon M transform") == \
           _grab(tsx, r"const TRANSFORM = '([^']+)'", "Mark.tsx M transform")


def test_viewboxes_match(sources) -> None:
    svg, tsx = sources
    assert _grab(svg, r'viewBox="([^"]+)"', "favicon viewBox") == \
           _grab(tsx, r'viewBox="([^"]+)"', "Mark.tsx viewBox")


def test_colours_are_allowed_to_differ(sources) -> None:
    """Pins the one difference that is intentional, so a future reader does not "fix" it.

    The favicon is a standalone file with no stylesheet, so it needs literal hex. The
    component is themeable and reads the same two colours from CSS variables.
    """
    svg, tsx = sources
    assert "#2e6da4" in svg and "#337ab7" in svg
    assert "var(--om-blue)" in tsx and "var(--om-blue-light)" in tsx

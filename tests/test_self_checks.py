"""Run every module's __main__ self-check under pytest.

Subprocess rather than import: the checks live under `if __name__ == "__main__"`, which an
import does not execute. Running each module as a module also proves it is executable
standalone, which the deploy path relies on.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

# Every module carrying an assert-based self-check. Add new ones here.
MODULES = [
    "build_info",
    "genetic_map",
    "panelbuilder",
    # Skips itself, exit 0, where primer3 is not installed: it is optional (GPLv2), so the
    # default image runs this check in its skipped form.
    "primers",
    "app.exports",
    "app.ispcr",
    "app.nl",
    "app.ldlink",
]


@pytest.mark.parametrize("mod", MODULES)
def test_module_self_check(mod: str) -> None:
    env = {
        **os.environ,
        # Without both of these the checks reach the network. TTL 0 is falsy and disables
        # expiry (panelbuilder reads `if CACHE_TTL_S and ...`); unpinned, the fixtures go
        # stale and the offline suite silently goes live.
        "PANELBUILDER_CACHE": str(ROOT / "tests" / "fixtures"),
        "PANELBUILDER_CACHE_TTL": "0",
    }
    p = subprocess.run(
        [sys.executable, "-m", mod],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=300,
    )
    assert p.returncode == 0, (
        f"{mod} self-check failed (exit {p.returncode})\n"
        f"--- stdout ---\n{p.stdout}\n--- stderr ---\n{p.stderr}"
    )


def test_every_self_check_is_registered() -> None:
    """A module that grows a self-check must be added to MODULES above, or nothing
    runs it.
    """
    found = set()
    for path in [*ROOT.glob("*.py"), *(ROOT / "app").glob("*.py")]:
        if path.name.startswith("test_"):
            continue
        src = path.read_text()
        if '__name__ == "__main__"' not in src or "assert" not in src:
            continue
        rel = path.relative_to(ROOT)
        found.add(str(rel.with_suffix("")).replace(os.sep, "."))
    missing = found - set(MODULES)
    assert not missing, f"modules with an unregistered self-check: {sorted(missing)}"

"""A key the app reads but compose never forwards is a feature that is silently off.

Compose reads .env for substitution and passes NOTHING into the container it was not asked
to. So a key can be correctly generated, correctly stored, and still never arrive: the app
sees an unset var, degrades gracefully by design, and reports the feature as unavailable.
Every layer behaves exactly as written and the feature does not work. There is no error to
find, which is what makes it worth a test.

This is not hypothetical. UCSC_API_KEY was added to the app, to .env on the NAS, and to the
deploy docs, and was missing from the compose environment block: `insilico_pcr_enabled`
would have stayed false with the key sitting one file away.

Only KEY-SHAPED vars are covered. The rest of what the app reads (ROOT_PATH, PORT, the
rate-limit knobs) carry working defaults, so an unforwarded one changes tuning, not
whether a feature exists.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COMPOSE = ROOT / "docker-compose.yml"

# What the app actually reads, found rather than listed: a hand-kept list is the thing that
# goes stale here.
_READS = re.compile(r"""environ(?:\.get\(|\.pop\(|\[)["']([A-Z][A-Z0-9_]*)["']""")

# A secret, by name. These gate a feature: unset means the feature is off.
_KEY_SHAPED = re.compile(r"(_API_KEY|_TOKEN|_SECRET|_PASSWORD)$")


def _vars_read() -> set:
    found = set()
    for path in [*ROOT.glob("*.py"), *(ROOT / "app").glob("*.py")]:
        found |= set(_READS.findall(path.read_text()))
    return found


def _vars_forwarded() -> set:
    # The environment: block's keys, i.e. what the container is handed.
    body = COMPOSE.read_text()
    block = re.search(r"^    environment:\n((?:^      .*\n)+)", body, re.M)
    assert block, "docker-compose.yml has no environment block, or its indentation changed"
    return set(re.findall(r"^      ([A-Z][A-Z0-9_]*):", block.group(1), re.M))


def test_every_key_the_app_reads_is_forwarded_by_compose() -> None:
    keys = {v for v in _vars_read() if _KEY_SHAPED.search(v)}
    assert keys, "found no key-shaped env vars: the regex has drifted from the code"
    missing = keys - _vars_forwarded()
    assert not missing, (
        f"the app reads {sorted(missing)} but docker-compose.yml does not pass "
        f"{'it' if len(missing) == 1 else 'them'} into the container. A key in .env would "
        f"be read for substitution and then dropped, and the feature would report itself "
        f"unavailable with the key correctly configured."
    )


def test_ucsc_key_is_documented_where_it_is_obtained() -> None:
    """The key is useless without the one page that issues it, which is not guessable."""
    doc = (ROOT / "deploy" / "README-deploy.md").read_text()
    assert "UCSC_API_KEY" in doc
    # The Hub Development page is where UCSC actually issues it. Nothing else does.
    assert "hgHubConnect" in doc, "the deploy doc must name the page that issues the key"


def test_every_root_module_the_app_imports_is_copied_into_the_image() -> None:
    """A root-level module the app imports but the Dockerfile does not COPY.

    The image is built from an explicit COPY list, not from the working tree, so a new
    module at the repo root is invisible to it. Nothing local notices: the tests import from
    the tree and pass, the build succeeds, and the container dies at startup on ImportError
    because `app/` is copied wholesale while its root-level imports are named one by one.
    primers.py shipped exactly this way.
    """
    dockerfile = (ROOT / "Dockerfile").read_text()
    copied = set()
    for line in dockerfile.splitlines():
        if line.startswith("COPY ") and "--from=" not in line:
            # "COPY a.py b.py ./" -> the names between COPY and the destination
            copied |= {p for p in line.split()[1:-1] if p.endswith(".py")}

    # What app/ actually imports from the repo root, found rather than listed.
    roots = {p.stem for p in ROOT.glob("*.py") if not p.name.startswith("test_")}
    needed = set()
    for path in (ROOT / "app").glob("*.py"):
        src = path.read_text()
        for mod in roots:
            if re.search(rf"^\s*(?:import {mod}\b|from {mod} import)", src, re.M):
                needed.add(f"{mod}.py")

    missing = needed - copied
    assert not missing, (
        f"app/ imports {sorted(missing)} at the repo root, and the Dockerfile never COPYs "
        f"{'it' if len(missing) == 1 else 'them'}. The image builds and the container dies "
        f"on ImportError at startup."
    )

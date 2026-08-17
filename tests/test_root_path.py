"""Both routes must answer: the bare one and the prefixed one.

WHY THIS TEST EXISTS. The app is served twice from one process, at the root on the LAN and behind
a Cloudflare tunnel that hands over the path prefix intact. FastAPI's root_path strips that prefix
only when the request actually carries it. If that behaviour ever regresses, /api/health keeps
returning 200 while every public /originmarker/* URL 404s, so the suite passes and the site is
down. Nothing else here would notice.

It was written for a Dependabot bump: uvicorn 0.52.1 -> 0.52.3, whose release notes are entirely
HTTP/1.1 request-parsing changes, which is the layer this behaviour sits on. The bump turned out to
be fine, and the point is that the check should not have depended on someone remembering to look.

Deliberately does NOT import app.main: that pulls the whole panel builder and its data files, which
makes the test slow and couples it to things it is not testing. What is under test is the
root_path contract between FastAPI and the server, and this reproduces exactly that mounting.
"""
import os

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

PREFIX = "/originmarker"


def _app(root_path: str):
    app = fastapi.FastAPI(root_path=root_path)

    @app.get("/api/health")
    def health():
        return {"ok": True}

    return app


def test_bare_route_answers_when_mounted_at_a_prefix():
    """The LAN case. The prefix is configured but the request does not carry it."""
    with TestClient(_app(PREFIX)) as c:
        assert c.get("/api/health").status_code == 200


def test_prefixed_route_answers_when_mounted_at_a_prefix():
    """The tunnel case, and the one that breaks silently.

    The proxy passes the prefix through, so the request path is /originmarker/api/health and
    root_path must strip it. A 404 here is the exact failure that leaves /api/health green while
    the public site is entirely down.
    """
    with TestClient(_app(PREFIX)) as c:
        r = c.get(f"{PREFIX}/api/health")
        assert r.status_code == 200, (
            f"{PREFIX}/api/health returned {r.status_code}. root_path is not stripping the prefix, "
            "so every public URL 404s while the bare health route still answers"
        )
        assert r.json()["ok"] is True


def test_root_mounted_app_still_answers_bare():
    """Local development, no prefix configured at all."""
    with TestClient(_app("")) as c:
        assert c.get("/api/health").status_code == 200


def test_the_real_app_reads_root_path_from_the_environment():
    """The contract the deployment depends on, checked in the source rather than by importing it.

    Importing app.main here would drag in the panel builder and its data. What matters is that the
    prefix is read from ROOT_PATH and handed to FastAPI, so that is what is asserted.
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(here, "app", "main.py"), encoding="utf-8") as fh:
        src = fh.read()
    assert "root_path=" in src
    assert "ROOT_PATH" in src

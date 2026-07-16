"""API self-check against the golden ABCC8 case, offline.

PANELBUILDER_CACHE is pointed at fixtures/ before app.main imports panelbuilder, so
/api/resolve and /api/panel replay recorded API responses and never touch the network.

    PANELBUILDER_CACHE=tests/fixtures .venv/bin/python -m pytest tests/test_api.py
"""
import json
import os
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
os.environ.setdefault("PANELBUILDER_CACHE", str(FIXTURES))
os.environ.setdefault("PANELBUILDER_CACHE_TTL", "0")   # fixtures are frozen, not stale

sys.path.insert(0, str(ROOT))

from fastapi.testclient import TestClient  # noqa: E402

from app import jobs  # noqa: E402
from app.main import app  # noqa: E402

GOLDEN = "NM_000352.6(ABCC8):c.3989-9G>A"
client = TestClient(app)


@pytest.fixture(scope="module", autouse=True)
def _need_fixtures():
    if not FIXTURES.exists():
        pytest.skip("tests/fixtures not present")


def test_health():
    j = client.get("/api/health").json()
    assert j["ok"] is True
    assert j["gnomad_dataset"] == "gnomad_r4"
    assert j["build"] == "GRCh38"
    assert set(j) >= {"version", "ensembl_release", "map_source",
                      "ldlink_enabled", "nl_enabled"}


def test_resolve_golden():
    r = client.post("/api/resolve", json={"variant": GOLDEN})
    assert r.status_code == 200
    j = r.json()
    assert j["variant"]["rsid"] == "rs151344623"
    assert j["variant"]["pos_grch38"] == 17_397_055
    assert j["variant"]["build"] == "GRCh38"
    assert j["transcript_sense"].startswith("G>A")      # minus strand, genomic C>T
    assert j["rarity"]["population_LD_usable"] is False
    assert "phasing" in j["ld_banner"].lower()
    # One verdict, not two: ld_banner must be rarity.reason verbatim, not a second wording
    # that can drift out of agreement with it. Only string identity pins that.
    assert j["ld_banner"] == j["rarity"]["reason"]
    assert "VCV000009088" in j["clinvar_url"]


def test_bad_variant_fails_loudly(monkeypatch):
    """Unresolvable input must 400, never fall back to a guessed coordinate."""
    # Monkeypatched because pb-raises -> 400 is this layer's whole contract, and the real
    # call costs ~90s of Ensembl retries. test_bad_variant_live does it for real.
    import panelbuilder as pb

    def boom(*a, **kw):
        raise pb.ApiError("cannot resolve rsNOPE on GRCh38")

    monkeypatch.setattr("panelbuilder.resolve_variant", boom)
    r = client.post("/api/resolve", json={"variant": "rsNOPE"})
    assert r.status_code == 400
    d = r.json()["detail"].lower()
    assert "resolve" in d
    assert "r1" not in d, "internal rule IDs must not appear in user-facing errors"


@pytest.mark.skipif(os.environ.get("RUN_LIVE") != "1", reason="live APIs; set RUN_LIVE=1")
def test_bad_variant_live():
    # ~90s: pb retries Ensembl's 404 with backoff. Hence the gate.
    r = client.post("/api/resolve", json={"variant": "rs00000000000000000"})
    assert r.status_code == 400


def test_panel_job_reaches_done():
    r = client.post("/api/panel", json={"variant": GOLDEN})
    assert r.status_code == 202
    job_id = r.json()["job_id"]

    for _ in range(600):                                 # cached => seconds, not minutes
        j = client.get(f"/api/panel/{job_id}").json()
        if j["status"] != "running":
            break
        time.sleep(0.1)
    assert j["status"] == "done", j.get("error")

    res = j["result"]
    assert len(res["candidates"]) == 1202
    assert len(res["recommended"]) == 20
    nearest = min(res["candidates"], key=lambda m: abs(m["dist"]))
    assert nearest["rsid"] == "rs757110" and nearest["dist"] == -125
    assert res["coverage"]["lower_count"] and res["coverage"]["higher_count"]
    assert res["provenance"]["disclaimer"] == pb_disclaimer()
    assert res["provenance"]["build"] == "GRCh38"


def test_sse_stream_ends_with_done():
    job_id = client.post("/api/panel", json={"variant": GOLDEN}).json()["job_id"]
    events, stages = [], []
    with client.stream("GET", f"/api/panel/{job_id}/stream") as r:
        assert r.status_code == 200
        assert r.headers["content-type"].startswith("text/event-stream")
        for line in r.iter_lines():
            if line.startswith("event:"):
                events.append(line.partition(":")[2].strip())
            elif line.startswith("data:") and events[-1:] == ["progress"]:
                stages.append(json.loads(line.partition(":")[2]))
    assert events[0] == "progress" and events[-1] == "done"
    assert stages[-1]["fraction"] == 1.0
    assert any("pulling variants" in s["stage"] for s in stages)


def test_panel_rate_limit_429(monkeypatch):
    monkeypatch.setattr("app.main.PER_IP_MAX", 0)
    r = client.post("/api/panel", json={"variant": GOLDEN})
    assert r.status_code == 429
    assert "rate limit" in r.json()["detail"].lower()


def test_invalid_query_rejected_before_any_network():
    # StructuredQuery.__post_init__ owns validation; the API must surface it as 400.
    r = client.post("/api/panel", json={"variant": GOLDEN, "window_bp": 99})
    assert r.status_code == 400
    r = client.post("/api/panel", json={"variant": GOLDEN, "ancestry": "KLINGON"})
    assert r.status_code == 400


def test_exports_round_trip():
    # exports.py consumes the PanelResult dataclass, so this also pins jobs.py to storing
    # the object rather than to_dict().
    job_id = client.post("/api/panel", json={"variant": GOLDEN}).json()["job_id"]
    for _ in range(600):
        if client.get(f"/api/panel/{job_id}").json()["status"] != "running":
            break
        time.sleep(0.1)

    for ext, magic in [("csv", b"rs757110"), ("json", b"rs757110"),
                       ("xlsx", b"PK"), ("pdf", b"%PDF")]:
        r = client.get(f"/api/export/{job_id}.{ext}")
        assert r.status_code == 200, f"{ext}: {r.text[:200]}"
        assert magic in r.content, f"{ext} body looks wrong"
        cd = r.headers["content-disposition"]
        assert cd.startswith("attachment;") and f".{ext}" in cd
        assert "GRCh38" in cd
    # The disclaimer rides every export, verbatim.
    body = client.get(f"/api/export/{job_id}.csv").content.decode()
    assert pb_disclaimer() in body

    assert client.get(f"/api/export/{job_id}.docx").status_code == 404


def test_ld_503_without_token(monkeypatch):
    # LD is optional annotation. No token -> 503, never a fabricated number.
    monkeypatch.delenv("LDLINK_TOKEN", raising=False)
    r = client.get("/api/ld", params={"a": "rs757110", "b": "rs2237982", "pop": "CEU"})
    assert r.status_code == 503
    assert "LDLINK_TOKEN" in r.json()["detail"]


def test_nl_parses_intent_without_coordinates():
    r = client.post("/api/nl", json={"text": f"build a panel for {GOLDEN} in EAS"})
    assert r.status_code == 200
    j = r.json()
    assert j["query"]["variant"] == GOLDEN
    assert j["query"]["ancestry"] == "EAS"
    # StructuredQuery has no coordinate fields, so the NL layer cannot emit one.
    assert not {"chrom", "pos", "pos_grch38", "ref", "alt"} & set(j["query"])
    assert client.post("/api/nl", json={"text": "hello"}).status_code == 400


def test_nl_refuses_to_choose_between_two_identifiers():
    """Two identifiers in one sentence is ambiguous: refuse, do not pick one."""
    r = client.post("/api/nl", json={"text": "not rs1801133, I mean rs151344623"})
    assert r.status_code == 400
    d = r.json()["detail"]
    assert "rs1801133" in d and "rs151344623" in d, "say which identifiers clash"


def test_nl_global_cap_is_checked_before_the_model_is_called(monkeypatch):
    """The cap bounds spend, so it must gate the call, not merely the response."""
    import app.nl

    billed = []
    monkeypatch.setattr(app.nl, "_llm_intent",
                        lambda text: billed.append(text) or {"variant": "rs151344623"})
    monkeypatch.setattr("app.main.GLOBAL_NL_MAX", 0)

    r = client.post("/api/nl", json={"text": "markers near the ABCC8 splice mutation"})
    assert r.status_code == 429
    assert billed == [], "capped requests must not reach the model at all"


def test_nl_free_path_is_never_metered(monkeypatch):
    """Text carrying an rsID is answered by regex, so it costs nothing and must not be
    charged against a budget that exists to bound spend, even with the wallet shut."""
    monkeypatch.setattr("app.main.GLOBAL_NL_MAX", 0)
    monkeypatch.setattr("app.main.PER_IP_NL", 0)
    r = client.post("/api/nl", json={"text": "panel for rs151344623 please"})
    assert r.status_code == 200
    assert r.json()["used_llm"] is False


def test_job_errors_do_not_leak_internals(monkeypatch):
    """job.error is served to the browser, so an unforeseen exception must not arrive
    there as a raw traceback string."""
    import panelbuilder as pb

    def boom(*a, **kw):
        raise OSError("[Errno 8] nodename nor servname provided")

    monkeypatch.setattr(pb, "build", boom)
    job_id = client.post("/api/panel", json={"variant": GOLDEN}).json()["job_id"]
    for _ in range(600):
        j = client.get(f"/api/panel/{job_id}").json()
        if j["status"] != "running":
            break
        time.sleep(0.1)
    assert j["status"] == "error"
    assert "OSError" not in j["error"] and "Errno" not in j["error"], j["error"]


def test_sse_stream_replays_to_a_client_that_subscribed_late():
    """A job can finish before the browser opens the stream. A late subscriber must get
    the buffered log and its terminal event, not hang on a 'done' it missed."""
    job_id = client.post("/api/panel", json={"variant": GOLDEN}).json()["job_id"]
    for _ in range(600):
        if jobs.get(job_id).status != "running":
            break
        time.sleep(0.1)
    assert jobs.get(job_id).status == "done"

    events = []
    with client.stream("GET", f"/api/panel/{job_id}/stream") as r:
        for line in r.iter_lines():
            if line.startswith("event:"):
                events.append(line.partition(":")[2].strip())
    assert events[-1] == "done", events[-3:]
    assert events.count("progress") > 1, "the buffered log should replay, not just status"


def test_finished_jobs_expire_and_the_registry_is_capped():
    """TTL and MAX_JOBS stop the registry pinning every PanelResult it ever built."""
    # The registry is poked directly: submitting MAX_JOBS+ real builds would take hours.
    saved = dict(jobs._jobs)
    try:
        now = time.time()
        jobs._jobs.clear()
        jobs._jobs["old"] = jobs.Job(id="old", status="done", created=now - jobs.TTL_S - 1)
        jobs._jobs["new"] = jobs.Job(id="new", status="done", created=now)
        jobs._jobs["busy"] = jobs.Job(id="busy", status="running",
                                      created=now - jobs.TTL_S - 1)
        jobs._evict_locked()
        assert set(jobs._jobs) == {"new", "busy"}, "TTL must drop finished jobs only"

        jobs._jobs.clear()
        for i in range(jobs.MAX_JOBS + 25):          # under TTL: the cap is what bites
            jobs._jobs[str(i)] = jobs.Job(id=str(i), status="done", created=now + i)
        jobs._evict_locked()
        assert len(jobs._jobs) == jobs.MAX_JOBS
        assert "0" not in jobs._jobs, "the cap must shed the oldest first"
        assert str(jobs.MAX_JOBS + 24) in jobs._jobs
    finally:
        jobs._jobs.clear()
        jobs._jobs.update(saved)


def test_unknown_job_404():
    assert client.get("/api/panel/deadbeef").status_code == 404
    assert client.get("/api/export/deadbeef.csv").status_code == 404


def pb_disclaimer():
    import panelbuilder as pb
    return pb.DISCLAIMER

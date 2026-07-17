"""API self-check against the golden ABCC8 case, offline.

PANELBUILDER_CACHE is pointed at fixtures/ before app.main imports panelbuilder, so
/api/resolve and /api/panel replay recorded API responses and never touch the network.

    PANELBUILDER_CACHE=tests/fixtures .venv/bin/python -m pytest tests/test_api.py
"""
import dataclasses
import json
import os
import re
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

import primers  # noqa: E402
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
    # R3 reaches this endpoint on DISCLAIMER, verbatim (R8), not duplicated into the LD
    # verdict: reason renders on a card with room for one line.
    assert "phasing" in j["disclaimer"].lower()
    assert j["disclaimer"] == pb_disclaimer()
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


def test_health_states_the_primer_defaults_the_form_draws_from(monkeypatch):
    """The manual primer form is gated on this field and draws nothing without it.

    It shipped absent, so a whole section of the form was invisible on the live site while
    every test passed: the flag said primers were enabled and the numbers never arrived.
    """
    monkeypatch.setattr("primers.available", lambda: True)
    j = client.get("/api/health").json()
    assert j["primers_enabled"] is True
    d = j["primer_defaults"]
    assert isinstance(d, dict) and d, "the form has no numbers to draw"
    # The engine's own, not a copy: every knob it names, at the values it will use.
    assert set(d) == set(dataclasses.asdict(primers.DEFAULTS))
    assert d["opt_tm"] == primers.DEFAULTS.opt_tm

    # No primer3, no design to configure, so no form: null rather than numbers nothing reads.
    monkeypatch.setattr("primers.available", lambda: False)
    j = client.get("/api/health").json()
    assert j["primers_enabled"] is False
    assert j["primer_defaults"] is None


def test_bundling_the_ucsc_check_spends_the_verify_budget_not_the_build_one(monkeypatch):
    """The checkbox must not be the cheaper door to UCSC's quota.

    Both doors reach one published daily limit. The build budget is 20 a window and the
    verify budget is 4, so a bundled check charged only to the build budget would be a
    5x rate-limit bypass on someone else's server.
    """
    monkeypatch.setattr("app.main.PER_IP_VERIFIES", 0)
    r = client.post("/api/panel", json={"variant": GOLDEN, "verify_primers": True})
    assert r.status_code == 429, "a bundled check ignored the verification budget"
    assert "verification" in r.json()["detail"].lower()

    # And it is only spent when actually asked for: an ordinary build is untouched by it.
    called = {}
    monkeypatch.setattr("app.jobs.submit", lambda q, verify_primers=False: called.setdefault(
        "verify", verify_primers) or "job1")
    r = client.post("/api/panel", json={"variant": GOLDEN})
    assert r.status_code == 202, r.text
    assert called["verify"] is False


def test_a_bundled_check_does_not_hold_the_build_slot_while_it_waits(monkeypatch):
    """A build slot bounds BUILDS, and the bundled check runs for minutes.

    Held through the UCSC waits, two ticked boxes refused every visitor a panel at the
    default limit of 2: the whole app blocked on someone else's rate limit, by two people
    who were not building anything any more. Verification is bounded by the per-IP budget
    and by ispcr's process-wide gate; it needs no build slot to stay polite.
    """
    import threading

    verifying = threading.Event()
    let_go = threading.Event()

    def slow_verify(job, on_log):
        verifying.set()
        let_go.wait(timeout=10)          # stand in for the 15s-per-pair UCSC gate

    monkeypatch.setattr("app.jobs._verify_after_build", slow_verify)
    before = jobs.active()
    jid = jobs.submit(pb_query(), verify_primers=True)
    try:
        assert verifying.wait(timeout=60), "the bundled check never started"
        # The build is done and the check is still running: the slot must already be back.
        assert jobs.active() == before, \
            "a bundled check is holding a build slot while it waits on UCSC"
    finally:
        let_go.set()
    for _ in range(600):
        if jobs.get(jid).status != "running":
            break
        time.sleep(0.05)
    assert jobs.get(jid).status == "done"
    assert jobs.active() == before, "the slot was released twice or not at all"


def test_verify_primers_is_not_a_query_field(monkeypatch):
    """It says what to do after the build, so pb.StructuredQuery must never see it.

    That shape is the typed front door to build() and carries only what to build. A field
    it does not define reaches it as a TypeError at submit time, which is a 500 on a request
    that was perfectly valid.
    """
    seen = {}

    def fake_submit(q, verify_primers=False):
        seen["q"] = q
        seen["verify"] = verify_primers
        return "job1"

    monkeypatch.setattr("app.jobs.submit", fake_submit)
    r = client.post("/api/panel", json={"variant": GOLDEN, "verify_primers": True})
    assert r.status_code == 202, r.text
    assert seen["verify"] is True
    assert not hasattr(seen["q"], "verify_primers"), \
        "the flag leaked into the query shape build() reads"


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


def _await_job(job_id):
    for _ in range(600):
        s = client.get(f"/api/panel/{job_id}").json()
        if s["status"] != "running":
            return s
        time.sleep(0.1)
    raise AssertionError(f"job {job_id} never finished")


def test_model_chosen_variant_is_disclosed_from_api_nl_to_the_filed_export(monkeypatch,
                                                                          request):
    """The disclosure has to survive the whole hop, not just render when handed the flag.

    /api/nl names the model, the client echoes it into /api/panel, build() records it and
    every export says so. Each side landing separately is what makes this worth pinning:
    any one of them going quiet leaves a model-chosen panel that reads as hand-typed, and
    nothing else in the app would notice.
    """
    import app.nl

    monkeypatch.setattr(app.nl, "_llm_intent", lambda text: {"variant": GOLDEN})
    # Both ends: a warm entry would hide the stub, and the stub's answer must not outlive
    # the test in a cache that another test would read as a real reply.
    app.nl._cached_intent.cache_clear()
    request.addfinalizer(app.nl._cached_intent.cache_clear)

    j = client.post("/api/nl", json={"text": "the splice mutation we discussed"}).json()
    assert j["used_llm"] is True, "prose naming no identifier must take the model path"
    assert j["model"], "/api/nl must name the model the client has to echo back"

    # Exactly the body the client builds from that response.
    body = dict(j["query"], nl_text=j["text"], nl_model=j["model"])
    assert not {"chrom", "pos", "pos_grch38", "ref", "alt"} & set(body), "R1"
    job_id = client.post("/api/panel", json=body).json()["job_id"]
    assert _await_job(job_id)["status"] == "done"

    model = j["model"]
    # Read the strings a human sees, not the raw bytes: reportlab compresses the streams,
    # so a substring check against the file passes or fails for the wrong reason.
    from app.exports import _pdf_text
    pdf = _pdf_text(client.get(f"/api/export/{job_id}.pdf").content)
    assert model in pdf, "the PDF does not say a model chose the variant"
    assert "the splice mutation we discussed" in pdf, "the PDF drops the text it was given"
    for ext in ("csv", "json"):
        assert model in client.get(f"/api/export/{job_id}.{ext}").text, f"{ext} is silent"
    # A reader stripping '#' comments must still see it, so it rides a column, not a note.
    csv_body = client.get(f"/api/export/{job_id}.csv").text
    rows = [ln for ln in csv_body.splitlines() if ln and not ln.startswith("#")]
    assert model in rows[0] or model in rows[1], "comment-stripped CSV loses the model"


def test_a_typed_identifier_is_never_reported_as_a_model_choice():
    """Silence is the honest rendering of an absent model: no caveat, and no 'none'."""
    j = client.post("/api/nl", json={"text": f"panel for {GOLDEN}"}).json()
    assert j["used_llm"] is False and j["model"] is None

    job_id = client.post("/api/panel", json=dict(j["query"], nl_text=j["text"],
                                                 nl_model=j["model"])).json()["job_id"]
    assert _await_job(job_id)["status"] == "done"

    for ext in ("csv", "json", "pdf"):
        body = client.get(f"/api/export/{job_id}.{ext}").content.decode("latin-1").lower()
        assert "language model" not in body, f"{ext} caveats a panel the user typed"
        assert "nl_model: none" not in body, f"{ext} renders a None nl_model"


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


def test_build_log_rides_both_transports_and_leaks_nothing():
    """The log has to survive the poll path, not just SSE, and a late subscriber has to get
    all of it. The leak needles are the last gate before a line reaches a browser: jobs.py
    writes the content, but this is where it goes out.

    The log assertions arm themselves once jobs.py buffers lines; until then the shape
    contract (a list, absent field degrades to empty) is what holds.
    """
    job_id = client.post("/api/panel", json={"variant": GOLDEN}).json()["job_id"]
    for _ in range(600):
        if jobs.get(job_id).status != "running":
            break
        time.sleep(0.1)
    polled = client.get(f"/api/panel/{job_id}").json()
    assert polled["status"] == "done", polled.get("error")
    log = polled["log"]
    assert isinstance(log, list), "the poller must carry the log: proxies buffer SSE"

    # Subscribed AFTER the job finished: replay must hand over every line, not the tail.
    streamed, ev = [], None
    with client.stream("GET", f"/api/panel/{job_id}/stream") as r:
        for line in r.iter_lines():
            if line.startswith("event:"):
                ev = line.partition(":")[2].strip()
            elif line.startswith("data:") and ev == "log":
                streamed.append(json.loads(line.partition(":")[2]))
    assert streamed == log, "a late subscriber must replay the whole log, not join at the end"

    blob = json.dumps(log).lower()
    for needle in ("traceback", "app.jobs", "panelbuilder.py", "/users/",
                   "key=", "token=", "authorization"):
        assert needle not in blob, f"build log leaks {needle!r} to the browser: {blob[:300]}"


def test_the_browser_knows_every_tag_the_engine_emits():
    """LOG_TAGS in web/src/api.ts must list exactly panelbuilder.TAGS.

    The UI renders a tag it does not know as INFO, so a tag the engine adds here reaches
    the reader disguised as routine chatter rather than as itself, and nothing else fails.
    Neither side can see the other, which is why the check has to sit across them.
    """
    import panelbuilder as pb

    src = (ROOT / "web" / "src" / "api.ts").read_text()
    m = re.search(r"export const LOG_TAGS = \[(.*?)\] as const", src, re.S)
    assert m, "LOG_TAGS is not where this check expects it in web/src/api.ts"
    assert set(re.findall(r"'([A-Z]+)'", m.group(1))) == set(pb.TAGS), (
        "web/src/api.ts LOG_TAGS and panelbuilder.TAGS disagree: a tag the engine emits "
        "and the UI does not list renders as INFO"
    )


def test_the_stream_delivers_a_line_that_lands_as_the_job_finishes(monkeypatch):
    """A build that fails writes its last line and flips status with no event between.

    The stream reads status and the log at different moments, so a line landing between
    those two reads is the one a reader most needs: on the error path it names the
    failure. The poller cannot lose it (it re-reads the whole log), which is exactly why
    the two transports have to be compared, not just the surviving one.
    """
    import panelbuilder as pb

    class FlipsWhenStatusIsRead:
        """Appends its final line at the instant status is read: the writer landing in the
        window between the stream's two reads, made deterministic.

        It fires on the SECOND read, not the first, so the seeded event and line are
        already drained by then. A job with anything still undrained keeps the loop going
        for another pass, which hides the very race this is here to catch.
        """
        def __init__(self):
            self.id, self.stage, self.fraction = "raced", "resolve", 0.1
            self.events = [{"stage": "resolve", "fraction": 0.1}]
            self.log = [{"tag": pb.Tag.INFO, "text": "the build started"}]
            self.result, self.error, self._status = None, None, "running"
            self.reads = 0

        @property
        def status(self):
            self.reads += 1
            if self.reads > 1 and self._status == "running":
                self.error = "ApiError: a source refused"
                self.log.append({"tag": pb.Tag.WARN, "text": self.error})   # jobs.py order
                self._status = "error"
            return self._status

    job = FlipsWhenStatusIsRead()
    monkeypatch.setattr(jobs, "get", lambda _id: job)
    streamed, ev = [], None
    with client.stream("GET", "/api/panel/raced/stream") as r:
        for line in r.iter_lines():
            if line.startswith("event:"):
                ev = line.partition(":")[2].strip()
            elif line.startswith("data:") and ev == "log":
                streamed.append(json.loads(line.partition(":")[2]))
    assert streamed == job.log, "the stream dropped a line the poller would have kept"


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


def pb_query(**kw):
    """The golden case as a StructuredQuery, for tests that drive jobs.py directly."""
    import panelbuilder as pb
    return pb.StructuredQuery(variant=GOLDEN, **kw)


def test_log_tags_match_the_frontend():
    """The tag set is a wire contract, written twice: once in panelbuilder and once in
    web/src/api.ts, which colours them. A tag the engine emits and the UI has never heard
    of renders unstyled, so the drift is invisible until someone squints at a log.
    """
    import re
    from pathlib import Path

    import panelbuilder as pb

    api_ts = Path(__file__).resolve().parent.parent / "web" / "src" / "api.ts"
    if not api_ts.exists():
        pytest.skip("web/src/api.ts not present")
    m = re.search(r"export const LOG_TAGS = \[([^\]]+)\]", api_ts.read_text())
    assert m, "could not read LOG_TAGS out of api.ts; this check is stale"
    ui = {t.strip().strip("'\"") for t in m.group(1).split(",") if t.strip()}
    assert ui == set(pb.TAGS), (
        f"engine and UI disagree about the log tags.\n"
        f"  engine only: {sorted(set(pb.TAGS) - ui)}\n"
        f"  ui only    : {sorted(ui - set(pb.TAGS))}"
    )


def test_index_html_is_never_cached():
    """index.html names the content-hashed bundles, so a browser holding an old copy loads
    an old app: a deploy live on the server and invisible in the browser. The assets
    themselves are immutable by construction and must stay cacheable.
    """
    dist = Path(__file__).resolve().parent.parent / "web" / "dist"
    if not (dist / "index.html").exists():
        pytest.skip("web/dist not built")

    # "/deep-link" stands for any client route: it 404s in StaticFiles and falls back to
    # index.html, which must carry the header too. Not "/docs": FastAPI's Swagger owns that
    # path, and the app's own docs are the hash route "#/docs".
    for path in ("/", "/index.html", "/deep-link"):
        r = client.get(path)
        assert r.status_code == 200, (path, r.status_code)
        assert "no-store" in r.headers.get("cache-control", ""), (
            f"{path} served without no-store: {r.headers.get('cache-control')!r}")

    asset = next((p for p in (dist / "assets").glob("index-*.js")), None)
    if asset:
        r = client.get(f"/assets/{asset.name}")
        assert r.status_code == 200
        assert "no-store" not in r.headers.get("cache-control", ""), (
            "a content-hashed asset must stay cacheable; only the html is no-store")

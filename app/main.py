"""OriginMarker HTTP API.

Thin transport over panelbuilder. This module owns NO genetics: every coordinate, rsID,
strand and ref/alt on the wire comes out of pb.* (R1). Served as a subpage at
https://ezrakruger.cc/originmarker/ via ROOT_PATH, and at '/' for local dev.

    ROOT_PATH=/originmarker uvicorn app.main:app
    uvicorn app.main:app --reload            # local dev
"""
from __future__ import annotations

import asyncio
import json
import logging
import mimetypes
import os
import threading
import time
from dataclasses import asdict
from functools import lru_cache
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from starlette.exceptions import HTTPException as StarletteHTTPException

import build_info
import genetic_map
import panelbuilder as pb
from app import jobs

VERSION = build_info.VERSION
DIST = Path(__file__).resolve().parent.parent / "web" / "dist"

_log = logging.getLogger(__name__)

MEDIA = {
    "csv": "text/csv",
    "json": "application/json",
    "pdf": "application/pdf",
    "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
}

app = FastAPI(title="OriginMarker", version=VERSION,
              root_path=os.environ.get("ROOT_PATH", ""),
              description="Candidate linkage-marker panels for PGT-M. " + pb.DISCLAIMER)

# Same-origin needs no CORS header; ALLOWED_ORIGINS exists only so the Vite dev server
# on :5173 can talk to a local API.
if _origins := [o.strip() for o in os.environ.get("ALLOWED_ORIGINS", "").split(",") if o.strip()]:
    from fastapi.middleware.cors import CORSMiddleware
    app.add_middleware(CORSMiddleware, allow_origins=_origins,
                       allow_methods=["GET", "POST"], allow_headers=["*"])


@app.on_event("startup")
def _startup() -> None:
    threading.Thread(target=_warm_ensembl_release, daemon=True).start()


# --------------------------------------------------------------------------- #
# Request models (validation at the trust boundary; pb.StructuredQuery re-checks)
# --------------------------------------------------------------------------- #

class ResolveIn(BaseModel):
    variant: str = Field(min_length=1, max_length=200)
    build: str = "GRCh38"


class PanelIn(BaseModel):
    variant: str = Field(min_length=1, max_length=200)
    gene: Optional[str] = Field(default=None, max_length=40)
    window_bp: int = pb.DEFAULT_WINDOW
    build: str = "GRCh38"
    ancestry: Optional[str] = None
    common_maf: float = pb.COMMON_MAF
    cross_check: bool = True
    # Echoed back by the client from /api/nl so the exports can disclose that a model, not
    # the user, chose the variant. Provenance only: the build reads neither. Capped here
    # because both are rendered into the PDF; nl_text matches NLIn.text's ceiling.
    nl_text: Optional[str] = Field(default=None, max_length=2000)
    nl_model: Optional[str] = Field(default=None, max_length=60)


class NLIn(BaseModel):
    text: str = Field(min_length=1, max_length=2000)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _lazy(name: str):
    """Import an optional sibling module at first use.

    These carry dependencies the core API does not need, so an import failure costs one
    feature rather than the whole API. Raises HTTPException(503); the cause goes to the
    log, never to the visitor.
    """
    try:
        mod = __import__(f"app.{name}", fromlist=[name])
    except Exception:  # noqa: BLE001
        _log.exception("optional module app.%s failed to import", name)
        raise HTTPException(503, "That part of the app is unavailable right now.")
    return mod


def _ensembl_release() -> Optional[int]:
    """Ensembl release for the provenance stamp (R6). Never blocks the request.

    None until the background warm lands; callers must render that as unknown. Reads
    panelbuilder's copy: a second copy here would let a frozen panel render the live
    server's release rather than its own.
    """
    return pb.ensembl_release()


def _warm_ensembl_release() -> None:
    """Fetch the Ensembl release once, retrying on a widening interval until it lands.

    /info/data is slow and rate-limited, so a single attempt at startup usually fails.
    Never cached to disk: the release changes about quarterly, and a stale stamped version
    is worse than none.
    """
    for delay in (0, 30, 120, 300):
        if delay:
            time.sleep(delay)
        try:
            j = json.loads(pb._get(pb.ENSEMBL, "/info/data", {"content-type": "application/json"},
                                   tries=1, timeout=25, use_cache=False))
            pb.set_ensembl_release(int(j["releases"][0]))
            return
        except Exception:  # noqa: BLE001 - provenance nicety, never fatal
            continue


def _client_ip(request: Request) -> str:
    """Best available identity for rate limiting.

    Order matters. A proxy APPENDS to X-Forwarded-For, so the first hop is
    attacker-controlled and XFF[0] lets a client rotate its own bucket. CF-Connecting-IP is
    written by Cloudflare and overwrites what the client sends; failing that, the LAST XFF
    hop is the one our nearest trusted proxy added.
    """
    cf = request.headers.get("cf-connecting-ip", "").strip()
    if cf:
        return cf
    xff = request.headers.get("x-forwarded-for", "")
    hops = [h.strip() for h in xff.split(",") if h.strip()]
    if hops:
        return hops[-1]
    return request.client.host if request.client else "?"


PER_IP_MAX = int(os.environ.get("PER_IP_BUILDS", "20"))
PER_IP_RESOLVES = int(os.environ.get("PER_IP_RESOLVES", "60"))
# Free-text is the only path that spends money: tighter per-IP bucket, plus a global
# ceiling that does not depend on identifying the caller at all.
PER_IP_NL = int(os.environ.get("PER_IP_NL", "10"))
GLOBAL_NL_MAX = int(os.environ.get("GLOBAL_NL_MAX", "120"))
PER_IP_WINDOW = 600.0
_ip_hits: dict[str, list[float]] = {}


def _rate_ok(key: str, limit: int) -> bool:
    now = time.time()
    if len(_ip_hits) > 1000:
        # Evict only what has aged out, and never the global key: clearing the table here
        # would let anyone reset every bucket by sending 1000 requests with distinct keys.
        for k in [k for k, v in _ip_hits.items()
                  if k != "nl:global" and not any(now - t < PER_IP_WINDOW for t in v)]:
            del _ip_hits[k]
        if len(_ip_hits) > 5000:             # still pathological: shed all but the global
            keep = _ip_hits.get("nl:global", [])
            _ip_hits.clear()
            _ip_hits["nl:global"] = keep
    hits = [t for t in _ip_hits.get(key, []) if now - t < PER_IP_WINDOW]
    _ip_hits[key] = hits
    if len(hits) >= limit:
        return False
    hits.append(now)
    return True


def _sse(event: str, data: dict) -> str:
    return f"event: {event}\ndata: {json.dumps(data)}\n\n"


# --------------------------------------------------------------------------- #
# API
# --------------------------------------------------------------------------- #

@app.get("/api/health")
def health():
    try:
        ldlink_enabled = _lazy("ldlink").available()   # ldlink owns what "configured" means
    except Exception:  # noqa: BLE001 - health must never be the thing that's down
        ldlink_enabled = False
    return {"ok": True, "version": VERSION,
            # "release" is the app version; "build" below is the GENOME build. Two senses
            # of the word, two keys: never merge them.
            "release": build_info.BUILD,
            "release_codename": build_info.CODENAME,
            "release_gloss": build_info.current().gloss,
            "gnomad_dataset": "gnomad_r4",
            "build": "GRCh38", "ensembl_release": _ensembl_release(),
            "map_source": genetic_map.MAP_SOURCE,
            "ldlink_enabled": ldlink_enabled,
            # R8: canonical wording ships from here; the frontend must never paraphrase it.
            "disclaimer": pb.DISCLAIMER,
            "layer_b_steps": pb.LAYER_B_STEPS,
            # Flags the LLM fallback only: the rsID/HGVS path works without a key.
            "nl_enabled": bool(os.environ.get("ANTHROPIC_API_KEY"))}


@app.post("/api/resolve")
def resolve(body: ResolveIn, request: Request):
    """Cheap pre-flight: ~3 cached calls. Fails loudly rather than guessing (R1).

    Rate-limited despite being cheap: an unresolvable variant is not. pb._http retries
    Ensembl's 404 with backoff, so a typo'd rsID holds a threadpool thread for ~90s and a
    few dozen of them wedge the API.
    """
    if not _rate_ok(f"resolve:{_client_ip(request)}", PER_IP_RESOLVES):
        raise HTTPException(429, f"Rate limit: {PER_IP_RESOLVES} lookups per "
                                 f"{int(PER_IP_WINDOW / 60)} minutes per client.")
    try:
        v = pb.resolve_variant(body.variant.strip(), build=body.build)
    except pb.ApiError as e:
        # pb.resolve_variant owns this wording: it says what was wrong and what to type
        # instead, so pass it through rather than wrapping it.
        raise HTTPException(400, str(e))
    except Exception:  # noqa: BLE001 - unexpected: keep it generic, leak nothing
        _log.exception("unexpected failure resolving %r on %s", body.variant, body.build)
        raise HTTPException(400, f"Could not resolve {body.variant!r} on {body.build}.")
    r = pb.assess_rarity(v)
    acc = v.clinvar_accession
    return {"variant": asdict(v), "rarity": asdict(r),
            "transcript_sense": v.transcript_sense_change(),          # R7
            # The R2/R3 verdict is worded once, in assess_rarity. Never hand-write a second
            # copy here: the UI renders {ld_banner || rarity.reason}, so a copy always wins.
            "ld_banner": r.reason,
            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{acc}/" if acc else ""}


@app.post("/api/panel", status_code=202)
async def panel(body: PanelIn, request: Request):
    # async: submit() is instant; the blocking pb.build() runs in jobs.py's own pool.
    try:
        q = pb.StructuredQuery(**body.model_dump())
    except ValueError as e:
        raise HTTPException(400, str(e))
    if not _rate_ok(f"panel:{_client_ip(request)}", PER_IP_MAX):
        raise HTTPException(429, f"Rate limit: {PER_IP_MAX} panel builds per "
                                 f"{int(PER_IP_WINDOW / 60)} minutes per client. A build "
                                 f"issues ~25 gnomAD queries; please pace them.")
    try:
        return {"job_id": jobs.submit(q)}
    except jobs.Busy:
        raise HTTPException(429, f"The server is already building its maximum of "
                                 f"{jobs.MAX_CONCURRENT} panels at once. Retry in a "
                                 f"minute.")


@app.get("/api/panel/{job_id}")
def panel_status(job_id: str):
    if (job := jobs.get(job_id)) is None:
        raise HTTPException(404, f"unknown or expired job {job_id!r}")
    return {"status": job.status, "stage": job.stage, "fraction": job.fraction,
            "result": job.result.to_dict() if job.result else None,   # serialise at the edge
            "error": job.error}


@app.get("/api/panel/{job_id}/stream")
async def panel_stream(job_id: str, request: Request):
    if (job := jobs.get(job_id)) is None:
        raise HTTPException(404, f"unknown or expired job {job_id!r}")

    async def gen():
        yield _sse("progress", {"stage": job.stage, "fraction": job.fraction})
        sent, last = 0, time.monotonic()
        while True:
            if sent < len(job.events):                  # drain before testing status, so a
                for ev in job.events[sent:]:            # trailing event is never lost
                    sent += 1
                    yield _sse("progress", ev)
                last = time.monotonic()
            elif job.status != "running":
                break
            else:
                if await request.is_disconnected():
                    return
                if time.monotonic() - last > 15:
                    # Cloudflare Tunnel drops idle streams; an SSE comment line is a no-op
                    # to the client but keeps the socket warm.
                    last = time.monotonic()
                    yield ": keepalive\n\n"
                await asyncio.sleep(0.1)
        if job.status == "done":
            yield _sse("done", {"job_id": job.id})
        else:
            yield _sse("error", {"message": job.error or "build failed"})

    return StreamingResponse(gen(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache",
                                      "X-Accel-Buffering": "no",
                                      "Connection": "keep-alive"})


@app.get("/api/export/{job_id}.{ext}")
def export(job_id: str, ext: str):
    if ext not in MEDIA:
        raise HTTPException(404, f"unsupported export format {ext!r}; expected one of "
                                 f"{sorted(MEDIA)}")
    if (job := jobs.get(job_id)) is None:
        raise HTTPException(404, f"unknown or expired job {job_id!r}")
    if job.status != "done":
        raise HTTPException(409, f"job {job_id} is {job.status}; nothing to export yet")

    exports = _lazy("exports")
    body = getattr(exports, f"to_{ext}")(job.result)      # exports takes the PanelResult
    name = exports.FILENAME(job.result, ext)              # stamps build + date (R6)
    return Response(body, media_type=MEDIA[ext],
                    headers={"Content-Disposition": f'attachment; filename="{name}"'})


@app.post("/api/nl")
def nl(body: NLIn, request: Request):
    """Intent only. The NL layer cannot emit a coordinate: StructuredQuery has no field
    for one, so R1 holds by construction rather than by trust.

    Rate limited per-IP and globally, being the only endpoint that spends money. The global
    cap is the one that matters: per-IP limiting rests on identifying the client, and
    identity on a public URL is only as good as the proxy in front of it.
    """
    mod = _lazy("nl")
    ip = _client_ip(request)
    # Both meters fire BEFORE parse(), because parse() is where the money leaves. Only text
    # that would actually reach the model is metered; anything carrying an rsID or HGVS is
    # answered by regex and costs nothing.
    if mod.needs_llm(body.text):
        if not _rate_ok(f"nl:{ip}", PER_IP_NL):
            raise HTTPException(429, f"Rate limit: {PER_IP_NL} free-text queries per "
                                     f"{int(PER_IP_WINDOW / 60)} minutes per client. Enter "
                                     f"an rsID or HGVS expression instead; those are not "
                                     f"limited.")
        if not _rate_ok("nl:global", GLOBAL_NL_MAX):
            raise HTTPException(429, "Free-text parsing is temporarily at capacity. Enter "
                                     "an rsID or HGVS expression; those are always "
                                     "available.")
    try:
        # By attribute, never by unpacking: nl.parse owns its return shape, and an exact
        # unpack turns a field being added there into a 400 whose detail is a Python unpack
        # error, since the ValueError below reads as the user having sent a bad query.
        p = mod.parse(body.text.strip())
    except HTTPException:
        raise
    except ValueError as e:
        # nl.parse owns this wording: it says what was wrong and what to send instead.
        raise HTTPException(400, str(e))
    except Exception:  # noqa: BLE001 - unexpected: log the cause, show the way out
        _log.exception("unexpected failure parsing free text")
        raise HTTPException(400, "That text could not be read as a variant request. Enter "
                                 "an rsID or HGVS expression, for example rs151344623.")
    # text/model are what the client echoes into /api/panel so an export can disclose that
    # a model chose the variant. used_llm is the sole authority on that: naming the model
    # off the regex path would caveat a panel the user typed word for word. Neither field
    # is a coordinate, and there is still no field here that could carry one.
    return {"query": asdict(p.query), "used_llm": bool(p.used_llm), "note": p.note or "",
            "text": body.text.strip(),
            "model": mod.MODEL if p.used_llm else None,
            "named_genes": list(getattr(p, "named_genes", []))}


@app.get("/api/genes")
def genes(q: str = ""):
    if len(q.strip()) < 2:
        return []
    return _gene_lookup(q.strip().upper())


@lru_cache(maxsize=512)
def _gene_lookup(symbol: str) -> list[dict]:
    # Ensembl REST has no prefix search (/xrefs/name/ABC returns []), so this is
    # exact-symbol only.
    try:
        g = json.loads(pb._get(pb.ENSEMBL, f"/lookup/symbol/homo_sapiens/{symbol}",
                               {"content-type": "application/json"}, tries=1, timeout=5))
    except Exception:  # noqa: BLE001 - 404 for anything that isn't a symbol
        return []
    return [{"symbol": g.get("display_name") or symbol, "description": g.get("description") or ""}]


@app.get("/api/ld")
def ld(a: str, b: str, pop: str = "CEU"):
    """Optional annotation BETWEEN TWO COMMON SNPS only. Never a ranking key, never an
    origin call (R2). The token stays in the ldlink module's env, never on the wire.

    ld_between() degrades rather than raising, so the {available: False} shape, not an
    exception, is what becomes the 503. allow_rare is deliberately not plumbed through:
    it is ldlink's R2 tripwire and has no correct value other than the default False.
    """
    try:
        out = _lazy("ldlink").ld_between(a, b, pop)
    except ValueError as e:                  # the R2 guard
        raise HTTPException(400, str(e))
    if not out.get("available"):
        raise HTTPException(503, out.get("note") or "LD annotation unavailable")
    return {k: out[k] for k in ("r2", "dprime", "pop", "note", "caveat")}


# --------------------------------------------------------------------------- #
# SPA (mounted last: /api wins every collision)
# --------------------------------------------------------------------------- #

# Python only learned .woff2 in 3.13; the runtime image is 3.12, where StaticFiles would
# otherwise guess application/octet-stream.
mimetypes.add_type("font/woff2", ".woff2")
mimetypes.add_type("font/woff", ".woff")


class SPAStatic(StaticFiles):
    """Serve index.html for unknown *routes* so client-side deep links work.

    StaticFiles raises HTTPException(404) rather than returning one, so this must catch
    rather than inspect a status. Catch Starlette's, not FastAPI's subclass.

    Only extensionless paths get the fallback: a blanket "any 404 -> index.html" answers a
    missing .js with 200 + HTML, which looks healthy to curl while the browser fails to
    parse it.
    """

    async def get_response(self, path, scope):
        try:
            return await super().get_response(path, scope)
        except StarletteHTTPException as e:
            if e.status_code != 404:
                raise
            if Path(path).suffix and Path(path).suffix != ".html":
                raise
            return await super().get_response("index.html", scope)


if (DIST / "index.html").exists():
    app.mount("/", SPAStatic(directory=DIST, html=True), name="spa")
else:
    # Checked at import, not per request: restart the server once web/dist is built.
    @app.get("/", response_class=HTMLResponse)
    def placeholder():
        return ("<h1>OriginMarker API</h1><p>The frontend bundle is not built yet "
                "(<code>web/dist/index.html</code> missing). The API is live at "
                "<a href='api/health'>api/health</a> and <code>api/docs</code>.</p>"
                f"<hr><footer><small>{pb.DISCLAIMER}</small></footer>")   # R8


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=int(os.environ.get("PORT", "8000")))

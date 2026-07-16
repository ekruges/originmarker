"""In-process panel job registry.

pb.build() is a long, blocking, network-bound call. This runs it in a thread pool and
records its on_progress and on_log callbacks so the SSE endpoint can stream them. Both
buffers are append-only, so a subscriber arriving mid-build replays from the start. Jobs
are TTL-evicted and capped so a public URL cannot grow memory forever.

Registry state is per-process, so this assumes a single uvicorn worker.
"""
from __future__ import annotations

import logging
import os
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Optional

import panelbuilder as pb

_log = logging.getLogger(__name__)

TTL_S = 3600                                                    # evict finished jobs after ~1h
MAX_JOBS = 200                                                  # hard cap on registry size
MAX_CONCURRENT = int(os.environ.get("MAX_CONCURRENT_BUILDS", "2"))
# ^ protects gnomAD's API (and the owner's IP): a build fans out ~25 concurrent region
#   queries, so 2 in flight is already 50 sockets at Broad.
MAX_LOG_LINES = 2000
# ^ a build emits about two lines per region chunk plus a handful per stage, so a default
#   window stays well under this and a wide one still fits. The cap is what bounds the
#   registry: MAX_JOBS x MAX_LOG_LINES short strings is the worst it can hold.


class Busy(RuntimeError):
    """Raised by submit() when MAX_CONCURRENT builds are already running."""


@dataclass
class Job:
    id: str
    # Append-only progress log, NOT an asyncio.Queue: a list is loop-agnostic
    # (list.append is atomic under the GIL) and lets several clients tail the same job by
    # index. A Queue would bind the job to whichever loop called submit().
    events: list = field(default_factory=list)
    # The build log, buffered like events and tailed the same way: append-only, so a late
    # subscriber replays from index 0 and gets the whole account. Capped by REFUSING further
    # lines rather than by dropping the oldest, since a subscriber tracks its position by
    # index and dropping from the front would silently renumber every line it has not read.
    log: list = field(default_factory=list)
    # pb.build() calls on_log from its region worker pool, so several threads append here.
    # list.append is atomic, but the cap is a read-then-append and the truncation notice
    # must be written exactly once.
    log_lock: threading.Lock = field(default_factory=threading.Lock)
    status: str = "running"                 # 'running' | 'done' | 'error'
    stage: str = "queued"
    fraction: float = 0.0
    # The live PanelResult, NOT to_dict(): app/exports.py walks Marker objects and calls
    # result.variant.transcript_sense_change(). The API serialises at the edge instead.
    result: Optional[pb.PanelResult] = None
    error: Optional[str] = None
    created: float = field(default_factory=time.time)


_jobs: dict[str, Job] = {}
_lock = threading.Lock()
_active = 0
_pool = ThreadPoolExecutor(max_workers=MAX_CONCURRENT, thread_name_prefix="panel")


def get(job_id: str) -> Optional[Job]:
    return _jobs.get(job_id)


def active() -> int:
    return _active


def submit(q: pb.StructuredQuery) -> str:
    """Start a build. Raises Busy when MAX_CONCURRENT builds are already in flight."""
    global _active
    with _lock:
        _evict_locked()
        if _active >= MAX_CONCURRENT:
            raise Busy(f"{_active} panel builds already running (limit {MAX_CONCURRENT})")
        _active += 1
        job = Job(id=uuid.uuid4().hex[:12])
        _jobs[job.id] = job
    _pool.submit(_run, job, q)
    return job.id


def _run(job: Job, q: pb.StructuredQuery) -> None:
    """Worker thread. Never raises: every outcome lands on the Job."""
    global _active

    def on_progress(stage: str, frac: float) -> None:
        job.stage, job.fraction = stage, frac
        job.events.append({"stage": stage, "fraction": round(frac, 3)})

    def on_log(tag: str, text: str) -> None:
        with job.log_lock:
            if len(job.log) < MAX_LOG_LINES:
                job.log.append({"tag": tag, "text": text})
            elif len(job.log) == MAX_LOG_LINES:
                # Says so rather than trailing off: a log that stops mid-build reads as a
                # build that stopped mid-build.
                job.log.append({"tag": pb.Tag.WARN,
                                "text": f"Build log truncated at {MAX_LOG_LINES} lines. "
                                        "The build itself continues."})

    try:
        job.result = pb.build(q, on_progress=on_progress, on_log=on_log)
        # status flips LAST, after every event is appended, so a reader that drains events
        # before testing status cannot miss a trailing one.
        job.status = "done"
    except (pb.ApiError, ValueError) as e:
        # pb writes these for the person who typed the query, so they go through verbatim.
        _log.warning("panel build %s failed: %s", job.id, e)
        job.error = str(e)
        # A failed build's log otherwise just stops, which reads like a lost connection
        # rather than a refusal: pb only emits DONE on the way out. Appended BEFORE the
        # status flip, for the same reason the events are.
        on_log(pb.Tag.WARN, job.error)
        job.status = "error"
    except Exception:  # noqa: BLE001
        # job.error reaches the browser: keep it generic, and leave the plumbing in the log.
        _log.exception("panel build %s failed", job.id)
        job.error = "The panel build did not finish. Retry in a moment."
        on_log(pb.Tag.WARN, job.error)
        job.status = "error"
    finally:
        with _lock:
            _active -= 1


def _evict_locked() -> None:
    now = time.time()
    dead = [jid for jid, j in _jobs.items()
            if j.status != "running" and now - j.created > TTL_S]
    over = len(_jobs) - len(dead) - MAX_JOBS
    if over > 0:
        oldest = sorted((j for j in _jobs.values()
                         if j.status != "running" and j.id not in dead),
                        key=lambda j: j.created)
        dead += [j.id for j in oldest[:over]]
    for jid in dead:
        _jobs.pop(jid, None)

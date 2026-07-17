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
import dataclasses
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Optional

import panelbuilder as pb
import primers
from app import ispcr

_log = logging.getLogger(__name__)

TTL_S = 3600                                                    # evict finished jobs after ~1h
MAX_JOBS = 200                                                  # hard cap on registry size
MAX_CONCURRENT = int(os.environ.get("MAX_CONCURRENT_BUILDS", "2"))
# ^ protects gnomAD's API (and the owner's IP): a build fans out ~25 concurrent region
#   queries, so 2 in flight is already 50 sockets at Broad.
MAX_LOG_LINES = 2000
# Verify jobs run one at a time, so the rest queue. Queued jobs are 'running' and eviction
# never touches those, so without a cap a public URL grows the registry without bound and
# builds a queue whose tail waits an hour. It is also someone else's budget: UCSC publishes
# 5,000 requests a day for the whole server.
MAX_QUEUED_VERIFIES = 4
# ^ a build emits about two lines per region chunk plus a handful per stage, so a default
#   window stays well under this and a wide one still fits. The cap is what bounds the
#   registry: MAX_JOBS x MAX_LOG_LINES short strings is the worst it can hold.


class Busy(RuntimeError):
    """Raised by submit() when MAX_CONCURRENT builds are already running."""


@dataclass
class Job:
    id: str
    # 'panel' | 'verify'. One registry holds both, so the TTL, the cap and the eviction
    # bound them together rather than twice. The two kinds fill DIFFERENT result fields,
    # and a reader must check this before reading either: handing a verify id to the panel
    # endpoint otherwise reads a None result as a build that produced nothing.
    kind: str = "panel"
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
    # A verify job's own results, keyed by rsID: {rsid: ispcr.verify() dict}. Its own field
    # rather than a widened `result`, so neither kind is typed as the other.
    #
    # These are ALSO written onto the panel's markers, by the loop in _run_verify. This once
    # said they never were, on the grounds that mutating a panel a user has already exported
    # turns one job id into two different documents. The concern was real and the conclusion
    # was wrong: kept here alone, a verdict reached this endpoint and the build log and never
    # the PDF anyone files, so the document said NOT CHECKED about a pair UCSC had already
    # called dangerous. Two exports of one job id do differ, and that is a pending answer
    # arriving rather than the panel contradicting itself: each export stamps its own
    # built_utc, and the earlier one's NOT CHECKED was true when it was taken.
    verdicts: dict = field(default_factory=dict)
    error: Optional[str] = None
    created: float = field(default_factory=time.time)


_jobs: dict[str, Job] = {}
_lock = threading.Lock()
_active = 0
_pool = ThreadPoolExecutor(max_workers=MAX_CONCURRENT, thread_name_prefix="panel")
# ONE worker, and not for this pool's sake: UCSC's published limit of one request per 15s
# belongs to UCSC, so it is global to this process. Any concurrency here breaches it no
# matter how the pool is sized. It is also why verification never runs inside a build: a
# 3-minute hold would occupy half of MAX_CONCURRENT and make the server look dead.
_verify_pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="verify")


def _log_fn(job: Job):
    """The capped, locked appender for a job's build log. Shared by both kinds, so the cap
    that bounds the registry is enforced in one place."""
    def on_log(tag: str, text: str) -> None:
        with job.log_lock:
            if len(job.log) < MAX_LOG_LINES:
                job.log.append({"tag": tag, "text": text})
            elif len(job.log) == MAX_LOG_LINES:
                # Says so rather than trailing off: a log that stops mid-build reads as a
                # build that stopped mid-build.
                job.log.append({"tag": pb.Tag.WARN,
                                "text": f"Build log truncated at {MAX_LOG_LINES} lines. "
                                        f"The build itself continues."})
    return on_log


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

    on_log = _log_fn(job)

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


def verifiable(panel: Job) -> list:
    """The finished panel's markers that have a pair to check, in shortlist order."""
    if panel.result is None:
        return []
    return [m for m in panel.result.recommended if m.primer and m.primer.ok]


def submit_verify(panel_job_id: str, rsids: Optional[list] = None) -> str:
    """Check a finished panel's primer pairs against the whole genome. Returns a job id.

    Its own job, never part of the build: at UCSC's published one request per 15 seconds a
    shortlist takes minutes, and a panel is useful before that. `rsids` selects markers;
    None means every pair the panel designed.
    """
    panel = _jobs.get(panel_job_id)
    if panel is None or panel.kind != "panel":
        raise KeyError(panel_job_id)
    if panel.status != "done":
        raise ValueError(f"panel {panel_job_id} is {panel.status}; nothing to verify yet")
    targets = verifiable(panel)
    if rsids is not None:
        want = set(rsids)
        targets = [m for m in targets if m.rsid in want]
        if unknown := want - {m.rsid for m in targets}:
            raise ValueError(f"no designed primer pair on this panel for "
                             f"{', '.join(sorted(unknown))}")
    if not targets:
        raise ValueError("this panel has no designed primer pair to verify")
    with _lock:
        _evict_locked()
        queued = sum(1 for j in _jobs.values()
                     if j.kind == "verify" and j.status == "running")
        if queued >= MAX_QUEUED_VERIFIES:
            raise Busy(f"{queued} primer verifications already queued (limit "
                       f"{MAX_QUEUED_VERIFIES}). They run one at a time at UCSC's published "
                       f"rate of one request every {ispcr.MIN_INTERVAL_S:.0f}s.")
        job = Job(id=uuid.uuid4().hex[:12], kind="verify")
        _jobs[job.id] = job
    _verify_pool.submit(_run_verify, job, targets)
    return job.id


def _run_verify(job: Job, targets: list) -> None:
    """Worker thread. Never raises: every outcome lands on the Job."""
    on_log = _log_fn(job)
    # ispcr raises its lines through pb._emit, which reads a ContextVar, so the sink must be
    # set on THIS thread or every line it writes lands in a context nobody reads.
    pb._log_sink.set(on_log)
    try:
        job.stage = "verifying"
        on_log(pb.Tag.INFO, f"UCSC In-Silico PCR: checking {len(targets)} pairs against "
                            f"{ispcr.DB} at one request per "
                            f"{ispcr.MIN_INTERVAL_S:.0f}s, about "
                            f"{len(targets) * ispcr.MIN_INTERVAL_S / 60:.0f} min")
        for i, m in enumerate(targets, 1):
            # "chr11", not "11": ispcr compares through its own _norm either way, but it
            # prints this string back at the reader. max_product is left at ispcr's own
            # default rather than lowered to the design's: it is the largest product UCSC
            # will REPORT, so a lower one hides a long spurious amplicon and reports the
            # pair as clean.
            out = ispcr.verify(m.primer.fwd.seq, m.primer.rev.seq, f"chr{m.chrom}",
                               m.primer.product_size)
            job.verdicts[m.rsid] = out
            # Written onto the panel too, or the verdict reaches the log and this endpoint
            # and never the document anyone files. An export taken later carries more than
            # one taken earlier, which is a pending answer arriving, not the same panel
            # saying two things: every export stamps its own built_utc, and an export taken
            # before this lands says NOT CHECKED, which is true at that moment.
            # replace(), not assignment: PrimerResult is frozen, which is what keeps a
            # primer from being prised apart from its warnings. A verdict is new information
            # about the pair, so it makes a new result rather than mutating one.
            #
            # The WARNINGS are replaced, not just the state. They opened with "NOT CHECKED
            # AGAINST THE GENOME", which this call is the checking of: left alone they would
            # contradict the verdict on the same row, and the reader is being asked to trust
            # exactly one of them. The verdict's own words go on instead, and the caveat with
            # them, because a pass here still is not a wet-lab result and the pair must not
            # arrive with nothing said about it.
            m.primer = dataclasses.replace(
                m.primer,
                insilico_pcr=out["state"],
                warnings=(
                    primers.Note(code=out["state"], short=out["short"], long=out["note"]),
                    primers.Note(code="ispcr_caveat", short=ispcr.CAVEAT_SHORT,
                                 long=ispcr.CAVEAT),
                ),
            )
            # ispcr owns every word of this verdict; nothing here restates it. A pass is one
            # line, and everything else is a WARN, because only one of its four states is a
            # pass and the other three all mean NOT VERIFIED.
            on_log(pb.Tag.INFO if out["state"] == ispcr.ONE_PRODUCT else pb.Tag.WARN,
                   f"{m.rsid}: {out['note']}")
            job.stage, job.fraction = f"verified {i}/{len(targets)}", i / len(targets)
            job.events.append({"stage": job.stage, "fraction": round(job.fraction, 3)})
        danger = sum(1 for v in job.verdicts.values() if v["state"] == ispcr.DANGER)
        clean = sum(1 for v in job.verdicts.values() if v["state"] == ispcr.ONE_PRODUCT)
        on_log(pb.Tag.DONE, f"{clean} of {len(targets)} pairs gave one product in "
                            f"{ispcr.DB}; {danger} are dangerous; "
                            f"{len(targets) - clean - danger} could not be checked. "
                            f"{ispcr.CAVEAT}")
        job.status = "done"
    except Exception:  # noqa: BLE001
        # job.error reaches the browser: keep it generic, and leave the plumbing in the log.
        # ispcr degrades rather than raising, so reaching here is a bug, not UCSC being down.
        _log.exception("primer verification %s failed", job.id)
        job.error = "Primer verification did not finish. Retry in a moment."
        on_log(pb.Tag.WARN, job.error)
        job.status = "error"


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

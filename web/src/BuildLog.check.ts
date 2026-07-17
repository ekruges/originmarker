// Self-check for the build-log wire contract. Run: node src/BuildLog.check.ts
// The component needs a DOM and is not imported here; this covers the normaliser that
// decides what reaches it, and the pure text composer that builds the .txt export.
import assert from 'node:assert/strict'
import { asLogLine, LOG_TAGS } from './api.ts'
import { buildLogText, type ExportCtx, type LogMeta } from './logfile.ts'

// app/jobs.py owns the shape and app/main.py passes it through unreshaped, so the frames
// below are what the wire is expected to carry.
assert.deepEqual(asLogLine({ tag: 'FETCH', text: 'gnomAD region chr11 chunk 3/12' }),
                 { tag: 'FETCH', text: 'gnomAD region chr11 chunk 3/12' })
for (const tag of LOG_TAGS) assert.equal(asLogLine({ tag, text: 'x' })?.tag, tag)

// A progress frame is not a log line. Both ride `event: progress` on some paths, and one
// read as the other would either blank the stage or invent a log entry.
assert.equal(asLogLine({ stage: 'fetching gnomAD', fraction: 0.4 }), null)
assert.equal(asLogLine(null), null)
assert.equal(asLogLine('[INFO] a bare string'), null)
assert.equal(asLogLine({ tag: 'INFO' }), null)

// A tag is a hint about a line, never a licence to drop one: an unknown tag, or main.py's
// {"text": ...} fallback for a non-mapping entry, must still render. INFO is the floor.
assert.equal(asLogLine({ tag: 'RETRY', text: 'gnomAD 429, backing off' })?.tag, 'INFO')
assert.equal(asLogLine({ text: 'no tag at all' })?.tag, 'INFO')
assert.equal(asLogLine({ tag: '', text: 'empty tag' })?.tag, 'INFO')
assert.equal(asLogLine({ tag: 'warn', text: 'lowercased on the wire' })?.tag, 'WARN')
assert.equal(asLogLine({ tag: 'FETCH', message: 'message, not text' })?.text, 'message, not text')

// The text is the line: never truncated, never re-worded, and an empty one is still a line.
const long = 'chr11:17,147,055-17,647,055 '.repeat(40)
assert.equal(asLogLine({ tag: 'FETCH', text: long })?.text, long)
assert.deepEqual(asLogLine({ tag: 'INFO', text: '' }), { tag: 'INFO', text: '' })

// --- the .txt export composer ---------------------------------------------------------
// Every line, verbatim, tagged. Then ONE footer line identifying where and when it came
// from, so a log in a bug report is not anonymous.
const CTX: ExportCtx = {
  now: '2026-07-17T20:00:00.000Z',
  url: 'https://ezrakruger.cc/originmarker/#/',
  agent: 'TestAgent/1.0',
}
const META: LogMeta = {
  release: 'Build 1.3.3.1 "Diakinesis"',
  jobId: 'abc123',
  provenance: {
    build: 'GRCh38', ensembl_release: 116, queried_utc: '2026-07-17T19:59:00Z',
    window_bp: 250000, common_maf: 0.05, candidate_n: 1202, elapsed_s: 1.3,
    sources: { gnomad: 'gnomad_r4', genetic_map: 'deCODE 2019' },
  },
}
const SAMPLE: import('./api.ts').LogLine[] = [
  { tag: 'INFO', text: '20 markers shortlisted' },
  { tag: 'WARN', text: 'rs886288: no primer pair' },
]
const out = buildLogText(SAMPLE, META, CTX)

// The body is the lines, tagged and in order, never truncated.
assert.ok(out.startsWith('[INFO] 20 markers shortlisted\n[WARN] rs886288: no primer pair'),
  'the body must be every line, tagged and in order')

// Exactly one footer line, `#`-prefixed, carrying the build/debug provenance.
const footer = out.trimEnd().split('\n').at(-1) ?? ''
assert.ok(footer.startsWith('# originmarker build/debug'), 'the footer identifies the file')
assert.equal(out.trimEnd().split('\n').filter((l) => l.startsWith('#')).length, 1,
  'the debug info is one line, not scattered through the file')
for (const fact of ['release=Build 1.3.3.1', 'job=abc123', 'exported=2026-07-17T20:00:00.000Z',
                    'url=https://ezrakruger.cc/originmarker', 'genome=GRCh38', 'gnomad=gnomad_r4',
                    'ensembl=116', 'agent=TestAgent/1.0']) {
  assert.ok(footer.includes(fact), `the footer must carry ${fact}`)
}
// The footer says WHERE it happened: the deploy (release), the run (job), the instance (url).
assert.ok(footer.includes('release=') && footer.includes('job=') && footer.includes('url='),
  'a log with no release/job/url could have come from anywhere')

// Degrades without a panel: an export taken before a build still identifies the instance
// rather than throwing on missing provenance.
const bare = buildLogText([{ tag: 'INFO', text: 'x' }], {}, CTX)
assert.ok(bare.includes('release=unknown') && bare.includes('job=none'))
assert.ok(!bare.includes('genome='), 'absent provenance drops its keys rather than printing undefined')

console.log('BuildLog.check.ts: all assertions passed')

// Self-check for the build-log wire contract. Run: node src/BuildLog.check.ts
// The component needs a DOM and is not imported here; this covers the normaliser that
// decides what reaches it, which is where the shape assumptions live.
import assert from 'node:assert/strict'
import { asLogLine, LOG_TAGS } from './api.ts'

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

console.log('BuildLog.check.ts: all assertions passed')

// Self-check for the query history store. Run: node src/history.check.ts
// The point of these assertions is the failure modes, not the happy path: localStorage is
// shared with every other version this origin ever served, and it can refuse a write.
import assert from 'node:assert/strict'
import {
  HISTORY_CAP, clearHistory, forgetQuery, noteCount, readHistory, recordQuery,
} from './history.ts'

// Node has no localStorage. This is the whole surface history.ts touches, plus a switch for
// the two failures a browser really produces: a full quota and a refused store.
let mem = new Map<string, string>()
let quotaFull = false
let refused = false
const fake: Storage = {
  get length() { return mem.size },
  key: (i: number) => [...mem.keys()][i] ?? null,
  getItem: (k: string) => mem.get(k) ?? null,
  setItem: (k: string, v: string) => {
    if (quotaFull) throw new DOMException('exceeded the quota', 'QuotaExceededError')
    mem.set(k, v)
  },
  removeItem: (k: string) => { mem.delete(k) },
  clear: () => mem.clear(),
}
const win = { get localStorage(): Storage {
  if (refused) throw new DOMException('access denied', 'SecurityError')
  return fake
} }
;(globalThis as unknown as { window: unknown }).window = win

const KEY = 'originmarker.history'
const reset = () => { mem = new Map(); quotaFull = false; refused = false }
const texts = (es: { text: string }[]) => es.map((e) => e.text)

// --- the basics: most recent first, one entry per distinct text ---------------
reset()
recordQuery('rs334')
recordQuery('rs6025 in Europeans')
assert.deepEqual(texts(readHistory()), ['rs6025 in Europeans', 'rs334'])
recordQuery('rs334')
assert.deepEqual(texts(readHistory()), ['rs334', 'rs6025 in Europeans'], 'a re-typed query is one entry, moved to the top')
// Verbatim: HGVS is case-sensitive, so these are two queries and not one.
recordQuery('NM_000518.5(HBB):c.20A>T')
recordQuery('nm_000518.5(hbb):c.20a>t')
assert.equal(readHistory().length, 4)

reset()
for (let i = 0; i < HISTORY_CAP + 5; i++) recordQuery(`rs${i}`)
assert.equal(readHistory().length, HISTORY_CAP, 'the cap holds')
assert.equal(readHistory()[0].text, `rs${HISTORY_CAP + 4}`, 'the cap drops the oldest, not the newest')

// --- the count: only ever a fact a build reported ----------------------------
reset()
recordQuery('rs334', 'rs334')
recordQuery('never built')
assert.equal(readHistory()[0].count, undefined)
noteCount('rs334', 148)
assert.equal(readHistory().find((e) => e.text === 'rs334')?.count, 148)
// A query that never built stays in the list rather than vanishing for want of a count.
assert.deepEqual(texts(readHistory()), ['never built', 'rs334'])
assert.equal(readHistory().find((e) => e.text === 'never built')?.count, undefined)
// A count for a variant nothing here asked for is dropped, not filed against a neighbour.
noteCount('rs9999999', 12)
assert.deepEqual(readHistory().map((e) => e.count), [undefined, 148])
// Only the entry that named the variant, and only the most recent one to name it.
recordQuery('the sickle cell mutation', 'rs334')
noteCount('rs334', 150)
assert.equal(readHistory().find((e) => e.text === 'the sickle cell mutation')?.count, 150)
assert.equal(readHistory().find((e) => e.text === 'rs334')?.count, 148)

// --- garbage: another version's, or nobody's ---------------------------------
reset()
mem.set(KEY, '{not json at all')
assert.deepEqual(readHistory(), [], 'a parse failure is an empty history, not a throw')
mem.set(KEY, '{"queries":["rs334"]}')
assert.deepEqual(readHistory(), [], 'a schema this version does not know is not read as one it does')
mem.set(KEY, JSON.stringify([
  1, null, 'rs334', { nope: true }, { text: '   ' },
  { text: 'ok', at: 5, count: -3 },
  { text: 'fine', at: 7, variant: 'rs334', count: 9 },
]))
assert.deepEqual(texts(readHistory()), ['ok', 'fine'], 'unreadable entries drop, readable ones survive')
assert.equal(readHistory()[0].count, undefined, 'a count that cannot have come from a panel is not shown as one')
assert.equal(readHistory()[1].count, 9)
// Recording against garbage must not throw, and must not resurrect it.
assert.deepEqual(texts(recordQuery('rs6025')), ['rs6025', 'ok', 'fine'])

// --- a store that refuses -----------------------------------------------------
reset()
recordQuery('rs334')
quotaFull = true
assert.deepEqual(texts(recordQuery('rs113993960')), ['rs113993960', 'rs334'], 'a full quota still returns the list the caller must render')
assert.deepEqual(texts(readHistory()), ['rs334'], 'and the write simply did not land')
assert.doesNotThrow(() => forgetQuery('rs334'))
quotaFull = false

reset()
refused = true
assert.deepEqual(readHistory(), [], 'no storage is an empty history')
assert.doesNotThrow(() => recordQuery('rs334'))
assert.doesNotThrow(() => noteCount('rs334', 5))
assert.doesNotThrow(() => forgetQuery('rs334'))
assert.doesNotThrow(() => clearHistory())

// --- clearing is the feature, not a nicety ------------------------------------
reset()
recordQuery('rs334', 'rs334')
recordQuery('rs6025')
assert.deepEqual(texts(forgetQuery('rs334')), ['rs6025'])
assert.deepEqual(texts(readHistory()), ['rs6025'], 'per-row clearing reaches the store, not just the screen')
assert.deepEqual(clearHistory(), [])
assert.deepEqual(readHistory(), [])
// An empty array left behind is still a record that someone was here.
assert.equal(mem.has(KEY), false, 'clearing all removes the key itself')

// --- what is stored is only what was declared ---------------------------------
reset()
recordQuery('rs334', 'rs334')
noteCount('rs334', 148)
const stored = JSON.parse(mem.get(KEY)!)
assert.deepEqual(Object.keys(stored[0]).sort(), ['at', 'count', 'text', 'variant'])

console.log('history.check.ts: all assertions passed')

// Self-check for the display helpers that make a claim about precision or absence.
// Run: node src/fmt.check.ts
import assert from 'node:assert/strict'
import { sig2, num, sci, orUnknown } from './fmt.ts'

// A fact the panel did not record renders as an admission: absence must be a value this
// helper CAN return, never confusable with a real release number or a formatting dash.
assert.equal(orUnknown(null), 'unknown')
assert.equal(orUnknown(undefined), 'unknown')
assert.equal(orUnknown(''), 'unknown')
assert.equal(orUnknown(116), '116')
assert.equal(orUnknown('deCODE 2019 sex-averaged'), 'deCODE 2019 sex-averaged')
assert.notEqual(orUnknown(null), orUnknown(116))
// 0 is a recorded value, not a gap: never swallow a real datum as "unknown".
assert.equal(orUnknown(0), '0')

// Absence and zero are different claims and must never render the same. This is why these
// helpers exist rather than a bare toFixed at each call site.
assert.equal(sig2(null), '-')
assert.equal(sig2(undefined), '-')
assert.equal(sig2(0), '0')
assert.notEqual(sig2(0), sig2(null))
assert.notEqual(num(0), num(null))
assert.notEqual(sci(0), sci(null))

// Two significant figures, held across the range the bundled map actually spans on a
// ±250 kb window (~5e-5 to ~0.8 cM), where any fixed decimal count is wrong at one end.
assert.equal(sig2(0.75758), '0.76')
assert.equal(sig2(0.0332), '0.033')
assert.equal(sig2(0.00005), '0.00005')
assert.equal(sig2(0.007519), '0.0075')
assert.equal(sig2(1.2), '1.2')

// The precision claim is the point: no rendering may imply digits the map cannot support.
for (const v of [0.75758, 0.0332, 0.00699, 0.007519, 0.00433]) {
  const digits = sig2(v).replace(/^0\.0*/, '').replace('.', '')
  assert.ok(digits.length <= 2, `${v} rendered as ${sig2(v)}: more than two significant figures`)
}

console.log('fmt.check.ts: all assertions passed')

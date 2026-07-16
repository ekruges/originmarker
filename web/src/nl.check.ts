// Self-check for the two decisions the free-text lane hands the UI: whether a model chose
// the variant, and whether the variant that resolved is in the gene the user asked about.
// Run: node src/nl.check.ts
import assert from 'node:assert/strict'
import { geneMismatch, withNlProvenance, type NLResponse, type StructuredQuery } from './api.ts'

const q: StructuredQuery = { variant: 'rs334', build: 'GRCh38' }
const nl = (o: Partial<NLResponse>): NLResponse =>
  ({ query: q, used_llm: false, note: '', ...o })

// --- withNlProvenance ------------------------------------------------------
// A human typing an identifier is not a parse: nothing to record, no keys invented.
assert.deepEqual(withNlProvenance(q, null), q)
assert.equal('nl_text' in withNlProvenance(q, null), false)
assert.equal('nl_model' in withNlProvenance(q, null), false)

// Path B: the model named the variant, and the build request has to say so.
const billed = withNlProvenance(q, nl({
  used_llm: true, model: 'claude-haiku-4-5-20251001', text: 'the sickle cell variant',
}))
assert.equal(billed.nl_model, 'claude-haiku-4-5-20251001')
assert.equal(billed.nl_text, 'the sickle cell variant')

// Path A: free, regex-matched. `used_llm` is the only authority, so a model name riding
// along on a free parse must NOT be recorded as having chosen anything. null, not the
// name, and not undefined: the build request states that no model chose this variant.
const free = withNlProvenance(q, nl({
  used_llm: false, model: 'claude-haiku-4-5-20251001', text: 'panel around rs334',
}))
assert.equal(free.nl_model, null)
assert.equal(free.nl_text, 'panel around rs334')
assert.notEqual(free.nl_model, billed.nl_model) // the two paths can never record the same

// R1 by construction, restated at this seam: folding provenance in adds no coordinate, and
// cannot, since neither shape has anywhere to put one.
for (const k of ['chrom', 'pos', 'pos_grch38', 'ref', 'alt', 'strand']) {
  assert.equal(k in billed, false, `provenance leaked a coordinate field: ${k}`)
}
assert.equal(billed.variant, q.variant) // and never rewrites the identifier itself

// --- geneMismatch ----------------------------------------------------------
// The refusal case: user said ABCC8, the identifier is HBB's.
assert.equal(geneMismatch('HBB', ['ABCC8']), true)
// Agreement, and agreement regardless of how the user cased or spaced it.
assert.equal(geneMismatch('ABCC8', ['ABCC8']), false)
assert.equal(geneMismatch('ABCC8', ['abcc8']), false)
assert.equal(geneMismatch('ABCC8', [' ABCC8 ']), false)
// Named several: matching any one of them is agreement, not a mismatch.
assert.equal(geneMismatch('HBB', ['ABCC8', 'HBB']), false)

// Silence is not disagreement. Naming no gene asks nothing, so there is no warning to
// raise; a record with no gene of its own cannot contradict anyone either. Both directions
// matter: invent a warning here and it fires on the ordinary rsID path forever.
assert.equal(geneMismatch('HBB', []), false)
assert.equal(geneMismatch('HBB', undefined), false)
assert.equal(geneMismatch(null, ['ABCC8']), false)
assert.equal(geneMismatch(null, []), false)
// Empty/whitespace symbols are absence, not a symbol that matches nothing.
assert.equal(geneMismatch('HBB', ['']), false)
assert.equal(geneMismatch('HBB', ['   ']), false)
assert.equal(geneMismatch('', ['ABCC8']), false)

console.log('nl.check.ts: all assertions passed')

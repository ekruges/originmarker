// Self-check for the Progenitor docs' structure, on the same terms as SyngamyDocs.check.ts.
//
// Section numbers are derived from SECTIONS order, so nothing can drift between a heading and
// the nav. What remains checkable is that every nav entry has a section and every section a nav
// entry, that no cross-reference points at an id that does not exist, and that every citation
// marker names a real citation.
//
// One check here has no counterpart on the other pages. This page quotes thresholds that decide
// whether a reconstructed reference is trusted, and a page asserting a floor the code does not
// enforce is worse than a page that says nothing: it is a claim about the tool's behaviour that
// a reader has no way to check. So the constants are imported and compared, not retyped.
// Run: node src/ProgenitorDocs.check.ts
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { CITATIONS } from './citations.ts'
import {
  MIN_ASCERTAINMENT, MIN_PRODUCTS, SAME_PARENT_MAX, DIFFERENT_PARENT_MIN,
} from './inferredReference.ts'
import { pct } from './parentage.ts'
import { PROGENITOR_EXAMPLES } from './progenitorExamples.ts'

const src = readFileSync(new URL('./ProgenitorDocs.tsx', import.meta.url), 'utf8')

const nav = [...src.matchAll(/\{ id: '([\w-]+)', label: (?:'[^']+'|"[^"]+") \}/g)].map((m) => m[1])
const heads = [...src.matchAll(/<Section id="([\w-]+)" title="/g)].map((m) => m[1])

assert.ok(nav.length > 0 && heads.length > 0, 'parsed nothing: the regexes have drifted')
assert.deepEqual(
  [...heads].sort(), [...nav].sort(),
  'every nav entry needs exactly one Section, and every Section an entry',
)
assert.ok(nav.length > 15, `a smaller page than intended: ${nav.length} sections`)
assert.equal(new Set(nav).size, nav.length, 'a duplicated section id would break its anchor')

// A heading that typed its own number would print it twice: "8 · 8 · Stage 3: how deep to build".
for (const [, id, title] of src.matchAll(/<Section id="([\w-]+)" title="([^"]*)"/g)) {
  assert.ok(!/^\d+\s*·/.test(title), `${id}'s title still hard-codes its number: ${title}`)
}

for (const [, id] of src.matchAll(/<SecRef id="([\w-]+)" \/>/g)) {
  assert.ok(nav.includes(id), `a cross-reference points at '${id}', which is not a section`)
}

const cited = [...src.matchAll(/<Ref id="([\w-]+)" \/>/g)].map((m) => m[1])
for (const id of cited) {
  assert.ok(id in CITATIONS, `a citation marker names '${id}', which is not in CITATIONS`)
}
assert.ok(cited.includes('zuccaro_2020'),
  'the series every measurement on this page came from is no longer cited')

// The sections a reader needs before trusting a reconstructed reference. Losing one silently
// would leave the page describing a method without its boundary.
for (const id of ['method', 'membership', 'depth', 'contamination', 'verify', 'floor',
  'refusals', 'validation', 'limits', 'privacy']) {
  assert.ok(nav.includes(id), `the '${id}' section is gone, and it is not an optional one`)
}

// Thresholds, compared against the constants rather than against the source text. Asserting that
// `pct(MIN_ASCERTAINMENT, 0)` appears only proves the page interpolates SOMETHING; it does not
// prove the reader sees the number the code enforces. So render each one the way the page does
// and require the rendered string to be present, which fails if either side moves.
assert.ok(src.includes('{MIN_PRODUCTS}') || src.includes(`${MIN_PRODUCTS} products`),
  'the product floor is not stated')
for (const [name, rendered] of [
  ['MIN_ASCERTAINMENT', pct(MIN_ASCERTAINMENT, 0)],
  ['SAME_PARENT_MAX', pct(SAME_PARENT_MAX, 1)],
  ['DIFFERENT_PARENT_MIN', pct(DIFFERENT_PARENT_MIN, 1)],
] as const) {
  assert.ok(src.includes(`pct(${name}, `),
    `${name} must be rendered from the constant, not typed`)
  assert.ok(!new RegExp(`[^\\w.]${rendered.replace(/[.%]/g, '\\$&')}`).test(src),
    `${name} renders as ${rendered}, and that literal is typed into the page as well: `
    + 'one of the two will drift')
}
// Newlines excluded on purpose: a string literal cannot contain one, and without that the
// character class runs past the closing quote and matches entities in the JSX below it.
assert.ok(!/'[^'\n]*&[a-z]+;[^'\n]*'/.test(src),
  'an HTML entity inside a plain string renders literally: use the character')
assert.equal(MIN_PRODUCTS, PROGENITOR_EXAMPLES.length,
  'the bundled examples must reach the floor, or the positive case cannot complete')

// Every documentation page links to the other two.
for (const other of ['#/syngamy-docs', '#/docs']) {
  assert.ok(src.includes(`href: '${other}'`), `this page must link to ${other}`)
}
assert.ok(/const PREFIX = 'progenitor-docs'/.test(src), 'the anchor route has changed')
assert.ok(src.includes('docsHelpers(PREFIX, SECTIONS)'),
  'section numbering must come from the shell, bound to this page\'s own route')

console.log(`ProgenitorDocs.check OK (${nav.length} sections, ${cited.length} citation markers, `
  + 'nav == headings, cross-refs and citations resolve, thresholds match the code)')

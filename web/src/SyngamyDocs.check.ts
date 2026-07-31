// Self-check for the Syngamy docs' structure, on the same terms as DocsPage.check.ts.
//
// Section numbers are derived from SECTIONS order, so nothing can drift between a heading and
// the nav. What remains checkable is that every nav entry has a section and every section a nav
// entry, that no cross-reference points at an id that does not exist, and that every citation
// marker names a real citation. That last one matters more here than on the panel docs: a
// reference list is the part a reviewer checks first, and a [n] pointing at nothing renders as
// nothing at all rather than as an error.
// Run: node src/SyngamyDocs.check.ts
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { CITATIONS } from './citations.ts'

const src = readFileSync(new URL('./SyngamyDocs.tsx', import.meta.url), 'utf8')

const nav = [...src.matchAll(/\{ id: '([\w-]+)', label: (?:'[^']+'|"[^"]+") \}/g)].map((m) => m[1])
const heads = [...src.matchAll(/<Section id="([\w-]+)" title="/g)].map((m) => m[1])

assert.ok(nav.length > 0 && heads.length > 0, 'parsed nothing: the regexes have drifted')
assert.deepEqual(
  [...heads].sort(), [...nav].sort(),
  'every nav entry needs exactly one Section, and every Section an entry',
)
assert.ok(nav.length > 20, `a smaller page than intended: ${nav.length} sections`)
assert.equal(new Set(nav).size, nav.length, 'a duplicated section id would break its anchor')

// A heading that typed its own number would print it twice: "5 · 5 · The three measurements".
for (const [, id, title] of src.matchAll(/<Section id="([\w-]+)" title="([^"]*)"/g)) {
  assert.ok(!/^\d+\s*·/.test(title), `${id}'s title still hard-codes its number: ${title}`)
}

for (const [, id] of src.matchAll(/<SecRef id="([\w-]+)" \/>/g)) {
  assert.ok(nav.includes(id), `a cross-reference points at '${id}', which is not a section`)
}

// Every inline [n] must resolve. Ref renders null for an unknown id, so a typo would silently
// drop the marker and leave a sentence citing nothing.
const cited = [...src.matchAll(/<Ref id="([\w-]+)" \/>/g)].map((m) => m[1])
assert.ok(cited.length >= 5, `too few citations for a page making this many claims: ${cited.length}`)
for (const id of cited) {
  assert.ok(id in CITATIONS, `a citation marker names '${id}', which is not in CITATIONS`)
}

// The five works this page's argument rests on. Losing one silently would leave the claim
// standing without its source, which is the specific failure this project cares about.
for (const id of ['ma_2017', 'egli_2018', 'zuccaro_2020', 'natesan_2014', 'ado']) {
  assert.ok(cited.includes(id), `${id} is no longer cited anywhere on the page`)
}

// The route the header links to, and the one the panel docs point across at.
assert.ok(src.includes("`#/syngamy-docs/${id}`"), 'the anchor route has changed')

// Both directions of the cross-link between the two documentation pages.
const panel = readFileSync(new URL('./DocsPage.tsx', import.meta.url), 'utf8')
assert.ok(src.includes('href="#/docs"'), 'this page must link back to the panel docs')
assert.ok(panel.includes('#/syngamy-docs'), 'the panel docs must link here')

console.log(`SyngamyDocs.check OK (${nav.length} sections, ${cited.length} citation markers, `
  + 'nav == headings, cross-refs and citations resolve)')

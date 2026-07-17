// Self-check for the docs' structure.
//
// Section numbers used to be typed into each heading by hand, and this check existed to
// catch them drifting from the nav. They are derived from SECTIONS order now, the way
// citation numbers always were, so there is nothing left to drift: what remains is that
// every nav entry has a section, every section is in the nav, and no cross-reference points
// at an id that does not exist.
// Run: node src/DocsPage.check.ts
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'

const src = readFileSync(new URL('./DocsPage.tsx', import.meta.url), 'utf8')

const nav = [...src.matchAll(/\{ id: '([\w-]+)', label: '([^']+)' \}/g)].map((m) => m[1])
const heads = [...src.matchAll(/<Section id="([\w-]+)" title="/g)].map((m) => m[1])

assert.ok(nav.length > 0 && heads.length > 0, 'parsed nothing: the regexes have drifted')
assert.deepEqual(
  [...heads].sort(), [...nav].sort(),
  'every nav entry needs exactly one Section, and every Section an entry',
)

// A heading that still typed its own number would print it twice: "5 · 5 · The star".
for (const [, id, title] of src.matchAll(/<Section id="([\w-]+)" title="([^"]*)"/g)) {
  assert.ok(!/^\d+\s*·/.test(title), `${id}'s title still hard-codes its number: ${title}`)
}

// Cross-references render whatever number their target holds, so they can no longer say the
// wrong one. They can still name a section that does not exist, which renders as "section 0".
for (const [, id] of src.matchAll(/<SecRef id="([\w-]+)" \/>/g)) {
  assert.ok(nav.includes(id), `a cross-reference points at '${id}', which is not a section`)
}

// Every primer warning on the panel, and the primer form, link here by this exact route
// (primers.py's PRIMER_DOCS, mirrored in api.ts). A rename that missed it would leave every
// note on the page pointing at nothing.
assert.ok(nav.includes('primers'), "the primer notes link to '#/docs/primers'")

console.log(`DocsPage.check OK (${nav.length} sections, nav == headings, cross-refs resolve)`)

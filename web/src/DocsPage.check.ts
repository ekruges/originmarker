// Self-check for the docs numbering, which is coupled in two places that cannot see each
// other: the nav numbers a section by its position in SECTIONS, while the heading carries
// the number as text. Drift is silent and the nav lies.
// Run: node src/DocsPage.check.ts
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'

const src = readFileSync(new URL('./DocsPage.tsx', import.meta.url), 'utf8')

const nav = [...src.matchAll(/\{ id: '([\w-]+)', label: '([^']+)' \}/g)].map((m) => m[1])
const heads = new Map(
  [...src.matchAll(/<Section id="([\w-]+)" title="(\d+) · /g)].map((m) => [m[1], Number(m[2])]),
)

assert.ok(nav.length > 0 && heads.size > 0, 'parsed nothing: the regexes have drifted')
assert.equal(heads.size, nav.length, 'every nav entry needs exactly one numbered Section')

// The nav numbers by array order, so position i is section i+1 and the heading must agree.
nav.forEach((id, i) => {
  assert.equal(heads.get(id), i + 1, `nav shows ${id} as ${i + 1}, its heading disagrees`)
})

// Cross-references are written as prose ("section 13"), so they cannot follow a renumber
// on their own. Each must name the number its target's heading actually carries.
for (const [, id, n] of src.matchAll(/docHref\('([\w-]+)'\)}>section (\d+)</g)) {
  assert.equal(heads.get(id), Number(n), `a link says section ${n} but ${id} is ${heads.get(id)}`)
}

console.log(`DocsPage.check OK (${nav.length} sections, nav == headings == cross-refs)`)

// Self-check for the disputed-position alarm. Run: node src/PanelTable.check.ts
// Node strips types but not JSX, so PanelTable.tsx cannot be imported directly and vite
// loads it instead. Keep this file plain .ts (createElement, no JSX) or node cannot run it.
import assert from 'node:assert/strict'
import { createElement as h } from 'react'
import { renderToStaticMarkup } from 'react-dom/server'
import { MantineProvider } from '@mantine/core'
import { createServer } from 'vite'
import type { Marker, PanelResult } from './api.ts'

const vite = await createServer({ server: { middlewareMode: true }, appType: 'custom' })
const { PanelTable } = await vite.ssrLoadModule('/src/PanelTable.tsx')
await vite.close()

// rank.check.ts covers the posMismatch predicate; this covers the banner it feeds, which
// the predicate check cannot see. Markers are hand-built because the golden ABCC8 fixture
// has no mismatches and never renders the banner.

const mk = (o: Partial<Marker> & { rsid: string }): Marker => ({
  variant_id: o.rsid, chrom: '11', pos: 17_397_055 + (o.dist ?? 0), ref: 'A', alt: 'G',
  af: 0.5, maf: 0.5, het: 0.5, het_max_pop: 0.5, dist: 0, side: 'lower coord',
  tier: 'A_core(<2kb)', per_pop_maf: {}, ensembl_pos_check: null, cm: null,
  recomb_fraction: null, hotspot_between: null, map_approx: null, ...o,
})

const panel = (ms: Marker[]): PanelResult => ({
  variant: {} as PanelResult['variant'],
  rarity: {} as PanelResult['rarity'],
  candidates: ms,
  recommended: ms,
  coverage: { lower_count: 0, higher_count: 0, lower_core_near: 0, higher_core_near: 0, flags: [] },
  params: {} as PanelResult['params'],
  provenance: {} as PanelResult['provenance'],
})

const render = (ms: Marker[]) =>
  renderToStaticMarkup(
    h(MantineProvider, null, h(PanelTable, { result: panel(ms), ancestry: null })),
  )

const text = (html: string) =>
  html.replace(/<[^>]+>/g, ' ').replace(/&#x27;/g, "'").replace(/\s+/g, ' ')

/** Text of the alarm ONLY, sliced off above the table. Must stay scoped: the table's coord
 *  column prints the same position, so a whole-document match would pass on a banner that
 *  dropped it. */
const banner = (html: string) => {
  const i = html.indexOf('role="alert"')
  return i < 0 ? '' : text(html.slice(i, html.indexOf('<table', i)))
}

// --- 1. A disputed marker raises the alarm, with both positions and the gap. --------
// Every figure on the row is computed from the gnomAD position, so the reader needs both
// numbers and the gap to judge whether the dispute matters.
const disputed = render([
  mk({ rsid: 'rs757110', pos: 17_396_930, dist: -125, ensembl_pos_check: 'MISMATCH:17401130' }),
  mk({ rsid: 'rs739689', pos: 17_395_957, dist: -1098, ensembl_pos_check: 'ok' }),
])
assert.match(disputed, /role="alert"/)
const dt = banner(disputed)
assert.match(dt, /disputed position/)
assert.match(dt, /rs757110/)
assert.match(dt, /gnomAD 17,396,930/)   // gnomAD's claim, named as gnomAD's
assert.match(dt, /Ensembl 17,401,130/)  // Ensembl's claim, not swallowed
assert.match(dt, /4,200 bp apart/)      // the gap, computed not asserted
assert.match(dt, /in the shortlist/)    // it is still being recommended, and says so
// The row itself is tinted: the banner scrolls away, the row does not.
assert.equal((disputed.match(/background:rgba\(224,\s*49,\s*49/g) ?? []).length, 1)

// --- 2. No alarm when nothing is disputed. -----------------------------------------
// 'ok' is agreement and null is silence; neither is a dispute. A banner that cries wolf on
// unchecked rows is the same failure as no banner.
const quiet = render([
  mk({ rsid: 'rs739689', ensembl_pos_check: 'ok' }),
  mk({ rsid: 'rs_unchecked', ensembl_pos_check: null }),
])
assert.doesNotMatch(quiet, /role="alert"/)
assert.doesNotMatch(text(quiet), /disputed/)
assert.equal((quiet.match(/background:rgba\(224,\s*49,\s*49/g) ?? []).length, 0)

// --- 3. The alarm is not silenceable by pagination. --------------------------------
// The banner must derive from result.candidates, not the current page, filter or scope
// toggle, or a reader hunting for a marker hides it by looking. 60 markers at PER_PAGE=50,
// the disputed one sorted last by het_max_pop, so it lands on page 2.
const many = [
  ...Array.from({ length: 59 }, (_, i) =>
    mk({ rsid: `rs_ok_${i}`, dist: i + 1, het_max_pop: 0.5, ensembl_pos_check: 'ok' })),
  mk({
    rsid: 'rs_offpage', pos: 17_396_930, dist: 900, het_max_pop: 0.01,
    ensembl_pos_check: 'MISMATCH:17401130',
  }),
]
const paged = render(many)
// Precondition: it really is off page 1, or the assertion below proves nothing.
const table = paged.slice(paged.indexOf('<table'))
assert.doesNotMatch(text(table), /rs_offpage/, 'setup: expected rs_offpage off page 1')
assert.match(paged, /role="alert"/)
assert.match(banner(paged), /rs_offpage/)     // named though its row is not rendered
assert.match(banner(paged), /4,200 bp apart/)

console.log('PanelTable.check.ts: all assertions passed')

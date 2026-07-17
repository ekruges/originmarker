// Self-check for the primer UI. Run: node src/PrimerOptions.check.ts
// Node strips types but not JSX, so the components cannot be imported directly and vite
// loads them instead. Keep this file plain .ts (createElement, no JSX) or node cannot run it.
import assert from 'node:assert/strict'
import { createElement as h } from 'react'
import { renderToStaticMarkup } from 'react-dom/server'
import { MantineProvider, Table } from '@mantine/core'
import { createServer } from 'vite'
import {
  PRIMER_SETTING_KEYS, PRIMER_SIZE_CAP, pcrDanger, pcrUnchecked, primerBuild,
  type Marker, type PanelResult, type PrimerBuild, type PrimerResult, type PrimerSettings,
  type PrimerWarning,
} from './api.ts'

const vite = await createServer({ server: { middlewareMode: true }, appType: 'custom' })
const { PanelTable, PrimerDetail } = await vite.ssrLoadModule('/src/PanelTable.tsx')
const { SearchPanel } = await vite.ssrLoadModule('/src/SearchPanel.tsx')
const { PrimerFields, PRIMER_FIELDS } = await vite.ssrLoadModule('/src/PrimerOptions.tsx')
await vite.close()

// primers.PrimerSettings' own defaults, as a panel records them.
const SETTINGS: PrimerSettings = {
  opt_tm: 69, min_tm: 67, max_tm: 71, max_pair_diff_tm: 5,
  min_size: 20, opt_size: 26, max_size: 35,
  min_gc: 40, max_gc: 60, gc_clamp: 1,
  salt_monovalent: 50, salt_divalent: 1.5, dntp_conc: 0.6, dna_conc: 50,
  min_product: 250, max_product: 600,
  max_poly_x: 4, max_self_any_th: 47, max_self_end_th: 47, max_hairpin_th: 47,
  max_ns_accepted: 0, tm_formula: 1, salt_corrections: 1,
  mask_maf: 0.01, target_pad: 50,
}
const BUILD: PrimerBuild = { scope: 'starred', settings: SETTINGS }
const PARAMS = { primer_scope: 'starred', primer_settings: SETTINGS }
// primers.NOT_CHECKED_WARNING is the engine's, not this file's: a sentinel proves the UI
// renders whatever it is handed rather than a copy it keeps. Both lengths are sentinelled
// apart, because the table must take `short` and only `short`: `long` is a paragraph and is
// the same paragraph on every row.
const WARN: PrimerWarning = {
  code: 'not_checked',
  short: 'SENTINEL-SHORT: the engine words this and the table prints it.',
  long: 'SENTINEL-LONG: the whole of it, which belongs in the exports and not in a table.',
  docs: '#/docs/primers',
}

const SEQ_L = 'GAGGGGGTCTCACTGTGTTGACCTA'
const SEQ_R = 'TTGACCATGGCTTAGGCATTCAGGA'
const PAIR: PrimerResult = {
  fwd: { seq: SEQ_L, pos: 17_396_800, idx: 100, length: SEQ_L.length, tm: 69.1, gc: 56 },
  rev: { seq: SEQ_R, pos: 17_397_212, idx: 512, length: SEQ_R.length, tm: 68.8, gc: 48 },
  product_size: 412,
  mask_note: {
    code: 'mask',
    short: 'SENTINEL-MASKSHORT: primer sites clear of 4 variants at MAF >= 1.000%.',
    long: 'SENTINEL-MASKLONG: the whole paragraph, which is the same on every row.',
    docs: '#/docs/primers',
  },
  warnings: [WARN],
  insilico_pcr: 'not_checked',
}

const mk = (o: Partial<Marker> & { rsid: string }): Marker => ({
  variant_id: o.rsid, chrom: '11', pos: 17_397_055 + (o.dist ?? 0), ref: 'A', alt: 'G',
  af: 0.5, maf: 0.5, het: 0.5, het_max_pop: 0.5, dist: 1, side: 'lower coord',
  tier: 'A_core(<2kb)', per_pop_maf: {}, ensembl_pos_check: null, cm: null,
  recomb_fraction: null, hotspot_between: null, map_approx: null, ...o,
})

const panel = (
  primer?: PrimerResult,
  params: Record<string, unknown> = PARAMS,
): PanelResult => {
  const ms = [mk({ rsid: 'rs757110', primer })]
  return {
    variant: {} as PanelResult['variant'],
    rarity: {} as PanelResult['rarity'],
    candidates: ms,
    recommended: ms,
    coverage: { lower_count: 0, higher_count: 0, lower_core_near: 0, higher_core_near: 0, flags: [] },
    params,
    provenance: {} as PanelResult['provenance'],
  }
}

const render = (r: PanelResult) =>
  renderToStaticMarkup(h(MantineProvider, null, h(PanelTable, { result: r, ancestry: null })))

/** The primer box forced open. SSR cannot click, and the box is collapsed by default, so
 *  what is inside it needs rendering directly. Wrapped in a Mantine Table because it is a
 *  row and Table.Tr reads that context. */
const renderOpen = (d: PrimerResult) =>
  renderToStaticMarkup(h(MantineProvider, null,
    h(Table, null, h(Table.Tbody, null, h(PrimerDetail, { d, defaultOpen: true })))))

const text = (html: string) => html.replace(/<[^>]+>/g, ' ').replace(/&#x27;/g, "'").replace(/\s+/g, ' ')

// --- 1. The three-state read, and an unknown state is not a pass. -------------------
// `insilico_pcr` is a STATE CODE off the wire, never prose: app/jobs.py writes ispcr's state
// there and welds the words on as the warning carrying that same code. This check used to
// feed it 'DANGEROUS: 3 products' and assert the prose came back, which the server never
// sent: it agreed with itself while the page rendered the bare word "danger" at the reader.
// Only the two known-safe tokens may render clean. A state this UI cannot read is a state
// nobody checked, and rendering it as clean is the one failure that turns an unverified
// primer into a verified-looking one.
assert.equal(pcrDanger({ insilico_pcr: 'not_checked' }), null)
assert.equal(pcrDanger({ insilico_pcr: 'one_product' }), null)
assert.equal(pcrDanger({}), null)                          // absent is silence, not danger
assert.equal(pcrDanger({ insilico_pcr: null }), null)
assert.ok(pcrDanger({ insilico_pcr: 'ok' }), 'an unreadable verdict must read as dangerous')
assert.ok(pcrDanger({ insilico_pcr: 'PASSED' }), 'an unreadable verdict must read as dangerous')

// 'unknown' is ispcr's FOURTH state and the one this file used to have no name for: asked,
// and the answer could not be read. It is neither of the other two things. UCSC never
// answered, so "do not order this pair without redesigning it" would be a finding invented
// here; and nothing was ruled out, so it is not clean either. Both predicates are asserted,
// because narrowing only pcrDangerous drops it through the gap and renders it green.
assert.equal(pcrDanger({ insilico_pcr: 'unknown' }), null,
  'a timeout is not UCSC contradicting the pair')
assert.equal(pcrUnchecked({ insilico_pcr: 'unknown' }), true,
  'a timeout is not a pass either: only a parsed single product is')

// The verdict's words are found by code, not by position: the caveat rides alongside it.
const DANGER_NOTE: PrimerWarning = {
  code: 'danger', short: 'SENTINEL: UCSC finds 3 products in hg38, not one.',
  long: 'the long form', docs: '#/docs/primers',
}
const CAVEAT_NOTE: PrimerWarning = {
  code: 'ispcr_caveat', short: 'in silico only', long: 'the long form', docs: '#/docs/primers',
}
assert.equal(
  pcrDanger({ insilico_pcr: 'danger', warnings: [CAVEAT_NOTE, DANGER_NOTE] }), DANGER_NOTE,
  'the verdict is the note whose code is the state, whatever order they arrive in',
)

// A dangerous state whose words did not arrive still has to shout. A red row that says
// nothing is read as a clean one.
assert.match(pcrDanger({ insilico_pcr: 'danger', warnings: [] })!.short, /unverified/)
assert.match(pcrDanger({ insilico_pcr: 'danger' })!.short, /unverified/)

// Unchecked and danger are separate questions: absent, null and the token all mean the
// genome was never consulted, and none of them is a pass.
assert.equal(pcrUnchecked({}), true)
assert.equal(pcrUnchecked({ insilico_pcr: null }), true)
assert.equal(pcrUnchecked({ insilico_pcr: 'not_checked' }), true)
assert.equal(pcrUnchecked({ insilico_pcr: 'one_product' }), false)

// --- 2. Settings are read off the panel, and a gap is not a setting. ----------------
// The form seeds from these numbers and rebuilds on them, so a params block missing a knob
// must yield null rather than a form showing a value nobody chose.
assert.deepEqual(primerBuild(panel(undefined, PARAMS)), BUILD)
assert.equal(primerBuild(panel(undefined, {})), null)
assert.equal(primerBuild(panel(undefined, { primer_scope: 'starred' })), null)
assert.equal(primerBuild(panel(undefined, { ...PARAMS, primer_scope: 'everything' })), null)
for (const k of PRIMER_SETTING_KEYS) {
  const short = { ...SETTINGS } as Record<string, unknown>
  delete short[k]
  assert.equal(primerBuild(panel(undefined, { ...PARAMS, primer_settings: short })), null,
    `${k} missing must void the set`)
  assert.equal(
    primerBuild(panel(undefined, { ...PARAMS, primer_settings: { ...SETTINGS, [k]: String(SETTINGS[k]) } })),
    null, `${k} as a string must void the set`,
  )
}

// --- 3. No primer3 on the server: nothing about primers renders. --------------------
// Not an empty section and not a form over an absent feature. Matched on primer-specific
// wording: the x-check column says "not checked" about a position, which is a different
// claim about a different thing and must not be what makes this pass.
const bare = render(panel(undefined, {}))
assert.doesNotMatch(text(bare), /Primers ·|Design primers for|primer pair|Tm |mask MAF/i)

// --- 4. A finding is never behind a click. ------------------------------------------
// The box is collapsed by default, because four lines of detail per row is a table nobody
// reads. What collapsing must never do is hide a FINDING: a dangerous pair and a failed
// design open themselves, so the two states worth interrupting the reader for are in the
// markup before anyone touches anything.
//
// jobs.py writes what a verified pair carries: the verdict as the note whose code IS the
// state, and the caveat beside it. Built that way here rather than hand-waved, so this
// agrees with the server or fails.
const verdict = (code: string, short: string): PrimerWarning[] => [
  { code, short, long: 'the long form, which belongs in the exports', docs: '#/docs/primers' },
  { code: 'ispcr_caveat', short: 'SENTINEL-CAVEAT: in silico only, not a wet-lab validation.',
    long: 'the long form', docs: '#/docs/primers' },
]

const danger = render(panel({
  ...PAIR, insilico_pcr: 'danger',
  warnings: verdict('danger', 'SENTINEL: UCSC finds 3 products in hg38, not one.'),
}))
assert.match(danger, /role="alert"/)
assert.match(text(danger), /3 products in hg38/)        // the engine's words, not a summary
assert.match(text(danger), new RegExp(SEQ_L), 'a dangerous pair must still show its primer')
assert.match(text(danger), new RegExp(SEQ_R), 'a dangerous pair must still show its primer')
assert.match(text(danger), /Hide primer design/, 'a dangerous pair opens itself')
// Named above the table too, or a reader hides the finding by paging past its row.
assert.match(text(danger).slice(0, danger.indexOf('<table')), /rs757110/)

// A verdict UCSC could not give is reported as a verdict UCSC could not give: no DANGER
// badge, no "In-silico PCR contradicts" banner, no green pass. The page said UCSC
// contradicted N pairs while the PDF from the same job id said NOT VERIFIED for all N: one
// job, two documents, two different instructions, and the screen's was invented.
const timedOut = render(panel({
  ...PAIR, insilico_pcr: 'unknown',
  warnings: verdict('unknown', 'SENTINEL: Still unverified: UCSC could not be reached.'),
}))
assert.doesNotMatch(text(timedOut), /DANGER/, 'a timeout is not a finding against the pair')
assert.doesNotMatch(text(timedOut), /contradicts/, 'UCSC contradicted nothing: it never answered')
assert.doesNotMatch(text(timedOut), /one product/, 'nor is a timeout a pass')
assert.doesNotMatch(timedOut, /role="alert"/, 'and it is not an alarm')
// Opened, it still carries ispcr's own sentence rather than silence.
assert.match(text(renderOpen({
  ...PAIR, insilico_pcr: 'unknown',
  warnings: verdict('unknown', 'SENTINEL: Still unverified: UCSC could not be reached.'),
})), /Still unverified/)

// A token this UI has never heard of takes the loud path, and still shows the pair.
const unknown = render(panel({ ...PAIR, insilico_pcr: 'looks fine' }))
assert.match(unknown, /role="alert"/)
assert.match(text(unknown), new RegExp(SEQ_L))

// --- 5. Not-checked is its own state, and it is not verified-clean. -----------------
// The benign states start shut. The line is what remains, and it names the state on itself:
// nothing here reads as a verdict, because none was reached.
const unchecked = render(panel(PAIR))
assert.match(text(unchecked), /Open primer design/)
assert.doesNotMatch(text(unchecked), new RegExp(SEQ_L), 'a benign pair starts collapsed')
assert.doesNotMatch(unchecked, /role="alert"/, 'unchecked is the normal state, not an alarm')
assert.doesNotMatch(text(unchecked), /DANGER|one product/, 'nothing was checked: say nothing')

// Opened, it is the engine's short form and a route to the rest. Never `long`: that is a
// paragraph, it is the same paragraph on every row, and the docs and the PDF both hold it.
const openUnchecked = renderOpen(PAIR)
assert.match(text(openUnchecked), /WARNING/)
assert.match(text(openUnchecked), new RegExp(WARN.short.slice(0, 30)), 'the engine words it')
assert.doesNotMatch(text(openUnchecked), /SENTINEL-LONG/, 'the table takes the short form')
assert.match(openUnchecked, /href="#\/docs\/primers"/, 'and links to the rest of it')
assert.match(text(openUnchecked), new RegExp(SEQ_L))
// The mask note is the longest line under a pair and the same paragraph on every row. Its
// counts are per-marker and stay; the paragraph after them is the docs' to carry.
assert.match(text(openUnchecked), /SENTINEL-MASKSHORT/)
assert.doesNotMatch(text(openUnchecked), /SENTINEL-MASKLONG/, 'the mask note has two lengths too')

// A note whose short form did not survive the wire still says something. A WARNING badge
// beside an empty line reads as a pair with nothing wrong with it, which is the one thing
// this module refuses. An older server's shape rendered exactly that.
const shapeless = renderOpen({ ...PAIR, warnings: [{ code: 'x', short: '', long: '', docs: '' }] })
assert.match(text(shapeless), /could not read/)

// A pass says so on the line, and the caveat rides with it: never stronger than what was
// done. The two states must not borrow each other's wording.
const clean = renderOpen({
  ...PAIR, insilico_pcr: 'one_product',
  warnings: verdict('one_product', 'SENTINEL-PASS: UCSC finds one product, as designed.'),
})
assert.match(text(clean), /one product/)
assert.match(text(clean), /SENTINEL-CAVEAT/, 'a pass still carries the caveat')
assert.doesNotMatch(text(clean), /NOT VERIFIED|DANGER/)

// The engine promises a pair never arrives without warnings. If one does, the promise is
// broken, and the pair must not render as the one clean primer on the page.
const silent = renderOpen({ ...PAIR, warnings: [] })
assert.match(text(silent), /which the design does not do/)
assert.match(text(silent), new RegExp(SEQ_L), 'even a contract breach still shows the primer')

// --- 6. A design that failed is a result, and says why. -----------------------------
// Silence here would read as a marker nobody tried, which is a different fact. It opens
// itself for the same reason a danger does: it is an answer, not detail.
const failed = render(panel({
  error: 'No 20-35 bp primer reaches 67 C in this window at 40-60% GC; the window is 31.2% GC.',
}))
assert.match(failed, /role="alert"/)
assert.match(text(failed), /No primer pair for this marker/)
assert.match(text(failed), /window is 31.2% GC/)        // the engine's words, and its knob
assert.match(text(failed), /Hide primer design/, 'a failed design opens itself')
// No pair, so nothing may claim anything about one.
assert.doesNotMatch(text(failed), /NOT VERIFIED|one product/)

// --- 7. The chip states the settings the panel was built under. ---------------------
const chipped = render(panel(PAIR))
assert.match(text(chipped), /Primers · starred · Tm 69/)

// --- 8. The form both mounts share offers every knob, seeded from the given numbers. -
// PrimerFields directly, because it is the unit: the chip and manual input render this
// same component, which is what stops the two paths offering different knobs.
const fields = (readOnly?: boolean) =>
  renderToStaticMarkup(h(MantineProvider, null, h(PrimerFields, {
    value: BUILD, onChange: () => {}, readOnly,
  })))

const form = fields()
assert.match(text(form), /Design primers for/)
// Every knob reaches the form. A knob on the wire with no field is a setting the panel
// reports and the user cannot touch, and PRIMER_NUM_KEYS is the list both sides agree on.
for (const k of Object.keys(PRIMER_FIELDS) as (keyof typeof PRIMER_FIELDS)[]) {
  assert.match(form, new RegExp(`value="${String(SETTINGS[k]).replace('.', '\\.')}"`),
    `${k} must render its stated value`)
}
// The cap is the primer3 Tm model's, not a preference: no length field may offer past it,
// because primer3 accepts a longer oligo and returns a Tm the model does not define there.
// Asserted on the bounds rather than the markup: Mantine clamps in JS and renders no max.
for (const k of ['min_size', 'opt_size', 'max_size'] as const) {
  assert.equal(PRIMER_FIELDS[k].max, PRIMER_SIZE_CAP, `${k} must stop at the Tm model's limit`)
}
// The salt must be on the form, not implied: the same oligo reads a different Tm under
// primer3's two default sets, so a Tm without its conditions is not reproducible.
assert.match(text(form), /divalent mM/)
assert.match(text(form), /dNTP mM/)

// Read-only is the no-rebuild state, and it must still state the numbers: they are the
// panel's provenance whether or not this page can act on them.
const ro = fields(true)
assert.match(ro, /readonly/i)
assert.match(ro, new RegExp(`value="${SETTINGS.opt_tm}"`))

// Scope 'none' asks for no primers, so the knobs that describe one have nothing to say.
const off = renderToStaticMarkup(h(MantineProvider, null, h(PrimerFields, {
  value: { ...BUILD, scope: 'none' }, onChange: () => {},
})))
assert.match(text(off), /Design primers for/)
assert.doesNotMatch(text(off), /mask MAF floor|min Tm/)

// --- 9. Manual input holds the fields back until the server states its defaults. ----
// Rendered in its default mode, which offers no manual knobs at all: the assertion is that
// primer fields never leak onto the search path, whatever the server sends.
const search = (health: unknown) =>
  renderToStaticMarkup(h(MantineProvider, null, h(SearchPanel, {
    health, busy: false, onResolve: () => {},
  })))
for (const health of [
  { nl_enabled: false, disclaimer: '', layer_b_steps: [] },
  { nl_enabled: false, disclaimer: '', layer_b_steps: [], primer_defaults: SETTINGS },
]) {
  assert.doesNotMatch(text(search(health)), /Design primers for|mask MAF floor/)
}

console.log('PrimerOptions.check.ts: all assertions passed')

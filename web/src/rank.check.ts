// Self-check for the client-side ranking mirror. Run: node src/rank.check.ts
import assert from 'node:assert/strict'
import { posMismatch, type Marker } from './api.ts'
import { twoPQ, ancestryHet, rankHet, shownHet, rankMarkers, minPopHet, applyPreset } from './rank.ts'

// het and het_max_pop must stay DIFFERENT numbers: where they coincide, a check cannot
// tell which one the code under test read. het_max_pop is a max over populations, so it
// is >= het.
const mk = (o: Partial<Marker> & { rsid: string }): Marker => ({
  variant_id: o.rsid, chrom: '11', pos: 17_397_055 + (o.dist ?? 0), ref: 'A', alt: 'G',
  af: 0.5, maf: 0.5, het: 0.32, het_max_pop: 0.48, dist: 0, side: 'lower coord',
  tier: 'A_core(<2kb)', per_pop_maf: {}, ensembl_pos_check: null, cm: null,
  recomb_fraction: null, hotspot_between: null, map_approx: null, ...o,
})

// 2pq is symmetric about p=0.5 and peaks there.
const close = (a: number, b: number) => assert.ok(Math.abs(a - b) < 1e-12, `${a} != ${b}`)
assert.equal(twoPQ(0.5), 0.5)
close(twoPQ(0.1), twoPQ(0.9)) // symmetric, so MAF vs AF both work
assert.ok(twoPQ(0.5) > twoPQ(0.3))

// ancestryHet reads per_pop_maf; null when the pop is absent or no ancestry selected.
const m = mk({ rsid: 'rs1', per_pop_maf: { AFR: 0.4, NFE: 0.1 } })
assert.equal(ancestryHet(m, 'AFR'), twoPQ(0.4))
assert.equal(ancestryHet(m, 'EAS'), null)
assert.equal(ancestryHet(m, null), null)

// rankHet falls back to the GLOBAL 2pq like panelbuilder._rank_key, never to het_max_pop:
// the PDF stamps provenance["ranking_key"] = "global 2pq (het)", and a drifted mirror
// orders the same panel differently on screen and in the export.
// het, het_max_pop and the AFR figure are three DIFFERENT numbers here on purpose.
const fb = mk({ rsid: 'rs_fb', het: 0.2952, het_max_pop: 0.5, per_pop_maf: { AFR: 0.4 } })
assert.equal(rankHet(fb, 'EAS'), fb.het) // unsampled population: global, not AFR's 0.5
assert.notEqual(rankHet(fb, 'EAS'), fb.het_max_pop)
assert.equal(rankHet(fb, null), fb.het) // no ancestry (the DEFAULT): global, not het_max_pop
assert.notEqual(rankHet(fb, null), fb.het_max_pop)
assert.equal(rankHet(fb, 'AFR'), twoPQ(0.4)) // sampled: that population's own figure

// shownHet is what a chart may DRAW, and it tracks the axis title or it is nothing.
const drawn = mk({ rsid: 'rs_draw', het: 0.32, het_max_pop: 0.5, per_pop_maf: { AFR: 0.4 } })
assert.equal(shownHet(drawn, null), drawn.het) // no ancestry: the global figure, as titled
assert.equal(shownHet(drawn, 'AFR'), twoPQ(0.4)) // ancestry: the ancestry's own figure
assert.equal(shownHet(drawn, 'EAS'), null) // unsampled population: no height, not zero
assert.notEqual(shownHet(drawn, 'EAS'), 0) // ...and never confusable with 2pq = 0
// Ranking falls back to the global figure so an unsampled population still gets an order;
// drawing under a titled axis has no such licence and shows nothing. Spelled as drawn.het
// rather than a literal: 0.5 is also het_max_pop and the AFR figure, so a literal here
// would not tell the three apart.
assert.equal(rankHet(drawn, 'EAS'), drawn.het)
assert.notEqual(rankHet(drawn, 'EAS'), drawn.het_max_pop)
assert.notEqual(shownHet(drawn, 'EAS'), rankHet(drawn, 'EAS'))
assert.notEqual(shownHet(drawn, 'AFR'), shownHet(drawn, null)) // selecting ancestry moves the bars

// Ancestry selection actually re-orders: rs_afr wins under AFR, loses under NFE.
const afr = mk({ rsid: 'rs_afr', het: 0.3, het_max_pop: 0.5, per_pop_maf: { AFR: 0.5, NFE: 0.05 } })
const nfe = mk({ rsid: 'rs_nfe', het: 0.3, het_max_pop: 0.5, per_pop_maf: { AFR: 0.05, NFE: 0.5 } })
assert.equal(rankMarkers([nfe, afr], 'AFR')[0].rsid, 'rs_afr')
assert.equal(rankMarkers([afr, nfe], 'NFE')[0].rsid, 'rs_nfe')

// Tie on het -> nearer marker wins.
const far = mk({ rsid: 'rs_far', dist: -50_000, per_pop_maf: { AFR: 0.5 } })
const near = mk({ rsid: 'rs_near', dist: -125, per_pop_maf: { AFR: 0.5 } })
assert.equal(rankMarkers([far, near], 'AFR')[0].rsid, 'rs_near')

// Global het breaks an ancestry tie before distance does.
const hiGlobal = mk({ rsid: 'rs_hi', het: 0.49, dist: -9_000, per_pop_maf: { AFR: 0.5 } })
const loGlobal = mk({ rsid: 'rs_lo', het: 0.2, dist: -100, per_pop_maf: { AFR: 0.5 } })
assert.equal(rankMarkers([loGlobal, hiGlobal], 'AFR')[0].rsid, 'rs_hi')

// minPopHet = worst-case population prior; drives 'cross-ancestry robust'.
assert.equal(minPopHet(mk({ rsid: 'r', per_pop_maf: { AFR: 0.5, NFE: 0.1 } })), twoPQ(0.1))
assert.equal(minPopHet(mk({ rsid: 'r' })), 0)
const robust = mk({ rsid: 'rs_robust', per_pop_maf: { AFR: 0.4, NFE: 0.4 } })
const lopsided = mk({ rsid: 'rs_lop', per_pop_maf: { AFR: 0.5, NFE: 0.01 } })
assert.equal(applyPreset([lopsided, robust], 'robust', null)[0].rsid, 'rs_robust')
assert.equal(applyPreset([far, near], 'closest', null)[0].rsid, 'rs_near')

// Sorting must not mutate the caller's array.
const orig = [far, near]
rankMarkers(orig, 'AFR')
assert.equal(orig[0].rsid, 'rs_far')

// posMismatch drives the panel-level alarm in PanelTable. Three states, three distinct
// claims, and collapsing any two of them is the failure this guards. Hand-built: the
// golden ABCC8 fixture cross-checks 'ok' or null on every marker, so no fixture reaches
// the mismatch branch.
assert.equal(posMismatch(mk({ rsid: 'rs_agreed', ensembl_pos_check: 'ok' })), false)
assert.equal(posMismatch(mk({ rsid: 'rs_disputed', ensembl_pos_check: 'MISMATCH:17401130' })), true)
// null is silence, not a pass: only the nearest few shortlisted markers are sent to
// Ensembl, and the check is skipped when cross_check is off or the call throws. Both
// directions must hold, or the alarm cries wolf on ~1200 unchecked rows.
assert.equal(posMismatch(mk({ rsid: 'rs_unchecked', ensembl_pos_check: null })), false)
assert.notEqual(
  posMismatch(mk({ rsid: 'rs_x', ensembl_pos_check: 'MISMATCH:17401130' })),
  posMismatch(mk({ rsid: 'rs_y', ensembl_pos_check: null })),
)

console.log('rank.check.ts: all assertions passed')

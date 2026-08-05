// Self-check for the reconstructed-parent layer. Run: node src/inferredReference.check.ts
//
// The arithmetic is duplicated from `tools/inferred_reference.py` deliberately, so it is checked
// against measurements rather than against itself. Every number quoted below came off a real
// array during validation, and `tools/parity_check.py` proves the two implementations agree on
// the public products to 1e-6 on every quantity a user is shown.
import assert from 'node:assert/strict'
import {
  ProductSet, groupByParent, kinship, pAgree, MIN_ASCERTAINMENT, MIN_PRODUCTS,
  SAME_PARENT_MAX, DIFFERENT_PARENT_MIN,
} from './inferredReference.ts'
import type { ProbeRow } from './ingest.ts'
import type { AB } from './informativity.ts'

const row = (probesetId: string, chrom: string, genotype: AB, baf: number | null = null):
ProbeRow => ({ probesetId, chrom, pos: 1000, log2R: null, baf, copyNumber: null, genotype,
  bestProbeset: true })

/** A parent with known heterozygosity, and haploid products drawn from him. */
function parent(nHom: number, nHet: number): { hom: string[]; het: string[] } {
  return {
    hom: Array.from({ length: nHom }, (_, i) => `hom${i}`),
    het: Array.from({ length: nHet }, (_, i) => `het${i}`),
  }
}

function products(
  p: { hom: string[]; het: string[] }, n: number, draw: (marker: number, product: number) => 0 | 1,
  callRate = 1,
): ProductSet {
  const ps = new ProductSet()
  for (let k = 0; k < n; k += 1) {
    const slot = ps.begin(`p${k}`)
    const baf = { inBand: 0, total: 0 }
    p.hom.forEach((id, i) => {
      if ((i * 7919 + k * 104729) % 1000 >= callRate * 1000) return
      ps.add(slot, row(id, '1', 'AA', 0), baf)
    })
    p.het.forEach((id, i) => {
      if ((i * 7919 + k * 104729) % 1000 >= callRate * 1000) return
      ps.add(slot, row(id, '1', draw(i, k) ? 'BB' : 'AA', 0), baf)
    })
    ps.end(slot, baf)
  }
  return ps
}

// --- 1. the parent's homozygous sites are recovered exactly ------------------------------------
{
  const p = parent(1000, 0)
  const ps = products(p, 5, () => 0)
  const ref = ps.build(4)
  assert.equal(ref.markers, 1000, 'every homozygous site should survive')
  assert.equal(ref.contamination, 0, 'nothing heterozygous went in, so nothing is contaminated')
  assert.equal([...ref.genotype.values()].every((g) => g === 'AA'), true)
}

// --- 2. a heterozygous site enters wrongly only when every product happens to agree -------------
//
// This is the whole error model. With m products called at a marker the chance they all drew the
// same allele is 2^(1-m), inflated at low m because dropout is probe-correlated rather than
// independent.
{
  const p = parent(0, 4000)
  // Deterministic alternating draw: products 0..3 differ at every marker, so nothing agrees.
  const disagreeing = products(p, 4, (i, k) => ((i + k) % 2) as 0 | 1)
  assert.equal(disagreeing.build(4).markers, 0,
    'if the products never agree, no heterozygous site can enter')

  // All products draw identically, so every heterozygous site enters wrongly.
  const agreeing = products(p, 4, () => 0)
  const blind = agreeing.build(4)
  assert.equal(blind.markers, 4000, 'every heterozygous site entered, all of it contamination')

  // And here is the honest limit of the estimator, pinned rather than papered over.
  //
  // Contamination is derived from how often the products DISAGREE. A set that never disagrees
  // is indistinguishable, from the inside, from a parent who is homozygous everywhere, so the
  // recovered heterozygosity is zero and the reported contamination is zero even though this
  // reference is nothing but contamination. Real product sets disagree at thousands of markers,
  // which is why this is a degenerate case and not a live hazard, but the tool cannot detect it
  // and does not pretend to.
  assert.equal(blind.hRetained, 0)
  assert.equal(blind.contamination, 0)

  // Given a heterozygosity to work from, contamination is the POSTERIOR that a marker which
  // agreed is nonetheless heterozygous, so it carries the prior that most markers are genuinely
  // homozygous. At h=0.17 with four products agreeing, only 3.5% of the reference is expected
  // to be contaminated, and that is the right answer for a real parent even though every marker
  // in this synthetic set happens to be heterozygous.
  for (const h of [0.17, 0.5, 0.9]) {
    const told = agreeing.build(4, [], h)
    const odds = h * pAgree(4)
    assert.ok(Math.abs(told.contamination - odds / (1 - h + odds)) < 1e-9,
      `contamination at h=${h} must equal h*pAgree(m)/(1-h+h*pAgree(m))`)
  }
}

// --- 3. m is chosen by ascertainment, not by a fixed offset from n ------------------------------
//
// The rule that `m = n - 1` is right at five products and wrong at eight was measured: at n=8 it
// demanded seven agreeing calls and left the retained set at roughly 62% of the genome's
// heterozygosity.
{
  const p = parent(8000, 2000)
  // Poor call rates, so deep thresholds bite hard and the ladder actually falls.
  const ps = products(p, 8, (i, k) => ((i * 3 + k) % 2) as 0 | 1, 0.55)
  const { mMin, ratios } = ps.chooseM()
  assert.ok(mMin >= 2 && mMin <= 8)
  assert.ok((ratios.get(mMin) ?? 0) >= MIN_ASCERTAINMENT,
    'the chosen depth must hold the ascertainment floor')
  for (const [m, r] of ratios) {
    if (m <= mMin) continue
    // A depth deeper than the chosen one is either too narrow, or no marker reached it at all
    // and its ascertainment is undefined. Neither may be selected, and an undefined one must
    // never be treated as passing.
    assert.ok(!Number.isFinite(r) || r < MIN_ASCERTAINMENT,
      `m=${m} holds ${r} yet was not chosen, so chooseM did not take the deepest`)
  }
  assert.ok(Number.isFinite(ratios.get(mMin) ?? NaN),
    'the chosen depth must have a measurable ascertainment, not an undefined one')
}

// --- 4. leave-one-out genuinely removes a product ----------------------------------------------
//
// A product scored against a reference containing itself reads exactly zero absence, a bias the
// size of the whole signal, so exclusion has to be real rather than cosmetic.
{
  const p = parent(500, 500)
  const ps = products(p, 5, (i, k) => (k === 0 ? 1 : 0))
  const all = ps.build(4)
  const without = ps.build(4, ['p0'])
  assert.equal(without.nProducts, 4)
  assert.ok(without.markers !== all.markers,
    'dropping the one product that drew differently must change the retained set')
  assert.equal(ps.build(4, ['p0', 'p1']).nProducts, 3, 'excluding several must remove several')
}

// --- 5. the agreement probability is the measured one, not the independent one ------------------
{
  assert.ok(Math.abs(pAgree(2) - 0.62) < 1e-9)
  assert.ok(Math.abs(pAgree(3) - 0.3525) < 1e-9)
  assert.ok(Math.abs(pAgree(4) - 0.1775) < 1e-9)
  assert.ok(Math.abs(pAgree(5) - 0.08) < 1e-9)
  // The excess vanishes by six, so a reference deep enough to be usable needs no correction.
  assert.equal(pAgree(6), 2 ** -5)
  assert.equal(pAgree(10), 2 ** -9)
}

// --- 6. membership needs every pair, because one pairwise misread is not hypothetical -----------
//
// Measured on real arrays: one genuine cross-parent pair reads 9.88%, under the same-parent cut,
// and grouping by chained links merges two parents through that single edge. The same product
// reads 12.1% to 12.4% against the rest of that group, so requiring all pairs rejects it.
{
  const A = [0, 1, 2]
  const B = [3, 4, 5]
  const r = (a: number, b: number): number => {
    if (a === b) return 0
    const sameSide = (A.includes(a) && A.includes(b)) || (B.includes(a) && B.includes(b))
    if (sameSide) return 0.07
    // The one bridging pair, exactly as measured, under the cut.
    if ((a === 2 && b === 3) || (a === 3 && b === 2)) return 0.0988
    return 0.13
  }
  assert.equal(kinship(0.0988), 'same parent', 'the bridging pair does read same parent')
  const groups = groupByParent(6, r)
  assert.equal(groups.length, 2, 'all-pairs must not chain through the single bridging edge')
  assert.deepEqual(groups.map((g) => g.length).sort(), [3, 3])
  for (const g of groups) {
    for (const a of g) for (const b of g) {
      if (a < b) assert.ok(r(a, b) < SAME_PARENT_MAX, 'every pair inside a group reads same parent')
    }
  }
}

// --- 7. one parent stays one group ---------------------------------------------------------------
{
  const groups = groupByParent(5, (a, b) => (a === b ? 0 : 0.07))
  assert.equal(groups.length, 1)
  assert.equal(groups[0].length, 5)
}

// --- 8. the thresholds and floors are the measured ones -----------------------------------------
assert.equal(MIN_PRODUCTS, 5)
assert.equal(SAME_PARENT_MAX, 0.105)
assert.equal(DIFFERENT_PARENT_MIN, 0.125)
assert.equal(MIN_ASCERTAINMENT, 0.9)
assert.equal(kinship(0.0468), 'same parent')   // tightest real within-parent pair
assert.equal(kinship(0.097), 'same parent')    // widest real within-parent pair
assert.equal(kinship(0.161), 'different parents')
assert.equal(kinship(0.115), 'ambiguous')      // the gap is 0.18 points wide, so this exists

// --- 9. concordance counts opposite homozygotes and nothing else ---------------------------------
{
  const ps = new ProductSet()
  const baf = { inBand: 0, total: 0 }
  const a = ps.begin('a')
  const b = ps.begin('b')
  for (let i = 0; i < 100; i += 1) ps.add(a, row(`m${i}`, '1', 'AA'), baf)
  for (let i = 0; i < 100; i += 1) ps.add(b, row(`m${i}`, '1', i < 25 ? 'BB' : 'AA'), baf)
  const o = ps.opposite(a, b)
  assert.equal(o.shared, 100)
  assert.equal(o.opposite, 25)
  assert.ok(Math.abs(o.rate - 0.25) < 1e-12)
}

// --- 10. only autosomal homozygous calls become evidence ------------------------------------------
{
  const ps = new ProductSet()
  const baf = { inBand: 0, total: 0 }
  const s = ps.begin('x')
  ps.add(s, row('auto', '1', 'AA'), baf)
  ps.add(s, row('sexX', 'X', 'AA'), baf)
  ps.add(s, row('sexY', 'Y', 'AA'), baf)
  ps.add(s, row('het', '1', 'AB'), baf)
  ps.add(s, row('nc', '1', 'NC'), baf)
  ps.end(s, baf)
  const ref = ps.build(1)
  assert.equal(ref.markers, 1, 'only the autosomal homozygous call is evidence')
  assert.equal(ref.genotype.has('auto'), true)
  assert.equal(ps.hetCalls[s], 1, 'a heterozygous call in a haploid is counted as the error it is')
  assert.equal(ps.called[s], 2, 'called counts genotyped autosomal markers, het included')
}

// --- 11. an empty or unbuildable reference reports nothing, never zero ---------------------------
{
  const ps = new ProductSet()
  const baf = { inBand: 0, total: 0 }
  const s = ps.begin('only')
  ps.add(s, row('m1', '1', 'AA'), baf)
  ps.end(s, baf)
  const impossible = ps.build(5)  // more agreeing calls than there are products
  assert.equal(impossible.markers, 0)
  assert.ok(Number.isNaN(impossible.contamination),
    'contamination of an empty reference is unknown, and reporting 0 would read as perfect')
}

console.log('inferredReference.check.ts OK')

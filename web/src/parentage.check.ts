// Self-check for the browser parentage layer. Run: node src/parentage.check.ts
//
// The arithmetic here is duplicated from `origin.py` deliberately, so it is checked against the
// same real measurements the Python side was validated on rather than against itself. Every
// number in section 1 came off a real file during that validation.
import assert from 'node:assert/strict'
import {
  ABSENCE_MARGIN, absenceExplainable, agreement, classify, emptyTally, isAutosome, pair, pct,
  secondParentSignal, tallyRow, type Tally,
} from './parentage.ts'
import type { ProbeRow } from './ingest.ts'
import type { AB } from './informativity.ts'

const row = (
  chrom: string, pos: number, genotype: AB, baf: number | null = null,
): ProbeRow => ({
  probesetId: `m${pos}`, chrom, pos, log2R: null, baf, copyNumber: null,
  genotype, bestProbeset: true,
})

// --- 1. the noise bound, against pairs of known relationship ----------------------------------
//
// Dropout only manufactures absence by turning a heterozygous call homozygous and discarding the
// parent's allele, so the bound is the product and neither term alone. A homozygous genome is
// immune however much drops out, which is why an androgenote holds at 0.16% absence despite a
// 13.6% no-call rate while a diploid embryo at 46.9% reaches 9.69% against its own father.
assert.equal(absenceExplainable(0.5, 0), 0, 'homozygous: immune to dropout')
assert.equal(absenceExplainable(0, 0.5), 0, 'no dropout: nothing to inflate')
for (const [nc, het, observed] of [
  [0.004, 0.318, 0.0005], // PennCNV offspring
  [0.136, 0.013, 0.0016], // 52461 androgenote 02
  [0.137, 0.020, 0.0022], // 52461 androgenote 04
  [0.469, 0.222, 0.0969], // Zuccaro A8, a TRUE pair scoring above unrelated ones
] as const) {
  assert.ok(absenceExplainable(nc, het) >= observed,
    `bound ${absenceExplainable(nc, het)} must cover observed ${observed}`)
}
assert.ok(Number.isNaN(absenceExplainable(NaN, 0.2)))

// --- 2. the second-parent signal is derived, not fitted ---------------------------------------
assert.ok(Math.abs(secondParentSignal(0.170) - 0.085) < 1e-9)
assert.ok(Math.abs(secondParentSignal(0.322) - 0.161) < 1e-9)
assert.ok(Number.isNaN(secondParentSignal(NaN)))

// --- 3. the three classes, on synthetic genomes with the real shapes ---------------------------
function build(opts: {
  absent: number; het: number; nonParental: number; n?: number; chrom?: string
}): { t: Tally } {
  const t = emptyTally()
  const n = opts.n ?? 4000
  for (let i = 0; i < n; i += 1) {
    const s = (i * 7919) % 1000
    const parent: AB = 'AA'
    let gt: AB = 'AA'
    if (s < opts.absent * 1000) gt = 'BB'
    else if (s >= 1000 - opts.nonParental * 1000) gt = 'AB'
    tallyRow(parent, row(opts.chrom ?? '1', 1000 + i * 1000, gt,
      s < opts.het * 1000 ? 0.5 : 0.0), t)
  }
  return { t }
}

const andro = classify(build({ absent: 0.002, het: 0.02, nonParental: 0.04 }).t, 0.17)
assert.equal(andro.verdict, 'parent_genome_present')
assert.equal(andro.zygosity, 'uniparental_homozygous')
assert.equal(andro.originClass, 'androgenetic')

const gyno = classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
assert.equal(gyno.verdict, 'no_parental_contribution')
assert.equal(gyno.originClass, 'gynogenetic')

const bip = classify(build({ absent: 0.002, het: 0.30, nonParental: 0.30 }).t, 0.17)
assert.equal(bip.zygosity, 'diploid')
assert.equal(bip.originClass, 'biparental')

// A homozygous genome never consults the narrow axis: one allele per locus, and if it is the
// parent's there is no room for a second complement.
const loud = classify(build({ absent: 0.002, het: 0.02, nonParental: 0.30 }).t, 0.17)
assert.ok(loud.nonParentalRate > loud.secondParentExpected, 'the narrow axis would say biparental')
assert.equal(loud.originClass, 'androgenetic', 'but zygosity settles it')

// --- 4. the uncalled band exists and says what would settle it --------------------------------
const mid = classify(build({ absent: 0.012, het: 0.30, nonParental: 0.30 }).t, 0.17)
assert.equal(mid.verdict, 'unclear')
assert.equal(mid.originClass, 'unclear')
assert.ok(mid.limits.some((l) => l.includes('array would measure dropout directly')))

// --- 5. a male offspring loses the paternal X by biology, not by defect ------------------------
{
  const t = build({ absent: 0.002, het: 0.02, nonParental: 0.04 }).t
  for (let i = 0; i < 4000; i += 1) {
    const s = (i * 7919) % 1000
    tallyRow('AA', row('X', 1000 + i * 1000, s < 315 ? 'BB' : 'AA', 0.0), t)
  }
  const pat = classify(t, 0.17, { role: 'paternal' })
  const x = pat.chroms.find((c) => c.chrom === 'X')
  assert.equal(x?.verdict, 'expected_absent', 'the father sent a Y instead')
  assert.equal(pat.spermType, 'Y_bearing')

  // A mother transmits an X to a child of either sex, so the same exemption must not apply.
  const mat = classify(t, 0.17, { role: 'maternal' })
  assert.equal(mat.chroms.find((c) => c.chrom === 'X')?.verdict, 'absent')
}

// --- 6. the X-bearing case a chrY test cannot reach -------------------------------------------
{
  const t = build({ absent: 0.002, het: 0.02, nonParental: 0.04 }).t
  for (let i = 0; i < 4000; i += 1) tallyRow('AA', row('X', 1000 + i * 1000, 'AA', 0.0), t)
  const r = classify(t, 0.17)
  assert.equal(r.spermType, 'X_bearing')
  assert.ok(r.notes.some((n) => n.includes('a chrY test cannot call')))
}

// --- 7. housekeeping ---------------------------------------------------------------------------
assert.ok(isAutosome('1') && isAutosome('22') && !isAutosome('X') && !isAutosome('23'))
assert.equal(pct(0.0016), '0.16%')
assert.equal(pct(NaN), 'n/a')
assert.equal(ABSENCE_MARGIN, 3)

// No-calls are excluded rather than counted as evidence in either direction.
{
  const t = emptyTally()
  for (let i = 0; i < 500; i += 1) tallyRow('AA', row('1', i * 1000, 'NC'), t)
  assert.equal(t.byChrom.size, 0, 'a no-call is not an informative marker')
  assert.equal(t.called, 0)
  assert.ok(t.markers === 500)
}


// --- 8. both parents, mirroring origin.both_parents --------------------------------------------
{
  const present = () => classify(build({ absent: 0.002, het: 0.30, nonParental: 0.30 }).t, 0.17)
  const absent = () => classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
  const middle = () => classify(build({ absent: 0.012, het: 0.30, nonParental: 0.30 }).t, 0.17)

  assert.equal(pair(present(), present(), NaN).originClass, 'biparental')
  assert.equal(pair(present(), absent(), NaN).originClass, 'androgenetic')
  assert.equal(pair(absent(), present(), NaN).originClass, 'gynogenetic')

  // The outcome a father-only run cannot reach: his absence alone reads as gynogenetic.
  const neither = pair(absent(), absent(), NaN)
  assert.equal(neither.originClass, 'neither_parent')
  assert.equal(classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
    .originClass, 'gynogenetic', 'which is exactly what one parent would have called it')
  assert.ok(neither.notes.some((n) => n.includes('check the pairing')))

  // An unresolved parent is never resolved by the other one.
  assert.equal(pair(middle(), present(), NaN).originClass, 'unclear')
  assert.equal(pair(present(), middle(), NaN).originClass, 'unclear')

  assert.ok(pair(present(), present(), 0.998).notes.some((n) => n.includes('relabelled file')))
  assert.ok(!pair(present(), present(), 0.62).notes.some((n) => n.includes('relabelled file')))
  assert.ok(!pair(present(), present(), NaN).notes.some((n) => n.includes('relabelled file')),
    'too few shared markers must not raise a duplicate warning')
}

// --- 9. parent agreement -----------------------------------------------------------------------
{
  const gt = (n: number, f: (i: number) => AB): Map<string, AB> =>
    new Map(Array.from({ length: n }, (_, i) => [`m${i}`, f(i)]))
  const a = gt(2000, (i) => (i % 3 === 0 ? 'AB' : 'AA'))
  assert.equal(agreement(a, a), 1, 'a file against itself')
  const third = gt(2000, (i) => (i % 3 === 0 ? 'AB' : 'BB'))
  assert.equal(agreement(a, third), 667 / 2000, 'agrees only where both say AB')
  assert.ok(Number.isNaN(agreement(a, gt(50, () => 'AA'))), 'too few shared markers')
  // A no-call on either side is not a disagreement, it is no observation.
  assert.equal(agreement(a, new Map([...a].map(([k, v], i) => [k, i % 4 ? v : 'NC' as AB]))), 1)
}

console.log('parentage.check.ts OK')

// The aggregation that answers the question this tool exists for: do breakpoints and losses fall
// on the two parental genomes differentially, or equally.
import assert from 'node:assert/strict'
import {
  parentalBalance, pairedWithinSample, mergeEvents, minAchievableP, MAX_POWER_SKEW,
  MIN_PER_GROUP, REPORTING_PER_GROUP,
} from './parentalBalance.ts'
import type { BalanceSample } from './parentalBalance.ts'

// Spread across chromosomes and far enough apart that merging keeps them distinct, so a test of
// the statistic is not accidentally a test of the merge.
const ev = (n: number, cls = 'loss') =>
  Array.from({ length: n }, (_, i) => ({
    cls, chrom: String((i % 22) + 1),
    startBp: Math.floor(i / 22) * 20e6, endBp: Math.floor(i / 22) * 20e6 + 5e5,
  }))

const sample = (
  name: string, parent: 'maternal' | 'paternal' | null, events: number,
  informative = 500_000, cls = 'loss', explainable = 0.01,
): BalanceSample => ({
  name, parent, informative, explainable, material: 'esc-single',
  originClass: parent === 'maternal' ? 'gynogenetic' : parent === 'paternal' ? 'androgenetic' : 'unclear',
  events: ev(events, cls),
})

// --- 1. A PLANTED DIFFERENCE IS FOUND ----------------------------------------------------------
{
  const rows = [
    ...Array.from({ length: 6 }, (_, i) => sample(`m${i}`, 'maternal', 20)),
    ...Array.from({ length: 6 }, (_, i) => sample(`p${i}`, 'paternal', 4)),
  ]
  const r = parentalBalance(rows, { permutations: 2000 })
  assert.equal(r.verdict, 'differential', `a 5x planted difference must be found: ${r.headline}`)
  assert.ok(r.headline.includes('UNEQUALLY'), r.headline)
  assert.ok(r.headline.includes('carry any change'), 'the per-genome binary must be in the headline')
  console.log(`  planted 5x maternal excess -> ${r.verdict}, p ${r.p.toFixed(4)}`)
}

// --- 2. NO DIFFERENCE IS NOT DRESSED UP AS ONE -------------------------------------------------
{
  const counts = [11, 9, 10, 12, 8, 10, 9, 11, 10, 10]
  const rows = counts.map((c, i) => sample(`s${i}`, i % 2 ? 'paternal' : 'maternal', c))
  const r = parentalBalance(rows, { permutations: 2000 })
  assert.equal(r.verdict, 'equal', `matched groups must read equal: ${r.headline}`)
  assert.ok(r.p > 0.05, `p should not be significant, got ${r.p}`)
  assert.ok(r.headline.includes('not the same as showing there is no difference'))
  console.log(`  matched groups -> ${r.verdict}, p ${r.p.toFixed(4)}`)
}

// --- 3. THE STATISTIC MUST HOLD ITS SIZE UNDER A CONFOUNDED DENOMINATOR ------------------------
//
// THE DEFECT THIS REPLACED. The pooled rate sum(events)/sum(markers) is not a symmetric function
// of the labels, because each sample's denominator travels with it, so the permutation null does
// not have the right size. Measured under a TRUE null: the pooled form ran at 0.164 at a marker
// skew of 1.5, more than three times nominal. The rank of per-sample counts does not use a
// denominator at all, so a marker skew cannot inflate it.
{
  let seed = 12345
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  let rejected = 0
  const TRIALS = 200
  for (let t = 0; t < TRIALS; t += 1) {
    // Both groups drawn from the SAME distribution: any rejection is a false positive. Markers are
    // deliberately skewed 1.5x between the groups, which is what broke the pooled statistic.
    const rows = [
      ...Array.from({ length: 6 }, (_, i) =>
        sample(`m${t}_${i}`, 'maternal', 3 + Math.floor(rnd() * 10), 600_000)),
      ...Array.from({ length: 6 }, (_, i) =>
        sample(`p${t}_${i}`, 'paternal', 3 + Math.floor(rnd() * 10), 400_000)),
    ]
    if (parentalBalance(rows, { permutations: 400, seed: t + 1 }).p < 0.05) rejected += 1
  }
  const rate = rejected / TRIALS
  console.log(`  type-I under a true null with 1.5x marker skew: ${rate.toFixed(3)}`
    + ' (pooled rate measured 0.164 here)')
  assert.ok(rate <= 0.12,
    `false positives ran at ${rate.toFixed(3)} against a nominal 0.05, which is the size failure `
    + 'the rank statistic replaced the pooled rate to fix')
}

// --- 4. UNEQUAL ARTEFACT PROPENSITY IS REFUSED -------------------------------------------------
// And it is measured on the explainable-noise ceiling, not the marker count. Marker count was
// measured to AMPLIFY this confounder rather than remove it.
{
  const rows = [
    ...Array.from({ length: 6 }, (_, i) => sample(`m${i}`, 'maternal', 20, 500_000, 'loss', 0.008)),
    ...Array.from({ length: 6 }, (_, i) => sample(`p${i}`, 'paternal', 20, 500_000, 'loss', 0.080)),
  ]
  const r = parentalBalance(rows, { permutations: 500 })
  assert.equal(r.verdict, 'not-comparable', `a 10x noise skew must refuse: ${r.headline}`)
  assert.ok(!r.power.withinTolerance && r.power.skew > MAX_POWER_SKEW)
  assert.ok(r.headline.includes('artefact'), r.headline)
  console.log(`  ${r.power.skew.toFixed(1)}x artefact-propensity skew -> ${r.verdict}`)
}

// --- 4b. AND SO IS A QUALITY-CORRELATED EXCLUSION ----------------------------------------------
// The skew above measures only the samples that SURVIVED. In this project's own run the three
// arrays that fell out as unclear were the three worst and all three were paternal, so the
// exclusion silently repaired the measured balance by dropping one group's weakest members.
{
  const rows = [
    ...Array.from({ length: 6 }, (_, i) => sample(`m${i}`, 'maternal', 10)),
    ...Array.from({ length: 6 }, (_, i) => sample(`p${i}`, 'paternal', 10)),
    ...Array.from({ length: 3 }, (_, i) => ({
      ...sample(`drop${i}`, null, 2), declaredParent: 'paternal' as const })),
  ]
  const r = parentalBalance(rows, { permutations: 500 })
  assert.equal(r.verdict, 'not-comparable', `a one-sided exclusion must refuse: ${r.headline}`)
  assert.equal(r.power.exclusion.paternal, 3)
  assert.equal(r.power.exclusion.maternal, 0)
  assert.ok(r.headline.includes('excluded unevenly'), r.headline)
  console.log(`  ${r.power.exclusion.paternal} paternal vs ${r.power.exclusion.maternal} maternal`
    + ` exclusions -> ${r.verdict}`)
}

// --- 5. TOO FEW GENOMES ------------------------------------------------------------------------
// Three per group cannot clear alpha 0.05 at any effect size, so the old floor of 3 permitted a
// test incapable of producing a finding.
{
  assert.equal(MIN_PER_GROUP, 5, 'the hard floor is where the achievable p first clears a '
    + 'Bonferroni threshold for three classes')
  const rows = [
    ...Array.from({ length: 3 }, (_, i) => sample(`m${i}`, 'maternal', 20)),
    ...Array.from({ length: 3 }, (_, i) => sample(`p${i}`, 'paternal', 2)),
  ]
  const r = parentalBalance(rows, { permutations: 500 })
  assert.equal(r.verdict, 'underpowered')
  assert.ok(r.headline.includes('No conclusion either way'))
  assert.ok(Math.abs(minAchievableP(3, 3, 500) - 0.100) < 1e-9,
    'the exact two-sided floor at 3 vs 3 is 2/C(6,3) = 0.100')
  assert.ok(Math.abs(minAchievableP(4, 4, 5000) - 2 / 70) < 1e-9,
    'and 0.0286 at 4 vs 4, not the 0.0423 the previous form returned')
  console.log(`  3 vs 3 -> ${r.verdict}; exact floors 0.100 at 3v3 and 0.0286 at 4v4`)
}

// --- 6. BETWEEN THE FLOORS, A RESULT IS LABELLED EXPLORATORY -----------------------------------
{
  assert.ok(REPORTING_PER_GROUP > MIN_PER_GROUP)
  const rows = [
    ...Array.from({ length: 6 }, (_, i) => sample(`m${i}`, 'maternal', 20)),
    ...Array.from({ length: 6 }, (_, i) => sample(`p${i}`, 'paternal', 4)),
  ]
  const r = parentalBalance(rows, { permutations: 2000 })
  assert.ok(r.headline.includes('EXPLORATORY'),
    `below ${REPORTING_PER_GROUP} per group the result must say so: ${r.headline}`)
  console.log(`  6 vs 6 is reported but labelled exploratory`)
}

// --- 7. THE UNIT IS DISTINCT EVENTS, NOT WINDOWS -----------------------------------------------
// A sliding detector cuts one change into as many rows as its windows, so a raw region count is
// defensible only as "windows".
{
  const sliced = Array.from({ length: 12 }, (_, i) => ({
    cls: 'loss', chrom: '4', startBp: 10e6 + i * 5e5, endBp: 10e6 + i * 5e5 + 5e5 }))
  assert.equal(mergeEvents(sliced).length, 1,
    '12 abutting windows on one chromosome are one distinct event')
  const apart = [
    { cls: 'loss', chrom: '4', startBp: 10e6, endBp: 11e6 },
    { cls: 'loss', chrom: '4', startBp: 90e6, endBp: 91e6 },
  ]
  assert.equal(mergeEvents(apart).length, 2, 'and two genuinely separate ones stay two')
  const diffClass = [
    { cls: 'loss', chrom: '4', startBp: 10e6, endBp: 11e6 },
    { cls: 'gain', chrom: '4', startBp: 10e6, endBp: 11e6 },
  ]
  assert.equal(mergeEvents(diffClass).length, 2, 'a loss and a gain are never merged together')
  console.log('  12 abutting windows merge to 1 event; separate and differently-classed stay apart')
}

// --- 8. REPRODUCIBLE ---------------------------------------------------------------------------
{
  const rows = [
    ...Array.from({ length: 6 }, (_, i) => sample(`m${i}`, 'maternal', 14)),
    ...Array.from({ length: 6 }, (_, i) => sample(`p${i}`, 'paternal', 6)),
  ]
  const a = parentalBalance(rows, { permutations: 1000, seed: 42 })
  const b = parentalBalance(rows, { permutations: 1000, seed: 42 })
  assert.equal(a.p, b.p, 'the same seed must give the same p')
  console.log(`  same seed, same p (${a.p.toFixed(4)})`)
}


// --- 9. THE WITHIN-EMBRYO DESIGN ---------------------------------------------------------------
// Both counts come from the SAME array, so array quality, marker count and amplification are
// identical on the two sides by construction and cannot explain a difference. That is what the
// between-group comparison cannot claim.
{
  // Ten embryos, nine of them carrying more maternal events.
  const rows = Array.from({ length: 10 }, (_, i) => ({
    name: `e${i}`, maternalEvents: i === 0 ? 1 : 5, paternalEvents: i === 0 ? 4 : 1,
  }))
  const r = pairedWithinSample(rows)
  assert.equal(r.verdict, 'differential', r.headline)
  assert.equal(r.informative, 10)
  assert.equal(r.maternalHigher, 9)
  // Exact two-sided binomial, 9 of 10: 2 * (C(10,0)+C(10,1)) / 2^10 = 22/1024.
  assert.ok(Math.abs(r.p - 22 / 1024) < 1e-12, `exact sign-test p, got ${r.p}`)
  console.log(`  9 of 10 embryos maternal-heavy -> ${r.verdict}, exact p ${r.p.toFixed(4)}`)
}

// A balanced set must not read as a difference.
{
  const rows = Array.from({ length: 10 }, (_, i) => ({
    name: `e${i}`, maternalEvents: i % 2 ? 5 : 1, paternalEvents: i % 2 ? 1 : 5 }))
  const r = pairedWithinSample(rows)
  assert.equal(r.verdict, 'equal')
  assert.equal(r.p, 1, 'five each way is the least extreme outcome there is')
  console.log(`  5 each way -> ${r.verdict}, p ${r.p.toFixed(4)}`)
}

// Ties carry no direction and are excluded, and below six informative samples even perfect
// agreement cannot clear alpha, so nothing is reported.
{
  const tied = Array.from({ length: 20 }, (_, i) => ({
    name: `t${i}`, maternalEvents: 3, paternalEvents: 3 }))
  const r = pairedWithinSample(tied)
  assert.equal(r.informative, 0, 'ties carry no direction')
  assert.equal(r.verdict, 'underpowered')

  const five = Array.from({ length: 5 }, (_, i) => ({
    name: `f${i}`, maternalEvents: 9, paternalEvents: 1 }))
  const r5 = pairedWithinSample(five)
  assert.equal(r5.verdict, 'underpowered',
    'five informative samples cannot reach alpha 0.05 even in perfect agreement, since '
    + '2 x 0.5^4 = 0.125')
  console.log('  ties excluded; below 6 informative samples nothing is reported')
}

console.log('parentalBalance: rank statistic holds its size, artefact propensity is the matched '
  + 'quantity, and the unit is distinct events')

// The aggregation that answers the question this tool exists for: do breakpoints and losses fall
// on the two parental genomes differentially, or equally.
import assert from 'node:assert/strict'
import { parentalBalance, MAX_POWER_SKEW, MIN_PER_GROUP } from './parentalBalance.ts'
import type { BalanceSample } from './parentalBalance.ts'

const ev = (n: number, cls = 'cnn-loh') =>
  Array.from({ length: n }, (_, i) => ({ cls, chrom: String((i % 22) + 1), startBp: i * 1e6, endBp: i * 1e6 + 5e5 }))

const sample = (
  name: string, parent: 'maternal' | 'paternal' | null, events: number,
  informative = 500_000, cls = 'cnn-loh',
): BalanceSample => ({
  name, parent, informative, material: 'esc-single',
  originClass: parent === 'maternal' ? 'gynogenetic' : parent === 'paternal' ? 'androgenetic' : 'unclear',
  events: ev(events, cls),
})

// --- 1. A PLANTED DIFFERENCE IS FOUND ----------------------------------------------------------
{
  const rows = [
    ...[0, 1, 2, 3].map((i) => sample(`m${i}`, 'maternal', 20)),
    ...[0, 1, 2, 3].map((i) => sample(`p${i}`, 'paternal', 4)),
  ]
  const r = parentalBalance(rows, { permutations: 2000 })
  assert.equal(r.verdict, 'differential', `5x planted difference must be found: ${r.headline}`)
  assert.ok(r.fold > 3, `maternal should carry more, fold ${r.fold}`)
  assert.ok(r.headline.includes('UNEQUALLY'), r.headline)
  console.log(`  planted 5x maternal excess -> ${r.verdict}, fold ${r.fold.toFixed(2)}, p ${r.p.toFixed(4)}`)
}

// --- 2. NO DIFFERENCE IS NOT DRESSED UP AS ONE -------------------------------------------------
// The failure that matters most here: a test that finds a difference in noise would put a false
// headline on a lab's report.
{
  const counts = [11, 9, 10, 12, 8, 10, 9, 11]
  const rows = counts.map((c, i) =>
    sample(`s${i}`, i % 2 ? 'paternal' : 'maternal', c))
  const r = parentalBalance(rows, { permutations: 2000 })
  assert.equal(r.verdict, 'equal', `matched groups must read equal: ${r.headline}`)
  assert.ok(r.p > 0.05, `p should not be significant, got ${r.p}`)
  assert.ok(r.headline.includes('not the same as showing there is no difference'),
    'an equal result must say what it does and does not establish')
  console.log(`  matched groups -> ${r.verdict}, p ${r.p.toFixed(4)}, and it states its own limit`)
}

// --- 3. UNEQUAL POWER IS REFUSED, NOT REPORTED -------------------------------------------------
// THE ONE THAT WOULD HAVE PRODUCED A CONFIDENT WRONG ANSWER. Detection scales with how many
// markers could have carried evidence. Give one group half the markers and it will show fewer
// events whatever the biology; reporting that as a parental difference reports the arrays.
{
  const rows = [
    ...[0, 1, 2, 3].map((i) => sample(`m${i}`, 'maternal', 20, 580_000)),
    ...[0, 1, 2, 3].map((i) => sample(`p${i}`, 'paternal', 20, 240_000)),
  ]
  const r = parentalBalance(rows, { permutations: 2000 })
  assert.equal(r.verdict, 'not-comparable',
    `a ${(580 / 240).toFixed(1)}x power skew must refuse: ${r.headline}`)
  assert.ok(r.power.skew > MAX_POWER_SKEW && !r.power.withinTolerance)
  assert.ok(r.headline.includes('reporting the arrays'), r.headline)
  console.log(`  ${r.power.skew.toFixed(2)}x power skew -> ${r.verdict}, refused rather than reported`)
}

// --- 4. TOO FEW GENOMES ------------------------------------------------------------------------
{
  const rows = [sample('m0', 'maternal', 20), sample('m1', 'maternal', 18), sample('p0', 'paternal', 3)]
  const r = parentalBalance(rows, { permutations: 500 })
  assert.equal(r.verdict, 'underpowered')
  assert.ok(r.headline.includes('No conclusion either way'), r.headline)
  assert.ok(MIN_PER_GROUP >= 3, 'the floor is a stated constant')
  console.log(`  ${r.groups[0].samples} vs ${r.groups[1].samples} genomes -> ${r.verdict}`)
}

// --- 5. BIPARENTAL AND UNCLEAR SAMPLES ARE EXCLUDED, AND SAID SO -------------------------------
// Dropping them silently would let a reader think the run had fewer samples than it did.
{
  const rows = [
    ...[0, 1, 2].map((i) => sample(`m${i}`, 'maternal', 10)),
    ...[0, 1, 2].map((i) => sample(`p${i}`, 'paternal', 10)),
    sample('unclear1', null, 5),
    { ...sample('nomarkers', 'maternal', 5), informative: 0 },
  ]
  const r = parentalBalance(rows, { permutations: 500 })
  assert.equal(r.excluded.length, 2, `both untestable samples must be named: ${JSON.stringify(r.excluded)}`)
  assert.ok(r.excluded.some((e) => e.name === 'unclear1' && /two candidate parents/.test(e.why)))
  assert.ok(r.excluded.some((e) => e.name === 'nomarkers' && /denominator/.test(e.why)))
  console.log(`  ${r.excluded.length} samples excluded and each says why`)
}

// --- 6. PER CLASS, WHICH IS THE FORM THE QUESTION WAS ASKED IN ---------------------------------
// "breakpoints/losses" is per class, not one lumped number, and the correction is for how many
// classes were actually tested. Group sizes are the ones a real run had: 8 gynogenetic and 7
// androgenetic.
{
  const rows = [
    ...Array.from({ length: 8 }, (_, i) => ({
      ...sample(`m${i}`, 'maternal', 0),
      events: [...ev(15, 'loss'), ...ev(5, 'cnn-loh')],
    })),
    ...Array.from({ length: 7 }, (_, i) => ({
      ...sample(`p${i}`, 'paternal', 0),
      events: [...ev(2, 'loss'), ...ev(5, 'cnn-loh')],
    })),
  ]
  const r = parentalBalance(rows, { permutations: 4000 })
  const loss = r.byClass.find((c) => c.cls === 'loss')
  const cnn = r.byClass.find((c) => c.cls === 'cnn-loh')
  assert.ok(loss && cnn, 'both classes reported')
  assert.ok(loss.fold > 3, `loss is planted maternal-heavy, fold ${loss.fold}`)
  assert.ok(Math.abs(cnn.fold - 1) < 0.01, `cnn-loh is planted equal, fold ${cnn.fold}`)
  assert.ok(loss.significant && !cnn.significant,
    'the planted class is flagged and the matched one is not')
  assert.ok(!loss.underpowered && !cnn.underpowered,
    '8 vs 7 genomes can reach a corrected threshold, so neither is underpowered')
  console.log(`  8 vs 7 genomes, per class: loss fold ${loss.fold.toFixed(2)} (flagged), `
    + `cnn-loh fold ${cnn.fold.toFixed(2)} (not flagged)`)
}

// --- 6b. A NULL THAT COULD NOT HAVE BEEN ANYTHING ELSE ------------------------------------------
// THE MISREADING THIS PREVENTS. Four genomes against four admit only 70 label assignments, so no
// p below about 0.029 exists and a corrected per-class threshold sits under it. The test cannot
// flag anything at any effect size. Reported as "no difference" that is a false negative dressed
// as a result, and a lab would take it as evidence the parents are alike.
{
  const rows = [
    ...Array.from({ length: 4 }, (_, i) => ({
      ...sample(`m${i}`, 'maternal', 0), events: [...ev(15, 'loss'), ...ev(5, 'cnn-loh')],
    })),
    ...Array.from({ length: 4 }, (_, i) => ({
      ...sample(`p${i}`, 'paternal', 0), events: [...ev(2, 'loss'), ...ev(5, 'cnn-loh')],
    })),
  ]
  const r = parentalBalance(rows, { permutations: 4000 })
  const loss = r.byClass.find((c) => c.cls === 'loss')
  assert.ok(loss, 'the class is still reported')
  assert.ok(loss.fold > 3, 'the effect is still large and still shown')
  assert.ok(loss.underpowered,
    'and it must be marked as unflaggable at this group size rather than read as no difference')
  assert.ok(r.minAchievableP > 0.02,
    `4 vs 4 admits no p below ~0.029, got ${r.minAchievableP}`)
  console.log(`  4 vs 4 genomes: a ${loss.fold.toFixed(1)}x effect cannot be flagged at all`
    + ` (floor p ${r.minAchievableP.toFixed(4)}), and it says so instead of reading as equal`)
}

// --- 7. REPRODUCIBLE ---------------------------------------------------------------------------
// A report a lab keeps must give the same number when it is regenerated.
{
  const rows = [
    ...[0, 1, 2, 3].map((i) => sample(`m${i}`, 'maternal', 14)),
    ...[0, 1, 2, 3].map((i) => sample(`p${i}`, 'paternal', 6)),
  ]
  const a = parentalBalance(rows, { permutations: 1000, seed: 42 })
  const b = parentalBalance(rows, { permutations: 1000, seed: 42 })
  assert.equal(a.p, b.p, 'the same seed must give the same p')
  console.log(`  same seed, same p (${a.p.toFixed(4)})`)
}

console.log('parentalBalance: finds a real difference, refuses an unequal comparison, and is reproducible')

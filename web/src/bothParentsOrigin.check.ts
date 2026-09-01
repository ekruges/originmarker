// Self-check for the two-parent origin rule. Run: node src/bothParentsOrigin.check.ts
//
// The rule exists because one measured asymmetry is large: the same channel, on the same arrays and
// the same events, is right 0.8539 of the time about the loaded parent's own copy and 0.2890 about
// the other parent's. So the ONE property that matters here is that a named parent can only ever
// come from a known-parent-lost side. Everything else in this file is a way of trying to break
// that.
import assert from 'node:assert/strict'
import { callBothParentsOrigin, TWO_PARENT_ACCURACY } from './bothParentsOrigin.ts'
import { OBVIOUS_EVENT_ACCURACY, type OneParentVerdict, type OneParentCall } from './oneParentOrigin.ts'
import type { Band } from './originPosterior.ts'

const call = (verdict: OneParentVerdict, posterior = 0.79, band: Band = 'C'): OneParentCall => ({
  verdict, posterior, band, markers: 4_000, exclusive: verdict === 'known-parent-lost' ? 400 : 4,
  heterozygous: verdict === 'both-present' ? 480 : 6, q: 0.1217,
  why: `fixture ${verdict}`,
})

const VERDICTS: OneParentVerdict[]
  = ['known-parent-lost', 'other-parent-lost', 'both-present', 'refused']

// --- 1. THE ASYMMETRY THIS RULE IS BUILT ON MUST STILL BE THERE ---------------------------------
//
// If the two directions ever measure the same, the rule below is arbitrary and should be revisited
// rather than left standing on a number that moved.
{
  const m = OBVIOUS_EVENT_ACCURACY
  assert.ok(m.knownParentLost.accuracy > m.otherParentLost.accuracy * 2,
    'this module selects one direction over the other because they were measured far apart: '
    + `${m.knownParentLost.accuracy} against ${m.otherParentLost.accuracy}. If that gap closes, `
    + 'the selection rule needs re-deriving, not re-running.')
}

// --- 2. A NAMED PARENT COMES ONLY FROM A known-parent-lost SIDE ---------------------------------
//
// Exhaustive over all 16 combinations. No input that lacks a known-parent-lost may name anyone.
{
  let namedFromWeak = 0
  let namedTotal = 0
  for (const p of VERDICTS) {
    for (const m of VERDICTS) {
      const r = callBothParentsOrigin(call(p), call(m))
      if (r.parent === null) continue
      namedTotal += 1
      if (p !== 'known-parent-lost' && m !== 'known-parent-lost') namedFromWeak += 1
      // And the side it came from must be the side that said it.
      assert.equal(r.parent === 'paternal' ? p : m, 'known-parent-lost',
        `parent ${r.parent} named from a ${r.parent === 'paternal' ? p : m} side`)
    }
  }
  assert.equal(namedFromWeak, 0,
    'a parent was named without either side reporting its own copy absent, which is the one thing '
    + 'this rule exists to prevent')
  assert.ok(namedTotal > 0, 'and the sweep must actually reach a named parent, or it asserts nothing')
}

// --- 3. THE DIRECTIONS ARE THE RIGHT WAY ROUND --------------------------------------------------
//
// Getting this backwards would name the wrong parent at high confidence on every event, which is
// the exact failure this tool is built to avoid, so it is pinned rather than assumed.
{
  const patGone = callBothParentsOrigin(call('known-parent-lost'), call('both-present'))
  assert.equal(patGone.verdict, 'parent-lost')
  assert.equal(patGone.parent, 'paternal',
    'the PATERNAL array reporting its own copy absent means the PATERNAL copy is gone')
  const matGone = callBothParentsOrigin(call('both-present'), call('known-parent-lost'))
  assert.equal(matGone.parent, 'maternal')
  // The confidence and band are the reporting side's own, not a combined quantity.
  const side = call('known-parent-lost', 0.83, 'C')
  const out = callBothParentsOrigin(side, call('refused'))
  assert.equal(out.posterior, side.posterior, 'no likelihood is combined here')
  assert.equal(out.band, side.band)
}

// --- 4. BOTH SIDES CLAIMING THEIR OWN COPY IS GONE NAMES NOBODY ---------------------------------
//
// That is nullisomy. It happens, and on amplified material it is far more often the array. Either
// way a rule that picked one of the two would be picking at random.
{
  const c = callBothParentsOrigin(call('known-parent-lost'), call('known-parent-lost'))
  assert.equal(c.verdict, 'contradiction')
  assert.equal(c.parent, null, 'two confident opposite answers must not resolve to one of them')
  assert.ok(Number.isNaN(c.posterior), 'and it must carry no confidence')
  assert.equal(c.band, 'F')
  assert.ok(/nullisomy/.test(c.why), 'the reason has to be legible to whoever reads the row')
}

// --- 5. THE WEAK DIRECTION IS RECORDED AND NEVER ACTED ON ---------------------------------------
{
  const agreeing = callBothParentsOrigin(call('known-parent-lost'), call('other-parent-lost'))
  assert.equal(agreeing.parent, 'paternal')
  assert.equal(agreeing.corroborated, true,
    'the maternal array pointing at the same parent is worth reporting')
  const silent = callBothParentsOrigin(call('known-parent-lost'), call('refused'))
  assert.equal(silent.parent, 'paternal', 'and a silent other side does not withdraw the call')
  assert.equal(silent.corroborated, false)
  // Corroboration must not be able to CREATE a call: two weak sides pointing at each other is the
  // 0.289 direction twice over, which is not evidence.
  const weakOnly = callBothParentsOrigin(call('other-parent-lost'), call('other-parent-lost'))
  assert.equal(weakOnly.parent, null,
    'two other-parent-lost calls are the unreliable direction twice, not a corroborated answer')
}

// --- 6. AND SILENCE IS DISTINGUISHED FROM AGREEMENT THAT NOTHING IS WRONG -----------------------
//
// A row that says "both copies present" and a row that says "neither array could tell" must not
// look the same to a reader: one is a result and the other is a gap in the data.
{
  const fine = callBothParentsOrigin(call('both-present'), call('both-present'))
  assert.equal(fine.verdict, 'both-present')
  assert.equal(fine.corroborated, true)
  const blind = callBothParentsOrigin(call('refused'), call('refused'))
  assert.equal(blind.verdict, 'refused')
  assert.ok(Number.isNaN(blind.posterior))
  assert.notEqual(fine.verdict, blind.verdict)
}

// --- 7. every string an operator reads is clean --------------------------------------------------
{
  for (const p of VERDICTS) {
    for (const m of VERDICTS) {
      const r = callBothParentsOrigin(call(p), call(m))
      assert.ok(r.why.length > 0)
      assert.ok(!r.why.includes('—'), `em dash in "${r.why.slice(0, 60)}"`)
    }
  }
}


// --- 8. THE MEASURED RESULT, PINNED ------------------------------------------------------------
//
// These are not decoration. The rule above is only worth its cost if the direction it selects is
// far better than the one it drops, and if the extra false calls it buys stay where they were
// measured. A future re-run that loses either property should fail here rather than ship quietly.
{
  const m = TWO_PARENT_ACCURACY
  assert.ok(m.direction.accuracy > m.directionOneArray.accuracy * 2,
    'the whole point is that the selected direction is far better than the dropped one: '
    + `${m.direction.accuracy} against ${m.directionOneArray.accuracy}`)
  assert.equal(m.direction.wrongParent, 0,
    'naming the WRONG parent is the failure this tool exists to avoid, and the measured count is '
    + 'zero. A re-measurement that is not zero is a different tool and must be looked at, not '
    + 'shipped.')
  assert.ok(m.direction.named > 200, 'and it must rest on a real number of calls')

  // THE COST IS REAL AND IS NOT HIDDEN. Loading a second array asks two questions where one was
  // asked, so it can only add false calls. The check asserts the direction of that, not its
  // absence, because asserting its absence would be asserting something untrue.
  assert.ok(m.specificity.falseCallRate >= 0.0084,
    'two arrays cannot produce FEWER false calls than one; if a re-measurement says otherwise, '
    + 'the arms are not being compared on the same material')

  // AND IT IS CONCENTRATED ON ONE MATERIAL, which is what makes it actionable rather than a
  // blanket caveat.
  const bm = m.specificityByMaterial
  for (const [name, row] of Object.entries(bm)) {
    if (name === 'blastomere') continue
    assert.equal(row.falseCallRate, 0,
      `${name} was measured at no false calls with both arrays loaded; if that moves, the claim `
      + 'that the second array is free on multi-cell material no longer holds')
  }
  assert.ok(bm.blastomere.falseCallRate > bm.trophectoderm.falseCallRate,
    'a single cell and a five-to-ten-cell biopsy do not share an error rate here either')
  assert.ok(bm.blastomere.falseCallRate > bm.blastomere.oneArray,
    'and on a single cell the second array does cost something, which the report must say')
  let arrays = 0
  for (const row of Object.values(bm)) arrays += row.arrays
  assert.ok(arrays >= 70 && arrays <= m.specificity.arrays,
    `the split must account for the arrays behind the pooled figure: ${arrays}`)
}


// --- 8. CORROBORATION IS REPORTED AND NOT REQUIRED, and the measurement says why ----------------
//
// The obvious next move on a 0.0813 blastomere false-call rate is to demand that both arrays agree.
// It was measured before it was rejected, and this pins the arithmetic so nobody has to re-derive
// it from an intuition that the filter "must" be worth it.
{
  const c = TWO_PARENT_ACCURACY.corroboration
  const inTrue = c.trueCalls.agreed / c.trueCalls.total
  const inFalse = c.falseCalls.agreed / c.falseCalls.total
  assert.ok(inTrue > inFalse * 3,
    `corroboration must carry signal to be worth showing: ${inTrue.toFixed(3)} of true calls `
    + `against ${inFalse.toFixed(3)} of false ones`)
  // And the reason it is not a gate: it would cost more true calls than it saves false ones.
  const trueLost = c.trueCalls.total - c.trueCalls.agreed
  const falseSaved = c.falseCalls.total - c.falseCalls.agreed
  assert.ok(trueLost > falseSaved * 4,
    `requiring corroboration would drop ${trueLost} real calls to remove ${falseSaved} false ones, `
    + 'which is the trade this constant exists to record')
  console.log(`  corroboration: ${(inTrue * 100).toFixed(1)}% of true calls, `
    + `${(inFalse * 100).toFixed(1)}% of false; as a gate it would cost ${trueLost} real calls `
    + `to remove ${falseSaved} false ones`)
}

console.log('bothParentsOrigin.check.ts: a parent is named only from the direction measured at '
  + '0.854, both sides claiming their own loss names nobody, and the weak direction never creates '
  + 'a call')

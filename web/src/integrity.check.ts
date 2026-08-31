// Self-check for the damaged-sample flag. Run: node src/integrity.check.ts
//
// The failure this guards against is the quiet one: an array that failed produces an EMPTY result,
// which on the page looks exactly like a clean sample. A flag that said nothing in that case would
// be worse than no flag, because it would add a reassurance the data does not support.
import assert from 'node:assert/strict'
import { assessIntegrity } from './integrity.ts'
import { COMPLEX_DEVIANT_FRACTION, COMPLEX_CALL_RATE } from './abnormalities.ts'

const base = { stage: 'bulk' as const, callRate: 0.98, affectedAutosomes: 0, totalAutosomes: 22, events: 0 }

// --- 1. AN EMPTY RESULT FROM A FAILED ARRAY IS NOT A CLEAN SAMPLE -------------------------------
{
  const clean = assessIntegrity(base)
  const failed = assessIntegrity({ ...base, stage: 'failed' })
  assert.equal(clean.level, 'intact')
  assert.equal(failed.level, 'uninterpretable',
    'a failed array with zero events must NOT read as intact: that is the whole point of the flag')
  assert.notEqual(failed.headline, clean.headline,
    'the two must not present identically, since they are identical in every other field')
  assert.ok(failed.limits.length > 0, 'an uninterpretable sample must say what not to do with it')
  assert.equal(clean.limits, '', 'and a clean one must not manufacture a caveat')
}

// --- 2. THE BOUNDS ARE THE ONES ALREADY MEASURED, exercised at their edges ----------------------
//
// No number in this module is new. If one of these ever fails, it is because a shared constant
// moved, and the flag should move with it rather than drift out of step.
{
  const justUnder = assessIntegrity({ ...base, callRate: COMPLEX_CALL_RATE - 0.001, events: 1,
    affectedAutosomes: 1 })
  const justOver = assessIntegrity({ ...base, callRate: COMPLEX_CALL_RATE + 0.001, events: 1,
    affectedAutosomes: 1 })
  assert.equal(justUnder.level, 'uninterpretable', 'below the measured call-rate bound')
  assert.equal(justOver.level, 'disturbed', 'above it, the sample is usable')

  const atBound = Math.ceil(COMPLEX_DEVIANT_FRACTION * 22)
  const wrecked = assessIntegrity({ ...base, affectedAutosomes: atBound, events: atBound })
  const ok = assessIntegrity({ ...base, affectedAutosomes: atBound - 1, events: atBound - 1 })
  assert.equal(wrecked.level, 'uninterpretable',
    `${atBound} of 22 autosomes is at or over the ${COMPLEX_DEVIANT_FRACTION} deviant-fraction bound`)
  assert.equal(ok.level, 'disturbed', 'one fewer is under it')
  // The genome, not the array, is the thing at fault here, and the wording has to say so or an
  // operator will re-run a perfectly good array.
  assert.ok(/GENOME/.test(wrecked.why), 'must distinguish a disturbed genome from a bad reaction')
  assert.ok(/positions stand|events themselves are real/i.test(wrecked.limits),
    'position needs no genotype, so it survives even here and the limits must say what does not')
}

// --- 3. A FAILED ARRAY OUTRANKS ITS OWN EVENT COUNT ---------------------------------------------
//
// The events on a failed array were found by the same broken measurement, so a low count is not
// evidence of health.
{
  const failedWithEvents = assessIntegrity({ ...base, stage: 'failed', events: 1, affectedAutosomes: 1 })
  assert.equal(failedWithEvents.level, 'uninterpretable')
  assert.ok(/not measuring a genome/.test(failedWithEvents.why))
}

// --- 4. the ordinary case reads as ordinary -----------------------------------------------------
{
  const two = assessIntegrity({ ...base, affectedAutosomes: 2, events: 2 })
  assert.equal(two.level, 'disturbed')
  assert.ok(two.headline.includes('2 changes'))
  assert.equal(two.limits, '', 'a usable sample must not carry a limits paragraph')
  assert.equal(assessIntegrity({ ...base, affectedAutosomes: 1, events: 1 }).headline
    .includes('1 change'), true, 'singular where it should be singular')
}

// --- 5. degenerate inputs do not produce a confident answer -------------------------------------
{
  const noChroms = assessIntegrity({ ...base, totalAutosomes: 0 })
  assert.ok(Number.isNaN(noChroms.measures.deviantFraction),
    'no autosomes means no fraction, and it must not be reported as 0')
  const noCallRate = assessIntegrity({ ...base, callRate: NaN, events: 1, affectedAutosomes: 1 })
  assert.equal(noCallRate.level, 'disturbed',
    'an unknown call rate must not silently trip the damaged branch')
}

// --- 6. every string an operator reads is clean --------------------------------------------------
{
  for (const c of [
    assessIntegrity(base),
    assessIntegrity({ ...base, stage: 'failed' }),
    assessIntegrity({ ...base, callRate: 0.5, events: 1, affectedAutosomes: 1 }),
    assessIntegrity({ ...base, affectedAutosomes: 20, events: 20 }),
    assessIntegrity({ ...base, affectedAutosomes: 2, events: 2 }),
  ]) {
    for (const s of [c.headline, c.why, c.limits]) {
      assert.ok(!s.includes('—'), `em dash in "${s.slice(0, 50)}"`)
    }
    assert.ok(c.headline.length > 0 && c.why.length > 0)
  }
}

console.log('integrity.check.ts: a failed array does not read as a clean one, and every bound is '
  + 'one the taxonomy already measured')

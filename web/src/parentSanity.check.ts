// Self-check for the parental slot check. Run: node src/parentSanity.check.ts
//
// The failure this guards against leaves no other trace. A valid array of a real person, put in the
// wrong slot, produces arithmetic that is entirely correct under a label that is entirely wrong. So
// the properties pinned here are: a decisive disagreement is a conflict, an ARRAY THAT CANNOT
// ANSWER IS NOT, and the direction is never reversed.
import assert from 'node:assert/strict'
import { reconcileParentSex } from './parentSanity.ts'
import { Y_CALL_MIN } from './sexing.ts'
import type { SexCall } from './sexing.ts'

/** The two shapes real donor arrays take, from the measurement in the module header. */
const MALE: SexCall = { yBearing: true, callRatio: 1.034, lrrShift: 0.12 }
const FEMALE: SexCall = { yBearing: false, callRatio: 0, lrrShift: -0.11 }
/** A panel with too few Y probes to ask. Not a finding about the sample. */
const SILENT: SexCall = { yBearing: null, callRatio: NaN, lrrShift: NaN }

// --- 1. THE FOUR ORDINARY COMBINATIONS, and only two of them are conflicts --------------------
{
  assert.equal(reconcileParentSex('paternal', MALE).conflict, undefined,
    'a male array in the paternal slot is what is expected')
  assert.equal(reconcileParentSex('maternal', FEMALE).conflict, undefined,
    'a female array in the maternal slot is what is expected')
  assert.equal(reconcileParentSex('paternal', FEMALE).conflict, 'paternal-array-has-no-y',
    'an egg donor in the paternal slot is THE mistake this module exists for')
  assert.equal(reconcileParentSex('maternal', MALE).conflict, 'maternal-array-has-y',
    'and the sperm donor in the maternal slot is the same mistake mirrored')
}

// --- 2. SILENCE IS NOT A CONFLICT ---------------------------------------------------------------
//
// yBearing is null when the panel carries too few Y probes to read one way or the other. Treating
// that as "no Y" would refuse every run on every panel without Y coverage, which is a large class
// of real arrays, and it would do so while claiming to have measured something.
{
  for (const role of ['paternal', 'maternal'] as const) {
    const r = reconcileParentSex(role, SILENT)
    assert.equal(r.conflict, undefined,
      `a panel that cannot answer must not be treated as an answer (${role})`)
    assert.equal(r.yBearing, null)
    assert.ok(/not checked/i.test(r.headline),
      'and it must say the check did not run rather than implying it passed')
    assert.ok(/Nothing is withheld/i.test(r.why))
  }
}

// --- 3. THE DIRECTION IS NEVER REVERSED ---------------------------------------------------------
//
// Getting this backwards would refuse every correct run and pass every swapped one, which is worse
// than not having the check at all.
{
  const swapped = reconcileParentSex('paternal', FEMALE)
  assert.ok(/PATERNAL ARRAY CARRIES NO CHROMOSOME Y/.test(swapped.headline))
  assert.ok(/egg donor/i.test(swapped.why),
    'the message has to name the likely cause, or an operator cannot act on it')
  const fine = reconcileParentSex('paternal', MALE)
  assert.ok(!/no chromosome y/i.test(fine.headline))
  assert.ok(/matches the slot/i.test(fine.headline))
}

// --- 4. A CONFLICT SAYS WHAT IT COSTS, AND WHAT IT DOES NOT --------------------------------------
//
// The events survive a wrong slot: they were found by channels that never read the parental
// array's sex. Only the naming is corrupted. If the message does not say that, an operator throws
// away a run that still carries real findings.
{
  for (const r of [reconcileParentSex('paternal', FEMALE),
    reconcileParentSex('maternal', MALE)]) {
    assert.ok(/No parent is named/i.test(r.headline), 'the cost has to be in the headline')
    assert.ok(/still real|position/i.test(r.why),
      'and what SURVIVES has to be there too, or a usable run gets discarded')
    assert.ok(/which file is in which slot/i.test(r.why), 'with the action to take')
  }
}

// --- 5. THE MEASUREMENT IS QUOTED, so the claim is checkable by hand ----------------------------
{
  const r = reconcileParentSex('paternal', FEMALE)
  assert.ok(r.why.includes('0.000'), `the measured Y call ratio has to appear: ${r.why}`)
  assert.ok(/1\.03/.test(r.why),
    'and the figure it is being compared against, which is what a sperm donor reads')
  // A conflict on an array whose ratio was never measured must still read sensibly rather than
  // printing NaN at an operator.
  const noNumbers = reconcileParentSex('paternal', { yBearing: false, callRatio: NaN, lrrShift: NaN })
  assert.ok(!/NaN/.test(noNumbers.why), `no NaN may reach an operator: ${noNumbers.why}`)
  assert.ok(/not measured/.test(noNumbers.why))
}

// --- 6. IT AGREES WITH THE THRESHOLD IT RESTS ON -------------------------------------------------
//
// This module does not re-decide what a Y looks like; sexing.ts does. Pinned so the two cannot
// drift apart silently.
{
  assert.ok(Y_CALL_MIN > 0 && Y_CALL_MIN < 1,
    'the separator sits between an egg donor at 0.000 and a sperm donor at 1.03')
  assert.ok(FEMALE.callRatio < Y_CALL_MIN && MALE.callRatio > Y_CALL_MIN,
    'the two real populations must fall either side of it')
}

// --- 7. every string an operator reads is clean --------------------------------------------------
{
  for (const role of ['paternal', 'maternal'] as const) {
    for (const sex of [MALE, FEMALE, SILENT]) {
      const r = reconcileParentSex(role, sex)
      for (const s of [r.headline, r.why]) {
        assert.ok(s.length > 0)
        assert.ok(!s.includes('—'), `em dash in "${s.slice(0, 50)}"`)
      }
    }
  }
}

console.log('parentSanity.check.ts: an egg donor in the paternal slot is a conflict, a panel that '
  + 'cannot answer is not, and a conflict says what it costs and what survives it')

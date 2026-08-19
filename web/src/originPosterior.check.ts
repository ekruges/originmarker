// Self-check for the class-marginal origin posterior. Run: node src/originPosterior.check.ts
//
// One property carries this module and it is the reason the module exists.
//
//   THE SIGN OF THE SHIFT NAMES A PARENT ONLY GIVEN THE CLASS. Gain inverts the map that loss and
//   copy-neutral LOH share. A version that resolved that ambiguity by assuming a class would be
//   confidently wrong on every event of the other kind, which is exactly the bug this replaces:
//   measured wrong-parent rate 0.551 to 0.580 on true low-fraction gains scored as losses.
//
// The tests below are written so that reintroducing that assumption fails loudly rather than
// quietly returning plausible numbers.
import assert from 'node:assert/strict'
import {
  originPosterior, shiftMean, shiftMagnitude, logRMean, fractionAt, applyCalibration,
  bandOf, classInvertedRisk, CLASS_SIGN, CLASSES, BAND_ACCURACY, BAND_A_MIN,
  VETO_MAX_F, bandObligateHet, SATURATING_CHANNEL_NOTE, type CalibrationMap,
} from './originPosterior.ts'

// --- 1. THE INVERSION IS REAL AND IT IS IN THE ALGEBRA, NOT IN A CONSTANT ------------------------
//
// Derived rather than transcribed: the expected shift when the LOADED parent's copy is affected
// must be positive under loss and copy-neutral LOH and NEGATIVE under gain, at every fraction. If
// this ever agrees across all three classes, the whole module is unnecessary and something has
// been flattened.
{
  for (const f of [0.02, 0.05, 0.10, 0.20, 0.40, 0.70]) {
    assert.ok(shiftMean('loss', f, 'loaded') > 0, `loss must read positive at f=${f}`)
    assert.ok(shiftMean('cnn-loh', f, 'loaded') > 0, `cnn-loh must read positive at f=${f}`)
    assert.ok(shiftMean('gain', f, 'loaded') < 0, `GAIN MUST INVERT at f=${f}`)
    // And the two hypotheses must be exact mirrors within a class, or the posterior acquires a
    // directional null of the kind this project has already had to remove once.
    assert.ok(Math.abs(shiftMean('loss', f, 'loaded') + shiftMean('loss', f, 'other')) < 1e-15)
    assert.ok(Math.abs(shiftMean('gain', f, 'loaded') + shiftMean('gain', f, 'other')) < 1e-15)
  }
  assert.equal(CLASS_SIGN.gain, -CLASS_SIGN.loss, 'gain opposes loss by construction')
  assert.equal(CLASS_SIGN['cnn-loh'], CLASS_SIGN.loss, 'copy-neutral LOH shares the loss direction')

  // Rank order of magnitudes: copy-neutral is the LARGEST signal, gain the smallest. This inverts
  // the framing the project was built on and is worth pinning.
  for (const f of [0.05, 0.2, 0.5]) {
    assert.ok(shiftMagnitude('cnn-loh', f) > shiftMagnitude('loss', f))
    assert.ok(shiftMagnitude('loss', f) > shiftMagnitude('gain', f))
  }
}

// --- 2. THE OLD BUG, PINNED AS A REGRESSION -------------------------------------------------------
//
// Take a shift produced by a true GAIN affecting the un-genotyped parent's copy. Under the old
// code that shift was inverted through the loss formula and named the loaded parent. The posterior
// must not do that: given intensity that resolves the class as a gain, it must name the OTHER
// parent, which is the truth.
{
  const f = 0.20
  const trueShift = shiftMean('gain', f, 'other')      // positive, and a loss-assumer reads it as "loaded"
  assert.ok(trueShift > 0, 'a gain on the other parent produces a POSITIVE shift')

  const withIntensity = originPosterior({
    shift: trueShift, shiftSd: 0.004,
    logR: logRMean('gain', f), logRSd: 0.01,        // intensity says gain, clearly
    material: 'bulk', markers: 800,
  })
  assert.equal(withIntensity.classResolved, 'gain', 'intensity must resolve this class')
  assert.equal(withIntensity.parent, 'other',
    'THE REGRESSION: a true gain on the other parent must name the OTHER parent, not the loaded '
    + 'one. Naming "loaded" here is precisely the inversion this module was written to remove')
  assert.ok(withIntensity.confidence > 0.9, 'and with the class resolved it should be confident')
}

// --- 3. WITHOUT THE CLASS, THE POSTERIOR HEDGES RATHER THAN GUESSING ------------------------------
//
// The same shift with no intensity channel must NOT produce a confident call, because the
// observation is equally consistent with a loss on the loaded parent and a gain on the other. An
// implementation that returns a confident answer here has assumed a class.
{
  const f = 0.20
  const trueShift = shiftMean('gain', f, 'other')
  const blind = originPosterior({
    shift: trueShift, shiftSd: 0.004, material: 'trophectoderm', markers: 800,
  })
  assert.equal(blind.classResolved, 'unresolved', 'no intensity means no class')
  assert.ok(blind.confidence < 0.90,
    `hedging across an inverting map must not reach band B, got ${blind.confidence.toFixed(4)}`)
  assert.ok(blind.why.includes('hedges across that inversion'))
  assert.equal(blind.limitedBy, 'class-unresolved',
    'and the reason offered must be the class, since that is the remedy the user can act on')
}

// --- 4. INTENSITY NAMES NO PARENT OF ITS OWN ------------------------------------------------------
//
// The invariant is narrower than "intensity never moves the origin", and getting it wrong once
// while writing this check is worth recording. At a FIXED class and a FIXED fraction the intensity
// term is identical under both hypotheses, so it cancels exactly. But intensity also informs the
// FRACTION, and a better-known fraction sharpens the expected shift magnitude, which legitimately
// sharpens the likelihood ratio between the two parents. That is proper inference, not a leak: the
// channel still never distinguishes the parents, it only says how big the event is.
//
// So the assertion is made with the fraction pinned as well, which is the only place exact
// cancellation is a property of the model.
{
  const shift = 0.02
  const only = { loss: 1, gain: 0, 'cnn-loh': 0 }
  const oneF = { classPrior: only, fGrid: [0.30] }
  const a = originPosterior({ shift, shiftSd: 0.004, material: 'bulk' }, oneF)
  const b = originPosterior(
    { shift, shiftSd: 0.004, logR: -0.08, logRSd: 0.01, material: 'bulk' }, oneF)
  assert.ok(Math.abs(a.pOther - b.pOther) < 1e-12,
    'with class AND fraction pinned, the intensity term must cancel exactly between the two '
    + 'parental hypotheses')

  // And with the fraction free, intensity SHOULD move the number, because it has learned the size
  // of the event. Pinning this stops a future refactor from "fixing" the behaviour above by
  // severing intensity from the fraction entirely.
  //
  // Measured at a deliberately WEAK shift. At the 5-sigma shift used above the posterior is already
  // pegged at 1e-8 and there is nothing left for intensity to sharpen, so testing it there would
  // have asserted a property of saturation rather than of the model.
  const free = { classPrior: only }
  const weak = 0.002
  const c = originPosterior({ shift: weak, shiftSd: 0.004, material: 'bulk' }, free)
  const d = originPosterior(
    { shift: weak, shiftSd: 0.004, logR: logRMean('loss', 0.5), logRSd: 0.01, material: 'bulk' },
    free)
  assert.ok(Math.abs(c.pOther - d.pOther) > 0.1,
    'with the fraction free, intensity must inform it and thereby sharpen the origin')
  assert.equal(c.parent, d.parent, 'but it must not FLIP the parent within one class')

  // Across classes it flips the answer, which is the whole mechanism.
  const asLoss = originPosterior({
    shift, shiftSd: 0.004, logR: logRMean('loss', 0.3), logRSd: 0.01, material: 'bulk' })
  const asGain = originPosterior({
    shift, shiftSd: 0.004, logR: logRMean('gain', 0.3), logRSd: 0.01, material: 'bulk' })
  assert.notEqual(asLoss.parent, asGain.parent,
    'the SAME shift must name different parents under a resolved loss and a resolved gain. If '
    + 'these agree, the inversion has been lost and the module is back to its old bug')
}

// --- 5. THE CALIBRATION MAP IS MONOTONE AND CANNOT CHANGE A CALL ----------------------------------
{
  const map: CalibrationMap = { raw: [0.5, 0.7, 0.9, 1.0], calibrated: [0.5, 0.82, 0.97, 1.0] }
  assert.ok(Math.abs(applyCalibration(map, 0.7) - 0.82) < 1e-12)
  assert.ok(Math.abs(applyCalibration(map, 0.8) - 0.895) < 1e-12, 'linear between knots')
  // Clamped, never extrapolated: a fit says nothing beyond its own support.
  assert.equal(applyCalibration(map, 0.1), 0.5)
  assert.equal(applyCalibration(map, 2.0), 1.0)
  // Monotone in, monotone out, which is what makes it unable to reorder two rows of a table.
  let prev = -Infinity
  for (let x = 0.5; x <= 1.0; x += 0.01) {
    const y = applyCalibration(map, x)
    assert.ok(y >= prev - 1e-12, 'the map must never decrease')
    prev = y
  }
  // And an absent map must be announced rather than silently passing the raw number off as final.
  const noMap = originPosterior({ shift: 0.03, shiftSd: 0.004, material: 'bulk' })
  assert.equal(noMap.uncalibrated, true)
  assert.ok(noMap.why.includes('RAW posterior'))
}

// --- 6. BANDS, AND THE MEASURED ACCURACIES THAT TRAVEL WITH THEM ----------------------------------
{
  assert.equal(bandOf(0.999), 'A')
  assert.equal(bandOf(BAND_A_MIN), 'A')
  assert.equal(bandOf(0.95), 'B')
  assert.equal(bandOf(0.80), 'C')
  assert.equal(bandOf(0.60), 'D')
  // An absent number is NOT band D. D is a measured band with a measured accuracy of about 0.62,
  // and handing an ungraded row that accuracy is exactly the borrowing this grade exists to stop.
  assert.equal(bandOf(NaN), 'F', 'an absent number is ungraded, not the weakest MEASURED band')
  assert.equal(bandOf(0.50), 'F', 'a coin flip is ungraded')
  assert.equal(bandOf(0.56), 'D', 'and just above the floor it is the weakest measured band')
  for (const m of ['bulk', 'esc-single', 'trophectoderm', 'blastomere'] as const) {
    assert.ok(!('F' in BAND_ACCURACY[m]),
      `F must carry no measured accuracy on ${m}: there is no cell to fill, and inventing one `
      + 'would let an ungraded row borrow a number from a graded one')
  }

  // Band D is weak but it is NOT a coin flip, and that is what lets it carry its number. Every
  // measured accuracy must sit clear of 0.5, or the decision to display it stops being honest.
  for (const m of ['bulk', 'esc-single', 'trophectoderm', 'blastomere'] as const) {
    assert.ok(BAND_ACCURACY[m].D > 0.55, `band D on ${m} must beat chance by a real margin`)
    assert.ok(BAND_ACCURACY[m].A > 0.99, `band A on ${m} must justify its label`)
    // Monotone across bands, or the labels mean nothing.
    assert.ok(BAND_ACCURACY[m].A > BAND_ACCURACY[m].B)
    assert.ok(BAND_ACCURACY[m].B > BAND_ACCURACY[m].C)
    assert.ok(BAND_ACCURACY[m].C > BAND_ACCURACY[m].D)
  }
}

// --- 7. THE OBLIGATE-HET CAP IS STRUCTURAL --------------------------------------------------------
//
// Not a power argument. The channel's failure mode is dropout, and dropout is the same event as
// the observation, so no amount of data separates them. It must never reach the top band.
{
  const strong = { shift: 0.14, shiftSd: 0.002, logR: logRMean('loss', 0.6), logRSd: 0.005,
    material: 'bulk' as const, markers: 900 }
  const uncapped = originPosterior(strong)
  assert.equal(uncapped.band, 'A', 'this evidence would otherwise be top band')
  const capped = originPosterior(strong, { obligateHetOnly: true })
  assert.equal(capped.band, 'B', 'resting on the obligate-het channel alone must cap at B')
  assert.equal(capped.confidence, uncapped.confidence,
    'the cap changes the BAND, not the number, so the calibration is not falsified by it')
}

// --- 8. THE ONE VETO, AND ITS SHAPE ---------------------------------------------------------------
//
// Amplified material, possible small gain, intensity unable to resolve a direction. Measured 0/34
// and 0/18 correct on TE and blastomere. No gate on observables separates these, so the parent is
// withheld rather than scored.
{
  assert.equal(classInvertedRisk('trophectoderm', 0.10, undefined), true,
    'amplified, small implied gain, no intensity: this is the cell')
  assert.equal(classInvertedRisk('blastomere', 0.10, 1.0), true,
    'intensity present but under 99% is still an unresolved direction')
  assert.equal(classInvertedRisk('trophectoderm', 0.10, 3.0), false,
    'intensity that resolves the direction removes the risk')
  assert.equal(classInvertedRisk('trophectoderm', 0.40, undefined), false,
    `above f=${VETO_MAX_F} the classes separate and the veto lifts`)
  assert.equal(classInvertedRisk('bulk', 0.10, undefined), false,
    'bulk resolves its own class and was never in this cell')
}

// --- 9. THE FRACTION INVERSION IS STATE-SPECIFIC --------------------------------------------------
{
  const d = 0.04
  assert.ok(Math.abs(fractionAt('loss', d) - (4 * d) / (1 + 2 * d)) < 1e-12)
  assert.ok(Math.abs(fractionAt('gain', d) - (4 * d) / (1 - 2 * d)) < 1e-12)
  assert.ok(Math.abs(fractionAt('cnn-loh', d) - 2 * d) < 1e-12)
  // Same deviation, three different fractions, and the gain reading is the largest.
  assert.ok(fractionAt('gain', d) > fractionAt('loss', d))
  assert.ok(fractionAt('loss', d) > fractionAt('cnn-loh', d))
  assert.ok(Number.isNaN(fractionAt('loss', 0)), 'a zero deviation implies no fraction')
}

// --- 10. A DEGENERATE INPUT REFUSES RATHER THAN INVENTING -----------------------------------------
{
  const noScale = originPosterior({ shift: 0.02, shiftSd: 0, material: 'bulk' })
  assert.equal(noScale.parent, 'withheld')
  assert.ok(Number.isNaN(noScale.confidence))
  // A withheld parent with no confidence is ungraded, which is now F rather than the weakest
  // measured band. The two must not share a symbol: one was measured at 0.62 and this was not
  // measured at all.
  assert.equal(bandOf(noScale.confidence), 'F')
  assert.ok(CLASSES.length === 3)
}

// --- 11. THE SATURATING CHANNELS ARE CAPPED, AND THE CAP IS NOT COSMETIC -------------------------
//
// Measured on the shipped genotype channels: the obligate-het posterior returns 1.0000 at 0, 2, 8,
// 40, 120 and 300 exclusive markers out of 400, flipping its verdict while never moving its number.
// A number that does not vary is not a confidence. Both genotype channels are therefore capped
// below the top band, and the note below is what the display uses to say so in words.
{
  assert.equal(bandObligateHet(1.0), 'B', 'a saturated posterior must not reach band A')
  assert.equal(bandObligateHet(0.999), 'B')
  assert.equal(bandObligateHet(0.95), 'B', 'below the cap it is unchanged')
  assert.equal(bandObligateHet(0.80), 'C')
  assert.equal(bandObligateHet(0.60), 'D')
  assert.ok(SATURATING_CHANNEL_NOTE.includes('band is the meaningful output'))
  assert.ok(SATURATING_CHANNEL_NOTE.includes('uncalibrated'))

  // The dosage posterior must NOT be capped: it is the one channel whose bands were measured.
  const strong = originPosterior({ shift: 0.14, shiftSd: 0.002, logR: logRMean('loss', 0.6),
    logRSd: 0.005, material: 'bulk', markers: 900 })
  assert.equal(strong.band, 'A',
    'the measured channel keeps its top band; capping it would discard the one calibration that '
    + 'was actually earned')
}

// --- 12. THE TWO ARMS ARE EXACT MIRRORS, WHICH IS WHAT KEEPS CALL-M-BY-ABSENCE HONEST ------------
//
// A one-parent run names a maternal event because the father's contribution is intact while a
// change is present. An external review measured that arm as far weaker than the paternal one on
// amplified material, a maternal gain on trophectoderm correct 0.402 of the time against 0.727,
// and its mechanism was that the loaded parent's markers ANCHOR THE FRAME so the arm needing
// displacement away from the anchor loses.
//
// This project removed that anchor before the review was written: self-referencing each chromosome
// against the array's own genome brought a one-parent null of -0.031 to -0.001..+0.006. Measured on
// 210,000 injections here the gap is 1.3 points and runs the other way, so no arm correction ships.
// These assertions are what stop the anchoring returning unnoticed, which would reintroduce a bias
// against exactly the calls a one-parent cohort depends on.
{
  for (const c of CLASSES) {
    for (let f = 0.02; f <= 0.70; f += 0.02) {
      assert.equal(shiftMean(c, f, 'loaded') + shiftMean(c, f, 'other'), 0,
        `the arms must be exact mirrors: ${c} at f=${f.toFixed(2)} is anchored`)
    }
  }
  // And a mirrored OBSERVATION must give a mirrored posterior, which is the property a reader of a
  // one-parent run is relying on without knowing it.
  for (const d of [0.004, 0.02, 0.08, 0.2]) {
    const a = originPosterior({ shift: +d, shiftSd: 0.01, material: 'trophectoderm', markers: 800 })
    const b = originPosterior({ shift: -d, shiftSd: 0.01, material: 'trophectoderm', markers: 800 })
    assert.ok(Math.abs((a.pOther + b.pOther) - 1) < 1e-9,
      `mirrored observations must give mirrored posteriors at d=${d}: got ${a.pOther + b.pOther}`)
  }
}

console.log('originPosterior.check.ts: all assertions passed, including the regression that a true '
  + 'gain on the un-genotyped parent is no longer inverted into a confident call for the loaded one')

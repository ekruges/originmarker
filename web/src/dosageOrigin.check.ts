// Self-check for the dosage channel. Run: node src/dosageOrigin.check.ts
//
// This module was rewritten after a methods review measured its predecessor's null on real
// material. The assertions below are organised around what that review found, because each finding
// is a way this statistic can look right and be wrong.
//
//   SELF-REFERENCING IS THE WHOLE THING. The raw one-parent null sits at -0.031 on trophectoderm
//   and -0.023 on blastomeres under NO event, pointing at the parent that was not genotyped. The
//   old version asserted symmetry between the two DIRECTIONS of the statistic, which was true and
//   which did not help: a symmetric statistic with an off-centre null is a confidently wrong answer
//   in one direction. So the test here is that a biased array produces no call.
//
//   THE ORDER OF THE QUESTIONS IS LOAD-BEARING. Array quality, then whether anything of this kind
//   could answer at this width, then whether there is an imbalance, then whether the sign is
//   secure. Each earlier question must short-circuit the later ones, or a bad array gets a
//   confident parent.
import assert from 'node:assert/strict'
import {
  callDosageOrigin, centroid, oriented, fractionFromShift, logitShift, materialOf,
  WINDOW_LO, WINDOW_HI, SIGN_SECURE_F, MAX_HET_BAF_SD, F80_CHROMOSOME, F80_SEGMENT,
  DRIFT_TAU, VIF_CHROMOSOME, RESIDUAL_R, floorFor, ORIGIN_WITHOUT_CLASS,
} from './dosageOrigin.ts'
import type { DosageState } from './dosageOrigin.ts'
import type { AB } from './informativity.ts'

const N = 4_000

/**
 * A region of `n` markers whose oriented dosage is centred on `mu`, with deterministic spread so
 * the tests do not depend on a random seed. Half the markers are parent-AA and half parent-BB with
 * their dosages mirrored, which exercises the orientation on every call.
 */
function region(n: number, mu: number, spread = 0.10): [AB, number][] {
  const out: [AB, number][] = []
  for (let i = 0; i < n; i += 1) {
    const jitter = spread * Math.sin(i * 2.399963) // deterministic, mean ~0 over many i
    const b = Math.min(0.999, Math.max(0.001, mu + jitter))
    out.push(i % 2 === 0 ? ['AA', b] : ['BB', 1 - b])
  }
  return out
}

// --- 1. the corrected expectation and its inversion ----------------------------------------------
//
// E[BAF] = 1/(2-f), not 0.5 + 0.5f. The old form is the mean of per-cell frequencies and assumes a
// monosomic cell contributes as much DNA as a disomic one; an array reads pooled dosage.
{
  const pooled = (f: number) => 1 / (2 - f)
  for (const f of [0.05, 0.10, 0.20, 0.30, 0.50, 0.70]) {
    const d = pooled(f) - 0.5
    assert.ok(Math.abs(fractionFromShift(d) - f) < 1e-12,
      `f = 4d/(1+2d) must invert the pooled expectation, failed at ${f}`)
    assert.ok(0.5 * f > d, `0.5+0.5f must overstate the shift at f=${f}`)
  }
  const at = (f: number) => (0.5 * f) / (1 / (2 - f) - 0.5)
  assert.ok(Math.abs(at(0.05) - 1.95) < 0.01, 'overstated 1.95x at low f')
  assert.ok(Math.abs(at(1.00) - 1.00) < 1e-9, 'and correct only at complete loss')
  // The 0.65 band edge is crossed at f = 0.46, not 0.30.
  assert.ok(pooled(0.30) < 0.65 && pooled(0.47) > 0.65)
  // Logit shift is additive and independent of baseline.
  assert.ok(Math.abs(logitShift(0.5) - Math.log(2)) < 1e-12)
  assert.ok(logitShift(0) === 0, 'no fraction is no shift (negative zero is still zero)')
  assert.ok(Number.isNaN(fractionFromShift(0)), 'a non-positive shift implies no fraction')
}

// --- 1b. THE INVERSION IS STATE-SPECIFIC AND THE THREE MUST NOT BE SHARED -------------------------
//
// Measured by brute-force copy counting to 1.1e-16. An earlier version applied the LOSS inversion
// to every state: right for a loss, and wrong in the OPPOSITE direction for a gain, which is the
// worse way to be wrong because the error inflates the statistic rather than shrinking it.
{
  const at = (state: DosageState) => fractionFromShift(0.04, state)
  assert.ok(Math.abs(at('loss') - 0.148) < 0.001, 'a loss inverts as 4d/(1+2d)')
  assert.ok(Math.abs(at('gain') - 0.174) < 0.001, 'a gain as 4d/(1-2d), which is LARGER not smaller')
  assert.ok(Math.abs(at('cnn-loh') - 0.080) < 0.001, 'and copy-neutral LOH as 2d')
  assert.ok(at('gain') > at('loss'), 'the two loss-like inversions diverge, so sharing one is an error')

  // Rank order by DEVIATION at fixed f is CNN-LOH > loss > gain, which is the opposite of the
  // intuition that a state with no copy-number signal must be hardest. Read through the
  // inversions: at the same deviation, CNN-LOH implies the SMALLEST fraction, i.e. the largest
  // signal per unit of f.
  assert.ok(at('cnn-loh') < at('loss') && at('loss') < at('gain'),
    'copy-neutral LOH is the largest-signal state, not the hardest')

  // "A gain shifts about a third as much as a loss" is the f = 1 endpoint only. The pooled ratio
  // is (2-f)/(2+f), so in the mosaic range the penalty is 1.1 to 1.7x.
  const ratio = (f: number) => (2 - f) / (2 + f)
  assert.ok(Math.abs(ratio(0.1) - 0.905) < 0.001)
  assert.ok(Math.abs(ratio(0.2) - 0.818) < 0.001)
  assert.ok(Math.abs(ratio(1.0) - 0.333) < 0.001, 'a third is where that claim comes from')

  // A gain's inversion has its own ceiling: no gain fraction reaches a deviation of 0.5.
  assert.ok(Number.isNaN(fractionFromShift(0.5, 'gain')))
  assert.ok(Number.isFinite(fractionFromShift(0.5, 'loss')))
}

// --- 2. orientation and the central window --------------------------------------------------------
{
  assert.equal(oriented('AA', 0.9), 0.9)
  assert.ok(Math.abs(oriented('BB', 0.9) - 0.1) < 1e-12)
  const c = centroid([['AA', 0.5], ['BB', 0.5], ['AA', 0.02], ['BB', 0.98],
    ['AB', 0.5], ['AA', null]])
  assert.equal(c.n, 2, 'only the central window counts toward the centroid')
  assert.equal(c.seen, 4, 'but every readable parent-homozygous marker is counted as seen')
  assert.ok(Math.abs(c.mean - 0.5) < 1e-12, 'and a mirrored pair orients to the same value')
  assert.ok(WINDOW_LO < 0.35 && WINDOW_HI > 0.65,
    'the window must be wider than the old middle band, or it selects for exactly the markers a '
    + 'shift has already moved out of')
}

// --- 3. SELF-REFERENCING REMOVES A DIRECTIONAL NULL ------------------------------------------------
//
// The failure this rewrite exists for. An array whose whole genome sits off-centre by the measured
// trophectoderm bias must produce NO call, because the region is not different from the rest of
// that array. Before self-referencing this returned a confident parent.
{
  const bias = -0.031
  const c = callDosageOrigin(region(N, 0.5 + bias), region(N, 0.5 + bias), 'bulk',
    { wholeChromosome: true })
  assert.equal(c.verdict, 'no-imbalance',
    'a uniformly biased array must produce no call, since nothing distinguishes the region')
  assert.ok(Math.abs(c.shift) < 1e-9)

  // On record: that bias is the shift a real fraction of 0.117 would produce, which is the size of
  // what self-referencing removes.
  assert.ok(Math.abs(fractionFromShift(Math.abs(bias)) - 0.117) < 0.002)
}

// --- 4. THE QUESTIONS ARE ASKED IN ORDER ----------------------------------------------------------
{
  const big = region(N, 0.62)
  const bg = region(N, 0.50)

  // A bad array short-circuits everything, including an interval that would otherwise call.
  const excluded = callDosageOrigin(big, bg, 'bulk', {
    wholeChromosome: true, hetBafSd: MAX_HET_BAF_SD + 0.01,
  })
  assert.equal(excluded.verdict, 'array-excluded')
  assert.ok(excluded.why.includes('independently of this interval'))

  // Not-evaluable is decided before the data is looked at, so a huge shift cannot override it.
  // Only where NO floor exists: a blastomere loss against one parent, at any fraction to 0.70.
  {
    const c = callDosageOrigin(big, bg, 'blastomere', { wholeChromosome: true, state: 'loss' })
    assert.equal(c.verdict, 'not-evaluable', 'a blastomere loss with one parent has no floor')
    assert.ok(Number.isNaN(c.z), 'and no statistic is computed for it')
    assert.ok(c.why.includes('no array of this kind'))
  }
  // Trophectoderm loss DOES have a floor, 0.625. Too high to be sign-secure, so it reports an
  // unassigned imbalance rather than refusing: strictly more informative than not-evaluable, and
  // the distinction is the point of separating the two questions.
  {
    const c = callDosageOrigin(big, bg, 'trophectoderm', { wholeChromosome: true, state: 'loss' })
    assert.notEqual(c.verdict, 'not-evaluable', 'a TE loss is evaluable, just not sign-secure')
    assert.ok(Number.isFinite(c.floor))
  }

  // Bulk is evaluable at both widths, which is what makes the above a measurement not a refusal.
  assert.notEqual(callDosageOrigin(big, bg, 'bulk', { wholeChromosome: false }).verdict,
    'not-evaluable')
  assert.ok(Number.isFinite(F80_CHROMOSOME.bulk) && Number.isFinite(F80_SEGMENT.bulk))
  assert.ok(Number.isNaN(F80_CHROMOSOME.blastomere) && Number.isNaN(F80_SEGMENT.trophectoderm))
}

// --- 4b. THE FLOOR IS STATE-AWARE, AND COPY-NEUTRAL LOH IS THE EASIEST STATE ----------------------
//
// The module shipped with loss-only floors and therefore refused the class it handles BEST. On a
// trophectoderm biopsy with one parent a loss floor is 0.625 and a copy-neutral one is 0.186,
// because copy number stays at 2 and the genotype caller never degrades.
{
  assert.ok(floorFor('trophectoderm', 'cnn-loh', true, 1) < 0.35,
    'a copy-neutral event on TE with ONE parent is callable, and was being refused')
  assert.ok(floorFor('trophectoderm', 'loss', true, 1) > 0.35,
    'while a loss on the same material is not')
  assert.ok(floorFor('trophectoderm', 'cnn-loh', true, 1)
    < floorFor('trophectoderm', 'loss', true, 1),
  'copy-neutral is the LARGER-signal state at every material')

  // A second parent is worth more than any file of the same kind: 0.625 -> 0.186 on TE.
  const one = floorFor('trophectoderm', 'loss', true, 1)
  const two = floorFor('trophectoderm', 'loss', true, 2)
  assert.ok(one / two > 3, `the second parent must buy over 3x, got ${(one / two).toFixed(2)}`)

  // And the one combination no input moves stays unreachable in both.
  assert.ok(Number.isNaN(floorFor('blastomere', 'loss', true, 1)))
  assert.ok(Number.isNaN(floorFor('blastomere', 'loss', false, 2)),
    'a blastomere segmental loss is a material limit, not an input one')
}

// --- 4c. THE CLASS IS A SEPARATE QUESTION AND MAY FAIL WITHOUT TAKING THE ORIGIN DOWN -------------
//
// The largest single source of refusals before this. On TE and blastomere material every detected
// event resolves an origin and none resolves a class at 400 markers, so a caller emitting one
// verdict had to refuse them all, discarding an origin it could support.
{
  const bg2 = region(N, 0.50)
  const strong = region(N, 0.5 + 0.20)

  // Wide log2R, which is every amplified material: origin named, class withheld.
  const wide = callDosageOrigin(strong, bg2, 'bulk',
    { wholeChromosome: true, state: 'loss', windowLogRSd: 0.20 })
  assert.equal(wide.verdict, 'known-parent-lost', 'the origin still resolves')
  assert.equal(wide.classVerdict, 'unresolved', 'and the class does not')
  assert.ok(wide.classWhy.includes('too wide'))

  // Narrow log2R, which only bulk reaches: both resolve.
  const narrow = callDosageOrigin(strong, bg2, 'bulk',
    { wholeChromosome: true, state: 'loss', windowLogRSd: 0.02 })
  assert.equal(narrow.verdict, 'known-parent-lost')
  assert.equal(narrow.classVerdict, 'loss')

  // Omitting it withholds the class rather than guessing one.
  const none = callDosageOrigin(strong, bg2, 'bulk', { wholeChromosome: true, state: 'loss' })
  assert.equal(none.classVerdict, 'unresolved')
  assert.ok(none.classWhy.includes('no window log2R spread'))

  // The measured share is carried in the message, so a reader sees it is the norm, not a defect.
  assert.equal(ORIGIN_WITHOUT_CLASS.trophectoderm, 1)
  assert.ok(ORIGIN_WITHOUT_CLASS.bulk < 0.5)
}

// --- 5. AN IMBALANCE IS REPORTED WITHOUT A PARENT WHEN THE SIGN IS NOT SECURE ----------------------
//
// Below f = 0.3 between half and all detections name the wrong parent, so the class is withheld.
// The MoChA precedent: report the event, leave it unassigned.
{
  const bg = region(N, 0.50)
  const small = callDosageOrigin(region(N, 0.5 + 0.024), bg, 'bulk', { wholeChromosome: true })
  assert.equal(small.verdict, 'imbalance-unassigned')
  assert.ok(small.impliedF < SIGN_SECURE_F)
  assert.ok(Math.abs(small.z) >= 2.576, 'the imbalance itself is still detected and reported')
  assert.ok(small.why.includes('WRONG parent'))

  // Above the bound the parent is named, and the direction follows the sign of the shift.
  const up = callDosageOrigin(region(N, 0.5 + 0.20), bg, 'bulk', { wholeChromosome: true })
  assert.equal(up.verdict, 'known-parent-lost', 'a positive shift is the loaded parent falling short')
  assert.ok(up.impliedF >= SIGN_SECURE_F)

  const down = callDosageOrigin(region(N, 0.5 - 0.20), bg, 'bulk', { wholeChromosome: true })
  assert.equal(down.verdict, 'other-parent-lost')

  // The two directions are reached on equal evidence. Necessary but, as the review showed, not
  // sufficient on its own: section 3 is what covers the null being off-centre.
  assert.ok(Math.abs(Math.abs(up.z) - Math.abs(down.z)) < 0.5)
  assert.ok(Math.abs(up.impliedF - down.impliedF) < 1e-9)
}

// --- 5b. THE JOINT TERM USES THE MEASURED CHANNEL CORRELATION ------------------------------------
//
// Intensity informs the STATE and never the ORIGIN: on haploid pronuclei, a complete-loss
// experiment with known parent, log2R cannot tell maternal from paternal (p = 0.54) while oriented
// dosage separates them at p = 1.3e-15. So intensity enters the detection step only, and the
// direction still comes from the dosage sign.
//
// The correlation is measured on THIS statistic, over 81 arrays: -0.058 on bulk but +0.49 to +0.63
// on amplified material, because a poorly amplified chromosome reads low in both channels.
// Treating them as independent would overstate the joint evidence by sqrt(2+2r)/sqrt(2).
{
  const bg = region(N, 0.50)
  // Just under the bar on dosage alone; supporting intensity carries it over.
  const sub = region(N, 0.5 + 0.0030)
  const alone = callDosageOrigin(sub, bg, 'bulk', { wholeChromosome: true })
  assert.equal(alone.verdict, 'no-imbalance')
  const joint = callDosageOrigin(sub, bg, 'bulk', { wholeChromosome: true, intensityZ: -6 })
  assert.notEqual(joint.verdict, 'no-imbalance', 'intensity must be able to carry a detection')

  // Intensity pointing the WRONG way (a gain) adds nothing: the term is one-sided on reduction.
  const wrongWay = callDosageOrigin(sub, bg, 'bulk', { wholeChromosome: true, intensityZ: +6 })
  assert.equal(wrongWay.verdict, 'no-imbalance')

  // Intensity NEVER changes the direction, only whether an event is declared.
  const up = callDosageOrigin(region(N, 0.5 + 0.20), bg, 'bulk',
    { wholeChromosome: true, intensityZ: -8 })
  assert.equal(up.verdict, 'known-parent-lost')
  const down = callDosageOrigin(region(N, 0.5 - 0.20), bg, 'bulk',
    { wholeChromosome: true, intensityZ: -8 })
  assert.equal(down.verdict, 'other-parent-lost',
    'the same intensity evidence must not push both directions to the same parent')

  // The correlation is applied, and it is material-specific. A correlated material must need MORE
  // joint evidence than an uncorrelated one for the same pair of channel z values.
  assert.ok(RESIDUAL_R.trophectoderm > 0.5 && Math.abs(RESIDUAL_R.bulk) < 0.1,
    'measured: amplified material is strongly correlated, bulk is not')
  const denom = (m: keyof typeof RESIDUAL_R) => Math.sqrt(2 + 2 * RESIDUAL_R[m])
  assert.ok(denom('trophectoderm') / denom('bulk') > 1.2,
    'ignoring the correlation would overstate the joint z by about a quarter on TE')
}

// --- 6. no event is no call ------------------------------------------------------------------------
{
  const c = callDosageOrigin(region(N, 0.50), region(N, 0.50), 'bulk', { wholeChromosome: true })
  assert.equal(c.verdict, 'no-imbalance')
  assert.ok(Math.abs(c.shift) < 1e-9)
}

// --- 7. DRIFT IS THE FLOOR ON THE STANDARD ERROR, NOT SAMPLING -------------------------------------
//
// The reason a z of 15 arrived on a verified-diploid chromosome. Adding markers must not drive the
// standard error to zero, because within-array drift does not average down.
{
  const shift = 0.02
  const small = callDosageOrigin(region(400, 0.5 + shift), region(N, 0.5), 'bulk',
    { wholeChromosome: true })
  const large = callDosageOrigin(region(40_000, 0.5 + shift), region(N, 0.5), 'bulk',
    { wholeChromosome: true })
  const ratio = Math.abs(large.z) / Math.abs(small.z)
  assert.ok(ratio < 10,
    `a hundredfold in markers must not buy tenfold in z, drift floors it; got ${ratio.toFixed(1)}`)
  assert.ok(DRIFT_TAU.blastomere > DRIFT_TAU.bulk * 10,
    'amplified material drifts far more than bulk, which is why the floors differ by material')
  assert.ok(VIF_CHROMOSOME.blastomere > 3 && VIF_CHROMOSOME.bulk < 1.1,
    'and its readings are spatially correlated where bulk is essentially white')
}

// --- 8. degenerate input is not-evaluable rather than a call ---------------------------------------
{
  assert.equal(callDosageOrigin([], region(N, 0.5), 'bulk', { wholeChromosome: true }).verdict,
    'not-evaluable')
  assert.equal(callDosageOrigin(region(N, 0.5), [], 'bulk', { wholeChromosome: true }).verdict,
    'not-evaluable')
  // A parent heterozygous everywhere carries nothing: every hypothesis predicts the same dosage.
  const het: [AB, number][] = Array.from({ length: N }, () => ['AB', 0.5])
  assert.equal(callDosageOrigin(het, region(N, 0.5), 'bulk', { wholeChromosome: true }).markers, 0)
}

// --- 9. material mapping is conservative for anything unrecognised ---------------------------------
{
  assert.equal(materialOf('bulk'), 'bulk')
  assert.equal(materialOf('trophectoderm'), 'trophectoderm')
  assert.equal(materialOf('blastomere'), 'blastomere')
  for (const s of ['single-cell', 'haploid', 'failed', 'unknown', 'nonsense']) {
    assert.equal(materialOf(s), 'esc-single',
      'an unrecognised stage must take amplified constants rather than bulk ones')
  }
}

console.log('dosageOrigin.check.ts: all assertions passed, including a uniformly biased array '
  + 'producing no call, not-evaluable decided before the data, and the class withheld under the '
  + 'sign-security bound')

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
  DRIFT_TAU, VIF_CHROMOSOME, RESIDUAL_R,
} from './dosageOrigin.ts'
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
  for (const m of ['trophectoderm', 'blastomere'] as const) {
    const c = callDosageOrigin(big, bg, m, { wholeChromosome: true })
    assert.equal(c.verdict, 'not-evaluable', `${m} has no whole-chromosome floor with one parent`)
    assert.ok(Number.isNaN(c.z), 'and no statistic is computed for it')
    assert.ok(c.why.includes('no array of this kind'))
  }
  // 12 Mb intervals are out of reach on every amplified material.
  for (const m of ['trophectoderm', 'blastomere', 'esc-single'] as const) {
    assert.equal(callDosageOrigin(big, bg, m, { wholeChromosome: false }).verdict, 'not-evaluable')
  }
  // Bulk is evaluable at both widths, which is what makes the above a measurement not a refusal.
  assert.notEqual(callDosageOrigin(big, bg, 'bulk', { wholeChromosome: false }).verdict,
    'not-evaluable')
  assert.ok(Number.isFinite(F80_CHROMOSOME.bulk) && Number.isFinite(F80_SEGMENT.bulk))
  assert.ok(Number.isNaN(F80_CHROMOSOME.blastomere) && Number.isNaN(F80_SEGMENT.trophectoderm))
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

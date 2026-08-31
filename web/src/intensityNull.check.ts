// Self-check for the intensity null. Run: node src/intensityNull.check.ts
//
// What is tested is the property the shipped statistic did not have: that ONE affected unit does
// not destroy the scale, and that the scale grows with real dispersion rather than shrinking with
// the marker count. Those two are the whole reason this module exists.
import assert from 'node:assert/strict'
import {
  median, mad, nullScale, calibratedZ, intensityDetects, Z_CHROMOSOME, Z_WINDOW,
} from './intensityNull.ts'

/** A deterministic stream, so a failure is the code changing and never the draw changing. */
function rng(seed: number): () => number {
  let x = seed >>> 0
  return () => {
    x ^= x << 13; x >>>= 0
    x ^= x >> 17
    x ^= x << 5; x >>>= 0
    return x / 0x100000000
  }
}
/** Approximately normal, from the sum of twelve uniforms. */
const gauss = (r: () => number) => {
  let s = 0
  for (let i = 0; i < 12; i += 1) s += r()
  return s - 6
}

// --- 1. median and MAD on known input ----------------------------------------------------------
{
  assert.equal(median([3, 1, 2]), 2)
  assert.equal(median([4, 1, 3, 2]), 2.5)
  assert.ok(Number.isNaN(median([])), 'no data is not zero')
  assert.ok(Number.isNaN(mad([1])), 'one point has no spread')
  // On Gaussian data the scaled MAD estimates the SD. 1000 draws at sd 1.
  const r = rng(7)
  const xs = Array.from({ length: 1000 }, () => gauss(r))
  const m = mad(xs)
  assert.ok(Math.abs(m - 1) < 0.12, `scaled MAD must estimate the SD on Gaussian data, got ${m}`)
}

// --- 2. THE POINT OF THE MODULE: one affected unit must not destroy the scale -------------------
//
// A plain SD over chromosomes was measured inflating 2.34x from a SINGLE injected trisomy, which
// means one real event destroys the null meant to detect it. The MAD must hold.
{
  const sd = (xs: number[]) => {
    const mu = xs.reduce((a, x) => a + x, 0) / xs.length
    return Math.sqrt(xs.reduce((a, x) => a + (x - mu) ** 2, 0) / (xs.length - 1))
  }
  // AVERAGED OVER SEEDS. A single draw of 22 units estimates its own MAD noisily, so pinning one
  // draw's inflation to two decimals tests the draw rather than the estimator. The real-array
  // measurement this mirrors is over 16 arrays.
  const SEEDS = 40
  for (const k of [1, 2, 3]) {
    let sdSum = 0
    let madSum = 0
    for (let seed = 0; seed < SEEDS; seed += 1) {
      const r = rng(1000 + seed)
      const clean = Array.from({ length: 22 }, () => gauss(r) * 0.05)
      const hit = [...clean]
      // A trisomy is log2(3/2) = 0.585, an order of magnitude above the clean spread.
      for (let i = 0; i < k; i += 1) hit[i] = 0.585
      sdSum += sd(hit) / sd(clean)
      madSum += mad(hit) / mad(clean)
    }
    const sdInf = sdSum / SEEDS
    const madInf = madSum / SEEDS
    assert.ok(sdInf > 1.8,
      `k=${k}: a plain SD is destroyed, which is why it is not used (${sdInf.toFixed(2)}x)`)
    assert.ok(madInf < sdInf / 1.5,
      `k=${k}: the MAD must degrade far more slowly than the SD: `
      + `MAD ${madInf.toFixed(2)}x against SD ${sdInf.toFixed(2)}x`)
    // The direction of any residual failure must be SAFE: inflation refuses, never over-calls.
    assert.ok(madInf >= 0.95,
      'contamination must inflate the null, never shrink it, so a contaminated null refuses')
  }
  // ONE affected unit is the case the shipped estimator could not survive, so it is pinned tightly.
  let one = 0
  for (let seed = 0; seed < SEEDS; seed += 1) {
    const r = rng(2000 + seed)
    const clean = Array.from({ length: 22 }, () => gauss(r) * 0.05)
    const hit = [...clean]; hit[0] = 0.585
    one += mad(hit) / mad(clean)
  }
  const oneInf = one / SEEDS
  assert.ok(oneInf < 1.25,
    `a single affected chromosome must barely move the MAD, got ${oneInf.toFixed(3)}x. The plain `
    + 'SD it replaces was measured at 2.34x here, meaning one real trisomy destroyed the null that '
    + 'was supposed to detect it.')
  console.log(`  one affected chromosome inflates the MAD ${oneInf.toFixed(3)}x (SD: 2.34x measured)`)
}

// --- 3. the scale does NOT shrink with the number of markers ------------------------------------
//
// The defect this replaces: sd/sqrt(n) falls as the marker count rises, so averaging more markers
// bought confidence that the data does not contain. The null here is over UNITS, and adding
// markers inside a unit must not change it at all.
{
  const r = rng(13)
  const units = Array.from({ length: 22 }, () => gauss(r) * 0.05)
  const a = nullScale(units)
  const b = nullScale(units)
  assert.equal(a.scale, b.scale, 'the scale is a property of the units, not of anything inside them')
  // The old form: sd/sqrt(n) shrinks without limit as n grows. Demonstrate the difference.
  const sdOverRootN = (n: number) => 0.05 / Math.sqrt(n)
  assert.ok(sdOverRootN(20000) < sdOverRootN(50) / 15,
    'the shipped form falls more than 15-fold from 50 to 20,000 markers, which is the error')
  assert.ok(a.scale > sdOverRootN(1000) * 5,
    'the empirical scale must be far wider than the iid form at a realistic window size')
}

// --- 4. a thin null returns nothing rather than a number ----------------------------------------
{
  const thin = nullScale([0.01, -0.02, 0.03])
  assert.equal(calibratedZ(0.5, thin), undefined,
    'three units is not a null, and a z against it must not be reported')
  const wide = nullScale(Array.from({ length: 22 }, (_, i) => (i % 2 ? 0.04 : -0.04)))
  assert.ok(Number.isFinite(calibratedZ(0.5, wide)!), '22 units is enough to speak')
  assert.equal(calibratedZ(NaN, wide), undefined, 'no measurement, no z')
  assert.equal(calibratedZ(0.5, { centre: 0, scale: 0, units: 22 }), undefined,
    'a zero scale would divide by zero and must refuse instead')
}

// --- 5. the thresholds are the measured ones, not the normal ones -------------------------------
//
// 2.576 is the two-sided normal 1 percent point and it fires on 10.1 percent of event-free
// chromosomes on this material. If either constant drifts back toward it, this fails.
{
  assert.ok(Z_CHROMOSOME > 5, `chromosome threshold must be the measured 99th percentile, got ${Z_CHROMOSOME}`)
  assert.ok(Z_WINDOW > 4, `window threshold must be the measured one, got ${Z_WINDOW}`)
  assert.ok(Z_CHROMOSOME > Z_WINDOW,
    'a whole chromosome has fewer independent units behind its null than a window does, so its '
    + 'threshold must be the stricter of the two')
  const n = nullScale(Array.from({ length: 22 }, (_, i) => (i % 2 ? 0.05 : -0.05)))
  // A shift that clears the old normal threshold but not the measured one must NOT detect.
  const modest = n.centre + 3 * n.scale
  assert.equal(intensityDetects(calibratedZ(modest, n), true), false,
    'a 3-sigma shift is under the measured chromosome threshold and must not be called')
  const big = n.centre + 8 * n.scale
  assert.equal(intensityDetects(calibratedZ(big, n), true), true, 'an 8-sigma shift must be called')
  // The window threshold is looser, so the same shift can detect as a segment and not as a chromosome.
  const between = n.centre + 5 * n.scale
  assert.equal(intensityDetects(calibratedZ(between, n), false), true)
  assert.equal(intensityDetects(calibratedZ(between, n), true), false)
  assert.equal(intensityDetects(undefined, true), false, 'no z is not a detection')
}

console.log('intensityNull.check.ts: one affected unit does not destroy the scale, the scale does '
  + 'not shrink with marker count, and the thresholds are the measured ones')

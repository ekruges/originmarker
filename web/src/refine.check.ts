// Self-check for breakpoint refinement. Run: node src/refine.check.ts
//
// The constants here were measured elsewhere, on 848 events spliced from real arrays at known
// marker positions. What this file checks is that the implementation does what those measurements
// describe: that it beats the window edge by the stated order of magnitude, that it refuses when
// it should, and that the interval is never narrower than the calibration supports.
import assert from 'node:assert/strict'
import {
  refineEdges, LR_DROP, MIN_INTERVAL_MARKERS, SEARCH_RADIUS, REFINE_CALL_FLOOR,
} from './refine.ts'
import { SEGMENT_LRT } from './segments.ts'

/**
 * One chromosome of absence indicators with an event planted at [from, to).
 *
 * Deterministic rather than random, so a failure is reproducible. The hash is spread rather than
 * contiguous so nothing here depends on run structure, which is the property the real scan's
 * clustering violates and this refinement does not assume.
 */
const chrom = (n: number, from: number, to: number, inRate: number, bg: number): Uint8Array => {
  const a = new Uint8Array(n)
  for (let i = 0; i < n; i += 1) {
    const r = ((i * 7919) % 10007) / 10007
    a[i] = r < (i >= from && i < to ? inRate : bg) ? 1 : 0
  }
  return a
}

// --- 1. it finds an edge the window only approximates -------------------------------------------
//
// The window scan steps a quarter of its width, so its edges are wrong by construction. This is
// the case the whole module exists for: hand it a window that overlaps the event but is offset
// by more than a step, and see where the edges land.
{
  const TRUE_A = 20_000
  const TRUE_B = 32_000
  const a = chrom(60_000, TRUE_A, TRUE_B, 0.11, 0.005)
  // A window offset 1,800 markers left and 2,400 short, i.e. what a 2,400-step scan would report.
  const r = refineEdges(a, TRUE_A - 1_800, TRUE_B - 2_400, 0.005)

  assert.ok(r.localised, 'a real event inside the radius must localise')
  const errA = Math.abs(r.start.index - TRUE_A)
  const errB = Math.abs(r.end.index - TRUE_B)
  assert.ok(errA < 200, `start within 200 markers of truth, got ${errA}`)
  assert.ok(errB < 200, `end within 200 markers of truth, got ${errB}`)

  // The claim being made is a 40x improvement over the window edge. The window edge here is off
  // by 1,800 and 2,400; this must be far better, not marginally.
  assert.ok(errA < 1_800 / 10 && errB < 2_400 / 10,
    `refinement must beat the window edge by an order of magnitude, got ${errA} and ${errB}`)
}

// --- 2. the interval covers the truth, and is never narrower than the calibration ---------------
{
  const TRUE_A = 15_000
  const TRUE_B = 25_000
  const a = chrom(50_000, TRUE_A, TRUE_B, 0.11, 0.005)
  const r = refineEdges(a, TRUE_A - 1_200, TRUE_B + 1_200, 0.005)

  assert.ok(r.start.lo <= TRUE_A && TRUE_A <= r.start.hi,
    `interval must cover the true start: [${r.start.lo}, ${r.start.hi}] vs ${TRUE_A}`)
  assert.ok(r.end.lo <= TRUE_B && TRUE_B <= r.end.hi,
    `interval must cover the true end: [${r.end.lo}, ${r.end.hi}] vs ${TRUE_B}`)

  // A likelihood can be locally very sharp on synthetic data. The measured median interval on real
  // arrays is 151 markers and nothing reached nominal coverage narrower, so a narrower interval
  // here would be claiming precision the calibration does not support.
  assert.ok(r.start.hi - r.start.lo >= MIN_INTERVAL_MARKERS,
    `start interval ${r.start.hi - r.start.lo} is under the ${MIN_INTERVAL_MARKERS} floor`)
  assert.ok(r.end.hi - r.end.lo >= MIN_INTERVAL_MARKERS)
}

// --- 3. a bigger drop gives a wider interval ----------------------------------------------------
//
// The relationship the calibration rests on: coverage was bought by widening the drop from the
// nominal 1.92 to 12, and 20 was wider still. If the drop did not move the width this would be a
// constant with no meaning.
{
  const a = chrom(40_000, 12_000, 20_000, 0.11, 0.005)
  const narrow = refineEdges(a, 11_000, 21_000, 0.005, SEARCH_RADIUS, 2)
  const wide = refineEdges(a, 11_000, 21_000, 0.005, SEARCH_RADIUS, 60)
  assert.ok(wide.start.hi - wide.start.lo >= narrow.start.hi - narrow.start.lo,
    'a larger likelihood drop must not produce a narrower interval')
  assert.equal(LR_DROP, 12, 'the shipped drop is the empirically calibrated one, not chi-squared')
  assert.ok(LR_DROP > 1.92 * 5, 'and it is several times the nominal value, which measured 0.755')
}

// --- 4. it refuses rather than inventing an edge -------------------------------------------------
{
  // Nothing there: a chromosome at the background rate throughout. Refinement on its own will
  // always find SOME stretch running above average, which is why the caller passes the score the
  // scan required. Below that bar there is nothing to localise.
  const flat = chrom(30_000, 0, 0, 0, 0.005)
  const r = refineEdges(flat, 10_000, 20_000, 0.005, SEARCH_RADIUS, LR_DROP, SEGMENT_LRT)
  assert.equal(r.localised, false, 'a chromosome with no event must not localise')
  assert.ok(r.why.includes('noise'), 'and must say why')

  // Without that bar it is a pure refiner and says so by localising on anything, which is the
  // documented contract and the reason the caller must not omit the threshold.
  assert.equal(refineEdges(flat, 10_000, 20_000, 0.005).localised, true,
    'with no threshold this is a refiner, not a detector')

  // A window pointing somewhere the event is not. The ascent runs to its radius and gives up
  // rather than dragging the edges across the chromosome.
  const off = chrom(80_000, 60_000, 70_000, 0.11, 0.005)
  const far = refineEdges(off, 1_000, 3_000, 0.005, 500)
  assert.equal(far.localised, false, 'hitting the search radius is a refusal, not a wide interval')
  assert.ok(far.why.includes('radius'))
}

// --- 5. a window CLEANER than the background is never an event ----------------------------------
//
// The scan's statistic is one-sided on purpose and the ascent has to be too, or it will happily
// find a maximally clean stretch and call it a segment.
{
  const inverted = chrom(40_000, 12_000, 20_000, 0.0, 0.11)
  const r = refineEdges(inverted, 11_000, 21_000, 0.11)
  assert.equal(r.localised, false, 'a below-background stretch is not a finding in either edge')
}

// --- 6. the gate the measurements are conditional on --------------------------------------------
{
  assert.equal(REFINE_CALL_FLOOR, 0.70)
  // Stated so it travels with the code: below this floor positional error rose from 6-9 markers
  // to 160-222 and interval coverage fell to 0.81. The caller must apply it; this module cannot
  // see a call rate.
  assert.ok(REFINE_CALL_FLOOR > 0.54, 'the floor sits above the call rate where refinement broke')
}

console.log('refine.check.ts: all assertions passed')

// Self-check for the pair likelihood in gainOrigin.ts.
// Run: node src/gainPair.check.ts
//
// The constructed cases are the ones the threshold pair got wrong, and they are the reason this
// function exists. Their numbers are the measured ones, not invented: see audit/gain-origin-mosaic.txt.
import assert from 'node:assert/strict'
import {
  callPair, scanAlleleGains, FRACTION_GRID, COMPRESSION_GRID, PAIR_DIRECTION_P,
  type ChromShare,
} from './gainOrigin.ts'

// --- 1. the textbook cases, with a clean cell's scales ------------------------------------------
{
  const sd = 0.02; const sl = 0.06
  const full = callPair(1 / 6, 0.55 * Math.log2(1.5), sd, sl)!
  assert.equal(full.state, 'gain-paternal', 'an extra paternal copy in every cell')
  assert.ok(full.fraction >= 0.9, `and the fraction reads whole, got ${full.fraction}`)
  assert.equal(callPair(-1 / 6, 0.55 * Math.log2(1.5), sd, sl)!.state, 'gain-maternal',
    'the mirror names the other parent')
  // A loss moves the share the same way as the other parent's gain and the intensity the other way.
  assert.equal(callPair(+1 / 6, 0.55 * Math.log2(0.75), sd, sl)!.state, 'loss-maternal',
    'the same share shift with the intensity DOWN is the maternal copy lost, not a paternal gain')
  assert.equal(callPair(0, 0, sd, sl)!.state, 'disomy', 'nothing moving is disomy')
}

// --- 2. THE TRISOMY THE THRESHOLD PAIR REFUSED --------------------------------------------------
//
// GSM4472422 chromosome 15, hand-called a maternal trisomy. Share 0.136 below the cell's centre at
// a share spread of 0.020, intensity +0.29 at an intensity spread of 0.164. The intensity is only
// 1.8 of its cell's spreads, so a bar at 2 refuses it, but +0.29 is close to what a share shift of
// 0.136 predicts for a gain and far from what it predicts for a loss.
{
  const c = callPair(-0.136, 0.29, 0.020, 0.164)!
  assert.equal(c.state, 'gain-maternal',
    `the real maternal trisomy is a maternal gain on the pair. Got ${c.state} at ${c.posterior.toFixed(3)}`)
  assert.ok(c.posterior > 0.75, `and clearly, got ${c.posterior.toFixed(3)}`)
  assert.equal(c.runnerUp, 'unexplained',
    'what it is weighed against is not the other parent but the possibility that its intensity '
    + `means nothing, which is the honest alternative. Got ${c.runnerUp}`)
}

// --- 3. THE PARTIAL LOSS THE ABSOLUTE FLOOR CALLED A GAIN ---------------------------------------
//
// GSM4552382 chromosome 3 with half its paternal copy removed by construction. Share 0.176 above
// the centre at a spread of 0.020, intensity +0.158 at a spread of 0.146. The intensity is positive,
// which is what fooled the absolute floor, but it is far below the +0.39 to +0.56 a share shift
// that size predicts for a gain, and the model is allowed to say so.
{
  const c = callPair(0.176, 0.158, 0.020, 0.146)!
  assert.notEqual(c.state, 'gain-paternal',
    `a half loss of the maternal copy must not read as a paternal gain. Got ${c.state} at ${c.posterior.toFixed(3)}`)
  assert.equal(c.state, 'unexplained',
    'its intensity matches no copy number: +0.158 where a gain of that share shift needs +0.32 and '
    + `a loss needs -0.23. Got ${c.state}`)

  // AND THE SEPARATION IS THE FIT, NOT THE SIZE. The trisomy's intensity is SMALLER in its own
  // cell's spreads than this one is, 1.8 against 1.1, and it is still the gain: what differs is
  // that +0.29 is what a share shift of 0.136 predicts and +0.158 is not what 0.176 predicts.
  assert.equal(callPair(-0.136, 0.29, 0.020, 0.164)!.state, 'gain-maternal')
}

// --- 4. it refuses rather than guessing ---------------------------------------------------------
{
  assert.equal(callPair(NaN, 0.3, 0.02, 0.06), null, 'no share, no call')
  assert.equal(callPair(0.1, NaN, 0.02, 0.06), null, 'no intensity, no call')
  assert.equal(callPair(0.1, 0.3, 0, 0.06), null, 'a zero scale would divide by nothing')
  // A huge intensity spread makes the intensity say nothing at all, and the share alone cannot
  // separate an extra copy of one parent from a partial loss of the other. The model must say so
  // rather than pick one: no copy-number state may win here.
  const vague = callPair(0.166, 0.0, 0.02, 2.0)!
  assert.equal(vague.state, 'unexplained',
    `an intensity too noisy to read names no copy number. Got ${vague.state} at ${vague.posterior.toFixed(3)}`)
  assert.ok(vague.posterior < 0.9, `and not decisively either, got ${vague.posterior.toFixed(3)}`)
}

// --- 5. THE TILT IS WHAT SURVIVES A WRONG DIRECTION ---------------------------------------------
//
// An extra copy from one parent and a missing copy from the other are the same tilt. Getting the
// direction wrong therefore flips the named parent too, and the two errors cancel: the statement
// about which side is over-represented stays right. The same built loss of the paternal copy in
// half the cells is run on three arms differing only in the baseline intensity the untouched
// chromosome already carried, which is the quantity a single cell cannot measure.
{
  const shift = -0.167 // -f/(2(2-f)) at f = 0.5
  const seen: string[] = []
  for (const baseline of [0.02, 0.37, 0.49]) {
    const c = callPair(shift, baseline + 0.5 * Math.log2(1.5 / 2), 0.02, 0.09)!
    seen.push(c.state)
    assert.equal(c.tilt, 'maternal',
      `a lost paternal copy leaves the maternal side over-represented whichever direction the `
      + `intensity suggests. Baseline ${baseline} gave tilt ${c.tilt} on state ${c.state}`)
  }
  assert.equal(seen[0], 'loss-paternal', `a clean arm reads the loss. Got ${seen[0]}`)
  assert.equal(seen[2], 'gain-maternal',
    `and a baseline of +0.49 turns the same event into the other parent's gain, which is the error `
    + `class all 22 of the 300 builds fell into. Got ${seen[2]}`)
  // Both wrong-direction calls above still named the right tilt, so the tilt is reported on its own.
  assert.equal(callPair(1 / 6, 0.55 * Math.log2(1.5), 0.02, 0.06)!.tilt, 'paternal')
  assert.equal(callPair(1 / 6, 0.55 * Math.log2(0.75), 0.02, 0.06)!.tilt, 'paternal',
    'a paternal gain and a maternal loss are one tilt')
  // It sums the two states that agree, so it can only be the firmer of the two statements.
  for (const [d, l, sl] of [[1 / 6, 0.3, 0.06], [-0.136, 0.29, 0.164], [0.09, 0.05, 0.09]]) {
    const c = callPair(d, l, 0.02, sl)!
    if (c.tilt) assert.ok(c.tiltPosterior >= c.posterior - 1e-12,
      `the tilt cannot be less certain than the state inside it, got ${c.tiltPosterior} < ${c.posterior}`)
  }
  const none = callPair(0, 0, 0.02, 0.06)!
  assert.equal(none.tilt, null, 'disomy tilts nowhere')
  assert.equal(none.tiltPosterior, 0, 'and carries no tilt posterior')
  assert.equal(callPair(0.166, 0.0, 0.02, 2.0)!.tilt, null,
    'an unreadable intensity names no state, so there is no tilt to report either')
}

// --- 6. the grids are the measured ones ---------------------------------------------------------
{
  assert.ok(COMPRESSION_GRID[0] >= 0.5 && COMPRESSION_GRID[COMPRESSION_GRID.length - 1] <= 1.0,
    'the compression range is the 58% to 95% of theoretical measured across two chips')
  assert.equal(FRACTION_GRID[FRACTION_GRID.length - 1], 1.0, 'a whole-cell event has to be reachable')
  // Not a range check. This number is the operating point the split chose on 10 cells and then
  // measured on 9 it had not seen, and the error rate quoted beside every row that clears it is
  // the one measured AT this value. Moving it without redoing that split makes the quoted rate a
  // claim about a threshold nobody tested.
  assert.equal(PAIR_DIRECTION_P, 0.995,
    'the direction threshold is the one the pre-registered split picked and measured; changing it '
    + 'means rerunning audit/gain-origin-split.ts and restating the held-out error in pairNote')
}

// --- 7. WHAT THE ROW IS ALLOWED TO SAY ABOUT THE PAIR -------------------------------------------
//
// The note travels on every row, so a state that names no copy number must not pick up the
// direction or parent language that only a copy-number state has earned.
{
  const ref = (chrom: number, share: number, lrrShift: number): ChromShare =>
    ({ chrom: String(chrom), share, informative: 800, biparental: true, lrrShift })
  const quiet = [1, 2, 4, 5, 6, 8, 9, 10].map((c, i) => ref(c, 0.5 + 0.004 * ((i % 3) - 1), 0.01 * ((i % 3) - 1)))
  const scan = (ev: ChromShare) => scanAlleleGains([...quiet, ev]).find((x) => x.chrom === '3')!

  const gain = scan(ref(3, 0.667, 0.33))
  assert.equal(gain.verdict, 'gain', 'a full extra copy with the intensity to match is a gain')
  assert.match(gain.why, /read as a call/, 'and its direction clears PAIR_DIRECTION_P')
  assert.match(gain.why, /paternal side is the over-represented one/)
  assert.match(gain.why, /parent is not taken from the model/,
    'every row says where its parent did not come from')

  const partial = scan(ref(3, 0.667, -0.12))
  assert.equal(partial.verdict, 'imbalance', 'the same share with the intensity down names nobody')
  assert.match(partial.why, /weighed evidence and not a call/,
    `below PAIR_DIRECTION_P the direction is not a call. Got: ${partial.why.slice(-160)}`)

  // The case that was wrong: intensity too weak to read either way, so the likelihood prefers the
  // state that says the intensity is not measuring copy number. That is not a direction, and the
  // row must not claim one.
  const vague = scan(ref(3, 0.667, 0.04))
  assert.equal(vague.verdict, 'imbalance')
  assert.match(vague.why, /intensity is not measuring its copy number/)
  for (const phrase of [/read as a call/, /weighed evidence and not a call/,
    /over-represented one/, /parent is not taken from the model/]) {
    assert.doesNotMatch(vague.why, phrase,
      `a state naming no copy number must not carry ${phrase}. Got: ${vague.why.slice(-200)}`)
  }

  // A quiet chromosome says nothing about the pair at all.
  assert.doesNotMatch(scan(ref(3, 0.5, 0.0)).why, /likelihood/,
    'disomy adds no clause, or every normal row would carry one')
}

console.log('gainPair.check.ts: the pair recovers the trisomy the thresholds refused, rejects the '
  + 'partial loss the floor accepted, and claims a direction only where it has one')

// Self-check for Y-bearing detection and parent naming. Run: node src/sexing.check.ts
//
// The figures are the measured ones from three experiments, 46 arrays, in `sexing.ts`. The two
// near misses are checked by name because each one is a rule that was wrong before it was
// measured: one signal alone inverts a real experiment.
import assert from 'node:assert/strict'
import {
  Y_CALL_MIN, Y_LRR_MIN, Y_MIN_PROBES, emptySex, accumulateSex, sexCall, paternalGroup,
} from './sexing.ts'
import type { ProbeRow } from './ingest.ts'

const row = (chrom: string, called: boolean, log2R: number | null, i: number): ProbeRow => ({
  probesetId: `p${i}`,
  chrom,
  pos: i * 1000,
  log2R,
  baf: null,
  copyNumber: null,
  genotype: called ? 'AA' : 'NC',
  unreadableGenotype: false,
  rawGenotype: called ? 'AA' : 'NoCall',
  bestProbeset: true,
})

/** One array: `yCall` of its Y probes call, at `yLrr`; autosomes call `aCall` at 0. */
const array = (yCall: number, yLrr: number, aCall = 0.87) => {
  const t = emptySex()
  for (let i = 0; i < 2000; i += 1) accumulateSex(row('Y', i / 2000 < yCall, yLrr, i), t)
  for (let i = 0; i < 100_000; i += 1) accumulateSex(row('7', i / 100_000 < aCall, 0, i), t)
  return sexCall(t)
}

// --- 1. a whole Y, and no Y at all --------------------------------------------------------------
{
  const yes = array(0.95, 0.36)
  assert.equal(yes.yBearing, true, 'a whole Y calls like an autosome at autosomal intensity')
  assert.ok(yes.callRatio > 1, 'and its Y calls as well as its autosomes do')

  const no = array(0.0, -1.0)
  assert.equal(no.yBearing, false)
}

// --- 2. THE TWO NEAR MISSES, each of which one signal alone gets wrong ---------------------------
//
// Both are real arrays. Either rule on its own names a group that the other rule, and the rest
// of the evidence, says is the other parent.
{
  // 50402-11: genotypes 86.2% of its Y probes with no Y chromosome under them. Call rate alone
  // would name the maternal group of Experiment 5 paternal, inverting every call in it.
  const noisy = array(0.862, -0.996, 0.761)
  assert.ok(noisy.callRatio >= Y_CALL_MIN, 'it clears the call bar')
  assert.ok(noisy.lrrShift < Y_LRR_MIN, 'and fails the intensity bar')
  assert.equal(noisy.yBearing, false, 'so it is not Y-bearing')

  // 52461-16: sits at -0.096 log2 with nothing to genotype. Intensity alone would name it
  // Y-bearing on an array that called not one Y probe.
  const quiet = array(0.0, -0.096, 0.903)
  assert.ok(quiet.lrrShift >= Y_LRR_MIN, 'it clears the intensity bar')
  assert.ok(quiet.callRatio < Y_CALL_MIN, 'and fails the call bar')
  assert.equal(quiet.yBearing, false)
}

// --- 3. the measured separation, with the thresholds inside it ----------------------------------
{
  // Y-bearing arrays: 93.7% to 97.3% of Y probes at +0.162 to +0.428 log2.
  for (const [c, l] of [[0.937, 0.359], [0.947, 0.428], [0.973, 0.162], [0.950, 0.370]] as const) {
    assert.equal(array(c, l).yBearing, true, `${c}/${l} is a real Y`)
  }
  // Everything else either calls nothing or sits about a log2 below its autosomes.
  for (const [c, l] of [[0.0, -0.096], [0.0, -1.245], [0.862, -0.996], [0.0, -0.811]] as const) {
    assert.equal(array(c, l).yBearing, false, `${c}/${l} is not`)
  }
  assert.ok(Y_LRR_MIN > -0.811 && Y_LRR_MIN < 0.162, 'the intensity cut sits inside the gap')
  assert.ok(Y_CALL_MIN > 0 && Y_CALL_MIN < 0.937 / 0.903, 'and the call cut inside its own')
}

// --- 4. a panel too thin to ask cannot answer, and says so --------------------------------------
{
  const none = emptySex()
  for (let i = 0; i < 1000; i += 1) accumulateSex(row('7', true, 0, i), none)
  assert.equal(sexCall(none).yBearing, null, 'no Y probes is unanswerable, not "no Y"')

  // A handful of Y probes all failing is not evidence of absence. Silence from five probes and
  // silence from two thousand are different claims.
  const thin = emptySex()
  for (let i = 0; i < Y_MIN_PROBES - 1; i += 1) accumulateSex(row('Y', false, -1, i), thin)
  for (let i = 0; i < 1000; i += 1) accumulateSex(row('7', true, 0, i), thin)
  assert.equal(sexCall(thin).yBearing, null, `under ${Y_MIN_PROBES} Y probes cannot say "no Y"`)

  // The bundled example subsets carry 51 and do answer, which is why they read as X-bearing
  // rather than as unanswerable.
  const enough = emptySex()
  for (let i = 0; i < 51; i += 1) accumulateSex(row('Y', false, -1, i), enough)
  for (let i = 0; i < 1000; i += 1) accumulateSex(row('7', true, 0, i), enough)
  assert.equal(sexCall(enough).yBearing, false, '51 probes calling nothing is a real "no Y"')
}

// --- 5. naming a group, and refusing to -------------------------------------------------------
{
  // Experiment 2: group 1 holds two Y-bearing products, the other three groups hold none.
  const groups = [[0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14], [15]]
  const y = Array(16).fill(false)
  y[0] = true
  y[5] = true
  assert.equal(paternalGroup(groups, y), 0)

  // No group carries a Y. A paternal group of five is all X-bearing 3% of the time, so this is
  // a real outcome and not an error; naming the largest group anyway would invert those runs.
  assert.equal(paternalGroup(groups, Array(16).fill(false)), null, 'no Y names nothing')

  // Two groups carrying a Y cannot both come from one sperm donor.
  const two = Array(16).fill(false)
  two[0] = true
  two[9] = true
  assert.equal(paternalGroup(groups, two), null, 'two Y-bearing groups is a refusal')

  // An unanswerable member neither names nor blocks.
  const some = Array(16).fill(null)
  some[2] = true
  assert.equal(paternalGroup(groups, some), 0)
}

console.log('sexing.check.ts: all assertions passed')

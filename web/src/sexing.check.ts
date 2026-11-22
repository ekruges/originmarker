// Self-check for Y-bearing detection and parent naming. Run: node src/sexing.check.ts
//
// The figures are the measured ones from three experiments, 46 arrays, in `sexing.ts`. The two
// near misses are checked by name because each one is a rule that was wrong before it was
// measured: one signal alone inverts a real experiment.
import assert from 'node:assert/strict'
import {
  Y_CALL_MIN, Y_LRR_MIN, Y_MIN_PROBES, Y_BAF_READABLE, Y_BAF_PRESENT, Y_BAF_ABSENT, emptySex,
  accumulateSex, sexCall, paternalGroup,
} from './sexing.ts'
import type { ProbeRow } from './ingest.ts'
import { callSexChromosomes } from './sexChromosome.ts'
import { sexMaterial, shownSex } from './sexMaterial.ts'
import type { StageCall } from './stage.ts'

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
  // A real array: genotypes 86.2% of its Y probes with no Y chromosome under them. Call rate
  // alone would name a maternal group paternal, inverting every call in it.
  const noisy = array(0.862, -0.996, 0.761)
  assert.ok(noisy.callRatio >= Y_CALL_MIN, 'it clears the call bar')
  assert.ok(noisy.lrrShift < Y_LRR_MIN, 'and fails the intensity bar')
  assert.equal(noisy.yBearing, null, 'so, with no allele fractions to settle it, it is not read')

  // A real array: sits at -0.096 log2 with nothing to genotype. Intensity alone would name it
  // Y-bearing on an array that called not one Y probe.
  const quiet = array(0.0, -0.096, 0.903)
  assert.ok(quiet.lrrShift >= Y_LRR_MIN, 'it clears the intensity bar')
  assert.ok(quiet.callRatio < Y_CALL_MIN, 'and fails the call bar')
  assert.equal(quiet.yBearing, null, 'no call is not no Y on a chip that blanks chrY by its X')
}

// --- 3. the measured separation, with the thresholds inside it ----------------------------------
{
  // Y-bearing arrays: 93.7% to 97.3% of Y probes at +0.162 to +0.428 log2.
  for (const [c, l] of [[0.937, 0.359], [0.947, 0.428], [0.973, 0.162], [0.950, 0.370]] as const) {
    assert.equal(array(c, l).yBearing, true, `${c}/${l} is a real Y`)
  }
  // Calling nothing a log2 below its autosomes is no Y; one signal alone is not read.
  for (const [c, l] of [[0.0, -1.245], [0.0, -0.811]] as const) {
    assert.equal(array(c, l).yBearing, false, `${c}/${l} is not`)
  }
  for (const [c, l] of [[0.0, -0.096], [0.862, -0.996]] as const) {
    assert.equal(array(c, l).yBearing, null, `${c}/${l} is not read`)
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

// --- 6. THE Y FROM ITS OWN MARKERS: allele fractions first, intensity as the check ----------------
//
// On this chip a Y beside two X is left uncalled, so the calls say nothing a Y's allele fractions do
// not. `share` of the chrY markers sit at a pole; `floor` of the array's own homozygous calls do.
const alleles = (o: { yCall: number; yLrr: number; share: number; floor: number }) => {
  const t = emptySex()
  for (let i = 0; i < 800; i += 1) {
    accumulateSex({ ...row('Y', i / 800 < o.yCall, o.yLrr, i), baf: i / 800 < o.share ? 0.02 : 0.5 }, t)
  }
  for (let i = 0; i < 100_000; i += 1) {
    const pole = ((i * 7919) % 1000) / 1000 < o.floor
    accumulateSex({ ...row('7', true, 0, i), baf: pole ? 0.98 : 0.7 }, t)
  }
  return sexCall(t)
}
{
  assert.deepEqual([Y_BAF_READABLE, Y_BAF_PRESENT, Y_BAF_ABSENT], [0.90, 0.75, 0.60], 'the bars as measured')
  const clean = { floor: 0.97 }
  // R2: a Y the chip left uncalled, read from its markers, and the reason says so with the numbers.
  const uncalled = alleles({ ...clean, yCall: 0, yLrr: -0.08, share: 0.88 })
  assert.equal(uncalled.yBearing, true)
  assert.ok(Math.abs(uncalled.bafFloor! - 0.97) < 0.01, `floor ${uncalled.bafFloor}`)
  assert.equal(uncalled.why, 'chrY is present: its markers read one allele each (88% of them, where '
    + "this array's own homozygous calls read 97%) at male intensity, though the chip left them uncalled")
  const called = alleles({ ...clean, yCall: 0.95, yLrr: 0.1, share: 0.88 })
  assert.equal(called.yBearing, true)
  assert.equal(called.why, undefined, 'a called Y keeps the plain reason')
  // R3: a woman's markers on an export whose intensity reads her chrY near her autosomes.
  assert.equal(alleles({ ...clean, yCall: 0, yLrr: -0.1, share: 0.5 }).yBearing, false)
  // R3: calls and intensity both say present against the markers: not read.
  const against = alleles({ ...clean, yCall: 0.95, yLrr: 0.2, share: 0.45 })
  assert.equal(against.yBearing, null)
  assert.ok(against.why!.startsWith('chrY was not read: only 45% of its markers'), against.why)
  // R4: between the bars, calls and intensity decide only where they agree.
  assert.equal(alleles({ ...clean, yCall: 0, yLrr: -0.98, share: 0.68 }).yBearing, false)
  assert.equal(alleles({ ...clean, yCall: 0.95, yLrr: 0.1, share: 0.68 }).yBearing, true)
  const between = alleles({ ...clean, yCall: 0, yLrr: -0.08, share: 0.62 })
  assert.equal(between.yBearing, null)
  assert.ok(between.why!.includes('between the 60% that is no Y and the 75% that is one'), between.why)
  assert.ok(between.why!.includes('says present'), between.why)
  // R4: a high share at an intensity under the bar is read only where calls and intensity agree.
  assert.equal(alleles({ ...clean, yCall: 0, yLrr: -0.7, share: 0.8 }).yBearing, false)
  assert.equal(alleles({ ...clean, yCall: 0.95, yLrr: -0.7, share: 0.8 }).yBearing, null)
  // R5: allele fractions not readable on a noisy array: calls and intensity, only where they agree.
  const noisy = { floor: 0.8, share: 0.9 }
  assert.equal(alleles({ ...noisy, yCall: 0.95, yLrr: 0.1 }).yBearing, true)
  assert.equal(alleles({ ...noisy, yCall: 0, yLrr: -1.0 }).yBearing, false)
  const loud = alleles({ ...noisy, yCall: 0, yLrr: -0.1 })
  assert.equal(loud.yBearing, null, 'no call is no longer no Y')
  assert.ok(loud.why!.includes('read only 80% at one allele, under the 90%'), loud.why)

  // R6: a Y not read prints no constitution of X alone, with a parent or without.
  const stage = { stage: 'single-cell', ploidy: 'diploid' } as StageCall
  const shown = sexMaterial({ parentLoaded: false, stage, hetRate: 0.16, dropIn: 0.02 })
  const notShown = sexMaterial({ parentLoaded: false, stage, hetRate: 0.15, dropIn: 0.14 })
  const unread = { ...emptySex(), yTotal: 800, yCalled: 0, autoCalled: 87_000, autoTotal: 100_000,
    yLrr: Array(800).fill(-0.08), autoLrr: Array(1000).fill(0), seen: 100_000,
    yBaf: 800, yBafPole: 496, homBaf: 8700, homBafPole: 8439 }
  const np = { material: shown, diploid: true, twoInsideNoise: true, chipGenotypesY: true }
  for (const c of [
    callSexChromosomes({ ...np, hetSex: 'female', hetRatio: 1.0, xShift: 0.02, xZ: 0.3, y: unread }),
    callSexChromosomes({ ...np, material: notShown, diploid: false, xShift: 0.02, xZ: 0.3, y: unread }),
    callSexChromosomes({ diploid: true, chipGenotypesY: true, hetSex: 'female', xShift: 0.02, xZ: 0.3, y: unread }),
  ]) {
    assert.equal(c.constitution, 'unresolved', c.why)
    assert.ok(c.why.includes('chrY was not read') && c.why.includes('no constitution of X alone'), c.why)
  }

  // R7: a Y with a heterozygous X is its own finding, with a parent or without; bulk DNA with chrX
  // at the autosomes' level keeps 47,XXY.
  const yTwoX = { ...unread, yBafPole: 704 }
  const r7 = 'chrY is present: its markers read one allele each (88% of them, where this array\'s own '
    + 'homozygous calls read 97%) at male intensity, though the chip left them uncalled; chrX reads two '
    + 'different X by its heterozygous calls. A Y with two different X is 47,XXY or a mix of cells, '
    + 'which one array does not separate'
  const noParent = callSexChromosomes({ ...np, hetSex: 'female', hetRatio: 1.0, xShift: 0.02, xZ: 0.3, y: yTwoX })
  assert.equal(noParent.constitution, '47,XXY or a mix of cells')
  assert.equal(noParent.aneuploid, false, 'a finding, not a constitution')
  assert.equal(noParent.why, r7)
  const withParent = callSexChromosomes({ diploid: true, chipGenotypesY: true, hetSex: 'female', hetRatio: 1.0,
    xShift: -0.1, xZ: -1, y: yTwoX })
  assert.equal(withParent.constitution, '47,XXY or a mix of cells')
  assert.equal(withParent.why, r7)
  assert.equal(callSexChromosomes({ diploid: true, bulk: true, chipGenotypesY: true, hetSex: 'female',
    hetRatio: 1.0, xShift: 0.05, xZ: 0.5, y: yTwoX }).constitution, '47,XXY')
  // The sex word: never female beside a chrY, and nothing else moves.
  assert.equal(shownSex('female', shown, true), 'not read: chrY present with a heterozygous X')
  assert.equal(shownSex('female', undefined, true), 'not read: chrY present with a heterozygous X')
  assert.equal(shownSex('female', undefined, false), 'female')
  assert.equal(shownSex('female', undefined, null), 'female')
  assert.equal(shownSex('male', undefined, true), 'male')
}

console.log('sexing.check.ts: all assertions passed')

// Self-check for the one-parent origin call. Run: node src/oneParentOrigin.check.ts
//
// This module exists because two previous attempts at a one-parent call failed in ways that were
// invisible in their output: one could only ever move in a single direction, and one returned the
// same answer to every input. So the assertions below are aimed at those two failures first, and
// at the arithmetic second.
//
//   SYMMETRY   both loss directions must be reachable, and neither may be privileged.
//   NOT CONSTANT   different inputs must produce different answers.
//   MENDELIAN   the allele the loaded parent does not carry must dominate the evidence, since it
//               is the one observation that no amount of dropout can manufacture.
import assert from 'node:assert/strict'
import {
  SYSTEMATIC_ERROR_BOUND, VALIDATION_UNITS,
  callOneParentOrigin, informative, inferDropout, CALL_POSTERIOR, MIN_MARKERS, DEFAULT_Q,
  MAX_REGION_HET, DROP_IN, OBVIOUS_EVENT_ACCURACY,
} from './oneParentOrigin.ts'
import { parentNamed } from './defects.ts'
import type { AB } from './informativity.ts'

type Pair = readonly [AB, AB]
/** A region where the loaded parent is AA throughout and the sample reads as given. */
const region = (sample: AB[], parent: AB = 'AA'): Pair[] => sample.map((s) => [parent, s] as Pair)
/** Deterministic draw, so a failure is reproducible. */
const draw = (n: number, f: (u: number, i: number) => AB): AB[] => {
  let seed = 20260813
  return Array.from({ length: n }, (_, i) => {
    seed = (seed * 1664525 + 1013904223) >>> 0
    return f(seed / 4294967296, i)
  })
}

// --- 1. the known parent's copy is gone: its allele is absent, the other one appears -------------
//
// Parent AA, so the child carries only the other parent's transmission: A at rate 1-q, B at rate q.
// The B observations are Mendelian proof the parent's copy is not there.
{
  const s = draw(400, (u) => (u < DEFAULT_Q ? 'BB' : 'AA'))
  const c = callOneParentOrigin(region(s), 0.308)
  assert.equal(c.verdict, 'known-parent-lost', `got ${c.verdict}: ${c.why}`)
  // The CALLING GATE still runs on the raw likelihood, which is what CALL_POSTERIOR is for. The
  // REPORTED number is bounded by what the validation supports, so it cannot exceed 1 - 0.206 and
  // must not be compared against a threshold the raw quantity is measured on.
  assert.ok(c.posterior <= 1 - SYSTEMATIC_ERROR_BOUND + 1e-9,
    `reported ${c.posterior} exceeds what ${VALIDATION_UNITS} validation units support`)
  assert.ok(c.posterior > 0.5, 'and a made call must still sit clearly above a coin flip')
}

// --- 8. dropout cannot manufacture the Mendelian evidence -----------------------------------------
//
// Dropout removes an allele; it never invents one. So raising the dropout rate on an intact region
// must not turn it into a loss of the loaded parent's copy.
{
  const s = draw(400, (u) => (u < DEFAULT_Q * 0.7 ? 'AB' : 'AA'))
  for (const ado of [0.05, 0.3, 0.6]) {
    const c = callOneParentOrigin(region(s), ado)
    assert.notEqual(c.verdict, 'known-parent-lost',
      `dropout at ${ado} must not produce the loaded parent's absence`)
  }
}

console.log('oneParentOrigin.check.ts: all assertions passed, including symmetry, '
  + 'not-a-constant, and the Mendelian evidence surviving dropout')

// --- 9. dropout is inferred from the sample, not declared ----------------------------------------
//
// The drag-and-drop contract: a trophectoderm biopsy and a blastomere differ six-fold in dropout
// and the user is never asked which they dropped. Heterozygosity is the readout, since dropout
// removes one allele of a heterozygote and so depresses the observed rate below the bulk figure.
{
  const cases: [number, number, number, string][] = [
    // sample het, expected dropout floor, ceiling, label
    [0.168, 0.00, 0.05, 'bulk-quality reads near zero dropout'],
    [0.160, 0.01, 0.10, 'trophectoderm, mild'],
    [0.135, 0.15, 0.25, 'single ESC'],
    [0.116, 0.26, 0.38, 'blastomere'],
  ]
  for (const [het, lo, hi, label] of cases) {
    const d = inferDropout(het)
    assert.ok(d >= lo && d <= hi, `${label}: het ${het} gave dropout ${d.toFixed(3)}, want ${lo}-${hi}`)
  }
  // Monotone: worse heterozygosity must never imply less dropout.
  assert.ok(inferDropout(0.116) > inferDropout(0.160), 'more dropout must read as more dropout')
  // Degenerate input falls back to the worst stage rather than to zero, which would be optimistic.
  assert.equal(inferDropout(NaN), 0.308)
  assert.equal(inferDropout(0), 0.308)
  // Bounded, so a pathological array cannot drive the likelihood to a degenerate value.
  assert.ok(inferDropout(0.0001) <= 0.6 && inferDropout(0.5) >= 0.01)
}

// --- an impossible heterozygosity is refused, not scored ------------------------------------------
//
// Found on a real array. A blastomere with a call-rate collapse on chromosome 1, 69% no-call,
// read 59.7% heterozygous across the markers where the father is homozygous. No genome does that:
// where the loaded parent is homozygous a biparental sample is heterozygous only when the other
// parent transmitted the allele this one lacks. The likelihood took the count at face value,
// because a heterozygote is near-impossible under either deletion hypothesis, and returned
// both-copies-present at posterior 1.000 from genotypes carrying no information.
{
  const junk: [AB, AB][] = []
  for (let i = 0; i < 2_000; i += 1) junk.push(['AA', i % 10 < 6 ? 'AB' : 'AA'])
  const c = callOneParentOrigin(junk, 0.30)
  assert.equal(c.verdict, 'refused', 'a 60% heterozygous region must not receive a verdict')
  assert.ok(Number.isNaN(c.posterior))
  assert.ok(c.why.includes('not measuring'))

  // A real biparental region sits near the allele frequency and must still be scored.
  //
  // BUILT FROM DEFAULT_Q RATHER THAN WRITTEN DOWN. This fixture used to sit at 30 percent
  // heterozygous, which was "expected" only because q was 0.30. At the measured q a real biparental
  // region sits near 12 percent, and a fixture pinned to the old number asserts that a region twice
  // the plausible heterozygosity must be scored, which is the opposite of what this guard is for.
  const pct = Math.max(1, Math.round(DEFAULT_Q * 10))
  const real: [AB, AB][] = []
  for (let i = 0; i < 2_000; i += 1) real.push(['AA', i % 10 < pct ? 'AB' : 'AA'])
  assert.notEqual(callOneParentOrigin(real, 0.30).verdict, 'refused',
    `a region at the expected heterozygosity (${pct * 10}%, from q=${DEFAULT_Q}) must still be called`)

  // And the ceiling must sit ABOVE that expectation with room, or it rejects real regions.
  assert.ok(MAX_REGION_HET > DEFAULT_Q,
    'the refusal ceiling must exceed the heterozygosity a real region shows, or every real region '
    + `is refused: ceiling ${MAX_REGION_HET} against q ${DEFAULT_Q}`)
  // DERIVED, NOT WRITTEN DOWN. This read `0.30 + DROP_IN`, which pinned the ceiling to a q that
  // has since been measured at less than half that. A check that hardcodes the constant it is
  // guarding cannot notice the constant being wrong.
  assert.ok(MAX_REGION_HET > DEFAULT_Q + DROP_IN,
    `the ceiling must clear q plus drop-in, or real regions are refused: `
    + `${MAX_REGION_HET} against ${DEFAULT_Q} + ${DROP_IN}`)
}

console.log('oneParentOrigin.check.ts: dropout inference pinned across four stages')

// --- THE OBVIOUS EVENT IS WHERE THIS CHANNEL EARNS ITS PLACE -----------------------------------
//
// A blastomere carries both parental genomes, so on a clear whole-chromosome loss this channel
// needs no detection floor: losing a parent's copy leaves an allele that parent does not have, and
// dropout removes alleles without inventing one. Measured by construction on real biparental
// arrays, removing one parent's copy across a chromosome in both directions from the same array.
{
  const m = OBVIOUS_EVENT_ACCURACY

  // THE WRONG-PARENT ARM MUST BE WORSE. This is the assertion the old figure could never have
  // passed: its harness scored 44 of 44 with the true father, the true mother, and two arrays that
  // are not parents at all. A parental measurement whose control does not move is not a parental
  // measurement, and this check exists so that can never ship again.
  assert.ok(m.specificityWrongParent.falseCallRate > m.specificity.falseCallRate * 3,
    'loading an unrelated adult must produce materially more false calls than loading the real '
    + `parent: ${m.specificityWrongParent.falseCallRate} against ${m.specificity.falseCallRate}. `
    + 'If these are close, the number is measuring the construction rather than the parentage.')
  // And the intervals must not overlap, or the separation is not established.
  assert.ok(m.specificityWrongParent.lo > m.specificity.hi,
    'the two arms\' intervals must be disjoint for the separation to be claimed')

  // THE DIRECTIONS ARE NOT EQUALLY RELIABLE, and the reporting layer depends on knowing which.
  // Absence of the loaded parent's allele is Mendelian and dropout cannot manufacture it. Absence
  // of the OTHER parent's copy is read from heterozygosity that is not there, which is exactly what
  // dropout produces, so it is far weaker on amplified material.
  assert.ok(m.knownParentLost.accuracy > 0.8,
    `known-parent-lost is the reportable direction: ${m.knownParentLost.accuracy}`)
  assert.ok(m.otherParentLost.accuracy < 0.5,
    'other-parent-lost is NOT reportable on this material, and a figure above 0.5 here would mean '
    + `the asymmetry had changed and the reporting rules must be revisited: ${m.otherParentLost.accuracy}`)
  assert.ok(m.knownParentLost.accuracy > m.otherParentLost.accuracy * 2,
    'the asymmetry between the two directions is the finding, not an artefact')

  // Every figure carries the array count it was clustered over, since many calls from one array
  // are one array. The per-material split is one level deeper and its rows are smaller by
  // construction: 76 arrays across five materials cannot each rest on 12.
  for (const [k, v] of Object.entries(m)) {
    if (k === 'specificityByMaterial') continue
    assert.ok((v as { arrays: number }).arrays >= 12,
      `${k} must state how many independent arrays it rests on`)
  }
  let split = 0
  for (const [k, v] of Object.entries(m.specificityByMaterial)) {
    assert.ok(v.arrays >= 1, `${k} must state its array count`)
    split += v.arrays
    // Each row is a rate with an interval, and its wrong-parent arm alongside it, or it cannot be
    // read without going back to the pooled figure it exists to replace.
    assert.ok(v.falseCallRate >= 0 && v.falseCallRate <= 1 && v.hi >= v.lo)
    assert.ok(v.wrong >= v.falseCallRate,
      `${k}: the unrelated-adult arm must not be BETTER than the true-parent arm`)
  }
  assert.ok(split >= 70 && split <= m.specificity.arrays,
    `the split must account for the arrays the pooled figure rests on: ${split} of `
    + `${m.specificity.arrays}`)
  // AND THE MATERIALS ARE NOT INTERCHANGEABLE, which is the reason the split exists. One cell
  // against five to ten is an order of magnitude in the error rate, and quoting the pooled number
  // on a blastomere understates it by four times.
  const blast = m.specificityByMaterial.blastomere
  const te = m.specificityByMaterial.trophectoderm
  assert.ok(blast.falseCallRate > te.falseCallRate,
    'a single cell cannot have the same false-call rate as a biopsy of several; if these ever '
    + 'match, the measurement has stopped separating them')
  assert.ok(blast.falseCallRate > m.specificity.falseCallRate * 2,
    'the pooled figure materially understates a blastomere, which is why callers must quote the '
    + `row: pooled ${m.specificity.falseCallRate} against blastomere ${blast.falseCallRate}`)
}

// --- BAND F NAMES NOBODY IN THIS CHANNEL EITHER -------------------------------------------------
//
// The dosage channel has withheld the parent in band F since a gynogenetic sample, a genome with no
// paternal contribution at all, was told its paternal copy was the one lost. This channel shares the
// band ladder and did not share the guard: it refused only when its RAW posterior missed the calling
// bar, so a run decisive under the model and knocked to chance by the systematic bound still printed
// a parent. audit/requirements.check.ts printed two of them unprompted, at 0.4824 and 0.5291.
{
  let named = 0
  let bandF = 0
  // A sweep, not a fixture: every combination the audit walks, plus the neighbourhood around the
  // two rows that were wrong. If any of them names a parent in band F this fails.
  for (const n of [40, 60, 80, 100, 150, 200, 300, 400, 800]) {
    for (const exclusive of [0, 1, 2, 3, 5, 8, 15, 40, 120]) {
      if (exclusive > n) continue
      const pairs: Pair[] = []
      for (let i = 0; i < n; i += 1) pairs.push(['AA', i < exclusive ? 'BB' : 'AA'])
      const r = callOneParentOrigin(pairs, 0.20)
      if (r.band !== 'F') continue
      bandF += 1
      if (parentNamed(r.verdict, 'paternal') !== null) {
        named += 1
        console.error(`band F named a parent: n=${n} exclusive=${exclusive} `
          + `verdict=${r.verdict} posterior=${r.posterior.toFixed(4)}`)
      }
    }
  }
  assert.ok(bandF > 0, 'the sweep must actually reach band F, or it is asserting nothing')
  assert.equal(named, 0,
    `${named} of ${bandF} band-F rows named a parent. Band F is where the tool says it cannot `
    + 'grade its own answer, and a parent attached to one is a coin flip wearing a name.')
  console.log(`  band F: ${bandF} rows in the sweep, 0 naming a parent`)
}

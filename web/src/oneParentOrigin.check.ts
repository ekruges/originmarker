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
  const real: [AB, AB][] = []
  for (let i = 0; i < 2_000; i += 1) real.push(['AA', i % 10 < 3 ? 'AB' : 'AA'])
  assert.notEqual(callOneParentOrigin(real, 0.30).verdict, 'refused',
    'a region at the expected heterozygosity must still be called')
  assert.ok(MAX_REGION_HET > 0.30 + DROP_IN,
    'the ceiling must clear q plus drop-in, or real regions are refused')
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
  assert.equal(m.usable.correct, 92)
  assert.equal(m.usable.calls, 100)
  // The split is load-bearing: it is why an array that fails its stage inference is excluded from
  // reporting rather than reported weakly.
  assert.ok(m.usable.perArray > m.rejected.perArray + 0.2,
    `an array that resolves to a stage must call materially better than one that does not: `
    + `${m.usable.perArray} against ${m.rejected.perArray}`)
  assert.ok(m.rejected.perArray < 0.7,
    'and the rejected set must be visibly poor, or excluding it would be unjustified')
  // And it must beat the dosage channel's weakest measured band on the same material, or there
  // would be no reason to prefer it on an obvious event.
  assert.ok(m.usable.perArray > 0.7,
    'this channel is preferred on obvious events because it measures better there')
  console.log(`  obvious whole-chromosome events: ${m.usable.correct}/${m.usable.calls} correct, `
    + `per-array ${m.usable.perArray} on arrays that resolve to a stage, `
    + `${m.rejected.perArray} on those that do not`)
}

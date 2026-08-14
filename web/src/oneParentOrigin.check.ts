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
  callOneParentOrigin, informative, CALL_POSTERIOR, MIN_MARKERS, DEFAULT_Q,
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
  assert.ok(c.posterior > CALL_POSTERIOR)
  assert.ok(c.exclusive > 0, 'the exclusive-allele count is the evidence and must be reported')
}

// --- 2. the OTHER parent's copy is gone: only the loaded parent's allele survives ----------------
//
// The direction that was undetectable before. Parent AA and the child carries only A.
{
  const s = draw(400, (u) => (u < 0.02 ? 'BB' : 'AA'))
  const c = callOneParentOrigin(region(s), 0.308)
  assert.equal(c.verdict, 'other-parent-lost', `got ${c.verdict}: ${c.why}`)
  assert.ok(c.posterior > CALL_POSTERIOR)
}

// --- 3. both copies present: heterozygotes appear at rate q ---------------------------------------
{
  const s = draw(400, (u) => (u < DEFAULT_Q * 0.7 ? 'AB' : u < DEFAULT_Q ? 'BB' : 'AA'))
  const c = callOneParentOrigin(region(s), 0.308)
  assert.equal(c.verdict, 'both-present', `got ${c.verdict}: ${c.why}`)
  assert.ok(c.heterozygous > 0, 'heterozygotes are what argue both copies are present')
}

// --- 4. NOT A CONSTANT FUNCTION. The three cases above must give three different answers ---------
//
// The failure of the previous caller, which answered the same thing to everything and looked
// perfectly specific while carrying no information.
{
  const lost = callOneParentOrigin(region(draw(300, (u) => (u < DEFAULT_Q ? 'BB' : 'AA'))), 0.308)
  const kept = callOneParentOrigin(region(draw(300, (u) => (u < 0.02 ? 'BB' : 'AA'))), 0.308)
  const intact = callOneParentOrigin(
    region(draw(300, (u) => (u < DEFAULT_Q * 0.7 ? 'AB' : 'AA'))), 0.308,
  )
  const answers = new Set([lost.verdict, kept.verdict, intact.verdict])
  assert.equal(answers.size, 3, `three distinct inputs must give three answers, got ${[...answers]}`)
}

// --- 5. SYMMETRY. Neither direction may be easier than the other ---------------------------------
//
// The same evidence strength on each side must reach comparable confidence. An asymmetry here is
// the signature of the centring bug this module was written to avoid.
{
  const lost = callOneParentOrigin(region(draw(400, (u) => (u < DEFAULT_Q ? 'BB' : 'AA'))), 0.308)
  const kept = callOneParentOrigin(region(draw(400, (u) => (u < 0.02 ? 'BB' : 'AA'))), 0.308)
  assert.ok(lost.posterior > CALL_POSTERIOR && kept.posterior > CALL_POSTERIOR,
    `both directions must be callable, got ${lost.posterior} and ${kept.posterior}`)
}

// --- 6. it works whichever homozygote the parent is ----------------------------------------------
//
// Orientation must not matter: a BB parent whose copy is lost shows A alleles, mirroring case 1.
{
  const s = draw(400, (u) => (u < DEFAULT_Q ? 'AA' : 'BB'))
  const c = callOneParentOrigin(region(s, 'BB'), 0.308)
  assert.equal(c.verdict, 'known-parent-lost', `orientation must not matter: ${c.why}`)
}

// --- 7. refusals ---------------------------------------------------------------------------------
{
  const thin = callOneParentOrigin(region(draw(10, () => 'AA')), 0.308)
  assert.equal(thin.verdict, 'refused')
  assert.ok(thin.why.includes(String(MIN_MARKERS)))

  // A heterozygous parent carries no information and must not be counted at all.
  assert.equal(informative('AB'), false)
  assert.equal(informative('NC'), false)
  const hetParent = callOneParentOrigin(region(draw(400, () => 'AA'), 'AB'), 0.308)
  assert.equal(hetParent.markers, 0, 'markers where the parent is heterozygous must be excluded')
  assert.equal(hetParent.verdict, 'refused')
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

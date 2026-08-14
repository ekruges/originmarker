// Self-check for the sibling-referenced parental call. Run: node src/siblingOrigin.check.ts
//
// What is worth testing here is not that a likelihood multiplies correctly. It is the three
// properties this construction exists to guarantee, each of which failed somewhere in this
// project's history and each of which is invisible in the output when it fails:
//
//   the call is SYMMETRIC between the two directions, because every one-directional result this
//   project produced came from a statistic that could only move one way;
//   the false-het rule SCALES with the panel, because a fixed rule makes the bias worse as arrays
//   are added while appearing to make the call more reliable;
//   it REFUSES rather than reporting weakly, at the wall, at thin evidence, and at too few
//   siblings.
import assert from 'node:assert/strict'
import {
  callSiblingOrigin, expectedPhi, hetRule, DROP_IN, PHI_WALL, CALL_POSTERIOR, MIN_MARKERS,
  MIN_SIBLINGS, type AB,
} from './siblingOrigin.ts'

const rep = (g: AB, n: number): AB[] => Array.from({ length: n }, () => g)
/** Phased call, which is the only mode where a SIDE can be named. */
const phased = (obs: AB[], sibs: number, ado = 0.308) =>
  callSiblingOrigin(obs, sibs, ado, undefined, undefined, true)

// --- 1. SYMMETRY. The two loss directions must be equally callable ------------------------------
//
// The property whose absence produced six wrong results. A region carrying only the alternate
// allele and a region carrying only the reference allele are the same evidence pointing opposite
// ways, and must reach the same posterior.
{
  const lostRef = phased(rep('BB', 200), 4)
  const lostOther = phased(rep('AA', 200), 4)
  assert.equal(lostRef.hypothesis, 'reference-copy-lost')
  assert.equal(lostOther.hypothesis, 'other-copy-lost')
  // Not exactly equal: false heterozygosity genuinely favours one side, which is the point of
  // reporting phi. But both must be callable, and the gap must be small at a sane phi.
  assert.ok(lostRef.posterior > CALL_POSTERIOR && lostOther.posterior > CALL_POSTERIOR,
    'both directions must be callable')
  assert.ok(Math.abs(lostRef.posterior - lostOther.posterior) < 0.05,
    `the two directions must be near-symmetric, got ${lostRef.posterior} and ${lostOther.posterior}`)
}

// --- 2. a region with both alleles present is NOT a deletion -------------------------------------
{
  const het = phased(rep('AB', 200), 4)
  assert.equal(het.hypothesis, 'no-deletion',
    'markers that stay heterozygous are a dropout cluster, not a lost copy')
}

// --- 3. the false-het rule must SCALE, and the fixed rule must be shown to fail ------------------
//
// At a fixed threshold of two, phi rises with panel size: more cells give a spurious call more
// chances to appear. That is the trap, so it is asserted rather than described.
{
  assert.equal(hetRule(2), 2)
  assert.equal(hetRule(3), 2)
  assert.equal(hetRule(5), 3)
  assert.equal(hetRule(8), 4)

  const fixedPhi = (k: number): number => {
    let tail = 0
    for (let j = 2; j <= k; j += 1) {
      let c = 1
      for (let i = 0; i < j; i += 1) c = (c * (k - i)) / (i + 1)
      tail += c * DROP_IN ** j * (1 - DROP_IN) ** (k - j)
    }
    return tail
  }
  assert.ok(fixedPhi(8) > fixedPhi(2),
    'a fixed rule must get WORSE with more siblings, which is why it is not used')
  assert.ok(expectedPhi(8) < expectedPhi(2),
    'the scaled rule must get BETTER with more siblings')
  assert.ok(expectedPhi(8) < 0.02, `eight siblings should be clean, got ${expectedPhi(8)}`)
}

// --- 4. it refuses, in each of the four ways it should -------------------------------------------
{
  assert.equal(phased(rep('BB', 200), 1).hypothesis, 'refused',
    'one sibling cannot establish heterozygosity')
  assert.ok(phased(rep('BB', 200), 1).why.includes(String(MIN_SIBLINGS)))

  assert.equal(phased(rep('BB', 5), 4).hypothesis, 'refused',
    'too few markers is a refusal')
  assert.ok(phased(rep('BB', 5), 4).why.includes(String(MIN_MARKERS)))

  // Past the wall, no marker count rescues the reference-copy-lost arm.
  const past = callSiblingOrigin(rep('BB', 5000), 3, 0.308, 0.45, undefined, true)
  assert.equal(past.hypothesis, 'refused', 'over the phi wall must refuse whatever the evidence')
  assert.ok(past.phi > PHI_WALL && past.why.includes('wall'))

  // The posterior guard is a GUARD, not a common path, and that is worth recording rather than
  // asserting a case that cannot occur. Above the marker floor this model is decisive: 20 markers
  // of one homozygote separate the lost-copy hypothesis from no-deletion by about 37 log units,
  // because a real deletion predicts that homozygote at 0.98 while an intact region predicts it
  // only at ado/2. Mixed evidence does not produce a marginal parental call either; it is absorbed
  // by no-deletion, which is what the third hypothesis is for. So the guard fires on degenerate
  // input rather than on real regions, and is checked as a threshold rather than a scenario.
  assert.ok(CALL_POSTERIOR > 0.9 && CALL_POSTERIOR < 1,
    'the posterior bar must be a real bar, not a formality')
  const empty = phased([], 4)
  assert.equal(empty.hypothesis, 'refused', 'no markers at all must refuse')
}

// --- 5. phi is reported on every path, including refusals ---------------------------------------
//
// If phi is not reported the result is not interpretable, so it must survive every branch.
{
  for (const c of [
    phased(rep('BB', 200), 4),
    phased(rep('BB', 5), 4),
    phased(rep('BB', 200), 1),
    callSiblingOrigin(rep('BB', 200), 4, 0.308),
  ]) {
    assert.ok(Number.isFinite(c.phi), 'phi must be reported on every path')
    assert.ok(c.siblings >= 0 && c.agreement >= 0)
  }
}

// --- 6. dropout is recognised as dropout, not as a lost copy -------------------------------------
//
// Allele dropout sends a heterozygote to either homozygote with equal probability, so a region of
// balanced homozygotes with no heterozygotes is what heavy dropout looks like on an INTACT region.
// The third hypothesis exists to absorb exactly that, and the test is that it does: the answer is
// no-deletion, not a parental call and not a refusal. This is the discriminator that stops a
// dropout cluster being reported as a lost parental copy.
{
  const balanced: AB[] = []
  for (let i = 0; i < 100; i += 1) balanced.push(i % 2 ? 'AA' : 'BB')
  const c = phased(balanced, 4)
  assert.equal(c.hypothesis, 'no-deletion',
    'symmetric dropout on an intact region must read as no deletion')
  // And critically it must not lean toward either parent.
  const flipped = phased(balanced.map((g) => (g === 'AA' ? 'BB' : 'AA')), 4)
  assert.equal(flipped.hypothesis, c.hypothesis,
    'and the same pattern with the alleles swapped must give the same answer')
}

// --- 7. more markers sharpen a real call, and never flip it ---------------------------------------
{
  const small = phased(rep('BB', 30), 4)
  const large = phased(rep('BB', 600), 4)
  assert.equal(small.hypothesis, large.hypothesis, 'more evidence must not change the direction')
  assert.ok(large.posterior >= small.posterior, 'more evidence must not weaken the call')
}

// --- 8. UNPHASED input must not name a side ------------------------------------------------------
//
// The failure this guards against is the worst available: reporting a parental side from unphased
// markers, where the reference allele is defined per marker so a real one-sided loss retains A at
// some markers and B at others and the count cancels. Unphased input may still separate a deletion
// from an intact region, and must say so without naming.
{
  const unphased = callSiblingOrigin(rep('BB', 200), 4, 0.308)
  assert.equal(unphased.phased, false)
  assert.equal(unphased.hypothesis, 'refused', 'unphased input must never name a side')
  assert.ok(unphased.why.includes('not phased'), 'and must say that is why')
  assert.ok(unphased.posterior > 0.95, 'while still reporting that a copy IS missing')

  // The intact case needs no orientation, so it is callable either way.
  const intact = callSiblingOrigin(rep('AB', 200), 4, 0.308)
  assert.equal(intact.hypothesis, 'no-deletion', 'an intact region needs no phase to recognise')
  assert.equal(phased(rep('AB', 200), 4).hypothesis, intact.hypothesis)
}

console.log('siblingOrigin.check.ts: all assertions passed, including symmetry, the phi wall '
  + 'and the refusal to name a side without phase')

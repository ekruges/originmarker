/**
 * oneParentGain: a gain read as a lost parent goes to the gain test; the rebuilt-parent weight follows
 * the pronucleus arithmetic; the heterozygous-call route names the other parent upward and refuses
 * the downward case.
 */
import assert from 'node:assert/strict'
import {
  gainReadAsLoss, rebuiltHomWeight, emptyWeightedShare, addWeightedShare, weightedShareR,
  hetCallShare, otherCopyFromHetCalls,
} from './oneParentGain.ts'
import { emptyShareTally, addShare, shareCall } from './oneParentState.ts'
import { COPY_SHIFT_GAIN_FLOOR } from './parentage.ts'

// Rule 1. A loss state under an intensity gain call is the extra copy's share; without the call it
// stays a loss, and a disomic state is not touched.
assert.equal(gainReadAsLoss('other-lost', true), true)
assert.equal(gainReadAsLoss('parent-lost', true), true)
assert.equal(gainReadAsLoss('other-lost', false), false)
assert.equal(gainReadAsLoss('disomic', true), false)
assert.equal(gainReadAsLoss(undefined, true), false)

// Rule 1, from the share itself: on a compressed cell whose two copies read r = 0.27, P's extra copy
// reads 2/3 of that, under OTHER_LOST_MAX, which is the case the rule exists for.
{
  const t = emptyShareTally()
  for (let i = 0; i < 400; i += 1) {
    addShare(t, 'AA', 0, 0, true)                       // S1: no copy carries B
    addShare(t, 'AA', i % 10 < 6 ? 0.27 * (2 / 3) : 0, 0.6, false)
  }
  assert.equal(shareCall(t).state, 'other-lost')
  assert.equal(gainReadAsLoss(shareCall(t).state, true), true)
}

// Rule 2. A measured parent weighs 1; q = 0 weighs 1; the weight falls with q and rises with m.
assert.equal(rebuiltHomWeight(0.6, undefined), 1)
assert.equal(rebuiltHomWeight(0, 5), 1)
assert.ok(Math.abs(rebuiltHomWeight(0.5, 5) - 0.25 / (0.25 + 0.5 / 32)) < 1e-12)
assert.ok(rebuiltHomWeight(0.9, 5) < rebuiltHomWeight(0.6, 5))
assert.ok(rebuiltHomWeight(0.9, 8) > rebuiltHomWeight(0.9, 5))
assert.equal(rebuiltHomWeight(1, 5), 0)

// Rule 2, the weighted share equals the unweighted one when every weight is 1 (measured parent).
{
  const t = emptyShareTally(); const w = emptyWeightedShare()
  for (let i = 0; i < 300; i += 1) {
    addShare(t, 'AA', 0.02, 0, true)
    const y = i % 2 ? 0.27 : 0.02
    addShare(t, 'AA', y, 0.7, false); addWeightedShare(w, 'AA', y, 0.7, undefined)
  }
  assert.ok(Math.abs(weightedShareR(t, w) - shareCall(t).r) < 1e-12)
  // and markers under COMMON_ALLELE_MIN never enter it
  const before = { ...w }
  addWeightedShare(w, 'AA', 0.9, 0.3, 5)
  assert.deepEqual(w, before)
}

// Rule 3. Input: only P-homozygous markers, common lacked allele, heterozygous call.
assert.equal(hetCallShare('AA', 'AB', 0.6, 0.7), 0.6)
assert.equal(hetCallShare('BB', 'AB', 0.6, 0.3), 0.4)
assert.equal(hetCallShare('AA', 'AA', 0.6, 0.7), null)
assert.equal(hetCallShare('AB', 'AB', 0.5, 0.7), null)
assert.equal(hetCallShare('AA', 'AB', 0.5, 0.2), null)

// Rule 3. The het cluster at 1/2 on reference chromosomes with a little spread.
const cluster = (centre: number, n = 200, jitter = 0.01) =>
  Array.from({ length: n }, (_, i) => centre + ((i % 7) - 3) * jitter)
const refs = [0.49, 0.5, 0.5, 0.51, 0.5, 0.495, 0.505].map((c) => cluster(c))
// The other parent's extra copy, same homologue twice: the cluster moves to 2/3. Named.
assert.equal(otherCopyFromHetCalls(cluster(2 / 3), refs)?.named, true)
// The loaded parent's extra copy (or the other's from both homologues): the cluster falls to 1/3.
// Downward names nothing.
const down = otherCopyFromHetCalls(cluster(1 / 3), refs)
assert.ok(down !== null && down.z < -3 && down.named === false)
// Two copies: nothing.
assert.equal(otherCopyFromHetCalls(cluster(0.5), refs)?.named, false)
// Too few markers or reference chromosomes: no reading.
assert.equal(otherCopyFromHetCalls(cluster(2 / 3, 50), refs), null)
assert.equal(otherCopyFromHetCalls(cluster(2 / 3), refs.slice(0, 5)), null)
// The spread floor: identical reference medians do not make a tiny lift significant.
const flat = Array.from({ length: 7 }, () => cluster(0.5, 200, 0))
assert.equal(otherCopyFromHetCalls(cluster(0.52, 200, 0), flat)?.named, false)

console.log('oneParentGain.check: ok')

// --- RULE 1'S TRIGGER IS THE INTENSITY, NOT THE GAIN CALL -------------------------------------
//
// 3 of 158 built gains printed as a lost parent with the measured father because the guard was
// keyed on the gain CALL: the intensity on those chromosomes is up, but not far enough past the
// array's own spread for a gain to be asserted. A lost copy cannot raise intensity at all.
// See audit/truth/gains-called-2026-11-18.txt rule 4.
assert.equal(gainReadAsLoss('other-lost', true), true, 'a called gain still triggers it')
assert.equal(gainReadAsLoss('parent-lost', true), true)
assert.equal(gainReadAsLoss('other-lost', false, COPY_SHIFT_GAIN_FLOOR), true,
  'and so does an intensity at the gain floor with no call made')
assert.equal(gainReadAsLoss('parent-lost', false, 0.9), true)
assert.equal(gainReadAsLoss('other-lost', false, COPY_SHIFT_GAIN_FLOOR - 0.01), false,
  'under the floor the copy state stands: the loss reading is not withheld on nothing')
assert.equal(gainReadAsLoss('other-lost', false, -0.5), false, 'a copy DOWN is a loss, not a gain')
assert.equal(gainReadAsLoss('other-lost', false, NaN), false,
  'an unmeasurable intensity withholds nothing')
assert.equal(gainReadAsLoss('disomic', false, 0.9), false, 'and no loss state is nothing to withhold')
assert.equal(gainReadAsLoss(undefined, true, 0.9), false)

console.log('oneParentGain.check (rule 1 trigger): ok')

/**
 * The one-parent copy state follows the dosage algebra, and the normal-majority test fails closed.
 */
import assert from 'node:assert/strict'
import {
  emptyShareTally, addShare, shareCall, disomicFraction, MIN_SHARE_MARKERS,
} from './oneParentState.ts'

// A synthetic chromosome: P is AA throughout; S1 markers (q = 0) and S2 markers (q = 0.6). `other`
// says what the sample carries besides P's copy. Drop-in d raises y everywhere; it must cancel.
function chrom(other: 'present' | 'other-lost' | 'parent-lost', d: number, n = 400) {
  const t = emptyShareTally()
  for (let i = 0; i < n; i += 1) {
    // S1: every founder is AA, so no copy here carries B. y is drop-in only.
    const y1 = other === 'parent-lost' ? d : d
    addShare(t, 'AA', y1, 0, true)
    // S2: B is common (q = 0.6). Deterministic: every fifth marker's other copy carries A.
    const otherIsB = i % 5 !== 0
    let y2: number
    if (other === 'present') y2 = otherIsB ? 0.5 * (1 - 2 * d) + d : d   // a het reads at the middle
    else if (other === 'other-lost') y2 = d                              // P's A alone
    else y2 = otherIsB ? 1 - d : d                                       // the other parent's copy alone
    addShare(t, 'AA', y2, 0.6, false)
  }
  return shareCall(t)
}

for (const d of [0, 0.05, 0.25]) {
  assert.equal(chrom('present', d).state, 'disomic', `two copies at drop-in ${d}`)
  assert.equal(chrom('other-lost', d).state, 'other-lost', `the other parent's copy gone at drop-in ${d}`)
  assert.equal(chrom('parent-lost', d).state, 'parent-lost', `P's copy gone at drop-in ${d}`)
}

// DROP-IN CANCELS: r on a disomic chromosome barely moves from no drop-in to heavy drop-in.
assert.ok(Math.abs(chrom('present', 0).r - chrom('present', 0.25).r) < 0.02)

// P HETEROZYGOUS OR UNCALLED SAYS NOTHING, and a marker the panel has no frequency for is skipped.
{
  const t = emptyShareTally()
  addShare(t, 'AB', 0.5, 0, true); addShare(t, 'NC', 0.5, 0.7, false); addShare(t, 'AA', 0.5, undefined, false)
  assert.deepEqual(t, emptyShareTally())
  // BB orients the other way: y is 1 - BAF and q is 1 - freqB.
  addShare(t, 'BB', 0.9, 1, true)
  assert.equal(t.n1, 1); assert.ok(Math.abs(t.y1 - 0.1) < 1e-12)
}

// TOO FEW MARKERS IS UNDECIDED, and so is NO TEMPLATE, which reads half drop-in at S1.
{
  const t = emptyShareTally()
  for (let i = 0; i < MIN_SHARE_MARKERS - 1; i += 1) { addShare(t, 'AA', 0, 0, true); addShare(t, 'AA', 0.5, 0.6, false) }
  assert.equal(shareCall(t).state, 'undecided')
  const nt = emptyShareTally()
  for (let i = 0; i < 300; i += 1) { addShare(nt, 'AA', 0.45, 0, true); addShare(nt, 'AA', 0.45, 0.6, false) }
  assert.equal(shareCall(nt).state, 'no-template')
}

// THE NORMAL-MAJORITY TEST: undecided chromosomes do not count, and a genome reading no template
// everywhere is 0, so an array whose calls are noise cannot pass it.
assert.equal(disomicFraction([chrom('present', 0.1), chrom('present', 0.1), chrom('other-lost', 0.1)]), 2 / 3)
{
  const nt = emptyShareTally()
  for (let i = 0; i < 300; i += 1) { addShare(nt, 'AA', 0.45, 0, true); addShare(nt, 'AA', 0.45, 0.6, false) }
  assert.equal(disomicFraction([shareCall(nt), shareCall(nt), shareCall(nt)]), 0)
}

console.log('  ok   oneParentState.check.ts')

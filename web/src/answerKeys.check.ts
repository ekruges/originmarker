/**
 * The answer keys count what they say they count, and the allele state follows the cell's own keys.
 */
import assert from 'node:assert/strict'
import {
  emptySameHom, addSameHom, dropIn, absentAlleleFloor, alleleState, ALLELE_STATE_EVIDENCE,
} from './answerKeys.ts'

// 1. ONLY MARKERS WHERE BOTH PARENTS ARE THE SAME HOMOZYGOTE COUNT.
{
  const t = emptySameHom()
  addSameHom(t, 'AA', 'AA', 'AA')   // right
  addSameHom(t, 'AB', 'AA', 'AA')   // drop-in
  addSameHom(t, 'BB', 'AA', 'AA')   // error
  addSameHom(t, 'BB', 'BB', 'BB')   // right
  addSameHom(t, 'AA', 'BB', 'BB')   // error
  addSameHom(t, 'NC', 'AA', 'AA')   // not called: not a reading
  addSameHom(t, 'AB', 'AA', 'BB')   // opposite homozygotes: informative, not an answer key
  addSameHom(t, 'AB', 'AB', 'AB')   // heterozygous parents: anything is possible
  addSameHom(t, 'AA', 'AA', 'NC')   // one parent uncalled: unknown answer
  assert.deepEqual(t, { n: 5, het: 1, wrong: 2 })
  assert.equal(dropIn(t), 1 / 5)
  assert.equal(absentAlleleFloor(t), 2 / 5)
  assert.ok(Number.isNaN(dropIn(emptySameHom())), 'no markers is no measurement, not zero')
}

// 2. THE STATES on a cell with drop-in: a 2% floor, two-copy chromosomes reading 15% one-parent.
{
  const s = (x: number, y: number, n: number) => alleleState(x, y, n, 0.02, 0.15).state
  assert.equal(s(15, 700, 1000), 'lost', 'the rarer allele at the floor, the other far over two copies')
  assert.equal(s(3, 700, 1000), 'lost', 'under the floor is a loss too')
  assert.equal(s(60, 650, 1000), 'partial', 'dominant far over two copies, rarer well over the floor')
  assert.equal(s(140, 160, 1000), 'none', 'both sides at the two-copy level')
  assert.equal(s(40, 180, 1000), 'none',
    'an extra copy barely moves the dominant side, so it is not read as a partial loss')
}

// 3. ERROR-FREE MATERIAL, which is where a rule on the rarer side alone could not fire: two-copy
// chromosomes carry no one-parent calls, so the floor and the two-copy level are both zero.
assert.equal(alleleState(0, 980, 1000, 0, 0).state, 'lost')
assert.equal(alleleState(0, 0, 1000, 0, 0).state, 'none')

// 4. NOTHING MEASURED, NOTHING CALLED.
assert.equal(alleleState(10, 700, 1000, NaN, 0.15).state, 'undecided')
assert.equal(alleleState(10, 700, 1000, 0.02, NaN).state, 'undecided')

// 5. THE CASE THAT SETS THE RULE'S SHAPE. One trophectoderm biopsy, published 45,XY,-6, carries its
// chr6 monosomy at 6.0% rarer and 69.6% dominant over 1,123 markers and a chr21 published normal at
// 6.8% and 70.7% over 263. Both are partial: the data do not separate them, and a rule that named a
// parent on one would name it on both.
assert.equal(alleleState(67, 782, 1123, 0.02, 0.15).state, 'partial')
assert.equal(alleleState(18, 186, 263, 0.02, 0.15).state, 'partial')

// 6. A LIKELIHOOD RATIO SCALES WITH THE EVIDENCE: a dominant 40% over 1,000 markers is decisive,
// over ten it is four markers and proves nothing.
assert.equal(alleleState(5, 400, 1000, 0.02, 0.15).state, 'lost')
assert.equal(alleleState(0, 4, 10, 0.02, 0.15).state, 'none')
assert.ok(Math.abs(ALLELE_STATE_EVIDENCE - Math.log(1000)) < 1e-12)

console.log('  ok   answerKeys.check.ts')

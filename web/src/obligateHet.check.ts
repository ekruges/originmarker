// Self-check for the obligate-heterozygote ploidy discriminator.
// Run: node src/obligateHet.check.ts
//
// The figures are the measured ones an external review produced on 135 GSE148488 arrays. They are
// pinned here because the whole point of this module is that it separates where the BAF band does
// not, and a change that closes the gap would be silent otherwise.
import assert from 'node:assert/strict'
import {
  HAPLOID_MAX, DIPLOID_MIN, ONE_PARENT_HAPLOID_MAX, ONE_PARENT_DIPLOID_MIN, MIN_INFORMATIVE,
  emptyHet, addOneParent, addTwoParent, hetCall,
} from './obligateHet.ts'
import { ploidyOf } from './productTriage.ts'
import type { AB } from './informativity.ts'

/** A tally at a given heterozygous fraction over `n` informative markers. */
const at = (fraction: number, n = 1000) => {
  const t = emptyHet()
  t.informative = n
  t.het = Math.round(fraction * n)
  return t
}

// --- 1. the measured per-stage fractions land on the right side --------------------------------
{
  // Haploid meiotic products: 0.0705 maternal, 0.0734 paternal, full range 0.047-0.100.
  for (const f of [0.047, 0.0705, 0.0734, 0.100]) {
    assert.equal(hetCall(at(f), 2).ploidy, 'uniparental', `${f} is a haploid product`)
  }
  // Everything post-zygotic. The measured per-stage medians all clear the boundary.
  for (const f of [0.5244, 0.7818, 0.8950, 0.9868]) {
    assert.equal(hetCall(at(f), 2).ploidy, 'biparental', `${f} carries two parents`)
  }

  // The boundary is deliberately conservative and this is the cost of that: the LOWEST diploid
  // array measured, at 0.4213, falls in the refusal band rather than being called. That is the
  // intended trade and it is pinned so nobody "fixes" it by moving the boundary down to 0.42,
  // which would put the boundary on the most extreme observation in the set.
  assert.equal(hetCall(at(0.4213), 2).ploidy, 'uncalled',
    'the lowest observed diploid is refused, not called, and that is the conservative choice')
  assert.ok(DIPLOID_MIN > 0.4213, 'the boundary sits above it on purpose')
}

// --- 2. THE GAP, which is the reason this module exists ----------------------------------------
//
// The BAF band has no gap on post-zygotic material: 91 of 120 non-haploid arrays fall inside the
// haploid band range, and 22 blastomeres sit below the 10% diploid-exclusion floor. This
// statistic's gap is four times the width of the entire haploid range.
{
  const highestHaploid = 0.100
  const lowestDiploid = 0.4213
  // Both boundaries sit inside the empty space between the two classes, not on either edge of it.
  assert.ok(HAPLOID_MAX > highestHaploid, 'the lower boundary clears the haploid range')
  assert.ok(DIPLOID_MIN < 1 && HAPLOID_MAX < DIPLOID_MIN, 'and the two do not cross')
  assert.ok(lowestDiploid - highestHaploid > 4 * (highestHaploid - 0.047),
    'the gap is several times the spread of the class below it')

  // The failure this replaces, stated as an executable fact: a blastomere whose BAF band and
  // heterozygosity sit in the haploid range is called haploid by the old gate and biparental by
  // this one. Both readings are of the same cell.
  assert.equal(ploidyOf(0.07, 0.10), 'haploid', 'the band gate reads a WGA blastomere as haploid')
  assert.equal(hetCall(at(0.5244), 2).ploidy, 'biparental', 'this one does not')
}

// --- 3. one parent widens the boundaries and is marked provisional ------------------------------
{
  assert.ok(ONE_PARENT_HAPLOID_MAX < HAPLOID_MAX, 'one parent dilutes the diploid signal')
  assert.ok(ONE_PARENT_DIPLOID_MIN < DIPLOID_MIN)
  const one = hetCall(at(0.15), 1)
  assert.equal(one.provisional, true, 'and every call under it says so')
  assert.equal(hetCall(at(0.0705), 2).provisional, false)

  // The same number reads differently depending on how the informative set was defined, which is
  // why the flag has to travel with the call. At 0.15: a clean uniparental call with both parents
  // known, and a refusal with only one, because one parent cannot tell a diluted diploid signal
  // from a haploid's error rate at that level.
  assert.equal(hetCall(at(0.15), 2).ploidy, 'uniparental')
  assert.equal(hetCall(at(0.15), 1).ploidy, 'uncalled')

  // The boundaries are exclusive: a fraction sitting exactly on one is not called by it.
  assert.equal(hetCall(at(HAPLOID_MAX), 2).ploidy, 'uncalled')
}

// --- 4. it refuses rather than guessing --------------------------------------------------------
{
  assert.equal(hetCall(at(0.5, MIN_INFORMATIVE - 1), 2).ploidy, 'uncalled',
    'too few informative markers is a refusal')
  assert.equal(hetCall(emptyHet(), 2).ploidy, 'uncalled', 'and none at all certainly is')
  // Between the boundaries is neither, not the nearer one.
  const mid = hetCall(at(0.30), 2)
  assert.equal(mid.ploidy, 'uncalled')
  assert.ok(mid.why.includes('neither'))
}

// --- 5. what counts as informative -------------------------------------------------------------
{
  const t = emptyHet()
  // One parent: that parent homozygous. A parent heterozygote says nothing.
  addOneParent('AB', 'AB', t)
  assert.equal(t.informative, 0, 'a heterozygous parent makes no marker informative')
  addOneParent('AA', 'NC', t)
  assert.equal(t.informative, 0, 'and a sample no-call is not evidence of homozygosity')
  addOneParent('AA', 'AB', t)
  addOneParent('BB', 'BB', t)
  assert.deepEqual([t.informative, t.het], [2, 1])

  const u = emptyHet()
  // Two parents: both homozygous AND opposite. Same homozygote is not obligate.
  addTwoParent('AA', 'AA', 'AB', u)
  assert.equal(u.informative, 0, 'two parents homozygous for the SAME allele is not obligate')
  addTwoParent('AA', 'BB', 'AB', u)
  addTwoParent('BB', 'AA', 'AA', u)
  assert.deepEqual([u.informative, u.het], [2, 1])
}

console.log('obligateHet.check.ts: all assertions passed')

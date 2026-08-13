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
  // why the flag has to travel with the call. At 0.15, with BOTH parents the informative set is
  // obligate-het and a diploid must read near 0.42, so 0.15 is uniparental. With ONE parent the
  // set is diluted to markers where the unseen parent differs and a diploid reads about 0.09, so
  // the same 0.15 is comfortably biparental.
  //
  // This second line previously expected a refusal, on the assumption that one parent could not
  // separate a diluted diploid from a haploid's error rate here. Measured against a bulk parent it
  // separates by a factor of two, one-genome products at 0.043-0.059 against children at 0.089 up.
  assert.equal(hetCall(at(0.15), 2).ploidy, 'uniparental')
  assert.equal(hetCall(at(0.15), 1).ploidy, 'biparental')

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

// --- the reference parent must be bulk, measured externally on GSE19247 ------------------------
//
// The numbers this pins are the ones that would otherwise be forgotten: scored against a
// SINGLE-CELL reference, haploid sperm and biparental blastomeres overlap almost completely, so
// no boundary placed anywhere separates them. Kept here so a future tightening of the one-parent
// band cannot quietly claim to have fixed a problem that is not in the band.
{
  const SPERM_HAPLOID_RANGE = [0.0021, 0.5132]   // n=23, every one haploid by construction
  const BLASTOMERE_RANGE = [0.0177, 0.3004]      // n=28, every one biparental

  const overlaps = SPERM_HAPLOID_RANGE[0] < BLASTOMERE_RANGE[1]
    && BLASTOMERE_RANGE[0] < SPERM_HAPLOID_RANGE[1]
  assert.ok(overlaps, 'the measured single-cell-reference distributions overlap; that is the point')

  // Whatever the boundaries are, they cannot separate that material. Assert it rather than trust
  // the prose: pick the shipped one-parent band and show it misclassifies in both directions.
  const worstSperm = SPERM_HAPLOID_RANGE[1]
  const bestBlastomere = BLASTOMERE_RANGE[0]
  assert.ok(worstSperm > ONE_PARENT_DIPLOID_MIN,
    'a haploid sperm reaches the diploid boundary, so the band cannot exclude it')
  assert.ok(bestBlastomere < ONE_PARENT_HAPLOID_MAX,
    'a biparental blastomere falls under the haploid boundary, so the band cannot admit it')
}

// --- the one-parent boundaries, against a BULK reference where they do work ---------------------
//
// The same statistic separates cleanly once the reference is bulk. These are the measured values,
// pinned so the boundaries cannot drift back to the scaled guesses that called every biparental
// child of a single genotyped parent uniparental.
{
  const ONE_GENOME = [0.0428, 0.0514, 0.0570, 0.0587]   // products of a bulk-genotyped parent
  const BIPARENTAL = [0.0894, 0.0973]                   // his children, two parental genomes

  for (const v of ONE_GENOME) {
    assert.equal(hetCall({ informative: 50_000, het: Math.round(50_000 * v) }, 1).ploidy,
      'uniparental', `${v} is a one-genome product and must call uniparental`)
  }
  for (const v of BIPARENTAL) {
    assert.equal(hetCall({ informative: 50_000, het: Math.round(50_000 * v) }, 1).ploidy,
      'biparental', `${v} is a biparental child and must call biparental`)
  }
  // The regression this replaces: at the old 0.12 bound every one of those biparental children
  // fell in the uniparental class.
  assert.ok(Math.min(...BIPARENTAL) < 0.12,
    'the old bound sat above the biparental class, which is why it was wrong')
  assert.ok(ONE_PARENT_HAPLOID_MAX < Math.min(...BIPARENTAL)
    && ONE_PARENT_DIPLOID_MIN <= Math.min(...BIPARENTAL),
    'the boundaries must bracket the measured separation')
}

console.log('obligateHet.check.ts: bulk-reference precondition pinned')

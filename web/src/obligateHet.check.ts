// Self-check for the obligate-heterozygote ploidy discriminator.
// Run: node src/obligateHet.check.ts
//
// The figures are the measured ones an external review produced on 135 GSE148488 arrays. They are
// pinned here because the whole point of this module is that it separates where the BAF band does
// not, and a change that closes the gap would be silent otherwise.
import assert from 'node:assert/strict'
import {
  informativeEnough, hetCallOneParent,
  HAPLOID_MAX, DIPLOID_MIN, ONE_PARENT_HAPLOID_MAX, ONE_PARENT_DIPLOID_MIN, MIN_INFORMATIVE,
  SURVIVOR_MIN, MIXED_SHARE, emptyHet, addOneParent, addTwoParent, hetCall,
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

// --- 6. whose allele survives names the copy that went, two parents only -----------------------
{
  // Father AA, mother BB. A cell with both copies reads AB; one carrying the father's copy alone
  // reads AA; the mother's alone, BB.
  const t = emptyHet()
  for (let i = 0; i < 300; i += 1) addTwoParent('AA', 'BB', i < 20 ? 'AB' : 'AA', t)
  assert.deepEqual([t.informative, t.het, t.onlyA, t.onlyB], [300, 20, 280, 0])
  const c = hetCall(t, 2)
  assert.equal(c.ploidy, 'uniparental')
  assert.equal(c.survivor, 'a', 'the first parent\'s allele alone makes that parent the survivor')
  assert.ok(Math.abs(c.onlyA - 280 / 300) < 1e-9 && c.onlyB === 0)

  // The measured margins on GSE148488 against Table S3: the survivor's share is 0.735 and up on
  // every lost chromosome, and neither parent's share reaches 0.627 on any intact one. The floor
  // sits inside that gap rather than on either edge of it.
  assert.ok(SURVIVOR_MIN > 0.627 && SURVIVOR_MIN < 0.735, 'the floor sits in the measured gap')

  // Uniparental with no survivor is what a chromosome the arrays are not measuring reads as, and
  // it names nobody. Here the two homozygotes split 140 to 140.
  const w = emptyHet()
  for (let i = 0; i < 300; i += 1) addTwoParent('AA', 'BB', i < 20 ? 'AB' : i < 160 ? 'AA' : 'BB', w)
  assert.equal(hetCall(w, 2).ploidy, 'uniparental')
  assert.equal(hetCall(w, 2).survivor, undefined, 'neither allele clears the floor')

  // A mixed signal sits between the boundaries, so its share floor must sit above the two-parent
  // diploid boundary or it would flag chromosomes the boundaries already call biparental.
  assert.ok(MIXED_SHARE > DIPLOID_MIN && MIXED_SHARE < SURVIVOR_MIN, 'the mixed floor sits between the boundaries and the survivor floor')

  // One parent loaded: the other allele has no owner, so nothing is attributed.
  const u = emptyHet()
  for (let i = 0; i < 300; i += 1) addOneParent('AA', 'AA', u)
  assert.deepEqual([u.onlyA, u.onlyB], [0, 0])
  assert.equal(hetCall(u, 1).survivor, undefined)
}

// --- A CHROMOSOME THAT WAS BARELY GENOTYPED MUST NOT REPORT A FRACTION --------------------------
//
// Measured on a laboratory series of 176 chromosome-observations. One chromosome was called at
// 0.175 of the markers the array holds for it and reported 55.1 percent heterozygous where its
// siblings read 0.07. Its informative count, 1,605, and its informative share of called markers,
// 0.789, are both unremarkable: the count cannot catch it and the call rate can. The series median
// call rate is 0.864, the fifth percentile 0.807, and nothing falls between 0.175 and 0.813.
{
  // The chromosome that must be refused, with its real numbers.
  assert.ok(!informativeEnough(1_605, 0.175), 'called at 0.175 of the array is not measured')
  // Its siblings, which must all be admitted.
  for (const cr of [0.813, 0.824, 0.848, 0.849, 0.869, 0.879, 0.905]) {
    assert.ok(informativeEnough(7_490, cr), `a call rate of ${cr} is fine`)
  }
  // The absolute floor still bites on its own, whatever the call rate.
  assert.ok(!informativeEnough(199, 0.95), 'under the absolute floor is refused')
  // No call rate available leaves the absolute floor as the only test, which is the old behaviour.
  assert.ok(informativeEnough(1_605, NaN), 'an unmeasurable call rate falls back to the floor')
  // AND IT MUST NOT PENALISE A SMALL CHROMOSOME. A first attempt compared the informative count to
  // the sample's median autosome and refused chr22 of a public cell for carrying 17 percent of it.
  // chr22 is small; its call rate is normal.
  assert.ok(informativeEnough(3_632, 0.90), 'chr22 is small, not unmeasured')
}

// --- THE MARKER CLASS THE ONE-PARENT CHANNEL DISCARDS, AND WHY IT IS A RATIO ---------------------
//
// Measured over 1,021 chromosome-observations of GSE148488 against the authors' published
// karyotype, with the same two-hypothesis rule applied to both classes at the shipped posterior:
//
//   parent-homozygous    one 50 of 102 = 0.490   two 874 of 917 = 0.953   overall 0.9068
//   parent-heterozygous  one 80 of 103 = 0.777   two 907 of 918 = 0.988   overall 0.9667
//
// The excluded class is better on both, and the calls on chromosomes the hand reads as normal fall
// from 43 in 918 to 11. The rate distributions OVERLAP, one contribution running to 0.4320 and two
// down to 0.1334, which is why this is a likelihood ratio and not a pair of boundaries.
//
// NOT WIRED INTO THE VERDICT. Substituting it changed calls on the synthetic fixtures in ways two
// other self-checks depend on, so the consumers need tracing first. The function and its
// measurement live here until that is done.
{
  // A chromosome carrying two contributions at a typical dropout reads about half the markers
  // heterozygous, and one contribution reads drop-in alone.
  const two = hetCallOneParent(Math.round(0.45 * 4_000), 4_000, 0.15)
  assert.equal(two.ploidy, 'biparental', two.why)
  const one = hetCallOneParent(Math.round(0.045 * 4_000), 4_000, 0.15)
  assert.equal(one.ploidy, 'uniparental', one.why)

  // IT REFUSES BETWEEN THE TWO PREDICTIONS rather than guessing, which is what keeps the false
  // calls at 11 in 918.
  const mid = hetCallOneParent(Math.round(0.20 * 300), 300, 0.15)
  assert.equal(mid.ploidy, 'uncalled', mid.why)

  // Too few markers is a refusal, and the count is in the reason so a reader can check it.
  const thin = hetCallOneParent(10, 20, 0.15)
  assert.equal(thin.ploidy, 'uncalled')
  assert.ok(thin.why.includes('20'), thin.why)

  // THE SURVIVOR COUNTS ARE NOT DEFINED ON THIS CLASS and must not read as measured zeroes: where
  // the loaded parent is heterozygous it carries both alleles, so "that parent's allele alone" has
  // no meaning. NaN rather than 0.
  assert.ok(Number.isNaN(two.onlyA) && Number.isNaN(two.onlyB))
  assert.equal(two.provisional, true, 'one parent defined the informative set')

  // AND DROPOUT MOVES THE BIPARENTAL PREDICTION, which is the whole reason this class is better:
  // at heavy dropout a biparental chromosome reads far less heterozygous and must still be called.
  // A 25 percent rate is biparental either way, because it is far closer to 0.5(1-ado) than to
  // drop-in at any dropout, but it sits ON the prediction at heavy dropout and well below it on
  // clean material, so the evidence is stronger there. That is the property, not the verdict.
  const heavy = hetCallOneParent(Math.round(0.25 * 4_000), 4_000, 0.55)
  const clean = hetCallOneParent(Math.round(0.25 * 4_000), 4_000, 0.02)
  assert.equal(heavy.ploidy, 'biparental')
  assert.equal(clean.ploidy, 'biparental')
  // The prediction it was weighed against moves with the dropout, and the reason says so, which is
  // what a reader checks by hand. 0.5(1-0.55) + drop-in is 26.8%; 0.5(1-0.02) + drop-in is 53.3%.
  assert.ok(heavy.why.includes('26.8%'), heavy.why)
  assert.ok(clean.why.includes('53.3%'), clean.why)
  // And a rate at drop-in is uniparental whatever the dropout, since one contribution predicts it.
  for (const ado of [0.02, 0.30, 0.55]) {
    assert.equal(hetCallOneParent(Math.round(0.045 * 4_000), 4_000, ado).ploidy, 'uniparental')
  }
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
    assert.equal(hetCall({ ...emptyHet(), informative: 50_000, het: Math.round(50_000 * v) }, 1).ploidy,
      'uniparental', `${v} is a one-genome product and must call uniparental`)
  }
  for (const v of BIPARENTAL) {
    assert.equal(hetCall({ ...emptyHet(), informative: 50_000, het: Math.round(50_000 * v) }, 1).ploidy,
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

/**
 * The parent-genome reading keeps its sign under any dropout, and noise is not read as a genome.
 */
import assert from 'node:assert/strict'
import {
  emptyCarryTally, addCarry, carryCall, parentGenome, RARE_MAX, PRESENT_MIN, ENRICHMENT_MIN,
  type CarryCall, type CarryTally,
} from './parentGenome.ts'

// A deterministic stream of uniform numbers, so the check is the same on every run.
let seed = 12345
const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 }

type Kind = 'child' | 'stranger' | 'haploid-stranger' | 'product' | 'noise'

/**
 * One chromosome against an adult. `f` is the frequency, in the sample's population, of the allele the
 * adult is homozygous for at rare markers; `dropout` the chance a heterozygote loses one allele. Where
 * the adult is the common homozygote the rare allele stays at the founders' 0.1: shared ancestry
 * raises an allele's frequency only where the adult is known to hold it.
 */
function tally(kind: Kind, f: number, dropout: number, n = 600, losesAdult = 0.5): CarryTally {
  const t = emptyCarryTally()
  // `losesAdult` is the chance that the allele a dropout removes is the one inherited from the adult.
  const observe = (fromAdult: string, other: string): string => {
    if (fromAdult === other) return fromAdult + fromAdult
    if (rnd() < dropout) return rnd() < losesAdult ? other + other : fromAdult + fromAdult
    return 'AB'
  }
  for (let i = 0; i < n; i += 1) {
    // Rare class: the adult is BB where B is rare (founder frequency 0.1).
    // Control class: the adult is AA at markers of the same class.
    for (const adult of ['BB', 'AA'] as const) {
      const fromPopulation = () => (rnd() < (adult === 'BB' ? f : 0.1) ? 'B' : 'A')
      let g: string
      if (kind === 'noise') g = ['AA', 'AB', 'BB'][Math.floor(rnd() * 3)]
      else if (kind === 'product') g = adult
      else if (kind === 'child') g = observe(adult[0], fromPopulation())
      else if (kind === 'haploid-stranger') { const a = fromPopulation(); g = a + a }
      else g = observe(fromPopulation(), fromPopulation())
      addCarry(t, adult, g, 0.1)
    }
  }
  return t
}
// 2,000 markers a class, a large chromosome, so that one chromosome's reading is not the draw's.
const chrom = (kind: Kind, f: number, dropout: number, n = 2000): CarryCall => carryCall(tally(kind, f, dropout, n))

// The sign holds at every dropout, which is what the opposite-homozygote rate cannot do.
for (const dropout of [0, 0.3, 0.6, 0.9, 1]) {
  assert.equal(chrom('child', 0.12, dropout).state, 'present', `a child at dropout ${dropout}`)
  assert.equal(chrom('stranger', 0.12, dropout).state, 'absent', `a stranger at dropout ${dropout}`)
}
// A stranger of the adult's own ancestry holds the allele far more often and still reads without it.
assert.equal(chrom('stranger', 0.4, 0.5).state, 'absent', 'a stranger of the same ancestry')
assert.equal(chrom('haploid-stranger', 0.12, 0).state, 'absent', "the other parent's pronucleus")
assert.ok(chrom('product', 0.12, 0).d > 0.95, "the adult's own product reads near one")
assert.equal(chrom('noise', 0.12, 0).state, 'noise', 'calls that are not a genome are not read')
assert.equal(chrom('child', 0.12, 0.5, 20).state, 'thin', 'too few markers are not read')

// Only a homozygous adult at a marker of known frequency is tallied, and only a called sample.
const t = emptyCarryTally()
addCarry(t, 'AB', 'AA', 0.1); addCarry(t, 'NC', 'AA', 0.1); addCarry(t, 'BB', 'NC', 0.1)
addCarry(t, 'BB', 'BB', undefined); addCarry(t, 'BB', 'BB', 0.5)
assert.deepEqual(t, emptyCarryTally())
addCarry(t, 'BB', 'BB', RARE_MAX); addCarry(t, 'AA', 'AB', 1 - RARE_MAX); addCarry(t, 'AA', 'BB', 0.05)
assert.deepEqual(t, { same: 1, het: 1, opp: 0, ctlSame: 0, ctlHet: 0, ctlOpp: 1 })

// The genome verdict.
const many = (kind: Kind, k: number, f = 0.12) => Array.from({ length: k }, () => tally(kind, f, 0.5))
assert.equal(parentGenome(many('child', 22)).verdict, 'carries')
assert.equal(parentGenome(many('stranger', 22)).verdict, 'none')
assert.equal(parentGenome(many('noise', 22)).verdict, 'unshown')
// A child that lost this parent's copy of most chromosomes still carries the genome on the rest.
assert.equal(parentGenome([...many('child', PRESENT_MIN), ...many('stranger', 19)]).verdict, 'carries')
assert.equal(parentGenome([...many('child', PRESENT_MIN - 1), ...many('stranger', 20)]).verdict, 'unshown')


// A child whose dropout takes the adult's allele four times in five reads below zero on every
// chromosome, and its heterozygotes still say whose it is. A stranger's never do, of any ancestry.
const lopsided = Array.from({ length: 22 }, () => tally('child', 0.12, 0.5, 600, 0.8))
assert.ok(lopsided.every((x) => carryCall(x).state === 'absent'), 'D alone reads the child as a stranger')
assert.equal(parentGenome(lopsided).verdict, 'carries')
assert.ok(parentGenome(lopsided).enrichment >= ENRICHMENT_MIN)
assert.ok(parentGenome(many('stranger', 22)).enrichment < 3)
assert.ok(parentGenome(many('stranger', 22, 0.4)).enrichment < ENRICHMENT_MIN, 'a stranger of the same ancestry')
assert.equal(parentGenome(many('stranger', 22, 0.4)).verdict, 'none')

// A poor array of a true child: 60 markers a chromosome and nine heterozygotes in ten dropped. Few
// chromosomes clear the bar alone; none reads without the genome, and the genome as a whole carries it.
// A stranger on the same thin chromosomes reads without it throughout and is not admitted.
const thinChild = Array.from({ length: 22 }, () => tally('child', 0.12, 0.9, 60))
assert.ok(thinChild.filter((x) => carryCall(x).state === 'present').length < PRESENT_MIN)
assert.ok(thinChild.every((x) => carryCall(x).state !== 'absent'))
assert.equal(parentGenome(thinChild).verdict, 'carries')
assert.notEqual(parentGenome(Array.from({ length: 22 }, () => tally('stranger', 0.12, 0.9, 60))).verdict, 'carries')
assert.notEqual(parentGenome(Array.from({ length: 22 }, () => tally('stranger', 0.4, 0.9, 60))).verdict, 'carries')

console.log('parentGenome.check ok')

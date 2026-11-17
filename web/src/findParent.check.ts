// Run: node --experimental-strip-types src/findParent.check.ts
import assert from 'node:assert/strict'
import {
  addCandidateRow, emptyCandidateTallies, candidateGenomes, carriers, findParent, fatherSearchText,
  FATHER_SEARCH_LIMIT,
} from './findParent.ts'
import type { ParentGenome } from './parentGenome.ts'

let x = 20261116
const rnd = (): number => { x ^= x << 13; x >>>= 0; x ^= x >> 17; x ^= x << 5; x >>>= 0; return x / 2 ** 32 }
const MARKERS = 5000
type Person = Map<string, string>
const ids = Array.from({ length: 22 }, (_, c) => Array.from({ length: MARKERS }, (_, i) => [`${c + 1}`, `p${c + 1}-${i}`]))
  .flat()
// Rare alleles at founder frequencies across the parent check's range, on either allele.
const freqB = new Map(ids.map(([, id], i) => {
  const q = 0.06 + 0.12 * rnd()
  return [id, i % 2 ? q : 1 - q]
}))
const allele = (id: string): 0 | 1 => (rnd() < freqB.get(id)! ? 1 : 0)
const gt = (a: number, b: number): string => (a + b === 0 ? 'AA' : a + b === 2 ? 'BB' : 'AB')
const adult = (): Person => new Map(ids.map(([, id]) => [id, gt(allele(id), allele(id))]))
const from = (p: Person, id: string): number => {
  const g = p.get(id)
  return g === 'AA' ? 0 : g === 'BB' ? 1 : rnd() < 0.5 ? 0 : 1
}
const child = (father: Person, mother: Person): Person =>
  new Map(ids.map(([, id]) => [id, gt(from(father, id), from(mother, id))]))
const pronucleus = (mother: Person): Person => new Map(ids.map(([, id]) => {
  const a = from(mother, id)
  return [id, gt(a, a)]
}))
const read = (sample: Person, candidates: Person[]): ParentGenome[] => {
  const t = emptyCandidateTallies(candidates.length)
  for (const [chrom, probesetId] of ids) {
    addCandidateRow(t, candidates, { chrom, probesetId, genotype: sample.get(probesetId)! }, freqB)
  }
  addCandidateRow(t, candidates, { chrom: 'X', probesetId: ids[0][1], genotype: 'BB' }, freqB)
  addCandidateRow(t, candidates, { chrom: '1', probesetId: 'no-frequency', genotype: 'BB' }, freqB)
  return candidateGenomes(t)
}

const [father, mother, stranger, other] = [adult(), adult(), adult(), adult()]
const cells = [child(father, mother), child(father, mother), pronucleus(mother)]

// One cell: the father alone carries, a stranger and the mother's pronucleus do not.
const one = read(cells[0], [stranger, father, other])
assert.deepEqual(one.map((g) => g.verdict), ['none', 'carries', 'none'])
assert.deepEqual(carriers(one), [1])
assert.deepEqual(read(cells[2], [stranger, father]).map((g) => g.verdict), ['none', 'none'])

// A run: named, with the pronucleus abstaining.
const run = findParent(cells.map((c) => read(c, [stranger, father, other])), 3)
assert.equal(run.outcome, 'named'); assert.equal(run.named, 1)
assert.deepEqual(run.counts.map((c) => c.carries), [0, 2, 0])
assert.equal(run.carryNone, 1); assert.equal(run.carrySeveral, 0)
const text = fatherSearchText(run, ['s.txt', 'f.txt', 'b.txt'])
assert.match(text.line, /^f\.txt is the father: its genome is in 2 of 3 arrays/)
assert.match(text.abstain, /^1 of 3 arrays carry none/)
assert.equal(text.limit, FATHER_SEARCH_LIMIT)
assert.deepEqual(text.rows[1], { name: 'f.txt', carries: 2, none: 1, unshown: 0 })

// The father left out: none named, every cell.
const out = findParent(cells.map((c) => read(c, [stranger, other])), 2)
assert.equal(out.outcome, 'none'); assert.equal(out.named, null); assert.equal(out.carryNone, 3)
assert.match(fatherSearchText(out, ['s', 'b']).line, /^None of the 2 candidates is shown to be the father/)

// Two candidates carry, and neither is named. The mother is the sure case. Close relatives share
// long stretches, which markers drawn one by one do not model, so that case is not built here.
const both = findParent(cells.map((c) => read(c, [father, mother, stranger])), 3)
assert.equal(both.outcome, 'several'); assert.equal(both.named, null); assert.equal(both.carrySeveral, 2)
assert.deepEqual(both.counts.map((c) => c.carries), [2, 3, 0])
assert.match(fatherSearchText(both, ['f', 'm', 's']).line, /^No father is named: f is in 2, m is in 3 of 3 arrays/)

console.log('findParent.check ok')

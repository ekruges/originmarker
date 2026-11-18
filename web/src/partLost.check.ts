/**
 * A chromosome called lost whole is read again in pieces: a deep loss over a third of it is reported
 * as that piece, a chromosome lost end to end stays a whole-chromosome call, a normal one prints
 * nothing, and a deeper piece inside a chromosome down one copy is reported beside the whole call.
 */
import assert from 'node:assert/strict'
import { readLostChromosome, wholeTestChanges, type CopyNull } from './partLost.ts'

let seed = 7
const rand = (): number => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 }
const gauss = (): number => Math.sqrt(-2 * Math.log(rand() || 1e-12)) * Math.cos(2 * Math.PI * rand())

/** 9,000 markers over 90 Mb. `level(pos)` is the true log R; a level at or under -1.5 is no copy left. */
function chromosome(level: (pos: number) => number) {
  return Array.from({ length: 9000 }, (_, i) => {
    const pos = i * 10_000 + 1
    const l = level(pos)
    return { chrom: '5', pos, called: l > -1.5 ? rand() < 0.9 : rand() < 0.05, log2R: l + 0.25 * gauss() }
  })
}
const n: CopyNull = { centre: 0, scale: 0.07, units: 20, floor: 0.40, callMedian: 0.9 }
const NULL_NO_CALL = 0.1

// A deep loss over the first 30 Mb, the rest at two copies: the piece, not the whole chromosome.
{
  const r = readLostChromosome(chromosome((p) => (p < 30e6 ? -2 : 0)), NULL_NO_CALL, 0, n)
  assert.equal(r.whole, false, 'the rest reads at the cell level, so the whole call goes')
  assert.equal(r.pieces.length >= 1, true)
  const end = Math.max(...r.pieces.map((s) => s.endBp))
  assert.ok(Math.abs(end - 30e6) < 1e6, `the piece ends at the cut (${end})`)
  assert.ok(Math.abs(r.restLrr) < 0.1)
}

// Lost end to end at one copy: no piece is deep enough, the whole call stands.
{
  const r = readLostChromosome(chromosome(() => -0.6), NULL_NO_CALL, 0, n)
  assert.equal(r.whole, true)
  assert.equal(r.pieces.length, 0)
}

// No copy left end to end: the pieces cover it, nothing is left to read, the whole call stands alone.
{
  const r = readLostChromosome(chromosome(() => -2), NULL_NO_CALL, 0, n)
  assert.equal(r.whole, true)
  assert.equal(r.pieces.length, 0, 'the whole loss is not described twice')
}

// Down one copy, and no copy left on the first 30 Mb: the whole call stands and the deeper piece is
// reported beside it.
{
  const r = readLostChromosome(chromosome((p) => (p < 30e6 ? -2 : -0.6)), NULL_NO_CALL, 0, n)
  assert.equal(r.whole, true)
  assert.equal(r.pieces.length >= 1, true, 'a further copy lost is a second event')
}

// Down one copy on a noisy array, where the rest fails the z: still not at the cell's level, so the
// whole call stands and the deeper piece is reported beside it.
{
  const wide: CopyNull = { ...n, scale: 0.15 }
  const r = readLostChromosome(chromosome((p) => (p > 60e6 ? -2 : -0.6)), NULL_NO_CALL, 0, wide)
  assert.equal(r.whole, true, 'a rest a copy down is not the cell level, significant or not')
}

// Pieces over all but a sliver near the end are the whole loss: whole call, no piece.
{
  const r = readLostChromosome(chromosome((p) => (p < 88e6 ? -2 : -0.9)), NULL_NO_CALL, 0, n)
  assert.equal(r.whole, true)
  assert.equal(r.pieces.length, 0)
}

// A normal chromosome: the whole test reads no change and the scan finds no piece.
{
  const ms = chromosome(() => 0)
  assert.equal(wholeTestChanges(ms, n), false)
  assert.equal(readLostChromosome(ms, NULL_NO_CALL, 0, n).pieces.length, 0)
}

console.log('partLost ok')

/**
 * motherFromRun: draws of one mother group and a stranger's do not; a noisy draw and one with no
 * mother's genome are not draws; the vote follows each draw's own drop-in; fewer than 2 draws do not
 * run; the route's reference honours the withheld cell, counts a lost state at the cell's intensity
 * level as both parents, and reads intensity without an allele-marker floor.
 */
import assert from 'node:assert/strict'
import {
  motherFromRun, oppositeHomPosterior, routeGains, sharedCopyShare, showsOf, ratesOfShows,
  SAME_MOTHER_MEASURED, MIN_DRAWS, type RunMarker, type RouteChrom,
} from './motherFromRun.ts'

let seed = 11
const rand = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 }

// Six chromosomes of 60 Mb, a marker every 5 kb, the loaded parent homozygous AA everywhere.
const markers: RunMarker[] = []
for (let c = 1; c <= 6; c += 1) {
  for (let pos = 0; pos < 60e6; pos += 5e3) markers.push({ id: `m${c}:${pos}`, chrom: String(c), pos, parent: 'AA', q: 0.02 + 0.96 * rand() })
}
type Mother = [boolean[], boolean[]]
const mother = (): Mother => [markers.map((m) => rand() < m.q), markers.map((m) => rand() < m.q)]
/** One copy of hers per chromosome, switching at one crossover. */
const copyOf = (h: Mother) => {
  const cross = new Map<string, [number, boolean]>()
  return markers.map((m, i) => {
    if (!cross.has(m.chrom)) cross.set(m.chrom, [rand() * 60e6, rand() < 0.5])
    const [at, first] = cross.get(m.chrom)!
    return h[(m.pos < at) === first ? 0 : 1][i]
  })
}
/** An embryo's calls: the copy from her and the loaded parent's A; drop-in b, dropout d. */
const embryo = (carried: boolean[], b = 0.01, d = 0.05) => {
  const m = new Map<string, string>()
  markers.forEach((x, i) => m.set(x.id, carried[i] ? (rand() < d ? 'AA' : 'AB') : rand() < b ? 'AB' : 'AA'))
  return m
}
const M = mother(); const S = mother()
const e1 = embryo(copyOf(M)); const e2 = embryo(copyOf(M)); const e3 = embryo(copyOf(M))
const stranger = embryo(copyOf(S))
const noisy = embryo(copyOf(M), 0.2)
const paternalOnly = embryo(markers.map(() => false))
const baf = (calls: Map<string, string>) => (id: string) => (calls.get(id) === 'AB' ? 0.5 : 0)

// Same mother: draws of M share her copies on about half the windows, M and S on none.
{
  const s = [e1, e2, stranger].map((c) => showsOf(markers, c))
  const r = s.map((x) => ratesOfShows(markers, x, true)!)
  const sib = sharedCopyShare(markers, s[0], s[1], r[0], r[1])
  const str = sharedCopyShare(markers, s[0], s[2], r[0], r[2])
  assert.ok(sib > 0.25 && sib < 0.85, `siblings share ${sib}`)
  assert.ok(str < SAME_MOTHER_MEASURED / 2, `strangers share ${str}`)
}

// The run: two clean draws of M vote; the stranger, the noisy draw and the paternal-only array do not.
{
  const run = motherFromRun({ markers, sample: e1, baf: baf(e1), others: [e2, e3, stranger, noisy, paternalOnly],
    parentKey: {}, rebuilt: false })
  assert.equal(run.status, 'used')
  assert.equal(run.draws, 2)
  // the paternal-only array shows no other parent at all and has no rate fit, so it is in no unit
  assert.equal(run.units, 5)
  assert.equal(run.readable, 4)
  // where she is the opposite homozygote the sample, carrying B from her and A from the loaded parent,
  // reads the loaded parent's share 1/2
  for (const [, v] of run.shares) assert.ok(Math.abs(v.share - 0.5) < 1e-9 && v.weight > 0)
}

// Fewer than MIN_DRAWS draws: the run is not used.
{
  const run = motherFromRun({ markers, sample: e1, baf: baf(e1), others: [e2, stranger], parentKey: {}, rebuilt: false })
  assert.equal(MIN_DRAWS, 2)
  assert.equal(run.status, 'too-few-draws')
  assert.equal(run.draws, 1)
  assert.equal(run.shares.size, 0)
}

// The vote. At the clean limit (b small, a = 1) every draw showing gives experiment 3's weight,
// q^2 / (q^2 + 2q(1-q) 2^-k); a draw with drop-in argues less, so the same calls carry less weight.
{
  const one: RunMarker[] = [{ id: 'x', chrom: '1', pos: 0, parent: 'AA', q: 0.6 }]
  const shows = Int8Array.of(1)
  const clean = { b: 0.005, a: 0.995, h: NaN, inflation: 1 }
  const w2 = oppositeHomPosterior(one, [{ s: shows, r: clean }, { s: shows, r: clean }])[0]
  const q = 0.6; const k = 2
  const limit = q * q / (q * q + 2 * q * (1 - q) * 2 ** -k)
  assert.ok(Math.abs(w2 - limit) < 0.01, `clean vote ${w2} against ${limit}`)
  const dirty = { b: 0.2, a: 0.6, h: NaN, inflation: 1 }
  const wd = oppositeHomPosterior(one, [{ s: shows, r: dirty }, { s: shows, r: dirty }])[0]
  assert.ok(wd < w2, `noisy vote ${wd} under clean ${w2}`)
  // one clean draw that does not show the allele rules her out as the opposite homozygote
  const no = oppositeHomPosterior(one, [{ s: shows, r: clean }, { s: Int8Array.of(0), r: clean }])[0]
  assert.equal(no, 0)
  // a rebuilt parent's marker also carries rule 2's weight
  const wr = oppositeHomPosterior(one, [{ s: shows, r: clean }, { s: shows, r: clean }], 5)[0]
  assert.ok(wr < w2 && wr > 0)
}

// The route's reference.
{
  const base = (chrom: string, share: number, lrrShift: number, extra: Partial<RouteChrom> = {}): RouteChrom =>
    ({ chrom, share, weight: 500, state: 'disomic', measured: true, intensityGain: false, lrrShift, ...extra })
  const ref = ['1', '2', '3', '4', '5', '6', '7', '8'].map((c, i) => base(c, 0.5 + (i % 3 - 1) * 0.01, (i % 3 - 1) * 0.03))
  const target = base('9', 0.75, 0.10)
  assert.deepEqual(routeGains([...ref, target], new Set()).map((g) => [g.chrom, g.loadedParent]), [['9', true]])
  // (a) a cell whose allele gains are withheld for want of a two-copy majority: nothing
  assert.deepEqual(routeGains([...ref, target], new Set(), true), [])
  // a chromosome the first build named is not tested again
  assert.deepEqual(routeGains([...ref, target], new Set(['9'])), [])
  // (b) a lost state with intensity at or above the cell's level reads both parents: four such chromosomes
  // widen the intensity spread and the target no longer clears it; below the level they stay out
  const up = [0.10, 0.12, 0.14, 0.16].map((l, i) => base(String(10 + i), 0.5, l, { state: 'parent-lost' }))
  assert.deepEqual(routeGains([...ref, ...up, target], new Set()).map((g) => g.chrom), [])
  const down = up.map((c) => ({ ...c, lrrShift: -c.lrrShift }))
  assert.deepEqual(routeGains([...ref, ...down, target], new Set()).map((g) => g.chrom), ['9'])
  // (c) chromosomes under the share's weight floor still count in the intensity reference, not the share's
  const thin = up.map((c) => ({ ...c, state: 'disomic', weight: 50, share: 0.9 }))
  assert.deepEqual(routeGains([...ref, ...thin, target], new Set()).map((g) => g.chrom), [])
  assert.equal(routeGains([...ref, ...thin, base('9', 0.75, 0.9)], new Set())[0].references, ref.length)
  // a lost state is tested only where the intensity scan calls the gain
  assert.deepEqual(routeGains([...ref, { ...target, state: 'other-lost' }], new Set()), [])
  assert.deepEqual(routeGains([...ref, { ...target, state: 'other-lost', intensityGain: true }], new Set()).map((g) => g.chrom), ['9'])
}

console.log('motherFromRun: all checks passed')

/**
 * The store: a value comes back as it went in, a key moves with every input byte, and an input not
 * on the open list is refused before it is read.
 *
 * Run: node --experimental-strip-types audit/gate/store.check.ts
 */
import assert from 'node:assert/strict'
import { mkdtempSync, writeFileSync, mkdirSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { createHash } from 'node:crypto'
import { serialize, deserialize } from 'node:v8'
import { gzipSync } from 'node:zlib'
import { pack, unpack, openStore, only, treeKey } from './store.ts'

const roundTrip = (x: unknown) => unpack(deserialize(serialize(pack(x))))

// Values: long integer arrays held as deltas, and everything else as v8 keeps it.
{
  const pos = Array.from({ length: 5000 }, (_, i) => 1_000_000 + i * 37)
  const shared = { n: 1 }
  const x: any = {
    pos, big: [2 ** 31 - 1, -(2 ** 31), ...Array(70).fill(0)], floats: Array.from({ length: 100 }, (_, i) => i / 3),
    negz: [-0, ...Array(70).fill(1)], mixed: [1, 'a', ...Array(70).fill(2)], short: [1, 2, 3],
    holes: Object.assign(new Array(80).fill(1), { 5: undefined }), m: new Map([['a', pos]]), s: new Set([1, 2]),
    t: new Float64Array([NaN, 1.5]), nan: NaN, u: undefined, a: shared, b: shared, e: new Error('kept'),
  }
  delete x.holes[7]
  x.self = x
  const y: any = roundTrip(x)
  assert.deepEqual(y.pos, pos)
  assert.deepEqual(y.big, x.big)
  assert.deepEqual(y.floats, x.floats)
  assert.ok(Object.is(y.negz[0], -0))
  assert.deepEqual(y.mixed, x.mixed)
  assert.ok(!(7 in y.holes) && y.holes.length === 80, 'a hole stays a hole')
  assert.deepEqual(y.m.get('a'), pos)
  assert.equal(y.m.get('a'), y.pos, 'one array named twice stays one array')
  assert.equal(y.a, y.b, 'a shared object stays shared')
  assert.equal(y.self, y, 'a cycle stays a cycle')
  assert.ok(Number.isNaN(y.nan) && 'u' in y && y.u === undefined)
  assert.ok(y.t instanceof Float64Array && Number.isNaN(y.t[0]))
  assert.equal(y.e.message, 'kept')
  const gz = (v: unknown) => gzipSync(serialize(v)).length
  assert.ok(gz(pack({ pos })) < gz({ pos }) / 2, 'positions are held small')
}

// Records held as columns come back key for key, value for value, in order; shared or uneven ones
// are kept as they are.
{
  const cn = Array.from({ length: 300 }, (_, i) => ({ chrom: String(1 + (i % 3)), pos: i * 10, called: i % 4 !== 0,
    log2R: i % 7 === 0 ? null : i % 11 === 0 ? NaN : i % 13 === 0 ? -0 : (i - 150) / 37, baf: i % 5 === 0 ? undefined : i / 300 }))
  const y: any = roundTrip({ m: new Map([['1', cn]]) })
  const back = y.m.get('1')
  assert.equal(back.length, cn.length)
  for (let i = 0; i < cn.length; i++) {
    assert.deepEqual(Object.keys(back[i]), Object.keys(cn[i]))
    for (const k of Object.keys(cn[i])) assert.ok(Object.is(back[i][k], (cn[i] as any)[k]), `${i} ${k}`)
  }
  const sh = { a: 1 }; const withShared = [sh, ...Array.from({ length: 80 }, (_, i) => ({ a: i }))]
  const z: any = roundTrip({ xs: withShared, again: sh })
  assert.equal(z.xs[0], z.again, 'an element reached twice stays one object')
  const uneven = Array.from({ length: 80 }, (_, i) => (i === 40 ? { a: 1, b: 2 } : { a: i }))
  assert.deepEqual(roundTrip(uneven), uneven)
  const nested = Array.from({ length: 80 }, (_, i) => ({ a: i, inner: { q: i } }))
  assert.deepEqual(roundTrip(nested), nested)
  const gz = (v: unknown) => gzipSync(serialize(v)).length
  assert.ok(gz(pack({ cn })) < gz({ cn }), 'columns are smaller')
}

// Keys, the open list, put and get.
const dir = mkdtempSync(join(tmpdir(), 'om-store-check-'))
try {
  const a = join(dir, 'a.probes'), b = join(dir, 'b.probes'), sealed = join(dir, 'sealed.probes')
  writeFileSync(a, 'one'); writeFileSync(b, 'two'); writeFileSync(sealed, 'never')
  const st = openStore(join(dir, 'store'), 'h1', new Set([a, b]))
  const k = st.key(['x'], [a, '-'])
  assert.equal(st.key(['x'], [a, '-']), k)
  assert.notEqual(st.key(['y'], [a, '-']), k, 'arguments enter the key')
  assert.notEqual(openStore(join(dir, 'store'), 'h2', new Set([a])).key(['x'], [a, '-']), k, 'the closure enters the key')
  assert.notEqual(st.key(['x'], [a, b]), k)
  assert.ok(!k.includes('probes'))
  assert.equal(st.allows([a, '-', sealed]), false)
  assert.throws(() => st.key(['x'], [sealed]), /refused/)
  writeFileSync(a, 'one, changed')
  assert.equal(openStore(join(dir, 'store'), 'h1', new Set([a])).key(['x'], [a, '-']) === k, false, 'the bytes enter the key')
  assert.equal(st.get(k), undefined)
  st.put(k, { res: { v: [1, 2] } })
  assert.deepEqual(st.get(k), { res: { v: [1, 2] } })
  assert.deepEqual(st.stats.hits, 1)

  // A held object answers what it holds and refuses the rest; a lazy field is read once.
  let n = 0
  const o: any = only({ sex: 'female', none: undefined }, { cn: () => ++n })
  assert.equal(o.sex, 'female'); assert.equal(o.none, undefined)
  assert.equal(o.cn, 1); assert.equal(o.cn, 1)
  assert.throws(() => o.heterozygosity, /does not hold/)
  assert.equal(o.then, undefined)

  // The closure's hash is trusted only for the tree it names.
  const tree = join(dir, 'tree'); mkdirSync(tree)
  writeFileSync(join(tree, 'x.ts'), 'export const x = 1')
  const sha = createHash('sha256').update('export const x = 1').digest('hex')
  writeFileSync(join(tree, 'GATE-KEYS.json'), JSON.stringify({ rows: { score: { hash: 'H', files: { 'x.ts': sha } } } }))
  assert.equal(treeKey(tree, 'rows'), 'H')
  writeFileSync(join(tree, 'x.ts'), 'export const x = 2')
  assert.throws(() => treeKey(tree, 'rows'), /does not match/)
} finally {
  rmSync(dir, { recursive: true, force: true })
}
console.log('store checks pass')

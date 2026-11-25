/**
 * ProbeMap behaves as a Map in every way the scorer uses one, and readerOf gives what Map.get gives,
 * for two chips' probes in one session.
 *
 * Run (from web/): node --experimental-strip-types src/probeIndex.check.ts
 */
import assert from 'node:assert/strict'
import { ProbeMap, probeNumber, readerOf } from './probeIndex.ts'

let seed = 3
const rnd = () => { seed = (seed * 1103515245 + 12345) % 2 ** 31; return seed / 2 ** 31 }
const ids = (p: string, n: number) => Array.from({ length: n }, (_, i) => `${p}${i}`)
// Two chips, partly the same probes, read in one session.
for (const chip of [ids('AX-', 3000), [...ids('AX-', 1500), ...ids('SNP_A-', 2000)]]) {
  const ref = new Map<string, unknown>(), pm = new ProbeMap<unknown>()
  for (let step = 0; step < 6000; step++) {
    const id = chip[Math.floor(rnd() * chip.length)]
    const op = rnd()
    if (op < 0.6) { const v = rnd() < 0.1 ? undefined : rnd() < 0.5 ? 'AB' : rnd(); ref.set(id, v); if (rnd() < 0.5) pm.set(id, v); else pm.setAt(probeNumber(id), v) }
    else if (op < 0.65) assert.equal(pm.delete(id), ref.delete(id))
    else assert.equal(pm.get(id), ref.get(id))
    assert.equal(pm.has(id), ref.has(id))
  }
  assert.equal(pm.size, ref.size)
  assert.deepEqual([...pm], [...ref])
  assert.deepEqual([...pm.keys()], [...ref.keys()])
  assert.deepEqual([...pm.values()], [...ref.values()])
  const seen: unknown[] = []; pm.forEach((v, k) => seen.push([k, v])); assert.deepEqual(seen, [...ref])
  assert.deepEqual([...new Map(pm)], [...ref])
  assert.equal(pm.get('never-seen'), undefined)
  // A reader gives what get gives, absent ids and stored undefined included, on each chip's own Map.
  const m = new Map<string, string | undefined>(chip.filter(() => rnd() < 0.5).map((id) => [id, rnd() < 0.1 ? undefined : 'AA']))
  const read = readerOf(m)
  for (let r = 0; r < 2; r++) for (const id of chip) assert.equal(read(probeNumber(id), id), m.get(id), id)
}
// Two Maps over the same ids keep separate values.
const a = new Map([['AX-1', 'AA']]), b = new Map([['AX-1', 'BB']])
assert.equal(readerOf(a)(probeNumber('AX-1'), 'AX-1'), 'AA')
assert.equal(readerOf(b)(probeNumber('AX-1'), 'AX-1'), 'BB')
console.log('probeIndex checks pass')

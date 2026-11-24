/**
 * sortNumbers gives what the comparator sort gives, element for element (Object.is), on ordinary
 * arrays and on the ones it hands back to the comparator.
 *
 * Run (from web/): node --experimental-strip-types src/sortNumbers.check.ts
 */
import assert from 'node:assert/strict'
import { sortNumbers } from './sortNumbers.ts'

let seed = 7
const rnd = () => { seed = (seed * 1103515245 + 12345) % 2 ** 31; return seed / 2 ** 31 }
const same = (xs: number[]) => {
  const a = sortNumbers([...xs]), b = [...xs].sort((p, q) => p - q)
  assert.equal(a.length, b.length)
  for (let i = 0; i < a.length; i++) assert.ok(Object.is(a[i], b[i]), `at ${i}: ${a[i]} against ${b[i]}`)
}
for (let r = 0; r < 200; r++) {
  const n = Math.floor(rnd() * 3000)
  const xs = Array.from({ length: n }, () => (rnd() < 0.1 ? 0 : (rnd() - 0.5) * 10 ** Math.floor(rnd() * 8 - 4)))
  if (rnd() < 0.3) for (let i = 0; i < 20 && n; i++) xs[Math.floor(rnd() * n)] = xs[0]
  same(xs)
}
same([])
same([1])
same([3, -0, 0, -0, 2])
same([NaN, 1, NaN, 0.5, -2])
same([Infinity, -Infinity, 1e308, -1e-308, 5e-324])
same([2, 1, 0.1 + 0.2, 0.3])
const arr = [3, 1, 2]
assert.equal(sortNumbers(arr), arr, 'sorts in place')
console.log('sortNumbers checks pass')

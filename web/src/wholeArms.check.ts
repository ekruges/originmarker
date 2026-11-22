/**
 * A whole call is read arm by arm: it stands where both usable arms carry it and on a chromosome with
 * one usable arm; it is withdrawn where one arm does not, naming the arm that does; a raised arm does
 * not carry a loss on its call fraction alone; and the arm a withdrawn call rests on prints as an
 * arm-level segment whose edges are not placed.
 */
import assert from 'node:assert/strict'
import { armSegment, P_ARM_END_GRCH37, readArms } from './wholeArms.ts'
import type { CopyNull } from './partLost.ts'
import { segmentCoords } from './segments.ts'
import { defectsFrom } from './defects.ts'

let seed = 11
const rand = (): number => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 }
const gauss = (): number => Math.sqrt(-2 * Math.log(rand() || 1e-12)) * Math.cos(2 * Math.PI * rand())

/** chr17 at one marker per 10 kb over 81 Mb, the p/q boundary at 24 Mb. `p` and `q` are each arm's
 *  true log R; `call` each arm's call rate. */
function chr17(p: number, q: number, call: [number, number] = [0.9, 0.9], chrom = '17', end = 81e6) {
  const split = P_ARM_END_GRCH37[chrom]
  return Array.from({ length: Math.floor(end / 10_000) }, (_, i) => {
    const pos = i * 10_000 + 1
    const onP = pos < split
    return { chrom, pos, called: rand() < (onP ? call[0] : call[1]), log2R: (onP ? p : q) + 0.2 * gauss() }
  })
}
const n: CopyNull = { centre: 0, scale: 0.05, units: 20, floor: 0.40, callMedian: 0.9 }

// A whole loss on both arms stands.
assert.equal(readArms(chr17(-0.6, -0.6), 'loss', -0.6, n).whole, true)
// A whole gain on both arms stands.
assert.equal(readArms(chr17(0.4, 0.4), 'gain', 0.4, n).whole, true)

// A q-arm loss moving the median: withdrawn, the q arm named.
{
  const r = readArms(chr17(0, -0.8), 'loss', -0.6, n)
  assert.equal(r.whole, false)
  assert.equal(r.arm, 'q')
  assert.ok(Math.abs(r.p!.shift) < 0.05 && r.q!.shift < -0.7)
}
// A q-arm gain: withdrawn, the q arm named. Half the gain floor (0.15) is the bar.
{
  const r = readArms(chr17(0.05, 0.5), 'gain', 0.4, n)
  assert.equal(r.whole, false)
  assert.equal(r.arm, 'q')
  assert.equal(readArms(chr17(0.2, 0.5), 'gain', 0.4, n).whole, true, 'p at +0.20 is over +0.15')
}
// A loss arm carried by its call fraction alone, at the cell's level: it carries the loss.
assert.equal(readArms(chr17(0, -0.8, [0.4, 0.9]), 'loss', -0.6, n).whole, true)
// One arm gone, the other raised and calling under CALL_COLLAPSE (a cell with a low call rate): the
// raised arm does not carry the loss, so the call is withdrawn and the lost arm named.
{
  const r = readArms(chr17(-1.8, 0.78, [0.25, 0.5]), 'loss', 0.38, n)
  assert.ok(r.q!.callFraction < 0.6)
  assert.equal(r.whole, false)
  assert.equal(r.arm, 'p')
}
// Both arms raised, a loss printed through the call fraction: neither carries it, no arm named.
{
  const r = readArms(chr17(0.37, 0.43, [0.5, 0.5]), 'loss', 0.4, n)
  assert.equal(r.whole, false)
  assert.equal(r.arm, undefined)
}
// An acrocentric: under 200 markers with a log R before the boundary, so the whole call stands.
{
  const ms = chr17(0, -0.8, [0.9, 0.9], '13', 115e6).filter((m) => m.pos >= 19.0e6 || m.pos % 1e6 < 10_000)
  const r = readArms(ms, 'loss', -0.8, n)
  assert.equal(r.p, null)
  assert.equal(r.whole, true)
}

// The arm call: the q arm from the boundary to the last marker, edges stated as not placed.
{
  const ms = chr17(0, -0.8)
  const sg = armSegment(ms, 'q', 'copy-loss')
  assert.equal(sg.startBp, 24_000_000)
  assert.equal(sg.endBp, ms[ms.length - 1].pos)
  assert.equal(sg.arm, 'q')
  assert.match(segmentCoords(sg).interval, /^not placed: an arm-level call on 17q/)
  const p = armSegment(ms, 'p', 'copy-gain')
  assert.equal(p.startBp, 1)
  assert.equal(p.endBp, 24_000_000)
  const [d] = defectsFrom([sg], [], [])
  assert.equal(d.kind, 'copy-loss')
  assert.equal(d.arm, 'q')
  assert.equal(d.interval, segmentCoords(sg).interval)
  assert.match(d.why, /^Arm-level call\. .* only the q arm carries the change/)
}

console.log('wholeArms: ok')

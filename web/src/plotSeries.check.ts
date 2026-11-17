// Series selection and box geometry for the call plots.
import assert from 'node:assert'
import {
  boxGeometry, boxLrr, CAP_INFORMATIVE, CAP_OTHER, callSeries, cellLrr, median, nearestIndex, segmentBox, tickStepMb, whoseLetter,
} from './plotSeries.ts'

// ---------------------------------------------------------------- whose letter
assert.equal(whoseLetter('BB', 'AA', 'BB'), 2, "father AA has no B: mother's letter")
assert.equal(whoseLetter('AA', 'AA', 'BB'), 1, "mother BB has no A: father's letter")
assert.equal(whoseLetter('AB', 'AA', 'BB'), 0, 'both letters')
assert.equal(whoseLetter('BB', 'AA'), 2, 'one parent: the other letter is the other parent')
assert.equal(whoseLetter('AA', 'AA'), 0, "one parent: the parent's own letter is not readable")
assert.equal(whoseLetter('BB', 'AA', 'AA'), 0, 'not Mendelian: neither')
assert.equal(whoseLetter('BB'), 0, 'no parent')
assert.equal(whoseLetter('NC', 'AA', 'BB'), 0, 'no call')

// ---------------------------------------------------------------- selection
let seed = 7
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x80000000 }
const build = (perChrom: Record<string, number>) => {
  const markerPos = new Map<string, { chrom: string; pos: number }>()
  const myBaf = new Map<string, number>()
  const myGt = new Map<string, string>()
  const father = new Map<string, string>()
  const cnByChrom = new Map<string, { pos: number; log2R: number | null }[]>()
  let k = 0
  for (const [chrom, n] of Object.entries(perChrom)) {
    const cn: { pos: number; log2R: number | null }[] = []
    // File order is not position order on these arrays.
    const order = [...Array(n).keys()].sort(() => rnd() - 0.5)
    for (const i of order) {
      const id = `p${k += 1}`
      const pos = 1000 + i * 100
      markerPos.set(id, { chrom, pos })
      myBaf.set(id, i / n)
      myGt.set(id, i % 3 === 0 ? 'AB' : 'BB')
      father.set(id, i % 20 === 0 ? 'AA' : 'AB')
      cn.push({ pos, log2R: -i / n })
    }
    cnByChrom.set(chrom, cn)
  }
  return { src: { markerPos, myBaf, myGt, cnByChrom }, father }
}

{
  const { src, father } = build({ 1: 20000, 21: 300 })
  const caps = { informative: 200, other: 500 }
  const s = callSeries(src, father, null, caps)
  const c1 = s.get('1')!
  assert.equal(c1.total, 20000)
  assert.equal(c1.pos.length, 700, 'a long chromosome keeps exactly the budget')
  for (let i = 1; i < c1.pos.length; i += 1) assert.ok(c1.pos[i] > c1.pos[i - 1], 'sorted by position')
  assert.equal(c1.endBp, 1000 + 19999 * 100, 'end is the last marker, kept or not')
  // Informative markers: father AA at every 20th, 1,000 in all; 200 kept, not 700 * 1/20.
  let inf = 0
  for (let i = 0; i < c1.pos.length; i += 1) if (((c1.pos[i] - 1000) / 100) % 20 === 0) inf += 1
  assert.equal(inf, 200, 'informative markers keep their own budget')
  // Every kept point carries its own marker's values, joined by position.
  for (let i = 0; i < c1.pos.length; i += 1) {
    const j = (c1.pos[i] - 1000) / 100
    assert.ok(Math.abs(c1.baf[i] - j / 20000) < 1e-6)
    assert.ok(Math.abs(c1.lrr[i] + j / 20000) < 1e-6)
    assert.equal(c1.who[i], j % 20 === 0 && j % 3 !== 0 ? 2 : 0)
  }
  // The kept positions spread over the chromosome rather than bunching at the start of the file.
  const half = [...c1.pos].filter((p) => p < 1000 + 10000 * 100).length
  assert.ok(half > 250 && half < 450, `half the chromosome holds ${half} of 700`)
  const c21 = s.get('21')!
  assert.equal(c21.pos.length, 300, 'a short chromosome keeps every marker')
  // Same input, same points.
  const again = callSeries(src, father, null, caps).get('1')!
  assert.deepEqual([...again.pos], [...c1.pos], 'deterministic')
}

// No parent: nothing is attributed, nothing is informative.
{
  const { src } = build({ 2: 1000 })
  const s = callSeries(src, null, null, { informative: 10, other: 50 }).get('2')!
  assert.equal(s.pos.length, 50)
  assert.ok([...s.who].every((w) => w === 0))
}

// ---------------------------------------------------------------- nearest
{
  const pos = new Uint32Array([10, 20, 40])
  assert.equal(nearestIndex(pos, 0), 0)
  assert.equal(nearestIndex(pos, 14), 0)
  assert.equal(nearestIndex(pos, 16), 1)
  assert.equal(nearestIndex(pos, 31), 2)
  assert.equal(nearestIndex(pos, 99), 2)
  assert.equal(nearestIndex(new Uint32Array(0), 5), -1)
}

// ---------------------------------------------------------------- box geometry
{
  // Refined edges: the box is the point estimates, the intervals are drawn either side.
  const b = segmentBox({ startBp: 9e6, endBp: 31e6, refined: {
    startBp: 10e6, startLoBp: 9.5e6, startHiBp: 10.5e6, endBp: 30e6, endLoBp: 29e6, endHiBp: 31e6,
  } }, 'LOSS', [])
  const g = boxGeometry(b, 100e6, 50, 1000)
  assert.equal(g.x0, 150)
  assert.equal(g.x1, 350)
  assert.deepEqual(g.left, [145, 155])
  assert.deepEqual(g.right, [340, 360])
  // Window edges only: no intervals.
  const w = boxGeometry(segmentBox({ startBp: 9e6, endBp: 31e6 }, 'LOSS', []), 100e6, 0, 100)
  assert.equal(w.x0, 9)
  assert.equal(w.x1, 31)
  assert.equal(w.left, undefined)
  // Whole chromosome: the whole plot.
  const whole = boxGeometry({ startBp: 0, endBp: 100e6, label: '', lines: [] }, 100e6, 50, 1000)
  assert.deepEqual([whole.x0, whole.x1], [50, 1050])
  // A run shorter than a pixel is still 2 px wide, centred on it.
  const tiny = boxGeometry({ startBp: 50e6, endBp: 50.001e6, label: '', lines: [] }, 100e6, 0, 1000)
  assert.ok(Math.abs(tiny.x1 - tiny.x0 - 2) < 1e-9 && Math.abs((tiny.x0 + tiny.x1) / 2 - 500.005) < 1e-6)
  // Past the end is clamped to the plot.
  const over = boxGeometry({ startBp: 90e6, endBp: 120e6, label: '', lines: [] }, 100e6, 0, 100)
  assert.equal(over.x1, 100)
}

assert.equal(median([3, NaN, 1, 2]), 2)
assert.equal(median([4, 1, 2, 3]), 2.5)
assert.ok(Number.isNaN(median([NaN])))
{
  const mk = (pos: number[], lrr: number[]) => ({
    pos: Uint32Array.from(pos), lrr: Float32Array.from(lrr), baf: new Float32Array(pos.length),
    gt: new Uint8Array(pos.length), who: new Uint8Array(pos.length), total: pos.length, endBp: pos.at(-1)!,
  })
  const a = mk([10, 20, 30, 40], [0, -1, -1, 0])
  assert.equal(boxLrr(a, 15, 35), -1, 'inside the box only')
  assert.equal(cellLrr(new Map([['1', a], ['2', mk([5, 6], [0, NaN])]])), 0)
}

assert.equal(tickStepMb(249e6), 50)
assert.equal(tickStepMb(48e6), 10)
assert.equal(tickStepMb(3e6), 1)

// ---------------------------------------------------------------- size and growth
// A real array is about 825,000 markers. Kept points are bounded per chromosome whatever the size,
// and the pass grows linearly: four times the markers costs about four times the time, where a
// pass that rescanned per marker would cost sixteen. Best of three draws, so one slow draw on a
// busy machine does not read as growth.
{
  const shares = Array.from({ length: 22 }, (_, i) => [String(i + 1), 0.08 - i * 0.003] as const)
  const sum = shares.reduce((a, [, x]) => a + x, 0)
  const sized = (total: number) => build(Object.fromEntries(shares.map(([c, x]) => [c, Math.round((x / sum) * total)])))
  const small = sized(206_000)
  const large = sized(825_000)
  const best = (a: ReturnType<typeof build>) => Math.min(...[0, 1, 2].map(() => {
    const t0 = performance.now()
    callSeries(a.src, a.father, null)
    return performance.now() - t0
  }))
  const out = callSeries(large.src, large.father, null)
  let kept = 0
  for (const s of out.values()) kept += s.pos.length
  assert.ok(kept <= 22 * (CAP_INFORMATIVE + CAP_OTHER), `${kept} points kept at 825,000 markers`)
  const ts = best(small)
  const tl = best(large)
  if (process.env.OM_TIMING) console.log(`callSeries: ${ts.toFixed(0)} ms at 206k, ${tl.toFixed(0)} ms at 825k`)
  assert.ok(tl / ts < 8, `four times the markers took ${(tl / ts).toFixed(1)}x the time`)
}

console.log('plotSeries: ok')

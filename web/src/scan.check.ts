// THE CRASH CHECK.
//
// A real run parses about 825,000 markers. The scan that shipped in 5.10.0 rescanned every marker
// of a chromosome inside its own window loop, so its cost was markers TIMES windows, and the tab
// locked up on every run. It passed every test in this repo, because the cost only appears at a
// size no fixture had and the code only existed inside a React component.
//
// So this check does two things a unit test does not. It runs at the real size, and it runs at two
// sizes and compares. A wall-clock ceiling alone is a machine-dependent flake; the RATIO between
// two sizes is not, and quadratic growth is exactly what the ratio catches: double the markers and
// a linear pass takes about twice as long, while a quadratic one takes about four times.
import assert from 'node:assert'
import { breathe, buildScanIndex, copyNeutralWindows, gatherInterval } from './scan.ts'
import type { CnMarker, SelfMarker } from './scan.ts'

// ---------------------------------------------------------------- synthetic array
// Marker spacing and per-chromosome counts in the proportions a genome-wide array has, so the
// window count and the chromosome count are both realistic.
const CHROM_SHARE = [8.1, 7.9, 6.5, 6.2, 5.9, 5.6, 5.2, 4.7, 4.5, 4.4, 4.4, 4.3,
  3.1, 2.9, 2.7, 2.6, 2.4, 2.3, 1.9, 1.9, 1.2, 1.3]

let seed = 20260818
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }

const buildArray = (total: number) => {
  const cnByChrom = new Map<string, CnMarker[]>()
  const selfMarkers: SelfMarker[] = []
  const markerPos = new Map<string, { chrom: string; pos: number }>()
  const parentGt = new Map<string, string>()
  const myBaf = new Map<string, number>()
  const myGt = new Map<string, string>()
  const share = CHROM_SHARE.reduce((a, x) => a + x, 0)
  let id = 0
  for (let ci = 0; ci < CHROM_SHARE.length; ci += 1) {
    const chrom = String(ci + 1)
    const n = Math.round((CHROM_SHARE[ci] / share) * total)
    const ms: CnMarker[] = []
    for (let i = 0; i < n; i += 1) {
      const probe = `p${id += 1}`
      const pos = i * 3000 + Math.floor(rnd() * 2000)
      const called = rnd() > 0.03
      const gt = rnd() < 0.33 ? 'AB' : rnd() < 0.5 ? 'AA' : 'BB'
      ms.push({ chrom, pos, called, log2R: rnd() < 0.02 ? null : (rnd() - 0.5) * 0.4 })
      markerPos.set(probe, { chrom, pos })
      parentGt.set(probe, rnd() < 0.33 ? 'AB' : rnd() < 0.5 ? 'AA' : 'BB')
      myBaf.set(probe, rnd())
      myGt.set(probe, gt)
      if (called) selfMarkers.push({ chrom, pos, het: gt === 'AB' })
    }
    cnByChrom.set(chrom, ms)
  }
  return { cnByChrom, selfMarkers, src: { markerPos, parentGt, myBaf, myGt, cnByChrom } }
}

const time = async (f: () => Promise<unknown> | unknown): Promise<number> => {
  const t0 = performance.now()
  await f()
  return performance.now() - t0
}

// ---------------------------------------------------------------- correctness first
// A fast wrong answer is not the goal. The prefix sum and the binary search have to agree with the
// obvious version they replaced, on an array small enough to compute the obvious version.
{
  const { cnByChrom, selfMarkers } = buildArray(40_000)
  const SEG = 600
  const { windows, chromEnd } = await copyNeutralWindows(cnByChrom, selfMarkers, 0.01, SEG)
  assert.deepStrictEqual([...chromEnd.keys()].sort(), [...cnByChrom.keys()].sort(),
    'every chromosome gets an end')

  const selfBy = new Map<string, SelfMarker[]>()
  for (const m of selfMarkers) {
    const cur = selfBy.get(m.chrom); if (cur) cur.push(m); else selfBy.set(m.chrom, [m])
  }
  let checked = 0
  for (const w of windows) {
    if (w.wholeChromosome) continue
    // The naive count: every called marker in the interval, no prefix sum, no binary search.
    const naive = (selfBy.get(w.chrom) ?? [])
      .filter((m) => m.pos >= w.startBp && m.pos <= w.endBp && m.het).length
    assert.strictEqual(w.het, naive, `window het count on chr${w.chrom} ${w.startBp}`)
    checked += 1
  }
  assert.ok(checked > 0, 'the scan produced sliding windows to check')
  console.log(`  windows agree with the naive count on ${checked} sliding windows`)
}

// gatherInterval used to walk the marker map three times per finding and build a fresh row for
// every marker it saw. It now slices a prebuilt index and SHARES the rows, so the rows arrive
// grouped by chromosome and sorted by position rather than in the order the file was parsed.
//
// That reordering is only safe because of what happens downstream: callDosageOrigin reduces both
// row sets to a centroid, a mean, which does not depend on order. So this checks the multiset of
// rows, and then checks the centroid itself, which is the number the origin call is actually made
// from. If a consumer ever starts reading these in order, this is the check that has to change
// with it.
{
  const { src } = buildArray(40_000)
  // Half of chr3, taken from where its markers actually are rather than a fixed base pair, so the
  // interval stays non-empty whatever size the synthetic array is built at.
  const on3 = [...src.markerPos.values()].filter((p) => p.chrom === '3').map((p) => p.pos)
  const mid = on3.sort((a, b) => a - b)[on3.length >> 1]
  const inside = (c: string, pos: number) => c === '3' && pos > mid
  const g = gatherInterval(buildScanIndex(src), { chrom: '3', startBp: mid + 1 })

  const naiveRegion: [string, number | null][] = []
  const naiveBackground: [string, number | null][] = []
  const naiveUnt: [string, string, number | null][] = []
  const naiveIn: number[] = []
  const naiveOut: number[] = []
  for (const [probe, p] of src.markerPos) {
    const pg = src.parentGt.get(probe)
    if (!pg) continue
    const b = src.myBaf.get(probe) ?? null
    ;(inside(p.chrom, p.pos) ? naiveRegion : naiveBackground).push([pg, b])
    if (inside(p.chrom, p.pos) && pg === 'AB') {
      naiveUnt.push([pg, src.myGt.get(probe) ?? 'NC', b])
    }
  }
  for (const [ch, ms] of src.cnByChrom) {
    for (const m of ms) {
      if (m.log2R === null || !Number.isFinite(m.log2R)) continue
      ;(inside(ch, m.pos) ? naiveIn : naiveOut).push(m.log2R)
    }
  }

  const bag = (rows: readonly (readonly unknown[])[]) => rows.map((r) => JSON.stringify(r)).sort()
  assert.deepStrictEqual(bag(g.region), bag(naiveRegion), 'same region rows')
  assert.deepStrictEqual(bag(g.background), bag(naiveBackground), 'same background rows')
  assert.deepStrictEqual(bag(g.untRows), bag(naiveUnt), 'same untransmitted rows')

  // The numbers the call is made from, not just the rows behind them.
  const mean = (xs: number[]) => xs.reduce((a, x) => a + x, 0) / xs.length
  const bafMean = (rows: readonly (readonly [string, number | null])[]) =>
    mean(rows.map((r) => r[1]).filter((x): x is number => x !== null))
  assert.ok(Math.abs(bafMean(g.region) - bafMean(naiveRegion)) < 1e-12, 'region centroid unchanged')
  assert.ok(Math.abs(bafMean(g.background) - bafMean(naiveBackground)) < 1e-12,
    'background centroid unchanged')
  assert.ok(Math.abs(mean(g.inL) - mean(naiveIn)) < 1e-12, 'in-region intensity mean unchanged')
  assert.ok(Math.abs(mean(g.outL) - mean(naiveOut)) < 1e-12, 'background intensity mean unchanged')
  assert.ok(g.region.length > 0 && g.inL.length > 0, 'the interval was not empty')
  console.log(`  indexed gather matches the per-finding walk on ${g.region.length} region and`
    + ` ${g.background.length} background rows, centroids identical`)
}

// A whole chromosome and the same chromosome given as an unbounded interval are the same request,
// and the fast path that skips the binary search has to agree with the one that does not.
{
  const { src } = buildArray(20_000)
  const idx = buildScanIndex(src)
  const whole = gatherInterval(idx, { chrom: '5' })
  const bounded = gatherInterval(idx, { chrom: '5', startBp: 0, endBp: Number.MAX_SAFE_INTEGER })
  assert.strictEqual(whole.region.length, bounded.region.length, 'whole chromosome, both ways')
  assert.strictEqual(whole.background.length, bounded.background.length, 'background, both ways')
  assert.strictEqual(whole.inL.length, bounded.inL.length, 'intensity, both ways')
  console.log(`  whole-chromosome fast path agrees with the bounded path`
    + ` (${whole.region.length} rows)`)
}

// MARKERS OUT OF ORDER. The index skips building a permutation when the positions already
// ascend, which is the normal case for an array export and most of what the build used to cost.
// The synthetic array above is built in order, so it only ever exercises the fast path. This
// shuffles the input and requires the same answer, because a file that arrives out of order and is
// silently treated as sorted would return the wrong markers for every interval on that chromosome.
{
  const { src } = buildArray(20_000)
  const shuffled = new Map<string, { chrom: string; pos: number }>()
  const rows = [...src.markerPos.entries()]
  // Deterministic reversal per chromosome, which is as out-of-order as input gets.
  for (const [probe, p] of rows.reverse()) shuffled.set(probe, p)
  const shuffledSrc = { ...src, markerPos: shuffled }

  const inOrder = gatherInterval(buildScanIndex(src), { chrom: '4' })
  const outOfOrder = gatherInterval(buildScanIndex(shuffledSrc), { chrom: '4' })
  const bag = (rows: readonly (readonly unknown[])[]) => rows.map((r) => JSON.stringify(r)).sort()
  assert.deepStrictEqual(bag(outOfOrder.region), bag(inOrder.region),
    'a shuffled input must yield the same region rows')
  assert.strictEqual(outOfOrder.background.length, inOrder.background.length,
    'and the same background size')

  // And a bounded interval, where a wrong permutation shows up as the wrong slice rather than the
  // wrong order.
  const on4 = [...src.markerPos.values()].filter((p) => p.chrom === '4').map((p) => p.pos)
    .sort((a, b) => a - b)
  const lo = on4[Math.floor(on4.length * 0.25)]
  const hi = on4[Math.floor(on4.length * 0.75)]
  const naive = rows.filter(([, p]) => p.chrom === '4' && p.pos >= lo && p.pos <= hi).length
  const sliced = gatherInterval(buildScanIndex(shuffledSrc),
    { chrom: '4', startBp: lo, endBp: hi })
  assert.strictEqual(sliced.region.length, naive,
    `a bounded slice of shuffled input must hold ${naive} markers, got ${sliced.region.length}`)
  console.log(`  shuffled input yields the same ${inOrder.region.length} rows, and a bounded slice`
    + ` of ${naive}`)
}

// ---------------------------------------------------------------- the growth curve
// COST PER MARKER, not total time. A linear pass costs the same per marker whatever the size of
// the array; a quadratic one costs four times as much per marker when the array is four times as
// big. That framing is what makes this robust on a shared runner: absolute times move with the
// machine, the SHAPE of the curve does not.
//
// Each size is built, warmed and timed on its own and then dropped. Holding two large arrays alive
// at once measures the garbage collector instead of the algorithm, which read as 4.3x growth on a
// scan that is provably flat.
const SMALL = 100_000
const LARGE = 400_000
// Four times the markers. Quadratic lands at 4x the per-marker cost. 2x leaves room for the sort
// terms and for a noisy runner, and still fails long before anything user-visible.
const GROWTH_LIMIT = 2.0

const usPerMarker = async (
  n: number,
  work: (a: ReturnType<typeof buildArray>) => Promise<unknown> | unknown,
): Promise<number> => {
  const arr = buildArray(n)
  await work(arr) // warm, so the first call does not pay for compilation
  const ms = await time(() => work(arr))
  return (ms * 1000) / n
}

const scanWork = (a: ReturnType<typeof buildArray>) =>
  copyNeutralWindows(a.cnByChrom, a.selfMarkers, 0.01, 600)

{
  const small = await usPerMarker(SMALL, scanWork)
  const large = await usPerMarker(LARGE, scanWork)
  const growth = large / small
  console.log(`  window scan: ${small.toFixed(3)}us/marker at ${SMALL / 1000}k,`
    + ` ${large.toFixed(3)}us/marker at ${LARGE / 1000}k (${growth.toFixed(2)}x,`
    + ' quadratic would be ~4x)')
  assert.ok(growth < GROWTH_LIMIT,
    `the window scan costs ${growth.toFixed(2)}x more per marker on an array four times as big,`
    + ` limit ${GROWTH_LIMIT}x. Something in copyNeutralWindows is rescanning, which is what`
    + ' locked the tab in 5.10.0.')
}

{
  // The index is built once per sample, so the per-finding cost is the walk and not the build.
  const gatherWork = (a: ReturnType<typeof buildArray>) => {
    const idx = buildScanIndex(a.src)
    for (let i = 0; i < 5; i += 1) gatherInterval(idx, { chrom: String(i + 1) })
  }
  const small = await usPerMarker(SMALL, gatherWork)
  const large = await usPerMarker(LARGE, gatherWork)
  const growth = large / small
  console.log(`  interval gather: ${small.toFixed(3)}us/marker at ${SMALL / 1000}k,`
    + ` ${large.toFixed(3)}us/marker at ${LARGE / 1000}k (${growth.toFixed(2)}x)`)
  assert.ok(growth < GROWTH_LIMIT,
    `gatherInterval costs ${growth.toFixed(2)}x more per marker on an array four times as big,`
    + ` limit ${GROWTH_LIMIT}x`)
}

// ---------------------------------------------------------------- the real size
// A whole run at the size a real array actually is, with the number of findings a disturbed sample
// actually produces. The ceiling is deliberately loose: it is here to catch a return to something
// that locks the tab, not to police a few milliseconds. A 2 second budget on this machine leaves
// room for a slow shared CI runner and still fails long before a user would.
const REAL_MARKERS = 825_000
// Deliberately loose. This is here to catch a return to something that locks the tab, not to
// police a few milliseconds, and a shared runner is slower than a laptop. The version that made
// the page unusable measured 9 seconds on this same shape, so anything like it still fails.
const RUN_BUDGET_MS = 4000
const FINDINGS = 25
{
  const { cnByChrom, selfMarkers, src } = buildArray(REAL_MARKERS)
  const counted = [...cnByChrom.values()].reduce((a, ms) => a + ms.length, 0)
  const tScan = await time(() => copyNeutralWindows(cnByChrom, selfMarkers, 0.01, 600))
  const tIndex = await time(() => buildScanIndex(src))
  const idx = buildScanIndex(src)
  // One gather per finding, which is what a run with this many events does.
  const tGather = await time(() => {
    for (let i = 0; i < FINDINGS; i += 1) {
      gatherInterval(idx, { chrom: String((i % 22) + 1), startBp: i * 100_000 })
    }
  })
  const total = tScan + tIndex + tGather
  console.log(`  a run at ${counted.toLocaleString()} markers with ${FINDINGS} findings:`
    + ` ${total.toFixed(0)}ms total`)
  console.log(`    window scan ${tScan.toFixed(0)}ms, index ${tIndex.toFixed(0)}ms,`
    + ` ${FINDINGS} gathers ${tGather.toFixed(0)}ms`)
  assert.ok(total < RUN_BUDGET_MS,
    `a run at ${counted} markers took ${total.toFixed(0)}ms, budget ${RUN_BUDGET_MS}ms`
    + ` (scan ${tScan.toFixed(0)}, index ${tIndex.toFixed(0)}, gathers ${tGather.toFixed(0)}).`
    + ' This is the check that stands between a slow path and a locked tab.')
}

// ---------------------------------------------------------------- the yield
// Fast is not the same as smooth. Whatever the arithmetic costs, the scan has to hand the page back
// between chromosomes, or the whole genome is one task and the tab is frozen for its duration.
{
  const { cnByChrom, selfMarkers } = buildArray(50_000)
  const yielded: string[] = []
  await copyNeutralWindows(cnByChrom, selfMarkers, 0.01, 600, (c) => { yielded.push(c) })
  assert.strictEqual(yielded.length, cnByChrom.size,
    'the scan yields once per chromosome so the page can paint')
  console.log(`  yields to the page ${yielded.length} times, once per chromosome`)
}

// THE BACKGROUND TAB. The first version of the yield used requestAnimationFrame, which does not
// fire at all while a tab is hidden: a run the user tabbed away from waited forever, which is a
// worse failure than the slow scan the yield was added to fix. Node has no rAF and no document, so
// this environment is the hidden tab, and a yield that resolves here resolves there.
{
  const t0 = performance.now()
  for (let i = 0; i < 200; i += 1) await breathe()
  const ms = performance.now() - t0
  console.log(`  200 yields with no requestAnimationFrame and no document: ${ms.toFixed(0)}ms`)
  assert.ok(ms < 1000, `200 yields took ${ms.toFixed(0)}ms, which means they are being throttled`)
}

console.log('scan: growth, budget, yielding and background-tab progress all hold')

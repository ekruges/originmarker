// Tracks, bins, copies from each parent, the level bands, the guard on a level of none, and the
// panels the figure draws: one for every change, each reachable, laid out so nothing collides.
import assert from 'node:assert'
import { readFileSync } from 'node:fs'
import {
  badgeX, BASE_BP, binBase, CELL_PER, copyBars, COPY_NONE_BELOW, COPY_TWO_ABOVE, copyWord, drawable,
  guardedBars, MIN_BIN_MARKERS, mergeBins, normalRange, perForStretch, plotData, sheet, tickStepMb,
  totalBars, valueOf, viewOf, wrap, zoomOrder,
  type Bin, type PlotData, type PlotSource, type Track,
} from './plotSeries.ts'

// ---------------------------------------------------------------- binning
{
  // A base bin is a stretch of the genome: markers inside one 2.5 Mb stretch make one bin.
  const pos = Float64Array.from({ length: 260 }, (_, i) => (259 - i) * 40_000)
  const x = Float64Array.from({ length: 260 }, (_, i) => ((259 - i) % 2 ? 1 : 0))
  const b = binBase(pos, x)
  assert.ok(b.every((y) => y.endBp - y.startBp >= BASE_BP), 'a bin covers at least one stretch')
  assert.equal(b[0].startBp, 0)
  assert.equal(b[0].endBp, BASE_BP)
  assert.equal(b.reduce((s, y) => s + y.n, 0), 260, 'every marker lands in one bin')
  assert.ok(b.every((y) => y.n >= MIN_BIN_MARKERS), 'no bin stands on almost nothing')
  const m = mergeBins(b, CELL_PER)
  assert.ok(m.length <= b.length, 'bars merge whole base bins')
  assert.equal(m.reduce((s, y) => s + y.n, 0), 260, 'merging keeps every marker')
  // A thin stretch widens its bin instead of making a bar out of two markers.
  const sparse = binBase(Float64Array.from([0, 1e6, 2e6, 9e6, 9.1e6]), Float64Array.from([1, 1, 1, 1, 1]))
  assert.ok(sparse.every((y) => y.n >= MIN_BIN_MARKERS || sparse.length === 1), 'a thin stretch widens')
  assert.equal(sparse.reduce((s, y) => s + y.n, 0), 5)
  const lv = binBase(Float64Array.from([1, 2, 3]), Float64Array.from([5, -1, 2]), true)
  assert.equal(lv[0].k, 2, 'a level bin is a median')
}

// ---------------------------------------------------------------- a cell to read
/**
 * One synthetic cell, built marker by marker so every case the figure must draw is in it.
 *
 * `plan` gives each chromosome what it carries: how many copies came from each parent, and whether
 * the array still reads its markers. The alleles follow from the copies, the intensity follows from
 * their sum, and the call rate follows from `read`.
 */
type Plan = { fat: number; mot: number; read?: number; from?: number }
const MARKERS = 2400
const PITCH = 50_000

function cell(plan: Record<string, Plan>): {
  src: PlotSource; F: Map<string, string>; M: Map<string, string>
} {
  const markerPos = new Map<string, { chrom: string; pos: number }>()
  const myGt = new Map<string, string>()
  const myBaf = new Map<string, number>()
  const cnByChrom = new Map<string, { chrom: string; pos: number; called: boolean; log2R: number | null }[]>()
  const F = new Map<string, string>(); const M = new Map<string, string>()
  let seed = 7
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x80000000 }
  for (const [chrom, p] of Object.entries(plan)) {
    const cn: { chrom: string; pos: number; called: boolean; log2R: number | null }[] = []
    for (let i = 0; i < MARKERS; i += 1) {
      const id = `${chrom}:${i}`
      const pos = i * PITCH
      const after = i >= (p.from ?? 0)
      const fat = after ? p.fat : 1
      const mot = after ? p.mot : 1
      const read = after ? p.read ?? 1 : 1
      markerPos.set(id, { chrom, pos })
      // Half the markers are opposite homozygotes, which is where whose DNA this is can be read.
      const f = i % 2 === 0 ? 'AA' : 'BB'
      const m = i % 2 === 0 ? 'BB' : 'AA'
      F.set(id, f); M.set(id, m)
      const total = fat + mot
      const called = rnd() < read
      // The father's allele fraction is his copies over the copies in all; with nothing left the
      // array reads background, which lands near a half and still returns calls.
      const bFromFather = f === 'BB' ? fat : mot
      const baf = total ? bFromFather / total : 0.5
      myBaf.set(id, Math.max(0, Math.min(1, baf + (rnd() - 0.5) * 0.04)))
      if (called) myGt.set(id, baf > 0.8 ? 'BB' : baf < 0.2 ? 'AA' : 'AB')
      cn.push({ chrom, pos, called, log2R: Math.log2(Math.max(0.05, total) / 2) + (rnd() - 0.5) * 0.04 })
    }
    cnByChrom.set(chrom, cn)
  }
  return { src: { markerPos, myGt, myBaf, cnByChrom }, F, M }
}

const PLAN: Record<string, Plan> = {
  1: { fat: 1, mot: 1 },
  2: { fat: 1, mot: 1 },
  3: { fat: 1, mot: 1 },
  4: { fat: 0, mot: 1 },                      // a whole chromosome, the father's copy gone
  5: { fat: 1, mot: 2 },                      // the mother's copy doubled
  6: { fat: 0, mot: 1, from: MARKERS / 2 },   // half of it, the father's copy gone
  7: { fat: 0, mot: 0, read: 0.45 },          // both copies gone, the markers fall silent
  8: { fat: 0, mot: 0, read: 1 },             // reads none with the markers still calling: guarded
}
const { src, F, M } = cell(PLAN)
const both = plotData(src, { loaded: { role: 'paternal', gt: F }, second: M })

const barsOf = (d: PlotData, key: Track['key'], chrom: string, per = 1): Bin[] => {
  const t = d.tracks.find((x) => x.key === key)!
  const level = d.tracks.find((x) => x.key === 'intensity')!.bins.get(chrom) ?? []
  return t.kind === 'level' ? totalBars(mergeBins(t.bins.get(chrom) ?? [], per, true))
    : copyBars(mergeBins(t.bins.get(chrom) ?? [], per), level)
}
const levelOf = (d: PlotData, key: Track['key'], chrom: string): string => {
  const v = barsOf(d, key, chrom).map((b) => b.k).filter(Number.isFinite).sort((a, b) => a - b)
  return copyWord(v[v.length >> 1])
}

// ---------------------------------------------------------------- both parents: the two tracks
{
  assert.deepEqual(both.tracks.filter((t) => t.shown).map((t) => t.key), ['father', 'mother'],
    'both parents measured: one track each')
  assert.equal(both.split, true)
  assert.deepEqual(both.tracks.filter((t) => t.shown).map((t) => t.label),
    ['Copies from the father', 'Copies from the mother'])
  for (const c of ['1', '2', '3']) {
    assert.equal(levelOf(both, 'father', c), 'one copy', `chr${c} father`)
    assert.equal(levelOf(both, 'mother', c), 'one copy', `chr${c} mother`)
  }
  const fa = both.tracks[0]
  const c1 = fa.bins.get('1') ?? []
  assert.ok(c1.every((b) => b.endBp - b.startBp >= BASE_BP), 'bars are stretches of the genome')
  // The run of bars along a chromosome is what the reader is asked to see, so there must be a run.
  assert.ok(mergeBins(c1, CELL_PER).length >= 8, `chromosome 1 draws ${mergeBins(c1, CELL_PER).length} bars`)
  assert.equal(both.endBp.get('1'), (MARKERS - 1) * PITCH)
}

// ---------------------------------------------------------------- a loss, a gain, a whole chromosome
{
  // A whole chromosome whose paternal copy is gone: his track at none, hers still at one copy.
  assert.equal(levelOf(both, 'father', '4'), 'none', 'whole chromosome, the father\'s copy gone')
  assert.equal(levelOf(both, 'mother', '4'), 'one copy', "the mother's copy stays")
  // A gain: the doubled track reaches two copies, and the figure's caption says to read that one.
  assert.equal(levelOf(both, 'mother', '5'), 'two copies', "the mother's copy doubled")
  // Half a chromosome: the called half reads none, the other half one copy.
  const half = (MARKERS / 2) * PITCH
  const inside = barsOf(both, 'father', '6').filter((b) => b.startBp >= half)
  const before = barsOf(both, 'father', '6').filter((b) => b.endBp < half)
  const mid = (xs: Bin[]) => xs.map((b) => b.k).sort((a, b) => a - b)[xs.length >> 1]
  assert.ok(mid(inside) < COPY_NONE_BELOW, `the called half reads none (${mid(inside).toFixed(2)})`)
  assert.ok(mid(before) > COPY_NONE_BELOW, `the rest of it does not (${mid(before).toFixed(2)})`)
  // A chromosome lost outright: both tracks at none.
  assert.equal(levelOf(both, 'father', '7'), 'none', 'both copies gone, the father\'s track')
  assert.equal(levelOf(both, 'mother', '7'), 'none', 'both copies gone, the mother\'s track')
}

// ---------------------------------------------------------------- the level bands
{
  assert.equal(copyWord(0), 'none')
  assert.equal(copyWord(COPY_NONE_BELOW - 0.01), 'none')
  assert.equal(copyWord(COPY_NONE_BELOW), 'one copy')
  assert.equal(copyWord(1), 'one copy')
  assert.equal(copyWord(COPY_TWO_ABOVE), 'one copy')
  assert.equal(copyWord(COPY_TWO_ABOVE + 0.01), 'two copies')
  assert.ok(COPY_NONE_BELOW === 0.68 && COPY_TWO_ABOVE === 1.65,
    'the measured bands, from audit/truth/gain-figure-levers-2026-11-17.txt')
}

// ---------------------------------------------------------------- the guard on a level of none
{
  const rate = both.tracks.find((t) => t.key === 'call-rate')!
  assert.ok(rate && !rate.shown, 'the call rate is kept in every setup and drawn as no track')
  const floor = normalRange(rate, 1, new Set(['7', '8'])).lo
  const guardSilent = guardedBars(barsOf(both, 'father', '7'), barsOf(both, 'mother', '7'),
    rate.bins.get('7') ?? [], floor)
  const guardCalling = guardedBars(barsOf(both, 'father', '8'), barsOf(both, 'mother', '8'),
    rate.bins.get('8') ?? [], floor)
  assert.equal(guardSilent.size, 0, 'a chromosome whose markers fell silent is drawn as none')
  assert.ok(guardCalling.size > barsOf(both, 'father', '8').length * 0.8,
    `both tracks read none while the markers keep calling: ${guardCalling.size} bars held back`)
  const kept = guardedBars(barsOf(both, 'father', '4'), barsOf(both, 'mother', '4'),
    rate.bins.get('4') ?? [], floor)
  assert.equal(kept.size, 0, 'one parent at none with the other at one copy is never guarded')
}

// ---------------------------------------------------------------- one parent, and none
{
  const alone = plotData(src, { loaded: { role: 'paternal', gt: F }, second: null })
  assert.equal(alone.split, false)
  assert.deepEqual(alone.tracks.filter((t) => t.shown).map((t) => t.key), ['intensity'],
    'one parent: the copies in all, and no parent named')
  assert.deepEqual(alone.tracks.filter((t) => t.shown).map((t) => t.label), ['Copies in all'])
  const none = plotData(src, { loaded: null, second: null })
  assert.deepEqual(none.tracks.filter((t) => t.shown).map((t) => t.key), ['intensity'], 'no parent: the same')
  assert.equal(levelOf(alone, 'intensity', '1'), 'two copies', 'two in all is a normal chromosome')
  assert.equal(levelOf(alone, 'intensity', '4'), 'one copy', 'a lost copy leaves one in all')
  assert.equal(levelOf(alone, 'intensity', '7'), 'none', 'both copies gone')
  assert.ok(none.tracks.find((t) => t.key === 'call-rate'), 'the guard is kept with no parent too')
  // Which setup it is decides nothing but how many tracks: two where both parents' alleles are
  // known, one otherwise. The mother alone and a rebuilt father behave as any single parent does,
  // and a rebuilt father beside a measured mother splits like any measured pair.
  const motherAlone = plotData(src, { loaded: { role: 'maternal', gt: M }, second: null })
  assert.deepEqual(motherAlone.tracks.filter((t) => t.shown).map((t) => t.key), ['intensity'],
    'the mother alone: the copies in all')
  const rebuilt = new Map([...F].map(([id, g]) => [id, g]))   // a reconstruction is homozygous throughout
  const rebuiltPair = plotData(src, { loaded: { role: 'paternal', gt: rebuilt }, second: M })
  assert.deepEqual(rebuiltPair.tracks.filter((t) => t.shown).map((t) => t.key), ['father', 'mother'],
    'a rebuilt father beside a measured mother: two tracks')
  assert.equal(plotData(src, { loaded: { role: 'paternal', gt: rebuilt }, second: null }).split, false,
    'a rebuilt father alone: one track')
}

// ---------------------------------------------------------------- every change gets its own panel
{
  const call = (chrom: string, mbSpan: number, status?: string) =>
    ({ chrom, startBp: 0, endBp: mbSpan * 1e6, status })
  const many = [call('1', 5), call('2', 90), call('3', 20), call('4', 60), call('5', 7), call('6', 40),
    call('7', 11), call('8', 120), call('9', 3), call('11', 8), call('14', 30), call('16', 4, 'also seen'),
    call('17', 60), call('19', 59)]
  const shown = zoomOrder(many)
  assert.equal(shown.length, many.length, 'a cell with fourteen changes draws fourteen panels')
  assert.deepEqual(shown.map((m) => m.chrom),
    ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '14', '16', '17', '19'], 'in chromosome order')
  assert.ok(shown.some((m) => m.status === 'also seen'),
    'a reading the card does not count keeps its panel, so a chip is never left with nowhere to go')
  const two = zoomOrder([{ chrom: '2', startBp: 5e6 }, { chrom: '2', startBp: 0 }])
  assert.deepEqual(two.map((m) => m.startBp), [0, 5e6], 'two changes on one chromosome keep their order')
  assert.equal(zoomOrder([]).length, 0, 'a cell with no change draws no panel')
}

// ---------------------------------------------------------------- the normal range and bar widths
{
  const fa = both.tracks[0]
  const r = normalRange(fa, 1, new Set(['4', '5', '6', '7', '8']))
  assert.ok(r.lo <= r.mid && r.mid <= r.hi, 'the range is ordered')
  const bins = Array.from({ length: 40 }, (_, i) => ({ startBp: i * 1e6, endBp: (i + 1) * 1e6, n: 20, k: 1 }))
  assert.equal(perForStretch(bins, 0, 40e6), CELL_PER, 'a long stretch takes the whole-cell bar')
  assert.equal(perForStretch(bins, 0, 12e6), 2, 'a shorter one keeps six bars')
  assert.equal(perForStretch(bins, 0, 3e6), 1)
  assert.equal(valueOf({ kind: 'share' }, { startBp: 0, endBp: 1, n: 4, k: 1 }), 0.25)
}

// ---------------------------------------------------------------- the window a zoom shows
assert.deepEqual(viewOf(10e6, 12e6, 135e6), { loBp: 4e6, hiBp: 18e6 })
assert.deepEqual(viewOf(0, 2e6, 135e6), { loBp: 0, hiBp: 14e6 })
assert.deepEqual(viewOf(80e6, 117e6, 135e6), { loBp: 0, hiBp: 135e6 })
assert.equal(tickStepMb(249e6), 50)
assert.deepEqual(wrap('one two three four', 9), ['one two', 'three', 'four'])

// ---------------------------------------------------------------- size and growth
// A real array is about 825,000 markers. Four times the markers costs about four times the time; a
// per-marker rescan would cost sixteen. Best of three draws.
{
  let seed = 11
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x80000000 }
  const sized = (total: number) => {
    const mp = new Map<string, { chrom: string; pos: number }>(); const g = new Map<string, string>()
    const baf = new Map<string, number>()
    const fa = new Map<string, string>(); const mo = new Map<string, string>()
    const cn = new Map<string, { chrom: string; pos: number; called: boolean; log2R: number | null }[]>()
    for (let i = 0; i < total; i += 1) {
      const id = `p${i}`; const chrom = String(1 + (i % 22))
      mp.set(id, { chrom, pos: i * 30 }); g.set(id, rnd() < 0.3 ? 'AB' : 'AA'); baf.set(id, rnd())
      fa.set(id, rnd() < 0.5 ? 'AA' : 'BB'); mo.set(id, rnd() < 0.5 ? 'AA' : 'BB')
      const c = cn.get(chrom) ?? []
      c.push({ chrom, pos: i * 30, called: true, log2R: 0 })
      cn.set(chrom, c)
    }
    return { src: { markerPos: mp, myGt: g, myBaf: baf, cnByChrom: cn }, fa, mo }
  }
  const small = sized(206_000); const large = sized(825_000)
  const best = (s: ReturnType<typeof sized>) => Math.min(...[0, 1, 2].map(() => {
    const t0 = performance.now(); plotData(s.src, { loaded: { role: 'paternal', gt: s.fa }, second: s.mo })
    return performance.now() - t0
  }))
  const ts = best(small); const tl = best(large)
  const d = plotData(large.src, { loaded: { role: 'paternal', gt: large.fa }, second: large.mo })
  let bars = 0
  for (const t of d.tracks) for (const b of t.bins.values()) bars += b.length
  assert.ok(bars <= 22 * 4 * 3.1e9 / (22 * BASE_BP), `${bars} bins kept`)
  if (process.env.OM_TIMING) console.log(`plotData: ${ts.toFixed(0)} ms at 206k, ${tl.toFixed(0)} ms at 825k, ${bars} bins`)
  assert.ok(tl / ts < 8, `four times the markers took ${(tl / ts).toFixed(1)}x the time`)
}


// ---------------------------------------------------------------- what the tracks are called (P8)
// A track label says what the track is, never how it is computed, and no reader-facing string of the
// figure says "letter" for an allele.
{
  const CALC = /share|median|log ?R|homozyg|\bAB\b|against|%|intensity|letter/i
  const alone = plotData(src, { loaded: { role: 'paternal', gt: F }, second: null })
  const none = plotData(src, { loaded: null, second: null })
  for (const d of [both, alone, none]) {
    for (const t of d.tracks) {
      if (t.shown) assert.ok(!CALC.test(t.label), `a track label names a calculation: ${t.label}`)
      assert.ok(!/letter/i.test(t.tech), `"letter" in a hover: ${t.tech}`)
    }
  }
  for (const f of ['CallPlot.tsx', 'plotSeries.ts', 'DefectCallout.tsx', 'pageSays.ts']) {
    const code = readFileSync(new URL(`./${f}`, import.meta.url), 'utf8')
      .split('\n').filter((l) => !/^\s*(\/\/|\*|\/\*\*)/.test(l)).join('\n')
    const said = [...code.matchAll(/(['"`])((?:(?!\1).)*)\1|>([^<>{}]+)</g)].map((m) => (m[2] ?? m[3] ?? '').replace(/\$\{[^}]*\}/g, ''))
    const bad = said.filter((x) => /letter/i.test(x) && !/letterSpacing/.test(x))
    assert.deepEqual(bad, [], `"letter" in a reader-facing string of ${f}`)
  }
  console.log('  track labels say what each track is; no figure string says "letter"')
}

// ---------------------------------------------------------------- the saved file carries every panel
// A saved figure has no control to click, so the one file it writes has to hold the cell and a panel
// for every change. The layout is checked here rather than in the browser: every panel placed, none
// overlapping, none outside the box, and every zoom below the cell.
{
  const sizes = [{ w: 1122, h: 520 }, ...Array.from({ length: 14 }, (_, i) => ({ w: 362, h: 200 + (i % 3) * 12 }))]
  const { width, height, at } = sheet(sizes, 3, 18, 24)
  assert.equal(at.length, sizes.length, 'every panel gets a place in the saved file')
  assert.ok(width >= 1122 && height > 520, 'the box holds the cell and the zooms under it')
  for (const [i, p] of at.entries()) {
    assert.ok(p.x >= 0 && p.y >= 0 && p.x + sizes[i].w <= width && p.y + sizes[i].h <= height,
      `panel ${i} sits inside the saved file`)
    if (i > 0) assert.ok(p.y >= sizes[0].h, `zoom ${i} is below the cell, not over it`)
  }
  for (let i = 0; i < at.length; i += 1) {
    for (let j = i + 1; j < at.length; j += 1) {
      const a = { ...at[i], ...sizes[i] }; const b = { ...at[j], ...sizes[j] }
      const over = a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h
      assert.ok(!over, `panels ${i} and ${j} do not overlap in the saved file`)
    }
  }
  // One change, one panel, and the grid still holds at a width that carries a single column.
  const one = sheet([{ w: 642, h: 480 }, { w: 642, h: 210 }], 1, 18, 24)
  assert.equal(one.at.length, 2, 'a cell with one change saves the cell and that one zoom')
  assert.equal(one.width, 642, 'a single column is as wide as the cell')
  assert.ok(one.at[1].y >= 480 + 24, 'the zoom clears the cell and the line that names the zooms')
}
// ---------------------------------------------------------------- no two badges sit on each other
// Every change carries a number under its own chromosome, so a cell with many changes crowds the
// short chromosomes at the right. A badge may move sideways; it may not land on its neighbour.
{
  const want = [900, 940, 948, 952, 1000, 1004]
  const got = badgeX(want, 17, 9, 1113)
  for (let i = 1; i < got.length; i += 1) {
    assert.ok(got[i] - got[i - 1] >= 17 - 1e-9, `badge ${i} clears the one before it`)
  }
  assert.ok(got[0] >= 9 && got[got.length - 1] <= 1113, 'the badges stay inside the figure')
  assert.deepEqual(badgeX([100, 300], 17, 9, 1113), [100, 300], 'badges that already clear are left alone')
  // A row packed to the right edge is pulled back rather than drawn off the figure.
  const tight = badgeX([1100, 1106, 1112, 1118], 17, 9, 1113)
  assert.ok(tight[3] <= 1113 && tight[0] >= 9, 'a row that will not fit at the edge is pulled back')
  for (let i = 1; i < tight.length; i += 1) assert.ok(tight[i] - tight[i - 1] >= 17 - 1e-9, 'and still clears')
}

// ---------------------------------------------------------------- a change with nothing to draw
// A chromosome the array left with no bins still gets its panel: the caption says it cannot be
// drawn, rather than the change going silently missing from a numbered list.
{
  const track = (bins: Map<string, Bin[]>): Track => ({
    key: 'father', label: 'Copies from the father', tech: '', parent: 'paternal', kind: 'share',
    shown: true, bins,
  } as Track)
  const bin: Bin = { startBp: 0, endBp: 1e6, value: 1, n: 10 } as Bin
  const tracks = [track(new Map([['2', [bin]]]))]
  assert.equal(drawable(tracks, '2'), true, 'a chromosome with bins is drawn')
  assert.equal(drawable(tracks, '3'), false, 'a chromosome with no bins is not')
  assert.equal(drawable([], '2'), false, 'and a figure with no tracks draws nothing')
}

console.log('plotSeries: ok')

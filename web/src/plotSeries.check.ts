// Tracks, bins, the normal range, the running score and which calls are drawn.
import assert from 'node:assert'
import { readFileSync } from 'node:fs'
import {
  BASE, binBase, CELL_BAR, drawable, mergeBins, MIN_DRAWN_SPAN_BP, normalRange, perForStretch, plotData,
  risesInside, runningScore, shareIn, shareTop, tickStepMb, valueOf, viewOf, wrap, type PlotData, type Track,
} from './plotSeries.ts'

// ---------------------------------------------------------------- binning
{
  const pos = Float64Array.from({ length: 260 }, (_, i) => (259 - i) * 10)
  const x = Float64Array.from({ length: 260 }, (_, i) => ((259 - i) % 2 ? 1 : 0))
  const b = binBase(pos, x)
  assert.deepEqual(b.map((y) => y.n), [50, 50, 50, 50, 60], 'a remainder of half a bin or less joins the last')
  assert.equal(b[0].startBp, 0)
  assert.equal(b[0].k, 25)
  const m = mergeBins(b, 4)
  assert.deepEqual(m.map((y) => y.n), [260], 'a short tail joins the bar before it')
  assert.equal(m[0].k, 130, 'counts add exactly')
  assert.deepEqual(mergeBins(b, 2).map((y) => y.n), [100, 100, 60])
  const lv = binBase(Float64Array.from([1, 2, 3]), Float64Array.from([5, -1, 2]), true)
  assert.equal(lv[0].k, 2, 'a level bin is a median')
}

// ---------------------------------------------------------------- tracks
// Two chromosomes; on chromosome 1 from marker 2,000 the father's copy is gone.
const N = 4000
const markerPos = new Map<string, { chrom: string; pos: number }>()
const myGt = new Map<string, string>()
const F = new Map<string, string>()
const M = new Map<string, string>()
let seed = 3
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x80000000 }
for (const chrom of ['1', '2']) {
  for (let i = 0; i < N; i += 1) {
    const id = `${chrom}:${i}`
    markerPos.set(id, { chrom, pos: i * 1000 })
    const f = i % 3 === 0 ? 'AA' : 'BB'
    const m = i % 2 === 0 ? 'BB' : f
    F.set(id, f); M.set(id, m)
    const lost = chrom === '1' && i >= 2000
    const het = f !== m && !lost
    myGt.set(id, het ? (rnd() < 0.15 ? f : 'AB') : lost ? m : f)
  }
}
const src = { markerPos, myGt, cnByChrom: new Map() }
{
  const d = plotData(src, { loaded: { role: 'paternal', gt: F }, second: M })
  assert.deepEqual(d.tracks.filter((t) => t.shown).map((t) => t.key), ['father', 'mother'], 'both parents: one track each')
  const fa = d.tracks[0]
  const c1 = fa.bins.get('1')!
  assert.ok(c1.every((b) => b.n >= BASE && b.n <= BASE * 1.5), 'equal-count bars')
  const inside = shareIn(c1, 2_000_000, 4_000_000)
  const outside = shareIn(fa.bins.get('2')!, 0, 4e6)
  assert.ok(inside.share < 0.05 && outside.share > 0.9, `father inside ${inside.share} outside ${outside.share}`)
  assert.ok(shareIn(d.tracks[1].bins.get('1')!, 2e6, 4e6).share > 0.95, "the mother's alleles stay")
  assert.equal(d.endBp.get('1'), (N - 1) * 1000)
  const r = normalRange(fa, 1, new Set(['1']))
  assert.ok(r.lo <= r.mid && r.mid <= r.hi && r.lo > 0.7, 'the cell\'s normal range from chromosomes with no call')
}
{
  // Father alone: his allele found where the other shows; heterozygous calls stand for the mother.
  const d = plotData(src, { loaded: { role: 'paternal', gt: F }, second: null })
  assert.deepEqual(d.tracks.filter((t) => t.shown).map((t) => t.key), ['father', 'two-letter'])
  assert.ok(Number.isNaN(d.dropIn))
  // With the run reading her: her track replaces the heterozygous calls.
  const opp = new Set([...markerPos.keys()].filter((id) => F.get(id) !== M.get(id)))
  const r = plotData(src, { loaded: { role: 'paternal', gt: F }, second: null, otherOpposite: opp })
  assert.deepEqual(r.tracks.filter((t) => t.shown).map((t) => t.key), ['father', 'mother'])
}
{
  const cn = new Map([['1', Array.from({ length: 400 }, (_, i) => ({ pos: i * 1000, called: i < 200, log2R: i < 200 ? 0 : -1 }))]])
  const d = plotData({ markerPos: new Map(), myGt: new Map(), cnByChrom: cn }, { loaded: null, second: null })
  assert.deepEqual(d.tracks.map((t) => t.key), ['call-rate', 'intensity'])
  const it = d.tracks[1]
  assert.ok(valueOf(it, it.bins.get('1')!.at(-1)!) < -0.5, 'intensity against the cell\'s own level')
}

// ---------------------------------------------------------------- running score
{
  const bins = [...Array(10)].map((_, i) => ({ startBp: i, endBp: i, n: 100, k: i >= 5 ? 5 : 85 }))
  const s = runningScore(bins, 0.85, 0.03)
  assert.ok(s.slice(0, 5).every((v) => v === 0), 'flat where the copy is there')
  assert.ok(s[9] > s[6] && s[6] > s[5] && s[5] > 0, 'climbs where it is not')
  assert.ok(runningScore(bins, 0.01, 0.03).every(Number.isNaN), 'no score where present is not above missing')
}

// ---------------------------------------------------------------- axis top and when the score is drawn
assert.equal(shareTop(0.12), 0.24, 'a low row is drawn against twice its normal level')
assert.equal(shareTop(0.85), 1, 'a high row keeps 0 to 100%')
assert.equal(shareTop(NaN), 1)
{
  const bins = [...Array(10)].map((_, i) => ({ startBp: i * 10, endBp: i * 10 + 9, n: 100, k: 0 }))
  const climbsInside = [0, 0, 0, 2, 4, 6, 8, 8, 7, 6]
  const climbsAfter = [0, 0, 0, 0, 0, 0, 0, 3, 6, 9]
  assert.equal(risesInside(climbsInside, bins, 30, 69), true)
  assert.equal(risesInside(climbsAfter, bins, 30, 69), false, 'a rise after the stretch ends is not drawn')
}

// ---------------------------------------------------------------- which calls are drawn
{
  const t = (key: Track['key'], parent: Track['parent']): Track => ({ key, parent, label: '', tech: '', kind: 'share', shown: true, bins: new Map() })
  const both: PlotData = { tracks: [t('father', 'paternal'), t('mother', 'maternal')], dropIn: 0.02, endBp: new Map() }
  const alone: PlotData = { tracks: [t('father', 'paternal'), t('two-letter', null)], dropIn: NaN, endBp: new Map() }
  const none: PlotData = { tracks: [t('call-rate', null), t('intensity', null)], dropIn: NaN, endBp: new Map() }
  const piece = (parent: 'paternal' | 'maternal' | null, gain = false, spanBp = 20e6) => ({ whole: false, spanBp, gain, parent })
  assert.equal(drawable(piece('maternal'), both).track?.key, 'mother')
  assert.equal(drawable(piece('maternal'), alone).draw, false, "her copy on a piece, father alone, no run: not drawn")
  assert.equal(drawable(piece('paternal'), alone).track?.key, 'father')
  assert.equal(drawable(piece('paternal', true), both).draw, false, 'a gain on a piece: not drawn')
  assert.equal(drawable(piece(null), none).draw, false, 'a piece with no parent: not drawn')
  assert.equal(drawable(piece('paternal', false, MIN_DRAWN_SPAN_BP - 1), both).draw, false, 'under one band: not drawn')
  assert.equal(drawable({ whole: true, spanBp: 1e8, gain: false, parent: null }, none).track?.key, 'intensity')
  assert.equal(drawable({ whole: true, spanBp: 1e8, gain: false, parent: 'maternal' }, alone).track?.key, 'two-letter')
}

// ---------------------------------------------------------------- geometry and text
{
  const bins = [...Array(40)].map((_, i) => ({ startBp: i * 1e6, endBp: i * 1e6 + 9e5, n: 50, k: 25 }))
  assert.equal(perForStretch(bins, 0, 40e6), CELL_BAR / BASE, 'a long stretch takes the whole-cell bar')
  assert.equal(perForStretch(bins, 0, 12e6), 2, 'a shorter one keeps six bars')
  assert.equal(perForStretch(bins, 0, 3e6), 1)
}
assert.deepEqual(viewOf(143.1e6, 144.7e6, 159e6), { loBp: 137.9e6, hiBp: 149.9e6 })
assert.deepEqual(viewOf(80e6, 117e6, 135e6), { loBp: 0, hiBp: 135e6 })
assert.equal(tickStepMb(249e6), 50)
assert.deepEqual(wrap('one two three four', 9), ['one two', 'three', 'four'])

// ---------------------------------------------------------------- size and growth
// A real array is about 825,000 markers. Four times the markers costs about four times the time; a
// per-marker rescan would cost sixteen. Best of three draws.
{
  const sized = (total: number) => {
    const mp = new Map<string, { chrom: string; pos: number }>(); const g = new Map<string, string>()
    const fa = new Map<string, string>(); const mo = new Map<string, string>()
    for (let i = 0; i < total; i += 1) {
      const id = `p${i}`; const chrom = String(1 + (i % 22))
      mp.set(id, { chrom, pos: i * 30 }); g.set(id, rnd() < 0.3 ? 'AB' : 'AA')
      fa.set(id, rnd() < 0.5 ? 'AA' : 'BB'); mo.set(id, rnd() < 0.5 ? 'AA' : 'BB')
    }
    return { src: { markerPos: mp, myGt: g, cnByChrom: new Map() }, fa, mo }
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
  assert.ok(bars <= 2 * 825_000 / BASE, `${bars} bins kept`)
  if (process.env.OM_TIMING) console.log(`plotData: ${ts.toFixed(0)} ms at 206k, ${tl.toFixed(0)} ms at 825k, ${bars} bins`)
  assert.ok(tl / ts < 8, `four times the markers took ${(tl / ts).toFixed(1)}x the time`)
}

console.log('plotSeries: ok')

// ---------------------------------------------------------------- what the rows are called (P8)
// A row label says what the row is, never how it is computed, and no reader-facing string of the
// figures says "letter" for an allele.
{
  const CALC = /share|median|log ?R|homozyg|\bAB\b|against|%|intensity|letter/i
  const labels = (d: PlotData) => d.tracks.filter((t) => t.shown).map((t) => t.label)
  const both = plotData(src, { loaded: { role: 'paternal', gt: F }, second: M })
  const alone = plotData(src, { loaded: { role: 'paternal', gt: F }, second: null })
  const opp = new Set([...markerPos.keys()].filter((id) => F.get(id) !== M.get(id)))
  const run = plotData(src, { loaded: { role: 'paternal', gt: F }, second: null, otherOpposite: opp })
  const cn = new Map([['1', Array.from({ length: 400 }, (_, i) => ({ pos: i * 1000, called: true, log2R: 0 }))]])
  const none = plotData({ markerPos: new Map(), myGt: new Map(), cnByChrom: cn }, { loaded: null, second: null })
  const noneAll = none.tracks.map((t) => t.label)
  assert.deepEqual(labels(both), ["Father's copy", "Mother's copy"])
  assert.deepEqual(labels(alone), ["Father's copy", "Mother's copy, indirect"])
  assert.deepEqual(labels(run), ["Father's copy", "Mother's copy"])
  assert.deepEqual(noneAll, ['DNA present', 'Copy number'])
  for (const d of [both, alone, run, none]) {
    for (const t of d.tracks) {
      assert.ok(!CALC.test(t.label), `a row label names a calculation: ${t.label}`)
      assert.ok(!/letter/i.test(t.tech), `"letter" in a hover: ${t.tech}`)
    }
  }
  // Every reader-facing string literal and JSX text of the figure files.
  for (const f of ['CallPlot.tsx', 'plotSeries.ts', 'DefectCallout.tsx', 'pageSays.ts']) {
    const code = readFileSync(new URL(`./${f}`, import.meta.url), 'utf8')
      .split('\n').filter((l) => !/^\s*(\/\/|\*|\/\*\*)/.test(l)).join('\n')
    const said = [...code.matchAll(/(['"`])((?:(?!\1).)*)\1|>([^<>{}]+)</g)].map((m) => (m[2] ?? m[3] ?? '').replace(/\$\{[^}]*\}/g, ''))
    const bad = said.filter((x) => /letter/i.test(x) && !/^two-letter$/.test(x) && !/letterSpacing/.test(x))
    assert.deepEqual(bad, [], `"letter" in a reader-facing string of ${f}`)
  }
  console.log('  row labels say what each row is; no figure string says "letter"')
}

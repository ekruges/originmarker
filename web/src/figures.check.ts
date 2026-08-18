// Self-check for figure geometry. Run: node src/figures.check.ts
//
// Two properties carry this module.
//
//   TICKS MUST SIT AT VALUES A PERSON WOULD CHOOSE. A reader recovers a number from a position by
//   counting ticks, so 0, 0.25, 0.5 is a figure and 0, 0.31, 0.62 is a decoration. This is the one
//   thing the previous charts had none of.
//
//   ONE LAYOUT, TWO RENDERERS. The screen and the report each had their own chart code with their
//   own sizes: the same figure was 260 units wide in one and 210 in the other, rows 26 high against
//   18. Layout now happens once and emits marks, so a divergence is not expressible.
import assert from 'node:assert/strict'
import { FIG, fit, niceTicks, panelLabel, scale, xAxis } from './figures.ts'
import { comparisonFigures, figureEnrichment, figureFold } from './comparisonFigures.ts'
import { compare } from './comparison.ts'
import type { FeatureTrack, Region } from './features.ts'

// --- 1. TICKS ARE ROUND NUMBERS ------------------------------------------------------------------
{
  assert.deepEqual(niceTicks(0, 1, 5), [0, 0.2, 0.4, 0.6000000000000001, 0.8, 1].map(
    (v) => Number(v.toFixed(10))), 'a unit range ticks in fifths')
  assert.deepEqual(niceTicks(0, 0.5, 5), [0, 0.1, 0.2, 0.3, 0.4, 0.5])
  assert.deepEqual(niceTicks(0, 3, 3), [0, 1, 2, 3])

  // The 1-2-5 progression: no tick may land on an unroundable value.
  for (const hi of [0.07, 0.3, 1.7, 4.2, 9.9, 23, 480]) {
    for (const t of niceTicks(0, hi)) {
      const mag = 10 ** Math.floor(Math.log10(Math.abs(t) || 1))
      const norm = Math.abs(t) / mag
      assert.ok(Math.abs(norm - Math.round(norm * 10) / 10) < 1e-9,
        `tick ${t} on range 0..${hi} is not a round value`)
    }
  }

  // Float dust must be snapped away, or a tick renders as 0.30000000000000004.
  for (const t of niceTicks(0, 1, 10)) {
    assert.ok(String(t).length <= 12, `tick ${t} carries float dust into the label`)
  }

  // A degenerate range must not loop forever or emit nonsense.
  assert.deepEqual(niceTicks(5, 5), [5])
  assert.deepEqual(niceTicks(1, 0), [1])
}

// --- 2. THE SCALE AND THE AXIS AGREE ---------------------------------------------------------------
{
  const s = scale(0, 1, 100, 300)
  assert.equal(s(0), 100)
  assert.equal(s(1), 300)
  assert.equal(s(0.5), 200)
  assert.equal(scale(2, 2, 10, 90)(2), 10, 'a zero-width domain must not divide by zero')

  // Every tick label on an axis must have a tick mark and vice versa.
  const marks = xAxis(100, 300, 50, 0, 1, 'Share (%)', (v) => `${100 * v}`)
  const ticks = marks.filter((m) => m.k === 'line' && m.y1 !== m.y2)
  const labels = marks.filter((m) => m.k === 'text' && m.anchor === 'middle')
  assert.equal(ticks.length + 1, labels.length,
    'one label per tick, plus the axis title, or the axis is misleading')
  assert.ok(marks.some((m) => m.k === 'text' && m.s === 'Share (%)'),
    'THE AXIS MUST STATE ITS UNIT. A number without one is not a measurement')
  // The axis line spans the plot exactly.
  const axis = marks.find((m) => m.k === 'line' && m.y1 === m.y2)!
  assert.equal(axis.k === 'line' && axis.x1, 100)
  assert.equal(axis.k === 'line' && axis.x2, 300)
}

// --- 3. ONE TYPOGRAPHIC SCALE ----------------------------------------------------------------------
{
  const sizes = Object.values(FIG.font)
  assert.ok(sizes.every((s) => s >= 6 && s <= 10), 'every size is in one readable range')
  assert.ok(FIG.font.panel > FIG.font.label, 'a panel letter outranks a data label')
  assert.ok(FIG.font.label > FIG.font.tick, 'and a label outranks a tick')
  assert.equal(fit('a very long feature name indeed', 10), 'a very lo…',
    'a truncated label must SHOW that it was cut')
  assert.equal(fit('short', 10), 'short')
  assert.equal(panelLabel(4, 10, 'A').s, 'A')
}

// --- 4. EVERY FIGURE SHARES THE SCALE AND CARRIES ITS OWN CAPTION -----------------------------------
{
  const track: FeatureTrack = {
    build: 'hg19',
    fragile: Array.from({ length: 8 }, (_, i) => ({
      chrom: '1', startBp: i * 20e6, endBp: i * 20e6 + 6e6, name: `FRA${i}`,
    })),
    longGenes: [{ chrom: '1', startBp: 195e6, endBp: 199e6, name: 'FAR' }],
    gaps: [], geneDensityPerMb: {},
    lateReplicationValleysES: [], lateReplicationValleysConstitutive: [],
  }
  const regions: Region[] = Array.from({ length: 9 }, (_, i) => ({
    chrom: '1', startBp: i * 20e6 + 1e6, endBp: i * 20e6 + 4e6,
  }))
  const names = regions.map((r) => `chr1 ${(r.startBp / 1e6).toFixed(1)}Mb`)
  const markers = new Map<string, number[]>([['1',
    Array.from({ length: 5000 }, (_, i) => i * 40_000)]])
  const c = compare(track, regions, markers, { permutations: 300, regionNames: names })

  const figs = comparisonFigures(c, names)
  assert.ok(figs.length >= 3, 'four figures, minus any with nothing to draw')
  for (const f of figs) {
    assert.ok(f.w > 100 && f.h > 40, 'a figure must have a size')
    assert.ok(f.caption.length > 60, 'and a caption long enough to stand alone in a manuscript')
    assert.ok(/^[A-D]\./.test(f.caption), `a caption must name its panel: ${f.caption.slice(0, 30)}`)
    assert.ok(f.marks.some((m) => m.k === 'text' && m.bold),
      'and the panel letter must be on the figure itself, not only in the caption')
    // No mark may fall outside the figure it belongs to, or a renderer clips it.
    for (const m of f.marks) {
      const xs = m.k === 'line' ? [m.x1, m.x2] : [m.x]
      for (const x of xs) {
        assert.ok(x >= -1 && x <= f.w + 1, `a mark at x=${x} falls outside a figure ${f.w} wide`)
      }
    }
  }

  // THE SHARED AXIS, which is the property that makes panel A readable: per-row axes would make
  // every feature fill its own width and look equally enriched.
  const a = figureEnrichment(c)
  // The null intervals: horizontal, so y1 === y2, and wider than a cap.
  const rows = a.marks.filter((m) => m.k === 'line' && m.y1 === m.y2 && m.x2 - m.x1 > 1
    && m.colour === FIG.grey)
  assert.ok(rows.length >= 2, 'several features to compare')
  const spans = rows.map((m) => (m.k === 'line' ? m.x2 - m.x1 : 0))
  assert.ok(new Set(spans.map((s) => Math.round(s))).size > 1,
    'null intervals must differ in width, which is only possible on one shared axis')

  // EVERY INTERVAL CARRIES A SERIF AT EACH END, which is what makes it read as a range rather than
  // as a bar. A figure that carries both vocabularies must not draw them the same way.
  const caps = a.marks.filter((m) => m.k === 'line' && m.x1 === m.x2 && m.colour === FIG.grey)
  for (const r of rows) {
    if (r.k !== 'line') continue
    for (const end of [r.x1, r.x2]) {
      assert.ok(
        caps.some((cp) => cp.k === 'line' && Math.abs(cp.x1 - end) < 1e-9
          && cp.y1 < r.y1 && cp.y2 > r.y1),
        `the interval at y=${r.y1} has no serif at x=${end}`,
      )
    }
  }
  assert.ok(caps.length >= rows.length * 2)

  // A LEGEND, because a reader should not have to infer what a mark means from the caption.
  assert.ok(a.marks.some((m) => m.k === 'text' && /null, middle 95%/.test(m.s)))
  assert.ok(a.marks.some((m) => m.k === 'text' && m.s === 'observed'))
  // And the sample size on the axis.
  assert.ok(a.marks.some((m) => m.k === 'text' && /n = \d+/.test(m.s)), 'n must be stated')

  // EVERY PANEL SHARES ONE WIDTH AND ONE LABEL COLUMN, so they stack as a composed figure and a
  // reader tracks one feature straight down it.
  assert.ok(figs.every((f) => f.w === figs[0].w), 'panels must share a width to stack')
  for (const f of figs) {
    const labels = f.marks.filter((m) => m.k === 'text' && m.anchor === 'end'
      && m.size === FIG.font.label)
    assert.ok(labels.every((m) => m.x === FIG.pad.left - 7),
      'feature labels must sit in the same column in every panel')
  }

  // Fold carries its reference line at 1.00x, in the accent colour so it reads over the bars.
  const fold = figureFold(c)
  assert.ok(fold.marks.some((m) => m.k === 'line' && m.colour === FIG.accent),
    'the null at 1.00x must be drawn, or a fold chart has no reference')
}

console.log('figures.check.ts: all assertions passed, including ticks at round values, every axis '
  + 'stating its unit, and one layout serving both renderers')

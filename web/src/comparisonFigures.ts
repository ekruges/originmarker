/**
 * The four comparison figures, laid out once.
 *
 * Each returns a Figure: a width, a height, a list of marks and a caption generated from the data
 * that produced it. Nothing here draws; the screen and the report each render the same marks, which
 * is what stops the two versions of a figure drifting apart.
 *
 * FOUR FIGURES BECAUSE THERE ARE FOUR QUESTIONS, and no more than that. Does the observation sit
 * outside what the matched null produces (A). What does that null actually look like, since a p
 * value is only a summary of it (B). How large is the effect, on one scale so features can be
 * ranked (C). And which specific region touched which feature, since an enrichment can rest on two
 * regions out of twenty and a p value hides that (D).
 */
import {
  FIG, fit, niceTicks, panelLabel, scale, xAxis, type Figure, type Mark,
} from './figures.ts'
import { enrichmentBars, foldChart, nullHistograms, type ComparisonResult } from './comparison.ts'

const pct = (v: number) => `${Math.round(100 * v)}%`

/**
 * A. Observed overlap against the matched null, one row per feature, one shared axis.
 *
 * The bar is the middle 95% of the null the permutation actually drew and the dot is the
 * observation. A SHARED AXIS ACROSS ROWS is the only reason this is readable: per-row axes make
 * every feature fill its own width and look equally enriched, which is the easiest way to draw a
 * chart that says the opposite of its data.
 */
export function figureEnrichment(c: ComparisonResult): Figure {
  const bars = enrichmentBars(c)
  const plotW = 210
  const w = FIG.pad.left + plotW + FIG.pad.right
  const h = FIG.pad.top + bars.length * FIG.rowH + FIG.pad.bottom
  const x0 = FIG.pad.left
  const x1 = x0 + plotW
  const top = FIG.pad.top
  const dMax = Math.max(0.05, ...bars.map((b) => Math.max(b.observed, b.hi)))
  const sx = scale(0, dMax, x0, x1)
  const marks: Mark[] = [panelLabel(4, 10, 'A')]

  // Gridlines first, so nothing sits on top of the data.
  for (const t of niceTicks(0, dMax)) {
    marks.push({ k: 'line', x1: sx(t), y1: top - 4, x2: sx(t), y2: top + bars.length * FIG.rowH,
      colour: FIG.line, width: 0.4, dash: true })
  }
  bars.forEach((b, i) => {
    const y = top + i * FIG.rowH + FIG.rowH / 2
    marks.push({ k: 'text', x: x0 - 6, y: y + 2.4, s: fit(b.label, 30), size: FIG.font.label,
      colour: FIG.ink, anchor: 'end' })
    marks.push({ k: 'line', x1: sx(b.lo), y1: y, x2: sx(b.hi), y2: y, colour: FIG.grey, width: 3.2 })
    marks.push({ k: 'dot', x: sx(b.observed), y, r: 2.6,
      colour: b.significant ? FIG.accent : FIG.ink })
    // Significance as a SYMBOL as well as colour, so the figure survives grey printing.
    marks.push({ k: 'text', x: x1 + 6, y: y + 2.4,
      s: `${pct(b.observed)} vs ${pct(b.expected)}${b.significant ? '  *' : ''}`,
      size: FIG.font.tick, colour: b.significant ? FIG.accent : FIG.grey, anchor: 'start' })
  })
  marks.push(...xAxis(x0, x1, top + bars.length * FIG.rowH, 0, dMax,
    'Regions overlapping the feature (%)', (v) => `${Math.round(100 * v)}`))

  const hit = bars.filter((b) => b.significant).length
  return {
    w,
    h,
    marks,
    caption: `A. Share of the ${c.regions} detected region${c.regions === 1 ? '' : 's'} overlapping `
      + `each feature (dot) against the middle 95% of a matched null (bar), over `
      + `${c.permutations.toLocaleString()} permutations. `
      + (hit ? `${hit} feature${hit === 1 ? '' : 's'} marked * clear${hit === 1 ? 's' : ''} the `
        + 'corrected threshold.' : 'No feature clears the corrected threshold.'),
  }
}

/**
 * B. The permutation null itself, with the observation marked in it.
 *
 * The figure a permutation test is obliged to show. A p value is a summary of this distribution and
 * summaries mislead in both directions: a ratio near one reaches a small p when the null is tight,
 * and a large ratio reaches nothing when it is broad. An observation the permutation never produced
 * is drawn at the edge and said to be outside, which is the strongest result available here.
 */
export function figureNulls(c: ComparisonResult): Figure {
  const hs = nullHistograms(c)
  const cellW = 116
  const cellH = 42
  const perRow = 3
  const rows = Math.ceil(hs.length / perRow)
  const w = 24 + perRow * cellW
  const h = FIG.pad.top + rows * (cellH + 30) + 10
  const marks: Mark[] = [panelLabel(4, 10, 'B')]

  hs.forEach((hst, i) => {
    const cx = 24 + (i % perRow) * cellW
    const cy = FIG.pad.top + Math.floor(i / perRow) * (cellH + 30)
    const plotW = cellW - 18
    const max = Math.max(1, ...hst.bins)
    const bw = plotW / hst.bins.length
    marks.push({ k: 'text', x: cx, y: cy - 3, s: fit(hst.label, 24), size: FIG.font.tick,
      colour: FIG.ink, anchor: 'start' })
    hst.bins.forEach((v, k) => {
      const bh = (v / max) * cellH
      if (bh > 0.3) {
        marks.push({ k: 'rect', x: cx + k * bw, y: cy + cellH - bh,
          w: Math.max(0.6, bw - 0.5), h: bh, colour: FIG.line })
      }
    })
    marks.push({ k: 'line', x1: cx + hst.observedAt * plotW, y1: cy, x2: cx + hst.observedAt * plotW,
      y2: cy + cellH, colour: hst.significant ? FIG.accent : FIG.ink, width: 1.2 })
    marks.push(...xAxis(cx, cx + plotW, cy + cellH, hst.lo, hst.hi, '',
      (v) => `${Math.round(100 * v)}`, 2))
    marks.push({ k: 'text', x: cx, y: cy + cellH + 20,
      s: hst.observedBin === -1 ? 'observed outside the null' : `p = ${hst.p.toExponential(1)}`,
      size: FIG.font.tick, colour: hst.significant ? FIG.accent : FIG.grey, anchor: 'start' })
  })
  const outside = hs.filter((x) => x.observedBin === -1).length
  return {
    w,
    h,
    marks,
    caption: 'B. The permutation null for each feature (bars) with the observation marked (line), '
      + 'on the same axis as A. A p value is a summary of this distribution and is only '
      + 'interpretable against it: a small p can come from a tight null rather than a large effect.'
      + (outside ? ` ${outside} observation${outside === 1 ? '' : 's'} fell outside the null `
        + 'entirely.' : ''),
  }
}

/** C. Fold enrichment, so features can be ranked on one scale against a reference line at 1. */
export function figureFold(c: ComparisonResult): Figure {
  const rows = foldChart(c)
  const plotW = 210
  const w = FIG.pad.left + plotW + FIG.pad.right
  const h = FIG.pad.top + rows.length * FIG.rowH + FIG.pad.bottom
  const x0 = FIG.pad.left
  const x1 = x0 + plotW
  const top = FIG.pad.top
  const dMax = Math.max(2, ...rows.map((r) => r.fold))
  const sx = scale(0, dMax, x0, x1)
  const marks: Mark[] = [panelLabel(4, 10, 'C')]

  for (const t of niceTicks(0, dMax)) {
    marks.push({ k: 'line', x1: sx(t), y1: top - 4, x2: sx(t), y2: top + rows.length * FIG.rowH,
      colour: FIG.line, width: 0.4, dash: true })
  }
  rows.forEach((r, i) => {
    const y = top + i * FIG.rowH + FIG.rowH / 2
    marks.push({ k: 'text', x: x0 - 6, y: y + 2.4, s: fit(r.label, 30), size: FIG.font.label,
      colour: FIG.ink, anchor: 'end' })
    marks.push({ k: 'rect', x: x0, y: y - 3.4, w: Math.max(0.6, sx(r.fold) - x0), h: 6.8,
      colour: r.fold >= 1 ? FIG.ink : FIG.grey })
    marks.push({ k: 'text', x: x1 + 6, y: y + 2.4, s: `${r.fold.toFixed(2)}x`,
      size: FIG.font.tick, colour: FIG.grey, anchor: 'start' })
  })
  // The reference line last, so it reads over the bars rather than under them.
  marks.push({ k: 'line', x1: sx(1), y1: top - 4, x2: sx(1), y2: top + rows.length * FIG.rowH,
    colour: FIG.accent, width: 0.9 })
  marks.push(...xAxis(x0, x1, top + rows.length * FIG.rowH, 0, dMax,
    'Fold enrichment over the matched null', (v) => v.toFixed(v < 1 ? 1 : 0)))
  return {
    w,
    h,
    marks,
    caption: 'C. Observed overlap divided by what the matched null expects. The vertical line is '
      + 'the null at 1.00x; bars left of it overlap the feature LESS than chance.',
  }
}

/** D. Which region touched which feature, which is what an enrichment is made of. */
export function figureGrid(c: ComparisonResult, names: readonly string[]): Figure {
  const colW = 15
  const labelW = 150
  const w = labelW + Math.max(1, c.features.length) * colW + 20
  const headH = 62
  const h = FIG.pad.top + headH + names.length * 11 + 26
  const top = FIG.pad.top + headH
  const marks: Mark[] = [panelLabel(4, 10, 'D')]

  c.features.forEach((f, i) => {
    // Rotated headers are not available to both renderers, so the header is stepped instead: each
    // label sits on its own line with a leader down to its column.
    const x = labelW + i * colW + colW / 2
    marks.push({ k: 'text', x: labelW - 6, y: FIG.pad.top + 8 + i * 10, s: fit(f.label, 30),
      size: FIG.font.tick, colour: FIG.grey, anchor: 'end' })
    marks.push({ k: 'line', x1: labelW - 3, y1: FIG.pad.top + 5 + i * 10, x2: x, y2: top - 3,
      colour: FIG.line, width: 0.4 })
  })
  names.forEach((n, r) => {
    const y = top + r * 11 + 7
    marks.push({ k: 'text', x: labelW - 6, y, s: fit(n, 30), size: FIG.font.tick,
      colour: FIG.ink, anchor: 'end' })
    c.features.forEach((f, i) => {
      const x = labelW + i * colW + colW / 2
      const hit = !!f.regionHits[r]
      marks.push({ k: 'dot', x, y: y - 2.4, r: hit ? 2.6 : 1, colour: hit ? FIG.accent : FIG.line,
        hollow: false })
    })
  })
  const touched = names.filter((_, i) => c.features.some((f) => f.regionHits[i])).length
  return {
    w,
    h,
    marks,
    caption: `D. Which region overlaps which feature. ${touched} of ${names.length} region`
      + `${names.length === 1 ? '' : 's'} touch at least one feature. An enrichment can rest on a `
      + 'handful of regions, and this is where that becomes visible rather than hiding in a p value.',
  }
}

/** Every figure, in order, so both renderers draw the same set. */
export const comparisonFigures = (c: ComparisonResult, names: readonly string[]): Figure[] => [
  figureEnrichment(c), figureNulls(c), figureFold(c), figureGrid(c, names),
].filter((f) => f.marks.length > 1)

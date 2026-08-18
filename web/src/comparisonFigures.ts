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
  FIG, FIG_W, fit, interval, legend, niceTicks, panelLabel, scale, xAxis, yAxis,
  type Figure, type Mark,
} from './figures.ts'
import {
  enrichmentBars, foldChart, nullHistograms, regionFlags, scannedTracks,
  type ComparisonResult,
} from './comparison.ts'

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
  const w = FIG_W
  const legendH = 12
  const untestedRows = c.features.filter((f) => !f.testable).length ? 12 : 0
  const h = FIG.pad.top + legendH + bars.length * FIG.rowH + FIG.pad.bottom + untestedRows
  const x0 = FIG.pad.left
  const x1 = x0 + FIG.plotW
  const top = FIG.pad.top + legendH
  const dMax = Math.max(0.05, ...bars.map((b) => Math.max(b.observed, b.hi)))
  const sx = scale(0, dMax, x0, x1)
  const marks: Mark[] = [panelLabel(4, 11, 'A')]

  marks.push(...legend(x0, FIG.pad.top + 2, [
    { kind: 'interval', colour: FIG.grey, s: 'null, middle 95%' },
    { kind: 'dot', colour: FIG.ink, s: 'observed' },
  ]))
  for (const t of niceTicks(0, dMax)) {
    marks.push({ k: 'line', x1: sx(t), y1: top, x2: sx(t), y2: top + bars.length * FIG.rowH,
      colour: FIG.line, width: 0.4, dash: true })
  }
  bars.forEach((b, i) => {
    const y = top + i * FIG.rowH + FIG.rowH / 2
    marks.push({ k: 'text', x: x0 - 7, y: y + 2.6, s: fit(b.label, 32), size: FIG.font.label,
      colour: FIG.ink, anchor: 'end' })
    marks.push(...interval(sx(b.lo), sx(b.hi), y))
    marks.push({ k: 'dot', x: sx(b.observed), y, r: 2.8,
      colour: b.significant ? FIG.accent : FIG.ink })
    marks.push({ k: 'text', x: x1 + 6, y: y + 2.6,
      s: `${pct(b.observed)}${b.significant ? ' *' : ''}`,
      size: FIG.font.tick, colour: b.significant ? FIG.accent : FIG.grey, anchor: 'start' })
  })
  marks.push(...xAxis(x0, x1, top + bars.length * FIG.rowH, 0, dMax,
    `Regions overlapping the feature (%), n = ${c.regions}`, (v) => `${Math.round(100 * v)}`))

  const untested = c.features.filter((f) => !f.testable)
  if (untested.length) {
    marks.push({ k: 'text', x: x0 - 7,
      y: top + bars.length * FIG.rowH + FIG.tick.len + FIG.font.tick + 22,
      s: `not testable: ${untested.map((f) => f.label).join(', ')}`,
      size: FIG.font.tick, colour: FIG.grey, anchor: 'end' })
  }

  const hit = bars.filter((b) => b.significant).length
  return {
    w,
    h,
    marks,
    caption: `A. Share of the ${c.regions} detected region${c.regions === 1 ? '' : 's'} overlapping `
      + 'each feature (filled circle) against the middle 95% of a matched null (interval), over '
      + `${c.permutations.toLocaleString()} permutations. The null draws intervals on the same `
      + 'chromosome carrying the same number of informative markers as the region tested. '
      + (hit ? `Asterisk: clears the Bonferroni-corrected threshold (${hit} of ${bars.length}).`
        : 'No feature clears the Bonferroni-corrected threshold.'),
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
  const cellH = 40
  const gap = 34
  const w = FIG_W
  const h = FIG.pad.top + hs.length * (cellH + gap) + 6
  const x0 = FIG.pad.left
  const plotW = FIG.plotW
  const marks: Mark[] = [panelLabel(4, 11, 'B')]

  hs.forEach((hst, i) => {
    const top = FIG.pad.top + i * (cellH + gap)
    const base = top + cellH
    const max = Math.max(1, ...hst.bins)
    const bw = plotW / hst.bins.length
    // The feature name sits in the SAME left column as panels A, C and D, so one feature reads
    // straight down the composed figure.
    marks.push({ k: 'text', x: x0 - 7, y: top + cellH / 2, s: fit(hst.label, 32),
      size: FIG.font.label, colour: FIG.ink, anchor: 'end' })
    hst.bins.forEach((v, k) => {
      const bh = (v / max) * cellH
      if (bh > 0.3) {
        marks.push({ k: 'rect', x: x0 + k * bw, y: base - bh, w: Math.max(0.6, bw - 0.5), h: bh,
          colour: FIG.line })
      }
    })
    // A histogram with no frequency axis asks a reader to compare heights they have no scale for.
    marks.push(...yAxis(x0, top, base, 0, max, 'permutations', (v) => String(Math.round(v)), 2))
    marks.push({ k: 'line', x1: x0 + hst.observedAt * plotW, y1: top - 3,
      x2: x0 + hst.observedAt * plotW, y2: base,
      colour: hst.significant ? FIG.accent : FIG.ink, width: 1.4 })
    marks.push(...xAxis(x0, x0 + plotW, base, hst.lo, hst.hi,
      i === hs.length - 1 ? 'Regions overlapping the feature (%)' : '',
      (v) => `${Math.round(100 * v)}`, 3))
    marks.push({ k: 'text', x: x0 + plotW + 6, y: top + cellH / 2,
      s: hst.observedBin === -1 ? 'outside null' : `p = ${hst.p.toExponential(1)}`,
      size: FIG.font.tick, colour: hst.significant ? FIG.accent : FIG.grey, anchor: 'start' })
  })
  const outside = hs.filter((x) => x.observedBin === -1).length
  return {
    w,
    h,
    marks,
    caption: 'B. The permutation null for each feature (grey bars, frequency on the left) with the '
      + 'observation marked (vertical line), on the same axis as A. A p value is a summary of this '
      + 'distribution and is only interpretable against it: a small p can come from a tight null '
      + 'rather than a large effect.'
      + (outside ? ` ${outside} observation${outside === 1 ? '' : 's'} fell outside the null `
        + 'entirely.' : ''),
  }
}

/** C. Fold enrichment, so features can be ranked on one scale against a reference line at 1. */
export function figureFold(c: ComparisonResult): Figure {
  const rows = foldChart(c)
  const w = FIG_W
  const h = FIG.pad.top + rows.length * FIG.rowH + FIG.pad.bottom
  const x0 = FIG.pad.left
  const x1 = x0 + FIG.plotW
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
    marks.push({ k: 'text', x: x0 - 7, y: y + 2.6, s: fit(r.label, 32), size: FIG.font.label,
      colour: FIG.ink, anchor: 'end' })
    marks.push({ k: 'rect', x: x0, y: y - 3.4, w: Math.max(0.6, sx(r.fold) - x0), h: 6.8,
      colour: r.fold >= 1 ? FIG.ink : FIG.line })
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

/**
 * D. Which region touched which feature, transposed so features stay rows.
 *
 * TRANSPOSED FOR TWO REASONS. It puts feature names in the same left column as A, B and C, so a
 * reader tracks one feature straight down the composed figure instead of re-reading four different
 * label columns. And it avoids rotated text entirely, which the report's renderer cannot draw, so
 * the screen and the printed figure stay identical rather than one degrading.
 */
export function figureGrid(c: ComparisonResult, names: readonly string[]): Figure {
  const w = FIG_W
  const cols = Math.max(1, names.length)
  const colW = Math.min(16, FIG.plotW / cols)
  const gridW = colW * cols
  const h = FIG.pad.top + c.features.length * FIG.rowH + FIG.pad.bottom
  const x0 = FIG.pad.left
  const top = FIG.pad.top
  const marks: Mark[] = [panelLabel(4, 11, 'D')]

  c.features.forEach((f, r) => {
    const y = top + r * FIG.rowH + FIG.rowH / 2
    marks.push({ k: 'text', x: x0 - 7, y: y + 2.6, s: fit(f.label, 32), size: FIG.font.label,
      colour: FIG.ink, anchor: 'end' })
    names.forEach((_, i) => {
      const x = x0 + i * colW + colW / 2
      const hit = !!f.regionHits[i]
      marks.push({ k: 'dot', x, y, r: hit ? 3 : 1.1, colour: hit ? FIG.accent : FIG.line })
    })
    const n = f.regionHits.filter(Boolean).length
    marks.push({ k: 'text', x: x0 + gridW + 8, y: y + 2.6, s: `${n}/${cols}`,
      size: FIG.font.tick, colour: FIG.grey, anchor: 'start' })
  })

  // Region indices along the bottom, every fifth labelled so the axis stays legible at any count.
  const base = top + c.features.length * FIG.rowH
  marks.push({ k: 'line', x1: x0, y1: base, x2: x0 + gridW, y2: base,
    colour: FIG.axis, width: 0.8 })
  names.forEach((_, i) => {
    if (i % 5 !== 0 && i !== names.length - 1) return
    const x = x0 + i * colW + colW / 2
    marks.push({ k: 'line', x1: x, y1: base, x2: x, y2: base + FIG.tick.len,
      colour: FIG.axis, width: 0.8 })
    marks.push({ k: 'text', x, y: base + FIG.tick.len + FIG.tick.gap + FIG.font.tick,
      s: String(i + 1), size: FIG.font.tick, colour: FIG.grey, anchor: 'middle' })
  })
  marks.push({ k: 'text', x: x0 + gridW / 2,
    y: base + FIG.tick.len + FIG.tick.gap + FIG.font.tick + 10,
    s: `Region, in the order listed (n = ${cols})`, size: FIG.font.axisTitle,
    colour: FIG.ink, anchor: 'middle' })

  const touched = names.filter((_, i) => c.features.some((f) => f.regionHits[i])).length
  return {
    w,
    h,
    marks,
    caption: `D. Which region overlaps which feature; filled circles are overlaps, and the count at `
      + `the right is per feature, over every site class scanned. ${touched} of ${cols} region`
      + `${cols === 1 ? '' : 's'} touch at `
      + 'least one feature. An enrichment can rest on a handful of regions, and this is where that '
      + 'becomes visible rather than hiding inside a p value.',
  }
}

/**
 * The analytical panels, in order. These go to the report.
 *
 * Kept out of the screen panel deliberately: four quantitative figures are what a reader wants when
 * they are checking a result, and a map is what they want when they are meeting one. The screen
 * shows the map; the report carries both, since the map without the panels beneath it is an
 * illustration rather than evidence.
 */
export const comparisonFigures = (c: ComparisonResult, names: readonly string[]): Figure[] => [
  figureEnrichment(c), figureNulls(c), figureFold(c), figureGrid(c, names),
].filter((f) => f.marks.length > 1)

/**
 * The headline: every chromosome to scale, with fragile sites and the detected regions on it.
 *
 * THIS IS THE FIGURE THAT ANSWERS "WHERE", and the other four answer "how much". A reader looking
 * at a run wants to see the genome and what happened to it before they are asked to read a p value,
 * and a coincidence between a change and a fragile site is something the eye settles in a moment
 * and a table takes a paragraph to say.
 *
 * DRAWN FROM THE SHIPPED TRACK, NOT FROM HARDCODED LENGTHS. Chromosome extents and centromere
 * positions both come from the gap annotations already in the file, so the ideogram cannot drift out
 * of step with the features drawn on it: chr1's telomere gap ends at 249,250,621, which is chr1.
 *
 * ONE SHARED BASE-PAIR AXIS across all chromosomes, so position is comparable between rows. Drawing
 * each chromosome to its own width would make chr21 and chr1 look the same size, which is the same
 * mistake as a per-row axis on the enrichment panel.
 */
export function figureGenomeMap(c: ComparisonResult): Figure {
  const track = c.track
  const tracks = scannedTracks(track)
  // Order fixed rather than taken from object order, so the lane a colour means never moves.
  const laneOrder = Object.keys(tracks).sort()
  const ends = new Map<string, number>()
  const centro = new Map<string, { startBp: number; endBp: number }>()
  for (const g of track.gaps ?? []) {
    ends.set(g.chrom, Math.max(ends.get(g.chrom) ?? 0, g.endBp))
    if ((g as { kind?: string }).kind === 'centromere') {
      centro.set(g.chrom, { startBp: g.startBp, endBp: g.endBp })
    }
  }
  const order = [...ends.keys()].filter((k) => /^\d+$/.test(k)).sort((a, b) => Number(a) - Number(b))
  if (!order.length) return { w: 10, h: 10, marks: [], caption: '' }

  const genomeMax = Math.max(...order.map((k) => ends.get(k) ?? 0))
  const labelW = 30
  const plotW = 452
  const laneH = 2.6
  const barH = 5.5
  const gap = 5
  const blockH = barH + laneOrder.length * laneH + 4 + gap
  const legendRows = 2
  const legendH = legendRows * 11 + 6
  const w = labelW + plotW + 14
  const h = 18 + legendH + order.length * blockH + 30
  const x0 = labelW
  const sx = scale(0, genomeMax, x0, x0 + plotW)
  const top = 18 + legendH
  const marks: Mark[] = []

  // A key over two rows, because six entries on one line runs past the figure.
  const half = Math.ceil(laneOrder.length / 2)
  marks.push(...legend(x0, 18 + 8, laneOrder.slice(0, half).map((k) => ({
    kind: 'swatch' as const, colour: LANE_COLOUR[k] ?? FIG.line, s: shortLane(k),
  }))))
  marks.push(...legend(x0, 18 + 19, [
    ...laneOrder.slice(half).map((k) => ({
      kind: 'swatch' as const, colour: LANE_COLOUR[k] ?? FIG.line, s: shortLane(k),
    })),
    { kind: 'swatch' as const, colour: FIG.accent, s: 'detected region' },
  ]))

  const flags = regionFlags(c)
  order.forEach((chrom, i) => {
    const blockTop = top + i * blockH
    const end = ends.get(chrom) ?? 0
    marks.push({ k: 'text', x: x0 - 4, y: blockTop + barH, s: chrom, size: FIG.font.tick,
      colour: FIG.ink, anchor: 'end' })
    // The chromosome outline.
    marks.push({ k: 'rect', x: sx(0), y: blockTop, w: sx(end) - sx(0), h: barH, colour: '#eceff2' })
    const cen = centro.get(chrom)
    if (cen) {
      const cx = sx((cen.startBp + cen.endBp) / 2)
      marks.push({ k: 'dot', x: cx, y: blockTop + barH / 2, r: 1.7, colour: FIG.axis, hollow: true })
    }
    // DETECTED REGIONS ON THE CHROMOSOME ITSELF, not in a lane, because they are the thing the
    // reader came for and the feature lanes are context underneath them.
    flags.forEach((fl) => {
      if (fl.chrom !== chrom) return
      marks.push({ k: 'rect', x: sx(fl.startBp), y: blockTop,
        w: Math.max(1.3, sx(fl.endBp) - sx(fl.startBp)), h: barH, colour: FIG.accent })
      // A STAR WHERE A MUNDANE EXPLANATION IS AVAILABLE: the region sits in a site class whose
      // coincidence with this set of regions is more than the matched null produces. Not a verdict
      // on the region, a prior on how to read it.
      if (fl.related.length) {
        marks.push({ k: 'text', x: (sx(fl.startBp) + sx(fl.endBp)) / 2, y: blockTop - 1.2,
          s: '*', size: 8, colour: FIG.accent, anchor: 'middle', bold: true })
      }
    })
    // One lane per scanned site class, beneath, in a fixed order.
    laneOrder.forEach((name, li) => {
      const ly = blockTop + barH + 2 + li * laneH
      marks.push({ k: 'rect', x: sx(0), y: ly, w: sx(end) - sx(0), h: laneH - 0.7,
        colour: '#f4f6f8' })
      for (const f of tracks[name]) {
        if (f.chrom !== chrom) continue
        marks.push({ k: 'rect', x: sx(f.startBp), y: ly,
          w: Math.max(0.6, sx(f.endBp) - sx(f.startBp)), h: laneH - 0.7,
          colour: LANE_COLOUR[name] ?? FIG.line })
      }
    })
  })

  marks.push(...xAxis(x0, x0 + plotW, top + order.length * blockH - gap + 2, 0, genomeMax / 1e6,
    'Position (Mb)', (v) => String(Math.round(v))))

  const onChrom = new Set(c.regionSpans.map((r) => r.chrom)).size
  const anyHit = c.features.filter((f) => f.testable && f.significant).map((f) => f.label)
  return {
    w,
    h,
    marks,
    caption: `Every autosome to scale on one axis. The ${c.regionSpans.length} detected region`
      + `${c.regionSpans.length === 1 ? '' : 's'} are marked on the chromosome itself across `
      + `${onChrom} chromosome${onChrom === 1 ? '' : 's'}; beneath each is one lane per site class `
      + `this tool scans for (${laneOrder.length} of them), in the key order. Open circles are `
      + 'centromeres. '
      + (anyHit.length ? `Coincidence clears the corrected threshold for ${anyHit.join(' and ')}; `
        + `the ${flags.filter((f) => f.related.length).length} region`
        + `${flags.filter((f) => f.related.length).length === 1 ? '' : 's'} marked * sit in one of `
        + 'those classes and so have an alternative explanation available.'
        : 'No site class clears the corrected threshold, so no region is marked.')
      + ' Positional only: the fragile compartment is established on both parental genomes, so'
      + ' nothing here indicates which parent a change came from.',
  }
}

/**
 * One colour per site class, fixed so a colour never means two things.
 *
 * Chosen to differ in LIGHTNESS as well as hue, because a figure that is only readable in colour is
 * not readable in a printed paper.
 */
export const LANE_COLOUR: Record<string, string> = {
  'common fragile site': '#2f6f9f',
  'gene over 500 kb': '#7fa8c9',
  'centromere or telomere': '#9a9a9a',
  'late-replication valley (ES)': '#c2a25a',
  'late-replication valley (constitutive)': '#6b5a2e',
}

/** Short names, because a key with five full labels runs past the figure. */
const SHORT: Record<string, string> = {
  'common fragile site': 'fragile site',
  'gene over 500 kb': 'long gene',
  'centromere or telomere': 'centromere/telomere',
  'late-replication valley (ES)': 'late-repl (ES)',
  'late-replication valley (constitutive)': 'late-repl (const.)',
}
const shortLane = (k: string) => SHORT[k] ?? k

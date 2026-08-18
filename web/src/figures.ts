/**
 * Figure geometry, computed once and drawn by two renderers that cannot disagree.
 *
 * WHY THIS EXISTS. The screen and the report each had their own chart code with their own hardcoded
 * sizes: the same figure was 260 units wide in one and 210 in the other, rows 26 high against 18,
 * and neither drew an axis at all. Two implementations of one picture drift, and a reader who
 * notices the drift has to distrust both. So layout happens ONCE and emits a list of marks; each
 * renderer only knows how to draw a line, a rectangle, a dot and a string.
 *
 * WHAT MAKES A FIGURE READABLE, and the previous ones had none of it. An axis with tick marks at
 * round values, so a reader can recover a number from a position. A stated unit, so they know what
 * the number is. One typographic scale rather than whatever size fitted. Panel letters, so a
 * caption can refer to a part. A caption generated from the data, so it cannot describe a figure it
 * did not draw. And significance carried by a symbol as well as by colour, so the figure survives
 * being printed in grey.
 */

/** One drawing primitive. Deliberately tiny: anything a renderer cannot do is not in a figure. */
export type Mark =
  | { k: 'line'; x1: number; y1: number; x2: number; y2: number; colour: string; width: number
    dash?: boolean }
  | { k: 'rect'; x: number; y: number; w: number; h: number; colour: string }
  | { k: 'dot'; x: number; y: number; r: number; colour: string; hollow?: boolean }
  | { k: 'text'; x: number; y: number; s: string; size: number; colour: string
    anchor: 'start' | 'middle' | 'end'; bold?: boolean; italic?: boolean }

export interface Figure {
  w: number
  h: number
  marks: Mark[]
  /** Generated from the data, so it cannot describe a figure it did not draw. */
  caption: string
}

/**
 * One typographic and spatial scale for every figure.
 *
 * Sizes are in points and shared by both renderers, which is what makes a screen figure and a
 * printed one the same figure rather than two drawings of one idea.
 */
export const FIG = {
  font: { panel: 8.5, axisTitle: 7.5, label: 7, tick: 6.5, caption: 6.8 },
  ink: '#1a1a1a',
  grey: '#6b6b6b',
  line: '#cfcfcf',
  accent: '#a4243b',
  /** Space around the plot area. Left is wide because feature names live there. */
  pad: { top: 16, right: 58, bottom: 30, left: 150 },
  tick: { len: 3, gap: 2 },
  rowH: 15,
} as const

/**
 * Round tick values covering a range, which is the single thing most missing from the old charts.
 *
 * A reader recovers a number from a position by counting ticks, so the ticks must sit at values a
 * person would choose: 0, 0.25, 0.5 rather than 0, 0.31, 0.62. The 1-2-5 progression is the
 * standard choice and is what every plotting library uses, reproduced here because pulling in a
 * charting dependency for eleven lines of arithmetic would be the larger cost.
 */
export function niceTicks(min: number, max: number, target = 5): number[] {
  if (!(max > min)) return [min]
  const raw = (max - min) / Math.max(1, target)
  const mag = 10 ** Math.floor(Math.log10(raw))
  const norm = raw / mag
  const step = (norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 5 ? 5 : 10) * mag
  const out: number[] = []
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) {
    // Snap away the float dust a repeated addition accumulates, or a tick reads 0.30000000000000004.
    out.push(Math.abs(v) < step * 1e-9 ? 0 : Number(v.toFixed(10)))
  }
  return out
}

/** A linear map from data units to figure units. */
export const scale = (d0: number, d1: number, r0: number, r1: number) =>
  (v: number) => (d1 === d0 ? r0 : r0 + ((v - d0) / (d1 - d0)) * (r1 - r0))

/**
 * A horizontal axis with outward ticks, tick labels and a title.
 *
 * Returned as marks rather than drawn, so the same axis appears identically in both renderers and
 * so its geometry can be checked without rendering anything.
 */
export function xAxis(
  x0: number, x1: number, y: number, d0: number, d1: number, title: string,
  fmt: (v: number) => string, target = 5,
): Mark[] {
  const sx = scale(d0, d1, x0, x1)
  const ticks = niceTicks(d0, d1, target)
  const marks: Mark[] = [
    { k: 'line', x1: x0, y1: y, x2: x1, y2: y, colour: FIG.line, width: 0.6 },
  ]
  for (const t of ticks) {
    const x = sx(t)
    marks.push({ k: 'line', x1: x, y1: y, x2: x, y2: y + FIG.tick.len, colour: FIG.line, width: 0.6 })
    marks.push({
      k: 'text', x, y: y + FIG.tick.len + FIG.tick.gap + FIG.font.tick, s: fmt(t),
      size: FIG.font.tick, colour: FIG.grey, anchor: 'middle',
    })
  }
  marks.push({
    k: 'text', x: (x0 + x1) / 2, y: y + FIG.tick.len + FIG.tick.gap + FIG.font.tick + 10,
    s: title, size: FIG.font.axisTitle, colour: FIG.ink, anchor: 'middle',
  })
  return marks
}

/** A panel letter, so a caption can refer to a part of a figure. */
export const panelLabel = (x: number, y: number, letter: string): Mark => ({
  k: 'text', x, y, s: letter, size: FIG.font.panel, colour: FIG.ink, anchor: 'start', bold: true,
})

/** Truncate a label to fit the space reserved for it, with an ellipsis so the cut is visible. */
export const fit = (s: string, chars: number): string =>
  (s.length <= chars ? s : `${s.slice(0, Math.max(1, chars - 1))}…`)

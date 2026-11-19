/**
 * The cell figure: how many copies each parent gave, at each place in the genome.
 *
 * One figure at two zoom levels. The top deck is the whole cell, every chromosome, one track per
 * parent where both are measured and the copies in all otherwise, with the chromosomes that carry no
 * change muted so the eye lands on the ones that do. Under it, one zoom per change in chromosome
 * order, each the same figure along that one chromosome with the called stretch marked, numbered to
 * match the badge above. The same zoom is what a single call's expander shows.
 */
import { useLayoutEffect, useRef, useState, type PointerEvent, type ReactNode } from 'react'
import {
  cellPer, copyBars, COPY_MAX, copyWord, guardedBars, mergeBins, normalRange, perForStretch, tickStepMb,
  totalBars, viewOf, wrap, zoomPicks, type Bin, type PlotData, type Role, type Track,
} from './plotSeries.ts'
import { headline, KIND_WORD, type Defect } from './defects.ts'
import { saidOf } from './pageSays.ts'
import { segmentCoords, type Segment } from './segments.ts'
import { OTHER_MIN_SPAN_BP } from './otherParentLoss.ts'
import type { Refusal } from './refusals.ts'

const FATHER = 'var(--om-blue)'
const MOTHER = 'var(--om-higher)'
const NEUTRAL = '#6b727b'
const CALL = 'var(--om-defect)'
const GAINC = 'var(--om-gain)'
const DIM = 'var(--om-text-dim)'
const INK = '#1d2126'
const RANGE = '#e9ecef'
const LINE = '#e3e6ea'
const colourOf = (p: Role | null | undefined) => (p === 'paternal' ? FATHER : p === 'maternal' ? MOTHER : NEUTRAL)
/** A loss, a gain, or a reading the card does not count: one rule, used by every part of the figure. */
const markColour = (m: { gain: boolean; withheld?: string; status?: string }) =>
  (m.withheld || (m.status !== undefined && m.status !== 'counted') ? DIM : m.gain ? GAINC : CALL)
const AUTOSOMES = Array.from({ length: 22 }, (_, i) => String(i + 1))
const mb = (bp: number) => (bp / 1e6).toFixed(1)

/** One entry of the list under the deck: the badge it carries under its chromosome, then its words. */
export interface Chip {
  n: number
  chrom: string
  /** Whole chromosome, or the stretch in Mb. */
  where: string
  /** What happened and whose copy, in the figure's own short words. */
  what: string
  /** Why the card does not count this reading, where it does not. */
  note?: string
  zoomed?: boolean
  colour: string
  bg: string
  muted?: boolean
}

/** The width a chip needs before the list is worth setting in two columns. */
const CHIP_MIN_W = 430

/** One call as the figures mark it. */
export interface Mark {
  chrom: string
  startBp: number
  endBp: number
  /** Edge intervals reaching past the call's edges. */
  outLoBp?: number
  outHiBp?: number
  whole: boolean
  gain: boolean
  parent: Role | null
  withheld?: string
  /** What the call is, in one plain sentence. */
  sentence: string
  /** The same finding in two or three words, for the list under the deck and the zoom's title. */
  word: string
  /** What the card's counting rule made of it: a counted change, or a reading it does not count. */
  status?: 'counted' | 'one reading' | 'also seen'
  /** What the row shows, for the hover. */
  lines: string[]
  /** The call states its own numbers beside the figure, so the caption quotes none of its own. */
  quoted?: boolean
  /** Intensity across the call as its scan measured it, quoted in place of the figure's own. */
  lrr?: number
  /** The call's own numbers as a plain sentence, for the line under its bars. */
  said?: string
}

function useWidth(): [React.RefObject<HTMLDivElement | null>, number] {
  const ref = useRef<HTMLDivElement>(null)
  const [w, setW] = useState(720)
  useLayoutEffect(() => {
    const el = ref.current
    if (!el) return
    const ro = new ResizeObserver(() => setW(Math.max(300, Math.floor(el.clientWidth))))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])
  return [ref, w]
}

/** The page's colour variables written out, on a white ground, so the file stands on its own. */
export function svgText(svg: SVGSVGElement): string {
  const css = getComputedStyle(document.documentElement)
  return svg.outerHTML.replace(/var\((--[a-z0-9-]+)\)/g, (_, v: string) => css.getPropertyValue(v).trim() || '#888')
}

function SaveSvg({ svg, name }: { svg: React.RefObject<SVGSVGElement | null>; name: string }) {
  return (
    <button
      type="button"
      onClick={() => {
        if (!svg.current) return
        const a = document.createElement('a')
        a.href = URL.createObjectURL(new Blob([svgText(svg.current)], { type: 'image/svg+xml' }))
        a.download = `${name}.svg`
        a.click()
        setTimeout(() => URL.revokeObjectURL(a.href), 1000)
      }}
      style={{
        border: 0, background: 'none', font: 'inherit', fontSize: 11, padding: 0, cursor: 'pointer',
        marginLeft: 'auto', paddingLeft: 12, whiteSpace: 'nowrap', color: 'var(--om-blue)', textDecoration: 'underline',
      }}
    >
      Save SVG
    </button>
  )
}

function Tip({ tip, w }: { tip: { x: number; y: number; lines: string[] } | null; w: number }) {
  if (!tip) return null
  return (
    <div style={{
      position: 'absolute', left: Math.max(0, Math.min(tip.x + 12, w - 300)), top: tip.y + 14, width: 'max-content',
      maxWidth: Math.min(380, w), pointerEvents: 'none', zIndex: 2, background: '#fff',
      border: '1px solid var(--om-border-strong)', padding: '4px 7px', fontSize: 12, lineHeight: 1.45,
    }}
    >
      {tip.lines.map((l, i) => <div key={i} style={{ fontWeight: i === 0 ? 700 : 400 }}>{l}</div>)}
    </div>
  )
}

/* --- the figure ------------------------------------------------------------------------------ */

/** The drawn tracks: both parents where both are measured, the copies in all otherwise. */
const drawnTracks = (data: PlotData): Track[] => data.tracks.filter((t) => t.shown)

/** One track's bars along a chromosome, already in copies. */
function barsOf(data: PlotData, t: Track, chrom: string, per: number): Bin[] {
  const own = t.bins.get(chrom) ?? []
  if (t.kind === 'level') return totalBars(mergeBins(own, per, true))
  const level = data.tracks.find((x) => x.key === 'intensity')?.bins.get(chrom) ?? []
  return copyBars(mergeBins(own, per), level)
}

/** The range this cell itself holds on the chromosomes with no change, in copies. */
function copyRange(data: PlotData, t: Track, per: number, skip: ReadonlySet<string>) {
  const v: number[] = []
  for (const chrom of t.bins.keys()) {
    if (skip.has(chrom)) continue
    for (const b of barsOf(data, t, chrom, per)) if (Number.isFinite(b.k)) v.push(b.k)
  }
  v.sort((a, b) => a - b)
  const at = (q: number) => (v.length ? v[Math.min(v.length - 1, Math.floor(q * v.length))] : NaN)
  return { lo: at(0.1), mid: at(0.5), hi: at(0.9) }
}

/** The call rate this cell holds on the chromosomes with no change: the guard's floor. */
function rateFloor(data: PlotData, per: number, skip: ReadonlySet<string>): number {
  const rate = data.tracks.find((t) => t.key === 'call-rate')
  return rate ? normalRange(rate, per, skip).lo : NaN
}

/** Bars both parents read as none on this chromosome while the markers are still calling. */
function guardOf(data: PlotData, chrom: string, per: number, floor: number): Set<number> {
  const f = data.tracks.find((t) => t.key === 'father')
  const m = data.tracks.find((t) => t.key === 'mother')
  const rate = data.tracks.find((t) => t.key === 'call-rate')
  if (!f || !m || !rate) return new Set()
  return guardedBars(barsOf(data, f, chrom, per), barsOf(data, m, chrom, per), rate.bins.get(chrom) ?? [], floor)
}

const copyLine = (t: Track, chrom: string, b: Bin, guarded: boolean): string[] => [
  `chr${chrom} ${mb(b.startBp)} to ${mb(b.endBp)} Mb`,
  guarded ? 'Both parents read none here while the markers are still calling, so no level is drawn.'
    : `${b.n} markers, ${copyWord(b.k)} (${b.k.toFixed(2)})`,
  t.tech,
]

/** The bars of one track, drawn from none upward; a guarded bar is an empty outline. */
function Bars({ bars, guard, x, y, colour, muted, inside }: {
  bars: readonly Bin[]; guard: ReadonlySet<number>; x: (bp: number) => number; y: (v: number) => number
  colour: string; muted?: boolean; inside?: (b: Bin) => boolean
}) {
  return (
    <g opacity={muted ? 0.3 : 1}>
      {bars.map((b, i) => {
        if (!Number.isFinite(b.k)) return null
        const x0 = x(b.startBp); const x1 = Math.max(x0 + 0.8, x(b.endBp))
        if (guard.has(i)) {
          return (
            <rect key={i} x={x0 + 0.5} y={y(1)} width={Math.max(1, x1 - x0 - 1)} height={Math.max(1, y(0) - y(1))}
              fill="none" stroke={colour} strokeWidth={1} strokeDasharray="2,2" opacity={0.7} />
          )
        }
        const top = y(Math.min(b.k, COPY_MAX))
        return (
          <rect key={i} x={x0 + 0.4} y={top} width={Math.max(0.8, x1 - x0 - 0.8)} height={Math.max(0.8, y(0) - top)}
            fill={colour} opacity={inside && !inside(b) ? 0.45 : 1} />
        )
      })}
    </g>
  )
}

/** The level lines of a track: none, one copy, two copies. */
function Levels({ x0, w, y, words, dim }: {
  x0: number; w: number; y: (v: number) => number; words: boolean; dim: boolean
}) {
  return (
    <g>
      {[0, 1, 2].map((v) => (
        <g key={v}>
          <line x1={x0} x2={x0 + w} y1={y(v)} y2={y(v)} stroke={v === 1 ? (dim ? LINE : DIM) : LINE} />
          {words ? <text x={x0 - 6} y={y(v) + 4} textAnchor="end" fontSize={10} fill={DIM}>{LEVEL_WORD[v]}</text> : null}
        </g>
      ))}
    </g>
  )
}

const LEVEL_WORD: Record<number, string> = { 0: 'none', 1: 'one copy', 2: 'two copies' }

/* --- the whole cell, the top deck -------------------------------------------------------------- */

export function WholeCell({ data, marks, name, numbers, focus, title, chips, caption }: {
  data: PlotData; marks: readonly Mark[]; name: string
  /** Each mark's number in the list below, keyed `chrom:startBp`. */
  numbers?: ReadonlyMap<string, number>
  /** A chromosome to bring forward on its own, for a single call's view. */
  focus?: string
  /** What the figure shows, every change it marks, and what the reader must know to read it. They
   *  are drawn inside the picture: the card offers the file for saving, and a figure that travels
   *  without its caption claims more than it measured. */
  title?: string
  /** One chip per mark, already numbered, in chromosome order. */
  chips?: readonly Chip[]
  caption?: string
}) {
  const [ref, W] = useWidth()
  const svg = useRef<SVGSVGElement>(null)
  const [tip, setTip] = useState<{ x: number; y: number; lines: string[] } | null>(null)
  const tracks = drawnTracks(data)
  if (!tracks.length) return null
  const called = new Set(marks.filter((m) => !m.withheld).map((m) => m.chrom))
  const changed = new Set(marks.map((m) => m.chrom))
  const pers = tracks.map(cellPer)
  const ranges = tracks.map((t, i) => copyRange(data, t, pers[i], called))
  const floor = rateFloor(data, pers[0], called)
  const GUT = 70
  const rows = W >= 760 ? [AUTOSOMES] : [AUTOSOMES.slice(0, 11), AUTOSOMES.slice(11)]
  const G = 3
  const TH = 74
  const len = (c: string) => data.endBp.get(c) ?? 1
  type Panel = { c: string; x: number; w: number; top: number; row: number }
  const panels: Panel[] = []
  // Below the tracks: the chromosome number, the mark, and the mark's number in the list. One lane,
  // so a reader never has to work out which chromosome a label belongs to.
  const blockOf = (row: readonly string[]) => tracks.length * (TH + 20) + 21
    + (row.some((c) => changed.has(c)) ? 24 : 6)
  const headLines = title ? wrap(title, Math.max(24, Math.floor(W / 7.6))) : []
  const headH = headLines.length ? 8 + headLines.length * 19 + 6 : 0
  // THE LIST READS AS CHIPS, NOT AS A PARAGRAPH. Twelve lines of one weight read as prose and the
  // eye cannot land on one; each entry gets the badge it carries under its chromosome and its own
  // line. Two columns where the width carries the longest chip, one where it does not.
  const chipList = chips ?? []
  const CHIP_H = 19
  const chipCols = W >= 2 * CHIP_MIN_W + 18 ? 2 : 1
  const chipW = (W - 18 * (chipCols - 1)) / chipCols
  const chipRows = Math.ceil(chipList.length / chipCols)
  const capLines = caption ? wrap(caption, Math.max(30, Math.floor(W / 5.6))) : []
  const footH = (chipList.length ? chipRows * CHIP_H + 14 : 0) + (capLines.length ? capLines.length * 14 + 10 : 0)
  const bodyTop = 6 + headH
  rows.forEach((row, r) => {
    const total = row.reduce((s, c) => s + len(c), 0)
    const usable = W - 2 - GUT - G * (row.length - 1)
    let x = 1 + GUT
    const top = bodyTop + rows.slice(0, r).reduce((s, q) => s + blockOf(q), 0)
    for (const c of row) {
      const w = usable * len(c) / total
      panels.push({ c, x, w, top, row: r })
      x += w + G
    }
  })
  const H = bodyTop + rows.reduce((s, row) => s + blockOf(row), 0) + 4 + footH
  const yOf = (top: number, v: number) => top + TH * (1 - Math.max(0, Math.min(COPY_MAX, v)) / COPY_MAX)
  const onMove = (e: PointerEvent<SVGSVGElement>) => {
    const r = e.currentTarget.getBoundingClientRect()
    const px = e.clientX - r.left; const py = e.clientY - r.top
    const p = panels.find((q) => px >= q.x && px <= q.x + q.w && py >= q.top && py <= q.top + tracks.length * (TH + 20))
    if (!p) { setTip(null); return }
    const ti = Math.floor((py - p.top) / (TH + 20))
    const t = tracks[ti]
    if (!t) { setTip(null); return }
    const bp = (px - p.x) / p.w * len(p.c)
    const m = marks.find((q) => q.chrom === p.c && bp >= q.startBp && bp <= q.endBp)
    if (m) { setTip({ x: px, y: py, lines: [m.sentence, ...m.lines.slice(0, 3)] }); return }
    const bars = barsOf(data, t, p.c, pers[ti])
    const i = bars.findIndex((q) => bp >= q.startBp && bp <= q.endBp)
    if (i < 0) { setTip(null); return }
    setTip({ x: px, y: py, lines: copyLine(t, p.c, bars[i], guardOf(data, p.c, pers[ti], floor).has(i)) })
  }
  return (
    <div ref={ref} style={{ position: 'relative', margin: '8px 0 2px', maxWidth: '100%' }}>
      <svg ref={svg} xmlns="http://www.w3.org/2000/svg" width={W} height={H} viewBox={`0 0 ${W} ${H}`}
        fontFamily="system-ui, -apple-system, sans-serif" fontSize={11} role="img"
        aria-label="Every chromosome of this cell, how many copies came from each parent"
        style={{ display: 'block', touchAction: 'pan-y' }}
        onPointerMove={onMove} onPointerDown={onMove} onPointerLeave={() => setTip(null)}
      >
        <rect width={W} height={H} fill="#fff" />
        {headLines.map((l, i) => (
          <text key={`h${i}`} x={1} y={20 + i * 19} fontSize={15} fontWeight={700} fill={INK}>{l}</text>
        ))}
        {rows.map((_row, r) => tracks.map((t, ti) => {
          const top = bodyTop + rows.slice(0, r).reduce((s, q) => s + blockOf(q), 0) + ti * (TH + 20) + 16
          // One cell has two tracks however many rows of chromosomes it takes, so the title is
          // written once and a later row carries a quiet continuation mark instead.
          return r === 0 ? (
            <text key={`l${r}${ti}`} x={1} y={top - 5} fontWeight={700} fill={colourOf(t.parent)}>{t.label}</text>
          ) : (
            <text key={`l${r}${ti}`} x={1} y={top - 5} fontSize={10} fill={DIM}>{t.label}, continued</text>
          )
        }))}
        {panels.map((p) => {
          const on = changed.has(p.c) && (!focus || p.c === focus)
          const guards = tracks.map((t, ti) => (t.kind === 'share' && ti === 0 ? guardOf(data, p.c, pers[ti], floor) : new Set<number>()))
          return (
            <g key={p.c}>
              {tracks.map((t, ti) => {
                const top = p.top + ti * (TH + 20) + 16
                const x = (bp: number) => p.x + p.w * bp / len(p.c)
                const rg = ranges[ti]
                const y = (v: number) => yOf(top, v)
                return (
                  <g key={t.key}>
                    <rect x={p.x} y={top} width={p.w} height={TH} fill="#fafbfc" />
                    {Number.isFinite(rg.lo) ? (
                      <rect x={p.x} y={y(rg.hi)} width={p.w} height={Math.max(1, y(rg.lo) - y(rg.hi))}
                        fill={RANGE} opacity={on ? 0.9 : 0.4} />
                    ) : null}
                    <Levels x0={p.x} w={p.w} y={y} words={p.c === AUTOSOMES[0] || (p.row > 0 && p.c === '12')} dim={!on} />
                    <Bars bars={barsOf(data, t, p.c, pers[ti])} guard={guards[ti]} x={x} y={y}
                      colour={colourOf(t.parent)} muted={!on} />
                  </g>
                )
              })}
              {(() => {
                const axis = p.top + tracks.length * (TH + 20) + 4
                const ms = marks.filter((m) => m.chrom === p.c)
                  .sort((a, b) => (numbers?.get(`${a.chrom}:${a.startBp}`) ?? 0) - (numbers?.get(`${b.chrom}:${b.startBp}`) ?? 0))
                // A chromosome carrying a loss and a gain says so: each mark draws its own piece of
                // the line in its own colour, and the number badge beneath it takes the same colour.
                const seg = ms.length ? p.w / ms.length : p.w
                return (
                  <g>
                    <text x={p.x + p.w / 2} y={axis + 12} textAnchor="middle" fontWeight={ms.length ? 700 : 400}
                      fill={ms.length ? markColour(ms[0]) : '#aeb4bb'}>{p.c}</text>
                    {ms.map((m, i) => (
                      <rect key={`b${i}`} x={p.x + i * seg + (i ? 1 : 0)} y={axis + 17}
                        width={Math.max(2, seg - (ms.length > 1 ? 2 : 0))} height={3} fill={markColour(m)} />
                    ))}
                    {ms.map((m, i) => {
                      const n = numbers?.get(`${m.chrom}:${m.startBp}`)
                      if (n === undefined) return null
                      const cx = p.x + p.w / 2 + (i - (ms.length - 1) / 2) * 17
                      return (
                        <g key={n}>
                          <circle cx={cx} cy={axis + 32} r={7.5} fill={markColour(m)} />
                          <text x={cx} y={axis + 36} textAnchor="middle" fontSize={10} fontWeight={700}
                            fill="#fff">{n}</text>
                        </g>
                      )
                    })}
                  </g>
                )
              })()}
            </g>
          )
        })}
        {chipList.map((c, i) => {
          const col = i % chipCols
          const row = Math.floor(i / chipCols)
          const x0 = col * (chipW + 18)
          const y0 = bodyTop + rows.reduce((s, q) => s + blockOf(q), 0) + 12 + row * CHIP_H
          return (
            <g key={`chip${c.n}`}>
              <rect x={x0} y={y0 - 1} width={chipW} height={CHIP_H - 3} rx={3}
                fill={c.muted ? '#f6f7f8' : c.bg} />
              <circle cx={x0 + 11} cy={y0 + 7} r={7.5} fill={c.colour} />
              <text x={x0 + 11} y={y0 + 11} textAnchor="middle" fontSize={10} fontWeight={700} fill="#fff">{c.n}</text>
              {/* One text run, so the chromosome and its words cannot overlap at any width. */}
              <text x={x0 + 24} y={y0 + 11} fontSize={11} fill={c.muted ? DIM : INK}>
                <tspan fontWeight={700}>{`chr${c.chrom}`}</tspan>
                {` ${c.where}: ${c.what}${c.note ? `, ${c.note}` : ''}`}
              </text>
              {c.zoomed ? (
                <text x={x0 + chipW - 4} y={y0 + 11} textAnchor="end" fontSize={9} fill={DIM}>shown below</text>
              ) : null}
            </g>
          )
        })}
        {capLines.map((l, i) => (
          <text key={`c${i}`} x={1}
            y={bodyTop + rows.reduce((s, q) => s + blockOf(q), 0) + 12 + (chipList.length ? chipRows * CHIP_H + 12 : 0)
              + i * 14}
            fontSize={11} fill={DIM}>{l}</text>
        ))}
      </svg>
      <Tip tip={tip} w={W} />
      <div style={{ display: 'flex', fontSize: 11, color: DIM, marginTop: 2 }}>
        <SaveSvg svg={svg} name={`${name}-cell`} />
      </div>
    </div>
  )
}

/* --- one change, the zoom ---------------------------------------------------------------------- */

export function Zoom({ data, mark, marks, number, width, title, caption, saveName }: {
  data: PlotData; mark: Mark; marks: readonly Mark[]; number?: number; width?: number
  /** The zoom's own words, drawn inside it so the saved file says what it shows and what it does not. */
  title?: string
  caption?: string
  /** Offer this zoom as a file. The button sits beside the picture, never inside it. */
  saveName?: string
}) {
  const [ref, auto] = useWidth()
  const svg = useRef<SVGSVGElement>(null)
  const [tip, setTip] = useState<{ x: number; y: number; lines: string[] } | null>(null)
  const W = width ?? auto
  const tracks = drawnTracks(data)
  const called = new Set(marks.filter((m) => !m.withheld).map((m) => m.chrom))
  const chromEnd = Math.max(data.endBp.get(mark.chrom) ?? 0, mark.outHiBp ?? mark.endBp)
  const { loBp, hiBp } = mark.whole ? { loBp: 0, hiBp: chromEnd } : viewOf(mark.startBp, mark.endBp, chromEnd)
  const first = tracks[0]
  const per = mark.whole ? cellPer(first) : perForStretch(first.bins.get(mark.chrom) ?? [], mark.startBp, mark.endBp)
  const floor = rateFloor(data, per, called)
  const guard = guardOf(data, mark.chrom, per, floor)
  const L = 66; const R = 8; const PW = Math.max(40, W - L - R)
  const x = (bp: number) => L + PW * Math.max(0, Math.min(1, (bp - loBp) / (hiBp - loBp)))
  const head = wrap(title ?? `Chromosome ${mark.chrom}: ${mark.word}`, Math.max(22, Math.floor(W / 7)))
  const capLines = caption ? wrap(caption, Math.max(28, Math.floor(W / 5.3))) : []
  const TH = W < 420 ? 66 : 78
  const top0 = (number ? 6 : 0) + head.length * 15 + 12
  const colour = markColour(mark)
  const H = top0 + tracks.length * (TH + 20) + 44 + (capLines.length ? capLines.length * 13 + 6 : 0)
  const span = hiBp - loBp
  const step = tickStepMb(span)
  const ticks: number[] = []
  for (let v = Math.ceil(loBp / 1e6 / step) * step; v * 1e6 <= hiBp; v += step) ticks.push(v)
  const axisY = top0 + (tracks.length - 1) * (TH + 20) + TH
  const onMove = (e: PointerEvent<SVGSVGElement>) => {
    const r = e.currentTarget.getBoundingClientRect()
    const px = e.clientX - r.left; const py = e.clientY - r.top
    const ti = Math.floor((py - top0) / (TH + 20))
    const t = tracks[ti]
    if (!t || px < L || px > L + PW) { setTip(null); return }
    const bp = loBp + (px - L) / PW * span
    const bars = barsOf(data, t, mark.chrom, per)
    const i = bars.findIndex((q) => bp >= q.startBp && bp <= q.endBp)
    if (i < 0) { setTip(null); return }
    setTip({ x: px, y: py, lines: copyLine(t, mark.chrom, bars[i], ti === 0 && guard.has(i)) })
  }
  return (
    <div ref={ref} style={{ position: 'relative', maxWidth: '100%' }}>
      <svg ref={svg} xmlns="http://www.w3.org/2000/svg" width={W} height={H} viewBox={`0 0 ${W} ${H}`}
        fontFamily="system-ui, -apple-system, sans-serif" fontSize={11} role="img" aria-label={mark.sentence}
        style={{ display: 'block', touchAction: 'pan-y' }}
        onPointerMove={onMove} onPointerDown={onMove} onPointerLeave={() => setTip(null)}
      >
        <rect width={W} height={H} fill="#fff" />
        {number ? (
          <g>
            <circle cx={9} cy={9} r={8} fill={colour} />
            <text x={9} y={13} textAnchor="middle" fontSize={10} fontWeight={700} fill="#fff">{number}</text>
          </g>
        ) : null}
        {head.map((l, i) => (
          <text key={i} x={number ? 23 : 1} y={13 + i * 15} fontSize={12} fontWeight={700} fill={INK}>{l}</text>
        ))}
        {tracks.map((t, ti) => {
          const top = top0 + ti * (TH + 20)
          const y = (v: number) => top + TH * (1 - Math.max(0, Math.min(COPY_MAX, v)) / COPY_MAX)
          const rg = copyRange(data, t, per, called)
          const bars = barsOf(data, t, mark.chrom, per).filter((b) => b.endBp >= loBp && b.startBp <= hiBp)
          return (
            <g key={t.key}>
              <text x={L} y={top - 6} fontWeight={700} fill={colourOf(t.parent)}>{t.label}</text>
              {Number.isFinite(rg.lo) ? (
                <rect x={L} y={y(rg.hi)} width={PW} height={Math.max(1, y(rg.lo) - y(rg.hi))} fill={RANGE} />
              ) : null}
              <rect x={x(mark.startBp)} y={top} width={Math.max(3, x(mark.endBp) - x(mark.startBp))} height={TH}
                fill={colour} opacity={0.07} />
              <Levels x0={L} w={PW} y={y} words dim={false} />
              <Bars bars={bars} guard={ti === 0 ? guard : new Set()} x={x} y={y} colour={colourOf(t.parent)}
                inside={(b) => b.startBp >= mark.startBp && b.endBp <= mark.endBp} />
              {[mark.startBp, mark.endBp].map((v, i) => (
                <line key={i} x1={x(v)} x2={x(v)} y1={top} y2={top + TH} stroke={colour} strokeWidth={1}
                  strokeDasharray="4,3" />
              ))}
              <line x1={L} x2={L + PW} y1={top + TH} y2={top + TH} stroke={DIM} />
            </g>
          )
        })}
        {ticks.map((v) => (
          <text key={v} x={x(v * 1e6)} y={axisY + 14} textAnchor="middle" fontSize={10} fill={DIM}>{v}</text>
        ))}
        <text x={L - 6} y={axisY + 14} textAnchor="end" fontSize={10} fill={DIM}>Mb</text>
        <rect x={x(mark.startBp)} y={axisY + 22} width={Math.max(3, x(mark.endBp) - x(mark.startBp))} height={3} fill={colour} />
        <text x={x(mark.startBp)} y={axisY + 38} fontSize={10} fontWeight={600} fill={colour}>
          {mark.withheld ? `found, not reported: ${mark.withheld}` : 'called here'}
        </text>
        {capLines.map((l, i) => (
          <text key={`c${i}`} x={1} y={axisY + 56 + i * 13} fontSize={10} fill={DIM}>{l}</text>
        ))}
      </svg>
      <Tip tip={tip} w={W} />
      {saveName ? (
        <div style={{ display: 'flex', fontSize: 11, color: DIM, marginTop: 2 }}>
          <SaveSvg svg={svg} name={saveName} />
        </div>
      ) : null}
    </div>
  )
}

/** One change on its own, for the expanders beside the change list and the segment table. */
export function CallFigure({ mark, data, marks, name }: {
  mark: Mark; data: PlotData | undefined; marks: readonly Mark[]; name: string
}) {
  if (!data) return null
  return (
    <div style={{ margin: '8px 0 2px', maxWidth: '100%' }}>
      <Zoom data={data} mark={mark} marks={marks} title={mark.sentence}
        caption={`${mark.said ? `${mark.said} ` : ''}${zoomCaption(data, mark)}`}
        saveName={`${name}-chr${mark.chrom}`} />
    </div>
  )
}

const CAPTION_SCALE = 'Each bar is a stretch of DNA and its height is how many copies of it came from that '
  + 'parent. One from each is a normal chromosome, and the grey band is the range this cell itself holds '
  + 'on the chromosomes with no change.'
const CAPTION_GUARD = 'Where both parents read none while the markers are still calling at this cell\'s own '
  + 'rate, the bar is an empty outline rather than a confident zero: a stretch with no DNA left still '
  + 'returns calls, so the call rate is the guard on a level of none.'
const CAPTION_LEVELS = 'None against one copy is read from the alleles, which need no reference level; one '
  + 'copy against two is read from the amount of DNA, which needs one.'
const CAPTION_GAIN = 'Read the doubled track: in a single cell a doubled copy takes more than its share of '
  + 'what is read, so the other parent\'s track dips below one copy although that copy is present. A '
  + 'surviving copy beside a lost one reads a little high for the same reason.'
const CAPTION_SOLO = 'With one parent or none the figure draws the copies in all and names nobody: the '
  + 'markers that say whose DNA is here are the ones where the parents carry different alleles, and they '
  + 'cannot be built from one parent.'

/** What a zoom must carry on its own: the scale, the band, and the caveat that bears on it. */
function zoomCaption(data: PlotData, mark: Mark): string {
  const scale = data.split
    ? 'Each bar is a stretch of DNA and its height is how many copies of it came from that parent. '
      + "One from each is a normal chromosome; the grey band is this cell's own range on the "
      + 'chromosomes with no change.'
    : 'Each bar is a stretch of DNA and its height is how many copies are here in all, against the '
      + "level this cell holds on its chromosomes with no change, which is the grey band. No parent "
      + 'is named.'
  const gain = data.split && mark.gain
    ? ' Read the doubled track: in a single cell a doubled copy takes more than its share of what is '
      + "read, so the other parent's track dips below one copy although that copy is present." : ''
  const guard = data.split
    ? ' A bar drawn as an empty outline is a stretch both parents read as none while the markers are '
      + 'still calling, which the call rate does not support as a zero.' : ''
  return scale + gain + guard
}

/** The figure for a result card: the whole cell, then a zoom for each change, numbered to match. */
export function CellFigure({ data, marks, name, title }: {
  data: PlotData | undefined; marks: readonly Mark[]; name: string
  /** What the figure shows, in one sentence: pageSays.figureTitle. */
  title: string
}): ReactNode {
  const [ref, W] = useWidth()
  if (!data) return null
  const { shown } = zoomPicks(marks.filter((m) => !m.withheld))
  // EVERY mark is numbered, in chromosome order, and the list below holds every number: a number
  // under its own chromosome and a line that names it cannot be read onto the wrong chromosome,
  // which three lanes of floating labels could. Whether a change is zoomed is a quiet note on its
  // line, never the presence or absence of a badge.
  const ordered = [...marks].sort((a, b) => +a.chrom - +b.chrom || a.startBp - b.startBp)
  const numbers = new Map<string, number>()
  ordered.forEach((m, i) => numbers.set(`${m.chrom}:${m.startBp}`, i + 1))
  const zoomed = new Set(shown.map((m) => `${m.chrom}:${m.startBp}`))
  // The card's own three kinds, in the card's own words: a counted change at full strength, and a
  // reading it does not count set quietly, so the two are told apart at a glance without a legend.
  const NOTE: Record<string, string> = {
    'one reading': 'seen once, not counted as a change',
    'also seen': 'read but not counted as a change',
  }
  const chips: Chip[] = ordered.map((m) => {
    const key = `${m.chrom}:${m.startBp}`
    const quiet = !!m.withheld || (m.status !== undefined && m.status !== 'counted')
    return {
      n: numbers.get(key)!,
      chrom: m.chrom,
      where: m.whole ? 'whole chromosome' : `${mb(m.startBp)} to ${mb(m.endBp)} Mb`,
      what: m.word,
      note: m.withheld ? 'found but not reported' : NOTE[m.status ?? 'counted'],
      zoomed: zoomed.has(key),
      colour: markColour(m),
      bg: quiet ? '#f6f7f8' : m.gain ? 'var(--om-gain-bg)' : 'var(--om-defect-bg)',
      muted: quiet,
    }
  })
  const cols = W >= 980 ? 3 : W >= 640 ? 2 : 1
  const colW = Math.floor((W - 18 * (cols - 1)) / cols)
  const hasGain = marks.some((m) => m.gain)
  // One source of words, and it lives inside the picture: the card offers the figure for saving and
  // a file that travels without its caption claims more than the measurement supports.
  const caption = `${CAPTION_SCALE} ${data.split ? CAPTION_GUARD : CAPTION_SOLO}`
    + `${data.split ? ` ${CAPTION_LEVELS}` : ''}${data.split && hasGain ? ` ${CAPTION_GAIN}` : ''}`
  return (
    <div ref={ref} style={{ marginTop: 14 }}>
      <WholeCell data={data} marks={marks} name={name} numbers={numbers}
        title={title} chips={chips} caption={caption} />
      {shown.length ? (
        <>
          <div style={{ fontSize: 13, fontWeight: 700, color: INK, marginTop: 12 }}>
            Each change along its own chromosome, numbered as above
          </div>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 18, marginTop: 6 }}>
            {shown.map((m) => (
              <Zoom key={`${m.chrom}:${m.startBp}`} data={data} mark={m} marks={marks}
                number={numbers.get(`${m.chrom}:${m.startBp}`)} width={colW} caption={zoomCaption(data, m)} />
            ))}
          </div>
        </>
      ) : null}
    </div>
  )
}

/* --- marks from rows ------------------------------------------------------------------------ */

const WHOSE: Record<Role, string> = { paternal: "the father's", maternal: "the mother's" }
const where = (chrom: string, s: number, e: number) =>
  `Chromosome ${chrom}, ${mb(s)} to ${mb(e)} Mb (${((e - s) / 1e6).toFixed(1)} Mb)`

/**
 * The finding in two or three words, for the list under the deck and the zoom's title.
 *
 * A class that is neither a plain loss nor a plain gain keeps its own words: a row the card calls
 * an imbalance must not read in the figure as a copy gone.
 */
function wordOf(gain: boolean, parent: Role | null, both: boolean, withheld?: string, kind?: string): string {
  if (withheld) return gain ? 'extra copy, not reported' : 'copy gone, not reported'
  const plain = !kind || kind === 'copy-loss' || kind === 'copy-gain' || kind === 'parental-absence'
    || kind === 'monosomy' || kind === 'trisomy' || kind === 'segmental-deletion'
    || kind === 'segmental-duplication'
  if (!plain) return KIND_WORD[kind] ?? (gain ? 'extra copy' : 'copy gone')
  if (parent) return `${parent === 'paternal' ? "father's" : "mother's"} ${gain ? 'doubled' : 'gone'}`
  return both ? (gain ? 'both doubled' : 'both gone') : gain ? 'extra copy' : 'copy gone'
}

/** A segment row as a mark. */
export function segmentMark(
  sg: Segment, loaded: Role | null, lines: string[], origin?: Role | 'unclear',
  /** The list calls this piece a second copy lost on a chromosome already lost whole, so the title
   *  says the same: a card must not say one thing in the row and another in the title. */
  deeper?: boolean,
): Mark {
  const co = segmentCoords(sg)
  const r = sg.refined
  const gain = sg.kind === 'copy-gain'
  let parent: Role | null
  let what: string
  if (sg.kind === 'copy-loss' || gain) {
    parent = origin === 'paternal' || origin === 'maternal' ? origin : null
    what = `${gain ? 'there is an extra copy' : deeper ? 'a second copy is lost' : 'one copy is lost'}`
      + `${parent ? `, ${WHOSE[parent]}` : ', parent not named'}`
  } else if (sg.contested) {
    parent = null; what = "both parents' alleles read missing, so no parent is named"
  } else {
    parent = sg.absentParent ?? loaded
    what = parent ? `${WHOSE[parent]} copy is missing` : "a parent's alleles are missing"
  }
  return {
    chrom: sg.chrom, startBp: co.start, endBp: co.end, outLoBp: r?.startLoBp, outHiBp: r?.endHiBp, whole: false, gain,
    parent, sentence: `${where(sg.chrom, co.start, co.end)}: ${what}.`, lines,
    word: wordOf(gain, parent, !!sg.contested),
  }
}

export const withheldRuns = (refusals: readonly Refusal[] | undefined): Segment[] =>
  (refusals ?? []).filter((x) => x.gate === 'other-parent.under-span' && x.withheld).map((x) => x.withheld as Segment)

export function withheldMark(sg: Segment, loaded: Role | null): Mark {
  const m = segmentMark(sg, loaded, [])
  const reason = `under ${OTHER_MIN_SPAN_BP / 1e6} Mb, no second line`
  return {
    ...m, withheld: reason, word: wordOf(m.gain, m.parent, false, reason),
    sentence: `${where(sg.chrom, m.startBp, m.endBp)}: ${m.parent ? WHOSE[m.parent] : "a parent's"} copy looks missing, but this is not reported: it is ${reason} of evidence.`,
    lines: [`${sg.markers.toLocaleString('en-US')} markers, score ${sg.score.toFixed(0)}`
      + (sg.lrrZ !== undefined ? `, intensity z ${sg.lrrZ.toFixed(2)}` : '')],
  }
}

/** A row of the list of changes as a mark; a whole-chromosome row is located by its chromosome. */
export function defectMark(
  d: Defect, segments: readonly Segment[], loaded: Role | null, endBp: number | undefined, deeper?: boolean,
): Mark | null {
  if (d.chrom === 'genome') return null
  const lines = [`${d.locus}${d.interval ? ` ${d.interval}` : ''}`,
    ...(d.informative !== undefined ? [`${d.informative.toLocaleString('en-US')} markers`] : []), d.why]
  const origin = d.origin
  const sg = segments.find((s) => s.chrom === d.chrom && segmentCoords(s).start === d.startBp && segmentCoords(s).end === d.endBp)
  if (sg) return { ...segmentMark(sg, loaded, [headline(d), ...lines], origin, deeper), quoted: true, lrr: d.lrr, said: saidOf(d) }
  const whole = d.locus === `chr${d.chrom}`
  const end = whole ? endBp ?? d.endBp : d.endBp
  const named = origin === 'paternal' || origin === 'maternal' ? origin : null
  const gain = d.kind === 'copy-gain' || d.kind === 'trisomy' || d.kind === 'segmental-duplication'
  const what = gain ? 'there is an extra copy' : d.kind === 'copy-loss' || d.kind === 'monosomy' ? 'one copy is lost'
    : d.kind === 'parental-absence' ? "one parent's copy is missing" : headline(d)
  return {
    chrom: d.chrom, startBp: whole ? 0 : d.startBp, endBp: end, whole, gain, parent: named,
    word: wordOf(gain, named, !!deeper, undefined, d.kind),
    sentence: whole ? `Whole chromosome ${d.chrom} (${(end / 1e6).toFixed(0)} Mb): ${what}${named ? `, ${WHOSE[named]}` : ', parent not named'}.`
      : `${where(d.chrom, d.startBp, d.endBp)}: ${what}${named ? `, ${WHOSE[named]}` : ''}.`,
    lines: [headline(d), ...lines],
    quoted: true,
    lrr: d.lrr,
  }
}

/** A joined row of the change list as one mark spanning its pieces. */
export function rowMark(
  pieces: readonly Defect[], segments: readonly Segment[], loaded: Role | null, endBp: number | undefined,
  deeper?: boolean,
): Mark | null {
  const first = defectMark(pieces[0], segments, loaded, endBp, deeper)
  if (!first || pieces.length < 2) return first
  const last = defectMark(pieces[pieces.length - 1], segments, loaded, endBp, deeper) ?? first
  const startBp = Math.min(...pieces.map((p) => p.startBp))
  const end = Math.max(...pieces.map((p) => p.endBp))
  const tail = first.sentence.slice(first.sentence.indexOf('):') + 1)
  return {
    ...first, startBp, endBp: end, outHiBp: last.outHiBp, lrr: undefined,
    said: pieces.map((p) => saidOf(p)).filter(Boolean).join(' ') || undefined,
    sentence: `${where(first.chrom, startBp, end)}${tail}`,
    lines: pieces.flatMap((p) => defectMark(p, segments, loaded, endBp, deeper)?.lines ?? []),
  }
}

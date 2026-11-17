/**
 * One chromosome as the lab reads it: allele fraction above, log R below, position in Mb across,
 * with a box over the called stretch. Hover a dot for that marker, the box for the call.
 */
import { useLayoutEffect, useRef, useState, type PointerEvent } from 'react'
import {
  boxGeometry, boxLrr, cellLrr, GT_TEXT, nearestIndex, segmentBox, tickStepMb, type ChromSeries, type PlotBox,
} from './plotSeries.ts'
import { headline, type Defect } from './defects.ts'
import { segmentCoords, type Segment } from './segments.ts'
import { OTHER_MIN_SPAN_BP } from './otherParentLoss.ts'
import type { Refusal } from './refusals.ts'

const L = 40
const R = 8
const TOP = 18
const H_BAF = 130
const GAP = 16
const H_LRR = 110
const AXIS = 30
const HEIGHT = TOP + H_BAF + GAP + H_LRR + AXIS
/** log R past this is drawn on the edge of its panel. */
const LRR_LIM = 1.5

const FATHER = 'var(--om-blue)'
const MOTHER = 'var(--om-higher)'
const OTHER = 'var(--om-text-dim)'
const BOX = 'var(--om-defect)'
const DIM = 'var(--om-text-dim)'

const yBaf = (v: number) => TOP + H_BAF * (1 - v)
const yLrr = (v: number) => TOP + H_BAF + GAP + H_LRR * (1 - (Math.max(-LRR_LIM, Math.min(LRR_LIM, v)) + LRR_LIM) / (2 * LRR_LIM))

/** Dots as one path per colour, so a few thousand points stay one element each. */
const dots = (xs: number[], ys: number[], r: number): string => {
  let d = ''
  for (let i = 0; i < xs.length; i += 1) {
    d += `M${(xs[i] - r).toFixed(1)} ${ys[i].toFixed(1)}a${r} ${r} 0 1 0 ${2 * r} 0a${r} ${r} 0 1 0 ${-2 * r} 0`
  }
  return d
}


/** The page's colour variables, written out, so the saved file looks the same outside the page. */
function saveSvg(svg: SVGSVGElement, filename: string): void {
  const css = getComputedStyle(document.documentElement)
  const text = svg.outerHTML.replace(/var\((--[a-z0-9-]+)\)/g, (_, v: string) => css.getPropertyValue(v).trim() || '#888')
  const a = document.createElement('a')
  a.href = URL.createObjectURL(new Blob([text], { type: 'image/svg+xml' }))
  a.download = filename
  a.click()
  setTimeout(() => URL.revokeObjectURL(a.href), 1000)
}

export function CallPlot({ chrom, series, boxes, name, cellLrr }: {
  chrom: string
  series: ChromSeries | undefined
  /** This cell's median log R over every chromosome, drawn as the line to read a box against. */
  cellLrr: number
  /** The first box is the row that opened the plot; the rest are drawn for context. */
  boxes: PlotBox[]
  /** For the saved file's name. */
  name: string
}) {
  const wrap = useRef<HTMLDivElement>(null)
  const svgRef = useRef<SVGSVGElement>(null)
  const [W, setW] = useState(640)
  const [tip, setTip] = useState<{ x: number; y: number; lines: string[] } | null>(null)
  useLayoutEffect(() => {
    const el = wrap.current
    if (!el) return
    const ro = new ResizeObserver(() => setW(Math.max(280, Math.floor(el.clientWidth))))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  if (!series || !series.pos.length) {
    return <div ref={wrap} style={{ fontSize: 12, color: DIM, margin: '6px 0' }}>No markers on chr{chrom} to plot.</div>
  }
  const endBp = Math.max(series.endBp, ...boxes.map((b) => b.endHiBp ?? b.endBp))
  const PW = W - L - R
  const x = (bp: number) => L + PW * (bp / endBp)
  const n = series.pos.length
  const xs: number[][] = [[], [], []]
  const yb: number[][] = [[], [], []]
  const xl: number[][] = [[], [], []]
  const yl: number[][] = [[], [], []]
  for (let i = 0; i < n; i += 1) {
    const w = series.who[i]
    const px = x(series.pos[i])
    if (Number.isFinite(series.baf[i])) { xs[w].push(px); yb[w].push(yBaf(series.baf[i])) }
    if (Number.isFinite(series.lrr[i])) { xl[w].push(px); yl[w].push(yLrr(series.lrr[i])) }
  }
  const colour = [OTHER, FATHER, MOTHER]
  const step = tickStepMb(endBp)
  const ticks: number[] = []
  for (let t = 0; t * step * 1e6 <= endBp; t += 1) ticks.push(t * step)
  const yBottom = TOP + H_BAF + GAP + H_LRR

  const onMove = (e: PointerEvent<SVGSVGElement>) => {
    const r = e.currentTarget.getBoundingClientRect()
    const px = e.clientX - r.left
    const py = e.clientY - r.top
    const inBaf = py >= TOP - 4 && py <= TOP + H_BAF + 4
    const inLrr = py >= TOP + H_BAF + GAP - 4 && py <= yBottom + 4
    if (px < L || px > W - R || (!inBaf && !inLrr)) { setTip(null); return }
    const bp = ((px - L) / PW) * endBp
    // The nearest dot within a few pixels, in the panel under the pointer.
    let best = -1
    let bestD = 36
    const i0 = nearestIndex(series.pos, bp)
    for (const dir of [-1, 1]) {
      for (let i = dir < 0 ? i0 : i0 + 1; i >= 0 && i < n; i += dir) {
        const dx = x(series.pos[i]) - px
        if (Math.abs(dx) > 6) break
        const v = inBaf ? series.baf[i] : series.lrr[i]
        if (!Number.isFinite(v)) continue
        const dy = (inBaf ? yBaf(v) : yLrr(v)) - py
        const d = dx * dx + dy * dy
        if (d < bestD) { bestD = d; best = i }
      }
    }
    if (best >= 0) {
      const w = series.who[best]
      setTip({ x: px, y: py, lines: [
        `chr${chrom}:${series.pos[best].toLocaleString('en-US')}`,
        `call ${GT_TEXT[series.gt[best]]}${w === 1 ? ", the father's letter only" : w === 2 ? ", the mother's letter only" : ''}`,
        `allele fraction ${Number.isFinite(series.baf[best]) ? series.baf[best].toFixed(3) : 'none'}`,
        `log R ${Number.isFinite(series.lrr[best]) ? series.lrr[best].toFixed(3) : 'none'}`,
      ] })
      return
    }
    const b = boxes.find((bx) => bp >= Math.min(bx.startBp, bx.startLoBp ?? Infinity)
      && bp <= Math.max(bx.endBp, bx.endHiBp ?? -Infinity))
    setTip(b ? { x: px, y: py, lines: b.lines } : null)
  }

  return (
    <div ref={wrap} style={{ position: 'relative', margin: '8px 0 2px', maxWidth: '100%' }}>
      <svg
        ref={svgRef} xmlns="http://www.w3.org/2000/svg" width={W} height={HEIGHT}
        viewBox={`0 0 ${W} ${HEIGHT}`} fontFamily="system-ui, sans-serif" fontSize={10}
        role="img" aria-label={`chr${chrom}: allele fraction and log R along the chromosome, ${boxes[0]?.label ?? ''}`}
        style={{ display: 'block', touchAction: 'pan-y' }}
        onPointerMove={onMove} onPointerDown={onMove} onPointerLeave={() => setTip(null)}
      >
        <rect x={0} y={0} width={W} height={HEIGHT} fill="#fff" />
        {/* Boxes under the dots, so the dots inside stay readable. */}
        {boxes.map((b, i) => {
          const g = boxGeometry(b, endBp, L, PW)
          const c = b.withheld ? DIM : BOX
          return (
            <g key={`f${i}`}>
              {[g.left, g.right].map((iv, k) => (iv ? (
                <rect key={k} x={iv[0]} y={TOP} width={Math.max(1, iv[1] - iv[0])} height={yBottom - TOP} fill={c} opacity={0.16} />
              ) : null))}
              <rect x={g.x0} y={TOP} width={g.x1 - g.x0} height={yBottom - TOP} fill={c} opacity={0.07} />
            </g>
          )
        })}
        {/* Panel frames and reference lines. */}
        <rect x={L} y={TOP} width={PW} height={H_BAF} fill="none" stroke="var(--om-border)" />
        <rect x={L} y={TOP + H_BAF + GAP} width={PW} height={H_LRR} fill="none" stroke="var(--om-border)" />
        {[0.5].map((v) => <line key={v} x1={L} x2={L + PW} y1={yBaf(v)} y2={yBaf(v)} stroke="var(--om-border)" strokeDasharray="3,3" />)}
        <line x1={L} x2={L + PW} y1={yLrr(0)} y2={yLrr(0)} stroke="var(--om-border)" strokeDasharray="3,3" />
        {Number.isFinite(cellLrr) ? (
          <g>
            <line x1={L} x2={L + PW} y1={yLrr(cellLrr)} y2={yLrr(cellLrr)} stroke={DIM} strokeDasharray="6,3" />
            <text x={L + PW - 3} y={yLrr(cellLrr) - 3} textAnchor="end" fill={DIM}>cell median</text>
          </g>
        ) : null}
        {[0, 0.5, 1].map((v) => (
          <text key={v} x={L - 4} y={yBaf(v) + 3} textAnchor="end" fill={DIM}>{v}</text>
        ))}
        {[-LRR_LIM, 0, LRR_LIM].map((v) => (
          <text key={v} x={L - 4} y={yLrr(v) + 3} textAnchor="end" fill={DIM}>{v > 0 ? `+${v}` : v}</text>
        ))}
        <text x={L + 4} y={TOP + 11} fill={DIM}>allele fraction</text>
        <text x={L + 4} y={TOP + H_BAF + GAP + 11} fill={DIM}>log R</text>
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t * 1e6)} x2={x(t * 1e6)} y1={yBottom} y2={yBottom + 4} stroke={DIM} />
            <text x={x(t * 1e6)} y={yBottom + 14} textAnchor="middle" fill={DIM}>{t}</text>
          </g>
        ))}
        <text x={L + PW} y={yBottom + 26} textAnchor="end" fill={DIM}>chr{chrom}, Mb</text>
        {/* Grey first, then each parent's letter on top. */}
        {[0, 1, 2].map((w) => (
          <g key={w} fill={colour[w]} opacity={w ? 0.9 : 0.35}>
            <path d={dots(xs[w], yb[w], w ? 1.6 : 1.2)} />
            <path d={dots(xl[w], yl[w], w ? 1.6 : 1.2)} />
          </g>
        ))}
        {boxes.map((b, i) => {
          const g = boxGeometry(b, endBp, L, PW)
          const c = b.withheld ? DIM : BOX
          const right = g.x0 > L + PW * 0.6
          const mid = boxLrr(series, b.startBp, b.endBp)
          return (
            <g key={`o${i}`}>
              {Number.isFinite(mid) ? (
                <line x1={g.x0} x2={g.x1} y1={yLrr(mid)} y2={yLrr(mid)} stroke={c} strokeWidth={2.5} />
              ) : null}
              <rect
                x={g.x0} y={TOP} width={g.x1 - g.x0} height={yBottom - TOP} fill="none" stroke={c}
                strokeWidth={i === 0 ? 2 : 1.2} strokeDasharray={b.withheld ? '5,3' : undefined}
              />
              <text
                x={right ? g.x1 : g.x0} y={TOP - 5} textAnchor={right ? 'end' : 'start'}
                fill={c} fontWeight={700}
              >
                {b.label}
              </text>
            </g>
          )
        })}
      </svg>
      {tip ? (
        <div style={{
          position: 'absolute', left: Math.max(0, Math.min(tip.x + 12, W - 280)), top: tip.y + 14,
          width: 'max-content', maxWidth: Math.min(360, W), pointerEvents: 'none', zIndex: 2,
          background: '#fff', border: '1px solid var(--om-border-strong)', padding: '4px 7px',
          fontSize: 11, lineHeight: 1.45,
        }}
        >
          {tip.lines.map((l, i) => <div key={i} style={{ fontWeight: i === 0 ? 700 : 400 }}>{l}</div>)}
        </div>
      ) : null}
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '2px 14px', fontSize: 11, color: DIM, marginTop: 3 }}>
        {xs[1].length + xl[1].length ? (
          <span><span style={{ color: FATHER }}>&#9679;</span> the father&rsquo;s letter only</span>
        ) : null}
        {xs[2].length + xl[2].length ? (
          <span><span style={{ color: MOTHER }}>&#9679;</span> the mother&rsquo;s letter only</span>
        ) : null}
        <span><span style={{ color: BOX }}>&#9472;</span> median log R in the box</span>
        <span>{n.toLocaleString('en-US')} of {series.total.toLocaleString('en-US')} markers shown</span>
        {boxes.some((b) => b.startLoBp !== undefined) ? <span>pale ends: where the edge may sit</span> : null}
        <button
          type="button"
          onClick={() => svgRef.current && saveSvg(svgRef.current, `${name}-chr${chrom}.svg`)}
          style={{
            border: '1px solid var(--om-border-strong)', background: '#fff', font: 'inherit',
            fontSize: 11, padding: '0 6px', cursor: 'pointer', marginLeft: 'auto',
          }}
        >
          Save SVG
        </button>
      </div>
    </div>
  )
}

const mbSpan = (a: number, b: number) => `${((b - a) / 1e6).toFixed(2)} Mb`

/** Runs under the span floor that were found and withheld, with their coordinates. */
export const withheldRuns = (refusals: readonly Refusal[] | undefined): Segment[] =>
  (refusals ?? []).filter((x) => x.gate === 'other-parent.under-span' && x.withheld)
    .map((x) => x.withheld as Segment)

/** The box for a withheld run, with what the limit says about it. */
export function withheldBox(sg: Segment): PlotBox {
  const co = segmentCoords(sg)
  return segmentBox(sg, 'withheld', [
    `withheld: ${sg.absentParent ?? 'other parent'} alleles absent, chr${sg.chrom}`,
    `${co.start.toLocaleString('en-US')}-${co.end.toLocaleString('en-US')}, ${co.interval}`,
    mbSpan(co.start, co.end),
    `${sg.markers.toLocaleString('en-US')} markers, score ${sg.score.toFixed(0)}`
      + (sg.lrrZ !== undefined ? `, intensity z ${sg.lrrZ.toFixed(2)}` : ''),
    `Under ${OTHER_MIN_SPAN_BP / 1e6} Mb a stretch both parents share reads the same in the alleles, `
      + 'and here neither intensity nor another array of the embryo shows the loss.',
  ], true)
}

/** The plot for one row of the list of changes: its own edges, and any withheld run on the
 *  same chromosome. */
export function DefectPlot({ d, segments, refusals, series, name }: {
  d: Defect
  segments: readonly Segment[]
  refusals: readonly Refusal[] | undefined
  series: ReadonlyMap<string, ChromSeries> | undefined
  name: string
}) {
  if (d.chrom === 'genome' || !series) return null
  const s = series.get(d.chrom)
  // A whole-chromosome row is located by its chromosome alone; its box is the whole chromosome,
  // whatever end the row carries.
  const whole = d.locus === `chr${d.chrom}`
  const endBp = whole ? s?.endBp ?? d.endBp : d.endBp
  const lines = [
    headline(d),
    `${d.locus}${d.interval ? ` ${d.interval}` : ''}`,
    whole ? `whole chromosome, ${mbSpan(0, endBp)}` : mbSpan(d.startBp, d.endBp),
    ...(d.informative !== undefined ? [`${d.informative.toLocaleString('en-US')} markers`] : []),
    d.why,
  ]
  const label = `${d.kind.replace(/-/g, ' ')}, ${d.origin === 'unclear' ? 'no parent named' : d.origin}`
  const sg = segments.find((s) => s.chrom === d.chrom && segmentCoords(s).start === d.startBp
    && segmentCoords(s).end === d.endBp)
  const own = sg ? segmentBox(sg, label, lines)
    : { startBp: whole ? 0 : d.startBp, endBp, label, lines }
  const held = withheldRuns(refusals).filter((w) => w.chrom === d.chrom && w !== sg).map(withheldBox)
  return (
    <CallPlot
      chrom={d.chrom} series={s} boxes={[own, ...held]} name={name} cellLrr={cellLrr(series)}
    />
  )
}

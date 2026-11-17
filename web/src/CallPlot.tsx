/**
 * The cell figures. The whole cell: every chromosome, one row per parent whose letters can be read,
 * with every printed call marked. One call: the same quantity closer in, with a plain sentence on top,
 * the stretch as a band and the cell's own normal range behind the bars; or, where no true quantity
 * separates such a stretch by eye, one line with the two numbers and the whole cell with it marked.
 */
import { useId, useLayoutEffect, useRef, useState, type PointerEvent, type ReactNode } from 'react'
import {
  BASE, cellPer, drawable, risesInside, shareTop, mergeBins, normalRange, perForStretch, runningScore, shareIn, tickStepMb, valueOf,
  viewOf, wrap, type Bin, type PlotData, type Role, type Track,
} from './plotSeries.ts'
import { headline, type Defect } from './defects.ts'
import { segmentCoords, type Segment } from './segments.ts'
import { OTHER_MIN_SPAN_BP } from './otherParentLoss.ts'
import type { Refusal } from './refusals.ts'

const FATHER = 'var(--om-blue)'
const MOTHER = 'var(--om-higher)'
const NEUTRAL = '#6b727b'
const CALL = 'var(--om-defect)'
const DIM = 'var(--om-text-dim)'
const INK = '#1d2126'
const RANGE = '#e9ecef'
const LINE = '#e3e6ea'
const colourOf = (p: Role | null | undefined) => (p === 'paternal' ? FATHER : p === 'maternal' ? MOTHER : NEUTRAL)
const AUTOSOMES = Array.from({ length: 22 }, (_, i) => String(i + 1))
const LEVEL_LO = -1.5
const LEVEL_HI = 1
const mb = (bp: number) => (bp / 1e6).toFixed(1)
const pct = (x: number) => `${Math.round(x * 100)}%`
const WORD: Record<string, string> = {
  father: "show the father's letter", mother: "show the mother's letter",
  'two-letter': 'read two letters', 'call-rate': 'were read',
}

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
  /** What the row shows, for the hover. */
  lines: string[]
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
        border: '1px solid var(--om-border-strong)', background: '#fff', font: 'inherit', fontSize: 11,
        padding: '0 6px', cursor: 'pointer', marginLeft: 'auto',
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

const binLine = (t: Track, chrom: string, b: Bin): string[] => [
  `chr${chrom} ${mb(b.startBp)} to ${mb(b.endBp)} Mb`,
  t.kind === 'share' ? `${b.n} markers, ${pct(b.k / b.n)} ${WORD[t.key] ?? ''}` : `${b.n} markers, ${b.k > 0 ? '+' : ''}${b.k.toFixed(2)} from the cell's level`,
  t.tech,
]

/* --- the whole cell ------------------------------------------------------------------------- */

export function WholeCell({ data, marks, name, focus }: {
  data: PlotData; marks: readonly Mark[]; name: string
  /** A chromosome to bring forward; the others are muted. */
  focus?: string
}) {
  const [ref, W] = useWidth()
  const svg = useRef<SVGSVGElement>(null)
  const [tip, setTip] = useState<{ x: number; y: number; lines: string[] } | null>(null)
  const tracks = data.tracks.filter((t) => t.shown)
  if (!tracks.length) return null
  const pers = tracks.map(cellPer)
  const called = new Set(marks.filter((m) => !m.withheld).map((m) => m.chrom))
  const ranges = tracks.map((t, i) => normalRange(t, pers[i], called))
  const tops = tracks.map((t, i) => (t.kind === 'share' ? shareTop(ranges[i].mid) : 1))
  const rowLabel = (t: Track, i: number) => (tops[i] < 1 ? `${t.label}, against this cell's own level (row top ${pct(tops[i])})` : t.label)
  const rows = W >= 760 ? [AUTOSOMES] : [AUTOSOMES.slice(0, 11), AUTOSOMES.slice(11)]
  const G = 3
  const TH = 58
  const block = tracks.length * (TH + 18) + 20
  const H = 6 + rows.length * block + 4
  const len = (c: string) => data.endBp.get(c) ?? 1
  type Panel = { c: string; x: number; w: number; top: number }
  const panels: Panel[] = []
  rows.forEach((row, r) => {
    const total = row.reduce((s, c) => s + len(c), 0)
    const usable = W - 2 - G * (row.length - 1)
    let x = 1
    for (const c of row) {
      const w = usable * len(c) / total
      panels.push({ c, x, w, top: 6 + r * block })
      x += w + G
    }
  })
  const yOf = (t: Track, top: number, v: number) => (t.kind === 'share' ? top + TH * (1 - Math.max(0, Math.min(1, v / tops[tracks.indexOf(t)])))
    : top + TH * (LEVEL_HI - Math.max(LEVEL_LO, Math.min(LEVEL_HI, v))) / (LEVEL_HI - LEVEL_LO))
  const marksFor = (c: string, t: Track) => marks.filter((m) => m.chrom === c && (!m.parent || !t.parent || m.parent === t.parent))
  const onMove = (e: PointerEvent<SVGSVGElement>) => {
    const r = e.currentTarget.getBoundingClientRect()
    const px = e.clientX - r.left; const py = e.clientY - r.top
    const p = panels.find((q) => px >= q.x && px <= q.x + q.w && py >= q.top && py <= q.top + block)
    if (!p) { setTip(null); return }
    const ti = Math.floor((py - p.top) / (TH + 18))
    const t = tracks[ti]
    if (!t) { setTip(null); return }
    const bp = (px - p.x) / p.w * len(p.c)
    const m = marksFor(p.c, t).find((q) => bp >= q.startBp && bp <= q.endBp)
    if (m) { setTip({ x: px, y: py, lines: [m.sentence, ...m.lines.slice(0, 3)] }); return }
    const b = mergeBins(t.bins.get(p.c) ?? [], pers[ti], t.kind === 'level').find((q) => bp >= q.startBp && bp <= q.endBp)
    if (!b) { setTip(null); return }
    // A bar under the cell's normal range with no call over it: say it is not a call, and why.
    const below = t.kind === 'share' ? valueOf(t, b) < ranges[ti].lo : valueOf(t, b) < -0.3
    setTip({ x: px, y: py, lines: [...binLine(t, p.c, b), ...(below ? [
      "Below this cell's normal range, but not a call: the stretch did not pass the tool's tests (too short, too few markers, or a score under the bar)."] : [])] })
  }
  return (
    <div ref={ref} style={{ position: 'relative', margin: '8px 0 2px', maxWidth: '100%' }}>
      <svg ref={svg} xmlns="http://www.w3.org/2000/svg" width={W} height={H} viewBox={`0 0 ${W} ${H}`}
        fontFamily="system-ui, -apple-system, sans-serif" fontSize={11} role="img"
        aria-label="Every chromosome of this cell, each parent's letters present"
        style={{ display: 'block', touchAction: 'pan-y' }}
        onPointerMove={onMove} onPointerDown={onMove} onPointerLeave={() => setTip(null)}
      >
        <rect width={W} height={H} fill="#fff" />
        {rows.map((_, r) => tracks.map((t, ti) => (
          <text key={`l${r}${ti}`} x={1} y={6 + r * block + ti * (TH + 18) + 10} fontWeight={600} fill={colourOf(t.parent)}>
            {rowLabel(t, ti)}
          </text>
        )))}
        {panels.map((p) => {
          const muted = focus && p.c !== focus
          return (
            <g key={p.c} opacity={muted ? 0.35 : 1}>
              {tracks.map((t, ti) => {
                const top = p.top + ti * (TH + 18) + 14
                const x = (bp: number) => p.x + p.w * bp / len(p.c)
                const rg = ranges[ti]
                return (
                  <g key={t.key}>
                    <rect x={p.x} y={top} width={p.w} height={TH} fill="none" stroke={LINE} />
                    {Number.isFinite(rg.lo) ? (
                      <rect x={p.x} y={yOf(t, top, rg.hi)} width={p.w} height={Math.max(1, yOf(t, top, rg.lo) - yOf(t, top, rg.hi))} fill={RANGE} />
                    ) : null}
                    {marksFor(p.c, t).map((m, i) => (
                      <rect key={i} x={x(m.startBp)} y={top} width={Math.max(2, x(m.endBp) - x(m.startBp))} height={TH}
                        fill={m.withheld ? 'none' : CALL} opacity={m.withheld ? 1 : 0.1}
                        stroke={m.withheld ? DIM : 'none'} strokeDasharray={m.withheld ? '3,2' : undefined} />
                    ))}
                    <path fill={colourOf(t.parent)} d={mergeBins(t.bins.get(p.c) ?? [], pers[ti], t.kind === 'level').map((b) => {
                      const y0 = t.kind === 'share' ? top + TH : yOf(t, top, 0)
                      const y1 = yOf(t, top, valueOf(t, b))
                      const x0 = x(b.startBp); const x1 = Math.max(x0 + 0.6, x(b.endBp))
                      const ya = Math.min(y0, y1); const yb = Math.max(ya + 0.6, Math.max(y0, y1))
                      return `M${x0.toFixed(1)} ${ya.toFixed(1)}H${x1.toFixed(1)}V${yb.toFixed(1)}H${x0.toFixed(1)}Z`
                    }).join('')} />
                  </g>
                )
              })}
              {marks.filter((m) => m.chrom === p.c).map((m, i) => {
                const x0 = p.x + p.w * m.startBp / len(p.c)
                const x1 = p.x + p.w * m.endBp / len(p.c)
                const y = p.top + tracks.length * (TH + 18) + 4
                return (
                  <line key={i} x1={x0} x2={Math.max(x0 + 2, x1)} y1={y} y2={y} strokeWidth={4}
                    stroke={m.withheld ? DIM : colourOf(m.parent) === NEUTRAL ? CALL : colourOf(m.parent)}
                    strokeDasharray={m.withheld ? '3,2' : undefined} />
                )
              })}
              <text x={p.x + p.w / 2} y={p.top + tracks.length * (TH + 18) + 18} textAnchor="middle"
                fontWeight={p.c === focus || called.has(p.c) ? 700 : 400} fill={called.has(p.c) ? INK : DIM}
              >
                {p.c}
              </text>
            </g>
          )
        })}
      </svg>
      <Tip tip={tip} w={W} />
      <div style={{ display: 'flex', fontSize: 11, color: DIM, marginTop: 2 }}>
        <span>
          Grey band: this cell&apos;s normal range. A line under a chromosome marks a call; a dip with no line under it is not a call.{' '}
          {marks.some((m) => m.withheld) ? 'Dashed: found, not reported. ' : ''}
          {tracks.some((t) => t.parent) ? 'A parent\'s bar is the share of markers where that parent\'s letter could be seen at which it is. ' : ''}
          Every bar of a row holds as many markers.
        </span>
        <SaveSvg svg={svg} name={`${name}-cell`} />
      </div>
    </div>
  )
}

/* --- one call ------------------------------------------------------------------------------- */

export function CallFigure({ mark, data, marks, name }: {
  mark: Mark; data: PlotData | undefined; marks: readonly Mark[]; name: string
}) {
  const [ref, W] = useWidth()
  const svg = useRef<SVGSVGElement>(null)
  const clip = useId().replace(/:/g, '')
  const [tip, setTip] = useState<{ x: number; y: number; lines: string[] } | null>(null)
  if (!data) return null
  const pick = drawable({ whole: mark.whole, spanBp: mark.endBp - mark.startBp, gain: mark.gain, parent: mark.parent }, data)
  const called = new Set(marks.filter((m) => !m.withheld).map((m) => m.chrom))
  const quote = (t: Track | null): string => {
    if (!t || t.kind !== 'share') return ''
    const ins = shareIn(t.bins.get(mark.chrom) ?? [], mark.startBp, mark.endBp)
    const rg = normalRange(t, cellPer(t), called)
    return Number.isFinite(ins.share)
      ? `Inside the stretch ${pct(ins.share)} of these markers ${WORD[t.key]}; on this cell's chromosomes with no call, ${pct(rg.mid)} (normal range ${pct(rg.lo)} to ${pct(rg.hi)}).`
      : ''
  }
  if (!pick.draw || !pick.track) {
    const t = pick.track ?? data.tracks.find((x) => x.key === (mark.parent ? (mark.parent === 'paternal' ? 'father' : 'mother') : ''))
      ?? data.tracks.find((x) => x.key === 'two-letter' || x.key === 'call-rate') ?? null
    const intensity = data.tracks.find((x) => x.key === 'intensity')
    const lvl = intensity ? mergeBins(intensity.bins.get(mark.chrom) ?? [], 1, true)
      .filter((b) => b.startBp >= mark.startBp && b.endBp <= mark.endBp).map((b) => b.k) : []
    const lvlText = lvl.length ? ` Intensity inside it sits ${lvl.sort((a, b) => a - b)[lvl.length >> 1].toFixed(2)} from the cell's level.` : ''
    return (
      <div ref={ref} style={{ margin: '8px 0 2px' }}>
        <div style={{ fontSize: 14, fontWeight: 700, color: INK, lineHeight: 1.35 }}>{mark.sentence}</div>
        <div style={{ fontSize: 12, color: INK, lineHeight: 1.5, marginTop: 3, maxWidth: 780 }}>
          {quote(t)}{lvlText} A single cell cannot show this stretch by eye
          {mark.endBp - mark.startBp < 6e6 ? ': it is shorter than one band' : ''}, so it is not drawn on its own.
        </div>
        <WholeCell data={data} marks={marks.filter((m) => m.chrom === mark.chrom)} name={name} focus={mark.chrom} />
      </div>
    )
  }
  const t = pick.track
  const bins = t.bins.get(mark.chrom) ?? []
  const chromEnd = Math.max(data.endBp.get(mark.chrom) ?? 0, mark.outHiBp ?? mark.endBp)
  const { loBp, hiBp } = mark.whole ? { loBp: 0, hiBp: chromEnd } : viewOf(mark.startBp, mark.endBp, chromEnd)
  const per = mark.whole ? cellPer(t) : perForStretch(bins, mark.startBp, mark.endBp)
  const bars = mergeBins(bins, per, t.kind === 'level').filter((b) => b.endBp >= loBp && b.startBp <= hiBp)
  const rg = normalRange(t, per, called)
  const colour = mark.withheld ? DIM : colourOf(t.parent ?? mark.parent)
  const L = 46; const R = 8; const PW = W - L - R
  const x = (bp: number) => L + PW * Math.max(0, Math.min(1, (bp - loBp) / (hiBp - loBp)))
  const head = wrap(mark.sentence, Math.max(30, Math.floor(W / 8)))
  let y = 8 + head.length * 19
  // A whole chromosome first, against the cell's others.
  const strip = mark.whole ? { top: y + 20, h: 70 } : null
  const pooled = strip ? AUTOSOMES.map((c) => {
    const bs = t.bins.get(c) ?? []
    return t.kind === 'share' ? shareIn(bs, 0, Infinity).share
      : (() => { const v = bs.map((b) => b.k).sort((a, b) => a - b); return v.length ? v[v.length >> 1] : NaN })()
  }) : []
  if (strip) y = strip.top + strip.h + 26
  const top = y + 44
  const TH = 120
  const vTop = t.kind === 'share' ? shareTop(rg.mid) : 1
  const yv = (v: number) => (t.kind === 'share' ? top + TH * (1 - Math.max(0, Math.min(1, v / vTop)))
    : top + TH * (LEVEL_HI - Math.max(LEVEL_LO, Math.min(LEVEL_HI, v))) / (LEVEL_HI - LEVEL_LO))
  const score = !mark.whole && t.kind === 'share' && Number.isFinite(data.dropIn) && t.parent
    ? runningScore(bars, rg.mid, data.dropIn) : null
  const scoreOk = !!score && risesInside(score, bars, mark.startBp, mark.endBp)
  const sTop = top + TH + 34
  const SH = 44
  const smax = scoreOk ? Math.max(...score!.filter(Number.isFinite), 1) : 1
  const words = quote(t)
  const under = wrap(words, Math.max(40, Math.floor((W - L) / 6.6)))
  const uTop = (scoreOk ? sTop + SH + 8 : top + TH + 8) + 16
  const axisY = uTop + under.length * 15 + 4
  const H = axisY + 32
  const span = hiBp - loBp
  const step = tickStepMb(span)
  const ticks: number[] = []
  for (let v = Math.ceil(loBp / 1e6 / step) * step; v * 1e6 <= hiBp; v += step) ticks.push(v)
  const xs = (i: number) => L + PW * (i + 0.5) / 22
  const ys = (v: number) => strip!.top + strip!.h * (1 - Math.max(0, Math.min(1, v / vTop)))
  const ysl = (v: number) => strip!.top + strip!.h * (LEVEL_HI - Math.max(LEVEL_LO, Math.min(LEVEL_HI, v))) / (LEVEL_HI - LEVEL_LO)
  const sy = (v: number) => (t.kind === 'share' ? ys(v) : ysl(v))
  const stripOthers = pooled.filter((v, i) => Number.isFinite(v) && !called.has(AUTOSOMES[i])).sort((a, b) => a - b)
  const sLo = stripOthers[Math.floor(stripOthers.length * 0.1)]
  const sHi = stripOthers[Math.min(stripOthers.length - 1, Math.floor(stripOthers.length * 0.9))]
  const onMove = (e: PointerEvent<SVGSVGElement>) => {
    const r = e.currentTarget.getBoundingClientRect()
    const px = e.clientX - r.left; const py = e.clientY - r.top
    if (strip && py >= strip.top && py <= strip.top + strip.h + 14) {
      const i = Math.floor((px - L) / (PW / 22))
      const v = pooled[i]
      setTip(i >= 0 && i < 22 && Number.isFinite(v) ? { x: px, y: py, lines: [`chr${AUTOSOMES[i]}`,
        t.kind === 'share' ? `${pct(v)} ${WORD[t.key]}` : `${v.toFixed(2)} from the cell's level`, t.tech] } : null)
      return
    }
    if (px < L || px > L + PW || py < top - 26 || py > top + TH) { setTip(null); return }
    const bp = loBp + (px - L) / PW * span
    if (bp >= mark.startBp && bp <= mark.endBp && py < top + 4) { setTip({ x: px, y: py, lines: [mark.sentence, ...mark.lines] }); return }
    const b = bars.find((q) => bp >= q.startBp && bp <= q.endBp)
    setTip(b ? { x: px, y: py, lines: binLine(t, mark.chrom, b) } : null)
  }
  const label = mark.withheld ? `found, not reported: ${mark.withheld}` : 'called stretch'
  const labelRight = x(mark.startBp) > L + PW * 0.55
  return (
    <div ref={ref} style={{ position: 'relative', margin: '8px 0 2px', maxWidth: '100%' }}>
      <svg ref={svg} xmlns="http://www.w3.org/2000/svg" width={W} height={H} viewBox={`0 0 ${W} ${H}`}
        fontFamily="system-ui, -apple-system, sans-serif" fontSize={12} role="img" aria-label={mark.sentence}
        style={{ display: 'block', touchAction: 'pan-y' }}
        onPointerMove={onMove} onPointerDown={onMove} onPointerLeave={() => setTip(null)}
      >
        <defs>
          <pattern id={`h${clip}`} width={6} height={6} patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
            <line x1={0} y1={0} x2={0} y2={6} stroke={DIM} strokeWidth={1.5} opacity={0.35} />
          </pattern>
        </defs>
        <rect width={W} height={H} fill="#fff" />
        {head.map((l, i) => <text key={i} x={2} y={18 + i * 19} fontSize={15} fontWeight={700} fill={INK}>{l}</text>)}
        {strip ? (
          <g>
            <text x={L} y={strip.top - 6} fontWeight={600} fill={INK}>{t.label}, each chromosome of this cell</text>
            {Number.isFinite(sLo) ? <rect x={L} y={sy(sHi)} width={PW} height={Math.max(1, sy(sLo) - sy(sHi))} fill={RANGE} /> : null}
            {pooled.map((v, i) => {
              if (!Number.isFinite(v)) return null
              const me = AUTOSOMES[i] === mark.chrom
              const y0 = t.kind === 'share' ? strip.top + strip.h : ysl(0)
              return (
                <g key={i}>
                  <rect x={xs(i) - PW / 22 * 0.3} y={Math.min(y0, sy(v))} width={PW / 22 * 0.6} height={Math.max(1, Math.abs(sy(v) - y0))}
                    fill={me ? colour : '#c4c9cf'} />
                  <text x={xs(i)} y={strip.top + strip.h + 13} textAnchor="middle" fontSize={me ? 12 : 10}
                    fontWeight={me ? 700 : 400} fill={me ? INK : DIM}>{AUTOSOMES[i]}</text>
                </g>
              )
            })}
          </g>
        ) : null}
        <text x={L} y={top - 30} fontWeight={600} fill={INK}>{t.label}{vTop < 1 ? ", against this cell's own level" : ''}{mark.whole ? `, along chromosome ${mark.chrom}` : ''}</text>
        {/* The cell's normal range, then the call, then the bars. */}
        {Number.isFinite(rg.lo) ? <rect x={L} y={yv(rg.hi)} width={PW} height={Math.max(1, yv(rg.lo) - yv(rg.hi))} fill={RANGE} /> : null}
        {mark.outLoBp !== undefined ? <rect x={x(mark.outLoBp)} y={top - 4} width={Math.max(0, x(mark.startBp) - x(mark.outLoBp))} height={TH + 4} fill={colour} opacity={0.06} /> : null}
        {mark.outHiBp !== undefined ? <rect x={x(mark.endBp)} y={top - 4} width={Math.max(0, x(mark.outHiBp) - x(mark.endBp))} height={TH + 4} fill={colour} opacity={0.06} /> : null}
        <rect x={x(mark.startBp)} y={top - 4} width={Math.max(3, x(mark.endBp) - x(mark.startBp))} height={TH + 4}
          fill={mark.withheld ? `url(#h${clip})` : colour} opacity={mark.withheld ? 1 : 0.12} />
        {[mark.startBp, mark.endBp].map((v, i) => (
          <line key={i} x1={x(v)} x2={x(v)} y1={top - 18} y2={top + TH} stroke={colour} strokeWidth={1.5}
            strokeDasharray={mark.withheld ? '4,3' : undefined} />
        ))}
        <text x={labelRight ? x(mark.endBp) : x(mark.startBp)} y={top - 7} textAnchor={labelRight ? 'end' : 'start'}
          fontWeight={700} fill={colour}>{label}</text>
        {[0, vTop / 2, vTop].map((v) => (t.kind === 'share' ? (
          <text key={v} x={L - 5} y={yv(v) + 4} textAnchor="end" fill={DIM}>{pct(v)}</text>
        ) : null))}
        {t.kind === 'level' ? <text x={L - 5} y={yv(0) + 4} textAnchor="end" fill={DIM}>level</text> : null}
        {bars.map((b, i) => {
          const v = valueOf(t, b)
          const inside = b.startBp >= mark.startBp && b.endBp <= mark.endBp
          const y0 = t.kind === 'share' ? top + TH : yv(0)
          return (
            <rect key={i} x={x(b.startBp) + 0.4} y={Math.min(y0, yv(v))} width={Math.max(0.8, x(b.endBp) - x(b.startBp) - 0.8)}
              height={Math.max(0.8, Math.abs(yv(v) - y0))} fill={colourOf(t.parent)} opacity={inside || mark.whole ? 1 : 0.5} />
          )
        })}
        {scoreOk ? (
          <g>
            <text x={L} y={sTop - 6} fontWeight={600} fill={INK}>
              {W < 560 ? 'Evidence the copy is missing, added up' : 'Evidence the copy is missing, added along the chromosome'}
            </text>
            <line x1={L} x2={L + PW} y1={sTop + SH} y2={sTop + SH} stroke={LINE} />
            <path fill="none" stroke={colour} strokeWidth={2}
              d={bars.map((b, i) => `${i ? 'L' : 'M'}${x(b.startBp).toFixed(1)} ${(sTop + SH * (1 - (score![i] || 0) / smax)).toFixed(1)}H${x(b.endBp).toFixed(1)}`).join('')} />
          </g>
        ) : null}
        {under.map((l, i) => <text key={i} x={L} y={uTop + i * 15} fill={INK}>{l}</text>)}
        <line x1={L} x2={L + PW} y1={axisY} y2={axisY} stroke="#8a9099" />
        {ticks.map((v) => (
          <g key={v}>
            <line x1={x(v * 1e6)} x2={x(v * 1e6)} y1={axisY} y2={axisY + 4} stroke={DIM} />
            <text x={x(v * 1e6)} y={axisY + 16} textAnchor="middle" fill={DIM}>{v}</text>
          </g>
        ))}
        <text x={L + PW} y={axisY + 29} textAnchor="end" fill={DIM}>
          {loBp > 0 || hiBp < chromEnd ? `${mb(loBp)} to ${mb(hiBp)} Mb of chromosome ${mark.chrom}` : `chromosome ${mark.chrom}, Mb`}
        </text>
      </svg>
      <Tip tip={tip} w={W} />
      <div style={{ display: 'flex', fontSize: 11, color: DIM, marginTop: 2 }}>
        <span>Bars of {per * BASE} markers. Grey band: this cell&apos;s normal range, the middle 80% of its chromosomes with no call.</span>
        <SaveSvg svg={svg} name={`${name}-chr${mark.chrom}`} />
      </div>
    </div>
  )
}

/* --- marks from rows ------------------------------------------------------------------------ */

const WHOSE: Record<Role, string> = { paternal: "the father's", maternal: "the mother's" }
const where = (chrom: string, s: number, e: number) =>
  `Chromosome ${chrom}, ${mb(s)} to ${mb(e)} Mb (${((e - s) / 1e6).toFixed(1)} Mb)`

/** A segment row as a mark. */
export function segmentMark(sg: Segment, loaded: Role | null, lines: string[], origin?: Role | 'unclear'): Mark {
  const co = segmentCoords(sg)
  const r = sg.refined
  const gain = sg.kind === 'copy-gain'
  let parent: Role | null
  let what: string
  if (sg.kind === 'copy-loss' || gain) {
    parent = origin === 'paternal' || origin === 'maternal' ? origin : null
    what = `${gain ? 'there is an extra copy' : 'one copy is lost'}${parent ? `, ${WHOSE[parent]}` : ', parent not named'}`
  } else if (sg.contested) {
    parent = null; what = "both parents' letters read missing, so no parent is named"
  } else {
    parent = sg.absentParent ?? loaded
    what = parent ? `${WHOSE[parent]} copy is missing` : "a parent's letters are missing"
  }
  return {
    chrom: sg.chrom, startBp: co.start, endBp: co.end, outLoBp: r?.startLoBp, outHiBp: r?.endHiBp, whole: false, gain,
    parent, sentence: `${where(sg.chrom, co.start, co.end)}: ${what}.`, lines,
  }
}

export const withheldRuns = (refusals: readonly Refusal[] | undefined): Segment[] =>
  (refusals ?? []).filter((x) => x.gate === 'other-parent.under-span' && x.withheld).map((x) => x.withheld as Segment)

export function withheldMark(sg: Segment, loaded: Role | null): Mark {
  const m = segmentMark(sg, loaded, [])
  const reason = `under ${OTHER_MIN_SPAN_BP / 1e6} Mb, no second line`
  return {
    ...m, withheld: reason,
    sentence: `${where(sg.chrom, m.startBp, m.endBp)}: ${m.parent ? WHOSE[m.parent] : "a parent's"} copy looks missing, but this is not reported: it is ${reason} of evidence.`,
    lines: [`${sg.markers.toLocaleString('en-US')} markers, score ${sg.score.toFixed(0)}`
      + (sg.lrrZ !== undefined ? `, intensity z ${sg.lrrZ.toFixed(2)}` : '')],
  }
}

/** A row of the list of changes as a mark; a whole-chromosome row is located by its chromosome. */
export function defectMark(d: Defect, segments: readonly Segment[], loaded: Role | null, endBp: number | undefined): Mark | null {
  if (d.chrom === 'genome') return null
  const lines = [`${d.locus}${d.interval ? ` ${d.interval}` : ''}`,
    ...(d.informative !== undefined ? [`${d.informative.toLocaleString('en-US')} markers`] : []), d.why]
  const origin = d.origin
  const sg = segments.find((s) => s.chrom === d.chrom && segmentCoords(s).start === d.startBp && segmentCoords(s).end === d.endBp)
  if (sg) return segmentMark(sg, loaded, [headline(d), ...lines], origin)
  const whole = d.locus === `chr${d.chrom}`
  const end = whole ? endBp ?? d.endBp : d.endBp
  const named = origin === 'paternal' || origin === 'maternal' ? origin : null
  const gain = d.kind === 'copy-gain' || d.kind === 'trisomy' || d.kind === 'segmental-duplication'
  const what = gain ? 'there is an extra copy' : d.kind === 'copy-loss' || d.kind === 'monosomy' ? 'one copy is lost'
    : d.kind === 'parental-absence' ? "one parent's copy is missing" : headline(d)
  return {
    chrom: d.chrom, startBp: whole ? 0 : d.startBp, endBp: end, whole, gain, parent: named,
    sentence: whole ? `Whole chromosome ${d.chrom} (${(end / 1e6).toFixed(0)} Mb): ${what}${named ? `, ${WHOSE[named]}` : ', parent not named'}.`
      : `${where(d.chrom, d.startBp, d.endBp)}: ${what}${named ? `, ${WHOSE[named]}` : ''}.`,
    lines: [headline(d), ...lines],
  }
}

/** The whole-cell figure with its heading, for a result card. */
export function CellFigure({ data, marks, name }: { data: PlotData | undefined; marks: readonly Mark[]; name: string }): ReactNode {
  if (!data) return null
  return (
    <div style={{ marginTop: 10 }}>
      <div style={{ fontSize: 12, fontWeight: 600, textTransform: 'uppercase', letterSpacing: '0.04em', color: DIM }}>
        Every chromosome
      </div>
      <WholeCell data={data} marks={marks} name={name} />
    </div>
  )
}

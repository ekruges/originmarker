/**
 * A figure the page draws, as a vector PDF.
 *
 * Reads the SVG markup CallPlot emits (svg, g, rect, line, circle, text, tspan) and draws each
 * element with the report's PDF builder: fills and strokes as path operators, text as Helvetica
 * text runs, so it stays selectable. Opacity is folded onto the white ground the figure is drawn
 * on, so the file holds no transparency group and no image. An element outside that list is not
 * drawn. The file's Title is the figure's own title, and a footer line under the figure carries
 * the version.
 */
import { Pdf, stringWidth, type FontName } from './pdf.ts'

/** The version line every figure file prints. The page sets it once the server has answered. */
export const figureFooter = { tool: 'OriginMarker, version not reported' }
// ponytail: one module-level value instead of threading the version through every figure prop.

const FOOT = 16

export interface FigurePage { svg: string; title: string }

const attrsOf = (s: string): Record<string, string> => {
  const out: Record<string, string> = {}
  for (const m of s.matchAll(/([\w:-]+)="([^"]*)"/g)) out[m[1]] = m[2]
  return out
}

const decode = (t: string): string => t.replace(/&(#x[\da-f]+|#\d+|amp|lt|gt|quot|apos);/gi, (_, e: string) => (
  e[0] === '#' ? String.fromCodePoint(e[1] === 'x' || e[1] === 'X' ? parseInt(e.slice(2), 16) : +e.slice(1))
    : ({ amp: '&', lt: '<', gt: '>', quot: '"', apos: "'" } as Record<string, string>)[e.toLowerCase()]))

/** A colour as [r, g, b] in 0..1, or null for none. */
function rgbOf(c: string | undefined, vars: (name: string) => string): [number, number, number] | null {
  if (!c || c === 'none' || c === 'transparent') return null
  const v = /^var\((--[\w-]+)\)$/.exec(c.trim())
  if (v) return rgbOf(vars(v[1]) || '#888888', vars)
  if (c === 'white') return [1, 1, 1]
  if (c === 'black') return [0, 0, 0]
  const rgb = /^rgba?\(([^)]+)\)$/.exec(c.trim())
  if (rgb) { const [r, g, b] = rgb[1].split(',').map((x) => +x.trim() / 255); return [r, g, b] }
  let h = c.trim().replace('#', '')
  if (h.length === 3) h = h.replace(/./g, (x) => x + x)
  const n = parseInt(h.slice(0, 6), 16)
  if (!Number.isFinite(n)) return [0.53, 0.53, 0.53]
  return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]
}

/** A colour seen through an opacity on white, as hex for the PDF builder. */
function hexOn(c: [number, number, number], a: number): string {
  return `#${c.map((x) => Math.round(255 * (1 - a * (1 - x))).toString(16).padStart(2, '0')).join('')}`
}

interface Ctx { dx: number; dy: number; alpha: number; fill?: string; size: number; weight: string; anchor: string }

/** Draw one SVG onto a fresh page of the builder. Returns the page size. */
function drawSvg(pdf: Pdf, svg: string, vars: (name: string) => string, footer: string): [number, number] {
  const root = attrsOf(/<svg\b([^>]*)>/.exec(svg)?.[1] ?? '')
  const W = +(root.width ?? 0)
  const H = +(root.height ?? 0)
  const PH = H + FOOT
  const Y = (y: number) => PH - y
  const stack: Ctx[] = [{ dx: 0, dy: 0, alpha: 1, size: +(root['font-size'] ?? 11), weight: '400', anchor: 'start' }]
  let text: { x: number; y: number; ctx: Ctx; runs: { s: string; bold: boolean }[] } | null = null
  let tspan: Ctx | null = null
  const paint = (a: Record<string, string>, ctx: Ctx, key: 'fill' | 'stroke'): string | null => {
    const c = rgbOf(a[key] ?? (key === 'fill' ? ctx.fill ?? '#000000' : undefined), vars)
    if (!c) return null
    const op = ctx.alpha * +(a.opacity ?? 1) * +(a[`${key}-opacity`] ?? 1)
    return op <= 0 ? null : hexOn(c, op)
  }
  const stroke = (a: Record<string, string>, ctx: Ctx): boolean => {
    const s = paint(a, ctx, 'stroke')
    if (!s) return false
    pdf.setStrokeColor(s)
    pdf.setLineWidth(+(a['stroke-width'] ?? 1))
    const d = (a['stroke-dasharray'] ?? '').split(/[ ,]+/).filter(Boolean).map(Number)
    pdf.setDash(d[0] ?? 0, d[1] ?? d[0] ?? 0)
    return true
  }
  for (const m of svg.matchAll(/<(\/?)([\w:]+)([^>]*?)(\/?)>|([^<]+)/g)) {
    const [, close, tag, rawAttrs, selfClose, chars] = m
    const ctx = stack[stack.length - 1]
    if (chars !== undefined) {
      if (text) text.runs.push({ s: decode(chars), bold: +(tspan ?? text.ctx).weight >= 600 })
      continue
    }
    if (close) {
      if (tag === 'g' || tag === 'svg') stack.pop()
      else if (tag === 'tspan') tspan = null
      else if (tag === 'text' && text) {
        const t = text
        const font = (b: boolean): FontName => (b ? 'Helvetica-Bold' : 'Helvetica')
        const width = t.runs.reduce((w, r) => w + stringWidth(r.s, font(r.bold), t.ctx.size), 0)
        let x = t.x - (t.ctx.anchor === 'middle' ? width / 2 : t.ctx.anchor === 'end' ? width : 0)
        const fill = paint({}, t.ctx, 'fill')
        if (fill) {
          pdf.setFillColor(fill)
          for (const r of t.runs) {
            if (!r.s) continue
            pdf.setFont(font(r.bold), t.ctx.size)
            pdf.drawString(x, Y(t.y), r.s)
            x += stringWidth(r.s, font(r.bold), t.ctx.size)
          }
        }
        text = null
      }
      continue
    }
    const a = attrsOf(rawAttrs)
    const inherit = (extra: Partial<Ctx> = {}): Ctx => ({
      ...ctx, alpha: ctx.alpha * +(a.opacity ?? 1), fill: a.fill ?? ctx.fill,
      size: a['font-size'] ? +a['font-size'] : ctx.size, weight: a['font-weight'] ?? ctx.weight,
      anchor: a['text-anchor'] ?? ctx.anchor, ...extra,
    })
    if (tag === 'svg') {
      // The root draws at the origin; a nested svg is placed by its own x and y.
      const nested = stack.length > 1 || m.index! > 0
      stack.push(inherit(nested ? { dx: ctx.dx + +(a.x ?? 0), dy: ctx.dy + +(a.y ?? 0) } : {}))
      if (selfClose) stack.pop()
      continue
    }
    if (tag === 'g') { stack.push(inherit()); if (selfClose) stack.pop(); continue }
    if (tag === 'tspan') { tspan = { ...(text?.ctx ?? ctx), weight: a['font-weight'] ?? (text?.ctx ?? ctx).weight }; continue }
    if (tag === 'text') {
      text = { x: ctx.dx + +(a.x ?? 0), y: ctx.dy + +(a.y ?? 0), ctx: inherit({ alpha: ctx.alpha * +(a.opacity ?? 1) }), runs: [] }
      if (selfClose) text = null
      continue
    }
    if (tag === 'rect') {
      const x = ctx.dx + +(a.x ?? 0); const y = ctx.dy + +(a.y ?? 0)
      const w = +(a.width ?? 0); const h = +(a.height ?? 0)
      const fill = paint(a, { ...ctx, fill: a.fill ?? ctx.fill }, 'fill')
      if (fill) { pdf.setFillColor(fill); pdf.rect(x, Y(y + h), w, h, true) }
      if (stroke(a, ctx)) pdf.rect(x, Y(y + h), w, h)
    } else if (tag === 'line') {
      if (stroke({ stroke: '#000000', ...a }, ctx)) {
        pdf.line(ctx.dx + +(a.x1 ?? 0), Y(ctx.dy + +(a.y1 ?? 0)), ctx.dx + +(a.x2 ?? 0), Y(ctx.dy + +(a.y2 ?? 0)))
      }
    } else if (tag === 'circle') {
      const cx = ctx.dx + +(a.cx ?? 0); const cy = Y(ctx.dy + +(a.cy ?? 0)); const r = +(a.r ?? 0)
      const fill = paint(a, { ...ctx, fill: a.fill ?? ctx.fill }, 'fill')
      if (fill) { pdf.setFillColor(fill); pdf.circle(cx, cy, r, true) }
      if (+(a['stroke-width'] ?? 1) > 0 && stroke(a, ctx)) pdf.circle(cx, cy, r)
    }
  }
  pdf.setDash()
  pdf.setFillColor('#666666')
  pdf.setFont('Helvetica', 6.5)
  pdf.drawString(2, 5, footer)
  return [Math.max(W, 1), PH]
}

/** One file of one or more figures, one page each, titled with the first figure's title. */
export function figurePdf(pages: readonly FigurePage[], opts: { title?: string; vars?: (name: string) => string; footer?: string } = {}): Blob {
  const pdf = new Pdf()
  pdf.setTitle(opts.title ?? pages[0]?.title ?? 'figure')
  for (const p of pages) {
    const size = drawSvg(pdf, p.svg, opts.vars ?? (() => ''), opts.footer ?? figureFooter.tool)
    pdf.showPage(size)
  }
  return pdf.save()
}

/** The page's colour variables, read from the document, for a figure taken off the page. */
export const pageVars = (name: string): string =>
  (typeof document === 'undefined' ? '' : getComputedStyle(document.documentElement).getPropertyValue(name).trim())

/** Hand the reader a PDF of these figures. */
export function savePdf(pages: readonly FigurePage[], name: string, title?: string): void {
  const a = document.createElement('a')
  a.href = URL.createObjectURL(figurePdf(pages, { title, vars: pageVars }))
  a.download = `${name}.pdf`
  a.click()
  setTimeout(() => URL.revokeObjectURL(a.href), 1000)
}

/**
 * What a PDF file holds, read from its bytes: whether every xref offset points at its object, how
 * many image objects or inline images it carries, its Title, its page count and every string it
 * draws. Used by the checks and the gate on every file this module and the report write.
 */
export function pdfFacts(txt: string): { xrefOk: boolean; images: number; title: string; pages: number; strings: string[] } {
  const unescape = (x: string) => x.replace(/\\([0-7]{3})|\\(.)/g, (_, o, c) => (o ? String.fromCharCode(parseInt(o, 8)) : c))
  const start = Number(/startxref\n(\d+)/.exec(txt)?.[1] ?? NaN)
  const xref = txt.slice(start).split('\n')
  const n = Number(/^0 (\d+)$/.exec(xref[1] ?? '')?.[1] ?? 0)
  let xrefOk = txt.startsWith('%PDF-') && txt.trimEnd().endsWith('%%EOF') && xref[0] === 'xref' && n > 1
  for (let i = 1; i < n && xrefOk; i += 1) {
    const off = Number(xref[2 + i].slice(0, 10))
    xrefOk = txt.startsWith(`${i} 0 obj`, off)
  }
  return {
    xrefOk,
    images: (txt.match(/\/Subtype\s*\/Image|\/XObject|\nBI\s|\sBI\s.*\sID\s/g) ?? []).length,
    title: unescape(/\/Title \(((?:\\.|[^\\)])*)\)/.exec(txt)?.[1] ?? ''),
    pages: (txt.match(/\/Type \/Page /g) ?? []).length,
    strings: [...txt.matchAll(/\(((?:\\.|[^\\)])*)\) Tj/g)].map((m) => unescape(m[1])),
  }
}

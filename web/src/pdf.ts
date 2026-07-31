/**
 * A PDF canvas, in the browser.
 *
 * Syngamy reads a family's arrays in this tab and has no endpoint to send them to, so the report
 * has to be typeset here rather than by the reportlab code that renders the 1.x panel. The API
 * below is deliberately the same shape as reportlab's canvas so the two layouts stay comparable
 * line for line, and both emit the same page furniture from the same numbers.
 *
 * Only the base-14 fonts are used, which is why nothing is embedded: their metrics are part of
 * the format. Every byte written is Latin-1, so string length is byte length and the xref offsets
 * are simply running character counts.
 */

export const LETTER: [number, number] = [612, 792]

export type FontName = 'Helvetica' | 'Helvetica-Bold' | 'Helvetica-Oblique' | 'Courier'
  | 'Courier-Bold'

const FONT_KEY: Record<FontName, string> = {
  'Helvetica': 'F1',
  'Helvetica-Bold': 'F2',
  'Helvetica-Oblique': 'F3',
  'Courier': 'F4',
  'Courier-Bold': 'F5',
}

/** Adobe's published widths for the base-14 fonts, per 1000 units, for characters 32 to 126. */
const W_HELV = [
  278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278,
  556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556,
  1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778,
  667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556,
  333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556,
  556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584,
]
const W_BOLD = [
  278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278,
  556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584, 611,
  975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722, 778,
  667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556,
  333, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611, 611,
  611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584,
]

/** Width of one string at one size, in points. The measure every wrap and right-edge uses. */
export function stringWidth(s: string, font: FontName, size: number): number {
  if (font === 'Courier' || font === 'Courier-Bold') return s.length * 0.6 * size
  const w = font === 'Helvetica-Bold' ? W_BOLD : W_HELV
  let total = 0
  for (let i = 0; i < s.length; i += 1) {
    const c = s.charCodeAt(i)
    // Anything outside the measured range is charged as a lowercase n, which is the safe way
    // to be wrong: the line is never narrower than it measures, so text cannot run off the page.
    total += c >= 32 && c <= 126 ? w[c - 32] : w[78]
  }
  return (total * size) / 1000
}

/** Greedy wrap to a width, breaking inside an over-long word rather than overflowing. */
export function wrap(s: string, font: FontName, size: number, width: number): string[] {
  const out: string[] = []
  for (const para of String(s).split('\n')) {
    let line = ''
    for (const word of para.split(/\s+/)) {
      if (!word) continue
      const next = line ? `${line} ${word}` : word
      if (stringWidth(next, font, size) <= width) { line = next; continue }
      if (line) { out.push(line); line = '' }
      let rest = word
      while (stringWidth(rest, font, size) > width) {
        let cut = rest.length
        while (cut > 1 && stringWidth(rest.slice(0, cut), font, size) > width) cut -= 1
        out.push(rest.slice(0, cut))
        rest = rest.slice(cut)
      }
      line = rest
    }
    out.push(line)
  }
  return out.length ? out : ['']
}

/** Escape for a PDF literal string, and fold anything Latin-1 cannot carry onto ASCII. */
function lit(s: string): string {
  let out = '('
  for (const ch of s) {
    const c = ch.codePointAt(0) ?? 63
    if (ch === '(' || ch === ')' || ch === '\\') out += `\\${ch}`
    else if (c < 32) out += ' '
    else if (c < 127) out += ch
    else if (c === 0x2013 || c === 0x2014 || c === 0x2212) out += '-'
    else if (c === 0x2018 || c === 0x2019) out += "'"
    else if (c === 0x201c || c === 0x201d) out += '"'
    else if (c === 0x2265) out += '>='
    else if (c === 0x2264) out += '<='
    else if (c === 0x00d7 || c === 0x2715) out += 'x'
    else if (c < 256) out += `\\${c.toString(8).padStart(3, '0')}`
    else out += '?'
  }
  return `${out})`
}

const f2 = (n: number): string => (Math.round(n * 100) / 100).toString()

const rgb = (hex: string): string => {
  const h = hex.replace('#', '')
  const n = parseInt(h.length === 3 ? h.replace(/./g, (c) => c + c) : h, 16)
  // eslint-disable-next-line no-bitwise
  return `${f2(((n >> 16) & 255) / 255)} ${f2(((n >> 8) & 255) / 255)} ${f2((n & 255) / 255)}`
}

export class Pdf {
  readonly width: number
  readonly height: number
  private pages: string[] = []
  private buf: string[] = []
  private font: FontName = 'Helvetica'
  private size = 10
  private title = 'report'

  constructor([w, h]: [number, number] = LETTER) {
    this.width = w
    this.height = h
  }

  setTitle(t: string): void { this.title = t }

  setFont(name: FontName, size: number): void {
    this.font = name
    this.size = size
  }

  setFillColor(hex: string): void { this.buf.push(`${rgb(hex)} rg`) }
  setStrokeColor(hex: string): void { this.buf.push(`${rgb(hex)} RG`) }
  setLineWidth(w: number): void { this.buf.push(`${f2(w)} w`) }

  drawString(x: number, y: number, s: string): void {
    this.buf.push(`BT /${FONT_KEY[this.font]} ${f2(this.size)} Tf ${f2(x)} ${f2(y)} Td `
      + `${lit(s)} Tj ET`)
  }

  drawRightString(x: number, y: number, s: string): void {
    this.drawString(x - stringWidth(s, this.font, this.size), y, s)
  }

  drawCentredString(x: number, y: number, s: string): void {
    this.drawString(x - stringWidth(s, this.font, this.size) / 2, y, s)
  }

  rect(x: number, y: number, w: number, h: number, fill = false): void {
    this.buf.push(`${f2(x)} ${f2(y)} ${f2(w)} ${f2(h)} re ${fill ? 'f' : 'S'}`)
  }

  line(x1: number, y1: number, x2: number, y2: number): void {
    this.buf.push(`${f2(x1)} ${f2(y1)} m ${f2(x2)} ${f2(y2)} l S`)
  }

  circle(cx: number, cy: number, r: number, fill = false): void {
    const k = r * 0.5522847498
    this.buf.push(
      `${f2(cx + r)} ${f2(cy)} m`,
      `${f2(cx + r)} ${f2(cy + k)} ${f2(cx + k)} ${f2(cy + r)} ${f2(cx)} ${f2(cy + r)} c`,
      `${f2(cx - k)} ${f2(cy + r)} ${f2(cx - r)} ${f2(cy + k)} ${f2(cx - r)} ${f2(cy)} c`,
      `${f2(cx - r)} ${f2(cy - k)} ${f2(cx - k)} ${f2(cy - r)} ${f2(cx)} ${f2(cy - r)} c`,
      `${f2(cx + k)} ${f2(cy - r)} ${f2(cx + r)} ${f2(cy - k)} ${f2(cx + r)} ${f2(cy)} c`,
      fill ? 'f' : 'S',
    )
  }

  /** One or more closed polygons, filled as a single path. */
  polys(list: [number, number][][], fill = true): void {
    for (const poly of list) {
      poly.forEach(([x, y], i) => this.buf.push(`${f2(x)} ${f2(y)} ${i ? 'l' : 'm'}`))
      this.buf.push('h')
    }
    this.buf.push(fill ? 'f' : 'S')
  }

  showPage(): void {
    this.pages.push(this.buf.join('\n'))
    this.buf = []
  }

  get pageCount(): number { return this.pages.length }

  save(): Blob {
    if (this.buf.length) this.showPage()
    const objs: string[] = []
    const push = (body: string): number => { objs.push(body); return objs.length }

    const catalog = push('')                                   // 1, filled once Pages is known
    const pagesObj = push('')                                  // 2, same
    const fonts = (Object.keys(FONT_KEY) as FontName[]).map((name) =>
      push(`<< /Type /Font /Subtype /Type1 /BaseFont /${name} /Encoding /WinAnsiEncoding >>`))
    const res = `<< /Font << ${(Object.keys(FONT_KEY) as FontName[])
      .map((n, i) => `/${FONT_KEY[n]} ${fonts[i]} 0 R`).join(' ')} >> >>`

    const kids: number[] = []
    for (const content of this.pages) {
      const stream = push(`<< /Length ${content.length} >>\nstream\n${content}\nendstream`)
      kids.push(push(`<< /Type /Page /Parent ${pagesObj} 0 R /MediaBox `
        + `[0 0 ${this.width} ${this.height}] /Resources ${res} /Contents ${stream} 0 R >>`))
    }
    objs[catalog - 1] = `<< /Type /Catalog /Pages ${pagesObj} 0 R >>`
    objs[pagesObj - 1] = `<< /Type /Pages /Kids [${kids.map((k) => `${k} 0 R`).join(' ')}] `
      + `/Count ${kids.length} >>`
    const info = push(`<< /Title ${lit(this.title)} /Producer (Syngamy) >>`)

    let out = '%PDF-1.4\n'
    const offsets: number[] = []
    objs.forEach((body, i) => {
      offsets.push(out.length)
      out += `${i + 1} 0 obj\n${body}\nendobj\n`
    })
    const startxref = out.length
    out += `xref\n0 ${objs.length + 1}\n0000000000 65535 f \n`
    for (const off of offsets) out += `${String(off).padStart(10, '0')} 00000 n \n`
    out += `trailer\n<< /Size ${objs.length + 1} /Root ${catalog} 0 R /Info ${info} 0 R >>\n`
      + `startxref\n${startxref}\n%%EOF\n`

    const bytes = new Uint8Array(out.length)
    for (let i = 0; i < out.length; i += 1) bytes[i] = out.charCodeAt(i) & 0xff
    return new Blob([bytes], { type: 'application/pdf' })
  }
}

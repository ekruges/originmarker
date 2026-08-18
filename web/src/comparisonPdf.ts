/**
 * The feature comparison as its own report, and as a section foldable into the main one.
 *
 * ONE DRAWING FUNCTION, TWO DESTINATIONS. `drawComparison` renders into whatever page it is given,
 * so the standalone report and the section bundled into the main report are the same marks from the
 * same numbers. A second implementation for the bundled version is how the two come to disagree,
 * and a reader who notices that has to distrust both.
 *
 * WHAT THE CHARTS ARE FOR. The reader's question is "do these regions overlap these features or
 * not", so there are three plain answers and nothing else: the observation against its own null,
 * the same thing as a fold so features can be ranked, and a grid of which region touched which
 * feature. The grid exists because an enrichment can clear its threshold on two regions out of
 * twenty, and that is invisible in a p value.
 *
 * THE NULL IS DRAWN, NOT SUMMARISED. A permutation p is uninterpretable without the distribution it
 * came from: a ratio near 1 reaches a small p when the null is tight, and a reader shown only the p
 * would call that a finding.
 */
import { LETTER, Pdf, wrap, type FontName } from './pdf.ts'
import type { Figure } from './figures.ts'
import { FIG } from './figures.ts'
import { comparisonFigures } from './comparisonFigures.ts'
import type { ComparisonResult } from './comparison.ts'

const INK = '#1a1a1a'
const GREY = '#6b6b6b'
const LINE = '#cfcfcf'
const FIG_CAPTION = FIG.font.caption

export interface Cursor { y: number }

/**
 * Render a Figure into a PDF page. Four primitives, exactly as the screen renderer draws them.
 *
 * The figure's own coordinates have y increasing DOWNWARD, as SVG does, because that is how the
 * layout was written. PDF has y increasing upward, so the flip happens here, once, rather than
 * being threaded through every layout calculation where it would eventually be got wrong.
 */
export function drawFigure(
  pdf: Pdf, fig: Figure, left: number, topY: number,
): number {
  const Y = (y: number) => topY - y
  for (const m of fig.marks) {
    if (m.k === 'line') {
      pdf.setStrokeColor(m.colour); pdf.setLineWidth(m.width)
      pdf.line(left + m.x1, Y(m.y1), left + m.x2, Y(m.y2))
    } else if (m.k === 'rect') {
      pdf.setFillColor(m.colour)
      pdf.rect(left + m.x, Y(m.y + m.h), m.w, m.h, true)
    } else if (m.k === 'dot') {
      pdf.setFillColor(m.colour)
      pdf.circle(left + m.x, Y(m.y), m.r, true)
    } else {
      pdf.setFont(m.bold ? 'Helvetica-Bold' : m.italic ? 'Helvetica-Oblique' : 'Helvetica', m.size)
      pdf.setFillColor(m.colour)
      if (m.anchor === 'end') pdf.drawRightString(left + m.x, Y(m.y), m.s)
      else if (m.anchor === 'middle') pdf.drawCentredString(left + m.x, Y(m.y), m.s)
      else pdf.drawString(left + m.x, Y(m.y), m.s)
    }
  }
  return topY - fig.h
}


/**
 * Draw the whole comparison into a page the caller owns.
 *
 * Returns the new y. The caller decides pagination, because the standalone report starts on a fresh
 * page while the bundled section continues one already in progress.
 */
export function drawComparison(
  pdf: Pdf,
  c: ComparisonResult,
  regionNames: readonly string[],
  opts: { left?: number; width?: number; y: number; newPage: () => number },
): number {
  const L = opts.left ?? 54
  const W = opts.width ?? 504
  let y = opts.y

  const need = (h: number) => { if (y - h < 60) y = opts.newPage() }
  const text = (
    s: string, size = 8, font: FontName = 'Helvetica', gap = 3, colour = INK, x = L,
  ) => {
    pdf.setFont(font, size); pdf.setFillColor(colour)
    for (const ln of wrap(s, font, size, W - (x - L))) {
      need(size + gap); pdf.drawString(x, y - size, ln); y -= size + gap
    }
    y -= 2
  }
  const heading = (s: string, size = 10) => {
    need(size + 12)
    pdf.setFont('Helvetica-Bold', size); pdf.setFillColor(INK)
    pdf.drawString(L, y - size, s); y -= size + 7
  }

  heading('Where these regions sit in the genome', 12)
  text(c.headline, 9, 'Helvetica-Bold', 3.4)
  text(`${c.regions} region${c.regions === 1 ? '' : 's'} tested against ${c.features.length} `
    + `feature set${c.features.length === 1 ? '' : 's'}, `
    + `${c.permutations.toLocaleString()} permutations.`, 8, 'Helvetica', 3, GREY)

  // THE FIGURES, laid out once in comparisonFigures.ts and drawn here by the same four primitives
  // the screen uses. Two implementations of one picture drift, and a reader who notices the drift
  // has to distrust both, so there is only one layout.
  for (const fig of comparisonFigures(c, regionNames)) {
    need(fig.h + 30)
    // A figure is never split across a page: half a chart is worse than a page break before it.
    y = drawFigure(pdf, fig, L, y)
    y -= 4
    text(fig.caption, FIG_CAPTION, 'Helvetica-Oblique', 2.6, GREY)
    y -= 6
  }

  // ---- what each feature would mean --------------------------------------------------------------
  heading('What a coincidence with each feature would mean', 9.4)
  for (const f of c.features) {
    text(`${f.label}: ${f.means}.`, 7.6, 'Helvetica', 2.6, GREY)
  }

  // ---- generated methods -------------------------------------------------------------------------
  heading('Methods', 9.4)
  text(c.methods, 7.6, 'Helvetica', 2.8, INK)

  // ---- the caveat, last and unavoidable ----------------------------------------------------------
  need(34)
  pdf.setStrokeColor(LINE); pdf.setLineWidth(2)
  const capTop = y
  text(c.caveat, 8, 'Helvetica-Bold', 3, INK, L + 10)
  pdf.line(L + 3, capTop - 2, L + 3, y + 4)
  return y
}

/** The standalone report, for a reader who ran the comparison and wants only that. */
export function comparisonPdf(
  c: ComparisonResult, regionNames: readonly string[], sampleName: string, build: string,
): Blob {
  const pdf = new Pdf(LETTER)
  pdf.setTitle(`OriginMarker feature comparison: ${sampleName}`)
  const H = LETTER[1]
  const newPage = () => { pdf.showPage(); return H - 54 }
  let y = H - 54

  pdf.setFont('Helvetica-Bold', 16); pdf.setFillColor(INK)
  pdf.drawString(54, y - 16, 'Feature comparison')
  y -= 24
  pdf.setFont('Helvetica', 9); pdf.setFillColor(GREY)
  pdf.drawString(54, y - 9, `${sampleName}   ${build}`)
  y -= 20
  pdf.setStrokeColor(LINE); pdf.setLineWidth(0.6); pdf.line(54, y, 558, y)
  y -= 14

  y = drawComparison(pdf, c, regionNames, { y, newPage })

  pdf.showPage()
  return pdf.save()
}

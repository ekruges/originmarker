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
import {
  enrichmentBars, foldChart, nullHistograms, type ComparisonResult,
} from './comparison.ts'

const INK = '#1a1a1a'
const GREY = '#6b6b6b'
const HIT = '#a4243b'
const LINE = '#cfcfcf'

export interface Cursor { y: number }

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

  // ---- chart 1: observed against its own null, one shared axis ----------------------------------
  const bars = enrichmentBars(c)
  if (bars.length) {
    heading('Overlap against a matched null', 9.4)
    need(bars.length * 18 + 26)
    const labelW = 168
    const plotW = 210
    const x0 = L + labelW
    const scale = (v: number) => x0 + (v / bars[0].axisMax) * plotW
    for (const b of bars) {
      need(18)
      const mid = y - 6
      pdf.setFont('Helvetica', 7.4); pdf.setFillColor(GREY)
      pdf.drawString(L, mid - 2, b.label.slice(0, 40))
      // the axis
      pdf.setStrokeColor(LINE); pdf.setLineWidth(0.5)
      pdf.line(x0, mid, x0 + plotW, mid)
      // the middle 95% of the null
      pdf.setStrokeColor(GREY); pdf.setLineWidth(3.4)
      pdf.line(scale(b.lo), mid, scale(b.hi), mid)
      // the observation
      pdf.setFillColor(b.significant ? HIT : INK)
      pdf.circle(scale(b.observed), mid, 2.6, true)
      pdf.setFont('Helvetica', 7); pdf.setFillColor(b.significant ? HIT : GREY)
      pdf.drawString(x0 + plotW + 6, mid - 2,
        `${(100 * b.observed).toFixed(0)}% vs ${(100 * b.expected).toFixed(0)}%`
        + (b.significant ? `  p=${b.p.toExponential(1)}` : ''))
      y -= 18
    }
    text('Bar is the middle 95% of the matched null, dot is what was observed. One axis across '
      + 'every row: a per-row axis would make each feature fill its own width and look equally '
      + 'enriched.', 7, 'Helvetica-Oblique', 2.6, GREY)
  }

  // ---- chart 2: the permutation null itself, with the observation in it --------------------------
  //
  // The most defensible figure here and it costs nothing, because the enrichment already computed
  // the distribution. A p value is a summary of this picture, and summaries of it mislead in both
  // directions: a ratio near one reaches a small p when the null is tight, a large ratio reaches
  // nothing when it is broad. Drawing it lets a reader judge instead of taking the p on trust.
  const hists = nullHistograms(c)
  if (hists.length) {
    heading('The permutation null, with the observation marked', 9.4)
    const cellW = 118
    const cellH = 34
    const perRow = Math.max(1, Math.floor(W / cellW))
    for (let i = 0; i < hists.length; i += perRow) {
      const row = hists.slice(i, i + perRow)
      need(cellH + 22)
      row.forEach((h, j) => {
        const x0 = L + j * cellW
        const top = y - 10
        pdf.setFont('Helvetica', 6.4); pdf.setFillColor(GREY)
        pdf.drawString(x0, y - 4, h.label.slice(0, 26))
        const max = Math.max(1, ...h.bins)
        const bw = (cellW - 14) / h.bins.length
        h.bins.forEach((v, k) => {
          pdf.setFillColor(k === h.observedBin ? HIT : LINE)
          const bh = (v / max) * cellH
          if (bh > 0.2) pdf.rect(x0 + k * bw, top - cellH + (cellH - bh), Math.max(0.6, bw - 0.4), bh, true)
        })
        pdf.setStrokeColor(h.significant ? HIT : INK); pdf.setLineWidth(0.9)
        pdf.line(x0 + h.observedAt * (cellW - 14), top - cellH, x0 + h.observedAt * (cellW - 14), top)
        pdf.setFont('Helvetica', 6); pdf.setFillColor(GREY)
        pdf.drawString(x0, top - cellH - 7,
          h.observedBin === -1 ? 'observed OUTSIDE the null' : `p=${h.p.toExponential(1)}`)
      })
      y -= cellH + 22
    }
    text('Bars are the matched null the permutation actually drew; the line is the observation. An '
      + 'observation outside the null entirely is the strongest result this analysis produces and '
      + 'is labelled as such.', 7, 'Helvetica-Oblique', 2.6, GREY)
  }

  // ---- chart 3: fold, so features can be ranked on one scale ------------------------------------
  const folds = foldChart(c)
  if (folds.length) {
    heading('Fold enrichment', 9.4)
    need(folds.length * 15 + 24)
    const labelW = 168
    const plotW = 210
    const x0 = L + labelW
    const max = Math.max(2, ...folds.map((f) => f.fold))
    // the null, at 1.00x
    const nullX = x0 + (1 / max) * plotW
    for (const f of folds) {
      need(15)
      const mid = y - 5
      pdf.setFont('Helvetica', 7.4); pdf.setFillColor(GREY)
      pdf.drawString(L, mid - 2, f.label.slice(0, 40))
      pdf.setFillColor(f.fold >= 1 ? INK : GREY)
      pdf.rect(x0, mid - 3.4, Math.max(0.6, (f.fold / max) * plotW), 6.8, true)
      pdf.setFont('Helvetica', 7); pdf.setFillColor(GREY)
      pdf.drawString(x0 + plotW + 6, mid - 2, `${f.fold.toFixed(2)}x`)
      y -= 15
    }
    pdf.setStrokeColor(LINE); pdf.setLineWidth(0.5)
    pdf.line(nullX, y + folds.length * 15 - 2, nullX, y + 2)
    text('The vertical line is the null at 1.00x. Left of it is LESS overlap than the matched null '
      + 'expects.', 7, 'Helvetica-Oblique', 2.6, GREY)
  }

  // ---- chart 4: the grid, which is the plain answer ---------------------------------------------
  if (regionNames.length && c.features.length) {
    heading('Which region touches which feature', 9.4)
    const colW = Math.min(64, (W - 150) / Math.max(1, c.features.length))
    need(regionNames.length * 11 + 40)
    // header
    pdf.setFont('Helvetica', 6.2); pdf.setFillColor(GREY)
    c.features.forEach((f, i) => {
      pdf.drawString(L + 150 + i * colW, y - 6, f.label.slice(0, 11))
    })
    y -= 13
    for (const n of regionNames) {
      need(11)
      pdf.setFont('Courier', 7); pdf.setFillColor(INK)
      pdf.drawString(L, y - 6, n.slice(0, 26))
      c.features.forEach((f, i) => {
        const touched = f.hits.includes(n)
        pdf.setFillColor(touched ? HIT : LINE)
        pdf.circle(L + 150 + i * colW + 8, y - 4, touched ? 2.4 : 1, true)
      })
      y -= 11
    }
    y -= 3
    text('A filled dot is an overlap. This is where an enrichment that rests on two regions out of '
      + 'twenty becomes visible, which a p value hides.', 7, 'Helvetica-Oblique', 2.6, GREY)
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

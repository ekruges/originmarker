import { Pdf, LETTER, wrap, type FontName } from './pdf.ts'
import { DISCLAIMER, parseMark, type MarkGeometry } from './syngamyPdf.ts'
import { pct } from './parentage.ts'
import type { RunProvenance, PairRate, SampleResult } from './inferredExports.ts'
import type { Limit } from './InferredLimits.ts'

/**
 * The report for a run whose parental reference was reconstructed rather than measured.
 *
 * Same document as the Syngamy report: Letter, the same margins, the same base-14 fonts, the
 * same masthead, the same disclaimer footer, typeset in the browser with no dependency and no
 * genotype leaving the page.
 *
 * One thing is deliberately different. Every page carries a line saying the reference was
 * inferred, not just page one, because a single page of a report is routinely photocopied,
 * emailed or pasted into a slide on its own, and a reader holding that page has to be able to
 * see that no array of the parent was ever used.
 */

const INK = '#1d2126'
const GREY = '#6b727b'
const WARN = '#a6572e'

export interface InferredReportInput {
  provenance: RunProvenance
  groups: string[][]
  samples: SampleResult[]
  pairs: PairRate[]
  matrixNames: string[]
  matrixRate: (a: number, b: number) => number
  ascertainment: Record<string, number>
  members: { name: string; absence: number; ceiling: number; ratio: number; verdict: string }[]
  controls: { name: string; role: string; absence: number; ratio: number; verdict: string }[]
  limits: Limit[]
  markSvg: string
}

async function mark(svg: string): Promise<MarkGeometry | null> {
  if (!svg) return null
  if (svg.trim().startsWith('<')) return parseMark(svg)
  try {
    const res = await fetch(svg)
    return res.ok ? parseMark(await res.text()) : null
  } catch { return null }
}

const ratio = (r: number): string => (Number.isFinite(r) ? `${r.toFixed(2)}x` : 'n/a')
const reads = (v: string): string =>
  ({ parent_genome_present: 'present', no_parental_contribution: 'absent' }[v] ?? 'uncalled')

export async function buildInferredPdf(input: InferredReportInput): Promise<Blob> {
  const [W, H] = LETTER
  const L = 40
  const R = W - 40
  const TOP = H - 42
  const BOT = 88
  const p = input.provenance
  const pdf = new Pdf(LETTER)
  pdf.setTitle('Progenitor')
  const logo = await mark(input.markSvg)

  let y = TOP
  let page = 0

  const footer = (): void => {
    // The inferred-reference line sits ABOVE the disclaimer and on every page. One page lifted
    // out of a report must still say the parent was never measured.
    pdf.setFont('Helvetica-Bold', 7)
    pdf.setFillColor(WARN)
    pdf.drawString(L, 54, 'Reference reconstructed from haploid products. No array of the '
      + 'parent was used.')
    pdf.setFont('Helvetica-Oblique', 7)
    pdf.setFillColor(INK)
    let yy = 44
    for (const line of wrap(DISCLAIMER, 'Helvetica-Oblique', 7, R - L)) {
      pdf.drawString(L, yy, line)
      yy -= 8.5
    }
    pdf.setFont('Helvetica', 6.5)
    pdf.setFillColor(GREY)
    pdf.drawString(L, 26, `Progenitor | ${p.tool} | reconstructed reference `
      + `| report ${p.reportId} `
      + `| generated ${p.generatedAt} | page ${page}`)
    pdf.setFillColor(INK)
  }

  const newPage = (): void => {
    if (page) { footer(); pdf.showPage() }
    page += 1
    y = TOP
  }
  const need = (h: number): void => { if (y - h < BOT) newPage() }

  const text = (
    s: string, size = 8, font: FontName = 'Helvetica', gap = 2, colour = INK, x0 = L,
  ): void => {
    for (const line of wrap(s, font, size, R - x0)) {
      need(size + gap)
      pdf.setFont(font, size)
      pdf.setFillColor(colour)
      pdf.drawString(x0, y - size, line)
      y -= size + gap
    }
    pdf.setFillColor(INK)
  }

  const heading = (s: string, size = 10): void => {
    need(size + 12)
    y -= 4
    text(s, size, 'Helvetica-Bold', 4)
  }

  const card = (rows: [string, string][], kw = 118): void => {
    for (const [k, v] of rows) {
      need(11)
      pdf.setFont('Helvetica-Bold', 8)
      pdf.setFillColor(INK)
      pdf.drawString(L, y - 8, k)
      pdf.setFont('Helvetica', 8)
      wrap(String(v), 'Helvetica', 8, R - L - kw).forEach((line, i) => {
        if (i) { need(10); pdf.setFont('Helvetica', 8) }
        pdf.drawString(L + kw, y - 8, line)
        y -= 10
      })
    }
  }

  interface Col { head: string; w: number; right?: boolean; font?: FontName; size?: number }

  const table = (cols: Col[], rows: (string | { v: string; colour: string })[][]): void => {
    const head = (): void => {
      need(14)
      pdf.setFont('Helvetica-Bold', 6.8)
      pdf.setFillColor(INK)
      let x = L
      for (const c of cols) {
        if (c.right) pdf.drawRightString(x + c.w - 6, y - 7, c.head)
        else pdf.drawString(x, y - 7, c.head)
        x += c.w
      }
      y -= 9
      pdf.setStrokeColor(INK)
      pdf.setLineWidth(0.4)
      pdf.line(L, y, R, y)
      y -= 2
    }
    head()
    for (const row of rows) {
      const lines = row.map((c, i) => wrap(
        typeof c === 'string' ? c : c.v, cols[i].font ?? 'Helvetica',
        cols[i].size ?? 6.8, cols[i].w - 6,
      ))
      const h = Math.max(...lines.map((l) => l.length)) * 8.4 + 0.6
      if (y - h < BOT) { newPage(); head() }
      let x = L
      row.forEach((c, i) => {
        pdf.setFont(cols[i].font ?? 'Helvetica', cols[i].size ?? 6.8)
        pdf.setFillColor(typeof c === 'string' ? INK : c.colour)
        lines[i].forEach((line, j) => {
          if (cols[i].right) pdf.drawRightString(x + cols[i].w - 6, y - 7 - j * 8.4, line)
          else pdf.drawString(x, y - 7 - j * 8.4, line)
        })
        x += cols[i].w
      })
      pdf.setFillColor(INK)
      y -= h
    }
  }

  newPage()

  // --- masthead, page 1 only ------------------------------------------------------------
  const MARK = 26
  if (logo) {
    const k = MARK / logo.view
    const at = (ux: number, uy: number): [number, number] => [L + ux * k, y - uy * k]
    const [cx, cy, r, sw] = logo.ring
    pdf.setStrokeColor(logo.ringColour)
    pdf.setLineWidth(sw * k)
    pdf.circle(...at(cx, cy), r * k)
    const [tx, ty, sx, sy] = logo.transform
    pdf.setFillColor(logo.fill)
    pdf.polys(logo.polys.map((poly: [number, number][]) =>
      poly.map(([gx, gy]: [number, number]) => at(tx + sx * gx, ty + sy * gy))))
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica-Bold', 14)
    pdf.drawString(L + MARK + 8, y - 18, 'Progenitor')
    y -= MARK + 4
  } else {
    text('Progenitor', 14, 'Helvetica-Bold', 4)
  }

  text('Parental genotype reconstructed from haploid products. '
    + `Experiment ${p.experiment}.`, 8.6)
  y -= 2

  // --- what this reference is, before any result -----------------------------------------
  heading('The reference was inferred, not measured')
  text('No array of this parent exists in this run. The genotype every rate below was measured '
    + `against was reconstructed from ${p.products.length} haploid meiotic products. At a marker `
    + 'where the parent is homozygous every product carries that allele, so agreement across '
    + 'products recovers it; where the parent is heterozygous each product carries one allele at '
    + 'random, and a site where they happen to agree enters the reference wrongly as homozygous.',
  7.4, 'Helvetica', 2.4, GREY)
  y -= 2
  card([
    ['Products', p.products.join(', ')],
    ['Agreement rule', `a marker enters the reference when at least ${p.mMin} products were `
      + 'called there and all of them agreed'],
    ['Markers', p.markers.toLocaleString()],
    ['Mean products', `${p.meanM.toFixed(2)} per retained marker`],
    ['Ascertainment', `${pct(p.ascertainment, 1)} of the parent's heterozygosity survives in the `
      + 'retained marker set'],
    ['Contamination', `${pct(p.contamination, 2)} of the reference is a heterozygous site read `
      + `as homozygous, which adds ${pct(p.spuriousAbsence, 2)} absence to a true haploid offspring`],
  ])

  // --- membership, which comes before the reference exists --------------------------------
  heading('Membership, established before the reference was built')
  text('Two haploid products of one parent differ only where that parent is heterozygous and the '
    + 'two drew differently. Measured across two experiments: 4.68% to 9.70% within one parent '
    + 'over 46 pairs, and 9.88% to 16.10% between parents over 45. Every pair inside a group '
    + 'reads same parent; a chain of pairwise links is not sufficient, because one genuine '
    + `cross-parent pair measured 9.88%, under the ${pct(p.sameParentMax, 1)} cut.`,
  7.4, 'Helvetica', 2.4, GREY)
  y -= 2
  table(
    [{ head: 'Group', w: 60 }, { head: 'Products', w: 300 },
      { head: 'Weakest pair inside', w: 112, right: true }],
    input.groups.map((g, i) => {
      const idx = g.map((n) => input.matrixNames.indexOf(n)).filter((x) => x >= 0)
      let worst = 0
      for (const a of idx) for (const b of idx) if (a < b) {
        worst = Math.max(worst, input.matrixRate(a, b))
      }
      return [`${i + 1}`, g.join(', '), idx.length > 1 ? pct(worst, 2) : 'n/a']
    }),
  )
  if (input.groups.length > 1) {
    y -= 2
    text(`These products are not all from one parent. ${input.groups.length} groups are `
      + 'internally consistent and inconsistent with each other, so only one of them was '
      + 'carried forward.', 7.6, 'Helvetica-Bold', 2.4, WARN)
  }

  // --- how m was chosen -------------------------------------------------------------------
  heading('How the agreement threshold was chosen')
  text('Requiring more products to agree lowers contamination and narrows the marker set toward '
    + 'probes nearly every product called. Those are enriched for sites where the parent is '
    + 'homozygous because the minor allele is rare, and unrelated people carry the common allele '
    + 'there too, so past a point the reference stops excluding them. The threshold is the '
    + 'deepest setting whose retained marker set still holds 90% of the heterozygosity measured '
    + 'at the shallowest.', 7.4, 'Helvetica', 2.4, GREY)
  y -= 2
  table(
    [{ head: 'Products agreeing', w: 110 }, { head: 'Ascertainment', w: 90, right: true },
      { head: '', w: 260 }],
    Object.entries(input.ascertainment).sort((a, b) => Number(a[0]) - Number(b[0]))
      .map(([m, r]) => [
        `at least ${m}`,
        pct(r, 0),
        Number(m) === p.mMin ? 'chosen' : r < 0.9 ? 'too narrow' : '',
      ]),
  )

  // --- the products, each judged without itself -------------------------------------------
  heading('Every product, scored against a reference built without it')
  text('A product scored against a reference containing itself reads exactly zero absence, a '
    + 'bias the size of the whole signal, so each is judged by the others alone.',
  7.4, 'Helvetica', 2.4, GREY)
  y -= 2
  table(
    [{ head: 'Product', w: 130 }, { head: 'Absent', w: 70, right: true },
      { head: 'Ceiling', w: 70, right: true }, { head: 'Ratio', w: 60, right: true },
      { head: 'Reads', w: 90 }],
    input.members.map((m) => [
      m.name, pct(m.absence, 2), pct(m.ceiling, 2), ratio(m.ratio),
      { v: reads(m.verdict), colour: m.verdict === 'parent_genome_present' ? INK : WARN },
    ]),
  )

  if (input.controls.length) {
    heading('Controls')
    table(
      [{ head: 'Sample', w: 130 }, { head: 'Expected', w: 150 },
        { head: 'Absent', w: 70, right: true }, { head: 'Ratio', w: 60, right: true },
        { head: 'Reads', w: 90 }],
      input.controls.map((c) => [
        c.name, c.role, pct(c.absence, 2), ratio(c.ratio),
        { v: reads(c.verdict), colour: c.verdict === 'no_parental_contribution' ? INK : WARN },
      ]),
    )
  }

  // --- every array that went in, including the ones that did not qualify -------------------
  heading('Every array submitted')
  table(
    [{ head: 'Sample', w: 120 }, { head: 'Role', w: 118 }, { head: 'Group', w: 56 },
      { head: 'Call', w: 48, right: true }, { head: 'BAF band', w: 58, right: true },
      { head: 'Note', w: 132 }],
    input.samples.map((s) => [
      s.sample, s.role, s.group || '',
      s.callRate === undefined ? '' : pct(s.callRate, 1),
      s.hetBand === undefined ? '' : pct(s.hetBand, 1),
      s.excludedBecause ? { v: s.excludedBecause, colour: WARN } : '',
    ]),
  )

  // --- the boundary of the evidence --------------------------------------------------------
  heading('Not reported against this reference')
  text('Each of the following was measured coming out wrong rather than assumed to be a risk. A '
    + 'refusal means the evidence does not reach, which is a different claim from a negative '
    + 'result.', 7.4, 'Helvetica', 2.4, GREY)
  y -= 2
  table(
    [{ head: 'Withheld', w: 210 }, { head: 'Because', w: 322 }],
    input.limits.map((l) => [{ v: l.what, colour: WARN }, l.why]),
  )

  heading('Reproducing this')
  text('The reconstruction is deterministic: no randomness, and no threshold fitted on the '
    + 'samples being scored. The product files listed above and their checksums are the '
    + 'reference, and re-running this tool at the same version on the same files reproduces it '
    + 'exactly. The genotype itself is deliberately not distributed as a file, because it would '
    + 'be several hundred thousand homozygous calls belonging to an identifiable person in a '
    + 'format that could be re-imported as though it had been measured.',
  7.4, 'Helvetica', 2.4, GREY)

  footer()
  return pdf.save()
}

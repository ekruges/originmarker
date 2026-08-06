/**
 * The Syngamy report: everything the run measured, on paper, for a supplementary file.
 *
 * Typeset from the same numbers the page displays, with the same page furniture as the 1.x panel
 * PDF, and assembled in this tab because the genotypes never leave it. What it prints is meant to
 * be checkable rather than merely readable: every call carries the figure it was made on and the
 * reference it was made against, every constant carries where it came from, and every file carries
 * its SHA-256.
 */
import { LETTER, Pdf, wrap, type FontName } from './pdf.ts'
import {
  GLOSS, pct,
  type ChromResult, type PairClass, type PairResult, type ParentageResult,
} from './parentage.ts'
import type { Gate, SampleProfile } from './ingest.ts'
import type { RunResult } from './runlength.ts'
import { int, utc } from './fmt.ts'
import { EXAMPLE_CITATION } from './examples.ts'

export const DISCLAIMER = 'Research use only. Parent of origin here is a statistical call from '
  + 'array genotypes and requires confirmation by an independent method in a qualified genetics '
  + 'laboratory. Not a clinical diagnostic.'

const INK = '#000000'
const GREY = '#666666'
const BLUE = '#337ab7'
const RULE = '#cccccc'
const WARN = '#b12222'

export interface ReportFile {
  name: string
  size: number
  sha256: string
  role: 'donor' | 'oocyte' | 'sample'
  markers: number
  profile?: SampleProfile
  gates?: Gate[]
  result?: ParentageResult
  maternal?: ParentageResult
  paired?: PairResult
  error?: string
}

export interface ReportInput {
  files: ReportFile[]
  donorHeterozygosity: number
  startedAt: string | null
  generatedAt: string
  tool: string
  reportId: string
  fromExamples: boolean
  /** favicon.svg's text. The browser fetches it; a CLI run passes it in. */
  markSvg?: string
  /** A per-locus deletion test, if one was run after the genome-wide call. */
  locus?: {
    chrom: string
    pos: number
    eventSizeBp?: number
    bySample: { name: string; results: RunResult[] }[]
  } | null
}

interface Col { head: string; w: number; right?: boolean; font?: FontName; size?: number }

/** SHA-256 of a file, hex. Read at export time only, so dropping a file stays instant. */
export async function sha256(file: File): Promise<string> {
  const d = await crypto.subtle.digest('SHA-256', await file.arrayBuffer())
  return [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, '0')).join('')
}

/** A short run identifier, from the moment and the files, so two runs never collide on paper. */
export function reportId(seed: string): string {
  let h = 0x811c9dc5
  for (let i = 0; i < seed.length; i += 1) {
    h ^= seed.charCodeAt(i)
    h = Math.imul(h, 0x01000193) >>> 0
  }
  return h.toString(16).padStart(8, '0')
}

export interface MarkGeometry {
  view: number
  ring: [number, number, number, number]
  ringColour: string
  fill: string
  transform: [number, number, number, number]
  polys: [number, number][][]
}

/**
 * The monogram's geometry, from the text of favicon.svg.
 *
 * Read with regexes rather than a DOM parser so this module works wherever it is asked to
 * typeset: the browser fetches the file, and the audit harness reads it off disk. Not a general
 * SVG reader, and it returns null rather than guessing at anything it does not recognise, since
 * a report without a logo is still a report.
 */
export function parseMark(svg: string): MarkGeometry | null {
  const view = /viewBox\s*=\s*"([-\d.\s]+)"/.exec(svg)
  const circle = /<circle\b([^>]*)>/.exec(svg)
  const path = /<path\b([^>]*)>/.exec(svg)
  if (!view || !circle || !path) return null
  const attr = (tag: string, name: string): string | null => {
    const m = new RegExp(`${name}\\s*=\\s*"([^"]*)"`).exec(tag)
    return m ? m[1] : null
  }
  const vw = Number(view[1].trim().split(/\s+/)[2])
  const ring = ['cx', 'cy', 'r', 'stroke-width'].map((k) => Number(attr(circle[1], k)))
  const t = /translate\(([-\d.]+),([-\d.]+)\)\s*scale\(([-\d.]+),([-\d.]+)\)/
    .exec(attr(path[1], 'transform') ?? '')
  const d = attr(path[1], 'd')
  if (!vw || !t || !d || ring.some((n) => !Number.isFinite(n))) return null
  return {
    view: vw,
    ring: ring as [number, number, number, number],
    ringColour: attr(circle[1], 'stroke') ?? BLUE,
    fill: attr(path[1], 'fill') ?? BLUE,
    transform: t.slice(1).map(Number) as [number, number, number, number],
    polys: svgPolys(d),
  }
}

/** The monogram as the site ships it. Null if it cannot be read. */
async function mark(svgSource?: string): Promise<MarkGeometry | null> {
  try {
    if (svgSource !== undefined) return parseMark(svgSource)
    const res = await fetch('favicon.svg')
    return res.ok ? parseMark(await res.text()) : null
  } catch { return null }
}

/** Sub-paths built from straight lines. Anything curved is dropped rather than approximated. */
function svgPolys(d: string): [number, number][][] {
  const toks = d.match(/[A-Za-z]|-?[\d.]+/g) ?? []
  const polys: [number, number][][] = []
  let cur: [number, number][] = []
  let x = 0
  let y = 0
  let op = ''
  let i = 0
  while (i < toks.length) {
    if (/[A-Za-z]/.test(toks[i])) {
      op = toks[i]
      i += 1
      if ((op === 'Z' || op === 'z' || op === 'M') && cur.length) { polys.push(cur); cur = [] }
      if (op === 'Z' || op === 'z') continue
    }
    if (op === 'M' || op === 'L') {
      x = Number(toks[i]); y = Number(toks[i + 1]); i += 2; op = 'L'
    } else if (op === 'H') { x = Number(toks[i]); i += 1 } else if (op === 'V') {
      y = Number(toks[i]); i += 1
    } else return []
    cur.push([x, y])
  }
  if (cur.length) polys.push(cur)
  return polys
}

export async function buildReportPdf(input: ReportInput): Promise<Blob> {
  const [W, H] = LETTER
  const L = 40
  const R = W - 40
  const TOP = H - 42
  const BOT = 78
  const pdf = new Pdf(LETTER)
  pdf.setTitle('Syngamy')
  const logo = await mark(input.markSvg)

  let y = TOP
  let page = 0

  const footer = (): void => {
    pdf.setFont('Helvetica-Oblique', 7)
    pdf.setFillColor(INK)
    let yy = 44
    for (const line of wrap(DISCLAIMER, 'Helvetica-Oblique', 7, R - L)) {
      pdf.drawString(L, yy, line)
      yy -= 8.5
    }
    pdf.setFont('Helvetica', 6.5)
    pdf.setFillColor(GREY)
    pdf.drawString(L, 26, `Syngamy | ${input.tool} | report ${input.reportId} | generated `
      + `${input.generatedAt} | page ${page}`)
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

  /** Label on the left at a fixed column, value wrapped in the remainder. */
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
      // A cell that has to wrap sets the row's height, so a long note never overprints the
      // row beneath it.
      const lines = row.map((cell, i) => wrap(
        typeof cell === 'string' ? cell : cell.v, cols[i].font ?? 'Helvetica',
        cols[i].size ?? 6.8, cols[i].w - 6,
      ))
      const h = Math.max(...lines.map((l) => l.length)) * 8.4 + 0.6
      if (y - h < BOT) { newPage(); head() }
      let x = L
      row.forEach((cell, i) => {
        pdf.setFont(cols[i].font ?? 'Helvetica', cols[i].size ?? 6.8)
        pdf.setFillColor(typeof cell === 'string' ? INK : cell.colour)
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

  // --- masthead ---------------------------------------------------------------------------
  // Page 1 only. A mark at the top of a report is a masthead; on every page it is noise.
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
    pdf.polys(logo.polys.map((poly) => poly.map(([gx, gy]) => at(tx + sx * gx, ty + sy * gy))))
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica-Bold', 14)
    pdf.drawString(L + MARK + 8, y - 18, 'Syngamy')
    y -= MARK + 4
  } else {
    text('Syngamy', 14, 'Helvetica-Bold', 4)
  }

  const results = input.files.filter((f) => f.result)
  const donor = input.files.find((f) => f.role === 'donor')
  const builds = [...new Set(input.files.map((f) => f.profile?.build.build ?? 'undetermined'))]
  const oocyte = input.files.find((f) => f.role === 'oocyte')
  const paired = results.some((f) => f.maternal)
  const classOf = (f: ReportFile): PairClass => f.paired?.originClass ?? f.result!.originClass
  const counts = new Map<PairClass, number>()
  for (const f of results) counts.set(classOf(f), (counts.get(classOf(f)) ?? 0) + 1)

  card([
    ['Report:', input.reportId],
    ['Generated:', input.generatedAt],
    ['Run started:', input.startedAt ? utc(input.startedAt) : 'not recorded'],
    ['Software:', `${input.tool}, Syngamy module`],
    ['Sperm donor:', donor ? `${donor.name}  (${int(donor.markers)} markers, SHA-256 `
      + `${donor.sha256.slice(0, 16)})` : 'none'],
    ['Oocyte donor:', oocyte ? `${oocyte.name}  (${int(oocyte.markers)} markers, SHA-256 `
      + `${oocyte.sha256.slice(0, 16)})`
      : 'none supplied, so a maternal origin is inferred from paternal absence rather than '
        + 'measured'],
    ['Donor heterozygosity:', Number.isFinite(input.donorHeterozygosity)
      ? `${pct(input.donorHeterozygosity, 3)} autosomal, of called markers` : 'not measured'],
    ['Samples:', `${results.length} analysed  (`
      + `${[...counts].map(([k, v]) => `${v} ${k.replace(/_/g, ' ')}`).join(', ')
        || 'none called'})`],
    ['Assembly:', builds.join(', ')],
    ['Data handling:', 'every file was read in the browser that produced this report. No '
      + 'genotype was transmitted, and none is retained.'],
    ...(input.fromExamples
      ? [['Example data:', `this run used the bundled public example files. ${EXAMPLE_CITATION}`
          + ' Every eighth marker, values unaltered.'] as [string, string]]
      : []),
  ])

  // --- summary ------------------------------------------------------------------------------
  heading('Result summary')
  const times = (a: number, b: number): string =>
    (Number.isFinite(a / b) ? `${(a / b).toFixed(1)}x` : '-')
  table([
    { head: 'Sample', w: 150 },
    { head: 'Origin', w: 68 },
    ...(paired
      ? [{ head: 'Pat absence', w: 50, right: true }, { head: 'Pat x ceil', w: 44, right: true },
        { head: 'Mat absence', w: 50, right: true }, { head: 'Mat x ceil', w: 44, right: true }]
      : [{ head: 'Absence', w: 46, right: true }, { head: 'Ceiling', w: 46, right: true },
        { head: 'x ceiling', w: 46, right: true }]),
    { head: 'Zygosity', w: 62 },
    { head: 'Sperm', w: 46 },
    { head: 'Informative', w: 46, right: true },
  ], results.map((f) => {
    const r = f.result!
    const m = f.maternal
    const cls = classOf(f)
    return [
      f.name,
      {
        v: cls.replace(/_/g, ' '),
        colour: cls === 'unclear' || cls === 'neither_parent' ? WARN : INK,
      },
      ...(paired
        ? [pct(r.genomeRate), times(r.genomeRate, r.explainable),
          m ? pct(m.genomeRate) : '-', m ? times(m.genomeRate, m.explainable) : '-']
        : [pct(r.genomeRate), pct(r.explainable), times(r.genomeRate, r.explainable)]),
      r.zygosity.replace(/_/g, ' '),
      r.spermType === 'unknown' ? '-' : r.spermType.replace('_', '-'),
      int(r.informative),
    ]
  }))
  text('Absence is the rate at which the donor\'s obligate allele is missing where he is '
    + 'homozygous. Ceiling is what this sample\'s own noise can manufacture: its no-call rate '
    + 'times its heterozygous fraction. A sample is called as lacking a paternal contribution '
    + 'only above 3x that ceiling, and as carrying one only at or below it; in between it is '
    + 'left uncalled. No confidence percentage is reported, deliberately: see Interpretation.',
  6.5, 'Helvetica-Oblique', 2, GREY)

  // --- one block per sample -------------------------------------------------------------------
  for (const f of results) {
    const r = f.result!
    const m = f.maternal
    const cls = classOf(f)
    need(96)
    y -= 8
    pdf.setStrokeColor(RULE)
    pdf.setLineWidth(0.6)
    pdf.line(L, y, R, y)
    y -= 6

    text(`${f.name}   vs ${donor?.name ?? 'sperm donor'}`
      + `${m && oocyte ? ` and ${oocyte.name}` : ''}`, 8, 'Helvetica', 3, GREY)

    // The verdict reversed out of a filled bar. This page gets printed in black and white, so
    // the call has to survive losing its colour.
    const bar = 22
    need(bar + 4)
    pdf.setFillColor(cls === 'unclear' || cls === 'neither_parent' ? WARN : BLUE)
    pdf.rect(L, y - bar, R - L, bar, true)
    pdf.setFillColor('#ffffff')
    pdf.setFont('Helvetica-Bold', 15)
    pdf.drawString(L + 7, y - bar + 6.5, cls.replace(/_/g, ' ').toUpperCase())
    pdf.setFont('Helvetica', 8)
    pdf.drawRightString(R - 7, y - bar + 8, cls === 'neither_parent'
      ? 'Neither declared parent accounts for this genome.' : GLOSS[cls])
    pdf.setFillColor(INK)
    y -= bar + 5

    if (r.segments.length) {
      // A banded callout rather than a line of text: a partly lost chromosome reads as "unclear"
      // in the whole-chromosome table, which looks like missing information rather than a located
      // event, and a reader skimming the report will not go looking for it.
      const chroms = [...new Set(r.segments.map((sg) => `chr${sg.chrom}`))].join(', ')
      const total = r.segments.reduce((a, sg) => a + sg.spanBp, 0)
      need(20)
      pdf.setFillColor(WARN)
      pdf.rect(L, y - 15, 3, 15, true)
      pdf.setFillColor(INK)
      pdf.setFont('Helvetica-Bold', 9.5)
      pdf.drawString(L + 8, y - 11, `Chromosomal change on ${chroms}`)
      y -= 17
      text(`The paternal genome is missing across ${r.segments.length === 1 ? 'one region'
        : `${r.segments.length} regions`} totalling ${(total / 1e6).toFixed(1)} Mb, at a rate the `
        + 'rest of this genome does not reach. The whole-chromosome verdict cannot show this: a '
        + 'chromosome that is partly lost reads as neither present nor absent. This is a LOSS of '
        + 'the paternal contribution over that region, not a statement about physical copy '
        + 'number, and not a gain: gains are not called on this platform at all. The regions are '
        + 'listed in full below.', 7.6, 'Helvetica', 2.6, INK, L + 8)
      y -= 2
    }

    text(`${m ? 'paternal ' : ''}absent ${pct(r.genomeRate)}  |  ceiling `
      + `${pct(r.explainable)}  |  ${times(r.genomeRate, r.explainable)}  |  `
      + `${int(r.informative)} informative markers  |  zygosity ${r.zygosity.replace(/_/g, ' ')}`
      + `${r.spermType === 'unknown' ? '' : `  |  ${r.spermType.replace('_', '-')} sperm`}`,
    8, 'Helvetica-Bold', 3)
    if (m) {
      text(`maternal absent ${pct(m.genomeRate)}  |  ceiling ${pct(m.explainable)}  |  `
        + `${times(m.genomeRate, m.explainable)}  |  ${int(m.informative)} informative markers`,
      8, 'Helvetica-Bold', 3)
    }
    if (f.paired?.notes.length) {
      for (const n of f.paired.notes) text(n, 7.4, 'Helvetica', 2.2)
      if (Number.isFinite(f.paired.agreement)) {
        text(`The two declared parents agree at ${pct(f.paired.agreement, 1)} of shared markers.`,
          7.4, 'Helvetica', 2.2, GREY)
      }
    }

    heading(m ? 'Evidence, sperm donor' : 'Evidence', 8.6)
    table([
      { head: 'Measurement', w: 152 },
      { head: 'Observed', w: 52, right: true },
      { head: 'Reference', w: 52, right: true },
      { head: 'The reference is', w: 246 },
    ], [
      ['Donor alleles absent', pct(r.genomeRate), pct(r.explainable),
        `this sample's no-call rate (${pct(r.noCallRate, 1)}) times its heterozygous fraction `
        + `(${pct(r.hetFraction, 1)}), plus a floor. Dropout manufactures `
        + `absence only by turning a heterozygous call homozygous and discarding the donor's `
        + `allele, so the bound is the product and neither term alone. Both factors are given `
        + `because a ratio can fall from either one moving.`],
      ['Alleles the donor lacks', pct(r.nonParentalRate), pct(r.secondParentExpected),
        'half the donor\'s heterozygosity. A second parent supplies an allele he lacks only '
        + 'where he is homozygous, at sum(pq), and his heterozygosity is sum(2pq).'],
      ['Heterozygous BAF band', pct(r.hetBand), '8.00%',
        'measured at 15-16% for diploid genomes and 1.3-3.4% for uniparental ones. Dropout does '
        + 'not mimic this: a biparental embryo at 33% dropout still reads 22.2%.'],
      ['Spread between chromosomes',
        Number.isFinite(r.dispersion) ? r.dispersion.toFixed(2) : '-', '0.35',
        'coefficient of variation of the per-chromosome rate. Measured at 1.11 for two arrays '
        + 'of one man and 0.76 for a degraded true offspring, against 0.10 for an unrelated '
        + 'adult and 0.10 for two genomes blended in one tube. Below the reference the '
        + 'difference is genome-wide rather than confined to part of the genome; it does not by '
        + 'itself separate an unrelated genome from a blended one.'],
      ['Cleanest chromosome',
        Number.isFinite(r.minChromRate) ? pct(r.minChromRate) : '-', pct(r.explainable),
        'the lowest per-chromosome absence rate. A genome that lost a segment leaves the rest '
        + 'untouched and reaches the ceiling somewhere; one that never reaches it anywhere is '
        + 'different everywhere.'],
    ])

    if (m) {
      heading('Evidence, oocyte donor', 8.6)
      table([
        { head: 'Measurement', w: 152 },
        { head: 'Observed', w: 52, right: true },
        { head: 'Reference', w: 52, right: true },
        { head: 'The reference is', w: 246 },
      ], [
        ['Donor alleles absent', pct(m.genomeRate), pct(m.explainable),
          'this sample\'s no-call rate times its heterozygous fraction, computed the same way '
          + 'as for the sperm donor. The two rates are comparable only because the instrument '
          + 'is the same.'],
        ['Alleles the donor lacks', pct(m.nonParentalRate), pct(m.secondParentExpected),
          'half the oocyte donor\'s heterozygosity, which is what a second parent contributes.'],
      ])
    }

    const sex = r.chroms.filter((c) => c.chrom === 'X' || c.chrom === 'Y' || c.chrom === '23')
    if (sex.length) {
      heading('Sex chromosomes', 8.6)
      chromTable(sex, r.explainable)
    }
    if (r.segments.length) {
      heading('Segments where the paternal genome is missing', 8.6)
      text('A partly lost chromosome reads as neither present nor absent, so the whole chromosome '
        + 'comes back unclear and the missing part goes unreported. Each region below is scored '
        + "against the median rate of this sample's OTHER chromosomes, which one event cannot "
        + 'move. The span is the resolution rather than the event: the smallest callable region '
        + 'is 2,400 informative markers, and a real event smaller than that is reported at the '
        + 'size of the window that found it.', 7.4, 'Helvetica', 2.4, GREY)
      table(
        [{ head: 'Chromosome', w: 150 }, { head: 'Span', w: 62, right: true },
          { head: 'Markers', w: 56, right: true }, { head: 'Absent', w: 54, right: true },
          { head: 'Against', w: 54, right: true }, { head: 'Score', w: 52, right: true }],
        r.segments.map((sg) => [
          `chr${sg.chrom} ${int(sg.startBp)}-${int(sg.endBp)}`,
          `${(sg.spanBp / 1e6).toFixed(1)} Mb`,
          int(sg.markers), pct(sg.rate, 2), pct(sg.nullRate, 2), sg.score.toFixed(0),
        ]),
      )
    }
    heading(m ? 'Chromosomes, sperm donor' : 'Chromosomes', 8.6)
    chromTable(r.chroms, r.explainable)
    if (m) {
      heading('Chromosomes, oocyte donor', 8.6)
      chromTable(m.chroms, m.explainable)
    }

    if (f.profile) {
      heading('Sample quality', 8.6)
      const p = f.profile
      const baf = [p.coding.meanBafHom0, p.coding.meanBafHet, p.coding.meanBafHom2]
        .map((x) => (x === null ? '-' : x.toFixed(3))).join('  /  ')
      card([
        ['Markers:', `${int(p.markers)} read, ${int(p.called)} called (${pct(p.callRate, 1)})`],
        ['No-call:', pct(p.nocallRate, 2)],
        ['Heterozygous:', `${pct(p.hetRate, 2)} of called`],
        ['chrX / autosomal het:', p.chrXHetRatio === null ? '-' : p.chrXHetRatio.toFixed(3)],
        ['Sex call:', p.sex],
        ['Product:', p.product],
        ['Assembly:', p.build.build
          ? `${p.build.build}  (${int(p.build.illegal.GRCh37)} placements illegal under GRCh37, `
            + `${int(p.build.illegal.GRCh38)} under GRCh38, ${int(p.build.tested)} tested)`
          : `undetermined  (${p.build.note})`],
        ['Allele coding:', `${p.coding.verdict}  ${p.coding.note}`],
        ['Mean BAF AA/AB/BB:', baf],
        ['SHA-256:', f.sha256],
        ['File size:', `${int(f.size)} bytes`],
      ], 108)
    }

    if (f.gates?.length) {
      heading('Quality gates', 8.6)
      table([
        { head: 'Gate', w: 108 },
        { head: 'Value', w: 46, right: true },
        { head: 'Verdict', w: 56 },
        { head: 'Basis', w: 292 },
      ], f.gates.map((g) => [
        g.name,
        g.value === null ? '-' : pct(g.value, 1),
        { v: g.verdict.replace(/_/g, ' '), colour: g.verdict === 'usable' ? INK : WARN },
        g.detail,
      ]))
    }

    const notes = [...new Set([...r.notes, ...m?.notes ?? []])]
    const limits = [...new Set([...r.limits, ...m?.limits ?? []])]
    if (notes.length) {
      heading('Findings', 8.6)
      for (const n of notes) text(n, 7.4, 'Helvetica', 2.2)
    }
    if (limits.length) {
      heading('Limits on this call', 8.6)
      for (const l of limits) text(l, 7.4, 'Helvetica-Bold', 2.2, WARN)
    }
  }

  function chromTable(chroms: ChromResult[], ceiling: number): void {
    table([
      { head: 'Chr', w: 34 },
      { head: 'Informative', w: 52, right: true },
      { head: 'Absent', w: 46, right: true },
      { head: 'Rate', w: 46, right: true },
      { head: 'x ceiling', w: 44, right: true },
      { head: 'Verdict', w: 76 },
      { head: 'Note', w: 234 },
    ], chroms.map((c) => [
      `chr${c.chrom}`, int(c.informative), int(c.absent), pct(c.rate),
      // A chromosome whose calls are not measuring it has no ratio to report; printing one
      // invites the reader to take the number and ignore the verdict beside it.
      c.verdict === 'not_measured' || !Number.isFinite(c.rate / ceiling)
        ? '-' : (c.rate / ceiling).toFixed(1),
      { v: c.verdict.replace(/_/g, ' '), colour: c.verdict === 'absent' ? WARN : INK },
      c.note ?? '',
    ]))
  }

  // --- per-locus test -------------------------------------------------------------------------
  if (input.locus?.bySample.length) {
    const L = input.locus
    newPage()
    heading(`Per-locus deletion test at chr${L.chrom}:${int(L.pos)}`, 11)
    text('Whether the paternal contribution is specifically absent around one site, rather than '
      + 'across the genome. Each sample was re-read on this chromosome and scored marker by '
      + 'marker: 1 where the sperm donor\'s allele is provably present, 0 where it is provably '
      + 'absent, and nothing where the marker cannot say. The statistic is the length of the '
      + 'longest run of consecutive absences, against the longest run independent genotyping '
      + 'error would produce at the same marker count. Contiguity is what separates a deletion '
      + 'from scattered error.', 8, 'Helvetica', 2.6)
    y -= 3
    text('This test requires the oocyte donor and was not run without her: absence is Mendelian '
      + 'and needs nothing from her, but presence is an identity claim and needs her homozygous '
      + 'for the other allele, or nothing scores as present and no run can be broken.',
    8, 'Helvetica', 2.6)
    if (L.eventSizeBp !== undefined) {
      y -= 3
      text(`An event size of ${int(L.eventSizeBp)} bp was declared, so a window whose resolution `
        + 'floor exceeds it reports below resolution rather than absence.', 8, 'Helvetica', 2.6)
    }

    for (const s of L.bySample) {
      need(70)
      y -= 6
      text(s.name, 8.6, 'Helvetica-Bold', 3)
      table([
        { head: 'Window', w: 84 },
        { head: 'L3', w: 40, right: true },
        { head: 'Scored', w: 40, right: true },
        { head: 'Absent', w: 40, right: true },
        { head: 'Run', w: 32, right: true },
        { head: 'r_min', w: 36, right: true },
        { head: 'p', w: 48, right: true },
        { head: 'Mat run', w: 40, right: true },
        { head: 'Floor bp', w: 52, right: true },
        { head: 'Verdict', w: 120 },
      ], s.results.map((w) => [
        w.window.replace(/_/g, ' '),
        int(w.nL3), int(w.nScored), int(w.zSum), int(w.longestRun),
        w.rMin === null ? '-' : String(w.rMin),
        w.nL3 ? w.runP.toExponential(1) : '-',
        w.nMat ? int(w.longestRunMaternal) : '-',
        w.resolutionFloorBp === null ? '-' : int(Math.round(w.resolutionFloorBp)),
        {
          v: w.verdict.replace(/_/g, ' '),
          colour: w.verdict === 'significant_run' ? WARN : INK,
        },
      ]))
      text('L3 is the markers where the sperm donor is homozygous, which is what enters r_min; '
        + 'scored is how many of those the sample was called at. Mat run is the same statistic '
        + 'against the oocyte donor over the same marker set: a genuine paternal deletion cannot '
        + 'produce a maternal run in the same place, and dropout produces both.',
      6.5, 'Helvetica-Oblique', 2, GREY)
      for (const w of s.results) {
        text(`${w.window.replace(/_/g, ' ')}: ${w.note}`, 7.4, 'Helvetica', 2.2)
      }
    }
    text('A run does not separate copy loss from copy-neutral loss of heterozygosity: both remove '
      + 'paternal alleles contiguously, and only the intensity channel distinguishes them. A '
      + 'window holding very few informative markers can reach significance on a short run, which '
      + 'is why the marker count is printed beside every verdict.', 7, 'Helvetica-Oblique', 2.2,
    GREY)
  }

  // --- methods --------------------------------------------------------------------------------
  newPage()
  heading('Methods', 11)
  const assembly = builds.length === 1 && builds[0] !== 'undetermined' ? builds[0] : null
  text(`Parent of origin was determined from SNP array genotypes using ${input.tool}. For each `
    + `of ${results.length} sample(s), the rate at which the sperm donor's obligate allele was `
    + 'absent was computed across autosomal markers where he is homozygous and the sample is '
    + 'called, and compared against the rate that sample\'s own genotyping noise can produce, '
    + 'taken as the product of its no-call rate and its heterozygous fraction. A sample was '
    + 'called as lacking a paternal contribution only where the observed rate exceeded that '
    + 'bound by 3-fold, and as carrying one only where it fell at or below the bound; rates '
    + 'between the two were left uncalled. Whether a second parent contributed was assessed from '
    + 'the rate at which the sample carries alleles the sperm donor does not possess, against '
    + 'half his heterozygosity, which is the contribution a second parent makes. Zygosity was '
    + 'taken from the fraction of B-allele frequencies falling in the heterozygous band, 0.35 to '
    + '0.65. Sperm type was inferred from the chrX rate rather than from chrY, so that an '
    + 'X-bearing sperm carrying a full paternal genome is distinguishable from an absent '
    + 'paternal contribution; neither case leaves a Y.', 8, 'Helvetica', 2.6)
  y -= 3
  text(paired
    ? 'An oocyte donor array was supplied, so each sample was tested against both parents by the '
      + 'same measurement and the class read off the pair. This separates a sample lacking a '
      + 'paternal contribution from one belonging to neither declared parent, which a single '
      + 'parent cannot do: with the sperm donor alone a maternal origin is the shape left behind '
      + 'by his absence rather than a presence that was measured. The two parents were also '
      + 'compared against each other, since two arrays of one person pass both tests and yield a '
      + 'confident biparental call.'
    : 'No oocyte donor array was supplied. A maternal origin below is therefore inferred from '
      + 'the absence of the paternal contribution rather than measured, and a sample belonging '
      + 'to neither declared parent is indistinguishable from one that is maternal in origin.',
  8, 'Helvetica', 2.6)
  y -= 3
  text(assembly
    ? `Assembly was determined from marker positions against the UCSC chromInfo and gap tables `
      + `and found to be ${assembly}. No liftOver was performed.`
    : 'Assembly could not be determined from marker positions against the UCSC chromInfo and gap '
      + 'tables, which carry GRCh37 and GRCh38 only. No liftOver was performed.',
  8, 'Helvetica', 2.6)
  if (input.locus?.bySample.length) {
    y -= 3
    text(`A per-locus deletion test was additionally run at chr${input.locus.chrom}:`
      + `${int(input.locus.pos)} over three windows of 25 kb, 10 Mb and the whole chromosome, `
      + 'using a run-length statistic over consecutive markers where the sperm donor is '
      + 'homozygous, with the significance threshold derived from the per-marker '
      + 'spurious-violation rate. Results are reported above.', 8, 'Helvetica', 2.6)
  }
  y -= 3
  const anySeg = results.some((f) => f.result!.segments.length)
  if (anySeg) {
    text('Segmental loss was scanned for within each chromosome, because a chromosome that is '
      + 'only partly lost reads as neither present nor absent and its whole-chromosome verdict '
      + 'is therefore uninformative. A multiscale likelihood-ratio scan over windows of at least '
      + '2,400 called informative markers was scored against a robust external null, the median '
      + "per-chromosome absence rate of that sample's OTHER chromosomes, so that one event cannot "
      + 'inflate the rate it is tested against. The calling threshold is empirical rather than a '
      + 'closed-form tail: absence artefact on amplified material is spatially clustered, so a '
      + 'tail computed under independence reports artefact as significant. It was set at 250 '
      + 'against a maximum of 139 reached by five genomes known to carry no event over 110 '
      + 'chromosome scans, and the weakest real event at the marker floor scores 431. A '
      + 'chromosome withheld by the allelic-ratio gate is not scanned. Reported spans are the '
      + 'resolution of the window that found the event, not the extent of the event itself.',
    8, 'Helvetica', 2.6)
    y -= 3
  }
  if (results.some((f) => f.result!.zygosity === 'diploid')) {
    text('Mosaicism was assessed from the continuous B-allele frequency rather than from the '
      + 'genotype derived from it. A genotype call is a threshold on that ratio and the '
      + 'heterozygous-to-homozygous boundary sits near 0.917, while a chromosome present in a '
      + 'fraction f of cells has a ratio of 1/(2-f), which does not cross that boundary until f '
      + 'is about 0.91; below that the call is unchanged and carries no information, so no number '
      + "of markers helps. Each chromosome's mean deviation from 0.5 over its heterozygous sites "
      + "was contrasted against the same sample's other chromosomes, in standard deviations, "
      + 'because a globally shifted array would make any absolute measure read as mosaic while a '
      + 'contrast cannot. Detection is reliable from about half the cells and partial from a '
      + 'third. No mosaic FRACTION is reported: inverting the statistic is biased low by roughly '
      + 'half, since the truncation that hides a mosaic from the genotype also removes the '
      + 'most-shifted markers from the mean. The axis is withheld on a uniparental genome, which '
      + 'is homozygous by construction and has no heterozygous sites for a mixture to shift.',
    8, 'Helvetica', 2.6)
    y -= 3
  }
  const fired = [...new Set(results.flatMap((f) => f.result!.limits))]
  if (fired.length) {
    text('The following limits applied to this run and are reported rather than resolved: '
      + fired.map((l, i) => `(${i + 1}) ${l.replace(/\.$/, '')}.`).join(' '), 8, 'Helvetica', 2.6)
    y -= 3
  }

  heading('Interpretation', 10)
  text('No confidence percentage is attached to a call, and this is deliberate. The per-sample '
    + 'likelihood ratio between "paternal genome present" and "absent" runs past 10^10000 on a '
    + 'genome-wide array, which is not a probability any reader should be handed. The honest '
    + 'figure is the empirical accuracy of the method, and on the pairs of known relationship it '
    + 'has been run against that is 9 of 9 correct, whose 95% lower bound is roughly 72%. About '
    + '300 consecutive correct calls would be needed to claim 99%. What is reported instead is '
    + 'the margin: the observed rate, the reference it was measured against, and the ratio '
    + 'between them, so a reader can see how far from the boundary the call sits.', 8,
  'Helvetica', 2.6)
  y -= 3
  text('The absence axis separates by roughly thirty-fold and carries the call. The '
    + 'alleles-the-donor-lacks axis separates by about 1.6-fold and only distinguishes a '
    + 'biparental genome from a paternal-only one that is heterozygous; treat a near-boundary '
    + 'call on that axis as provisional. A sample scoring in the uncalled band is not a weak '
    + 'positive, it is an absence of evidence: an array from the second parent measures dropout '
    + 'directly and settles it.', 8, 'Helvetica', 2.6)

  heading('Constants and where they come from', 10)
  table([
    { head: 'Constant', w: 132 },
    { head: 'Value', w: 118 },
    { head: 'Provenance', w: 282 },
  ], [
    ['Mendelian inconsistency floor', '0.0063',
      'Kothiyal 2019, 0.63% of variant sites across 1,314 nuclear families. A lower bound on the '
      + 'spurious-violation rate, never the value itself.'],
    ['Paternal-absence noise bound', 'no-call rate x heterozygous fraction',
      'Derived per sample, not fitted. Bounds every related pair measured, within 2x: 0.13% '
      + 'predicted against 0.05% observed, 0.18% against 0.16%, 10.4% against 9.69%.'],
    ['Absence calling margin', '3x',
      'How far above the noise bound an absence must sit before it is called.'],
    ['Residual absence floor', '0.005',
      'Absence on clean data from genotyping error alone, measured at 0.03% and 0.05%. Added to '
      + 'the bound so a near-perfect array still has a non-zero reference.'],
    ['Second-parent signal', 'half the donor\'s heterozygosity',
      'Derived. A second parent supplies an allele the donor lacks only where he is homozygous, '
      + 'at sum(pq); his heterozygosity is sum(2pq) over the same markers.'],
    ['Diploid heterozygous BAF band', '0.08 of markers in 0.35-0.65',
      'Measured at 15-16% for diploid genomes and 1.3-3.4% for uniparental ones. Dropout does '
      + 'not mimic this: a biparental embryo at 33% dropout still reads 22.2%.'],
    ['Segment floor', '2,400 called informative markers',
      'Titrated on real spliced genomes. At 1,200 markers the weakest construction scores 1.09x '
      + 'the null maximum, which is not a detection; at 2,400 it scores 3.1x. A marker count is '
      + 'not a resolution: array spacing runs from 1 bp to 21 kb.'],
    ['Segment threshold', '250 log-likelihood ratio',
      'Empirical, and fitted rather than validated. Five genomes carrying no event reached a '
      + 'maximum of 139 over 110 chromosome scans; the weakest real event at the floor scores '
      + '431. Awaits an out-of-sample clean cohort.'],
    ['Mosaic contrast', '8 standard deviations',
      'Per chromosome, on the heterozygous allelic ratio against the same sample\'s other '
      + 'chromosomes. Four bulk diploid arrays with no mosaic ran -1.65 to 5.18 over 88 '
      + 'chromosome observations; a titrated mosaic reads 9.7 to 22.7 at a third of cells and '
      + '11.7 to 48.5 at half. The genotype cannot see one below nine tenths at all.'],
    ['Allelic-ratio floor', '0.40 of BAFs outside 0.15-0.85',
      'Below this a chromosome mis-clustered and its genotypes are not measuring it. Correctly '
      + 'clustered chromosomes run 0.752 to 0.976; one reported ABSENT at 18.36% ran 0.130.'],
    ['Per-chromosome minimum', '200 informative markers',
      'Below this a chromosome is not reported: the rate is too noisy to place against the '
      + 'ceiling.'],
    ['Assembly tables', 'UCSC chromInfo and gap, hg19 and hg38',
      'Freely redistributable. No liftOver chain file is used or shipped: those carry a '
      + 'non-commercial field-of-use restriction Apache 2.0 cannot sublicense.'],
  ])

  // --- inventory ------------------------------------------------------------------------------
  heading('Files read', 10)
  table([
    { head: 'File', w: 152 },
    { head: 'Role', w: 34 },
    { head: 'Bytes', w: 54, right: true },
    { head: 'Markers', w: 46, right: true },
    { head: 'SHA-256', w: 246, font: 'Courier', size: 6 },
  ], input.files.map((f) => [
    f.name, f.role, int(f.size), int(f.markers), f.sha256,
  ]))
  text('The hash is of the file as read, so a reviewer can confirm the input without being sent '
    + 'it.', 6.5, 'Helvetica-Oblique', 2, GREY)

  heading('Citation', 10)
  text(`${input.tool}. Syngamy: parent of origin from SNP array genotypes. `
    + `Report ${input.reportId}, generated ${input.generatedAt}.`, 8, 'Helvetica', 2.6)
  if (input.fromExamples) {
    text(`Example data: ${EXAMPLE_CITATION}`, 8, 'Helvetica', 2.6)
  }
  text('Kothiyal P et al. Mendelian inconsistent signatures from 1,314 ancestrally diverse '
    + 'family trios. J Comput Biol 2019;26(5):405-419. doi:10.1089/cmb.2018.0253. Cited for the '
    + 'inconsistency floor above.', 8, 'Helvetica', 2.6)

  footer()
  pdf.showPage()
  return pdf.save()
}

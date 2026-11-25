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
import { BAND_WORD, KIND_WORD, runAlerts, type Defect } from './defects.ts'
import { unanswerable } from './abnormalities.ts'
import { directionOf, split, stretchWhere, stretchWord, type Mechanism, type Stretch } from './timing.ts'
import {
  pct,
  type ChromResult, type PairResult, type ParentageResult,
} from './parentage.ts'
import { segmentCoords } from './segments.ts'
import {
  aneuploidyCells, cardSays, chromCells, countText, disturbedText, gainCells, kindSeen, notReportable, NOT_REPORTED_SAYS,
  readAgainst, seenHeading, seenLocus, segmentKindText, stateWord, summaryParts, verdictText, KIND_SAYS, MEASURED, MEASURED_SOURCE, pdfCells,
  PARENT_SENTENCE, sexLine, WORD,
  type CardSays, type ChangeRow,
} from './pageSays.ts'
import { locus } from './stage.ts'
import { sexGates, shownSex } from './sexMaterial.ts'
import { DEFAULT_PERMUTATIONS, MIN_REGIONS } from './features.ts'
import { drawComparison } from './comparisonPdf.ts'
import type { ComparisonResult } from './comparison.ts'
import type { ParentSearchText } from './findParent.ts'
import type { Gate, SampleProfile } from './ingest.ts'
import type { RunResult } from './runlength.ts'
import { int, utc } from './fmt.ts'
import { EXAMPLE_CITATION, EXAMPLE_DENSITY } from './examples.ts'
import { barIn, barsApplied, barValue, type Setup } from './bars.ts'
import { withheldHeading, withheldRows } from './withheld.ts'
import { methodsParagraph, type Methods } from './methods.ts'
import { PROVEN_CHIP } from './platform.ts'

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
  /** This file is a reconstructed genotype rather than a measured array. */
  inferred?: boolean
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
  /**
   * The run-wide feature comparison, present only once someone chose to run it.
   *
   * ONE FOR THE RUN, not one per sample. A single chip contributes a handful of regions and the
   * matched null has no shape under five, so a per-sample comparison answers "no conclusion" on
   * almost every card. The question is a cohort question in any case.
   */
  comparison?: ComparisonResult
  /** The search among several files labelled sperm, when there was one. See findParent.ts. */
  fatherSearch?: ParentSearchText | null
  /** A per-locus deletion test, if one was run after the genome-wide call. */
  locus?: {
    chrom: string
    pos: number
    eventSizeBp?: number
    bySample: { name: string; results: RunResult[] }[]
  } | null
  /** Other arrays read beside these in the run, where the scoring had more than the report holds. */
  runArrays?: number
}

/** The setup the registry's bars are chosen by: parents loaded, other arrays in the run, comparison. */
export function reportSetup(input: Pick<ReportInput, 'files' | 'runArrays' | 'comparison'>): Setup {
  const loaded = (input.files.some((f) => f.role === 'donor') ? 1 : 0) + (input.files.some((f) => f.role === 'oocyte') ? 1 : 0)
  // A result read against a parent counts that parent, whatever files the report was handed.
  const read = Math.max(0, ...input.files.filter((f) => f.result).map((f) => (f.result!.noParent ? 0 : f.maternal ? 2 : 1)))
  const samples = input.files.filter((f) => f.result).length
  // A result grouped with other cells was read beside them, whatever the report holds.
  const grouped = input.files.some((f) => (f.result?.units ?? 0) > 1)
  return { parents: Math.max(loaded, read) as 0 | 1 | 2, run: grouped || (input.runArrays ?? samples - 1) > 0, comparison: !!input.comparison }
}

/** The methods paragraph of a report's run: the same function the page's copy block reads. */
export function reportMethods(input: Pick<ReportInput, 'files' | 'runArrays' | 'comparison' | 'tool'>): Methods {
  const donor = input.files.find((f) => f.role === 'donor')
  const results = input.files.filter((f) => f.result)
  return methodsParagraph({
    tool: input.tool, chip: PROVEN_CHIP,
    builds: [...new Set(results.map((f) => f.profile?.build.build ?? 'undetermined'))],
    samples: results.length,
    father: donor ? (donor.inferred ? 'rebuilt' : 'measured')
      : results.some((f) => !f.result!.noParent && f.result!.role === 'paternal') ? 'measured' : 'none',
    mother: input.files.some((f) => f.role === 'oocyte')
      || results.some((f) => !f.result!.noParent && (!!f.maternal || f.result!.role === 'maternal')),
    grouped: Math.max(0, ...results.map((f) => f.result!.units ?? 0)),
    sexRows: results.filter((f) => {
      const l = cardSays(f.result!, f.maternal, f.paired?.originClass ?? f.result!.originClass, !!donor?.inferred).list
      return [...l.counted, ...l.once, ...l.seen, ...l.notReported?.rows ?? []].some((g) => g.pieces[0]?.step === 'sex')
    }).length,
    applied: barsApplied(reportSetup(input)),
  })
}

/** The timing section of a sample not reportable: one line, no rows. */
export const NO_TIMING_SAYS = 'No timing is read on a sample that is not reportable.'

/** The tool string, or a plain statement that the server gave no version. */
export const toolLine = (tool: string): string => (/\d/.test(tool) ? tool : `${tool}, version not reported`)

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

/**
 * The Where cell of the timing table. A whole-chromosome row carries no coordinates: its end is the
 * chromosome's length only where a parent was loaded, so printed it would read as a measured span.
 */
export const timingWhere = (u: { chrom: string; startBp: number; endBp: number; whole?: boolean }): string =>
  (u.whole ? `chr${u.chrom}, whole chromosome`
    : `chr${u.chrom} ${(u.startBp / 1e6).toFixed(1)}-${(u.endBp / 1e6).toFixed(1)}Mb`)

/**
 * The timing table's lines: one per gain or loss, or one per stretch where its stretches read
 * differently. A finding has no direction and so no origin to time; it is listed in its own table,
 * and listing it here too printed each finding of a sample twice.
 */
export function timingRows(uniformity: readonly {
  chrom: string; kind: string; startBp: number; endBp: number; whole?: boolean; mechanism: Mechanism; why: string; stretches?: Stretch[]
}[]): { where: string; mechanism: Mechanism; why: string }[] {
  return uniformity.filter((u) => directionOf(u.kind)).flatMap((u) => (u.stretches && split(u.stretches)
    ? u.stretches.map((s) => ({ where: stretchWhere(u, s), mechanism: s.mechanism, why: s.why }))
    : [{ where: timingWhere(u), mechanism: u.mechanism, why: u.why }]))
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
    pdf.drawString(L, 26, `Syngamy | ${toolLine(input.tool)} | report ${input.reportId} | generated `
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
  // NO PARENTAL ARRAY, NO PARENTAL SENTENCE. Every text below that describes a donor's alleles, an
  // absence ratio or an inference from one is printed only where a parent was loaded.
  const noParent = !donor && !oocyte
  // The bars this run was decided by, from the registry: the table near the end and every value
  // quoted in the text read them, so a printed bar cannot differ from the one the scoring used.
  const applied = barsApplied(reportSetup(input))
  // The headline and the change list, from the module the card reads, so the two cannot differ.
  const fatherRebuilt = !!donor?.inferred
  const saysOf = new Map<ReportFile, CardSays>()
  for (const f of results) {
    saysOf.set(f, cardSays(f.result!, f.maternal, f.paired?.originClass ?? f.result!.originClass, fatherRebuilt))
  }
  const says = (f: ReportFile) => saysOf.get(f)!
  const counts = new Map<string, number>()
  for (const f of results) counts.set(says(f).head.word, (counts.get(says(f).head.word) ?? 0) + 1)

  card([
    ['Report:', input.reportId],
    ['Generated:', input.generatedAt],
    ['Run started:', input.startedAt ? utc(input.startedAt) : 'not recorded'],
    ['Software:', `${toolLine(input.tool)}, Syngamy module`],
    ['Sperm donor:', donor ? `${donor.name}  (${int(donor.markers)} markers, SHA-256 `
      + `${donor.sha256.slice(0, 16)})` : 'none'],
    // The inference sentence holds only where a father was loaded; with no parent nothing is inferred.
    ['Oocyte donor:', oocyte ? `${oocyte.name}  (${int(oocyte.markers)} markers, SHA-256 `
      + `${oocyte.sha256.slice(0, 16)})`
      : donor ? 'none supplied, so a maternal origin is inferred from paternal absence rather than '
        + 'measured' : 'none'],
    ['Donor heterozygosity:', Number.isFinite(input.donorHeterozygosity)
      ? `${pct(input.donorHeterozygosity, 3)} autosomal, of called markers` : 'not measured'],
    ['Samples:', `${results.length} analysed  (`
      + `${[...counts].map(([k, v]) => `${v} ${k.toLowerCase()}`).join(', ')
        || 'none called'})`],
    ['Assembly:', builds.join(', ')],
    ['Data handling:', 'every file was read in the browser that produced this report. No '
      + 'genotype was transmitted, and none is retained.'],
    ...(input.fromExamples
      ? [['Example data:', `this run used the bundled public example files. ${EXAMPLE_CITATION} `
          + EXAMPLE_DENSITY] as [string, string]]
      : []),
  ])

  if (input.fatherSearch) {
    heading('Which file is the father')
    card(input.fatherSearch.rows.map((r) => [r.name, `carries in ${r.carries}, none in ${r.none}, `
      + `not shown in ${r.unshown}`] as [string, string]))
    for (const s of [input.fatherSearch.line, input.fatherSearch.abstain, input.fatherSearch.limit]) {
      if (s) text(s, 8, 'Helvetica', 3)
    }
  }

  // --- summary ------------------------------------------------------------------------------
  heading('Result summary')

  // Said before anything else in the section. Every number below inherits it, and a reader who
  // takes only the first section must not come away thinking a person was measured.
  const inferredIn = input.files.filter((f) => f.inferred)
  if (inferredIn.length) {
    const msg = `Measured against a RECONSTRUCTED genotype, not an array of anyone: `
      + `${inferredIn.map((f) => `${f.name} (${f.role})`).join(', ')}. It was rebuilt from the `
      + 'haploid cells that parent produced, so it is homozygous everywhere by construction and '
      + 'covers only the markers those cells could establish. Absence rates against it are not '
      + 'the same quantity as absence against a measured array and must not be compared with one.'
    const h = wrap(msg, 'Helvetica-Bold', 9, R - L - 10).length * 11 + 3
    need(h + 4)
    pdf.setFillColor(WARN)
    pdf.rect(L, y - h, 3, h, true)
    text(msg, 9, 'Helvetica-Bold', 2, WARN, L + 10)
    y -= 4
  }
  const times = (a: number, b: number): string =>
    (Number.isFinite(a / b) ? `${(a / b).toFixed(1)}x` : '-')

  /** The counted changes, then how many more were seen once, in the fewest words; findings apart,
   *  by the card's own split (pageSays.rowDirection). */
  const changeOf = (f: ReportFile): string => summaryParts(says(f).list).changes
  const findingOf = (f: ReportFile): string => summaryParts(says(f).list).findings
  const found = results.filter((f) => findingOf(f))
  const flagged = results.filter((f) => changeOf(f))
  if (flagged.length) {
    // Said before the table as well as inside it. A reader who takes only the first section must
    // not miss a chromosome-scale finding because it sat in a column they skimmed past.
    const msg = `Chromosomal change in ${flagged.length} of ${results.length} sample(s): `
      + flagged.map((f) => `${f.name} ${changeOf(f)}`).join('; ')
      + '. Coordinates and evidence per sample below.'
    const h = wrap(msg, 'Helvetica-Bold', 9, R - L - 10).length * 11 + 3
    need(h + 4)
    pdf.setFillColor(WARN)
    pdf.rect(L, y - h, 3, h, true)
    text(msg, 9, 'Helvetica-Bold', 2, WARN, L + 10)
    y -= 4
  }
  if (found.length) {
    text(`Findings, neither a gain nor a loss, in ${found.length} of ${results.length} sample(s): `
      + found.map((f) => `${f.name} ${findingOf(f)}`).join('; ') + '.', 8.4, 'Helvetica-Bold', 3)
    y -= 2
  }

  const wide = paired ? 0 : 1
  table([
    { head: 'Sample', w: 100 + wide * 12 },
    { head: 'Origin', w: 58 },
    { head: 'Chromosomal change', w: 86 + wide * 12 },
    ...(paired
      ? [{ head: 'Pat abs', w: 42, right: true }, { head: 'Pat x', w: 36, right: true },
        { head: 'Mat abs', w: 42, right: true }, { head: 'Mat x', w: 36, right: true }]
      : [{ head: 'Absence', w: 44, right: true }, { head: 'Ceiling', w: 44, right: true },
        { head: 'x ceiling', w: 44, right: true }]),
    { head: 'Zygosity', w: 54 },
    { head: 'Sperm', w: 38 },
    { head: 'Informative', w: 40, right: true },
  ], results.map((f) => {
    const r = f.result!
    const m = f.maternal
    const head = says(f).head
    return [
      f.name,
      { v: head.word.toLowerCase(), colour: head.decided ? INK : WARN },
      says(f).list.notReported ? { v: 'not reportable', colour: GREY }
        : { v: [changeOf(f), findingOf(f) ? `findings: ${findingOf(f)}` : ''].filter(Boolean).join('; ') || 'none detected',
          colour: changeOf(f) ? WARN : GREY },
      ...(paired
        ? [pct(r.genomeRate), times(r.genomeRate, r.explainable),
          m ? pct(m.genomeRate) : '-', m ? times(m.genomeRate, m.explainable) : '-']
        : [pct(r.genomeRate), pct(r.explainable), times(r.genomeRate, r.explainable)]),
      r.zygosity.replace(/_/g, ' '),
      r.spermType === 'unknown' ? '-' : r.spermType.replace('_', '-'),
      int(r.informative),
    ]
  }))
  text('Chromosomal change is a whole chromosome or a region whose DNA is absent or present in '
    + 'extra copies, read from the array no longer genotyping it with the intensity agreeing, and '
    + 'reported per sample below with its coordinates. "none detected" is not a clean bill: a '
    + 'region under the size floor, or on a chromosome whose calls are not measuring it, is not '
    + 'reported either way.', 6.5, 'Helvetica-Oblique', 2, GREY)
  if (!noParent) text('Absence is the rate at which the donor\'s obligate allele is missing where he is '
    + 'homozygous. Ceiling is what this sample\'s own noise can manufacture: its no-call rate '
    + 'times its heterozygous fraction. A sample is called as lacking a paternal contribution '
    + `only above ${barValue(barIn(applied, 'ABSENCE_MARGIN'))} that ceiling,`
    + ' and as carrying one only at or below it; in between it is '
    + 'left uncalled. No confidence percentage is reported, deliberately: see Interpretation.',
  6.5, 'Helvetica-Oblique', 2, GREY)

  // --- how sure each word is: the measured count behind every word printed below ---------------
  heading('How sure each word is', 8.6)
  if (noParent) {
    text('With no parent loaded each change carries the word "position only": it is read from copy '
      + 'number, with no parent named. A row with neither a gain nor a loss is a finding.',
    7.4, 'Helvetica', 2.4, GREY)
  } else {
    text('The headline reads each loaded parent\'s genome chromosome by chromosome, and the cell\'s '
      + 'copy number; the genome-wide absence ratio above is pooled over every autosome, the changed '
      + 'ones included, so it decides the headline only where it calls a parent absent. Each change '
      + 'carries one word for its kind of call. A change found two ways or consistent is counted; a '
      + 'stretch seen by one reading is listed, not counted. Measured on the 31 published cells of '
      + 'Zuccaro 2020 Table S3 (GSE148488): right of scored, by kind of call and setup.',
    7.4, 'Helvetica', 2.4, GREY)
    {
      const kinds = ['two-ways', 'whole', 'second-test', 'run-only', 'allele-scan', 'unnamed', 'intensity',
        'gain-parent', 'imbalance', 'doubled'] as const
      const cellOf = (c: readonly [number, number] | undefined): string => (c ? `${c[0]} of ${c[1]}` : '-')
      table(
        [{ head: 'Word', w: 70 }, { head: 'Kind of call', w: 236 },
          { head: 'Both parents', w: 72, right: true }, { head: 'Father measured', w: 76, right: true },
          { head: 'Father rebuilt', w: 76, right: true }],
        kinds.map((k) => [WORD[k], KIND_SAYS[k], cellOf(MEASURED[k]?.two), cellOf(MEASURED[k]?.measured),
          cellOf(MEASURED[k]?.rebuilt)]),
      )
    }
    text(PARENT_SENTENCE, 7.6, 'Helvetica-Bold', 2.4)
    text(`Counts from ${MEASURED_SOURCE}. Right is the right place by the counting rule (calls and the `
      + 'stated region joined per chromosome overlap by half of the longer) and, where both name one, '
      + 'the right parent. Wrong means the paper states no change there; some may be the paper\'s misses.',
    6.5, 'Helvetica-Oblique', 2, GREY)
  }

  // --- one block per sample -------------------------------------------------------------------
  for (const f of results) {
    const r = f.result!
    const m = f.maternal
    const { setup, list, head, defects } = says(f)
    // Every state word below comes from pageSays.stateWord through its helpers.
    const seen = notReportable(r)
    need(96)
    y -= 8
    pdf.setStrokeColor(RULE)
    pdf.setLineWidth(0.6)
    pdf.line(L, y, R, y)
    y -= 6

    // The card's own line: with no parent loaded it names none.
    text(`${f.name}   ${readAgainst(r, !!m, donor?.name ?? 'sperm donor', oocyte?.name ?? 'oocyte donor')}`,
      8, 'Helvetica', 3, GREY)

    // The verdict reversed out of a filled bar. This page gets printed in black and white, so
    // the call has to survive losing its colour.
    const bar = 22
    need(bar + 4)
    pdf.setFillColor(head.decided ? BLUE : WARN)
    pdf.rect(L, y - bar, R - L, bar, true)
    pdf.setFillColor('#ffffff')
    pdf.setFont('Helvetica-Bold', 15)
    pdf.drawString(L + 7, y - bar + 6.5, `${head.word}${countText(list)}`)
    pdf.setFillColor(INK)
    y -= bar + 5
    text(head.line, 8.4, 'Helvetica-Bold', 3)

    // The changes as the card lists them: one row per change, each piece with its own refined
    // coordinates and the scan's own pair of numbers.
    const changeTable = (rows: ChangeRow[]) => table(
      [{ head: 'Where', w: 118 }, { head: 'Whose', w: 92 }, { head: 'Word', w: 62 },
        { head: 'Coordinates and signal', w: 260 }],
      rows.map((g) => {
        const [where, whose, word, detail] = pdfCells(g, setup)
        return [where, whose, { v: word, colour: word === 'one reading' || word === 'also seen' || word === 'not reported' ? GREY : INK },
          { v: detail, colour: GREY }]
      }),
    )
    if (list.counted.length) changeTable(list.counted)
    if (list.once.length) {
      text('Seen once, not confirmed', 7.6, 'Helvetica-Bold', 2, GREY)
      changeTable(list.once)
    }
    if (list.seen.length) {
      text('Also read, not counted as a change', 7.6, 'Helvetica-Bold', 2, GREY)
      changeTable(list.seen)
    }
    if (list.disturbed) text(disturbedText(list), 7.6, 'Helvetica', 2, GREY)
    if (list.notReported && (list.notReported.rows.length || list.notReported.disturbed)) {
      text('Seen on this array. The sample is not reportable, so none of this is reported as a gain or a loss.', 7.6, 'Helvetica-Bold', 2, GREY)
      if (list.notReported.rows.length) changeTable(list.notReported.rows)
      if (list.notReported.disturbed) text(disturbedText(list.notReported, false), 7.6, 'Helvetica', 2, GREY)
    }
    // Every refusal and every row not reported, with its cause, as the page lists them.
    {
      const held = withheldRows(r, list)
      text(withheldHeading(held.length), 7.6, 'Helvetica-Bold', 2, GREY)
      if (held.length) {
        table([{ head: 'What', w: 190 }, { head: 'Why', w: 220 }, { head: 'Measured', w: 122 }],
          held.map((w) => [w.what, { v: w.why, colour: GREY }, { v: w.measured || '-', colour: GREY }]))
      }
    }

    // THE FEATURE COMPARISON IS NO LONGER PRINTED HERE. It answers a different question from the
    // rest of this report: everything else asks whose a change is, that asks whether the change
    // sits where the genome breaks anyway. Printed side by side, the second reads as evidence about
    // the first, and it is not: the fragile compartment is established on BOTH parental genomes
    // from the first cell cycle. It is now an addon with its own report, folded in below only once
    // someone has chosen to run it.

    text(r.noParent ? `read with no parent loaded  |  zygosity ${r.zygosity.replace(/_/g, ' ')}`
      : `${m ? 'paternal ' : ''}absent ${pct(r.genomeRate)}  |  ceiling `
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

    // WHAT WOULD STOP A READER, BEFORE ANY EVIDENCE. Three checks the run performs reached no
    // surface at all before this: the parental array not matching the slot it was loaded under,
    // the sample being too damaged to interpret, and the declared material disagreeing with what
    // the array reads. Each one changes how every number below should be read, so each is printed
    // above them in the same banner the reconstructed-genotype notice uses.
    for (const a of runAlerts(r)) {
      const msg = `${a.headline} ${a.body}`
      const h = wrap(msg, 'Helvetica-Bold', 8.4, R - L - 10).length * 10.4 + 3
      need(h + 4)
      pdf.setFillColor(WARN)
      pdf.rect(L, y - h, 3, h, true)
      text(msg, 8.4, 'Helvetica-Bold', 3, WARN, L + 10)
    }

    if (!r.noParent) heading(m ? 'Evidence, sperm donor' : 'Evidence', 8.6)
    if (!r.noParent) table([
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
      ['Heterozygous BAF band', pct(r.hetBand), pct(barIn(applied, 'HET_BAND_DIPLOID').value),
        'measured at 15-16% for diploid genomes and 1.3-3.4% for uniparental ones. Dropout does '
        + 'not mimic this: a biparental embryo at 33% dropout still reads 22.2%.'],
      ['Spread between chromosomes',
        Number.isFinite(r.dispersion) ? r.dispersion.toFixed(2) : '-', barValue(barIn(applied, 'UNIFORM_CV')),
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
      heading(seenHeading(seen, 'Sex chromosomes'), 8.6)
      if (r.sexChromosome) text(sexLine(r.sexChromosome, r.sexMaterial), 7.6, 'Helvetica', 2)
      chromTable(sex, r.explainable, seen, defects)
    }
    const aneu = r.chroms.filter((c) => c.aneuploidy || c.lostParent || c.mixed
      || c.alleleGain?.verdict === 'gain' || c.alleleGain?.verdict === 'imbalance')
    if (aneu.length) {
      heading(seenHeading(seen, 'Aneuploidy'), 8.6)
      text('A chromosome that is not there cannot be genotyped, so a whole-chromosome change '
        + 'shows up as a collapse in call rate rather than as an allelic mismatch; the intensity '
        + 'then says which way it went. Measured over 1,012 chromosomes, an intact one calls at '
        + '0.78x to 1.16x of its genome median and never leaves -0.79 to +0.42 log2. A parent is '
        + 'attached only where something survives on the chromosome to attribute and the sample '
        + 'carries that parent elsewhere. With both parents loaded, a chromosome reading one '
        + 'parent\'s allele alone where the parents are opposite homozygotes is missing the other '
        + 'parent\'s copy; that row says PARENT ABSENT, because one remaining copy and two read the '
        + 'same at those markers and only the intensity can separate them. An extra copy is read '
        + 'from the paternal share against the cell\'s own other chromosomes, with intensity as a '
        + 'gate: a shift with intensity below the cell\'s level is an IMBALANCE and names no parent, '
        + 'because a partial loss of the other parent reads the same way.',
      7.4, 'Helvetica', 2.4, GREY)
      table(
        [{ head: 'Chromosome', w: 92 }, { head: 'Change', w: 56 },
          { head: 'Call rate', w: 62, right: true }, { head: 'Intensity', w: 62, right: true },
          { head: 'Whose copy', w: 156 }],
        aneu.map((c) => {
          const [change, whose] = aneuploidyCells(seen, c)
          return [`chr${c.chrom}`, { v: change, colour: seen ? GREY : WARN },
            `${c.callFraction.toFixed(2)}x`, `${c.lrrShift > 0 ? '+' : ''}${c.lrrShift.toFixed(2)} log2`, whose]
        }),
      )
    }
    if (r.gains.length) {
      if (r.findings?.length) {
        heading('Other chromosomal classes found', 8.6)
        text('These need no parental array at all: copy-neutral loss of heterozygosity, runs of '
          + 'homozygosity indicating uniparental isodisomy, ploidy, and a genome too disturbed to '
          + 'reference against. They are listed with the events above rather than apart from them, '
          + 'because a reader compares rows. Where a class carries no parental origin AT ALL, that '
          + 'is stated as a property of the class rather than left as an empty field.',
        7.4, 'Helvetica', 2.4, GREY)
        table(
          [{ head: 'Class', w: 132 }, { head: 'Where', w: 110 },
            { head: 'Origin', w: 92 }, { head: 'Evidence', w: 198 }],
          r.findings.map((f) => [
            { v: stateWord(seen, KIND_WORD[f.cls] ?? f.cls, kindSeen(f.cls)), colour: INK },
            f.chrom === 'genome' ? 'whole genome' : `chr${f.chrom}`,
            { v: f.originBlocked ? 'none for this class' : 'scored below',
              colour: f.originBlocked ? GREY : INK },
            { v: seen ? NOT_REPORTED_SAYS : [f.evidence, f.originBlocked, f.flag].filter(Boolean).join('. '), colour: GREY },
          ]),
        )
      }

      // WHAT WAS LOOKED FOR AND CANNOT BE ANSWERED. Without this table a report that lists what was
      // found says nothing about what was never checked, and on this platform the difference
      // between "no heterodisomy present" and "heterodisomy is unreachable from an embryo alone" is
      // the difference between a result and a misunderstanding.
      {
        const walls = unanswerable(!!r.twoParents, r.units ?? 1)
        if (walls.length) {
          heading('Classes this array cannot answer', 8.6)
          table(
            [{ head: 'Class', w: 150 }, { head: 'Why not', w: 382 }],
            walls.map((t) => [{ v: t.label, colour: INK },
              { v: t.limit ?? t.by, colour: GREY }]),
          )
        }
      }

      heading(seenHeading(seen, 'Extra copies, and where they came from'), 8.6)
      table(
        [{ head: 'Where', w: 116 }, { head: 'Kind', w: 72 }, { head: 'Origin', w: 84 },
          { head: 'Confidence', w: 108 },
          { head: 'Evidence', w: 152 }],
        r.gains.map((g) => [
          g.where,
          gainCells(seen, g).kind,
          { v: gainCells(seen, g).origin, colour: g.called ? INK : WARN },
          { v: g.confidence !== undefined && Number.isFinite(g.confidence)
            ? `${g.confidence.toFixed(4)}  ${BAND_WORD[g.band ?? ''] ?? ''}`.trim()
            : 'not stated',
          colour: g.band === 'A' || g.band === 'B' ? INK : GREY },
          { v: gainCells(seen, g).why, colour: GREY },
        ]),
      )
      text('An extra copy in a cell carrying ONE parent\'s genome is that parent\'s by '
        + 'construction, and the question worth asking is whether it is that parent\'s other '
        + 'homologue, which is meiotic and shows as heterozygosity in a genome that should have '
        + 'none, or a duplicate of the same homologue, which is bit-identical to a single copy '
        + 'and cannot be detected in any channel this array carries. Naming a PARENT for a gain '
        + 'needs a cell carrying both parents, and both parents loaded. NOT VALIDATED ON A TRUE '
        + 'POSITIVE: no confirmed gain with a known origin exists in the material behind these '
        + 'thresholds, so this has been shown to refuse correctly and never to fire correctly.',
      6.5, 'Helvetica-Oblique', 2, GREY)
    }
    // TIMING, WHICH ONLY EXISTS WHEN TWO UNITS OF ONE EMBRYO WERE ARRAYED, whether or not a gain was read.
    if (timingRows(r.uniformity ?? []).length) {
      heading('Gamete-borne or post-fertilisation', 8.6)
      text('Whether a change was carried by the gamete or arose after fertilisation cannot be '
        + 'read on one cell. Across the units of one embryo, read on each stretch between the '
        + 'breakpoints of all its cells: a change every unit carries, same direction and parent, '
        + 'is from the gamete with three units or more and only consistent with it with two; one '
        + 'unit gaining what another loses, or a change some units carry while another reads that '
        + 'stretch normal, arose after fertilisation. A change over stretches that read '
        + 'differently has a line for each. The units were grouped by genotype concordance rather '
        + 'than by their labels.',
      7.4, 'Helvetica', 2.4, GREY)
      if (seen) text(NO_TIMING_SAYS, 7.6, 'Helvetica', 2, GREY)
      else table(
        [{ head: 'Where', w: 132 }, { head: 'Timing', w: 118 }, { head: 'Evidence', w: 282 }],
        timingRows(r.uniformity!).map((u) => [
          u.where,
          { v: stretchWord(u), colour: u.mechanism === 'none' ? GREY : INK },
          { v: u.why, colour: GREY },
        ]),
      )
    }
    if (r.segments.length) {
      heading(seenHeading(seen, 'Segmental change'), 8.6)
      text('A partly lost chromosome reads as neither present nor absent, so the whole chromosome '
        + 'comes back unclear and the missing part goes unreported. Each region below is scored '
        + "against the median rate of this sample's OTHER chromosomes, which one event cannot "
        + 'move. The span is the resolution rather than the event: the smallest callable region '
        + 'is 2,400 informative markers, and a real event smaller than that is reported at the '
        + 'size of the window that found it.', 7.4, 'Helvetica', 2.4, GREY)
      table(
        [{ head: 'Chromosome', w: 118 }, { head: 'What', w: 82 },
          { head: 'Breakpoint interval', w: 78 },
          { head: 'Span', w: 48, right: true },
          { head: 'Markers', w: 44, right: true }, { head: 'Rate', w: 44, right: true },
          { head: 'Against', w: 44, right: true }, { head: 'Score', w: 40, right: true }],
        r.segments.map((sg) => [
          locus(sg.chrom, segmentCoords(sg).start, segmentCoords(sg).end),
          { v: segmentKindText(seen, sg, 'paternal'), colour: sg.kind === 'parental-absence' || seen ? INK : WARN },
          { v: seen && sg.arm ? seenLocus({ locus: '', chrom: sg.chrom, arm: sg.arm }).replace(/^, /, '') : segmentCoords(sg).interval,
            colour: segmentCoords(sg).localised ? GREY : WARN },
          `${(segmentCoords(sg).spanBp / 1e6).toFixed(1)} Mb`,
          int(sg.markers), pct(sg.rate, 1), pct(sg.nullRate, 1), sg.arm ? 'n/a' : sg.score.toFixed(0),
        ]),
      )
    }
    // PARENT OF ORIGIN, GENOTYPE CHANNEL.
    if (r.oneParent?.length) {
      const two = r.oneParent.some((o) => o.twoParents)
      heading(seenHeading(seen, 'Parent of origin, from genotypes'), 8.6)
      text('A marker where a loaded parent is homozygous and the sample carries the allele that '
        + "parent does not have is Mendelian evidence that the parent's copy is absent: "
        + 'amplification removes alleles and cannot invent one. The exclusive count is that '
        + 'evidence directly, and is checkable without trusting the model.',
      7.4, 'Helvetica', 2.4, GREY)
      text(two
        // BOTH ARRAYS LOADED. The asymmetry is the reason this reads the way it does, and a
        // reader deciding whether to act on a row needs the number that applies to it.
        ? 'Both parental arrays were loaded, so each was asked only the question it can answer. '
          + "This channel is right 0.8539 [0.6784, 1.0000] of the time about a loaded parent's OWN "
          + 'copy being absent, and 0.2890 [0.1273, 0.4506] about the other parent\'s, measured on '
          + 'material where the answer is known from dissection rather than from this tool. The '
          + 'parent named below is taken only from the array that reported its own copy missing. '
          + 'Where both arrays reported that, nothing is named: it asserts no parental copy at all, '
          + 'which on amplified material is more often the array than the genome.'
        : 'One parental array was loaded. This channel is right 0.8539 [0.6784, 1.0000] of the time '
          + "about that array's OWN copy being absent and 0.2890 [0.1273, 0.4506] about the other "
          + "parent's, so a row reading \u201cother parent lost\u201d is close to a coin flip. "
          + 'Loading the second parental array puts every event on the reliable direction.',
      7.4, 'Helvetica', 2.4, GREY)
      table(
        [{ head: 'Region', w: 138 }, { head: 'Verdict', w: 96 },
          { head: 'Whose copy', w: 82 },
          { head: 'Informative', w: 58, right: true },
          { head: 'Exclusive', w: 54, right: true },
          { head: 'Confidence', w: 92 }],
        r.oneParent.map((o) => [
          o.where,
          { v: seen ? verdictText(true, o.verdict) : o.verdict.replace(/-/g, ' '),
            colour: o.verdict === 'refused' ? GREY : INK },
          { v: o.twoParents
            ? (o.parent ? `${o.parent}${o.corroborated ? ', agreed' : ''}` : 'not named')
            : 'see verdict',
          colour: o.twoParents && o.parent ? INK : GREY },
          int(o.markers), int(o.exclusive),
          // This channel is capped at band B by construction, so it never prints "very confident".
          { v: Number.isFinite(o.posterior)
            ? `${o.posterior.toFixed(4)}  ${BAND_WORD[o.band ?? ''] ?? ''}`.trim()
            : 'not assigned',
          colour: o.band === 'B' ? INK : GREY },
        ]),
      )
    }

    // PARENT OF ORIGIN, DOSAGE CHANNEL. A separate table on purpose: it answers the same question
    // from a different measurement, on the events the genotype channel structurally cannot reach.
    if (r.dosageCalls?.length) {
      heading(seenHeading(seen, 'Whole chromosomes, from allele dosage'), 8.6)
      text('A whole chromosome is DETECTED by the collapse of its genotype call rate, so on exactly '
        + 'those events the genotype channel above has no evidence left. Allele dosage is read '
        + 'whether or not a genotype is emitted. Each chromosome is measured against the rest of '
        + "THIS array's own genome, which is what removes a directional bias that otherwise points "
        + 'at the parent that was not genotyped. The fraction column is the mosaic fraction the '
        + 'shift implies. Which parent that shift names depends on the copy-number class, because '
        + 'a gain inverts the direction a loss and a copy-neutral event share, so the class is '
        + 'marginalised into a calibrated probability rather than assumed, and the confidence '
        + 'column carries its band. Where a '
        + 'material and interval width have no detection floor at all, the row reads not evaluable '
        + 'rather than refused: no array of that kind could have answered, whatever its quality.',
      7.4, 'Helvetica', 2.4, GREY)
      table(
        [{ head: 'Chromosome', w: 64 }, { head: 'Verdict', w: 104 },
          { head: 'Confidence', w: 116 },
          { head: 'Material', w: 66 },
          { head: 'Shift', w: 44, right: true }, { head: 'z', w: 36, right: true },
          { head: 'Fraction', w: 40, right: true },
          { head: 'Window', w: 42, right: true }],
        r.dosageCalls.map((d) => [
          d.where,
          { v: seen ? verdictText(true, d.verdict) : d.verdict.replace(/-/g, ' '),
            colour: (d.verdict === 'loaded-parent' || d.verdict === 'other-parent')
              ? INK : GREY },
          // Bands C and D print grey, so a weak number cannot be read as a strong one.
          { v: d.confidence !== undefined && Number.isFinite(d.confidence)
            ? `${d.confidence.toFixed(3)}  ${BAND_WORD[d.band ?? ''] ?? ''}`.trim()
            : 'not stated',
          colour: d.band === 'A' || d.band === 'B' ? INK : GREY },
          d.material,
          Number.isFinite(d.shift) ? d.shift.toFixed(4) : '-',
          Number.isFinite(d.z) ? d.z.toFixed(2) : '-',
          Number.isFinite(d.impliedF) ? d.impliedF.toFixed(3) : '-',
          int(d.window),
        ]),
      )
    }

    heading(seenHeading(notReportable(r), m ? 'Chromosomes, sperm donor' : 'Chromosomes'), 8.6)
    chromTable(r.chroms, r.explainable, notReportable(r), defects)
    if (m) {
      heading(seenHeading(notReportable(r), 'Chromosomes, oocyte donor'), 8.6)
      chromTable(m.chroms, m.explainable, notReportable(r), defects)
    }

    if (f.profile) {
      if (r.stage) {
        heading('Material and amplification', 8.6)
        const d = Number.isFinite(r.stage.dropout) ? r.stage.dropout.toFixed(3) : 'not assigned'
        const sd = Number.isFinite(f.profile.hetBafSd)
          ? f.profile.hetBafSd.toFixed(4) : 'not measured'
        text(`Inferred stage: ${r.stage.stage}. ${r.stage.why}. Template copies per locus `
          + `${r.stage.templates}, reported for context and not used to predict dropout. Allele `
          + `dropout ${d} (basis: ${r.stage.basis}) parameterises every directional call, and `
          + `${r.stage.markerFloor} informative markers are required for one.`, 7.6)
        text(`The stage rests on TWO measurements. Heterozygosity, `
          + `${pct(f.profile.hetRate, 2)} of called markers, is how many heterozygous markers `
          + `SURVIVED, which tracks how much template this sample had. Allele-fraction spread at `
          + `those markers, ${sd}, is how far the survivors SCATTERED, which tracks how much `
          + 'amplification it went through. Unamplified genomic DNA scatters at about 0.088 and '
          + 'whole-genome amplified material at 0.21 to 0.30. Both are required for the bulk rung, '
          + 'because heterozygosity alone placed amplified samples there and handed them a drift '
          + 'constant seventeen times tighter than a single cell\'s.', 7.4)
        text(`Limits of this figure: ${r.stage.caveat}.`, 7.4)
        text('Stage is inferred from the array rather than declared. Dropout is a sampling failure '
          + 'during amplification: a heterozygous locus in a single cell is one molecule per '
          + 'allele before amplification, and an allele that fails to prime early is absent from '
          + 'everything downstream. The per-stage figures are calibration constants measured on '
          + 'this material rather than predictions from template count, which cannot reproduce '
          + 'them, and the low rungs are dominated by genotyping error rather than by loss. Call '
          + 'rate is gated before ploidy, because a failed amplification drives heterozygosity '
          + 'toward zero and would otherwise be read as one genome.', 7.4)
      }
      heading('Sample quality', 8.6)
      const p = f.profile
      const baf = [p.coding.meanBafHom0, p.coding.meanBafHet, p.coding.meanBafHom2]
        .map((x) => (x === null ? '-' : x.toFixed(3))).join('  /  ')
      card([
        ['Markers:', `${int(p.markers)} read, ${int(p.called)} called (${pct(p.callRate, 1)})`],
        ['No-call:', pct(p.nocallRate, 2)],
        ['Heterozygous:', `${pct(p.hetRate, 2)} of called`],
        ['Allele-fraction spread:', Number.isFinite(p.hetBafSd)
          ? `${p.hetBafSd.toFixed(4)} at heterozygous calls `
            + `(unamplified DNA ~0.088, amplified 0.21-0.30)`
          : 'not measured'],
        ['chrX / autosomal het:', p.chrXHetRatio === null ? '-' : p.chrXHetRatio.toFixed(3)],
        ['Sex call:', shownSex(p.sex, r.sexMaterial, r.sexChromosome)],
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
      ], sexGates(f.gates, r.sexMaterial).map((g) => [
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

  function chromTable(chroms: ChromResult[], ceiling: number, seen = false, rows: readonly Defect[] = []): void {
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
      { v: chromCells(seen, c, rows).verdict, colour: c.verdict === 'absent' && !seen ? WARN : INK },
      chromCells(seen, c, rows).note,
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
  {
    const methods = reportMethods(input)
    text(methods.text, 8, 'Helvetica', 2.6)
    for (const c of methods.cites) text(c.ref, 7.4, 'Helvetica', 2.2, GREY)
    y -= 3
  }
  heading('How each reading works', 9.5)
  const assembly = builds.length === 1 && builds[0] !== 'undetermined' ? builds[0] : null
  if (noParent) {
    text(`No parental array was supplied. Each of ${results.length} sample(s) was read from its own SNP `
      + `array with ${input.tool}: its stage, its sex chromosomes and its copy number, from its call `
      + 'rate, intensity and B-allele frequencies against its own chromosomes. No parental origin is '
      + 'named, and nothing is inferred about either parent.', 8, 'Helvetica', 2.6)
  } else text(`Parent of origin was determined from SNP array genotypes using ${input.tool}. For each `
    + `of ${results.length} sample(s), the rate at which the sperm donor's obligate allele was `
    + 'absent was computed across autosomal markers where he is homozygous and the sample is '
    + 'called, and compared against the rate that sample\'s own genotyping noise can produce, '
    + 'taken as the product of its no-call rate and its heterozygous fraction. A sample was '
    + 'called as lacking a paternal contribution only where the observed rate exceeded that '
    + `bound by ${barValue(barIn(applied, 'ABSENCE_MARGIN'))},`
    + ' and as carrying one only where it fell at or below the bound; rates '
    + 'between the two were left uncalled. Whether a second parent contributed was assessed from '
    + 'the rate at which the sample carries alleles the sperm donor does not possess, against '
    + 'half his heterozygosity, which is the contribution a second parent makes. Zygosity was '
    + 'taken from the fraction of B-allele frequencies falling in the heterozygous band, 0.35 to '
    + '0.65. Sperm type was inferred from the chrX rate rather than from chrY, so that an '
    + 'X-bearing sperm carrying a full paternal genome is distinguishable from an absent '
    + 'paternal contribution; neither case leaves a Y.', 8, 'Helvetica', 2.6)
  y -= 3
  if (!noParent) text(paired
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
      + `${barValue(barIn(applied, 'MIN_SEGMENT_MARKERS'))}, called and informative, was scored against a robust external null, the median `
      + "per-chromosome absence rate of that sample's OTHER chromosomes, so that one event cannot "
      + 'inflate the rate it is tested against. The calling threshold is empirical rather than a '
      + 'closed-form tail: absence artefact on amplified material is spatially clustered, so a '
      + `tail computed under independence reports artefact as significant. It was set at ${barValue(barIn(applied, 'SEGMENT_LRT'))} `
      + 'against a maximum of 139 reached by five genomes known to carry no event over 110 '
      + 'chromosome scans, and the weakest real event at the marker floor scores 431. A '
      + 'chromosome withheld by the allelic-ratio gate is not scanned. Reported spans are the '
      + 'resolution of the window that found the event, not the extent of the event itself.',
    8, 'Helvetica', 2.6)
    y -= 3
    text('The same scan was run a second time on a DIFFERENT indicator, because the one above '
      + "asks where a parent's alleles are missing, which a region can be while its DNA is "
      + 'entirely present. The copy-number channel asks where the DNA itself is gone, read from '
      + 'the array no longer calling at all, and a region is reported only where that no-call '
      + `rate rises AND the median intensity departs from the rest of the genome by at least ${barValue(barIn(applied, 'SEGMENT_LRR_SHIFT'))}`
      + '. That second requirement is load-bearing: across 46 arrays the scan without it '
      + 'returns 31 regions including a 44.9 Mb loss on a bulk diploid adult who cannot carry '
      + 'one, and with it 17 survive with none on any of the six bulk diploid arrays. The sign of '
      + 'the shift then gives the direction. Extra copies use the same machinery reversed; no '
      + 'segmental gain occurs in those 46 arrays, so that direction was validated on a '
      + 'constructed positive class, a block of a whole-chromosome-gain array spliced into a '
      + 'clean array of the same series, called copy-gain at all four sizes tested while the same '
      + 'splice from a euploid array returns nothing. A chromosome already called whole is not '
      + 'scanned for a region inside it, since that is one event described twice.',
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

  // Written to be pasted into a methods section, so it names versions, counts and definitions
  // rather than describing them. Only printed when a run actually scored placement.
  const placed = results.find((r) => (r.result?.placement ?? [])
    .some((e) => Number.isFinite(e.p)))
  if (placed) {
    const pl = (placed.result?.placement ?? []).filter((e) => Number.isFinite(e.p))
    heading('Methods: positional enrichment', 9.5)
    text('Called regions were tested for positional enrichment against reference feature sets on '
      + 'GRCh37/hg19. Common fragile sites were the classical aphidicolin-induced set, given '
      + 'coordinates from their published cytogenetic bands through the UCSC cytoBand table. '
      + 'Genes were UCSC refGene, from which intervals of genes exceeding 500 kb were merged. '
      + 'Centromeric and telomeric intervals were taken from the UCSC gap table. Replication '
      + 'timing was represented by Repli-seq valleys, the local minima of the ENCODE '
      + 'wgEncodeUwRepliSeq wavelet-smoothed signal, taken both from the BG02ES embryonic stem '
      + 'line and as a constitutive set called in at least 8 of 10 ENCODE cell lines; valleys '
      + 'mark the latest-replicating positions rather than late-replicating domains.', 7.6)
    text('Significance was assessed by permutation. Because a region can only be called where the '
      + 'array carries markers and where amplification produced calls, and because marker density '
      + 'covaries with gene density, intervals drawn uniformly along the genome are not a valid '
      + 'null. Each null interval was therefore drawn on the same chromosome as the observed '
      + 'region and constrained to contain the same number of informative markers, so that a '
      + 'region of n markers was compared against other stretches of n markers rather than '
      + `against an equivalent physical length. ${DEFAULT_PERMUTATIONS} permutations were drawn `
      + 'per feature set from a fixed seed, making each reported p reproducible from the same '
      + 'inputs. The reported p is the proportion of permutations at or beyond the observed value '
      + 'in the tail the observation falls in, so depletion is reported as well as enrichment, '
      + `with a Laplace correction giving a floor of 1/${DEFAULT_PERMUTATIONS + 1}. Enrichment was `
      + `not computed for fewer than ${MIN_REGIONS} regions.`, 7.6)
    text(`In this report ${pl.length} feature sets were tested against `
      + `${pl[0].regions} region(s). Effect size should govern interpretation rather than the p `
      + 'value alone: with a narrow null distribution a difference of a few per cent can reach a '
      + 'small p while carrying no biological meaning, which is why the null interval and its '
      + 'distribution are plotted beside every result.', 7.6)
    text('This analysis is positional and carries no information about parent of origin. It uses '
      + 'breakpoint position only and requires no parental genotype, and the late-replicating '
      + 'fragile compartment is reported to be established on both parental genomes from the '
      + 'first cell cycle, so a positional result here does not support a parental one.', 7.6)
  }

  heading('Interpretation', 10)
  if (noParent) {
    text('No parental call is made in this report, since no parental array was supplied. Copy '
      + 'number is read against each cell\'s own chromosomes, and each change is printed with its '
      + 'position and the reading it rests on.', 8, 'Helvetica', 2.6)
  } else text('No confidence percentage is attached to a call, and this is deliberate. The per-sample '
    + 'likelihood ratio between "paternal genome present" and "absent" runs past 10^10000 on a '
    + 'genome-wide array, which is not a probability any reader should be handed. The honest '
    + 'figure is the empirical accuracy of the method, and on the pairs of known relationship it '
    + 'has been run against that is 9 of 9 correct, whose 95% lower bound is roughly 72%. About '
    + '300 consecutive correct calls would be needed to claim 99%. What is reported instead is '
    + 'the margin: the observed rate, the reference it was measured against, and the ratio '
    + 'between them, so a reader can see how far from the boundary the call sits.', 8,
  'Helvetica', 2.6)
  y -= 3
  if (!noParent) text('The absence axis separates by roughly thirty-fold and carries the call. The '
    + 'alleles-the-donor-lacks axis separates by about 1.6-fold and only distinguishes a '
    + 'biparental genome from a paternal-only one that is heterozygous; treat a near-boundary '
    + 'call on that axis as provisional. A sample scoring in the uncalled band is not a weak '
    + 'positive, it is an absence of evidence: an array from the second parent measures dropout '
    + 'directly and settles it.', 8, 'Helvetica', 2.6)

  // THE RUN-WIDE FEATURE COMPARISON, once, after every sample and before the constants. Printed
  // only if someone ran it: an addon nobody chose has no result, and a heading with nothing under
  // it reads as a negative, which this analysis has a genuine one of.
  if (input.comparison) {
    newPage()
    y = drawComparison(pdf, input.comparison, input.comparison.regionNames, {
      y, newPage: () => { newPage(); return y },
    })
    y -= 6
  }

  heading('Bars applied in this run', 10)
  text('Every named bar the scoring compared a measurement against in this run, read from the code '
    + 'itself. Bars of a step this run did not take (a parent not loaded, no other arrays, no '
    + 'comparison) are not listed. Assembly was read against the UCSC chromInfo and gap tables, hg19 '
    + 'and hg38; no liftOver chain file is used or shipped.', 7.4, 'Helvetica', 2.4, GREY)
  table([
    { head: 'Bar', w: 150 },
    { head: 'Value', w: 62, right: true },
    { head: 'What it decides', w: 320 },
  ], applied.map((b) => [{ v: b.name, colour: INK }, barValue(b), { v: b.says, colour: GREY }]))

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

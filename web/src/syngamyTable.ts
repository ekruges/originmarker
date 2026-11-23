/**
 * Every change in a run as one tab-separated table.
 *
 * WHY THIS IS THE PRIMARY OUTPUT AND THE PDF IS THE RECORD. The reader of this tool is a
 * bioinformatician or research-track embryologist re-analysing arrays after the fact, and their
 * next step is joining these calls against their own sample metadata. A 12-page PDF is the wrong
 * artefact for that: it signals a clinical deliverable while the disclaimer says it is not one,
 * and it cannot be joined to anything. The table can.
 *
 * ONE ROW PER CHANGE, and the columns say where the parent came from rather than only what it is.
 * A parent MEASURED on the interval, one INHERITED from a genome-level call, and one that is a
 * bare direction are three different claims, and a table that flattened them into one "parent"
 * column would invite exactly the counting error the panel warns about.
 */
import type { ParentageResult } from './parentage.ts'
import { kindOf, NOT_REPORTED_SAYS, SEEN, WORD, type CardSays } from './pageSays.ts'
import { directionOf } from './timing.ts'

export type TableChange = {
  chrom: string
  startBp: number
  endBp: number
  kind: string
  origin?: string
  band?: string
  confidence?: number
  inheritedMargin?: number
  stage?: string
  informative?: number
  why?: string
  /** The word the result card puts on this change (pageSays.WORD). */
  word?: string
}

export type TableSample = {
  name: string
  originClass?: string
  /** The result card's headline word (pageSays.headlineOf). */
  headline?: string
  zygosity?: string
  role?: string
  genomeRate?: number
  explainable?: number
  informative?: number
  stage?: string
  changes: readonly TableChange[]
}

/**
 * One sample of the table, from the card's own reading. With no parent loaded `role` is empty: the
 * result's role is a placeholder there. On a sample not reportable every row's word is "not
 * reported", as on the card, and its class column states only how it read against the array.
 */
export function tableSampleOf(name: string, r: ParentageResult, says: CardSays): TableSample {
  return {
    name,
    headline: says.head.word,
    originClass: r.originClass,
    zygosity: r.zygosity,
    role: r.noParent ? undefined : r.role,
    genomeRate: r.genomeRate,
    explainable: r.explainable,
    informative: r.informative,
    stage: r.stage?.stage,
    // On a sample not reportable a row states what was seen against the rest of the array, with no
    // class, parent or evidence sentence that would state a copy lost or gained.
    changes: says.defects.map((d) => (says.list.notReported
      ? { chrom: d.chrom, startBp: d.startBp, endBp: d.endBp, kind: SEEN[directionOf(d.kind) ?? 'none'],
        stage: d.stage, informative: d.informative, why: NOT_REPORTED_SAYS, word: 'not reported' }
      : {
        chrom: d.chrom, startBp: d.startBp, endBp: d.endBp, kind: d.kind,
        origin: d.origin, band: d.band, confidence: d.confidence,
        inheritedMargin: d.inheritedMargin, stage: d.stage, informative: d.informative,
        why: d.why, word: WORD[kindOf(d, says.setup)],
      })),
  }
}

export const COLUMNS = [
  'sample', 'sample_headline', 'sample_origin_class', 'sample_zygosity', 'loaded_parent', 'material',
  'chrom', 'start_bp', 'end_bp', 'span_mb', 'class', 'word',
  'parent', 'parent_basis', 'band', 'confidence', 'inherited_margin',
  'independent_determination', 'informative_markers', 'evidence',
] as const

/**
 * How the parent on this row was arrived at, which decides whether the row is its own evidence.
 *
 * `inherited` rows in one sample are ONE determination however many rows carry it, so the
 * `independent_determination` column is 1 on the first and 0 on the rest. Anything summing that
 * column gets the number of real determinations rather than the number of rows.
 */
const basisOf = (c: TableChange): string => {
  if (c.band === 'inherited') return 'inherited-from-genome-call'
  if (c.band === 'F') return 'direction-only'
  if (!c.origin || c.origin === 'unclear') return 'none'
  return 'measured-on-interval'
}

const cell = (v: unknown): string => {
  if (v === undefined || v === null) return ''
  if (typeof v === 'number') return Number.isFinite(v) ? String(v) : ''
  // Tabs and newlines would break the format; the evidence text legitimately contains neither
  // after this, and a quote is not special in TSV.
  return String(v).replace(/[\t\r\n]+/g, ' ').trim()
}

export function syngamyTable(samples: readonly TableSample[]): string {
  const lines: string[] = [COLUMNS.join('\t')]
  for (const s of samples) {
    let inheritedSeen = false
    for (const c of s.changes) {
      const basis = basisOf(c)
      // The first inherited row in a sample carries the determination; the rest carry none.
      let independent = '1'
      if (basis === 'inherited-from-genome-call') {
        independent = inheritedSeen ? '0' : '1'
        inheritedSeen = true
      } else if (basis === 'none') independent = '0'
      lines.push([
        cell(s.name), cell(s.headline), cell(s.originClass), cell(s.zygosity), cell(s.role), cell(c.stage ?? s.stage),
        cell(c.chrom), cell(c.startBp), cell(c.endBp),
        cell(Number.isFinite(c.endBp - c.startBp) ? ((c.endBp - c.startBp) / 1e6).toFixed(3) : ''),
        cell(c.kind), cell(c.word),
        cell(c.origin === 'unclear' ? '' : c.origin), cell(basis), cell(c.band),
        // A CONFIDENCE ONLY WHERE ONE IS REPORTABLE. Bands C and D suppress their digits in the
        // interface for the reason the mosaicism episode established, and a table that carried
        // them anyway would be the channel through which they are quoted.
        cell(c.band === 'A' || c.band === 'B' ? c.confidence : undefined),
        cell(c.inheritedMargin), independent,
        cell(c.informative ?? s.informative), cell(c.why),
      ].join('\t'))
    }
  }
  return `${lines.join('\n')}\n`
}

/**
 * The run in one flat table: samples by stage and reportable class, changes by direction, extent
 * and parent, embryos, and when each change arose. Plan item 3.9, rule C.
 *
 * BUILT ONLY FROM PER-SAMPLE RECORDS. A change is a row the result card counts
 * (pageSays.changeList(...).counted), never a raw defect row: one event printed as a copy loss and
 * a parental absence is one change. A sample not reportable is counted among the samples and
 * contributes nothing else. audit/truth/overview_check.py recounts every number from the same
 * records without this file.
 */
import { reportOf, rowDirection, type ChangeList, type Report } from './pageSays.ts'
import { ORIGIN_WORD, stretchWord } from './timing.ts'
import { armBoundary } from './wholeArms.ts'

/** One counted row of a result card, as the overview reads it. */
export interface OverviewRow {
  chrom: string
  direction: 'gain' | 'loss' | 'none'
  whole: boolean
  startBp: number
  endBp: number
  origin: string
  /** The origin words printed on the row's changes, one entry per distinct word; empty for none. A
   *  change read by stretch adds each stretch's word, "not called" for a stretch with none. */
  when: string[]
  /** The row's pieces, so a harness can attach the timing of each as the page does. */
  pieces: { startBp: number; endBp: number; kind: string }[]
}

/** What the overview reads of one sample. */
export interface OverviewSample {
  stage: string
  report: Report
  rows: OverviewRow[]
  /** Rows the card lists as seen once, not confirmed. */
  once: number
  /** The genome is too disturbed to reference against (one finding). */
  disturbed: boolean
  /** The page's embryo group, where this sample is a unit of one with another. */
  embryo?: number
}

/** A sample's overview reading, from its card's change list. */
export function overviewSampleOf(
  list: ChangeList,
  r: { stage?: { stage?: string }; refusals?: readonly { claim: string; gate: string }[] },
  embryo?: number,
): OverviewSample {
  return {
    stage: r.stage?.stage ?? '',
    report: reportOf(r),
    rows: list.counted.map((g) => ({
      chrom: g.chrom,
      direction: rowDirection(g),
      whole: g.whole,
      startBp: g.startBp,
      endBp: g.endBp,
      origin: g.origin,
      // A change whose stretches read differently carries each stretch's word, "not called" included,
      // so it is counted as read differently and never under one word.
      when: [...new Set(g.pieces.flatMap((p) => (p.mechanismStretches ? p.mechanismStretches.map(stretchWord)
        : p.mechanism && ORIGIN_WORD[p.mechanism] ? [ORIGIN_WORD[p.mechanism]] : [])))],
      pieces: g.pieces.map((p) => ({ startBp: p.startBp, endBp: p.endBp, kind: p.kind })),
    })),
    once: list.once.length,
    disturbed: !!list.disturbed,
    embryo,
  }
}

/**
 * Length of each chromosome, from the gap table the tool ships (hg19_features.json `gaps`, as
 * features.normaliseTrack reads it). Where an arm ends at the centromere is wholeArms.armBoundary,
 * the boundary the whole-call rule reads.
 */
export type Arms = ReadonlyMap<string, { end: number }>

export function armsFrom(gaps: readonly { chrom: string; startBp: number; endBp: number }[]): Arms {
  const m = new Map<string, { end: number }>()
  for (const g of gaps) m.set(g.chrom, { end: Math.max(m.get(g.chrom)?.end ?? 0, g.endBp) })
  return m
}

export const ARM_MARGIN_BP = 5e6

/** Whole, arm or piece. An arm runs from one end of the chromosome to the p/q boundary, 5 Mb either way. */
export function extentOf(r: Pick<OverviewRow, 'chrom' | 'whole' | 'startBp' | 'endBp'>, arms: Arms): 'whole' | 'arm' | 'piece' {
  if (r.whole) return 'whole'
  const a = arms.get(r.chrom)
  const cen = armBoundary(r.chrom)
  if (!a || cen === undefined) return 'piece'
  const nearCen = (x: number) => Math.abs(x - cen) <= ARM_MARGIN_BP
  const p = r.startBp <= ARM_MARGIN_BP && nearCen(r.endBp)
  const q = nearCen(r.startBp) && r.endBp >= a.end - ARM_MARGIN_BP
  return p || q ? 'arm' : 'piece'
}

const STAGES = ['single-cell', 'blastomere', 'trophectoderm', 'haploid', 'bulk', 'failed'] as const
const STAGE_WORD: Record<string, string> = {
  'single-cell': 'single cell', blastomere: 'blastomere', trophectoderm: 'trophectoderm',
  haploid: 'haploid', bulk: 'bulk', failed: 'failed', '': 'no stage called',
}
const WHEN = [ORIGIN_WORD.meiotic, ORIGIN_WORD['consistent-with-gamete'], ORIGIN_WORD['mitotic-not-all']]

/** The table, one [label, count] per line, in print order. */
export function overviewTable(samples: readonly OverviewSample[], arms: Arms): [string, number][] {
  const out: [string, number][] = []
  const n = (label: string, v: number) => out.push([label, v])
  n('samples', samples.length)
  const stages = [...STAGES, ...new Set(samples.map((s) => s.stage).filter((s) => !(STAGES as readonly string[]).includes(s)))]
  for (const st of stages) n(`stage: ${STAGE_WORD[st] ?? st}`, samples.filter((s) => s.stage === st).length)
  for (const c of ['reportable', 'limited', 'not reportable'] as const) {
    n(`${c}${c === 'limited' ? ' (integrity disturbed)' : c === 'not reportable' ? ' (integrity uninterpretable)' : ''}`,
      samples.filter((s) => s.report === c).length)
  }
  const read = samples.filter((s) => s.report !== 'not reportable')
  const rows = read.flatMap((s) => s.rows)
  const changes = rows.filter((r) => r.direction !== 'none')
  n('changes with a gain or a loss', changes.length)
  for (const d of ['loss', 'gain'] as const) {
    const of = changes.filter((r) => r.direction === d)
    n(d === 'loss' ? 'losses' : 'gains', of.length)
    for (const e of ['whole', 'arm', 'piece'] as const) n(`${d}, ${e}`, of.filter((r) => extentOf(r, arms) === e).length)
    for (const p of ['maternal', 'paternal'] as const) n(`${d}, ${p}`, of.filter((r) => r.origin === p).length)
    n(`${d}, parent not named`, of.filter((r) => r.origin !== 'maternal' && r.origin !== 'paternal').length)
  }
  n('findings (copy-neutral, uneven shares, complex genome)',
    rows.filter((r) => r.direction === 'none').length + read.filter((s) => s.disturbed).length)
  n('stretches seen once, not counted', read.reduce((a, s) => a + s.once, 0))
  const groups = new Map<number, number>()
  for (const s of samples) if (s.embryo !== undefined) groups.set(s.embryo, (groups.get(s.embryo) ?? 0) + 1)
  n('embryos of two or more units', groups.size)
  n('samples in them', [...groups.values()].reduce((a, b) => a + b, 0))
  n('samples in no embryo group', samples.filter((s) => s.embryo === undefined).length)
  for (const k of [...new Set(groups.values())].sort((a, b) => a - b)) {
    n(`embryos of ${k} units`, [...groups.values()].filter((v) => v === k).length)
  }
  // A row whose pieces read two different origins is counted on its own line, never under either.
  for (const w of WHEN) n(`origin: ${w}`, changes.filter((r) => r.when.length === 1 && r.when[0] === w).length)
  n('origin: pieces read differently', changes.filter((r) => r.when.length > 1).length)
  n('origin: no origin word', changes.filter((r) => r.when.length === 0).length)
  return out
}

/** The same table as tab-separated text, for the download. */
export const overviewTsv = (t: readonly [string, number][]): string =>
  `line\tcount\n${t.map(([l, v]) => `${l}\t${v}`).join('\n')}\n`

/**
 * What a result card says: the headline word, the changes it counts, and the word on each change.
 *
 * One module for the card, the PDF, the summary table and the figure, so the four cannot say
 * different things about one sample. Pure, so `pageSays.check.ts` runs it in node.
 *
 * THE HEADLINE reads the per-chromosome parent check (parentGenome.ts) and the cell's copy number,
 * not the genome-wide absence ratio alone. That ratio is pooled over every autosome, so a cell's own
 * losses push it toward "unclear"; on the 31 published cells it left 21 of 31 unclear with both
 * parents loaded although the parent check found both genomes on 17 of them. The ratio still
 * decides where it says a parent is ABSENT: that is kept as a disagreement and the card says so.
 *
 * THE WORD ON A CHANGE follows the kind of call, because on the answer key the kind of call tracks
 * right and wrong and the band does not. Each word carries a measured count, from
 * audit/truth/page-says-2026-11-17.txt.
 *
 * Self-check: node --experimental-strip-types src/pageSays.check.ts
 */
import { int } from './fmt.ts'
import { pct, type ParentageResult } from './parentage.ts'
import type { DosageVerdict } from './dosageOrigin.ts'
import {
  defectsFrom, findingToDefect, KIND_WORD, mendelParent, parentNamed, withMechanism, type Defect,
} from './defects.ts'
import { segmentCoords } from './segments.ts'
import { directionOf, split, stretchWhere, stretchWord } from './timing.ts'

export type Parent = 'paternal' | 'maternal'

/* --- the headline ---------------------------------------------------------------------------- */

/** What the headline reads. Small on purpose, so a run record can supply it as well as a result. */
export interface HeadInput {
  noParent?: boolean
  /** The stage word, read with no parent. */
  stage?: string
  /** The array failed: its calls are noise. */
  failed: boolean
  zygosity: string
  /** The class the genome-wide ratio gives: the pair's with both parents, else the loaded one's. */
  cls: string
  /**
   * The sample's own damaged-sample flag, where it reads uninterpretable (integrity.ts).
   *
   * Carried as the words it printed, because the headline has to say which of the three
   * conditions fired. A sample this flag is set on must not have a parent read from it: each of
   * the three branches says so in its own limits, from "nothing on this sample is reportable" to
   * "do not read parental calls from this sample". The card used to print the flag in its notices
   * and a parental verdict above them at the same time.
   */
  damaged?: string
  /** True where the loaded parent is the mother alone. */
  motherAlone?: boolean
  loaded: {
    who: 'father' | 'mother'
    /** The per-chromosome parent check: carries, none, unshown or unread. */
    check: string
    present: number
    read: number
    /** The genome-wide ratio calls this parent's genome absent. */
    absent: boolean
    /** Absence over the noise ceiling, where measured. */
    ratio?: number
  }[]
}

export interface Headline {
  word: string
  /** Blue where the card states a result, orange where it does not. */
  decided: boolean
  /** One plain line under the word. */
  line: string
  /** Set where the word is "unclear": which rule kept it. */
  reason?: 'failed' | 'not-shown' | 'disagree' | 'copy-number' | 'class' | 'damaged'
}

const GLOSS_WORD: Record<string, string> = {
  androgenetic: 'Paternal genome only. No maternal complement.',
  gynogenetic: 'No paternal contribution. Maternal in origin.',
  biparental: 'Both parents contributed.',
  unclear: 'Not separated by the available evidence.',
}

const his = (who: 'father' | 'mother') => (who === 'father' ? 'his' : 'her')

/** The inputs from a scored result. `cls` is the pair's class where both parents were loaded. */
export function headInput(r: ParentageResult, m: ParentageResult | undefined, cls: string): HeadInput {
  const loadedWho = r.role === 'paternal' ? 'father' as const : 'mother' as const
  const ratio = (x: ParentageResult) => x.genomeRate / x.explainable
  const one = (x: ParentageResult, who: 'father' | 'mother', g: ParentageResult['parentGenome']) => ({
    who,
    check: g?.verdict ?? 'unread',
    present: g?.present ?? 0,
    read: g?.read ?? 0,
    absent: x.verdict === 'no_parental_contribution',
    ratio: Number.isFinite(ratio(x)) ? ratio(x) : undefined,
  })
  return {
    noParent: r.noParent,
    stage: r.stage?.stage,
    failed: r.stage?.stage === 'failed',
    zygosity: r.zygosity,
    cls,
    damaged: r.integrity?.level === 'uninterpretable' ? r.integrity.headline : undefined,
    motherAlone: !m && r.role === 'maternal',
    loaded: r.noParent ? [] : [
      one(r, loadedWho, r.parentGenome),
      ...(m ? [one(m, 'mother' as const, r.parentGenomeOther)] : []),
    ],
  }
}

/**
 * The headline rule of plan item 3.14 P1.
 *
 * Both parents (or the loaded parent and a second copy) where every loaded parent's check carries
 * that parent's genome, the cell reads two copies, and no loaded parent is called absent
 * genome-wide. Unclear, with its reason, on a failed array, a sample the run calls uninterpretable,
 * a parent not shown, a parent check and an absence ratio that disagree, or a copy number not read.
 * A genome read as one copy keeps the class the ratio gives, which is where that ratio is the right
 * instrument, unless that class is biparental, which a one-genome cell cannot be.
 */
export function headlineOf(h: HeadInput, changes: string): Headline {
  if (h.noParent) {
    return {
      word: (h.stage ?? 'unread').replace(/[_-]/g, ' ').toUpperCase(),
      decided: h.stage !== 'failed',
      line: `Read with no parent: the stage, the sex chromosomes and copy number. ${changes}`,
    }
  }
  const unclear = (reason: Headline['reason'], line: string): Headline => ({
    word: 'UNCLEAR', decided: false, reason, line: `${line} ${changes}`,
  })
  if (h.failed) {
    return unclear('failed', 'The array failed: its calls are noise, so no parent\'s genome can be measured in it.')
  }
  // A SAMPLE THE RUN ITSELF CALLS UNINTERPRETABLE GETS NO PARENTAL VERDICT. The flag and the
  // verdict used to be printed on one card: a cell at a 55.7% call rate read BIPARENTAL above its
  // own notice saying "do not read parental calls from this sample". The notice is the one that
  // is measured, so the headline says it instead of speaking over it. The changes keep their
  // positions, which is what that flag's limits leave standing.
  if (h.damaged) return unclear('damaged', `${h.damaged}, so no parent is read from this sample.`)
  if (h.zygosity === 'uniparental_homozygous' || h.cls === 'neither_parent') {
    // ONE GENOME AND BOTH PARENTS IN IT ARE TWO READINGS OF ONE CELL THAT CANNOT BOTH HOLD, and
    // the card said both: the pair's class printed BIPARENTAL beside the cell's own badge reading
    // uniparental homozygous. Neither reading is dropped; the headline says they disagree.
    if (h.zygosity === 'uniparental_homozygous' && h.cls === 'biparental') {
      return unclear('disagree', 'Two readings disagree. This cell reads one genome, and the '
        + 'genome-wide rates place both parents in it.')
    }
    const notIn = h.loaded.length === 1 && h.loaded[0].check === 'none'
    const word = h.cls.replace(/_/g, ' ').toUpperCase()
    return {
      word,
      decided: h.cls !== 'unclear' && !notIn,
      reason: h.cls === 'unclear' ? 'class' : undefined,
      line: (h.cls === 'neither_parent' ? 'Neither declared parent accounts for this genome.'
        : `${GLOSS_WORD[h.cls] ?? ''}${notIn ? ' That holds if the loaded array is this sample\'s '
          + 'parent, which the alleles cannot show.' : ''}`) + ` ${changes}`,
    }
  }
  const missing = h.loaded.find((p) => p.check !== 'carries')
  if (missing) {
    return unclear('not-shown', missing.check === 'none'
      ? `The ${missing.who}'s genome is not in this cell, so the loaded array may not be ${his(missing.who) === 'his' ? 'his' : 'hers'}.`
      : missing.check === 'unread'
        ? `Too few chromosomes could be read to show the ${missing.who}'s genome in this cell.`
        : `The ${missing.who}'s genome is not shown in this cell: found on ${missing.present} of `
          + `${missing.read} chromosomes read.`)
  }
  const absent = h.loaded.find((p) => p.absent)
  if (absent) {
    return unclear('disagree', `Two readings disagree. The ${absent.who}'s genome is found on `
      + `${absent.present} of ${absent.read} chromosomes read, but ${his(absent.who)} alleles are missing `
      + `across the genome${absent.ratio !== undefined ? ` at ${absent.ratio.toFixed(1)}x the noise ceiling` : ''}.`)
  }
  if (h.zygosity !== 'diploid') {
    return unclear('copy-number', 'The copy number could not be read, so one genome cannot be told from two.')
  }
  if (h.loaded.length === 2) {
    return { word: 'BOTH PARENTS', decided: true, line: `Both parents' genomes are in this cell. ${changes}` }
  }
  const who = h.loaded[0].who
  return {
    word: `${who.toUpperCase()} + SECOND PARENT`,
    decided: true,
    line: `The ${who}'s genome is in this cell, with a second parent's copy. ${changes}`,
  }
}

/* --- the kind of call, its word and its measured count ----------------------------------------- */

/** Which setup a count was measured in. The mother alone was not measured. */
export type Setup = 'two' | 'measured' | 'rebuilt' | 'mother' | 'none'

export const setupOf = (r: ParentageResult, m: ParentageResult | undefined, rebuilt: boolean): Setup =>
  (r.noParent ? 'none' : m ? 'two' : r.role === 'maternal' ? 'mother' : rebuilt ? 'rebuilt' : 'measured')

export type CallKind = 'whole' | 'two-ways' | 'run-only' | 'second-test' | 'allele-scan' | 'unnamed' | 'doubled'
  | 'intensity' | 'gain-parent' | 'imbalance' | 'finding' | 'inherited' | 'no-parent'

export type Word = 'found two ways' | 'consistent' | 'one reading' | 'also seen' | 'position only'
  | 'not reported'

const named = (d: Defect): d is Defect & { origin: Parent } => d.origin === 'paternal' || d.origin === 'maternal'

/** The kind of a printed call. Mirrors audit/truth/page_says.py, which measured each kind. */
export function kindOf(d: Defect, setup: Setup): CallKind {
  if (setup === 'none') return 'no-parent'
  if (d.band === 'inherited') return 'inherited'
  if (d.step === 'imbalance') return 'imbalance'
  if (d.step === 'finding') return 'finding'
  if (d.kind === 'copy-gain' && named(d)) return 'gain-parent'
  if (d.step === 'whole-alleles') return named(d) ? 'whole' : 'intensity'
  if (d.step === 'allele-scan') return 'allele-scan'
  if (d.step === 'run') {
    if (!named(d)) return 'unnamed'
    return d.second?.verdict === 'parent-lost' && d.second.parent === d.origin ? 'two-ways' : 'run-only'
  }
  if (d.step === 'window') return named(d) ? 'second-test' : 'unnamed'
  return 'intensity'
}

export const WORD: Record<CallKind, Word> = {
  whole: 'consistent',
  'two-ways': 'found two ways',
  'second-test': 'consistent',
  inherited: 'consistent',
  'run-only': 'one reading',
  'allele-scan': 'one reading',
  unnamed: 'one reading',
  intensity: 'one reading',
  'gain-parent': 'one reading',
  finding: 'one reading',
  imbalance: 'also seen',
  'no-parent': 'position only',
  doubled: 'one reading',
}

/** A change is counted in the headline when its word is one of these. */
export const COUNTED: ReadonlySet<Word> = new Set<Word>(['found two ways', 'consistent', 'position only'])

/** Where every count below was measured. */
export const MEASURED_SOURCE = 'audit/truth/page-says-2026-11-17.txt'

/**
 * Right of scored, per kind of call and setup, on the 31 published cells of Zuccaro 2020 Table S3
 * (GSE148488). Right is the right place by the plan's counting rule and, where the call and the
 * paper both name a parent, the right parent; for a row naming no parent, the right place.
 */
export const MEASURED: Partial<Record<CallKind, Partial<Record<'two' | 'measured' | 'rebuilt', readonly [number, number]>>>> = {
  whole: { two: [89, 96], measured: [93, 103], rebuilt: [92, 102] },
  'two-ways': { two: [19, 26] },
  'run-only': { two: [7, 20] },
  'second-test': { measured: [13, 15], rebuilt: [7, 11] },
  'allele-scan': { two: [3, 17], measured: [5, 21], rebuilt: [7, 24] },
  unnamed: { two: [22, 38], measured: [16, 36], rebuilt: [12, 21] },
  intensity: { two: [3, 7], measured: [4, 8], rebuilt: [4, 9] },
  'gain-parent': { two: [1, 4] },
  imbalance: { two: [10, 24], measured: [5, 10], rebuilt: [3, 4] },
  doubled: { two: [0, 3], measured: [0, 3], rebuilt: [0, 1] },
}

/** On a stretch in the place the paper states, the parent named: right, of named, three setups. */
export const PARENT_ON_STRETCH = [59, 59] as const

export const PARENT_SENTENCE = `Where a stretch sits where the published karyotype states a change, the parent this page `
  + `names was right on ${PARENT_ON_STRETCH[0]} of ${PARENT_ON_STRETCH[1]} (31 published cells, three setups). `
  + 'What is uncertain is whether a stretch is there, not whose it is.'

const SETUP_WORDS: Record<'two' | 'measured' | 'rebuilt', string> = {
  two: 'both parents loaded', measured: 'the father measured', rebuilt: 'the father rebuilt from his haploid cells',
}

export const KIND_SAYS: Record<CallKind, string> = {
  whole: 'The alleles show one parent\'s copy missing along the whole chromosome.',
  'two-ways': 'Two readings agree: the absence scan, and the Mendelian test over the same stretch.',
  'run-only': 'One reading: the absence scan. The Mendelian test over the same stretch did not confirm it.',
  'second-test': 'The Mendelian test against the loaded parent names whose copy is missing.',
  'allele-scan': 'One reading: the scan of one parent\'s alleles along the chromosome.',
  unnamed: 'One reading: a stretch with a parent\'s alleles missing, no parent named.',
  intensity: 'One reading: copy number from intensity.',
  'gain-parent': 'One reading: an extra copy, with its parent named.',
  finding: 'One reading.',
  inherited: 'The parent follows from the genome-level call: 12 of 12 right against dissection.',
  imbalance: 'The parents\' alleles read unbalanced along the chromosome. Not counted as a change.',
  'no-parent': 'Read from copy number, with no parent loaded.',
  doubled: 'One reading: a copy lost here, reported by the copy scan and by an allele scan over the same stretch.',
}

/** The measured count behind a word, in one sentence, or what stands in for it. */
export function measuredSentence(kind: CallKind, setup: Setup): string {
  if (kind === 'inherited' || kind === 'no-parent') return ''
  const s = setup === 'mother' ? 'measured' : setup
  const c = s === 'none' ? undefined : MEASURED[kind]?.[s]
  if (!c) return 'Not measured on the answer key in this setup.'
  const what = kind === 'unnamed' || kind === 'intensity' ? 'sit where the published karyotype states a change'
    : kind === 'imbalance' ? 'sit on a chromosome the published karyotype states changed'
      : 'agree with the published karyotype'
  return `${c[0]} of ${c[1]} calls of this kind ${what} (31 published cells, `
    + `${SETUP_WORDS[s as 'two' | 'measured' | 'rebuilt']}${setup === 'mother' ? '; not measured with the mother alone' : ''}).`
}

/* --- the change list --------------------------------------------------------------------------- */

/** One change as the card lists it: the pieces of one chromosome, kind and parent, joined. */
export interface ChangeRow {
  chrom: string
  gain: boolean
  origin: Parent | 'unclear'
  call: CallKind
  word: Word
  /** The printed calls it joins, in position order. Each keeps its own refined coordinates. */
  pieces: Defect[]
  whole: boolean
  startBp: number
  endBp: number
  /**
   * A copy-loss piece from the copy scan on a chromosome that also carries a whole-chromosome loss:
   * partLost.ts prints such a piece only where it sits a further copy below the rest, so it is a
   * second copy lost there, not the whole loss again.
   */
  deeper?: boolean
}

export interface ChangeList {
  /** Counted in the headline. */
  counted: ChangeRow[]
  /** Listed under them: seen by one reading, not confirmed. */
  once: ChangeRow[]
  /** Read, not counted as a change. */
  seen: ChangeRow[]
  /**
   * Set on a sample not reportable (reportOf): every row read on it, each with the word
   * "not reported", and the other three lists empty. Counted nowhere.
   */
  notReported?: { rows: ChangeRow[]; disturbed?: { on: number; of: number } }
  /**
   * Set where changes sit on so many autosomes that the genome has no undisturbed remainder
   * (abnormalities.detectComplex). A property of the whole genome, so it is a sentence, not a row.
   */
  disturbed?: { on: number; of: number }
}

const chromOrder = (c: string): number => (/^\d+$/.test(c) ? Number(c) : c === 'X' ? 23 : c === 'Y' ? 24 : 99)

/**
 * The rows, joined per chromosome, gain or loss, parent and kind of call: the unit the plan's
 * counting rule scores, so one loss printed as two adjacent pieces is one change.
 */
export function changeList(defects: readonly Defect[], setup: Setup, reportable = true): ChangeList {
  if (!reportable) {
    const l = changeList(defects, setup)
    const rows = [...l.counted, ...l.once, ...l.seen].map((g) => ({ ...g, word: 'not reported' as const }))
      .sort((a, b) => chromOrder(a.chrom) - chromOrder(b.chrom) || a.startBp - b.startBp)
    return { counted: [], once: [], seen: [], notReported: { rows, disturbed: l.disturbed } }
  }
  const groups = new Map<string, ChangeRow>()
  let disturbed: ChangeList['disturbed']
  for (const d of defects) {
    if (d.kind === 'complex' && d.chrom === 'genome') {
      const m = (d.why ?? '').match(/(\d+) of (\d+) autosomes/)
      disturbed = { on: m ? Number(m[1]) : NaN, of: m ? Number(m[2]) : NaN }
      continue
    }
    const call = kindOf(d, setup)
    const gain = d.kind === 'copy-gain' || d.kind === 'trisomy' || d.kind === 'segmental-duplication'
    const whole = d.locus === `chr${d.chrom}` || d.chrom === 'genome'
    const key = `${d.chrom}|${gain}|${d.origin}|${call}|${whole}`
    const g = groups.get(key)
    if (g) {
      g.pieces.push(d)
      g.startBp = Math.min(g.startBp, d.startBp)
      g.endBp = Math.max(g.endBp, d.endBp)
    } else {
      groups.set(key, {
        chrom: d.chrom, gain, origin: named(d) ? d.origin : 'unclear', call, word: WORD[call],
        pieces: [d], whole, startBp: d.startBp, endBp: d.endBp,
      })
    }
  }
  const rows = mergeDoubled([...groups.values()])
  for (const g of rows) g.pieces.sort((a, b) => a.startBp - b.startBp)
  for (const g of rows) {
    g.deeper = !g.whole && g.pieces.every((p) => p.kind === 'copy-loss' && p.step === 'intensity')
      && rows.some((w) => w.whole && w.chrom === g.chrom && !w.gain && w.call !== 'imbalance')
  }
  rows.sort((a, b) => chromOrder(a.chrom) - chromOrder(b.chrom) || a.startBp - b.startBp)
  return {
    counted: rows.filter((g) => COUNTED.has(g.word)),
    once: rows.filter((g) => g.word === 'one reading'),
    seen: rows.filter((g) => g.word === 'also seen'),
    disturbed,
  }
}

export type Report = 'reportable' | 'limited' | 'not reportable'

/** Whether a sample is reportable, from its own sample refusals: the test the overview counts by. */
export const reportOf = (r: { refusals?: readonly { claim: string; gate: string }[] }): Report => {
  const gates = (r.refusals ?? []).filter((x) => x.claim === 'sample').map((x) => x.gate)
  return gates.includes('integrity.uninterpretable') ? 'not reportable'
    : gates.includes('integrity.disturbed') ? 'limited' : 'reportable'
}

/**
 * A row's direction: the one definition of a change. A change is a gain or a loss; a counted row
 * with neither is a finding. The card, the PDF and the overview all read this.
 */
export const rowDirection = (g: Pick<ChangeRow, 'gain' | 'pieces'>): 'gain' | 'loss' | 'none' =>
  (g.gain ? 'gain' : g.pieces.some((p) => directionOf(p.kind) === 'loss') ? 'loss' : 'none')

/** What a card counts: its changes, and its findings with the disturbed genome as one of them. */
export function cardCounts(l: ChangeList): { changes: number; findings: number } {
  const changes = l.counted.filter((g) => rowDirection(g) !== 'none').length
  return { changes, findings: l.counted.length - changes + (l.disturbed ? 1 : 0) }
}

const plural = (n: number, one: string): string => `${n} ${one}${n === 1 ? '' : 's'}`

/** The count after the headline word: " · 7 changes · 1 finding", a zero part left out. */
export function countText(l: ChangeList): string {
  const { changes, findings } = cardCounts(l)
  return `${changes ? ` · ${plural(changes, 'change')}` : ''}${findings ? ` · ${plural(findings, 'finding')}` : ''}`
}

/** What a row seen on a sample not reportable states: how it read against the rest of the array. */
export const SEEN: Record<'gain' | 'loss' | 'none', string> = {
  gain: 'higher than the rest of this array',
  loss: 'lower than the rest of this array',
  none: 'read unlike the rest of this array',
}

export const NOT_REPORTED_SAYS = 'Not reported: this sample is not reportable.'

/** The sentence under an opened row: the kind of call and its measured count, or what stands in. */
export function rowSays(g: ChangeRow, setup: Setup): string {
  if (g.word === 'not reported') return NOT_REPORTED_SAYS
  if (rowDirection(g) === 'none' && g.call === 'no-parent') {
    return 'Read with no parent loaded. Neither a gain nor a loss, so it is a finding, not a change.'
  }
  return `${KIND_SAYS[g.call]} ${measuredSentence(g.call, setup)}`.trim()
}

/** Where a piece seen on a sample not reportable lies; an arm piece names its arm, with no call. */
export const seenLocus = (p: Pick<Defect, 'locus' | 'chrom' | 'arm'>): string =>
  (p.arm ? `${p.locus}, the ${p.chrom}${p.arm} arm: its edge and its extent inside the arm are not placed` : p.locus)

/** One piece's detail line; on a row not reported, its place and intensity only. */
export const pieceText = (p: Defect, g: ChangeRow): string =>
  (g.word === 'not reported' ? [seenLocus(p), lrrText(p)] : [locusText(p), pairText(p), lrrText(p)])
    .filter(Boolean).join(' · ')

/** A figure mark's words on a sample not reportable: what was seen, and that it is not reported. */
export function seenMarkText(m: { chrom: string; startBp: number; endBp: number; whole: boolean; gain: boolean; dir?: string }):
{ word: string; sentence: string } {
  const seen = SEEN[m.gain ? 'gain' : m.dir === 'none' ? 'none' : 'loss']
  const where = m.whole ? `Whole chromosome ${m.chrom}` : `Chromosome ${m.chrom}, ${mb1(m.startBp)} to ${mb1(m.endBp)} Mb`
  return { word: seen, sentence: `${where}: ${seen}. ${NOT_REPORTED_SAYS}` }
}

/** What the whole-cell figure draws of a result: one reading per printed row piece, plus each run withheld. */
export function drawnOf(r: { refusals?: readonly { gate: string; withheld?: unknown }[] }, defects: readonly Defect[]):
{ listed: number; withheld: number } {
  return {
    listed: defects.filter((d) => d.chrom !== 'genome').length,
    withheld: (r.refusals ?? []).filter((x) => x.gate === 'other-parent.under-span' && x.withheld).length,
  }
}

/**
 * The figure's note under its title, from the card's own counts: how many readings are drawn, what
 * the ones beyond the list's rows are, and what the card counts. None where every drawn reading is
 * one counted row.
 */
export function figureNote(l: ChangeList, drawn: { listed: number; withheld: number }): string | undefined {
  const rows = [...l.counted, ...l.once, ...l.seen, ...l.notReported?.rows ?? []]
  const marked = rows.filter((g) => g.chrom !== 'genome').length
  const parts = drawn.listed - marked
  const total = drawn.listed + drawn.withheld
  if (!total || (!l.notReported && !drawn.withheld && !parts && !l.once.length && !l.seen.length)) return undefined
  const extra = [
    parts > 0 ? `${parts} more ${parts === 1 ? 'piece' : 'pieces'} of those rows, since a row joins the pieces of one stretch` : '',
    drawn.withheld ? `${drawn.withheld} found and withheld, not in the list` : '',
  ].filter(Boolean)
  const head = `${total} ${total === 1 ? 'reading is' : 'readings are'} drawn and numbered below: `
    + `${plural(marked, 'row')} of the list${extra.length ? `, ${extra.join(', and ')}` : ''}.`
  if (l.notReported) return `${head} This sample is not reportable, so none of them is reported.`
  const { changes, findings } = cardCounts(l)
  // How many drawn readings the counted rows are, so the sentence adds up with what is drawn.
  const counted = l.counted.filter((g) => g.chrom !== 'genome').reduce((n, g) => n + g.pieces.length, 0)
  const rest = total - counted
  return `${head} The card counts ${plural(changes, 'change')} and ${plural(findings, 'finding')}, drawn as `
    + `${plural(counted, 'reading')}${!rest ? '' : rest === 1
      ? '; the other 1 is drawn as an outline and says on its own line why it is not counted'
      : `; the other ${rest} are drawn as outlines and say on their own line why they are not counted`}.`
}

/**
 * A sample in the fewest words, as the report's summary prints it: its changes and the stretches
 * seen once, and apart from them its findings, by the card's own split (rowDirection).
 */
export function summaryParts(l: ChangeList): { changes: string; findings: string } {
  const short = (g: ChangeRow): string => `chr${g.chrom} ${g.whole ? 'whole' : `${spanMb(g)}Mb`} ${whoseText(g)}`
  return {
    changes: [...l.counted.filter((g) => rowDirection(g) !== 'none').map(short),
      ...(l.once.length ? [`${l.once.length} more seen once`] : [])].join(', '),
    findings: [...l.counted.filter((g) => rowDirection(g) === 'none').map(short),
      ...(l.disturbed ? ['genome widely disturbed'] : [])].join(', '),
  }
}

/* --- what a sample not reportable may print of a chromosome or a row -------------------------- */

/** How a reading sat against the rest of the array; `same` is a reading at the array's own level. */
export type Seen = 'gain' | 'loss' | 'none' | 'same'

const SEEN_WORD: Record<Seen, string> = { ...SEEN, same: 'at the level of the rest of this array' }

/** The sample is one the overview counts as not reportable. */
export const notReportable = (r: { refusals?: readonly { claim: string; gate: string }[] }): boolean =>
  reportOf(r) === 'not reportable'

/**
 * THE ONE GUARD for every state word a chromosome, a stretch or a bar prints: on a sample not
 * reportable it returns only how the reading sat against the rest of the array. Every table, note,
 * badge and tooltip that names a loss or a gain takes its word from here; the numbers beside it
 * (rates, log2) are what was seen and stay.
 */
export const stateWord = (seen: boolean, word: string, dir: Seen): string => (seen ? SEEN_WORD[dir] : word)

/**
 * The chromosomes a list of interval labels covers, collapsed.
 *
 * A refusal that applies to three hundred intervals is one fact, and a reader wants to know how
 * far it reached rather than to scroll three hundred coordinates. The full list stays in the
 * report and the export.
 */
export const listChroms = (wheres: readonly string[]): string => {
  const seen: string[] = []
  for (const w of wheres) {
    const c = /^chr([\dXY]+)/.exec(w)?.[1]
    if (c && !seen.includes(c)) seen.push(c)
  }
  if (!seen.length) return `${wheres.length} intervals`
  const n = seen.length
  return n <= 6 ? `chr${seen.join(', chr')}` : `${n} chromosomes, chr${seen[0]} to chr${seen[n - 1]}`
}

/**
 * The run log's lines for a sample's gains, losses and absent stretches. On a sample not reportable
 * each states only how the reading sat against the rest of the array (stateWord), with no class
 * or parent, as every other surface of that sample does.
 */
export function changeLogLines(
  r: Pick<ParentageResult, 'gains' | 'losses' | 'segments' | 'refusals'> & { dosageCalls?: ParentageResult['dosageCalls'] },
  sampleName: string, soloRole: string,
): ['DONE' | 'WARN', string][] {
  const seen = notReportable(r)
  const out: ['DONE' | 'WARN', string][] = []
  // A NAMED PARENT IS ALWAYS PRINTED IN FULL; A REFUSAL ONCE PER DISTINCT REASON. A run on
  // single-cell material produced 319 of these, every one the same sentence; the intervals stay
  // listed, and the full table is in the report and the export. Grouped by class as well as
  // reason, so the summary says what it stood for. On a sample not reportable no class, parent or
  // reason is named: the origin was not read.
  const refusals = new Map<string, { cls: string; why: string; at: string[] }>()
  for (const c of r.dosageCalls ?? []) {
    const named = c.verdict === 'loaded-parent' || c.verdict === 'other-parent'
    const cls = (c as { cls?: string }).cls ?? 'dosage'
    if (named && !seen) { out.push(['DONE', `${cls} origin ${c.where}: ${c.why}`]); continue }
    const key = seen ? '' : `${cls}|${c.why}`
    const at = refusals.get(key)
    if (at) at.at.push(c.where)
    else refusals.set(key, { cls, why: c.why, at: [c.where] })
  }
  for (const { cls, why, at } of refusals.values()) {
    if (seen) {
      out.push(['WARN', `${at.length === 1 ? at[0] : `${at.length} intervals, on ${listChroms(at)}`}: origin ${verdictText(true, '')}`])
    } else if (at.length === 1) out.push(['WARN', `${cls} origin ${at[0]}: ${why}`])
    else out.push(['WARN', `${at.length} ${cls} intervals not evaluable for the same reason, on ${listChroms(at)}: ${why}`])
  }
  for (const g of r.gains ?? []) {
    out.push(seen ? ['WARN', `${g.where}: ${stateWord(true, '', 'gain')}. ${NOT_REPORTED_SAYS}`]
      : [g.called ? 'DONE' : 'WARN', `gain ${g.where}: ${g.origin}. ${g.why}`])
  }
  for (const l of r.losses ?? []) {
    out.push(seen ? ['WARN', `${l.where}: ${stateWord(true, '', 'loss')}. ${NOT_REPORTED_SAYS}`]
      : [l.called ? 'DONE' : 'WARN', `loss ${l.where}: ${l.origin} copy missing. ${l.why}`])
  }
  const segs = r.segments ?? []
  if (segs.length) {
    out.push(['WARN', seen
      ? `${sampleName}: ${segs.length} stretch(es) ${stateWord(true, '', 'none')}: `
        + `${segs.map((x) => `chr${x.chrom} ${(segmentCoords(x).spanBp / 1e6).toFixed(1)}Mb`).join(', ')}. ${NOT_REPORTED_SAYS}`
      : `${sampleName}: ${segs.length} segment(s) where a parental `
        + `genome is missing: ${segs.map((x) => `chr${x.chrom} ${x.absentParent ?? soloRole} `
          + `${(segmentCoords(x).spanBp / 1e6).toFixed(1)}Mb `
          + `${segmentCoords(x).localised ? segmentCoords(x).interval : '(not localised)'} `
          + `at ${pct(x.rate, 1)}`).join(', ')}`])
  }
  return out
}

/** A section heading over readings of a sample not reportable says they are seen and not reported. */
export const seenHeading = (seen: boolean, title: string): string =>
  (seen ? `${title}: readings seen on this array, not reported` : title)

/** Which way a class reads, for the guard. */
export const kindSeen = (kind: string): Seen => directionOf(kind) ?? 'none'

interface ChromState {
  chrom: string
  verdict: string
  note?: string
  aneuploidy?: 'loss' | 'gain'
  lostParent?: string
  aneuploidyParent?: string
  mixed?: boolean
  alleleGain?: { verdict: string; parent?: string | null }
  callFraction: number
  lrrShift: number
}

/** Which way a chromosome's state reads, for the guard. */
export function chromSeen(c: ChromState): Seen {
  if (c.aneuploidy) return c.aneuploidy
  if (c.alleleGain?.verdict === 'gain') return 'gain'
  if (c.lostParent) return 'loss'
  if (c.mixed || c.alleleGain?.verdict === 'imbalance') return 'none'
  return 'same'
}

const signed = (x: number): string => `${x > 0 ? '+' : ''}${x.toFixed(2)}`

/**
 * The whole call's note ("LOST: ...", "GAINED: ...") is written when the chromosome is first read and
 * outlives a whole call withdrawn later (arm rule, stage or normal-majority refusal, part-lost
 * chromosome). Such a chromosome has `aneuploidy` unset, which is also what keeps it out of the
 * report's aneuploidy table, and its final call is the card's change rows on it.
 */
const WHOLE_NOTE = /^(LOST|GAINED): /

function withdrawnNote(c: ChromState, rows: readonly { chrom: string; kind: string; startBp: number; endBp: number }[]): string {
  const on = rows.filter((d) => d.chrom === c.chrom)
  return on.length
    ? `not called whole: read in pieces, ${on.map((d) => `${directionOf(d.kind) ?? 'finding'} ${mb1(d.startBp)} to ${mb1(d.endBp)} Mb`).join(', ')}`
    : 'not called whole: no change stands on it'
}

/** A chromosome's verdict and note, as every chromosome table prints them. `rows` is the sample's
 *  change list (defectsForResult), read where a whole call was withdrawn. */
export function chromCells(
  seen: boolean, c: ChromState, rows: readonly { chrom: string; kind: string; startBp: number; endBp: number }[] = [],
): { verdict: string; note: string } {
  const verdict = stateWord(seen && c.verdict === 'absent', c.verdict.replace(/_/g, ' '), 'none')
  if (!seen && c.note && !c.aneuploidy && WHOLE_NOTE.test(c.note)) return { verdict, note: withdrawnNote(c, rows) }
  if (!seen || !c.note) return { verdict, note: c.note ?? '' }
  const nums = [Number.isFinite(c.callFraction) ? `calls at ${c.callFraction.toFixed(2)}x of the genome median` : '',
    Number.isFinite(c.lrrShift) ? `intensity ${signed(c.lrrShift)} log2` : ''].filter(Boolean).join(', ')
  return { verdict, note: [stateWord(true, '', chromSeen(c)), nums].filter(Boolean).join(': ') }
}

/** The change and whose-copy cells of the report's aneuploidy table. */
export function aneuploidyCells(seen: boolean, c: ChromState): [string, string] {
  const change = c.aneuploidy?.toUpperCase() ?? (c.mixed ? 'MIXED'
    : c.lostParent ? 'PARENT ABSENT'
      : c.alleleGain?.verdict === 'gain' ? 'GAIN (ALLELES)' : 'IMBALANCE')
  const whose = c.mixed ? 'not called: between the boundaries, one parent dominant'
    : !c.lostParent && !c.aneuploidy && c.alleleGain?.verdict === 'gain'
      ? `${c.alleleGain.parent} supplied the extra copy, from the alleles`
      : !c.lostParent && !c.aneuploidy && c.alleleGain?.verdict === 'imbalance'
        ? 'not named: intensity below the cell\'s level'
        : c.lostParent ? (c.aneuploidy === 'gain'
          ? `no ${c.lostParent} copy: every copy is the other parent's, from the alleles`
          : `${c.lostParent} copy absent, from the alleles`)
          : c.aneuploidyParent === 'this' ? 'the parent scored against'
            : c.aneuploidyParent === 'other' ? 'the other parent'
              : 'not determined'
  return [stateWord(seen, change, chromSeen(c)), seen ? 'not named: the sample is not reportable' : whose]
}

/** What a segment row is, in one line, on the page and in the report. The parent comes from the
 *  run that found the segment, falling back to the loaded role. */
export function segmentKindText(seen: boolean,
  sg: { kind: string; contested?: boolean; absentParent?: string }, role: string): string {
  const word = sg.kind === 'copy-loss' ? 'DNA absent'
    : sg.kind === 'copy-gain' ? 'extra copies'
      : sg.contested ? 'alleles absent, both parents'
        : `${sg.absentParent ?? role} alleles absent`
  return stateWord(seen, word, sg.kind === 'copy-gain' ? 'gain' : sg.kind === 'copy-loss' ? 'loss' : 'none')
}

/** An extra-copy annotation as printed: what, whose and why. */
export function gainCells(seen: boolean, g: { kind: string; origin: string; why: string }):
{ kind: string; origin: string; why: string } {
  return seen ? { kind: stateWord(true, g.kind, 'gain'), origin: 'not named', why: NOT_REPORTED_SAYS } : g
}

/** A verdict word of an origin table (genotypes, allele dosage), or what stands in for it. */
export const verdictText = (seen: boolean, verdict: string): string =>
  (seen ? 'not read: the sample is not reportable' : verdict.replace(/[-_]/g, ' '))

/** A copy level under a figure bar: its word, or on a sample not reportable how it sat. */
export const levelWord = (seen: boolean, word: string, copies: number, none: number, two: number): string =>
  stateWord(seen, word, copies < none ? 'loss' : copies > two ? 'gain' : 'same')

/** The disturbed-genome finding in one plain sentence, or nothing. */
export function disturbedText(l: Pick<ChangeList, 'disturbed'>, reported = true): string {
  if (!l.disturbed) return ''
  const { on, of } = l.disturbed
  const n = Number.isFinite(on) ? `${on} of ${of}` : 'Many'
  // On a sample not reportable nothing read is a change, so the sentence names none.
  return reported ? `${n} chromosomes carry a change: this cell's genome is widely disturbed.`
    : `${n} chromosomes read away from the rest: this cell's genome is widely disturbed.`
}

const ABSENCE_STEPS = new Set(['run', 'window', 'allele-scan'])

/**
 * One loss reported by the copy scan and by an absence or allele scan over the same stretch is one
 * row: each covers at least half of the longer of the two, the counting rule's overlap. Both calls
 * stay in `pieces`. The row names the parent where one call does; two calls naming different parents,
 * or an absence call already confirmed, stay two rows.
 */
function mergeDoubled(rows: ChangeRow[]): ChangeRow[] {
  const piece = (g: ChangeRow) => !g.whole && !g.gain && g.chrom !== 'genome'
  const copy = (g: ChangeRow) => piece(g) && g.pieces.every((p) => p.kind === 'copy-loss' && p.step === 'intensity')
  const absence = (g: ChangeRow) => piece(g) && g.pieces.every((p) => p.kind === 'parental-absence' && ABSENCE_STEPS.has(p.step ?? ''))
  const gone = new Set<ChangeRow>()
  for (const c of rows.filter(copy)) {
    const a = rows.find((x) => !gone.has(x) && absence(x) && x.chrom === c.chrom
      && Math.max(0, Math.min(x.endBp, c.endBp) - Math.max(x.startBp, c.startBp))
        >= 0.5 * Math.max(x.endBp - x.startBp, c.endBp - c.startBp)
      && x.word === 'one reading'
      && (c.origin === 'unclear' || x.origin === 'unclear' || c.origin === x.origin))
    if (!a) continue
    gone.add(c)
    a.pieces.push(...c.pieces)
    a.origin = a.origin !== 'unclear' ? a.origin : c.origin
    a.startBp = Math.min(a.startBp, c.startBp)
    a.endBp = Math.max(a.endBp, c.endBp)
    // Never found two ways: the two scans read one signal, a lost copy takes its alleles with it.
    a.call = 'doubled'
    a.word = WORD.doubled
  }
  return rows.filter((g) => !gone.has(g))
}

/**
 * What the change list reads: the autosomes. The sex chromosomes are read on their own line
 * (sexLine); a chrX row reaches the list only where a loaded parent's X reads absent with no chrY
 * called (parentage.ts), which no record of the answer key or the laboratory rows shows.
 */
export const LIST_SCOPE = 'chromosomes 1 to 22'

/** The count line used under every headline. */
export function changesText(l: ChangeList): string {
  if (l.notReported) return 'Not reportable, so no change is reported.'
  const { changes: n, findings } = cardCounts(l)
  const once = l.once.length
  // Scoped to what the list reads, so it cannot contradict the sex-chromosome line beside it.
  return `${n ? `${plural(n, 'change')} found` : 'No gain or loss found'} on ${LIST_SCOPE}${findings ? `, ${plural(findings, 'finding')}` : ''}`
    + `${once ? `; ${once} more stretch${once === 1 ? '' : 'es'} seen once, not confirmed` : ''}.`
}

const mb1 = (bp: number): string => (bp / 1e6).toFixed(1)

/**
 * When it arose, with the units it was read on, where units of the same embryo were read. A change
 * whose stretches read differently gets one line per stretch, a stretch with no word included, so
 * no word reads as covering the whole change.
 */
export const whenLines = (g: ChangeRow): string[] => [...new Set(g.pieces.flatMap((p) => (p.mechanismStretches
  ? p.mechanismStretches.map((s) => `When it arose, ${stretchWhere(p, s)}: ${stretchWord(s)}. ${cap(s.why)}.`)
  : p.mechanism ? [`When it arose: ${stretchWord({ mechanism: p.mechanism })}. ${cap(p.mechanismWhy ?? '')}.`] : [])))]

const cap = (x: string): string => x.charAt(0).toUpperCase() + x.slice(1)

/** Whose, in words. */
export function whoseText(g: ChangeRow): string {
  // A row on a sample not reportable states what was seen against the rest of the array, never a
  // copy lost or gained: the notice on such a card says the cell's own level is not established.
  if (g.word === 'not reported') return SEEN[rowDirection(g)]
  const who = g.origin === 'paternal' ? "father's" : g.origin === 'maternal' ? "mother's" : ''
  if (g.call === 'imbalance') return 'parents unbalanced'
  // THE TAXONOMY'S PLAIN PHRASE, NOT ITS CODE WORD. A row from the taxonomy printed its class as
  // the code spells it, so a copy-neutral stretch read "segmental upd" in a column whose other
  // rows read "the father's copy missing". defects.KIND_WORD already holds a phrase for every
  // class and is what the rest of the page prints.
  // A row with neither a gain nor a loss is named by its class, with no word for a copy lost.
  if (g.call === 'finding' || rowDirection(g) === 'none') {
    const kind = g.pieces[0].kind ?? ''
    // With the parent where the row has one, so the row does not withhold what its own figure
    // prints: the figure's heading is defects.headline, which names it.
    return `${KIND_WORD[kind] ?? kind.replace(/-/g, ' ')}${who ? `, the ${who}` : ''}`
  }
  if (g.gain) return who ? `extra copy, the ${who}` : 'extra copy, parent not named'
  if (g.deeper) return 'second copy lost too'
  if ((g.call === 'doubled' || g.pieces[0].kind === 'copy-loss') && !who) return 'copy lost, parent not named'
  return who ? `${who} copy missing` : "a parent's copy missing, parent not named"
}

/**
 * The title over one change's own figure: where it is and whose copy, in one sentence.
 *
 * IT READS THE ROW, so the two cannot name different parents. The figure used to compose its own
 * sentence and fell back to whichever parent was loaded where the call named none, so a row reading
 * "a parent's copy missing, parent not named" drew a figure under it headed "the father's copy is
 * missing", on the strength of the father being the array in the slot. Which track that figure
 * draws is a different question and still follows the loaded array: the track is the evidence the
 * call was made on, whoever turns out to own the missing copy, and it carries its own label.
 */
export function callTitle(g: ChangeRow): string {
  // A whole row's end is the last marker only where a parent was loaded (0 with none), so a whole
  // chromosome prints no span. A row read by stretch says so instead of carrying one origin word.
  const stretched = g.pieces.some((p) => p.mechanismStretches && split(p.mechanismStretches))
  const where = g.chrom === 'genome' ? 'The whole genome'
    : g.whole ? `Whole chromosome ${g.chrom}`
      : `Chromosome ${g.chrom}, ${mb1(g.startBp)} to ${mb1(g.endBp)} Mb `
        + `(${((g.endBp - g.startBp) / 1e6).toFixed(1)} Mb)`
  return `${where}: ${whoseText(g)}${stretched ? ', origin read by stretch' : ''}.`
}

/** Where, in Mb, one entry per joined piece; the refined coordinates every surface prints. */
export function whereText(g: ChangeRow): string {
  if (g.chrom === 'genome') return 'whole genome'
  if (g.whole) return `chr${g.chrom}, whole chromosome`
  return `chr${g.chrom} ${[...new Set(g.pieces.map((p) => `${mb1(p.startBp)}-${mb1(p.endBp)}`))].join(', ')} Mb`
}

/** The span of a row in Mb: the length its pieces cover, overlaps counted once. */
export const spanMb = (g: ChangeRow): string => {
  let total = 0; let end = -Infinity
  for (const p of [...g.pieces].sort((a, b) => a.startBp - b.startBp)) {
    total += Math.max(0, p.endBp - Math.max(p.startBp, end))
    end = Math.max(end, p.endBp)
  }
  return (total / 1e6).toFixed(1)
}

/** The scan's own pair of numbers for one piece, or nothing where the call has no pair. */
export function pairText(d: Defect): string {
  if (d.rate === undefined || d.nullRate === undefined) return ''
  return `missing-copy signal ${pct(d.rate, 1)} inside, ${pct(d.nullRate, 1)} with both copies present`
}

/** Intensity across one piece, as the scan measured it. */
export const lrrText = (d: Defect): string => (d.lrr === undefined ? ''
  : `intensity ${d.lrr > 0 ? '+' : ''}${d.lrr.toFixed(2)} log2`)

/** The full coordinates of one piece: locus and breakpoint intervals. */
export const locusText = (d: Defect): string => `${d.locus}${d.interval ? `, edges ${d.interval}` : ''}`

/** A stretch's own pair of numbers in plain words, for the line under its figure. */
export function saidOf(d: Defect): string | undefined {
  if (d.kind !== 'parental-absence' || d.rate === undefined || d.nullRate === undefined) return undefined
  const who = d.origin === 'paternal' ? "The father's" : d.origin === 'maternal' ? "The mother's" : "A parent's"
  return `${who} copy reads missing at ${pct(d.rate, 1)} of informative markers here, against ${pct(d.nullRate, 1)} with both copies present.`
}

/** The cells of one row of the PDF's change table: where, whose, word, then coordinates and signal. */
export function pdfCells(g: ChangeRow, setup: Setup): [string, string, string, string] {
  const seen = g.word === 'not reported'
  const pieces = g.whole ? '' : g.pieces.map((p) => [seen ? seenLocus(p) : `${p.locus}${p.interval ? ` (edges ${p.interval})` : ''}`,
    seen ? '' : pairText(p), lrrText(p)].filter(Boolean).join(', ')).join('; ')
  return [whereText(g), whoseText(g), g.word,
    [pieces, seen ? NOT_REPORTED_SAYS : measuredSentence(g.call, setup)].filter(Boolean).join('. ')]
}

/**
 * The line above the headline: which loaded arrays this sample was read against.
 *
 * THE NAME IS THE ONE THAT WAS LOADED. It used to be the sperm donor's slot whatever was in it, so
 * a run with the mother alone printed "<file> vs " and named nobody. Either parent alone is enough
 * everywhere else in this tool, and this line is the only place a reader can see which array the
 * card's words are held against, which is what the one-parent caveat on a headline depends on.
 */
export function readAgainst(
  r: { noParent?: boolean; role: Parent }, bothLoaded: boolean, donorName: string, oocyteName: string,
): string {
  if (r.noParent) return 'no parent loaded'
  if (bothLoaded) return `vs ${donorName} and ${oocyteName}`
  return `vs ${r.role === 'maternal' ? oocyteName : donorName}`
}

/** The sex-chromosome line the page and the PDF print: the constitution and its reason. */
export const sexLine = (sc: { constitution: string; why: string }): string =>
  (sc.constitution === 'unresolved' ? `Not read: ${sc.why}.` : `${sc.constitution}: ${sc.why}.`)

/** One line per row, the same text on the card, in the PDF and in the check. */
export function rowLine(g: ChangeRow): string {
  return `${whereText(g)} · ${whoseText(g)} · ${g.word}`
}

const NUM = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten']
const count = (n: number, one: string, many: string): string => `${n <= 10 ? NUM[n] : n} ${n === 1 ? one : many}`

/**
 * The whole-cell figure's title: what the figure shows, from the same headline and change list as
 * the card, never how a bar is computed.
 */
export function figureTitle(head: Headline, l: ChangeList, h: HeadInput): string {
  if (l.notReported) {
    return `${h.noParent ? 'Copy number along every chromosome' : head.line.split('. ')[0].replace(/\.$/, '')}. `
      + 'This sample is not reportable, so nothing marked is reported.'
  }
  const { changes: n, findings: f } = cardCounts(l)
  const allPieces = n > 0 && l.counted.every((g) => !g.whole || rowDirection(g) === 'none')
  const findings = f ? count(f, 'finding', 'findings') : ''
  const join = (a: string) => [a, findings].filter(Boolean).join(' and ')
  const what = join(n ? count(n, allPieces ? 'stretch' : 'change', allPieces ? 'stretches' : 'changes') : '')
  const except = what ? `, except ${what}` : ''
  const none = n ? '' : ` No gain or loss was found on ${LIST_SCOPE}.`
  const marked = join(n ? count(n, 'change', 'changes') : '')
  if (h.noParent) return `Copy number along every chromosome${marked ? `, with ${marked} marked` : ''}.${none}`
  if (head.word === 'BOTH PARENTS') return `Both parents' copies are present on every chromosome${except}.${none}`
  if (head.word.endsWith('+ SECOND PARENT')) {
    const who = head.word.startsWith('MOTHER') ? "mother's" : "father's"
    return `The ${who} copy and a second parent's copy are present on every chromosome${except}.${none}`
  }
  if (head.word === 'UNCLEAR') {
    return `Whose copies this cell carries is unclear${n ? `; ${count(n, 'change is', 'changes are')} marked` : ''}.`
  }
  return `${head.line.split('. ')[0].replace(/\.$/, '')}${n ? `; ${count(n, 'change is', 'changes are')} marked` : ''}.`
}

/* --- one sample, as every surface says it --------------------------------------------------- */

export interface CardSays { setup: Setup; list: ChangeList; head: Headline; defects: Defect[]; title: string }

/**
 * The headline and the change list of one scored sample: what the card, the PDF and the summary
 * table print. `cls` is the pair's class with both parents loaded; `rebuilt` that the loaded father
 * is a reconstruction.
 */
export function cardSays(r: ParentageResult, m: ParentageResult | undefined, cls: string, rebuilt: boolean): CardSays {
  const setup = setupOf(r, m, rebuilt)
  const defects = defectsForResult(r)
  const list = changeList(defects, setup, reportOf(r) !== 'not reportable')
  const h = headInput(r, m, cls)
  const head = headlineOf(h, changesText(list))
  return { setup, list, head, defects, title: figureTitle(head, list, h) }
}

/**
 * Where the card, the PDF and the figure could say different things about one scored sample: an
 * empty list means one count, one set of coordinates and one pair of numbers on all of them.
 * Run over real arrays by audit/truth/run.ts with OM_SURFACES=1.
 */
export function surfaceProblems(r: ParentageResult, m: ParentageResult | undefined, cls: string, rebuilt: boolean): string[] {
  const out: string[] = []
  const { list, head, setup, title } = cardSays(r, m, cls, rebuilt)
  const n = cardCounts(list).changes
  if (n && !head.line.includes(`${n} change`)) out.push(`headline line does not carry the count ${n}: ${head.line}`)
  if (!n && head.line.includes(' changes found')) out.push(`headline counts changes the list does not: ${head.line}`)
  if (n && (head.word === 'BOTH PARENTS' || head.word.endsWith('SECOND PARENT')) && !title.includes(n <= 10 ? NUM[n] : String(n))) {
    out.push(`figure title does not carry the count ${n}: ${title}`)
  }
  for (const g of [...list.counted, ...list.once, ...list.seen, ...list.notReported?.rows ?? []]) {
    const [where, , word, detail] = pdfCells(g, setup)
    if (where !== whereText(g)) out.push(`PDF where differs: ${where} / ${whereText(g)}`)
    if (word !== g.word || !['found two ways', 'consistent', 'one reading', 'also seen', 'position only', 'not reported'].includes(word)) {
      out.push(`word: ${word}`)
    }
    for (const p of g.whole ? [] : g.pieces) {
      if (!detail.includes(p.locus)) out.push(`PDF lacks the locus ${p.locus}`)
      const pair = pairText(p).match(/[\d.]+%/g) ?? []
      if (!pair.every((x) => detail.includes(x))) out.push(`PDF pair differs at ${p.locus}`)
      const said = saidOf(p)
      if (said && !pair.every((x) => said.includes(x))) out.push(`figure pair differs at ${p.locus}`)
      const sg = r.segments.find((x) => x.chrom === p.chrom && segmentCoords(x).start === p.startBp && segmentCoords(x).end === p.endBp)
      if (p.step && p.step !== 'finding' && !p.step.startsWith('whole') && p.step !== 'imbalance' && !sg) {
        out.push(`piece ${p.locus} is not a segment's refined coordinates`)
      }
    }
  }
  return out
}

/* --- the change list of a scored result ------------------------------------------------------ */

/**
 * An F-grade call for a finding that covers the whole genome rather than an interval.
 *
 * There is no interval to measure, so the evidence is the array-wide absence of the loaded parent,
 * which is the right evidence for a genome-scoped statement and the wrong evidence for a located
 * one. Graded F because it is a direction and not a measurement of anything in particular.
 */
const genomeGrade = (r: {
  genomeRate?: number; explainable?: number
}): { verdict: string; shift: number; z: number; impliedF: number; why: string
  confidence: number; band: string } | undefined => {
  const { genomeRate, explainable } = r
  if (!Number.isFinite(genomeRate) || !Number.isFinite(explainable)) return undefined
  return {
    verdict: (genomeRate as number) > (explainable as number) ? 'other-parent' : 'loaded-parent',
    // No interval, so there is no interval statistic. Left absent rather than filled with the
    // array-wide numbers, which would read as a measurement of this finding.
    shift: NaN,
    z: NaN,
    impliedF: NaN,
    confidence: 0.5,
    band: 'F',
    why: 'Graded F and named anyway: this finding covers the whole genome and has no interval to '
      + `measure, so the direction is the array's own, from absence of the loaded parent at `
      + `${((genomeRate as number) * 100).toFixed(2)}% against a `
      + `${((explainable as number) * 100).toFixed(2)}% ceiling. Do not report or count this row`,
  }
}

/**
 * The change list for one sample, built once so the panel and the export cannot disagree.
 *
 * Previously assembled inline in the render, which meant any machine-readable output would have
 * had to reassemble it and could have drifted from what the reader was looking at.
 */
/** The last marker position seen on a chromosome, for a whole-chromosome event's extent. */
function chromEndOf(r: ParentageResult, chrom: string): number {
  const pos = r.markerPositions?.get(chrom)
  return pos && pos.length ? pos[pos.length - 1] : 0
}

export function defectsForResult(r: ParentageResult): Defect[] {
  // WHOLE-CHROMOSOME ANEUPLOIDY IS A CHROMOSOMAL CHANGE, and belongs in the one list with the rest.
  // It used to have a panel of its own, so a reader comparing a whole-chromosome loss against a
  // segmental one was comparing two boxes with different vocabularies and no shared ordering.
  // AND A CHROMOSOME THE ALLELES READ AS MISSING ONE PARENT, with no intensity verdict: that is a
  // parental absence rather than a copy loss, because one parent's allele alone is what a
  // monosomy and an isodisomy both read. See ChromResult.lostParent.
  // AND AN EXTRA COPY OR AN IMBALANCE THE ALLELES FOUND: see ChromResult.alleleGain. A gain the
  // intensity detector also flagged keeps its row, and the alleles supply the parent where they
  // name one the dosage channel does not or name it better; where the two name different parents
  // the row names nobody.
  const alleleEvent = (c: { alleleGain?: { verdict: string } }) =>
    c.alleleGain?.verdict === 'gain' || c.alleleGain?.verdict === 'imbalance'
  const aneu = (r.chroms ?? []).filter((c) => c.aneuploidy || c.lostParent || c.mixed
    || alleleEvent(c)).map((c) => {
    const scored = (r.dosageCalls ?? []).find((d: { where: string }) => d.where === `chr${c.chrom}`)
    const other = r.role === 'paternal' ? 'maternal' : 'paternal'
    // THE MENDELIAN CHANNEL FIRST, measured at 0.920 on obvious whole-chromosome events against
    // the dosage channel's band D on the same material. It asks whether an allele is there at all
    // rather than whether a mean has moved, so it needs no detection floor.
    const mendel = (r.oneParent ?? []).find((o: { where: string }) => o.where === `chr${c.chrom}`)
    const named = mendel ? mendelParent(mendel, r.role, c.aneuploidy === 'gain') : null
    // THE ALLELES ALONE, with two bulk parents loaded: see ChromResult.lostParent. Where the
    // Mendelian channel also names a parent the two must agree, or the row names nobody. Two
    // readings of one genome that disagree are a fact about the genome, not a tie to break.
    // On a GAIN the origin is the parent who supplied the extra copy. One parent's allele alone on
    // a three-copy chromosome means every copy is that parent's, so the absent parent is inverted.
    const letters = !c.lostParent ? null
      : c.aneuploidy === 'gain' ? (c.lostParent === 'paternal' ? 'maternal' : 'paternal')
        : c.lostParent
    const ag = c.alleleGain
    const gainParent = ag?.verdict === 'gain' ? ag.parent ?? null : null
    const dosageParent = scored ? parentNamed(scored.verdict as DosageVerdict, r.role) : null
    const gainConflict = !!(gainParent && c.aneuploidy === 'gain' && dosageParent
      && dosageParent !== gainParent)
    const conflict = !!(named && letters && named !== letters)
    const origin = conflict || c.mixed || gainConflict
      || (ag?.verdict === 'imbalance' && !c.aneuploidy && !c.lostParent) ? 'unclear'
      : gainParent && !letters ? gainParent
      : named ?? letters
      ?? (c.aneuploidyParent === 'this' ? r.role
        : c.aneuploidyParent === 'other' ? other
          : (scored && parentNamed(scored.verdict as DosageVerdict, r.role)) || 'unclear')
    const k = c.contribution
    const [father, mother] = k
      ? (r.role === 'paternal' ? [k.onlyA, k.onlyB] : [k.onlyB, k.onlyA]) : [NaN, NaN]
    const notOneOrTwo = c.aneuploidy ? '' : ' Whether one copy remains or two is a question for the intensity, '
      + 'which did not reach a verdict on this chromosome: a monosomy and an isodisomy read '
      + 'the same at these markers.'
    // The opposite-homozygote count needs both parents' arrays; with one parent the call rests on the
    // copy-state score (oneParentState.ts), so that is what the row quotes, or nothing.
    const loadedWho = r.role === 'paternal' ? 'father' : 'mother'
    const otherWho = r.role === 'paternal' ? 'mother' : 'father'
    const scoreHere = c.copyState?.r
    const scoreRest = (r.chroms ?? []).filter((x) => x.chrom !== c.chrom && /^\d+$/.test(x.chrom))
      .map((x) => x.copyState?.r).filter((x): x is number => x !== undefined && Number.isFinite(x)).sort((a, b) => a - b)
    const lettersWhy = !letters ? ''
      : r.twoParents && k
        ? ` Of ${int(k.informative)} markers where the parents are opposite homozygotes, `
          + `${pct(k.fraction, 1)} read both alleles, ${pct(father, 1)} the father's alone and `
          + `${pct(mother, 1)} the mother's alone, so `
          + (c.aneuploidy === 'gain' ? `every copy present is ${letters}.`
            : `the ${c.lostParent} copy is the one missing.`)
          + notOneOrTwo
        : !r.twoParents && c.aneuploidy !== 'gain' && scoreHere !== undefined && Number.isFinite(scoreHere) && scoreRest.length
          ? ` Read with the ${loadedWho} alone: the alleles the ${loadedWho} lacks score `
            + `${scoreHere.toFixed(2)} on this chromosome, against a median of `
            + `${scoreRest[scoreRest.length >> 1].toFixed(2)} on this cell's other autosomes; near 1 means `
            + `the ${loadedWho}'s copy is gone and near 0 the ${otherWho}'s, so the ${c.lostParent === 'paternal' ? 'father' : 'mother'}'s `
            + 'copy is the one missing.' + notOneOrTwo
          : ''
    // WHY THE ALLELE CHANNEL IS SILENT ON A ROW THAT EXISTS. The intensity channel can reach a
    // verdict on a chromosome the alleles declined, and until now the row simply carried no allele
    // sentence: a reader comparing the two channels could not tell a chromosome with nothing to
    // read from one the array failed on. The refusal states which.
    const withheldWhy = !k && c.contributionWithheld
      ? ` The alleles reached no verdict here: ${c.contributionWithheld}.`
      : ''
    const conflictWhy = conflict
      ? ` The allele count and the Mendelian channel name different parents, so no parent is named.`
      : ''
    const mixedWhy = c.mixed && k
      ? ` Mixed signal, not called. Of ${int(k.informative)} markers where the parents are opposite `
        + `homozygotes, ${pct(k.fraction, 1)} read both alleles, ${pct(father, 1)} the father's alone `
        + `and ${pct(mother, 1)} the mother's alone: between the boundaries with one parent dominant. `
        + 'A mosaic, a partial loss or noise, not separable here.'
      : ''
    return {
      chrom: c.chrom,
      startBp: 0,
      endBp: chromEndOf(r, c.chrom),
      kind: c.aneuploidy === 'gain' || (!c.aneuploidy && !c.lostParent && gainParent) ? 'copy-gain'
        : c.aneuploidy === 'loss' ? 'copy-loss'
          : !c.aneuploidy && !c.lostParent && ag?.verdict === 'imbalance' ? 'whole-chromosome'
            : 'parental-absence',
      locus: `chr${c.chrom}`,
      origin,
      band: gainParent && !named ? 'C' : named ? mendel?.band : scored?.band,
      confidence: named ? mendel?.posterior : gainParent ? undefined : scored?.confidence,
      inheritedMargin: (scored as { inheritedMargin?: number } | undefined)?.inheritedMargin,
      stage: r.stage?.stage,
      step: !c.aneuploidy && !c.lostParent && !gainParent ? 'imbalance'
        : c.lostParent || named || gainParent ? 'whole-alleles' : 'whole-intensity',
      why: (!c.aneuploidy && !c.lostParent && !c.mixed && ag
        ? `${ag.why[0].toUpperCase()}${ag.why.slice(1)}.`
        : c.aneuploidy
        ? `calls at ${c.callFraction.toFixed(2)}x the genome rate with intensity `
          + `${c.lrrShift > 0 ? '+' : ''}${c.lrrShift.toFixed(2)} log2 from the rest. An intact `
          + 'chromosome calls at 0.78x to 1.16x of its genome median and never leaves -0.79 to '
          + '+0.42 log2, measured over 1,012 chromosomes.'
        : `calls at ${c.callFraction.toFixed(2)}x the genome rate with intensity `
          + `${c.lrrShift > 0 ? '+' : ''}${c.lrrShift.toFixed(2)} log2 from the rest, neither of `
          + 'which reached a whole-chromosome verdict; the alleles did.')
        + `${lettersWhy}${withheldWhy}${conflictWhy}${mixedWhy}`
        + (c.aneuploidy === 'gain' && ag && alleleEvent(c) ? ` By the alleles: ${ag.why}.` : '')
        + (gainConflict ? ' The alleles and the dosage channel name different parents, so no parent is named.' : '')
        + `${named ? ` ${mendel?.why}` : ''}`
        + (scored?.why && !named && !letters && !gainParent ? ` ${scored.why}` : ''),
    } as never
  })

  return withMechanism([
    ...aneu,
      ...defectsFrom(r.segments, r.gains, r.losses ?? [], r.oneParent ?? [], r.role,
        r.stage, r.dosageCalls ?? []),
      // One list, not two. A copy-neutral event and a deletion are different measurements
      // of the same kind of thing and a reader compares them against each other.
      ...(r.findings ?? []).map((f) => findingToDefect(f, r.stage,
        // The call scored over THIS finding's own interval, matched by the same label the
        // scorer wrote, so a finding shows the origin that was actually measured for it
        // rather than one borrowed from its chromosome.
        // A GENOME-SCOPED FINDING HAS NO INTERVAL TO SCORE, so it never reached the
        // scorer and was the one row still leaving without a grade. Its evidence is the
        // whole array's, which is exactly what a genome-scoped finding should be judged
        // on, so it is graded F from the same array-wide absence measurement the other
        // ungradable rows fall back to.
        f.chrom === 'genome'
          ? genomeGrade(r)
          : (r.dosageCalls ?? []).find((d: { where: string }) => d.where
            === `chr${f.chrom} ${(f.startBp / 1e6).toFixed(1)}-${(f.endBp / 1e6).toFixed(1)}Mb`),
        r.role)),
    ], r.uniformity)
}

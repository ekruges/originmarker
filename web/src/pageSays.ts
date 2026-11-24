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
import { ORIGIN_WORD, split, stretchWhere, stretchWord } from './timing.ts'

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
export function changeList(defects: readonly Defect[], setup: Setup): ChangeList {
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

/** The disturbed-genome finding in one plain sentence, or nothing. */
export function disturbedText(l: ChangeList): string {
  if (!l.disturbed) return ''
  const { on, of } = l.disturbed
  return `${Number.isFinite(on) ? `${on} of ${of}` : 'Many'} chromosomes carry a change: this cell's genome is widely disturbed.`
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

/** The count line used under every headline. */
export function changesText(l: ChangeList): string {
  const n = l.counted.length
  const once = l.once.length
  return `${n ? `${n} change${n === 1 ? '' : 's'} found` : 'No change found'}`
    + `${once ? `; ${once} more stretch${once === 1 ? '' : 'es'} seen once, not confirmed` : ''}.`
}

const mb1 = (bp: number): string => (bp / 1e6).toFixed(1)

/**
 * When it arose, with the units it was read on, where units of the same embryo were read. A change
 * whose stretches read differently gets one line per stretch, a stretch with no word included, so
 * no word reads as covering the whole change.
 */
export const whenLines = (g: ChangeRow): string[] => [...new Set(g.pieces.flatMap((p) => (p.mechanismStretches
  ? p.mechanismStretches.map((s) => `When it arose, ${stretchWhere(p, s)}: ${stretchWord(s)}. ${s.why}.`)
  : p.mechanism && ORIGIN_WORD[p.mechanism] ? [`When it arose: ${ORIGIN_WORD[p.mechanism]}. ${p.mechanismWhy ?? ''}.`] : [])))]

/** Whose, in words. */
export function whoseText(g: ChangeRow): string {
  const who = g.origin === 'paternal' ? "father's" : g.origin === 'maternal' ? "mother's" : ''
  if (g.call === 'imbalance') return 'parents unbalanced'
  // THE TAXONOMY'S PLAIN PHRASE, NOT ITS CODE WORD. A row from the taxonomy printed its class as
  // the code spells it, so a copy-neutral stretch read "segmental upd" in a column whose other
  // rows read "the father's copy missing". defects.KIND_WORD already holds a phrase for every
  // class and is what the rest of the page prints.
  if (g.call === 'finding') {
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
  const pieces = g.whole ? '' : g.pieces.map((p) => [`${p.locus}${p.interval ? ` (edges ${p.interval})` : ''}`,
    pairText(p), lrrText(p)].filter(Boolean).join(', ')).join('; ')
  return [whereText(g), whoseText(g), g.word, [pieces, measuredSentence(g.call, setup)].filter(Boolean).join('. ')]
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
  const n = l.counted.length
  const allPieces = n > 0 && l.counted.every((g) => !g.whole)
  const except = n ? `, except ${count(n, allPieces ? 'stretch' : 'change', allPieces ? 'stretches' : 'changes')}` : ''
  const none = n ? '' : ' No change was found.'
  if (h.noParent) return `Copy number along every chromosome${n ? `, with ${count(n, 'change', 'changes')} marked` : ''}.${none}`
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
  const list = changeList(defects, setup)
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
  const n = list.counted.length
  if (n && !head.line.includes(`${n} change`)) out.push(`headline line does not carry the count ${n}: ${head.line}`)
  if (!n && head.line.includes(' changes found')) out.push(`headline counts changes the list does not: ${head.line}`)
  if (n && (head.word === 'BOTH PARENTS' || head.word.endsWith('SECOND PARENT')) && !title.includes(n <= 10 ? NUM[n] : String(n))) {
    out.push(`figure title does not carry the count ${n}: ${title}`)
  }
  for (const g of [...list.counted, ...list.once, ...list.seen]) {
    const [where, , word, detail] = pdfCells(g, setup)
    if (where !== whereText(g)) out.push(`PDF where differs: ${where} / ${whereText(g)}`)
    if (word !== g.word || !['found two ways', 'consistent', 'one reading', 'also seen', 'position only'].includes(word)) {
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

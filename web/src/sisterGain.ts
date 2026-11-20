/**
 * A gain read where a sister cell shows the loss.
 *
 * WHY. A copy that one daughter cell lacks went to the other, so the two changes sit over one interval:
 * a whole chromosome that failed to separate, or the piece beyond a break. The loss is the easy half.
 * It stops the array genotyping, and with a parent the alleles name it. The gain is the hard half: a
 * third copy moves intensity by a fraction of what a lost one does, and a scan of the whole genome for
 * it has to clear a threshold set for every place it might be. Told where to look, it is one test.
 *
 * THE READING. Intensity is taken as the median log R of each 5 Mb window. For a cell and an interval
 * shown lost by another cell of the same embryo:
 *
 *   shift = median of the cell's windows in the interval, less the median of its windows on the
 *           chromosomes that share its commonest copy state (parentage.copyReference)
 *   z     = shift over the spread of runs of as many windows on those chromosomes
 *
 * A run of windows is the unit because amplification noise is correlated along a chromosome: single
 * windows vary far more than their count suggests (intensityNull.ts). The tested chromosome is left
 * out of its own null.
 *
 * MEASURED by audit/truth/sisters.ts on 257 open cells of 72 embryos, no parent loaded, on chromosomes
 * the cell did not already report a change on. Where a sister lost the whole chromosome: 5 of 11
 * chromosomes stated gained are read, 0 of 38 stated normal and 0 of 11 stated lost. The same intervals
 * put on cells of OTHER embryos, where nothing reciprocal can be, are read on 4 of 917 chromosomes
 * stated normal. Losses with a breakpoint are few without a parent, 6 stated-normal tests and none
 * read, so the segment case is not yet measured.
 *
 * WHICH CELLS COUNT AS SISTERS. The embryo label is not evidence: of 431 pairs sharing a label on the
 * laboratory's own arrays only 129 read as one genome, and 29 of the draft's 30 calls came from a
 * "sister" that was a different genome. A pair counts only where the two arrays read as ONE genome, by
 * relatedness.oppositeHom at SISTER_OPP_HOM_MAX, with SISTER_MIN_SHARED shared homozygous calls needed
 * before the rate is read at all. An unmeasurable pair is not a sister.
 *
 * WHY NOT parentage.agreement, WHICH IS THE PAGE'S OWN DEFINITION OF A UNIT OF ONE EMBRYO. Identity of
 * calls is FIRST order in dropout: one cell dropping an allele breaks it. An opposite homozygote is
 * SECOND order: it needs both cells to drop at the same heterozygous marker in opposite directions. On
 * these single cells, where dropout runs at a median of 0.324, that difference decides the measure.
 * Measured over 14,365 pairs of 170 open arrays (audit/truth/sister-gain-2026-11-19.txt): identity
 * correlates +0.772 with the pair's lower call rate and separates same-label from different-label pairs
 * by 1.12x in the median, while the opposite-homozygote rate correlates -0.445 and separates them by
 * 3.3x. abnormalities.SAME_EMBRYO_AGREEMENT is 0.90 because two TROPHECTODERM BIOPSIES concord at
 * 95.8%; two blastomeres of one embryo sit at 0.68 to 0.85, under that threshold, so identity cannot
 * carry this rule on the material it is for.
 *
 * SO THIS RULE AND THE PAGE'S UNIT GROUPING USE DIFFERENT MEASURES, which is a thing to know rather
 * than a thing to like. The grouping drives the gamete-borne against post-zygotic call and fails
 * closed on blastomeres: it reports fewer units and withholds, rather than grouping wrongly.
 *
 * WHAT IT IS NOT. Evidence about a cell with no sister on the array, and a claim about which division
 * made the pair. A cell is only asked about intervals a sister shows lost; the gain it reports stands
 * on its own intensity.
 */
import { oppositeHom } from './relatedness.ts'
import type { AB } from './informativity.ts'

export const WINDOW_BP = 5_000_000
/** Markers with an intensity a window needs before its median is kept. */
export const MIN_WINDOW_MARKERS = 100

/** One window's median log R: chromosome, window index (position / WINDOW_BP, floored), median. */
export type LrrWindow = [chrom: string, index: number, median: number]

const medianOf = (xs: number[]): number => {
  if (!xs.length) return NaN
  xs.sort((a, b) => a - b)
  const m = xs.length >> 1
  return xs.length % 2 ? xs[m] : (xs[m - 1] + xs[m]) / 2
}

export function lrrWindows(
  byChrom: ReadonlyMap<string, readonly { pos: number; log2R: number | null }[]>,
): LrrWindow[] {
  const out: LrrWindow[] = []
  for (const [chrom, ms] of byChrom) {
    const bins = new Map<number, number[]>()
    for (const m of ms) {
      if (m.log2R === null || !Number.isFinite(m.log2R)) continue
      const i = Math.floor(m.pos / WINDOW_BP)
      const b = bins.get(i)
      if (b) b.push(m.log2R)
      else bins.set(i, [m.log2R])
    }
    for (const [i, xs] of [...bins].sort((a, b) => a[0] - b[0])) {
      if (xs.length >= MIN_WINDOW_MARKERS) out.push([chrom, i, medianOf(xs)])
    }
  }
  return out
}

export interface IntervalReading { shift: number; z: number; windows: number }

/** Reference chromosomes and runs needed before an interval is read at all. */
export const MIN_REFERENCE_CHROMS = 8
export const MIN_NULL_RUNS = 20
/** Windows an interval must cover: two, 10 Mb, so that one window's noise is not the reading. */
export const MIN_INTERVAL_WINDOWS = 2
/**
 * The longest run the null is made of. A whole chromosome 1 is fifty windows and no other chromosome
 * holds twenty runs that long, so a longer interval is read against runs of this many. A longer run
 * varies less, so the null is then too wide, never too narrow: it errs toward not reading.
 */
export const MAX_RUN_WINDOWS = 12

/**
 * One interval of one cell against that cell's own reference chromosomes. Null where the interval
 * covers too few windows or the cell has too little reference to be read against.
 */
export function readInterval(
  windows: readonly LrrWindow[], reference: readonly string[],
  chrom: string, startBp: number, endBp: number,
): IntervalReading | null {
  const by = new Map<string, Map<number, number>>()
  for (const [c, i, v] of windows) {
    let m = by.get(c)
    if (!m) by.set(c, m = new Map())
    m.set(i, v)
  }
  const own = by.get(chrom)
  if (!own || reference.length < MIN_REFERENCE_CHROMS) return null
  // A window belongs to the interval when all of it lies inside, give or take half a window.
  const inside = [...own.keys()].sort((a, b) => a - b)
    .filter((i) => i * WINDOW_BP >= startBp - WINDOW_BP / 2 && (i + 1) * WINDOW_BP <= endBp + WINDOW_BP / 2)
  const k = inside.length
  if (k < MIN_INTERVAL_WINDOWS) return null
  const run = Math.min(k, MAX_RUN_WINDOWS)
  const refValues: number[] = []
  const runs: number[] = []
  for (const c of reference) {
    const m = by.get(c)
    if (!m) continue
    for (const v of m.values()) refValues.push(v)
    if (c === chrom) continue
    const keys = [...m.keys()].sort((a, b) => a - b)
    const step = Math.max(1, run >> 1)
    for (let j = 0; j + run <= keys.length; j += step) {
      runs.push(medianOf(keys.slice(j, j + run).map((i) => m.get(i) as number)))
    }
  }
  if (runs.length < MIN_NULL_RUNS) return null
  const centre = medianOf(refValues)
  const runCentre = medianOf([...runs])
  const scale = 1.4826 * medianOf(runs.map((x) => Math.abs(x - runCentre)))
  const shift = medianOf(inside.map((i) => own.get(i) as number)) - centre
  return { shift, z: scale > 0 ? shift / scale : 0, windows: k }
}

/** Spreads of the cell's own runs, and log2 above its reference, before a sister's loss is read as a gain. */
export const SISTER_Z = 3
export const SISTER_SHIFT = 0.25


export interface Interval { chrom: string; startBp: number; endBp: number; whole: boolean }

export interface SisterUnit {
  windows: readonly LrrWindow[]
  /** Genotype calls by probe id. Required: a pair whose genomes cannot be compared is not a pair. */
  gt: ReadonlyMap<string, AB>
  /** Chromosomes sharing the cell's commonest copy state, with no copy change reported on them. */
  reference: readonly string[]
  /** What this cell shows lost: whole chromosomes and segments. */
  losses: readonly Interval[]
  /** Chromosomes this cell already reports a copy change on, which are not asked again. */
  changed: ReadonlySet<string>
}

export interface SisterGain extends Interval, IntervalReading {
  /** Index, in the group, of the cell that shows the loss. */
  from: number
}

/**
 * For each cell of one embryo, the intervals a sister shows lost that this cell reads as gained.
 * Where several sisters show a loss on one chromosome the strongest reading is kept.
 */
/**
 * Opposite-homozygote rate at or under which two arrays read as one genome.
 *
 * Measured on these single cells: pairs sharing a laboratory label sit at a median of 0.0366 and pairs
 * that do not at 0.1208, and the groups the plan's own earlier measurement read as one genome sat at
 * 0.000 to 0.012 against 0.027 for the closest pair of different genomes.
 *
 * This rule's own number. relatedness.OPPOSITE_HOM_MAX happens to be the same 0.020 and answers a
 * different question, parent against stranger; retuning that one must not move this one.
 */
export const SISTER_OPP_HOM_MAX = 0.020
/**
 * Shared homozygous calls a pair needs before its rate is read at all, as relatedness.relate requires.
 *
 * A pair under it is NOT a sister. The fault this rule exists to stop is reading a gain off a cell of
 * another genome, so an unmeasurable pair must not supply an interval.
 */
export const SISTER_MIN_SHARED = 10_000

/** Whether two arrays read as one genome. Unmeasurable is not a sister. */
export function oneGenome(a: SisterUnit, b: SisterUnit): boolean {
  const { rate, n } = oppositeHom(a.gt as Map<string, AB>, b.gt as Map<string, AB>)
  return n >= SISTER_MIN_SHARED && Number.isFinite(rate) && rate <= SISTER_OPP_HOM_MAX
}

/** The parts of a scored cell this rule reads, so the page and audit/truth/sisters.ts derive a unit
 *  the one way. `changed` is every chromosome already carrying a reported change, which is not asked
 *  again, and the reference drops those chromosomes. */
export interface ScoredCell {
  chroms: readonly { chrom: string; aneuploidy?: 'loss' | 'gain' | null }[]
  segments: readonly { chrom: string; kind: string; start: number; end: number }[]
  copyReference?: readonly string[] | null
}

const AUTOSOMES = Array.from({ length: 22 }, (_, i) => String(i + 1))

export function sisterUnitOf(
  cell: ScoredCell, windows: readonly LrrWindow[], gt: ReadonlyMap<string, AB>,
): SisterUnit {
  const autosome = (c: string) => /^\d+$/.test(c)
  const whole = cell.chroms.filter((c) => c.aneuploidy && autosome(c.chrom))
  const segs = cell.segments.filter(
    (sg) => autosome(sg.chrom) && (sg.kind === 'copy-loss' || sg.kind === 'copy-gain'))
  const changed = new Set<string>([...whole.map((c) => c.chrom), ...segs.map((sg) => sg.chrom)])
  return {
    windows, gt, changed,
    reference: (cell.copyReference ?? AUTOSOMES).filter((c) => !changed.has(c)),
    losses: [
      ...whole.filter((c) => c.aneuploidy === 'loss').map((c) => WHOLE(c.chrom)),
      ...segs.filter((sg) => sg.kind === 'copy-loss')
        .map((sg) => ({ chrom: sg.chrom, startBp: sg.start, endBp: sg.end, whole: false })),
    ],
  }
}

export function sisterGains(units: readonly SisterUnit[]): SisterGain[][] {
  // Symmetric by construction, so each unordered pair is measured once.
  const pairs = new Map<string, boolean>()
  const same = (i: number, j: number): boolean => {
    const key = i < j ? `${i}:${j}` : `${j}:${i}`
    let v = pairs.get(key)
    if (v === undefined) pairs.set(key, v = oneGenome(units[i], units[j]))
    return v
  }
  return units.map((u, b) => {
    const best = new Map<string, SisterGain>()
    units.forEach((other, a) => {
      if (a === b || !same(a, b)) return
      for (const loss of other.losses) {
        if (u.changed.has(loss.chrom)) continue
        const r = readInterval(u.windows, u.reference, loss.chrom, loss.startBp, loss.endBp)
        if (!r || !(r.shift >= SISTER_SHIFT && r.z >= SISTER_Z)) continue
        const had = best.get(loss.chrom)
        if (!had || r.z > had.z) best.set(loss.chrom, { ...loss, ...r, from: a })
      }
    })
    return [...best.values()].sort((x, y) => Number(x.chrom) - Number(y.chrom))
  })
}

/** A whole chromosome as an interval. */
export const WHOLE = (chrom: string): Interval => ({ chrom, startBp: 0, endBp: Number.MAX_SAFE_INTEGER, whole: true })

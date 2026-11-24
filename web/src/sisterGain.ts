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
 * parentage.agreement at SISTER_AGREE_MIN, taken OFF the chromosomes either cell reports a change on.
 * agreement returns NaN under its own floor of shared calls, so an unmeasurable pair is not a sister.
 *
 * THE EXCLUSION IS THE POINT, NOT A DETAIL. Two cells of one embryo that lost opposite parental copies
 * of a region read as opposite homozygotes across that whole region (relatedness.ts: four of 122
 * same-embryo blastomere pairs, at 0.053 to 0.056). That is FIRST order and it is the very reciprocal
 * event this rule hunts, so a genome-wide comparison pushes exactly the pairs worth reading above any
 * gate, while two ordinary siblings with quiet genomes stay under it. Measured on the published cells,
 * the worst false pair genome-wide is two trophectoderm biopsies of different embryos of one egg donor.
 *
 * WHY NOT relatedness.oppositeHom, WHICH THE PLAN FIRST NAMED. It admits siblings at every threshold.
 * Measured against the 31 published cells' own embryo labels, 51 true sister pairs against 414 others:
 * an opposite-homozygote gate of 0.02 admits 23 sisters and 7 non-sisters, and its best threshold still
 * admits 5. relatedness.ts says why: over 106 published trios sibling pairs sit at 0.0114, under the
 * 0.020 gate, so the rate alone cannot make that split.
 *
 * WHY THE THRESHOLD IS NOT abnormalities.SAME_EMBRYO_AGREEMENT. That is 0.90 because two TROPHECTODERM
 * BIOPSIES concord at 95.8%. Blastomeres of one embryo sit lower, and on the published truth 0.90
 * admits 4 of 51 true sisters while 0.86 admits 16 with no false pair of the 414. The page keeps 0.90
 * for its unit grouping, which drives the gamete-borne against post-zygotic call; that grouping fails
 * closed on blastomeres, reporting fewer units rather than grouping wrongly, and its calibration is
 * its own item.
 *
 * WHAT WAS TRIED AND DOES NOT WORK. 8*oo*hh/hx^2 is one under a dropout model of a single genome and
 * is measurably free of dropout (correlation -0.064 with the implied rate, against -0.383 for
 * identity), but it carries no genome-sharing signal at all: 0.380 on true sisters against 0.345 on
 * the rest. It measures how the arrays behave, not whether they are one genome.
 *
 * WHAT IT IS NOT. Evidence about a cell with no sister on the array, and a claim about which division
 * made the pair. A cell is only asked about intervals a sister shows lost; the gain it reports stands
 * on its own intensity.
 */
import { agreement } from './parentage.ts'
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
  /** Where each probe sits, so the same-genome comparison can leave out the changed chromosomes.
   *  The shape parentage and relatedness already use, so no second map is built. A probe missing
   *  here is left out rather than compared. */
  pos: ReadonlyMap<string, { chrom: string }>
  /** Chromosomes sharing the cell's commonest copy state, with no copy change reported on them. */
  reference: readonly string[]
  /** What this cell shows lost: whole chromosomes and segments. */
  losses: readonly Interval[]
  /** Chromosomes this cell already reports a copy change on, which are not asked again. */
  changed: ReadonlySet<string>
  /** Staged haploid. A haploid cell is only ever the sister of a haploid cell (samePloidy). */
  haploid?: boolean
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
 * Share of shared called markers, off the changed chromosomes, above which two arrays read as one
 * genome.
 *
 * Calibrated on the 31 published cells of one series against their own embryo labels, 51 true sister
 * pairs and 414 others: no pair of different embryos reaches 0.8577, and this admits 16 of the 51
 * sisters with none of the 414. The cost is sensitivity, which is the right way to be wrong here,
 * because a pair that is not read leaves every call untouched while a pair wrongly admitted can print
 * a gain that is not there.
 *
 * This rule's own number, deliberately not abnormalities.SAME_EMBRYO_AGREEMENT.
 */
export const SISTER_AGREE_MIN = 0.86

/**
 * Whether two cells may be sisters at all by their stage: a cell staged haploid only with another
 * staged haploid, any other cell only with a cell not staged haploid. A pronucleus or polar body is
 * not a unit of a diploid embryo, and two haploid products of one parent agree at most markers
 * because both are homozygous everywhere, so agreement alone cannot refuse such a pair.
 */
export const samePloidy = (a: Pick<SisterUnit, 'haploid'>, b: Pick<SisterUnit, 'haploid'>): boolean =>
  !!a.haploid === !!b.haploid

/** Whether two arrays read as one genome, off the chromosomes either reports a change on, and are
 *  of one ploidy stage. Unmeasurable is not a sister: NaN fails the comparison. */
export function oneGenome(a: SisterUnit, b: SisterUnit): boolean {
  if (!samePloidy(a, b)) return false
  const skip = (probe: string): boolean => {
    const c = a.pos.get(probe)?.chrom
    return c === undefined || a.changed.has(c) || b.changed.has(c)
  }
  return agreement(a.gt as Map<string, AB>, b.gt as Map<string, AB>, skip) > SISTER_AGREE_MIN
}

/** The parts of a scored cell this rule reads, so the page and audit/truth/sisters.ts derive a unit
 *  the one way. `changed` is every chromosome already carrying a reported change, which is not asked
 *  again, and the reference drops those chromosomes. */
export interface ScoredCell {
  chroms: readonly { chrom: string; aneuploidy?: 'loss' | 'gain' | null }[]
  segments: readonly { chrom: string; kind: string; start: number; end: number }[]
  copyReference?: readonly string[] | null
  stage?: { stage?: string }
}

const AUTOSOMES = Array.from({ length: 22 }, (_, i) => String(i + 1))

export function sisterUnitOf(
  cell: ScoredCell, windows: readonly LrrWindow[], gt: ReadonlyMap<string, AB>,
  pos: ReadonlyMap<string, { chrom: string }>,
): SisterUnit {
  const autosome = (c: string) => /^\d+$/.test(c)
  const whole = cell.chroms.filter((c) => c.aneuploidy && autosome(c.chrom))
  const segs = cell.segments.filter(
    (sg) => autosome(sg.chrom) && (sg.kind === 'copy-loss' || sg.kind === 'copy-gain'))
  const changed = new Set<string>([...whole.map((c) => c.chrom), ...segs.map((sg) => sg.chrom)])
  return {
    windows, gt, pos, changed, haploid: cell.stage?.stage === 'haploid',
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
  return units.map((_, b) => sisterGainsFor(units, b, same))
}

/** The readings for cell `b` of the group alone, which asks only the pairs that include it. */
export function sisterGainsFor(
  units: readonly SisterUnit[], b: number,
  same: (i: number, j: number) => boolean = (i, j) => oneGenome(units[i], units[j]),
): SisterGain[] {
  const u = units[b]
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
}

/** A whole chromosome as an interval. */
export const WHOLE = (chrom: string): Interval => ({ chrom, startBp: 0, endBp: Number.MAX_SAFE_INTEGER, whole: true })

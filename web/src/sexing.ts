/**
 * Does this haploid product carry a Y chromosome?
 *
 * This is the only signal that names a reconstructed parent without having an array of that
 * parent. Grouping products by shared parent is reliable and says nothing about WHICH parent: a
 * group of siblings looks the same whether the shared parent is the mother or the father. A
 * maternal pronucleus cannot carry a Y, so one member with a whole Y names its group paternal.
 *
 * TWO signals are required and neither is sufficient. Measured over the 46 arrays of three
 * experiments:
 *
 *   - Call rate alone is wrong on 50402-11, which genotypes 86.2% of its Y probes while its Y
 *     intensity sits at -1.00 against its own autosomes, exactly where arrays with no Y sit.
 *     An absent chromosome still produces calls; that is noise on nothing, not a chromosome.
 *   - Intensity alone is wrong on 52461-16, which reads -0.10 while calling 0.0% of its Y
 *     probes. Nothing is there to genotype.
 *
 * Requiring both separates cleanly with room on either side. Y-bearing arrays call 93.7% to
 * 97.3% of Y probes at +0.16 to +0.43 log2; every other array either calls 0.0% or sits at
 * -0.81 to -1.25 log2. The nearest miss on each axis is the array named above.
 */
import type { ProbeRow } from './ingest.ts'
import { isAutosome } from './parentage.ts'

/** Y call rate, as a fraction of the array's own autosomal call rate. A whole Y genotypes about
 *  as well as an autosome does; the arrays that clear this reach 1.03 to 1.11. */
export const Y_CALL_MIN = 0.5

/** Y intensity, relative to the array's own autosomal median. A single-copy Y in a haploid cell
 *  sits at its autosomal level; an absent one sits about a log2 below. */
export const Y_LRR_MIN = -0.5

/**
 * Fewest Y probes that can say "no Y". A whole Y genotypes 93.7% to 97.3% of its probes, so a
 * handful of them all failing is not evidence of absence, it is a panel too thin to ask. The
 * bundled example subsets carry 51 and answer; a file carrying five does not.
 */
export const Y_MIN_PROBES = 20

/** Every hundredth autosomal marker, which is 8,000 or so values, fixes the autosomal median to
 *  far inside the 1.3 log2 gap this decides. Holding all 800,000 would cost megabytes to move
 *  the answer by nothing. */
const LRR_STRIDE = 100

export interface SexTally {
  yCalled: number
  yTotal: number
  autoCalled: number
  autoTotal: number
  yLrr: number[]
  autoLrr: number[]
  /** Counts markers seen, so the stride samples the same way regardless of file order. */
  seen: number
}

export const emptySex = (): SexTally => ({
  yCalled: 0, yTotal: 0, autoCalled: 0, autoTotal: 0, yLrr: [], autoLrr: [], seen: 0,
})

/** Fed one row at a time, in the pass that reads the file anyway. */
export function accumulateSex(r: ProbeRow, t: SexTally): void {
  if (r.chrom === 'Y') {
    t.yTotal += 1
    if (r.genotype !== 'NC') t.yCalled += 1
    if (r.log2R !== null) t.yLrr.push(r.log2R)
    return
  }
  if (!isAutosome(r.chrom)) return
  t.autoTotal += 1
  if (r.genotype !== 'NC') t.autoCalled += 1
  t.seen += 1
  if (t.seen % LRR_STRIDE === 0 && r.log2R !== null) t.autoLrr.push(r.log2R)
}

const median = (xs: number[]): number =>
  (xs.length ? [...xs].sort((a, b) => a - b)[xs.length >> 1] : NaN)

export interface SexCall {
  /** True only when both signals agree. Null when the file carries too few Y probes to ask,
   *  which is a panel that cannot answer rather than an array with no Y. */
  yBearing: boolean | null
  /** Y call rate over the array's own autosomal call rate. */
  callRatio: number
  /** Y median log2 ratio over the array's own autosomal median. */
  lrrShift: number
}

export function sexCall(t: SexTally): SexCall {
  if (t.yTotal < Y_MIN_PROBES || !t.autoTotal) {
    return { yBearing: null, callRatio: NaN, lrrShift: NaN }
  }
  const auto = t.autoCalled / t.autoTotal
  const callRatio = auto > 0 ? (t.yCalled / t.yTotal) / auto : NaN
  const lrrShift = median(t.yLrr) - median(t.autoLrr)
  // An array with Y probes but no intensity column can still be judged on calls alone, which is
  // the weaker of the two and is said to be so rather than silently treated as equal.
  if (!Number.isFinite(lrrShift)) {
    return { yBearing: callRatio >= Y_CALL_MIN, callRatio, lrrShift }
  }
  return {
    yBearing: callRatio >= Y_CALL_MIN && lrrShift >= Y_LRR_MIN,
    callRatio,
    lrrShift,
  }
}

/**
 * Which group is the father's.
 *
 * Returns the index of the one group carrying a Y, or null when that cannot be decided. Null is
 * returned rather than a guess in both directions that matter:
 *
 *   - No group carries a Y. A paternal group of n products is all X-bearing with probability
 *     2^-n, which at n=5 is 3%. Silence here is a real possibility, not an error.
 *   - More than one group carries a Y. One sperm donor cannot produce two unrelated groups of
 *     products, so the input is not what it was said to be and naming either would be a guess.
 */
export function paternalGroup(
  groups: readonly (readonly number[])[], yBearing: readonly (boolean | null)[],
): number | null {
  const withY = groups
    .map((g, i) => (g.some((m) => yBearing[m] === true) ? i : -1))
    .filter((i) => i >= 0)
  return withY.length === 1 ? withY[0] : null
}

/**
 * Does this array carry a Y chromosome?
 *
 * This is the only signal that names a reconstructed parent without having an array of that
 * parent. Grouping products by shared parent is reliable and says nothing about WHICH parent: a
 * group of siblings looks the same whether the shared parent is the mother or the father. A
 * maternal pronucleus cannot carry a Y, so one member with a whole Y names its group paternal.
 *
 * ON THIS CHIP THE CALLS ARE NOT A SECOND SIGNAL. The UK Biobank Axiom chip leaves chrY uncalled
 * on any sample it reads as female from its X heterozygosity, so "chrY has calls" restates the X
 * and says nothing of a Y that sits beside two X. A Y that is there is hemizygous: its markers
 * read one allele each, as the array's own homozygous calls do. A Y that is not there is noise.
 * So the reading rests on the Y markers' allele fractions, with intensity as the check, and the
 * calls may only block a reading where they and the intensity both say the opposite:
 *
 *   - allele fractions are read only where this array's own homozygous autosomal calls sit at a
 *     pole at least Y_BAF_READABLE of the time; under that the Y share reads noise whatever the Y;
 *   - there, a Y share of Y_BAF_PRESENT or more at an intensity of Y_LRR_MIN or more is a Y, and a
 *     share of Y_BAF_ABSENT or less is none, unless calls and intensity both say present;
 *   - everywhere else the calls and intensity are read together, and only where they agree.
 *
 * Neither calls nor intensity alone suffices, which is why they are read together where allele
 * fractions do not decide. Call rate alone is wrong on an array that genotypes 86.2% of its Y
 * probes while its Y intensity sits at -1.00 against its own autosomes, where arrays with no Y
 * sit: an absent chromosome still produces calls on a cell whose X reads as one copy. Intensity
 * alone is wrong where its scale moves with the export and the ploidy: a woman's bulk array on
 * the normalized-intensity export reads chrY at -0.16 to -0.02 against her autosomes, and a
 * haploid maternal cell at +0.04 to +0.26.
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

/** Share of this array's own homozygous autosomal calls at a pole for its allele fractions to read
 *  a Y. Under it the arrays are noisy enough that a chrY share reads noise whatever the Y. */
export const Y_BAF_READABLE = 0.90

/** Share of chrY markers at a pole that is a Y, where allele fractions are readable: a present Y
 *  reads 0.75 or more on 80 of 85 arrays with a known Y. */
export const Y_BAF_PRESENT = 0.75

/** Share of chrY markers at a pole that is no Y: arrays with no Y read 0.60 or less on 118 of 122,
 *  and every known no-Y parent and maternal haploid cell 0.38 to 0.57. */
export const Y_BAF_ABSENT = 0.60

/** A B allele fraction under this, or over one less it, sits at a pole: one allele. */
const POLE = 0.15

/** Every tenth autosomal marker measures how polarised this array's homozygous calls are. */
const FLOOR_STRIDE = 10

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
  /** chrY markers carrying an allele fraction, and those of them at a pole. */
  yBaf: number
  yBafPole: number
  /** Sampled homozygous autosomal calls carrying an allele fraction, and those of them at a pole. */
  homBaf: number
  homBafPole: number
}

export const emptySex = (): SexTally => ({
  yCalled: 0, yTotal: 0, autoCalled: 0, autoTotal: 0, yLrr: [], autoLrr: [], seen: 0,
  yBaf: 0, yBafPole: 0, homBaf: 0, homBafPole: 0,
})

const atPole = (b: number): boolean => b < POLE || b > 1 - POLE

/** Fed one row at a time, in the pass that reads the file anyway. */
export function accumulateSex(r: ProbeRow, t: SexTally): void {
  if (r.chrom === 'Y') {
    t.yTotal += 1
    if (r.genotype !== 'NC') t.yCalled += 1
    if (r.log2R !== null) t.yLrr.push(r.log2R)
    if (r.baf !== null && Number.isFinite(r.baf)) { t.yBaf += 1; if (atPole(r.baf)) t.yBafPole += 1 }
    return
  }
  if (!isAutosome(r.chrom)) return
  t.autoTotal += 1
  if (r.genotype !== 'NC') t.autoCalled += 1
  t.seen += 1
  if (t.seen % LRR_STRIDE === 0 && r.log2R !== null) t.autoLrr.push(r.log2R)
  if (t.seen % FLOOR_STRIDE === 0 && r.baf !== null && Number.isFinite(r.baf)
    && (r.genotype === 'AA' || r.genotype === 'BB')) {
    t.homBaf += 1
    if (atPole(r.baf)) t.homBafPole += 1
  }
}

const median = (xs: number[]): number =>
  (xs.length ? [...xs].sort((a, b) => a - b)[xs.length >> 1] : NaN)

export interface SexCall {
  /** Null when the file carries too few Y probes to ask, which is a panel that cannot answer
   *  rather than an array with no Y, or where the Y markers do not settle it (`why` says how). */
  yBearing: boolean | null
  /** Y call rate over the array's own autosomal call rate. */
  callRatio: number
  /** Y median log2 ratio over the array's own autosomal median. */
  lrrShift: number
  /** Share of chrY markers whose allele fraction sits at a pole. */
  yBafShare?: number
  /** The same share over this array's own homozygous autosomal calls. */
  bafFloor?: number
  /** True where the Y was read from its allele fractions. */
  byAlleles?: boolean
  /** Why, where the reading is not the plain one: a Y the chip left uncalled, or a Y not read. */
  why?: string
}

const pct = (x: number): string => `${Math.round(x * 100)}%`

/** The allele-fraction clause of a Y read from its markers. */
export const yAlleles = (c: SexCall): string => `its markers read one allele each `
  + `(${pct(c.yBafShare!)} of them, where this array's own homozygous calls read ${pct(c.bafFloor!)}) `
  + `at male intensity${c.callRatio >= Y_CALL_MIN ? '' : ', though the chip left them uncalled'}`

export function sexCall(t: SexTally): SexCall {
  if (t.yTotal < Y_MIN_PROBES || !t.autoTotal) {
    return { yBearing: null, callRatio: NaN, lrrShift: NaN }
  }
  const auto = t.autoCalled / t.autoTotal
  const callRatio = auto > 0 ? (t.yCalled / t.yTotal) / auto : NaN
  const lrrShift = median(t.yLrr) - median(t.autoLrr)
  const yBafShare = t.yBaf ? t.yBafPole / t.yBaf : NaN
  const bafFloor = t.homBaf ? t.homBafPole / t.homBaf : NaN
  const base = { callRatio, lrrShift, yBafShare, bafFloor }
  const calls = callRatio >= Y_CALL_MIN
  // No intensity column: the calls have nothing to agree with.
  const inten = Number.isFinite(lrrShift) ? lrrShift >= Y_LRR_MIN : null
  const readable = bafFloor >= Y_BAF_READABLE
  if (readable && yBafShare >= Y_BAF_PRESENT && inten === true) {
    const c = { ...base, yBearing: true, byAlleles: true }
    return calls ? c : { ...c, why: `chrY is present: ${yAlleles(c)}` }
  }
  const shares = `${pct(yBafShare)} of its markers read one allele each, where this array's own `
    + `homozygous calls read ${pct(bafFloor)}`
  if (readable && yBafShare <= Y_BAF_ABSENT) {
    if (!(calls && inten === true)) return { ...base, yBearing: false }
    return { ...base, yBearing: null, why: `chrY was not read: only ${shares}, at or under the `
      + `${pct(Y_BAF_ABSENT)} that is no Y, while its calls and intensity both say present` }
  }
  if (calls && inten === true) return { ...base, yBearing: true }
  if (!calls && inten === false) return { ...base, yBearing: false }
  const alleles = !Number.isFinite(bafFloor)
    ? 'this file carries no allele fractions at its homozygous calls'
    : !readable ? `this array's own homozygous calls read only ${pct(bafFloor)} at one allele, under `
      + `the ${pct(Y_BAF_READABLE)} at which allele fractions can show a Y`
      : !Number.isFinite(yBafShare) ? 'its markers carry no allele fractions'
        : yBafShare < Y_BAF_PRESENT ? `${shares}, between the ${pct(Y_BAF_ABSENT)} that is no Y and the `
          + `${pct(Y_BAF_PRESENT)} that is one`
          : `${shares}, at an intensity under the male level`
  const agree = inten === null ? 'this file carries no chrY intensity to check its calls against'
    : `its calls (${callRatio.toFixed(2)} of the autosomes' rate, against ${Y_CALL_MIN.toFixed(2)}) say `
      + `${calls ? 'present' : 'absent'} while its intensity (${lrrShift.toFixed(2)} log2 from the `
      + `autosomes, against ${Y_LRR_MIN.toFixed(2)}) says ${inten ? 'present' : 'absent'}`
  return { ...base, yBearing: null, why: `chrY was not read: ${alleles}, and ${agree}` }
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

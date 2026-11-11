/**
 * Answer keys: markers whose right answer is known before the sample is read, measured per cell.
 *
 * WHY. Every quality gate this tool had compared a number to a line chosen once, on other material:
 * a heterozygosity ceiling derived from bulk drop-in, a heterozygous fraction under which a chromosome
 * counts as one copy. On a cell whose error rates differ from the material the line was drawn on, the
 * line refuses real calls or admits false ones, and nothing in the output says which. An answer key
 * measures the error on THIS cell, at markers where the parents have already said what the right
 * answer is, so a claim can be judged against what this array actually does.
 *
 * SAME-HOMOZYGOTE MARKERS. Where both measured parents carry the same homozygote, the sample can only
 * be that homozygote: no parent, no relative and no recombination can put another allele there. So a
 * heterozygous call is DROP-IN, and the opposite homozygote is an outright ERROR, the rate at which an
 * allele absent from the template appears as a homozygous call. That second rate is what a lost
 * parent's allele should read on a chromosome that has lost that parent's copy.
 *
 * Measured by the over-refusal consult on GSE148488: drop-in 0.220 to 0.307 on the seven cells refused
 * by the heterozygosity ceiling against a 0.020 median on admitted cells, wrong homozygote 0.007 to
 * 0.039 on the same seven. Both need two MEASURED bulk parents; a reconstructed parent is homozygous at
 * markers where the real one is not, and those markers would read as errors that are not.
 */
import type { AB } from './informativity.ts'

export interface SameHomTally {
  /** Called sample genotypes at markers where both parents are the same homozygote. */
  n: number
  /** Of those, heterozygous: drop-in. */
  het: number
  /** Of those, the opposite homozygote: an allele absent from the template, called homozygous. */
  wrong: number
}

export const emptySameHom = (): SameHomTally => ({ n: 0, het: 0, wrong: 0 })

export const addSameHom = (t: SameHomTally, sample: AB | null, fa: AB | null, mo: AB | null): void => {
  if (fa !== mo || (fa !== 'AA' && fa !== 'BB')) return
  if (sample !== 'AA' && sample !== 'AB' && sample !== 'BB') return
  t.n += 1
  if (sample === 'AB') t.het += 1
  else if (sample !== fa) t.wrong += 1
}

/** Heterozygous calls where there can be none. */
export const dropIn = (t: SameHomTally): number => (t.n ? t.het / t.n : NaN)

/** The opposite homozygote where there can be none: this cell's floor for an absent allele. */
export const absentAlleleFloor = (t: SameHomTally): number => (t.n ? t.wrong / t.n : NaN)

/** Same-homozygote markers needed before either rate is a measurement. A cell carries hundreds of
 *  thousands; under this the parents are barely overlapping and the rates are not reported. */
export const MIN_SAME_HOM = 10_000

/**
 * The evidence standard for a single chromosome's allele state: 1000 to 1. Same order as SEGMENT_LRT,
 * and stated as a likelihood ratio rather than a fraction so that it means the same thing on a cell
 * with 300 informative markers and one with 3,000.
 */
export const ALLELE_STATE_EVIDENCE = Math.log(1000)

/** Binomial log-likelihood, without the constant, with probabilities held off 0 and 1. */
const ll = (x: number, n: number, p: number): number => {
  const q = Math.min(1 - 1e-12, Math.max(1e-12, p))
  return x * Math.log(q) + (n - x) * Math.log(1 - q)
}

/** How strongly the observed rate x/n rejects a hypothesised rate p: the log likelihood ratio of the
 *  best fit against p. Zero when x/n equals p. */
const rejects = (x: number, n: number, p: number): number => ll(x, n, x / n) - ll(x, n, p)

export type AlleleState = 'lost' | 'partial' | 'none' | 'undecided'

export interface AlleleStateCall {
  state: AlleleState
  /** The two answer keys the call was made against, carried so a reader sees the comparison. */
  floor: number
  twoCopyMajority: number
  /** Log likelihood ratios: how strongly the dominant parent's share exceeds what a two-copy
   *  chromosome reads on this cell, and how strongly the rarer parent's exceeds the error floor. */
  majorityAboveTwoCopy: number
  minorityAboveFloor: number
}

/**
 * Whether one parent's copy is gone from a chromosome, at markers where the two parents are opposite
 * homozygotes, judged against this cell's own answer keys.
 *
 * Two readings, one per parent: the share of calls carrying the DOMINANT parent's allele alone, and
 * the share carrying the RARER parent's alone.
 *
 *   lost      the dominant share is clearly above what this cell's two-copy chromosomes read, and the
 *             rarer share is at this cell's error floor: the rarer parent's copy is gone.
 *   partial   the dominant share is clearly above the two-copy level and the rarer share clearly above
 *             the floor: the copy is gone from part of the sample, or the sample mixes states. No
 *             parent is named.
 *   none      the dominant share is not clearly above the two-copy level: no loss on this channel.
 *
 * WHY THE DOMINANT SIDE DECIDES WHETHER ANYTHING HAPPENED. On a two-copy chromosome the one-parent
 * calls come from dropout and fall on both sides alike; a lost copy pushes one side to most of the
 * calls. On error-free material two-copy chromosomes carry no one-parent calls at all, so a test made
 * on the rarer side could not tell a loss from nothing there. And an EXTRA copy barely moves the
 * dominant side, because dropping the single minority allele is about as likely as on two copies, so
 * a trisomy reads `none` here rather than as a partial loss; gains are gainOrigin's question.
 *
 * NO FRACTION THRESHOLD. The lines are this cell's floor and its own two-copy level, and the decision
 * is a likelihood ratio at ALLELE_STATE_EVIDENCE. On a cell whose calls are noise the two-copy level
 * rises toward the dominant share of everything, and nothing separates.
 */
export function alleleState(
  minority: number, majority: number, n: number, floor: number, twoCopyMajority: number,
): AlleleStateCall {
  const base = { floor, twoCopyMajority, majorityAboveTwoCopy: NaN, minorityAboveFloor: NaN }
  if (!(n > 0) || !Number.isFinite(floor) || !Number.isFinite(twoCopyMajority)) {
    return { ...base, state: 'undecided' }
  }
  const majorityAboveTwoCopy = majority / n <= twoCopyMajority ? 0 : rejects(majority, n, twoCopyMajority)
  const minorityAboveFloor = minority / n <= floor ? 0 : rejects(minority, n, floor)
  const call = { floor, twoCopyMajority, majorityAboveTwoCopy, minorityAboveFloor }
  const E = ALLELE_STATE_EVIDENCE
  if (majorityAboveTwoCopy < E) return { ...call, state: 'none' }
  return { ...call, state: minorityAboveFloor < E ? 'lost' : 'partial' }
}

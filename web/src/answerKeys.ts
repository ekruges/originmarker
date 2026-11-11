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

/** Same-homozygote markers needed before drop-in is a measurement. A cell carries hundreds of
 *  thousands; under this the parents barely overlap and the rate is not used. */
export const MIN_SAME_HOM = 10_000

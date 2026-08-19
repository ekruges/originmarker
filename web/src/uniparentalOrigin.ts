/**
 * Parental origin from ZYGOSITY, for the events the dosage channel cannot reach.
 *
 * THE CASE THIS EXISTS FOR. A real run: 18 samples, 8 gynogenetic and 7 androgenetic, 1,272
 * detected regions, and not one parental origin called. Every row read "not evaluable", because a
 * sub-chromosomal interval on single-cell material has no detection floor at any mosaic fraction,
 * with one genotyped parent or two. That is true, and it is the wrong question to have asked.
 *
 * The dosage channel asks WHICH OF TWO PARENTS a shift points at. That question only arises when
 * both parents' genomes are present. In a uniparental sample they are not: a gynogenetic embryo
 * carries maternal DNA and no paternal complement, so a chromosome lost from it was maternal,
 * because there was no other chromosome there to lose. An androgenetic one is the mirror. The
 * answer is Mendelian and needs no detection floor, exactly as the obligate-het channel needs none.
 *
 * This is the same error that 5.10.0 removed from triploidy and complex genomes: a limit that is
 * real for ONE channel, applied to every channel, so a question already answered elsewhere was
 * reported unanswerable.
 *
 * WHAT IT IS NOT. It is not a per-event measurement. Every event in one sample gets the SAME
 * answer and the same confidence, because they all rest on the single genome-level call. The
 * report must say so rather than let a reader count 1,272 independent confirmations of one fact.
 */
import { SYSTEMATIC_ERROR_BOUND, reportedConfidence } from './oneParentOrigin.ts'

export type UniparentalInput = {
  /** The genome-level call. Only the two uniparental classes carry an answer. */
  originClass: 'androgenetic' | 'gynogenetic' | 'biparental' | 'unclear' | string
  /**
   * Whether the sample carries one parental contribution or two.
   *
   * The shipped values are 'diploid', 'uniparental_homozygous' and 'unknown'. Matched on the
   * prefix rather than the whole string, because the first version of this compared against
   * 'uniparental', which is not one of them, and the channel silently answered nothing on every
   * real sample while its own tests passed on a value the application never produces.
   */
  zygosity: string
  /** Which parent the LOADED array is. A run may load either, and assuming names the wrong one. */
  role: 'paternal' | 'maternal'
  /** Absence measured across the genome, and the level that dropout and error alone explain. */
  genomeRate: number
  explainable: number
}

export type UniparentalCall = {
  verdict: 'loaded-parent' | 'other-parent'
  /** Which parent in absolute terms, so a reader never has to resolve `loaded` themselves. */
  parent: 'paternal' | 'maternal'
  confidence: number
  band: 'A' | 'B' | 'C' | 'D'
  why: string
  /** How far above the explainable ceiling the genome-level call sits. */
  foldOverCeiling: number
}

/**
 * How fast confidence rises with the genome-level margin.
 *
 * CHOSEN, NOT FITTED, and it is the one number here that is. There is no injection series for this
 * channel, so nothing would justify claiming a calibration. What bounds it instead is the same
 * systematic error bound every other channel is bounded by, from 13 validation units: the cap is
 * 1 - 0.206, so a call from this channel cannot reach band B however decisive the zygosity is.
 * That ceiling is the honest statement. A margin of 2x is the least that gets called at all and a
 * margin of 13x is the most seen in practice, so a half-evidence of 1.0 on the log spreads that
 * range across the reportable interval rather than pinning it at either end.
 */
export const MARGIN_HALF_EVIDENCE = 1.0

/** Bands, shared with every other channel so a reader compares rows rather than scales. */
const bandOf = (p: number): 'A' | 'B' | 'C' | 'D' =>
  (p >= 0.985 ? 'A' : p >= 0.90 ? 'B' : p >= 0.75 ? 'C' : 'D')

/**
 * The parent every event in this sample belongs to, or null where the question is still open.
 *
 * Returns null for a biparental or unclear sample, which is the case the dosage channel is for,
 * and for a genome-level call that did not clear its own ceiling.
 */
export function uniparentalOrigin(input: UniparentalInput): UniparentalCall | null {
  const genomeIs = input.originClass === 'gynogenetic' ? 'maternal'
    : input.originClass === 'androgenetic' ? 'paternal' : null
  // A biparental or unclear genome has two candidate parents, which is the dosage channel's job.
  if (!genomeIs) return null
  // The class is read off the absence measurement, so a sample whose zygosity was not resolved
  // has no basis for it however the class was labelled.
  if (!input.zygosity?.startsWith('uniparental')) return null
  const { genomeRate, explainable } = input
  if (!Number.isFinite(genomeRate) || !Number.isFinite(explainable) || explainable <= 0) return null
  const fold = genomeRate / explainable
  // Below its own ceiling the genome-level call is not made, and inheriting from it would be
  // inventing evidence.
  if (!(fold > 1)) return null

  const confidence = reportedConfidence(Math.log(fold), {
    bound: SYSTEMATIC_ERROR_BOUND, halfEvidence: MARGIN_HALF_EVIDENCE,
  })
  const verdict = genomeIs === input.role ? 'loaded-parent' : 'other-parent'
  return {
    verdict,
    parent: genomeIs,
    confidence,
    band: bandOf(confidence),
    foldOverCeiling: fold,
    why: `this sample carries ONE parental genome, called ${input.originClass} at `
      + `${fold.toFixed(1)}x the level dropout and error alone explain, so every change in it is `
      + `${genomeIs}: there was no ${genomeIs === 'maternal' ? 'paternal' : 'maternal'} copy `
      + 'present to lose, gain or rearrange. This is inherited from the genome-level call rather '
      + 'than measured on this interval, so it is the same answer and the same confidence for '
      + 'every change in this sample, and it is capped below band B because the channel has no '
      + 'injection series of its own',
  }
}

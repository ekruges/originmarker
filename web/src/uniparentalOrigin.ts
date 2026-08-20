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
import { HET_BAND_DIPLOID } from './parentage.ts'

/**
 * What the inherited grade is worth, measured against DISSECTION rather than inference.
 *
 * Fifteen single pronuclei were dissected from 2PN zygotes, so which parent each carries is known
 * physically and not from any statistic. Running the genome-level call on them:
 *
 *   called and correct   12
 *   called and wrong      0
 *   declined as unclear   3
 *
 * Zero errors. That is the accuracy every row carrying this grade rests on, because the grade IS
 * that call inherited.
 *
 * CLUSTERED ON THE ZYGOTE, which is the sampling unit: the two pronuclei of one zygote come from
 * one fertilisation, and every paternal pronucleus in the series comes from a single sperm donor.
 * Eight zygotes with zero events give a one-sided 95% error bound of 0.3123, so the reportable
 * floor is 0.688. Treating the 12 calls as independent would give 0.779 and would be the same
 * cluster error a methods review already measured at 3.5 to 5.8x too narrow on this material.
 *
 * WHAT IT CANNOT SUPPORT. A dissected pronucleus has been through neither syngamy nor a mitosis, so
 * mosaicism, post-zygotic loss, lineage segregation and chimaerism are structurally absent from
 * this series. Those are exactly the mechanisms that would break the inheritance step in real
 * embryo material, so this is an upper bound on embryo performance rather than an estimate of it.
 */
export const INHERITED_VALIDATION = {
  correct: 12,
  wrong: 0,
  declined: 3,
  /** Zygotes, not arrays. See above. */
  clusters: 8,
  /** One-sided 95% floor from zero events over the clusters: 1 - 0.05^(1/8). */
  accuracyFloor: 0.688,
} as const

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
  /**
   * Fraction of calls in the heterozygous band, which is what says ONE parental contribution or
   * two independently of whose it is.
   *
   * Load-bearing for the case where the genome carries the LOADED parent. Absence cannot measure
   * that case at all: a parent that is present is not absent, so the absence route is silent by
   * construction and this is the only evidence left.
   */
  hetBand: number
}

export type UniparentalCall = {
  verdict: 'loaded-parent' | 'other-parent'
  /** Which parent in absolute terms, so a reader never has to resolve `loaded` themselves. */
  parent: 'paternal' | 'maternal'
  /**
   * NO CONFIDENCE NUMBER, DELIBERATELY.
   *
   * The number this channel used to emit was `reportedConfidence(log(margin))`, which is monotone
   * in the margin and in nothing else: a strictly increasing reparameterisation of a quantity
   * already reported beside it. A reader holding the margin learned nothing from it, and a reader
   * holding only the number had lost the interpretable quantity.
   *
   * Its range made it worse than uninformative. Every value the channel could emit was band D
   * until the margin passed roughly 13,000x, and the observed range of 4.8x to 17.4x compressed
   * into 0.65 to 0.68: three decimal places moving by 0.03 across every input the channel will
   * ever see, carrying the same band the dosage channel uses for its weakest measured guesses,
   * where measured accuracy is 0.60 to 0.64. A gynogenetic inference is not 60% accurate.
   * Conditional on the genome-level call it is deductive, and the band said the opposite.
   *
   * The bound that produced the cap is a rule of three on 13 units, which bounds the GENOME-LEVEL
   * call. It is not a bound on the inheritance step, which has no error of its own: given the
   * class, the parent follows. The cap was being applied to the wrong link in the chain.
   *
   * What is emitted instead is the verdict, the margin, and the basis, which cannot be quoted as a
   * per-event accuracy because it is not one.
   */
  band: 'inherited'
  why: string
  /** How far above the explainable ceiling the genome-level call sits. */
  foldOverCeiling: number
}


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
  // WHICH ROUTE ESTABLISHED THE CLASS DECIDES WHICH MEASUREMENT SUPPORTS IT.
  //
  // The first version tested absence of the loaded parent for every sample, which can only ever
  // be large when the loaded parent is MISSING. Measured on this project's own audit set, that
  // gate fired on 14 of 14 gynogenetic samples and 0 of 8 androgenetic ones, because an
  // androgenetic genome loaded against the sperm donor contains the loaded parent and its absence
  // sits at 0.09 to 0.62 of the explainable ceiling. Half the samples were dropped in silence, and
  // any aggregation fed from these rows would have found a maternal excess whatever the biology.
  //
  // A genome that CARRIES the loaded parent is established by zygosity instead: one contribution
  // rather than two, from the heterozygous band. That correctly refuses the three androgenetic
  // arrays sitting at hetBand 0.15 to 0.16, above the diploid threshold, which the absence route
  // had been admitting on evidence that says the opposite.
  const matchesLoaded = genomeIs === input.role
  const { hetBand } = input
  if (matchesLoaded && (!Number.isFinite(hetBand) || hetBand <= 0)) return null
  const margin = matchesLoaded ? HET_BAND_DIPLOID / hetBand : genomeRate / explainable
  const basis = matchesLoaded
    ? `${(hetBand * 100).toFixed(2)}% of calls in the heterozygous band against the `
      + `${(HET_BAND_DIPLOID * 100).toFixed(0)}% a second parental contribution would produce, `
      + `${margin.toFixed(1)}x below it`
    : `absence of the loaded parent at ${margin.toFixed(1)}x the level dropout and error alone `
      + 'explain'
  // Below its own threshold the genome-level call is not supported by the route that has to carry
  // it, and inheriting from it would be inventing the evidence.
  if (!(margin > 1)) return null
  const fold = margin
  const verdict = matchesLoaded ? 'loaded-parent' : 'other-parent'
  return {
    verdict,
    parent: genomeIs,
    band: 'inherited',
    foldOverCeiling: fold,
    why: `This sample carries ONE parental genome, called ${input.originClass} on ${basis}, so `
      + `every change in it is ${genomeIs}: there was no `
      + `${genomeIs === 'maternal' ? 'paternal' : 'maternal'} copy present to lose, gain or `
      + 'rearrange. This is inherited from the genome-level call rather than measured on this '
      + 'interval, so it is the same answer and the same confidence for every change in this '
      + 'sample, and it is capped below band B because the channel has no injection series of '
      + 'its own',
  }
}

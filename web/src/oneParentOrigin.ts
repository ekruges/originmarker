/**
 * Whose copy is missing, from ONE genotyped parent.
 *
 * THE MENDELIAN FACT THIS RESTS ON. At a marker where the known parent is homozygous AA, a child
 * reading BB cannot be carrying that parent's copy: the parent has no B to give. So a single BB
 * observation is direct evidence the KNOWN parent's copy is gone, and it is evidence no amount of
 * dropout can manufacture, because dropout removes an allele and never invents one. The estimator
 * below is that observation counted properly, with the two ways it can be faked priced in.
 *
 * The three hypotheses at an informative marker, with q the frequency of the allele the known
 * parent does NOT carry:
 *
 *                              child reads AA        child reads AB      child reads BB
 *   both copies present        1 - q                 q                   0
 *   KNOWN parent's copy lost   1 - q                 0                   q
 *   OTHER parent's copy lost   1                     0                   0
 *
 * Read the columns rather than the rows. AB says both copies are present, because the known parent
 * contributes A and something contributed B. BB says the known parent's copy is absent. AA alone
 * says little, since every hypothesis allows it, which is why a region of pure AA is the hard case
 * and is separated by the RATE at which BB fails to appear rather than by any single marker.
 *
 * WHY THIS IS NOT THE STATISTIC THAT FAILED BEFORE. An earlier attempt formed an allele share and
 * centred it on the sample's own median. For every q below 0.5 the majority of these markers sit
 * at the top of the range, so the median was the ceiling and the upward direction had no room:
 * "the other parent's copy is lost" could not be detected at all, and a hundred percent of calls
 * came back the same way. Here there is no share and no centre. Each marker contributes a
 * likelihood and the region's evidence is their product, so neither direction is privileged and
 * the asymmetry that produced those results cannot arise.
 *
 * WHAT IT STILL CANNOT DO. It names WHICH PARENT only in the sense of known against unknown. If
 * the loaded parent is the father, "known parent's copy lost" means paternal; the caller supplies
 * that mapping. It cannot say anything at all about a marker where the known parent is
 * heterozygous, and those are simply excluded.
 */
import type { AB } from './informativity.ts'
import { bandObligateHet, type Band } from './originPosterior.ts'

/** Frequency of the allele the known parent lacks, when it cannot be estimated per marker. */
export const DEFAULT_Q = 0.30

/** Genotyping error: a homozygote read as the other homozygote, or a spurious heterozygote. */
export const GENOTYPE_ERROR = 0.02

/**
 * Drop-in, a heterozygous call at a truly homozygous marker. Measured on this platform across 113
 * same-genome array pairs, medians 0.0390 to 0.0525 by experiment. It matters here because a
 * spurious AB argues for both copies being present and so masks a real deletion.
 */
export const DROP_IN = 0.0435

/** Posterior a hypothesis must reach before the region is called. */
export const CALL_POSTERIOR = 0.95

/** Fewest informative markers before a region is scored. */
export const MIN_MARKERS = 50

/**
 * Heterozygosity above which the region's genotypes are not measuring it, and no verdict is issued.
 *
 * At a marker where the known parent is HOMOZYGOUS, a biparental child is heterozygous exactly
 * when the other parent transmitted the allele this one lacks. That is the panel's own allele
 * frequency, so with `q` at its default 0.30 and drop-in at the top of the range measured here,
 * the most a real region can show is 0.30 + 0.0525, and dropout only lowers it. The ceiling sits
 * above that with margin.
 *
 * WHY IT IS NEEDED, from a real array. A blastomere carrying a call-rate collapse on chromosome 1,
 * 69% no-call, returned 59.7% heterozygous across the markers where the father is homozygous. No
 * genome reads that way. The likelihood took the heterozygote count at face value, since a
 * heterozygote is near-impossible under either deletion hypothesis, and returned both-copies-
 * present at posterior 1.000 on genotypes that carried no information at all. A confident verdict
 * from meaningless input is worse than a refusal, because only the refusal is visible.
 */
export const MAX_REGION_HET = 0.40

/**
 * Allele dropout for a sample, inferred from the sample itself rather than declared.
 *
 * The drag-and-drop contract means nothing may be asked of the user that the data can answer, and
 * developmental stage is one of those things: a trophectoderm biopsy and a blastomere differ by
 * six-fold in dropout (0.050 against 0.308) and the difference is visible in the arrays.
 *
 * Heterozygosity is the readout. Dropout removes one allele of a heterozygote, so a stage that
 * drops heavily reads fewer heterozygotes than its genome contains. Against the bulk figure for
 * this platform, the shortfall IS the dropout rate. Measured medians by stage: bulk 0.013,
 * trophectoderm 0.050, single ESC 0.199, blastomere 0.308.
 *
 * Bounded at both ends: a sample cannot have negative dropout, and above the blastomere rate the
 * estimate stops being informative and the caller should be refusing on marker count anyway.
 */
export const BULK_HETEROZYGOSITY = 0.168

export function inferDropout(sampleHeterozygosity: number): number {
  if (!Number.isFinite(sampleHeterozygosity) || sampleHeterozygosity <= 0) return 0.308
  const shortfall = 1 - sampleHeterozygosity / BULK_HETEROZYGOSITY
  return Math.min(0.6, Math.max(0.01, shortfall))
}

export type OneParentVerdict = 'known-parent-lost' | 'other-parent-lost' | 'both-present' | 'refused'

export interface OneParentCall {
  verdict: OneParentVerdict
  posterior: number
  /**
   * The band this posterior lands in, CAPPED AT B because the channel can never reach A.
   *
   * Structural rather than a power argument. This channel asks whether a forbidden allele appeared
   * where the homozygous parent had none to give, so its failure mode is dropout, and dropout is
   * the same event as the observation: "she did not transmit it" and "amplification lost it" look
   * identical. No quantity of markers separates them, so no amount of evidence earns the top band.
   */
  band: Band
  /** Markers where the known parent is homozygous and the sample is called. */
  markers: number
  /** Of those, how many carry the allele the known parent does not have. Mendelian evidence. */
  exclusive: number
  /** Of those, how many are heterozygous, which argues both copies are present. */
  heterozygous: number
  q: number
  why: string
}

/**
 * A marker is informative only where the known parent is HOMOZYGOUS.
 *
 * Where that parent is heterozygous every hypothesis predicts the same distribution of child
 * genotypes, so the marker carries nothing and including it would dilute the evidence with noise.
 */
export const informative = (parent: AB): boolean => parent === 'AA' || parent === 'BB'

/**
 * Per-marker likelihood under the three hypotheses.
 *
 * `carriesOther` is true when the sample shows the allele the known parent lacks; `het` when it is
 * heterozygous. Orientation is handled by the caller so that this function is the same whichever
 * homozygote the parent is.
 */
function likelihood(
  het: boolean, carriesOther: boolean, q: number, ado: number, eps: number, dropIn: number,
): [number, number, number] {
  if (het) {
    // Both copies present predicts heterozygosity at rate q, surviving dropout. A deletion cannot
    // produce it except by drop-in, which is why drop-in bounds this estimator.
    return [q * (1 - ado) + dropIn, dropIn, dropIn]
  }
  if (carriesOther) {
    // The known parent has no such allele, so this is impossible unless its copy is absent.
    // Under both-copies-present it requires the heterozygote to have dropped the parent's allele.
    return [q * (ado / 2), q, eps / 2]
  }
  // Homozygous for the known parent's own allele: every hypothesis allows it.
  return [(1 - q) + q * (ado / 2), 1 - q, 1 - eps]
}

/**
 * Call which copy is missing across a region, from one genotyped parent.
 *
 * `pairs` are (parent genotype, sample genotype) at each marker in the region. `q` is the frequency
 * of the allele the parent lacks; passing a per-marker value is better than the default, and
 * getting it wrong biases the two loss hypotheses against each other, so it is reported.
 */
export function callOneParentOrigin(
  pairs: readonly (readonly [AB, AB])[],
  ado: number,
  q = DEFAULT_Q,
  eps = GENOTYPE_ERROR,
  dropIn = DROP_IN,
  /**
   * Thresholds, so a caller that knows what it is doing can move them.
   *
   * The web tool never passes this: it offers no knob whose correct setting the user would have
   * to know. The command line does, because the person there is usually asking what happens at a
   * different threshold, and the honest way to answer is to let them move it and say that they
   * did. Defaults are the measured values, so omitting this is the shipped configuration.
   */
  opts: {
    minMarkers?: number, maxRegionHet?: number, callPosterior?: number,
  } = {},
): OneParentCall {
  const minMarkers = opts.minMarkers ?? MIN_MARKERS
  const maxRegionHet = opts.maxRegionHet ?? MAX_REGION_HET
  const callPosterior = opts.callPosterior ?? CALL_POSTERIOR
  let n = 0
  let exclusive = 0
  let het = 0
  const logs: [number, number, number] = [0, 0, 0]
  for (const [p, s] of pairs) {
    if (!informative(p) || s === 'NC') continue
    n += 1
    const isHet = s === 'AB'
    // The allele the parent lacks: B when the parent is AA, A when BB.
    const carriesOther = !isHet && s !== p
    if (isHet) het += 1
    if (carriesOther) exclusive += 1
    const l = likelihood(isHet, carriesOther, q, ado, eps, dropIn)
    logs[0] += Math.log(l[0])
    logs[1] += Math.log(l[1])
    logs[2] += Math.log(l[2])
  }

  const base = { markers: n, exclusive, heterozygous: het, q }
  if (n < minMarkers) {
    return {
      ...base, verdict: 'refused', posterior: NaN, band: 'D' as Band,
      why: `${n} informative markers is under the ${minMarkers} this needs; a marker only counts `
        + 'where the loaded parent is homozygous',
    }
  }

  // The genotypes must be measuring the region before their likelihood means anything. This is
  // the region-level form of the array-level ceiling in stage.ts, and it exists for the same
  // reason: an impossible rate is evidence about the reaction, not about the genome.
  if (het / n > maxRegionHet) {
    return {
      ...base, verdict: 'refused', posterior: NaN, band: 'D' as Band,
      why: `${(100 * het / n).toFixed(1)}% of the informative markers read heterozygous, over the `
        + `${(100 * maxRegionHet).toFixed(0)}% ceiling. Where the loaded parent is homozygous a `
        + 'biparental sample can only be heterozygous when the other parent transmitted the allele '
        + 'this one lacks, so a rate this high means the genotypes are not measuring the region. '
        + 'No origin is called from them',
    }
  }

  const top = Math.max(...logs)
  const w = logs.map((x) => Math.exp(x - top))
  const sum = w[0] + w[1] + w[2]
  const post = w.map((x) => x / sum)
  const names: OneParentVerdict[] = ['both-present', 'known-parent-lost', 'other-parent-lost']
  let best = 0
  for (let i = 1; i < 3; i += 1) if (post[i] > post[best]) best = i

  if (post[best] < callPosterior) {
    return {
      ...base, verdict: 'refused', posterior: post[best], band: 'D' as Band,
      why: `best hypothesis reaches ${post[best].toFixed(3)}, under the ${callPosterior} needed`,
    }
  }
  const rate = n ? exclusive / n : 0
  return {
    ...base,
    verdict: names[best],
    posterior: post[best],
    band: bandObligateHet(post[best]),
    why: `${names[best]} at posterior ${post[best].toFixed(4)} over ${n} markers where the loaded `
      + `parent is homozygous. ${exclusive} carry the allele that parent does not have `
      + `(${(rate * 100).toFixed(1)}%), which only its absence explains, and ${het} are `
      + `heterozygous, which only both copies being present explains`,
  }
}

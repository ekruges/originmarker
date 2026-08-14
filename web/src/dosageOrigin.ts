/**
 * Which parent's copy is missing, read from ALLELE DOSAGE rather than from genotype calls.
 *
 * WHY THIS EXISTS. `oneParentOrigin.ts` answers the same question from genotypes and is the better
 * instrument where genotypes exist. It cannot answer it for a whole-chromosome loss, and the
 * reason is structural rather than a limitation of the code: a whole chromosome is DETECTED by the
 * collapse of its genotype call rate, and the origin is CALLED from genotypes, so the signal that
 * identifies the event is the destruction of the evidence that would assign it. Measured on
 * GSE148488: all four segmental losses scored, all three whole-chromosome losses refused, and on
 * one of them chromosome 1 was 69% no-call with the survivors reading 59.7% heterozygous where the
 * father is homozygous, which is not a genome.
 *
 * The B-allele frequency is measured at every marker whether or not the caller emits a genotype.
 * It survives the collapse. So this reads the same Mendelian fact out of the channel that is
 * still there.
 *
 * THE MENDELIAN FACT, and it is the same one the genotype version rests on. Take a marker where
 * the loaded parent is HOMOZYGOUS, and orient so that parent's own allele is A and the allele it
 * lacks is B. The parent can transmit only A. So:
 *
 *     both copies present     the other parent gave B (rate q) -> AB, dosage in the MIDDLE band
 *                             the other parent gave A          -> AA, dosage LOW
 *     the loaded parent's     only the other parent's allele is there, one copy:
 *     copy is absent            it gave B (rate q)             -> B alone, dosage HIGH
 *                               it gave A                      -> A alone, dosage LOW
 *     the other parent's      only this parent's allele is there, which is A:
 *     copy is absent                                           -> dosage LOW, always
 *
 * Two of those cells are impossible under the others, which is where all the power is. A HIGH
 * dosage at a marker where the loaded parent is homozygous requires that parent's copy to be
 * ABSENT, because if it were present it would contribute an A and hold the dosage at or below the
 * middle. A MIDDLE dosage requires BOTH copies, because one copy cannot be heterozygous.
 *
 * THIS IS NOT AN ALLELE SHARE AND IT IS NOT CENTRED. Two earlier attempts in this project formed
 * a paternal share and centred it on the sample's own median, and both produced one-directional
 * results. For a one-parent informative set the expectation is 1 - q/2 under two copies and 1
 * under loss-of-unknown, so for every q below one half the majority of markers sit at the ceiling:
 * the median IS the ceiling, upward headroom is zero by construction, and every call came back the
 * same way. What is counted here is the RATE OF AN EVENT THAT IS IMPOSSIBLE UNDER THE ALTERNATIVE,
 * in each direction separately. A quantity that is never centred cannot be mis-centred, and the
 * two directions are counted by the same code with the orientation flipped, so neither can be
 * given a headroom the other lacks.
 *
 * WHERE IT IS WEAKER THAN THE GENOTYPE VERSION, and it is weaker. Dosage is a continuous
 * measurement and amplification distorts it: preferential amplification of one allele drags a
 * true middle toward an extreme, and a poorly amplified locus drags an extreme toward the middle.
 * The second direction is the safe one, since it moves evidence from the decisive bands into the
 * band every hypothesis allows, so noise costs power rather than buying a wrong answer. The first
 * is not safe, which is why the extreme bands here are narrow and why `MIN_MARKERS_DOSAGE` is far
 * above the genotype caller's floor: this is intended for whole chromosomes and large segments,
 * where there are thousands of markers, not for the small regions the genotype channel handles.
 */
import type { AB } from './informativity.ts'

/**
 * Dosage bands, taken from the ones this project already uses on the same channel.
 *
 * `tallyRow` counts a BAF as extreme below 0.15 or above 0.85, and as in-band between 0.35 and
 * 0.65. Reusing those means a marker is classified the same way here as everywhere else, rather
 * than by a second set of thresholds that would drift from the first.
 */
export const BAND_EXTREME = 0.15
export const BAND_MID_LO = 0.35
export const BAND_MID_HI = 0.65

/**
 * Probability a marker lands outside the band its true state predicts.
 *
 * Covers everything that moves a dosage: amplification bias, probe cross-hybridisation, and the
 * clustering that produced the reading. Set at the drop-in rate measured on this platform, since
 * a false heterozygote and a dosage that has wandered into the middle band are the same physical
 * event seen through two channels.
 */
export const DOSAGE_NOISE = 0.0435

/**
 * Informative markers a region needs before a dosage verdict.
 *
 * Far above the genotype caller's 50. A dosage band is a noisier observation than a genotype call
 * and the decisive bands are entered at rate q rather than at rate 1, so the evidence per marker
 * is lower. This channel exists for whole chromosomes and multi-megabase segments, which carry
 * thousands of informative markers; a region that cannot reach this floor is a region the genotype
 * channel should be answering.
 */
export const MIN_MARKERS_DOSAGE = 400

/** Posterior a hypothesis must reach before it is named. Same bar as the genotype channel. */
export const DOSAGE_POSTERIOR = 0.95

/**
 * How far the region's unresolved fraction may exceed the sample's own genome before it is refused.
 *
 * THIS IS THE GUARD THAT ALMOST DID NOT GET WRITTEN, and the reason it exists is worth keeping.
 * The middle band is impossible with one copy, so middle-band markers are enormously informative
 * for two copies being present. That makes this channel vulnerable in exactly the way the genotype
 * channel was: a region whose intensity is noise scatters its dosages across the range, a large
 * share land in the middle band, and the likelihood reads that as overwhelming evidence of two
 * copies. Measured on a real blastomere, chromosome 1 against the rest of its own genome:
 *
 *     chr1              call 30.7%   own  7.7%   middle 49.7%   excluded 5.4%   between 37.2%
 *     its other chroms  call 83.7%   own 88.1%   middle  4.0%   excluded 0.9%   between  7.1%
 *
 * Twelve times the middle band and five times the between band. It returned both-copies-present at
 * posterior 1.0000. The tell is the BETWEEN band, the one no hypothesis puts mass in: a genome
 * resolves its alleles into clusters, and a region that does not is not being measured. Compared
 * against the sample's OWN genome rather than an absolute, because amplification sets that
 * baseline per array and it varies enormously between a blastomere and bulk DNA.
 *
 * WHAT THIS GUARD CANNOT SEPARATE, and it matters for which samples it should be trusted on. An
 * unresolved dosage has two causes. In a SINGLE cell it is noise: one genome cannot be partly
 * anything, so a reading between the bands is the measurement failing. In a MULTI-CELL biopsy it
 * can be real. A mosaic loss carried by some cells and not others puts the biopsy's average dosage
 * genuinely between the clean bands, which is the signal a mosaicism caller would want. This guard
 * treats both as unmeasurable and refuses. That is right for a blastomere and conservative for a
 * trophectoderm biopsy, where it will refuse regions that are mosaic rather than broken. Measured:
 * two trophectoderm biopsies of one embryo showed 29.9% and 23.7% unresolved over a chr16 loss
 * against a 7% genome background, and both were refused here while the genotype channel called
 * them both-present. Separating the two needs the per-cell dosage distribution rather than the
 * biopsy mean, which this function does not have.
 */
/**
 * Frequency of the allele the loaded parent lacks, among alleles the other parent transmits.
 *
 * Shared with the genotype channel so the two cannot disagree about the same population quantity
 * while claiming to answer the same question.
 */
export const DEFAULT_DOSAGE_Q = 0.30

export const MAX_BETWEEN_RATIO = 3
export const MAX_BETWEEN_FLOOR = 0.15

export type DosageVerdict =
  | 'known-parent-lost'
  | 'other-parent-lost'
  | 'both-present'
  | 'refused'

export interface DosageCall {
  verdict: DosageVerdict
  posterior: number
  /** Markers where the loaded parent is homozygous AND a dosage was read. */
  markers: number
  /** Dosage at the extreme the loaded parent CANNOT produce. Only its absence explains these. */
  excluded: number
  /** Dosage in the middle band. Only two copies explain these. */
  middle: number
  /** Dosage at the loaded parent's own extreme. Every hypothesis allows these. */
  own: number
  /** Dosage in no band at all. No hypothesis predicts these, so a high rate means the intensity
   *  is not resolving alleles in this region rather than that the region is unusual. */
  between: number
  q: number
  why: string
}

/** LOW is the loaded parent's own allele, MID is two copies, HIGH is the allele it cannot give. */
type Band = 'own' | 'middle' | 'excluded' | 'between'

/**
 * Which band a dosage falls in, ORIENTED so the loaded parent's own allele is always low.
 *
 * The flip is the whole reason the two directions cannot be given different sensitivity: a
 * BB parent is handled by the identical code path as an AA parent with the axis reversed, so any
 * asymmetry would have to be written twice to survive.
 */
export function band(parent: AB, baf: number): Band {
  if (!Number.isFinite(baf)) return 'between'
  const b = parent === 'BB' ? 1 - baf : baf
  if (b < BAND_EXTREME) return 'own'
  if (b > 1 - BAND_EXTREME) return 'excluded'
  if (b >= BAND_MID_LO && b <= BAND_MID_HI) return 'middle'
  return 'between'
}

/**
 * Per-marker likelihood under the three hypotheses, in the order
 * [both copies present, the loaded parent's copy lost, the other parent's copy lost].
 *
 * The zeros are the point. A middle dosage is impossible with one copy, and an excluded-extreme
 * dosage is impossible while the loaded parent's copy is there. Noise keeps them from being
 * literally zero, so a single stray marker cannot veto a hypothesis outright, but the ratio is
 * what carries the call.
 */
function likelihood(b: Band, q: number, noise: number): [number, number, number] {
  const raw: Record<Band, [number, number, number]> = {
    // Two copies: middle at rate q, own-extreme otherwise. Never the excluded extreme.
    // One copy from the other parent: excluded extreme at rate q, own otherwise. Never middle.
    // One copy from this parent: its own allele, always.
    own: [1 - q, 1 - q, 1],
    middle: [q, 0, 0],
    excluded: [0, q, 0],
    // Between the bands. Every hypothesis produces these at the noise rate and none is favoured.
    between: [0, 0, 0],
  }
  const p = raw[b]
  // Spread `noise` of the mass uniformly over the four bands so nothing is impossible and a
  // single distorted marker cannot decide anything.
  return [
    p[0] * (1 - noise) + noise / 4,
    p[1] * (1 - noise) + noise / 4,
    p[2] * (1 - noise) + noise / 4,
  ]
}

/**
 * Call which copy is missing across a region, from one genotyped parent and the sample's dosage.
 *
 * `pairs` are (parent genotype, sample B-allele frequency) at each marker. Markers where the
 * parent is heterozygous carry nothing and are dropped: every hypothesis predicts the same dosage
 * distribution there, so including them would dilute the evidence with noise.
 */
export function callDosageOrigin(
  pairs: readonly (readonly [AB, number | null])[],
  q: number = DEFAULT_DOSAGE_Q,
  noise = DOSAGE_NOISE,
  opts: {
    minMarkers?: number
    posterior?: number
    /**
     * The same band profile computed over the REST of this sample's genome, excluding the region
     * under test. Without it the unresolved-dosage guard cannot run, because there is nothing to
     * compare against, and the call is made without it. Supply it wherever it can be computed.
     */
    background?: { between: number }
    maxBetweenRatio?: number
  } = {},
): DosageCall {
  const minMarkers = opts.minMarkers ?? MIN_MARKERS_DOSAGE
  const need = opts.posterior ?? DOSAGE_POSTERIOR
  const betweenRatio = opts.maxBetweenRatio ?? MAX_BETWEEN_RATIO

  const logs: [number, number, number] = [0, 0, 0]
  let n = 0
  let excluded = 0
  let middle = 0
  let own = 0
  let between = 0
  for (const [parent, baf] of pairs) {
    if (parent !== 'AA' && parent !== 'BB') continue
    if (baf === null || !Number.isFinite(baf)) continue
    n += 1
    const b = band(parent, baf)
    if (b === 'excluded') excluded += 1
    else if (b === 'middle') middle += 1
    else if (b === 'own') own += 1
    else between += 1
    const l = likelihood(b, q, noise)
    logs[0] += Math.log(l[0])
    logs[1] += Math.log(l[1])
    logs[2] += Math.log(l[2])
  }

  const base = { markers: n, excluded, middle, own, between, q }
  if (n < minMarkers) {
    return {
      ...base, verdict: 'refused', posterior: NaN,
      why: `${n} markers carry a dosage where the loaded parent is homozygous, under the `
        + `${minMarkers} this channel needs. Dosage is a noisier observation than a genotype call, `
        + 'so its floor is higher; a region this small should be read from genotypes',
    }
  }

  // The intensity must be resolving alleles before its likelihood means anything. Compared
  // against this sample's own genome, since amplification sets that baseline per array.
  const betweenRate = between / n
  const bg = opts.background?.between
  const ceiling = bg === undefined ? Infinity
    : Math.max(MAX_BETWEEN_FLOOR, bg * betweenRatio)
  if (betweenRate > ceiling) {
    return {
      ...base, verdict: 'refused', posterior: NaN,
      why: `${(100 * betweenRate).toFixed(1)}% of the dosages sit in no band at all, against `
        + `${(100 * (bg as number)).toFixed(1)}% over the rest of this sample's genome. No `
        + 'hypothesis predicts an unresolved dosage, so a rate this far above the sample\'s own '
        + 'background means the intensity is not resolving alleles here rather than that the '
        + 'region is unusual. A middle-band reading is impossible with one copy, so scattered '
        + 'noise would otherwise read as strong evidence of two',
    }
  }

  const top = Math.max(...logs)
  const w = logs.map((l) => Math.exp(l - top))
  const sum = w[0] + w[1] + w[2]
  const post = w.map((x) => x / sum)
  const best = post.indexOf(Math.max(...post))
  const names: DosageVerdict[] = ['both-present', 'known-parent-lost', 'other-parent-lost']
  const pctOf = (x: number) => `${((100 * x) / n).toFixed(1)}%`

  if (post[best] < need) {
    return {
      ...base, verdict: 'refused', posterior: post[best],
      why: `best hypothesis reaches ${post[best].toFixed(3)}, under the ${need} needed. `
        + `${excluded} markers (${pctOf(excluded)}) sit at the dosage the loaded parent cannot `
        + `produce and ${middle} (${pctOf(middle)}) sit in the two-copy band`,
    }
  }
  return {
    ...base,
    verdict: names[best],
    posterior: post[best],
    why: `${names[best]} at posterior ${post[best].toFixed(4)} from allele dosage over ${n} `
      + `markers where the loaded parent is homozygous. ${excluded} (${pctOf(excluded)}) sit at `
      + 'the dosage that parent cannot produce, which only its absence explains, and '
      + `${middle} (${pctOf(middle)}) sit in the band that needs two copies. Read from intensity `
      + 'rather than genotype calls, so a collapsed call rate does not remove the evidence',
  }
}

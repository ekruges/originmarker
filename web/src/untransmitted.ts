/**
 * The untransmitted parental haplotype, and the mechanism of a trisomy, from markers the
 * obligate-het channel throws away.
 *
 * TWO CHANNELS, ONE MARKER CLASS. Everything this project has called origin from until now uses
 * markers where the loaded parent is HOMOZYGOUS. This module uses the complement, markers where
 * that parent is HETEROZYGOUS, which is a disjoint set of comparable size that was simply unused:
 * 27,302 to 58,571 markers per array, 86.2 to 99.8% of what is informative.
 *
 * WHY THE UNUSED SET IS THE BETTER ONE, and it is not what I assumed. The gain is not drift
 * cancellation and not a larger denominator. It is that the marker set is CLEAN BY CONSTRUCTION.
 * In the parent-homozygous window only 32 to 90% of markers are truly heterozygous in the child,
 * because the other parent frequently transmits the same allele and those markers carry no dosage
 * information at all while still counting toward the denominator. At a marker where the loaded
 * parent is heterozygous and the child reads homozygous, the transmission is determined: the
 * parent gave the allele the child shows, and did not give the other. Every such marker carries
 * information. Measured signal-to-noise is 1.40 to 1.98x above the obligate-het channel, on every
 * amplified material at every mosaic fraction tested.
 *
 * AND IT IS THE ONLY CHANNEL THAT GIVES A SINGLE BLASTOMERE A FLOOR AT ALL. Against one parent, a
 * blastomere loss has no detection floor at any fraction to 0.70 on the obligate-het channel and
 * still none at eight cells. Through this channel it has one, 0.628. That is still too high to
 * call, so the honest statement is that this converts a material impossibility into a measurable
 * quantity rather than into an answer. It is the difference between "no array of this kind can
 * answer" and "this needs a fraction above 0.63", which is a different sentence to a reader.
 *
 * WHAT IS MEASURED HERE AND WHAT IS INHERITED. The marker selection, the orientation and the
 * impossible-reading rate below are computed from the arrays. The floors and the SNR ratio come
 * from an external methods review (audit/ORIGIN-CONSULT.md) which measured them on 135 arrays of
 * this platform, and are inherited rather than re-derived. That distinction is kept explicit
 * because the two carry different warranties.
 */
import type { AB } from './informativity.ts'

/**
 * Floors for this channel, from the review. NaN means no floor at any fraction to 0.70.
 *
 * The blastomere entry is the reason this module exists: the obligate-het channel has no floor
 * there at all, in either parent configuration.
 */
export const F80_UNTRANSMITTED = {
  onePparent: { bulk: 0.040, 'esc-single': 0.313, trophectoderm: 0.442, blastomere: NaN },
  twoParents: { bulk: 0.036, 'esc-single': 0.232, trophectoderm: 0.186, blastomere: 0.628 },
} as const

/** Measured advantage over the obligate-het channel, for reporting rather than for computing. */
export const SNR_GAIN_RANGE = [1.40, 1.98] as const

export interface UntransmittedPair {
  /** The allele the loaded parent did NOT transmit, at this marker. */
  untransmitted: 'A' | 'B'
  /** The sample's B-allele frequency there. */
  baf: number
}

/**
 * Select the markers this channel runs on, and record which allele went untransmitted.
 *
 * Informative only where the loaded parent is HETEROZYGOUS and the sample reads HOMOZYGOUS. A
 * heterozygous sample is ambiguous: it could hold the parent's A with the other parent's B, or the
 * parent's B with the other parent's A, and nothing distinguishes them without the second parent.
 * Those markers are dropped rather than guessed at.
 */
export function untransmittedPairs(
  pairs: readonly (readonly [AB, AB, number | null])[],
): { pairs: UntransmittedPair[]; considered: number; ambiguous: number } {
  const out: UntransmittedPair[] = []
  let considered = 0
  let ambiguous = 0
  for (const [parent, sample, baf] of pairs) {
    if (parent !== 'AB') continue
    considered += 1
    if (sample === 'AB') { ambiguous += 1; continue }
    if (sample !== 'AA' && sample !== 'BB') continue
    if (baf === null || !Number.isFinite(baf)) continue
    // The sample shows one allele homozygously, so the parent transmitted THAT one and withheld
    // the other. This is the whole inference and it needs no second parent.
    out.push({ untransmitted: sample === 'AA' ? 'B' : 'A', baf })
  }
  return { pairs: out, considered, ambiguous }
}

/**
 * Fraction of readings that no genotype can produce, which is a direct per-array error rate.
 *
 * At a marker where the loaded parent is heterozygous, a sample homozygous for one allele should
 * carry dosage at that allele's extreme. A reading at the OPPOSITE extreme contradicts its own
 * genotype call, so the rate of those measures how far the two channels disagree on this array
 * without needing any external truth. Measured range on this platform is 1.7 to 13.8%.
 */
export function impossibleRate(ps: readonly UntransmittedPair[]): number {
  if (!ps.length) return NaN
  let bad = 0
  for (const p of ps) {
    // Untransmitted B means the sample called AA, so its dosage should sit low.
    const contradicts = p.untransmitted === 'B' ? p.baf > 0.85 : p.baf < 0.15
    if (contradicts) bad += 1
  }
  return bad / ps.length
}

/**
 * Oriented dosage on this channel, so that the UNTRANSMITTED allele always reads high.
 *
 * The untransmitted allele is the one the loaded parent held back. It reaches the sample only
 * through the OTHER parent, so its dosage rising is the same evidence a rising parent-allele share
 * is on the obligate-het channel: the loaded parent's contribution is short. Orienting here rather
 * than at the call site is what keeps the two directions symmetric.
 */
export const orientUntransmitted = (p: UntransmittedPair): number =>
  (p.untransmitted === 'B' ? p.baf : 1 - p.baf)

// ---------------------------------------------------------------------------------------------
// Trisomy mechanism: both parental homologues, or one duplicated.

/**
 * What the occupancy statistic can actually name, which is less than the usual vocabulary implies.
 *
 * THIS TYPE USED TO SAY 'SPH' AND THAT WAS A CLAIM IT COULD NOT SUPPORT. "Single parental
 * homologue" is read by every embryologist as "post-zygotic mitotic", and the rule behind that
 * reading is wrong. Meiosis II non-disjunction WITH a crossover produces both-homologue tracts
 * distal to the centromere, so BPH does not separate meiotic from mitotic; only MII WITHOUT
 * recombination gives a chromosome-wide single-homologue pattern. The primary literature
 * acknowledges the collapse rather than resolving it, and nobody separates the two on genotype
 * alone at any material quality.
 *
 *   APCAD2 2025: "SPH trisomies are expected to be predominantly mitotic in origin ... even
 *   though a meiotic segregation error in meiosis II without recombination cannot be excluded."
 *
 * So the honest second category is NOT-BPH, and what it pools is stated wherever it is reported:
 * non-recombinant MII together with post-zygotic mitotic duplication. Naming it 'SPH' invited a
 * reader to collapse that pooling into a mechanism the data does not choose.
 */
export type Mechanism = 'BPH' | 'not-BPH' | 'unresolved'

/**
 * Sex-specific recombination rate, cM/Mb, from the deCODE 1-Mb map (UCSC hg19 recombRate).
 *
 * Carried because it says WHERE this class of reasoning is safest, and the asymmetry is large
 * enough to change how a result should be read. Suppression near the centromere is far stronger in
 * male meiosis: the female-to-male ratio rises from 1.63 genome-wide to 4.79 pericentromerically.
 * Any argument resting on retained pericentromeric heterozygosity is therefore substantially safer
 * for a paternal-origin event than a maternal one, and the informative window is wider on the male
 * side. Reported alongside the call rather than used to compute it.
 */
export const RECOMB_CM_PER_MB = {
  genomeWide: { female: 1.598, male: 0.982 },
  pericentromeric2Mb: { female: 0.660, male: 0.138 },
} as const

/**
 * Band centres for a trisomy at markers where the loaded parent is heterozygous.
 *
 * BOTH-PARENTAL-HOMOLOGUES means the extra copy is the parent's OTHER homologue, so the child
 * holds one of each of that parent's alleles plus one from the other parent: the other parent's
 * allele share is always 1/3, and the distribution has ONE band.
 *
 * SINGLE-PARENTAL-HOMOLOGUE means one homologue was duplicated, so the child holds two copies of
 * one parental allele plus one from the other parent: the share is 1/3 or 1.0 depending on which
 * allele was doubled, and the distribution has TWO bands, NEITHER at 2/3.
 *
 * THE DIAGNOSTIC IS OCCUPANCY, NOT A MEAN, and that is what makes it work on this material. A
 * per-marker three-band assignment needs BAF spread below 0.053 and fails on every amplified class
 * with 13.3 to 21.8% misassignment. A fractional occupancy over hundreds of markers is insensitive
 * to per-marker dispersion, because it asks where mass sits rather than which band a marker is in.
 * AUC 1.000 from 400 markers with zero contamination between the two.
 */
export const BAND_BPH = 2 / 3
export const BAND_SPH_LOW = 1 / 3
export const BAND_SPH_HIGH = 1.0
/** Half-width around a band centre. Wide, since only the OCCUPANCY has to be right. */
export const BAND_HALF_WIDTH = 0.12
/** Markers below which the occupancy statistic is not attempted. AUC reaches 1.000 from here. */
export const MIN_MARKERS_MECHANISM = 400

export interface MechanismCall {
  mechanism: Mechanism
  /** Fraction of markers whose share sits at the both-homologues band. */
  atBph: number
  /** Fraction at either single-homologue band. */
  atSph: number
  markers: number
  why: string
}

/**
 * Call the mechanism of an established trisomy from band occupancy.
 *
 * REPORTS BPH POSITIVELY AND SPH ONLY BY EXCLUSION, which the measurements require rather than
 * being a stylistic choice. The single-homologue bands at 1/3 and 1.0 are both populated by an
 * ordinary euploid genome, so occupancy there is not evidence of anything on its own: its power
 * against euploid at a 1% false-positive rate is 0.134. The both-homologues band at 2/3 is not
 * populated by a euploid genome, so occupancy there is decisive. `pairs` must already be
 * established as copy number three by the dosage channel; this answers only which mechanism.
 */
export function callMechanism(
  ps: readonly UntransmittedPair[],
  opts: { minMarkers?: number; halfWidth?: number; copyNumberThree?: boolean } = {},
): MechanismCall {
  const need = opts.minMarkers ?? MIN_MARKERS_MECHANISM
  // COPY NUMBER THREE MUST ALREADY BE ESTABLISHED, and this gate is not a formality. The two
  // single-homologue bands are populated by an ordinary euploid genome, so run on a normal
  // chromosome this statistic reports SPH every time, by exclusion, and is wrong every time.
  // Caught exactly that way: a euploid chromosome of a real parent-child pair came back "SPH".
  if (opts.copyNumberThree === false) {
    return {
      mechanism: 'unresolved', atBph: NaN, atSph: NaN, markers: ps.length,
      why: 'copy number three is not established for this interval, and the mechanism question '
        + 'only exists once it is. The two single-homologue bands are also populated by a euploid '
        + 'genome, so asking this of a normal chromosome returns SPH by exclusion and means '
        + 'nothing',
    }
  }
  const w = opts.halfWidth ?? BAND_HALF_WIDTH
  let bph = 0
  let sph = 0
  let n = 0
  for (const p of ps) {
    // Share of the OTHER parent's contribution, which is what the band centres are stated in.
    const share = orientUntransmitted(p)
    n += 1
    if (Math.abs(share - BAND_BPH) <= w) bph += 1
    else if (Math.abs(share - BAND_SPH_LOW) <= w || Math.abs(share - BAND_SPH_HIGH) <= w) sph += 1
  }
  const atBph = n ? bph / n : NaN
  const atSph = n ? sph / n : NaN
  const base = { atBph, atSph, markers: n }
  if (n < need) {
    return {
      ...base, mechanism: 'unresolved',
      why: `${n} markers where the loaded parent is heterozygous, under the ${need} this needs. `
        + 'The statistic is an occupancy over many markers rather than a per-marker call, so it '
        + 'needs the count rather than the precision',
    }
  }
  // A euploid genome puts mass at 1/3 and 1.0 but none at 2/3, so occupancy at 2/3 is the only
  // positive evidence available here.
  if (atBph >= 0.25) {
    return {
      ...base, mechanism: 'BPH',
      why: `${(100 * atBph).toFixed(1)}% of markers sit at the both-homologues band, which a `
        + 'euploid genome does not populate and a duplicated single homologue cannot reach. Both '
        + "of the loaded parent's homologues are present",
    }
  }
  return {
    ...base, mechanism: 'not-BPH',
    why: `${(100 * atBph).toFixed(1)}% of markers sit at the both-homologues band, too few for `
      + `both homologues, and ${(100 * atSph).toFixed(1)}% sit at the two single-homologue bands. `
      + 'Reported by EXCLUSION rather than positively: those two bands are also populated by an '
      + 'ordinary euploid genome, so their occupancy has power 0.134 against euploid on its own '
      + 'and is only meaningful once copy number three is already established. This class POOLS '
      + 'non-recombinant meiosis II with post-zygotic mitotic duplication and does not choose '
      + 'between them: a meiosis II error with a crossover would have shown both-homologue tracts '
      + 'distal to the centromere, so the absence of them excludes only the recombinant case',
  }
}

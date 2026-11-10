/**
 * Whether an array's genotype calls carry any information about any genome.
 *
 * WHAT THIS IS FOR. A refusal can say that too few markers were informative, and it can say that
 * the call rate collapsed. Neither says that the calls themselves are noise. An array whose
 * genotypes are noise still returns a call rate, still returns markers, and reads above the
 * diploid heterozygosity ceiling, where the ceiling's own caveat is that a mixture, a
 * contaminated reaction and a destroyed array cannot be told apart. This separates them.
 *
 * THE STATISTIC. Over the markers where the sample and a reference are both homozygous, measure
 * how often they are homozygous for OPPOSITE alleles, and divide by the rate expected if the
 * sample's calls said nothing about the reference. The expected rate comes from the two arrays'
 * own AA fractions over that same marker set, so the null is carried by the data rather than
 * assumed, and R is dimensionless:
 *
 *     R = P(opposite) / [ s(1 - r) + (1 - s)r ]
 *
 *     R near 1          the calls say nothing about this reference, or about any genome
 *     R well below 1    the calls are informative
 *
 * THE REFERENCE NEED NOT BE A PARENT, and that is what makes this usable. Two unrelated people
 * still share population structure, so an informative array sits well below 1 against a stranger,
 * while an array whose calls are noise sits at 1 against everyone.
 *
 * A DNA MIXTURE READS LOW, NOT HIGH, by construction: a homozygous call in a mixture requires
 * every contributor to carry that allele, so the surviving homozygotes agree with each
 * contributor more often than chance, not less.
 *
 * THE REFERENCE MUST BE MEASURED, NEVER RECONSTRUCTED. A reconstructed reference is homozygous at
 * every marker by construction, and one built from this sample's own relatives shares its calls
 * for a reason that is not array quality.
 *
 * WHY NOT hetBafSd. It measures how much amplification a sample went through, not how good the
 * array is, and arrays with no genotype signal read BELOW most usable amplified ones: 0.158 to
 * 0.168 against 0.173 to 0.345 over the 112-array Axiom corpus. An upper threshold on it refuses
 * good single cells and passes every failed array.
 *
 * Bands and the corpus they were fixed on: audit/noise-information-ratio-2026-11-08.txt.
 */
import type { AB } from './informativity.ts'

/** At or below this, the calls are informative. Not shipped as a positive claim: see the run note. */
export const INFORMATIVE_MAX = 0.30
/** At or above this, the calls carry no information. The cut the shipped reason fires on. */
export const NO_INFORMATION_MIN = 0.60
/** A reference whose A/B coding is flipped against the sample reads near 2, not 1. */
export const FLIPPED_CODING_MIN = 1.5
/** Below this many shared homozygous markers the ratio is not estimated. */
export const MIN_SHARED_HOMOZYGOTES = 500

/** The four counts the ratio needs, over markers homozygous in BOTH arrays. */
export interface HomTally {
  both: number
  opp: number
  sampleAA: number
  refAA: number
}

export const emptyHomTally = (): HomTally => ({ both: 0, opp: 0, sampleAA: 0, refAA: 0 })

/**
 * Add one marker to a tally. Anything heterozygous or uncalled on either side carries no
 * opposite-homozygote information and is skipped.
 */
export const addHomozygous = (t: HomTally, sample: AB | null, reference: AB | null): void => {
  if (sample !== 'AA' && sample !== 'BB') return
  if (reference !== 'AA' && reference !== 'BB') return
  t.both += 1
  if (sample === 'AA') t.sampleAA += 1
  if (reference === 'AA') t.refAA += 1
  if (sample !== reference) t.opp += 1
}

/** The information ratio, or NaN where there are too few shared homozygotes to estimate it. */
export const informationRatio = (t: HomTally): number => {
  if (t.both < MIN_SHARED_HOMOZYGOTES) return NaN
  const s = t.sampleAA / t.both
  const r = t.refAA / t.both
  const expected = s * (1 - r) + (1 - s) * r
  return expected > 0 ? (t.opp / t.both) / expected : NaN
}

export type InformationBand =
  | 'informative' | 'unresolved' | 'no-information' | 'flipped-coding' | 'too-few-markers'

export const informationBand = (r: number): InformationBand =>
  !Number.isFinite(r) ? 'too-few-markers'
    : r >= FLIPPED_CODING_MIN ? 'flipped-coding'
      : r >= NO_INFORMATION_MIN ? 'no-information'
        : r <= INFORMATIVE_MAX ? 'informative' : 'unresolved'

/**
 * The array-level reason, or undefined where the calls still carry information.
 *
 * Says what the failure is NOT, because that is the part the heterozygosity ceiling cannot do,
 * and does not offer a rerun: across 29 groups of replicate arrays in public data on this
 * platform, replicates of one amplification product never landed on opposite sides of these
 * gates.
 */
export const whyNoInformation = (r: number, reference: string): string | undefined => {
  const band = informationBand(r)
  if (band === 'flipped-coding') {
    return `this array and ${reference} disagree at opposite homozygotes ${r.toFixed(2)} times `
      + 'as often as unrelated calls would, which is the reading an inverted A/B coding gives '
      + 'rather than a failure. The two files should be checked against one another before '
      + 'anything here is believed'
  }
  if (band !== 'no-information') return undefined
  return `this array's genotype calls carry no information about any genome: against `
    + `${reference} they disagree at opposite homozygotes at ${r.toFixed(2)} of the rate `
    + 'expected from calls that say nothing, where 1.0 is that rate exactly. This is not a '
    + 'mixture of two people, which reads below 1 rather than at it, and it is not a limit of '
    + 'the markers, which are there. Where the same amplified material has been run twice, the '
    + 'failure reproduced every time, so running that product again is unlikely to help'
}

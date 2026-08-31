/**
 * Whether a sample is intact enough for anything else on the page to mean what it says.
 *
 * A SAMPLE-LEVEL FLAG, NOT A REFUSAL, and the difference is the point. The tool already stops
 * looking for events on an array that failed its own stage inference, and already emits a
 * `complex` finding when a genome has no undisturbed remainder to reference against. Neither of
 * those reaches a reader as a statement about the SAMPLE: the first produces an empty result that
 * looks like a clean one, and the second is a single row among however many others the run found.
 * An operator scanning a run needs to see, per sample and before reading anything else, that this
 * one is wrecked.
 *
 * NO NEW THRESHOLDS. Every boundary here is one the project already measured for a different
 * purpose, reused rather than re-chosen:
 *   COMPLEX_DEVIANT_FRACTION 0.30  the share of deviant autosomes past which self-reference fails,
 *                                  measured on a series where all 70 whole-chromosome losses sat
 *                                  above it and none was assignable.
 *   COMPLEX_CALL_RATE 0.70         the call rate past which the same applies.
 *   stage 'failed'                 the array did not resolve to any material, which the stage
 *                                  inference decides on its own measured quality gates.
 * Adding a fourth number chosen by eye would put this flag in the class of thresholds this project
 * has repeatedly had to withdraw.
 */
import { COMPLEX_DEVIANT_FRACTION, COMPLEX_CALL_RATE } from './abnormalities.ts'
import type { Stage } from './stage.ts'

export type Integrity =
  /** No change found, and the array met its own quality gates. */
  | 'intact'
  /** Changes found, and the genome still has an undisturbed remainder to measure them against. */
  | 'disturbed'
  /**
   * Nothing on this array can be trusted. Either it did not resolve to a material at all, or so
   * much of the genome is deviant that every self-referenced statistic is measuring the reference.
   */
  | 'uninterpretable'

export interface IntegrityCall {
  level: Integrity
  /** The one-line flag an operator reads first. */
  headline: string
  why: string
  /** What must not be done with this sample's rows. Empty when the sample is usable. */
  limits: string
  measures: {
    affectedAutosomes: number
    totalAutosomes: number
    deviantFraction: number
    callRate: number
    stage: Stage
    events: number
  }
}

export interface IntegrityInput {
  stage: Stage
  callRate: number
  /** Autosomes carrying any change: aneuploidy, segment, or taxonomy finding. */
  affectedAutosomes: number
  totalAutosomes: number
  /** Total changes reported, which is what the reader will see listed. */
  events: number
}

/**
 * The flag.
 *
 * Order matters. A failed array is uninterpretable whatever its event count says, because the
 * event count was produced by the same broken measurement. A genome past the deviant-fraction
 * bound is uninterpretable even if the array itself is pristine, because the fault is in the
 * genome rather than the reaction, and those are different sentences to write.
 */
export function assessIntegrity(input: IntegrityInput): IntegrityCall {
  const { stage, callRate, affectedAutosomes, totalAutosomes, events } = input
  const deviantFraction = totalAutosomes > 0 ? affectedAutosomes / totalAutosomes : NaN
  const measures = {
    affectedAutosomes, totalAutosomes, deviantFraction, callRate, stage, events,
  }

  if (stage === 'failed') {
    return {
      level: 'uninterpretable',
      headline: 'DAMAGED SAMPLE: this array did not resolve to any material',
      why: 'The stage inference rejected this array on its own quality gates, so there is no '
        + 'material whose dropout could parameterise a likelihood and no genome the statistics '
        + 'could describe. Any change listed for this sample was found by a measurement that is '
        + 'not measuring a genome.',
      limits: 'Nothing on this sample is reportable. Do not read its changes, its parental calls '
        + 'or its counts, and do not include it in any comparison across samples.',
      measures,
    }
  }

  if (Number.isFinite(callRate) && callRate < COMPLEX_CALL_RATE) {
    return {
      level: 'uninterpretable',
      headline: `DAMAGED SAMPLE: call rate ${(callRate * 100).toFixed(1)}% is below the `
        + `${(COMPLEX_CALL_RATE * 100).toFixed(0)}% this tool can reference against`,
      why: 'Every origin statistic here is measured against the rest of this array\'s own genome. '
        + 'Below this call rate there are not enough calls left for that reference to be stable, '
        + 'so a measured difference describes the reference as much as the target.',
      limits: 'Do not read parental calls from this sample. Positions may still be indicative, '
        + 'because position comes from intensity and needs no genotype, but they are not reportable.',
      measures,
    }
  }

  if (Number.isFinite(deviantFraction) && deviantFraction >= COMPLEX_DEVIANT_FRACTION) {
    return {
      level: 'uninterpretable',
      headline: `DAMAGED SAMPLE: ${affectedAutosomes} of ${totalAutosomes} autosomes are affected`,
      why: `A share of ${deviantFraction.toFixed(2)} is over the ${COMPLEX_DEVIANT_FRACTION} at `
        + 'which self-reference fails. The array may be perfect: it is the GENOME that is '
        + 'disturbed, and there is no undisturbed part of it left to measure the disturbed parts '
        + 'against. In the series this bound was measured on, all 70 whole-chromosome losses sat '
        + 'here and none of them was assignable to a parent.',
      limits: 'No parent is named on any event of this sample, at any confidence. The events '
        + 'themselves are real and their positions stand; which parent they came from does not.',
      measures,
    }
  }

  if (events > 0) {
    return {
      level: 'disturbed',
      headline: `${events} change${events === 1 ? '' : 's'} on `
        + `${affectedAutosomes} of ${totalAutosomes} autosomes`,
      why: 'Changes were found and the rest of the genome is undisturbed enough to measure them '
        + 'against, which is the condition every origin statistic here depends on.',
      limits: '',
      measures,
    }
  }

  return {
    level: 'intact',
    headline: 'No chromosomal change found',
    why: `Call rate ${(callRate * 100).toFixed(1)}% and no autosome carries a detected change.`,
    limits: '',
    measures,
  }
}

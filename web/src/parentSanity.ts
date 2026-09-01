/**
 * Whether the array in a parental slot is the kind of person that slot says it is.
 *
 * THE FAILURE THIS EXISTS FOR PRODUCES NO ERROR ANYWHERE ELSE. An operator opens a folder holding
 * a sperm donor and four egg donors, picks the wrong one, and marks it paternal. Every genotype in
 * that file is real, the array is high quality, no reader can object and every calculation
 * downstream is arithmetically correct. What is wrong is only the LABEL, and the label is what the
 * report is for: every event the tool then attributes comes out named for the parent it did not
 * come from. A calibration cannot catch it, because a calibration always loads the right file.
 *
 * IT IS DECIDABLE, CHEAPLY, AND THE TOOL ALREADY HAD THE MEASUREMENT. A sperm donor is male and an
 * egg donor is female, so the paternal array carries a Y and the maternal one does not. `sexing.ts`
 * has read chromosome Y since long before this module existed; nothing ever asked it about the
 * PARENT. Measured over the 12 donor arrays of GSE148488, where the roles come from the GEO titles
 * and not from this tool:
 *
 *   4 sperm donor arrays    yBearing true    Y call ratio 1.029 to 1.039
 *   8 egg donor arrays      yBearing false   Y call ratio 0.000, all eight
 *
 * Twelve of twelve, with no overlap and nothing in between. The separator is not a threshold on a
 * noisy quantity, it is the presence or absence of a chromosome.
 *
 * SILENCE IS NOT A CONFLICT. `yBearing` is null when the panel carries too few Y probes to ask,
 * and a panel that cannot answer must not be treated as an answer of "no Y": that would refuse
 * every run on every array without Y coverage. Only a decisive measurement that disagrees with the
 * declaration is a conflict.
 *
 * WHAT A CONFLICT COSTS, and it is deliberately not the whole run. The events are still real: they
 * were found by channels that never read the parental array's sex, and their positions and classes
 * do not depend on which slot the parent went into. What is withheld is every NAMED PARENT, since
 * that is the one output the swap corrupts, and it is withheld rather than inverted. Guessing that
 * the operator meant the other role would be a second guess on top of the first.
 *
 * A 46,XX MALE OR A FAILED Y REACTION LANDS HERE TOO, and stopping is still right. Both are cases
 * where the array cannot support the parental attribution the report would print, and both deserve
 * a person looking at them rather than a silent answer.
 */
import type { SexCall } from './sexing.ts'

export type ParentSexConflict =
  /** Declared paternal, no Y. The likeliest cause is an egg donor in the father's slot. */
  | 'paternal-array-has-no-y'
  /** Declared maternal, carries a Y. The likeliest cause is the sperm donor in the mother's slot. */
  | 'maternal-array-has-y'

export interface ParentSanity {
  declared: 'paternal' | 'maternal'
  /** What the array measured. Null when the panel cannot answer, which is not a conflict. */
  yBearing: boolean | null
  conflict?: ParentSexConflict
  /** One line an operator can act on. */
  headline: string
  why: string
}

/**
 * Reconcile a parental array's measured sex with the role it was declared under.
 *
 * `declared` is the slot the operator chose. `sex` is what `sexCall` read off that same array.
 */
export function reconcileParentSex(
  declared: 'paternal' | 'maternal', sex: SexCall,
): ParentSanity {
  const base = { declared, yBearing: sex.yBearing }
  const measured = `Y call ratio ${Number.isFinite(sex.callRatio) ? sex.callRatio.toFixed(3) : 'not measured'}`
    + `, Y intensity ${Number.isFinite(sex.lrrShift) ? sex.lrrShift.toFixed(3) : 'not measured'}`
    + " against this array's own autosomes"

  if (sex.yBearing === null) {
    return {
      ...base,
      headline: `The ${declared} array's sex was not checked.`,
      why: 'this array carries too few chromosome Y probes to read one way or the other, so the '
        + 'check that the file matches the slot it was put in could not run. That is a limit of '
        + 'the panel and not a finding about the sample. Nothing is withheld for it',
    }
  }

  if (declared === 'paternal' && sex.yBearing === false) {
    return {
      ...base,
      conflict: 'paternal-array-has-no-y',
      headline: 'THE PATERNAL ARRAY CARRIES NO CHROMOSOME Y. No parent is named.',
      why: `the array loaded as the paternal parent has no Y (${measured}). Across the twelve `
        + 'donor arrays this was measured on, every sperm donor read a Y call ratio of 1.03 and '
        + 'every egg donor read exactly 0, so this is the presence of a chromosome rather than a '
        + 'threshold on a noisy number. The likeliest cause is that an egg donor was put in the '
        + "paternal slot. Every event below is still real and its position stands, because the "
        + 'channels that found them never read this array\'s sex. What is withheld is the parent '
        + 'of origin: with the slots crossed, every one of those labels would name the parent the '
        + 'event did NOT come from. Check which file is in which slot and run it again',
    }
  }

  if (declared === 'maternal' && sex.yBearing === true) {
    return {
      ...base,
      conflict: 'maternal-array-has-y',
      headline: 'THE MATERNAL ARRAY CARRIES A CHROMOSOME Y. No parent is named.',
      why: `the array loaded as the maternal parent carries a Y (${measured}). The likeliest cause `
        + 'is that the sperm donor was put in the maternal slot. Every event below is still real '
        + 'and its position stands. What is withheld is the parent of origin, which with the slots '
        + 'crossed would name the wrong parent every time. Check which file is in which slot and '
        + 'run it again',
    }
  }

  return {
    ...base,
    headline: `The ${declared} array's sex matches the slot it was loaded into.`,
    why: `${measured}, which is what a ${declared === 'paternal' ? 'male' : 'female'} array reads`,
  }
}

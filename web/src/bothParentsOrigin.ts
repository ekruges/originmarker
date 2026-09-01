/**
 * Which parent's copy is missing, when BOTH parents' arrays are loaded.
 *
 * THE ASYMMETRY IS THE WHOLE REASON THIS EXISTS. The Mendelian channel answers two questions from
 * one loaded parent, and they are not equally answerable. Measured on the 14 linkage-resolvable
 * pronuclei of GSE148488, where exactly one parental copy is missing and which one is known from
 * the dissection rather than from this tool:
 *
 *   loaded parent's own copy is gone      known-parent-lost correct   0.8539 [0.6784, 1.0000]
 *   the OTHER parent's copy is gone       other-parent-lost correct   0.2890 [0.1273, 0.4506]
 *
 * Same arrays, same events, same code. The only difference is which parent was loaded. The reason
 * is structural: "the loaded parent's allele is absent" rests on alleles that parent does not
 * carry, and dropout cannot manufacture an allele. "The other parent's copy is absent" has to be
 * read from heterozygosity that is not there, and missing heterozygosity is exactly what dropout
 * produces. In the second arm the channel returned both-copies-present on 172 of 308 chromosomes
 * that carry one parental contribution.
 *
 * SO WITH BOTH PARENTS LOADED, ASK EACH ONE ONLY THE QUESTION IT CAN ANSWER. A paternal loss is a
 * known-parent-lost when the father is loaded; a maternal loss is a known-parent-lost when the
 * mother is loaded. Take the parent from whichever side returns that verdict and ignore the weak
 * direction entirely. Every loss then lands on the 0.854 arm instead of half of them landing on
 * the 0.289 one.
 *
 * The tool previously loaded only the first parent for this channel however many were supplied, so
 * roughly half of all losses were answered by the direction that is wrong about seven times in ten.
 *
 * NOTHING HERE IS A NEW MODEL. No likelihood is combined, nothing is multiplied, no independence is
 * assumed between the two calls. They are two readings of one genome and are not independent. This
 * only selects which of the two existing calls to believe, on a rule fixed in advance by a
 * measurement that predates it.
 *
 * BOTH SIDES SAYING known-parent-lost IS A CONTRADICTION AND IS REPORTED AS ONE. It asserts neither
 * parent contributed a copy here, which is nullisomy, and on real material is far more often an
 * artefact of the array than a genome missing a chromosome outright. It names no parent. Silence
 * with a stated reason beats a coin flip between two confident answers.
 */
import type { Band } from './originPosterior.ts'
import type { OneParentCall } from './oneParentOrigin.ts'

/**
 * What loading the second array actually bought, measured after the rule below was fixed.
 *
 * THE RULE WAS WRITTEN FIRST AND MEASURED SECOND, deliberately. Tuning it to these numbers would
 * be fitting it to the set it is scored on, which is the failure this whole calibration exists to
 * undo. Nothing here was adjusted after the run.
 *
 * `audit/trio-calibration.ts` against GSE148488, dropout 0.199, q at the shipped default.
 *
 * DIRECTION, on the 14 linkage-resolvable pronuclei, 308 chromosomes. Each is missing exactly one
 * parental complement and which one is known from the dissection.
 *
 *   one array, its own copy present    other-parent-lost correct   0.2890   89 right, 172 missed
 *   one array, its own copy absent     known-parent-lost correct   0.8539   263 right
 *   BOTH ARRAYS                        correct parent named        0.8539   263 named, 45 refused
 *
 * and of those 263 named, the number naming the WRONG parent is 0. The 172 misses are the case
 * this fixes: with one array loaded the channel returned both-copies-present on 172 chromosomes
 * that carry no second parental complement at all.
 *
 * SPECIFICITY, on the 76 usable trios, where every autosome carries both copies by pedigree so any
 * named parent is false. This is the cost side and it is real.
 *
 *   material         one array   both arrays
 *   trophectoderm      0.0000      0.0000     36 arrays
 *   esc-line           0.0000      0.0000     16
 *   esc-single         0.0000      0.0000      4
 *   blastomere         0.0335      0.0813     19
 *   pooled             0.0084      0.0203     76
 *
 * EVERY FALSE CALL IN BOTH ARMS IS ON A BLASTOMERE. On a trophectoderm biopsy or an ESC the second
 * array costs nothing measurable and converts the weak direction into the strong one. On a single
 * cell it roughly doubles the false-call rate, from about one intact chromosome in thirty to about
 * one in twelve, because two arrays are being asked instead of one and each carries its own chance
 * of a spurious answer.
 *
 * WHAT THAT RATE APPLIES TO, since the arm above is a worst case by construction. It scores every
 * autosome. The shipped tool runs this channel only on loci where an event was already detected by
 * another channel, so the population it is exposed to is far smaller and already selected. The
 * figure is the channel's own error rate, not the rate of wrong rows in a report.
 */
export const TWO_PARENT_ACCURACY = {
  /** Correct parent named on material where the answer is known from dissection. */
  direction: { accuracy: 0.8539, lo: 0.6784, hi: 1.0, arrays: 14, named: 263, wrongParent: 0 },
  /** The same events with one array loaded, on the half that fell to the weak direction. */
  directionOneArray: { accuracy: 0.2890, lo: 0.1273, hi: 0.4506, arrays: 14 },
  /** A parent named where both copies are present. Pooled over all 76 trios. */
  specificity: { falseCallRate: 0.0203, lo: 0.0029, hi: 0.0378, arrays: 76 },
  /**
   * WHAT CORROBORATION IS WORTH, and it is why it is reported rather than required.
   *
   * The other array pointing at the same parent by the weak direction is genuinely enriched in
   * true calls: 89 of 263 real ones carry it against 2 of 34 false ones, 0.338 against 0.059,
   * an enrichment of 5.8x. So it is real signal and a reader should see it.
   *
   * AS A GATE IT IS A BAD TRADE AND WAS MEASURED BEFORE BEING REJECTED. Requiring it would drop
   * 32 of the 34 false calls and 174 of the 263 true ones: two thirds of every real answer, to
   * remove a false-call rate already under one in fifty. The channel would name a parent on 89
   * events instead of 263. Refusing is cheaper than a wrong answer, but not at that price, and the
   * weak direction is 0.2890 accurate on its own, so its silence is barely evidence at all.
   */
  corroboration: {
    trueCalls: { agreed: 89, total: 263 },
    falseCalls: { agreed: 2, total: 34 },
  },
  /** And by material, because every one of those false calls is on a single cell. */
  specificityByMaterial: {
    trophectoderm: { falseCallRate: 0.0000, oneArray: 0.0000, arrays: 36 },
    'esc-line': { falseCallRate: 0.0000, oneArray: 0.0000, arrays: 16 },
    'esc-single': { falseCallRate: 0.0000, oneArray: 0.0000, arrays: 4 },
    blastomere: { falseCallRate: 0.0813, lo: 0.0179, hi: 0.1447, oneArray: 0.0335, arrays: 19 },
  },
} as const

export type BothParentsVerdict =
  /** One side said its own parent's copy is gone and the other did not contradict it. */
  | 'parent-lost'
  /** Neither side found its own parent's copy missing, and at least one saw both copies. */
  | 'both-present'
  /** Both sides said their own parent's copy is gone, which is nullisomy or an artefact. */
  | 'contradiction'
  /** Neither side could answer. */
  | 'refused'

export interface BothParentsCall {
  verdict: BothParentsVerdict
  /** The parent whose copy is missing, or null. Only ever set from a known-parent-lost side. */
  parent: 'paternal' | 'maternal' | null
  /** The reported confidence of the side the parent came from. Not a combined quantity. */
  posterior: number
  band: Band
  /**
   * Whether the OTHER side independently agreed, by returning other-parent-lost.
   *
   * REPORTED, NEVER REQUIRED. It carries real signal, measured at 0.338 of true calls against
   * 0.059 of false ones, but requiring it would cost two thirds of every real answer to remove a
   * false-call rate already under one in fifty. See TWO_PARENT_ACCURACY.corroboration for the
   * counts. It never creates or changes a call.
   */
  corroborated: boolean
  why: string
}

/**
 * Combine the two single-parent calls.
 *
 * `pat` must be the call made with the PATERNAL array loaded and `mat` the one with the maternal
 * array loaded, both over the same interval of the same sample.
 */
export function callBothParentsOrigin(
  pat: OneParentCall, mat: OneParentCall,
): BothParentsCall {
  const patLost = pat.verdict === 'known-parent-lost'
  const matLost = mat.verdict === 'known-parent-lost'

  if (patLost && matLost) {
    return {
      verdict: 'contradiction',
      parent: null,
      posterior: NaN,
      band: 'F',
      corroborated: false,
      why: 'both parental arrays report their OWN copy absent here, which asserts no parental copy '
        + 'is present at all. That is nullisomy, and on amplified material it is far more often the '
        + 'array than the genome. No parent is named. '
        + `Paternal side: ${pat.why}. Maternal side: ${mat.why}`,
    }
  }

  if (patLost || matLost) {
    const side = patLost ? pat : mat
    const other = patLost ? mat : pat
    const parent = patLost ? 'paternal' as const : 'maternal' as const
    // The weak direction from the other array, pointing at the same parent. Reported, never relied
    // on: other-parent-lost is correct 0.2890 of the time on this material.
    const corroborated = other.verdict === 'other-parent-lost'
    return {
      verdict: 'parent-lost',
      parent,
      posterior: side.posterior,
      band: side.band,
      corroborated,
      why: `the ${parent} array's own copy is absent here, which is read from alleles that parent `
        + 'does not carry and which dropout cannot manufacture. That direction is correct 0.8539 '
        + '[0.6784, 1.0000] of the time on known material, against 0.2890 [0.1273, 0.4506] for the '
        + 'same channel asked about the parent that is NOT loaded, so the answer is taken from this '
        + `side alone. ${side.why}. The ${patLost ? 'maternal' : 'paternal'} array `
        + (corroborated
          ? 'independently reports the other parent\'s copy missing, which agrees; that direction '
            + 'is the weak one and is noted rather than counted'
          : `reports ${other.verdict}, which is not evidence against this call: only the loaded `
            + 'parent\'s own absence is reliable in this channel'),
    }
  }

  // Neither side found its own parent's copy missing.
  if (pat.verdict === 'both-present' || mat.verdict === 'both-present') {
    const side = pat.verdict === 'both-present' ? pat : mat
    return {
      verdict: 'both-present',
      parent: null,
      posterior: side.posterior,
      band: side.band,
      corroborated: pat.verdict === 'both-present' && mat.verdict === 'both-present',
      why: 'neither parental array reports its own copy absent, and at least one sees both copies '
        + `present. Paternal side: ${pat.verdict}. Maternal side: ${mat.verdict}`,
    }
  }

  return {
    verdict: 'refused',
    parent: null,
    posterior: NaN,
    band: 'F',
    corroborated: false,
    why: 'neither parental array could answer here. '
      + `Paternal side: ${pat.why}. Maternal side: ${mat.why}`,
  }
}

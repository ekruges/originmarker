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
import { bandObligateHet, bandOf, type Band } from './originPosterior.ts'

/**
 * Frequency of the allele the known parent lacks, when it cannot be estimated per marker.
 *
 * MEASURED MODEL-FREE, which is the only way it can be measured without assuming the thing it is
 * used to test. At a marker where the loaded parent is homozygous, the OTHER parent's own genotype
 * gives the transmission probability of the lacked allele exactly: 0 if that parent lacks it too,
 * 0.5 if heterozygous, 1 if homozygous for it. No Hardy-Weinberg, no dropout model, no population
 * reference. Over 28 cross-person donor orientations on GSE148488, 607,000 to 658,000 markers each:
 * min 0.1080, median 0.1217, mean 0.1197, max 0.1303. Technical replicates of one donor agree to
 * 0.0033, so the spread above is between people, not measurement.
 *
 * THE CLOSED FORM DOES NOT REPRODUCE IT HERE, and the reason is worth recording. q_cond =
 * (H/2)/(1-H) evaluated at each parent's own heterozygosity gives 0.0822 to 0.1053, about 15
 * percent below the measured value, because it uses the LOADED parent's allele spectrum while q is
 * a property of what the OTHER parent transmits. In a cross-population cross those differ, and
 * every trio in this series is one. Prefer the model-free figure; use the closed form only where
 * both parents are known to be drawn from the same population.
 *
 * WHY IT MATTERS. The previous value, 0.30, was about 2.5 times this. The error is one-sided:
 * too high manufactures "the other parent's copy is lost" on chromosomes that are intact.
 */
export const DEFAULT_Q = 0.1217

/** Genotyping error: a homozygote read as the other homozygote, or a spurious heterozygote. */
export const GENOTYPE_ERROR = 0.02

/**
 * Drop-in, a heterozygous call at a truly homozygous marker. Measured on this platform across 113
 * same-genome array pairs, medians 0.0390 to 0.0525 by experiment. It matters here because a
 * spurious AB argues for both copies being present and so masks a real deletion.
 */
export const DROP_IN = 0.0435

/** Top of the same measured drop-in range, which is what a ceiling must be sized against. */
export const DROPIN_CEILING = 0.0525

/** Posterior a hypothesis must reach before the region is called. */
export const CALL_POSTERIOR = 0.95

/**
 * The error this channel's validation cannot rule out, which is what caps a reported confidence.
 *
 * THE LIKELIHOOD SATURATES BECAUSE IT IS CONDITIONED ON THE MODEL BEING RIGHT, and over hundreds of
 * near-independent Mendelian markers it is right or wrong decisively. Measured on the shipped code:
 * 1.0000 at 0, 4, 8, 12, 16, 17, 18, 20, 24, 30 and 60 exclusive markers out of 400, including on
 * both sides of the point where the verdict flips. That is the likelihood answering the question it
 * was asked, not a defect in it.
 *
 * It is the wrong quantity to REPORT, because the model can be wrong in ways no marker count
 * touches: a contaminated reaction, a mislabelled file, an error rate mis-specified for a chemistry
 * not in the calibration set. What a reader needs is the chance the CALL is right, and that is
 * bounded by what the validation supports rather than by what the likelihood computes.
 *
 * Both anchors are this project's own. The trio validation is 12 of 12 with the mother hidden. The
 * methods review computed the cluster-robust form of exactly this kind of claim: 251 of 251 correct
 * from 13 INDEPENDENT UNITS, a one-sided zero-event error bound of 0.206 rather than 0.012, and it
 * says in terms that quoting the uncorrected figure "overstates by a factor of 17". So a reported
 * confidence above 1 - 0.206 claims validation nobody has.
 */
/**
 * What this channel is worth on an OBVIOUS whole-chromosome event, measured by construction.
 *
 * A blastomere carries both parental genomes, which is what makes this channel work there without
 * any detection floor: at a marker where the loaded parent is homozygous, losing that parent's copy
 * leaves an allele the parent does not have, and dropout removes alleles without inventing one. It
 * asks whether an allele is present, not whether a mean has moved.
 *
 * MEASURED ON REAL TRIOS WITH A WRONG-PARENT CONTROL, which is what the previous figure lacked.
 *
 * The number that used to sit here, 0.920, came from a harness that built the child's post-loss
 * genotype out of the LOADED parent's own alleles. The signal it scored therefore existed whoever
 * was loaded: 44 of 44 with the true father, 44 of 44 with the true mother, and 44 of 44 with
 * arrays that are certainly not parents. It was not a measurement of parent of origin.
 *
 * These come from `audit/trio-calibration.ts` against GSE148488, where the pedigree is external and
 * genetically verified, and every arm is reported twice: once with the real parent and once with an
 * unrelated adult. Dropout 0.199, q at the shipped default.
 *
 * SPECIFICITY, intact chromosomes of real children, where both-copies-present is correct by
 * pedigree. A false call is any loss verdict. OVER ALL 76 USABLE TRIOS, not the 24-trio slice the
 * first version of this figure ran on.
 *   true parent loaded    0.0084 [0.0006, 0.0161]   76 arrays, 1,672 chromosome calls
 *   unrelated adult       0.0598 [0.0230, 0.0966]   76 arrays
 * The arms separate and do not overlap, which is what makes this a measurement of parentage rather
 * than of the construction.
 *
 * AND IT IS NOT ONE NUMBER. A blastomere is ONE cell and a trophectoderm biopsy is five to ten, so
 * they do not drop out at the same rate and cannot share a false-call rate. The pooled figure above
 * is dominated by whichever material is most numerous, which is a property of this dataset and not
 * of anyone's samples. Split, true parent loaded against unrelated adult:
 *
 *   trophectoderm    36   0.0000 [0.0000, 0.0000]     0.0025 [0.0000, 0.0060]
 *   esc-line         16   0.0000 [0.0000, 0.0000]     0.0057 [0.0000, 0.0133]
 *   esc-single        4   0.0000 [0.0000, 0.0000]     0.3523 [0.0000, 0.7953]
 *   blastomere       19   0.0335 [0.0047, 0.0623]     0.1555 [0.0666, 0.2444]
 *   zygote            1   0.0000                      0.0000
 *
 * THE BLASTOMERE ROW IS THE ONE TO READ BEFORE USING THIS ON CLEAVAGE-STAGE MATERIAL. About one
 * chromosome in thirty comes back with a false loss on a genome that has both parental copies,
 * against none in 792 trophectoderm chromosome calls. The two arms still separate by a factor of
 * 4.6, so the call is reading parentage rather than noise, but the error rate is a different order
 * of magnitude from a biopsy of several cells and a single blastomere call should be read as one
 * observation rather than as a result.
 *
 * DIRECTION, on the 15 pronuclei, where exactly one parental copy is missing and which one is known
 * from the dissection. Parents resolved by linkage, 14 of 15 resolved.
 *   loaded parent is NOT this genome's   known-parent-lost correct   0.8539 [0.6784, 1.0000]
 *   loaded parent IS this genome's       other-parent-lost correct   0.2890 [0.1273, 0.4506]
 *
 * THE ASYMMETRY IS THE RESULT AND IT IS STRUCTURAL. "The loaded parent's copy is absent" is read
 * from alleles the parent does not carry, which dropout cannot manufacture, and it holds at 0.854.
 * "The OTHER parent's copy is absent" has to be inferred from heterozygosity that is not there, and
 * missing heterozygosity is exactly what dropout produces, so it holds at 0.289 and returned
 * both-copies-present on 173 chromosomes that carry one parental contribution.
 *
 * So `known-parent-lost` is reportable and `other-parent-lost` is not, on this material.
 */
export const OBVIOUS_EVENT_ACCURACY = {
  /** Intact chromosomes, real children, true parent loaded. A false call is any loss verdict. */
  specificity: { falseCallRate: 0.0084, lo: 0.0006, hi: 0.0161, arrays: 76 },
  /** The same, with an unrelated adult loaded. Must be worse, or nothing was measured. */
  specificityWrongParent: { falseCallRate: 0.0598, lo: 0.0230, hi: 0.0966, arrays: 76 },
  /**
   * The same specificity, BY MATERIAL, because one cell and ten cells do not share a dropout rate.
   * A caller that knows what it is holding should quote the row for that material and not the
   * pooled figure. `wrong` is the unrelated-adult arm on the same arrays.
   */
  specificityByMaterial: {
    trophectoderm: { falseCallRate: 0.0000, lo: 0.0000, hi: 0.0000, wrong: 0.0025, arrays: 36 },
    'esc-line': { falseCallRate: 0.0000, lo: 0.0000, hi: 0.0000, wrong: 0.0057, arrays: 16 },
    'esc-single': { falseCallRate: 0.0000, lo: 0.0000, hi: 0.0000, wrong: 0.3523, arrays: 4 },
    blastomere: { falseCallRate: 0.0335, lo: 0.0047, hi: 0.0623, wrong: 0.1555, arrays: 19 },
  },
  /** Pronuclei: the loaded parent's own copy is genuinely gone. Read from exclusive alleles. */
  knownParentLost: { accuracy: 0.8539, lo: 0.6784, hi: 1.0, arrays: 14 },
  /** Pronuclei: the OTHER parent's copy is gone. Inferred from absent heterozygosity. */
  otherParentLost: { accuracy: 0.2890, lo: 0.1273, hi: 0.4506, arrays: 14 },
} as const

export const VALIDATION_UNITS = 13
export const SYSTEMATIC_ERROR_BOUND = 0.206

/**
 * Log-odds at which half the distance from the floor to the cap is reached.
 *
 * Sets how quickly evidence buys confidence. A log-odds of 12 is a likelihood ratio of about
 * 160,000, which is where a Mendelian call over a few hundred markers stops being marginal.
 */
export const HALF_EVIDENCE = 12

/**
 * Turn the saturated likelihood into a confidence a reader can act on.
 *
 * TWO CHANGES, AND THE SECOND IS WHY THE NUMBER MOVES AT ALL. The cap comes first: nothing is
 * reported above what the truth set supports. Then the ORDERING is taken from the likelihood's
 * LOG-ODDS rather than its probability, because that is where the evidence varies. A raw posterior
 * runs 0.999869 to 1.0000 over the whole usable range, one part in ten thousand and unreadable; the
 * same evidence in log-odds runs from about 9 to over a thousand.
 *
 * Monotone, so it cannot reorder two calls, and a reparameterisation rather than new evidence: the
 * VERDICT is untouched. Only the number beside it changes, from one that claimed certainty the
 * validation cannot support to one bounded by it.
 *
 * IT TAKES THE LOG MARGIN, NOT THE PROBABILITY, and that distinction is the whole thing. Once the
 * log-likelihoods are exponentiated and normalised, the best hypothesis reaches exactly 1.0 in
 * floating point and every difference between a marginal call and an overwhelming one is gone: a
 * first attempt read the probability, clamped it away from 1, and ended up reporting the CLAMP
 * rather than the evidence. The margin between the best hypothesis and the next best is the same
 * information before it is destroyed, and it runs from about zero to over a thousand.
 */
export function reportedConfidence(
  logMargin: number,
  opts: { bound?: number; halfEvidence?: number; floor?: number } = {},
): number {
  if (!Number.isFinite(logMargin)) return NaN
  const bound = opts.bound ?? SYSTEMATIC_ERROR_BOUND
  const half = opts.halfEvidence ?? HALF_EVIDENCE
  /**
   * CHANCE, WHICH DEPENDS ON HOW MANY HYPOTHESES ARE ON THE TABLE.
   *
   * One third is right for THIS channel, which chooses among three: both copies present, the
   * loaded parent's copy absent, the other parent's absent. It is wrong for a channel asking a
   * two-way question, where an uninformative answer is a coin flip at one half and starting the
   * scale below that reports ignorance as worse than chance, dragging every number above it down
   * with it. A caller with a different hypothesis count passes its own.
   */
  const floor = opts.floor ?? 1 / 3
  const cap = 1 - bound
  if (!(logMargin > 0)) return floor
  return floor + (cap - floor) * (logMargin / (logMargin + half))
}

/** Fewest informative markers before a region is scored. */
export const MIN_MARKERS = 50

/**
 * Heterozygosity above which the region's genotypes are not measuring it, and no verdict is issued.
 *
 * At a marker where the known parent is HOMOZYGOUS, a biparental child is heterozygous exactly
 * when the other parent transmitted the allele this one lacks. That is the panel's own allele
 * frequency, so the most a real region can show is `q` plus drop-in at the top of the range
 * measured here, and dropout only lowers it.
 *
 * DERIVED FROM DEFAULT_Q RATHER THAN WRITTEN DOWN, because it is a function of it and the two had
 * already drifted apart. This was 0.40, sized against a q of 0.30. At the measured q the real
 * ceiling is 0.1217 + 0.0525 = 0.1742, so the old value admitted regions at 2.3 times the
 * heterozygosity it exists to reject, and correcting q alone would have left that gap open.
 *
 * WHY IT IS NEEDED, from a real array. A blastomere carrying a call-rate collapse on chromosome 1,
 * 69% no-call, returned 59.7% heterozygous across the markers where the father is homozygous. No
 * genome reads that way. The likelihood took the heterozygote count at face value, since a
 * heterozygote is near-impossible under either deletion hypothesis, and returned both-copies-
 * present at posterior 1.000 on genotypes that carried no information at all. A confident verdict
 * from meaningless input is worse than a refusal, because only the refusal is visible.
 */
export const MAX_REGION_HET = DEFAULT_Q + DROPIN_CEILING

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
    /** Error the validation cannot rule out, which caps the reported confidence. */
    systematicErrorBound?: number,
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
      ...base, verdict: 'refused', posterior: NaN, band: 'F' as Band,
      why: `${n} informative markers is under the ${minMarkers} this needs; a marker only counts `
        + 'where the loaded parent is homozygous',
    }
  }

  // The genotypes must be measuring the region before their likelihood means anything. This is
  // the region-level form of the array-level ceiling in stage.ts, and it exists for the same
  // reason: an impossible rate is evidence about the reaction, not about the genome.
  if (het / n > maxRegionHet) {
    return {
      ...base, verdict: 'refused', posterior: NaN, band: 'F' as Band,
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

  // The margin over the next-best hypothesis, in log space, which is the evidence before
  // normalising destroys it.
  const runnerUp = Math.max(...logs.filter((_, i) => i !== best))
  const reported = reportedConfidence(logs[best] - runnerUp, { bound: opts.systematicErrorBound })

  if (post[best] < callPosterior) {
    return {
      // This refusal DOES carry a number, so it is graded from that number rather than pinned to
      // the weakest measured band. A posterior of 0.51 and an absent one are not the same row.
      ...base, verdict: 'refused', posterior: post[best], band: bandOf(reported) as Band,
      why: `best hypothesis reaches ${post[best].toFixed(3)}, under the ${callPosterior} needed`,
    }
  }
  const rate = n ? exclusive / n : 0
  // BAND F NAMES NOBODY IN THIS CHANNEL EITHER. The dosage channel has withheld the parent in
  // band F since a gynogenetic sample was told its paternal copy was the one lost; this channel
  // never had the guard, and the two share the band ladder, so the ladder's own docstring read as
  // if one guard covered both.
  //
  // The gap is real rather than theoretical. `post[best]` is decisive under the model while
  // `reported` is what survives SYSTEMATIC_ERROR_BOUND, and on a thin panel the second falls to
  // chance while the first is still over `callPosterior`. Run through the shipped composition,
  // 100 markers with 2 exclusive and 200 with 5 give `other-parent-lost` at reported 0.4824 and
  // 0.5291, band F, each printing a parent's name. audit/requirements.check.ts prints those two
  // rows unprompted. A parent named at 0.48 is a coin flip wearing a name, which is the one thing
  // this tool must not produce.
  //
  // No row on real material moves: over the 76 usable trios of GSE148488, in both the true-parent
  // and the wrong-parent arm, 0 band-F rows named a parent. This closes the hole rather than
  // changing an answer anyone has seen.
  const band = bandObligateHet(reported)
  if (band === 'F') {
    return {
      ...base,
      verdict: 'refused',
      posterior: reported,
      band,
      why: `${names[best]} is the best hypothesis, at posterior ${post[best].toFixed(4)} over `
        + `${n} markers where the loaded parent is homozygous. Reported confidence after the `
        + `systematic bound is ${reported.toFixed(4)}, which is band F. No parent is named there: `
        + 'the evidence is reported and the attribution is withheld',
    }
  }
  return {
    ...base,
    verdict: names[best],
    // Reported rather than raw: the likelihood is decisive under its own model, and the model is
    // not the only thing that can be wrong. See SYSTEMATIC_ERROR_BOUND.
    posterior: reported,
    band,
    why: `${names[best]} at posterior ${post[best].toFixed(4)} over ${n} markers where the loaded `
      + `parent is homozygous. ${exclusive} carry the allele that parent does not have `
      + `(${(rate * 100).toFixed(1)}%), which only its absence explains, and ${het} are `
      + `heterozygous, which only both copies being present explains`,
  }
}

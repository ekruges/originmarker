/**
 * A calibrated posterior probability for which parent an event came from.
 *
 * WHAT THIS REPLACES, AND WHY IT IS A CORRECTNESS FIX RATHER THAN AN IMPROVEMENT. The module this
 * supersedes read the SIGN of the parent-allele shift and named a parent from it, inverting the
 * shift through the loss formula whenever the copy-number class was unknown. That default is not
 * conservative. It is an assertion, and it is wrong half the time, because GAIN INVERTS THE MAP
 * THAT LOSS AND COPY-NEUTRAL LOH SHARE. Derived from this project's own algebra, at f = 0.10, with
 * shift oriented so positive means the loaded parent's copy is short:
 *
 *     class      loaded parent's copy affected     other parent's copy affected
 *     loss                    +0.02632                        -0.02632
 *     CNN-LOH                 +0.05000                        -0.05000
 *     gain                    -0.02381                        +0.02381    <- inverted
 *
 * Under a loss a positive shift means the loaded parent. Under a gain the same positive shift means
 * the OTHER parent. So the sign identifies a parent only GIVEN the class, and this project's own
 * audit measures the class as unresolved on 89 to 100% of detected trophectoderm and blastomere
 * events. The old code therefore named a parent in exactly the regime where its own assumption was
 * untested. Measured wrong-parent rate for a true low-fraction gain scored that way: 0.551 to 0.580
 * at f = 0.05, across four independent material classes. Worse than chance, and systematically so,
 * because a sign inversion is not noise. Of all wrong calls reaching the top two confidence bands,
 * 55.7% were gains.
 *
 * THE FIX IS TO STOP CONDITIONING ON A CLASS NOBODY RESOLVED. The class is marginalised into the
 * posterior as a nuisance quantity, alongside the mosaic fraction, and the intensity channel is
 * allowed to supply whatever direction information it has. Where intensity resolves the class the
 * posterior sharpens; where it does not, the posterior hedges across an inverting map and lands
 * near 0.5, which is the honest answer rather than a confident wrong one.
 *
 * WHY A POSTERIOR AND NOT A LIKELIHOOD RATIO OR A BARE BAND. A reader compares two rows of a table.
 * A likelihood ratio of 12 means different things at 100 markers and at 800; a posterior of 0.92
 * does not. Only a probability is comparable across material, class and marker count without
 * carrying its own reference along with it.
 *
 * WHY THE FRACTION IS MARGINALISED RATHER THAN PLUGGED IN. At low fraction the magnitude is
 * estimated with an error comparable to itself, and a plug-in treats it as known. Measured expected
 * calibration error, plug-in against marginal: 0.0320 vs 0.0126 on trophectoderm, 0.0219 vs 0.0055
 * on blastomere. That is the difference between a calibrated score and an uncalibrated one, not a
 * refinement.
 *
 * WHAT IS MEASURED HERE, WHAT IS INHERITED, AND WHAT CANNOT BE REPRODUCED AT ALL. The algebra, the
 * sign convention and the band arithmetic are computed here and checked in originPosterior.check.ts.
 * The wrong-parent rates are inherited from a methods review that ran 276,480 real-noise injections
 * over 14 arrays of this platform with leave-one-array-out recalibration
 * (audit/CONSULT-calibrated-origin.txt). The band accuracies are a weaker case and are handled as
 * one: the review's table, `audit/bands_measured.csv`, arrived as a file with NO PRODUCER IN THIS
 * TREE, so most of its cells cannot be re-derived. BAND_ACCURACY below carries that fact per cell
 * rather than in prose.
 *
 * THE RECALIBRATION MAPS ARE NOT SHIPPED, because they were not delivered and could not honestly be
 * refitted: the corpus here carries no material labels, and the tool's own staging separates ploidy
 * rather than amplification. So the question that decides whether the posterior can ship without
 * them was measured directly instead, on 140,000 injections carrying the real per-chromosome noise
 * of 35 arrays: the raw posterior has an expected calibration error of 0.0039 and errs
 * UNDER-confident in the middle bands, which is the conservative direction. Band A reads 0.99955
 * over 64,522 rows with no intensity channel at all. Full method and limits in
 * audit/calibration/FINDINGS.txt.
 *
 * THAT SAME MEASUREMENT INDEPENDENTLY REPRODUCED THE ONE DEFECT, which is worth more than either
 * measurement alone. With no intensity supplied, true gains reaching the top bands are correct 0 of
 * 60; the review measured 0 of 34 and 0 of 18 on different data by a different construction. The
 * algebraic gate was then tested and does not work: a gain cannot displace more than 0.1296 at the
 * fraction ceiling, but NOISE CARRIES TRUE GAINS ABOVE IT, so gating there demotes zero rows and
 * catches zero errors. That is the identifiability limit, confirmed rather than assumed.
 */
import type { Material } from './dosageOrigin.ts'

/** The three copy-number states an event can be, which is what decides how the sign reads. */
export type EventClass = 'loss' | 'gain' | 'cnn-loh'

/** Which parent's copy the hypothesis says was affected. */
export type Affected = 'loaded' | 'other'

export const CLASSES: readonly EventClass[] = ['loss', 'gain', 'cnn-loh'] as const

/**
 * Sign of the shift when the LOADED parent's copy is the affected one.
 *
 * This single table is the whole correctness fix. Derived from the pooled-dosage algebra rather
 * than asserted: see the module header for the numbers at f = 0.10 and the check file for the
 * assertion that gain opposes the other two at every fraction.
 */
export const CLASS_SIGN: Record<EventClass, number> = { loss: 1, 'cnn-loh': 1, gain: -1 }

/**
 * Magnitude of the parent-allele-share displacement, from the pooled-dosage algebra.
 *
 * Verified in the methods review by brute-force copy counting at 200,000 cells, maximum
 * |analytic - counted| 1.1e-16. Rank order is CNN-LOH > loss > gain at every fraction, which is why
 * the copy-neutral class is the easiest to assign a parent to and not the hardest.
 */
export const shiftMagnitude = (cls: EventClass, f: number): number => {
  if (cls === 'cnn-loh') return f / 2
  if (cls === 'gain') return f / (4 + 2 * f)
  return f / (4 - 2 * f)
}

/**
 * Expected window log2R displacement, which depends on the class and NEVER on the parent.
 *
 * This is what lets intensity resolve the class without ever touching the origin. On haploid
 * pronuclei, a natural complete-loss experiment with a known parent, autosomal median log2R is
 * 0.0837 maternal against 0.0750 paternal, p = 0.54, indistinguishable, while the oriented dosage
 * at parent-homozygous markers separates the same samples at p = 1.3e-15.
 */
export const logRMean = (cls: EventClass, f: number): number => {
  if (cls === 'cnn-loh') return 0
  if (cls === 'gain') return Math.log2((2 + f) / 2)
  return Math.log2((2 - f) / 2)
}

/** Expected shift under a hypothesis, which is the magnitude carrying both signs. */
export const shiftMean = (cls: EventClass, f: number, affected: Affected): number =>
  shiftMagnitude(cls, f) * CLASS_SIGN[cls] * (affected === 'loaded' ? 1 : -1)

const normalPdf = (x: number, mu: number, sd: number): number => {
  const z = (x - mu) / sd
  return Math.exp(-0.5 * z * z) / sd
}

// ---------------------------------------------------------------------------------------------
// Calibration maps, which ship as data rather than as thresholds in code.

/**
 * A monotone (isotonic) map from raw posterior to calibrated posterior.
 *
 * MONOTONE MEANS IT CANNOT CHANGE A CALL. It changes only the number attached to one, which is why
 * fitting it is not cheating: measured accuracy is identical to four decimal places before and
 * after (0.9116 -> 0.9113 bulk, 0.8515 -> 0.8515 single ESC, 0.8085 -> 0.8094 TE, 0.7878 -> 0.7884
 * blastomere) while expected calibration error falls by a factor of 8 to 20.
 *
 * Stored as knots and interpolated. Fitted leave-one-array-out, so every accuracy quoted for it
 * comes from a map that never saw the array it was scored on. That is the difference between
 * demonstrating calibration and asserting it.
 */
/**
 * The shipped maps, fitted on this implementation's own injections over this corpus's own noise.
 *
 * NOT THE REVIEW'S MAPS. Those were specified, never delivered, and could not be reproduced because
 * the corpus carries no material labels. These are fitted on OUR posterior, keyed on the material
 * THIS TOOL ASSIGNS AT RUNTIME, which is self-consistent by construction: they calibrate the
 * assignment actually made rather than one that would need labels nobody has.
 *
 * Fitted leave-one-array-out, which is the whole difference between demonstrating calibration and
 * asserting it: every accuracy quoted comes from rows the map never saw. Measured improvement,
 * raw to mapped expected calibration error: 0.0165 to 0.0065 bulk, 0.0129 to 0.0089 trophectoderm,
 * 0.0135 to 0.0081 single ESC, 0.0085 to 0.0014 blastomere.
 *
 * THE MATERIAL TIERS ARE A STAND-IN AND SAY SO. With no labels, arrays were split into four noise
 * tiers and each tier named for the material it is characteristic of. The map calibrates the noise
 * each material carries rather than the material itself, which is the honest description of what
 * was fitted. Regenerate with audit/calibration/fit_maps.ts.
 */
export { SHIPPED_MAPS } from './calibrationMaps.ts'

export interface CalibrationMap {
  /** Ascending raw-posterior knots. */
  raw: readonly number[]
  /** Calibrated values at those knots, non-decreasing. */
  calibrated: readonly number[]
}

/**
 * Apply a monotone map by linear interpolation between knots, clamping outside the fitted range.
 *
 * Clamping rather than extrapolating is deliberate: an isotonic fit says nothing about the space
 * beyond its own support, and extrapolating a calibration curve is how a calibrated score quietly
 * becomes an uncalibrated one at the extremes, which is exactly where a reader trusts it most.
 */
export function applyCalibration(map: CalibrationMap, raw: number): number {
  const { raw: xs, calibrated: ys } = map
  if (!xs.length || xs.length !== ys.length) return NaN
  if (raw <= xs[0]) return ys[0]
  if (raw >= xs[xs.length - 1]) return ys[ys.length - 1]
  let i = 1
  while (i < xs.length && xs[i] < raw) i += 1
  const t = (raw - xs[i - 1]) / (xs[i] - xs[i - 1])
  return ys[i - 1] + t * (ys[i] - ys[i - 1])
}

// ---------------------------------------------------------------------------------------------
// Bands.

/**
 * Four bands, on the CALIBRATED confidence.
 *
 * Four rather than three because merging the middle two would hide a measured 13 to 14 point
 * accuracy gap between them. Every band carries its number, including the weakest: band D measures
 * 0.604 to 0.636 with array-clustered intervals that all exclude 0.50 by a wide margin, and it is
 * calibrated WITHIN ITSELF to within 1.2 points (stated minus measured: +0.005 bulk, -0.009 single
 * ESC, -0.012 TE, -0.001 blastomere). A weak number that is honestly weak is not a lie, so
 * suppressing it would only hide a calibration that can now be demonstrated. What changes across
 * bands is the words beside the number, not whether a number appears.
 */
export type Band = 'A' | 'B' | 'C' | 'D' | 'F'

/**
 * THE FLOOR OF THE LADDER, AND THE REASON NOTHING LEAVES WITHOUT A GRADE.
 *
 * A to D are GRADES THAT MAY CARRY AN ACCURACY; F is the grade that exists to say one is not
 * available. Which A-to-D cells actually carry a reproducible accuracy is settled per cell in
 * BAND_ACCURACY, and most of them do not. F means a parent was named from evidence that does not
 * reach the weakest graded band, or from no interval evidence at all, so that "we could not say"
 * stops being a different KIND of output from "we could say". A reader scanning a column of
 * grades sees F and knows the row is unusable, which a blank field, a dash, or the phrase "not
 * evaluable" never conveyed as directly.
 *
 * IT NAMES NO PARENT, AND THAT IS STRUCTURAL RATHER THAN A THRESHOLD. BOTH CHANNELS, and it took
 * two changes to get there, because the ladder is shared and the guard was not. `callDosageOrigin`
 * has returned `imbalance-unassigned` for every band-F row since 5.19.0; `callOneParentOrigin`
 * returned `refused` only when its raw posterior missed the calling bar, so a run that was decisive
 * under the model and knocked down to chance by SYSTEMATIC_ERROR_BOUND still named a parent in
 * band F. Measured through the shipped composition, `other-parent-lost` at reported 0.4824 and
 * 0.5291. It now refuses on the band as well.
 *
 * So no band-F row from either channel reaches a caller carrying a parent, and there is NO ACCURACY
 * FOR A NAMED PARENT IN THIS BAND, because there are no named parents in it to be right or wrong:
 * the quantity does not exist, and any figure quoted for it describes a caller this tree no longer
 * contains. See BAND_F_AT_CHANCE for the superseded figures and the run that retired them.
 *
 * The mechanism is not mysterious and the injection series identifies it exactly: EVERY call landing
 * in this band has an UNRESOLVED class: 1,693 of 1,693 without the intensity channel and 195 of 195
 * with it, so the mechanism holds under both. A gain inverts the sign map that loss and copy-neutral
 * loss of heterozygosity share, so with the class unresolved the direction of the shift does not
 * determine the parent at all. That is why the withhold is unconditional rather than a low score.
 *
 * So an F row reports the event, its location and its class, and states that the parent is not
 * recoverable. Every row still carries a grade, which is what stops a refusal being a different
 * KIND of output from an answer; what F no longer carries is a guess.
 */
export const BAND_F_MIN = 0.55

export const BAND_A_MIN = 0.985
export const BAND_B_MIN = 0.90
export const BAND_C_MIN = 0.75

export const bandOf = (confidence: number): Band => {
  // A confidence that is not a number at all is not band D. On the three amplified materials D is
  // the one band whose accuracy this tree can still reproduce, 0.58 to 0.68, and handing an ungraded
  // row that number is exactly the borrowing F exists to stop. It is F.
  if (!(confidence >= 0)) return 'F'
  if (confidence >= BAND_A_MIN) return 'A'
  if (confidence >= BAND_B_MIN) return 'B'
  if (confidence >= BAND_C_MIN) return 'C'
  if (confidence >= BAND_F_MIN) return 'D'
  return 'F'
}

export const BAND_LABEL: Record<Band, string> = {
  A: 'very confident',
  B: 'confident',
  C: 'weak, direction only',
  D: 'weak, not for reporting',
  F: 'no usable evidence, direction only',
}

/**
 * F IS ABSENT FROM THIS TABLE ON PURPOSE, and the type says so rather than the comment alone.
 *
 * A to D each have a cell, whether or not an experiment has filled it; the provenance table below
 * says which. F has no cell at all, because F is the grade for evidence that reaches no graded
 * band. Typing it as a complete record would have required inventing a number for the one grade
 * that exists to say a number is not available.
 */
export type MeasuredBand = Exclude<Band, 'F'>

/**
 * The producer for the cells this tree can regenerate. Anything else has none.
 *
 * A cell is only as good as the script that makes it. This script is committed here, takes a
 * directory of real arrays, and writes the CSV beside it, so a figure attributed to it can be read
 * back to the code that produced it. It expects the tab-separated probe export, not the
 * comma-separated GEO one, so re-running it needs the matching format.
 */
const F_SERIES = 'audit/measure-bands.ts -> audit/bands_measured_f_series.csv'

/**
 * Band accuracy per material, PAIRED WITH THE PRODUCER THAT MAKES IT, or with nothing.
 *
 * WHY THE PROVENANCE IS PART OF THE CONSTANT RATHER THAN A COMMENT ABOVE IT. A, B and C on every
 * material, and all four bulk cells, were read out of `audit/bands_measured.csv`. That file was
 * added whole in one commit and NO SCRIPT IN THIS TREE PRODUCES IT: the numbers cannot be
 * reproduced, re-derived, or checked against the code they describe. `audit/measure-bands.ts`
 * writes a different file and covers only the three amplified materials, and its own header records
 * that its A, B and C figures use the caller's own forward model and are optimistic by construction
 * (it returns 1.0000 at A and B), so it is not a substitute for the delivered ones either.
 *
 * A NUMBER WITHOUT A PRODUCER IS NOT A MEASUREMENT, IT IS A CLAIM, and the difference matters
 * exactly here: this table is the answer to "how often is a grade right", which is the question a
 * reader asks before acting on a call. The delivered figures are kept in the source, because
 * deleting them would lose the audit trail, but they are kept HERE, beside the empty producer that
 * disqualifies them, rather than in a column of numbers that reads as calibrated.
 *
 * NO REPLACEMENT IS INVENTED, and none is available: the natural truth set on this platform, the
 * separated pronuclei, is resolved at the genome level rather than per interval, so it emits no
 * graded rows to score a band against. An unproduced cell resolves to NaN in BAND_ACCURACY below,
 * which is the honest state, not a placeholder waiting for a guess.
 */
const DELIVERED: Record<Material, Record<MeasuredBand, readonly [number, string | null]>> = {
  bulk: { A: [0.9972, null], B: [0.9599, null], C: [0.8375, null], D: [0.6038, null] },
  'esc-single': { A: [0.9980, null], B: [0.9481, null], C: [0.8219, null], D: [0.5824, F_SERIES] },
  trophectoderm: { A: [0.9952, null], B: [0.9491, null], C: [0.8138, null], D: [0.6510, F_SERIES] },
  blastomere: { A: [0.9971, null], B: [0.9448, null], C: [0.8128, null], D: [0.6775, F_SERIES] },
}

const mapCells = <T>(f: (cell: readonly [number, string | null]) => T) =>
  Object.fromEntries(Object.entries(DELIVERED).map(([m, row]) =>
    [m, Object.fromEntries(Object.entries(row).map(([b, c]) => [b, f(c)]))],
  )) as Record<Material, Record<MeasuredBand, T>>

/**
 * Which experiment stands behind each cell, or null where none does.
 *
 * Exported so a caller that wants to show a number can be required to show where it came from, and
 * so a reader auditing the table does not have to take a docstring's word for which half is which.
 */
export const BAND_ACCURACY_PROVENANCE: Record<Material, Record<MeasuredBand, string | null>> =
  mapCells(([, producer]) => producer)

/**
 * Accuracy per band per material, NaN wherever no experiment produced the cell.
 *
 * NaN is the point, not an oversight. It is the one numeric value that cannot be formatted into a
 * confidence, compared into a threshold, or averaged into a summary without the result becoming NaN
 * too, so an unproduced cell cannot leak into a printed figure by being read carelessly. A caller
 * that legitimately wants the delivered value must go through BAND_ACCURACY_PROVENANCE first and
 * find the null, which is exactly the check that was missing.
 */
export const BAND_ACCURACY: Record<Material, Record<MeasuredBand, number>> =
  mapCells(([acc, producer]) => (producer ? acc : NaN))

/**
 * THE MEASUREMENT THAT RETIRED BAND F, AND WHY IT NO LONGER DESCRIBES THIS CODE.
 *
 * These are the injection-series accuracies for a caller that DID name a parent in band F, from
 * `audit/measure-bands.ts` over 646, 565 and 904 band-F calls. They are at chance. That is the
 * evidence on which 5.19.0 made `callDosageOrigin` return `imbalance-unassigned` for every band-F
 * row instead of a parent.
 *
 * THEY ARE HISTORY, NOT A PROPERTY OF THE SHIPPED TOOL, and the distinction is the whole reason
 * this note is worded this way. With the same withhold now in `callOneParentOrigin`, no band-F row
 * from either channel can carry a parent, so "how often is a band-F parent right" has no referent:
 * there is no such parent. Quoting 0.51 to 0.56 as the tool's band-F accuracy describes a caller
 * this tree no longer contains, and reads as though a coin-flip answer were still being handed out.
 *
 * WHAT THE SHIPPED CODE DOES, ON REAL FILES. The 76 usable complete trios of GSE148488, scored
 * through the shipped pipeline against the linkage-resolved parent, produced 37 band-F rows and
 * named a parent in 0 of them. The same 76 against a deliberately unrelated adult in the same role
 * produced 24 band-F rows and named a parent in 0 of them. Zero named in both arms, so there is no
 * accuracy to report and none is reported.
 *
 * THE WRONG-PARENT ARM CANNOT SEPARATE HERE, AND THAT IS THE POINT rather than a failed control:
 * the withhold is structural, so neither arm names anybody and there is nothing for parentage to
 * move. A band-F figure that DID separate would mean the withhold had been reintroduced as a
 * threshold, which is the state this note exists to prevent.
 *
 * Kept out of BAND_ACCURACY, as before, because a number meaning "this is a coin flip" does not
 * belong in a column of numbers meaning "this is how often it is right".
 */
export const BAND_F_AT_CHANCE: Partial<Record<Material, number>> = {
  'esc-single': 0.5257,
  trophectoderm: 0.5088,
  blastomere: 0.5557,
}

// ---------------------------------------------------------------------------------------------
// The one cell where a number is still withheld.

/** Materials on which the class-inversion veto applies. Bulk resolves its own class and is exempt. */
const AMPLIFIED: ReadonlySet<Material> = new Set<Material>(['esc-single', 'trophectoderm', 'blastomere'])

/** Implied fraction below which a possible gain cannot be told from a possible loss. */
export const VETO_MAX_F = 0.15
/** |z| the intensity channel must reach before it is treated as having resolved a direction. */
export const VETO_DIRECTION_Z = 2.576

/**
 * The one place the band structure is not honest, reported as a defect rather than shipped over.
 *
 * Amplified material, a true gain at f <= 0.15 affecting the un-genotyped parent's copy, on rows
 * that nonetheless reach band A or B: measured accuracy 0 of 34 on trophectoderm and 0 of 18 on
 * blastomere. Not weakly wrong, RELIABLY INVERTED, and the pooled band statistics hide it because
 * the cell is 430 rows out of 207,360.
 *
 * NO GATE ON OBSERVABLES FIXES THIS, and that was searched rather than assumed: nine threshold
 * combinations, the best of which catches 48% of the cell while touching 9.7% of amplified band A/B
 * rows and moving their accuracy only 0.9747 to 0.9779. Pushing to catch 63% costs a third of
 * top-band coverage. The reason no gate works is that these rows are inverted PRECISELY BECAUSE the
 * observables cannot separate a small gain from a small loss, so no function of those observables
 * can separate them either. That is an identifiability limit, not a tuning failure.
 *
 * So the parent is withheld here and the event is still reported. The remedy is categorical rather
 * than incremental: with a second parental array a gain shows a third parental haplotype and the
 * class stops being a threshold call at all.
 */
export function classInvertedRisk(
  material: Material,
  impliedGainF: number,
  intensityZ: number | undefined,
): boolean {
  if (!AMPLIFIED.has(material)) return false
  if (!(impliedGainF < VETO_MAX_F)) return false
  // Direction is resolved only if intensity says so at 99%. Absent intensity, it is unresolved.
  const resolved = intensityZ !== undefined && Number.isFinite(intensityZ)
    && Math.abs(intensityZ) >= VETO_DIRECTION_Z
  return !resolved
}

// ---------------------------------------------------------------------------------------------
// The posterior itself.

export interface PosteriorInput {
  /** Self-referenced displacement of the parent-allele share over the event window. */
  shift: number
  /** Standard error of that displacement, from the array's own genome at the same window size. */
  shiftSd: number
  /** Self-referenced displacement of window log2R. Omitted, the class is left to the prior. */
  logR?: number
  /** Standard error of the intensity displacement. */
  logRSd?: number
  material: Material
  /** Informative markers behind the shift, carried for the decomposition rather than the algebra. */
  markers?: number
}

export interface PosteriorOptions {
  /** Prior over the three classes. Defaults to uniform, which is what marginalising means here. */
  /**
   * A SOFT prior over the copy-number class. Never a point mass.
   *
   * THE TRAP THIS EXISTS TO DOCUMENT. The caller resolves the class from the same intensity and
   * allele data this function then scores, so handing that estimate back as a hard prior conditions
   * the posterior on a deterministic function of its own conditioning data. It adds no information
   * and discards the uncertainty in the estimate, which makes the result OVERCONFIDENT: the band
   * label becomes a lie rather than merely a wrong number.
   *
   * Measured over 6,000 simulated events, band A defined as confidence at or above 0.985:
   *
   *   prior                                band A n   band A acc   confidently wrong
   *   uniform, intensity supplied            2913       0.9997        0.0002
   *   SOFT, from the class likelihood        3607       0.9939        0.0037
   *   HARD, on the ESTIMATED class           5243       0.8741        0.1100
   *   HARD, on the TRUE class                5193       1.0000        0.0000
   *
   * A run that hard-conditions on the estimated class reports 0.8741 inside a band labelled 0.985
   * and gets 11 percent of ALL events confidently wrong. The class estimated from this data is
   * correct about 0.72 of the time; hard conditioning needs roughly 0.98 for band A to keep its
   * label, which is two to three orders of magnitude of evidence away.
   *
   * An earlier attempt here wired the resolved class straight in and measured 48 of 48 correct at
   * band A. That is the signature of the fourth row, not the third: it is what conditioning on the
   * TRUE class looks like, and it is unreachable in production.
   *
   * So: a distribution, never a point mass, and only when the evidence behind it is not also being
   * supplied to the likelihood. `assertSoft` below enforces the first half.
   */
  classPrior?: Partial<Record<EventClass, number>>
  /**
   * Permit a point mass, for verifying a model identity. PRODUCTION MUST NEVER SET THIS.
   *
   * Exact cancellation of the intensity term between the two parental hypotheses is a property of
   * the model at a fixed class and fraction, and demonstrating it needs the class actually fixed.
   * That is a statement about the arithmetic, not a way to score a real event.
   *
   * It is a separate, awkwardly named flag rather than a looser threshold because the failure being
   * guarded is a caller wiring in the class it just estimated. Nobody does that by writing
   * `pointMassClassPrior: true`. originPosterior.check.ts asserts that no file outside a check sets
   * it, so the backdoor cannot quietly become a door.
   */
  pointMassClassPrior?: boolean
  /** Fraction grid to marginalise over. Defaults to 0.01 to 0.70, the range the floors were measured on. */
  fGrid?: readonly number[]
  /** Per-material, per-class calibration maps. Without one the result is flagged uncalibrated. */
  calibration?: Partial<Record<Material, Partial<Record<EventClass | 'marginal', CalibrationMap>>>>
  /**
   * Whether the origin evidence rests on the obligate-het channel alone, which caps the band at B.
   *
   * That cap is structural rather than a power argument. The channel asks whether a forbidden
   * allele appeared where the homozygous parent had none to give, so its failure mode is dropout,
   * and dropout is the same event as the observation. No amount of data separates them. Measured
   * boundary: the pooled heterozygote share (1-f)/(2-f) does not cross the AA/AB decision boundary
   * at BAF 0.19 until f = 0.765, so the channel is structurally blind below roughly three quarters.
   */
  obligateHetOnly?: boolean
}

export interface OriginPosterior {
  /** Probability the affected copy came from the parent that was NOT loaded. */
  pOther: number
  /** The same after the isotonic map, which is the number a reader should see. */
  confidence: number
  /** Which parent the posterior names, and how sure it is of that rather than of "other". */
  parent: 'loaded' | 'other' | 'withheld'
  band: Band
  /** True where the number is the raw posterior because no map was supplied for this cell. */
  uncalibrated: boolean
  /** Posterior over the classes, which is what the intensity channel actually buys. */
  classPosterior: Record<EventClass, number>
  /** The most probable class, and whether it is resolved enough to state. */
  classResolved: EventClass | 'unresolved'
  /** Why the confidence is what it is, since the four causes have different remedies. */
  limitedBy: 'fraction' | 'amplification' | 'markers' | 'class-unresolved' | 'none'
  why: string
}

/** Default grid. Fine enough that the marginal is smooth, coarse enough to stay cheap in a browser. */
const DEFAULT_F_GRID: readonly number[] = Array.from({ length: 70 }, (_, i) => (i + 1) / 100)

/**
 * Marginalise class and fraction, and return a calibrated probability for the un-genotyped parent.
 *
 * The two hypotheses are which parent's copy was affected. Both are scored against the SAME
 * observations under every class and fraction, so a class that would invert the reading is not
 * discarded, it is weighed. That is the whole point: when the observations cannot tell a small gain
 * from a small loss the two hypotheses receive nearly equal support and the posterior sits near
 * 0.5, which is the honest statement. The old code resolved that ambiguity by fiat, in favour of
 * loss, and was therefore confidently wrong on every gain.
 */
/**
 * Most weight one class may carry in a supplied prior.
 *
 * Above this the prior is a point mass in all but name. Set at the accuracy hard conditioning
 * would need for band A to keep its label, 0.98, which the class estimate on this material does
 * not come close to reaching.
 */
export const MAX_CLASS_PRIOR = 0.98

export function originPosterior(
  input: PosteriorInput,
  opts: PosteriorOptions = {},
): OriginPosterior {
  const { shift, shiftSd, logR, logRSd, material } = input
  const grid = opts.fGrid ?? DEFAULT_F_GRID
  const prior = { loss: 1 / 3, gain: 1 / 3, 'cnn-loh': 1 / 3, ...opts.classPrior }

  const bad = (why: string): OriginPosterior => ({
    pOther: NaN, confidence: NaN, parent: 'withheld', band: 'F', uncalibrated: true,
    classPosterior: { loss: NaN, gain: NaN, 'cnn-loh': NaN }, classResolved: 'unresolved',
    limitedBy: 'none', why,
  })

  // A POINT MASS IS REFUSED RATHER THAN HONOURED. See the note on classPrior: conditioning on a
  // class estimated from this same data reports 0.8741 accuracy inside a band labelled 0.985. A
  // caller that means to supply evidence about the class must supply its uncertainty with it.
  const priorWeights: number[] = Object.values(prior)
  const priorTotal = priorWeights.reduce((a, x) => a + x, 0)
  const heaviest = priorTotal > 0
    ? Math.max(...priorWeights.map((x) => x / priorTotal)) : 0
  if (!opts.pointMassClassPrior && heaviest > MAX_CLASS_PRIOR) {
    return bad(`a class prior putting ${heaviest.toFixed(3)} on one class is a point mass in all `
      + `but name, over the ${MAX_CLASS_PRIOR} allowed. The class is estimated from the same data `
      + 'this scores, so conditioning that hard on it discards the uncertainty in the estimate and '
      + 'overstates the result. Supply a distribution, or nothing.')
  }

  if (!Number.isFinite(shift) || !(shiftSd > 0)) {
    return bad('no self-referenced shift or no scale to read it against, so no posterior exists')
  }

  // Joint weight for every (hypothesis, class, fraction) cell. The intensity term is shared by both
  // hypotheses within a class, which is exactly why it informs the class and never the origin.
  let wLoaded = 0
  let wOther = 0
  const byClass: Record<EventClass, number> = { loss: 0, gain: 0, 'cnn-loh': 0 }
  const useIntensity = logR !== undefined && Number.isFinite(logR) && (logRSd ?? 0) > 0

  for (const cls of CLASSES) {
    const pc = prior[cls] ?? 0
    if (!(pc > 0)) continue
    for (const f of grid) {
      const pf = 1 / grid.length
      const lIntensity = useIntensity
        ? normalPdf(logR as number, logRMean(cls, f), logRSd as number)
        : 1
      const base = pc * pf * lIntensity
      const a = base * normalPdf(shift, shiftMean(cls, f, 'loaded'), shiftSd)
      const b = base * normalPdf(shift, shiftMean(cls, f, 'other'), shiftSd)
      wLoaded += a
      wOther += b
      byClass[cls] += a + b
    }
  }

  const total = wLoaded + wOther
  if (!(total > 0)) return bad('every hypothesis received zero weight, which means the shift lies '
    + 'far outside the range any class and fraction could produce')

  const pOther = wOther / total
  const classPosterior = {
    loss: byClass.loss / total, gain: byClass.gain / total, 'cnn-loh': byClass['cnn-loh'] / total,
  }
  // A class is stated only when it carries most of the posterior mass. On amplified material this
  // is usually false, and that is the finding rather than a failure.
  let top: EventClass = 'loss'
  for (const c of CLASSES) if (classPosterior[c] > classPosterior[top]) top = c
  const classResolved: EventClass | 'unresolved' = classPosterior[top] >= 0.90 ? top : 'unresolved'

  // Confidence is in the named direction, not in "other": a posterior of 0.02 for the other parent
  // is a confident call for the loaded one.
  const rawConfidence = Math.max(pOther, 1 - pOther)
  const named: 'loaded' | 'other' = pOther >= 0.5 ? 'other' : 'loaded'

  const maps = opts.calibration?.[material]
  const map = maps?.[classResolved === 'unresolved' ? 'marginal' : classResolved] ?? maps?.marginal
  const calibrated = map ? applyCalibration(map, rawConfidence) : rawConfidence
  const uncalibrated = !map

  let band = bandOf(calibrated)
  // The obligate-het cap, applied after the band and before anything reads it.
  if (opts.obligateHetOnly && band === 'A') band = 'B'

  // The decomposition, because low confidence from four different causes has four different
  // remedies and a single number hides which one applies.
  const impliedF = fractionAt(classResolved === 'unresolved' ? 'loss' : classResolved, Math.abs(shift))
  let limitedBy: OriginPosterior['limitedBy'] = 'none'
  if (band === 'A') limitedBy = 'none'
  else if (classResolved === 'unresolved' && material !== 'bulk') limitedBy = 'class-unresolved'
  else if (impliedF <= 0.15) limitedBy = 'fraction'
  else if ((input.markers ?? Infinity) < 400) limitedBy = 'markers'
  else limitedBy = 'amplification'

  const pct = (x: number) => (100 * x).toFixed(1)
  const why = `posterior ${calibrated.toFixed(4)} for the `
    + `${named === 'other' ? 'un-genotyped' : 'loaded'} parent, band ${band} `
    + `(${BAND_LABEL[band]}), marginalised over the three copy-number classes and over mosaic `
    + `fraction. Class posterior: loss ${pct(classPosterior.loss)}%, gain `
    + `${pct(classPosterior.gain)}%, copy-neutral ${pct(classPosterior['cnn-loh'])}%`
    + (classResolved === 'unresolved'
      ? '. The class is NOT resolved, and because a gain inverts the sign that loss and '
        + 'copy-neutral LOH share, the posterior hedges across that inversion rather than assuming '
        + 'one. That hedging is why the number is what it is'
      : `, resolved as ${classResolved}`)
    + (uncalibrated
      ? '. No isotonic recalibration map was supplied for this cell, so this is the RAW posterior. '
        + 'Measured on 140,000 injections carrying the real noise of 35 arrays of this platform, '
        + 'the raw posterior has an expected calibration error of 0.0039 and errs UNDER-confident '
        + 'in the middle bands by 2 to 3 points, so the band is a lower bound rather than an '
        + 'overstatement (audit/calibration/FINDINGS.txt)'
      : '')

  return {
    pOther, confidence: calibrated, parent: named, band, uncalibrated,
    classPosterior, classResolved, limitedBy, why,
  }
}

/**
 * Fraction implied by a shift magnitude under an assumed class, kept here for the decomposition.
 *
 * This is the inverse of shiftMagnitude and it is deliberately NOT used to pick a class. It exists
 * so a reader can be told "the fraction looks small" as a reason for low confidence.
 */
export const fractionAt = (cls: EventClass, d: number): number => {
  if (!(d > 0)) return NaN
  if (cls === 'cnn-loh') return 2 * d
  if (cls === 'gain') return d >= 0.5 ? NaN : (4 * d) / (1 - 2 * d)
  return (4 * d) / (1 + 2 * d)
}

// ---------------------------------------------------------------------------------------------
// The other two origin channels.
//
// EVERY CHANNEL MUST EMIT A CONFIDENCE, not just the dosage one. Until this existed, a two-parent
// genotype call, which is the STRONGEST evidence this tool has, printed a bare parent name with no
// number, while a weak dosage call printed 0.62 beside it. A reader comparing those two rows was
// being told the opposite of the truth by the formatting alone.

/**
 * Posterior for a signed statistic whose two hypotheses sit symmetrically either side of zero.
 *
 * This is the shape both remaining channels have: a self-referenced deviation that should read
 * +separation under one parent and -separation under the other. With Gaussian noise the posterior
 * reduces to a logistic in the observed deviation, which is written out here rather than in each
 * caller so the two channels cannot drift apart.
 *
 * Returns the probability of the POSITIVE hypothesis.
 */
export function twoPointPosterior(
  deviation: number, se: number, separation: number,
): number {
  if (!Number.isFinite(deviation) || !(se > 0) || !(separation > 0)) return NaN
  // log LR = 2 * deviation * separation / se^2, clamped so a huge ratio does not overflow to NaN.
  const l = Math.max(-700, Math.min(700, (2 * deviation * separation) / (se * se)))
  return 1 / (1 + Math.exp(-l))
}

/**
 * Band for a call resting on the obligate-het channel ALONE, which can never be band A.
 *
 * The cap is structural rather than a power argument, and that distinction is the whole reason it
 * is enforced here rather than left to a threshold. The channel asks whether a forbidden allele
 * appeared where the homozygous parent had none to give, so its failure mode is DROPOUT, and
 * dropout is the same event as the observation: not seeing the allele is what both "she did not
 * transmit it" and "amplification lost it" look like. No quantity of markers separates them.
 *
 * Measured boundary for the same conclusion arriving from a different direction: the pooled
 * heterozygote share (1-f)/(2-f) does not cross the AA/AB genotype boundary at BAF 0.19 until
 * f = 0.765, so the channel is structurally blind below roughly three quarters.
 */
export const bandObligateHet = (posterior: number): Band => {
  const b = bandOf(posterior)
  return b === 'A' ? 'B' : b
}

/**
 * THE GENOTYPE CHANNELS' POSTERIORS SATURATE, so their BAND is the output and their digits are not.
 *
 * Measured on the shipped code rather than suspected. The obligate-het posterior returns 1.0000 at
 * 0, 2, 8, 40, 120 and 300 exclusive markers out of 400: it flips its VERDICT between 8 and 40 and
 * never moves its NUMBER. The two-parent share model behaves the same way, reading 0.9974 at the
 * margin and 1.0000 everywhere above it.
 *
 * That is not a bug in either model, it is what a likelihood ratio over hundreds of near-independent
 * Mendelian markers does when the model is taken at face value. The real error rate on those
 * channels is set by things no likelihood here represents: contamination, a mis-specified dropout
 * rate, a sample that is not the sample it is labelled. So the digits carry no information a reader
 * can use, and presenting them as a varying confidence would be the same laundering this rework
 * exists to stop.
 *
 * WHAT IS DONE ABOUT IT. Both channels are capped below the top band, both report themselves
 * uncalibrated, and this constant is what the display uses to say so in words. Only the dosage
 * posterior earned its bands from measurement, by marginalising its nuisance parameters instead of
 * plugging them in and by flooring its error with drift that does not average down. Applying that
 * same treatment to the genotype channels is the outstanding work, and it needs a parent-child
 * truth set that the corpus to hand does not contain.
 */
export const SATURATING_CHANNEL_NOTE =
  'this channel\'s posterior saturates, so the band is the meaningful output and the digits are '
  + 'not: the same number is returned across a wide range of evidence. It is capped below the top '
  + 'band and reported as uncalibrated for that reason'

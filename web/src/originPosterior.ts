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
 * WHAT IS MEASURED HERE AND WHAT IS INHERITED. The algebra, the sign convention and the band
 * arithmetic are computed here and checked in originPosterior.check.ts. The band accuracies and the
 * wrong-parent rates come from a methods review that ran 276,480 real-noise injections over 14
 * arrays of this platform with leave-one-array-out recalibration (audit/CONSULT-calibrated-origin.txt,
 * audit/bands_measured.csv), and are inherited.
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
export type Band = 'A' | 'B' | 'C' | 'D'

export const BAND_A_MIN = 0.985
export const BAND_B_MIN = 0.90
export const BAND_C_MIN = 0.75

export const bandOf = (confidence: number): Band => {
  if (!(confidence >= 0)) return 'D'
  if (confidence >= BAND_A_MIN) return 'A'
  if (confidence >= BAND_B_MIN) return 'B'
  if (confidence >= BAND_C_MIN) return 'C'
  return 'D'
}

export const BAND_LABEL: Record<Band, string> = {
  A: 'very confident',
  B: 'confident',
  C: 'weak, direction only',
  D: 'weak, not for reporting',
}

/**
 * Measured accuracy per band per material, from the injection experiment.
 *
 * Array-clustered, matching audit/bands_measured.csv. A naive Wilson interval on pooled rows is 3.5
 * to 5.8x too narrow on amplified material, because rows from one array are not independent.
 */
export const BAND_ACCURACY: Record<Material, Record<Band, number>> = {
  bulk: { A: 0.9972, B: 0.9599, C: 0.8375, D: 0.6038 },
  'esc-single': { A: 0.9980, B: 0.9481, C: 0.8219, D: 0.6360 },
  trophectoderm: { A: 0.9952, B: 0.9491, C: 0.8138, D: 0.6346 },
  blastomere: { A: 0.9971, B: 0.9448, C: 0.8128, D: 0.6185 },
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
  classPrior?: Partial<Record<EventClass, number>>
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
export function originPosterior(
  input: PosteriorInput,
  opts: PosteriorOptions = {},
): OriginPosterior {
  const { shift, shiftSd, logR, logRSd, material } = input
  const grid = opts.fGrid ?? DEFAULT_F_GRID
  const prior = { loss: 1 / 3, gain: 1 / 3, 'cnn-loh': 1 / 3, ...opts.classPrior }

  const bad = (why: string): OriginPosterior => ({
    pOther: NaN, confidence: NaN, parent: 'withheld', band: 'D', uncalibrated: true,
    classPosterior: { loss: NaN, gain: NaN, 'cnn-loh': NaN }, classResolved: 'unresolved',
    limitedBy: 'none', why,
  })
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

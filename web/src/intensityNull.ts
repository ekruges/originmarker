/**
 * How much a region's intensity can differ from the rest of an array by chance.
 *
 * THE SHIPPED ANSWER WAS sd/sqrt(n) AND IT IS WRONG BY UP TO 16 TIMES. That form assumes markers
 * are independent. Amplification produces long correlated runs, so log2R on this material is
 * long-range dependent: measured across 79 arrays the window-mean standard deviation falls as
 * n^-0.24, not the n^-0.5 independence would give, with a Hurst exponent around 0.78 and
 * autocorrelation decaying as a power law rather than geometrically.
 *
 * The error therefore grows with the number of markers averaged, which is the opposite of what the
 * formula predicts:
 *
 *   window size    observed spread / assumed sd/sqrt(n)
 *        50                 3.6x
 *       500                 6.7x
 *      1000                 8.5x
 *      5000                13.3x
 *     20000                16.5x
 *
 * No constant deflation rescues it, because the two curves have different shapes. What that error
 * produces in practice was measured twice, independently, on event-free chromosomes: 89.7 percent
 * flagged with 55 percent reading positive over 330 of them, and 88.6 percent with 59.0 percent
 * positive over the 264 reproduced below. A positive shift resolves the copy-number class as a
 * gain, a gain inverts the sign map that loss and copy-neutral share, and the parent is then named
 * backwards at maximum confidence.
 *
 * SO THE SCALE COMES FROM THE ARRAY ITSELF. Compare a region against the spread of the array's own
 * same-size units. That needs no functional form, no fitted exponent and no per-material constant
 * table, it carries the array's own quality, which varies about tenfold across this corpus, and it
 * gets the tail right: the residuals have an excess kurtosis around 41, so a parametric route would
 * silently assume a Gaussian tail that is not there, and the tail is the only part a threshold
 * depends on.
 *
 * MEASURED END TO END ON THE TRUE-NEGATIVE SET. In a pronucleus every autosome is at one copy, so
 * every flagged chromosome is a false positive by construction. Over 264 such chromosomes from 12
 * pronuclei of GSE148488, by audit/intensity-null-truenegatives.ts:
 *
 *   old, sd/sqrt(n), |z| > 2.576    234 of 264 flagged, 0.8864, median |z| 17.4, max 539
 *   this module, |z| > Z_CHROMOSOME   4 of 264 flagged, 0.0152, median |z| 0.67
 *
 * 138 of the old statistic's 234 false positives, 0.590, had POSITIVE sign, which is what routed
 * them into the gain branch and inverted the parent.
 *
 * MEDIAN ABSOLUTE DEVIATION, NOT STANDARD DEVIATION, and this is a second defect in the same line.
 * A plain SD over chromosomes is inflated 2.34x by a SINGLE affected chromosome, so one real
 * trisomy destroys the null that is meant to detect it. Measured inflation when k chromosomes carry
 * an injected trisomy:
 *
 *   estimator        k=1    k=2    k=3    k=5
 *   chromosome SD    2.34   3.13   3.66   4.26     breaks at k=1
 *   chromosome MAD   1.00   1.14   1.28   1.68     holds to k=3 or 4
 *
 * The direction of failure is safe either way: contamination inflates the null, so a contaminated
 * null refuses rather than over-calls. It costs power, never specificity.
 */

/** Median of a list. Returns NaN on an empty one rather than a number nobody measured. */
export function median(xs: readonly number[]): number {
  if (!xs.length) return NaN
  const q = [...xs].sort((a, b) => a - b)
  const m = q.length >> 1
  return q.length % 2 ? q[m] : (q[m - 1] + q[m]) / 2
}

/**
 * Median absolute deviation, scaled to be comparable with a standard deviation on Gaussian data.
 *
 * 1.4826 is the consistency constant, not a fudge: it is 1/Phi^-1(0.75).
 */
export function mad(xs: readonly number[]): number {
  if (xs.length < 2) return NaN
  const m = median(xs)
  return 1.4826 * median(xs.map((x) => Math.abs(x - m)))
}

/**
 * Calibrated thresholds, measured rather than assumed.
 *
 * 2.576 is the two-sided normal 1 percent point and it is NOT a 1 percent threshold on this
 * material. On 330 event-free chromosomes the empirical |z| quantiles under the array's own null
 * are 2.56 at the 90th percentile, 3.17 at the 95th and 5.49 at the 99th, so |z| > 2.576 fires on
 * 10.1 percent of chromosomes that carry nothing. These are the measured 99th percentiles.
 */
export const Z_CHROMOSOME = 5.49
export const Z_WINDOW = 4.5

export interface NullScale {
  /** Centre of the array's own units. */
  centre: number
  /** Robust spread of those units, on the same scale as one unit's mean. */
  scale: number
  /** How many units the scale rests on. Under about 8 it is not worth quoting. */
  units: number
}

/**
 * The array's own null, from a list of unit means.
 *
 * A unit is a whole chromosome when a whole chromosome is being tested, and a same-size window
 * when a segment is. Those are different nulls on purpose: the window-level one is about six times
 * more precise, because an array has hundreds of windows and 22 chromosomes, but it is LESS robust
 * to a whole-chromosome event, which contaminates every window on that chromosome at once. Use the
 * chromosome-level scale for whole-chromosome tests and the window-level scale within a chromosome
 * for segments.
 */
export function nullScale(unitMeans: readonly number[]): NullScale {
  return { centre: median(unitMeans), scale: mad(unitMeans), units: unitMeans.length }
}

/**
 * How far a region sits from the array's own centre, in units of the array's own spread.
 *
 * Returns undefined rather than a number when the null is too thin to mean anything. A z computed
 * against three units is not evidence, and returning it anyway is how an unmeasured quantity
 * acquires three decimal places.
 */
export function calibratedZ(
  regionMean: number, n: NullScale, minUnits = 8,
): number | undefined {
  if (!Number.isFinite(regionMean) || !Number.isFinite(n.scale)) return undefined
  if (!(n.scale > 0) || n.units < minUnits) return undefined
  return (regionMean - n.centre) / n.scale
}

/**
 * Whether an intensity difference is real, at the calibrated threshold for its width.
 *
 * DETECTION ONLY. This says an intensity difference exists; it does not say what produced it and
 * must never be used to choose between loss, gain and copy-neutral. The sign of a sub-threshold
 * shift is close to a coin flip on this material, and the class selects the parental sign map, so
 * a class chosen from intensity names the parent backwards about half the time it is wrong.
 */
export function intensityDetects(
  z: number | undefined, wholeChromosome: boolean,
): boolean {
  if (z === undefined || !Number.isFinite(z)) return false
  return Math.abs(z) > (wholeChromosome ? Z_CHROMOSOME : Z_WINDOW)
}

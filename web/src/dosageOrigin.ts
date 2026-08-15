/**
 * Copy-number imbalance and, where the evidence supports it, its parent, from ALLELE DOSAGE.
 *
 * WHY A DOSAGE CHANNEL AT ALL. A whole chromosome is DETECTED by the collapse of its genotype call
 * rate, and `oneParentOrigin.ts` assigns from genotypes, so on exactly those events the genotype
 * channel has no evidence left. Dosage is read whether or not a genotype is emitted. Measured on
 * GSE148488: all four segmental losses scored from genotypes and all three whole-chromosome losses
 * refused.
 *
 * THIS IS A REWRITE. The first version was a band-occupancy likelihood, and an external methods
 * review (audit/MOSAIC-AUDIT.txt) measured its null on real material. Four things it had wrong,
 * kept here because each is a trap worth naming rather than quietly deleting:
 *
 *   THE NULL WAS OFF-CENTRE AND POINTED AT THE WRONG PARENT. Raw one-parent centroid null medians
 *   are -0.031 on trophectoderm and -0.023 on blastomeres under NO event, and on TE that offset is
 *   the shift a mosaic fraction of 0.117 would produce. The orientation step causes it: the loaded
 *   parent's allele is the one always identifiable, so dropout of the other parent's allele moves
 *   readings toward it more often than the reverse. Symmetry between the two DIRECTIONS of a
 *   statistic, which the old check asserted and which was true, says nothing about where its null
 *   sits. Self-referencing against the array's own clean genome is therefore not optional; it
 *   brings the null median to -0.001..+0.006.
 *
 *   THE PARTIAL EXPECTATION WAS WRONG BY (2-f). An array reads POOLED dosage over all DNA in the
 *   well, and a loss removes template from the denominator too, so the mean of per-cell
 *   frequencies, 0.5 + 0.5f, is not it. Correct is 1/(2-f) for loss of the loaded parent's copy and
 *   (1-f)/(2-f) for loss of the other's, smaller by up to 1.95x at low f. In logit space the shift
 *   is exactly additive, -log(1-f), independent of baseline. Loss WITH reduplication of the
 *   remaining homologue is the exception: total DNA is unchanged and (1+f)/2 is right, so the two
 *   mechanisms differ about twofold in implied f and must be separated by log2R before a fraction
 *   is quoted.
 *
 *   MOSAICISM DOES NOT LIVE IN THE UNRESOLVED READINGS. That was the premise of the old guard. The
 *   unresolved zone's share of band redistribution never reaches a majority at any f at any stage,
 *   and below f = 0.3 in a blastomere it is 0.02 to 0.20 while the extremes take 0.57 to 0.92: the
 *   low-f signal is a change in the dropout BALANCE, readings tipping extreme to extreme, not
 *   intermediate readings accumulating. Counting them had power 0.01 to 0.02 at f <= 0.3.
 *
 *   BAND-THRESHOLD SELECTION DESTROYS THE MEAN. Taking a segment mean over markers passing a 0.15
 *   band threshold has power 0.00 at every f on blastomere and TE, because that selection keeps
 *   only 27-39% of true heterozygotes and the mean is then set by wherever the contaminating
 *   homozygous tail happens to sit. The central window is the fix.
 *
 * WHAT THIS DOES INSTEAD. Orient so the loaded parent's own allele is low; keep markers where that
 * parent is homozygous; take the central window; compute the window centroid MINUS the same
 * quantity on the array's own clean chromosomes; divide by a standard error whose floor is
 * systematic within-array drift rather than sampling. Then decide, in this order: could any array
 * of this kind at this width answer, is THIS array usable, is there an imbalance, and only last,
 * is the sign secure enough to name a parent.
 *
 * THAT ORDER IS DELIBERATE. Nearly every amplified array fails the MoChA quality gate, so asking
 * quality first would report a QC failure for what is really a study-design limit, and send a
 * reader to re-run a sample when what they need is a second genotyped parent or a wider interval.
 *
 * NAMING IS THE LAST QUESTION AND USUALLY THE ANSWER IS NO. Below f = 0.3 the proportion of
 * detections naming the WRONG parent runs 0.50 to 1.00, so an imbalance is reported without a
 * class unless the implied fraction clears that. This follows MoChA, which left 29% of its events
 * unassigned because power to detect an imbalance exceeded power to resolve which one it was.
 */
import type { AB } from './informativity.ts'

/** Material class, which sets every noise constant below. Amplification, not developmental age. */
export type Material = 'bulk' | 'esc-single' | 'trophectoderm' | 'blastomere'

/**
 * The central window, on the oriented scale.
 *
 * Wide on purpose. A narrow band around 0.5 would select for markers whose allelic balance
 * survived amplification intact, which is the subset a mosaic shift has already moved out of.
 */
export const WINDOW_LO = 0.20
export const WINDOW_HI = 0.80

/**
 * Variance inflation over independent sampling, measured per material.
 *
 * Bulk is essentially white. Every amplified stage is not: readings are spatially correlated with
 * an autocorrelation integral length of 1.0 to 2.2 Mb, so a 12 Mb interval holds only 5.5 to 12
 * independent blocks whatever its marker density, and effective independent markers saturate near
 * 250 per chromosome out of a nominal 900 to 1,100. Adding markers past that buys nothing.
 */
export const VIF_CHROMOSOME: Record<Material, number> = {
  bulk: 0.94, 'esc-single': 3.83, trophectoderm: 3.34, blastomere: 3.63,
}
export const VIF_SEGMENT: Record<Material, number> = {
  bulk: 1.04, 'esc-single': 2.14, trophectoderm: 1.74, blastomere: 2.31,
}

/**
 * Systematic within-array drift, in BAF units, which is the FLOOR on the standard error.
 *
 * Not sampling noise: it does not average down with more markers. Against a sampling term of
 * 0.0019 to 0.0102 it runs 1.44 to 1.90 times larger on amplified material, so an interval's
 * uncertainty is set by how much the array wanders rather than by how many markers it carries.
 * Omitting it is how a z of 15 arrives on a chromosome that is independently verified diploid.
 */
export const DRIFT_TAU: Record<Material, number> = {
  bulk: 0.0011, 'esc-single': 0.0153, trophectoderm: 0.0123, blastomere: 0.0193,
}

/**
 * Smallest mosaic fraction detectable at 80% power and a 1% false-positive rate, by material and
 * interval width. NaN means no floor was reached at any fraction up to 0.70: NOT EVALUABLE.
 *
 * These decide what this module will attempt. A 12 Mb interval is out of reach on every amplified
 * material, and a single blastomere against one parent has no floor at all, which is why both
 * return not-evaluable rather than a refusal dressed as a QC failure. The distinction matters to a
 * reader: a refusal says this array was bad, not-evaluable says no array of this kind at this
 * width could have answered.
 */
export const F80_CHROMOSOME: Record<Material, number> = {
  bulk: 0.040, 'esc-single': 0.443, trophectoderm: NaN, blastomere: NaN,
}
export const F80_SEGMENT: Record<Material, number> = {
  bulk: 0.056, 'esc-single': NaN, trophectoderm: NaN, blastomere: NaN,
}

/**
 * Implied mosaic fraction below which the SIGN is not secure, so no parent is named.
 *
 * Measured proportion of detections naming the wrong parent, with loss of the loaded parent's copy
 * as the truth: at f 0.05 it is 1.00; at 0.10, 0.86 on a TE 12 Mb segment and 1.00 on a blastomere
 * chromosome; at 0.20, 0.50 and 0.33; at 0.30, 0.15 and 0.17. It does not become tolerable until
 * 0.30, and that is what this constant is.
 */
export const SIGN_SECURE_F = 0.30

/** |z| a shift must reach before an imbalance is reported at all. Two-sided, 1% FPR. */
export const Z_DETECT = 2.576

/**
 * Array-level exclusion, adopted from MoChA rather than invented here: median SD of BAF at
 * heterozygous sites above this and the array is not analysed. It is the one directly reusable
 * gate in a literature that otherwise publishes no floor for this material.
 */
export const MAX_HET_BAF_SD = 0.11

export type DosageVerdict =
  | 'known-parent-lost'
  | 'other-parent-lost'
  | 'imbalance-unassigned'
  | 'no-imbalance'
  | 'not-evaluable'
  | 'array-excluded'

export interface DosageCall {
  verdict: DosageVerdict
  /** Signed, self-referenced centroid shift. Positive means the loaded parent's copy is short. */
  shift: number
  z: number
  /** Mosaic fraction implied by the shift, via f = 4d/(1+2d). NaN where the shift is not positive. */
  impliedF: number
  /** Markers in the central window, which is the only denominator that carries the signal. */
  window: number
  /** Markers where the loaded parent is homozygous and a dosage was read. */
  markers: number
  material: Material
  /** Smallest fraction this material and width could detect. NaN where none could. */
  floor: number
  why: string
}

/** Orient so the loaded parent's own allele reads low, whichever homozygote it is. */
export const oriented = (parent: AB, baf: number): number =>
  (parent === 'BB' ? 1 - baf : baf)

/**
 * Mosaic fraction implied by a pooled-dosage shift.
 *
 * From E[BAF] = 1/(2-f), so d = 1/(2-f) - 1/2 and f = 4d/(1+2d). NOT f = 2d, which inverts the
 * per-cell formula and understates the fraction about 1.8x.
 */
export const fractionFromShift = (d: number): number =>
  (d <= 0 ? NaN : (4 * d) / (1 + 2 * d))

/** The additive logit-space displacement a fraction f produces. Independent of baseline. */
export const logitShift = (f: number): number => -Math.log(1 - f)

/** Mean and SD of the central window, over markers where the loaded parent is homozygous. */
export function centroid(
  pairs: readonly (readonly [AB, number | null])[],
): { mean: number; sd: number; n: number; seen: number } {
  const xs: number[] = []
  let seen = 0
  for (const [parent, baf] of pairs) {
    if (parent !== 'AA' && parent !== 'BB') continue
    if (baf === null || !Number.isFinite(baf)) continue
    seen += 1
    const b = oriented(parent, baf)
    if (b >= WINDOW_LO && b <= WINDOW_HI) xs.push(b)
  }
  if (!xs.length) return { mean: NaN, sd: NaN, n: 0, seen }
  const mean = xs.reduce((a, b) => a + b, 0) / xs.length
  const sd = xs.length > 1
    ? Math.sqrt(xs.reduce((a, b) => a + (b - mean) ** 2, 0) / (xs.length - 1))
    : NaN
  return { mean, sd, n: xs.length, seen }
}

/**
 * Call an interval from allele dosage.
 *
 * `region` and `background` are (parent genotype, sample dosage) pairs. The background must come
 * from the SAME array and EXCLUDE the interval under test: it is what removes the directional null
 * bias, and a background including the region would subtract the signal along with the bias.
 */
export function callDosageOrigin(
  region: readonly (readonly [AB, number | null])[],
  background: readonly (readonly [AB, number | null])[],
  material: Material,
  opts: {
    /** True for a whole chromosome, which has its own noise scale and its own floor. */
    wholeChromosome?: boolean
    /** Median SD of BAF at the sample's heterozygous sites, for the array-level gate. */
    hetBafSd?: number
    signSecureF?: number
    zDetect?: number
  } = {},
): DosageCall {
  const whole = opts.wholeChromosome ?? false
  const signSecure = opts.signSecureF ?? SIGN_SECURE_F
  const zNeed = opts.zDetect ?? Z_DETECT
  const floor = (whole ? F80_CHROMOSOME : F80_SEGMENT)[material]

  const r = centroid(region)
  const b = centroid(background)
  const base = {
    shift: NaN, z: NaN, impliedF: NaN, window: r.n, markers: r.seen, material, floor,
  }

  // 1. COULD ANY ARRAY OF THIS KIND ANSWER AT THIS WIDTH. Asked FIRST, and before looking at the
  // data, for two reasons. A favourable draw cannot then talk the answer into existence. And it is
  // the more informative reply: telling someone their array failed quality control, when the truth
  // is that no array of its kind at this width could have answered, sends them to re-run a sample
  // instead of to change the design. Nearly every amplified array fails the MoChA gate below, so
  // asking that first would report a QC failure for what is really a study-design limit.
  if (!Number.isFinite(floor)) {
    return {
      ...base, verdict: 'not-evaluable',
      why: `a ${whole ? 'whole chromosome' : '12 Mb-scale interval'} on ${material} material has `
        + 'no detection floor at 1% false positives for any mosaic fraction up to 0.70 with one '
        + 'genotyped parent. This is not a refusal of this array: no array of this kind at this '
        + 'width could answer. Genotyping the second parent is worth about threefold in detectable '
        + 'fraction and is the largest single improvement available',
    }
  }
  // 2. IS THIS ARRAY USABLE. Kept separate from the interval question on purpose: this is a
  // property of the array, and conflating the two was exactly the old guard's fault.
  if (opts.hetBafSd !== undefined && opts.hetBafSd > MAX_HET_BAF_SD) {
    return {
      ...base, verdict: 'array-excluded',
      why: `BAF spread at heterozygous sites is ${opts.hetBafSd.toFixed(3)}, over the `
        + `${MAX_HET_BAF_SD} gate. The array is too noisy to analyse, independently of this `
        + 'interval',
    }
  }

  if (!r.n || !b.n || !Number.isFinite(r.sd)) {
    return {
      ...base, verdict: 'not-evaluable',
      why: `${r.n} markers in the central window and ${b.n} in the array's own background, which `
        + 'is not enough to form a self-referenced shift',
    }
  }

  // 3. IS THERE AN IMBALANCE. Self-referenced, with drift as the floor on the standard error.
  const shift = r.mean - b.mean
  const vif = (whole ? VIF_CHROMOSOME : VIF_SEGMENT)[material]
  const tau = DRIFT_TAU[material]
  const se = Math.sqrt((vif * r.sd * r.sd) / r.n + tau * tau)
  const z = shift / se
  const impliedF = fractionFromShift(Math.abs(shift))
  const out = { ...base, shift, z, impliedF }

  if (Math.abs(z) < zNeed) {
    return {
      ...out, verdict: 'no-imbalance',
      why: `centroid shift ${shift.toFixed(4)} against this array's own genome, z ${z.toFixed(2)}, `
        + `under the ${zNeed.toFixed(2)} needed. Standard error ${se.toFixed(4)} is floored by `
        + `within-array drift at ${tau}, which does not average down with more markers`,
    }
  }

  // 4. IS THE SIGN SECURE. Only now, and usually it is not.
  const secure = Number.isFinite(impliedF) && impliedF >= signSecure && impliedF >= floor
  if (!secure) {
    return {
      ...out, verdict: 'imbalance-unassigned',
      why: `an imbalance is present, shift ${shift.toFixed(4)} at z ${z.toFixed(2)}, but the `
        + `implied fraction ${Number.isFinite(impliedF) ? impliedF.toFixed(3) : 'n/a'} is under `
        + `the ${signSecure} at which the sign becomes secure, and below a fraction of 0.3 between `
        + 'half and all such detections name the WRONG parent. The event is reported and the '
        + 'parent is not',
    }
  }
  const lost = shift > 0
  return {
    ...out,
    verdict: lost ? 'known-parent-lost' : 'other-parent-lost',
    why: `${lost ? "the loaded parent's" : "the other parent's"} copy is short over this interval: `
      + `centroid shift ${shift.toFixed(4)} against this array's own genome at z ${z.toFixed(2)}, `
      + `implying a fraction of ${impliedF.toFixed(3)}, over both the ${signSecure} sign-security `
      + `bound and this material's ${floor} detection floor. Read from dosage, so a collapsed call `
      + 'rate does not remove the evidence',
  }
}

/**
 * The material class this module's constants are indexed by, from the stage the array was given.
 *
 * Deliberately coarse: what sets the noise is how much amplification happened, not how old the
 * embryo was. A trophectoderm biopsy and a single ES cell differ by more than a blastomere and a
 * single ES cell do.
 */
export function materialOf(stage: string): Material {
  if (stage === 'bulk') return 'bulk'
  if (stage === 'trophectoderm') return 'trophectoderm'
  if (stage === 'blastomere') return 'blastomere'
  // single-cell, haploid and anything unrecognised take the amplified single-cell constants,
  // which are the more conservative of the two remaining sets on every parameter.
  return 'esc-single'
}

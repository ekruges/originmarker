/**
 * Copy-number imbalance and, where the evidence supports it, its parent, from ALLELE DOSAGE.
 *
 * WHY A DOSAGE CHANNEL AT ALL. A whole chromosome is DETECTED by the collapse of its genotype call
 * rate, and `oneParentOrigin.ts` assigns from genotypes, so on exactly those events the genotype
 * channel has less to work with. Dosage is read whether or not a genotype is emitted. Measured on
 * GSE148488: all four segmental losses scored from genotypes and all three whole-chromosome losses
 * refused.
 *
 * THAT MECHANISM IS HALF WRONG AND THE CORRECTION MATTERS. A later review measured what is
 * actually left at a vendor no-call: the discrete genotype carries 0.0000 bits about origin, but
 * the continuous BAF carries 0.0471, and BAF is present at 88-91% of no-call markers. The cause is
 * low intensity rather than an off-cluster reading. So it is DISCRETISATION that throws the
 * information away, half of it at called markers and all of it at no-called ones, and the evidence
 * is not destroyed by the event. The refusal on whole-chromosome losses still stands, but for a
 * different reason: in this dataset a loss large enough to see in a single cell co-occurs with an
 * array too damaged to self-reference, and those two conditions are perfectly confounded. All 70
 * whole-chromosome losses in the series sit on arrays with 40 to 100% of autosomes deviant.
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
 * which parent.
 *
 * THAT ORDER IS DELIBERATE. The width question is asked before the array question so that a limit
 * of the study design is never reported as a fault of the sample, which would send a reader to
 * re-run an array when what they need is a second genotyped parent or a wider interval.
 *
 * AND THE ARRAY QUESTION IS NOW ONLY THE STRUCTURAL ONE. A BAF-spread gate adopted from MoChA used
 * to sit here at 0.11 and refused every amplified array: 9 of the 14 worked examples and 852 of 877
 * in the reference corpus, which is every single cell, blastomere and trophectoderm biopsy. On
 * exactly the material this tool exists for, no parent was ever named. It was redundant with the
 * bands rather than protective of them, and that was measured: over 140,000 injections carrying the
 * real per-chromosome noise of 35 arrays, band A accuracy on the noisiest third was 0.9993 against
 * 0.9999 on the cleanest. Array noise already reaches the answer through the self-referenced
 * standard error, so a noisy array earns a lower band rather than a refusal. What remains is the
 * refusal with a reason: a genome with no undisturbed remainder has nothing to self-reference
 * against, which is a property of the genome rather than a verdict on the array.
 *
 * NAMING NO LONGER READS THE SIGN, AND THAT IS THE CORRECTION THIS MODULE MOST NEEDED. It used to
 * ask whether the implied fraction cleared 0.30 and then take the sign of the shift, inverting
 * through the loss formula whenever nobody had resolved the class. Both halves of that were wrong.
 *
 *   THE SIGN DOES NOT NAME A PARENT ON ITS OWN. Gain inverts the map that loss and copy-neutral LOH
 *   share: at f = 0.10 the loaded parent's copy reads +0.0263 under a loss and -0.0238 under a
 *   gain. Defaulting to loss is not conservative, it is an assertion, and this project's own audit
 *   puts the class as unresolved on 89 to 100% of amplified detections. Measured wrong-parent rate
 *   for a true low-fraction gain scored that way: 0.551 to 0.580 at f = 0.05, across four material
 *   classes. Worse than chance and systematic, because a sign inversion is not noise.
 *
 *   THE 0.30 THRESHOLD WAS GUARDING THE WRONG FAILURE. The posterior is measurably UNDER-confident
 *   before recalibration, not over-confident, so the blanket withhold was discarding honest calls:
 *   on blastomeres it refused 24.1% of events that sit in a band measured at 0.9971 accuracy. A
 *   single threshold on fraction cannot express this in any case, since fraction explains only
 *   51.7% of the variance in achievable confidence and material, array identity, marker count and
 *   class carry the rest.
 *
 * So the class is marginalised into a posterior, the intensity channel supplies what direction
 * information exists, and the result is a calibrated probability placed in one of four bands. Every
 * band carries its number, including the weakest, because each is calibrated within itself. The one
 * cell where a number is still withheld is narrow and named: see classInvertedRisk.
 */
import type { AB } from './informativity.ts'
import {
  originPosterior, classInvertedRisk, BAND_LABEL, VETO_MAX_F, SHIPPED_MAPS,
  type OriginPosterior, type CalibrationMap, type EventClass,
} from './originPosterior.ts'

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
 * Detection floors by material, STATE and how many parents are loaded.
 *
 * The loss-only tables above were what this module shipped with, and they made it refuse work its
 * own evidence supports. Three things the fuller measurement changes:
 *
 *   COPY-NEUTRAL LOH IS THE EASIEST STATE, not the hardest. Copy number stays at 2 so the genotype
 *   caller never degrades: median call rate 0.867 at copy-neutral log2R against 0.287 below -1.5.
 *   Its deviation is f/2, the largest of the three. On a trophectoderm biopsy with ONE parent its
 *   origin floor is 0.186, comfortably callable, where a loss on the same material is 0.625.
 *
 *   A SECOND PARENT IS WORTH MORE THAN MORE CELLS. On the same single file it moves a TE whole
 *   chromosome from 0.625 to 0.186, a factor of 3.36, where going from one cell to five moves
 *   0.186 to 0.135, a factor of 1.38.
 *
 *   ONE COMBINATION IS A MATERIAL LIMIT AND NO INPUT MOVES IT. A single blastomere, one parent, a
 *   loss: no floor at any fraction to 0.70 on any channel, and still none at eight cells.
 *
 * NaN means no floor was reached, which is not-evaluable rather than a refusal of this array.
 */
export interface Floors { chromosome: number; segment: number }
const F: (c: number, sg: number) => Floors = (chromosome, segment) => ({ chromosome, segment })

const FLOOR_ONE_PARENT: Record<DosageState, Record<Material, Floors>> = {
  loss: {
    bulk: F(0.050, 0.056),
    'esc-single': F(0.348, NaN),
    trophectoderm: F(0.625, NaN),
    blastomere: F(NaN, NaN),
  },
  'cnn-loh': {
    bulk: F(0.040, 0.040),
    'esc-single': F(0.399, NaN),
    trophectoderm: F(0.186, 0.186),
    blastomere: F(0.399, NaN),
  },
  gain: {
    bulk: F(0.044, 0.044),
    'esc-single': F(NaN, NaN),
    trophectoderm: F(0.511, NaN),
    blastomere: F(0.620, NaN),
  },
}

const FLOOR_TWO_PARENTS: Record<DosageState, Record<Material, Floors>> = {
  loss: {
    bulk: F(0.040, 0.044),
    'esc-single': F(0.232, NaN),
    trophectoderm: F(0.186, 0.327),
    blastomere: F(0.628, NaN),
  },
  'cnn-loh': {
    bulk: F(0.040, 0.040),
    'esc-single': F(0.232, NaN),
    trophectoderm: F(0.135, 0.186),
    blastomere: F(0.399, NaN),
  },
  gain: {
    bulk: F(0.044, 0.044),
    'esc-single': F(0.443, NaN),
    trophectoderm: F(0.511, NaN),
    blastomere: F(0.620, NaN),
  },
}

/**
 * Whether ANY array of this kind, at this width, could name a parent at any mosaic fraction.
 *
 * Exported because the answer is a table lookup that needs no data, and a caller that assembles
 * the data first pays for a background it is about to be told is irrelevant. On single-cell
 * material every SEGMENT is undefined at every fraction, so a run whose findings are all segments
 * spends most of its time preparing evidence for a question already known to be unanswerable.
 *
 * THE CALLER AND STEP 1 OF `callDosageOrigin` MUST AGREE. They read this one predicate rather than
 * each testing the floor themselves, so a caller cannot skip the background for a call that then
 * turns out to need it. `dosageOrigin.check.ts` walks the full cross product to hold that.
 */
export const originUnreachable = (
  material: Material, state: DosageState, wholeChromosome: boolean, parents: 1 | 2,
): boolean => !Number.isFinite(floorFor(material, state, wholeChromosome, parents))

/** The floor for one combination. `parents` is how many parental arrays are loaded. */
export function floorFor(
  material: Material, state: DosageState, wholeChromosome: boolean, parents: 1 | 2,
): number {
  const t = (parents === 2 ? FLOOR_TWO_PARENTS : FLOOR_ONE_PARENT)[state][material]
  return wholeChromosome ? t.chromosome : t.segment
}

/**
 * Fraction of detected events whose ORIGIN resolves while the CLASS does not, by material.
 *
 * The measurement that most changes what this module should emit. Power to detect an allelic
 * imbalance exceeds power to resolve which state produced it, by a lot on exactly the material
 * this tool targets: at 400 informative markers it is 0.296 on bulk but 1.000 on trophectoderm and
 * blastomere. Lifting it needs a 3 to 5x reduction in window log2R spread, measured at 0.17-0.22
 * against the 0.029-0.081 required, and four times the markers buys only 1.2-1.4x. It is not
 * reachable by collecting more of the same.
 *
 * So a caller that must emit a class refuses most of its own detections. Emitting origin and class
 * as separate fields with separate confidence is what turns those into answers.
 */
export const ORIGIN_WITHOUT_CLASS: Record<Material, number> = {
  bulk: 0.296, 'esc-single': 0.759, trophectoderm: 1.000, blastomere: 1.000,
}

/** |z| a shift must reach before an imbalance is reported at all. Two-sided, 1% FPR. */
export const Z_DETECT = 2.576

/**
 * Array-level exclusion, adopted from MoChA rather than invented here: median SD of BAF at
 * heterozygous sites above this and the array is not analysed. It is the one directly reusable
 * gate in a literature that otherwise publishes no floor for this material.
 */
export const MAX_HET_BAF_SD = 0.11

/**
 * Residual correlation between the dosage and intensity channels, MEASURED ON THIS STATISTIC.
 *
 * A joint term has to know how much the two channels are already saying the same thing. The
 * methods review measured this two ways and they disagreed: technical-replicate DIFFERENCES gave
 * -0.055 to +0.036, replicate MEANS gave up to -0.46, the gap being a common amplification
 * artefact the differences cancel and the means do not.
 *
 * Neither is the number this code needs, because the shipped statistic is SELF-REFERENCED and that
 * subtraction removes exactly the artefact the means carry. So it was measured directly here, on
 * the shipped quantity: per array, the leave-one-out per-chromosome BAF centroid shift against the
 * leave-one-out per-chromosome mean log2R shift, over 81 arrays of GSE148488
 * (audit/asymmetry/residual_corr.ts, residual-corr.json).
 *
 *     material         arrays   median r
 *     bulk                 41     -0.058
 *     blastomere           10     +0.486
 *     esc-single           22     +0.508
 *     trophectoderm         8     +0.633
 *
 * Bulk is independent and quadrature addition would be fine there. Every amplified material is
 * strongly POSITIVELY correlated, which is the opposite of what quadrature assumes: ignoring it
 * inflates the joint z by sqrt(2+2r)/sqrt(2), which is 1.27x on trophectoderm. The mechanism is the
 * one behind the null bias: a chromosome that amplified poorly reads both lower in intensity and
 * lower in oriented dosage, because dropout moves readings toward the allele that can always be
 * identified.
 */
export const RESIDUAL_R: Record<Material, number> = {
  bulk: -0.058, 'esc-single': 0.508, trophectoderm: 0.633, blastomere: 0.486,
}

export type DosageVerdict =
  /** The loaded parent's copy is the affected one. Named from the posterior, not from the sign. */
  | 'loaded-parent'
  /** The un-genotyped parent's copy is the affected one. */
  | 'other-parent'
  /** Amplified material where a small gain and a small loss name opposite parents. See the veto. */
  | 'class-inverted-risk'
  | 'imbalance-unassigned'
  | 'no-imbalance'
  | 'not-evaluable'
  | 'array-excluded'

/** Whether the copy-number STATE could be separated from its nearest feasible alternative. */
export type ClassVerdict = 'loss' | 'gain' | 'cnn-loh' | 'unresolved'

export interface DosageCall {
  verdict: DosageVerdict
  /**
   * The copy-number class, resolved SEPARATELY from the origin and usually not resolved at all.
   *
   * Kept apart because the two have different power: on trophectoderm and blastomere material
   * every detected event resolves its origin and none resolves its class at 400 markers. A caller
   * that emitted one verdict had to refuse those, discarding an origin it could actually support.
   */
  classVerdict: ClassVerdict
  classWhy: string
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
  /**
   * The calibrated posterior, which is what names the parent. Absent where the call never got far
   * enough to compute one (array excluded, no imbalance, nothing to reference against).
   */
  posterior?: OriginPosterior
  why: string
}

/** Orient so the loaded parent's own allele reads low, whichever homozygote it is. */
export const oriented = (parent: AB, baf: number): number =>
  (parent === 'BB' ? 1 - baf : baf)

/**
 * The copy-number states a deviation can come from, which decide how it inverts.
 *
 * A methods review measured the full algebra by brute-force copy counting at 200,000 cells, max
 * |analytic - counted| 1.1e-16. Deviation from 0.5 at a truly heterozygous marker, oriented:
 *
 *     state                              E[parent-allele share]     deviation
 *     loss of the loaded parent's copy   (1-f)/(2-f)                -f/(4-2f)
 *     loss of the other parent's copy    1/(2-f)                    +f/(4-2f)
 *     CNN-LOH, loaded parent's lost      (1-f)/2                    -f/2
 *     CNN-LOH, other parent's lost       (1+f)/2                    +f/2
 *     gain of the loaded parent's copy   (1+f)/(2+f)                +f/(4+2f)
 *     gain of the other parent's copy    1/(2+f)                    -f/(4+2f)
 *
 * TWO THINGS THIS INVERTS ABOUT THE FRAMING THIS MODULE WAS BUILT ON. Rank order by deviation is
 * CNN-LOH > loss > gain at every f, so copy-neutral loss of heterozygosity is the LARGEST-signal
 * state, not the hardest: a class with no copy-number signal at all is the easiest one to assign a
 * parent to. And "a gain shifts about a third as much as a loss" is the f = 1 endpoint only; the
 * pooled ratio is (2-f)/(2+f), which is 0.905 at f = 0.1 and 0.818 at f = 0.2, so in the mosaic
 * range the gain penalty is 1.1 to 1.7x rather than 3x.
 */
export type DosageState = 'loss' | 'gain' | 'cnn-loh'

/**
 * Mosaic fraction implied by a pooled-dosage shift, FOR A GIVEN STATE.
 *
 * The three inversions are different and must not be shared. At d = 0.04 they give 0.148 for a
 * loss, 0.174 for a gain and 0.080 for copy-neutral LOH, against 0.240 for the per-cell 6d form.
 * An earlier version applied the loss inversion to everything: correct where the event was a loss
 * and wrong in the OPPOSITE direction on a gain, which is the worst way to be wrong, since the
 * error grows the statistic rather than shrinking it.
 */
export const fractionFromShift = (d: number, state: DosageState = 'loss'): number => {
  if (!(d > 0)) return NaN
  if (state === 'cnn-loh') return 2 * d
  // A gain's inversion diverges as d approaches 0.5, which is its own ceiling: no gain fraction
  // produces a deviation at or above it, so a larger one is not a gain.
  if (state === 'gain') return d >= 0.5 ? NaN : (4 * d) / (1 - 2 * d)
  return (4 * d) / (1 + 2 * d)
}

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
    /** Which state the deviation is being read as. Decides the inversion AND the floor. */
    state?: DosageState
    /** How many parental arrays are loaded. Two is worth 3.36x on a TE chromosome. */
    parents?: 1 | 2
    /**
     * Spread of the window log2R, which decides whether the CLASS can be separated from its
     * nearest feasible alternative. Two states are separable when they differ by more than
     * 2 x 2.576 of it. Measured 0.17-0.22 on every amplified class against the 0.029-0.081 needed,
     * so on that material the class is almost never resolvable and the origin usually is.
     */
    windowLogRSd?: number
    zDetect?: number
    /**
     * Self-referenced intensity evidence for this interval, as a signed z where NEGATIVE means
     * total dosage is reduced. Optional: omitted, the call rests on dosage alone.
     *
     * IT INFORMS THE STATE AND NEVER THE ORIGIN. On haploid pronuclei, a natural complete-loss
     * experiment with known parent, autosomal median log2R is 0.0837 maternal against 0.0750
     * paternal, p = 0.54, indistinguishable, while the oriented dosage at father-homozygous
     * markers separates them at p = 1.3e-15. So this is combined into the DETECTION step and the
     * direction still comes from the dosage sign alone.
     */
    intensityZ?: number
    /**
     * Set where the array has no undisturbed genome to reference against, which is the one
     * array-level condition that still refuses outright. See detectComplex: above 30% deviant
     * autosomes, or a call rate under 0.70, the self-reference this all depends on does not exist.
     */
    noSelfReference?: boolean
    /** Window log2R displacement in its natural units, if the caller has it directly. */
    logRShift?: number
    /** Standard error of that displacement. Derived from windowLogRSd and the window when absent. */
    logRShiftSe?: number
    /** Per-material, per-class isotonic maps. Without them the posterior is flagged uncalibrated. */
    calibration?: Partial<Record<Material, Partial<Record<EventClass | 'marginal', CalibrationMap>>>>
    /** Set where the origin evidence rests on the obligate-het channel alone, which caps at band B. */
    obligateHetOnly?: boolean
  } = {},
): DosageCall {
  const whole = opts.wholeChromosome ?? false
  const zNeed = opts.zDetect ?? Z_DETECT
  // The floor is state-aware and parent-count aware. Reading a copy-neutral event against a loss
  // floor is what made this module refuse the easiest class it has.
  const state = opts.state ?? 'loss'
  const parents = opts.parents ?? 1
  const floor = floorFor(material, state, whole, parents)

  const r = centroid(region)
  const b = centroid(background)
  const base = {
    shift: NaN, z: NaN, impliedF: NaN, window: r.n, markers: r.seen, material, floor,
    classVerdict: 'unresolved' as ClassVerdict, classWhy: '',
  }

  // 1. COULD ANY ARRAY OF THIS KIND ANSWER AT THIS WIDTH. Asked FIRST, and before looking at the
  // data, for two reasons. A favourable draw cannot then talk the answer into existence. And it is
  // the more informative reply: telling someone their array failed quality control, when the truth
  // is that no array of its kind at this width could have answered, sends them to re-run a sample
  // instead of to change the design. Nearly every amplified array fails the MoChA gate below, so
  // asking that first would report a QC failure for what is really a study-design limit.
  if (originUnreachable(material, state, whole, parents)) {
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
  // 2. IS THERE ANYTHING LEFT TO REFERENCE AGAINST. This is the structural refusal and the only
  // array-level one that survives: every statistic here is measured against the rest of this
  // array's own genome, so a genome with no undisturbed remainder has nothing to measure against.
  // A refusal here is a property of the GENOME rather than a verdict on the array's quality.
  if (opts.noSelfReference) {
    return {
      ...base, verdict: 'array-excluded',
      why: 'this array has no undisturbed genome left to reference against, so the self-referenced '
        + 'statistic every origin call here depends on would be measuring its own reference. No '
        + 'parent is named on any interval of this array, at any confidence',
    }
  }

  // THE BAF-SPREAD GATE NO LONGER BLOCKS A CALL, and removing it is the change that made this
  // feature work on the material it exists for. Adopted from MoChA as an array-level exclusion at
  // 0.11, it refused 9 of the 14 worked examples and 852 of 877 arrays in the reference corpus:
  // every amplified sample, which is every single cell, blastomere and trophectoderm biopsy. On
  // exactly the material this tool is for, no parent was ever named.
  //
  // It is redundant, which was measured rather than argued. Array noise already reaches the answer
  // through the self-referenced standard error, so a noisy array gets a wider error, a lower
  // posterior and a lower band, which is the banding the whole design rests on. Over 140,000
  // injections carrying the real per-chromosome noise of 35 arrays, band A accuracy on the NOISIEST
  // third was 0.9993 against 0.9999 on the cleanest, and overall accuracy 0.6875 against 0.7019.
  // The gate was not protecting anything the bands were not already handling; it was a blanket
  // threshold of exactly the kind the fraction threshold was, refusing wholesale where the design
  // says to band.
  //
  // The measurement is still made and still reported: a reader should know the array is noisy. It
  // informs the confidence rather than replacing it.

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
  // Inverted as a LOSS unless the caller establishes otherwise, because that is the state this
  // module's floors were measured on. A gain inverts through a different formula and a
  // copy-neutral event through a third; passing the wrong one is an error in the wrong direction.
  const impliedF = fractionFromShift(Math.abs(shift), opts.state ?? 'loss')
  const out = { ...base, shift, z, impliedF }

  // The joint detection statistic. Both terms point the same way by construction: |z| for the
  // dosage channel, which is two-sided because either parent's copy may be the one short, and the
  // reduction in intensity for the other, which is one-sided because a loss lowers total dosage
  // whichever parent it came from. Stouffer with the MEASURED correlation, not quadrature: on
  // amplified material the two channels correlate around +0.5 to +0.6 and treating them as
  // independent would overstate the combined evidence by about a quarter.
  const rho = RESIDUAL_R[material]
  const iz = opts.intensityZ
  const zEvent = iz !== undefined && Number.isFinite(iz)
    ? (Math.abs(z) + Math.max(0, -iz)) / Math.sqrt(2 + 2 * rho)
    : Math.abs(z)

  if (zEvent < zNeed) {
    return {
      ...out, verdict: 'no-imbalance',
      why: `centroid shift ${shift.toFixed(4)} against this array's own genome, z ${z.toFixed(2)}`
        + (iz !== undefined && Number.isFinite(iz)
          ? `, jointly with intensity ${zEvent.toFixed(2)} at a measured channel correlation of `
            + `${rho}`
          : '')
        + `, under the ${zNeed.toFixed(2)} needed. Standard error ${se.toFixed(4)} is floored by `
        + `within-array drift at ${tau}, which does not average down with more markers`,
    }
  }

  // THE DETECTION FLOOR DOES NOT GATE THE ASSIGNMENT, and putting it here was a mistake worth
  // recording because it looked obviously right. The f80 floor is the smallest fraction detectable
  // at 80% power; by this line a detection has ALREADY been made, at z >= 2.576 against a standard
  // error floored by drift. Conditioning on that, re-applying a power threshold refuses events the
  // array did in fact see, and it is the same single-threshold-on-fraction the rest of this rework
  // removes: fraction explains only 51.7% of the variance in achievable confidence.
  //
  // It also made the class-inversion veto below UNREACHABLE. Measured against our own constants:
  // the veto needs an implied gain fraction under 0.15, which is a shift under 0.0349, which is an
  // implied loss fraction under 0.1304, and the lowest amplified floor we carry is 0.135. So every
  // inverted cell was being caught by a gate that gave a less informative reason for it.
  //
  // The floor stays where it belongs: question 1, whether ANY array of this kind at this width
  // could answer, and as reported information on the call. Selection effects on a below-floor
  // detection are real, and they are handled where they belong, by marginalising the fraction
  // rather than by refusing the row.
  // THE CLASS IS A SEPARATE QUESTION with its own power, asked after the origin and allowed to
  // fail without taking the origin down with it.
  const sd = opts.windowLogRSd
  const separable = sd !== undefined && Number.isFinite(sd) && 2 * 2.576 * sd < 0.29
  const cls: ClassVerdict = separable ? state : 'unresolved'
  const share = `${(100 * ORIGIN_WITHOUT_CLASS[material]).toFixed(0)}%`
  const why = sd === undefined
    ? 'no window log2R spread was supplied'
    : `window log2R spread ${sd.toFixed(3)} is too wide to separate the two closest feasible states`
  const classWhy = separable
    ? `window log2R spread ${(sd as number).toFixed(3)} separates ${state} from its nearest `
      + 'feasible alternative'
    : `the class is not resolved: ${why}. On this material that is the usual outcome, since at 400 `
      + `informative markers ${share} of detected events resolve an origin without a class. `
      + 'Lifting it needs a three to fivefold narrower spread, which more markers do not deliver: '
      + 'four times as many buys 1.2 to 1.4x'

  // 4. WHICH PARENT. Marginalising the class rather than assuming one, which is the whole change.
  //
  // This used to be `const lost = shift > 0`, a naked sign read, with the class defaulted to loss
  // wherever nobody had resolved it. That default is not conservative: gain inverts the map that
  // loss and copy-neutral LOH share, so on a true gain the naked sign named the WRONG parent, at a
  // measured 0.551 to 0.580 across four material classes at f = 0.05. Since this project's own
  // audit puts the class as unresolved on 89 to 100% of amplified detections, the old line was
  // wrong in the majority regime rather than at the edges.
  //
  // The intensity channel is handed over in its natural units so it can inform the class. It still
  // names no parent of its own: within a class and a fraction its term is identical under both
  // parental hypotheses and cancels exactly.
  const seL = opts.windowLogRSd !== undefined && Number.isFinite(opts.windowLogRSd) && r.n > 0
    ? (opts.windowLogRSd as number) / Math.sqrt(r.n)
    : undefined
  const logRShift = opts.logRShift ?? (
    iz !== undefined && Number.isFinite(iz) && seL !== undefined ? iz * seL : undefined)
  const logRSe = opts.logRShiftSe ?? seL

  // The shipped maps are the default, so a caller gets a CALIBRATED number without having to know
  // they exist. Passing `calibration` overrides them; passing an empty object opts out and gets the
  // raw posterior, which announces itself as raw.
  const post = originPosterior(
    { shift, shiftSd: se, logR: logRShift, logRSd: logRSe, material, markers: r.n },
    {
      calibration: opts.calibration ?? (SHIPPED_MAPS as never),
      obligateHetOnly: opts.obligateHetOnly,
    },
  )

  // THE ONE VETO. Amplified material, a possible small gain, intensity unable to resolve a
  // direction: measured 0 of 34 correct on trophectoderm and 0 of 18 on blastomere among rows that
  // still reached the top two bands. Reliably inverted rather than weakly wrong, and no gate on
  // observables separates them, because the rows are inverted precisely because the observables
  // cannot tell a small gain from a small loss. The event is still reported; only the parent is
  // withheld.
  const gainF = fractionFromShift(Math.abs(shift), 'gain')
  if (classInvertedRisk(material, gainF, iz)) {
    return {
      ...out, classVerdict: cls, classWhy, posterior: post,
      verdict: 'class-inverted-risk',
      why: `an imbalance is present, shift ${shift.toFixed(4)} at z ${z.toFixed(2)}, but on `
        + `${material} material a shift this size is equally consistent with a small loss and a `
        + `small gain (implied gain fraction ${Number.isFinite(gainF) ? gainF.toFixed(3) : 'n/a'}, `
        + `under ${VETO_MAX_F}), and the intensity channel does not resolve the direction at 99%. `
        + 'Those two readings name OPPOSITE parents. In this exact cell the measured accuracy of a '
        + 'stated parent is 0 of 34 on trophectoderm and 0 of 18 on blastomere, so a number here '
        + 'would not be weak, it would be inverted. A second parental array resolves the class '
        + 'categorically rather than by threshold and removes this cell entirely',
    }
  }

  const named = post.parent === 'other' ? 'other-parent' : 'loaded-parent'
  return {
    ...out,
    classVerdict: cls,
    classWhy,
    posterior: post,
    verdict: named,
    why: `${post.parent === 'loaded' ? "the loaded parent's" : "the un-genotyped parent's"} copy is `
      + `the affected one, at a calibrated confidence of ${post.confidence.toFixed(4)} `
      + `(band ${post.band}, ${BAND_LABEL[post.band]}). Centroid shift ${shift.toFixed(4)} against `
      + `this array's own genome at z ${z.toFixed(2)}, over this material's ${floor} detection `
      + `floor. ${post.why}. Read from dosage, so a collapsed call rate does not remove the evidence`,
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

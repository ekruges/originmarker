/**
 * Which parent an extra copy came from, and when that cannot be said.
 *
 * A gain is not one question. It is two, and they are answered by different evidence:
 *
 *   In a BIPARENTAL cell, both parents are present and the question is which one contributed
 *   the extra copy. Allele dosage answers it: at a marker where the parents are opposite
 *   homozygotes a euploid cell carries one of each, and an extra copy tips the ratio toward
 *   whichever parent supplied it.
 *
 *   In a UNIPARENTAL cell - a haploid meiotic product - everything present came from one parent
 *   already, so "which parent" is answered by construction. The question worth asking is whether
 *   the extra copy is the OTHER homologue of that parent, which is a meiotic error, or a
 *   duplicate of the same homologue, which is mitotic. One of those is visible and one is not.
 *
 * ORIENTATION IS THE TRAP. Working in raw B-allele frequency, an extra paternal copy pushes BAF
 * DOWN where the father is AA and UP where he is BB. Averaging raw BAF across markers of both
 * kinds cancels the signal, and getting the sign wrong inverts every call. So nothing here works
 * in raw BAF: each informative marker is converted to the PATERNAL ALLELE SHARE first, which
 * points the same way at every marker whichever parent is A. Euploid sits at 0.5, an extra
 * paternal copy at 2/3, an extra maternal copy at 1/3.
 *
 * WINDOW STATISTIC, NEVER PER MARKER. The measured per-marker spread is 0.0585 in bulk ESC
 * material, 0.1410 in trophectoderm and 0.3079 in a WGA blastomere. The band separation is 0.167.
 * That is 1.9 sd in bulk and 0.34 sd in a single cell, so a per-marker call is worthless at every
 * stage and hopeless at the one that matters most. Only a window median is used.
 *
 * MEDIAN, NOT MEAN. On WGA single cells the band separation compresses 36%, from 0.167 to 0.106,
 * but the band MEDIANS do not move (0.3336 / 0.5003 / 0.6669). The compression is in the tails,
 * so a median is stable where a mean is not. It also compresses symmetrically - |paternal| 0.1061
 * against |maternal| 0.1062 - so the DIRECTION is unbiased even where the magnitude shrinks.
 * That is why this reports direction and refuses to report a copy count.
 *
 * RECENTRING IS REQUIRED, NOT OPTIONAL. The centre is the sample's own genome-wide median share,
 * never the theoretical 0.5. Against a reconstructed parent the theoretical centre is biased
 * toward maternal, because a marker where the father is truly heterozygous cannot be represented
 * in an all-homozygous reference and is silently promoted into the informative set carrying a
 * maternal-looking allele share. Measured offsets: +0.013 at five products, +0.032 at four,
 * +0.073 to +0.077 at three. That last is 46-72% of the entire band separation, which is enough
 * to invert calls outright. Recentring on the sample's own median removes it (residual 0.0001 to
 * 0.0130). Confirmed on a real array in the audit: A8's own centre is 0.5993 rather than 0.5000,
 * and against the theoretical centre its UNTOUCHED euploid genome reads +0.0993, over the margin,
 * i.e. a false paternal gain. Recentring puts it at +0.0000. It is valid only while less than about half the genome carries a gain, which is stated
 * rather than checked because a sample past that point is not one this tool can reason about.
 */
import type { AB } from './informativity.ts'
import { twoPointPosterior, bandOf, type Band } from './originPosterior.ts'
import { mad, median as medianOf } from './intensityNull.ts'

/**
 * Informative markers a window needs before a direction is reported.
 *
 * WAS 400, LOWERED TO 200 ON MEASUREMENT. The original figure came from the smallest count at
 * which the worst array in a stage stays under 1% two-way error, taken as 100 in a bulk ESC line,
 * 50 to 100 in a trophectoderm biopsy and 400 in a WGA blastomere, with the strictest shipped
 * because the tool cannot see which stage it has been handed.
 *
 * The blastomere figure was re-measured against a BULK-GENOTYPED FATHER on five of his biparental
 * children, real WGA arrays, using the same criterion: every window here is euploid, so any
 * directional call is an error.
 *
 *     markers    worst array, two-way error
 *        100         0.0051
 *        200         0.0030
 *        400         0.0023
 *
 * Every count tested is already under the 1% bar, and 200 is under it by a factor of three.
 *
 * SPECIFICITY ALONE WOULD NOT JUSTIFY THIS, since a caller that never calls anything has no
 * error. Sensitivity was measured separately, on a real HapMap CEPH trio with a trisomy built
 * from the parents' own alleles and dropout drawn from a Markov model fitted to real WGA arrays
 * (a marker beside a dropped one is 2.18x more likely to drop). At 200 markers the direction is
 * recovered at 1.00 with half the alleles dropping; at 50 it falls to 0.87.
 *
 * So both halves are measured on the material the floor governs, and 200 is where they meet.
 *
 * The cost, restated at the new value: 200 informative markers span a median 23.8 Mb rather than
 * 47.6, which halves the smallest region that can carry a parental label. A focal gain still
 * cannot be annotated; an arm, and now a large sub-arm region, can.
 */
export const MIN_INFORMATIVE_DEFAULT = 200
export const MIN_INFORMATIVE_TROPHECTODERM = 100
export const MIN_INFORMATIVE_BULK = 100

/**
 * Distance from the sample's own centre before a direction is called.
 *
 * The band positions are measured. THIS THRESHOLD IS NOT: it is placed at a third of the measured
 * separation, which is a choice this project is making rather than a figure anyone reported. The
 * reasoning is a noise calculation. At the default 400 markers with the worst measured per-marker
 * spread of 0.3079, the standard error of a median is about 1.25 * 0.3079 / sqrt(400) = 0.019, so
 * 0.056 is roughly three of those. Anything closer to the centre than this is reported as unclear.
 *
 * Because it is derived rather than measured, it errs toward refusing. A run that wants the
 * measured error rates behind it should quote the marker counts above, which are measured, and
 * treat the band between the thresholds as uncalled rather than as weak evidence.
 */
export const SHARE_MARGIN = 0.056

/** Where the bands sit, measured, relative to the centre. Reported for context, not used as a cut. */
export const EXPECTED_SEPARATION = 0.167

export interface DosageMarker {
  chrom: string
  pos: number
  /** Fraction of signal at this marker attributable to the father's allele. */
  patShare: number
}

/**
 * The paternal allele share at one marker, or null where the marker says nothing.
 *
 * Informative means the two parents are homozygous for DIFFERENT alleles. Same-homozygote markers
 * carry no dosage information at all: both parents give the same allele and every copy number
 * reads identically. A heterozygous parent is likewise uninformative, and against a reconstructed
 * reference cannot occur, which is the bias `recentre` exists to absorb.
 */
export function paternalShare(father: AB, mother: AB, baf: number | null): number | null {
  if (baf === null || !Number.isFinite(baf)) return null
  if (father === 'AA' && mother === 'BB') return 1 - baf
  if (father === 'BB' && mother === 'AA') return baf
  return null
}

/**
 * With only one parent measured, a marker is informative where that parent is homozygous and the
 * sample is heterozygous, which means the other parent contributed the opposite allele. The share
 * is oriented the same way, so downstream code cannot tell the two paths apart by accident.
 *
 * DO NOT USE THIS. It is retained only so that callers referencing it fail loudly rather than
 * silently, and it is not wired into any shipped path.
 *
 * The defect is the conditioning, and it is total rather than a matter of degree. Selecting on the
 * SAMPLE'S OWN heterozygosity selects exactly the markers where both alleles survived, which forces
 * the share to 0.5 at every retained marker: simulated mean 0.5000 with zero markers at either
 * tail. There is no information left to read. It also re-derives the region from heterozygosity,
 * which is the failure that produced every artefact in the audit.
 *
 * The correct one-parent construction establishes heterozygosity from EVENT-FREE SIBLING CELLS of
 * the same embryo, never from the cell being scored. Sibling-established het keeps both tails
 * (15.3% of markers at 1.0 and 15.3% at 0.0 in the same simulation) because the affected cell's own
 * calls never enter the selection.
 *
 * See audit/CONSULT-parental-origin.txt section 1A, including the false-het fraction phi that
 * governs the bias and the sibling rule that controls it.
 */
export function paternalShareOneParent(father: AB, sample: AB, baf: number | null): number | null {
  if (baf === null || !Number.isFinite(baf)) return null
  if (sample !== 'AB') return null
  if (father === 'AA') return 1 - baf
  if (father === 'BB') return baf
  return null
}

const median = (xs: readonly number[]): number => {
  if (!xs.length) return NaN
  const s = [...xs].sort((a, b) => a - b)
  const h = s.length >> 1
  return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2
}

/**
 * The sample's own centre: the median paternal share over every informative marker it has.
 *
 * CORRECT FOR THE TWO-PARENT CLASS, WRONG FOR THE ONE-PARENT CLASS, and the difference is not a
 * matter of degree. Read the caveat before reusing this.
 *
 * With BOTH parents the informative markers are obligate-het, the disomic expectation is 0.5, and
 * the two loss hypotheses sit at 0 and 1. The distribution is centred and its median is a sensible
 * reference point. Using the theoretical 0.5 instead is what inverts calls against a reconstructed
 * parent, which is why this exists.
 *
 * With ONE parent the informative markers are that parent's homozygous sites, and the disomic
 * expectation is 1 - q/2 where q is the alternate-allele frequency. The two loss hypotheses are
 * still SYMMETRIC about it, at 1 - q and 1. But for every q < 0.5 the majority of markers sit
 * exactly at 1.0, so the MEDIAN IS THE CEILING and the upward headroom is 0.000 by construction:
 *
 *     q      markers at share 1.0    median   mean     headroom above median / above mean
 *     0.10   0.916                   1.000    0.9502   +0.0000 / +0.0498
 *     0.25   0.788                   1.000    0.8749   +0.0000 / +0.1251
 *     0.40   0.661                   1.000    0.7995   +0.0000 / +0.2005
 *
 * A statistic whose centre is a boundary point of its own support can only move one way. That is
 * the generator of the one-directional results recorded in the audit, and it is a centring defect
 * rather than a dropout effect: allele dropout sends AB to either homozygote with equal
 * probability, so the expected share is invariant in it (0.8754 / 0.8748 / 0.8748 / 0.8754 at ADO
 * 0.000 / 0.199 / 0.308 / 0.600 against an analytic 0.8750). Dropout inflates variance only.
 *
 * So for a one-parent informative set, centre on the model expectation 1 - q/2, or drop the share
 * scalar for a per-marker likelihood, which needs no centring at all. Do not pass a one-parent
 * marker set to this function.
 */
export const recentre = (all: readonly DosageMarker[]): number =>
  median(all.map((m) => m.patShare))

/**
 * The disomic expectation for a ONE-PARENT informative set, from the alternate-allele frequency.
 *
 * `q` is estimated from the marker set itself rather than assumed: at a reference-homozygous
 * marker the sample carries the alternate allele at rate q/2 under disomy, so q = 2(1 - mean share).
 * The mean is used deliberately, since the median is the ceiling here.
 */
export const oneParentCentre = (all: readonly DosageMarker[]): number => {
  if (!all.length) return NaN
  const mean = all.reduce((a, m) => a + m.patShare, 0) / all.length
  return mean
}

export type GainOrigin = 'paternal' | 'maternal' | 'unclear'

export interface OriginCall {
  origin: GainOrigin
  /** Median paternal share inside the region. */
  share: number
  /** The sample's own genome-wide median, which the share is judged against. */
  centre: number
  /** share - centre. Positive is paternal, negative is maternal. */
  deviation: number
  informative: number
  /**
   * Calibrated probability that the named parent is right, and the band it lands in.
   *
   * This channel used to emit a bare direction at a hard 0.056 margin, so the strongest evidence
   * the tool has printed no number while a weak dosage call printed one. A reader comparing the
   * two rows was misled by the formatting alone. NaN where no direction was called.
   */
  confidence: number
  band: Band
  /**
   * Always true on this channel, and it is not a placeholder.
   *
   * The number below rests on a model of the share's noise, not on a measurement of how often the
   * model is right. The dosage channel's bands were earned from 140,000 injections carrying real
   * array noise; this channel has no equivalent, because no clean parent-child pair was available
   * to build one from. Until one is, the band is capped at B and this flag travels with every call
   * so a reader is never shown a modelled number dressed as a measured one.
   */
  uncalibrated?: boolean
  why: string
}

/**
 * Spread of the per-marker paternal share, which is what turns a deviation into a confidence.
 *
 * Median absolute deviation rather than a standard deviation, scaled to be comparable, because the
 * share distribution carries a contaminating homozygous tail that a plain SD would chase. That tail
 * is the same one that destroyed the old band-threshold estimator in the dosage channel.
 */
export const shareSpread = (region: readonly DosageMarker[]): number => {
  const xs = region.map((m) => m.patShare).filter((x): x is number => Number.isFinite(x))
  if (xs.length < 3) return NaN
  const med = median(xs)
  return 1.4826 * median(xs.map((x) => Math.abs(x - med)))
}

/**
 * Which parent supplied the extra copy in one region.
 *
 * `all` is every informative marker on the array, used only for the centre. `region` is the subset
 * inside the event. Both are the caller's job to assemble, so this stays independent of how
 * segments are found.
 */
export function callGainOrigin(
  region: readonly DosageMarker[],
  centre: number,
  minInformative = MIN_INFORMATIVE_DEFAULT,
  margin = SHARE_MARGIN,
): OriginCall {
  const informative = region.length
  const share = median(region.map((m) => m.patShare))
  const deviation = share - centre
  const base = { share, centre, deviation, informative }

  if (!Number.isFinite(centre)) {
    return {
      ...base, origin: 'unclear', confidence: NaN, band: 'F',
      why: 'the sample has no informative markers to centre on',
    }
  }
  if (informative < minInformative) {
    return {
      ...base,
      origin: 'unclear', confidence: NaN, band: 'F',
      why: `${informative} informative markers is under the ${minInformative} needed for a `
        + 'direction at this stage',
    }
  }
  if (Math.abs(deviation) < margin) {
    return {
      ...base,
      origin: 'unclear', confidence: NaN, band: 'F',
      why: `the share sits ${deviation >= 0 ? '+' : ''}${deviation.toFixed(3)} from this `
        + `sample's own centre, inside the ${margin} band that is not called either way`,
    }
  }
  // The confidence is in the DIRECTION, so a loss reading inherits it unchanged: the same
  // deviation supports "paternal gained" exactly as strongly as it supports "maternal lost".
  //
  // THE STANDARD ERROR IS FLOORED, AND THE FLOOR IS THE WHOLE STORY. Sampling error alone,
  // spread/sqrt(n), is around 0.002 over the marker counts this channel sees, which would make any
  // deviation clearing the 0.056 margin read as a certainty. The dosage channel has been here
  // before: within-array drift does not average down with more markers, and omitting it is how a z
  // of 15 arrives on a chromosome independently verified normal. So the floor is SHARE_MARGIN, the
  // project's own long-standing estimate of the scale at which this statistic stops being
  // trustworthy, and the sampling term is added in quadrature beneath it.
  //
  // AND THE RESULT IS CAPPED AT BAND B, because it is a MODEL rather than a measurement. The
  // dosage posterior earned its bands from 140,000 injections carrying real array noise
  // (audit/calibration/FINDINGS.txt). This one has had no equivalent test: the corpus to hand
  // holds no clean parent-child pair to build one from, and the candidates in it are not what
  // their filenames suggest, one nominal parent reading 55% heterozygous against a diploid's
  // 16.8%. An unvalidated model must not be allowed to print "very confident", so it cannot,
  // and every call it makes reports itself as uncalibrated.
  const spread = shareSpread(region)
  const sampling = Number.isFinite(spread) ? spread / Math.sqrt(informative) : NaN
  const se = Number.isFinite(sampling)
    ? Math.sqrt(sampling * sampling + margin * margin) : NaN
  const pPaternal = twoPointPosterior(deviation, se, EXPECTED_SEPARATION)
  const confidence = Number.isFinite(pPaternal) ? Math.max(pPaternal, 1 - pPaternal) : NaN
  const raw = bandOf(confidence)
  return {
    ...base,
    origin: deviation > 0 ? 'paternal' : 'maternal',
    confidence,
    band: raw === 'A' ? 'B' : raw,
    uncalibrated: true,
    why: `${deviation > 0 ? 'paternal' : 'maternal'} alleles are over-represented by `
      + `${Math.abs(deviation).toFixed(3)} against this sample's own centre of `
      + `${centre.toFixed(3)}, over the ${margin} margin`
      + (Number.isFinite(confidence)
        ? `, at a confidence of ${confidence.toFixed(4)} against a per-marker spread of `
          + `${spread.toFixed(3)} over ${informative} markers`
        : ''),
  }
}

/**
 * Which parent's copy is MISSING across a region, in a biparental cell.
 *
 * THE DIRECTION IS THE OPPOSITE OF A GAIN, and this function exists so that inversion is written
 * down once rather than rediscovered. An extra paternal copy raises the paternal share; a missing
 * paternal copy removes it. So the same deviation means different parents depending on which
 * event it belongs to:
 *
 *     paternal GAIN   share rises      deviation > 0
 *     paternal LOSS   share falls      deviation < 0
 *
 * Passing a loss to callGainOrigin therefore names the wrong parent every time.
 *
 * A loss is the larger signal of the two. One extra of two copies moves the share to 0.667, while
 * losing one of two removes that parent's alleles entirely and moves it to 0 or 1. The margin is
 * kept as a floor rather than raised to suit, because whole-genome amplification compresses the
 * observed value and a partial or mosaic loss lands between the two.
 */
export function callLossOrigin(
  region: readonly DosageMarker[],
  centre: number,
  minInformative = MIN_INFORMATIVE_DEFAULT,
  margin = SHARE_MARGIN,
): OriginCall {
  const call = callGainOrigin(region, centre, minInformative, margin)
  if (call.origin === 'unclear') return call
  // Same measurement, opposite reading: the parent under-represented here is the one lost.
  const lost = call.deviation < 0 ? 'paternal' : 'maternal'
  return {
    ...call,
    origin: lost,
    why: `${lost} alleles are under-represented by ${Math.abs(call.deviation).toFixed(3)} against `
      + `this sample's own centre of ${centre.toFixed(3)}, over the ${margin} margin, which is `
      + `the ${lost} copy missing rather than the other parent's gained`,
  }
}

// --- the uniparental case ----------------------------------------------------------------------

/**
 * Which homologue an extra copy carries, in a cell that has only one parent's genome.
 *
 * ONE OF THE TWO ANSWERS IS UNREACHABLE AND THE TOOL SAYS SO. A haploid genome cannot be
 * heterozygous, so every heterozygous call in one is error. If the extra copy is the parent's
 * OTHER homologue, the two differ wherever the parent is heterozygous and the region reads
 * heterozygous far above that error rate: measured AUC 1.000 against the real amplification
 * background from 50 markers upward, holding even against degraded arrays running 24-45% spurious
 * heterozygosity. If the extra copy is a duplicate of the SAME homologue it is bit-identical to
 * one copy at every marker, and no genotype channel can see it. Measured AUC 0.037 to 0.119,
 * which is at chance, and intensity does not rescue it either: 0.503 to 0.561 for one copy
 * against two, on a mean log2R difference of 0.005 to 0.013 inside a spread of 0.07 to 0.16.
 *
 * So the honest output is: report the meiotic case, and where the evidence is absent report that
 * an identical duplicate cannot be distinguished from the normal single copy. NOT "no gain".
 */
export const HETERO_MULTIPLE = 2.0

/**
 * And an ABSOLUTE floor, because a ratio alone is not enough.
 *
 * A region carrying the parent's other homologue reads heterozygous wherever that parent is
 * heterozygous, which is roughly the parent's own heterozygosity: about 17% on this material. A
 * ratio test alone fires on 7.4% against a 3.4% background, which is two clean chromosomes'
 * amplification noise differing from each other, not a second homologue. Measured on the audit
 * set, requiring both a doubling AND this floor removes every such call while keeping the ones
 * that read at 55-61%.
 */
export const HETERO_ABSOLUTE_MIN = 0.10

/**
 * The heterozygosity a region is judged against: the median over the OTHER chromosomes.
 *
 * NOT the genome-wide rate. A cell carrying several gained chromosomes raises the genome-wide
 * heterozygosity with the very events being tested, so each one is measured against a bar its
 * siblings have already lifted, and a real gain reads as ordinary. This is the same self-null
 * failure the segment scan solved with an external per-chromosome null, and it is solved the same
 * way here: exclude the chromosome under test, take a median so a handful of affected chromosomes
 * cannot drag it, and floor it so a spotless genome does not divide by nothing.
 *
 * Measured consequence of getting it wrong, on a real array carrying five whole-chromosome gains:
 * genome-wide heterozygosity 69.2% against per-chromosome values of 40-51%, so every one of the
 * five reads BELOW its own background and all five refuse.
 */
export const HET_BACKGROUND_FLOOR = 0.02

export function externalHetBackground(
  byChrom: ReadonlyMap<string, { informative: number; het: number }>,
  exclude: string,
  minInformative = 200,
): number {
  const rates: number[] = []
  for (const [c, h] of byChrom) {
    if (c === exclude) continue
    if (h.informative < minInformative) continue
    rates.push(h.het / h.informative)
  }
  if (!rates.length) return NaN
  return Math.max(HET_BACKGROUND_FLOOR, median(rates))
}

export interface HomologueCall {
  verdict: 'other homologue' | 'indistinguishable'
  regionHet: number
  backgroundHet: number
  informative: number
  why: string
}

export function callHomologue(
  regionHet: number,
  regionMarkers: number,
  backgroundHet: number,
  minMarkers = 50,
): HomologueCall {
  const base = { regionHet, backgroundHet, informative: regionMarkers }
  if (regionMarkers < minMarkers || !Number.isFinite(backgroundHet)) {
    return {
      ...base,
      verdict: 'indistinguishable',
      why: `${regionMarkers} markers is under the ${minMarkers} this needs`,
    }
  }
  // Order matters: each refusal has to give ITS OWN reason. A region sitting at the background
  // is the iso-duplication case and needs the bit-identical explanation; a region that doubles
  // the background but is still only a few percent is amplification noise and needs a different
  // one. Collapsing them into one message loses the finding that matters.
  if (regionHet >= HETERO_MULTIPLE * backgroundHet && regionHet < HETERO_ABSOLUTE_MIN) {
    return {
      ...base,
      verdict: 'indistinguishable',
      why: `${(regionHet * 100).toFixed(1)}% heterozygous doubles the ${
        (backgroundHet * 100).toFixed(1)}% background but is under the ${
        (HETERO_ABSOLUTE_MIN * 100).toFixed(0)}% a second homologue would produce. A region `
        + `carrying both homologues reads near the parent's own heterozygosity; two clean `
        + 'chromosomes whose amplification noise differs by a factor of two do not',
    }
  }
  if (regionHet >= HETERO_MULTIPLE * backgroundHet) {
    return {
      ...base,
      verdict: 'other homologue',
      why: `${(regionHet * 100).toFixed(1)}% heterozygous against a ${
        (backgroundHet * 100).toFixed(1)}% background for this array, which a single copy of one `
        + 'homologue cannot produce: the extra copy is the parent\'s other homologue, a meiotic '
        + 'error',
    }
  }
  return {
    ...base,
    verdict: 'indistinguishable',
    why: `${(regionHet * 100).toFixed(1)}% heterozygous against a ${
      (backgroundHet * 100).toFixed(1)}% background, which is what one copy looks like. An extra `
      + 'copy of the SAME homologue is bit-identical to one copy and cannot be detected in any '
      + 'channel this array carries, so this is not evidence of no gain',
  }
}

// --- whole-chromosome gain, from the alleles ---------------------------------------------------

/**
 * Which parent supplied an extra whole chromosome, read from the alleles at the markers where the
 * parents are opposite homozygotes, with intensity as a gate and never as the detector.
 *
 * WHY THIS EXISTS BESIDE callGainOrigin. That function names a parent for a gain somebody ELSE
 * detected, from call rate or intensity, and centres on the genome-wide median share. On an
 * amplified single cell neither half holds. The detectors rarely fire: 0 gain rows on the six
 * GSE148488 cells carrying hand-called gains, three of them on arrays the stage inference refuses.
 * And the genome-wide centre is the wrong reference in a cell with several events: the median
 * share of five of the 31 two-parent cells sits between 0.15 and 0.76, against 0.44 to 0.54 when
 * taken over the chromosomes the alleles call biparental.
 *
 * THE STATISTIC. The median paternal share over a chromosome's informative markers, judged against
 * the OTHER biparental chromosomes of the same cell: their median is the centre, their
 * median-absolute-deviation the scale, so every cell is its own null and a noisy cell earns a wide
 * scale and refuses rather than over-calls. The scale has a floor, because a very clean cell would
 * otherwise turn a 0.01 shift into a certainty.
 *
 * AN EXTRA PATERNAL COPY READS +, AN EXTRA MATERNAL ONE -, and a partial loss of the other parent
 * reads the SAME WAY: one parent over-represented. Allele ratios alone cannot separate them, since
 * an extra maternal copy in every cell and a paternal copy lost in half of them both put the
 * maternal share at two thirds. Total copy number separates them, three against one and a half,
 * which is why intensity is a gate here: a shift that clears the bar with intensity BELOW the
 * cell's own level is reported as an imbalance and names nobody.
 *
 * WHAT IS MEASURED AND WHAT IS NOT. On the 31 two-parent cells of GSE148488, the hand karyotype of
 * Zuccaro 2020 lists eight whole-chromosome or arm gains with a named parent. The sign of the
 * shift names the right parent on all of them. The whole-chromosome statistic clears the bar on
 * the one clean whole-chromosome trisomy (chromosome 15, maternal) and on the large-arm
 * duplication that dominates chromosome 6; the others are arm-sized and dilute over the
 * chromosome, or sit in cells too chaotic to scale. The same cells carry chromosomes that clear the
 * bar and that the hand karyotype lists as normal; sibling blastomeres share several, so some may
 * be real gains the hand karyotype did not tabulate, and the hand calls cannot say which. So this
 * channel reports a candidate with its numbers, capped at band C, and never a certainty.
 *
 * WHAT IT DOES NOT DO. It does not say meiotic or mitotic. It does not resolve arms. It needs two
 * bulk parental arrays and a genome that reads both parents on enough chromosomes, the same
 * conditions as the allele count.
 *
 * AND IT READS A FULL EXTRA COPY ONLY, measured rather than argued. Gains built into real cells at
 * a cell fraction of 0.3 and 0.5 are detected 0 times in 120: at half the share should move +0.10
 * and the measured median movement is +0.012, because dropout has already spread these markers
 * toward 0 and 1 and scaling one allele by 1.5 barely moves a marker reading 0.95. A full copy,
 * which scales by 2, does move it. So a mosaic gain is silent here, and that is the commoner
 * finding in an embryo biopsy.
 *
 * THE INTENSITY GATE IS AN ABSOLUTE FLOOR AND THE INTENSITY IT GATES IS NOT. On 60 partial losses
 * built at half, 8 keep an intensity over the floor because the chromosome already read high in
 * the untouched cell, and 3 of those also clear the share bar and are reported as an extra copy of
 * the other parent. That is the error this gate exists to prevent, at 3 in 60 on a construction
 * that does not model the heterozygosity a real partial loss also removes. See
 * audit/gain-origin-mosaic.txt. A row from this channel is a candidate for that reason as much as
 * for any other.
 */
export const ALLELE_GAIN_Z = 6
/** The smallest shift that is a gain however wide the cell's own scale is. A third of the 0.167
 *  an extra copy moves a clean median, about half of what WGA leaves of it. */
export const ALLELE_GAIN_SHIFT = 0.08
export const ALLELE_GAIN_MIN_INFORMATIVE = 200
/** Other biparental chromosomes a cell must have before it can be its own null. */
export const ALLELE_GAIN_MIN_REFERENCE = 6
/** Floor on the cell's own scale. The cleanest cells measure 0.005 and would otherwise read a
 *  shift of 0.03 as six sigma. */
export const ALLELE_GAIN_SCALE_FLOOR = 0.02
/** The same floor for the intensity spread. The cleanest cell measured here is 0.042, so this binds
 *  on nothing observed and exists so a cell with an implausibly tight intensity cannot manufacture a
 *  large z from a small shift. */
export const ALLELE_GAIN_LRR_SCALE_FLOOR = 0.04
/**
 * How far above the cell's own intensity level an extra copy must read, in units of that cell's own
 * intensity spread.
 *
 * AN ABSOLUTE FLOOR WAS WRONG AND THE REASON IS MEASURED. The rule this replaces required +0.15
 * log2 whatever the cell. But a chromosome's log2R is far noisier than its allele share: over the
 * 30 scorable GSE148488 cells the median absolute deviation of chromosome intensity WITHIN one
 * cell runs 0.042 to 0.246, median 0.092, which is 4.9 times the same cell's share spread. A fixed
 * floor is therefore several spreads in a clean cell and under one spread in a noisy one. Both
 * cells that produced a false gain on a constructed partial loss sit at the noisy end, 0.146 and
 * 0.093, and their chromosome 3 already read +0.29 and +0.46 before anything was built into it, so
 * removing half a copy still left it over +0.15 and it was reported as an extra copy of the other
 * parent (audit/gain-origin-mosaic.txt).
 *
 * SO THE BAR IS THE CELL'S OWN SPREAD, strict where the intensity is unreliable and sensitive where
 * it is not. Two spreads is the conventional one-sided bar and is NOT fitted to the test that found
 * the defect: in the median cell it is +0.18, close to the floor it replaces, while in the noisiest
 * it is +0.29 and in the cleanest +0.08. It sits between the two expectations it has to separate,
 * a full extra copy at +0.30 absolute, which is 3.3 spreads in a median cell, and a half loss at
 * -0.21, which is -2.3.
 */
export const ALLELE_GAIN_LRR_Z = 2

export interface ChromShare {
  chrom: string
  /** Median paternal share over this chromosome's informative markers. */
  share: number
  informative: number
  /** The heterozygous share at these markers read both parents, so the chromosome can be a
   *  reference for the others and can itself carry a gain. */
  biparental: boolean
  /** Median log2R against the cell's own null centre. */
  lrrShift: number
}

export interface AlleleGain {
  chrom: string
  /** 'gain' names a parent; 'imbalance' is a shift with intensity below the cell's level, which
   *  is a partial loss of the other parent and not an extra copy; 'none' cleared nothing;
   *  'unscored' had too little to be judged. */
  verdict: 'gain' | 'imbalance' | 'none' | 'unscored'
  /** The parent who supplied the extra copy. Only on 'gain'. */
  parent?: 'paternal' | 'maternal'
  shift: number
  z: number
  scale: number
  centre: number
  informative: number
  lrrShift: number
  /** The intensity against the same reference chromosomes, in units of their own spread. */
  lrrZ?: number
  /** The pair read as a likelihood, reported beside the verdict. See callPair. */
  pair?: PairCall
  why: string
}

/** One clause naming what the pair likelihood makes of the same two numbers, where it disagrees
 *  with the verdict or adds a cell fraction the thresholds cannot see. Evidence, never a verdict. */
const pairNote = (p: PairCall | undefined): string => {
  if (!p || p.state === 'disomy') return ''
  const what = p.state === 'unexplained'
    ? 'that this chromosome\'s intensity is not measuring its copy number'
    : `${p.state.replace('-', ' of the ')} copy`
      + (p.fraction < 1 ? ` in about ${(p.fraction * 100).toFixed(0)}% of the cells` : '')
  return `. Weighing the two numbers together as a likelihood prefers ${what}, at ${p.posterior.toFixed(2)} `
    + `against ${p.runnerUp.replace('-', ' of the ')}`
    + (p.tilt
      ? `. The direction of that call rests on intensity, which carries a per-chromosome baseline a `
        + `single cell cannot measure; the statement that the ${p.tilt} side is the over-represented `
        + `one rests on the alleles and holds at ${p.tiltPosterior.toFixed(2)}`
      : '')
    + `; that model's calling threshold has not been set on held-out data, so it is reported and `
    + 'does not decide this row'
}

export function scanAlleleGains(chroms: readonly ChromShare[]): AlleleGain[] {
  const out: AlleleGain[] = []
  for (const c of chroms) {
    const base = {
      chrom: c.chrom, shift: NaN, z: NaN, scale: NaN, centre: NaN,
      informative: c.informative, lrrShift: c.lrrShift,
    }
    if (c.informative < ALLELE_GAIN_MIN_INFORMATIVE || !Number.isFinite(c.share)) {
      out.push({ ...base, verdict: 'unscored',
        why: `${c.informative} informative markers is under the ${ALLELE_GAIN_MIN_INFORMATIVE} needed` })
      continue
    }
    if (!c.biparental) {
      out.push({ ...base, verdict: 'unscored',
        why: 'this chromosome does not read both parents, so it is a loss question and not a gain one' })
      continue
    }
    const ref = chroms.filter((o) => o.chrom !== c.chrom && o.biparental
      && o.informative >= ALLELE_GAIN_MIN_INFORMATIVE && Number.isFinite(o.share))
    if (ref.length < ALLELE_GAIN_MIN_REFERENCE) {
      out.push({ ...base, verdict: 'unscored',
        why: `${ref.length} other chromosomes read both parents, under the ${ALLELE_GAIN_MIN_REFERENCE} `
          + 'this cell needs to serve as its own reference' })
      continue
    }
    const shares = ref.map((o) => o.share)
    const centre = medianOf(shares)
    const scale = Math.max(mad(shares), ALLELE_GAIN_SCALE_FLOOR)
    const shift = c.share - centre
    const z = shift / scale
    // THE INTENSITY ON THE SAME TERMS AS THE SHARE: this chromosome against the same reference
    // chromosomes of the same cell, in units of their own spread. See ALLELE_GAIN_LRR_Z.
    const lrrs = ref.map((o) => o.lrrShift).filter((x) => Number.isFinite(x))
    const lrrCentre = lrrs.length >= ALLELE_GAIN_MIN_REFERENCE ? medianOf(lrrs) : NaN
    const lrrScale = lrrs.length >= ALLELE_GAIN_MIN_REFERENCE
      ? Math.max(mad(lrrs), ALLELE_GAIN_LRR_SCALE_FLOOR) : NaN
    const lrrZ = (c.lrrShift - lrrCentre) / lrrScale
    // THE PAIR, REPORTED BESIDE THE VERDICT AND NOT YET DECIDING IT. The thresholds below are
    // measured on constructed events; this is a model whose calling threshold has not been set on
    // a held-out half, so it travels with every row as evidence a reader can weigh and is not
    // allowed to overturn a verdict on its own. See callPair.
    const pair = callPair(shift, c.lrrShift, scale, lrrScale) ?? undefined
    const m = { ...base, shift, z, scale, centre, lrrZ, pair }
    if (Math.abs(z) < ALLELE_GAIN_Z || Math.abs(shift) < ALLELE_GAIN_SHIFT) {
      out.push({ ...m, verdict: 'none',
        why: `the share sits ${shift >= 0 ? '+' : ''}${shift.toFixed(3)} from the cell's own centre, `
          + `${z.toFixed(1)} of its own spreads` })
      continue
    }
    const over = shift > 0 ? 'paternal' : 'maternal'
    const other = shift > 0 ? 'maternal' : 'paternal'
    const nums = `${over} alleles are over-represented by ${Math.abs(shift).toFixed(3)} against the `
      + `${ref.length} other chromosomes of this cell (centre ${centre.toFixed(3)}, ${Math.abs(z).toFixed(1)} `
      + `of their own spreads of ${scale.toFixed(3)})`
    if (Number.isFinite(lrrZ) && lrrZ >= ALLELE_GAIN_LRR_Z) {
      out.push({ ...m, verdict: 'gain', parent: over,
        why: `${nums}, and the intensity reads ${c.lrrShift >= 0 ? '+' : ''}${c.lrrShift.toFixed(2)} log2 `
          + `from the cell's own level, ${lrrZ.toFixed(1)} of that cell's own intensity spreads of `
          + `${lrrScale.toFixed(2)}, so the ${over} parent supplied an extra copy. Candidate, not `
          + 'calibrated: a partial loss of the other parent in a fraction of cells cannot be excluded '
          + `from the alleles alone${pairNote(pair)}` })
    } else {
      const lrr = Number.isFinite(c.lrrShift) ? c.lrrShift : NaN
      out.push({ ...m, verdict: 'imbalance',
        why: !Number.isFinite(lrr) || !Number.isFinite(lrrZ)
          ? `${nums}, but this file carries no usable intensity, so an extra copy cannot be told from a `
            + `partial loss of the ${other} copy and no parent is named`
          : lrrZ < 0
          ? `${nums}, but the intensity reads ${lrr.toFixed(2)} log2 from the cell's own level, `
            + `${lrrZ.toFixed(1)} of that cell's own intensity spreads, below it. An extra copy would read `
            + `above. A partial loss of the ${other} copy fits, and no parent is named${pairNote(pair)}`
          : `${nums}, but the intensity reads ${lrr >= 0 ? '+' : ''}${lrr.toFixed(2)} log2, only `
            + `${lrrZ.toFixed(1)} of that cell's own intensity spreads of ${lrrScale.toFixed(2)}, under the `
            + `${ALLELE_GAIN_LRR_Z} an extra copy must clear. Neither an extra copy nor a loss is `
            + `supported, and no parent is named${pairNote(pair)}` })
    }
  }
  return out
}

// --- the pair, as a likelihood --------------------------------------------------------------

/**
 * Which of five states a chromosome is in, from the allele share and the intensity TOGETHER.
 *
 * WHY A PAIR AND NOT TWO THRESHOLDS. The share says a parent is over-represented and cannot say
 * why: an extra copy of one parent and a partial loss of the other move it the same way and by the
 * same amount. The intensity separates them by its SIGN. Deciding on each number alone throws away
 * what their RELATIONSHIP carries, and that relationship is the stable part. Measured on 300 gains
 * and partial losses built into real cells, the ratio of intensity to share shift is
 *
 *     gain    +2.04 at a cell fraction of 0.3, +1.96 at 0.5, +2.01 at 1.0
 *     loss    -0.75 at 0.3, -0.97 at 0.5
 *
 * so it barely moves with the fraction, which is what makes it usable when the fraction is unknown.
 * Both failures of the threshold pair come from ignoring it. A chromosome whose intensity sits a
 * little above its cell is called an extra copy however small its share shift, until the bar is
 * raised, and then a real trisomy whose intensity is only moderately raised is refused although its
 * intensity is exactly what its share predicts.
 *
 * THE STATES AND THEIR PREDICTIONS. Writing d for the paternal share shift and l for the intensity
 * against the cell's own level, at a cell fraction f:
 *
 *     disomy                 d = 0                     l = 0
 *     extra paternal copy    d = +f / (2(2+f))         l = k log2((2+f)/2)
 *     extra maternal copy    d = -f / (2(2+f))         l = k log2((2+f)/2)
 *     paternal copy lost     d = -f / (2(2-f))         l = k log2((2-f)/2)
 *     maternal copy lost     d = +f / (2(2-f))         l = k log2((2-f)/2)
 *
 * The share arithmetic is exact: a full extra copy puts the share at two thirds, which is a shift
 * of one sixth, and that is what a full build measures (0.166 against 0.1667). The intensity
 * carries a chemistry factor k, because one extra copy of two is log2(3/2) = 0.585 in theory and
 * reads 0.5547 on one chip of GSE19247 and 0.3420 on another, 95 and 58 percent. k is not known for
 * an unseen array, so it is marginalised over that measured range rather than assumed.
 *
 * AND A SIXTH STATE, WHICH IS WHAT MAKES THE OTHER FIVE SAFE. A chromosome whose intensity matches
 * nothing its share predicts is not evidence for whichever of the five fits least badly. That case
 * is real and common: a half loss of one parent on a chromosome whose intensity already read high
 * leaves the intensity positive but far below what a gain of that share shift requires, and the
 * five-state model then prefers the gain by default. `unexplained` lets the share move freely and
 * draws the intensity from the cell's own intensity distribution rather than from any copy number,
 * which is the hypothesis that this chromosome's intensity is not measuring its copy number. It
 * wins exactly when no copy-number state fits, and it loses to a state whose prediction the
 * intensity actually matches. On the two chromosomes that broke the threshold pair, the real
 * maternal trisomy's intensity sits 0.2 of a spread from what its share predicts for a gain and the
 * constructed half loss sits 1.1 away, which is the whole of the separation.
 *
 * WHICH OF THE TWO STATEMENTS TO TRUST. An extra copy from one parent and a missing copy from the
 * other are the same tilt and move the share the same way, so only intensity tells them apart. Over
 * 300 constructed builds the model names the state on 267 and the tilt is right on 262 of those,
 * including 94 of 98 at a third of the cells, where the shipped thresholds detect nothing at all.
 * Every one of the 22 calls that named the wrong one of gain and loss was on a chromosome whose own
 * baseline intensity already sat between +0.18 and +0.49, against a range of -0.05 to +0.10 on every
 * arm that never erred. A per-chromosome baseline offset of that size is not separable from a real
 * copy-number change by intensity alone within one cell, and the cell's intensity spread does not
 * reveal it, which is why a wrong call can carry a posterior of 1.00. The tilt comes from the
 * alleles and is immune to it, so it is reported as its own number rather than folded into the
 * state.
 *
 * WHAT THIS DOES NOT DO. It is a model, not a measurement of how often the model is right, and the
 * posterior it returns is therefore not a calibrated probability. It assumes the two numbers are
 * independent given the state, which they are not exactly. And it cannot see a third copy that
 * duplicates a homologue already present, which moves neither number the way an extra parental copy
 * does.
 */
export type PairState = 'disomy' | 'gain-paternal' | 'gain-maternal' | 'loss-paternal' | 'loss-maternal'
  | 'unexplained'

/** The chemistry factor on one copy of intensity, over the range measured across two chips. */
export const COMPRESSION_GRID = [0.55, 0.7, 0.85, 1.0]
/** Cell fractions considered. A whole-cell event is 1; below 0.1 nothing here can see anything. */
export const FRACTION_GRID = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

/**
 * Which parent ends up over-represented. An extra copy from one parent and a missing copy from the
 * other are the same tilt, and they move the allele share the same way: only intensity separates
 * them. Intensity carries a per-chromosome baseline that a single cell cannot measure, so the tilt
 * is the more reliable of the two statements and is reported separately.
 */
const TILT: Record<PairState, 'paternal' | 'maternal' | null> = {
  'disomy': null, 'unexplained': null,
  'gain-paternal': 'paternal', 'loss-maternal': 'paternal',
  'gain-maternal': 'maternal', 'loss-paternal': 'maternal',
}

export interface PairCall {
  state: PairState
  /** Posterior of the winning state under the model above. NOT calibrated: see the docstring. */
  posterior: number
  /** The runner-up, so a reader can see what it was weighed against. */
  runnerUp: PairState
  /** Cell fraction of the best fit, and the compression it assumed. */
  fraction: number
  compression: number
  /** The over-represented parent, and the posterior summed over both states that agree on it. */
  tilt: 'paternal' | 'maternal' | null
  tiltPosterior: number
}

const NORMAL_LOG = (x: number, mu: number, sd: number): number => {
  const z = (x - mu) / sd
  return -0.5 * z * z - Math.log(sd)
}

export function callPair(
  shift: number, lrr: number, shareScale: number, lrrScale: number,
): PairCall | null {
  if (![shift, lrr, shareScale, lrrScale].every((x) => Number.isFinite(x))) return null
  if (shareScale <= 0 || lrrScale <= 0) return null
  const best = new Map<PairState, { lp: number; f: number; k: number }>()
  const consider = (s: PairState, d: number, l: number, f: number, k: number) => {
    const lp = NORMAL_LOG(shift, d, shareScale) + NORMAL_LOG(lrr, l, lrrScale)
    const cur = best.get(s)
    if (!cur || lp > cur.lp) best.set(s, { lp, f, k })
  }
  consider('disomy', 0, 0, 0, 0)
  // The share is free here, so it costs nothing to explain, and the intensity is drawn from the
  // cell's own spread about its own centre rather than from a copy number.
  consider('unexplained', shift, 0, 0, 0)
  for (const f of FRACTION_GRID) {
    const dGain = f / (2 * (2 + f))
    const dLoss = f / (2 * (2 - f))
    for (const k of COMPRESSION_GRID) {
      const lGain = k * Math.log2((2 + f) / 2)
      const lLoss = k * Math.log2((2 - f) / 2)
      consider('gain-paternal', +dGain, lGain, f, k)
      consider('gain-maternal', -dGain, lGain, f, k)
      consider('loss-paternal', -dLoss, lLoss, f, k)
      consider('loss-maternal', +dLoss, lLoss, f, k)
    }
  }
  const rows = [...best.entries()]
  const top = Math.max(...rows.map(([, v]) => v.lp))
  const w = rows.map(([s, v]) => [s, Math.exp(v.lp - top)] as [PairState, number])
  const sum = w.reduce((a, [, x]) => a + x, 0)
  const sorted = w.map(([s, x]) => [s, x / sum] as [PairState, number]).sort((a, b) => b[1] - a[1])
  const win = best.get(sorted[0][0])!
  const tilt = TILT[sorted[0][0]]
  return {
    state: sorted[0][0], posterior: sorted[0][1], runnerUp: sorted[1][0],
    fraction: win.f, compression: win.k,
    tilt,
    tiltPosterior: tilt ? sorted.reduce((a, [s, p]) => a + (TILT[s] === tilt ? p : 0), 0) : 0,
  }
}

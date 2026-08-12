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
 * 0.0130). It is valid only while less than about half the genome carries a gain, which is stated
 * rather than checked because a sample past that point is not one this tool can reason about.
 */
import type { AB } from './informativity.ts'

/**
 * Informative markers a window needs before a direction is reported.
 *
 * Stage-dependent, measured as the smallest count at which the WORST array in the stage is under
 * 1% two-way error: 100 in a bulk ESC line, 50 to 100 in a trophectoderm biopsy, 400 in a WGA
 * blastomere. The default is the strictest of those, because the tool cannot see which stage it
 * has been handed and a wrong parent is worse than a refusal.
 *
 * The cost is stated plainly: 400 informative markers span a median 47.6 Mb, so at the default a
 * whole chromosome or a large arm can be annotated and a focal gain cannot. A 5 Mb region carries
 * a median 51 informative markers, which clears the trophectoderm bar and is an eighth of this.
 */
export const MIN_INFORMATIVE_DEFAULT = 400
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
 * Weaker than the two-parent form and deliberately marked so by the caller: the informative set is
 * smaller and it is conditioned on the sample being heterozygous, which a gain itself perturbs.
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
 * This is the reference point every window is judged against. Using the theoretical 0.5 instead is
 * the single mistake that inverts calls against a reconstructed parent.
 */
export const recentre = (all: readonly DosageMarker[]): number =>
  median(all.map((m) => m.patShare))

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
  why: string
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
    return { ...base, origin: 'unclear', why: 'the sample has no informative markers to centre on' }
  }
  if (informative < minInformative) {
    return {
      ...base,
      origin: 'unclear',
      why: `${informative} informative markers is under the ${minInformative} needed for a `
        + 'direction at this stage',
    }
  }
  if (Math.abs(deviation) < margin) {
    return {
      ...base,
      origin: 'unclear',
      why: `the share sits ${deviation >= 0 ? '+' : ''}${deviation.toFixed(3)} from this `
        + `sample's own centre, inside the ${margin} band that is not called either way`,
    }
  }
  return {
    ...base,
    origin: deviation > 0 ? 'paternal' : 'maternal',
    why: `${deviation > 0 ? 'paternal' : 'maternal'} alleles are over-represented by `
      + `${Math.abs(deviation).toFixed(3)} against this sample's own centre of `
      + `${centre.toFixed(3)}, over the ${margin} margin`,
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
 * which is below chance, and intensity does not rescue it either: 0.503 to 0.561 for one copy
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

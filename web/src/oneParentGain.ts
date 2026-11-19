/**
 * Which parent gave an extra chromosome with ONE parent loaded: the parts the copy-state share
 * (oneParentState.ts) needs before gainOrigin.scanAlleleGains can read it. Plan item 3.0d; scored in
 * audit/truth/extra-copy-one-parent-2026-11-16.txt.
 *
 * Notation: P the loaded parent, homozygous at a marker; the allele P lacks has founder frequency q;
 * y is the sample's BAF share of that allele. Two copies: y = 1/2 where the other parent passed the
 * allele (probability q), else 0. P's extra copy: 1/3. The other parent's extra copy: 2/3 where both
 * of its copies carry the allele, 1/3 where one does.
 *
 * 1. A GAIN READ AS A LOST PARENT. P's extra copy puts the copy-state share at 2/3 of the disomic
 *    level and the other parent's at 4/3, which on some cells crosses OTHER_LOST_MAX or
 *    PARENT_LOST_MIN. A lost copy cannot raise intensity, so where the intensity scan calls a gain
 *    the loss state is not reported and the chromosome goes to the gain test.
 * 2. A REBUILT P. A heterozygous parent all of whose m pronuclei showed one allele reads as that
 *    homozygote, with chance 2^-m per allele. Each marker of the gain test's share is weighted by the
 *    chance the homozygote is real, from q and m; a measured P has weight 1.
 * 3. THE HETEROZYGOUS CALLS, UPWARD ONLY. Where the sample is called heterozygous at P's homozygous
 *    markers, the median y is the position of the het cluster: 1/2 with two copies. Only the other
 *    parent's extra copy lifts it (two copies of an allele P does not carry); P's extra copy, and an
 *    extra copy made of both of the other parent's homologues, lower it alike, so downward names
 *    nothing.
 */
import { mad, median } from './intensityNull.ts'
import { COPY_SHIFT_GAIN_FLOOR } from './parentage.ts'
import { COMMON_ALLELE_MIN, MIN_SHARE_MARKERS, type ShareTally, type CopyState } from './oneParentState.ts'
import { ONE_PARENT_ORIENT_Z, ALLELE_GAIN_SCALE_FLOOR, ALLELE_GAIN_MIN_REFERENCE } from './gainOrigin.ts'

/**
 * Rule 1: a copy state that reads a lost parent where the intensity reads MORE DNA than the cell's
 * own level, whether or not a gain can be called there.
 *
 * THE TRIGGER IS THE INTENSITY, NOT THE CALL. Keyed on the gain CALL alone it left 3 of 158 built
 * gains printing as a lost parent with the measured father (the D2 line of
 * audit/truth/extra-copy-one-parent-2026-11-16.txt): the intensity on those chromosomes is up, but
 * not far enough past the array's own spread for a gain to be asserted, so the guard never fired
 * and the opposite direction was printed instead. A lost copy cannot raise intensity at all, so
 * evidence too weak to assert a gain is still enough to refuse its opposite. The magnitude asked
 * is the gain floor, parentage.COPY_SHIFT_GAIN_FLOOR.
 */
export const gainReadAsLoss = (
  state: CopyState | undefined, intensityGain: boolean, lrrShift = NaN,
): boolean => (state === 'other-lost' || state === 'parent-lost')
  && (intensityGain || lrrShift >= COPY_SHIFT_GAIN_FLOOR)

/** Rule 2: the chance a rebuilt homozygote is really homozygous. q: founder frequency of the allele
 *  it lacks. m: the fewest pronuclei the rebuild requires per marker. 1 where m is not known. */
export function rebuiltHomWeight(q: number, m: number | undefined): number {
  if (m === undefined || !Number.isFinite(m)) return 1
  const hom = (1 - q) ** 2
  const den = hom + 2 * q * (1 - q) * 2 ** -m
  return den > 0 ? hom / den : 0
}

/** The weighted S2 sums of the share; S1 (the base) stays in the ShareTally, where q = 0 weighs 1. */
export interface WeightedShare { w: number; wy: number; wq: number }
export const emptyWeightedShare = (): WeightedShare => ({ w: 0, wy: 0, wq: 0 })

/** One marker, under the same conditions oneParentState.addShare puts it into S2. */
export function addWeightedShare(
  t: WeightedShare, parent: string | null | undefined, baf: number | null | undefined,
  freqB: number | undefined, m: number | undefined,
): void {
  if (parent !== 'AA' && parent !== 'BB') return
  if (baf === null || baf === undefined || !Number.isFinite(baf) || freqB === undefined) return
  const q = parent === 'AA' ? freqB : 1 - freqB
  if (q < COMMON_ALLELE_MIN) return
  const w = rebuiltHomWeight(q, m)
  t.w += w; t.wy += w * (parent === 'AA' ? baf : 1 - baf); t.wq += w * q
}

/** The share r of oneParentState.shareCall with the weighted S2. NaN where shareCall reads none. */
export function weightedShareR(base: ShareTally, t: WeightedShare | undefined): number {
  if (!t || base.n1 < MIN_SHARE_MARKERS || base.n2 < MIN_SHARE_MARKERS || !(t.w > 0)) return NaN
  const b = base.y1 / base.n1
  const k = 1 - 2 * b
  return k > 0 && t.wq > 0 ? (t.wy / t.w - b) / (k * t.wq / t.w) : NaN
}

/** Rule 3's per-marker input: y where P is homozygous, q >= COMMON_ALLELE_MIN and the call is AB. */
export function hetCallShare(
  parent: string | null | undefined, genotype: string, baf: number | null | undefined,
  freqB: number | undefined,
): number | null {
  if ((parent !== 'AA' && parent !== 'BB') || genotype !== 'AB') return null
  if (baf === null || baf === undefined || !Number.isFinite(baf) || freqB === undefined) return null
  const q = parent === 'AA' ? freqB : 1 - freqB
  return q >= COMMON_ALLELE_MIN ? (parent === 'AA' ? baf : 1 - baf) : null
}

export interface OtherCopyCall { z: number; median: number; centre: number; scale: number; named: boolean }

/**
 * Rule 3: the het-call median of `target` against the same medians on `reference` (the cell's
 * other two-parent chromosomes). Null where either side has too few markers or chromosomes.
 */
export function otherCopyFromHetCalls(
  target: readonly number[], reference: readonly (readonly number[])[],
): OtherCopyCall | null {
  if (target.length < MIN_SHARE_MARKERS) return null
  const meds = reference.filter((r) => r.length >= MIN_SHARE_MARKERS).map((r) => median(r))
  if (meds.length < ALLELE_GAIN_MIN_REFERENCE) return null
  const m = median(target)
  const centre = median(meds)
  const scale = Math.max(mad(meds), ALLELE_GAIN_SCALE_FLOOR)
  const z = (m - centre) / scale
  return { z, median: m, centre, scale, named: z >= ONE_PARENT_ORIENT_Z }
}

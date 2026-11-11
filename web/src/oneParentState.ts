/**
 * Whether a chromosome carries both parents' copies, the loaded parent's alone, or the other's alone,
 * from ONE parent's array and a founder panel.
 *
 * WHY. The laboratory never has both parents arrayed: it has one, usually reconstructed from that
 * parent's haploid products and homozygous at every marker by construction. The two-parent reading of
 * a lost copy (obligateHet.hetCall) needs markers where the parents are opposite homozygotes, so with
 * one parent it never runs. This reads the same question from one parent P and the population.
 *
 * THE READING. At markers where P is homozygous, the sample's BAF share of the allele P lacks is y.
 * Split those markers by q, how common that allele is among unrelated founders:
 *
 *   S1  q = 0: every founder carries P's allele, so no copy can carry the other one. y here is drop-in
 *       and noise alone, `base`, and 1 - 2 x base is this chromosome's own compression.
 *   S2  q >= 0.5: the allele P lacks is common, so the other parent usually carries it.
 *
 *   r = (mean y over S2 - base) / ((1 - 2 base) x mean q over S2)
 *
 * A chromosome carrying the other parent's copy shows that allele in proportion to q, and r sits near
 * the share of the other parent's alleles that differ from P. One that has lost the other parent's copy
 * shows P's allele everywhere whatever q is, so r is near 0. One that has lost P's copy carries only the
 * other parent's, so r is near 1. Drop-in only raises `base`, which cancels in S2 - S1.
 *
 * Measured by the one-parent parity consult on GSE148488 against a father rebuilt from five of his
 * pronuclei: disomic r median 0.442 (0.343 to 0.611), the other parent's copy lost 0.027 (0.009 to
 * 0.151), P's copy lost 0.888 (0.675 to 0.962), and no template reads base 0.41 to 0.46. The rebuilt
 * father's calls matched the measured father's on 634 of 638 chromosomes.
 *
 * CHROMOSOME-LOCAL, NOT GENOME-WIDE. A genome-wide compression under-scales a chromosome whose base
 * differs: three published paternal monosomies read r 0.53 to 0.62 that way and were missed.
 *
 * WHAT IT CANNOT SEE. A gain: an extra copy of P's moves r toward two thirds of disomic and of the
 * other's toward four thirds, both inside the disomic range. Gains need intensity.
 */

export interface ShareTally {
  /** S1: markers where every founder carries P's allele. */
  n1: number; y1: number
  /** S2: markers where the allele P lacks has founder frequency 0.5 or more. */
  n2: number; y2: number; q2: number
}
export const emptyShareTally = (): ShareTally => ({ n1: 0, y1: 0, n2: 0, y2: 0, q2: 0 })

/** Founder frequency at or over which the allele P lacks is common enough to read a copy by. */
export const COMMON_ALLELE_MIN = 0.5

/**
 * One marker into the tally. `parent` is P's genotype, `baf` the sample's B allele frequency, which is
 * used at no-calls as well, since it is reported there and carries the signal. `freqB` is the founders'
 * B allele frequency and `monomorphic` whether the panel calls the marker one shared homozygote.
 */
export function addShare(
  t: ShareTally, parent: string | null | undefined, baf: number | null | undefined,
  freqB: number | undefined, monomorphic: boolean,
): void {
  if (parent !== 'AA' && parent !== 'BB') return
  if (baf === null || baf === undefined || !Number.isFinite(baf) || freqB === undefined) return
  const y = parent === 'AA' ? baf : 1 - baf
  const q = parent === 'AA' ? freqB : 1 - freqB
  if (q === 0 && monomorphic) { t.n1 += 1; t.y1 += y; return }
  if (q >= COMMON_ALLELE_MIN) { t.n2 += 1; t.y2 += y; t.q2 += q }
}

/** Markers needed in each set before a chromosome is read. Chromosome 21 carries a few thousand. */
export const MIN_SHARE_MARKERS = 100

/** `base` over which the chromosome carries no template: drop-in and noise fill half the calls. */
export const NO_TEMPLATE_BASE = 0.30
/** r under which the other parent's copy is gone, and over which P's is. */
export const OTHER_LOST_MAX = 0.20
export const PARENT_LOST_MIN = 0.75
/**
 * The low end of where a chromosome missing P's copy reads: the 5% quantile of r over the 59 such
 * chromosomes of Zuccaro 2020 Table S3 (median 0.888, 95% at 0.962), measured by the parity consult.
 * PARENT_LOST_MIN is the call, set high so a disomic chromosome is never called lost; this is the
 * range, for a test that asks whether ANY chromosome could be one.
 */
export const PARENT_LOST_LOW = 0.675

export type CopyState = 'disomic' | 'other-lost' | 'parent-lost' | 'no-template' | 'undecided'

export interface ShareCall {
  state: CopyState
  base: number
  r: number
  n1: number
  n2: number
}

export function shareCall(t: ShareTally): ShareCall {
  const base = t.n1 ? t.y1 / t.n1 : NaN
  if (t.n1 < MIN_SHARE_MARKERS || t.n2 < MIN_SHARE_MARKERS) {
    return { state: 'undecided', base, r: NaN, n1: t.n1, n2: t.n2 }
  }
  if (base > NO_TEMPLATE_BASE) return { state: 'no-template', base, r: NaN, n1: t.n1, n2: t.n2 }
  const k = 1 - 2 * base
  const dq = t.q2 / t.n2
  const r = k > 0 && dq > 0 ? (t.y2 / t.n2 - base) / (k * dq) : NaN
  const state: CopyState = !Number.isFinite(r) ? 'undecided'
    : r < OTHER_LOST_MAX ? 'other-lost' : r > PARENT_LOST_MIN ? 'parent-lost' : 'disomic'
  return { state, base, r, n1: t.n1, n2: t.n2 }
}

/**
 * The share of decided autosomes in the genome's most common copy state: the normal-majority test,
 * from the alleles. Intensity is read against the array's own per-chromosome median, which is a normal
 * chromosome only when most chromosomes share one state. That state is two copies in a diploid and one
 * in a haploid product, whichever parent it came from, so all three count and the commonest decides.
 * No template is never a state a majority can share: an array whose calls are noise reads it
 * throughout and scores 0, so the test fails closed.
 *
 * Counting disomy alone withheld copy number from every haploid product, which reads one copy on every
 * chromosome. On GSE148488 the cells published chaotic read 0.000 to 0.727 and the rest 0.909 or more
 * either way; one of them, staged haploid, is 15 chromosomes disomic and 8 single, 0.682.
 */
export function normalFraction(calls: readonly ShareCall[]): number {
  const decided = calls.filter((c) => c.state !== 'undecided')
  if (!decided.length) return NaN
  const n = (s: CopyState) => decided.filter((c) => c.state === s).length
  return Math.max(n('disomic'), n('other-lost'), n('parent-lost')) / decided.length
}

/** Normal fraction at or over which a normal majority is shown and copy number can be read. */
export const NORMAL_MAJORITY_MIN = 0.80

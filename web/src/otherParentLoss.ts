/**
 * Where the OTHER parent's copy is missing along a chromosome, with one parent P loaded.
 *
 * A marker is `seen` where the sample's call shows an allele P's copy would not give it: where P is
 * homozygous, the allele P lacks, which the other parent carries with founder frequency q; where P
 * is heterozygous, a heterozygous call. The rates are the array's own, read from its chromosomes
 * that carry both parents:
 *
 *   both copies          P homozygous   b + a x q      P heterozygous   h
 *   other parent's gone  either         b
 *
 * A run of markers scores the summed log-likelihood ratio, gone against present, divided by the
 * array's clustering. Calls on amplified material fail in clusters, so a sum over them is worth
 * fewer independent calls than it counts; `inflation` is how many fewer, measured on the same
 * reference chromosomes.
 *
 * WHAT IT CANNOT TELL APART. Where both parents carry the same haplotype the sample shows no allele
 * P lacks, exactly as where the other parent's copy is gone. Between unrelated parents such
 * stretches are short, which OTHER_MIN_SPAN_BP excludes. Between related parents they are not, and a
 * run is then true of the alleles and not of the copy number.
 *
 * Nor an extra copy of P's: it outnumbers the other parent's alleles and they drop out of the calls.
 * The caller drops a run that intensity reads as gained (`coveredShare`), and runs nothing on a
 * genome not read as two parents', where a haploid product shows P's other haplotype in blocks.
 *
 * Not run with both parents loaded: each parent's own homozygous markers answer this directly.
 * Measurements: audit/truth.
 *
 * Self-check:  node src/otherParentLoss.check.ts
 */
import { median } from './intensityNull.ts'
import { ABSENT_PARENT_OPP_MIN } from './obligateHet.ts'
import type { Segment } from './segments.ts'

export interface OtherMarker {
  pos: number
  /** Founder frequency of the allele P lacks, over 0; null where P is heterozygous. */
  q: number | null
  /** The call shows an allele P's copy would not give: the allele P lacks, or a heterozygote. */
  seen: boolean
  /** The call is the homozygote P cannot give. Implies `seen`; never set where P is heterozygous. */
  opp: boolean
}

export interface OtherRates {
  /** Rate at which a call shows an allele that is not there. */
  b: number
  /** Rise of the seen rate with q where P is homozygous and both copies are present. */
  a: number
  /** Heterozygous share where P is heterozygous and both copies are present; NaN with no such markers. */
  h: number
  /** Variance of block sums over what independent calls would give. Never under 1. */
  inflation: number
}

/** Divided score a run must reach. */
export const OTHER_LRT = 15
/**
 * Shortest run reported. Under it, a haplotype two unrelated parents share reads the same, and
 * nothing in one parent's alleles separates the two. A shorter run that passes every other test is
 * returned to the caller, which withholds it and says where it is.
 */
export const OTHER_MIN_SPAN_BP = 3_000_000
/** Markers per block when the clustering is measured: under the shortest run reported. */
export const BLOCK_MARKERS = 500
/** Blocks needed before the clustering is read. */
export const MIN_BLOCKS = 20
/** Lowest `b` taken from a fit. A fitted intercept at or under zero is the fit's error. */
export const NOISE_FLOOR = 0.005
/** Homozygous-class markers a reference chromosome needs before its line is fitted. */
const MIN_FIT_MARKERS = 500
/** Heterozygous-class markers a reference chromosome needs before its share is read. */
const MIN_HET_MARKERS = 200
/** Homozygous calls at P's homozygous markers a run needs before P's own copy is said to be there. */
const MIN_OPP_CALLS = 200
/** Median of a chi-squared variable with one degree of freedom. */
const CHI2_MEDIAN = 0.4549

/** The seen rate with both copies present. At or under `b` the marker says nothing. */
const present = (m: OtherMarker, r: { b: number; a: number; h: number }): number => (
  Math.min(0.999, m.q === null ? r.h : r.b + r.a * m.q))

/**
 * The array's rates and clustering, from chromosomes that carry both parents. Each is one
 * chromosome's markers in position order. Null where fewer than two can be fitted or the clustering
 * cannot be measured, and nothing is then scanned.
 */
export function otherRates(reference: readonly (readonly OtherMarker[])[]): OtherRates | null {
  const bs: number[] = []
  const as: number[] = []
  const hs: number[] = []
  for (const ms of reference) {
    let n = 0; let sq = 0; let sqq = 0; let ss = 0; let sqs = 0; let nh = 0; let sh = 0
    for (const m of ms) {
      if (m.q === null) { nh += 1; if (m.seen) sh += 1; continue }
      n += 1; sq += m.q; sqq += m.q * m.q
      if (m.seen) { ss += 1; sqs += m.q }
    }
    const varQ = n ? sqq / n - (sq / n) ** 2 : 0
    if (n >= MIN_FIT_MARKERS && varQ > 0) {
      const a = (sqs / n - (sq / n) * (ss / n)) / varQ
      as.push(a)
      bs.push(ss / n - a * (sq / n))
    }
    if (nh >= MIN_HET_MARKERS) hs.push(sh / nh)
  }
  if (bs.length < 2) return null
  const r = {
    b: Math.min(0.5, Math.max(NOISE_FLOOR, median(bs))),
    a: median(as),
    h: hs.length ? median(hs) : NaN,
  }
  const ratios: number[] = []
  for (const ms of reference) {
    let k = 0; let x = 0; let v = 0
    for (const m of ms) {
      const p = present(m, r)
      if (!(p > r.b)) continue
      x += (m.seen ? 1 : 0) - p
      v += p * (1 - p)
      k += 1
      if (k === BLOCK_MARKERS) { ratios.push((x * x) / v); k = 0; x = 0; v = 0 }
    }
  }
  if (ratios.length < MIN_BLOCKS) return null
  return { ...r, inflation: Math.max(1, median(ratios) / CHI2_MEDIAN) }
}

/**
 * Runs along one chromosome where the other parent's copy is missing. `markers` in position order.
 *
 * The best-scoring run is taken, then the best on either side of it, until none reaches the
 * threshold. A run is returned where it also spans `minSpanBp` and P's own copy is shown to be the
 * one that remains: under ABSENT_PARENT_OPP_MIN of the homozygous calls at P's homozygous markers
 * inside it are the homozygote P cannot give. Only a run of OTHER_MIN_SPAN_BP or more is a finding.
 */
export function scanOtherParent(
  chrom: string, markers: readonly OtherMarker[], rates: OtherRates, minSpanBp = OTHER_MIN_SPAN_BP,
): Segment[] {
  const ms = markers.filter((m) => present(m, rates) > rates.b)
  const n = ms.length
  const rate = new Float64Array(n)
  const cum = new Float64Array(n + 1)
  for (let i = 0; i < n; i += 1) {
    const p = present(ms[i], rates)
    rate[i] = p
    cum[i + 1] = cum[i] + (ms[i].seen ? Math.log(rates.b / p) : Math.log((1 - rates.b) / (1 - p)))
  }
  const need = OTHER_LRT * rates.inflation
  const out: Segment[] = []
  const todo: [number, number][] = [[0, n]]
  while (todo.length) {
    const [lo, hi] = todo.pop()!
    let low = lo; let from = lo; let to = lo; let best = 0
    for (let j = lo + 1; j <= hi; j += 1) {
      if (cum[j - 1] < cum[low]) low = j - 1
      if (cum[j] - cum[low] > best) { best = cum[j] - cum[low]; from = low; to = j }
    }
    if (best < need) continue
    todo.push([lo, from], [to, hi])
    const spanBp = ms[to - 1].pos - ms[from].pos
    if (spanBp < minSpanBp) continue
    let seen = 0; let expected = 0; let homCalls = 0; let opp = 0
    for (let k = from; k < to; k += 1) {
      const m = ms[k]
      expected += rate[k]
      if (m.seen) seen += 1
      if (m.q !== null && (m.opp || !m.seen)) { homCalls += 1; if (m.opp) opp += 1 }
    }
    if (homCalls < MIN_OPP_CALLS || opp / homCalls >= ABSENT_PARENT_OPP_MIN) continue
    const w = to - from
    out.push({
      kind: 'parental-absence', chrom, markers: w, absent: w - seen, rate: (w - seen) / w,
      startBp: ms[from].pos, endBp: ms[to - 1].pos, spanBp,
      score: best / rates.inflation, nullRate: 1 - expected / w,
    })
  }
  return out.sort((x, y) => x.startBp - y.startBp)
}

/** The share of [startBp, endBp] that `intervals` cover between them. */
export function coveredShare(
  startBp: number, endBp: number, intervals: readonly { startBp: number; endBp: number }[],
): number {
  let covered = 0
  let reach = startBp
  for (const iv of [...intervals].sort((x, y) => x.startBp - y.startBp)) {
    const from = Math.max(iv.startBp, reach)
    const to = Math.min(iv.endBp, endBp)
    if (to > from) { covered += to - from; reach = to }
  }
  return endBp > startBp ? covered / (endBp - startBp) : 0
}

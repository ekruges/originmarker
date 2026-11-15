/**
 * Where a parent's copy is missing along a chromosome, read against one loaded parent P: the OTHER
 * parent's copy (scanOtherParent), and P's own (scanOwnCopy).
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
 * Nor an extra copy: it outnumbers the other parent's alleles and they drop out of the calls, in
 * either direction. The caller refuses a run that intensity reads as gained, and runs nothing on a
 * genome not read as two parents', where a haploid product shows P's other haplotype in blocks.
 *
 * P'S OWN COPY is read at the markers where P is homozygous: with it gone the call is the homozygote
 * P cannot give, at the rate the other parent's allele shows alone (e + a x q), against e + d x q
 * with both present. The same runs, the same division, on the opposite-homozygote calls. It reaches
 * what the window scan in segments.ts cannot for want of markers, a short loss or a parent rebuilt
 * from few cells.
 *
 * The other parent's copy is not read this way with both parents loaded: each parent's own
 * homozygous markers answer it directly.
 * Measurements: audit/truth.
 *
 * Self-check:  node src/otherParentLoss.check.ts
 */
import { median } from './intensityNull.ts'
import { ABSENT_PARENT_OPP_MIN } from './obligateHet.ts'
import { LR_DROP } from './refine.ts'
import { NULL_FLOOR, type Segment } from './segments.ts'

export interface OtherMarker {
  /** The probe, so another array of the run can be read at the same marker. */
  id: string
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
/**
 * Fall in the divided score that bounds an edge: every place the run could start, or end, and score
 * within this of its best. The window scan's own bound, on this scan's scale.
 */
export const EDGE_DROP = LR_DROP
/**
 * Share of blocks at or over which another array reads as the same genome as this one: midway
 * between the most arrays of two embryos of one couple reach and the least arrays of one embryo give.
 */
export const SAME_GENOME_SHARE = 0.75
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
 * Every run of a cumulative score that gains `need` or more: the best, then the best on either side
 * of it, until none reaches it. Each is [first marker, one past the last, gain].
 */
function bestRuns(cum: Float64Array, need: number): [number, number, number][] {
  const out: [number, number, number][] = []
  const todo: [number, number][] = [[0, cum.length - 1]]
  while (todo.length) {
    const [lo, hi] = todo.pop()!
    let low = lo; let from = lo; let to = lo; let best = 0
    for (let j = lo + 1; j <= hi; j += 1) {
      if (cum[j - 1] < cum[low]) low = j - 1
      if (cum[j] - cum[low] > best) { best = cum[j] - cum[low]; from = low; to = j }
    }
    if (best < need) continue
    todo.push([lo, from], [to, hi])
    out.push([from, to, best])
  }
  return out
}

/**
 * The edges of a run [from, to) and the interval around each, as positions. An edge's interval is
 * the unbroken stretch of places it could move to with the run's score within `drop` of its best.
 */
function edgesOf(
  cum: Float64Array, pos: (i: number) => number, from: number, to: number, drop: number,
): NonNullable<Segment['refined']> {
  let startLo = from
  while (startLo > 0 && cum[startLo - 1] <= cum[from] + drop) startLo -= 1
  let startHi = from
  while (startHi < to - 1 && cum[startHi + 1] <= cum[from] + drop) startHi += 1
  let endLo = to
  while (endLo > from + 1 && cum[endLo - 1] >= cum[to] - drop) endLo -= 1
  let endHi = to
  while (endHi < cum.length - 1 && cum[endHi + 1] >= cum[to] - drop) endHi += 1
  return {
    startBp: pos(from), startLoBp: pos(startLo), startHiBp: pos(startHi),
    endBp: pos(to - 1), endLoBp: pos(endLo - 1), endHiBp: pos(endHi - 1),
    spanBp: pos(to - 1) - pos(from),
  }
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
  const out: Segment[] = []
  for (const [from, to, best] of bestRuns(cum, OTHER_LRT * rates.inflation)) {
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
      refined: edgesOf(cum, (i) => ms[i].pos, from, to, EDGE_DROP * rates.inflation),
    })
  }
  return out.sort((x, y) => x.startBp - y.startBp)
}

export interface OwnRates {
  /** Opposite-homozygote rate where the other parent cannot have given the allele P lacks. */
  e: number
  /** Rise of that rate with q with both copies present: P's allele dropped from a heterozygote. */
  d: number
  /** Its rise with q with P's copy gone: the other parent's allele alone, at the rate it shows. */
  a: number
  /** Variance of block sums over what independent calls would give. Never under 1. */
  inflation: number
}

/**
 * The array's rates and clustering for P's OWN copy, from chromosomes that carry both parents. `a`
 * is `OtherRates.a` of the same parent. Null where fewer than two chromosomes can be fitted, the
 * clustering cannot be measured, or a missing copy would read no higher than a present one.
 */
export function ownRates(reference: readonly (readonly OtherMarker[])[], a: number): OwnRates | null {
  const es: number[] = []
  const ds: number[] = []
  for (const ms of reference) {
    let n = 0; let sq = 0; let sqq = 0; let so = 0; let sqo = 0
    for (const m of ms) {
      if (m.q === null) continue
      n += 1; sq += m.q; sqq += m.q * m.q
      if (m.opp) { so += 1; sqo += m.q }
    }
    const varQ = n ? sqq / n - (sq / n) ** 2 : 0
    if (n < MIN_FIT_MARKERS || !(varQ > 0)) continue
    const d = (sqo / n - (sq / n) * (so / n)) / varQ
    ds.push(d)
    es.push(so / n - d * (sq / n))
  }
  if (es.length < 2) return null
  const e = Math.min(0.5, Math.max(NULL_FLOOR, median(es)))
  const d = Math.max(0, median(ds))
  if (!(a > d)) return null
  const ratios: number[] = []
  for (const ms of reference) {
    let k = 0; let x = 0; let v = 0
    for (const m of ms) {
      if (m.q === null) continue
      const p = Math.min(0.999, e + d * m.q)
      x += (m.opp ? 1 : 0) - p
      v += p * (1 - p)
      k += 1
      if (k === BLOCK_MARKERS) { ratios.push((x * x) / v); k = 0; x = 0; v = 0 }
    }
  }
  if (ratios.length < MIN_BLOCKS) return null
  return { e, d, a, inflation: Math.max(1, median(ratios) / CHI2_MEDIAN) }
}

/**
 * Runs along one chromosome where P's OWN copy is missing: the call is the homozygote P cannot
 * give, at the rate the other parent's allele shows alone. `markers` in position order.
 *
 * An opposite homozygote is P's allele seen missing, not merely unseen, so a haplotype the parents
 * share does not produce it and the floor is the one every segment has.
 */
export function scanOwnCopy(
  chrom: string, markers: readonly OtherMarker[], rates: OwnRates, minSpanBp: number,
): Segment[] {
  const ms = markers.filter((m) => m.q !== null)
  const n = ms.length
  const cum = new Float64Array(n + 1)
  const base = new Float64Array(n)
  for (let i = 0; i < n; i += 1) {
    const q = ms[i].q as number
    const p0 = Math.min(0.999, rates.e + rates.d * q)
    const p1 = Math.min(0.999, rates.e + rates.a * q)
    base[i] = p0
    cum[i + 1] = cum[i] + (ms[i].opp ? Math.log(p1 / p0) : Math.log((1 - p1) / (1 - p0)))
  }
  const out: Segment[] = []
  for (const [from, to, best] of bestRuns(cum, OTHER_LRT * rates.inflation)) {
    const spanBp = ms[to - 1].pos - ms[from].pos
    if (spanBp < minSpanBp) continue
    let opp = 0; let expected = 0
    for (let k = from; k < to; k += 1) { expected += base[k]; if (ms[k].opp) opp += 1 }
    const w = to - from
    out.push({
      kind: 'parental-absence', chrom, markers: w, absent: opp, rate: opp / w,
      startBp: ms[from].pos, endBp: ms[to - 1].pos, spanBp,
      score: best / rates.inflation, nullRate: expected / w,
      refined: edgesOf(cum, (i) => ms[i].pos, from, to, EDGE_DROP * rates.inflation),
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

/**
 * How much of the genome another array shares with this one, as the share of blocks that read the
 * same copy of the other parent. `other` holds the other array's reading, in genome order, at each
 * marker where P is homozygous and THIS array shows the allele P lacks: an array carrying the same
 * copy shows it too, at b + a, and one carrying an unrelated copy at b + a x q. NaN on too few blocks.
 *
 * Arrays of one embryo share both parents' copies throughout, so what the parents share reads as
 * missing in all of them and a loss in one does not. Arrays of two embryos share about half.
 */
export function sameGenomeShare(other: readonly { q: number; seen: boolean }[], rates: OtherRates): number {
  const same = Math.min(0.999, rates.b + rates.a)
  let blocks = 0; let agree = 0; let k = 0; let sum = 0
  for (const m of other) {
    const apart = Math.min(0.999, rates.b + rates.a * m.q)
    sum += m.seen ? Math.log(same / apart) : Math.log((1 - same) / (1 - apart))
    k += 1
    if (k === BLOCK_MARKERS) { blocks += 1; if (sum > 0) agree += 1; k = 0; sum = 0 }
  }
  return blocks >= MIN_BLOCKS ? agree / blocks : NaN
}

/**
 * How strongly an array shows the other parent's alleles across [startBp, endBp]: the summed
 * log-likelihood ratio, present against gone, divided by the array's clustering.
 */
export function presentScore(
  markers: readonly OtherMarker[], rates: OtherRates, startBp: number, endBp: number,
): number {
  let sum = 0
  for (const m of markers) {
    if (m.pos < startBp || m.pos > endBp) continue
    const p = present(m, rates)
    if (!(p > rates.b)) continue
    sum += m.seen ? Math.log(p / rates.b) : Math.log((1 - p) / (1 - rates.b))
  }
  return sum / rates.inflation
}

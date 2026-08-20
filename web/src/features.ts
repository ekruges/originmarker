/**
 * Where the called regions fall, against fragile sites, long genes and centromeres. hg19.
 *
 * This asks a different question from everything else in Syngamy: not WHICH PARENT a change came
 * from, but WHERE in the genome it landed. It therefore needs only breakpoint position and no
 * parental genotype, which is why it can run on every sample rather than the few with a usable
 * reference parent.
 *
 * THE NULL IS THE ENTIRE DIFFICULTY, and it is why this lives in the app rather than in a script
 * run afterwards. A region can only be called where the array carries markers AND where
 * amplification produced calls. Draw random intervals uniformly along the genome and you compare
 * callable genome against uncallable, which reports an enrichment for anything correlated with
 * marker density - gene density is, and so, therefore, are fragile sites, which sit in large genes.
 * The permutation here draws each null interval on the SAME CHROMOSOME with the SAME NUMBER OF
 * INFORMATIVE MARKERS as the observed region, so a 2,000-marker region is compared against other
 * 2,000-marker stretches rather than against 2,000 kilobases of anything.
 *
 * Only the sample's own marker positions can supply that, and only Syngamy has them.
 *
 * WHAT IT CANNOT SAY. Nothing about parent of origin. An external review notes the late-replicating
 * fragile compartment is established symmetrically on both parental genomes from the first cell
 * cycle, so a positional result here does not license a parental one.
 */

/** One interval on the reference. `name` is the feature, e.g. FRA3B or WWOX. */
export interface FeatureInterval {
  chrom: string
  startBp: number
  endBp: number
  name: string
}

export interface FeatureTrack {
  build: string
  fragile: FeatureInterval[]
  longGenes: FeatureInterval[]
  gaps: (FeatureInterval & { kind: string })[]
  /** genes per 1-Mb bin, chrom -> bin index -> count */
  geneDensityPerMb: Record<string, Record<string, number>>
  /**
   * Repli-seq VALLEYS, which are local minima of the wave signal and therefore the latest
   * replicating POINTS rather than late domains: 4,843 of them spanning 5 Mb in total, where the
   * genome's late-replicating fraction is around a third. They are the sharp marker associated
   * with common fragile sites, and must not be described as late-replicating regions.
   * ES is BG02ES, an embryonic stem line, which is the closest available cell type to this
   * material. Constitutive is late in at least 8 of 10 ENCODE lines.
   */
  lateReplicationValleysES: FeatureInterval[]
  lateReplicationValleysConstitutive: FeatureInterval[]
  /**
   * The WIDENED track set, name -> intervals, scored alongside the five above. Open-ended
   * deliberately: a track is added by putting it in this map and naming its scoring mode in
   * `modes`, not by editing scoreAll, so the shipped code path does not have to change every
   * time a factor is added.
   */
  extra?: Record<string, FeatureInterval[]>
  /**
   * Scoring mode per extra track. THIS IS NOT COSMETIC, and defaulting it is how the ES valley
   * track produced a meaningless 1.00. The observed regions have a MEDIAN SIZE OF 29 Mb, so a
   * track of many small intervals is touched by every region AND by every marker-matched null
   * interval; binary overlap then reports 1.00 by construction and carries no information.
   *   binary   - fraction of regions touching the track. Only valid when a typical region does
   *              NOT always contain one; the builder measures this and refuses otherwise.
   *   density  - intervals per Mb of region. For sparse point-like features (CTCF sites, CpG
   *              islands, SV calls) where the COUNT carries the signal.
   *   coverage - fraction of the region's base pairs covered. For broad domains (replication
   *              timing, LADs, heterochromatin) where the EXTENT carries the signal.
   */
  modes?: Record<string, 'binary' | 'density' | 'coverage'>
}

export interface Region {
  chrom: string
  startBp: number
  endBp: number
}

export interface Enrichment {
  feature: string
  /** fraction of observed regions touching the feature */
  observed: number
  nullMean: number
  ratio: number
  /** permutation p, two-sided by construction: the tail the observation falls in. */
  p: number
  regions: number
  /** Names hit, so a result can be read rather than only counted. */
  hits: string[]
  why: string
  /**
   * The null distribution itself, kept so a report can DRAW it. A permutation p is only
   * interpretable against the distribution it came from: a ratio of 0.97 reached p 0.002 on this
   * data because the null was tight, and a reader shown only the p would call that a finding.
   * `q` is the 2.5th, 50th and 97.5th percentiles; `hist` is 20 equal bins between `lo` and `hi`.
   */
  nullQuantiles: [number, number, number]
  nullHist: { lo: number, hi: number, counts: number[] }
}

/** Fewest permutations that can produce a p below 0.05 at all. */
export const MIN_PERMUTATIONS = 200
export const DEFAULT_PERMUTATIONS = 2000

/** Fewest regions before an enrichment is reported rather than refused. */
export const MIN_REGIONS = 5

const touches = (f: FeatureInterval, r: Region): boolean =>
  f.chrom === r.chrom && f.startBp < r.endBp && r.startBp < f.endBp

/** Deterministic generator, so a reported p is reproducible from the same inputs. */
const rng = (seed: number) => () => {
  seed = (seed * 1664525 + 1013904223) >>> 0
  return seed / 4294967296
}

/**
 * Markers strictly inside a region, by binary search on a sorted position array.
 * The count, not the span, is what the null matches on.
 */
export function markersIn(sorted: number[], from: number, to: number): number {
  let lo = 0
  let hi = sorted.length
  while (lo < hi) {
    const m = (lo + hi) >> 1
    if (sorted[m] < from) lo = m + 1
    else hi = m
  }
  const start = lo
  hi = sorted.length
  while (lo < hi) {
    const m = (lo + hi) >> 1
    if (sorted[m] <= to) lo = m + 1
    else hi = m
  }
  return lo - start
}

/**
 * Enrichment of `regions` in `features`, against a null matched on marker count per chromosome.
 *
 * `markersByChrom` must be the informative markers the regions were actually called from, sorted
 * ascending. Passing every marker on the array instead overstates the callable genome and returns
 * an enrichment that is a property of amplification.
 */
export function enrichment(
  label: string,
  features: FeatureInterval[],
  regions: Region[],
  markersByChrom: Map<string, number[]>,
  permutations = DEFAULT_PERMUTATIONS,
  seed = 7,
): Enrichment {
  const usable = regions.filter((r) => (markersByChrom.get(r.chrom)?.length ?? 0) > 10)
  const hits = features.filter((f) => usable.some((r) => touches(f, r))).map((f) => f.name)
  const observed = usable.length ? usable.filter((r) => features.some((f) => touches(f, r))).length
    / usable.length : NaN

  const emptyNull = {
    nullQuantiles: [NaN, NaN, NaN] as [number, number, number],
    nullHist: { lo: NaN, hi: NaN, counts: [] },
  }
  if (usable.length < MIN_REGIONS) {
    return {
      feature: label, observed, nullMean: NaN, ratio: NaN, p: NaN,
      regions: usable.length, hits, ...emptyNull,
      why: `${usable.length} regions is under the ${MIN_REGIONS} this needs for a rate`,
    }
  }

  const next = rng(seed)
  const draws: number[] = []
  for (let t = 0; t < permutations; t += 1) {
    let hit = 0
    for (const r of usable) {
      const m = markersByChrom.get(r.chrom)!
      const n = markersIn(m, r.startBp, r.endBp)
      if (n < 1 || n >= m.length) continue
      const i = Math.floor(next() * (m.length - n))
      if (features.some((f) => touches(f, { chrom: r.chrom, startBp: m[i], endBp: m[i + n] }))) {
        hit += 1
      }
    }
    draws.push(hit / usable.length)
  }
  const nullMean = draws.reduce((a, b) => a + b, 0) / draws.length
  // The tail the observation actually falls in, so depletion is reportable as well as enrichment.
  const tail = observed >= nullMean
    ? draws.filter((d) => d >= observed).length
    : draws.filter((d) => d <= observed).length
  const p = (1 + tail) / (permutations + 1)
  const sorted = [...draws].sort((a, b) => a - b)
  const q = (f: number): number => sorted[Math.min(sorted.length - 1, Math.floor(f * sorted.length))]
  const lo = Math.min(sorted[0], observed)
  const hi = Math.max(sorted[sorted.length - 1], observed)
  const BINS = 20
  const counts = new Array<number>(BINS).fill(0)
  const width = (hi - lo) || 1
  for (const d of draws) counts[Math.min(BINS - 1, Math.floor(((d - lo) / width) * BINS))] += 1
  return {
    nullQuantiles: [q(0.025), q(0.5), q(0.975)],
    nullHist: { lo, hi, counts },
    feature: label,
    observed,
    nullMean,
    ratio: nullMean > 0 ? observed / nullMean : NaN,
    p,
    regions: usable.length,
    hits,
    why: `${(observed * 100).toFixed(1)}% of ${usable.length} regions touch ${label}, against `
      + `${(nullMean * 100).toFixed(1)}% for intervals matched on marker count, `
      + `p ${p < 0.001 ? '< 0.001' : p.toFixed(3)}`,
  }
}

/**
 * A per-chromosome sorted interval index with a running maximum end.
 *
 * WHY THIS EXISTS. `enrichment` scans every feature for every region for every permutation, which
 * is fine for 24 fragile sites and intolerable for the 494,000 G4-seq intervals the widened track
 * set carries: 2,000 permutations x 417 regions x 494,000 intervals is 4 x 10^11 comparisons.
 * The index makes each query a binary search plus a short walk. `maxEnd[i]` is the largest end
 * among intervals 0..i, which is what lets a query stop early on intervals that start before the
 * region but end well before it too.
 */
export interface IntervalIndex {
  byChrom: Map<string, {
    /** raw interval bounds, each sorted independently, for counting */
    starts: number[]
    endsSorted: number[]
    names: string[]
    /** the same intervals MERGED (disjoint, ascending) plus a prefix sum of covered bp */
    mStart: number[]
    mEnd: number[]
    mCum: number[]
  }>
}

export function buildIndex(features: FeatureInterval[]): IntervalIndex {
  const byChrom: IntervalIndex['byChrom'] = new Map()
  const grouped = new Map<string, FeatureInterval[]>()
  for (const f of features) {
    const g = grouped.get(f.chrom)
    if (g) g.push(f)
    else grouped.set(f.chrom, [f])
  }
  for (const [c, fs] of grouped) {
    fs.sort((a, b) => a.startBp - b.startBp)
    const starts = fs.map((f) => f.startBp)
    const names = fs.map((f) => f.name)
    const endsSorted = fs.map((f) => f.endBp).sort((a, b) => a - b)
    const mStart: number[] = []
    const mEnd: number[] = []
    let s = fs[0].startBp
    let e = fs[0].endBp
    for (let i = 1; i < fs.length; i += 1) {
      if (fs[i].startBp <= e) e = Math.max(e, fs[i].endBp)
      else { mStart.push(s); mEnd.push(e); s = fs[i].startBp; e = fs[i].endBp }
    }
    mStart.push(s); mEnd.push(e)
    const mCum = new Array<number>(mStart.length + 1)
    mCum[0] = 0
    for (let i = 0; i < mStart.length; i += 1) mCum[i + 1] = mCum[i] + (mEnd[i] - mStart[i])
    byChrom.set(c, { starts, endsSorted, names, mStart, mEnd, mCum })
  }
  return { byChrom }
}

/** Number of entries strictly less than `p`. */
const countBelow = (a: number[], p: number): number => {
  let lo = 0
  let hi = a.length
  while (lo < hi) {
    const m = (lo + hi) >> 1
    if (a[m] < p) lo = m + 1
    else hi = m
  }
  return lo
}

/** Number of entries less than or equal to `p`. */
const countAtOrBelow = (a: number[], p: number): number => {
  let lo = 0
  let hi = a.length
  while (lo < hi) {
    const m = (lo + hi) >> 1
    if (a[m] <= p) lo = m + 1
    else hi = m
  }
  return lo
}

/**
 * Intervals overlapping [from, to), in O(log n).
 *
 * An interval overlaps iff it starts before `to` AND ends after `from`. Every interval that ends
 * at or before `from` also starts before `to` (start < end <= from < to), so the overlapping
 * count is exactly (# starting before `to`) minus (# ending at or before `from`) - which holds
 * whether or not the intervals overlap each other, and needs only two binary searches. The
 * alternative, walking the candidate list, is what made the 494,000-interval G4 track
 * unscoreable: 2,000 permutations x 417 regions x a walk of every interval in a 29-Mb window.
 */
export function countIn(idx: IntervalIndex, chrom: string, from: number, to: number): number {
  const c = idx.byChrom.get(chrom)
  if (!c) return 0
  return countBelow(c.starts, to) - countAtOrBelow(c.endsSorted, from)
}

export const overlapsIndexed = (idx: IntervalIndex, chrom: string, from: number, to: number):
  boolean => countIn(idx, chrom, from, to) > 0

/**
 * Base pairs of [from, to) covered by the track, counting a base once however many intervals
 * cover it. Read off the MERGED prefix sum, so it is two binary searches and two clipped edges
 * rather than a scan.
 */
export function coveredBp(idx: IntervalIndex, chrom: string, from: number, to: number): number {
  const c = idx.byChrom.get(chrom)
  if (!c || to <= from) return 0
  const lo = countAtOrBelow(c.mEnd, from)      // first merged interval that can reach `from`
  const hi = countBelow(c.mStart, to)          // one past the last that starts before `to`
  if (hi <= lo) return 0
  let total = c.mCum[hi] - c.mCum[lo]
  // clip the two boundary intervals to the query window
  if (c.mStart[lo] < from) total -= from - c.mStart[lo]
  if (c.mEnd[hi - 1] > to) total -= c.mEnd[hi - 1] - to
  return total
}

export type ScoreMode = 'binary' | 'density' | 'coverage'

/** Most names reported per track. The statistic is in `observed`; this list is for reading. */
export const HIT_NAMES = 200
/** No track interval here is longer than this, so the name scan can stop. */
const MAX_FEATURE_SPAN = 10_000_000

/**
 * Enrichment under any of the three statistics, against THE SAME NULL as `enrichment` above.
 *
 * The null is not re-derived here and must not be: each null interval is drawn on the SAME
 * chromosome with the SAME NUMBER OF INFORMATIVE MARKERS as the observed region, using the same
 * generator, the same seed and the same draw order, so the only thing that differs from the
 * validated path is WHICH STATISTIC is computed on each drawn interval. `features.check.ts`
 * asserts that mode 'binary' reproduces `enrichment` exactly on the same inputs; if that
 * assertion ever fails, this function has forked the null and the numbers it reports are not
 * comparable with the validation figures.
 *
 * For 'binary' the statistic is a RATE over regions (fraction touching), as before. For
 * 'density' and 'coverage' it is a MEAN over regions of a per-region quantity, so a ratio near
 * 1.00 means the regions look like the rest of the callable genome on that measure.
 */
export function enrichmentBy(
  label: string,
  features: FeatureInterval[],
  regions: Region[],
  markersByChrom: Map<string, number[]>,
  mode: ScoreMode = 'binary',
  permutations = DEFAULT_PERMUTATIONS,
  seed = 7,
  idxIn?: IntervalIndex,
): Enrichment {
  const idx = idxIn ?? buildIndex(features)
  const usable = regions.filter((r) => (markersByChrom.get(r.chrom)?.length ?? 0) > 10)

  const stat = (chrom: string, from: number, to: number): number => {
    if (mode === 'binary') return overlapsIndexed(idx, chrom, from, to) ? 1 : 0
    const span = Math.max(to - from, 1)
    if (mode === 'density') return countIn(idx, chrom, from, to) / (span / 1e6)
    return coveredBp(idx, chrom, from, to) / span
  }

  // Names hit, so a result can be read rather than only counted. Capped: the widened tracks run
  // to hundreds of thousands of intervals and an unbounded name list is not a readable output.
  // Names hit, so a result can be read rather than only counted. CAPPED at HIT_NAMES: the
  // widened tracks run to hundreds of thousands of intervals, and a name list that long is not a
  // readable output and costs more to build than the statistic itself. The count is in
  // `observed`; this list is for identifying WHICH features, which only reads usefully when
  // there are few of them.
  const hits: string[] = []
  {
    const seen = new Set<string>()
    for (const r of usable) {
      const c = idx.byChrom.get(r.chrom)
      // SKIP this region, do not stop: a region on a chromosome the track has no intervals on
      // says nothing about the regions after it. An earlier revision used `break` here, which
      // silently truncated the name list at the first such region.
      if (!c) continue
      const upper = countBelow(c.starts, r.endBp)
      for (let i = upper - 1; i >= 0 && hits.length < HIT_NAMES; i -= 1) {
        if (c.starts[i] < r.startBp - MAX_FEATURE_SPAN) break
        if (c.endsSorted.length && c.names[i] !== undefined && !seen.has(c.names[i])) {
          seen.add(c.names[i])
          hits.push(c.names[i])
        }
      }
      if (hits.length >= HIT_NAMES) break
    }
  }

  const observed = usable.length
    ? usable.reduce((a, r) => a + stat(r.chrom, r.startBp, r.endBp), 0) / usable.length
    : NaN

  const emptyNull = {
    nullQuantiles: [NaN, NaN, NaN] as [number, number, number],
    nullHist: { lo: NaN, hi: NaN, counts: [] },
  }
  if (usable.length < MIN_REGIONS) {
    return {
      feature: label, observed, nullMean: NaN, ratio: NaN, p: NaN,
      regions: usable.length, hits, ...emptyNull,
      why: `${usable.length} regions is under the ${MIN_REGIONS} this needs for a rate`,
    }
  }

  const next = rng(seed)
  const draws: number[] = []
  for (let t = 0; t < permutations; t += 1) {
    let acc = 0
    for (const r of usable) {
      const m = markersByChrom.get(r.chrom)!
      const n = markersIn(m, r.startBp, r.endBp)
      if (n < 1 || n >= m.length) continue
      const i = Math.floor(next() * (m.length - n))
      acc += stat(r.chrom, m[i], m[i + n])
    }
    draws.push(acc / usable.length)
  }
  const nullMean = draws.reduce((a, b) => a + b, 0) / draws.length
  const tail = observed >= nullMean
    ? draws.filter((d) => d >= observed).length
    : draws.filter((d) => d <= observed).length
  const p = (1 + tail) / (permutations + 1)
  const sorted = [...draws].sort((a, b) => a - b)
  const q = (f: number): number => sorted[Math.min(sorted.length - 1, Math.floor(f * sorted.length))]
  const lo = Math.min(sorted[0], observed)
  const hi = Math.max(sorted[sorted.length - 1], observed)
  const BINS = 20
  const counts = new Array<number>(BINS).fill(0)
  const width = (hi - lo) || 1
  for (const d of draws) counts[Math.min(BINS - 1, Math.floor(((d - lo) / width) * BINS))] += 1
  const unit = mode === 'binary' ? 'of regions touch'
    : mode === 'density' ? 'intervals/Mb in' : 'of region bp covered by'
  const fmt = (v: number) => (mode === 'density' ? v.toFixed(3) : `${(v * 100).toFixed(1)}%`)
  return {
    nullQuantiles: [q(0.025), q(0.5), q(0.975)],
    nullHist: { lo, hi, counts },
    feature: label,
    observed,
    nullMean,
    ratio: nullMean > 0 ? observed / nullMean : NaN,
    p,
    regions: usable.length,
    hits,
    why: `${fmt(observed)} ${unit} ${label} over ${usable.length} regions, against `
      + `${fmt(nullMean)} for intervals matched on marker count, `
      + `p ${p < 0.001 ? '< 0.001' : p.toFixed(3)} (${mode})`,
  }
}

/** Every track in one pass. Gaps are reported so artefact near a centromere is visible. */
export function scoreAll(
  track: FeatureTrack,
  regions: Region[],
  markersByChrom: Map<string, number[]>,
  permutations = DEFAULT_PERMUTATIONS,
): Enrichment[] {
  const five = [
    enrichment('common fragile site', track.fragile, regions, markersByChrom, permutations),
    enrichment('gene over 500 kb', track.longGenes, regions, markersByChrom, permutations),
    enrichment('centromere or telomere', track.gaps, regions, markersByChrom, permutations),
    enrichment('late-replication valley (ES)', track.lateReplicationValleysES ?? [],
               regions, markersByChrom, permutations),
    enrichment('late-replication valley (constitutive)',
               track.lateReplicationValleysConstitutive ?? [],
               regions, markersByChrom, permutations),
  ]
  // The five above are left on the ORIGINAL code path, unchanged and unindexed, so the figures
  // already reported for them stay reproducible byte for byte. The widened set is scored beside
  // them under the mode each track declares.
  const extra = track.extra ?? {}
  const modes = track.modes ?? {}
  const rest = Object.keys(extra).sort().map((name) => enrichmentBy(
    name, extra[name], regions, markersByChrom,
    modes[name] ?? 'binary', permutations,
  ))
  return [...five, ...rest]
}

/** Genes per Mb across a region, from the binned track. NaN when the region spans no bin. */
export function geneDensity(track: FeatureTrack, r: Region): number {
  const bins = track.geneDensityPerMb[r.chrom]
  if (!bins) return NaN
  let total = 0
  let n = 0
  for (let b = Math.floor(r.startBp / 1e6); b <= Math.floor(r.endBp / 1e6); b += 1) {
    total += bins[String(b)] ?? 0
    n += 1
  }
  return n ? total / n : NaN
}

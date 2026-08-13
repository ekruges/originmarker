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

/** Every track in one pass. Gaps are reported so artefact near a centromere is visible. */
export function scoreAll(
  track: FeatureTrack,
  regions: Region[],
  markersByChrom: Map<string, number[]>,
  permutations = DEFAULT_PERMUTATIONS,
): Enrichment[] {
  return [
    enrichment('common fragile site', track.fragile, regions, markersByChrom, permutations),
    enrichment('gene over 500 kb', track.longGenes, regions, markersByChrom, permutations),
    enrichment('centromere or telomere', track.gaps, regions, markersByChrom, permutations),
    enrichment('late-replication valley (ES)', track.lateReplicationValleysES ?? [],
               regions, markersByChrom, permutations),
    enrichment('late-replication valley (constitutive)',
               track.lateReplicationValleysConstitutive ?? [],
               regions, markersByChrom, permutations),
  ]
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

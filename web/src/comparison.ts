/**
 * The optional in-depth comparison: do these regions sit where the genome is fragile, or not?
 *
 * WHY THIS IS SEPARATE AND OPTIONAL. It answers a different question from the rest of the tool.
 * Everything else here asks WHOSE a change is; this asks whether the change is where the genome
 * breaks anyway. Those are not competing answers and one does not qualify the other, but a reader
 * who meets them in one undifferentiated report will treat the second as evidence about the first,
 * and it is not. So it runs on request, produces its own report, and is folded into the main one
 * only once someone has chosen to run it.
 *
 * WHAT IT IS FOR. A segmental loss at a common fragile site, a late-replication valley or an
 * assembly gap has a mundane explanation available: that part of the genome is hard to replicate
 * and hard to assay, and the same interval turns up across unrelated samples. A segmental loss with
 * no such coincidence has no such explanation and is more likely to be what it appears to be. That
 * is the whole contribution, and it is a prior on interpretation rather than a verdict.
 *
 * THE ONE THING IT CANNOT DO, STATED HERE BECAUSE IT IS THE MOST TEMPTING MISREADING. It says
 * nothing whatever about parental origin. The late-replicating fragile compartment is established
 * on BOTH parental genomes from the first cell cycle, so a region coinciding with it is no more
 * paternal than maternal. A result here must never be read as support for a parental call, and the
 * report says so in those words on every page that carries a number.
 *
 * THE NULL IS THE POINT. Regions can only be called where informative markers are, and marker
 * density tracks gene density, so a comparison against the genome at large reports an enrichment
 * for almost any feature. Every figure here is scored against intervals drawn on the same
 * chromosome carrying the same number of informative markers as the region being tested, which is
 * the null the shipped enrichment already uses; this module reads it, it does not reinvent it.
 */
import type { Enrichment, FeatureInterval, FeatureTrack, Region } from './features.ts'
import { DEFAULT_PERMUTATIONS, MIN_REGIONS, scoreAll } from './features.ts'

/** Human names and one line on what a coincidence with each feature would MEAN. */
export const FEATURE_MEANING: Record<string, { label: string; means: string }> = {
  'common fragile site': {
    label: 'Common fragile sites',
    means: 'these break under replication stress in normal cells, so a change here has a mundane '
      + 'explanation available that a change elsewhere does not',
  },
  'gene over 500 kb': {
    label: 'Long genes',
    means: 'genes over the length at which transcription and replication collide, which is the '
      + 'mechanism behind most common fragile sites',
  },
  'centromere or telomere': {
    label: 'Centromeres and telomeres',
    means: 'sequence the reference does not contain, so markers thin out and an interval can be '
      + 'called from the thinning rather than from the sample. Assembly gaps sit inside these',
  },
  'late-replication valley (ES)': {
    label: 'Late-replication valleys, embryonic stem',
    means: 'the latest-replicating POINTS in an embryonic stem line, the closest available cell '
      + 'type to this material. Points rather than domains: 4,843 of them spanning 5 Mb in total',
  },
  'late-replication valley (constitutive)': {
    label: 'Late-replication valleys, constitutive',
    means: 'late in at least 8 of 10 cell lines, so late irrespective of cell type',
  },
}

export const meaningFor = (feature: string): { label: string; means: string } =>
  FEATURE_MEANING[feature] ?? { label: feature, means: 'a genomic feature supplied with the track' }

/**
 * Significance at which a coincidence is called, and it is deliberately not 0.05.
 *
 * Every feature in the track is tested on the same regions, so the more features a track carries
 * the more likely one clears an uncorrected bound by chance. The threshold is divided by the number
 * of features actually tested, which is the plainest correction available and the one whose effect
 * a reader can check by counting the rows.
 */
export const ALPHA = 0.05

export interface FeatureComparison {
  feature: string
  label: string
  means: string
  /**
   * FRACTION of regions touching this feature, which is what the shipped enrichment reports.
   * Kept as a fraction rather than converted to a count, because the null is a fraction too and
   * mixing the two units is how a chart ends up comparing incomparable numbers.
   */
  observed: number
  /** Count, for a reader who wants the plain number beside the fraction. */
  observedCount: number
  /** What the matched null expects, and the middle 95% of that null. */
  expected: number
  nullLo: number
  nullHi: number
  p: number
  /** observed / expected. NaN where the null expects nothing. */
  fold: number
  /** True where p clears the corrected bound. */
  significant: boolean
  /** Named regions that touch this feature. */
  hits: string[]
  /**
   * The permutation distribution itself, 20 bins, so a report can DRAW it.
   *
   * The most honest single figure this analysis can produce, and it costs nothing because the
   * enrichment already computed it. A permutation p is uninterpretable without the distribution it
   * came from: a ratio near 1 reaches a small p when the null is tight, and a reader shown only the
   * p would call that a finding. Drawing the null makes that judgement available rather than
   * asking the reader to take the p on trust.
   */
  nullHist: { lo: number; hi: number; counts: number[] }
}

export type Verdict = 'feature-coincident' | 'independent' | 'underpowered'

export interface ComparisonResult {
  verdict: Verdict
  /** One sentence a reader can act on. */
  headline: string
  regions: number
  /**
   * The names of the regions compared, carried ON the result.
   *
   * So a report cannot draw a grid against a different set of regions from the one the enrichment
   * was computed on. Passing them alongside invites exactly that, and the grid is the figure a
   * reader checks the enrichment against, so a mismatch there is worse than no grid.
   */
  regionNames: string[]
  permutations: number
  features: FeatureComparison[]
  /** Generated, so the report's methods always describe the run that produced it. */
  methods: string
  /** Carried on every output. The single most tempting misreading of this analysis. */
  caveat: string
}

export const NOT_ABOUT_PARENTS =
  'Positional only. This says nothing about which parent a change came from: the late-replicating '
  + 'fragile compartment is established on both parental genomes from the first cell cycle, so a '
  + 'region coinciding with it is no more paternal than maternal. A result here must never be read '
  + 'as support for a parental call.'

/**
 * Run the comparison and reduce it to something a reader can act on.
 *
 * The verdict is deliberately coarse. The underlying numbers are per-feature and are all shown, but
 * a reader wants to know one thing first: is there a mundane explanation available for where these
 * regions landed. Three answers cover it, and the third is not a failure.
 */
/** Join a list the way a sentence does, so five coincidences do not read as one long chain. */
export const listOf = (xs: readonly string[]): string => (xs.length <= 1 ? (xs[0] ?? '')
  : `${xs.slice(0, -1).join(', ')} and ${xs[xs.length - 1]}`)

export function compare(
  track: FeatureTrack,
  regions: readonly Region[],
  markersByChrom: Map<string, number[]>,
  opts: {
    permutations?: number; alpha?: number; minRegions?: number
    /** Names in the same order as `regions`, carried onto the result for the grid. */
    regionNames?: readonly string[]
  } = {},
): ComparisonResult {
  const permutations = opts.permutations ?? DEFAULT_PERMUTATIONS
  const minRegions = opts.minRegions ?? MIN_REGIONS
  const raw: Enrichment[] = regions.length
    ? scoreAll(track, regions as Region[], markersByChrom, permutations)
    : []
  const usable = raw.filter((e) => Number.isFinite(e.p))
  // Corrected for the number of features actually tested, not the number the track could carry.
  const alpha = (opts.alpha ?? ALPHA) / Math.max(1, usable.length)

  const features: FeatureComparison[] = usable.map((e) => {
    const m = meaningFor(e.feature)
    return {
      feature: e.feature,
      label: m.label,
      means: m.means,
      observed: e.observed,
      observedCount: Math.round(e.observed * e.regions),
      expected: e.nullMean,
      // The 2.5th and 97.5th percentiles of the null the permutation actually drew, so the figure
      // is drawn against its own distribution rather than a normal approximation of it.
      nullLo: e.nullQuantiles[0],
      nullHi: e.nullQuantiles[2],
      p: e.p,
      fold: e.nullMean > 0 ? e.observed / e.nullMean : NaN,
      significant: e.p < alpha,
      hits: e.hits ?? [],
      nullHist: e.nullHist,
    }
  }).sort((a, b) => a.p - b.p)

  const hit = features.filter((f) => f.significant)
  const verdict: Verdict = regions.length < minRegions ? 'underpowered'
    : hit.length ? 'feature-coincident' : 'independent'

  const headline = verdict === 'underpowered'
    ? `${regions.length} region${regions.length === 1 ? '' : 's'} is too few to compare against `
      + `anything: the null needs at least ${minRegions} to have a shape. No conclusion either way.`
    : verdict === 'feature-coincident'
      ? `These regions coincide with ${listOf(hit.map((f) => f.label.toLowerCase()))}. A `
        + 'mundane explanation is available for where they landed, so treat them more cautiously '
        + 'than their detection statistics alone suggest.'
      : 'These regions do not coincide with any feature tested. Nothing here offers an alternative '
        + 'explanation for where they landed, which leaves their detection statistics standing.'

  return {
    verdict,
    headline,
    regions: regions.length,
    regionNames: [...(opts.regionNames
      ?? regions.map((r) => `chr${r.chrom} ${(r.startBp / 1e6).toFixed(1)}-${(r.endBp / 1e6).toFixed(1)}Mb`))],
    permutations,
    features,
    methods: methodsText(regions.length, permutations, usable.length, alpha),
    caveat: NOT_ABOUT_PARENTS,
  }
}

/**
 * The methods paragraph, generated from the run rather than written once and left to drift.
 *
 * Every number in it comes from the analysis that produced the figures beside it, so a report
 * cannot describe a permutation count or a correction it did not use.
 */
export function methodsText(
  regions: number, permutations: number, tested: number, alpha: number,
): string {
  return `Each detected region was tested against ${tested} genomic feature set`
    + `${tested === 1 ? '' : 's'} for positional coincidence. The null is the part that decides `
    + 'whether any of this means anything: a region can only be called where informative markers '
    + 'are, and marker density tracks gene density, so a comparison against the genome at large '
    + 'reports an enrichment for almost any feature. Instead, each region was compared against '
    + `intervals drawn on the SAME chromosome carrying the SAME number of informative markers as `
    + `the region itself, resampled ${permutations.toLocaleString()} times. Reported for each `
    + 'feature: the number of the '
    + `${regions} region${regions === 1 ? '' : 's'} overlapping it, what that matched null expects, `
    + 'the middle 95% of the null, and the resulting p value. Significance is called at '
    + `${alpha.toExponential(2)}, which is 0.05 divided by the ${tested} feature set`
    + `${tested === 1 ? '' : 's'} tested, since testing many features on one set of regions is `
    + 'otherwise a way of finding one that fits.'
}

// ---------------------------------------------------------------------------------------------
// Chart data. Shapes only, so the same numbers drive the screen and the report.

export interface EnrichmentBar {
  label: string
  observed: number
  expected: number
  lo: number
  hi: number
  /** Axis maximum this bar should be drawn against, shared across the set. */
  axisMax: number
  significant: boolean
  p: number
}

/**
 * Observed against its own null, one row per feature.
 *
 * A SHARED AXIS ACROSS ROWS, which is the only reason the picture is readable. Per-row axes make
 * every feature look equally enriched, because each fills its own width, and that is the single
 * easiest way to draw a chart that says the opposite of its data.
 */
export function enrichmentBars(c: ComparisonResult): EnrichmentBar[] {
  const axisMax = Math.max(
    1, ...c.features.map((f) => Math.max(f.observed, f.nullHi, f.expected)),
  )
  return c.features.map((f) => ({
    label: f.label,
    observed: f.observed,
    expected: f.expected,
    lo: f.nullLo,
    hi: f.nullHi,
    axisMax,
    significant: f.significant,
    p: f.p,
  }))
}

export interface RegionCell { region: string; feature: string; touched: boolean }

/**
 * Which region touches which feature, as a grid.
 *
 * The plainest possible answer to "overlap or lack thereof", and the one a reader checks the
 * enrichment against: a feature can clear its null on two regions out of twenty, and the grid is
 * where that becomes visible rather than being hidden inside a p value.
 */
export function regionGrid(c: ComparisonResult, regionNames: readonly string[]): RegionCell[] {
  const out: RegionCell[] = []
  for (const r of regionNames) {
    for (const f of c.features) out.push({ region: r, feature: f.label, touched: f.hits.includes(r) })
  }
  return out
}

/** Fold enrichment per feature, for the one chart that puts every feature on a single scale. */
export const foldChart = (c: ComparisonResult): { label: string; fold: number; p: number }[] =>
  c.features.filter((f) => Number.isFinite(f.fold))
    .map((f) => ({ label: f.label, fold: f.fold, p: f.p }))
    .sort((a, b) => b.fold - a.fold)

/**
 * Which of a run's changes this analysis can be asked about.
 *
 * SEGMENTS AND TAXONOMY FINDINGS BOTH. A copy-neutral event and an isodisomy are chromosomal
 * changes exactly as a deletion is, and the question is identical for all of them: does it sit
 * where the genome breaks anyway. Taking only the segments left the panel reporting "0 regions" on
 * runs whose changes were all findings.
 *
 * Genome-wide findings are dropped, because a triploidy or a chaotic genome has no interval to
 * compare against anything: they are properties of the whole array rather than places in it.
 */
export function comparisonRegions(run: {
  segments?: readonly { chrom: string }[]
  findings?: readonly { chrom: string; startBp: number; endBp: number }[]
  segmentCoords?: unknown
}, coordsOf?: (sg: unknown) => { start: number; end: number }): {
  region: Region; name: string
}[] {
  const out: { region: Region; name: string }[] = []
  const name = (chrom: string, a: number, b: number) =>
    `chr${chrom} ${(a / 1e6).toFixed(1)}-${(b / 1e6).toFixed(1)}Mb`
  for (const sg of run.segments ?? []) {
    const co = coordsOf ? coordsOf(sg) : null
    if (!co) continue
    out.push({
      region: { chrom: sg.chrom, startBp: co.start, endBp: co.end },
      name: name(sg.chrom, co.start, co.end),
    })
  }
  for (const f of run.findings ?? []) {
    if (f.chrom === 'genome' || !(f.endBp > f.startBp)) continue
    out.push({
      region: { chrom: f.chrom, startBp: f.startBp, endBp: f.endBp },
      name: name(f.chrom, f.startBp, f.endBp),
    })
  }
  return out
}

/**
 * Accept a track whose intervals are TUPLES as well as one whose intervals are objects.
 *
 * WHY THIS EXISTS, AND IT IS NOT DEFENSIVE PROGRAMMING. The shipped hg19_features.json stores each
 * interval as ["1", 61300000, 84900000, "FRA1B"], while FeatureInterval is an object and every
 * consumer reads `f.chrom` and `f.startBp`. On a tuple those are undefined, so `touches` compared
 * undefined against a chromosome name and returned false for every marker of every region. The
 * enrichment therefore reported ZERO overlap with every feature on every real run, whatever the
 * regions actually overlapped, and reported it with a p value.
 *
 * Nothing failed, which is why it survived: the check file builds objects, so the tests exercised a
 * shape the application never sees. The generator emits objects too, so the file was compacted
 * somewhere between the two.
 *
 * Both shapes are accepted rather than one being declared correct, because that makes the fix
 * independent of regenerating a data file, and a track from either path now works.
 */
export function normaliseTrack(raw: unknown): FeatureTrack | null {
  if (!raw || typeof raw !== 'object') return null
  const t = raw as Record<string, unknown>
  const list = (v: unknown): FeatureInterval[] => {
    if (!Array.isArray(v)) return []
    return v.map((f) => (Array.isArray(f)
      ? { chrom: String(f[0]), startBp: Number(f[1]), endBp: Number(f[2]), name: String(f[3] ?? '') }
      : f as FeatureInterval))
      .filter((f) => f && typeof f.chrom === 'string'
        && Number.isFinite(f.startBp) && Number.isFinite(f.endBp))
  }
  return {
    build: typeof t.build === 'string' ? t.build : 'unknown',
    fragile: list(t.fragile),
    longGenes: list(t.longGenes),
    gaps: list(t.gaps) as FeatureTrack['gaps'],
    geneDensityPerMb: (t.geneDensityPerMb ?? {}) as FeatureTrack['geneDensityPerMb'],
    lateReplicationValleysES: list(t.lateReplicationValleysES),
    lateReplicationValleysConstitutive: list(t.lateReplicationValleysConstitutive),
    ...(t.extra ? { extra: Object.fromEntries(
      Object.entries(t.extra as Record<string, unknown>).map(([k, v]) => [k, list(v)]),
    ) } : {}),
    ...(t.modes ? { modes: t.modes as never } : {}),
  } as FeatureTrack
}

/**
 * The permutation null as a drawable histogram, with the observation placed in it.
 *
 * Returns the bins, the bar the observation falls in, and where it sits along the axis, so the
 * screen and the report draw the same picture from the same numbers.
 */
export interface NullHistogram {
  label: string
  bins: number[]
  lo: number
  hi: number
  /** Index of the bin the observation falls in, or -1 where it lands outside the null entirely. */
  observedBin: number
  /** Position of the observation on [0, 1] across the axis, clamped so it stays drawable. */
  observedAt: number
  observed: number
  p: number
  significant: boolean
}

export function nullHistograms(c: ComparisonResult): NullHistogram[] {
  return c.features.filter((f) => f.nullHist?.counts?.length).map((f) => {
    const { lo, hi, counts } = f.nullHist
    const span = hi - lo
    const at = span > 0 ? (f.observed - lo) / span : 0
    // OUTSIDE THE NULL IS THE INTERESTING CASE, so it is reported rather than clamped away: an
    // observation the permutation never reached is the strongest result this analysis produces.
    const idx = span > 0 ? Math.floor(at * counts.length) : -1
    return {
      label: f.label,
      bins: counts,
      lo,
      hi,
      observedBin: idx >= 0 && idx < counts.length ? idx : -1,
      observedAt: Math.max(0, Math.min(1, at)),
      observed: f.observed,
      p: f.p,
      significant: f.significant,
    }
  })
}

/**
 * Every region a whole run found, pooled, with its sample name attached.
 *
 * POOLED BECAUSE A SINGLE SAMPLE IS NOT A DATASET. The matched null needs at least five regions to
 * have a shape, and one chip routinely contributes fewer, so a per-sample comparison answers "no
 * conclusion" on almost every card while the run as a whole carries plenty. The question is a
 * cohort question in any case: whether the changes THIS EXPERIMENT found sit where the genome
 * breaks anyway is not a property of one embryo.
 *
 * Names carry the sample so a coincidence can be traced back to the array that produced it.
 */
export function pooledRegions<E extends {
  result?: { segments?: readonly { chrom: string }[]
    findings?: readonly { chrom: string; startBp: number; endBp: number }[] }
  file?: { name: string }
}>(entries: readonly E[], coordsOf: (sg: unknown) => { start: number; end: number }): {
  region: Region; name: string
}[] {
  const out: { region: Region; name: string }[] = []
  for (const e of entries) {
    if (!e.result) continue
    const tag = (e.file?.name ?? '').replace(/\.[^.]+$/, '').slice(0, 14)
    for (const r of comparisonRegions(e.result, coordsOf)) {
      out.push({ region: r.region, name: tag ? `${tag} ${r.name}` : r.name })
    }
  }
  return out
}

/**
 * One marker map for the pooled run.
 *
 * The null is drawn from where informative markers are, so pooling regions across samples means
 * pooling the positions they could have been called at. Positions are unioned rather than
 * concatenated: two samples of the same platform share almost every marker, and counting each twice
 * would inflate the density the null is matched on.
 */
export function pooledMarkers<E extends { result?: { markerPositions?: Map<string, number[]> } }>(
  entries: readonly E[],
): Map<string, number[]> {
  const byChrom = new Map<string, Set<number>>()
  for (const e of entries) {
    for (const [c, ps] of e.result?.markerPositions ?? []) {
      if (!byChrom.has(c)) byChrom.set(c, new Set())
      const set = byChrom.get(c)!
      for (const p of ps) set.add(p)
    }
  }
  return new Map([...byChrom].map(([c, set]) => [c, [...set].sort((a, b) => a - b)]))
}

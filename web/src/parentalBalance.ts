/**
 * Whether chromosomal change falls on the two parental genomes DIFFERENTIALLY or EQUALLY.
 *
 * This is the question the tool exists to answer, and until now nothing in it aggregated: every
 * event got its own row and nothing ever counted the rows.
 *
 * WHY UNIPARENTAL SAMPLES ARE THE CLEAN COMPARISON. In a biparental sample, deciding whose a
 * change is means calling it per event, and detection power differs by event class and material,
 * so a raw count measures the assay as much as the biology. A gynogenetic embryo carries maternal
 * DNA and no paternal complement, and an androgenetic one is the mirror. Each genome has exactly
 * one parent in it, so every change in it is that parent's by construction and there is no
 * per-event assignment to be biased. Comparing the two groups compares the parents directly.
 *
 * WHAT IS STILL CONFOUNDED, and it is not optional to report. Power to see a change at all depends
 * on how many informative markers a sample has and how much of it was called, and those vary
 * enormously between samples: 580,133 informative on one array against 282,055 on another in the
 * same run. If one group is systematically better measured than the other, a difference in counts
 * is a difference in looking. So counts are normalised per informative marker, and the imbalance
 * in power between the groups is measured and reported beside the result rather than left for a
 * reader to wonder about.
 */

/**
 * Fewest genomes in EACH group before a difference can be reported at all.
 *
 * FIVE, NOT THREE. At 3 against 3 the exact two-sided minimum achievable p is 2/C(6,3) = 0.100,
 * which cannot clear alpha 0.05 by any margin, so the old floor permitted a test incapable of
 * producing a finding. Five is the smallest size whose floor clears a Bonferroni threshold for
 * three classes.
 */
export const MIN_PER_GROUP = 5

/**
 * Group size below which a result is exploratory rather than a headline.
 *
 * Measured on this corpus's own dispersion, negative binomial k = 0.81: power at a 3x effect is
 * 0.14 at n=4, 0.38 at n=8, and 0.82 at n=20. Twenty is the smallest group with 80% power at 3x.
 * Between the hard floor and this, the test runs and the report must say what it could have found.
 */
export const REPORTING_PER_GROUP = 20
/** Label permutations for the null. */
export const DEFAULT_PERMUTATIONS = 10_000
/** Significance for the difference. */
export const ALPHA = 0.05
/**
 * How lopsided the two groups' detection power may be before the comparison is not trusted.
 *
 * A ratio of median informative markers. Beyond this the groups are not comparable: the better
 * measured group would show more events whatever the biology, and reporting a difference from it
 * would be reporting the sequencing.
 */
export const MAX_POWER_SKEW = 1.4

export type BalanceSample = {
  name: string
  /** The single parent this genome carries, or null where the sample has two or is unclear. */
  parent: 'maternal' | 'paternal' | null
  /**
   * Which group this sample WOULD have joined had it been classified, so an exclusion that is
   * correlated with quality is visible. Set on excluded samples; ignored on included ones.
   */
  declaredParent?: 'maternal' | 'paternal'
  originClass: string
  /** Markers that could have carried evidence. Reported, but NOT used as a denominator. */
  informative: number
  /**
   * The array's own explainable-noise ceiling, which is what predicts how many false regions it
   * will produce. This is the quantity the groups must be matched on, not the marker count.
   */
  explainable: number
  material: string
  events: readonly { cls: string; chrom: string; startBp: number; endBp: number }[]
}

export type GroupStat = {
  parent: 'maternal' | 'paternal'
  samples: number
  /** Distinct events after merging, summed over the group. */
  events: number
  eventsPerSample: number
  /** Genomes carrying ANY change, which is the headline unit. */
  carryingAny: number
  /** The value the rank test actually compares. */
  medianEvents: number
  /** Median artefact propensity, the quantity the groups must be matched on. */
  medianNoiseCeiling: number
  /** Reported for context only. Deliberately not a denominator. */
  medianInformative: number
}

export type ClassStat = {
  cls: string
  maternal: number
  paternal: number
  maternalRate: number
  paternalRate: number
  /** Maternal rate over paternal rate. Above 1 means maternal genomes carry more. */
  fold: number
  p: number
  significant: boolean
  /**
   * True where this class could NOT have been flagged at any effect size, because the group sizes
   * do not admit a small enough p once corrected. Without this a reader takes a null for evidence
   * of no difference when it is evidence of not enough genomes.
   */
  underpowered: boolean
}

export type BalanceResult = {
  verdict: 'differential' | 'equal' | 'underpowered' | 'not-comparable'
  headline: string
  groups: GroupStat[]
  byClass: ClassStat[]
  fold: number
  p: number
  alpha: number
  permutations: number
  /**
   * The smallest p these group sizes can produce, from the number of distinct label assignments.
   *
   * Eight genomes split four and four give 70 assignments, so nothing can be reported below about
   * 0.029, and a corrected per-class threshold sits under that. The test is then incapable of
   * flagging anything however large the effect, and saying so is the difference between "no
   * difference was found" and "no difference could have been found".
   */
  minAchievableP: number
  power: {
    maternalMedianNoiseCeiling: number
    paternalMedianNoiseCeiling: number
    maternalMedianInformative: number
    paternalMedianInformative: number
    skew: number
    withinTolerance: boolean
    /** Whether the two groups lost samples at comparable rates before the comparison began. */
    exclusion: { maternal: number; paternal: number; balanced: boolean }
  }
  excluded: { name: string; why: string }[]
  methods: string
}

/**
 * The smallest two-sided p a label permutation can reach with these group sizes.
 *
 * One over the number of distinct assignments, doubled for the mirror, with the add-one the
 * estimator uses. Computed rather than approximated because it is the number that decides whether
 * a null is informative.
 */
export const minAchievableP = (nA: number, nB: number, permutations: number): number => {
  if (nA < 1 || nB < 1) return 1
  const n = nA + nB
  let choose = 1
  for (let i = 1; i <= Math.min(nA, nB); i += 1) choose = (choose * (n - i + 1)) / i
  // THE EXACT TWO-SIDED MINIMUM: the most extreme split and its mirror, over all distinct label
  // assignments. The previous form added the estimator's add-one to both parts and returned 0.1429
  // at 3 vs 3 where the exact value is 0.100, and 0.0423 at 4 vs 4 where it is 0.0286. It was
  // conservative rather than wrong, so it over-declared underpowered.
  const exact = 2 / choose
  // A permutation run cannot resolve below its own resolution either.
  return Math.min(1, Math.max(exact, 1 / (permutations + 1)))
}

const median = (xs: number[]): number => {
  if (!xs.length) return NaN
  const s = [...xs].sort((a, b) => a - b)
  const m = s.length >> 1
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}

/** Same generator the feature enrichment uses, so a report is reproducible from its seed. */
const rng = (seed: number) => () => {
  seed = (seed * 1664525 + 1013904223) >>> 0
  return seed / 4294967296
}

/**
 * Distinct events, by joining regions that sit within `joinBp` of each other on one chromosome.
 *
 * WINDOW COUNT IS NOT EVENT COUNT. A sliding detector cuts one biological change into as many rows
 * as its windows, so a raw region count is defensible only as "windows". Measured on this
 * project's own 214-cell corpus, merging removes 17.9% of rows and reorders samples barely at all
 * (Spearman 0.966 against the raw count), so this is a modest correction rather than a rescue. It
 * is made because the merged count is defensible as distinct events and the raw one is not.
 */
export const MERGE_JOIN_BP = 1_000_000

export function mergeEvents(
  events: readonly { cls: string; chrom: string; startBp: number; endBp: number }[],
  joinBp = MERGE_JOIN_BP,
): { cls: string; chrom: string; startBp: number; endBp: number }[] {
  const byKey = new Map<string, { cls: string; chrom: string; startBp: number; endBp: number }[]>()
  for (const e of events) {
    const k = `${e.cls}|${e.chrom}`
    const cur = byKey.get(k)
    if (cur) cur.push(e)
    else byKey.set(k, [e])
  }
  const out: { cls: string; chrom: string; startBp: number; endBp: number }[] = []
  for (const group of byKey.values()) {
    const sorted = [...group].sort((a, b) => a.startBp - b.startBp)
    let cur = { ...sorted[0] }
    for (let i = 1; i < sorted.length; i += 1) {
      const e = sorted[i]
      if (e.startBp - cur.endBp <= joinBp) cur.endBp = Math.max(cur.endBp, e.endBp)
      else { out.push(cur); cur = { ...e } }
    }
    out.push(cur)
  }
  return out
}

/**
 * The rank-sum of per-sample values, which is the statistic the permutation actually supports.
 *
 * A POOLED RATIO IS NOT EXCHANGEABLE-SAFE. `sum(events)/sum(markers)` over a group is not a
 * symmetric function of the labels, because each sample's denominator travels with it, so the
 * permutation null does not have the right size. Measured under a true null at 7 against 5: the
 * pooled statistic runs at 0.044 when the groups are matched but 0.164 at a marker skew of 1.5,
 * more than three times nominal, while the rank statistic holds at 0.049 and 0.077. The rank is
 * also robust to the one disturbed genome that dominates a pooled count.
 */
export const rankSum = (values: readonly number[], inGroupA: readonly boolean[]): number => {
  const idx = values.map((_, i) => i).sort((a, b) => values[a] - values[b])
  const rank = new Array<number>(values.length)
  // Midranks, so ties do not manufacture a difference.
  let i = 0
  while (i < idx.length) {
    let j = i
    while (j + 1 < idx.length && values[idx[j + 1]] === values[idx[i]]) j += 1
    const mid = (i + j) / 2 + 1
    for (let k = i; k <= j; k += 1) rank[idx[k]] = mid
    i = j + 1
  }
  let sum = 0
  for (let k = 0; k < values.length; k += 1) if (inGroupA[k]) sum += rank[k]
  return sum
}

/**
 * Distinct events in one sample, which is the per-sample value every test below ranks.
 *
 * NO MARKER DENOMINATOR. Dividing by informative markers was measured to AMPLIFY the confounder it
 * was meant to remove: across this project's arrays, informative markers span 1.82x while the
 * array's own artefact propensity spans 19.7x, and the two are NEGATIVELY correlated (Spearman
 * -0.533, p = 0.005) because cleaner arrays also carry more markers. Comparing the cleanest against
 * the noisiest arrays, the bias runs 5.89x on raw counts and 7.06x once divided by markers. The
 * quantity that predicts how many false regions an array produces is its own explainable-noise
 * ceiling, not its marker count, so counts stay raw and the noise term is reported beside them.
 */
export const eventCount = (s: BalanceSample): number => mergeEvents(s.events).length

/** Whether this genome carries any change at all, the unit immune to slicing and to sensitivity. */
export const carriesAny = (s: BalanceSample): boolean => s.events.length > 0

/**
 * Permutation p on the rank-sum of per-sample values.
 *
 * The labels are shuffled and each sample keeps its own value, so the null is exactly "the
 * parental label is unrelated to the value" and nothing about the values themselves is assumed.
 */
function rankTest(
  values: readonly number[], inA: readonly boolean[], permutations: number, seed: number,
): number {
  const nA = inA.filter(Boolean).length
  if (!nA || nA === values.length) return 1
  const observed = rankSum(values, inA)
  const centre = (nA * (values.length + 1)) / 2
  const target = Math.abs(observed - centre)
  const next = rng(seed)
  let atLeast = 0
  for (let i = 0; i < permutations; i += 1) {
    const shuffled = [...inA]
    for (let j = shuffled.length - 1; j > 0; j -= 1) {
      const k = Math.floor(next() * (j + 1))
      ;[shuffled[j], shuffled[k]] = [shuffled[k], shuffled[j]]
    }
    if (Math.abs(rankSum(values, shuffled) - centre) >= target - 1e-12) atLeast += 1
  }
  return (atLeast + 1) / (permutations + 1)
}

export function parentalBalance(
  samples: readonly BalanceSample[],
  opts: { permutations?: number; alpha?: number; seed?: number; minPerGroup?: number } = {},
): BalanceResult {
  const permutations = opts.permutations ?? DEFAULT_PERMUTATIONS
  const alpha = opts.alpha ?? ALPHA
  const minPerGroup = opts.minPerGroup ?? MIN_PER_GROUP
  const seed = opts.seed ?? 7

  const excluded: { name: string; why: string }[] = []
  const usable: BalanceSample[] = []
  for (const s of samples) {
    if (!s.parent) {
      excluded.push({ name: s.name, why: `${s.originClass}, so it carries two candidate parents `
        + 'and cannot be attributed to one without a per-event call' })
      continue
    }
    usable.push(s)
  }

  const mats = usable.filter((s) => s.parent === 'maternal')
  const pats = usable.filter((s) => s.parent === 'paternal')
  const mkGroup = (parent: 'maternal' | 'paternal', rows: BalanceSample[]): GroupStat => ({
    parent,
    samples: rows.length,
    events: rows.reduce((a, s) => a + eventCount(s), 0),
    eventsPerSample: rows.length
      ? rows.reduce((a, s) => a + eventCount(s), 0) / rows.length : NaN,
    carryingAny: rows.filter(carriesAny).length,
    medianEvents: median(rows.map(eventCount)),
    medianNoiseCeiling: median(rows.map((s) => s.explainable)),
    medianInformative: median(rows.map((s) => s.informative)),
  })
  const groups = [mkGroup('maternal', mats), mkGroup('paternal', pats)]

  // POWER IS THE ARTEFACT PROPENSITY, NOT THE MARKER COUNT. See eventCount.
  const mNoise = groups[0].medianNoiseCeiling
  const pNoise = groups[1].medianNoiseCeiling
  const skew = Number.isFinite(mNoise) && Number.isFinite(pNoise) && Math.min(mNoise, pNoise) > 0
    ? Math.max(mNoise, pNoise) / Math.min(mNoise, pNoise) : NaN
  const withinTolerance = Number.isFinite(skew) && skew <= MAX_POWER_SKEW

  // AND THE EXCLUSIONS THEMSELVES MAY BE QUALITY-CORRELATED, which the skew above cannot see
  // because it measures only the samples that survived. In this project's own run the three arrays
  // that fell out as unclear were the three worst, and all three were paternal, so the exclusion
  // silently repaired the measured skew by dropping one group's worst members.
  const droppedFrom = (parent: 'maternal' | 'paternal') =>
    samples.filter((s) => !s.parent && s.declaredParent === parent).length
  const exclusion = {
    maternal: droppedFrom('maternal'),
    paternal: droppedFrom('paternal'),
    balanced: Math.abs(droppedFrom('maternal') - droppedFrom('paternal')) <= 1,
  }

  const power = {
    maternalMedianNoiseCeiling: mNoise,
    paternalMedianNoiseCeiling: pNoise,
    maternalMedianInformative: groups[0].medianInformative,
    paternalMedianInformative: groups[1].medianInformative,
    skew,
    withinTolerance,
    exclusion,
  }

  const floorP = minAchievableP(mats.length, pats.length, permutations)
  const methods = `Uniparental samples only, ${usable.length} of ${samples.length}: a genome with `
    + 'one parent in it attributes every change in it to that parent by construction, so no '
    + 'per-event call is needed and no per-event detection bias can enter. The countable unit is '
    + `distinct events, regions merged within ${(MERGE_JOIN_BP / 1e6).toFixed(0)} Mb, because a `
    + 'sliding detector cuts one change into as many rows as it has windows. Counts are NOT '
    + 'normalised by marker count: marker count does not predict how many false regions an array '
    + "produces and dividing by it amplifies the bias, so the array's own explainable-noise "
    + 'ceiling is reported beside the counts instead. The test is the rank-sum of per-sample '
    + `counts under ${permutations.toLocaleString()} label permutations, two-sided, against alpha `
    + `${alpha}; a pooled ratio is not a symmetric function of the labels and does not hold its `
    + 'size. The headline is the fraction of genomes carrying any change, which is the only unit '
    + 'immune to both slicing and detector sensitivity.'
    + ' WHAT THIS CANNOT ANSWER: gynogenetic and androgenetic conceptuses are not a sample of '
    + 'embryos, they are a sample of fertilisation failures, and the two classes arise by '
    + 'different mechanisms with different replication histories, so a difference between them is '
    + 'confounded with the mechanism that produced the class. Where the material is dissected '
    + 'pronuclei rather than conceptuses, it is uniparental by dissection and has been through '
    + 'neither syngamy nor a mitosis, which makes a run on it a POSITIVE CONTROL for the '
    + 'detection step and not a result about parental balance in embryos at all.'

  const base = {
    groups, power, excluded, methods, alpha, permutations, minAchievableP: floorP,
    byClass: [] as ClassStat[], fold: NaN, p: NaN,
  }

  if (mats.length < minPerGroup || pats.length < minPerGroup) {
    return {
      ...base,
      verdict: 'underpowered',
      headline: `${mats.length} maternal and ${pats.length} paternal genome`
        + `${pats.length === 1 ? '' : 's'} is too few to compare: each group needs at least `
        + `${minPerGroup}, below which the smallest reachable p is ${floorP.toFixed(4)}. `
        + 'No conclusion either way.',
    }
  }
  if (!withinTolerance) {
    return {
      ...base,
      verdict: 'not-comparable',
      headline: 'The two groups are not equally prone to artefact: their median explainable-noise '
        + `ceilings differ by ${Number.isFinite(skew) ? skew.toFixed(2) : '?'}x, past the `
        + `${MAX_POWER_SKEW}x this comparison allows. The noisier group would carry more detected `
        + 'change whatever the biology. No conclusion either way.',
    }
  }
  if (!exclusion.balanced) {
    return {
      ...base,
      verdict: 'not-comparable',
      headline: `Samples were excluded unevenly: ${exclusion.maternal} maternal against `
        + `${exclusion.paternal} paternal. A quality-correlated exclusion repairs the measured `
        + 'balance of the survivors by removing one group\'s worst members, so the comparison is '
        + 'between groups that were filtered differently. No conclusion either way.',
    }
  }

  const inA = usable.map((s) => s.parent === 'maternal')
  const counts = usable.map(eventCount)
  const p = rankTest(counts, inA, permutations, seed)
  const mMed = groups[0].medianEvents
  const pMed = groups[1].medianEvents
  const fold = pMed > 0 ? mMed / pMed : NaN

  const classes = [...new Set(usable.flatMap((s) => s.events.map((e) => e.cls)))].sort()
  const classAlpha = alpha / Math.max(1, classes.length)
  const byClass: ClassStat[] = classes.map((cls) => {
    const per = usable.map((s) => mergeEvents(s.events.filter((e) => e.cls === cls)).length)
    const cp = rankTest(per, inA, permutations, seed + 1)
    const mm = median(per.filter((_, i) => inA[i]))
    const pp = median(per.filter((_, i) => !inA[i]))
    return {
      cls,
      maternal: per.filter((_, i) => inA[i]).reduce((a, x) => a + x, 0),
      paternal: per.filter((_, i) => !inA[i]).reduce((a, x) => a + x, 0),
      maternalRate: mm,
      paternalRate: pp,
      fold: pp > 0 ? mm / pp : NaN,
      p: cp,
      significant: cp < classAlpha,
      underpowered: floorP >= classAlpha,
    }
  }).sort((a, b) => a.p - b.p)

  const differential = p < alpha
  const hits = byClass.filter((c) => c.significant)
  const exploratory = mats.length < REPORTING_PER_GROUP || pats.length < REPORTING_PER_GROUP
  const anyLine = `${groups[0].carryingAny} of ${groups[0].samples} maternal genomes carry any `
    + `change, against ${groups[1].carryingAny} of ${groups[1].samples} paternal.`
  const scope = exploratory
    ? ` EXPLORATORY: with ${Math.min(mats.length, pats.length)} genomes in the smaller group this `
      + `has about 80% power only at very large effects; ${REPORTING_PER_GROUP} per group is where `
      + 'a 3x effect is found reliably.'
    : ''

  return {
    ...base,
    verdict: differential ? 'differential' : floorP >= alpha ? 'underpowered' : 'equal',
    byClass,
    fold,
    p,
    headline: differential
      ? `Chromosomal change falls on the two parental genomes UNEQUALLY: median `
        + `${mMed} distinct events on maternal genomes against ${pMed} on paternal, `
        + `p = ${p.toFixed(4)}.${hits.length ? ` Carried by ${hits.map((c) => c.cls).join(', ')}.`
          : ''} ${anyLine}${scope}`
      : floorP >= alpha
        ? `No difference could have been detected at this size: ${mats.length} maternal and `
          + `${pats.length} paternal genomes admit no p below ${floorP.toFixed(4)}. This is not `
          + `evidence that the genomes are alike. ${anyLine}`
        : 'Chromosomal change falls on the two parental genomes EQUALLY, within what this many '
          + `genomes can resolve: median ${mMed} distinct events on maternal against ${pMed} on `
          + `paternal, p = ${p.toFixed(4)}. An equal result is a result, but it is not the same as `
          + `showing there is no difference. ${anyLine}${scope}`,
  }
}

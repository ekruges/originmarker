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

/** Fewest samples in EACH group before a difference means anything. */
export const MIN_PER_GROUP = 3
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
export const MAX_POWER_SKEW = 1.5

export type BalanceSample = {
  name: string
  /** The single parent this genome carries, or null where the sample has two or is unclear. */
  parent: 'maternal' | 'paternal' | null
  originClass: string
  /** Markers that could have carried evidence. The denominator, and the confounder. */
  informative: number
  material: string
  events: readonly { cls: string; chrom: string; startBp: number; endBp: number }[]
}

export type GroupStat = {
  parent: 'maternal' | 'paternal'
  samples: number
  events: number
  eventsPerSample: number
  /** The comparable number: events per 100,000 informative markers. */
  ratePer100k: number
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
    maternalMedianInformative: number
    paternalMedianInformative: number
    skew: number
    withinTolerance: boolean
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
  const n = nA + nB
  if (nA < 1 || nB < 1) return 1
  let choose = 1
  for (let i = 1; i <= Math.min(nA, nB); i += 1) choose = (choose * (n - i + 1)) / i
  const distinct = Math.min(choose, permutations)
  return (Math.min(2, distinct) + 1) / (Math.min(permutations, distinct) + 1)
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

/** Events per 100,000 informative markers, pooled over a group. */
const rateOf = (rows: readonly BalanceSample[], pick?: (cls: string) => boolean): number => {
  let ev = 0
  let inf = 0
  for (const s of rows) {
    inf += s.informative
    ev += pick ? s.events.filter((e) => pick(e.cls)).length : s.events.length
  }
  return inf > 0 ? (ev / inf) * 100_000 : NaN
}

/**
 * The difference between the groups, against a null built by shuffling which genome is whose.
 *
 * Permuting the LABELS rather than assuming a distribution: counts per sample are neither normal
 * nor Poisson here, since one disturbed genome contributes many correlated regions, and a test
 * that assumed either would be confident for the wrong reason. Shuffling the labels keeps each
 * sample's own event count and marker count together and only asks whether attaching them to
 * maternal rather than paternal was what produced the difference.
 */
function permutedP(
  rows: readonly BalanceSample[], pick: ((cls: string) => boolean) | undefined,
  permutations: number, seed: number,
): { fold: number; p: number; mat: number; pat: number } {
  const mats = rows.filter((r) => r.parent === 'maternal')
  const pats = rows.filter((r) => r.parent === 'paternal')
  const mat = rateOf(mats, pick)
  const pat = rateOf(pats, pick)
  const observed = Math.abs(mat - pat)
  const nMat = mats.length
  const pool = [...mats, ...pats]
  const next = rng(seed)
  let atLeast = 0
  for (let i = 0; i < permutations; i += 1) {
    const shuffled = [...pool]
    for (let j = shuffled.length - 1; j > 0; j -= 1) {
      const k = Math.floor(next() * (j + 1))
      ;[shuffled[j], shuffled[k]] = [shuffled[k], shuffled[j]]
    }
    const a = rateOf(shuffled.slice(0, nMat), pick)
    const b = rateOf(shuffled.slice(nMat), pick)
    if (Math.abs(a - b) >= observed - 1e-12) atLeast += 1
  }
  return {
    fold: pat > 0 ? mat / pat : NaN,
    // Add-one, so a p of exactly zero is never reported off a finite number of permutations.
    p: (atLeast + 1) / (permutations + 1),
    mat,
    pat,
  }
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
    if (!Number.isFinite(s.informative) || s.informative <= 0) {
      excluded.push({ name: s.name, why: 'no informative markers, so it has no denominator' })
      continue
    }
    usable.push(s)
  }

  const mats = usable.filter((s) => s.parent === 'maternal')
  const pats = usable.filter((s) => s.parent === 'paternal')
  const mkGroup = (parent: 'maternal' | 'paternal', rows: BalanceSample[]): GroupStat => ({
    parent,
    samples: rows.length,
    events: rows.reduce((a, s) => a + s.events.length, 0),
    eventsPerSample: rows.length
      ? rows.reduce((a, s) => a + s.events.length, 0) / rows.length : NaN,
    ratePer100k: rateOf(rows),
    medianInformative: median(rows.map((s) => s.informative)),
  })
  const groups = [mkGroup('maternal', mats), mkGroup('paternal', pats)]

  const mMed = groups[0].medianInformative
  const pMed = groups[1].medianInformative
  const skew = Number.isFinite(mMed) && Number.isFinite(pMed) && Math.min(mMed, pMed) > 0
    ? Math.max(mMed, pMed) / Math.min(mMed, pMed) : NaN
  const withinTolerance = Number.isFinite(skew) && skew <= MAX_POWER_SKEW
  const power = {
    maternalMedianInformative: mMed,
    paternalMedianInformative: pMed,
    skew,
    withinTolerance,
  }

  const methods = `Uniparental samples only, ${usable.length} of ${samples.length}: a genome with `
    + 'one parent in it attributes every change in it to that parent by construction, so no '
    + 'per-event call is needed and no per-event detection bias can enter. Counts are expressed '
    + 'per 100,000 informative markers, because power to see a change at all scales with how many '
    + 'markers could have carried evidence. The null is built by shuffling which genome is whose '
    + `over ${permutations.toLocaleString()} permutations, keeping each sample's own event and `
    + 'marker counts together, so the test asks only whether the parental label produced the '
    + `difference. Two-sided, reported against alpha ${alpha}. Detection power between the groups `
    + `is measured as a ratio of median informative markers and must sit within ${MAX_POWER_SKEW}x `
    + 'for the comparison to be reported at all.'

  const tooFew = mats.length < minPerGroup || pats.length < minPerGroup
  if (tooFew) {
    return {
      verdict: 'underpowered',
      headline: `${mats.length} maternal and ${pats.length} paternal genome`
        + `${pats.length === 1 ? '' : 's'} is too few to compare: each group needs at least `
        + `${minPerGroup}. No conclusion either way.`,
      groups, byClass: [], fold: NaN, p: NaN, alpha, permutations, power, excluded, methods,
      minAchievableP: minAchievableP(mats.length, pats.length, permutations),
    }
  }
  if (!withinTolerance) {
    return {
      verdict: 'not-comparable',
      headline: 'The two groups are not measured equally well: median informative markers differ '
        + `by ${Number.isFinite(skew) ? skew.toFixed(2) : '?'}x, past the ${MAX_POWER_SKEW}x this `
        + 'comparison allows. The better measured group would carry more detected change whatever '
        + 'the biology, so any difference here would be reporting the arrays rather than the '
        + 'genomes. No conclusion either way.',
      groups, byClass: [], fold: NaN, p: NaN, alpha, permutations, power, excluded, methods,
      minAchievableP: minAchievableP(mats.length, pats.length, permutations),
    }
  }

  const overall = permutedP(usable, undefined, permutations, seed)
  const floorP = minAchievableP(mats.length, pats.length, permutations)
  const classes = [...new Set(usable.flatMap((s) => s.events.map((e) => e.cls)))].sort()
  // Corrected for the number of classes actually tested, the same plainest correction the feature
  // enrichment uses, and one a reader can check by counting the rows.
  const classAlpha = alpha / Math.max(1, classes.length)
  const byClass: ClassStat[] = classes.map((cls) => {
    const r = permutedP(usable, (c) => c === cls, permutations, seed + 1)
    return {
      cls,
      maternal: mats.reduce((a, s) => a + s.events.filter((e) => e.cls === cls).length, 0),
      paternal: pats.reduce((a, s) => a + s.events.filter((e) => e.cls === cls).length, 0),
      maternalRate: r.mat,
      paternalRate: r.pat,
      fold: r.fold,
      p: r.p,
      significant: r.p < classAlpha,
      underpowered: floorP >= classAlpha,
    }
  }).sort((a, b) => a.p - b.p)

  const differential = overall.p < alpha
  const dir = overall.fold > 1 ? 'maternal' : 'paternal'
  const hits = byClass.filter((c) => c.significant)
  const headline = differential
    ? `Chromosomal change falls on the two parental genomes UNEQUALLY: `
      + `${overall.mat.toFixed(2)} per 100,000 informative markers on maternal genomes against `
      + `${overall.pat.toFixed(2)} on paternal, ${Number.isFinite(overall.fold)
        ? `${(overall.fold > 1 ? overall.fold : 1 / overall.fold).toFixed(2)}x more on ${dir}` : ''}`
      + `, p = ${overall.p.toFixed(4)}`
      + (hits.length ? `. Carried by ${hits.map((c) => c.cls).join(', ')}.` : '.')
    : floorP >= alpha
      ? `No difference could have been detected at this size: ${mats.length} maternal and `
        + `${pats.length} paternal genomes admit no p below ${floorP.toFixed(4)}, which is not `
        + `under alpha ${alpha}. This is not evidence that the genomes are alike. More genomes, `
        + 'not more markers, is what changes it.'
      : 'Chromosomal change falls on the two parental genomes EQUALLY, within what this many genomes '
      + `can resolve: ${overall.mat.toFixed(2)} per 100,000 informative markers on maternal against `
      + `${overall.pat.toFixed(2)} on paternal, p = ${overall.p.toFixed(4)}. `
      + 'An equal result is a result, but it is not the same as showing there is no difference: '
      + 'it says none was detectable at this size.'

  return {
    verdict: differential ? 'differential' : 'equal',
    headline,
    groups,
    byClass,
    fold: overall.fold,
    p: overall.p,
    alpha,
    permutations,
    minAchievableP: floorP,
    power,
    excluded,
    methods,
  }
}

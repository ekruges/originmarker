/**
 * runlength - the H2-versus-H3 discriminator: do paternal-absence violations RUN?
 *
 * Phase 2 of OriginMarker 2.0 (build spec 5, 10.2). Adapted from haplarithmisis module 4
 * (Zamani Esteki 2015, 10.1016/j.ajhg.2015.04.011).
 *
 * The principle is the whole method: **artefacts are independent across markers, real structural
 * events are contiguous.** A dropout at marker i says nothing about marker i+1. A deleted segment
 * removes every marker inside it. So the discriminating quantity is not whether a violation
 * occurred but whether violations run - which is why phase 1's positive-H3 states are
 * uninterpretable one at a time. At the Kothiyal floor a single Mendelian violation among 100
 * markers is simply expected.
 *
 * This needs no signal model and no fitted distribution, only genotypes and one fitted rate, so
 * it runs client-side in the same instant pass as phase 1. It is also the first thing here that
 * produces a publishable sentence: "the paternal allele is contiguously undetected across four
 * informative markers, p = 1.1e-08".
 *
 * What it does NOT do (spec 5.5), and the second of these is the likeliest way to misuse it:
 *   - it cannot separate copy-loss from copy-neutral LOH, because both remove paternal alleles
 *     contiguously. That needs LRR, which is phase 3.
 *   - it cannot see an event smaller than r_min x marker spacing, and the honest report for one
 *     of those is "below resolution", never "absent".
 *
 * Self-check:  node src/runlength.check.ts
 */

import { WINDOWS, type AB, type Marker } from './informativity.ts'

/**
 * Kothiyal 2019's per-trio Mendelian-inconsistency floor: 0.63% of variant sites across 1,314
 * nuclear families (10.1089/cmb.2018.0253). A LOWER BOUND on q, never the value itself.
 */
export const KOTHIYAL_FLOOR = 0.0063

/** Per-marker parental scores. null means the marker cannot testify on that side. */
export interface MarkerScore {
  /** 1 = the parent's allele is represented, 0 = absent, null = this marker cannot say. */
  sPat: 0 | 1 | null
  sMat: 0 | 1 | null
}

const isHom = (g: AB): boolean => g === 'AA' || g === 'BB'
const only = (g: AB): string => g[0]
const carries = (g: AB, a: string): boolean => g !== 'NC' && g.includes(a)

/**
 * Score one marker (spec 5.1).
 *
 * s_pat is defined only where the father is HOMOZYGOUS, because a heterozygous father is
 * compatible with either observed allele and so cannot testify to presence or absence. That is
 * exactly the L3-capable class of phase 1, so the scoring rule and the informativity table are
 * one object seen from two directions rather than two rules that must be kept in step.
 */
export function scoreMarker(father: AB, mother: AB | null, embryo: AB): MarkerScore {
  if (embryo === 'NC' || father === 'NC') return { sPat: null, sMat: null }
  const sPat: 0 | 1 | null = isHom(father) ? (carries(embryo, only(father)) ? 1 : 0) : null
  const known = mother !== null && mother !== 'NC'
  const sMat: 0 | 1 | null = known && isHom(mother)
    ? (carries(embryo, only(mother)) ? 1 : 0)
    : null
  return { sPat, sMat }
}

/**
 * P(longest run >= r) among n independent Bernoulli(q) trials, Erdos-Renyi approximation
 * (spec 5.2): 1 - (1 - q^r)^(n - r + 1).
 *
 * Computed as -expm1(m * log1p(-x)) rather than literally, because the direct form cancels
 * catastrophically: at q = 0.0063 and r = 10 the true value is 9.8e-23, but `1 - x` rounds to
 * exactly 1 in float64 and the subtraction returns 0. Reporting p = 0 for a real event would be
 * both wrong and unfalsifiable. Same expression, evaluated where the precision lives.
 */
export function runLengthP(r: number, n: number, q: number): number {
  if (r <= 0) return 1
  if (r > n || q <= 0) return 0
  if (q >= 1) return 1
  const m = n - r + 1
  const x = q ** r
  if (x >= 1) return 1
  return -Math.expm1(m * Math.log1p(-x))
}

/**
 * Smallest run length whose p-value clears the window-wise error budget (spec 5.3).
 *
 * Never hardcoded: it falls out of the fitted q and the number of windows tested. Report it
 * alongside the result so the threshold is auditable instead of magic.
 */
export function rMin(n: number, q: number, alpha = 0.05, windows = WINDOWS.length): number | null {
  const budget = alpha / Math.max(1, windows)
  for (let r = 1; r <= n; r++) if (runLengthP(r, n, q) <= budget) return r
  return null
}

/**
 * Compose the two sources of q (spec 5.4).
 *
 * q is NOT the dropout rate. It is the probability a marker reads paternal-absent for any reason
 * other than a real absence: genotyping error, contamination, dropout, mapping artefact.
 *
 * The max, not the sum. Kothiyal is a genome-wide trio floor that already contains error classes
 * the empirical estimator also captures, so adding them double-counts. It functions as a lower
 * bound: an empirical fit below 0.0063 is optimistic and the floor is used instead.
 */
export function qHat(qEmpirical: number | null): { q: number; source: 'empirical' | 'kothiyal_floor' } {
  if (qEmpirical === null || !Number.isFinite(qEmpirical) || qEmpirical <= KOTHIYAL_FLOOR) {
    return { q: KOTHIYAL_FLOOR, source: 'kothiyal_floor' }
  }
  return { q: qEmpirical, source: 'empirical' }
}

export type RunVerdict =
  /** A run at or above r_min: contiguous paternal absence, not independent artefact. */
  | 'significant_run'
  /** A run exists but is shorter than r_min, so it is indistinguishable from artefact here. */
  | 'no_significant_run'
  /** The window's marker spacing cannot resolve an event of the size in question. */
  | 'below_resolution'
  /** No L3-capable markers: the father is heterozygous across the window, so s_pat is undefined. */
  | 'undefined_father_heterozygous'

export interface RunResult {
  window: string
  halfWidthBp: number
  /** L3-capable markers: father homozygous, counted whether or not the embryo was called, since
   *  capability is a parental property. This is the n that enters r_min. */
  nL3: number
  /** Of those, how many the embryo was actually called at. The gap between this and nL3 is
   *  observation lost to no-calls, and it is reported rather than folded into n. */
  nScored: number
  /** Total paternal-absence markers in the window, contiguous or not. */
  zSum: number
  longestRun: number
  runP: number
  /** The mirror side, reported separately and never pooled with the paternal run. */
  nMat: number
  zSumMaternal: number
  longestRunMaternal: number
  runPMaternal: number
  rMin: number | null
  q: number
  qSource: 'empirical' | 'kothiyal_floor'
  /** Median spacing between consecutive L3-capable markers. */
  medianSpacingL3Bp: number | null
  /** r_min x median L3 spacing: the smallest event this window could possibly resolve. */
  resolutionFloorBp: number | null
  verdict: RunVerdict
  note: string
}

const median = (xs: number[]): number | null => {
  if (xs.length === 0) return null
  const s = [...xs].sort((a, b) => a - b)
  const m = s.length >> 1
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}

const longestRunOf = (flags: readonly boolean[]): number => {
  let best = 0, cur = 0
  for (const f of flags) { cur = f ? cur + 1 : 0; if (cur > best) best = cur }
  return best
}

export interface RunOptions {
  /** Fitted per-marker spurious-violation rate. null falls back to the Kothiyal floor. */
  qEmpirical?: number | null
  alpha?: number
  /** An event size the caller cares about, in bp. Supplying it turns an absent run into an
   *  explicit "below resolution" rather than letting it read as "no event". */
  eventSizeOfInterestBp?: number
}

/**
 * Run the discriminator over the three spec windows.
 *
 * Markers need not arrive sorted; they are ordered by position here, because the runs are
 * physical and an unsorted input would silently fragment them.
 */
export function analyseRuns(
  markers: readonly Marker[],
  father: Map<string, AB>,
  mother: Map<string, AB> | undefined,
  embryo: Map<string, AB>,
  variantChrom: string,
  variantPos: number,
  opts: RunOptions = {},
): RunResult[] {
  const alpha = opts.alpha ?? 0.05
  const { q, source } = qHat(opts.qEmpirical ?? null)
  const want = variantChrom.replace(/^chr/i, '').toLowerCase()

  const ordered = markers
    .filter((m) => m.chrom.replace(/^chr/i, '').toLowerCase() === want && m.pos !== variantPos)
    .slice()
    .sort((a, b) => a.pos - b.pos)

  return WINDOWS.map(({ name, halfWidthBp }) => {
    const inWin = ordered.filter((m) => Math.abs(m.pos - variantPos) <= halfWidthBp)

    const patFlags: boolean[] = []
    const matFlags: boolean[] = []
    const l3Positions: number[] = []
    let scored = 0
    for (const m of inWin) {
      const fa = father.get(m.rsid) ?? 'NC'
      const mo = mother ? (mother.get(m.rsid) ?? 'NC') : null
      const em = embryo.get(m.rsid) ?? 'NC'

      // L3-capable is a property of the PARENTS alone, so it is counted whether or not the embryo
      // was called here. That is deliberate and it is the conservative choice: n enters r_min, and
      // a larger n demands a LONGER run before significance is declared. Counting only the markers
      // that happened to produce a score would shrink n on a noisy sample and make significance
      // easier to reach exactly where the data is worst.
      const capable = fa === 'AA' || fa === 'BB'
      if (!capable) continue
      l3Positions.push(m.pos)
      const { sPat, sMat } = scoreMarker(fa, mo, em)
      if (sPat !== null) scored++
      // A no-call BREAKS a run rather than bridging it: an unobserved marker is not a
      // demonstrated absence, and letting a run span one would manufacture contiguity.
      patFlags.push(sPat === 0)
      // The mirror is scored on the SAME marker set, not on every mother-homozygous marker.
      // Both runs then share one denominator and are directly comparable, which is the point of
      // reporting them side by side: a paternal run of 4 and a maternal run of 4 mean the same
      // thing about opposite parents. Scoring the mirror over a wider set would make the two
      // numbers incommensurable while looking like a pair.
      if (mo === 'AA' || mo === 'BB') matFlags.push(sMat === 0)
    }

    const nL3 = patFlags.length
    const zSum = patFlags.filter(Boolean).length
    const run = longestRunOf(patFlags)
    const nMat = matFlags.length
    const zSumMat = matFlags.filter(Boolean).length
    const runMat = longestRunOf(matFlags)

    const gaps = l3Positions.slice(1).map((p, i) => p - l3Positions[i])
    const spacing = median(gaps)
    const threshold = nL3 > 0 ? rMin(nL3, q, alpha) : null
    const floor = threshold !== null && spacing !== null ? threshold * spacing : null

    let verdict: RunVerdict
    let note: string
    if (nL3 === 0) {
      verdict = 'undefined_father_heterozygous'
      note = 'no L3-capable markers in this window: the father is heterozygous throughout, so '
        + 'paternal presence and absence are both undefined here. This is not a negative result.'
    } else if (threshold !== null && run >= threshold) {
      verdict = 'significant_run'
      note = `${run} consecutive L3-capable markers with the paternal allele undetected, against `
        + `r_min ${threshold} at q ${q.toFixed(4)}. p = ${runLengthP(run, nL3, q).toExponential(2)}. `
        + 'Contiguity is what separates this from independent artefact. It does NOT distinguish '
        + 'copy-loss from copy-neutral LOH: both remove paternal alleles contiguously, and only '
        + 'LRR separates them.'
    } else if (
      opts.eventSizeOfInterestBp !== undefined && floor !== null
      && opts.eventSizeOfInterestBp < floor
    ) {
      verdict = 'below_resolution'
      note = `an event of ${opts.eventSizeOfInterestBp} bp cannot produce a significant run here: `
        + `r_min ${threshold} x median L3 spacing ${spacing} bp puts this window's floor at `
        + `${floor} bp. Report this as BELOW RESOLUTION, never as absent.`
    } else {
      verdict = 'no_significant_run'
      note = `longest paternal-absence run ${run} of ${nL3} L3-capable markers, below r_min `
        + `${threshold ?? 'n/a'}. Indistinguishable from independent artefact at q ${q.toFixed(4)}`
        + (floor !== null ? `. This window cannot resolve events under ${floor} bp.` : '.')
    }

    return {
      window: name,
      halfWidthBp,
      nL3,
      nScored: scored,
      zSum,
      longestRun: run,
      runP: nL3 > 0 ? runLengthP(run, nL3, q) : 1,
      nMat,
      zSumMaternal: zSumMat,
      longestRunMaternal: runMat,
      runPMaternal: nMat > 0 ? runLengthP(runMat, nMat, q) : 1,
      rMin: threshold,
      q,
      qSource: source,
      medianSpacingL3Bp: spacing,
      resolutionFloorBp: floor,
      verdict,
      note,
    }
  })
}

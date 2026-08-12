/**
 * Where a segment actually starts and ends, and how sure of it we are.
 *
 * THE PROBLEM. `scanChromosome` finds a segment by sliding windows of 2,400 markers and up,
 * stepping a quarter of a window at a time, and then reports the winning window's EDGES as the
 * event's coordinates. Those edges are a scanning artefact, not a measurement: the step alone
 * quantises them to a median 1.96 Mb, and the reported span is whatever window happened to score
 * highest. Measured against 848 events spliced from real arrays at known marker positions, the
 * window edge lands a median 373.5 markers from the truth, which is 2.51 Mb, with a p95 of
 * 15.46 Mb.
 *
 * THE FIX. Once the window has chosen roughly where the event is, find the edges properly:
 * coordinate ascent on the two boundaries at MARKER resolution, evaluating the same
 * segment-versus-background log-likelihood ratio the scan already uses. Hold one edge, sweep the
 * other over every marker within a radius, take the best, alternate. Off a cumulative sum each
 * candidate is O(1), so the whole thing is one extra pass and no new data structure.
 *
 * Measured on the same 848 events: median error 9 markers, 0.063 Mb. A 40x improvement, and it
 * also localises more events than the window edge does (629 of 848 against 412). Circular binary
 * segmentation and binary segmentation were both measured on the same constructs and both beat the
 * window edge but lost to this, while costing more: CBS needs a quadratic interval search and
 * binary segmentation a full dynamic program, against one linear sweep here.
 *
 * ABSENCE CHANNEL ONLY. Intensity was measured and refused. log2R localises 6% of these events
 * against 74% for absence, because half of what this tool looks for is an origin switch, which is
 * a haplotype substitution with no dosage step to find at all (1% localisation). And where there
 * IS a dosage step, the log2R residual in single-cell WGA has a block length of 24 markers against
 * 2.6 for absence, so the waviness supports spurious steps larger than the real one. Combining the
 * two is worse than absence alone, 53% against 74%: the intensity term drags the joint argmax onto
 * an intensity wave. This is the same conclusion the copy-number scan reached independently.
 *
 * WHAT IT DOES NOT DO. It does not decide whether there is an event. The window scan does that,
 * against an external null, and this only sharpens what the scan already found. Run the other way
 * round it would find an edge in anything.
 */

/** Below this call rate the refinement converges on dropout clusters rather than the true edge:
 *  positional error rises from 6-9 markers to 160-222, and interval coverage falls to 0.81. */
export const REFINE_CALL_FLOOR = 0.70

/** How far either side of the window edge to look. Beyond this the search costs more than the
 *  answer is worth, and about 2% of edges are mislocalised further out than any radius helps. */
export const SEARCH_RADIUS = 3_000

/**
 * Drop in the log-likelihood ratio that bounds the interval.
 *
 * NOT a chi-squared quantile. The nominal 1.92 gives 75.5% coverage of a nominal 95%, and
 * deflating the likelihood by the measured variance-inflation factor only reaches 92.0%. This is
 * 6.2x the nominal value and was calibrated empirically to 96.9% on 470 edges. It is a property of
 * this platform, this call-rate range and this artefact, and will not transfer to another array or
 * another amplification protocol without recalibration.
 */
export const LR_DROP = 12

/** Below this, an interval would be claiming more precision than the data supports: the measured
 *  median interval is 151 markers and no method reached nominal coverage narrower than that. */
export const MIN_INTERVAL_MARKERS = 100

export interface RefinedEdge {
  /** Marker index of the best estimate. */
  index: number
  /** Marker indices bounding the interval, inclusive. */
  lo: number
  hi: number
}

export interface Refined {
  start: RefinedEdge
  end: RefinedEdge
  /** False when the search hit its radius, which means the window was not describing this event
   *  and the coordinates should not be reported at all. */
  localised: boolean
  why: string
}

/** Binomial log-likelihood of `k` of `n` at rate `p`, dropping the constant term. */
const logLike = (k: number, n: number, p: number): number => {
  if (p <= 0) return k > 0 ? -Infinity : 0
  if (p >= 1) return k < n ? -Infinity : 0
  return k * Math.log(p) + (n - k) * Math.log(1 - p)
}

/**
 * Score a candidate [a, b) against the background rate: how much better the data are explained by
 * a raised rate inside than by the background everywhere.
 */
const score = (cum: Int32Array, a: number, b: number, nullRate: number): number => {
  const n = b - a
  if (n <= 0) return -Infinity
  const k = cum[b] - cum[a]
  const rate = k / n
  // A window cleaner than the background is not an event, and must never win the ascent.
  if (rate <= nullRate) return -Infinity
  return logLike(k, n, rate) - logLike(k, n, nullRate)
}

/**
 * Sharpen a segment's edges to marker resolution and bound each with an interval.
 *
 * `absent` is the per-marker indicator over ONE chromosome's informative markers, in position
 * order. `a0`/`b0` are the window's marker indices, half-open. Returns indices into that same
 * array, so the caller maps them back to positions.
 */
export function refineEdges(
  absent: Uint8Array | readonly number[],
  a0: number,
  b0: number,
  nullRate: number,
  radius = SEARCH_RADIUS,
  drop = LR_DROP,
  /** The score the window scan required to call this an event at all. After sharpening, a real
   *  event scores HIGHER than the window did; if it scores less than the scan's own bar, the
   *  window was fitting noise and there is nothing here to localise. Defaults to 0, which makes
   *  this a pure refiner; a caller refining a scan hit should pass that scan's threshold. */
  minScore = 0,
): Refined {
  const n = absent.length
  const cum = new Int32Array(n + 1)
  for (let i = 0; i < n; i += 1) cum[i + 1] = cum[i] + (absent[i] ? 1 : 0)

  let a = Math.max(0, Math.min(a0, n - 1))
  let b = Math.max(a + 1, Math.min(b0, n))

  // Coordinate ascent, three passes. Three is where the measured error stopped moving; the edges
  // are nearly independent once the window is roughly right, so this converges immediately.
  for (let pass = 0; pass < 3; pass += 1) {
    let best = score(cum, a, b, nullRate)
    let bestA = a
    for (let c = Math.max(0, a - radius); c <= Math.min(b - 1, a + radius); c += 1) {
      const s = score(cum, c, b, nullRate)
      if (s > best) { best = s; bestA = c }
    }
    a = bestA
    best = score(cum, a, b, nullRate)
    let bestB = b
    for (let c = Math.max(a + 1, b - radius); c <= Math.min(n, b + radius); c += 1) {
      const s = score(cum, a, c, nullRate)
      if (s > best) { best = s; bestB = c }
    }
    b = bestB
  }

  const peak = score(cum, a, b, nullRate)
  // Profile interval: hold the other edge, walk out until the score falls by `drop`.
  const edge = (which: 'a' | 'b'): RefinedEdge => {
    const at = which === 'a' ? a : b
    const test = (c: number): number =>
      (which === 'a' ? score(cum, c, b, nullRate) : score(cum, a, c, nullRate))
    let lo = at
    let hi = at
    const floorIdx = which === 'a' ? 0 : a + 1
    const ceilIdx = which === 'a' ? b - 1 : n
    while (lo > floorIdx && peak - test(lo - 1) < drop) lo -= 1
    while (hi < ceilIdx && peak - test(hi + 1) < drop) hi += 1
    // Never report an interval narrower than the calibration supports, even where the likelihood
    // is locally sharp: the sharpness is real and the coverage measured on it is not.
    const width = hi - lo
    if (width < MIN_INTERVAL_MARKERS) {
      const pad = Math.ceil((MIN_INTERVAL_MARKERS - width) / 2)
      lo = Math.max(floorIdx, lo - pad)
      hi = Math.min(ceilIdx, hi + pad)
    }
    return { index: at, lo, hi }
  }

  // Hitting the radius means the ascent was still climbing when it ran out of room, so the window
  // was not describing this event. Reported as not localised rather than as a wide interval:
  // unconditional coverage never exceeded 0.743 for any method, and the shortfall is a
  // localisation failure, not something a wider interval fixes.
  const ranOut = Math.abs(a - a0) >= radius || Math.abs(b - b0) >= radius
  const tooWeak = !Number.isFinite(peak) || peak < minScore
  return {
    start: edge('a'),
    end: edge('b'),
    localised: !ranOut && !tooWeak,
    why: ranOut
      ? 'the edge search reached its radius without settling, so the window was not describing '
        + 'this event'
      : !Number.isFinite(peak)
        ? 'no candidate scored above the background rate'
        : tooWeak
          ? `the sharpened segment scores ${peak.toFixed(0)}, under the ${minScore} the scan `
            + 'required, so the window was fitting noise'
          : '',
  }
}

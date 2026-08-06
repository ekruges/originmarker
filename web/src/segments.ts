/**
 * segments - where along a chromosome is the parental genome missing, rather than whether.
 *
 * The per-chromosome verdict in `parentage` answers "is this whole chromosome absent". A real
 * event is often smaller than that, and a chromosome that is half lost reads as neither present
 * nor absent: the rate lands between the two references and the chromosome comes back unclear,
 * with nothing said about the half that IS missing.
 *
 * Three things had to be got right here, each of them a measurement rather than a preference.
 *
 * THE NULL IS EXTERNAL. Scoring a window against the sample's OWN rate misses the events that
 * matter most, because a large event inflates the rate it is being tested against. Measured on a
 * degraded sample carrying five true whole-chromosome losses: against the self null those five
 * score far below the threshold and come back negative, while against a robust median of the
 * other chromosomes they score orders of magnitude above it. So the null is the median
 * per-chromosome rate over the OTHER chromosomes, which one event cannot move.
 *
 * INDEPENDENCE IS NOT ASSUMED. A Bonferroni-corrected exact binomial produces 4 to 30 false
 * segments per known-clean genome, and an Erdos-Renyi run-length scan fails the same way, because
 * absence artefact on amplified single cells is spatially clustered: runs of two occur 6 to 10
 * times more often than independence predicts. Permuting the same indicator, which preserves the
 * rate and destroys the clustering, drops the false segments to zero. That is the diagnosis. So
 * the threshold here is EMPIRICAL, taken from genomes known to carry nothing, and never from a
 * closed-form tail.
 *
 * A FLOOR IN MARKERS IS NOT A RESOLUTION. Marker spacing on this array runs from 1 bp to 21 kb,
 * so a fixed 200-marker window spans anywhere from about 1 kb to 32 Mb. Every segment reported
 * here carries both its marker count and its physical span, because only one of those is a
 * resolution and it is not the one the statistic is computed in.
 *
 * Self-check:  node src/segments.check.ts
 */
import { isAutosome } from './parentage.ts'

/**
 * Smallest segment that can be called, in CALLED informative markers.
 *
 * Titrated on real genomes rather than simulated ones: a block of a MATERNAL pronucleus spliced
 * into a clean paternal pronucleus of the same series, so the segment is a genuine alternative
 * genome carrying real amplification artefact, arrayed on the same platform in the same lab.
 * Twelve constructions per size, four backgrounds by three donors:
 *
 *     spliced markers   score range   detected at SEGMENT_LRT
 *          1,200          152 -  197        no
 *          2,400          431 -  556        yes
 *          4,800          968 - 1319        yes
 *          9,600         1807 - 2450        yes
 *         19,200         3489 - 4582        yes
 *
 * Against a null maximum of 139 across five genomes carrying no event. At 1,200 markers the
 * weakest construction scores 1.09x that maximum, which is not a detection; at 2,400 the weakest
 * scores 3.1x it. The floor sits where separation becomes real rather than where it becomes
 * nonzero.
 *
 * A marker count is not a resolution. At this floor the reported span was 11 Mb on chromosome 2;
 * on a denser stretch the same 2,400 markers cover far less. Both numbers travel with every
 * segment for that reason.
 */
export const MIN_SEGMENT_MARKERS = 2_400

/**
 * Absence rate below which a chromosome contributes nothing to the null.
 *
 * A sample with a genuinely clean genome has per-chromosome rates near zero, and a median of
 * those would make any window look catastrophic. The floor is the residual absence genotyping
 * error alone produces, measured at 0.03% and 0.05% on clean data.
 */
export const NULL_FLOOR = 0.002

/**
 * Log-likelihood ratio a segment must clear.
 *
 * EMPIRICAL, and FITTED rather than validated. Measured on the five paternal pronuclei of the
 * bundled series, each a meiotic product of the sperm donor and therefore present on every
 * autosome, so the truth is zero segments: across 110 chromosome scans the maximum score reached
 * was 139, and per-chromosome maxima ran 58 to 139. The weakest real event at the marker floor
 * scores 431. This sits between them, 1.8x above the null maximum and 1.7x below the weakest
 * detection.
 *
 * What that does and does not buy: it bounds the false rate at zero per genome ON THOSE FIVE
 * GENOMES, which is not the same as bounding it in general. The threshold is tuned on the same
 * null it is evaluated against, so the honest reading is that it awaits an out-of-sample clean
 * cohort. A segment scoring just over it is a weak finding and the score is reported beside every
 * one so that reads as what it is.
 */
export const SEGMENT_LRT = 250

export interface Segment {
  chrom: string
  /** Called informative markers inside the segment. The statistic is computed in these. */
  markers: number
  absent: number
  rate: number
  startBp: number
  endBp: number
  /** Physical span. A marker count is not a resolution on an array with 1 bp to 21 kb spacing. */
  spanBp: number
  /** Log-likelihood ratio against the external null. Reported so a marginal call reads marginal. */
  score: number
  /** The rate this was scored against, and where it came from. */
  nullRate: number
}

export interface MarkerAbsence {
  chrom: string
  pos: number
  /** The parent's allele is not present here. Only markers where the parent is homozygous and
   *  the sample is called can say this; the caller supplies those and no others. */
  absent: boolean
}

/** log P(k absent of n) under a binomial with rate p, up to terms that cancel in the ratio. */
const logLike = (k: number, n: number, p: number): number => {
  if (p <= 0) return k > 0 ? -Infinity : 0
  if (p >= 1) return k < n ? -Infinity : 0
  return k * Math.log(p) + (n - k) * Math.log1p(-p)
}

/**
 * The robust null: the median per-chromosome absence rate over every chromosome EXCEPT the one
 * being scanned, floored.
 *
 * Excluding the chromosome under test is what stops a whole-chromosome event from raising the bar
 * it has to clear. Taking the median rather than the mean is what stops several events from doing
 * it collectively.
 */
export function externalNull(byChrom: Map<string, [number, number]>, exclude: string): number {
  const rates: number[] = []
  for (const [c, [n, a]] of byChrom) {
    if (c === exclude || !isAutosome(c) || n < 200) continue
    rates.push(a / n)
  }
  if (!rates.length) return NULL_FLOOR
  rates.sort((x, y) => x - y)
  const m = rates.length >> 1
  const med = rates.length % 2 ? rates[m] : (rates[m - 1] + rates[m]) / 2
  return Math.max(med, NULL_FLOOR)
}

/**
 * Scan one chromosome for segments where the parental genome is missing.
 *
 * Multiscale: windows are tried at several sizes because an event has no reason to match one, and
 * the best-scoring window at any size is the one that describes it. See the peak-picking below
 * for why that is not the same as merging everything that overlaps.
 */
export function scanChromosome(
  markers: readonly MarkerAbsence[],
  nullRate: number,
  minMarkers = MIN_SEGMENT_MARKERS,
  threshold = SEGMENT_LRT,
): Segment[] {
  const ms = [...markers].sort((a, b) => a.pos - b.pos)
  const n = ms.length
  if (n < minMarkers) return []

  // Prefix sums, so a window's absence count is O(1) and the scan stays linear per scale.
  const cum = new Int32Array(n + 1)
  for (let i = 0; i < n; i += 1) cum[i + 1] = cum[i] + (ms[i].absent ? 1 : 0)

  const hits: Segment[] = []
  for (let w = minMarkers; w <= n; w *= 2) {
    const step = Math.max(1, Math.floor(w / 4))
    for (let i = 0; i + w <= n; i += step) {
      const k = cum[i + w] - cum[i]
      const rate = k / w
      // Only a RAISED rate is an event. A window cleaner than the genome is not a finding.
      if (rate <= nullRate) continue
      const score = logLike(k, w, rate) - logLike(k, w, nullRate)
      if (score < threshold) continue
      hits.push({
        chrom: ms[i].chrom,
        markers: w,
        absent: k,
        rate,
        startBp: ms[i].pos,
        endBp: ms[i + w - 1].pos,
        spanBp: ms[i + w - 1].pos - ms[i].pos,
        score,
        nullRate,
      })
    }
  }
  if (!hits.length) return []

  // One event clearing the threshold at several scales is one event, and the question is which
  // window DESCRIBES it. Taking the union of everything that overlaps is wrong: a whole-chromosome
  // window carrying a diluted version of the same event overlaps the tight one, and unioning them
  // reported a 12 Mb loss as spanning 231 Mb. So this is peak-picking, not merging. The
  // highest-scoring window wins, everything overlapping it is discarded as another view of the
  // same event, and the scan repeats on what is left.
  hits.sort((a, b) => b.score - a.score)
  const picked: Segment[] = []
  for (const h of hits) {
    if (picked.some((p) => h.startBp <= p.endBp && h.endBp >= p.startBp)) continue
    picked.push(h)
  }
  return picked.sort((a, b) => a.startBp - b.startBp)
}

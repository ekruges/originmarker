/**
 * A chromosome the whole-chromosome test calls lost, read again in pieces.
 *
 * The whole test reads one number per chromosome, its median log R. A deep loss over a third of a
 * chromosome can move that median over both whole-chromosome gates on its own, and a whole call
 * then hides the piece and the edge where the chromosome was cut. So on a chromosome called lost
 * whole the copy-number scan runs as on any other, and the part of the chromosome outside its loss
 * pieces is put to the same whole-chromosome test, with the same null, floor and call-rate bar:
 *
 *  - the rest reads at the cell's own level (its shift under the floor the array earned, its call
 *    rate not collapsed): the pieces account for the shift, and they are reported instead of the
 *    whole call. The z is not asked here: it says a change is bigger than the array's noise, and a
 *    noisy rest a full copy down fails it without reading at the cell's level;
 *  - otherwise, or where the rest is too short to read: the whole call stands, and a piece is
 *    reported beside it only where it sits a further SEGMENT_LRR_SHIFT below the rest, a second
 *    copy lost inside a chromosome already down one. Any other piece there is the whole loss
 *    described twice;
 *  - pieces covering RESTATED_SHARE of the chromosome's span are the whole loss: the whole call
 *    stands and no piece is reported. What is left beside them is a sliver, not the chromosome.
 *
 * No constant here is new: the floor and the call-rate median are the ones parentage.ts judged the
 * whole chromosome by, passed through as `ParentageResult.copyNull`; RESTATED_SHARE is the share at
 * which the copy scan already takes a piece to restate an allele-named loss.
 */
import { calibratedZ, intensityDetects } from './intensityNull.ts'
import { CALL_COLLAPSE } from './parentage.ts'
import { scanCopyNumber, SEGMENT_LRR_SHIFT, type Segment } from './segments.ts'

/** Share of a chromosome's span past which its loss pieces are the whole loss. Same as scoreSample's restatement filter. */
export const RESTATED_SHARE = 0.8

export interface CopyNull { centre: number; scale: number; units: number; floor: number; callMedian: number }

/** The fewest markers a chromosome, or the rest of one, needs before the whole test reads it. */
export const REST_MIN_MARKERS = 200

type Marker = { chrom: string; pos: number; called: boolean; log2R: number | null }

const medianOf = (xs: number[]): number => {
  if (!xs.length) return NaN
  const q = [...xs].sort((x, y) => x - y)
  return q[q.length >> 1]
}

/** Whether a stretch of markers is a change by the whole-chromosome test of parentage.ts. */
export function wholeTestChanges(ms: readonly Marker[], n: CopyNull): boolean {
  const med = medianOf(ms.map((m) => m.log2R).filter((x): x is number => x !== null))
  const byIntensity = intensityDetects(calibratedZ(med, n), true) && Number.isFinite(med)
    && Math.abs(med - n.centre) >= n.floor
  return byIntensity || callFraction(ms, n) < CALL_COLLAPSE
}

/** Whether a stretch of markers reads at the cell's own level: shift under the floor, calls not collapsed. */
export function atCellLevel(ms: readonly Marker[], n: CopyNull): boolean {
  const med = medianOf(ms.map((m) => m.log2R).filter((x): x is number => x !== null))
  return Number.isFinite(med) && Math.abs(med - n.centre) < n.floor && !(callFraction(ms, n) < CALL_COLLAPSE)
}

const callFraction = (ms: readonly Marker[], n: CopyNull): number => (
  n.callMedian > 0 && ms.length ? ms.filter((m) => m.called).length / ms.length / n.callMedian : NaN)

export interface PartLost {
  /** The whole-chromosome loss call stands. */
  whole: boolean
  /** The loss pieces to report on this chromosome. */
  pieces: Segment[]
  /** Median log R of the rest, for the record. NaN when the rest is too short. */
  restLrr: number
}

/**
 * Read a chromosome the whole test called lost. `nullNoCall` and `genomeLrr` are what the copy scan
 * uses on every other chromosome.
 */
export function readLostChromosome(
  ms: readonly Marker[], nullNoCall: number, genomeLrr: number, n: CopyNull,
): PartLost {
  const pieces = scanCopyNumber(ms, nullNoCall, genomeLrr).filter((s) => s.kind === 'copy-loss')
  if (!pieces.length) return { whole: true, pieces: [], restLrr: NaN }
  const inside = (m: Marker) => pieces.some((s) => m.pos >= s.startBp && m.pos <= s.endBp)
  const rest = ms.filter((m) => !inside(m))
  const restLrr = medianOf(rest.map((m) => m.log2R).filter((x): x is number => x !== null))
  let lo = Infinity; let hi = -Infinity
  for (const m of ms) { if (m.pos < lo) lo = m.pos; if (m.pos > hi) hi = m.pos }
  const covered = pieces.reduce((t, s) => t + (s.endBp - s.startBp), 0)
  if (hi > lo && covered >= RESTATED_SHARE * (hi - lo)) return { whole: true, pieces: [], restLrr }
  if (rest.length >= REST_MIN_MARKERS && atCellLevel(rest, n)) {
    return { whole: false, pieces, restLrr }
  }
  const deeper = pieces.filter((s) => {
    const own = medianOf(ms.filter((m) => m.pos >= s.startBp && m.pos <= s.endBp)
      .map((m) => m.log2R).filter((x): x is number => x !== null))
    return rest.length >= REST_MIN_MARKERS && own <= restLrr - SEGMENT_LRR_SHIFT
  })
  return { whole: true, pieces: deeper, restLrr }
}

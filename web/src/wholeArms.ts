/**
 * A whole-chromosome call, read arm by arm.
 *
 * The whole test (parentage.ts, classify) reads one median log R and one call fraction per
 * chromosome. An arm holding most of a chromosome's markers moves that median by close to its own
 * shift, and a deep loss on one arm drops the whole chromosome's call fraction under CALL_COLLAPSE,
 * so a change on one arm passes for the whole chromosome. A whole call is a claim about each arm, so
 * each usable arm is read on its own against the centre and call-rate median the whole test used:
 *
 *  - an arm carries a loss where its shift is at or under -need/2, or where its call fraction is
 *    under CALL_COLLAPSE while its shift is under +need/2. A loss lowers intensity, so a raised arm
 *    does not carry a loss on its call fraction alone;
 *  - an arm carries a gain where its shift is at or over +need/2;
 *  - the whole call stands where both usable arms carry it, and where fewer than two arms are usable:
 *    on an acrocentric the long arm is all of the chromosome the chip reads;
 *  - otherwise it is withdrawn, and where exactly one arm carries it that arm is named, for the
 *    caller to print where the piece scans print nothing of that kind on it (armSegment).
 *
 * `need` is the floor the whole call itself had to clear, set from the sign of the chromosome's
 * shift as classify sets it, so need/2 is no new bar. An arm is usable at REST_MIN_MARKERS markers
 * carrying a log R. The arms are GRCh37's: on any other assembly the caller does not ask.
 */
import { CALL_COLLAPSE, COPY_SHIFT_GAIN_FLOOR } from './parentage.ts'
import { callFraction, medianOf, REST_MIN_MARKERS, type CopyNull } from './partLost.ts'
import type { Segment } from './segments.ts'

/** GRCh37 p/q boundary per autosome: the end of the last p-arm centromere band, UCSC hg19 cytoBand
 *  (audit/ucsc-hg19-cytoband.tsv). A marker before it is on p, at or after it on q. */
export const P_ARM_END_GRCH37: Readonly<Record<string, number>> = {
  1: 125_000_000, 2: 93_300_000, 3: 91_000_000, 4: 50_400_000, 5: 48_400_000, 6: 61_000_000,
  7: 59_900_000, 8: 45_600_000, 9: 49_000_000, 10: 40_200_000, 11: 53_700_000, 12: 35_800_000,
  13: 17_900_000, 14: 17_600_000, 15: 19_000_000, 16: 36_600_000, 17: 24_000_000, 18: 17_200_000,
  19: 26_500_000, 20: 27_500_000, 21: 13_200_000, 22: 14_700_000,
}

type Marker = { chrom: string; pos: number; called: boolean; log2R: number | null }

export interface ArmReading {
  /** The whole call stands. */
  whole: boolean
  /** Withdrawn, and this arm alone carries the change. */
  arm?: 'p' | 'q'
  /** Each arm's shift from the cell's centre and its call fraction; null where the arm is not usable. */
  p: { shift: number; callFraction: number } | null
  q: { shift: number; callFraction: number } | null
}

const readArm = (ms: readonly Marker[], n: CopyNull) => {
  const lrr = ms.map((m) => m.log2R).filter((x): x is number => x !== null)
  return lrr.length < REST_MIN_MARKERS ? null
    : { shift: medianOf(lrr) - n.centre, callFraction: callFraction(ms, n) }
}

/** Read a chromosome the whole test called `kind`, whose whole-chromosome shift is `shift`. */
export function readArms(
  ms: readonly Marker[], kind: 'loss' | 'gain', shift: number, n: CopyNull,
): ArmReading {
  const split = P_ARM_END_GRCH37[ms[0]?.chrom ?? '']
  if (split === undefined) return { whole: true, p: null, q: null }
  const p = readArm(ms.filter((m) => m.pos < split), n)
  const q = readArm(ms.filter((m) => m.pos >= split), n)
  if (!p || !q) return { whole: true, p, q }
  const need = shift > 0 ? Math.min(n.floor, COPY_SHIFT_GAIN_FLOOR) : n.floor
  const carries = (a: { shift: number; callFraction: number }): boolean => (kind === 'gain'
    ? a.shift >= need / 2
    : a.shift <= -need / 2 || (a.callFraction < CALL_COLLAPSE && a.shift < need / 2))
  const cp = carries(p); const cq = carries(q)
  if (cp && cq) return { whole: true, p, q }
  return { whole: false, arm: cp ? 'p' : cq ? 'q' : undefined, p, q }
}

/**
 * The arm a withdrawn whole call rests on, as a copy-number segment: the p arm from its first marker
 * to the boundary, the q arm from the boundary to its last marker. Its edges are not placed: the
 * change may start anywhere on the arm. No likelihood is computed for it, so `score` is 0 and the
 * rates are NaN; `arm` marks it.
 */
export function armSegment(ms: readonly Marker[], arm: 'p' | 'q', kind: 'copy-loss' | 'copy-gain'): Segment {
  const split = P_ARM_END_GRCH37[ms[0].chrom]
  const on = ms.filter((m) => (arm === 'p' ? m.pos < split : m.pos >= split))
  let lo = Infinity; let hi = -Infinity
  for (const m of on) { if (m.pos < lo) lo = m.pos; if (m.pos > hi) hi = m.pos }
  const startBp = arm === 'p' ? lo : split
  const endBp = arm === 'p' ? split : hi
  return {
    kind, chrom: ms[0].chrom, arm, markers: on.filter((m) => m.log2R !== null).length, absent: 0,
    rate: NaN, nullRate: NaN, score: 0, startBp, endBp, spanBp: endBp - startBp,
  }
}

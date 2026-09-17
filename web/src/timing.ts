import { callUniformity, unitsCarrying, type EventLocation } from './abnormalities.ts'
import { segmentCoords, type Segment } from './segments.ts'

/**
 * When a change happened: carried by the gamete, or arisen after fertilisation.
 *
 * NOT SEPARABLE ON ONE UNIT at any material quality, which is why this lives apart from the
 * per-sample pipeline: it reads two or more independently sampled units of the SAME embryo. Two
 * patterns separate the two cases, and one of them is decisive:
 *
 *   reciprocal     one unit gains what another loses over the same interval. Two daughter cells
 *                  split one chromosome between them, so the error happened at a division AFTER
 *                  fertilisation. A change carried by the gamete is in every cell and cannot read
 *                  this way.
 *   uniform        present in every unit. Measured on this corpus: meiotic 64 of 64 uniform,
 *                  post-zygotic 6 of 7 non-uniform.
 *
 * Whole chromosomes are included here as well as segments. A reciprocal whole-chromosome pair is
 * the commonest mitotic signature in cleavage-stage embryos, and reading only segments missed it.
 */

export type Mechanism = 'meiotic' | 'post-zygotic' | 'reciprocal' | 'unresolved'

/** An event with its direction, where the tool called one. Findings carry none. */
export interface DirectedEvent extends EventLocation {
  direction?: 'gain' | 'loss'
}

/** A scored sample, in the shape both the page and the harnesses hold it. */
interface ScoredLike {
  chroms?: readonly { chrom: string; aneuploidy?: 'loss' | 'gain' }[]
  segments?: readonly Segment[]
  findings?: readonly { chrom: string; startBp: number; endBp: number }[]
}

/**
 * Every event of one scored sample, as intervals with a direction.
 *
 * A whole chromosome is written as the whole interval rather than its called extent: the call is
 * about the chromosome, and an interval intersection is what the timing test does.
 */
export function eventsOf(r: ScoredLike): DirectedEvent[] {
  const out: DirectedEvent[] = []
  for (const c of r.chroms ?? []) {
    if (c.aneuploidy) {
      out.push({ chrom: c.chrom, startBp: 0, endBp: Number.MAX_SAFE_INTEGER, direction: c.aneuploidy })
    }
  }
  for (const sg of r.segments ?? []) {
    const co = segmentCoords(sg)
    const kind = String(sg.kind)
    out.push({
      chrom: sg.chrom,
      startBp: co.start,
      endBp: co.end,
      direction: /gain/i.test(kind) ? 'gain' : /loss|absence/i.test(kind) ? 'loss' : undefined,
    })
  }
  for (const f of r.findings ?? []) {
    if (f.chrom !== 'genome') out.push({ chrom: f.chrom, startBp: f.startBp, endBp: f.endBp })
  }
  return out
}

const overlapping = (a: DirectedEvent, b: DirectedEvent): boolean =>
  a.chrom === b.chrom && a.startBp <= b.endBp && b.startBp <= a.endBp

/**
 * When one event of one unit happened, against the other units of its embryo.
 *
 * `perUnit` holds every unit of the embryo INCLUDING the one the event came from, which is what
 * makes the uniform case a count of all of them.
 */
export function callMechanism(
  event: DirectedEvent, perUnit: readonly (readonly DirectedEvent[])[],
): { mechanism: Mechanism; why: string } {
  if (perUnit.length < 2) return callUniformity(1, 1)

  if (event.direction) {
    const opposite = event.direction === 'gain' ? 'loss' : 'gain'
    const mirrored = perUnit.filter(
      (evs) => evs.some((e) => e.direction === opposite && overlapping(e, event)),
    ).length
    if (mirrored > 0) {
      return {
        mechanism: 'reciprocal',
        why: `${mirrored} other unit${mirrored === 1 ? '' : 's'} of this embryo `
          + `${mirrored === 1 ? 'carries' : 'carry'} the opposite change over the same interval: a `
          + `${event.direction} here against a ${opposite} there. Two daughter cells splitting one `
          + 'chromosome between them is a division error after fertilisation. A change carried by '
          + 'the gamete is in every cell and cannot be reciprocal',
      }
    }
  }
  return callUniformity(unitsCarrying(event, perUnit), perUnit.length)
}

/**
 * The timing of every event of every unit of one embryo.
 *
 * Returned in the order the units were given, so a caller can attach each list to its own sample
 * without matching anything back up.
 */
export function timeGroup(units: readonly ScoredLike[]): (DirectedEvent & {
  mechanism: Mechanism; why: string
})[][] {
  const perUnit = units.map(eventsOf)
  return perUnit.map((evs) => evs.map((e) => ({ ...e, ...callMechanism(e, perUnit) })))
}

/**
 * When a change happened: carried by the gamete, or arisen after fertilisation.
 *
 * NOT SEPARABLE ON ONE UNIT, which is why this lives apart from the per-sample pipeline: it reads
 * two or more units of the SAME embryo (abnormalities.groupUnits) and reads results, changing none.
 * Plan item 3.10, rule B1:
 *
 *   from the gamete        every unit carries a change of the same direction on that chromosome,
 *                          overlapping, and where a parent is named on both, the same parent.
 *                          Three units or more.
 *   consistent with the    the same with two units. Two cells can share a mitotic error made at
 *   gamete                 the first division, so two units never read "from the gamete".
 *   after fertilisation,   another unit gains what this one loses (or the reverse) over an
 *   reciprocal             overlapping interval of the same chromosome.
 *   after fertilisation,   some units carry it and at least one other reads that chromosome
 *   not in every cell      measured normal: measured, no row of any kind on it, cell reportable.
 *   none                   anything else, and any single unit. No origin word is printed.
 */

export type Mechanism = 'meiotic' | 'consistent-with-gamete' | 'mitotic-reciprocal' | 'mitotic-not-all' | 'none'

/** The words the page prints after a change, or nothing. */
export const ORIGIN_WORD: Record<Mechanism, string> = {
  meiotic: 'from the gamete',
  'consistent-with-gamete': 'consistent with the gamete (two cells read)',
  'mitotic-reciprocal': 'after fertilisation',
  'mitotic-not-all': 'after fertilisation',
  none: '',
}

/** One change of one unit, as the result card lists it (pageSays.defectsForResult). */
export interface UnitChange {
  chrom: string
  startBp: number
  endBp: number
  kind: string
  origin: string
  /** A whole-chromosome row: it spans the chromosome for overlap. */
  whole: boolean
}

/** What one unit says, in the shape both the page and the record harness can build. */
export interface UnitReading {
  changes: readonly UnitChange[]
  /** Chromosomes the unit measured (`chroms`); an unmeasured chromosome never reads normal. */
  measured: readonly string[]
  /** Stage not failed and no sample-level integrity refusal. */
  reportable: boolean
}

const GAIN = new Set(['copy-gain', 'trisomy', 'segmental-duplication'])
const LOSS = new Set(['copy-loss', 'parental-absence', 'monosomy', 'segmental-deletion'])

/** Gain, loss, or none for a kind that has no direction (copy-neutral findings, the imbalance). */
export const directionOf = (kind: string): 'gain' | 'loss' | undefined =>
  (GAIN.has(kind) ? 'gain' : LOSS.has(kind) ? 'loss' : undefined)

const overlapping = (a: UnitChange, b: UnitChange): boolean =>
  a.chrom === b.chrom && (a.whole || b.whole || (a.startBp <= b.endBp && b.startBp <= a.endBp))

const named = (o: string): boolean => o === 'maternal' || o === 'paternal'

/** A unit's reading of a scored sample, from its card rows. */
export function readingOf(
  defects: readonly { chrom: string; startBp: number; endBp: number; kind: string; origin: string; locus?: string }[],
  r: {
    chroms?: readonly { chrom: string }[]
    stage?: { stage?: string }
    integrity?: { level?: string }
  },
): UnitReading {
  return {
    changes: defects.filter((d) => d.chrom !== 'genome').map((d) => ({
      chrom: d.chrom, startBp: d.startBp, endBp: d.endBp, kind: d.kind, origin: d.origin,
      whole: d.locus === `chr${d.chrom}`,
    })),
    measured: (r.chroms ?? []).map((c) => c.chrom),
    reportable: r.stage?.stage !== 'failed' && (r.integrity?.level ?? 'intact') === 'intact',
  }
}

const unitsWord = (n: number): string => `${n} unit${n === 1 ? '' : 's'}`

/**
 * The origin of one change of unit `self`, against every unit of its embryo (`self` included).
 * The reason states how many units carry it and how many read the chromosome normal.
 */
export function callMechanism(
  self: number, change: UnitChange, units: readonly UnitReading[],
): { mechanism: Mechanism; why: string } {
  const dir = directionOf(change.kind)
  const n = units.length
  if (n < 2) return { mechanism: 'none', why: 'one unit of this embryo was read; one cell cannot say when a change arose' }
  if (!dir) return { mechanism: 'none', why: 'a change with no direction carries no origin' }

  const opposite = dir === 'gain' ? 'loss' : 'gain'
  const mirrored = units.filter((u, i) => i !== self
    && u.changes.some((x) => directionOf(x.kind) === opposite && overlapping(x, change))).length
  const carrying = units.filter((u) => u.changes.some((x) => directionOf(x.kind) === dir
    && overlapping(x, change) && (!named(x.origin) || !named(change.origin) || x.origin === change.origin))).length
  const normal = units.filter((u, i) => i !== self && u.reportable && u.measured.includes(change.chrom)
    && !u.changes.some((x) => x.chrom === change.chrom)).length
  const read = `${carrying} of ${n} units carry this ${dir}, ${normal} read${normal === 1 ? 's' : ''} `
    + `chromosome ${change.chrom} normal`

  if (mirrored > 0) {
    return {
      mechanism: 'mitotic-reciprocal',
      why: `${read}, and ${unitsWord(mirrored)} of this embryo carr${mirrored === 1 ? 'ies' : 'y'} the `
        + `opposite change over the same interval. One daughter cell gaining what another loses is a `
        + 'division error after fertilisation',
    }
  }
  if (carrying === n) {
    return n >= 3
      ? { mechanism: 'meiotic', why: `${read}: every unit of this embryo carries it, as a change carried by the gamete is in every cell` }
      : { mechanism: 'consistent-with-gamete', why: `${read}: both cells read carry it. Two cells can share an error made at the first division, so two cells do not settle it` }
  }
  if (normal > 0) {
    return {
      mechanism: 'mitotic-not-all',
      why: `${read}: a cell of the same embryo without it means it arose after fertilisation, in a lineage`,
    }
  }
  return {
    mechanism: 'none',
    why: `${read}: not every unit carries it and none reads the chromosome measured normal, so no origin is called`,
  }
}

/** One change of a unit with its origin, located as the card row it came from. */
export interface TimedChange extends UnitChange { mechanism: Mechanism; why: string }

/** Every change of every unit of one embryo, in the order the units were given. */
export function timeGroup(units: readonly UnitReading[]): TimedChange[][] {
  return units.map((u, i) => u.changes.map((c) => ({ ...c, ...callMechanism(i, c, units) })))
}

/**
 * When a change happened: carried by the gamete, or arisen after fertilisation.
 *
 * NOT SEPARABLE ON ONE UNIT, which is why this lives apart from the per-sample pipeline: it reads
 * two or more units of the SAME embryo (abnormalities.groupUnits) and reads results, changing none.
 * Plan item 3.10, rule B1:
 *
 *   from the gamete        every unit carries a change of the same direction over the stretch,
 *                          and the parents named on them are not two different parents.
 *                          Three units or more.
 *   consistent with the    the same with two units. Two cells can share a mitotic error made at
 *   gamete                 the first division, so two units never read "from the gamete".
 *   after fertilisation,   one unit gains what another loses over the same stretch.
 *   reciprocal
 *   after fertilisation,   some units carry it and at least one other reads the stretch measured
 *   not in every cell      normal: chromosome measured, no row of any kind over it, cell reportable.
 *   none                   anything else, and any single unit. No origin word is printed.
 *
 * READ PER STRETCH, NOT PER CALL (decision D5): the same stretch prints the same word in every cell
 * of the embryo, whatever shape each cell's call takes. See stretchesOf.
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
  /** Stage not failed and no sample-level integrity refusal: may read a stretch normal. */
  reportable: boolean
  /**
   * The sample is reportable at all (pageSays.reportOf, decided by the caller). A unit not shown
   * takes no part in its embryo's reading and gets no origin word: a cell whose normal level cannot
   * be told cannot show that it shares a change or lacks one. Absent means shown.
   */
  shown?: boolean
}

const GAIN = new Set(['copy-gain', 'trisomy', 'segmental-duplication'])
const LOSS = new Set(['copy-loss', 'parental-absence', 'monosomy', 'segmental-deletion'])

/** Gain, loss, or none for a kind that has no direction (copy-neutral findings, the imbalance). */
export const directionOf = (kind: string): 'gain' | 'loss' | undefined =>
  (GAIN.has(kind) ? 'gain' : LOSS.has(kind) ? 'loss' : undefined)

const named = (o: string): boolean => o === 'maternal' || o === 'paternal'

/** A unit's reading of a scored sample, from its card rows. */
export function readingOf(
  defects: readonly { chrom: string; startBp: number; endBp: number; kind: string; origin: string; locus?: string }[],
  r: {
    chroms?: readonly { chrom: string }[]
    stage?: { stage?: string }
    integrity?: { level?: string }
  },
  shown = true,
): UnitReading {
  return {
    changes: defects.filter((d) => d.chrom !== 'genome').map((d) => ({
      chrom: d.chrom, startBp: d.startBp, endBp: d.endBp, kind: d.kind, origin: d.origin,
      whole: d.locus === `chr${d.chrom}`,
    })),
    measured: (r.chroms ?? []).map((c) => c.chrom),
    reportable: r.stage?.stage !== 'failed' && (r.integrity?.level ?? 'intact') === 'intact',
    shown,
  }
}

const unitsWord = (n: number): string => `${n} unit${n === 1 ? '' : 's'}`

/**
 * Breakpoints closer than this are one breakpoint. A refined edge carries an interval of a median
 * 151 markers (segments.segmentCoords) and an edge not localised is a window edge with a larger
 * error, while the placement bar takes an edge within 1 Mb as placed. Two edges nearer than that are
 * not shown to be two places, so the stretch between them gets no reading of its own.
 */
export const SAME_BREAK_BP = 1e6

/** The end of a stretch that runs to the end of the chromosome. */
export const TO_END = Number.MAX_SAFE_INTEGER

/** One stretch of a change, with the origin the embryo reads there. */
export interface Stretch { fromBp: number; toBp: number; mechanism: Mechanism; why: string }

/** Cuts of one chromosome: every edge of every row of every unit, near-equal ones joined. */
function cutsOf(chrom: string, units: readonly UnitReading[]): { snap: (x: number) => number; at: number[] } {
  const edges = [...new Set([0, ...units.flatMap((u) => u.changes)
    .filter((c) => c.chrom === chrom && !c.whole).flatMap((c) => [c.startBp, c.endBp])])].sort((a, b) => a - b)
  const rep = new Map<number, number>()
  const at: number[] = []
  // ponytail: single linkage, so a run of cuts each under 1 Mb from the last joins however long it
  // grows; cluster by the edges' own intervals if a run ever spans several Mb.
  edges.forEach((x, i) => {
    if (i === 0 || x - edges[i - 1] >= SAME_BREAK_BP) at.push(x)
    rep.set(x, at[at.length - 1])
  })
  return { snap: (x) => rep.get(x) ?? x, at }
}

type Covers = (c: UnitChange) => boolean

/** The embryo's reading of one stretch for one direction. The same for every unit that carries it. */
function readStretch(
  chrom: string, dir: 'gain' | 'loss', covers: Covers, units: readonly UnitReading[],
): { mechanism: Mechanism; why: string } {
  const n = units.length
  const on = (u: UnitReading) => u.changes.filter((x) => x.chrom === chrom && covers(x))
  const opposite = dir === 'gain' ? 'loss' : 'gain'
  const carry = units.map((u, i) => (on(u).some((x) => directionOf(x.kind) === dir) ? i : -1)).filter((i) => i >= 0)
  const mirror = units.map((u, i) => (on(u).some((x) => directionOf(x.kind) === opposite) ? i : -1)).filter((i) => i >= 0)
  const parents = new Set(units.flatMap(on).filter((x) => directionOf(x.kind) === dir && named(x.origin)).map((x) => x.origin))
  const normal = units.filter((u) => u.reportable && u.measured.includes(chrom) && !on(u).length).length
  const read = `${carry.length} of ${n} units carry this ${dir}, ${normal} read${normal === 1 ? 's' : ''} `
    + `chromosome ${chrom} normal`
    + (parents.size > 1 ? ', and they name different parents' : '')

  if (carry.some((i) => mirror.some((j) => j !== i))) {
    return {
      mechanism: 'mitotic-reciprocal',
      why: `${read}, and ${unitsWord(mirror.length)} of this embryo carr${mirror.length === 1 ? 'ies' : 'y'} the `
        + 'opposite change over the same stretch. One daughter cell gaining what another loses is a '
        + 'division error after fertilisation',
    }
  }
  if (carry.length === n && parents.size < 2) {
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
    why: `${read}: not every unit carries it and none reads it measured normal, so no origin is called`,
  }
}

/**
 * The origin of one change over each stretch it covers, against every unit of its embryo.
 *
 * A stretch is cut by the edges of every row of every unit on the chromosome (SAME_BREAK_BP joins
 * near-equal ones) and is read once for the embryo, so one stretch prints one word in every cell.
 * A change over stretches that read differently keeps each stretch with its own word: a copy
 * lacking its centromere end leaves a piece with no centromere, shared out unevenly at division, so
 * a whole loss in one cell and a part loss in its sister are both read as they are.
 */
export function stretchesOf(change: UnitChange, units: readonly UnitReading[]): Stretch[] {
  const dir = directionOf(change.kind)
  const own = { fromBp: change.whole ? 0 : change.startBp, toBp: change.whole ? TO_END : change.endBp }
  if (units.length < 2) return [{ ...own, mechanism: 'none', why: 'one unit of this embryo was read; one cell cannot say when a change arose' }]
  if (!dir) return [{ ...own, mechanism: 'none', why: 'a change with no direction carries no origin' }]
  const { snap, at } = cutsOf(change.chrom, units)
  const parts: Stretch[] = []
  if (!change.whole && snap(change.startBp) === snap(change.endBp)) {
    // Shorter than the join: read over its own extent.
    const covers: Covers = (x) => x.whole || (x.startBp <= change.endBp && change.startBp <= x.endBp)
    parts.push({ ...own, ...readStretch(change.chrom, dir, covers, units) })
  } else {
    at.forEach((a, k) => {
      const b = at[k + 1] ?? TO_END
      if (!change.whole && (a < snap(change.startBp) || b > snap(change.endBp))) return
      const covers: Covers = (x) => x.whole || (snap(x.startBp) <= a && b <= snap(x.endBp))
      parts.push({
        fromBp: change.whole ? a : Math.max(a, change.startBp),
        toBp: change.whole ? b : Math.min(b, change.endBp),
        ...readStretch(change.chrom, dir, covers, units),
      })
    })
  }
  // Neighbours with the same word and the same counts are one stretch.
  return parts.reduce<Stretch[]>((out, p) => {
    const last = out[out.length - 1]
    if (last && last.mechanism === p.mechanism && last.why === p.why) last.toBp = p.toBp
    else out.push({ ...p })
    return out
  }, [])
}

/** Whether a change's stretches print more than one word (no word counts as one). */
export const split = (stretches: readonly Stretch[]): boolean =>
  new Set(stretches.map((s) => ORIGIN_WORD[s.mechanism])).size > 1

/**
 * The origin of one change as one call: the word its stretches share, or none where they differ.
 * The reading is the embryo's, so `self` only names the unit the change is on.
 */
export function callMechanism(
  _self: number, change: UnitChange, units: readonly UnitReading[],
): { mechanism: Mechanism; why: string } {
  const st = stretchesOf(change, units)
  if (split(st)) return { mechanism: 'none', why: 'its stretches read differently' }
  return { mechanism: st[0].mechanism, why: [...new Set(st.map((s) => s.why))].join('. ') }
}

/** One change of a unit with its origin, located as the card row it came from. `stretches` always
 *  holds at least one; `mechanism` is the shared word, or none where the stretches differ. */
export interface TimedChange extends UnitChange { mechanism: Mechanism; why: string; stretches: Stretch[] }

/** Every change of every unit of one embryo, in the order the units were given, read among the units shown. */
export function timeGroup(units: readonly UnitReading[]): TimedChange[][] {
  const shown = units.filter((u) => u.shown !== false)
  return units.map((u) => u.changes.map((c) => {
    if (u.shown === false) {
      const why = 'the sample is not reportable, so no origin is read on it'
      return { ...c, mechanism: 'none' as const, why,
        stretches: [{ fromBp: c.whole ? 0 : c.startBp, toBp: c.whole ? TO_END : c.endBp, mechanism: 'none' as const, why }] }
    }
    const stretches = stretchesOf(c, shown)
    const one = !split(stretches)
    return {
      ...c, stretches,
      mechanism: one ? stretches[0].mechanism : 'none',
      why: one ? [...new Set(stretches.map((s) => s.why))].join('. ') : 'its stretches read differently',
    }
  }))
}

const mb = (bp: number): string => (bp / 1e6).toFixed(1)

/** A stretch's extent in words: "0.0 to 52.6 Mb", "52.6 Mb to the end". */
export const stretchText = (s: Pick<Stretch, 'fromBp' | 'toBp'>): string =>
  (s.toBp >= TO_END ? `${mb(s.fromBp)} Mb to the end` : `${mb(s.fromBp)} to ${mb(s.toBp)} Mb`)

/** The word for one stretch, or "not called" where the embryo reads none there. */
export const stretchWord = (s: Pick<Stretch, 'mechanism'>): string => ORIGIN_WORD[s.mechanism] || 'not called'

/** The call a stretch belongs to, as the page, the report and the download name it. */
export interface StretchOf { chrom: string; kind: string; startBp: number; endBp: number; whole?: boolean; locus?: string }

const wholeOf = (c: StretchOf): boolean => c.whole ?? c.locus === `chr${c.chrom}`

/**
 * Where a stretch is, the stretch first and the call named as what it belongs to, so a stretch of
 * a whole-chromosome call does not read as the whole chromosome: "chr17 0.0 to 52.6 Mb, part of the
 * whole-chromosome loss", "chr9 16.7 to 133.9 Mb, part of the loss over 7.2 to 139.4 Mb".
 */
export function stretchWhere(c: StretchOf, s: Pick<Stretch, 'fromBp' | 'toBp'>): string {
  const dir = directionOf(c.kind) ?? 'change'
  return `chr${c.chrom} ${stretchText(s)}, part of the `
    + (wholeOf(c) ? `whole-chromosome ${dir}` : `${dir} over ${mb(c.startBp)} to ${mb(c.endBp)} Mb`)
}

/** A change's origin as one line of text: its word, or each stretch with where it is and its word. */
export const originText = (d: StretchOf & { mechanism?: Mechanism; mechanismStretches?: readonly Stretch[] }): string =>
  (d.mechanismStretches ? d.mechanismStretches.map((s) => `${stretchWhere(d, s)}: ${stretchWord(s)}`).join('; ')
    : d.mechanism ? ORIGIN_WORD[d.mechanism] : '')

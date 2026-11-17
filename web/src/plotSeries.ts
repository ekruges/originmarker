/**
 * The per-marker series a call plot draws, and the geometry of the box over a call.
 *
 * Built from the accumulator the one scoring pass already filled, after it is scored: no second
 * read of the file and no scan. Bounded: per chromosome at most CAP_INFORMATIVE markers where a
 * loaded parent is homozygous (the ones that can show whose letter is missing) and CAP_OTHER of
 * the rest, each drawn by a seeded reservoir so a chromosome of any length keeps the same budget
 * and the same file always keeps the same points.
 */

/** Markers where a loaded parent is homozygous, per chromosome. With both parents loaded only
 *  about 3 in 100 markers are, so a shared budget would leave too few to show a parent's letter
 *  missing; they get a budget of their own. */
export const CAP_INFORMATIVE = 1500
/** Every other marker, per chromosome. A plot is at most about 1,000 px wide, so 4,000 points in
 *  all is about four per pixel column: enough for the bands to read, few enough to draw at once. */
export const CAP_OTHER = 2500

/** Sample genotype per marker: 0 no call, 1 AA, 2 AB, 3 BB. */
export const GT_TEXT = ['no call', 'AA', 'AB', 'BB'] as const
const GT_CODE: Record<string, number> = { AA: 1, AB: 2, BB: 3 }

/** Whose letter alone the sample shows at a marker: 0 not readable, 1 the father's, 2 the mother's. */
export type Who = 0 | 1 | 2

export interface ChromSeries {
  /** Sorted by position. */
  pos: Uint32Array
  /** B allele frequency; NaN where the file gives none. */
  baf: Float32Array
  /** log R; NaN where the file gives none. */
  lrr: Float32Array
  gt: Uint8Array
  who: Uint8Array
  /** Markers on the chromosome before the budget, so the plot can say what it shows. */
  total: number
  /** Last marker position on the chromosome, kept or not. */
  endBp: number
}

/** The fields of scoreSample's Collected this reads. */
export interface SeriesSource {
  markerPos: ReadonlyMap<string, { chrom: string; pos: number }>
  myBaf: ReadonlyMap<string, number>
  myGt: ReadonlyMap<string, string>
  cnByChrom: ReadonlyMap<string, readonly { pos: number; log2R: number | null }[]>
}

const hom = (g: string | undefined): boolean => g === 'AA' || g === 'BB'

/** Whose letter alone a homozygous call shows. A parent homozygous for the other letter has none
 *  of this one to give, so it came from the other parent. Both excluded is not Mendelian: 0. */
export function whoseLetter(sample: string | undefined, father?: string, mother?: string): Who {
  if (!hom(sample)) return 0
  const notFather = hom(father) && sample !== father
  const notMother = hom(mother) && sample !== mother
  if (notFather && notMother) return 0
  return notFather ? 2 : notMother ? 1 : 0
}

/** Seeded reservoir of indices into a stream. Deterministic per chromosome. */
class Reservoir {
  readonly keep: string[] = []
  seen = 0
  readonly cap: number
  private seed: number
  constructor(cap: number, seed: number) { this.cap = cap; this.seed = seed || 1 }
  private rnd(): number {
    this.seed = (this.seed * 1103515245 + 12345) & 0x7fffffff
    return this.seed / 0x80000000
  }
  offer(id: string): void {
    this.seen += 1
    if (this.keep.length < this.cap) { this.keep.push(id); return }
    const j = Math.floor(this.rnd() * this.seen)
    if (j < this.cap) this.keep[j] = id
  }
}

export function callSeries(
  src: SeriesSource,
  father: ReadonlyMap<string, string> | null,
  mother: ReadonlyMap<string, string> | null,
  caps: { informative: number; other: number } = { informative: CAP_INFORMATIVE, other: CAP_OTHER },
): Map<string, ChromSeries> {
  const pools = new Map<string, { inf: Reservoir; rest: Reservoir; end: number }>()
  let lastChrom = ''
  let p: { inf: Reservoir; rest: Reservoir; end: number } | undefined
  for (const [id, { chrom, pos }] of src.markerPos) {
    if (chrom !== lastChrom) { p = pools.get(chrom); lastChrom = chrom }
    if (!p) {
      // Seeded by the chromosome name, so the same file always keeps the same points.
      const s = [...chrom].reduce((a, c) => a * 31 + c.charCodeAt(0), 7)
      pools.set(chrom, p = { inf: new Reservoir(caps.informative, s), rest: new Reservoir(caps.other, s + 1), end: 0 })
    }
    if (pos > p.end) p.end = pos
    const informative = (father !== null && hom(father.get(id)) || mother !== null && hom(mother.get(id)))
      && src.myGt.has(id)
    ;(informative ? p.inf : p.rest).offer(id)
  }
  const out = new Map<string, ChromSeries>()
  for (const [chrom, p] of pools) {
    const ids = [...p.inf.keep, ...p.rest.keep]
      .map((id) => ({ id, pos: src.markerPos.get(id)!.pos }))
      .sort((a, b) => a.pos - b.pos)
    // log R is held per position on the copy-number channel; join the kept markers by position.
    const lrrAt = new Map<number, number>()
    for (const { pos } of ids) lrrAt.set(pos, NaN)
    for (const m of src.cnByChrom.get(chrom) ?? []) {
      if (m.log2R !== null && lrrAt.has(m.pos)) lrrAt.set(m.pos, m.log2R)
    }
    const n = ids.length
    const s: ChromSeries = {
      pos: new Uint32Array(n), baf: new Float32Array(n), lrr: new Float32Array(n),
      gt: new Uint8Array(n), who: new Uint8Array(n), total: p.inf.seen + p.rest.seen, endBp: p.end,
    }
    ids.forEach(({ id, pos }, i) => {
      const g = src.myGt.get(id)
      s.pos[i] = pos
      s.baf[i] = src.myBaf.get(id) ?? NaN
      s.lrr[i] = lrrAt.get(pos) ?? NaN
      s.gt[i] = g ? GT_CODE[g] ?? 0 : 0
      s.who[i] = whoseLetter(g, father?.get(id), mother?.get(id))
    })
    out.set(chrom, s)
  }
  return out
}

/** Median of the finite values in a list, or NaN. */
export function median(xs: readonly number[]): number {
  const v = xs.filter(Number.isFinite).sort((a, b) => a - b)
  if (!v.length) return NaN
  const m = v.length >> 1
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2
}

/** The cell's own log R level: the median over every kept marker of every chromosome. A line to
 *  read a box against, not the tool's reference, which is chosen by its own rule. */
export function cellLrr(series: ReadonlyMap<string, ChromSeries>): number {
  const all: number[] = []
  for (const s of series.values()) for (const v of s.lrr) all.push(v)
  return median(all)
}

/** Median log R of the kept markers inside [startBp, endBp]. */
export function boxLrr(s: ChromSeries, startBp: number, endBp: number): number {
  const v: number[] = []
  for (let i = 0; i < s.pos.length; i += 1) if (s.pos[i] >= startBp && s.pos[i] <= endBp) v.push(s.lrr[i])
  return median(v)
}

/** Index of the kept marker nearest a position, or -1 on an empty series. */
export function nearestIndex(pos: ArrayLike<number>, bp: number): number {
  const n = pos.length
  if (!n) return -1
  let lo = 0
  let hi = n - 1
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (pos[mid] < bp) lo = mid + 1
    else hi = mid
  }
  return lo > 0 && bp - pos[lo - 1] <= pos[lo] - bp ? lo - 1 : lo
}

/** One box over a chromosome plot. `lines` are what the row that opened it shows, for the hover. */
export interface PlotBox {
  startBp: number
  endBp: number
  /** Edge intervals, where the tool localised the edges. */
  startLoBp?: number
  startHiBp?: number
  endLoBp?: number
  endHiBp?: number
  withheld?: boolean
  label: string
  lines: string[]
}

/** The tool's own edges for a segment: refined edges with their intervals, or the window edges. */
export function segmentBox(sg: {
  startBp: number; endBp: number
  refined?: { startBp: number; startLoBp: number; startHiBp: number; endBp: number; endLoBp: number; endHiBp: number }
}, label: string, lines: string[], withheld = false): PlotBox {
  const r = sg.refined
  return r
    ? { startBp: r.startBp, endBp: r.endBp, startLoBp: r.startLoBp, startHiBp: r.startHiBp,
      endLoBp: r.endLoBp, endHiBp: r.endHiBp, label, lines, withheld }
    : { startBp: sg.startBp, endBp: sg.endBp, label, lines, withheld }
}

export interface BoxGeometry {
  /** The called stretch, in pixels. */
  x0: number
  x1: number
  /** The edge intervals, in pixels, drawn lighter. Absent where the edges were not localised. */
  left?: [number, number]
  right?: [number, number]
}

/** Pixels for a box on a plot spanning [0, endBp] over [left, left + width]. Clamped to the plot,
 *  and never narrower than 2 px, so a call shorter than a pixel is still seen. */
export function boxGeometry(b: PlotBox, endBp: number, left: number, width: number): BoxGeometry {
  const x = (bp: number) => left + width * Math.max(0, Math.min(1, bp / Math.max(endBp, 1)))
  let x0 = x(b.startBp)
  let x1 = x(b.endBp)
  if (x1 - x0 < 2) { const m = (x0 + x1) / 2; x0 = m - 1; x1 = m + 1 }
  const g: BoxGeometry = { x0, x1 }
  if (b.startLoBp !== undefined && b.startHiBp !== undefined) g.left = [x(b.startLoBp), x(b.startHiBp)]
  if (b.endLoBp !== undefined && b.endHiBp !== undefined) g.right = [x(b.endLoBp), x(b.endHiBp)]
  return g
}

/** Tick spacing in Mb for an axis of this length: about six ticks. */
export function tickStepMb(endBp: number): number {
  const mb = endBp / 1e6
  return [1, 2, 5, 10, 20, 25, 50, 100].find((s) => mb / s <= 7) ?? 100
}

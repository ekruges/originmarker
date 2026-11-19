/**
 * What the cell figure draws: how many copies each parent gave, at each place in the genome.
 *
 * One copy from each parent is a normal chromosome. A lost copy takes that parent's track to none, a
 * doubled copy takes it to two, a chromosome lost outright takes both to none, so losses and gains are
 * the same picture in opposite directions and whose copy it is is which track moves.
 *
 * A parent's copies are that parent's share of the DNA times the copies in all: the share comes from
 * the alleles at the markers where the parents carry different alleles and needs no reference level,
 * and only the total needs one. The product was measured against each part on the open public cells
 * with stated karyotypes (audit/truth/gain-figure-levers-2026-11-17.txt): it matches the better part
 * everywhere and beats the total alone in the cells carrying many changes, where the level the cell is
 * measured against is itself shifted.
 *
 * Built once per sample from the accumulator the one scoring pass filled, over every marker; only the
 * bins are kept, so memory is bounded by markers over bin size.
 *
 * Bins hold a fixed number of the track's own markers, so every bar is equally trustworthy and a thin
 * region never makes a ragged bar. The base bin is BASE markers; a figure merges whole base bins.
 */

/**
 * A base bin is a stretch of the genome, not a count of markers: the measurement this figure rests on
 * binned by position (audit/truth/gain-figure-levers-2026-11-17.txt), and a run of bars along a
 * chromosome is what the reader is asked to see. Markers where the parents carry different alleles are
 * a small share of the array, about 3 in 100 on the arrays shipped as examples, so bars of a fixed
 * marker count would each cover tens of megabases and the run would collapse into two or three blocks.
 */
export const BASE_BP = 2.5e6
/** Base bins per bar across the whole cell: the 10 Mb bar of the measured figure. */
export const CELL_PER = 4
/** Fewest markers a bar may hold. Under this the bar joins the one before it, so a thin stretch makes
 *  a wider bar rather than a bar that reads as a level on almost nothing. */
export const MIN_BIN_MARKERS = 8

/** Base bins per bar for a track on the whole-cell figure. */
export function cellPer(_t: Track): number {
  return CELL_PER
}

/** Fewest bars a stretch gets on its own figure. */
export const MIN_STRETCH_BARS = 6
export type Role = 'paternal' | 'maternal'

/**
 * The levels a track is read at, fixed from the labelled stretches before the figure was drawn and
 * never tuned afterwards: 117 stretches where a parent's copy is missing, 1,075 published-normal
 * tracks and 6 published gains, on 28 open cells with both parents measured
 * (audit/truth/gain-figure-levers-2026-11-17.txt, 17 November 2026). On these boundaries a missing
 * copy lands in "none" 96% of the time, a normal track in "one copy" 92%, and a doubled copy in
 * "two" five times of six.
 *
 * The two boundaries do not rest on the same signal. None against one copy rests on the alleles:
 * the share puts a lost copy outside the cell's own band on 86% of whole losses against 42% for
 * intensity alone. One copy against two rests on intensity, where a doubled copy shows as more DNA
 * and heterozygosity makes the comparison worse (audit/truth/het-and-copies-2026-11-17.txt).
 */
export const COPY_NONE_BELOW = 0.68
export const COPY_TWO_ABOVE = 1.65
/** The top of a track's axis: above four copies in one cell the reading is amplification, not count. */
export const COPY_MAX = 2.9

/** What a track's level says, in the figure's own words. */
export function copyWord(copies: number): 'none' | 'one copy' | 'two copies' {
  return copies < COPY_NONE_BELOW ? 'none' : copies > COPY_TWO_ABOVE ? 'two copies' : 'one copy'
}

export interface Bin {
  startBp: number
  endBp: number
  /** Markers in the bin. */
  n: number
  /** Share tracks: the parent's share of the DNA, added over the bin's markers. Level tracks: the
   *  median of the bin's readings (log R less the cell's own level, or copies). */
  k: number
}

/** One row of the figure. */
export interface Track {
  key: 'father' | 'mother' | 'call-rate' | 'intensity'
  /** In words, for the figure. */
  label: string
  /** How the bar is computed, in one sentence, for the hover and "How this was read". */
  tech: string
  /** The parent whose alleles it reads, for its colour; none for the parent-free tracks. */
  parent: Role | null
  /** 'share' bars from 0 to 1; 'level' bars that already hold their own value. */
  kind: 'share' | 'level'
  /** Drawn as a track of the figure: the two parents where both are measured, the total otherwise. */
  shown: boolean
  bins: Map<string, Bin[]>
}

export interface PlotData {
  tracks: Track[]
  /** Both parents measured, so copies can be split between them. */
  split: boolean
  /** Last marker position per chromosome. */
  endBp: Map<string, number>
}

/** The fields of scoreSample's Collected this reads. */
export interface PlotSource {
  markerPos: ReadonlyMap<string, { chrom: string; pos: number }>
  myGt: ReadonlyMap<string, string>
  myBaf: ReadonlyMap<string, number>
  cnByChrom: ReadonlyMap<string, readonly { pos: number; called: boolean; log2R: number | null }[]>
}

const hom = (g: string | undefined): g is 'AA' | 'BB' => g === 'AA' || g === 'BB'
const isAuto = (c: string) => /^\d+$/.test(c) && +c >= 1 && +c <= 22

export function median(xs: readonly number[]): number {
  const v = xs.filter(Number.isFinite).sort((a, b) => a - b)
  if (!v.length) return NaN
  const m = v.length >> 1
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2
}

/** Base bins of BASE markers in position order; a last remainder of half a bin or less joins the bin
 *  before it. `level` gives the bin a median instead of a count. */
export function binBase(pos: Float64Array, x: Float64Array, level = false, size = BASE_BP): Bin[] {
  const o = [...pos.keys()].sort((a, b) => pos[a] - pos[b])
  const out: Bin[] = []
  let i = 0
  while (i < o.length) {
    const startBp = Math.floor(pos[o[i]] / size) * size
    let edge = startBp + size
    let j = i
    // One bin is one stretch of the genome, widened a stretch at a time until it holds enough markers.
    while (j < o.length) {
      if (pos[o[j]] < edge) { j += 1; continue }
      if (j - i >= MIN_BIN_MARKERS) break
      edge += size
    }
    const part = o.slice(i, j)
    const k = level ? median(part.map((t) => x[t])) : part.reduce((s, t) => s + x[t], 0)
    out.push({ startBp, endBp: edge, n: part.length, k })
    i = j
  }
  // A last bin too thin to stand on its own joins the one before it.
  if (out.length > 1 && out[out.length - 1].n < MIN_BIN_MARKERS) {
    const last = out.pop()!
    const prev = out[out.length - 1]
    out[out.length - 1] = { startBp: prev.startBp, endBp: last.endBp, n: prev.n + last.n,
      k: level ? (prev.k * prev.n + last.k * last.n) / (prev.n + last.n) : prev.k + last.k }
  }
  return out
}


/** Bars of `per` base bins each. Counts add exactly; a level bar takes the median of its bins. */
export function mergeBins(bins: readonly Bin[], per: number, level = false): Bin[] {
  if (per <= 1) return [...bins]
  const out: Bin[] = []
  for (let i = 0; i < bins.length; i += per) {
    const g = bins.slice(i, i + per)
    if (out.length && g.length < per / 2) {
      const last = out[out.length - 1]
      const all = [last, ...g]
      out[out.length - 1] = { startBp: last.startBp, endBp: g[g.length - 1].endBp, n: all.reduce((s, b) => s + b.n, 0),
        k: level ? median(all.map((b) => b.k)) : all.reduce((s, b) => s + b.k, 0) }
      continue
    }
    out.push({ startBp: g[0].startBp, endBp: g[g.length - 1].endBp, n: g.reduce((s, b) => s + b.n, 0),
      k: level ? median(g.map((b) => b.k)) : g.reduce((s, b) => s + b.k, 0) })
  }
  return out
}

/** The bar's value: the share for a share track, the level for an intensity track. */
export const valueOf = (t: Pick<Track, 'kind'>, b: Bin): number => (t.kind === 'share' ? b.k / b.n : b.k)

/**
 * The tracks for one sample.
 *
 * Both parents measured: each parent's share of the DNA at the markers where the parents carry
 * different alleles, which the figure turns into that parent's copies. One parent or none: the
 * copies in all and no claim about whose. With one parent the share cannot be built from the
 * markers that matter, and his track then reads three quarters of a copy where a copy is missing
 * against none with both parents, so the figure draws no parent track rather than a level it cannot
 * measure (audit/truth/gain-figure-levers-2026-11-17.txt).
 */
export function plotData(src: PlotSource, x: {
  loaded: { role: Role; gt: ReadonlyMap<string, string> } | null
  second: ReadonlyMap<string, string> | null
  genomeLrr?: number
}): PlotData {
  type Acc = Map<string, { p: number[]; x: number[] }>
  const acc = new Map<string, Acc>()
  const put = (t: string, chrom: string, pos: number, v: number) => {
    let a = acc.get(t)
    if (!a) acc.set(t, a = new Map())
    let c = a.get(chrom)
    if (!c) a.set(chrom, c = { p: [], x: [] })
    c.p.push(pos); c.x.push(v)
  }
  const endBp = new Map<string, number>()
  const loaded = x.loaded
  const other: Role | null = loaded ? (loaded.role === 'paternal' ? 'maternal' : 'paternal') : null
  const split = !!(loaded && x.second)
  for (const [id, { chrom, pos }] of src.markerPos) {
    if (!isAuto(chrom)) continue
    if (pos > (endBp.get(chrom) ?? 0)) endBp.set(chrom, pos)
    if (!split) continue
    // Only a marker where the parents carry DIFFERENT alleles says whose DNA is here: there the
    // loaded parent's allele fraction is his share of what the array read.
    const p = loaded!.gt.get(id)
    const q = x.second!.get(id)
    if (!hom(p) || !hom(q) || p === q) continue
    const baf = src.myBaf.get(id)
    if (baf === undefined || !Number.isFinite(baf)) continue
    const mine = p === 'BB' ? baf : 1 - baf
    put(loaded!.role, chrom, pos, Math.max(0, Math.min(1, mine)))
    put(other!, chrom, pos, Math.max(0, Math.min(1, 1 - mine)))
  }
  {
    // The cell's own level: the median over every autosomal reading. The call rate is kept in every
    // setup, because it is the only guard on a level of none.
    const all: number[] = []
    for (const [chrom, ms] of src.cnByChrom) if (isAuto(chrom)) for (const m of ms) if (m.log2R !== null) all.push(m.log2R)
    const sorted = Float64Array.from(all).sort()
    const normal = x.genomeLrr ?? (sorted.length ? sorted[sorted.length >> 1] : 0)
    for (const [chrom, ms] of src.cnByChrom) {
      if (!isAuto(chrom)) continue
      for (const m of ms) {
        put('call-rate', chrom, m.pos, m.called ? 1 : 0)
        if (m.log2R !== null) put('intensity', chrom, m.pos, m.log2R - normal)
      }
    }
  }
  const NAME: Record<Role, string> = { paternal: 'father', maternal: 'mother' }
  const spec = (key: string): Omit<Track, 'bins' | 'shown'> | null => {
    if (key === 'paternal' || key === 'maternal') {
      const who = NAME[key as Role]
      return {
        key: who as 'father' | 'mother', parent: key as Role, kind: 'share',
        label: `Copies from the ${who}`,
        tech: `Each bar is the share of the DNA that is the ${who}'s, at the markers where the parents `
          + 'carry different alleles, times how many copies are here in all. One from each parent is a '
          + 'normal chromosome.',
      }
    }
    if (key === 'call-rate') {
      return { key, parent: null, kind: 'share', label: 'Markers read',
        tech: 'Each bar is the share of markers the array read at all. It is not drawn as a track: it is '
          + 'the guard on a level of none, because a stretch with no DNA left still returns calls.' }
    }
    if (key === 'intensity') {
      return { key, parent: null, kind: 'level', label: 'Copies in all',
        tech: 'Each bar is how many copies are here in all, from the signal strength against this cell\'s '
          + 'own level (median log R ratio). Two in all is a normal chromosome.' }
    }
    return null
  }
  const order = ['paternal', 'maternal', 'call-rate', 'intensity']
  const tracks: Track[] = []
  for (const key of order) {
    const a = acc.get(key)
    const sp = spec(key)
    if (!a || !sp) continue
    const bins = new Map<string, Bin[]>()
    for (const [chrom, c] of a) bins.set(chrom, binBase(Float64Array.from(c.p), Float64Array.from(c.x), sp.kind === 'level'))
    tracks.push({ ...sp, shown: sp.key === 'call-rate' ? false : split ? sp.key !== 'intensity' : sp.key === 'intensity', bins })
  }
  return { tracks, split, endBp }
}

/**
 * One track's copies along a chromosome: the parent's share of the DNA times the copies in all,
 * both over the same bars. `level` is the intensity track's own bins, in position order.
 */
export function copyBars(share: readonly Bin[], level: readonly Bin[]): Bin[] {
  const out: Bin[] = []
  let i = 0
  for (const b of share) {
    while (i < level.length && (level[i].startBp + level[i].endBp) / 2 < b.startBp) i += 1
    const v: number[] = []
    for (let j = i; j < level.length && (level[j].startBp + level[j].endBp) / 2 <= b.endBp; j += 1) v.push(level[j].k)
    const lvl = median(v)
    const copies = Number.isFinite(lvl) ? 2 * 2 ** lvl : NaN
    out.push({ ...b, k: (b.n ? b.k / b.n : NaN) * copies })
  }
  return out
}

/** The copies in all along a chromosome, for the setups that name no parent. */
export const totalBars = (level: readonly Bin[]): Bin[] =>
  level.map((b) => ({ ...b, k: 2 * 2 ** b.k }))

/**
 * Bars where both parents' tracks read none while the markers are still calling at the rate this
 * cell holds elsewhere. A stretch with no template reads background in both channels and the caller
 * still emits calls, at 0.665 heterozygous against 0.151 at two copies, so neither intensity nor the
 * alleles can say "none" on their own; what marks it is the call rate collapsing, 0.438 against
 * 0.852 (audit/truth/het-and-copies-2026-11-17.txt). The figure draws these bars as an empty
 * outline rather than a confident zero.
 */
export function guardedBars(
  father: readonly Bin[], mother: readonly Bin[], rate: readonly Bin[], floor: number,
): Set<number> {
  const out = new Set<number>()
  if (!Number.isFinite(floor)) return out
  for (let i = 0; i < father.length; i += 1) {
    const f = father[i]; const m = mother[i]
    if (!m || !(f.k < COPY_NONE_BELOW) || !(m.k < COPY_NONE_BELOW)) continue
    const over = rate.filter((b) => (b.startBp + b.endBp) / 2 >= f.startBp && (b.startBp + b.endBp) / 2 <= f.endBp)
    const r = median(over.map((b) => (b.n ? b.k / b.n : NaN)))
    if (Number.isFinite(r) && r >= floor) out.add(i)
  }
  return out
}

/** The middle 80% of a track's bars on the chromosomes with no call: the cell's own normal range. */
export function normalRange(t: Track, per: number, skip: ReadonlySet<string>): { lo: number; mid: number; hi: number } {
  const v: number[] = []
  for (const [c, bs] of t.bins) if (!skip.has(c)) for (const b of mergeBins(bs, per, t.kind === 'level')) v.push(valueOf(t, b))
  v.sort((a, b) => a - b)
  const at = (q: number) => (v.length ? v[Math.min(v.length - 1, Math.floor(q * v.length))] : NaN)
  return { lo: at(0.1), mid: median(v), hi: at(0.9) }
}

/** A share over the bars that lie inside [startBp, endBp], pooled by marker counts. */
export function shareIn(bins: readonly Bin[], startBp: number, endBp: number): { share: number; n: number } {
  let k = 0; let n = 0
  for (const b of bins) if (b.startBp >= startBp && b.endBp <= endBp) { k += b.k; n += b.n }
  return { share: n ? k / n : NaN, n }
}

/** Base bins per bar on a stretch's own figure: as many as keep MIN_STRETCH_BARS bars inside it, up to
 *  the whole-cell bar. */
export function perForStretch(bins: readonly Bin[], startBp: number, endBp: number): number {
  const inside = bins.filter((b) => b.startBp >= startBp && b.endBp <= endBp).length
  return Math.max(1, Math.min(CELL_PER, Math.floor(inside / MIN_STRETCH_BARS)))
}

/** The stretch a figure shows: a call under a sixth of its chromosome in a window around it. */
export function viewOf(startBp: number, endBp: number, chromEndBp: number): { loBp: number; hiBp: number } {
  const span = endBp - startBp
  if (span >= chromEndBp / 6) return { loBp: 0, hiBp: chromEndBp }
  const half = Math.max(span * 3.5, 6e6)
  const mid = (startBp + endBp) / 2
  let loBp = Math.max(0, mid - half)
  let hiBp = Math.min(chromEndBp, mid + half)
  if (hiBp - loBp < 2 * half) { if (loBp === 0) hiBp = Math.min(chromEndBp, 2 * half); else loBp = Math.max(0, chromEndBp - 2 * half) }
  return { loBp, hiBp }
}

/** Most changes the figure zooms. Three zooms fit across a page and stay readable, so six is two
 *  rows; every change past that is named in one line instead. */
export const MAX_ZOOMS = 6

/**
 * The changes the figure zooms: the biggest called stretches, laid out in chromosome order.
 *
 * A reading the card does not count as a change never takes a zoom from one it does: the panels
 * below are there to show the changes the headline reports.
 */
export function zoomPicks<T extends { chrom: string; startBp: number; endBp: number; status?: string }>(
  calls: readonly T[], max = MAX_ZOOMS,
): { shown: T[]; rest: T[] } {
  const counts = (c: T) => (c.status === undefined || c.status === 'counted' ? 0 : 1)
  const by = [...calls].sort((a, b) => counts(a) - counts(b) || (b.endBp - b.startBp) - (a.endBp - a.startBp))
  const shown = by.slice(0, max)
  const rest = by.slice(max).sort((a, b) => +a.chrom - +b.chrom)
  shown.sort((a, b) => +a.chrom - +b.chrom || a.startBp - b.startBp)
  return { shown, rest }
}

/** Tick spacing in Mb for an axis of this length: about six ticks. */
export function tickStepMb(spanBp: number): number {
  const mb = spanBp / 1e6
  return [1, 2, 5, 10, 20, 25, 50, 100].find((s) => mb / s <= 7) ?? 100
}

/** Greedy word wrap at about `chars` per line. */
export function wrap(text: string, chars: number): string[] {
  const out: string[] = []
  let line = ''
  for (const w of text.split(' ')) {
    if (line && line.length + 1 + w.length > chars) { out.push(line); line = w } else line = line ? `${line} ${w}` : w
  }
  if (line) out.push(line)
  return out
}

/**
 * What the cell figures draw: for each parent whose alleles can be read, the share of the markers where
 * that parent's allele could be seen at which it is, in bins along every chromosome.
 *
 * Built once per sample from the accumulator the one scoring pass filled, over every marker; only the
 * bins are kept, so memory is bounded by markers over bin size. Measured on the open public cells with
 * stated karyotypes (audit/truth/call-plot-levers-2026-11-17.txt).
 *
 * Bins hold a fixed number of the track's own markers, so every bar is equally trustworthy and a thin
 * region never makes a ragged bar. The base bin is BASE markers; a figure merges whole base bins.
 */

/** Markers in a base bin. Figures merge these into bars of BASE * k markers. */
export const BASE = 50
/** Markers per bar on the whole-cell figure: equal-count bins of 200 separated stated losses in noisy
 *  cells better than 5 Mb bins or smaller counts (0.47 against 0.30 of inside bars past every other). */
export const CELL_BAR = 200
/** Most bars a track draws across the whole cell. A dense track (heterozygous calls, call rate) gets bars of
 *  more markers, so the figure stays one bar per few pixels and each bar is no less trustworthy. */
export const CELL_MAX_BARS = 300

/** Base bins per bar for a track on the whole-cell figure. */
export function cellPer(t: Track): number {
  let n = 0
  for (const bs of t.bins.values()) for (const b of bs) n += b.n
  return Math.max(CELL_BAR / BASE, Math.ceil(n / (BASE * CELL_MAX_BARS)))
}

/** Fewest bars a stretch gets on its own figure; the bar is shrunk to BASE markers at most to reach it. */
export const MIN_STRETCH_BARS = 6
/** Under this a stretch is not drawn on its own: no measure separated one band (about 6 Mb) by eye. */
export const MIN_DRAWN_SPAN_BP = 6e6

export type Role = 'paternal' | 'maternal'

export interface Bin {
  startBp: number
  endBp: number
  /** Markers in the bin. */
  n: number
  /** Share tracks: markers where the allele is found. Intensity: the median log R less the normal. */
  k: number
}

/** One row of the figure. */
export interface Track {
  key: 'father' | 'mother' | 'two-letter' | 'call-rate' | 'intensity'
  /** In words, for the figure. */
  label: string
  /** How the bar is computed, in one sentence, for the hover and "How this was read". */
  tech: string
  /** The row stands in for a parent's copy without measuring it directly. */
  indirect?: boolean
  /** The parent whose alleles it counts, for its colour; none for the parent-free tracks. */
  parent: Role | null
  /** 'share' bars from 0 to 1; 'level' bars around 0. */
  kind: 'share' | 'level'
  /** Drawn as a row of the whole-cell figure. Intensity is kept for every run, for gains and the
   *  one-line notes, and drawn only with no parent loaded. */
  shown: boolean
  bins: Map<string, Bin[]>
}

export interface PlotData {
  tracks: Track[]
  /** The cell's rate of false heterozygous calls, where both parents are measured; NaN otherwise. */
  dropIn: number
  /** Last marker position per chromosome. */
  endBp: Map<string, number>
}

/** The fields of scoreSample's Collected this reads. */
export interface PlotSource {
  markerPos: ReadonlyMap<string, { chrom: string; pos: number }>
  myGt: ReadonlyMap<string, string>
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
export function binBase(pos: Float64Array, x: Float64Array, level = false, size = BASE): Bin[] {
  const o = [...pos.keys()].sort((a, b) => pos[a] - pos[b])
  const out: Bin[] = []
  for (let i = 0; i < o.length;) {
    let j = Math.min(o.length, i + size)
    if (o.length - j <= size / 2) j = o.length
    const part = o.slice(i, j)
    const k = level ? median(part.map((t) => x[t])) : part.reduce((s, t) => s + x[t], 0)
    out.push({ startBp: pos[part[0]], endBp: pos[part[part.length - 1]], n: part.length, k })
    i = j
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
 * Both parents measured: each parent's allele found where the parents are opposite homozygotes. One
 * parent loaded: that parent's allele found where the other allele shows; the other parent's allele
 * found where the run's other embryos read her as the opposite homozygote (`otherOpposite`), or else
 * the heterozygous calls. No parent: call rate, heterozygous calls and intensity.
 */
export function plotData(src: PlotSource, x: {
  loaded: { role: Role; gt: ReadonlyMap<string, string> } | null
  second: ReadonlyMap<string, string> | null
  /** Probe ids where the run's other draws read the other parent as the loaded one's opposite homozygote. */
  otherOpposite?: ReadonlySet<string> | null
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
  let sameHom = 0; let sameHomAB = 0
  for (const [id, { chrom, pos }] of src.markerPos) {
    if (!isAuto(chrom)) continue
    if (pos > (endBp.get(chrom) ?? 0)) endBp.set(chrom, pos)
    const s = src.myGt.get(id)
    if (!loaded) continue
    const p = loaded.gt.get(id)
    if (!s || !hom(p)) continue
    if (x.second) {
      const q = x.second.get(id)
      if (!hom(q)) continue
      if (q === p) { sameHom += 1; if (s === 'AB') sameHomAB += 1; continue }
      put(loaded.role, chrom, pos, s === 'AB' || s === p ? 1 : 0)
      put(other!, chrom, pos, s === 'AB' || s === q ? 1 : 0)
      continue
    }
    if (s !== p) put(loaded.role, chrom, pos, s === 'AB' ? 1 : 0)
    if (x.otherOpposite?.has(id)) put(other!, chrom, pos, s !== p ? 1 : 0)
  }
  if (!loaded || (!x.second && !x.otherOpposite)) {
    for (const [id, { chrom, pos }] of src.markerPos) {
      if (!isAuto(chrom)) continue
      const s = src.myGt.get(id)
      if (s) put('two-letter', chrom, pos, s === 'AB' ? 1 : 0)
    }
  }
  {
    // The cell's own level: the median over every autosomal reading.
    const all: number[] = []
    for (const [chrom, ms] of src.cnByChrom) if (isAuto(chrom)) for (const m of ms) if (m.log2R !== null) all.push(m.log2R)
    const sorted = Float64Array.from(all).sort()
    const normal = x.genomeLrr ?? (sorted.length ? sorted[sorted.length >> 1] : 0)
    for (const [chrom, ms] of src.cnByChrom) {
      if (!isAuto(chrom)) continue
      for (const m of ms) {
        if (!loaded) put('call-rate', chrom, m.pos, m.called ? 1 : 0)
        if (m.log2R !== null) put('intensity', chrom, m.pos, m.log2R - normal)
      }
    }
  }
  const NAME: Record<Role, string> = { paternal: 'father', maternal: 'mother' }
  const spec = (key: string): Omit<Track, 'bins' | 'shown'> | null => {
    if (key === 'paternal' || key === 'maternal') {
      const who = NAME[key as Role]
      const both = !!x.second
      const own = loaded?.role === key
      return {
        key: who as 'father' | 'mother', parent: key as Role, kind: 'share',
        label: `${who[0].toUpperCase()}${who.slice(1)}'s copy`,
        tech: both ? `Each bar is the share of markers where the parents carry different alleles at which the ${who}'s allele is found.`
          : own ? `Each bar is the share of markers where the ${who} is homozygous and the other parent's allele shows at which the ${who}'s allele is found too.`
            : `Each bar is the share of markers where the run's other cells read the ${who} as homozygous for the allele the loaded parent lacks at which that allele is found.`,
      }
    }
    if (key === 'two-letter') {
      // With one parent loaded and nothing reading the other, heterozygous calls stand for her copy.
      const stands = loaded && !x.second && !x.otherOpposite ? NAME[loaded.role === 'paternal' ? 'maternal' : 'paternal'] : null
      return {
        key, parent: null, kind: 'share', indirect: !!stands,
        label: stands ? `${stands[0].toUpperCase()}${stands.slice(1)}'s copy, indirect` : 'Heterozygous calls',
        tech: `Each bar is the share of read markers that carry two alleles${stands
          ? `; a heterozygous call needs the ${stands}'s allele, so with one parent loaded this row stands in for the ${stands}'s copy` : ''}.`,
      }
    }
    if (key === 'call-rate') return { key, parent: null, kind: 'share', label: 'DNA present', tech: 'Each bar is the share of markers the array read at all.' }
    if (key === 'intensity') return { key, parent: null, kind: 'level', label: 'Copy number', tech: 'Each bar is the signal strength against this cell\'s own level (median log R ratio).' }
    return null
  }
  const order = ['paternal', 'maternal', 'two-letter', 'call-rate', 'intensity']
  const tracks: Track[] = []
  for (const key of order) {
    const a = acc.get(key)
    const sp = spec(key)
    if (!a || !sp) continue
    const bins = new Map<string, Bin[]>()
    for (const [chrom, c] of a) bins.set(chrom, binBase(Float64Array.from(c.p), Float64Array.from(c.x), sp.kind === 'level'))
    tracks.push({ ...sp, shown: sp.key !== 'intensity' || !loaded, bins })
  }
  return { tracks, dropIn: sameHom > 200 ? sameHomAB / sameHom : NaN, endBp }
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
  return Math.max(1, Math.min(CELL_BAR / BASE, Math.floor(inside / MIN_STRETCH_BARS)))
}

/**
 * Evidence that the parent's copy is missing, added along the chromosome and held at zero from below,
 * so it is flat where the copy is there and climbs where it is not. Per marker the log ratio of the
 * marker's reading with the copy missing (present at `missing`, the cell's false heterozygous rate) to with
 * it there (present at `present`, the cell's normal level); summed per bar from its counts.
 */
export function runningScore(bins: readonly Bin[], present: number, missing: number): number[] {
  const r1 = Math.max(missing, 0.02)
  if (!(present > r1 && present < 1)) return bins.map(() => NaN)
  const up = Math.log(r1 / present)
  const dn = Math.log((1 - r1) / (1 - present))
  let s = 0
  return bins.map((b) => (s = Math.max(0, s + b.k * up + (b.n - b.k) * dn)))
}

/**
 * The top of a share row's axis. A row whose normal level sits under a quarter (heterozygous calls) is
 * drawn against twice that level, so its normal range sits mid-row instead of in a thin strip at the
 * floor; no value changes, only where the axis ends. A row whose normal level is high keeps 0 to 100%,
 * since its drop toward zero is the event and a cut axis would make a partial drop read as a full one.
 */
export function shareTop(normal: number): number {
  return Number.isFinite(normal) && normal < 0.25 ? Math.min(1, Math.max(0.05, 2 * normal)) : 1
}

/** Whether a running score climbs inside the called stretch: by half its height or more between the
 *  bar before the stretch and its last bar. A score that rises only outside the stretch is not drawn. */
export function risesInside(score: readonly number[], bins: readonly Bin[], startBp: number, endBp: number): boolean {
  const inside = bins.map((b, i) => (b.startBp >= startBp && b.endBp <= endBp ? i : -1)).filter((i) => i >= 0)
  if (!inside.length) return false
  const top = Math.max(0, ...score.filter(Number.isFinite))
  const from = inside[0] > 0 ? score[inside[0] - 1] : 0
  const rise = (score[inside[inside.length - 1]] || 0) - (from || 0)
  return top > 0 && rise >= top / 2
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

/**
 * Whether a call's own figure is drawn, from what separated on the open public cells. Not drawn: a
 * stretch under MIN_DRAWN_SPAN_BP; a gain on part of a chromosome; a partial loss with no parent; the
 * other parent's copy on part of a chromosome where only heterozygous calls stand for her.
 */
export function drawable(call: { whole: boolean; spanBp: number; gain: boolean; parent: Role | null },
  data: PlotData): { draw: boolean; track: Track | null } {
  const pick = (k: Track['key']) => data.tracks.find((t) => t.key === k) ?? null
  const own = call.parent === 'paternal' ? pick('father') : call.parent === 'maternal' ? pick('mother') : null
  if (call.whole) return { draw: true, track: own ?? pick('intensity') ?? pick('two-letter') }
  if (call.spanBp < MIN_DRAWN_SPAN_BP || call.gain) return { draw: false, track: null }
  if (own) return { draw: true, track: own }
  return { draw: false, track: null }
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

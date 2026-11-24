/**
 * Events built at known places into the parsed rows of a real array, for the known-answer list.
 *
 * An event is a specification, applied in memory to the rows of a base cell just before scoring.
 * Nothing is written per event. Rows outside [start, end) of the event's chromosome are returned as
 * the same objects, so a chromosome the event does not touch is identical to the base.
 *
 * INTENSITY. Every marker inside with a finite log2R moves by exactly `amount`.
 *
 * WHAT A MARKER AT A NEW INTENSITY LOOKS LIKE, with the cell's level table (levelTable): in this
 * material the no-call rate, and how far a call's B-allele frequency sits from where the call puts
 * it, both follow a marker's own log2R. So a moved marker is read as the cell's own reference
 * markers read at its new level (a marker that stays in its bin keeps its own): its B-allele frequency takes the value at the same rank among the
 * cell's markers of the same call in the new level's bin, and it becomes a no-call, or a no-call
 * becomes called, with the probability that turns the no-call rate of its old level into that of its
 * new one, taking the B-allele frequency of a no-call at the new level, on the side it was. A no-call
 * moved to a new level keeps its rank among the cell's no-calls there. A no-call that becomes called in a gain is heterozygous at the new level's heterozygous
 * share, else the homozygote its B-allele frequency is nearer; in a loss it is the kept allele.
 *
 * The B-allele frequencies are read at the level base + `bafLevel` x amount (default 1): in some
 * material they spread faster than log2R alone says, and that factor is set by measuring real losses.
 *
 * GAIN WITH A PARENT. Where `parent` names whose copy was added and a loaded parent's genotype is
 * given, the alleles are re-weighted as audit/validate/constructed-gains.ts does for whole
 * chromosomes: `alleleWeight` of one copy's worth, w (A + B) / 2, is added to the allele that parent
 * passed where that can be read, split evenly where it cannot, and B / (A + B) recomputed, which is
 * (BAF + w / 2) / (1 + w / 2) toward B and BAF / (1 + w / 2) toward A. With no parent no allele
 * ratio is built.
 *
 * LOSS. A heterozygous call stays heterozygous with probability `hetKept`, keeping its own B-allele
 * frequency or, with keptBaf 'kept', reading as the kept allele's homozygote does; else it becomes the
 * homozygote of the allele the kept copy carries (the parent not lost, where the loaded parent's
 * genotype shows it, otherwise drawn at random), with a B-allele frequency drawn from the cell's
 * homozygous markers at the new level. With both parents loaded, at their opposite homozygotes a
 * homozygous call of the lost parent's allele becomes the kept allele's homozygote. Where `noCallRatio` is given the no-call rate inside is
 * multiplied by it instead of following the level table.
 *
 * Which allele each parent passed, at one marker, with one parent's genotype `p` loaded and the cell's
 * own call `c`: where both parents are loaded and both homozygous, each passed its own allele; where p is homozygous the loaded parent passed p's allele, and the other parent passed
 * the other allele where c is heterozygous and the same allele where c is homozygous. Where p is
 * heterozygous and c homozygous both passed c's allele. Otherwise it cannot be read.
 */

export interface Row {
  probesetId: string; chrom: string; pos: number; log2R: number | null; baf: number | null
  copyNumber: number | null; genotype: 'AA' | 'AB' | 'BB' | 'NC'; unreadableGenotype: boolean
  rawGenotype: string; bestProbeset: boolean
}

export interface EventSpec {
  id: string
  chrom: string
  kind: 'gain' | 'loss'
  /** The interval changed, [start, end) in bp. */
  start: number
  end: number
  /** log2R added to every marker inside. */
  amount: number
  /** Gain: whose copy was added, against the loaded parent. Loss: whose copy was lost. Null: unknown. */
  parent?: 'loaded' | 'other' | null
  /** Gain with a parent: the share of one copy's worth added on the allele side. */
  alleleWeight?: number
  /** Loss: probability a heterozygous call stays heterozygous. */
  hetKept?: number
  /** Loss: factor on the no-call rate inside, in place of the level table. */
  noCallRatio?: number
  /** The B-allele frequency is read at the level base + bafLevel x amount (default 1). */
  bafLevel?: number
  /** Loss: a heterozygous call that stays heterozygous keeps its own B-allele frequency ('own',
   *  the default) or takes the kept allele's, as a hemizygous marker reads ('kept'). */
  keptBaf?: 'own' | 'kept'
  seed: number
}

/** What the builder reads from the base cell beyond its rows. */
export interface BaseContext {
  /** The loaded parent's genotype by probe id, or null with no parent. */
  parentGt: { get(id: string): string | undefined } | null
  /** The second parent's genotype, where both are loaded. */
  otherGt?: { get(id: string): string | undefined } | null
  /** B-allele frequencies of the cell's own homozygous calls on its reference chromosomes. */
  homBaf: { AA: Float64Array; BB: Float64Array }
  /** How the cell's own reference markers read at each level of log2R (levelTable). */
  levels?: LevelTable
}

/**
 * How a cell's own markers read at each level of log2R, from its reference chromosomes: per bin of
 * LEVEL_WIDTH, the no-call rate, the heterozygous share of calls, and the sorted distances of the
 * B-allele frequency from where its call puts it (from the pole for a homozygote, from 0.5 for a
 * heterozygote and for a no-call). A bin with under LEVEL_MIN markers reads as the nearest bin toward 0 that has them.
 */
export interface LevelTable { nc: Float64Array; ab: Float64Array; hom: Float64Array[]; het: Float64Array[]; noCall: Float64Array[]; use: Int32Array }
export const LEVEL_LO = -3
export const LEVEL_WIDTH = 0.1
export const LEVEL_BINS = 60
export const LEVEL_MIN = 200

export function levelTable(rows: readonly Row[], keep: (r: Row) => boolean): LevelTable {
  const n = new Float64Array(LEVEL_BINS), nc = new Float64Array(LEVEL_BINS), ab = new Float64Array(LEVEL_BINS)
  const hom: number[][] = Array.from({ length: LEVEL_BINS }, () => [])
  const het: number[][] = Array.from({ length: LEVEL_BINS }, () => [])
  const noCall: number[][] = Array.from({ length: LEVEL_BINS }, () => [])
  for (const r of rows) {
    if (!keep(r) || r.log2R === null || !Number.isFinite(r.log2R)) continue
    const b = rawBin(r.log2R)
    n[b]++
    if (r.genotype === 'NC') {
      nc[b]++
      if (r.baf !== null && Number.isFinite(r.baf)) noCall[b].push(Math.abs(r.baf - 0.5))
      continue
    }
    if (r.genotype === 'AB') ab[b]++
    if (r.baf === null || !Number.isFinite(r.baf)) continue
    if (r.genotype === 'AB') het[b].push(Math.abs(r.baf - 0.5))
    else hom[b].push(r.genotype === 'AA' ? r.baf : 1 - r.baf)
  }
  const ok = (b: number) => n[b] >= LEVEL_MIN && hom[b].length >= 20 && het[b].length >= 20 && noCall[b].length >= 20
  const zero = rawBin(0)
  const use = new Int32Array(LEVEL_BINS)
  for (let b = 0; b < LEVEL_BINS; b++) {
    let c = b
    while (!ok(c) && c !== zero) c += c < zero ? 1 : -1
    use[b] = c
  }
  const sorted = (xs: number[]) => Float64Array.from(xs).sort()
  return { nc: nc.map((x, b) => (n[b] ? x / n[b] : 0)), ab: ab.map((x, b) => (n[b] - nc[b] ? x / (n[b] - nc[b]) : 0)),
    hom: hom.map(sorted), het: het.map(sorted), noCall: noCall.map(sorted), use }
}
const rawBin = (l: number) => Math.min(LEVEL_BINS - 1, Math.max(0, Math.floor((l - LEVEL_LO) / LEVEL_WIDTH)))
/** Rank of v in a sorted array, as a share in (0, 1). */
function rankOf(a: Float64Array, v: number): number {
  let lo = 0, hi = a.length
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < v) lo = m + 1; else hi = m }
  let up = lo
  while (up < a.length && a[up] === v) up++
  return (lo + up + 1) / 2 / (a.length + 1)
}
const at = (a: Float64Array, q: number) => a[Math.min(a.length - 1, Math.max(0, Math.floor(q * a.length)))]

/** A small seeded generator, so every built event is reproducible from its seed. */
export function rng(seed: number): () => number {
  let a = seed >>> 0
  return () => {
    a = (a + 0x6d2b79f5) >>> 0
    let t = a
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

const other = (a: 'A' | 'B'): 'A' | 'B' => (a === 'A' ? 'B' : 'A')

/** The allele each parent passed at one marker, where it can be read. See the header. */
export function passed(p: string | undefined, c: string, o?: string): { loaded: 'A' | 'B'; other: 'A' | 'B' } | null {
  const hom = (g: string | undefined) => (g === 'AA' ? 'A' : g === 'BB' ? 'B' : null)
  const pa = hom(p); const ca = hom(c); const oa = hom(o)
  if (pa && oa) return { loaded: pa, other: oa }
  if (pa) {
    if (c === 'AB') return { loaded: pa, other: other(pa) }
    if (ca) return { loaded: pa, other: ca }
    return null
  }
  if (p === 'AB' && ca) return { loaded: ca, other: ca }
  return null
}

const RAW: Record<string, string> = { AA: '0', AB: '1', BB: '2', NC: '-1' }

/** The base rows with one event applied. The input is never changed. */
export function applyEvent(rows: readonly Row[], e: EventSpec, ctx: BaseContext): Row[] {
  const rand = rng(e.seed)
  const lv = ctx.levels
  const draw = (g: 'AA' | 'BB', bin: number) => {
    if (lv && bin >= 0) { const d = at(lv.hom[bin], rand()); return g === 'AA' ? d : 1 - d }
    const pool = ctx.homBaf[g]
    return pool.length ? pool[Math.floor(rand() * pool.length)] : (g === 'AA' ? 0 : 1)
  }
  const inside = (r: Row) => r.chrom === e.chrom && r.pos >= e.start && r.pos < e.end
  const ratio = e.noCallRatio
  let up = 0
  if (e.kind === 'loss' && ratio !== undefined && ratio > 1) {
    let n = 0, nc = 0
    for (const r of rows) if (inside(r)) { n++; if (r.genotype === 'NC') nc++ }
    up = noCallUp(n ? nc / n : 0, ratio)
  }
  const out = rows.slice()
  for (let i = 0; i < rows.length; i++) {
    const r = rows[i]
    if (!inside(r)) continue
    const finite = r.log2R !== null && Number.isFinite(r.log2R)
    const log2R = finite ? r.log2R! + e.amount : r.log2R
    const b0 = lv && finite ? lv.use[rawBin(r.log2R!)] : -1
    const b1 = lv && finite ? lv.use[rawBin(log2R!)] : -1
    const bb = lv && finite ? lv.use[rawBin(r.log2R! + (e.bafLevel ?? 1) * e.amount)] : -1
    const who = ctx.parentGt && e.parent ? passed(ctx.parentGt.get(r.probesetId), r.genotype, ctx.otherGt?.get(r.probesetId)) : null
    let genotype = r.genotype
    let baf = r.baf
    const u = rand(), v = rand(), w0 = rand()
    // The B-allele frequency reads as the cell's own markers of the same call read at the new level,
    // at the same rank.
    if (lv && b0 >= 0 && bb !== b0 && baf !== null && Number.isFinite(baf)) {
      if (genotype === 'NC') {
        const d = at(lv.noCall[bb], rankOf(lv.noCall[b0], Math.abs(baf - 0.5)))
        baf = baf >= 0.5 ? 0.5 + d : 0.5 - d
      } else if (genotype === 'AB') {
        // A heterozygous call keeps its own B-allele frequency: real lost copies in this material
        // keep their heterozygous calls in the middle band rather than spread as the cell's
        // low-intensity heterozygotes do.
      } else {
        const d = at(lv.hom[bb], rankOf(lv.hom[b0], genotype === 'AA' ? baf : 1 - baf))
        baf = genotype === 'AA' ? d : 1 - d
      }
    }
    if (e.kind === 'gain') {
      const w = e.alleleWeight ?? 0
      if (ctx.parentGt && e.parent && w > 0 && baf !== null && Number.isFinite(baf) && finite) {
        // One copy's worth, w (A + B) / 2, on the allele that parent passed: B / (A + B) recomputed.
        const a = who ? who[e.parent] : null
        const addB = a === 'B' ? w / 2 : a === 'A' ? 0 : w / 4
        baf = (baf + addB) / (1 + w / 2)
      }
    } else if (who && ctx.otherGt && e.parent && genotype !== 'NC' && genotype !== 'AB'
      && genotype !== ((e.parent === 'loaded' ? who.other : who.loaded) === 'A' ? 'AA' : 'BB')
      && ctx.parentGt!.get(r.probesetId) !== ctx.otherGt.get(r.probesetId)) {
      // Both parents loaded and opposite homozygotes: a homozygous call of the lost allele is the
      // kept allele's drop-out, and with the lost copy gone it reads as the kept allele.
      const k = e.parent === 'loaded' ? who.other : who.loaded
      genotype = k === 'A' ? 'AA' : 'BB'; baf = draw(genotype, bb)
    } else if (genotype === 'AB' && u < (e.hetKept ?? 0)) {
      if (e.keptBaf === 'kept') {
        const k = who && e.parent ? (e.parent === 'loaded' ? who.other : who.loaded) : (w0 < 0.5 ? 'A' : 'B')
        baf = draw(k === 'A' ? 'AA' : 'BB', bb)
      }
    } else if (genotype === 'AB') {
      // Loss. The allele the kept copy carries: the parent not lost, where it can be read.
      const k = who && e.parent ? (e.parent === 'loaded' ? who.other : who.loaded) : (w0 < 0.5 ? 'A' : 'B')
      genotype = k === 'A' ? 'AA' : 'BB'; baf = draw(genotype, bb)
    }
    // No-calls: a fixed ratio where the specification gives one, else as the cell's own markers
    // read at the new level.
    if (ratio !== undefined || !lv || b0 < 0) {
      if (e.kind === 'loss' && ratio !== undefined) {
        if (genotype !== 'NC' && v < up) genotype = 'NC'
        else if (genotype === 'NC' && r.genotype === 'NC' && ratio < 1 && v < 1 - ratio) {
          const k = who && e.parent ? (e.parent === 'loaded' ? who.other : who.loaded) : (w0 < 0.5 ? 'A' : 'B')
          genotype = k === 'A' ? 'AA' : 'BB'; baf = draw(genotype, bb)
        }
      }
    } else {
      const p0 = lv.nc[b0], p1 = lv.nc[b1]
      if (genotype !== 'NC' && p1 > p0 && v < (p1 - p0) / (1 - p0)) {
        genotype = 'NC'
        if (baf !== null && Number.isFinite(baf)) baf = 0.5 + (baf >= 0.5 ? 1 : -1) * at(lv.noCall[bb], rand())
      }
      else if (r.genotype === 'NC' && p1 < p0 && v < (p0 - p1) / p0) {
        const het = e.kind === 'gain' && w0 < lv.ab[b1]
        if (het) { genotype = 'AB'; baf = 0.5 + (u < 0.5 ? -1 : 1) * at(lv.het[bb], rand()) } else {
          const k = e.kind === 'loss' && who && e.parent ? (e.parent === 'loaded' ? who.other : who.loaded)
            : (baf !== null && Number.isFinite(baf) ? (baf < 0.5 ? 'A' : 'B') : (u < 0.5 ? 'A' : 'B'))
          genotype = k === 'A' ? 'AA' : 'BB'; baf = draw(genotype, bb)
        }
      }
    }
    out[i] = genotype === r.genotype ? { ...r, log2R, baf } : { ...r, log2R, baf, genotype, rawGenotype: RAW[genotype] }
  }
  return out
}

/**
 * The per-marker probability a called marker inside becomes a no-call, for a no-call rate inside of
 * `n0` multiplied by `ratio`: (ratio n0 - n0) / (1 - n0), capped to [0, 1].
 */
export function noCallUp(n0: number, ratio: number): number {
  if (!(ratio > 1) || !(n0 < 1)) return 0
  return Math.min(1, Math.max(0, (ratio * n0 - n0) / (1 - n0)))
}

// --- what a stretch of markers looks like, read the same way on real and built events ----------

/** B-allele-frequency bins: [0, 0.05), [0.05, 0.1), ... [0.95, 1]. */
export const BAF_BINS = 20

export interface Signature {
  markers: number
  /** Median log2R. */
  lrr: number
  /** Heterozygous calls over called markers. */
  abShare: number
  /** No-calls over markers. */
  noCall: number
  /** Share of markers with BAF in [0.2, 0.8], the tool's off-pole share. */
  offPole: number
  /** Share of markers with a finite BAF in each bin. */
  bafHist: number[]
}

export function signature(rows: readonly Row[], keep: (r: Row) => boolean): Signature {
  const lrr: number[] = []
  let n = 0, called = 0, ab = 0, nc = 0, withBaf = 0, off = 0
  const hist = new Array(BAF_BINS).fill(0)
  for (const r of rows) {
    if (!keep(r)) continue
    n++
    if (r.genotype === 'NC') nc++; else { called++; if (r.genotype === 'AB') ab++ }
    if (r.log2R !== null && Number.isFinite(r.log2R)) lrr.push(r.log2R)
    if (r.baf !== null && Number.isFinite(r.baf)) {
      withBaf++
      if (r.baf >= 0.2 && r.baf <= 0.8) off++
      hist[Math.min(BAF_BINS - 1, Math.floor(r.baf * BAF_BINS))]++
    }
  }
  lrr.sort((a, b) => a - b)
  return { markers: n, lrr: lrr.length ? lrr[lrr.length >> 1] : NaN, abShare: called ? ab / called : NaN,
    noCall: n ? nc / n : NaN, offPole: withBaf ? off / withBaf : NaN,
    bafHist: hist.map((h) => (withBaf ? h / withBaf : NaN)) }
}

// --- the allele reading with a parent (R11) ---------------------------------------------------

/** x above this is a marker where the cell carries the allele the loaded parent lacks. */
export const X_MIN = 0.15

/**
 * At markers where the loaded parent is homozygous: x = the cell's B-allele frequency where the
 * parent is AA, one minus it where BB. Returns the mean of x over markers with x above X_MIN, the
 * share of such markers, and the number of markers read.
 */
export function xStat(rows: readonly Row[], keep: (r: Row) => boolean, parentGt: { get(id: string): string | undefined }) {
  let n = 0, k = 0, s = 0
  for (const r of rows) {
    if (!keep(r) || r.baf === null || !Number.isFinite(r.baf)) continue
    const p = parentGt.get(r.probesetId)
    const x = p === 'AA' ? r.baf : p === 'BB' ? 1 - r.baf : NaN
    if (Number.isNaN(x)) continue
    n++
    if (x > X_MIN) { k++; s += x }
  }
  return { mean: k ? s / k : NaN, share: n ? k / n : NaN, n }
}

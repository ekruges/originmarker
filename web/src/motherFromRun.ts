/**
 * The other parent read off the other draws of her a run holds, with ONE parent loaded. Plan item 3.0d,
 * second build; rule and measurements in audit/truth/mother-from-the-run-2026-11-16.txt.
 *
 * A draw of the mother is any array carrying one of her copies: another embryo of hers, her pronucleus,
 * a polar body. At a marker where the loaded parent P is homozygous only she can give the allele P lacks
 * (founder frequency q). Two draws of her carry the same copy of hers over long stretches and different
 * ones elsewhere, so their showing of that allele is correlated (K = 2 where they share a copy, 0 where
 * not); draws of different women share nothing, whatever their ancestry or noise shifts every window by.
 * Each draw then votes on her genotype, weighted by its own drop-in (b) and dropout (a), and the
 * two-parent test reads the sample's share at the markers where she is the opposite homozygote.
 *
 * 1. EMBRYOS. Arrays reading as one genome (sameGenomeShare, lower direction, >= SAME_GENOME_SHARE) are
 *    one unit; only arrays with a rate fit and b < NO_TEMPLATE_BASE are joined. A unit's call is the
 *    majority of its called arrays.
 * 2. A DRAW shows the other parent's alleles (a >= SHOWS_OTHER_MIN) and is clean (b <= DROPIN_CEILING).
 *    The sample's own unit passes both or the run is not used for it.
 * 3. SAME MOTHER: the share of WINDOW_BP windows reading as one shared copy (K >= SHARED_COPY_K) is at or
 *    over the bar for the parent's kind.
 * 4. HER GENOTYPE: per marker, the posterior that she is the opposite homozygote, entered at
 *    >= POSTERIOR_MIN, with rule 2 of oneParentGain.ts for a rebuilt parent.
 * 5. THE TEST: gainOrigin's bars, on a reference of the chromosomes reading both parents.
 */
import { otherRates, sameGenomeShare, SAME_GENOME_SHARE, type OtherRates } from './otherParentLoss.ts'
import { NO_TEMPLATE_BASE, MIN_SHARE_MARKERS } from './oneParentState.ts'
import { DROPIN_CEILING } from './oneParentOrigin.ts'
import { RARE_MIN } from './parentGenome.ts'
import {
  ALLELE_GAIN_MIN_INFORMATIVE, ALLELE_GAIN_MIN_REFERENCE, ALLELE_GAIN_SCALE_FLOOR, ALLELE_GAIN_LRR_SCALE_FLOOR,
  ALLELE_GAIN_Z, ALLELE_GAIN_SHIFT, ALLELE_GAIN_LRR_Z,
} from './gainOrigin.ts'
import { rebuiltHomWeight } from './oneParentGain.ts'
import { mad, median } from './intensityNull.ts'

/** Fitted signal a at or over which a unit shows the other parent's alleles: the middle of the measured
 *  gap between units with no mother's genome (0.00 to 0.05) and units with one (0.32 and up). */
export const SHOWS_OTHER_MIN = 0.19
/** Same-mother bars: 1.9 times the highest window share any two embryos of different mothers reached
 *  (0.143 against a measured parent, 0.052 against a rebuilt one), the placing of 3.0b rule 6. */
export const SAME_MOTHER_MEASURED = 0.272
export const SAME_MOTHER_REBUILT = 0.100
/** The window the two bars were measured with. */
export const WINDOW_BP = 10e6
/** Window score reading as a shared copy: midway between the shared (2) and unshared (0) levels. */
export const SHARED_COPY_K = 1
/** A marker enters where she is more likely than not the opposite homozygote. */
export const POSTERIOR_MIN = 0.5
/** Fewest draws of her before the run is read for her genotype. */
export const MIN_DRAWS = 2
/** The cap on a probability sameGenomeShare uses. */
const P_MAX = 0.999

/** A marker where the loaded parent is homozygous; q the founder frequency of the allele it lacks. */
export interface RunMarker { id: string; chrom: string; pos: number; parent: 'AA' | 'BB'; q: number }
/** Genotype calls of one array, by probe id. */
export interface Calls { get(id: string): string | undefined }

/** 1 the call shows the allele the parent lacks, 0 it is the parent's homozygote, -1 no call. */
export function showsOf(markers: readonly RunMarker[], calls: Calls): Int8Array {
  const s = new Int8Array(markers.length)
  for (let i = 0; i < markers.length; i += 1) {
    const g = calls.get(markers[i].id)
    s[i] = g === 'AA' || g === 'AB' || g === 'BB' ? (g === markers[i].parent ? 0 : 1) : -1
  }
  return s
}

/** The rate fit of otherParentLoss.otherRates on one set of calls; `inside` bounds the q used. */
export function ratesOfShows(markers: readonly RunMarker[], s: Int8Array, open1 = false): OtherRates | null {
  const by = new Map<string, { id: string; pos: number; q: number; seen: boolean; opp: boolean }[]>()
  for (let i = 0; i < markers.length; i += 1) {
    const m = markers[i]
    if (s[i] < 0 || !(m.q > 0) || (open1 && !(m.q < 1))) continue
    let v = by.get(m.chrom)
    if (!v) by.set(m.chrom, v = [])
    v.push({ id: m.id, pos: m.pos, q: m.q, seen: s[i] === 1, opp: false })
  }
  return otherRates([...by.values()])
}

/** Same-genome share of y read at the markers where x shows the allele, y's rates. */
function shareAt(markers: readonly RunMarker[], x: Int8Array, y: Int8Array, ry: OtherRates): number {
  const where: { q: number; seen: boolean }[] = []
  for (let i = 0; i < markers.length; i += 1) {
    if (x[i] === 1 && y[i] >= 0 && markers[i].q > 0) where.push({ q: markers[i].q, seen: y[i] === 1 })
  }
  return sameGenomeShare(where, ry)
}

/** Joined through shared partners: `link(i, j)` for i < j, returns groups of indices. */
function groupsOf(n: number, link: (i: number, j: number) => boolean): number[][] {
  const par = Array.from({ length: n }, (_, i) => i)
  const find = (x: number): number => { while (par[x] !== x) x = par[x] = par[par[x]]; return x }
  for (let i = 0; i < n; i += 1) {
    for (let j = i + 1; j < n; j += 1) {
      if (find(i) !== find(j) && link(i, j)) par[find(i)] = find(j)
    }
  }
  const g = new Map<number, number[]>()
  for (let i = 0; i < n; i += 1) { const r = find(i); if (!g.has(r)) g.set(r, []); g.get(r)!.push(i) }
  return [...g.values()]
}

/** Majority call of a unit's arrays; a tie is no call. */
export function majority(arrays: readonly Int8Array[]): Int8Array {
  const n = arrays[0].length
  const s = new Int8Array(n)
  for (let i = 0; i < n; i += 1) {
    let k = 0; let c = 0
    for (const a of arrays) { if (a[i] < 0) continue; c += 1; k += a[i] }
    s[i] = c && 2 * k !== c ? (2 * k > c ? 1 : 0) : -1
  }
  return s
}

/** Share of windows in which two units read as carrying one copy of the same mother. NaN with none. */
export function sharedCopyShare(
  markers: readonly RunMarker[], x: Int8Array, y: Int8Array, rx: OtherRates, ry: OtherRates,
): number {
  const win = new Map<string, number[]>()
  for (let i = 0; i < markers.length; i += 1) {
    const m = markers[i]
    if (x[i] < 0 || y[i] < 0 || !(m.q >= RARE_MIN) || !(m.q < 1)) continue
    const p = (x[i] - (rx.b + rx.a * m.q)) * (y[i] - (ry.b + ry.a * m.q))
    const key = `${m.chrom}:${Math.floor(m.pos / WINDOW_BP)}`
    let w = win.get(key)
    if (!w) win.set(key, w = [0, 0, 0])
    w[0] += p; w[1] += rx.a * ry.a * m.q * (1 - m.q) / 2; w[2] += 1
  }
  let n = 0; let shared = 0
  for (const w of win.values()) {
    if (w[2] < MIN_SHARE_MARKERS) continue
    n += 1
    if (w[0] / w[1] >= SHARED_COPY_K) shared += 1
  }
  return n ? shared / n : NaN
}

/** A draw's own likelihood of showing the allele: mother opposite homozygote, heterozygous, P's homozygote. */
const showChance = (r: OtherRates): [number, number, number] =>
  [Math.min(P_MAX, r.b + r.a), Math.min(P_MAX, r.b + r.a / 2), Math.min(P_MAX, r.b)]

/** Per marker, the chance she is the opposite homozygote, 0 where under POSTERIOR_MIN or unvoted. */
export function oppositeHomPosterior(
  markers: readonly RunMarker[], draws: readonly { s: Int8Array; r: OtherRates }[], mMin?: number,
): Float64Array {
  const w = new Float64Array(markers.length)
  const p = draws.map((d) => showChance(d.r))
  for (let i = 0; i < markers.length; i += 1) {
    const q = markers[i].q
    if (!(q > 0 && q < 1)) continue
    let lBB = 2 * Math.log(q); let lAB = Math.log(2 * q * (1 - q)); let lAA = 2 * Math.log(1 - q); let k = 0
    for (let d = 0; d < draws.length; d += 1) {
      const v = draws[d].s[i]
      if (v < 0) continue
      k += 1
      const [bb, ab, aa] = p[d]
      if (v === 1) { lBB += Math.log(bb); lAB += Math.log(ab); lAA += Math.log(aa) }
      else { lBB += Math.log(1 - bb); lAB += Math.log(1 - ab); lAA += Math.log(1 - aa) }
    }
    if (!k) continue
    const top = Math.max(lBB, lAB, lAA)
    const post = Math.exp(lBB - top) / (Math.exp(lBB - top) + Math.exp(lAB - top) + Math.exp(lAA - top))
    if (post < POSTERIOR_MIN) continue
    w[i] = post * rebuiltHomWeight(q, mMin)
  }
  return w
}

/** Why the run's other arrays were or were not used for this sample. */
export type RunStatus = 'used' | 'sample-unfit' | 'sample-not-a-draw' | 'too-few-draws'

export interface MotherFromRun {
  status: RunStatus
  /** Draws of the mother that voted. */
  draws: number
  /** Units the run's arrays formed, and how many of them passed rule 2. */
  units: number
  readable: number
  /** Indices into `others` of the arrays whose units voted. */
  drawArrays: number[]
  /** Per chromosome: the weighted median of the loaded parent's allele share, and its summed weight. */
  shares: Map<string, { share: number; weight: number }>
  /** Probe ids where she is entered as the loaded parent's opposite homozygote. For figures; scores nothing. */
  opposite?: Set<string>
}

const unitOf = (markers: readonly RunMarker[], members: readonly Int8Array[], arrays: number[]) => {
  const s = members.length === 1 ? members[0] : majority(members)
  const r = ratesOfShows(markers, s, true)
  return { s, r, ok: !!r && r.a >= SHOWS_OTHER_MIN && r.b <= DROPIN_CEILING, arrays }
}
const genomeShare = (markers: readonly RunMarker[], x: Int8Array, rx: OtherRates, y: Int8Array, ry: OtherRates) =>
  Math.min(shareAt(markers, x, y, ry), shareAt(markers, y, x, rx))

/** What does not depend on the sample, kept per loaded parent across calls: each array's calls and fit,
 *  the same-genome share of each pair, each unit, and the window share of each pair of units. A sample
 *  then joins the run's arrays, which gives the same embryos and groups as reading all of them at once. */
interface Unit { s: Int8Array; r: OtherRates | null; ok: boolean; arrays: number[] }
interface ParentCache {
  markers: RunMarker[]
  arrays: WeakMap<Calls, { id: number; s: Int8Array; r: OtherRates | null }>
  genome: Map<number, number>
  units: Map<string, Unit>
  shared: Map<string, number>
  next: number
}
const caches = new WeakMap<object, ParentCache>()

/**
 * The mother read off the run, for one sample. `markers` are the loaded parent's homozygous markers with
 * a founder frequency; `sample` the sample's calls; `others` every other array of the run.
 */
export function motherFromRun(input: {
  markers: RunMarker[]; sample: Calls; baf: (id: string) => number | undefined
  others: readonly Calls[]; parentKey: object; rebuilt: boolean; mMin?: number
}): MotherFromRun {
  const { sample, others, rebuilt, mMin } = input
  let pc = caches.get(input.parentKey)
  if (!pc) caches.set(input.parentKey, pc = { markers: input.markers, arrays: new WeakMap(), genome: new Map(), units: new Map(), shared: new Map(), next: 0 })
  const c = pc
  const { markers } = c
  const none = (status: RunStatus, units = 0, readable = 0): MotherFromRun =>
    ({ status, draws: 0, units, readable, drawArrays: [], shares: new Map() })
  const own = showsOf(markers, sample)
  const ownRates = ratesOfShows(markers, own)
  if (!ownRates || !(ownRates.b < NO_TEMPLATE_BASE)) return none('sample-unfit')
  const arrs = others.map((o) => {
    let a = c.arrays.get(o)
    if (!a) { const sh = showsOf(markers, o); c.arrays.set(o, a = { id: c.next++, s: sh, r: ratesOfShows(markers, sh) }) }
    return a
  })
  // 1. Embryos: the sample (index 0) and the usable arrays, joined at the same-genome bar.
  const usable = arrs.map((a, i) => ({ a, i })).filter(({ a }) => a.r && a.r.b < NO_TEMPLATE_BASE)
  const genome = (i: number, j: number): number => {
    if (i === 0) return genomeShare(markers, own, ownRates, usable[j - 1].a.s, usable[j - 1].a.r!)
    const x = usable[i - 1].a; const y = usable[j - 1].a
    const key = Math.min(x.id, y.id) * 1048576 + Math.max(x.id, y.id)
    let v = c.genome.get(key)
    if (v === undefined) c.genome.set(key, v = genomeShare(markers, x.s, x.r!, y.s, y.r!))
    return v
  }
  const embryos = groupsOf(usable.length + 1, (i, j) => genome(i, j) >= SAME_GENOME_SHARE)
  // 2. Units: the sample's is read fresh, the others are kept by their arrays.
  const units = embryos.map((e) => {
    if (e.includes(0)) {
      const rest = e.filter((i) => i > 0)
      return { key: '', u: unitOf(markers, [own, ...rest.map((i) => usable[i - 1].a.s)], rest.map((i) => usable[i - 1].i)) }
    }
    const key = e.map((i) => usable[i - 1].a.id).sort((x, y) => x - y).join(',')
    let u = c.units.get(key)
    if (!u) c.units.set(key, u = unitOf(markers, e.map((i) => usable[i - 1].a.s), []))
    return { key, u: { ...u, arrays: e.map((i) => usable[i - 1].i) } }
  })
  const mine = units.find((x) => x.key === '')!
  const readable = units.filter((x) => x.key !== '' && x.u.ok)
  if (!mine.u.ok) return none('sample-not-a-draw', units.length, readable.length)
  // 3. Same mother.
  const bar = rebuilt ? SAME_MOTHER_REBUILT : SAME_MOTHER_MEASURED
  const all = [mine, ...readable]
  const groups = groupsOf(all.length, (i, j) => {
    const x = all[i]; const y = all[j]
    if (i === 0) return sharedCopyShare(markers, x.u.s, y.u.s, x.u.r!, y.u.r!) >= bar
    const key = `${x.key}|${y.key}`
    let v = c.shared.get(key)
    if (v === undefined) c.shared.set(key, v = sharedCopyShare(markers, x.u.s, y.u.s, x.u.r!, y.u.r!))
    return v >= bar
  })
  const voters = groups.find((g) => g.includes(0))!.filter((i) => i !== 0)
  const draws = voters.map((i) => ({ s: all[i].u.s, r: all[i].u.r! }))
  const drawArrays = voters.flatMap((i) => all[i].u.arrays)
  if (draws.length < MIN_DRAWS) {
    return { ...none('too-few-draws', units.length, readable.length + 1), draws: draws.length, drawArrays }
  }
  // 4. Her genotype, and the sample's share where she is the opposite homozygote.
  const w = oppositeHomPosterior(markers, draws, rebuilt ? mMin : undefined)
  const by = new Map<string, [number, number][]>()
  for (let i = 0; i < markers.length; i += 1) {
    if (!(w[i] > 0)) continue
    const b = input.baf(markers[i].id)
    if (b === undefined || !Number.isFinite(b)) continue
    let v = by.get(markers[i].chrom)
    if (!v) by.set(markers[i].chrom, v = [])
    v.push([markers[i].parent === 'AA' ? 1 - b : b, w[i]])
  }
  const shares = new Map<string, { share: number; weight: number }>()
  for (const [c, v] of by) {
    v.sort((x, y) => x[0] - y[0])
    const tot = v.reduce((s, x) => s + x[1], 0)
    let acc = 0; let med = NaN
    for (const [s, x] of v) { acc += x; if (acc >= tot / 2) { med = s; break } }
    shares.set(c, { share: med, weight: tot })
  }
  const opposite = new Set<string>()
  for (let i = 0; i < markers.length; i += 1) if (w[i] > 0) opposite.add(markers[i].id)
  return { status: 'used', draws: draws.length, units: units.length, readable: readable.length + 1, drawArrays, shares, opposite }
}

/** One chromosome as the route reads it. `share` is the LOADED parent's allele share. */
export interface RouteChrom {
  chrom: string; share: number; weight: number; state?: string; measured: boolean
  intensityGain: boolean; lrrShift: number
}
export interface RouteCall {
  chrom: string; z: number; shift: number; scale: number; centre: number; lrrZ: number; weight: number
  references: number
  /** True where the extra copy is the loaded parent's. */
  loadedParent: boolean
}

const lostState = (s?: string) => s === 'other-lost' || s === 'parent-lost'

/**
 * The two-parent test on the share read against her. Returns the chromosomes it names; `skip` holds the
 * ones the first build already named a parent on. Nothing where the cell's allele gains are withheld for
 * want of a two-copy majority (`withheld`).
 */
export function routeGains(chroms: readonly RouteChrom[], skip: ReadonlySet<string>, withheld = false): RouteCall[] {
  if (withheld) return []
  const disomic = (c: RouteChrom) => c.state === 'disomic' && c.measured
  const levels = chroms.filter((c) => disomic(c) && Number.isFinite(c.lrrShift)).map((c) => c.lrrShift)
  const level = levels.length ? median(levels) : Infinity
  // A lost copy cannot raise intensity: a lost state at or above the cell's level reads both parents.
  const both = (c: RouteChrom) => disomic(c)
    || (lostState(c.state) && c.measured && Number.isFinite(c.lrrShift) && c.lrrShift >= level)
  const out: RouteCall[] = []
  for (const t of chroms) {
    if (skip.has(t.chrom) || !t.measured || t.state === 'no-template' || t.state === 'undecided') continue
    if (lostState(t.state) && !t.intensityGain) continue
    if (!(t.weight >= ALLELE_GAIN_MIN_INFORMATIVE) || !Number.isFinite(t.share)) continue
    // The share reference needs ALLELE_GAIN_MIN_INFORMATIVE of summed weight (weight, not a count).
    const ref = chroms.filter((c) => c.chrom !== t.chrom && both(c) && c.weight >= ALLELE_GAIN_MIN_INFORMATIVE
      && Number.isFinite(c.share))
    if (ref.length < ALLELE_GAIN_MIN_REFERENCE) continue
    const shares = ref.map((c) => c.share)
    const centre = median(shares)
    const scale = Math.max(mad(shares), ALLELE_GAIN_SCALE_FLOOR)
    const shift = t.share - centre
    const z = shift / scale
    // The intensity reference needs no allele markers.
    const lrrs = chroms.filter((c) => c.chrom !== t.chrom && both(c) && Number.isFinite(c.lrrShift)).map((c) => c.lrrShift)
    if (lrrs.length < ALLELE_GAIN_MIN_REFERENCE || !Number.isFinite(t.lrrShift)) continue
    const lrrZ = (t.lrrShift - median(lrrs)) / Math.max(mad(lrrs), ALLELE_GAIN_LRR_SCALE_FLOOR)
    if (Math.abs(z) >= ALLELE_GAIN_Z && Math.abs(shift) >= ALLELE_GAIN_SHIFT && lrrZ >= ALLELE_GAIN_LRR_Z) {
      out.push({ chrom: t.chrom, z, shift, scale, centre, lrrZ, weight: t.weight, references: ref.length,
        loadedParent: shift > 0 })
    }
  }
  return out
}

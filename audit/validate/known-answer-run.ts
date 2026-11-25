/**
 * The known-answer list, scored by the shipped scoreSample: each base cell is read and parsed once,
 * scored untouched, then scored again with each of its built events applied in memory
 * (known-answer.ts). Nothing is written per event but its record.
 *
 * The job list is JSON lines, one base cell each:
 *   { key, path, parent (the loaded parent's array, '-' for none), other (the second parent's array,
 *     optional), role ('maternal' | 'paternal', the loaded parent's), exclude: [chromosomes a manual
 *     or published call names], measure: [{ name, chrom, start, end }], events: [EventSpec] }
 * `measure` names real events to read the same way as the built ones; an end of 0 is the
 * chromosome's end.
 *
 * TWO READINGS PER CELL, each in the cell's own units, beside the shipped tool. For a stretch:
 * intensity, the mean of log2R clipped to +-2; and with a parent x, the share of the cell's signal on
 * the allele the loaded parent lacks, at markers where it is homozygous (with both parents, only
 * where they are opposite homozygotes), as its mean over x above X_MIN and the share of such
 * markers. Each is a z against the cell's stretches of the same span class (5, 10, 20, 40, 80 or
 * 160 Mb, the nearest on a log scale; step a quarter span) on its reference chromosomes other than
 * the one read: (value - median) / (1.4826 MAD), at least 8 stretches. A stretch counts only if it
 * holds at least half the median marker count of the cell's stretches of that span.
 *
 * Records, JSON lines:
 *  - `base`: the untouched result (every printed change, each autosome's whole-chromosome reading,
 *    copyReference, the integrity level and the gates that withheld anything), the cell's reference
 *    chromosomes (autosomes neither excluded nor carrying a printed change), its noise, each
 *    autosome's marker extent, the signature and the two readings of each `measure` region, the
 *    signature of the reference, the time spent reading, unpacking, parsing and scoring, and the
 *    readings' null: [zI, zMean, zShare] of every reference stretch, per span class.
 *  - `event`: the result with the event applied (every printed change, withheld claims and the
 *    whole-chromosome reading on the event's chromosome), and the signature and the two readings of
 *    the built interval. A reading carries the reference stretches' medians beside its z.
 *
 * Noise is 1.4826 x MAD of the medians of 10 Mb windows (at least 100 markers) on the reference,
 * about their median.
 *
 * With OM_SIGONLY=1 nothing is scored: the records carry the signatures alone, the reference being
 * the autosomes the job does not exclude. That is how the builder is compared with real events.
 *
 * Run: node --experimental-strip-types audit/validate/known-answer-run.ts <jobs.jsonl> <shard> <shards> <out.jsonl>
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { performance } from 'node:perf_hooks'
import { applyEvent, signature, levelTable, X_MIN, type Row, type EventSpec } from './known-answer.ts'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const keys = await import(`${W}answerKeys.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)
const platform = await import(`${W}platform.ts`)
const pageSays = await import(`${W}pageSays.ts`)

const [jobFile, shard, shards, out] = [process.argv[2], Number(process.argv[3]), Number(process.argv[4]), process.argv[5]]
const SIG_ONLY = process.env.OM_SIGONLY === '1'
const AUTOSOMES = Array.from({ length: 22 }, (_, i) => String(i + 1))
/** Stretch spans for the allele reading's null, Mb; an event uses the nearest on a log scale. */
const X_SPANS = [5, 10, 20, 40, 80, 160]
const done = new Set<string>()
if (existsSync(out)) {
  for (const l of readFileSync(out, 'utf8').split('\n')) {
    if (l.startsWith('{')) { const r = JSON.parse(l); done.add(`${r.key}|${r.id ?? ''}`) }
  }
}
const jobs = readFileSync(jobFile, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l))
  .filter((_: any, i: number) => i % shards === shard)

const t = { read: 0, gunzip: 0, parse: 0 }
function text(p: string) {
  let a = performance.now(); const raw = readFileSync(p); t.read += performance.now() - a
  a = performance.now(); const s = (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8'); t.gunzip += performance.now() - a
  return s
}
/** The rows of an array as run.ts reads them, with the same chip recognition. */
function rows(path: string): Row[] {
  const s = text(path)
  const a = performance.now()
  const lines = s.split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = ingest.headerMap(lines[h])
  if (!map) throw new Error(`no header map: ${path.split('/').pop()}`)
  const o: Row[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) o.push(r) }
  const b = ingest.emptyBuildSums()
  for (const r of o) ingest.accumulateBuild(r, b)
  const call = platform.recogniseRows(o, ingest.buildVerdict(b), map, panel.freqB)
  if (!call.proven) throw new Error(`chip not proven: ${path.split('/').pop()}`)
  t.parse += performance.now() - a
  return o
}
const panel = keys.parsePanel(text(new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
const parents = new Map<string, any>()
function parentIndex(path: string) {
  if (path === '-') return score.noParent()
  if (!parents.has(path)) {
    const banner = inferredArray.parseInferredBanner(text(path).slice(0, inferredArray.MARK_BYTES)) ?? undefined
    const acc = score.emptyParent()
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rows(path)) {
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectParentRow(r, acc)
    }
    parents.set(path, score.finishParent(acc, ingest.finishProfile(path, byChrom, bafSums, first).build.build, banner))
  }
  return parents.get(path)
}

/** One array's rows scored as run.ts scores a row with no other arrays of its run (`mat` the other parent). */
async function scoreRows(rs: readonly Row[], pat: any, mat: any, role: string, name: string) {
  const acc = score.emptyCollected(pat, mat, panel as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); const builds = ingest.emptyBuildSums(); let first = ''
  for (const r of rs) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectRow(r, pat, mat, acc as never)
    ingest.accumulateBuild(r, builds)
  }
  const profile = ingest.finishProfile(name, byChrom, bafSums, first, builds)
  return await score.scoreSample({ acc, profile, pat, mat, soloRole: role === 'maternal' ? 'maternal' : 'paternal',
    sibs: [], sampleName: name, log: () => {} } as never) as any
}
const changesOf = (res: any) => pageSays.defectsForResult(res).map((d: any) => ({ chrom: d.chrom, start: d.startBp,
  end: d.endBp, kind: d.kind, origin: d.origin, step: d.step ?? null }))
const chromOf = (res: any, c: string) => {
  const x = (res.chroms ?? []).find((y: any) => y.chrom === c)
  return x ? { aneuploidy: x.aneuploidy ?? '', lostParent: x.lostParent ?? '', lrrShift: x.lrrShift ?? null,
    lrrZ: x.lrrZ ?? null, gain: x.alleleGain?.verdict ?? '', gainParent: x.alleleGain?.parent ?? '' } : null
}
const r4 = (x: number) => (Number.isFinite(x) ? Math.round(x * 1e4) / 1e4 : null)
const sig = (s: ReturnType<typeof signature>) => ({ markers: s.markers, lrr: r4(s.lrr), abShare: r4(s.abShare),
  noCall: r4(s.noCall), offPole: r4(s.offPole), bafHist: s.bafHist.map(r4) })
const median = (xs: number[]) => { const q = xs.filter(Number.isFinite).sort((a, b) => a - b); return q.length ? q[q.length >> 1] : NaN }
const mad = (xs: number[], m: number) => 1.4826 * median(xs.map((x) => Math.abs(x - m)))

function noise(rs: readonly Row[], ref: Set<string>) {
  const win = new Map<string, number[]>()
  for (const r of rs) {
    if (!ref.has(r.chrom) || r.log2R === null || !Number.isFinite(r.log2R)) continue
    const k = `${r.chrom}:${Math.floor(r.pos / 10e6)}`
    let w = win.get(k); if (!w) win.set(k, w = []); w.push(r.log2R)
  }
  const meds = [...win.values()].filter((w) => w.length >= 100).map(median)
  return mad(meds, median(meds))
}

/**
 * Per chromosome, sorted by position, prefix sums for the two per-cell readings: intensity (log2R
 * clipped to +-LRR_CLIP) at every marker with one, and x at the informative markers (see `xp` below).
 */
const LRR_CLIP = 2
type Track = { pos: Float64Array; n: Float64Array; s: Float64Array; xpos: Float64Array; xn: Float64Array; xk: Float64Array; xs: Float64Array }
function tracksOf(rs: readonly Row[], xp: Map<string, string> | null, only?: string) {
  const by = new Map<string, { l: [number, number][]; x: [number, number][] }>()
  for (const r of rs) {
    if (only ? r.chrom !== only : !AUTOSOMES.includes(r.chrom)) continue
    let b = by.get(r.chrom); if (!b) by.set(r.chrom, b = { l: [], x: [] })
    if (r.log2R !== null && Number.isFinite(r.log2R)) b.l.push([r.pos, Math.max(-LRR_CLIP, Math.min(LRR_CLIP, r.log2R))])
    if (xp && r.baf !== null && Number.isFinite(r.baf)) {
      const p = xp.get(r.probesetId)
      if (p === 'AA') b.x.push([r.pos, r.baf]); else if (p === 'BB') b.x.push([r.pos, 1 - r.baf])
    }
  }
  const out = new Map<string, Track>()
  for (const [c, b] of by) {
    b.l.sort((u, v) => u[0] - v[0]); b.x.sort((u, v) => u[0] - v[0])
    const n = new Float64Array(b.l.length + 1), s = new Float64Array(b.l.length + 1)
    b.l.forEach(([, v], j) => { n[j + 1] = n[j] + 1; s[j + 1] = s[j] + v })
    const xn = new Float64Array(b.x.length + 1), xk = new Float64Array(b.x.length + 1), xs = new Float64Array(b.x.length + 1)
    b.x.forEach(([, v], j) => { xn[j + 1] = xn[j] + 1; xk[j + 1] = xk[j] + (v > X_MIN ? 1 : 0); xs[j + 1] = xs[j] + (v > X_MIN ? v : 0) })
    out.set(c, { pos: Float64Array.from(b.l.map((u) => u[0])), n, s, xpos: Float64Array.from(b.x.map((u) => u[0])), xn, xk, xs })
  }
  return out
}
const lowerBound = (a: Float64Array, v: number) => { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < v) lo = m + 1; else hi = m } return lo }
interface Read { n: number; lrr: number; xn: number; mean: number; share: number }
function readAt(t: Track, a: number, b: number): Read {
  const i = lowerBound(t.pos, a), j = lowerBound(t.pos, b)
  const xi = lowerBound(t.xpos, a), xj = lowerBound(t.xpos, b)
  const n = t.n[j] - t.n[i], xn = t.xn[xj] - t.xn[xi], k = t.xk[xj] - t.xk[xi]
  return { n, lrr: n ? (t.s[j] - t.s[i]) / n : NaN, xn, mean: k >= 10 ? (t.xs[xj] - t.xs[xi]) / k : NaN, share: xn >= 50 ? k / xn : NaN }
}
/** Every stretch of `span` Mb on a chromosome, step span/4. */
function stretches(t: Track | undefined, span: number): Read[] {
  const out: Read[] = []
  if (!t || !t.pos.length) return out
  const L = span * 1e6, first = t.pos[0], last = t.pos[t.pos.length - 1]
  for (let a = first; a + L <= last + 1; a += L / 4) out.push(readAt(t, a, a + L))
  return out
}
const spanClass = (mb: number) => X_SPANS.reduce((b, s) => (Math.abs(Math.log(s / mb)) < Math.abs(Math.log(b / mb)) ? s : b))

for (const job of jobs) {
  const t0 = { ...t }
  let a = performance.now()
  const pat = parentIndex(job.parent)
  const mat = job.other ? parentIndex(job.other) : null
  const base = rows(job.path)
  for (const r of base) Object.freeze(r)
  Object.freeze(base)
  const pgt: Map<string, string> | null = job.parent === '-' ? null : pat.gt
  const ogt: Map<string, string> | null = mat ? mat.gt : null
  // x is read where the loaded parent is homozygous; with both parents, only where they are
  // opposite homozygotes.
  let xp: Map<string, string> | null = null
  if (pgt) {
    xp = new Map()
    for (const [id, g] of pgt) {
      if (g !== 'AA' && g !== 'BB') continue
      if (ogt && ogt.get(id) !== (g === 'AA' ? 'BB' : 'AA')) continue
      xp.set(id, g)
    }
  }
  const pathTimes = { read: t.read - t0.read, gunzip: t.gunzip - t0.gunzip, parse: t.parse - t0.parse }
  a = performance.now()
  const res0 = SIG_ONLY ? { chroms: [], refusals: [] } : await scoreRows(base, pat, mat, job.role, job.key)
  const scoreMs = performance.now() - a
  const changes0 = SIG_ONLY ? [] : changesOf(res0)
  const changed = new Set(changes0.map((c: any) => String(c.chrom)))
  const ref = new Set(AUTOSOMES.filter((c) => !job.exclude.includes(c) && !changed.has(c)))
  const extent: Record<string, [number, number, number]> = {}
  for (const r of base) {
    if (!AUTOSOMES.includes(r.chrom)) continue
    const e = extent[r.chrom] ?? (extent[r.chrom] = [Infinity, -Infinity, 0])
    if (r.pos < e[0]) e[0] = r.pos
    if (r.pos > e[1]) e[1] = r.pos
    e[2]++
  }
  // Homozygous B-allele frequencies on the reference, for a built loss to draw from.
  const AA: number[] = [], BB: number[] = []
  for (const r of base) {
    if (!ref.has(r.chrom) || r.baf === null || !Number.isFinite(r.baf)) continue
    if (r.genotype === 'AA') AA.push(r.baf); else if (r.genotype === 'BB') BB.push(r.baf)
  }
  // R25: a job may name the parents whose genotypes BUILD a gain's allele side (buildParent, buildOther) while
  // the cell is scored with the parent it names in \`parent\` (none for R25).
  const bpg = job.buildParent ? parentIndex(job.buildParent).gt : pgt
  const bog = job.buildOther ? parentIndex(job.buildOther).gt : job.buildParent ? null : ogt
  const ctx = { parentGt: bpg, otherGt: bog, homBaf: { AA: Float64Array.from(AA), BB: Float64Array.from(BB) },
    levels: levelTable(base, (r) => ref.has(r.chrom)) }
  const refSig = signature(base, (r) => ref.has(r.chrom))
  const tracks = tracksOf(base, xp)
  // The reference stretches of a span on the reference chromosomes other than `not`. A stretch is
  // read only where it holds at least half the median marker count of the cell's stretches of that
  // span, so tiles over centromeres and short arms do not set the spread.
  const byChrom = new Map<number, Map<string, Read[]>>()
  const stretchesBy = (span: number) => {
    if (!byChrom.has(span)) {
      const m = new Map([...ref].map((c) => [c, stretches(tracks.get(c), span)] as [string, Read[]]))
      const nMin = median([...m.values()].flat().map((r) => r.n)) / 2
      const xMin = median([...m.values()].flat().map((r) => r.xn)) / 2
      for (const [c, v] of m) m.set(c, v.filter((r) => r.n >= nMin && (!xp || r.xn >= xMin)))
      byChrom.set(span, m)
    }
    return byChrom.get(span)!
  }
  const nullCache = new Map<string, { lrr: number[]; mean: number[]; share: number[] }>()
  const refOf = (span: number, not: string) => {
    const k = `${span}|${not}`
    if (!nullCache.has(k)) {
      const all = [...stretchesBy(span)].filter(([c]) => c !== not).flatMap(([, v]) => v)
      nullCache.set(k, { lrr: all.map((x) => x.lrr), mean: all.map((x) => x.mean).filter(Number.isFinite),
        share: all.map((x) => x.share).filter(Number.isFinite) })
    }
    return nullCache.get(k)!
  }
  /** A stretch in the cell's own units: z against the reference stretches of its span class. */
  const zOf = (v: Read, span: number, not: string) => {
    const s = refOf(span, not)
    const z = (x: number, xs: number[]) => { if (xs.length < 8 || !Number.isFinite(x)) return null; const m = median(xs); return r4((x - m) / mad(xs, m)) }
    return { zI: z(v.lrr, s.lrr), zMean: xp ? z(v.mean, s.mean) : null, zShare: xp ? z(v.share, s.share) : null,
      refLrr: r4(median(s.lrr)), refMean: xp ? r4(median(s.mean)) : null, refShare: xp ? r4(median(s.share)) : null }
  }
  const readOut = (v: Read, span: number, not: string) => ({ n: v.n, lrr: r4(v.lrr), xn: v.xn, mean: r4(v.mean), share: r4(v.share),
    span, ...zOf(v, span, not) })
  if (!done.has(`${job.key}|`)) {
    // The per-cell readings' null: every reference stretch against the other reference chromosomes.
    const nul: Record<string, (number | null)[][]> = {}
    for (const span of X_SPANS) {
      nul[span] = []
      for (const [c, v] of stretchesBy(span)) for (const r of v) { const z = zOf(r, span, c); nul[span].push([z.zI, z.zMean, z.zShare]) }
    }
    const measures = (job.measure ?? []).map((m: any) => {
      const e = extent[m.chrom]
      const end = m.end || (e ? e[1] + 1 : 0)
      const keep = (r: Row) => r.chrom === m.chrom && r.pos >= m.start && r.pos < end
      const tr = tracks.get(m.chrom)
      return { ...m, end, sig: sig(signature(base, keep)),
        read: tr ? readOut(readAt(tr, m.start, end), spanClass((end - m.start) / 1e6), m.chrom) : null }
    })
    appendFileSync(out, JSON.stringify({ kind: 'base', key: job.key, parents: (pgt ? 1 : 0) + (ogt ? 1 : 0), role: job.role,
      ms: { ...Object.fromEntries(Object.entries(pathTimes).map(([k, v]) => [k, Math.round(v)])), score: Math.round(scoreMs) },
      stage: res0.stage?.stage ?? null, ploidy: res0.stage?.ploidy ?? null, integrity: res0.integrity?.level ?? null,
      gates: [...new Set((res0.refusals ?? []).map((y: any) => y.gate))],
      ref: [...ref], noise: r4(noise(base, ref)), extent, refSig: sig(refSig), measures,
      changes: changes0, chroms: Object.fromEntries(AUTOSOMES.map((c) => [c, chromOf(res0, c)])),
      copyReference: res0.copyReference ?? null, null: nul }) + '\n')
  }
  for (const e of (job.events ?? []) as EventSpec[]) {
    if (done.has(`${job.key}|${e.id}`)) continue
    a = performance.now()
    const built = applyEvent(base, e, ctx)
    const buildMs = performance.now() - a
    a = performance.now()
    const res = SIG_ONLY ? { chroms: [], refusals: [] } : await scoreRows(built, pat, mat, job.role, job.key)
    const ms = performance.now() - a
    const keep = (r: Row) => r.chrom === e.chrom && r.pos >= e.start && r.pos < e.end
    const tr = tracksOf(built, xp, e.chrom).get(e.chrom)
    appendFileSync(out, JSON.stringify({ kind: 'event', key: job.key, id: e.id, ms: Math.round(ms), buildMs: Math.round(buildMs),
      changes: SIG_ONLY ? [] : changesOf(res), chrom: chromOf(res, e.chrom),
      withheld: (res.refusals ?? []).filter((y: any) => y.chrom === e.chrom).map((y: any) => ({ claim: y.claim, gate: y.gate,
        cause: y.cause, kind: y.withheld?.kind ?? (typeof y.withheld === 'string' ? y.withheld : null),
        start: y.withheld?.startBp ?? null, end: y.withheld?.endBp ?? null })),
      sig: sig(signature(built, keep)),
      read: tr ? readOut(readAt(tr, e.start, e.end), spanClass((e.end - e.start) / 1e6), e.chrom) : null,
    }) + '\n')
  }
}
appendFileSync(out, `DONE ${shard}\n`)

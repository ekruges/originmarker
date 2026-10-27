/**
 * The het-only rebuild with a TWO-SIDED window, which the failure of the one-sided version
 * implies.
 *
 * One-sided, asserting heterozygosity wherever at least k embryos read heterozygous, reached a
 * precision of 0.6271 against a modelled 0.9920. The cause was measured: drop-in is a property of
 * the MARKER and not of the embryo, so a bad probe collects a vote from every sample and the count
 * at high votes flattens instead of dying, by a factor of 1.2e9 at fourteen of thirty.
 *
 * But a truly heterozygous mother transmits her allele to about HALF her children, so her vote is
 * mid-range by Mendel and never near-unanimous. So near-unanimity is evidence of an artefact rather
 * than of heterozygosity, and the upper bound follows from the binomial rather than from a search:
 * at thirty samples and this dropout the vote is Binomial(30, 0.428), and a bound of 21 discards
 * 0.07 percent of true heterozygous markers while cutting the flat tail.
 *
 * Both bounds are derived. k from the contamination target, u from the transmission distribution.
 *
 * Run: OM_TRIOS=<dir> OM_MOTHER=<gsm> node --experimental-strip-types audit/rebuild-two-sided.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const segs = await import(`${W}segments.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('contested-flag-tests: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

const fileOf = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const m = /^(GSM\d+)/.exec(f)
  if (m && (f.endsWith('.gz') || f.endsWith('.probes') || f.endsWith('.txt'))) fileOf.set(m[1], join(DIR, f))
}
function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r }
  }
}
const cache = new Map<string, unknown>()
function parentIndex(gsm: string) {
  const hit = cache.get(gsm); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  cache.set(gsm, idx)
  return idx
}
async function scoreWith(gsm: string, father: string, mother: string) {
  const pat = parentIndex(father); const mat = parentIndex(mother)
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
  } as never) as never
  return { res, profile } as { res: Res; profile: { callRate: number } }
}










const op = await import(`${W}oneParentOrigin.ts`)
const DROP_IN = op.DROP_IN as number
const H = 0.168  // the egg donor's heterozygosity on this platform

const MOTHER = process.env.OM_MOTHER
if (!MOTHER) { console.log('OM_MOTHER not set. NOT RUN.'); process.exit(0) }

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
let father = ''
const kids: string[] = []
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const m = at(r, 'mother_gsms')
  if (m !== MOTHER) continue
  const g = at(r, 'gsm')
  if (!fileOf.has(g)) continue
  father = father || at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  kids.push(g)
}
if (!father || kids.length < 2) { console.log('family not found. NOT RUN.'); process.exit(0) }

// The embryo each sample came from, read off GEO's own filename: the letter before the cell digit.
const embryoOf = new Map<string, string>()
for (const g of kids) {
  const f = (fileOf.get(g) ?? '').split('/').pop() ?? ''
  const m = /^GSM\d+_([A-Za-z]+)/.exec(f)
  embryoOf.set(g, m ? m[1].toUpperCase() : g)
}
const embryos = [...new Set(embryoOf.values())].sort()
console.log(`family: ${kids.length} samples over ${embryos.length} embryos, one father, `
  + 'the mother arrayed\n')

const gtOf = (g: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(fileOf.get(g)!)()) {
    if (!/^\d+$/.test(r.chrom)) continue
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const faGt = gtOf(father)
const moGt = gtOf(MOTHER)
console.log(`father ${faGt.size.toLocaleString()} called, mother ${moGt.size.toLocaleString()} called`)

// Every sample's genotype at the father-homozygous markers, plus its dropout.
const sample = new Map<string, { gt: Map<string, string>; ado: number }>()
for (const g of kids) {
  const gt = gtOf(g)
  let het = 0
  for (const v of gt.values()) if (v === 'AB') het += 1
  sample.set(g, { gt, ado: op.inferDropout(gt.size ? het / gt.size : NaN) as number })
  process.stdout.write('.')
}
console.log('\n')

// k derived from the contamination target the shipped rebuilt mother carries on its banner.
const TARGET = 0.02749
const lchoose = (n: number, k: number) => {
  let s = 0
  for (let i = 1; i <= k; i += 1) s += Math.log(n - k + i) - Math.log(i)
  return s
}
const tailGE = (n: number, p: number, k: number) => {
  let s = 0
  for (let i = k; i <= n; i += 1) s += Math.exp(lchoose(n, i) + i * Math.log(p) + (n - i) * Math.log(1 - p))
  return Math.min(1, Math.max(0, s))
}
function deriveK(n: number, ado: number) {
  const pSig = 0.5 * (1 - ado)
  for (let k = 1; k <= n; k += 1) {
    const d = tailGE(n, DROP_IN, k)
    const s = tailGE(n, pSig, k)
    const contam = (1 - H) * d / ((1 - H) * d + H * s)
    if (contam <= TARGET) return { k, contam, complete: s }
  }
  return { k: n, contam: NaN, complete: NaN }
}

/** The rebuild: assert the mother heterozygous where at least k of the given samples read AB. */
function rebuild(from: string[]) {
  const ados = from.map((g) => sample.get(g)!.ado).filter((x) => Number.isFinite(x))
  const ado = ados.length ? ados.reduce((a, b) => a + b, 0) / ados.length : 0.3
  const { k, contam, complete } = deriveK(from.length, ado)
  const votes = new Map<string, number>()
  for (const g of from) {
    const gt = sample.get(g)!.gt
    for (const [probe, fg] of faGt) {
      if (fg === 'AB') continue
      if (gt.get(probe) === 'AB') votes.set(probe, (votes.get(probe) ?? 0) + 1)
    }
  }
  const asserted = new Set<string>()
  for (const [probe, v] of votes) if (v >= k) asserted.add(probe)
  return { k, contam, complete, asserted, ado }
}

/** Score an assertion set against the mother's own array. */
function grade(asserted: Set<string>) {
  let hit = 0; let miss = 0
  for (const probe of asserted) {
    const t = moGt.get(probe)
    if (!t) continue
    if (t === 'AB') hit += 1; else miss += 1
  }
  // Completeness: of her true heterozygous markers at father-homozygous sites, how many found.
  let trueHet = 0; let found = 0
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    if (moGt.get(probe) !== 'AB') continue
    trueHet += 1
    if (asserted.has(probe)) found += 1
  }
  return { precision: hit + miss ? hit / (hit + miss) : NaN, scored: hit + miss,
    completeness: trueHet ? found / trueHet : NaN, trueHet }
}


/** The upper bound: the vote a truly heterozygous mother can reach, from the binomial alone. */
function deriveU(n: number, ado: number, discard = 0.001) {
  const p = 0.5 * (1 - ado)
  const lch = (N: number, k: number) => {
    let s = 0
    for (let i = 1; i <= k; i += 1) s += Math.log(N - k + i) - Math.log(i)
    return s
  }
  const pmf = (k: number) => Math.exp(lch(n, k) + k * Math.log(p) + (n - k) * Math.log(1 - p))
  for (let u = n; u >= 1; u -= 1) {
    let tail = 0
    for (let k = u + 1; k <= n; k += 1) tail += pmf(k)
    if (tail > discard) return { u: u + 1, discarded: tail }
  }
  return { u: n, discarded: 0 }
}

function rebuild2(from: string[], twoSided: boolean) {
  const ados = from.map((g) => sample.get(g)!.ado).filter((x) => Number.isFinite(x))
  const ado = ados.length ? ados.reduce((a, b) => a + b, 0) / ados.length : 0.3
  const { k, contam, complete } = deriveK(from.length, ado)
  const { u, discarded } = deriveU(from.length, ado)
  const votes = new Map<string, number>()
  for (const g of from) {
    const gt = sample.get(g)!.gt
    for (const [probe, fg] of faGt) {
      if (fg === 'AB') continue
      if (gt.get(probe) === 'AB') votes.set(probe, (votes.get(probe) ?? 0) + 1)
    }
  }
  const asserted = new Set<string>()
  for (const [probe, v] of votes) {
    if (v < k) continue
    if (twoSided && v > u) continue
    asserted.add(probe)
  }
  return { k, u, contam, complete, discarded, asserted, ado }
}

console.log('THE TWO-SIDED WINDOW AGAINST THE ONE-SIDED RULE, all samples in')
console.log(`${'rule'.padEnd(14)} ${'k'.padStart(3)} ${'u'.padStart(4)} ${'precision'.padStart(10)} `
  + `${'completeness'.padStart(13)} ${'asserted'.padStart(10)}`)
for (const [label, twoSided] of [['one-sided', false], ['two-sided', true]] as [string, boolean][]) {
  const r = rebuild2(kids, twoSided);
  const s = grade(r.asserted)
  console.log(`${label.padEnd(14)} ${String(r.k).padStart(3)} `
    + `${(twoSided ? String(r.u) : '-').padStart(4)} ${s.precision.toFixed(4).padStart(10)} `
    + `${s.completeness.toFixed(4).padStart(13)} ${r.asserted.size.toLocaleString().padStart(10)}`)
}
{
  const r = rebuild2(kids, true)
  console.log(`\n  the upper bound discards ${(100 * r.discarded).toFixed(4)}% of true `
    + 'heterozygous markers by the binomial, which is the price of it')
}

console.log('\nTHE WINDOW SWEPT, so the bound is seen against what it buys')
console.log(`${'u'.padStart(4)} ${'precision'.padStart(10)} ${'completeness'.padStart(13)} `
  + `${'asserted'.padStart(10)}`)
{
  const base = rebuild2(kids, false)
  const votes = new Map<string, number>()
  for (const g of kids) {
    const gt = sample.get(g)!.gt
    for (const [probe, fg] of faGt) {
      if (fg === 'AB') continue
      if (gt.get(probe) === 'AB') votes.set(probe, (votes.get(probe) ?? 0) + 1)
    }
  }
  for (const u of [12, 14, 16, 18, 20, 21, 24, 30]) {
    const a = new Set<string>()
    for (const [probe, v] of votes) if (v >= base.k && v <= u) a.add(probe)
    const s = grade(a)
    console.log(`${String(u).padStart(4)} ${s.precision.toFixed(4).padStart(10)} `
      + `${s.completeness.toFixed(4).padStart(13)} ${a.size.toLocaleString().padStart(10)}`)
  }
}

console.log('\nLEAVE ONE EMBRYO OUT, two-sided')
const ps: number[] = []
for (const e of embryos) {
  const from = kids.filter((g) => embryoOf.get(g) !== e)
  if (from.length < 2) continue
  const r = rebuild2(from, true)
  const s = grade(r.asserted)
  ps.push(s.precision)
  console.log(`  ${e.padEnd(4)} n ${String(from.length).padStart(3)} k ${String(r.k).padStart(2)} `
    + `u ${String(r.u).padStart(2)}  precision ${s.precision.toFixed(4)}  `
    + `completeness ${s.completeness.toFixed(4)}`)
}
if (ps.length) {
  console.log(`  mean precision ${(ps.reduce((a, b) => a + b, 0) / ps.length).toFixed(4)}`)
}

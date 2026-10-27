/**
 * The informative floor on the series whose reference parent is not at position 00.
 *
 * An external review reports one chromosome of that series at 312 parent-heterozygous markers
 * against a typical 4,000, passing MIN_INFORMATIVE and reporting a heterozygous rate outside any
 * biological range. This measures every autosome of every sample in that series against its
 * measured reference, and reports the informative counts, their fraction of the called markers,
 * and the rates, so the floor question can be answered rather than argued.
 *
 * No laboratory sample identifier is printed; samples are reported by position in the series.
 *
 * Run: OM_TRIOS=<dir> OM_REF=<reference file> node --experimental-strip-types audit/lab-floor.ts
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










const obl = await import(`${W}obligateHet.ts`)
const REF = process.env.OM_REF
if (!REF || !existsSync(join(DIR, REF))) {
  console.log('lab-floor: OM_REF not readable. NOT RUN.')
  process.exit(0)
}
const refGt = new Map<string, string>()
for (const r of rowsOf(join(DIR, REF))()) {
  if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') refGt.set(r.probesetId, r.genotype)
}
console.log(`reference array: ${refGt.size.toLocaleString()} called markers\n`)

const prefix = /^(\d{4,6})-/.exec(REF)?.[1] ?? ''
const samples = readdirSync(DIR)
  .filter((f) => f.startsWith(`${prefix}-`) && /\.probes/.test(f) && f !== REF).sort()

type Acc = { nHom: number; hetHom: number; nHet: number; hetHet: number; called: number }
const all: { pos: string; chrom: string; a: Acc }[] = []
for (const f of samples) {
  const pos = /-(\d{2})\.probes/.exec(f)?.[1] ?? '??'
  const by = new Map<string, Acc>()
  for (const r of rowsOf(join(DIR, f))()) {
    if (!/^\d+$/.test(r.chrom)) continue
    const s = r.genotype
    if (s !== 'AA' && s !== 'BB' && s !== 'AB') continue
    let a = by.get(r.chrom)
    if (!a) { a = { nHom: 0, hetHom: 0, nHet: 0, hetHet: 0, called: 0 }; by.set(r.chrom, a) }
    a.called += 1
    const rg = refGt.get(r.probesetId)
    if (!rg) continue
    if (rg === 'AB') { a.nHet += 1; if (s === 'AB') a.hetHet += 1 }
    else { a.nHom += 1; if (s === 'AB') a.hetHom += 1 }
  }
  for (const [chrom, a] of by) all.push({ pos, chrom, a })
  process.stdout.write('.')
}
console.log('\n')

const q = (xs: number[], p: number) => {
  const s = xs.slice().sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.max(0, Math.round(p * (s.length - 1))))] : NaN
}
const FLOOR = obl.MIN_INFORMATIVE as number
const medHet = q(all.map((x) => x.a.nHet), 0.5)
const medHom = q(all.map((x) => x.a.nHom), 0.5)
console.log(`${all.length} chromosome-observations over ${samples.length} samples`)
console.log(`median parent-HET count ${medHet.toLocaleString()}, parent-HOM ${medHom.toLocaleString()}`)
console.log(`\nTHE THINNEST CHROMOSOMES, which is where an absolute floor of ${FLOOR} lets things through`)
console.log(`${'pos'.padEnd(5)} ${'chr'.padStart(3)} ${'nHom'.padStart(8)} ${'nHet'.padStart(7)} `
  + `${'het frac'.padStart(9)} ${'hetHom'.padStart(7)} ${'hetHet'.padStart(7)}  passes the floor?`)
let admitted = 0
for (const x of all.slice().sort((a, b) => a.a.nHet - b.a.nHet).slice(0, 10)) {
  const frac = x.a.nHet / Math.max(1, medHet)
  const passes = x.a.nHom >= FLOOR || x.a.nHet >= FLOOR
  if (passes && frac < 0.25) admitted += 1
  console.log(`${x.pos.padEnd(5)} ${x.chrom.padStart(3)} ${x.a.nHom.toLocaleString().padStart(8)} `
    + `${x.a.nHet.toLocaleString().padStart(7)} ${(100 * frac).toFixed(1).padStart(8)}% `
    + `${(x.a.hetHom / Math.max(1, x.a.nHom)).toFixed(4).padStart(7)} `
    + `${(x.a.hetHet / Math.max(1, x.a.nHet)).toFixed(4).padStart(7)}  ${passes ? 'yes' : 'no'}`)
}
const pass = all.filter((x) => x.a.nHom >= FLOOR || x.a.nHet >= FLOOR)
const thin = pass.filter((x) => x.a.nHet < medHet * 0.25)
console.log(`\nchromosomes passing the absolute floor whose parent-HET count is under a quarter of `
  + `the median: ${thin.length} of ${pass.length}`)
console.log('a fractional floor would catch those; an absolute one does not.')

/**
 * Which reference catches a chromosome that is not being measured, without penalising a small one.
 *
 * A fraction of the sample's median AUTOSOME refuses chr22 for being chr22. Three candidates are
 * printed per chromosome so the right denominator can be chosen from the data:
 *   informative / called on that chromosome
 *   informative / the markers that chromosome carries on the array
 *   called / the markers it carries, which is the chromosome's own call rate
 *
 * Run: OM_TRIOS=<dir> OM_REF=<reference> node --experimental-strip-types audit/floor-reference.ts
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
if (!REF || !existsSync(join(DIR, REF))) { console.log('NOT RUN'); process.exit(0) }
const refGt = new Map<string, string>()
for (const r of rowsOf(join(DIR, REF))()) {
  if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') refGt.set(r.probesetId, r.genotype)
}
const prefix = /^(\d{4,6})-/.exec(REF)?.[1] ?? ''
const samples = readdirSync(DIR)
  .filter((f) => f.startsWith(`${prefix}-`) && /\.probes/.test(f) && f !== REF).sort()

type Acc = { onArray: number; called: number; inf: number; het: number }
const rows: { pos: string; chrom: string; a: Acc }[] = []
for (const f of samples) {
  const pos = /-(\d{2})\.probes/.exec(f)?.[1] ?? '??'
  const by = new Map<string, Acc>()
  for (const r of rowsOf(join(DIR, f))()) {
    if (!/^\d+$/.test(r.chrom)) continue
    let a = by.get(r.chrom)
    if (!a) { a = { onArray: 0, called: 0, inf: 0, het: 0 }; by.set(r.chrom, a) }
    a.onArray += 1
    const s = r.genotype
    if (s !== 'AA' && s !== 'BB' && s !== 'AB') continue
    a.called += 1
    const rg = refGt.get(r.probesetId)
    if (!rg || rg === 'AB') continue
    a.inf += 1
    if (s === 'AB') a.het += 1
  }
  for (const [chrom, a] of by) rows.push({ pos, chrom, a })
  process.stdout.write('.')
}
console.log('\n')
console.log(`${'pos'.padEnd(5)} ${'chr'.padStart(3)} ${'onArray'.padStart(8)} ${'called'.padStart(8)} `
  + `${'inf'.padStart(7)} ${'inf/called'.padStart(11)} ${'inf/array'.padStart(10)} `
  + `${'callRate'.padStart(9)} ${'het'.padStart(7)}`)
for (const x of rows.slice().sort((a, b) => (a.a.inf / Math.max(1, a.a.called)) - (b.a.inf / Math.max(1, b.a.called))).slice(0, 8)) {
  const a = x.a
  console.log(`${x.pos.padEnd(5)} ${x.chrom.padStart(3)} ${a.onArray.toLocaleString().padStart(8)} `
    + `${a.called.toLocaleString().padStart(8)} ${a.inf.toLocaleString().padStart(7)} `
    + `${(a.inf / Math.max(1, a.called)).toFixed(4).padStart(11)} `
    + `${(a.inf / Math.max(1, a.onArray)).toFixed(4).padStart(10)} `
    + `${(a.called / Math.max(1, a.onArray)).toFixed(3).padStart(9)} `
    + `${(a.het / Math.max(1, a.inf)).toFixed(4).padStart(7)}`)
}
const q = (xs: number[], p: number) => { const s = xs.slice().sort((a, b) => a - b)
  return s.length ? s[Math.round(p * (s.length - 1))] : NaN }
for (const [name, fn] of [['inf/called', (a: Acc) => a.inf / Math.max(1, a.called)],
  ['inf/onArray', (a: Acc) => a.inf / Math.max(1, a.onArray)],
  ['callRate', (a: Acc) => a.called / Math.max(1, a.onArray)]] as const) {
  const xs = rows.map((x) => fn(x.a))
  console.log(`\n  ${name.padEnd(12)} min ${Math.min(...xs).toFixed(4)}  p5 ${q(xs, 0.05).toFixed(4)}`
    + `  median ${q(xs, 0.5).toFixed(4)}`)
}
void obl

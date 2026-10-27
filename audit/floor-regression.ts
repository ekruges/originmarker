/**
 * Regression for the fractional informative floor: how many chromosomes lose their contribution
 * call to it, on the public corpus and on the laboratory series where the defect was found.
 *
 * The pre-registration's control is that nothing in the public corpus falls below 81 percent of
 * its own median, so the floor should refuse nothing there and refuse exactly one chromosome of
 * the laboratory series.
 *
 * Run: OM_TRIOS=<dir> [OM_REF=<reference>] node --experimental-strip-types audit/floor-regression.ts
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
const FLOOR = obl.MIN_INFORMATIVE as number

const refFile = process.env.OM_REF
const manifestPath = join(DIR, 'trio_manifest_full.csv')
type Pair = { gsm: string; father: string; mother?: string }
const pairs: Pair[] = []
if (existsSync(manifestPath)) {
  const manifest = readFileSync(manifestPath, 'utf8').trim().split('\n')
  const cols = manifest[0].split(',')
  const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
  const hand = new Set<string>()
  for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
    if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
    hand.add(l.split('\t')[0])
  }
  for (const line of manifest.slice(1)) {
    const r = line.split(',')
    const g = at(r, 'gsm')
    const f = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
    const m = at(r, 'mother_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
    if (g && f && m && fileOf.has(g) && hand.has(g)) pairs.push({ gsm: g, father: f, mother: m })
  }
} else if (refFile) {
  const prefix = /^(\d{4,6})-/.exec(refFile)?.[1] ?? ''
  for (const f of readdirSync(DIR)) {
    if (!f.startsWith(`${prefix}-`) || !/\.probes/.test(f) || f === refFile) continue
    pairs.push({ gsm: f, father: refFile })
  }
}
console.log(`${pairs.length} samples, floor ${FLOOR} absolute, plus the chromosome's own call rate\n`)

let refused = 0; let scored = 0
const hits: string[] = []
for (const p of pairs) {
  const refGt = new Map<string, string>()
  const src = fileOf.get(p.father) ?? join(DIR, p.father)
  for (const r of rowsOf(src)()) {
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') refGt.set(r.probesetId, r.genotype)
  }
  const sampleSrc = fileOf.get(p.gsm) ?? join(DIR, p.gsm)
  const by = new Map<string, { inf: number; called: number; onArray: number }>()
  for (const r of rowsOf(sampleSrc)()) {
    if (!/^\d+$/.test(r.chrom)) continue
    let a = by.get(r.chrom)
    if (!a) { a = { inf: 0, called: 0, onArray: 0 }; by.set(r.chrom, a) }
    a.onArray += 1
    const s = r.genotype
    if (s !== 'AA' && s !== 'BB' && s !== 'AB') continue
    a.called += 1
    const rg = refGt.get(r.probesetId)
    if (!rg || rg === 'AB') continue
    a.inf += 1
  }
  for (const [chrom, a] of by) {
    const was = a.inf >= FLOOR
    const cr = a.onArray > 0 ? a.called / a.onArray : NaN
    const now = obl.informativeEnough(a.inf, cr) as boolean
    if (was) scored += 1
    if (was && !now) {
      refused += 1
      hits.push(`${p.gsm.replace(/\.probes.*/, '')} chr${chrom}: ${a.inf} informative, `
        + `call rate ${cr.toFixed(3)}`)
    }
  }
}
console.log(`chromosomes that passed the absolute floor:        ${scored}`)
console.log(`of those, refused by the fractional floor:         ${refused}`)
for (const h of hits) console.log(`  ${h}`)

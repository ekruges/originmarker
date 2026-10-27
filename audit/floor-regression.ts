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

// The third arm is the DISCARDED rule: informative count against a fraction of the sample's own
// median autosome. Its reference is the sample itself, so a sample whose every chromosome has
// collapsed has an equally collapsed median and every chromosome passes. That is the failure mode,
// and it is the opposite of refusing small chromosomes for being small.
//
// Two criteria, both measured here:
//   1. the shipped rule refuses no chromosome whose own call rate is healthy
//   2. the shipped rule catches collapsed chromosomes that the discarded rule lets through
const MEDIAN_FRACTION = 0.20
const HEALTHY = 0.60

let refused = 0; let scored = 0
let healthyRefusedNow = 0
let refusedOld = 0; let healthyRefusedOld = 0
const smallOld = new Map<string, number>()
const hits: string[] = []
const oldHits: string[] = []
let caughtOnly = 0
let worstLet = 1
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
  const infs = [...by.values()].map((a) => a.inf).filter((n) => n >= FLOOR).sort((x, y) => x - y)
  const median = infs.length ? infs[Math.floor(infs.length / 2)] : 0
  for (const [chrom, a] of by) {
    const was = a.inf >= FLOOR
    const cr = a.onArray > 0 ? a.called / a.onArray : NaN
    const now = obl.informativeEnough(a.inf, cr) as boolean
    const oldRule = a.inf >= MEDIAN_FRACTION * median
    if (was) scored += 1
    if (was && !oldRule) {
      refusedOld += 1
      oldHits.push(`${p.gsm.replace(/\.probes.*/, '')} chr${chrom}: ${a.inf} informative, `
        + `${(100 * a.inf / Math.max(1, median)).toFixed(0)}% of the median ${median}, `
        + `call rate ${cr.toFixed(3)}`)
      if (cr >= HEALTHY) {
        healthyRefusedOld += 1
        smallOld.set(chrom, (smallOld.get(chrom) ?? 0) + 1)
      }
    }
    if (was && !now && cr >= HEALTHY) healthyRefusedNow += 1
    if (was && !now && oldRule) { caughtOnly += 1; if (cr < worstLet) worstLet = cr }
    if (was && !now) {
      refused += 1
      hits.push(`${p.gsm.replace(/\.probes.*/, '')} chr${chrom}: ${a.inf} informative, `
        + `call rate ${cr.toFixed(3)}`)
    }
  }
}
console.log(`chromosomes that passed the absolute floor:        ${scored}`)
console.log(`refused by the SHIPPED rule (floor + call rate):   ${refused}`)
console.log(`  of those, with a healthy call rate:             ${healthyRefusedNow}`)
console.log(`refused by the DISCARDED rule (median fraction):   ${refusedOld}`)
console.log(`  of those, with a healthy call rate:             ${healthyRefusedOld}`)
if (smallOld.size) {
  const by2 = [...smallOld.entries()].sort((a, b) => b[1] - a[1])
  console.log(`  which chromosomes those were: `
    + by2.map(([c, n]) => `chr${c} x${n}`).join(', '))
}
console.log()
console.log(`collapsed, caught by the shipped rule only:        ${caughtOnly}`)
console.log(`  worst call rate the discarded rule let through:  ${worstLet.toFixed(3)}`)
console.log()
if (healthyRefusedNow > 0) {
  console.log(`FAIL criterion 1: the shipped rule refuses ${healthyRefusedNow} chromosome(s) with a `
    + 'healthy call rate.')
} else if (caughtOnly > 0) {
  console.log('PASS: the shipped rule refuses nothing whose own call rate is healthy, and it catches')
  console.log(`      ${caughtOnly} collapsed chromosome(s) the discarded rule passes, down to a call`)
  console.log(`      rate of ${worstLet.toFixed(3)}. The discarded rule measures a collapsed sample`)
  console.log('      against its own collapsed median, so it cannot see this.')
} else {
  console.log('INCONCLUSIVE: the two rules agree on this corpus.')
}
console.log('\nrefused by the DISCARDED rule, in full:')
for (const h of oldHits) console.log(`  ${h}`)
console.log('\nrefused by the shipped rule, in full:')
for (const h of hits) console.log(`  ${h}`)

/**
 * WHICH GATE WITHHELD A ROW, in the gate's own words.
 *
 * A row that reports no class can be stopped at several points and they are not interchangeable:
 * an implausible heterozygous band says the array failed, a low call rate says zygosity cannot be
 * read, and an absence between the ceiling and 3x it says the sample is real but the measurement
 * does not separate. A sheet that prints one reason where another fired is telling a reader the
 * wrong thing about the material.
 *
 * So this prints, for named samples, the stage verdict and its reasoning, the zygosity, the
 * absence verdict, and every limit string the pipeline attached. Nothing is summarised.
 *
 * Run: OM_LAB=<dir> OM_PAIRS=<tsv> OM_FOCUS=<barcode,...> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/spot-check-gates.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_LAB
const PAIRS = process.env.OM_PAIRS
const FOCUS = new Set((process.env.OM_FOCUS ?? '').split(',').map((x) => x.trim()).filter(Boolean))
if (!DIR || !existsSync(DIR) || !PAIRS || !existsSync(PAIRS) || !FOCUS.size) {
  console.log('spot-check-gates: OM_LAB, OM_PAIRS and OM_FOCUS must all be set. NOT RUN.')
  process.exit(0)
}

const fileOf = new Map<string, string>()
const walk = (d: string) => {
  for (const e of readdirSync(d)) {
    const p = join(d, e)
    if (statSync(p).isDirectory()) { walk(p); continue }
    const m = /(\d{4,6}-\d{2})/.exec(e)
    if (m && /\.probes(\.txt)?(\.gz)?$/.test(e) && !fileOf.has(m[1])) fileOf.set(m[1], p)
  }
}
walk(DIR)

function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () { for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r } }
}
const cache = new Map<string, unknown>()
function parentIndex(bc: string) {
  const hit = cache.get(bc); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(bc)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(bc, byChrom as never, bafSums as never, first).build.build)
  cache.set(bc, idx); return idx
}

const pairs = readFileSync(PAIRS, 'utf8').trim().split('\n')
  .filter((l) => l && !l.startsWith('#'))
  .map((l) => { const [sample, father, mother] = l.split('\t').map((x) => x.trim()); return { sample, father, mother } })
  .filter((p) => FOCUS.has(p.sample))

const pct = (x: number, d = 2) => (Number.isFinite(x) ? `${(x * 100).toFixed(d)}%` : '-')

for (const p of pairs) {
  if (!fileOf.has(p.sample)) { console.log(`${p.sample}: no array\n`); continue }
  const hasF = !!p.father && fileOf.has(p.father)
  const hasM = !!p.mother && fileOf.has(p.mother)
  if (!hasF && !hasM) { console.log(`${p.sample}: no parental array\n`); continue }
  const soloRole = hasF ? 'paternal' : 'maternal'
  const pat = hasF ? parentIndex(p.father) : parentIndex(p.mother)
  const mat = hasF && hasM ? parentIndex(p.mother) : null
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(p.sample)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(p.sample, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat, soloRole, sibs: [], sampleName: p.sample, log: () => {},
  } as never) as {
    stage?: { stage?: string; why?: string; dropout?: number }
    zygosity?: string; verdict?: string; originClass?: string
    genomeRate?: number; explainable?: number; limits?: string[]; notes?: string[]
  }
  console.log(`=== ${p.sample}  (parent loaded: ${hasF ? p.father : p.mother}, role ${soloRole}) ===`)
  console.log(`  call rate            ${pct(profile.callRate, 1)}`)
  console.log(`  heterozygous band    ${pct(profile.hetBand, 1)}`)
  console.log(`  stage                ${res.stage?.stage} - ${res.stage?.why}`)
  console.log(`  zygosity             ${res.zygosity}`)
  console.log(`  absence              ${pct(res.genomeRate ?? NaN)} against a ceiling of ${pct(res.explainable ?? NaN)}`)
  console.log(`  absence verdict      ${res.verdict}`)
  console.log(`  origin class         ${res.originClass ?? '-'}`)
  for (const l of res.limits ?? []) console.log(`  LIMIT  ${l}`)
  console.log()
}

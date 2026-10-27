/**
 * Arm (b) end to end: the SHIPPED pipeline with only one parent array loaded, scored against the
 * published karyotype.
 *
 * audit/arm-b-naming.ts measures the statistic. This measures what the application emits, which is
 * the thing a reader sees: `oneParent[].parent` on a run where no mother was supplied.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=2400 audit/arm-b-pipeline.ts
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










// Table S3's parent labels, per chromosome.
const label = new Map<string, Map<string, 'maternal' | 'paternal'>>()
const abnormal = new Map<string, Set<string>>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const c = l.split('\t')
  const gsm = c[0]; const chrom = c[3]; const clause = (c[4] ?? '').toLowerCase()
  if (!gsm || !chrom) continue
  if (!abnormal.has(gsm)) abnormal.set(gsm, new Set())
  if (/monosomy|trisomy|del|dup|gain|loss|breakpoint/.test(clause)) abnormal.get(gsm)!.add(chrom)
  const par = /\(pat\)/.test(clause) ? 'paternal' as const
    : /\(mat\)/.test(clause) ? 'maternal' as const : undefined
  if (!par) continue
  // A gain names the parent whose copy is ABSENT in the table's own convention, and this channel
  // reports the absent copy, so both are on the same footing. Gains are excluded anyway: the
  // pre-registration scopes arm (b) to losses.
  if (/\bgain\b|trisomy|\bdup\b/.test(clause) && !/\bdel\b/.test(clause)) continue
  if (!label.has(gsm)) label.set(gsm, new Map())
  label.get(gsm)!.set(chrom, par)
}

const manifestPath = join(DIR, 'trio_manifest_full.csv')
if (!existsSync(manifestPath)) { console.log('no manifest. NOT RUN.'); process.exit(0) }
const manifest = readFileSync(manifestPath, 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''

async function scoreSolo(gsm: string, father: string) {
  const pat = parentIndex(father)
  const acc = score.emptyCollected(pat as never, undefined as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, undefined as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  return await score.scoreSample({
    acc, profile, pat, mat: undefined, soloRole: 'paternal', sibs: [], sampleName: gsm,
    log: () => {},
  } as never) as { oneParent?: { where: string; parent?: string; verdict: string
    oneArrayNamed?: boolean; why: string }[] }
}

let right = 0; let wrong = 0; let unnamed = 0; let noRow = 0; let falseNamed = 0
const wrongDetail: string[] = []
const falseDetail: string[] = []
let cells = 0
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const gsm = at(r, 'gsm')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  if (!gsm || !father || !fileOf.has(gsm) || !label.has(gsm)) continue
  let res
  try { res = await scoreSolo(gsm, father) } catch (e) {
    console.log(`  ${gsm}: ${(e as Error).message}`); continue
  }
  cells += 1
  const rows = res.oneParent ?? []
  const want = label.get(gsm)!
  for (const [chrom, par] of want) {
    const hit = rows.find((x) => new RegExp(`\\bchr${chrom}\\b`).test(x.where))
    if (!hit) { noRow += 1; continue }
    if (!hit.parent) { unnamed += 1; continue }
    if (hit.parent === par) right += 1
    else { wrong += 1; wrongDetail.push(`${gsm} chr${chrom}: hand ${par}, emitted ${hit.parent}`) }
  }
  // Rows the pipeline raised on chromosomes the table calls normal, and named.
  for (const x of rows) {
    const m = /\bchr(\d+)\b/.exec(x.where)
    if (!m || !x.parent) continue
    if (abnormal.get(gsm)?.has(m[1])) continue
    falseNamed += 1
    falseDetail.push(`${gsm} ${x.where}: named ${x.parent}, verdict ${x.verdict}`)
  }
}

console.log(`\n${cells} cells scored with ONLY THE FATHER LOADED, through the shipped pipeline\n`)
console.log('HAND-CALLED LOSSES CARRYING A PARENT LABEL')
console.log(`  named correctly:        ${right}`)
console.log(`  named wrongly:          ${wrong}`)
for (const d of wrongDetail) console.log(`    ${d}`)
console.log(`  row raised, not named:  ${unnamed}`)
console.log(`  no row raised at all:   ${noRow}`)
console.log(`\nNAMED ON A CHROMOSOME THE TABLE READS AS NORMAL: ${falseNamed}`)
for (const d of falseDetail.slice(0, 12)) console.log(`    ${d}`)
if (falseDetail.length > 12) console.log(`    and ${falseDetail.length - 12} more`)
console.log(right > 0 && wrong * 2 < right
  ? `\nPASS: the pipeline names ${right} of these with one array and gets ${wrong} wrong. `
    + 'Before this round it named none of them.'
  : `\nFAIL: ${right} right against ${wrong} wrong.`)

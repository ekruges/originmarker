/**
 * Detection re-measured POSITIONALLY, against coordinates rather than chromosome names.
 *
 * Every detection figure this directory has published joined a called segment to a hand-called
 * event by chromosome alone, so a segment anywhere on a chromosome counted as detecting an event
 * anywhere else on it. The regions come from audit/karyotype_regions.py, which resolves every one
 * of Table S3's 28 segmental clauses to GRCh37 coordinates out of the UCSC hg19 cytoBand track.
 *
 * Reports, per hand-called clause: whether any called segment overlaps it at all, the best
 * reciprocal overlap achieved, and which parent the overlapping segment names.
 *
 * Run: OM_TRIOS=<dir> OM_REGIONS=<regions.tsv> node --experimental-strip-types audit/positional-detection.ts
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




type Seg = { chrom: string; kind: string; startBp: number; endBp: number; score: number
  absentParent?: string; contested?: boolean }
type Res = { segments?: Seg[]; stage?: { stage: string } }

async function scoreCell(gsm: string, father: string, mother: string) {
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
  return await score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
  } as never) as Res
}

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const parentsOf = new Map<string, { father: string; mother: string }>()
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const f = at(r, 'father_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  const m = at(r, 'mother_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  if (f && m) parentsOf.set(at(r, 'gsm'), { father: f, mother: m })
}

// The hand-called clauses with coordinates, written by audit/karyotype_regions.py.
const REG = process.env.OM_REGIONS ?? '/root/regions.tsv'
type Ev = { gsm: string; chrom: string; start: number; end: number; clause: string
  parent: string; kind: 'del' | 'dup' }
const events: Ev[] = []
for (const l of readFileSync(REG, 'utf8').split('\n')) {
  if (!l.trim() || l.startsWith('#')) continue
  const [gsm, chrom, start, end, clause] = l.split('\t')
  const parent = /\(\s*pat/i.test(clause) ? 'paternal'
    : /\(\s*mat/i.test(clause) ? 'maternal' : ''
  if (!parent) continue
  events.push({ gsm, chrom, start: Number(start), end: Number(end), clause, parent,
    kind: /dup/i.test(clause) && !/del/i.test(clause) ? 'dup' : 'del' })
}
console.log(`${events.length} hand-called clauses with a named parent and real coordinates\n`)

const r0 = (a: [number, number], b: [number, number]) => {
  const inter = Math.min(a[1], b[1]) - Math.max(a[0], b[0])
  return inter <= 0 ? 0 : inter / (Math.max(a[1], b[1]) - Math.min(a[0], b[0]))
}
const cells = [...new Set(events.map((e) => e.gsm))].filter((g) => parentsOf.has(g)).sort()
const segsOf = new Map<string, Seg[]>()
for (const g of cells) {
  const p = parentsOf.get(g)!
  const res = await scoreCell(g, p.father, p.mother)
  segsOf.set(g, res.segments ?? [])
  process.stdout.write('.')
}
console.log('\n')

console.log('EVERY HAND-CALLED CLAUSE, AGAINST WHAT WAS CALLED IN THAT REGION')
console.log(`${'cell'.padEnd(12)} ${'chr'.padStart(3)} ${'hand region'.padEnd(16)} `
  + `${'kind'.padEnd(4)} ${'par'.padEnd(9)} ${'best overlap'.padEnd(12)} called there`)
let anyOverlap = 0; let half = 0; let chromOnly = 0
const byKind = { del: [0, 0, 0], dup: [0, 0, 0] }
for (const e of events.filter((x) => segsOf.has(x.gsm))) {
  const segs = segsOf.get(e.gsm)!.filter((s) => s.chrom === e.chrom)
  let best = 0; let bestSeg: Seg | undefined
  for (const s of segs) {
    const r = r0([e.start, e.end], [s.startBp, s.endBp])
    if (r > best) { best = r; bestSeg = s }
  }
  const hit = best > 0
  if (hit) anyOverlap += 1
  if (best >= 0.5) half += 1
  if (!hit && segs.length) chromOnly += 1
  const slot = byKind[e.kind]
  slot[0] += hit ? 1 : 0; slot[1] += best >= 0.5 ? 1 : 0; slot[2] += 1
  console.log(`${e.gsm.padEnd(12)} ${e.chrom.padStart(3)} `
    + `${`${(e.start / 1e6).toFixed(1)}-${(e.end / 1e6).toFixed(1)}Mb`.padEnd(16)} `
    + `${e.kind.padEnd(4)} ${e.parent.padEnd(9)} ${best.toFixed(3).padEnd(12)} `
    + (bestSeg
      ? `${bestSeg.kind} ${(bestSeg.startBp / 1e6).toFixed(1)}-${(bestSeg.endBp / 1e6).toFixed(1)}Mb`
        + ` ${bestSeg.absentParent ?? '-'}${bestSeg.absentParent && bestSeg.absentParent !== e.parent ? ' WRONG PARENT' : ''}`
      : (segs.length ? `NOTHING IN REGION, but ${segs.length} elsewhere on the chromosome` : 'nothing on the chromosome')))
}
console.log('\nDETECTION, POSITIONAL')
console.log(`  clauses with ANY overlap:            ${anyOverlap} of ${events.length}`)
console.log(`  clauses at reciprocal overlap >=0.5: ${half} of ${events.length}`)
console.log(`  clauses a CHROMOSOME join would have counted but position does not: ${chromOnly}`)
console.log(`  deletions: any overlap ${byKind.del[0]} of ${byKind.del[2]}, >=0.5 ${byKind.del[1]}`)
console.log(`  duplications: any overlap ${byKind.dup[0]} of ${byKind.dup[2]}, >=0.5 ${byKind.dup[1]}`)

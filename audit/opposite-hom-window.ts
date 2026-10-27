/**
 * The opposite-homozygote survivor share INSIDE each hand-called region, with the floor each one
 * implies. Corrects a run that computed the share over whole chromosomes and so scored a
 * segmental deletion covering 27 percent of a chromosome as though the whole chromosome were lost.
 *
 * Per hand-called clause: the opposite-homozygote markers inside its coordinates, the share
 * reading one parent's allele alone there, the same over the rest of that chromosome as the local
 * null, and N = threshold / KL(event || null) from those two.
 *
 * Run: OM_TRIOS=<dir> OM_REGIONS=<regions.tsv> node --experimental-strip-types audit/opposite-hom-window.ts
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







const segs2 = await import(`${W}segments.ts`)
const THRESHOLD = segs2.SEGMENT_LRT as number
const kl = (p: number, q: number) => (p <= 0 || p >= 1 || q <= 0 || q >= 1 ? NaN
  : p * Math.log(p / q) + (1 - p) * Math.log((1 - p) / (1 - q)))

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
type Ev = { gsm: string; chrom: string; start: number; end: number; clause: string; parent: string }
const events: Ev[] = []
for (const l of readFileSync(process.env.OM_REGIONS ?? '/root/regions.tsv', 'utf8').split('\n')) {
  if (!l.trim() || l.startsWith('#')) continue
  const [gsm, chrom, start, end, clause] = l.split('\t')
  const parent = /\(\s*pat/i.test(clause) ? 'paternal' : /\(\s*mat/i.test(clause) ? 'maternal' : ''
  // Deletion-bearing only: an absence statistic cannot see a copy present twice.
  if (!parent || !/\bdel/i.test(clause)) continue
  if (parentsOf.has(gsm)) events.push({ gsm, chrom, start: Number(start), end: Number(end), clause, parent })
}
console.log(`${events.length} deletion-bearing clauses with coordinates and both parents\n`)

console.log(`threshold ${THRESHOLD}`)
console.log(`${'cell'.padEnd(12)} ${'chr'.padStart(3)} ${'par'.padEnd(9)} ${'n in'.padStart(6)} `
  + `${'share in'.padStart(9)} ${'n out'.padStart(6)} ${'share out'.padStart(10)} `
  + `${'floor'.padStart(8)}  reachable?`)
let ok = 0
const byCell = new Map<string, Ev[]>()
for (const e of events) byCell.set(e.gsm, [...(byCell.get(e.gsm) ?? []), e])

for (const [gsm, evs] of [...byCell].sort()) {
  const p = parentsOf.get(gsm)!
  const fa = parentIndex(p.father) as unknown as { gt: Map<string, string> }
  const mo = parentIndex(p.mother) as unknown as { gt: Map<string, string> }
  // The opposite-homozygote markers, with the sample's genotype and position.
  type M = { chrom: string; pos: number; fatherAllele: string; motherAllele: string; gt: string }
  const ms: M[] = []
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    const fg = fa.gt.get(r.probesetId); const mg = mo.gt.get(r.probesetId)
    if (!fg || !mg) continue
    if (!((fg === 'AA' && mg === 'BB') || (fg === 'BB' && mg === 'AA'))) continue
    if (r.genotype !== 'AA' && r.genotype !== 'BB' && r.genotype !== 'AB') continue
    ms.push({ chrom: r.chrom, pos: r.pos, fatherAllele: fg[0], motherAllele: mg[0], gt: r.genotype })
  }
  for (const e of evs) {
    // The share reading the SURVIVING parent's allele alone, which is the other parent's.
    const want = e.parent === 'paternal' ? 'motherAllele' : 'fatherAllele'
    const share = (set: M[]) => {
      const n = set.length
      if (!n) return [NaN, 0] as const
      let alone = 0
      for (const m of set) {
        const surv = m[want]
        if (m.gt === `${surv}${surv}`) alone += 1
      }
      return [alone / n, n] as const
    }
    const here = ms.filter((m) => m.chrom === e.chrom && m.pos >= e.start && m.pos <= e.end)
    const rest = ms.filter((m) => m.chrom === e.chrom && (m.pos < e.start || m.pos > e.end))
    const [si, ni] = share(here)
    const [so, no] = share(rest)
    const d = kl(si, so)
    const floor = d > 0 ? THRESHOLD / d : Infinity
    const reach = Number.isFinite(floor) && floor <= ni
    if (reach) ok += 1
    console.log(`${gsm.padEnd(12)} ${e.chrom.padStart(3)} ${e.parent.padEnd(9)} `
      + `${String(ni).padStart(6)} ${(Number.isFinite(si) ? si.toFixed(3) : '-').padStart(9)} `
      + `${String(no).padStart(6)} ${(Number.isFinite(so) ? so.toFixed(3) : '-').padStart(10)} `
      + `${(Number.isFinite(floor) ? Math.ceil(floor).toLocaleString() : 'inf').padStart(8)}  `
      + `${reach ? 'YES' : 'no'}`)
  }
}
console.log(`\n  clauses whose derived floor fits in the markers the region carries: ${ok} of ${events.length}`)
console.log('  share in  = the surviving parent\'s allele alone, inside the hand-called region')
console.log('  share out = the same over the rest of that chromosome, which is the local null')

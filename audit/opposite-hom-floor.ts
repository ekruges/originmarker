/**
 * Step 1 of audit/opposite-hom-floor-preregistration.txt: the two rates on the
 * OPPOSITE-HOMOZYGOTE set, measured per cell and per chromosome against the hand karyotype,
 * before any floor is proposed.
 *
 * U1 intact chromosomes bounded above by 0.627 per parent, rising with the cell's dropout.
 * U2 lost chromosomes at or above 0.735 for the survivor.
 * U3 the two populations do not overlap.
 *
 * Then the floor each cell implies, N = threshold / KL(event || null), from its own rates.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/opposite-hom-floor.ts
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






type Contribution = { informative: number; fraction: number; onlyA: number; onlyB: number }
type Chrom = { chrom: string; verdict: string; contribution?: Contribution }
type Res = { chroms?: Chrom[]; stage?: { stage: string; dropout: number } }

const segs2 = await import(`${W}segments.ts`)
const obl = await import(`${W}obligateHet.ts`)
const THRESHOLD = segs2.SEGMENT_LRT as number

const kl = (p: number, q: number) => (p <= 0 || p >= 1 || q <= 0 || q >= 1 ? NaN
  : p * Math.log(p / q) + (1 - p) * Math.log((1 - p) / (1 - q)))

async function run(gsm: string, father: string, mother: string) {
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
  } as never) as Res
  return { res, cr: (profile as unknown as { callRate: number }).callRate }
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
const hand = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , , chrom, manual] = l.split('\t')
  hand.set(`${gsm}|${chrom}`, manual)
}
const cells = [...new Set([...hand.keys()].map((k) => k.split('|')[0]))]
  .filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()

type Row = { gsm: string; chrom: string; inf: number; het: number; share: number
  kind: 'intact' | 'lost' | 'other'; dropout: number; cr: number }
const rows: Row[] = []
for (const g of cells) {
  const p = parentsOf.get(g)!
  const { res, cr } = await run(g, p.father, p.mother)
  for (const c of (res.chroms ?? [])) {
    const k = c.contribution
    if (!/^\d+$/.test(c.chrom) || !k || k.informative < 200) continue
    const h = hand.get(`${g}|${c.chrom}`) ?? ''
    const kind = /^normal/i.test(h.trim()) ? 'intact' as const
      : /del|monosom|nullisom/i.test(h) ? 'lost' as const : 'other' as const
    rows.push({ gsm: g, chrom: c.chrom, inf: k.informative, het: k.fraction,
      share: Math.max(k.onlyA, k.onlyB), kind,
      dropout: res.stage?.dropout ?? NaN, cr })
  }
  process.stdout.write('.')
}
console.log('\n')

const q = (xs: number[], p: number) => {
  if (!xs.length) return NaN
  const s = xs.slice().sort((a, b) => a - b)
  return s[Math.min(s.length - 1, Math.max(0, Math.round(p * (s.length - 1))))]
}
console.log('U1 AND U2. THE SURVIVOR SHARE ON THE OPPOSITE-HOMOZYGOTE SET, BY HAND READING')
console.log(`SURVIVOR_MIN ${obl.SURVIVOR_MIN}; the calibration says intact < 0.627 and lost >= 0.735`)
console.log(`${'hand reads'.padEnd(14)} ${'n'.padStart(5)} ${'p10'.padStart(7)} `
  + `${'median'.padStart(7)} ${'p90'.padStart(7)} ${'max'.padStart(7)}`)
for (const k of ['intact', 'lost', 'other'] as const) {
  const xs = rows.filter((r) => r.kind === k).map((r) => r.share)
  if (!xs.length) continue
  console.log(`${k.padEnd(14)} ${String(xs.length).padStart(5)} ${q(xs, 0.1).toFixed(3).padStart(7)} `
    + `${q(xs, 0.5).toFixed(3).padStart(7)} ${q(xs, 0.9).toFixed(3).padStart(7)} `
    + `${Math.max(...xs).toFixed(3).padStart(7)}`)
}
const intact = rows.filter((r) => r.kind === 'intact').map((r) => r.share)
const lost = rows.filter((r) => r.kind === 'lost').map((r) => r.share)
console.log(`\nU3. OVERLAP: intact max ${Math.max(...intact).toFixed(3)} against `
  + `lost min ${Math.min(...lost).toFixed(3)}  `
  + `${Math.max(...intact) < Math.min(...lost) ? 'NO OVERLAP' : 'THEY OVERLAP'}`)

console.log('\nTHE INTACT SHARE AGAINST THE CELL CALL RATE, since dropout is what raises it')
console.log(`${'call rate'.padEnd(12)} ${'n'.padStart(5)} ${'median intact share'.padStart(20)}`)
for (const [lo, hi] of [[0, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.01]] as [number, number][]) {
  const xs = rows.filter((r) => r.kind === 'intact' && r.cr >= lo && r.cr < hi).map((r) => r.share)
  if (xs.length) {
    console.log(`${`${lo.toFixed(2)}-${hi.toFixed(2)}`.padEnd(12)} `
      + `${String(xs.length).padStart(5)} ${q(xs, 0.5).toFixed(3).padStart(20)}`)
  }
}

console.log('\nTHE FLOOR EACH CELL IMPLIES, from its own intact and lost shares')
console.log(`threshold ${THRESHOLD}. chr16 carries about 937 of these markers whole, 245 in 16q23.1-tel.`)
console.log(`${'cell'.padEnd(12)} ${'cr'.padStart(5)} ${'intact'.padStart(7)} ${'lost'.padStart(7)} `
  + `${'floor'.padStart(8)} ${'inf/chrom'.padStart(10)}  reaches a 245-marker window?`)
let reach = 0; let total = 0
for (const g of cells) {
  const mine = rows.filter((r) => r.gsm === g)
  const ii = mine.filter((r) => r.kind === 'intact').map((r) => r.share)
  const ll = mine.filter((r) => r.kind === 'lost').map((r) => r.share)
  if (!ii.length || !ll.length) continue
  const nul = q(ii, 0.5); const ev = q(ll, 0.5)
  const d = kl(ev, nul)
  const floor = d > 0 ? THRESHOLD / d : Infinity
  const med = q(mine.map((r) => r.inf), 0.5)
  total += 1
  if (floor <= 245) reach += 1
  console.log(`${g.padEnd(12)} ${mine[0].cr.toFixed(2).padStart(5)} ${nul.toFixed(3).padStart(7)} `
    + `${ev.toFixed(3).padStart(7)} `
    + `${(Number.isFinite(floor) ? Math.ceil(floor).toLocaleString() : 'inf').padStart(8)} `
    + `${med.toLocaleString().padStart(10)}  ${floor <= 245 ? 'YES' : 'no'}`)
}
console.log(`\n  cells whose derived floor fits inside a 245-marker window: ${reach} of ${total}`)

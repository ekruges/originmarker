/**
 * Which marker class carries one parental contribution against two, with ONE parent loaded.
 *
 * Per chromosome and per loaded parent: the heterozygous rate over markers where that parent is
 * HOMOZYGOUS, which is the class the one-parent caller keeps, and over markers where it is
 * HETEROZYGOUS, which it excludes. Truth is the authors' by-eye karyotype.
 *
 * V1-V4 fixed in audit/marker-class-preregistration.txt. Also carries the informative counts and
 * their fraction of each chromosome's called markers, for the floor question.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/marker-class.ts
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
const op = await import(`${W}oneParentOrigin.ts`)

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

type Obs = { gsm: string; chrom: string; loaded: 'father' | 'mother'
  nHom: number; hetHom: number; nHet: number; hetHet: number; called: number
  truth: 'one' | 'two' | 'other' }
const obs: Obs[] = []

for (const g of cells) {
  const p = parentsOf.get(g)!
  const fa = parentIndex(p.father) as unknown as { gt: Map<string, string> }
  const mo = parentIndex(p.mother) as unknown as { gt: Map<string, string> }
  // Per chromosome, per loaded parent: the two marker classes.
  type Acc = { nHom: number; hetHom: number; nHet: number; hetHet: number; called: number }
  const byChrom = new Map<string, { father: Acc; mother: Acc }>()
  const blank = (): Acc => ({ nHom: 0, hetHom: 0, nHet: 0, hetHet: 0, called: 0 })
  for (const r of rowsOf(fileOf.get(g)!)()) {
    if (!/^\d+$/.test(r.chrom)) continue
    const s = r.genotype
    if (s !== 'AA' && s !== 'BB' && s !== 'AB') continue
    let e = byChrom.get(r.chrom)
    if (!e) { e = { father: blank(), mother: blank() }; byChrom.set(r.chrom, e) }
    e.father.called += 1
    e.mother.called += 1
    for (const [who, gt] of [['father', fa.gt.get(r.probesetId)],
      ['mother', mo.gt.get(r.probesetId)]] as const) {
      if (!gt || gt === 'NC') continue
      const a = e[who]
      if (gt === 'AB') { a.nHet += 1; if (s === 'AB') a.hetHet += 1 }
      else { a.nHom += 1; if (s === 'AB') a.hetHom += 1 }
    }
  }
  for (const [chrom, e] of byChrom) {
    const h = hand.get(`${g}|${chrom}`) ?? ''
    const normal = /^normal/i.test(h.trim())
    const lostPat = /\(\s*pat/i.test(h) && /del|monosom|nullisom/i.test(h)
    const lostMat = /\(\s*mat/i.test(h) && /del|monosom|nullisom/i.test(h)
    const whole = /monosom|nullisom/i.test(h)
    for (const who of ['father', 'mother'] as const) {
      const a = e[who]
      // ONE contribution for this loaded parent means the OTHER parent's copy is gone, so the
      // sample carries only this parent's. A segmental loss is not a whole-chromosome statement
      // and is excluded from both truth classes.
      const otherGone = who === 'father' ? lostMat : lostPat
      const truth = normal ? 'two' as const
        : (otherGone && whole) ? 'one' as const : 'other' as const
      obs.push({ gsm: g, chrom, loaded: who, nHom: a.nHom, hetHom: a.hetHom,
        nHet: a.nHet, hetHet: a.hetHet, called: a.called, truth })
    }
  }
  process.stdout.write('.')
}
console.log('\n')

const q = (xs: number[], p: number) => {
  if (!xs.length) return NaN
  const s = xs.slice().sort((a, b) => a - b)
  return s[Math.min(s.length - 1, Math.max(0, Math.round(p * (s.length - 1))))]
}
const FLOOR = obl.MIN_INFORMATIVE as number
console.log(`V1 AND V2. THE TWO CLASSES, at the shipped informative floor of ${FLOOR}`)
for (const [label, floorHom, floorHet] of [['shipped floor', FLOOR, FLOOR],
  ['floor 1000 in both', 1000, 1000]] as [string, number, number][]) {
  console.log(`\n  ${label}`)
  console.log(`  ${'class'.padEnd(22)} ${'truth'.padEnd(6)} ${'n'.padStart(5)} `
    + `${'min'.padStart(7)} ${'median'.padStart(7)} ${'max'.padStart(7)}`)
  for (const [cls, nKey, hKey, fl] of [
    ['parent-HOM (shipped)', 'nHom', 'hetHom', floorHom],
    ['parent-HET (excluded)', 'nHet', 'hetHet', floorHet]] as const) {
    const vals: Record<string, number[]> = { one: [], two: [] }
    for (const o of obs) {
      if (o.truth === 'other') continue
      const n = o[nKey as 'nHom']; const h = o[hKey as 'hetHom']
      if (n < (fl as number)) continue
      vals[o.truth].push(h / n)
    }
    for (const t of ['one', 'two']) {
      const xs = vals[t]
      console.log(`  ${cls.padEnd(22)} ${t.padEnd(6)} ${String(xs.length).padStart(5)} `
        + `${(xs.length ? Math.min(...xs) : NaN).toFixed(4).padStart(7)} `
        + `${q(xs, 0.5).toFixed(4).padStart(7)} `
        + `${(xs.length ? Math.max(...xs) : NaN).toFixed(4).padStart(7)}`)
    }
    const one = vals.one; const two = vals.two
    if (one.length && two.length) {
      const gap = Math.min(...two) - Math.max(...one)
      console.log(`  ${''.padEnd(22)} ${'gap'.padEnd(6)} `
        + `${gap.toFixed(4).padStart(27)}  ${gap > 0 ? 'SEPARATES' : 'OVERLAPS'}`)
    }
  }
}

console.log(`\nV3. THE SHIPPED BOUNDARIES AGAINST THE PARENT-HOM CLASS`)
console.log(`  ONE_PARENT_HAPLOID_MAX ${obl.ONE_PARENT_HAPLOID_MAX}, `
  + `ONE_PARENT_DIPLOID_MIN ${obl.ONE_PARENT_DIPLOID_MIN}`)
{
  const two = obs.filter((o) => o.truth === 'two' && o.nHom >= FLOOR)
  const below = two.filter((o) => o.hetHom / o.nHom < (obl.ONE_PARENT_HAPLOID_MAX as number))
  console.log(`  chromosomes the hand reads NORMAL that fall below the haploid max: `
    + `${below.length} of ${two.length}  (${(100 * below.length / Math.max(1, two.length)).toFixed(1)}%)`)
  const twoHet = obs.filter((o) => o.truth === 'two' && o.nHet >= FLOOR)
  const belowHet = twoHet.filter((o) => o.hetHet / o.nHet < (obl.ONE_PARENT_HAPLOID_MAX as number))
  console.log(`  the same on the parent-HET class: ${belowHet.length} of ${twoHet.length}  `
    + `(${(100 * belowHet.length / Math.max(1, twoHet.length)).toFixed(1)}%)`)
}

console.log('\nV4. THE INFORMATIVE FLOOR: CHROMOSOMES THAT ARE NOT BEING MEASURED')
console.log(`  ${'cell'.padEnd(12)} ${'chr'.padStart(3)} ${'loaded'.padEnd(7)} ${'nHom'.padStart(7)} `
  + `${'nHet'.padStart(6)} ${'frac of called'.padStart(15)} ${'hetHom'.padStart(7)} ${'hetHet'.padStart(7)}`)
const fr = obs.filter((o) => o.called > 0).map((o) => ({ o, f: (o.nHom + o.nHet) / o.called }))
const medFrac = q(fr.map((x) => x.f), 0.5)
let flagged = 0
for (const { o, f } of fr.sort((a, b) => a.f - b.f).slice(0, 8)) {
  flagged += 1
  console.log(`  ${o.gsm.padEnd(12)} ${o.chrom.padStart(3)} ${o.loaded.padEnd(7)} `
    + `${o.nHom.toLocaleString().padStart(7)} ${o.nHet.toLocaleString().padStart(6)} `
    + `${(100 * f).toFixed(1).padStart(14)}% `
    + `${(o.hetHom / Math.max(1, o.nHom)).toFixed(4).padStart(7)} `
    + `${(o.hetHet / Math.max(1, o.nHet)).toFixed(4).padStart(7)}`)
}
console.log(`  the median chromosome has ${(100 * medFrac).toFixed(1)}% of its called markers informative`)
const pass200 = fr.filter((x) => x.o.nHom >= FLOOR || x.o.nHet >= FLOOR)
const under = pass200.filter((x) => x.f < medFrac / 4)
console.log(`  chromosomes passing the absolute floor of ${FLOOR} whose informative FRACTION is `
  + `under a quarter of the median: ${under.length} of ${pass200.length}`)
void op

/**
 * The het-only rebuild of a missing parent from diploid siblings, built and validated.
 *
 * At a father-homozygous marker an embryo reading heterozygous PROVES the mother carries the other
 * allele: dropout removes alleles and cannot invent one. An embryo reading homozygous proves
 * nothing and is discarded, so a sample that lost her copy writes nothing in either direction.
 * Only heterozygous assertions are made about her; never homozygous ones.
 *
 * The validation family is public and carries truth: the egg donor is herself arrayed, so the
 * rebuild is scored against her own genotypes. k is derived from a contamination target rather
 * than chosen. Hold-out is at the EMBRYO, since siblings share events.
 *
 * Y1-Y4 fixed in audit/arm-b-and-rebuild-preregistration.txt.
 *
 * Run: OM_TRIOS=<dir> OM_MOTHER=<gsm> node --experimental-strip-types audit/het-only-rebuild.ts
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










const op = await import(`${W}oneParentOrigin.ts`)
const DROP_IN = op.DROP_IN as number
const H = 0.168  // the egg donor's heterozygosity on this platform

const MOTHER = process.env.OM_MOTHER
if (!MOTHER) { console.log('OM_MOTHER not set. NOT RUN.'); process.exit(0) }

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
let father = ''
const kids: string[] = []
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const m = at(r, 'mother_gsms')
  if (m !== MOTHER) continue
  const g = at(r, 'gsm')
  if (!fileOf.has(g)) continue
  father = father || at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  kids.push(g)
}
if (!father || kids.length < 2) { console.log('family not found. NOT RUN.'); process.exit(0) }

// The embryo each sample came from, read off GEO's own filename: the letter before the cell digit.
const embryoOf = new Map<string, string>()
for (const g of kids) {
  const f = (fileOf.get(g) ?? '').split('/').pop() ?? ''
  const m = /^GSM\d+_([A-Za-z]+)/.exec(f)
  embryoOf.set(g, m ? m[1].toUpperCase() : g)
}
const embryos = [...new Set(embryoOf.values())].sort()
console.log(`family: ${kids.length} samples over ${embryos.length} embryos, one father, `
  + 'the mother arrayed\n')

const gtOf = (g: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(fileOf.get(g)!)()) {
    if (!/^\d+$/.test(r.chrom)) continue
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const faGt = gtOf(father)
const moGt = gtOf(MOTHER)
console.log(`father ${faGt.size.toLocaleString()} called, mother ${moGt.size.toLocaleString()} called`)

// Every sample's genotype at the father-homozygous markers, plus its dropout.
const sample = new Map<string, { gt: Map<string, string>; ado: number }>()
for (const g of kids) {
  const gt = gtOf(g)
  let het = 0
  for (const v of gt.values()) if (v === 'AB') het += 1
  sample.set(g, { gt, ado: op.inferDropout(gt.size ? het / gt.size : NaN) as number })
  process.stdout.write('.')
}
console.log('\n')

// k derived from the contamination target the shipped rebuilt mother carries on its banner.
const TARGET = 0.02749
const lchoose = (n: number, k: number) => {
  let s = 0
  for (let i = 1; i <= k; i += 1) s += Math.log(n - k + i) - Math.log(i)
  return s
}
const tailGE = (n: number, p: number, k: number) => {
  let s = 0
  for (let i = k; i <= n; i += 1) s += Math.exp(lchoose(n, i) + i * Math.log(p) + (n - i) * Math.log(1 - p))
  return Math.min(1, Math.max(0, s))
}
function deriveK(n: number, ado: number) {
  const pSig = 0.5 * (1 - ado)
  for (let k = 1; k <= n; k += 1) {
    const d = tailGE(n, DROP_IN, k)
    const s = tailGE(n, pSig, k)
    const contam = (1 - H) * d / ((1 - H) * d + H * s)
    if (contam <= TARGET) return { k, contam, complete: s }
  }
  return { k: n, contam: NaN, complete: NaN }
}

/** The rebuild: assert the mother heterozygous where at least k of the given samples read AB. */
function rebuild(from: string[]) {
  const ados = from.map((g) => sample.get(g)!.ado).filter((x) => Number.isFinite(x))
  const ado = ados.length ? ados.reduce((a, b) => a + b, 0) / ados.length : 0.3
  const { k, contam, complete } = deriveK(from.length, ado)
  const votes = new Map<string, number>()
  for (const g of from) {
    const gt = sample.get(g)!.gt
    for (const [probe, fg] of faGt) {
      if (fg === 'AB') continue
      if (gt.get(probe) === 'AB') votes.set(probe, (votes.get(probe) ?? 0) + 1)
    }
  }
  const asserted = new Set<string>()
  for (const [probe, v] of votes) if (v >= k) asserted.add(probe)
  return { k, contam, complete, asserted, ado }
}

/** Score an assertion set against the mother's own array. */
function grade(asserted: Set<string>) {
  let hit = 0; let miss = 0
  for (const probe of asserted) {
    const t = moGt.get(probe)
    if (!t) continue
    if (t === 'AB') hit += 1; else miss += 1
  }
  // Completeness: of her true heterozygous markers at father-homozygous sites, how many found.
  let trueHet = 0; let found = 0
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    if (moGt.get(probe) !== 'AB') continue
    trueHet += 1
    if (asserted.has(probe)) found += 1
  }
  return { precision: hit + miss ? hit / (hit + miss) : NaN, scored: hit + miss,
    completeness: trueHet ? found / trueHet : NaN, trueHet }
}

console.log('Y1 AND Y2. ALL SAMPLES IN, derived k and the contamination it targets')
{
  const r = rebuild(kids)
  const s = grade(r.asserted)
  console.log(`  n ${kids.length}, mean dropout ${r.ado.toFixed(3)}, derived k ${r.k}`)
  console.log(`  modelled contamination ${r.contam.toFixed(5)}, so predicted precision `
    + `${(1 - r.contam).toFixed(4)}`)
  console.log(`  MEASURED precision   ${s.precision.toFixed(4)} over ${s.scored.toLocaleString()} asserted markers`)
  console.log(`  modelled completeness ${r.complete.toFixed(4)}`)
  console.log(`  MEASURED completeness ${s.completeness.toFixed(4)} of ${s.trueHet.toLocaleString()} true het markers`)
  console.log(`  Y1 ${s.precision >= 1 - r.contam ? 'HOLDS' : 'FAILS'}`)
}

console.log('\nY3. LEAVE ONE EMBRYO OUT, scored on the held-out embryo\'s own contribution')
console.log(`  ${'held out'.padEnd(10)} ${'n in'.padStart(5)} ${'k'.padStart(3)} `
  + `${'precision'.padStart(10)} ${'completeness'.padStart(13)} ${'asserted'.padStart(10)}`)
const precs: number[] = []
for (const e of embryos) {
  const from = kids.filter((g) => embryoOf.get(g) !== e)
  if (from.length < 2) continue
  const r = rebuild(from)
  const s = grade(r.asserted)
  precs.push(s.precision)
  console.log(`  ${e.padEnd(10)} ${String(from.length).padStart(5)} ${String(r.k).padStart(3)} `
    + `${s.precision.toFixed(4).padStart(10)} ${s.completeness.toFixed(4).padStart(13)} `
    + `${r.asserted.size.toLocaleString().padStart(10)}`)
}
if (precs.length) {
  const mean = precs.reduce((a, b) => a + b, 0) / precs.length
  console.log(`  mean precision over ${precs.length} hold-outs: ${mean.toFixed(4)}`)
}

console.log('\nY4. DOES PRECISION FALL AS MORE SAMPLES ENTER? A majority vote would.')
console.log(`  ${'n in'.padStart(5)} ${'k'.padStart(3)} ${'precision'.padStart(10)} `
  + `${'completeness'.padStart(13)} ${'asserted'.padStart(10)}`)
for (const n of [4, 8, 12, 16, 20, kids.length].filter((x, i, a) => x <= kids.length && a.indexOf(x) === i)) {
  const r = rebuild(kids.slice(0, n))
  const s = grade(r.asserted)
  console.log(`  ${String(n).padStart(5)} ${String(r.k).padStart(3)} ${s.precision.toFixed(4).padStart(10)} `
    + `${s.completeness.toFixed(4).padStart(13)} ${r.asserted.size.toLocaleString().padStart(10)}`)
}

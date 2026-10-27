/**
 * Why the parent-heterozygous class cannot replace the shipped one-parent call.
 *
 * It scores better on the published karyotype because that corpus's one-contribution cases are
 * monosomies: one parent's copy gone, so the sample carries a single haplotype and reads near
 * drop-in in both classes. A GYNOGENETIC sample is different. It carries TWO maternal haplotypes,
 * which are heterozygous wherever a crossover separated them, so it reads high on the
 * parent-heterozygous class and looks biparental there while carrying no paternal contribution at
 * all. An external review warned of exactly this in words, and `crossSurface.check.ts` caught it
 * when the class was substituted.
 *
 * This measures it on real arrays rather than on a fixture: samples of a series with a measured
 * egg donor and no father, where the two carrying two maternal haplotypes are known from the
 * marker channels.
 *
 * No laboratory sample identifier is printed; samples are reported by position in the series.
 *
 * Run: OM_TRIOS=<dir> OM_REF=<reference> node --experimental-strip-types audit/gynogenetic-blind.ts
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
const REF = process.env.OM_REF
if (!REF || !existsSync(join(DIR, REF))) { console.log('OM_REF not readable. NOT RUN.'); process.exit(0) }

const refGt = new Map<string, string>()
for (const r of rowsOf(join(DIR, REF))()) {
  if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') refGt.set(r.probesetId, r.genotype)
}
const prefix = /^(\d{4,6})-/.exec(REF)?.[1] ?? ''
const samples = readdirSync(DIR)
  .filter((f) => f.startsWith(`${prefix}-`) && /\.probes/.test(f) && f !== REF).sort()

console.log(`reference ${refGt.size.toLocaleString()} called markers\n`)
console.log(`${'pos'.padEnd(5)} ${'oppHom'.padStart(8)} ${'het'.padStart(7)} ${'nHom'.padStart(8)} `
  + `${'hetHom'.padStart(7)} ${'nHet'.padStart(8)} ${'hetHet'.padStart(7)} `
  + `${'HOM says'.padEnd(13)} ${'HET says'.padEnd(13)}  truth from the channels`)
for (const f of samples) {
  const pos = /-(\d{2})\.probes/.exec(f)?.[1] ?? '??'
  let nHom = 0; let hHom = 0; let nHet = 0; let hHet = 0
  let opp = 0; let oppN = 0; let called = 0; let het = 0; let rows = 0
  for (const r of rowsOf(join(DIR, f))()) {
    rows += 1
    const s = r.genotype
    if (s !== 'AA' && s !== 'BB' && s !== 'AB') continue
    called += 1
    if (s === 'AB') het += 1
    if (!/^\d+$/.test(r.chrom)) continue
    const rg = refGt.get(r.probesetId)
    if (!rg) continue
    if (rg === 'AB') { nHet += 1; if (s === 'AB') hHet += 1 } else {
      nHom += 1; if (s === 'AB') hHom += 1
      if (s !== 'AB') { oppN += 1; if (s !== rg) opp += 1 }
    }
  }
  if (called / rows < 0.60) continue
  const hetRate = called ? het / called : NaN
  const ado = op.inferDropout(hetRate) as number
  const homCall = obl.hetCall({ informative: nHom, het: hHom, onlyA: 0, onlyB: 0 } as never, 1) as { ploidy: string }
  const hetCall2 = obl.hetCallOneParent(hHet, nHet, ado) as { ploidy: string }
  const oppRate = oppN ? opp / oppN : NaN
  // A product of the reference reads near zero; a sample carrying none of her genome reads high.
  const truth = oppRate < 0.01 ? 'carries her genome only'
    : oppRate > 0.04 ? 'carries no genome of hers' : 'intermediate'
  console.log(`${pos.padEnd(5)} ${(100 * oppRate).toFixed(2).padStart(7)}% ${hetRate.toFixed(3).padStart(7)} `
    + `${nHom.toLocaleString().padStart(8)} ${(hHom / Math.max(1, nHom)).toFixed(4).padStart(7)} `
    + `${nHet.toLocaleString().padStart(8)} ${(hHet / Math.max(1, nHet)).toFixed(4).padStart(7)} `
    + `${homCall.ploidy.padEnd(13)} ${hetCall2.ploidy.padEnd(13)}  ${truth}`)
}
console.log('\nA sample carrying her genome only has NO second contribution. Where the HET column')
console.log('reads biparental on such a sample, that class cannot replace the shipped one.')

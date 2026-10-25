/**
 * The four measurements of audit/refusal-split-preregistration.txt.
 *
 * H: the MZ pronuclei's parents, from the opposite-homozygote rate rather than from the paper.
 * I: the segment scan's false-segment baseline on clean haploid products of the sperm donor.
 * J: the same on POOR products, out of sample, which is the specificity cost of admitting them.
 * K: whether per-chromosome heterozygosity spread separates a chaotic genome from a mixture.
 *
 * Scored with ONE parent, the sperm donor. With the mother loaded every autosome of a paternal
 * pronucleus reads as having lost her, which puts every chromosome in `lostByAlleles`, and the
 * scan skips those by design: the false-positive count would be zero by construction.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/refusal-split-baseline.ts
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



type Res = {
  segments?: { chrom: string; kind: string; startBp: number; endBp: number; score: number
    absentParent?: string }[]
  stage?: { stage: string; why: string }
  chroms?: { chrom: string; verdict: string; lostParent?: string }[]
}

async function scoreSolo(gsm: string, parent: string, declared?: string) {
  const pat = parentIndex(parent)
  const acc = score.emptyCollected(pat as never, null as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, null as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat: null, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
    ...(declared ? { declaredStage: declared } : {}),
  } as never) as Res
  return { res, profile: profile as { callRate: number; hetRate: number
    byChrom: Map<string, { called: number; het: number }> } }
}

// The bulk parents and the MZ pronuclei, by GEO's own supplementary filenames.
const named = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const g = /^(GSM\d+)/.exec(f)?.[1]
  if (!g) continue
  const m = /_(sperm_DNA|donor_[A-D]|[pm]MII-\d+)_/.exec(f)
  if (m) named.set(m[1], g)
}
const sperm = named.get('sperm_DNA')!
const donors = [...named.keys()].filter((k) => k.startsWith('donor_')).sort()
const pn = [...named.keys()].filter((k) => /MII-/.test(k)).sort()

// --- H ----------------------------------------------------------------------------------------
// The opposite-homozygote rate, computed directly: markers where the parent is homozygous and the
// sample is homozygous for the other allele. A haploid product of that parent reads 0.2 to 0.8%.
const gtOf = (gsm: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const oppRate = (s: Map<string, string>, p: Map<string, string>) => {
  let n = 0; let opp = 0
  for (const [id, pg] of p) {
    if (pg === 'AB') continue
    const sg = s.get(id); if (!sg || sg === 'AB') continue
    n += 1; if (sg !== pg) opp += 1
  }
  return n ? opp / n : NaN
}
console.log('TEST H. THE PRONUCLEI AGAINST EVERY BULK PARENT, OPPOSITE-HOMOZYGOTE RATE')
const parentGt = new Map<string, Map<string, string>>()
for (const k of ['sperm_DNA', ...donors]) parentGt.set(k, gtOf(named.get(k)!))
console.log(`${'pronucleus'.padEnd(11)} ` + ['sperm_DNA', ...donors].map((k) => k.padStart(11)).join(' '))
const pairing = new Map<string, string>()
for (const p of pn) {
  const s = gtOf(named.get(p)!)
  const cells: string[] = []; let best = ''; let bestR = 1
  for (const k of ['sperm_DNA', ...donors]) {
    const r = oppRate(s, parentGt.get(k)!)
    cells.push(`${(r * 100).toFixed(2)}%`.padStart(11))
    if (r < bestR) { bestR = r; best = k }
  }
  if (bestR < 0.01) pairing.set(p, best)
  console.log(`${p.padEnd(11)} ${cells.join(' ')}   -> ${bestR < 0.01 ? best : 'no parent'}`)
}
console.log(`\n  pronuclei with exactly one parent under 1%: ${pairing.size} of ${pn.length}`)
const eggs = new Set([...pairing.values()].filter((v) => v.startsWith('donor_')))
console.log(`  egg donors identified: ${[...eggs].join(', ') || 'none'}`)

// --- I and J ----------------------------------------------------------------------------------
// Every parental-absence segment on a haploid product of the loaded parent is false: his alleles
// are present on every autosome by dissection.
const CLEAN = ['pMII-1', 'pMII-2', 'pMII-3', 'pMII-4', 'pMII-6']
const POOR = ['pMII-5', 'pMII-7', 'pMII-8']
for (const [label, set, declared] of [
  ['TEST I. BASELINE, THE FIVE CLEAN PATERNAL PRONUCLEI, ONE PARENT LOADED', CLEAN, undefined],
  ['TEST J. POOR PATERNAL PRONUCLEI, OUT OF SAMPLE, STAGE DECLARED', POOR, 'single-cell'],
] as [string, string[], string | undefined][]) {
  console.log(`\n${label}`)
  console.log(`${'array'.padEnd(9)} ${'stage'.padEnd(13)} ${'callRate'.padStart(8)} `
    + `${'het'.padStart(6)} ${'scanned'.padStart(7)} ${'FALSE segs'.padStart(10)}  where`)
  let segs = 0; let scanned = 0
  for (const p of set) {
    const { res, profile } = await scoreSolo(named.get(p)!, sperm, declared)
    const auto = (res.chroms ?? []).filter((c) => /^\d+$/.test(c.chrom))
    const n = auto.filter((c) => c.verdict !== 'not_measured').length
    const bad = (res.segments ?? []).filter((s) => s.kind === 'parental-absence'
      && /^\d+$/.test(s.chrom))
    segs += bad.length; scanned += n
    console.log(`${p.padEnd(9)} ${(res.stage?.stage ?? '?').padEnd(13)} `
      + `${profile.callRate.toFixed(2).padStart(8)} ${profile.hetRate.toFixed(3).padStart(6)} `
      + `${String(n).padStart(7)} ${String(bad.length).padStart(10)}  `
      + bad.map((s) => `chr${s.chrom} ${(s.startBp / 1e6).toFixed(0)}-${(s.endBp / 1e6).toFixed(0)}Mb LRT ${s.score.toFixed(0)}`).join('; '))
  }
  console.log(`  FALSE parental-absence segments: ${segs} over ${scanned} chromosome scans`)
}

// --- K ----------------------------------------------------------------------------------------
console.log('\nTEST K. PER-CHROMOSOME HETEROZYGOSITY SPREAD, EVERY CELL WITH A HAND KARYOTYPE')
const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const parentsOf = new Map<string, string>()
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const f = at(r, 'father_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  if (f) parentsOf.set(at(r, 'gsm'), f)
}
const hand = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , karyo] = l.split('\t')
  if (!hand.has(gsm)) hand.set(gsm, karyo)
}
const cells = [...hand.keys()].filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()
const out: { gsm: string; stage: string; cr: number; het: number; cv: number; k: string }[] = []
for (const g of cells) {
  const { res, profile } = await scoreSolo(g, parentsOf.get(g)!)
  const rates = [...profile.byChrom.entries()].filter(([c]) => /^\d+$/.test(c))
    .filter(([, s]) => s.called > 500).map(([, s]) => s.het / s.called)
  const mean = rates.reduce((a, b) => a + b, 0) / rates.length
  const sd = Math.sqrt(rates.reduce((a, b) => a + (b - mean) ** 2, 0) / rates.length)
  out.push({ gsm: g, stage: res.stage?.stage ?? '?', cr: profile.callRate,
    het: profile.hetRate, cv: sd / mean, k: hand.get(g) ?? '' })
}
console.log(`${'cell'.padEnd(12)} ${'stage'.padEnd(14)} ${'callRate'.padStart(8)} `
  + `${'het'.padStart(6)} ${'CV'.padStart(6)}  hand karyotype`)
for (const r of out.sort((a, b) => b.cv - a.cv)) {
  console.log(`${r.gsm.padEnd(12)} ${r.stage.padEnd(14)} ${r.cr.toFixed(2).padStart(8)} `
    + `${r.het.toFixed(3).padStart(6)} ${r.cv.toFixed(3).padStart(6)}  ${r.k.slice(0, 32)}`)
}
const refused = out.filter((r) => r.stage === 'failed')
console.log(`\n  cells the stage inference refuses: ${refused.length}`)
for (const r of refused) {
  console.log(`    ${r.gsm}  CV ${r.cv.toFixed(3)}  het ${r.het.toFixed(3)}  `
    + `callRate ${r.cr.toFixed(2)}  hand: ${r.k.slice(0, 30)}`)
}

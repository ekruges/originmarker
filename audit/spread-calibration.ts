/**
 * Calibration of the conditioned heterozygosity spread, on four real classes.
 *
 *   uniform diploid      the 5 bulk parental arrays. Truth is structural.
 *   uniform hemizygous   the 15 MII pronuclei. Truth is the aspiration.
 *   localised hemizygous the 31 hand-karyotyped cells, graded 0 to 12 monosomies. Truth is the
 *                        authors' by-eye karyotype.
 *   normal, poor array   a cell the hand reads 46,XX at a 0.67 call rate, refused by the ceiling.
 *
 * Predictions Q1 to Q6 are fixed in audit/refusal-split-preregistration.txt.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/spread-calibration.ts
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




type Stats = { markers: number; called: number; het: number }
type Res = { stage?: { stage: string }; chroms?: { chrom: string; verdict: string }[] }

const stg = await import(`${W}stage.ts`)
const LOC = stg.HET_SPREAD_LOCALISED as number

// The statistic as it stood, unconditioned, so the two read side by side.
const oldSpread = (b: Map<string, Stats>) => {
  const r = [...b.entries()].filter(([c, v]) => /^\d+$/.test(c) && v.called >= 500)
    .map(([, v]) => v.het / v.called)
  if (r.length < 5) return NaN
  const m = r.reduce((a, x) => a + x, 0) / r.length
  return Math.sqrt(r.reduce((a, x) => a + (x - m) ** 2, 0) / r.length) / m
}

async function profileOf(gsm: string, father?: string, mother?: string) {
  const pat = father ? parentIndex(father) : null
  const mat = mother ? parentIndex(mother) : null
  const acc = pat ? score.emptyCollected(pat as never, mat as never) : null
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    if (acc && pat) score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  let res: Res = {}
  if (acc && pat) {
    res = await score.scoreSample({
      acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
    } as never) as Res
  }
  const p = profile as unknown as { callRate: number; hetRate: number; byChrom: Map<string, Stats> }
  return { p, res, now: stg.hetSpread(p.byChrom as never) as number, was: oldSpread(p.byChrom) }
}

// The bulk parents and the pronuclei, by GEO's own supplementary filenames.
const named = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const g = /^(GSM\d+)/.exec(f)?.[1]; if (!g) continue
  const m = /_(sperm_DNA|donor_[A-D]|[pm]MII-\d+)_/.exec(f)
  if (m) named.set(m[1], g)
}
const sperm = named.get('sperm_DNA')!

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

// The graded truth: hand-called monosomies and nullisomies per cell.
const mono = new Map<string, number>(); const nulli = new Map<string, number>()
const norm = new Map<string, number>(); const karyo = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , k, , manual] = l.split('\t')
  if (!karyo.has(gsm)) karyo.set(gsm, k)
  const m = manual.toLowerCase()
  const bump = (x: Map<string, number>) => x.set(gsm, (x.get(gsm) ?? 0) + 1)
  if (m.includes('nullisom')) bump(nulli)
  else if (m.includes('monosom')) bump(mono)
  else if (m.startsWith('normal')) bump(norm)
}

const line = (label: string, gsm: string, r: Awaited<ReturnType<typeof profileOf>>, extra = '') => {
  const auto = (r.res.chroms ?? []).filter((c) => /^\d+$/.test(c.chrom))
  const nm = auto.filter((c) => c.verdict !== 'not_measured').length
  console.log(`${label.padEnd(22)} ${gsm.padEnd(12)} ${r.p.callRate.toFixed(2).padStart(8)} `
    + `${r.p.hetRate.toFixed(3).padStart(6)} `
    + `${(Number.isFinite(r.was) ? r.was.toFixed(3) : 'NaN').padStart(6)} `
    + `${(Number.isFinite(r.now) ? r.now.toFixed(3) : 'NaN').padStart(6)} `
    + `${(Number.isFinite(r.now) && r.now >= LOC ? 'LOCALISED' : 'uniform').padEnd(10)} `
    + `${String(nm || '-').padStart(3)}  ${extra}`)
}
const HEAD = `${'class'.padEnd(22)} ${'array'.padEnd(12)} ${'callRate'.padStart(8)} `
  + `${'het'.padStart(6)} ${'was'.padStart(6)} ${'now'.padStart(6)} ${'verdict'.padEnd(10)} `
  + `${'nMe'.padStart(3)}  truth`

console.log(`threshold HET_SPREAD_LOCALISED = ${LOC}\n`)
console.log('Q1 AND Q2. UNIFORM DIPLOID AND UNIFORM HEMIZYGOUS, TRUTH NOT FROM ANY ARRAY READING')
console.log(HEAD)
let q1 = 0; let q2 = 0
for (const k of ['sperm_DNA', 'donor_A', 'donor_B', 'donor_C', 'donor_D']) {
  if (!named.has(k)) continue
  const r = await profileOf(named.get(k)!)
  if (Number.isFinite(r.now) && r.now >= LOC) q1 += 1
  line('uniform diploid', named.get(k)!, r, `bulk ${k}`)
}
for (const k of [...named.keys()].filter((x) => /MII-/.test(x)).sort()) {
  const r = await profileOf(named.get(k)!, sperm)
  if (Number.isFinite(r.now) && r.now >= LOC) q2 += 1
  line('uniform hemizygous', named.get(k)!, r, `${k}, one copy of every autosome`)
}
console.log(`\n  Q1 bulk diploid arrays reading LOCALISED, must be 0:   ${q1}`)
console.log(`  Q2 pronuclei reading LOCALISED, must be 0:              ${q2}`)

console.log('\nQ3 TO Q6. THE HAND-KARYOTYPED CELLS, GRADED BY ABNORMAL AUTOSOME COUNT')
console.log(HEAD)
const pts: { gsm: string; now: number; cr: number; share: number; nm: number }[] = []
const cells = [...karyo.keys()].filter((g) => parentsOf.has(g) && fileOf.has(g))
for (const g of cells) {
  const p = parentsOf.get(g)!
  const r = await profileOf(g, p.father, p.mother)
  const auto = (r.res.chroms ?? []).filter((c) => /^\d+$/.test(c.chrom))
  const nm = auto.filter((c) => c.verdict !== 'not_measured').length
  const ab = (mono.get(g) ?? 0) + (nulli.get(g) ?? 0)
  if (Number.isFinite(r.now) && nm > 0) {
    pts.push({ gsm: g, now: r.now, cr: r.p.callRate, share: ab / 22, nm })
  }
  pts.sort((a, b) => b.now - a.now)
  line(r.res.stage?.stage ?? '?', g, r,
    `mono ${mono.get(g) ?? 0}, null ${nulli.get(g) ?? 0}, normal ${norm.get(g) ?? 0}  `
    + `[${(karyo.get(g) ?? '').slice(0, 24)}]`)
}
const corr = (xs: number[], ys: number[]) => {
  const n = xs.length; const mx = xs.reduce((a, b) => a + b, 0) / n
  const my = ys.reduce((a, b) => a + b, 0) / n
  let sxy = 0; let sx = 0; let sy = 0
  for (let i = 0; i < n; i += 1) {
    sxy += (xs[i] - mx) * (ys[i] - my); sx += (xs[i] - mx) ** 2; sy += (ys[i] - my) ** 2
  }
  return sxy / Math.sqrt(sx * sy)
}
console.log('\n  Q6. WHAT THE CONDITIONED SPREAD TRACKS, over the cells where it is measurable')
console.log(`    n = ${pts.length}`)
console.log(`    against the hand-called abnormal autosome share: r = ${corr(pts.map((p) => p.share), pts.map((p) => p.now)).toFixed(3)}`)
console.log(`    against the call rate:                           r = ${corr(pts.map((p) => p.cr), pts.map((p) => p.now)).toFixed(3)}`)
console.log('    the first must be the stronger association, and positive, or the statistic is')
console.log('    still reading dropout rather than hemizygosity.')
for (const [q, g, want] of [['Q3', 'GSM4472421', 'uniform'], ['Q4', 'GSM4472409', 'uniform'],
  ['Q5', 'GSM4472405', 'LOCALISED']] as [string, string, string][]) {
  const p = pts.find((x) => x.gsm === g)
  const got = p ? (p.now >= LOC ? 'LOCALISED' : 'uniform') : 'not measurable'
  console.log(`    ${q} ${g}: predicted ${want.padEnd(9)} got ${got.padEnd(14)} `
    + `${got === want ? 'HELD' : 'FAILED'}`)
}

/**
 * Step 1 of audit/null-replacement-preregistration.txt: the SHAPE of the per-chromosome absence
 * distribution, measured against the hand karyotype, before any estimator is named.
 *
 * T1 bimodal, lower mode near the floor on hand-intact chromosomes and upper near the measured
 *    event rate of 0.115 to 0.124 on hand-lost ones.
 * T2 on the cells where the median fails, the hand-normal chromosomes still sit at the lower mode.
 * T3 how many autosomes the hand reads as normal, per cell.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/null-shape.ts
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





type Chrom = { chrom: string; verdict: string }
type Res = { chroms?: Chrom[]; stage?: { stage: string } }

const segs2 = await import(`${W}segments.ts`)
const obl = await import(`${W}obligateHet.ts`)

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
  const a = acc as unknown as {
    t: { byChrom: Map<string, [number, number]> }
    matAbsenceTally: Map<string, [number, number]>
  }
  const lossRate = Math.min(
    obl.lossSignalRate((pat as never as { gt: Map<string, string> }).gt,
      (mat as never as { gt: Map<string, string> }).gt) as number,
    obl.lossSignalRate((mat as never as { gt: Map<string, string> }).gt,
      (pat as never as { gt: Map<string, string> }).gt) as number,
  )
  return { res, acc: a, lossRate, callRate: (profile as unknown as { callRate: number }).callRate }
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
// The hand reading of each chromosome, and whether it says this parent's copy went.
const hand = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , , chrom, manual] = l.split('\t')
  hand.set(`${gsm}|${chrom}`, manual)
}
const cells = [...new Set([...hand.keys()].map((k) => k.split('|')[0]))]
  .filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()

type Row = { gsm: string; chrom: string; run: 'pat' | 'mat'; rate: number; n: number
  handLost: boolean; handNormal: boolean; whose: string }
const rows: Row[] = []
const perCell: { gsm: string; stage: string; cr: number; normals: number; abnormal: number
  medPat: number; medMat: number; loss: number }[] = []

for (const g of cells) {
  const p = parentsOf.get(g)!
  const { res, acc, lossRate, callRate } = await run(g, p.father, p.mother)
  let normals = 0; let abnormal = 0
  for (const c of (res.chroms ?? [])) {
    if (!/^\d+$/.test(c.chrom)) continue
    const h = hand.get(`${g}|${c.chrom}`) ?? ''
    const isNormal = /^normal/i.test(h.trim())
    const lostPat = /\(\s*pat/i.test(h) && /del|monosom|nullisom/i.test(h)
    const lostMat = /\(\s*mat/i.test(h) && /del|monosom|nullisom/i.test(h)
    const anyLoss = /del|monosom|nullisom/i.test(h)
    if (isNormal) normals += 1; else if (h) abnormal += 1
    for (const [run, tally] of [['pat', acc.t.byChrom], ['mat', acc.matAbsenceTally]] as const) {
      const e = tally.get(c.chrom)
      if (!e || e[0] < 200) continue
      rows.push({
        gsm: g, chrom: c.chrom, run, rate: e[1] / e[0], n: e[0],
        // A chromosome that lost THIS parent's copy should read high on THIS run.
        handLost: run === 'pat' ? lostPat : lostMat,
        handNormal: isNormal, whose: anyLoss ? (lostPat ? 'pat' : lostMat ? 'mat' : '?') : '',
      })
    }
  }
  const med = (t: Map<string, [number, number]>) => {
    const xs = [...t.entries()].filter(([c, e]) => /^\d+$/.test(c) && e[0] >= 200)
      .map(([, e]) => e[1] / e[0]).sort((a, b) => a - b)
    return xs.length ? xs[xs.length >> 1] : NaN
  }
  perCell.push({ gsm: g, stage: res.stage?.stage ?? '?', cr: callRate, normals, abnormal,
    medPat: med(acc.t.byChrom), medMat: med(acc.matAbsenceTally), loss: lossRate })
  process.stdout.write('.')
}
console.log('\n')

const q = (xs: number[], p: number) => {
  if (!xs.length) return NaN
  const s = xs.slice().sort((a, b) => a - b)
  return s[Math.min(s.length - 1, Math.max(0, Math.round(p * (s.length - 1))))]
}
const show = (label: string, xs: number[]) => {
  if (!xs.length) { console.log(`  ${label.padEnd(34)} none`); return }
  console.log(`  ${label.padEnd(34)} n=${String(xs.length).padStart(4)}  `
    + `p10 ${q(xs, 0.10).toFixed(4)}  median ${q(xs, 0.50).toFixed(4)}  `
    + `p90 ${q(xs, 0.90).toFixed(4)}`)
}

console.log('T1. THE SHAPE. Per-chromosome absence rate, split by what the hand karyotype says.')
console.log(`  the measured rate a true single-parent loss produces: `
  + `${(perCell[0]?.loss ?? NaN).toFixed(4)}`)
show('hand reads this parent LOST', rows.filter((r) => r.handLost).map((r) => r.rate))
show('hand reads the chromosome NORMAL', rows.filter((r) => r.handNormal).map((r) => r.rate))
show('hand reads a loss of the OTHER parent', rows.filter((r) => !r.handLost && r.whose && r.whose !== r.run).map((r) => r.rate))
show('everything else', rows.filter((r) => !r.handLost && !r.handNormal).map((r) => r.rate))

console.log('\nT2 AND T3. PER CELL: DOES A LOW QUANTILE STILL FIND THE CLEAN CHROMOSOMES?')
console.log(`${'cell'.padEnd(12)} ${'stage'.padEnd(13)} ${'norm'.padStart(4)} ${'abn'.padStart(4)} `
  + `${'medPat'.padStart(7)} ${'medMat'.padStart(7)} ${'p10 pat'.padStart(8)} ${'p10 mat'.padStart(8)} `
  + `${'normOnly'.padStart(9)}  median blind?`)
const BOUND = 0.0253
for (const c of perCell.sort((a, b) => b.abnormal - a.abnormal)) {
  const mine = rows.filter((r) => r.gsm === c.gsm)
  const p10p = q(mine.filter((r) => r.run === 'pat').map((r) => r.rate), 0.10)
  const p10m = q(mine.filter((r) => r.run === 'mat').map((r) => r.rate), 0.10)
  // What the null WOULD be if it were taken only over the chromosomes the hand reads as normal.
  const normOnly = q(mine.filter((r) => r.handNormal).map((r) => r.rate), 0.50)
  const blind = !segs2.canDetectLoss(Math.max(c.medPat, c.medMat), c.loss)
  console.log(`${c.gsm.padEnd(12)} ${c.stage.padEnd(13)} ${String(c.normals).padStart(4)} `
    + `${String(c.abnormal).padStart(4)} ${c.medPat.toFixed(4).padStart(7)} `
    + `${c.medMat.toFixed(4).padStart(7)} ${p10p.toFixed(4).padStart(8)} `
    + `${p10m.toFixed(4).padStart(8)} `
    + `${(Number.isFinite(normOnly) ? normOnly.toFixed(4) : '-').padStart(9)}  `
    + `${blind ? 'BLIND' : 'ok'}`)
}
console.log(`\n  the blindness boundary is a null of ${BOUND}`)
console.log('  normOnly is the median over only the chromosomes the hand reads as normal: the')
console.log('  best a within-sample null could do if it knew which they were.')

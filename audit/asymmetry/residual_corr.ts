/**
 * Residual correlation between the BAF and intensity channels, measured on the SHIPPED statistic.
 *
 * A joint term needs to know how far the two channels are already saying the same thing. The
 * methods review measured this on technical-replicate differences (r = -0.055 bulk, -0.007 TE,
 * +0.036 blastomere) and on replicate means (up to -0.46), and those disagree because the means
 * carry a common amplification artefact the differences cancel.
 *
 * Neither is the number this code needs. The shipped statistic is SELF-REFERENCED: each
 * chromosome's centroid has the same array's other chromosomes subtracted, which removes exactly
 * the common-array artefact that drives the -0.46. So the residual correlation on the shipped
 * statistic has to be measured on the shipped statistic, per array, across chromosomes, on
 * material with no event.
 *
 *   usage: node --experimental-strip-types audit/asymmetry/residual_corr.ts <ref> <dir> [out.json]
 */
import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs'
import { join, basename } from 'node:path'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const { isAutosome } = await import(`${W}parentage.ts`)
const { inferStage } = await import(`${W}stage.ts`)
const dosage = await import(`${W}dosageOrigin.ts`)

const REF = process.argv[2]
const DIR = process.argv[3]
const OUT = process.argv[4] ?? 'audit/asymmetry/residual-corr.json'
const STRIDE = Number(process.argv[5] ?? 8)
if (!REF || !DIR) throw new Error('usage: residual_corr.ts <ref> <dir> [out.json]')

type AB = 'AA' | 'AB' | 'BB' | 'NC'
const readLines = (p: string) => {
  const raw = readFileSync(p)
  return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8').split('\n')
}
function header(lines: string[]) {
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  return { h, map: h < 0 ? null : ingest.headerMap(lines[h]) }
}

// The reference genotypes, once.
const refLines = readLines(REF)
const rh = header(refLines)
if (!rh.map) throw new Error('reference header not recognised')
const refGt = new Map<string, AB>()
for (let i = rh.h + 1; i < refLines.length; i += 1) {
  const r = ingest.parseRow(refLines[i], rh.map)
  if (r && r.genotype !== 'NC' && isAutosome(r.chrom)) refGt.set(r.probesetId, r.genotype as AB)
}
process.stderr.write(`reference ${basename(REF)}: ${refGt.size} autosomal calls\n`)

function files(dir: string, rel = ''): string[] {
  const out: string[] = []
  for (const e of readdirSync(join(dir, rel)).sort()) {
    const r = rel ? join(rel, e) : e
    if (statSync(join(dir, r)).isDirectory()) out.push(...files(dir, r))
    else if (/\.(probes|txt\.gz)$/.test(e)) out.push(r)
  }
  return out
}

const pearson = (xs: number[], ys: number[]): number => {
  const n = xs.length
  if (n < 3) return NaN
  const mx = xs.reduce((a, b) => a + b, 0) / n
  const my = ys.reduce((a, b) => a + b, 0) / n
  let sxy = 0
  let sxx = 0
  let syy = 0
  for (let i = 0; i < n; i += 1) {
    const dx = xs[i] - mx
    const dy = ys[i] - my
    sxy += dx * dy; sxx += dx * dx; syy += dy * dy
  }
  return sxx && syy ? sxy / Math.sqrt(sxx * syy) : NaN
}
const median = (xs: number[]) => {
  if (!xs.length) return NaN
  const q = [...xs].sort((a, b) => a - b)
  return q[q.length >> 1]
}

const rows: { id: string, material: string, chroms: number, r: number }[] = []
for (const f of files(DIR)) {
  const lines = readLines(join(DIR, f))
  const { h, map } = header(lines)
  if (!map) continue

  // Per chromosome: the central-window BAF readings, and the log2R readings.
  const baf = new Map<string, number[]>()
  const lrr = new Map<string, number[]>()
  const byChrom = new Map()
  const sums = ingest.emptyBafSums()
  let first = ''
  // Stride, and it costs nothing here. This measures a correlation across ~22 per-chromosome
  // means, and the review's own figure is that effective independent markers saturate near 250
  // per chromosome out of 900-1,100, so a full read buys no precision and costs eight times as
  // much. At stride 8 each chromosome still carries thousands of readings.
  for (let i = h + 1; i < lines.length; i += STRIDE) {
    const r = ingest.parseRow(lines[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, sums as never)
    if (!isAutosome(r.chrom)) continue
    if (r.log2R !== null && Number.isFinite(r.log2R)) {
      lrr.set(r.chrom, [...(lrr.get(r.chrom) ?? []), r.log2R])
    }
    const pg = refGt.get(r.probesetId)
    if (!pg || (pg !== 'AA' && pg !== 'BB')) continue
    if (r.baf === null || !Number.isFinite(r.baf)) continue
    const b = dosage.oriented(pg as never, r.baf)
    if (b < dosage.WINDOW_LO || b > dosage.WINDOW_HI) continue
    baf.set(r.chrom, [...(baf.get(r.chrom) ?? []), b])
  }
  const id = basename(f).replace(/\.(probes|CEL\.txt\.gz)$/, '').replace(/^(GSM\d+)_.*$/, '$1')
  const prof = ingest.finishProfile(id, byChrom as never, sums as never, first)
  const st = inferStage(prof)
  if (st.stage === 'failed') continue
  const material = dosage.materialOf(st.stage)

  // SELF-REFERENCED both channels: each chromosome against the array's other chromosomes, which
  // is what the shipped statistic does and what removes the common-array artefact.
  const chroms = [...baf.keys()].filter((c) => (baf.get(c)?.length ?? 0) >= 200
    && (lrr.get(c)?.length ?? 0) >= 200)
  // Leave-one-chromosome-out by SUMS rather than by rebuilding the background each time. The
  // concatenating version is quadratic in chromosomes times markers and took minutes per array.
  let totB = 0
  let nB = 0
  let totL = 0
  let nL = 0
  for (const c of chroms) {
    for (const v of baf.get(c) ?? []) { totB += v; nB += 1 }
    for (const v of lrr.get(c) ?? []) { totL += v; nL += 1 }
  }
  const xs: number[] = []
  const ys: number[] = []
  for (const c of chroms) {
    const B = baf.get(c) ?? []
    const L = lrr.get(c) ?? []
    const sB = B.reduce((a, b) => a + b, 0)
    const sL = L.reduce((a, b) => a + b, 0)
    const inB = sB / B.length
    const inL = sL / L.length
    const outB = (totB - sB) / (nB - B.length)
    const outL = (totL - sL) / (nL - L.length)
    if (Number.isFinite(inB - outB) && Number.isFinite(inL - outL)) {
      xs.push(inB - outB)
      ys.push(inL - outL)
    }
  }
  const r = pearson(xs, ys)
  if (Number.isFinite(r)) {
    rows.push({ id, material, chroms: xs.length, r })
    process.stderr.write(`${id.padEnd(14)} ${material.padEnd(14)} ${xs.length} chrom  r ${r.toFixed(3)}\n`)
  }
}

const byMaterial = new Map<string, number[]>()
for (const x of rows) byMaterial.set(x.material, [...(byMaterial.get(x.material) ?? []), x.r])
process.stderr.write('\nresidual correlation of the SHIPPED self-referenced statistic\n')
const summary: Record<string, { arrays: number, medianR: number, meanR: number }> = {}
for (const [m, rs] of byMaterial) {
  const mean = rs.reduce((a, b) => a + b, 0) / rs.length
  summary[m] = { arrays: rs.length, medianR: median(rs), meanR: mean }
  process.stderr.write(`  ${m.padEnd(14)} n=${String(rs.length).padStart(3)}  `
    + `median r ${median(rs).toFixed(3)}  mean r ${mean.toFixed(3)}\n`)
}
writeFileSync(OUT, JSON.stringify({ reference: basename(REF), dir: DIR, summary, rows }, null, 1))
process.stderr.write(`\n-> ${OUT}\n`)

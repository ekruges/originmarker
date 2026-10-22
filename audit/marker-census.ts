/**
 * How many markers on an array are informative, and how many per chromosome.
 *
 * A marker is informative for parent of origin when the two parents are homozygous for OPPOSITE
 * alleles, because only then does each parental copy carry a different letter and the share between
 * them can be read. A marker where either parent is heterozygous, or where both are homozygous for
 * the same allele, says nothing about which parent a copy came from however well it is measured.
 *
 * Counted on the sample's own called markers, since a marker the sample did not call is not
 * available to it whatever the parents read.
 *
 * Run: OM_TRIOS=<dir> [OM_CELLS=n] node --experimental-strip-types audit/marker-census.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('marker-census: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

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
const cache = new Map<string, { gt: Map<string, string> }>()
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
  cache.set(gsm, idx as never)
  return idx as unknown as { gt: Map<string, string> }
}

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const trios: { gsm: string; father: string; mother: string }[] = []
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  const mother = at(r, 'mother_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  const gsm = at(r, 'gsm')
  if (gsm && father && mother && fileOf.has(gsm)) trios.push({ gsm, father, mother })
}
const LIMIT = Number(process.env.OM_CELLS ?? trios.length)
const med = (xs: number[]) => xs.length ? xs.slice().sort((a, b) => a - b)[xs.length >> 1] : NaN

const perCell: { gsm: string; rows: number; called: number; informative: number; perChrom: number[]; aabb: number; bbaa: number }[] = []
for (const t of trios.slice(0, LIMIT)) {
  const pat = parentIndex(t.father); const mat = parentIndex(t.mother)
  let rows = 0; let called = 0; let aabb = 0; let bbaa = 0
  const byChrom = new Map<string, number>()
  for (const r of rowsOf(fileOf.get(t.gsm)!)()) {
    rows += 1
    if (r.genotype === 'NC' || r.genotype === null) continue
    called += 1
    const fa = pat.gt.get(r.probesetId); const mo = mat.gt.get(r.probesetId)
    // BOTH orientations are informative and both are counted. The share is flipped for the second
    // so that it always means the paternal fraction; see paternalShare.
    if (fa === 'AA' && mo === 'BB') { aabb += 1; byChrom.set(r.chrom, (byChrom.get(r.chrom) ?? 0) + 1) }
    else if (fa === 'BB' && mo === 'AA') { bbaa += 1; byChrom.set(r.chrom, (byChrom.get(r.chrom) ?? 0) + 1) }
  }
  const autosomes = [...byChrom.entries()].filter(([c]) => /^\d+$/.test(c)).map(([, n]) => n)
  const informative = autosomes.reduce((a, b) => a + b, 0)
  perCell.push({ gsm: t.gsm, rows, called, informative, perChrom: autosomes, aabb, bbaa })
  process.stdout.write('.')
}
console.log()
const pct = (a: number, b: number) => `${((100 * a) / b).toFixed(2)}%`
console.log(`\n${perCell.length} cells, each against its own two parents\n`)
console.log(`${'cell'.padEnd(12)} ${'markers on array'.padStart(17)} ${'called in cell'.padStart(15)} ${'informative'.padStart(12)} ${'of called'.padStart(10)} ${'median/chrom'.padStart(13)}`)
for (const c of perCell) {
  console.log(`${c.gsm.padEnd(12)} ${c.rows.toLocaleString().padStart(17)} ${c.called.toLocaleString().padStart(15)} `
    + `${c.informative.toLocaleString().padStart(12)} ${pct(c.informative, c.called).padStart(10)} ${med(c.perChrom).toLocaleString().padStart(13)}`)
}
const inf = perCell.map((c) => c.informative)
const allChrom = perCell.flatMap((c) => c.perChrom)
console.log(`\nACROSS ${perCell.length} CELLS`)
console.log(`  markers on the array            ${med(perCell.map((c) => c.rows)).toLocaleString()} (median)`)
console.log(`  called in the cell              ${med(perCell.map((c) => c.called)).toLocaleString()} (median), `
  + `${pct(med(perCell.map((c) => c.called)), med(perCell.map((c) => c.rows)))} of the array`)
console.log(`  informative, autosomes          ${med(inf).toLocaleString()} (median), range ${Math.min(...inf).toLocaleString()} to ${Math.max(...inf).toLocaleString()}`)
console.log(`  informative as a share of called ${pct(med(inf), med(perCell.map((c) => c.called)))}`)
console.log(`  informative per chromosome      ${med(allChrom).toLocaleString()} (median over ${allChrom.length} chromosomes), `
  + `range ${Math.min(...allChrom).toLocaleString()} to ${Math.max(...allChrom).toLocaleString()}`)
const under = allChrom.filter((n) => n < 200).length
console.log(`  chromosomes under the 200 floor ${under} of ${allChrom.length} (${pct(under, allChrom.length)})`)

const A = perCell.map((c) => c.aabb); const B = perCell.map((c) => c.bbaa)
console.log(`\n  father AA, mother BB            ${med(A).toLocaleString()} (median)`)
console.log(`  father BB, mother AA            ${med(B).toLocaleString()} (median)`)
console.log(`  both orientations are counted, and the second has its share flipped so that it`)
console.log(`  always reads as the paternal fraction. Neither is being left out.`)

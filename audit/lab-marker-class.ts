/**
 * The three specific figures an external review reported on the LABORATORY's arrays, reproduced.
 *
 * Its measurements were on the 32477 experiment and the 52461 children, where mine were on the
 * public corpus with a published karyotype, and the two disagree by a wide margin. These are the
 * review's own numbers on the review's own material:
 *
 *   W1  the opposite-homozygote rate of every 52461 sample against the sperm donor, which is the
 *       basis of its claim that only 7 of 18 are consistent with his paternity.
 *   W2  of the chromosomes of the confirmed children, how many fall below ONE_PARENT_HAPLOID_MAX
 *       in each marker class. The review reports 111 of 154 on the shipped class.
 *   W3  a sample of that series chr21, which the review reports at 312 parent-heterozygous markers against a
 *       typical 4,000, passing MIN_INFORMATIVE and reporting nonsense.
 *
 * No laboratory sample identifier is printed: samples are reported by position in the series.
 *
 * Run: OM_TRIOS=<dir with the lab arrays> node --experimental-strip-types audit/lab-marker-class.ts
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

// The series are discovered rather than named, and reported by position, so no identifier is
// printed or committed. The reference array of each series is the one the pairs file names.
const series = new Map<string, string[]>()
for (const f of readdirSync(DIR)) {
  const m = /^(\d{4,6})-(\d{2})\.probes/.exec(f)
  if (!m) continue
  series.set(m[1], [...(series.get(m[1]) ?? []), f].sort())
}
const label = new Map<string, string>()
let si = 0
for (const k of [...series.keys()].sort()) { si += 1; label.set(k, `series ${si}`) }

const gtOf = (file: string) => {
  const m = new Map<string, string>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (const r of rowsOf(join(DIR, file))()) {
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') {
      m.set(r.probesetId, r.genotype)
      pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
    }
  }
  return { gt: m, pos }
}

console.log('W1. THE OPPOSITE-HOMOZYGOTE RATE AGAINST EACH SERIES REFERENCE ARRAY')
console.log('    a true child reads 0.003 to 0.019 by this project\'s record; an unrelated adult near 0.07')
for (const [key, files] of [...series].sort()) {
  // The reference is the -00 array where one exists, which is the measured parent of that series.
  const ref = files.find((f) => /-00\.probes/.test(f))
  if (!ref) continue
  const R = gtOf(ref)
  console.log(`\n  ${label.get(key)}, reference array is position 00, ${files.length - 1} samples`)
  console.log(`  ${'position'.padEnd(10)} ${'oppHom'.padStart(8)} ${'n'.padStart(9)} `
    + `${'callRate'.padStart(9)} ${'het'.padStart(7)}  reading`)
  const rates: { pos: string; rate: number }[] = []
  for (const f of files) {
    if (f === ref) continue
    const S = gtOf(f)
    let n = 0; let opp = 0; let called = 0; let het = 0; let rows = 0
    for (const r of rowsOf(join(DIR, f))()) {
      rows += 1
      const sg = r.genotype
      if (sg !== 'AA' && sg !== 'BB' && sg !== 'AB') continue
      called += 1
      if (sg === 'AB') het += 1
      if (!/^\d+$/.test(r.chrom)) continue
      const rg = R.gt.get(r.probesetId)
      if (!rg || rg === 'AB' || sg === 'AB') continue
      n += 1
      if (sg !== rg) opp += 1
    }
    const rate = n ? opp / n : NaN
    const p = /-(\d{2})\.probes/.exec(f)?.[1] ?? '??'
    rates.push({ pos: p, rate })
    const read = !Number.isFinite(rate) ? '-'
      : called / rows < 0.60 ? 'no usable DNA'
      : rate < 0.010 ? 'consistent with this reference'
      : rate < 0.040 ? 'intermediate'
      : 'carries no genome of this reference'
    console.log(`  ${p.padEnd(10)} ${(100 * rate).toFixed(3).padStart(7)}% ${n.toLocaleString().padStart(9)} `
      + `${(called / rows).toFixed(2).padStart(9)} ${(het / called).toFixed(3).padStart(7)}  ${read}`)
  }
  const kids = rates.filter((r) => r.rate < 0.010).length
  const not = rates.filter((r) => r.rate >= 0.040).length
  console.log(`  consistent ${kids}, intermediate ${rates.length - kids - not}, not his ${not}`)
}

console.log('\n\nW2 AND W3. THE TWO MARKER CLASSES ON THE SAMPLES THAT ARE THIS REFERENCE\'S CHILDREN')
for (const [key, files] of [...series].sort()) {
  const ref = files.find((f) => /-00\.probes/.test(f))
  if (!ref) continue
  const R = gtOf(ref)
  // Children only, by the project's own linkage statistic.
  const kids: string[] = []
  for (const f of files) {
    if (f === ref) continue
    let n = 0; let opp = 0
    for (const r of rowsOf(join(DIR, f))()) {
      const sg = r.genotype
      if (!/^\d+$/.test(r.chrom) || (sg !== 'AA' && sg !== 'BB')) continue
      const rg = R.gt.get(r.probesetId)
      if (!rg || rg === 'AB') continue
      n += 1; if (sg !== rg) opp += 1
    }
    if (n && opp / n < 0.010) kids.push(f)
  }
  if (!kids.length) continue
  console.log(`\n  ${label.get(key)}: ${kids.length} samples consistent with the reference`)
  console.log(`  ${'pos'.padEnd(5)} ${'chr'.padStart(3)} ${'nHom'.padStart(7)} ${'nHet'.padStart(6)} `
    + `${'hetHom'.padStart(7)} ${'hetHet'.padStart(7)}  below the haploid max?`)
  let belowHom = 0; let belowHet = 0; let total = 0; const thin: string[] = []
  for (const f of kids) {
    const p = /-(\d{2})\.probes/.exec(f)?.[1] ?? '??'
    type Acc = { nHom: number; hetHom: number; nHet: number; hetHet: number }
    const by = new Map<string, Acc>()
    for (const r of rowsOf(join(DIR, f))()) {
      if (!/^\d+$/.test(r.chrom)) continue
      const sg = r.genotype
      if (sg !== 'AA' && sg !== 'BB' && sg !== 'AB') continue
      const rg = R.gt.get(r.probesetId)
      if (!rg) continue
      let a = by.get(r.chrom)
      if (!a) { a = { nHom: 0, hetHom: 0, nHet: 0, hetHet: 0 }; by.set(r.chrom, a) }
      if (rg === 'AB') { a.nHet += 1; if (sg === 'AB') a.hetHet += 1 }
      else { a.nHom += 1; if (sg === 'AB') a.hetHom += 1 }
    }
    for (const [chrom, a] of [...by].sort((x, y) => Number(x[0]) - Number(y[0]))) {
      if (a.nHom < (obl.MIN_INFORMATIVE as number)) continue
      total += 1
      const rh = a.hetHom / a.nHom; const rt = a.nHet ? a.hetHet / a.nHet : NaN
      const bh = rh < (obl.ONE_PARENT_HAPLOID_MAX as number)
      const bt = rt < (obl.ONE_PARENT_HAPLOID_MAX as number)
      if (bh) belowHom += 1
      if (bt) belowHet += 1
      if (a.nHet < 1000) thin.push(`pos ${p} chr${chrom}: ${a.nHom} hom, ${a.nHet} het, `
        + `hetHom ${rh.toFixed(4)}, hetHet ${Number.isFinite(rt) ? rt.toFixed(4) : '-'}`)
    }
  }
  console.log(`  chromosomes scored: ${total}`)
  console.log(`  below ONE_PARENT_HAPLOID_MAX on the shipped parent-HOM class: ${belowHom} of ${total}`
    + `  (${(100 * belowHom / Math.max(1, total)).toFixed(1)}%)`)
  console.log(`  below it on the parent-HET class:                             ${belowHet} of ${total}`
    + `  (${(100 * belowHet / Math.max(1, total)).toFixed(1)}%)`)
  if (thin.length) {
    console.log(`  chromosomes with under 1,000 parent-HET markers, which the floor of `
      + `${obl.MIN_INFORMATIVE} admits:`)
    for (const t of thin.slice(0, 6)) console.log(`    ${t}`)
  } else {
    console.log('  no chromosome has under 1,000 parent-HET markers')
  }
}

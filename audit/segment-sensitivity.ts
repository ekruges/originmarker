/**
 * What the missed segmental deletions fail, and what catching them would cost.
 *
 * Sweeps the two constants the absence scan turns on, the smallest window it will consider and the
 * likelihood-ratio threshold, over the cells whose karyotype was read by hand. Reports detection of
 * the hand-called deletions against calls on chromosomes the same hand reads as normal, so a gain
 * in sensitivity is always shown beside what it costs.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/segment-sensitivity.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const seg = await import(`${W}segments.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('segment-sensitivity: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

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
  if (!map) throw new Error(`${path}: no usable header`)
  return function* () { for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r } }
}
const cache = new Map<string, { gt: Map<string, string> }>()
function parentIndex(gsm: string) {
  const hit = cache.get(gsm); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  cache.set(gsm, idx as never); return idx as unknown as { gt: Map<string, string> }
}

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const parentsOf = new Map<string, { father: string; mother: string }>()
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  const mother = at(r, 'mother_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  if (father && mother) parentsOf.set(at(r, 'gsm'), { father, mother })
}
const hand = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , , chrom, manual] = l.split('\t')
  hand.set(`${gsm}|${chrom}`, manual)
}
const cells = [...new Set([...hand.keys()].map((k) => k.split('|')[0]))]
  .filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()

// Per cell: the two absence series and their nulls, computed once.
type Abs = { chrom: string; pos: number; absent: boolean }
// pat/mat: one parent homozygous, the diluted set the scan ships with.
// oppPat/oppMat: BOTH parents homozygous for opposite alleles, where a lost copy shows at every
// marker instead of half of them. Fewer markers, undiluted signal. Which wins is measured below.
interface Series { pat: Map<string, Abs[]>; mat: Map<string, Abs[]>
  oppPat: Map<string, Abs[]>; oppMat: Map<string, Abs[]> }
const byCell = new Map<string, Series>()
for (const gsm of cells) {
  const p = parentsOf.get(gsm)!
  const pa = parentIndex(p.father); const ma = parentIndex(p.mother)
  const s: Series = { pat: new Map(), mat: new Map(), oppPat: new Map(), oppMat: new Map() }
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!/^\d+$/.test(r.chrom) || r.genotype === 'NC') continue
    const fa = pa.gt.get(r.probesetId) ?? 'NC'; const mo = ma.gt.get(r.probesetId) ?? 'NC'
    if (fa === 'AA' || fa === 'BB') {
      const a = s.pat.get(r.chrom) ?? []; a.push({ chrom: r.chrom, pos: r.pos, absent: r.genotype !== 'AB' && r.genotype !== fa }); s.pat.set(r.chrom, a)
    }
    if (mo === 'AA' || mo === 'BB') {
      const a = s.mat.get(r.chrom) ?? []; a.push({ chrom: r.chrom, pos: r.pos, absent: r.genotype !== 'AB' && r.genotype !== mo }); s.mat.set(r.chrom, a)
    }
    if ((fa === 'AA' && mo === 'BB') || (fa === 'BB' && mo === 'AA')) {
      const op = s.oppPat.get(r.chrom) ?? []
      op.push({ chrom: r.chrom, pos: r.pos, absent: r.genotype === mo })
      s.oppPat.set(r.chrom, op)
      const om = s.oppMat.get(r.chrom) ?? []
      om.push({ chrom: r.chrom, pos: r.pos, absent: r.genotype === fa })
      s.oppMat.set(r.chrom, om)
    }
  }
  byCell.set(gsm, s); process.stdout.write('.')
}
console.log('\n')

const tally = (m: Map<string, { absent: boolean }[]>) => {
  const t = new Map<string, [number, number]>()
  for (const [c, ms] of m) t.set(c, [ms.length, ms.filter((x) => x.absent).length])
  return t
}
const isDel = (h: string) => /del|monosomy|nullisomy/i.test(h) && !/^normal/i.test(h)
const isNormal = (h: string) => /^normal/i.test(h.trim())

// --- the split and the rule, exactly as audit/segment-statistic-preregistration.txt fixes them ---
const tune = new Set(cells.filter((_, i) => i % 2 === 0))
const half = (gsm: string) => (tune.has(gsm) ? 'tune' : 'held')

function run(which: 'shipped' | 'opposite' | 'before', minW: number, thr: number, want: 'tune' | 'held' | 'all') {
  let found = 0; let n = 0; let fp = 0; let normals = 0; let segFound = 0; let segWant = 0
  for (const gsm of cells) {
    if (want !== 'all' && half(gsm) !== want) continue
    const s = byCell.get(gsm)!
    // 'before' is the one-parent scan exactly as it stood before the symmetric change.
    const A = which === 'opposite' ? s.oppPat : s.pat
    const B = which === 'before' ? new Map<string, Abs[]>() : (which === 'shipped' ? s.mat : s.oppMat)
    const nA = tally(A); const nB = tally(B)
    for (const chrom of new Set([...A.keys(), ...B.keys()])) {
      const h = hand.get(`${gsm}|${chrom}`) ?? ''
      const hits = [
        ...seg.scanChromosome(A.get(chrom) ?? [], seg.externalNull(nA as never, chrom), minW, thr),
        ...seg.scanChromosome(B.get(chrom) ?? [], seg.externalNull(nB as never, chrom), minW, thr),
      ]
      if (isDel(h)) {
        n += 1; if (hits.length) found += 1
        if (/del|dup/i.test(h)) { segWant += 1; if (hits.length) segFound += 1 }
      } else if (isNormal(h)) { normals += 1; if (hits.length) fp += 1 }
    }
  }
  return { found, n, fp, normals, segFound, segWant, rate: normals ? fp / normals : NaN }
}

const WIN = [2400, 1200, 600, 300, 150]
const LRT = [250, 150, 80]
console.log(`cells: ${cells.length} total, ${tune.size} tuning, ${cells.length - tune.size} held out`)
const tA = run('shipped', 2400, 250, 'tune'); const hA = run('shipped', 2400, 250, 'held')
console.log(`hand-called deletions: ${tA.n} on the tuning half, ${hA.n} on the held-out half`)
console.log(`hand-normal chromosomes: ${tA.normals} tuning, ${hA.normals} held out\n`)

console.log('TUNING HALF, the full grid. The cap is the shipped statistic at 2400/250 on this half.')
console.log(`  cap = ${tA.fp} of ${tA.normals} = ${(100 * tA.rate).toFixed(2)}%\n`)
console.log(`  ${'window'.padStart(7)} ${'LRT'.padStart(5)} ${'deletions found'.padStart(17)} ${'on hand-normal'.padStart(16)} ${'within cap'.padStart(11)}`)
let best: { w: number; t: number; found: number } | null = null
for (const w of WIN) {
  for (const t of LRT) {
    const r = run('opposite', w, t, 'tune')
    const ok = r.rate <= tA.rate + 1e-12
    console.log(`  ${w.toString().padStart(7)} ${t.toString().padStart(5)} ${`${r.found} of ${r.n}`.padStart(17)} `
      + `${`${r.fp} of ${r.normals}`.padStart(16)} ${(ok ? 'yes' : 'no').padStart(11)}`)
    if (!ok) continue
    // most detections; ties to the higher threshold, then the larger window
    if (!best || r.found > best.found
      || (r.found === best.found && t > best.t)
      || (r.found === best.found && t === best.t && w > best.w)) best = { w, t, found: r.found }
  }
}
if (!best) { console.log('\nNO POINT MET THE CAP. The shipped statistic stays.'); process.exit(0) }
console.log(`\nTHE RULE PICKS window ${best.w}, LRT ${best.t}\n`)

const hB = run('opposite', best.w, best.t, 'held')
console.log('HELD-OUT HALF, cells neither statistic was tuned on')
console.log(`  ${'statistic'.padEnd(26)} ${'deletions found'.padStart(17)} ${'segmental only'.padStart(16)} ${'on hand-normal'.padStart(16)}`)
console.log(`  ${'shipped, at 2400/250'.padEnd(26)} ${`${hA.found} of ${hA.n}`.padStart(17)} ${`${hA.segFound} of ${hA.segWant}`.padStart(16)} ${`${hA.fp} of ${hA.normals}`.padStart(16)}`)
console.log(`  ${`opposite, at ${best.w}/${best.t}`.padEnd(26)} ${`${hB.found} of ${hB.n}`.padStart(17)} ${`${hB.segFound} of ${hB.segWant}`.padStart(16)} ${`${hB.fp} of ${hB.normals}`.padStart(16)}`)

console.log('\nAGAINST THE BAR FIXED BEFORE THE SPLIT')
const moreFound = hB.found >= hA.found
const noMoreFp = hB.fp <= hA.fp
console.log(`  detects at least as many:        ${moreFound ? 'MET' : 'NOT MET'} (${hB.found} against ${hA.found})`)
console.log(`  calls no more hand-normal:       ${noMoreFp ? 'MET' : 'NOT MET'} (${hB.fp} against ${hA.fp})`)
console.log(`  VERDICT: ${moreFound && noMoreFp ? 'the opposite-homozygote statistic ships' : 'the shipped statistic stays'}`)

console.log('\nHELD-OUT CURVE, so the chosen point can be read against what was available')
console.log(`  ${'window'.padStart(7)} ${'LRT'.padStart(5)} ${'deletions found'.padStart(17)} ${'on hand-normal'.padStart(16)}`)
for (const w of WIN) for (const t of LRT) {
  const r = run('opposite', w, t, 'held')
  console.log(`  ${w.toString().padStart(7)} ${t.toString().padStart(5)} ${`${r.found} of ${r.n}`.padStart(17)} ${`${r.fp} of ${r.normals}`.padStart(16)}`)
}

console.log('\n\nBEFORE AND AFTER THE SYMMETRIC SCAN, all 31 cells, shipped setting 2400/250')
const bf = run('before', 2400, 250, 'all')
const af = run('shipped', 2400, 250, 'all')
console.log(`  ${'scan'.padEnd(26)} ${'deletions found'.padStart(17)} ${'segmental only'.padStart(16)} ${'on hand-normal'.padStart(16)}`)
console.log(`  ${'one parent (before)'.padEnd(26)} ${`${bf.found} of ${bf.n}`.padStart(17)} ${`${bf.segFound} of ${bf.segWant}`.padStart(16)} ${`${bf.fp} of ${bf.normals}`.padStart(16)}`)
console.log(`  ${'both parents (now)'.padEnd(26)} ${`${af.found} of ${af.n}`.padStart(17)} ${`${af.segFound} of ${af.segWant}`.padStart(16)} ${`${af.fp} of ${af.normals}`.padStart(16)}`)
console.log(`\n  detections gained: ${af.found - bf.found}`)
console.log(`  hand-normal calls gained: ${af.fp - bf.fp}`)

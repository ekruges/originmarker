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

function run(which: 'shipped' | 'opposite', minW: number, thr: number) {
  let found = 0; let want = 0; let fp = 0; let normals = 0; let segFound = 0; let segWant = 0
  for (const gsm of cells) {
    const s = byCell.get(gsm)!
    const A = which === 'shipped' ? s.pat : s.oppPat
    const B = which === 'shipped' ? s.mat : s.oppMat
    const nA = tally(A); const nB = tally(B)
    for (const chrom of new Set([...A.keys(), ...B.keys()])) {
      const h = hand.get(`${gsm}|${chrom}`) ?? ''
      const hits = [
        ...seg.scanChromosome(A.get(chrom) ?? [], seg.externalNull(nA as never, chrom), minW, thr),
        ...seg.scanChromosome(B.get(chrom) ?? [], seg.externalNull(nB as never, chrom), minW, thr),
      ]
      if (isDel(h)) {
        want += 1; if (hits.length) found += 1
        if (/del|dup/i.test(h)) { segWant += 1; if (hits.length) segFound += 1 }
      } else if (isNormal(h)) { normals += 1; if (hits.length) fp += 1 }
    }
  }
  return { found, want, fp, normals, segFound, segWant }
}
console.log(`${'statistic'.padEnd(10)} ${'minWin'.padStart(7)} ${'LRT'.padStart(5)} `
  + `${'all deletions'.padStart(15)} ${'segmental only'.padStart(16)} ${'on hand-normal'.padStart(16)}`)
for (const which of ['shipped', 'opposite'] as const) {
  for (const minW of which === 'opposite' ? [300, 150] : [2400]) {
    for (const thr of [250, 150, 80]) {
      const r = run(which, minW, thr)
      console.log(`${which.padEnd(10)} ${minW.toString().padStart(7)} ${thr.toString().padStart(5)} `
        + `${`${r.found} of ${r.want}`.padStart(15)} ${`${r.segFound} of ${r.segWant}`.padStart(16)} `
        + `${`${r.fp} of ${r.normals}`.padStart(16)}`)
    }
  }
}

/**
 * The two tests of the seven missed deletion-bearing events.
 *
 * F: the four on cells the stage inference refuses. Declaring the stage is the documented way
 *    past that gate, and the gain scan already recovers 2 of 2 events that way. Does the segment
 *    scan recover theirs, and with the right parent.
 * G: the three terminal chr16 deletions on cells the scan accepts. Either the whole-chromosome
 *    channel already claimed chr16, in which case the scan skips it by design and the parent IS
 *    attributed without the breakpoint, or the scan ran and the window minimum is the limit.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/missed-segment-tests.ts
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


type Seg = {
  chrom: string; startBp: number; endBp: number; kind: string; spanBp: number
  absentParent?: 'paternal' | 'maternal'; contested?: boolean
  alleleGain?: { shift: number; informative: number }
}
type Res = {
  segments?: Seg[]
  stage?: { stage: string; why: string }
  chroms?: { chrom: string; verdict: string; lostParent?: string; aneuploidy?: string }[]
}

async function scoreDeclared(gsm: string, father: string, mother: string, declared?: string) {
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
  return await score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
    ...(declared ? { declaredStage: declared } : {}),
  } as never) as Res
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

// TEST F. The four misses on refused cells, all del(6q)(pat), plus the other refused cells for
// context. Declared single-cell, which is what an operator would do with a cell they know is one.
const F = [
  { gsm: 'GSM4472405', chrom: '6', hand: 'del(6q)(pat)', parent: 'paternal' },
  { gsm: 'GSM4552387', chrom: '6', hand: 'del(6q) (pat)', parent: 'paternal' },
  { gsm: 'GSM4552388', chrom: '6', hand: 'del(6q)(pat)', parent: 'paternal' },
  { gsm: 'GSM4552389', chrom: '6', hand: 'del(6q) (pat)', parent: 'paternal' },
]
console.log('TEST F. THE FOUR MISSES ON CELLS THE STAGE INFERENCE REFUSES')
console.log(`${'cell'.padEnd(12)} ${'declared'.padEnd(12)} ${'inferred'.padEnd(10)} `
  + `${'segs'.padStart(5)} ${'on the chrom'.padEnd(34)} parent`)
let fRecovered = 0
for (const t of F) {
  const p = parentsOf.get(t.gsm)!
  for (const declared of [undefined, 'single-cell']) {
    const res = await scoreDeclared(t.gsm, p.father, p.mother, declared)
    const here = (res.segments ?? []).filter((s) => s.chrom === t.chrom)
    const run = here.map((s) => s.absentParent).find(Boolean)
    if (declared && here.length) fRecovered += 1
    console.log(`${t.gsm.padEnd(12)} ${(declared ?? '(inferred)').padEnd(12)} `
      + `${(res.stage?.stage ?? '?').padEnd(10)} ${String((res.segments ?? []).length).padStart(5)} `
      + `${(here.length
        ? here.map((s) => `${s.kind} ${(s.startBp / 1e6).toFixed(1)}-${(s.endBp / 1e6).toFixed(1)}Mb`).join('; ')
        : 'nothing').padEnd(34)} `
      + `${run ?? '-'}${run && run !== t.parent ? '  WRONG' : (run ? '  right' : '')}`)
  }
}
console.log(`\n  recovered with the stage declared: ${fRecovered} of ${F.length}`)

// TEST G. The three terminal chr16 maternal deletions on cells the scan accepts.
const G = ['GSM4472423', 'GSM4552385', 'GSM4552386']
console.log('\nTEST G. THE THREE TERMINAL chr16 MATERNAL DELETIONS, ON ACCEPTED CELLS')
console.log(`${'cell'.padEnd(12)} ${'stage'.padEnd(14)} ${'chr16 verdict'.padEnd(16)} `
  + `${'lostParent'.padEnd(11)} ${'aneuploidy'.padEnd(11)} chr16 segments`)
let gSkipped = 0
for (const gsm of G) {
  const p = parentsOf.get(gsm)!
  const res = await scoreDeclared(gsm, p.father, p.mother)
  const c16 = (res.chroms ?? []).find((c) => c.chrom === '16')
  const here = (res.segments ?? []).filter((s) => s.chrom === '16')
  // The scan skips a chromosome the alleles already read as having lost a parent, or one already
  // called whole by copy number. Either is the same absence described twice.
  if (c16?.lostParent || c16?.aneuploidy) gSkipped += 1
  console.log(`${gsm.padEnd(12)} ${(res.stage?.stage ?? '?').padEnd(14)} `
    + `${(c16?.verdict ?? 'none').padEnd(16)} ${(c16?.lostParent ?? '-').padEnd(11)} `
    + `${(c16?.aneuploidy ?? '-').padEnd(11)} `
    + `${here.length ? here.map((s) => s.kind).join(',') : 'nothing'}`)
}
console.log(`\n  chr16 already claimed whole, so the segment scan skips it: ${gSkipped} of ${G.length}`)
console.log('  If that is 3 of 3 the parent IS attributed on these events and the BREAKPOINT is')
console.log('  not, which is a different shortfall from a missed event.')

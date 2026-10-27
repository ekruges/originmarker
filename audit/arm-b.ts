/**
 * Arm B scored as a CALLER: the same two-hypothesis decision rule applied to both marker classes,
 * so the only difference between the arms is which markers they use.
 *
 *   parent HOMOZYGOUS    two contributions  het = q(1-ado) + d
 *   parent HETEROZYGOUS  two contributions  het = 0.5(1-ado) + d
 *   either               one contribution   het = d
 *
 * Called at the shipped posterior, refused below it, scored per chromosome against the authors'
 * published karyotype with each parent taken as the loaded one in turn.
 *
 * X1-X3 fixed in audit/arm-b-and-rebuild-preregistration.txt.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/arm-b.ts
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









const op = await import(`${W}oneParentOrigin.ts`)
const Q = op.DEFAULT_Q as number
const DROP_IN = op.DROP_IN as number
const POST = op.CALL_POSTERIOR as number
const MIN_MARKERS = op.MIN_MARKERS as number

/** Two hypotheses, one decision rule, whichever marker class supplies the counts. */
function callTwo(het: number, n: number, pTwo: number, dropIn = DROP_IN) {
  if (n < MIN_MARKERS) return { verdict: 'refused' as const, post: NaN }
  const pOne = dropIn
  const lp = (p: number) => het * Math.log(p) + (n - het) * Math.log(1 - p)
  const a = lp(Math.min(0.999, Math.max(1e-9, pTwo)))
  const b = lp(Math.min(0.999, Math.max(1e-9, pOne)))
  const m = Math.max(a, b)
  const post = Math.exp(a - m) / (Math.exp(a - m) + Math.exp(b - m))
  if (post >= POST) return { verdict: 'two' as const, post }
  if (1 - post >= POST) return { verdict: 'one' as const, post }
  return { verdict: 'refused' as const, post }
}

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const parentsOf = new Map<string, { father: string; mother: string }>()
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const f = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  const m = at(r, 'mother_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  if (f && m) parentsOf.set(at(r, 'gsm'), { father: f, mother: m })
}
const hand = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , , chrom, manual] = l.split('\t')
  hand.set(`${gsm}|${chrom}`, manual)
}
const cells = [...new Set([...hand.keys()].map((k) => k.split('|')[0]))]
  .filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()

type Score = { right: number; wrong: number; refused: number }
const blank = (): Score => ({ right: 0, wrong: 0, refused: 0 })
const tally: Record<string, Record<string, Score>> = {
  hom: { one: blank(), two: blank() }, het: { one: blank(), two: blank() },
}

for (const g of cells) {
  const p = parentsOf.get(g)!
  const fa = parentIndex(p.father) as unknown as { gt: Map<string, string> }
  const mo = parentIndex(p.mother) as unknown as { gt: Map<string, string> }
  type Acc = { nHom: number; hHom: number; nHet: number; hHet: number; called: number; markers: number }
  const by = new Map<string, { father: Acc; mother: Acc }>()
  const blankA = (): Acc => ({ nHom: 0, hHom: 0, nHet: 0, hHet: 0, called: 0, markers: 0 })
  for (const r of rowsOf(fileOf.get(g)!)()) {
    if (!/^\d+$/.test(r.chrom)) continue
    let e = by.get(r.chrom)
    if (!e) { e = { father: blankA(), mother: blankA() }; by.set(r.chrom, e) }
    e.father.markers += 1; e.mother.markers += 1
    const s = r.genotype
    if (s !== 'AA' && s !== 'BB' && s !== 'AB') continue
    e.father.called += 1; e.mother.called += 1
    for (const [who, gt] of [['father', fa.gt.get(r.probesetId)],
      ['mother', mo.gt.get(r.probesetId)]] as const) {
      if (!gt || gt === 'NC') continue
      const a = e[who]
      if (gt === 'AB') { a.nHet += 1; if (s === 'AB') a.hHet += 1 }
      else { a.nHom += 1; if (s === 'AB') a.hHom += 1 }
    }
  }
  // Dropout from the sample's own heterozygosity, as the shipped caller does.
  let called = 0; let het = 0
  for (const e of by.values()) { called += e.father.called; het += e.father.hHom + e.father.hHet }
  const ado = op.inferDropout(called ? het / called : NaN) as number
  for (const [chrom, e] of by) {
    const h = hand.get(`${g}|${chrom}`) ?? ''
    const normal = /^normal/i.test(h.trim())
    const lostPat = /\(\s*pat/i.test(h) && /monosom|nullisom/i.test(h)
    const lostMat = /\(\s*mat/i.test(h) && /monosom|nullisom/i.test(h)
    for (const who of ['father', 'mother'] as const) {
      const a = e[who]
      const otherGone = who === 'father' ? lostMat : lostPat
      const truth = normal ? 'two' : otherGone ? 'one' : '';
      if (!truth) continue
      // Refuse a chromosome that was barely genotyped, as the shipped pipeline now does.
      if (a.markers > 0 && a.called / a.markers < 0.60) continue
      for (const [cls, n, hh, pTwo] of [
        ['hom', a.nHom, a.hHom, Q * (1 - ado) + DROP_IN],
        ['het', a.nHet, a.hHet, 0.5 * (1 - ado) + DROP_IN]] as const) {
        const c = callTwo(hh, n, pTwo)
        const s2 = tally[cls][truth]
        if (c.verdict === 'refused') s2.refused += 1
        else if (c.verdict === truth) s2.right += 1
        else s2.wrong += 1
      }
    }
  }
  process.stdout.write('.')
}
console.log('\n')
console.log(`X1-X3. THE SAME DECISION RULE ON BOTH MARKER CLASSES, posterior ${POST}`)
console.log(`${'class'.padEnd(22)} ${'truth'.padEnd(6)} ${'right'.padStart(6)} ${'wrong'.padStart(6)} `
  + `${'refused'.padStart(8)} ${'accuracy of calls'.padStart(18)}`)
for (const cls of ['hom', 'het'] as const) {
  const name = cls === 'hom' ? 'parent-HOM (shipped)' : 'parent-HET (excluded)'
  for (const t of ['one', 'two'] as const) {
    const s = tally[cls][t]
    const made = s.right + s.wrong
    console.log(`${name.padEnd(22)} ${t.padEnd(6)} ${String(s.right).padStart(6)} `
      + `${String(s.wrong).padStart(6)} ${String(s.refused).padStart(8)} `
      + `${(made ? s.right / made : NaN).toFixed(4).padStart(18)}`)
  }
  const o = tally[cls].one; const w = tally[cls].two
  const right = o.right + w.right; const made = right + o.wrong + w.wrong
  const all = made + o.refused + w.refused
  console.log(`${name.padEnd(22)} ${'BOTH'.padEnd(6)} ${String(right).padStart(6)} `
    + `${String(o.wrong + w.wrong).padStart(6)} ${String(o.refused + w.refused).padStart(8)} `
    + `${(made ? right / made : NaN).toFixed(4).padStart(18)}`
    + `   answered ${((made / Math.max(1, all)) * 100).toFixed(1)}% of chromosomes`)
}

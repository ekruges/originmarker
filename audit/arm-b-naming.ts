/**
 * Arm (b): can the absent parent be named with ONE parent array loaded?
 *
 * Two statistics, both computable from one parent. `hetCallOneParent` on markers where the loaded
 * parent is HETEROZYGOUS says whether the chromosome carries one contribution or two. The
 * opposite-homozygote rate at markers where the loaded parent is HOMOZYGOUS says whether the
 * loaded parent's copy is the one present: at a father-AA marker a sample reads BB only if his
 * copy is absent, because he cannot transmit B.
 *
 * Whole-chromosome losses only. A whole-chromosome statistic applied to a segmental event reads as
 * near-intact.
 *
 * Criteria are in audit/arm-b-naming-preregistration.txt and were committed before this ran.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/arm-b-naming.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('arm-b-naming: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

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

const obl = await import(`${W}obligateHet.ts`)
const op = await import(`${W}oneParentOrigin.ts`)
const FLOOR = obl.MIN_INFORMATIVE as number

// Table S3, whole-chromosome losses with a parent label, plus every chromosome it reads as normal.
type Truth = { lost?: 'maternal' | 'paternal'; normal: boolean }
const truth = new Map<string, Map<string, Truth>>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const c = l.split('\t')
  const gsm = c[0]; const chrom = c[3]; const clause = (c[4] ?? '').toLowerCase()
  if (!gsm || !chrom) continue
  let m = truth.get(gsm)
  if (!m) { m = new Map(); truth.set(gsm, m) }
  const whole = /monosomy/.test(clause)
  const segmental = /\bdel\b|\bdup\b|breakpoint/.test(clause)
  const par = /\(pat\)/.test(clause) ? 'paternal' as const
    : /\(mat\)/.test(clause) ? 'maternal' as const : undefined
  // A chromosome the table does not call abnormal is normal; a segmental event is neither, and is
  // excluded from both arms rather than counted as normal.
  if (whole && par) m.set(chrom, { lost: par, normal: false })
  else if (segmental) m.set(chrom, { normal: false })
  else if (!m.has(chrom)) m.set(chrom, { normal: !/monosomy|trisomy|gain|loss/.test(clause) })
}

const manifestPath = join(DIR, 'trio_manifest_full.csv')
if (!existsSync(manifestPath)) { console.log('no manifest. NOT RUN.'); process.exit(0) }
const manifest = readFileSync(manifestPath, 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''

type Row = { gsm: string; chrom: string; nHom: number; hetHom: number; nHet: number
  hetHet: number; oppN: number; opp: number; callRate: number; ado: number; t: Truth }
const rows: Row[] = []
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const gsm = at(r, 'gsm')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  if (!gsm || !father || !fileOf.has(gsm) || !truth.has(gsm)) continue
  const fa = new Map<string, string>()
  for (const x of rowsOf(fileOf.get(father)!)()) {
    if (x.genotype === 'AA' || x.genotype === 'BB' || x.genotype === 'AB') fa.set(x.probesetId, x.genotype)
  }
  type Acc = { nHom: number; hetHom: number; nHet: number; hetHet: number
    oppN: number; opp: number; called: number; onArray: number }
  const by = new Map<string, Acc>()
  let allCalled = 0; let allRows = 0; let allHet = 0
  for (const x of rowsOf(fileOf.get(gsm)!)()) {
    if (!/^\d+$/.test(x.chrom)) continue
    let a = by.get(x.chrom)
    if (!a) { a = { nHom: 0, hetHom: 0, nHet: 0, hetHet: 0, oppN: 0, opp: 0, called: 0, onArray: 0 }
      by.set(x.chrom, a) }
    a.onArray += 1; allRows += 1
    const sg = x.genotype
    if (sg !== 'AA' && sg !== 'BB' && sg !== 'AB') continue
    a.called += 1; allCalled += 1
    if (sg === 'AB') allHet += 1
    const g = fa.get(x.probesetId)
    if (!g) continue
    if (g === 'AB') { a.nHet += 1; if (sg === 'AB') a.hetHet += 1 } else {
      a.nHom += 1
      if (sg === 'AB') a.hetHom += 1
      else { a.oppN += 1; if (sg !== g) a.opp += 1 }
    }
  }
  if (!allRows || allCalled / allRows < 0.60) continue
  const ado = op.inferDropout(allHet / allCalled) as number
  for (const [chrom, a] of by) {
    const t = truth.get(gsm)!.get(chrom)
    if (!t) continue
    const cr = a.onArray ? a.called / a.onArray : NaN
    if (!(obl.informativeEnough(a.nHom, cr) as boolean)) continue
    rows.push({ gsm, chrom, nHom: a.nHom, hetHom: a.hetHom, nHet: a.nHet, hetHet: a.hetHet,
      oppN: a.oppN, opp: a.opp, callRate: cr, ado, t })
  }
}

const rate = (r: Row) => (r.oppN ? r.opp / r.oppN : NaN)
const losses = rows.filter((r) => r.t.lost)
const normals = rows.filter((r) => r.t.normal)
console.log(`${rows.length} chromosome-observations scored: ${losses.length} hand-called `
  + `whole-chromosome losses with a parent label, ${normals.length} read as normal\n`)

// CRITERION 3 FIRST: the separation, before any threshold is applied.
console.log('CRITERION 3. THE OPPOSITE-HOMOZYGOTE RATE, BY WHAT THE HAND SAYS')
const show = (name: string, set: Row[]) => {
  const v = set.map(rate).filter(Number.isFinite).sort((a, b) => a - b)
  if (!v.length) { console.log(`  ${name.padEnd(34)} none`); return }
  const q = (p: number) => v[Math.min(v.length - 1, Math.floor(p * v.length))]
  console.log(`  ${name.padEnd(34)} n=${String(v.length).padStart(4)}  `
    + `min ${(100 * v[0]).toFixed(2)}%  median ${(100 * q(0.5)).toFixed(2)}%  `
    + `max ${(100 * v[v.length - 1]).toFixed(2)}%`)
}
show('normal, father present', normals)
show('maternal loss, father present', losses.filter((r) => r.t.lost === 'maternal'))
show('paternal loss, father ABSENT', losses.filter((r) => r.t.lost === 'paternal'))
const present = [...normals, ...losses.filter((r) => r.t.lost === 'maternal')].map(rate)
  .filter(Number.isFinite)
const absent = losses.filter((r) => r.t.lost === 'paternal').map(rate).filter(Number.isFinite)
const gapLo = present.length ? Math.max(...present) : NaN
const gapHi = absent.length ? Math.min(...absent) : NaN
console.log(`\n  highest rate where his copy is present: ${(100 * gapLo).toFixed(2)}%`)
console.log(`  lowest rate where his copy is absent:   ${(100 * gapHi).toFixed(2)}%`)
console.log(gapHi > gapLo
  ? `  PASS criterion 3: a gap of ${(100 * (gapHi - gapLo)).toFixed(2)} points, so a boundary exists.`
  : '  FAIL criterion 3 AS PRE-REGISTERED: the ranges overlap, so the boundary below is chosen and'
    + '\n  not derived. What stands in its place is the sweep at the end: whether the result moves'
    + '\n  when the boundary does. The outliers are listed first.')

// The overlap is one-sided, so name the chromosomes that cause it rather than leaving it at a max.
const high = normals.filter((r) => rate(r) > 0.05).sort((a, b) => rate(b) - rate(a))
console.log(`\n  chromosomes the hand reads as NORMAL with an opposite-homozygote rate over 5%: `
  + `${high.length} of ${normals.length}`)
for (const r of high.slice(0, 10)) {
  console.log(`    ${r.gsm} chr${r.chrom}: opp ${(100 * rate(r)).toFixed(2)}%, `
    + `hetHet ${(r.hetHet / Math.max(1, r.nHet)).toFixed(4)}, hetHom `
    + `${(r.hetHom / Math.max(1, r.nHom)).toFixed(4)}, call rate ${r.callRate.toFixed(3)}`)
}
const lowLoss = losses.filter((r) => r.t.lost === 'paternal' && rate(r) < 0.05)
console.log(`  hand-called PATERNAL losses with a rate under 5%: ${lowLoss.length} of `
  + `${losses.filter((r) => r.t.lost === 'paternal').length}`)
for (const r of lowLoss.slice(0, 6)) {
  console.log(`    ${r.gsm} chr${r.chrom}: opp ${(100 * rate(r)).toFixed(2)}%, `
    + `hetHet ${(r.hetHet / Math.max(1, r.nHet)).toFixed(4)}`)
}

const BOUND = Number.isFinite(gapLo) && Number.isFinite(gapHi) && gapHi > gapLo
  ? (gapLo + gapHi) / 2 : 0.03
console.log(`\n  boundary used below, the midpoint of the gap: ${(100 * BOUND).toFixed(2)}%\n`)

// THE TABLE.
const callAt = (r: Row, bound: number) => {
  const one = obl.hetCallOneParent(r.hetHet, r.nHet, r.ado) as { ploidy: string }
  if (one.ploidy !== 'uniparental') return undefined
  return rate(r) < bound ? 'maternal' as const : 'paternal' as const
}
const call = (r: Row) => callAt(r, BOUND)
let right = 0; let wrong = 0; const missed: string[] = []
const wrongDetail: string[] = []
for (const r of losses) {
  const c = call(r)
  if (!c) { missed.push(`${r.gsm} chr${r.chrom} (${r.t.lost})`); continue }
  if (c === r.t.lost) right += 1
  else { wrong += 1; wrongDetail.push(`${r.gsm} chr${r.chrom}: hand ${r.t.lost}, named ${c}, `
    + `opp ${(100 * rate(r)).toFixed(2)}%`) }
}
console.log('CRITERION 1. DIRECTION, ON HAND-CALLED WHOLE-CHROMOSOME LOSSES')
console.log(`  named correctly:              ${right} of ${losses.length}`)
console.log(`  named wrongly:                ${wrong}`)
console.log(`  not named, read two contributions: ${missed.length}`)
for (const d of wrongDetail) console.log(`    ${d}`)
console.log(wrong * 2 < right
  ? '  PASS criterion 1: the wrong calls are fewer than half the right ones.'
  : `  FAIL criterion 1: ${wrong} wrong against ${right} right.`)

let falseNamed = 0
const falseDetail: string[] = []
for (const r of normals) {
  const c = call(r)
  if (c) { falseNamed += 1; falseDetail.push(`${r.gsm} chr${r.chrom}: named ${c}, `
    + `opp ${(100 * rate(r)).toFixed(2)}%, hetHet ${(r.hetHet / Math.max(1, r.nHet)).toFixed(4)}`) }
}
const fp = normals.length ? falseNamed / normals.length : NaN
console.log('\nCRITERION 2. FALSE NAMING, ON CHROMOSOMES THE HAND READS AS NORMAL')
console.log(`  named as carrying an absent parent: ${falseNamed} of ${normals.length} `
  + `(${(100 * fp).toFixed(1)}%)`)
for (const d of falseDetail.slice(0, 15)) console.log(`    ${d}`)
if (falseDetail.length > 15) console.log(`    and ${falseDetail.length - 15} more`)
console.log(fp < 0.05
  ? '  PASS criterion 2: under 5 percent.'
  : `  FAIL criterion 2: ${(100 * fp).toFixed(1)} percent, over the 5 the pre-registration set.`)

// WHAT REPLACES THE FAILED GAP TEST: if the answer is the same across every boundary a reader
// might have picked, the choice of boundary is not carrying the result.
console.log('\nTHE BOUNDARY SWEPT')
console.log(`  ${'bound'.padStart(6)} ${'right'.padStart(6)} ${'wrong'.padStart(6)} `
  + `${'unnamed'.padStart(8)} ${'false'.padStart(6)} ${'false%'.padStart(7)}`)
for (const b of [0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]) {
  let rt = 0; let wr = 0; let un = 0
  for (const r of losses) {
    const c = callAt(r, b)
    if (!c) { un += 1; continue }
    if (c === r.t.lost) rt += 1; else wr += 1
  }
  const fn = normals.filter((r) => callAt(r, b)).length
  console.log(`  ${(100 * b).toFixed(1).padStart(5)}% ${String(rt).padStart(6)} `
    + `${String(wr).padStart(6)} ${String(un).padStart(8)} ${String(fn).padStart(6)} `
    + `${(100 * fn / Math.max(1, normals.length)).toFixed(1).padStart(6)}%`)
}

// LEAVE ONE CELL OUT. The boundary is chosen once for the corpus, so a per-cell split says whether
// any single cell is carrying the numbers.
const cells = [...new Set(losses.map((r) => r.gsm))]
let worst = 1; let worstCell = ''
for (const gsm of cells) {
  const set = losses.filter((r) => r.gsm === gsm)
  let rt = 0; let n = 0
  for (const r of set) { const c = call(r); if (c) { n += 1; if (c === r.t.lost) rt += 1 } }
  if (n >= 2) { const acc = rt / n; if (acc < worst) { worst = acc; worstCell = gsm } }
}
console.log(`\n  cells carrying a hand-called labelled loss: ${cells.length}`)
console.log(`  worst single cell with 2 or more named: ${worstCell} at `
  + `${(100 * worst).toFixed(0)}% correct`)

// ===========================================================================================
// THE SEGMENTAL ARM. The two statistics are per-marker, so applying them inside a region is the
// same measurement on fewer markers. That is the claim the wiring rests on, and it is tested here
// against Table S3's hand-called DELETIONS with a parent label, resolved to coordinates through
// the cytoband track. Duplications are excluded: a gain does not remove a parent's copy.
//
// Direction only. Arm (b) names the parent on a region something else has already called, so the
// risk it adds is measured by criterion 2 above, on regions nothing called.
// ===========================================================================================
type Band = { chrom: string; start: number; end: number; name: string }
const bands: Band[] = []
for (const l of readFileSync(new URL('./ucsc-hg19-cytoband.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#')) continue
  const c = l.split('\t')
  const chrom = (c[0] ?? '').replace(/^chr/, '')
  if (!/^\d+$/.test(chrom)) continue
  bands.push({ chrom, start: Number(c[1]), end: Number(c[2]), name: c[3] ?? '' })
}
const armOf = (chrom: string, arm: 'p' | 'q') => {
  const b = bands.filter((x) => x.chrom === chrom && x.name.startsWith(arm))
  if (!b.length) return undefined
  return { start: Math.min(...b.map((x) => x.start)), end: Math.max(...b.map((x) => x.end)) }
}
const chromEnd = (chrom: string) => Math.max(...bands.filter((x) => x.chrom === chrom).map((x) => x.end))
const bandStart = (chrom: string, name: string) => {
  const hit = bands.find((x) => x.chrom === chrom && x.name === name)
  return hit ? hit.start : undefined
}
// `del 16q23.1-tel`, `del12p`, `del (6q21-tel)`, `del8p23.1`, `del 6q`.
const resolve = (chrom: string, clause: string) => {
  const m = /del\s*\(?\s*(\d+)?\s*([pq])([\d.]+)?\s*(?:-\s*(tel|[\d.]+))?/.exec(clause)
  if (!m) return undefined
  const arm = m[2] as 'p' | 'q'
  const whole = armOf(chrom, arm)
  if (!whole) return undefined
  if (!m[3]) return whole
  const from = bandStart(chrom, `${arm}${m[3]}`)
  if (from === undefined) return undefined
  if (m[4] === 'tel') return { start: from, end: arm === 'q' ? chromEnd(chrom) : whole.end }
  if (m[4]) {
    const to = bands.find((x) => x.chrom === chrom && x.name === `${arm}${m[4]}`)
    if (to) return { start: Math.min(from, to.start), end: Math.max(from, to.end) }
    return undefined
  }
  const only = bands.find((x) => x.chrom === chrom && x.name === `${arm}${m[3]}`)
  return only ? { start: only.start, end: only.end } : undefined
}

type Seg = { gsm: string; chrom: string; start: number; end: number
  lost: 'maternal' | 'paternal'; clause: string }
const segs: Seg[] = []
const unresolved: string[] = []
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const c = l.split('\t')
  const gsm = c[0]; const chrom = c[3]; const clause = c[4] ?? ''
  if (!gsm || !chrom || !fileOf.has(gsm)) continue
  if (!/\bdel\b|\bdel\(|\bdel\d/.test(clause.toLowerCase())) continue
  const par = /\(pat\)/.test(clause) ? 'paternal' as const
    : /\(mat\)/.test(clause) ? 'maternal' as const : undefined
  if (!par) continue
  const co = resolve(chrom, clause)
  if (!co) { unresolved.push(`${gsm} chr${chrom}: ${clause.trim()}`); continue }
  segs.push({ gsm, chrom, start: co.start, end: co.end, lost: par, clause: clause.trim() })
}
console.log(`\n\nTHE SEGMENTAL ARM. ${segs.length} hand-called deletions with a parent label `
  + `resolved to coordinates, ${unresolved.length} not resolved`)
for (const u of unresolved) console.log(`  unresolved: ${u}`)

// Group by sample so each array is read once.
const bySample = new Map<string, Seg[]>()
for (const g of segs) {
  const k = bySample.get(g.gsm) ?? []
  k.push(g); bySample.set(g.gsm, k)
}
let sRight = 0; let sWrong = 0; let sUnnamed = 0; let sThin = 0
console.log(`\n  ${'gsm'.padEnd(12)} ${'region'.padEnd(22)} ${'Mb'.padStart(6)} ${'nHet'.padStart(6)} `
  + `${'hetHet'.padStart(7)} ${'opp'.padStart(7)} ${'hand'.padEnd(9)} named`)
for (const [gsm, list] of bySample) {
  const father = (() => {
    for (const line of manifest.slice(1)) {
      const r = line.split(',')
      if (at(r, 'gsm') !== gsm) continue
      return at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
    }
    return undefined
  })()
  if (!father) continue
  const fa = new Map<string, string>()
  for (const x of rowsOf(fileOf.get(father)!)()) {
    if (x.genotype === 'AA' || x.genotype === 'BB' || x.genotype === 'AB') fa.set(x.probesetId, x.genotype)
  }
  type A = { nHet: number; hetHet: number; oppN: number; opp: number }
  const accs = list.map(() => ({ nHet: 0, hetHet: 0, oppN: 0, opp: 0 } as A))
  let called = 0; let rowsN = 0; let hetN = 0
  for (const x of rowsOf(fileOf.get(gsm)!)()) {
    rowsN += 1
    const sg = x.genotype
    if (sg !== 'AA' && sg !== 'BB' && sg !== 'AB') continue
    called += 1
    if (sg === 'AB') hetN += 1
    const g = fa.get(x.probesetId)
    if (!g) continue
    for (let i = 0; i < list.length; i += 1) {
      const seg = list[i]
      if (x.chrom !== seg.chrom || x.pos < seg.start || x.pos > seg.end) continue
      const a = accs[i]
      if (g === 'AB') { a.nHet += 1; if (sg === 'AB') a.hetHet += 1 } else if (sg !== 'AB') {
        a.oppN += 1; if (sg !== g) a.opp += 1
      }
    }
  }
  if (!rowsN || called / rowsN < 0.60) continue
  const ado = op.inferDropout(hetN / called) as number
  for (let i = 0; i < list.length; i += 1) {
    const seg = list[i]; const a = accs[i]
    const one = obl.hetCallOneParent(a.hetHet, a.nHet, ado) as { ploidy: string; why: string }
    const r = a.oppN ? a.opp / a.oppN : NaN
    const named = one.ploidy !== 'uniparental' ? undefined
      : r < BOUND ? 'maternal' as const : 'paternal' as const
    if (a.nHet < 50) sThin += 1
    else if (!named) sUnnamed += 1
    else if (named === seg.lost) sRight += 1
    else sWrong += 1
    console.log(`  ${gsm.padEnd(12)} ${seg.clause.slice(0, 22).padEnd(22)} `
      + `${((seg.end - seg.start) / 1e6).toFixed(0).padStart(6)} ${a.nHet.toLocaleString().padStart(6)} `
      + `${(a.hetHet / Math.max(1, a.nHet)).toFixed(4).padStart(7)} `
      + `${(100 * r).toFixed(2).padStart(6)}% ${seg.lost.padEnd(9)} ${named ?? 'not named'}`)
  }
}
console.log(`\n  named correctly:              ${sRight} of ${sRight + sWrong + sUnnamed + sThin}`)
console.log(`  named wrongly:                ${sWrong}`)
console.log(`  not named, two contributions: ${sUnnamed}`)
console.log(`  under 50 parent-het markers:  ${sThin}`)
console.log(sWrong * 2 < sRight
  ? '  The table holds inside a region: the wrong calls are fewer than half the right ones.'
  : `  The table does NOT hold inside a region: ${sWrong} wrong against ${sRight} right.`)

console.log(`\nTHE SHIPPED ONE-PARENT CHANNEL NAMES NONE OF THESE ${losses.length} LOSSES TODAY.`)
void FLOOR

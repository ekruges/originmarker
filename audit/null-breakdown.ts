/**
 * The segment null's breakdown, and what the obligate set can name without it.
 *
 * P1-P3 are an external review's predictions about `externalNull`; S1-S4 are mine about the
 * obligate-heterozygote set. Fixed in audit/null-breakdown-preregistration.txt before this ran.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/null-breakdown.ts
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




type Seg = { chrom: string; kind: string; startBp: number; endBp: number; score: number
  absentParent?: string; contested?: boolean }
type Contribution = { ploidy: string; fraction: number; informative: number
  onlyA: number; onlyB: number; survivor?: string }
type Chrom = { chrom: string; verdict: string; contribution?: Contribution }
type Res = { segments?: Seg[]; stage?: { stage: string }; chroms?: Chrom[] }

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
    absenceByChrom: Map<string, { pos: number; absent: boolean }[]>
    matAbsenceByChrom: Map<string, { pos: number; absent: boolean }[]>
  }
  return { res, acc: a, profile: profile as unknown as { callRate: number } }
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
const FLOOR = segs2.NULL_FLOOR as number
const MININF = obl.MIN_INFORMATIVE as number
const SMIN = obl.SURVIVOR_MIN as number
// The review's normalisation, which the shipped code does not compute.
const theirs = (c?: Contribution) => (c && (c.onlyA + c.onlyB) > 0
  ? Math.max(c.onlyA, c.onlyB) / (c.onlyA + c.onlyB) : NaN)

// --- P1 and P2 ---------------------------------------------------------------------------------
console.log('P1 AND P2. THE chr16 NULL AND THE 16q23.1-tel WINDOW. 16q23.1 starts 74.1 Mb.')
console.log(`NULL_FLOOR ${FLOOR}\n`)
console.log(`${'cell'.padEnd(12)} ${'nullPat'.padStart(8)} ${'nullMat'.padStart(8)} `
  + `${'rateMat'.padStart(8)} ${'n'.padStart(6)} ${'oppHom'.padStart(7)} rejected before scoring?`)
for (const g of ['GSM4472423', 'GSM4552385', 'GSM4552386']) {
  const p = parentsOf.get(g)!
  const { res, acc } = await run(g, p.father, p.mother)
  const nullPat = segs2.externalNull(acc.t.byChrom as never, '16') as number
  const nullMat = segs2.externalNull(acc.matAbsenceTally as never, '16') as number
  // The window itself, 74.1 Mb to the telomere, on the maternal marker set.
  const ms = (acc.matAbsenceByChrom.get('16') ?? []).filter((m) => m.pos >= 74_100_000)
  const rate = ms.length ? ms.filter((m) => m.absent).length / ms.length : NaN
  const c16 = (res.chroms ?? []).find((c) => c.chrom === '16')
  console.log(`${g.padEnd(12)} ${nullPat.toFixed(4).padStart(8)} ${nullMat.toFixed(4).padStart(8)} `
    + `${rate.toFixed(4).padStart(8)} ${String(ms.length).padStart(6)} `
    + `${String(c16?.contribution?.informative ?? '-').padStart(7)} `
    + `${rate <= nullMat ? 'YES, rate <= null' : 'no, it would be scored'}`)
}

// --- P3 ----------------------------------------------------------------------------------------
console.log('\nP3. THE CONTESTED WINDOWS: ARE THE TWO ABSENCE RATES EQUAL OR ASYMMETRIC?')
console.log(`${'cell'.padEnd(12)} ${'chr'.padStart(3)} ${'nullPat'.padStart(8)} `
  + `${'nullMat'.padStart(8)} ${'ratePat'.padStart(8)} ${'rateMat'.padStart(8)} `
  + `${'ratio'.padStart(6)}  reading`)
let sym = 0; let asym = 0
for (const g of ['GSM4552378', 'GSM4552380', 'GSM4552382', 'GSM4552385', 'GSM4552386', 'GSM4552420']) {
  const p = parentsOf.get(g)!
  const { res, acc } = await run(g, p.father, p.mother)
  const done = new Set<string>()
  for (const sg of (res.segments ?? []).filter((s) => s.contested)) {
    if (done.has(sg.chrom)) continue
    done.add(sg.chrom)
    const win = (xs: { pos: number; absent: boolean }[]) => {
      const m = xs.filter((x) => x.pos >= sg.startBp && x.pos <= sg.endBp)
      return m.length ? m.filter((x) => x.absent).length / m.length : NaN
    }
    const rp = win(acc.absenceByChrom.get(sg.chrom) ?? [])
    const rm = win(acc.matAbsenceByChrom.get(sg.chrom) ?? [])
    const ratio = Math.max(rp, rm) / (rp + rm)
    if (ratio < 0.60) sym += 1; else asym += 1
    console.log(`${g.padEnd(12)} ${sg.chrom.padStart(3)} `
      + `${(segs2.externalNull(acc.t.byChrom as never, sg.chrom) as number).toFixed(4).padStart(8)} `
      + `${(segs2.externalNull(acc.matAbsenceTally as never, sg.chrom) as number).toFixed(4).padStart(8)} `
      + `${rp.toFixed(4).padStart(8)} ${rm.toFixed(4).padStart(8)} ${ratio.toFixed(3).padStart(6)}  `
      + `${ratio < 0.60 ? 'SYMMETRIC, dropout' : 'asymmetric, P3 FAILS here'}`)
  }
}
console.log(`\n  symmetric ${sym}, asymmetric ${asym}  (P3 holds if symmetric dominates)`)

// --- S1 ----------------------------------------------------------------------------------------
console.log('\nS1. chr6 ON THE FOUR REFUSED CELLS CARRYING A HAND-CALLED del(6q)(pat)')
console.log(`SURVIVOR_MIN ${SMIN} on the shipped denominator; the review normalises to onlyA+onlyB`)
console.log(`${'cell'.padEnd(12)} ${'stage'.padEnd(10)} ${'inf'.padStart(6)} ${'het'.padStart(7)} `
  + `${'onlyA'.padStart(7)} ${'onlyB'.padStart(7)} ${'shipped'.padStart(8)} ${'review'.padStart(7)} `
  + `${'survivor'.padStart(9)}  names?`)
for (const g of ['GSM4472405', 'GSM4552387', 'GSM4552388', 'GSM4552389']) {
  const p = parentsOf.get(g)!
  const { res } = await run(g, p.father, p.mother)
  const c = (res.chroms ?? []).find((x) => x.chrom === '6')?.contribution
  const shipped = c ? Math.max(c.onlyA, c.onlyB) : NaN
  const rev = theirs(c)
  const who = c?.survivor === 'a' ? 'paternal' : c?.survivor === 'b' ? 'maternal' : '-'
  console.log(`${g.padEnd(12)} ${(res.stage?.stage ?? '?').padEnd(10)} `
    + `${String(c?.informative ?? 0).padStart(6)} ${(c?.fraction ?? NaN).toFixed(4).padStart(7)} `
    + `${(c?.onlyA ?? NaN).toFixed(4).padStart(7)} ${(c?.onlyB ?? NaN).toFixed(4).padStart(7)} `
    + `${shipped.toFixed(3).padStart(8)} ${rev.toFixed(3).padStart(7)} ${who.padStart(9)}  `
    + `${shipped >= SMIN ? 'shipped names it' : rev >= 0.9 ? 'review normalisation names it' : 'neither'}`)
}

// --- S2 ----------------------------------------------------------------------------------------
// The rate at which a sample carrying ONE parent's allele alone reads an opposite homozygote
// against the OTHER parent. Measured on the parental arrays themselves: at markers where the
// mother is homozygous, how often is the father's transmitted allele the other one. With the
// father's genotype known, the two transmissible alleles are equally likely.
console.log('\nS2. THE RATE A TRUE SINGLE-PARENT LOSS PRODUCES, FROM THE PARENTAL ARRAYS ALONE')
{
  const someone = [...parentsOf.values()][0]
  const fa = parentIndex(someone.father) as unknown as { gt: Map<string, string> }
  const mo = parentIndex(someone.mother) as unknown as { gt: Map<string, string> }
  let n = 0; let opp = 0
  for (const [id, mg] of mo.gt) {
    if (mg !== 'AA' && mg !== 'BB') continue
    const fg = fa.gt.get(id)
    if (!fg || fg === 'NC') continue
    // the father transmits each of his alleles with probability 1/2
    const other = mg === 'AA' ? 'B' : 'A'
    const p = fg === 'AB' ? 0.5 : (fg[0] === other ? 1 : 0)
    n += 1; opp += p
  }
  console.log(`  maternal copy lost, read against the mother: ${(opp / n).toFixed(4)} over ${n.toLocaleString()} markers`)
  let n2 = 0; let opp2 = 0
  for (const [id, fg] of fa.gt) {
    if (fg !== 'AA' && fg !== 'BB') continue
    const mg = mo.gt.get(id)
    if (!mg || mg === 'NC') continue
    const other = fg === 'AA' ? 'B' : 'A'
    const p = mg === 'AB' ? 0.5 : (mg[0] === other ? 1 : 0)
    n2 += 1; opp2 += p
  }
  console.log(`  paternal copy lost, read against the father: ${(opp2 / n2).toFixed(4)} over ${n2.toLocaleString()} markers`)
  console.log('  the review derives 0.100 from a fitted spectrum; this assumes none')
}

// --- S3 ----------------------------------------------------------------------------------------
console.log('\nS3. OPPOSITE-HOMOZYGOTE MARKERS INSIDE 16q23.1-tel, AGAINST THE 200 FLOOR')
{
  const g = 'GSM4552385'; const p = parentsOf.get(g)!
  const fa = parentIndex(p.father) as unknown as { gt: Map<string, string> }
  const mo = parentIndex(p.mother) as unknown as { gt: Map<string, string> }
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (const r of rowsOf(fileOf.get(g)!)()) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  let whole = 0; let inWin = 0
  for (const [id, fg] of fa.gt) {
    const mg = mo.gt.get(id)
    if (!mg) continue
    if (!((fg === 'AA' && mg === 'BB') || (fg === 'BB' && mg === 'AA'))) continue
    const q = pos.get(id); if (!q || q.chrom !== '16') continue
    whole += 1
    if (q.pos >= 74_100_000) inWin += 1
  }
  console.log(`  chr16 whole:        ${whole.toLocaleString()} opposite-homozygote markers`)
  console.log(`  chr16 16q23.1-tel:  ${inWin.toLocaleString()}  against MIN_INFORMATIVE ${MININF}  `
    + `${inWin >= MININF ? 'ENOUGH' : 'TOO FEW, the obligate set cannot speak at window scale'}`)
  console.log('  the review scaled this to about 171 and flagged it as an estimate')
}

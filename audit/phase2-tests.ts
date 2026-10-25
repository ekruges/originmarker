/**
 * Phase 2 of the second review: four questions in one pass, plus the effect of conditioning the
 * heterozygosity spread on per-chromosome call rate.
 *
 * L: what the conditioned spread reads on every cell, against the unconditioned one it replaces.
 * M: how many chromosomes of an admitted cell then fail the per-chromosome mis-clustering gate.
 *    No-calls sit in that denominator, so the genome-level gate may be re-applying one chromosome
 *    at a time.
 * N: the MATERNAL-homozygous informative count on chr16 for the three terminal misses, which is
 *    the only candidate that explains three sibling cells failing identically.
 * O: `absentParent` on the duplication-only segment, where the prediction is MATERNAL: a paternal
 *    duplication puts BAF at 1/3 and miscalls the mother's allele to absent.
 * P: BAF band occupancy on the five uniformly-refused cells. A triploid has no 0.5 band at all;
 *    a mixture keeps one.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/phase2-tests.ts
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




type Res = {
  segments?: { chrom: string; kind: string; startBp: number; endBp: number; score: number
    absentParent?: string }[]
  stage?: { stage: string; why: string }
  chroms?: { chrom: string; verdict: string; lostParent?: string; informative?: number }[]
}
type Stats = { markers: number; called: number; het: number }

async function score2(gsm: string, father: string, mother: string | null, declared?: string) {
  const pat = parentIndex(father)
  const mat = mother ? parentIndex(mother) : null
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  // BAF band occupancy, per sample, at HETEROZYGOUS calls: a triploid's heterozygous sites sit at
  // 1/3 and 2/3 and leave the half band empty, a mixture keeps it. Counted here rather than taken
  // from the profile, because the profile's band is a different interval.
  const baf = { nHet: 0, half: 0, third: 0, twoThird: 0, other: 0 }
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
    const b = (r as { baf?: number }).baf
    if (r.genotype === 'AB' && typeof b === 'number' && Number.isFinite(b)) {
      baf.nHet += 1
      if (b >= 0.45 && b <= 0.55) baf.half += 1
      else if (b >= 0.27 && b < 0.40) baf.third += 1
      else if (b > 0.60 && b <= 0.73) baf.twoThird += 1
      else baf.other += 1
    }
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
    ...(declared ? { declaredStage: declared } : {}),
  } as never) as Res
  return { res, profile: profile as unknown as { callRate: number; hetRate: number
    byChrom: Map<string, Stats> }, acc: acc as never, baf }
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
const hand = new Map<string, string>(); const karyo = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , k, chrom, manual] = l.split('\t')
  hand.set(`${gsm}|${chrom}`, manual)
  if (!karyo.has(gsm)) karyo.set(gsm, k)
}
const cells = [...karyo.keys()].filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()

const stg = await import(`${W}stage.ts`)
const par = await import(`${W}parentage.ts`)

// The unconditioned statistic, as it stood, so the two can be read side by side.
const oldSpread = (b: Map<string, Stats>) => {
  const rates = [...b.entries()].filter(([c, v]) => /^\d+$/.test(c) && v.called >= 500)
    .map(([, v]) => v.het / v.called)
  if (rates.length < 5) return NaN
  const m = rates.reduce((a, x) => a + x, 0) / rates.length
  const sd = Math.sqrt(rates.reduce((a, x) => a + (x - m) ** 2, 0) / rates.length)
  return sd / m
}

const rows: { gsm: string; stage: string; cr: number; het: number; was: number; now: number
  nMeasured: number; nNotMeasured: number; segs: number; k: string
  baf: { nHet: number; half: number; third: number; twoThird: number } }[] = []
for (const g of cells) {
  const p = parentsOf.get(g)!
  const { res, profile, baf } = await score2(g, p.father, p.mother)
  const auto = (res.chroms ?? []).filter((c) => /^\d+$/.test(c.chrom))
  rows.push({
    gsm: g, stage: res.stage?.stage ?? '?', cr: profile.callRate, het: profile.hetRate,
    was: oldSpread(profile.byChrom), now: stg.hetSpread(profile.byChrom as never) as number,
    nMeasured: auto.filter((c) => c.verdict !== 'not_measured').length,
    nNotMeasured: auto.filter((c) => c.verdict === 'not_measured').length,
    segs: (res.segments ?? []).length, k: karyo.get(g) ?? '', baf,
  })
  process.stdout.write('.')
}
console.log('\n')

// --- L and M -----------------------------------------------------------------------------------
console.log('TEST L AND M. THE CONDITIONED SPREAD, AND WHAT THE PER-CHROMOSOME GATE THEN LEAVES')
console.log(`${'cell'.padEnd(12)} ${'stage'.padEnd(14)} ${'callRate'.padStart(8)} `
  + `${'was'.padStart(6)} ${'now'.padStart(6)} ${'measured'.padStart(8)} ${'notMeas'.padStart(7)} `
  + `${'segs'.padStart(5)}  hand`)
const LOC = stg.HET_SPREAD_LOCALISED as number
for (const r of rows.filter((x) => x.stage === 'failed').sort((a, b) => b.now - a.now)) {
  const wasIn = Number.isFinite(r.was) && r.was >= LOC
  const nowIn = Number.isFinite(r.now) && r.now >= LOC
  console.log(`${r.gsm.padEnd(12)} ${r.stage.padEnd(14)} ${r.cr.toFixed(2).padStart(8)} `
    + `${(Number.isFinite(r.was) ? r.was.toFixed(3) : 'NaN').padStart(6)} `
    + `${(Number.isFinite(r.now) ? r.now.toFixed(3) : 'NaN').padStart(6)} `
    + `${String(r.nMeasured).padStart(8)} ${String(r.nNotMeasured).padStart(7)} `
    + `${String(r.segs).padStart(5)}  ${r.k.slice(0, 26)}`
    + `${wasIn && !nowIn ? '   WAS ADMITTED, NOW REFUSED' : ''}`
    + `${!wasIn && nowIn ? '   NEWLY ADMITTED' : ''}${wasIn && nowIn ? '   still admitted' : ''}`)
}
const accepted = rows.filter((x) => x.stage !== 'failed')
const moved = accepted.filter((r) => (Number.isFinite(r.was) && r.was >= LOC)
  !== (Number.isFinite(r.now) && r.now >= LOC))
console.log(`\n  accepted cells whose classification moved, which must be 0: ${moved.length}`)
console.log('  (an accepted cell never reaches this gate, so this is a sanity check on the statistic)')

// --- N -----------------------------------------------------------------------------------------
console.log('\nTEST N. THE chr16 MARKER SETS ON THE THREE TERMINAL MATERNAL MISSES')
console.log('  pat set = markers where the FATHER is homozygous; mat set = where the MOTHER is.')
console.log(`${'cell'.padEnd(12)} ${'pat set'.padStart(8)} ${'mat set'.padStart(8)} `
  + `${'floor'.padStart(6)} ${'verdict'.padEnd(14)} reachable?`)
for (const g of ['GSM4472423', 'GSM4552385', 'GSM4552386']) {
  const p = parentsOf.get(g)!
  const { res, acc } = await score2(g, p.father, p.mother)
  const a = acc as unknown as {
    absenceByChrom: Map<string, unknown[]>; matAbsenceByChrom: Map<string, unknown[]>
  }
  const nPat = a.absenceByChrom.get('16')?.length ?? 0
  const nMat = a.matAbsenceByChrom.get('16')?.length ?? 0
  const v = (res.chroms ?? []).find((c) => c.chrom === '16')?.verdict ?? 'none'
  const floor = segs.MIN_SEGMENT_MARKERS as number
  console.log(`${g.padEnd(12)} ${String(nPat).padStart(8)} ${String(nMat).padStart(8)} `
    + `${String(floor).padStart(6)} ${v.padEnd(14)} `
    + `pat ${nPat >= floor ? 'yes' : 'NO'}, mat ${nMat >= floor ? 'yes' : 'NO'}`)
}

// --- O -----------------------------------------------------------------------------------------
console.log('\nTEST O. THE DUPLICATION-ONLY EVENT REPORTED BY THE ABSENCE CHANNEL')
{
  const g = 'GSM4552365'; const p = parentsOf.get(g)!
  const { res } = await score2(g, p.father, p.mother)
  console.log(`  hand: ${hand.get(`${g}|16`)}`)
  for (const sg of (res.segments ?? []).filter((s) => s.chrom === '16')) {
    console.log(`  ${sg.kind.padEnd(18)} ${(sg.startBp / 1e6).toFixed(1)}-${(sg.endBp / 1e6).toFixed(1)}Mb `
      + `LRT ${sg.score.toFixed(0).padStart(6)}  absentParent ${sg.absentParent ?? '-'}`)
  }
  console.log('  PREDICTION: a parental-absence segment here reads MATERNAL and spans 16q21-tel,')
  console.log('  57.4-90.4Mb. If it reads paternal, the review\'s mechanism is wrong.')
}

// --- P -----------------------------------------------------------------------------------------
console.log('\nTEST P. BAF BANDS AT HETEROZYGOUS CALLS. A TRIPLOID HAS NO HALF BAND.')
console.log(`${'cell'.padEnd(12)} ${'het'.padStart(6)} ${'nHet'.padStart(8)} `
  + `${'0.5 band'.padStart(9)} ${'1/3'.padStart(7)} ${'2/3'.padStart(7)}  hand`)
const uniform = ['GSM4472422', 'GSM4552387', 'GSM4552388', 'GSM4552389', 'GSM4472421',
  'GSM4552418', 'GSM4552422', 'GSM4472409', 'GSM4472405']
for (const r of rows.filter((x) => uniform.includes(x.gsm))) {
  const pc = (n: number) => `${((100 * n) / Math.max(1, r.baf.nHet)).toFixed(1)}%`
  console.log(`${r.gsm.padEnd(12)} ${r.het.toFixed(3).padStart(6)} `
    + `${r.baf.nHet.toLocaleString().padStart(8)} ${pc(r.baf.half).padStart(9)} `
    + `${pc(r.baf.third).padStart(7)} ${pc(r.baf.twoThird).padStart(7)}  ${r.k.slice(0, 24)}`)
}
const ctrl = rows.filter((x) => x.stage !== 'failed' && /normal|46,X/.test(x.k)).slice(0, 4)
console.log('  controls, cells the stage inference accepts:')
for (const r of ctrl) {
  const pc = (n: number) => `${((100 * n) / Math.max(1, r.baf.nHet)).toFixed(1)}%`
  console.log(`${r.gsm.padEnd(12)} ${r.het.toFixed(3).padStart(6)} `
    + `${r.baf.nHet.toLocaleString().padStart(8)} ${pc(r.baf.half).padStart(9)} `
    + `${pc(r.baf.third).padStart(7)} ${pc(r.baf.twoThird).padStart(7)}  ${r.k.slice(0, 24)}`)
}
void par

/**
 * Which parent supplied an extra whole chromosome, on real cells and on gains built into real cells.
 *
 * TWO ARMS, because neither alone can carry the claim.
 *
 * REAL EVENTS. Zuccaro 2020 Table S3 (audit/zuccaro-table-s3.tsv) records, for the 31 GSE148488
 * cells that have both parents arrayed, every chromosome the authors read as a gain or duplication
 * and whose parent it was. This runs the shipped pipeline on those cells with both parents loaded
 * and reports what the allele gain scan says on each such chromosome. It is eight events, five of
 * them arm-sized, so it can show that the direction is right and cannot show how often.
 *
 * AND THE HAND KARYOTYPE IS NOT A COMPLETE NEGATIVE SET FOR GAINS. The authors tabulated the
 * chromosomes they were looking for, which were losses. A chromosome the scan calls a gain that the
 * table lists as normal is therefore NOT counted as a false positive here. It is listed, with
 * whether sibling cells of the same embryo agree, and the reader decides. Counting them as errors
 * would score the scan against a truth that was never built to hold it.
 *
 * CONSTRUCTED GAINS. A real cell with both parents, one chromosome's allele intensities
 * re-weighted so the named parent's allele is doubled (BAF' from A' = 2A or B' = 2B, the
 * construction audit/gain-positive-control.ts uses), the log2R raised by 0.3, and the rest of the
 * cell untouched. The shipped pipeline is run on the result, so the noise, the dropout, the
 * reference chromosomes and every gate are the real ones. The truth is the construction. Parent
 * accuracy and detection rate are reported by how clean the cell's own scale is, since a noisy
 * cell is supposed to refuse, and a refusal is not an error.
 *
 * AND THE ARRAYS THE STAGE INFERENCE REFUSES. A cell with several gains reads far more heterozygous
 * than one diploid can, and the inference calls it "no genome reads this way" and withholds every
 * reading on it, the allele scan included. Four of the 31 cells are refused, and they hold the two
 * clearest real gains and 23 of the hand-called losses. A second pass reads those cells with the
 * stage DECLARED single-cell, the operator override the tool already has, so the table shows what
 * would be reported if the refusal were lifted. That is the decision the chaotic-cell question is
 * about, and it is made here with numbers instead of by assertion. Nothing in the shipped tool
 * changes by it.
 *
 * Run: OM_TRIOS=<dir with GSM*.probes.gz and trio_manifest_full.csv> [OM_SHARD=i/n] [OM_CONSTRUCT=1]
 *        [OM_REAL=0] [OM_CHROMS=3,9,15] [OM_OUT=out.json] \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/gain-origin-truth.ts
 */
import { readFileSync, existsSync, readdirSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

// MERGE MODE: combine the constructed arm of several shards into one table and stop. The shards
// each build into a quarter or a third of the cells, and a table per shard says little.
if (process.env.OM_MERGE) {
  const all: BuiltRow[] = []
  for (const f of process.env.OM_MERGE.split(',')) all.push(...(JSON.parse(readFileSync(f, 'utf8')).built as BuiltRow[]))
  console.log(`merged ${process.env.OM_MERGE.split(',').length} shards: ${all.length} constructed gains over ${new Set(all.map((b) => b.gsm)).size} cells`)
  summariseBuilt(all)
  process.exit(0)
}

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('gain-origin-truth: OM_TRIOS is not a readable directory. NOT RUN.')
  process.exit(0)
}
const [SHARD_I, SHARD_N] = (process.env.OM_SHARD ?? '0/1').split('/').map(Number)
const CONSTRUCT = process.env.OM_CONSTRUCT === '1'
const CHROMS = (process.env.OM_CHROMS ?? '3,9,15').split(',')
const REAL = process.env.OM_REAL !== '0'

const fileOf = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const m = /^(GSM\d+)/.exec(f)
  if (m && (f.endsWith('.gz') || f.endsWith('.probes') || f.endsWith('.txt'))) fileOf.set(m[1], join(DIR, f))
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
interface Truth { gsm: string; chrom: string; manual: string }
const truth: Truth[] = readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8')
  .split('\n').filter((l) => l && !l.startsWith('#') && !l.startsWith('gsm'))
  .map((l) => { const [gsm, , , chrom, manual] = l.split('\t'); return { gsm, chrom, manual } })
const handOf = new Map(truth.map((t) => [`${t.gsm}|${t.chrom}`, t.manual]))
const cells = [...new Set(truth.map((t) => t.gsm))].filter((g) => parentsOf.has(g) && fileOf.has(g))
  .sort().filter((_, i) => i % SHARD_N === SHARD_I)

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
const parentCache = new Map<string, { gt: Map<string, string> }>()
function parentIndex(gsm: string) {
  const hit = parentCache.get(gsm)
  if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  parentCache.set(gsm, idx as never)
  return idx as unknown as { gt: Map<string, string> }
}

type Build = { chrom: string; extra: 'paternal' | 'maternal' } | null
async function scoreWith(gsm: string, father: string, mother: string, build: Build = null, declared?: string) {
  const pat = parentIndex(father); const mat = parentIndex(mother)
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (build && r.chrom === build.chrom) {
      // THE CONSTRUCTION. At a marker where the parents are opposite homozygotes the named parent's
      // allele intensity is doubled, which is what a third copy of that parent's chromosome does.
      // Elsewhere on the chromosome nothing changes the share, so only those markers are touched.
      const fa = pat.gt.get(r.probesetId); const mo = mat.gt.get(r.probesetId)
      if (r.baf !== null && ((fa === 'AA' && mo === 'BB') || (fa === 'BB' && mo === 'AA'))) {
        let A = 1 - r.baf; let B = r.baf
        const fatherAllele = fa === 'AA' ? 'A' : 'B'
        const doubled = build.extra === 'paternal' ? fatherAllele : (fatherAllele === 'A' ? 'B' : 'A')
        if (doubled === 'A') A *= 2; else B *= 2
        r.baf = B / (A + B)
      }
      if (r.log2R !== null) r.log2R += 0.3
    }
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  return score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
    ...(declared ? { declaredStage: declared } : {}),
  }) as Promise<Record<string, unknown>>
}

interface Chrom {
  chrom: string; aneuploidy?: string; lostParent?: string
  alleleGain?: { verdict: string; parent?: string; z: number; shift: number; scale: number; lrrShift: number }
}
const gainOf = (c: Chrom | undefined) => c?.alleleGain

interface BuiltRow { gsm: string; chrom: string; extra: string; verdict: string; parent?: string; scale: number }
function summariseBuilt(built: BuiltRow[]) {
  const classes: [string, (s: number) => boolean][] = [
    ['clean cells, own scale under 0.03', (s) => s < 0.03],
    ['noisy cells, own scale 0.03 to 0.10', (s) => s >= 0.03 && s < 0.10],
    ['chaotic cells, own scale over 0.10', (s) => s >= 0.10],
    ['no scale (too few biparental chromosomes)', (s) => !Number.isFinite(s)],
  ]
  console.log('\nCONSTRUCTED GAINS: one chromosome, one parent, built into a real two-parent cell')
  for (const [name, test] of classes) {
    const b = built.filter((x) => test(x.scale))
    if (!b.length) continue
    const g = b.filter((x) => x.verdict === 'gain')
    const right = g.filter((x) => x.parent === x.extra).length
    console.log(`  ${name.padEnd(44)} ${String(b.length).padStart(4)} built: gain called ${String(g.length).padStart(4)} (${(100 * g.length / b.length).toFixed(0)}%), right parent ${right} of ${g.length}, `
      + `imbalance ${b.filter((x) => x.verdict === 'imbalance').length}, none ${b.filter((x) => x.verdict === 'none').length}, unscored ${b.filter((x) => x.verdict === 'unscored').length}`)
  }
  const all = built.filter((x) => x.verdict === 'gain')
  console.log(`  ALL ${built.length} built: gain called ${all.length}, right parent ${all.filter((x) => x.parent === x.extra).length}, WRONG parent ${all.filter((x) => x.parent !== x.extra).length}`)
}

// --- real events ---------------------------------------------------------------------------------
interface RealRow {
  gsm: string; chrom: string; hand: string; handParent?: string; arm: boolean
  verdict: string; parent?: string; z?: number; shift?: number; lrr?: number; scale?: number; stage: string
}
interface Loss { gsm: string; chrom: string; hand: string; handParent: string; lost?: string }
const scaleOf = new Map<string, number>()
const refusedCells: string[] = []

async function realArm(title: string, declared?: string, only?: string[]) {
  const real: RealRow[] = []
  const losses: Loss[] = []
  const unlabelled: { gsm: string; chrom: string; verdict: string; parent?: string; z: number; shift: number; lrr: number }[] = []
  for (const gsm of (only ?? cells)) {
    const p = parentsOf.get(gsm)!
    const res = await scoreWith(gsm, p.father, p.mother, null, declared)
    const chroms = (res.chroms ?? []) as Chrom[]
    const stageName = String((res.stage as { stage?: string })?.stage ?? '')
    if (!declared) {
      if (stageName === 'failed') refusedCells.push(gsm)
      const scales = chroms.map((c) => gainOf(c)?.scale).filter((x): x is number => Number.isFinite(x as number))
      scaleOf.set(gsm, scales.length ? scales.sort((a, b) => a - b)[scales.length >> 1] : NaN)
    }
    for (const c of chroms) {
      const hand = handOf.get(`${gsm}|${c.chrom}`) ?? ''
      const g = gainOf(c)
      const isGain = /tris|gain|dup/i.test(hand)
      const parent = /\(pat\)/.test(hand) ? 'paternal' : /\(mat\)/.test(hand) ? 'maternal' : undefined
      const mono = /^monosomy \d+ \((pat|mat)\)/.exec(hand)
      if (mono) losses.push({ gsm, chrom: c.chrom, hand, handParent: mono[1] === 'pat' ? 'paternal' : 'maternal', lost: c.lostParent })
      if (isGain) {
        real.push({ gsm, chrom: c.chrom, hand, handParent: parent, arm: /dup|del/i.test(hand) && !/gain|tris/i.test(hand),
          verdict: g?.verdict ?? 'unscored', parent: g?.parent, z: g?.z, shift: g?.shift, lrr: g?.lrrShift, scale: g?.scale, stage: stageName })
      } else if (g && (g.verdict === 'gain' || g.verdict === 'imbalance') && hand === 'normal') {
        unlabelled.push({ gsm, chrom: c.chrom, verdict: g.verdict, parent: g.parent, z: g.z, shift: g.shift, lrr: g.lrrShift })
      }
    }
    process.stdout.write('.')
  }
  console.log()
  console.log(`\n${title}: chromosomes the hand karyotype reads as a gain or duplication, ${(only ?? cells).length} cells`)
  for (const r of real) {
    const call = r.verdict === 'gain' ? `GAIN ${r.parent}` : r.verdict
    const ok = r.verdict === 'gain' ? (r.parent === r.handParent ? 'right parent' : 'WRONG parent') : 'not called'
    console.log(`  ${r.gsm} chr${r.chrom.padStart(2)} ${r.arm ? 'arm-sized ' : 'whole     '} hand ${(r.handParent ?? '?').padEnd(8)} | ${call.padEnd(16)} ${ok.padEnd(13)} `
      + `z ${Number.isFinite(r.z as number) ? (r.z as number).toFixed(1).padStart(5) : '    -'} shift ${Number.isFinite(r.shift as number) ? (r.shift as number).toFixed(3).padStart(6) : '     -'} `
      + `lrr ${Number.isFinite(r.lrr as number) ? (r.lrr as number).toFixed(2).padStart(5) : '    -'} scale ${Number.isFinite(r.scale as number) ? (r.scale as number).toFixed(3) : '-'} stage ${r.stage}`
      + `  [${r.hand.slice(0, 34)}]`)
  }
  const named = real.filter((r) => r.verdict === 'gain')
  console.log(`  summary: ${real.length} events, called ${named.length}, right parent ${named.filter((r) => r.parent === r.handParent).length}, wrong ${named.filter((r) => r.parent !== r.handParent).length}`)
  const lc = losses.filter((l) => l.lost)
  console.log(`  hand-called monosomies with a parent: ${losses.length}, named by the alleles ${lc.length}, right parent ${lc.filter((l) => l.lost === l.handParent).length}, wrong ${lc.filter((l) => l.lost !== l.handParent).length}`)
  console.log(`  called gains and imbalances on chromosomes the hand karyotype lists as normal (not counted as errors, see header): ${unlabelled.length}`)
  for (const u of unlabelled.sort((a, b) => a.chrom.localeCompare(b.chrom, undefined, { numeric: true }))) {
    console.log(`    ${u.gsm} chr${u.chrom.padStart(2)} ${u.verdict.padEnd(9)} ${(u.parent ?? '-').padEnd(9)} z ${u.z.toFixed(1).padStart(6)} shift ${u.shift.toFixed(3).padStart(6)} lrr ${u.lrr.toFixed(2).padStart(5)}`)
  }
  return { real, losses, unlabelled }
}
let realOut: Awaited<ReturnType<typeof realArm>> | null = null
let declaredOut: Awaited<ReturnType<typeof realArm>> | null = null
if (REAL) {
  realOut = await realArm('REAL EVENTS, AS SHIPPED')
  if (refusedCells.length) {
    console.log(`\nthe stage inference refuses ${refusedCells.length} of ${cells.length} cells: ${refusedCells.join(', ')}`)
    declaredOut = await realArm('REAL EVENTS ON THE REFUSED ARRAYS, STAGE DECLARED single-cell', 'single-cell', refusedCells)
  }
}

// --- constructed gains ---------------------------------------------------------------------------
interface Built { gsm: string; chrom: string; extra: string; verdict: string; parent?: string; scale: number }
const built: Built[] = []
const skipped = { handEvent: 0, alreadyShifted: 0, unscorable: 0 }
if (CONSTRUCT) {
  for (const gsm of cells) {
    const p = parentsOf.get(gsm)!
    // THE UNTOUCHED CELL FIRST. A gain built into a chromosome that already carries an event is not
    // a gain built into a euploid one, and what the scan then names depends on both. So a chromosome
    // is built into only if the hand karyotype lists it as normal AND the untouched scan reads it as
    // none. One that cannot be scored at all (too few markers, too few reference chromosomes) is
    // counted and skipped: building into it can only return the same refusal.
    const natural = await scoreWith(gsm, p.father, p.mother)
    const nat = (natural.chroms ?? []) as Chrom[]
    const ss = nat.map((x) => gainOf(x)?.scale).filter((x): x is number => Number.isFinite(x as number)).sort((a, b) => a - b)
    const scale = ss.length ? ss[ss.length >> 1] : NaN
    for (const chrom of CHROMS) {
      const g0 = gainOf(nat.find((x) => x.chrom === chrom))
      if ((handOf.get(`${gsm}|${chrom}`) ?? '') !== 'normal') { skipped.handEvent += 1; continue }
      if (g0 && (g0.verdict === 'gain' || g0.verdict === 'imbalance')) { skipped.alreadyShifted += 1; continue }
      if (!g0 || g0.verdict === 'unscored') { skipped.unscorable += 1; continue }
      for (const extra of ['paternal', 'maternal'] as const) {
        const res = await scoreWith(gsm, p.father, p.mother, { chrom, extra })
        const c = ((res.chroms ?? []) as Chrom[]).find((x) => x.chrom === chrom)
        const g = gainOf(c)
        built.push({ gsm, chrom, extra, verdict: g?.verdict ?? 'unscored', parent: g?.parent, scale })
        process.stdout.write('+')
      }
    }
  }
  console.log()
  console.log(`  chromosomes not built into: ${skipped.handEvent} the hand karyotype lists an event on, ${skipped.alreadyShifted} the untouched scan already calls, ${skipped.unscorable} that cannot be scored in that cell`)
  summariseBuilt(built)
}
if (process.env.OM_OUT) writeFileSync(process.env.OM_OUT, JSON.stringify({ realOut, declaredOut, built, scaleOf: [...scaleOf] }, null, 1))

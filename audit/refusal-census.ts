/**
 * The refusal census: every claim the shipped tool withheld, scored against a published karyotype.
 *
 * WHY THIS EXISTS. Refusal in this tool was "fixed" about eight times, each time by checking that a
 * threshold was right on its own terms or by giving it a better sentence, and the over-refusal kept
 * returning because nothing ever counted what the refusals threw away. This counts it. Every
 * withheld claim is read from `result.refusals`, where the shipped `scoreSample` records it with the
 * value it would have reported, and scored against Zuccaro 2020 Table S3, a per-chromosome karyotype
 * published and reviewed for every arrayed cell of GSE148488.
 *
 * IT RUNS THE PRODUCT, NOT A COPY OF IT. A harness that re-implemented part of the pipeline scanned
 * arrays the product refuses, and produced every gains figure this project recorded until
 * 10 November. This one calls `scoreSample` exactly as the application does, on the same rows.
 *
 * SCORED BY CLAIM TYPE AND DIRECTION, as the over-refusal consult specified, because normal
 * chromosomes outnumber abnormal ones about four to one and a pooled count rewards saying "normal".
 * Wrong parent is its own column. A published "normal" contradicting a claim on a cell published
 * chaotic is reported apart, since the consult found most such chromosomes carry the signature of
 * no template and the published reading may be the wrong one; those are to be adjudicated, not
 * counted either way.
 *
 * Run: OM_TRIOS=<dir with trio_manifest_full.csv and the arrays> \
 *        node --experimental-strip-types --max-old-space-size=4096 audit/refusal-census.ts
 */
import { readFileSync, existsSync, readdirSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'
import { lineReader } from './joy-common.ts'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('refusal-census: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

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

/** gsm -> chromosome -> the table's word for it; and gsm -> whole-cell karyotype. */
const truth = new Map<string, Map<string, string>>()
const karyotype = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , kary, chrom, manual] = l.split('\t')
  if (!gsm) continue
  if (!truth.has(gsm)) truth.set(gsm, new Map())
  truth.get(gsm)!.set(chrom.replace(/^chr/i, ''), (manual ?? '').trim().toLowerCase())
  karyotype.set(gsm, (kary ?? '').trim())
}

function* rowsOf(path: string) {
  const raw = readFileSync(path)
  const next = lineReader(path.endsWith('.gz') ? gunzipSync(raw) : raw)
  const head: string[] = []
  for (let i = 0; i < 60; i += 1) { const l = next(); if (l === null) break; head.push(l) }
  const h = head.findIndex((l) => l && !l.startsWith('#'))
  const map = ingest.headerMap(head[h])
  if (!map) throw new Error(`${path}: header not recognised`)
  for (const l of head.slice(h + 1)) { const r = ingest.parseRow(l, map); if (r) yield r }
  for (let l = next(); l !== null; l = next()) { const r = ingest.parseRow(l, map); if (r) yield r }
}
const parentCache = new Map<string, unknown>()
function parentIndex(gsm: string) {
  if (parentCache.has(gsm)) return parentCache.get(gsm)
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never,
    ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  parentCache.set(gsm, idx)
  return idx
}
async function scoreCell(gsm: string, father: string, mother: string) {
  const pat = parentIndex(father); const mat = parentIndex(mother)
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  return score.scoreSample({ acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm,
    log: () => {} }) as Promise<Record<string, any>>
}

type Score = 'confirmed' | 'contradicted' | 'wrong-parent' | 'adjudicate' | 'unscorable'
type Direction = 'abnormal' | 'normal' | 'none'

/** A claim against the table's word for that chromosome. Returns direction and verdict. */
function judge(claim: string, value: any, t: string | undefined): [Direction, Score] {
  if (t === undefined || t === '' || t === 'no call') return ['none', 'unscorable']
  const normal = t === 'normal'
  const mono = /^monosomy \d+ \((pat|mat)\)/.exec(t)
  const lostPar = mono ? (mono[1] === 'pat' ? 'paternal' : 'maternal') : undefined
  const lossLike = !!mono || /nullisom|del|only/.test(t)
  const gainLike = /trisom|dup/.test(t)
  const loss = (v: boolean): Score => (normal ? 'contradicted' : lossLike ? (v ? 'confirmed' : 'contradicted')
    : gainLike ? (v ? 'contradicted' : 'confirmed') : 'adjudicate')
  switch (claim) {
    case 'aneuploidy':
      return ['abnormal', value === 'loss' ? loss(true) : loss(false)]
    case 'lostParent':
    case 'aneuploidyParent': {
      const par = claim === 'lostParent' ? value : (value === 'this' ? 'paternal' : 'maternal')
      if (normal) return ['abnormal', 'contradicted']
      if (lostPar) return ['abnormal', par === lostPar ? 'confirmed' : 'wrong-parent']
      return ['abnormal', 'adjudicate']
    }
    case 'alleleGain': {
      const v = value?.verdict
      if (v === 'gain') return ['abnormal', normal ? 'contradicted' : gainLike ? 'confirmed' : 'adjudicate']
      if (v === 'none') return ['normal', normal ? 'confirmed' : (lossLike || gainLike) ? 'contradicted' : 'adjudicate']
      return ['none', 'unscorable']
    }
    case 'segment': {
      const k = value?.kind as string | undefined
      if (!k) return ['none', 'unscorable']
      return ['abnormal', loss(k !== 'copy-gain')]
    }
    default:
      return ['none', 'unscorable']
  }
}
const brief = (v: any): string => (v === undefined ? '' : typeof v === 'object'
  ? (v.verdict ?? v.kind ?? JSON.stringify(v).slice(0, 40)) : String(v))

const SHARDS = Number(process.env.OM_SHARDS ?? 1)
const SHARD = Number(process.env.OM_SHARD ?? 0)
const cells = [...truth.keys()].filter((g) => fileOf.has(g) && parentsOf.has(g)).sort()
  .filter((_, i) => i % SHARDS === SHARD)
console.log(`${cells.length} published cells with both parents on disk (shard ${SHARD} of ${SHARDS})`)

const out: string[] = [['kind', 'gsm', 'stage', 'refused_by', 'gate', 'claim', 'chrom', 'value',
  'published', 'cell_karyotype', 'direction', 'score', 'heterozygosity', 'call_rate'].join('\t')]
for (const gsm of cells) {
  const p = parentsOf.get(gsm)!
  let res
  try { res = await scoreCell(gsm, p.father, p.mother) } catch (e) {
    console.log(`  ${gsm}: ${(e as Error).message}`); continue
  }
  const t = truth.get(gsm)!
  const kary = karyotype.get(gsm) ?? ''
  const stg = res.stage?.stage ?? 'unknown'
  const by = res.stage?.refusedBy ?? ''
  const emit = (kind: string, gate: string, claim: string, chrom: string, value: any, measured: any) => {
    const [dir, sc] = judge(claim, value, t.get(chrom))
    out.push([kind, gsm, stg, by, gate, claim, chrom, brief(value), t.get(chrom) ?? '', kary, dir, sc,
      Number(measured?.heterozygosity ?? NaN).toFixed(3), Number(measured?.callRate ?? NaN).toFixed(3)].join('\t'))
  }
  for (const r of (res.refusals ?? [])) emit(r.soft ? 'soft' : 'withheld', r.gate, r.claim, r.chrom ?? '', r.withheld, r.measured)
  // AND WHAT IT MADE, scored the same way, so over-refusal and over-calling sit side by side.
  const m = { heterozygosity: res.stage ? NaN : NaN, callRate: NaN }
  for (const c of (res.chroms ?? [])) {
    for (const claim of ['aneuploidy', 'lostParent', 'aneuploidyParent', 'alleleGain']) {
      if (c[claim] !== undefined) emit('made', '', claim, c.chrom, c[claim], m)
    }
  }
  for (const sg of (res.segments ?? [])) emit('made', '', 'segment', sg.chrom, sg, m)
}
if (process.env.OM_TSV) writeFileSync(process.env.OM_TSV, out.join('\n') + '\n')
console.log(`${out.length - 1} claims recorded${process.env.OM_TSV ? `, written to ${process.env.OM_TSV}` : ''}`)

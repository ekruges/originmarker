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
const inferredArray = await import(`${W}inferredArray.ts`)
const answerKeys = await import(`${W}answerKeys.ts`)

/**
 * THE LABORATORY'S SETTING. OM_FATHER names a parent array, usually one reconstructed from the
 * children, and the cell is scored against it ALONE, with no mother: the configuration every
 * laboratory run actually has. Unset, both measured parents are loaded, which is the reference the
 * one-parent result is compared against.
 *
 * OM_PANEL_DIR holds founder panels, panel-val-<mother>.tsv, built WITHOUT the cell's own parents,
 * so the parent-free answer keys read against it owe nothing to the family.
 */
const FATHER = process.env.OM_FATHER
const PANEL_DIR = process.env.OM_PANEL_DIR
const MOTHER_LETTER: Record<string, string> = { GSM4472407: 'A', GSM4472417: 'B', GSM4472418: 'B',
  GSM4472419: 'B', GSM4472420: 'B', GSM4472415: 'C', GSM4472416: 'C', GSM4472408: 'D' }
const panelCache = new Map<string, Map<string, 'AA' | 'BB'>>()
/** OM_PANEL_FILE, where set, is one panel for every cell: the shipped one, used to check that the
 *  decisions it makes match those of the per-family panels it was validated with. */
const PANEL_FILE = process.env.OM_PANEL_FILE
function panelFor(mother: string): Map<string, 'AA' | 'BB'> | null {
  if (PANEL_FILE) {
    // Read by the product's own parser, so the census reads a panel as the app does.
    if (!panelCache.has('file')) {
      const raw = readFileSync(PANEL_FILE)
      panelCache.set('file', answerKeys.parsePanel(
        (PANEL_FILE.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')) as never)
    }
    return panelCache.get('file')!
  }
  const letter = MOTHER_LETTER[mother]
  if (!PANEL_DIR || !letter) return null
  if (!panelCache.has(letter)) {
    const m = new Map<string, 'AA' | 'BB'>()
    for (const l of readFileSync(join(PANEL_DIR, `panel-val-${letter}.tsv`), 'utf8').split('\n').slice(1)) {
      const [id, , allele] = l.split('\t')
      if (id) m.set(id, allele === 'A' ? 'AA' : 'BB')
    }
    panelCache.set(letter, m)
  }
  return panelCache.get(letter)!
}

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
  return parentAt(fileOf.get(gsm)!)
}
/** A parent array from a path, reconstructed or measured, read as cli/om.ts loadParent reads it:
 *  the banner at the head of a reconstructed file is what marks it as one. */
function parentAt(path: string) {
  if (parentCache.has(path)) return parentCache.get(path)
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  const raw = readFileSync(path)
  const banner = inferredArray.parseInferredBanner(
    (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8', 0, inferredArray.MARK_BYTES)) ?? undefined
  for (const r of rowsOf(path)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never,
    ingest.finishProfile(path, byChrom as never, bafSums as never, first).build.build, banner)
  parentCache.set(path, idx)
  return idx
}
async function scoreCell(gsm: string, father: string, mother: string) {
  const pat = FATHER ? parentAt(FATHER) : parentIndex(father)
  const mat = FATHER ? null : parentIndex(mother)
  const acc = score.emptyCollected(pat as never, mat as never, panelFor(mother) as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); const builds = ingest.emptyBuildSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    ingest.accumulateBuild(r as never, builds)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first, builds)
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
const only = process.env.OM_ONLY ? new Set(process.env.OM_ONLY.split(',')) : null
const cells = [...truth.keys()].filter((g) => fileOf.has(g) && parentsOf.has(g) && (!only || only.has(g))).sort()
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
      Number(measured?.heterozygosity ?? NaN).toFixed(3), Number(measured?.callRate ?? measured?.autosomalCallRate ?? NaN).toFixed(3)].join('\t'))
  }
  // ONE ROW PER CELL with what gates the allele channel, so a hidden gate shows up as a number.
  const k = res.answerKeys
  const f = (x: number) => (Number.isFinite(x) ? x.toFixed(4) : 'none')
  for (const c of (res.chroms ?? [])) {
    if (!c.copyState) continue
    const cs = c.copyState
    out.push(['state', gsm, stg, by, '', 'copyState', c.chrom, cs.state, t.get(c.chrom) ?? '', kary,
      'none', 'unscorable', Number.isFinite(cs.r) ? cs.r.toFixed(3) : 'nan', Number.isFinite(cs.base) ? cs.base.toFixed(3) : 'nan'].join('\t'))
  }
  out.push(['cell', gsm, stg, by, `zygosity=${res.zygosity}`, 'cell', '',
    `dropIn parents=${f(k?.dropIn.parents)} panel=${f(k?.dropIn.panel)} `
      + `information parents=${f(k?.information.parents)} panel=${f(k?.information.panel)} `
      + `spread=${f(k?.heterozygositySpread)} normal=${f(k?.normalFraction ?? NaN)}`,
    '', kary, 'none', 'unscorable', '', ''].join('\t'))
  for (const r of (res.refusals ?? [])) emit(r.soft ? 'soft' : 'withheld', r.gate, r.claim, r.chrom ?? '', r.withheld, r.measured)
  // AND WHAT IT MADE, scored the same way, so over-refusal and over-calling sit side by side.
  const m = { heterozygosity: res.stage ? NaN : NaN, callRate: NaN }
  for (const c of (res.chroms ?? [])) {
    for (const claim of ['aneuploidy', 'lostParent', 'aneuploidyParent', 'alleleGain', 'mixed']) {
      // For a made claim the measurement columns carry the chromosome's own log R shift, so an
      // allele call can be read against the intensity channel without a second run.
      if (c[claim] !== undefined) emit('made', '', claim, c.chrom, c[claim], { heterozygosity: c.lrrShift, callRate: c.callFraction })
    }
  }
  for (const sg of (res.segments ?? [])) emit('made', '', 'segment', sg.chrom, sg, m)
  // AND WHAT IT LEFT UNDECIDED. A chromosome whose contribution is 'uncalled' or withheld names no
  // parent and is recorded nowhere as refused, so it is the refusal no gate admits to. Read off the
  // product's own result, with the numbers its verdict was made on.
  for (const c of (res.chroms ?? [])) {
    if (!/^\d+$/.test(c.chrom)) continue
    const k = c.contribution
    const why = c.contributionWithheld
    if (!why && k?.ploidy !== 'uncalled') continue
    const tv = t.get(c.chrom) ?? ''
    out.push(['undecided', gsm, stg, by, why ? 'contribution.withheld' : 'contribution.uncalled',
      'contribution', c.chrom, why ? 'withheld' : 'uncalled', tv, kary,
      /^monosomy/.test(tv) ? 'abnormal' : tv === 'normal' ? 'normal' : 'none',
      tv === '' || tv === 'no call' ? 'unscorable' : 'undecided',
      Number(k?.fraction ?? NaN).toFixed(3), Number(k?.informative ?? NaN).toFixed(0),
      Number(k?.onlyA ?? NaN).toFixed(3), Number(k?.onlyB ?? NaN).toFixed(3)].join('\t'))
  }
}
if (process.env.OM_TSV) writeFileSync(process.env.OM_TSV, out.join('\n') + '\n')
console.log(`${out.length - 1} claims recorded${process.env.OM_TSV ? `, written to ${process.env.OM_TSV}` : ''}`)

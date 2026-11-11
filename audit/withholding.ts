/**
 * Withhold arrays, rebuild, and compare with the measured parent: parts A and B of
 * audit/withholding-preregistration.txt. Part C, where each laboratory cell is withheld and the
 * shipped command rebuilds from the rest, is audit/withholding-c.sh.
 *
 * Every statistic is the shipped one: judgeProducts and ProductSet.build for the rebuild, addShare
 * and shareCall for the copy state, the shipped founder panel for both. Arrays are read once into
 * typed arrays over one probe index, so the hundreds of rebuilds and the cells scored against each
 * cost one pass over the data, not one each.
 *
 * Run: OM_TRIOS=<dir> OM_S23=<dir> OM_S23_SERIES=<series number> OM_OUT=<dir> \
 *        node --experimental-strip-types --max-old-space-size=6000 audit/withholding.ts
 */
import { readFileSync, readdirSync, existsSync, writeFileSync, mkdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'
import { lineReader } from './joy-common.ts'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const inferred = await import(`${W}inferredReference.ts`)
const keys = await import(`${W}answerKeys.ts`)
const ops = await import(`${W}oneParentState.ts`)
const { isAutosome } = await import(`${W}parentage.ts`)

const TRIOS = process.env.OM_TRIOS
const S23 = process.env.OM_S23
const OUT = process.env.OM_OUT ?? '.'
if (!TRIOS || !existsSync(TRIOS)) { console.log('OM_TRIOS not readable. NOT RUN.'); process.exit(0) }
mkdirSync(OUT, { recursive: true })

const panel = keys.parsePanel(gunzipSync(readFileSync(
  new URL(`../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url))).toString('utf8'))
const freqB = panel.freqB!
const mono = (id: string, g: 'AA' | 'BB'): boolean => panel.get(id) === g

// --- one probe index for every array ------------------------------------------------------------
const CAP = 1_400_000
const idx = new Map<string, number>()
const probe: string[] = []
const chromNum = new Uint8Array(CAP)
const NONE = 0, AA = 1, BB = 2, AB = 3
interface Arr { id: string; gt: Uint8Array; baf: Float32Array }

function* rows(path: string) {
  const raw = readFileSync(path)
  const next = lineReader(path.endsWith('.gz') ? gunzipSync(raw) : raw)
  let map: ReturnType<typeof ingest.headerMap> = null
  for (let l = next(); l !== null; l = next()) {
    if (!map) { if (l && !l.startsWith('#')) map = ingest.headerMap(l); continue }
    const r = ingest.parseRow(l, map)
    if (r) yield r
  }
}

/** Read an array, into the product set when one is given (a truth array never is). */
function load(path: string, id: string, ps?: InstanceType<typeof inferred.ProductSet>): Arr {
  const gt = new Uint8Array(CAP)
  const baf = new Float32Array(CAP).fill(NaN)
  const slot = ps?.begin(id)
  const band = { inBand: 0, total: 0 }
  for (const r of rows(path)) {
    if (ps) ps.add(slot!, r as never, band)
    if (!isAutosome(r.chrom)) continue
    let i = idx.get(r.probesetId)
    if (i === undefined) {
      i = idx.size
      if (i >= CAP) throw new Error('probe index over capacity')
      idx.set(r.probesetId, i); probe.push(r.probesetId); chromNum[i] = Number(r.chrom)
    }
    gt[i] = r.genotype === 'AA' ? AA : r.genotype === 'BB' ? BB : r.genotype === 'AB' ? AB : NONE
    if (r.baf !== null && Number.isFinite(r.baf)) baf[i] = r.baf
  }
  if (ps) ps.end(slot!, band)
  return { id, gt, baf }
}

const fileIn = (dir: string, prefix: string): string => {
  const f = readdirSync(dir).find((x) => x.startsWith(prefix) && /\.probes(\.gz)?$/.test(x))
  if (!f) throw new Error(`${prefix} not in ${dir}`)
  return join(dir, f)
}

/** The rebuild's homozygous calls, against the truth array of the same person. */
function measured(genotype: Map<string, string>, truth: Arr) {
  let n = 0, het = 0, opp = 0
  for (const [id, g] of genotype) {
    const t = truth.gt[idx.get(id) ?? -1]
    if (!t) continue
    n += 1
    if (t === AB) het += 1
    else if ((t === AA) !== (g === 'AA')) opp += 1
  }
  return { n, contamination: het / n, opposite: opp / n }
}

function* subsets<T>(xs: T[], min: number): Generator<T[]> {
  for (let mask = 1; mask < 1 << xs.length; mask += 1) {
    const s = xs.filter((_, i) => mask & (1 << i))
    if (s.length >= min) yield s
  }
}

const lines: string[] = []
const say = (s: string) => { console.log(s); lines.push(s) }

// --- A: products withheld -----------------------------------------------------------------------
function partA(label: string, ps: InstanceType<typeof inferred.ProductSet>, paternal: string[],
  others: string[], kindOf: (id: string) => string, truth: Arr) {
  const rowsA: string[] = ['subset\tsize\tmember_fail\twithheld_fail\tfalse_admit\tmasked_markers\t'
    + 'masked_n\tmodel_contam\tmeasured_masked\tmeasured_plain\topp_masked']
  let memberFail = 0, withheldFail = 0, members = 0, withheldN = 0, maskWorse = 0, built = 0
  const falseAdmit = new Map<string, number>()
  const verdicts = new Map<string, Map<string, number>>()
  const t0 = Date.now()
  for (const S of subsets(paternal, inferred.MIN_SEED)) {
    const withheld = paternal.filter((p) => !S.includes(p))
    const nonS = ps.ids.filter((id) => !S.includes(id))
    const mMin = ps.chooseM(nonS).mMin
    const judged = inferred.judgeProducts(ps, S, [...withheld, ...others], mMin, freqB, mono)
    const mf = judged.filter((j) => S.includes(j.id) && j.verdict !== 'product')
    const wf = judged.filter((j) => withheld.includes(j.id) && j.verdict !== 'product')
    const fa = judged.filter((j) => others.includes(j.id) && j.verdict === 'product')
    for (const j of judged) {
      if (!others.includes(j.id)) continue
      const k = kindOf(j.id)
      if (!verdicts.has(k)) verdicts.set(k, new Map())
      verdicts.get(k)!.set(j.verdict, (verdicts.get(k)!.get(j.verdict) ?? 0) + 1)
    }
    members += S.length; withheldN += withheld.length
    memberFail += mf.length; withheldFail += wf.length
    for (const j of fa) falseAdmit.set(j.id, (falseAdmit.get(j.id) ?? 0) + 1)
    let tail = '\t\t\t\t\t'
    if (S.length >= inferred.MIN_PRODUCTS) {
      const masked = ps.build(mMin, nonS, undefined, freqB)
      const plain = ps.build(mMin, nonS)
      const m = measured(masked.genotype, truth)
      const p = measured(plain.genotype, truth)
      built += 1
      if (!(m.contamination < p.contamination)) maskWorse += 1
      tail = `\t${masked.masked}\t${masked.markers}\t${masked.contamination.toFixed(5)}\t`
        + `${m.contamination.toFixed(5)}\t${p.contamination.toFixed(5)}\t${m.opposite.toFixed(5)}`
    }
    rowsA.push(`${S.join(',')}\t${S.length}\t${mf.map((j) => `${j.id}:${j.verdict}:${j.call.r.toFixed(3)}`).join(',')}\t`
      + `${wf.map((j) => `${j.id}:${j.verdict}:${j.call.r.toFixed(3)}`).join(',')}\t`
      + `${fa.map((j) => `${j.id}:${j.call.r.toFixed(3)}`).join(',')}${tail}`)
  }
  writeFileSync(join(OUT, `withholding-A-${label}.tsv`), rowsA.join('\n') + '\n')
  say(`A ${label}: ${rowsA.length - 1} subsets in ${((Date.now() - t0) / 1000).toFixed(0)}s`)
  say(`  A1 members not 'product': ${memberFail} of ${members}`)
  say(`  A2 withheld paternal not 'product': ${withheldFail} of ${withheldN}`)
  say(`  A3 non-paternal arrays read 'product': ${falseAdmit.size === 0 ? 'none'
    : [...falseAdmit].map(([id, n]) => `${id} (${kindOf(id)}) in ${n} subsets`).join('; ')}`)
  for (const [k, m] of verdicts) say(`     ${k}: ${[...m].map(([v, n]) => `${v} ${n}`).join(', ')}`)
  say(`  A4 masked rebuild NOT below unmasked in measured contamination: ${maskWorse} of ${built}`)
}

// --- GSE148488 ----------------------------------------------------------------------------------
const PN = ['GSM4774680', 'GSM4774681', 'GSM4774682', 'GSM4774683',
  'GSM4774684', 'GSM4774685', 'GSM4774686', 'GSM4774687']
const MN = ['GSM4774673', 'GSM4774674', 'GSM4774675', 'GSM4774676',
  'GSM4774677', 'GSM4774678', 'GSM4774679']
const EGG = ['GSM4472407', 'GSM4472415']
const truthWord = new Map<string, Map<string, string>>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , , chrom, manual] = l.split('\t')
  if (!gsm) continue
  if (!truthWord.has(gsm)) truthWord.set(gsm, new Map())
  truthWord.get(gsm)!.set(chrom.replace(/^chr/i, ''), (manual ?? '').trim().toLowerCase())
}
const S3 = [...truthWord.keys()].filter((g) => readdirSync(TRIOS).some((f) => f.startsWith(g)))

const ps = new inferred.ProductSet()
const arrays = new Map<string, Arr>()
const t0 = Date.now()
for (const g of [...PN, ...MN, ...EGG, ...S3]) arrays.set(g, load(fileIn(TRIOS, g), g, ps))
const father = load(fileIn(TRIOS, 'GSM4472397'), 'GSM4472397')
say(`loaded ${arrays.size} GSE148488 arrays and the father's own in ${((Date.now() - t0) / 1000).toFixed(0)}s, `
  + `${idx.size} autosomal probes`)
const kindG = (id: string) => MN.includes(id) ? 'maternal pronucleus' : EGG.includes(id) ? 'egg donor'
  : S3.includes(id) ? 'Table S3 cell' : PN.includes(id) ? 'paternal pronucleus' : 'other'
partA('GSE148488', ps, PN, [...MN, ...EGG, ...S3], kindG, father)

// --- B: the one-parent copy state against withheld rebuilds -------------------------------------
if (!process.env.OM_SKIP_B) {
  let seed = 20261111
  const rand = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296)
  const pick = (size: number, n: number): string[][] => {
    const all = [...subsets(PN, size)].filter((s) => s.length === size)
    const out: string[][] = []
    while (out.length < n && all.length) out.push(all.splice(Math.floor(rand() * all.length), 1)[0])
    return out
  }
  type Variant = { name: string; products: string[]; masked: boolean; gt: Uint8Array }
  const variants: Variant[] = [{ name: 'measured', products: [], masked: false, gt: father.gt }]
  const add = (products: string[], masked: boolean) => {
    const nonS = ps.ids.filter((id) => !products.includes(id))
    const r = ps.build(ps.chooseM(nonS).mMin, nonS, undefined, masked ? freqB : undefined)
    const gt = new Uint8Array(CAP)
    for (const [id, g] of r.genotype) gt[idx.get(id)!] = g === 'AA' ? AA : BB
    variants.push({ name: `${masked ? 'masked' : 'plain'}-${products.length}:${products.map((p) => p.slice(-2)).join('')}`,
      products, masked, gt })
  }
  add(PN, true); add(PN, false)
  for (const p of PN) add(PN.filter((x) => x !== p), true)
  for (const s of pick(6, 6)) add(s, true)
  for (const s of pick(5, 6)) add(s, true)
  const OLD5 = ['GSM4774680', 'GSM4774681', 'GSM4774682', 'GSM4774683', 'GSM4774685']
  add(OLD5, true); add(OLD5, false)
  say(`B: ${variants.length} parents (the measured father and ${variants.length - 1} rebuilds)`)

  const freqArr = new Float64Array(CAP).fill(NaN)
  const panelArr = new Uint8Array(CAP)
  for (const [id, i] of idx) {
    const q = freqB.get(id); if (q !== undefined) freqArr[i] = q
    const g = panel.get(id); if (g) panelArr[i] = g === 'AA' ? AA : BB
  }
  const cells = [...S3, ...PN, ...MN]
  // calls[cell][variant] = state per autosome 1..22 (undefined where no tally)
  const calls = new Map<string, (ReturnType<typeof ops.shareCall> | undefined)[][]>()
  const V = variants.length
  for (const cell of cells) {
    const a = arrays.get(cell)!
    const T = variants.map(() => Array.from({ length: 23 }, () => ops.emptyShareTally()))
    const touched = variants.map(() => new Uint8Array(23))
    for (let i = 0; i < idx.size; i += 1) {
      const b = a.baf[i]
      if (Number.isNaN(b)) continue
      const q = freqArr[i]
      if (Number.isNaN(q)) continue
      const c = chromNum[i]
      const pg = panelArr[i]
      for (let v = 0; v < V; v += 1) {
        const g = variants[v].gt[i]
        if (g !== AA && g !== BB) continue
        ops.addShare(T[v][c], g === AA ? 'AA' : 'BB', b, q, pg === g)
        touched[v][c] = 1
      }
    }
    calls.set(cell, T.map((t, v) => t.map((x, c) => (c >= 1 && touched[v][c] ? ops.shareCall(x) : undefined))))
  }

  const lostTag = (w: string): 'pat' | 'mat' | null =>
    /^monosomy \d+ \((pat|mat)\)$/.exec(w)?.[1] as 'pat' | 'mat' | undefined ?? null
  const out: string[] = ['variant\tagree\tinformative\twrong_parent\tlost_on_normal_measured_disomic\t'
    + 'copy_decision_same\tpn_other_lost\tmn_parent_lost\tmonosomy_confirmed']
  for (let v = 0; v < V; v += 1) {
    let agree = 0, inf = 0, onNormal = 0, decisionSame = 0, confirmed = 0
    const wrong: string[] = []
    for (const cell of S3) {
      const cv = calls.get(cell)![v], cm = calls.get(cell)![0]
      for (let c = 1; c <= 22; c += 1) {
        const x = cv[c], m = cm[c]
        const w = truthWord.get(cell)?.get(String(c)) ?? ''
        if (x && m && x.state !== 'undecided' && m.state !== 'undecided') {
          inf += 1
          if (x.state === m.state) agree += 1
        }
        if (x && w === 'normal' && (x.state === 'other-lost' || x.state === 'parent-lost') && m?.state === 'disomic') onNormal += 1
        const tag = lostTag(w)
        if (x && tag) {
          const said = x.state === 'parent-lost' ? 'pat' : x.state === 'other-lost' ? 'mat' : null
          if (said === tag) confirmed += 1
          else if (said) wrong.push(`${cell}:${c}`)
        }
      }
      const frac = (k: number) => ops.normalFraction(calls.get(cell)![k].filter((x): x is NonNullable<typeof x> => !!x))
      const dv = frac(v), dm = frac(0)
      const decide = (d: number) => (Number.isFinite(d) ? d >= ops.NORMAL_MAJORITY_MIN : 'none')
      if (decide(dv) === decide(dm)) decisionSame += 1
    }
    const share = (list: string[], want: string) => {
      let hit = 0, n = 0
      for (const cell of list) {
        if (variants[v].products.includes(cell)) continue
        for (const x of calls.get(cell)![v]) {
          if (!x || x.state === 'undecided') continue
          n += 1; if (x.state === want) hit += 1
        }
      }
      return n ? hit / n : NaN
    }
    out.push(`${variants[v].name}\t${agree}\t${inf}\t${wrong.join(',')}\t${onNormal}\t${decisionSame}\t`
      + `${share(PN, 'other-lost').toFixed(4)}\t${share(MN, 'parent-lost').toFixed(4)}\t${confirmed}`)
  }
  writeFileSync(join(OUT, 'withholding-B.tsv'), out.join('\n') + '\n')
  for (const l of out) say(`  ${l}`)
}

// --- S23 ----------------------------------------------------------------------------------------
// The laboratory's series number stays out of the repository: see audit/local/barcode-map.tsv.
const SERIES = process.env.OM_S23_SERIES
if (S23 && existsSync(S23) && SERIES) {
  const P = ['02', '04', '05', '09', '11', '13']
  const M = ['01', '03', '06', '07', '08', '10', '14', '18']
  const U = ['12', '15', '16', '17']
  const name = (n: string) => `${SERIES}-${n}`
  const ps23 = new inferred.ProductSet()
  for (const n of [...P, ...M, ...U]) load(fileIn(S23, name(n)), name(n), ps23)
  const truth = load(fileIn(S23, name('00')), name('00'))
  const kindS = (id: string) => M.map(name).includes(id) ? 'sheet maternal' : 'unsettled'
  partA('S23', ps23, P.map(name), [...M, ...U].map(name), kindS, truth)
}
writeFileSync(join(OUT, 'withholding-AB-summary.txt'), lines.join('\n') + '\n')

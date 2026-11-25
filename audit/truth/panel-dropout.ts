/**
 * Allele drop-out and call rate by material, on the lab chip's public trios: the first measurement
 * for the marker panel, set aside because the estimator is pulled down by drop-in and the mean over
 * arrays is pulled up by arrays missing a parent's copy. The panel reads audit/truth/panel-reads.ts.
 *
 * d is the tool's own estimator (ingest.estimateDropout): 1 - 2h at markers where the father is
 * heterozygous and the mother homozygous, so it needs no population frequency. The call rate is
 * the share of chromosome 1 to 22 markers called. A marker reads with chance callRate x (1 - d).
 * A parent over the diploid ceiling is left out by the rule the tool applies
 * (stage.parentOverCeiling), as in audit/trio-calibration.ts.
 *
 * Run: OM_TRIOS=<dir of GSE148488 probes.gz + trio_manifest_full.csv> \
 *        node --experimental-strip-types audit/truth/panel-dropout.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const stage = await import(`${W}stage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('panel-dropout: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Arr { gt: Map<string, AB>; het: number; callRate: number }

const cache = new Map<string, Arr>()
function load(gsm: string, keep: boolean): Arr {
  const hit = cache.get(gsm)
  if (hit) return hit
  const lines = gunzipSync(readFileSync(`${DIR}/${gsm}.probes.gz`)).toString('utf8').split('\n')
  let h = 0
  while (lines[h].startsWith('#')) h += 1
  const map = ingest.headerMap(lines[h])
  const gt = new Map<string, AB>()
  let het = 0; let called = 0; let total = 0
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || !/^\d+$/.test(r.chrom)) continue
    total += 1
    if (r.genotype === 'NC') continue
    gt.set(r.probesetId, r.genotype as AB)
    called += 1
    if (r.genotype === 'AB') het += 1
  }
  const out = { gt, het: called ? het / called : NaN, callRate: total ? called / total : NaN }
  if (keep) cache.set(gsm, out)
  return out
}

const rows = readFileSync(`${DIR}/trio_manifest_full.csv`, 'utf8').trim().split('\n')
const head = rows[0].split(',')
const recs = rows.slice(1).map((l) => {
  const f = l.split(',')
  return Object.fromEntries(head.map((k, i) => [k, f[i] ?? ''])) as Record<string, string>
})
const present = (g: string) => existsSync(`${DIR}/${g}.probes.gz`)
const parents = recs.filter((r) => r.role === 'parent' && present(r.gsm)).map((r) => r.gsm)
const hetOf = new Map(parents.map((g) => [g, load(g, true).het]))
const panelHet = Math.min(...[...hetOf.values()].filter((x) => Number.isFinite(x)))
const over = new Set(parents.filter((g) => stage.parentOverCeiling(hetOf.get(g) as number, panelHet)))
console.log(`parents read ${parents.length}, over the diploid ceiling and left out ${over.size}`)
const first = (csv: string) =>
  csv.split(/[;| ]/).map((x) => x.trim()).filter((x) => x && !over.has(x) && present(x))[0]

const LABEL: Record<string, string> = {
  'esc-line': 'bulk DNA (ES cell line)', trophectoderm: 'trophectoderm biopsy', blastomere: 'single cell (blastomere)',
}
const by = new Map<string, { d: number; cr: number; p: number; n: number }[]>()
for (const r of recs) {
  if (r.complete_trio.toLowerCase() !== 'true' || !present(r.gsm) || !LABEL[r.material]) continue
  const fa = first(r.father_gsms); const mo = first(r.mother_gsms)
  if (!fa || !mo) continue
  const child = load(r.gsm, false)
  const est = ingest.estimateDropout(load(fa, true).gt, load(mo, true).gt, child.gt)
  if (!est) continue
  const list = by.get(r.material) ?? []
  list.push({ d: est.d, cr: child.callRate, p: child.callRate * (1 - est.d), n: est.n })
  by.set(r.material, list)
}

const mean = (x: number[]) => x.reduce((s, v) => s + v, 0) / x.length
const median = (x: number[]) => { const s = [...x].sort((a, b) => a - b); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2 }
const f4 = (v: number) => v.toFixed(4)
console.log('material\tarrays\tmean d\tmedian d\tmin d\tmax d\tmean call rate\tmean read chance\tmin read chance\tmax read chance\tmedian markers in d')
for (const m of Object.keys(LABEL)) {
  const l = by.get(m) ?? []
  if (!l.length) { console.log(`${LABEL[m]}\t0`); continue }
  const d = l.map((x) => x.d); const p = l.map((x) => x.p)
  console.log([LABEL[m], l.length, f4(mean(d)), f4(median(d)), f4(Math.min(...d)), f4(Math.max(...d)),
    f4(mean(l.map((x) => x.cr))), f4(mean(p)), f4(Math.min(...p)), f4(Math.max(...p)),
    median(l.map((x) => x.n))].join('\t'))
}

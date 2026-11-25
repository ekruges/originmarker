/**
 * The chance one marker reads, by material, on the lab chip's public trios. The marker panel's
 * count of markers per side is read from these numbers (originmarker/panelbuilder.py, MATERIALS).
 *
 * At a site where the two parents are opposite homozygotes, a child that carries both parents'
 * copies must be heterozygous. The read chance is the share of those sites called heterozygous; it
 * fails by a no-call or by a homozygous call. Counted only on arrays the tool reads as carrying both
 * parents, and on those only on chromosomes with no change on the result card, both from the
 * shipped scoreSample with both parents loaded (audit/truth/run.ts records). Beside it, for each
 * array, the drop-out estimator ingest.estimateDropout (1 - 2h at father-het by mother-hom sites).
 * A parent over the diploid ceiling is left out by stage.parentOverCeiling.
 *
 * Arrays are named by material and index only.
 *
 *   OM_TRIOS=<dir> node --experimental-strip-types audit/truth/panel-reads.ts list > two.tsv
 *   OM_TRIOS=<dir> node --experimental-strip-types audit/truth/panel-reads.ts measure <records.jsonl>
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const stage = await import(`${W}stage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('panel-reads: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}
const mode = process.argv[2]

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Arr { gt: Map<string, AB>; chrom: Map<string, string>; het: number; callRate: number }
const cache = new Map<string, Arr>()
function load(gsm: string, keep: boolean): Arr {
  const hit = cache.get(gsm)
  if (hit) return hit
  const lines = gunzipSync(readFileSync(`${DIR}/${gsm}.probes.gz`)).toString('utf8').split('\n')
  let h = 0
  while (lines[h].startsWith('#')) h += 1
  const map = ingest.headerMap(lines[h])
  const gt = new Map<string, AB>(); const chrom = new Map<string, string>()
  let het = 0; let called = 0; let total = 0
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || !/^\d+$/.test(r.chrom)) continue
    total += 1
    if (!keep) chrom.set(r.probesetId, r.chrom)   // only a child's sites are placed
    gt.set(r.probesetId, r.genotype as AB)
    if (r.genotype === 'NC') continue
    called += 1
    if (r.genotype === 'AB') het += 1
  }
  const out = { gt, chrom, het: called ? het / called : NaN, callRate: total ? called / total : NaN }
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
const first = (csv: string) =>
  csv.split(/[;| ]/).map((x) => x.trim()).filter((x) => x && !over.has(x) && present(x))[0]

const LABEL: Record<string, string> = {
  'esc-line': 'bulk DNA (ES cell line)', trophectoderm: 'trophectoderm biopsy', blastomere: 'single cell (blastomere)',
}
const trios = recs
  .filter((r) => r.complete_trio.toLowerCase() === 'true' && present(r.gsm) && LABEL[r.material])
  .map((r) => ({ gsm: r.gsm, material: r.material, father: first(r.father_gsms), mother: first(r.mother_gsms) }))
  .filter((t) => t.father && t.mother)

if (mode === 'list') {
  for (const t of trios) {
    console.log([t.gsm, `${DIR}/${t.gsm}.probes.gz`, `${DIR}/${t.father}.probes.gz`, 'two',
      `${DIR}/${t.mother}.probes.gz`].join('\t'))
  }
  process.exit(0)
}

// The tool's own reading of each child, both parents loaded.
const verdicts = new Map<string, { both: boolean; changed: Set<string>; why: string }>()
for (const l of readFileSync(process.argv[3], 'utf8').split('\n')) {
  if (!l.startsWith('{')) continue   // run.ts ends each file with a DONE line
  const d = JSON.parse(l)
  if (d.tag !== 'two') continue
  const pv = d.parentGenome?.verdict; const mv = d.parentGenomeOther?.verdict
  const changed = new Set<string>((d.card?.counted ?? []).map((c: { chrom: string }) => String(c.chrom)))
  verdicts.set(d.id, { both: pv === 'carries' && mv === 'carries', changed, why: `father ${pv}, mother ${mv}` })
}

const mean = (x: number[]) => x.reduce((s, v) => s + v, 0) / x.length
const median = (x: number[]) => { const s = [...x].sort((a, b) => a - b); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2 }
const f4 = (v: number) => (Number.isFinite(v) ? v.toFixed(4) : 'n/a')

console.log(`parents read ${parents.length}, over the diploid ceiling and left out ${over.size}`)
console.log('per array: material, index, estimator d, call rate, both parents by the tool, chromosomes with a change, '
  + 'opposite-homozygote sites counted, read, not called, called homozygous (of those called), counted')
const by = new Map<string, { read: number; lost: number; nc: number }[]>()
const idx = new Map<string, number>()
for (const t of trios) {
  const i = (idx.get(t.material) ?? 0) + 1
  idx.set(t.material, i)
  const fa = load(t.father, true); const mo = load(t.mother, true); const ch = load(t.gsm, false)
  // No-calls are in the maps; the estimator skips them, as it skips a missing call.
  const est = ingest.estimateDropout(fa.gt, mo.gt, ch.gt)
  const v = verdicts.get(t.gsm)
  let n = 0; let het = 0; let nc = 0; let hom = 0
  if (v?.both) {
    for (const [id, f] of fa.gt) {
      const m = mo.gt.get(id)
      if (!((f === 'AA' && m === 'BB') || (f === 'BB' && m === 'AA'))) continue
      if (v.changed.has(ch.chrom.get(id) ?? '')) continue
      const c = ch.gt.get(id)
      if (c === undefined) continue
      n += 1
      if (c === 'AB') het += 1
      else if (c === 'NC') nc += 1
      else hom += 1
    }
  }
  const counted = !!v?.both && n > 0
  if (counted) {
    const l = by.get(t.material) ?? []
    l.push({ read: het / n, nc: nc / n, lost: hom / (n - nc) })
    by.set(t.material, l)
  }
  console.log([LABEL[t.material], i, f4(est?.d ?? NaN), f4(ch.callRate), v ? (v.both ? 'yes' : `no (${v.why})`) : 'not scored',
    v ? ([...v.changed].sort((a, b) => +a - +b).join(' ') || 'none') : 'n/a', n, f4(n ? het / n : NaN), f4(n ? nc / n : NaN),
    f4(n - nc ? hom / (n - nc) : NaN), counted ? 'yes' : 'no'].join('\t'))
}

console.log('\nmaterial\tarrays counted\tmean read\tmedian read\tmin read\tmax read\tmean not called\tmean called homozygous\tmedian called homozygous')
for (const m of Object.keys(LABEL)) {
  const l = by.get(m) ?? []
  if (!l.length) { console.log(`${LABEL[m]}\t0`); continue }
  const r = l.map((x) => x.read)
  console.log([LABEL[m], l.length, f4(mean(r)), f4(median(r)), f4(Math.min(...r)), f4(Math.max(...r)),
    f4(mean(l.map((x) => x.nc))), f4(mean(l.map((x) => x.lost))), f4(median(l.map((x) => x.lost)))].join('\t'))
}

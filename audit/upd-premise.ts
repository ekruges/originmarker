/**
 * Whether the heterozygosity-loss detectors hold on the data alone, with the ploidy call withheld.
 *
 * `detectLoh` and `detectUpd` exclude a genome the tool has called uniparental. On a genome whose
 * ploidy was refused, `unknown`, that exclusion does not apply, so what keeps them honest there is
 * the array's own background heterozygosity. This drives both detectors with zygosity `unknown` on
 * genomes whose answer is known from the source:
 *   one complement   GSE12713 complete moles, stated monospermic by the series; GSE148488 pronuclei,
 *                    separated by micromanipulation; GSE19247 single sperm. Homozygous by biology,
 *                    so every finding is false.
 *   diploid, bulk    GSE148488 adult gamete donors; GSE60909 parents and grandparents; GSE18932 pure
 *                    fibroblast lines. A whole-chromosome finding is false. Segmental runs are
 *                    counted and not scored, because long homozygous stretches occur in healthy
 *                    adults.
 *   isodisomy        GSE18932 GM11496 (chr7) and GM15603 (chr8), karyotype stated per sample.
 *
 * Run: any of OM_MOLES_GT OM_TRIOS OM_GSE OM_GSE_B OM_PGD OM_CELLS set to their directories, then
 *        node --experimental-strip-types --max-old-space-size=6144 audit/upd-premise.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const ab = await import(`${W}abnormalities.ts`)
const scan = await import(`${W}scan.ts`)

type Kind = 'one' | 'bulk' | 'upd'
type Item = { path: string; set: string; kind: Kind; id: string; chrom?: string }
const items: Item[] = []
const json = (p: string) => JSON.parse(readFileSync(p, 'utf8'))
const dir = (k: string) => { const d = process.env[k]; return d && existsSync(d) ? d : null }
const push = (d: string, file: string, set: string, kind: Kind, id: string, chrom?: string) => {
  if (existsSync(join(d, file))) items.push({ path: join(d, file), set, kind, id, chrom })
}

let d = dir('OM_MOLES_GT')
if (d) for (const t of json(join(d, 'truth.json'))) push(d, `${t.gsm}.probes`, 'GSE12713 monospermic moles', 'one', t.gsm)
d = dir('OM_TRIOS')
if (d) {
  const [head, ...rows] = readFileSync(join(d, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
  const h = head.split(',')
  const gi = h.indexOf('gsm'); const mi = h.indexOf('material'); const ri = h.indexOf('role')
  for (const l of rows) {
    const f = l.split(',')
    if (f[mi] === 'pronucleus') push(d, `${f[gi]}.probes.gz`, 'GSE148488 pronuclei', 'one', f[gi])
    else if (f[mi] === 'bulk-gdna' && f[ri] === 'parent') push(d, `${f[gi]}.probes.gz`, 'GSE148488 adult donors', 'bulk', f[gi])
  }
}
for (const k of ['OM_GSE', 'OM_GSE_B']) {
  d = dir(k)
  if (d) for (const t of json(join(d, 'truth.json'))) {
    if (/^sperm/i.test(t.src ?? '')) push(d, `${t.gsm}.probes`, 'GSE19247 single sperm', 'one', t.gsm)
  }
}
d = dir('OM_PGD')
if (d) for (const t of json(join(d, 'truth.json'))) {
  if (/^(Mother|Father|Grandmother|Grandfather)_PGD/i.test(t.title)) push(d, `${t.gsm}.probes`, 'GSE60909 parents and grandparents', 'bulk', t.gsm)
}
d = dir('OM_CELLS')
if (d) for (const t of json(join(d, 'truth.json'))) {
  const c = /UPID\s*chr(\w+)/i.exec(t.karyotype ?? '')?.[1]
  if (c) push(d, `${t.gsm}.probes`, 'GSE18932 isodisomy lines', 'upd', t.gsm, c)
  else if (/mixture (100%GM\d+_0%|0%GM\d+_100%)/.test(t.title)) push(d, `${t.gsm}.probes`, 'GSE18932 pure lines', 'bulk', t.gsm)
}

function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}
const isAutosome = (c: string) => /^\d+$/.test(c) && +c >= 1 && +c <= 22
const q = (xs: number[], p: number) => {
  const s = [...xs].sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}

interface Tally { n: number; threw: number; bg: number[]; wholeArrays: number; segArrays: number; whole: number; seg: number
  stated: number; examples: string[] }
const sets = new Map<string, Tally & { kind: Kind }>()
for (const it of items) {
  const t = sets.get(it.set) ?? { kind: it.kind, n: 0, threw: 0, bg: [], wholeArrays: 0, segArrays: 0, whole: 0, seg: 0, stated: 0, examples: [] }
  sets.set(it.set, t)
  try {
    const cnByChrom = new Map<string, { chrom: string; pos: number; called: boolean; log2R: number | null }[]>()
    const self: { chrom: string; pos: number; het: boolean }[] = []
    const lrr: number[] = []
    for (const r of rowsOf(it.path)()) {
      if (!isAutosome(r.chrom)) continue
      const cn = cnByChrom.get(r.chrom) ?? []
      cn.push({ chrom: r.chrom, pos: r.pos, called: r.genotype !== 'NC', log2R: r.log2R })
      cnByChrom.set(r.chrom, cn)
      if (r.log2R !== null && Number.isFinite(r.log2R)) lrr.push(r.log2R)
      if (r.genotype !== 'NC') self.push({ chrom: r.chrom, pos: r.pos, het: r.genotype === 'AB' })
    }
    lrr.sort((a, b) => a - b)
    const genomeLrr = lrr.length ? lrr[lrr.length >> 1] : 0
    const background = self.length ? self.filter((m) => m.het).length / self.length : NaN
    const { windows, chromEnd } = await scan.copyNeutralWindows(cnByChrom, self, genomeLrr, ab.LOH_SEGMENT_MARKERS)
    const found = [
      ...ab.mergeLoh(ab.detectLoh(windows, { zygosity: 'unknown' })),
      ...ab.detectUpd(ab.runsOfHomozygosity(self, { chromEndBp: chromEnd }),
        { zygosity: 'unknown', backgroundHet: background, intensity: lrr.length > 0 }),
    ] as { cls: string; chrom: string; wholeChromosome: boolean }[]
    t.n += 1
    t.bg.push(background)
    const whole = found.filter((f) => f.wholeChromosome && f.chrom !== it.chrom)
    const seg = found.filter((f) => !f.wholeChromosome)
    if (it.chrom && found.some((f) => f.wholeChromosome && f.chrom === it.chrom)) t.stated += 1
    t.whole += whole.length
    t.seg += seg.length
    if (whole.length) t.wholeArrays += 1
    if (seg.length) t.segArrays += 1
    if ((whole.length || (it.kind === 'one' && seg.length)) && t.examples.length < 4) {
      t.examples.push(`${it.id} bg ${background.toFixed(4)}: ${[...whole, ...seg].slice(0, 5).map((f) => `${f.cls}:${f.chrom}`).join(' ')}`)
    }
  } catch (e) {
    t.threw += 1
    if (t.examples.length < 4) t.examples.push(`THREW ${it.id}: ${String((e as Error).message ?? e).slice(0, 80)}`)
  }
}

console.log('detectLoh and detectUpd driven with zygosity `unknown`, on genomes with a known answer')
console.log('')
for (const [name, t] of sets) {
  console.log(`${name}  (${t.kind === 'one' ? 'one complement: every finding false' : t.kind === 'bulk' ? 'diploid bulk: whole-chromosome false, segmental counted' : 'stated isodisomy'})`)
  console.log(`  arrays ${t.n}, threw ${t.threw}; background heterozygosity min ${q(t.bg, 0).toFixed(4)} `
    + `p50 ${q(t.bg, 0.5).toFixed(4)} max ${q(t.bg, 1).toFixed(4)}`)
  if (t.kind === 'upd') console.log(`  stated chromosome flagged whole-chromosome: ${t.stated}/${t.n}`)
  console.log(`  whole-chromosome findings${t.kind === 'upd' ? ' on other chromosomes' : ''}: ${t.whole} on ${t.wholeArrays} arrays`)
  console.log(`  segmental findings: ${t.seg} on ${t.segArrays} arrays`)
  for (const e of t.examples) console.log(`    ${e}`)
}

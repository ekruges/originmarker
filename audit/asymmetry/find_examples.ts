/**
 * Find bundle-ready examples in a public series: a parent, and the children that best show what
 * the tool does.
 *
 * Everything is judged at the SUBSET density the examples ship at, because an example that only
 * works on the full array demonstrates nothing once bundled. A confirmed segmental loss was
 * generated at three densities during an earlier attempt and found at none of them; whole
 * chromosomes survive because they use every marker on the chromosome.
 *
 *   usage: node --experimental-strip-types audit/asymmetry/find_examples.ts <dir> [stride]
 */
import { readFileSync, readdirSync, writeFileSync } from 'node:fs'
import { join, basename } from 'node:path'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const { isAutosome, emptyTally, tallyRow, classify } = await import(`${W}parentage.ts`)
const { inferStage } = await import(`${W}stage.ts`)
const ob = await import(`${W}obligateHet.ts`)
const dosage = await import(`${W}dosageOrigin.ts`)

const DIR = process.argv[2]
const STRIDE = Number(process.argv[3] ?? 8)
if (!DIR) throw new Error('usage: find_examples.ts <dir> [stride]')

type AB = 'AA' | 'AB' | 'BB' | 'NC'
const rd = (p: string) => {
  const b = readFileSync(p)
  return (p.endsWith('.gz') ? gunzipSync(b) : b).toString('utf8').split('\n')
}
const hdr = (L: string[]) => {
  let h = -1
  for (let i = 0; i < 60; i += 1) if (L[i] && !L[i].startsWith('#')) { h = i; break }
  return { h, map: h < 0 ? null : ingest.headerMap(L[h]) }
}

interface Arr {
  file: string
  id: string
  gt: Map<string, AB>
  callRate: number
  hetRate: number
  stage: string
  bafSd: number
  aneu: { chrom: string, kind: string }[]
}

/** One pass per array at the subset density: genotypes, stage, quality, and its own aneuploidy. */
function read(file: string): Arr | null {
  const L = rd(join(DIR, file))
  const { h, map } = hdr(L)
  if (!map) return null
  const gt = new Map<string, AB>()
  const byChrom = new Map()
  const sums = ingest.emptyBafSums()
  const t = emptyTally()
  const het: number[] = []
  let first = ''
  for (let i = h + 1; i < L.length; i += STRIDE) {
    const r = ingest.parseRow(L[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, sums as never)
    // Self-tally, so a chromosome's call-rate collapse is visible without any reference.
    tallyRow((r.genotype ?? 'NC') as never, r as never, t as never)
    if (!isAutosome(r.chrom)) continue
    if (r.genotype !== 'NC') gt.set(r.probesetId, r.genotype as AB)
    if (r.genotype === 'AB' && r.baf !== null && Number.isFinite(r.baf)) het.push(r.baf)
  }
  const p = ingest.finishProfile('x', byChrom as never, sums as never, first)
  const st = inferStage(p)
  const mu = het.length ? het.reduce((a, b) => a + b, 0) / het.length : NaN
  const bafSd = het.length > 1
    ? Math.sqrt(het.reduce((a, x) => a + (x - mu) ** 2, 0) / (het.length - 1)) : NaN
  const cls = classify(t as never, 0.1626, { role: 'paternal' }) as {
    chroms: { chrom: string, aneuploidy?: string }[]
  }
  return {
    file,
    id: basename(file).replace(/\.(probes|CEL\.txt\.gz|txt\.gz)$/, '').replace(/^(GSM\d+)_.*$/, '$1'),
    gt,
    callRate: p.callRate,
    hetRate: p.hetRate,
    stage: st.stage,
    bafSd,
    aneu: cls.chroms.filter((c) => c.aneuploidy).map((c) => ({ chrom: c.chrom, kind: c.aneuploidy! })),
  }
}

const files = readdirSync(DIR).filter((f) => /\.(probes|probes\.gz|txt\.gz)$/.test(f)).sort()
process.stderr.write(`${files.length} arrays, stride ${STRIDE}\n`)
const arrays: Arr[] = []
for (const f of files) {
  const a = read(f)
  if (a && a.stage !== 'failed') arrays.push(a)
}
process.stderr.write(`${arrays.length} usable\n\n`)

// A PARENT IS A BULK ARRAY. Judged on call rate rather than the stage label, which is derived from
// heterozygosity against a population anchor and has mislabelled real bulk gDNA as single-cell.
// Every array carrying an event, with the quality that decides whether dosage can speak for it.
process.stderr.write('arrays with a whole-chromosome event:\n')
for (const a of arrays.filter((x) => x.aneu.length)) {
  process.stderr.write(`  ${a.id.padEnd(12)} ${a.stage.padEnd(14)} bafSd ${a.bafSd.toFixed(3)} `
    + `call ${a.callRate.toFixed(3)}  ${a.aneu.map((e) => `chr${e.chrom}:${e.kind}`).join(' ')}\n`)
}
process.stderr.write('\n')

const parents = arrays
  .filter((a) => a.callRate >= 0.93 && a.hetRate >= 0.12 && a.hetRate <= 0.22 && a.bafSd < 0.13)
  .sort((a, b) => a.bafSd - b.bafSd)
process.stderr.write(`bulk-like candidates (${parents.length}):\n`)
for (const p of parents.slice(0, 12)) {
  process.stderr.write(`  ${p.id.padEnd(12)} call ${p.callRate.toFixed(3)} het ${p.hetRate.toFixed(3)} `
    + `bafSd ${p.bafSd.toFixed(3)} ${p.stage} ${p.aneu.length ? `[${p.aneu.length} aneuploid]` : ''}\n`)
}

const opp = (a: Map<string, AB>, b: Map<string, AB>) => {
  let n = 0
  let o = 0
  for (const [probe, ga] of a) {
    if (ga !== 'AA' && ga !== 'BB') continue
    const gb = b.get(probe)
    if (gb !== 'AA' && gb !== 'BB') continue
    n += 1
    if (ga !== gb) o += 1
  }
  return { rate: n ? o / n : NaN, n }
}

const found: Record<string, unknown>[] = []
for (const p of parents.slice(0, 6)) {
  const kids: { a: Arr, opp: number, frac: number }[] = []
  for (const c of arrays) {
    if (c.id === p.id) continue
    const o = opp(p.gt, c.gt)
    if (!(o.rate <= 0.02) || o.n < 5_000) continue
    const h = ob.emptyHet()
    for (const [probe, pg] of p.gt) {
      const cg = c.gt.get(probe)
      if (cg) ob.addOneParent(pg as never, cg as never, h as never)
    }
    const link = ob.hetCall(h as never, 1) as { ploidy: string, fraction: number }
    if (link.ploidy !== 'biparental') continue
    kids.push({ a: c, opp: o.rate, frac: link.fraction })
  }
  if (!kids.length) continue
  process.stderr.write(`\n=== ${p.id}: ${kids.length} confirmed child(ren)\n`)
  for (const k of kids) {
    const ev = k.a.aneu.map((x) => `chr${x.chrom}:${x.kind}`).join(' ')
    process.stderr.write(`  ${k.a.id.padEnd(12)} opp ${k.opp.toFixed(4)} second ${(100 * k.frac).toFixed(1)}% `
      + `${k.a.stage.padEnd(14)} bafSd ${k.a.bafSd.toFixed(3)} ${ev || 'euploid'}\n`)
    found.push({
      parent: p.id, child: k.a.id, opp: k.opp, second: k.frac, stage: k.a.stage,
      bafSd: k.a.bafSd, callRate: k.a.callRate, events: k.a.aneu,
      dosageUsable: k.a.bafSd <= dosage.MAX_HET_BAF_SD,
    })
  }
}
writeFileSync('audit/asymmetry/example-candidates.json',
  JSON.stringify({ dir: DIR, stride: STRIDE, parents: parents.map((p) => p.id), found }, null, 1))
process.stderr.write(`\n${found.length} parent-child pairs -> audit/asymmetry/example-candidates.json\n`)

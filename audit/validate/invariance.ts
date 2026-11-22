/**
 * Metamorphic tests: transformations of the input that must not change the answer, so no truth is
 * needed to catch a fault. Every Table S3 cell of GSE148488 is scored by the shipped scoreSample
 * against ONE parent, the laboratory's configuration: the father reconstructed from his pronuclei
 * (OM_FATHER, a file om reconstruct wrote), or his measured array where none is given. Four ways:
 *
 *   orig     as read
 *   shuffle  every file's rows in a different order: nothing may change
 *   flip     A and B swapped everywhere: genotypes AA and BB exchanged and B allele frequency read
 *            as 1 - BAF in the cell and the parent, the founder panel's alleles and frequencies with
 *            them. Allele coding is arbitrary, so nothing may change
 *   half     every other marker, the same ones in every file: claims may be lost, never contradicted
 *
 * One JSON line per cell and variant; audit/validate/scorecard.py compares them.
 *
 * Run: OM_TRIOS=<dir> [OM_FATHER=<file>] node --experimental-strip-types audit/validate/invariance.ts <shard> <shards> <out>
 */
import { readFileSync, appendFileSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const keys = await import(`${W}answerKeys.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)

const TRIOS = process.env.OM_TRIOS!
const [shard, shards, out] = [Number(process.argv[2]), Number(process.argv[3]), process.argv[4]]
const FATHER = 'GSM4472397'

const file = (gsm: string) => join(TRIOS, readdirSync(TRIOS).find((f) => f.startsWith(gsm) && /probes/.test(f))!)
const cells = [...new Set(readFileSync(new URL('../zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')
  .filter((l) => l.startsWith('GSM')).map((l) => l.split('\t')[0]))].sort()
  .filter((g) => readdirSync(TRIOS).some((f) => f.startsWith(g))).filter((_, i) => i % shards === shard)

function rows(path: string): any[] {
  const raw = readFileSync(path)
  const lines = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8').split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = ingest.headerMap(lines[h])
  const outRows: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) outRows.push(r) }
  return outRows
}

const hash = (s: string) => { let x = 2166136261; for (let i = 0; i < s.length; i++) x = Math.imul(x ^ s.charCodeAt(i), 16777619); return x >>> 0 }
const swap = (g: string) => (g === 'AA' ? 'BB' : g === 'BB' ? 'AA' : g)
type Variant = 'orig' | 'shuffle' | 'flip' | 'half'
function transform(rs: any[], v: Variant): any[] {
  if (v === 'shuffle') return [...rs].sort((a, b) => hash(a.probesetId + 'x') - hash(b.probesetId + 'x'))
  if (v === 'half') return rs.filter((r) => hash(r.probesetId) % 2 === 0)
  if (v === 'flip') return rs.map((r) => ({ ...r, genotype: swap(r.genotype), baf: r.baf === null ? null : 1 - r.baf }))
  return rs
}

const panel0 = keys.parsePanel(gunzipSync(readFileSync(
  new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url))).toString('utf8'))
function panelFor(v: Variant) {
  if (v !== 'flip') return panel0
  const p: any = new Map([...panel0].map(([id, g]) => [id, swap(g)]))
  p.calibration = panel0.calibration
  p.freqB = new Map([...panel0.freqB!].map(([id, q]) => [id, 1 - q]))
  return p
}

// A reconstruction states its own heterozygosity and absence floor in its banner, read here as the
// application reads it.
const fatherPath = process.env.OM_FATHER ?? file(FATHER)
const banner = (() => {
  const raw = readFileSync(fatherPath)
  return inferredArray.parseInferredBanner((fatherPath.endsWith('.gz') ? gunzipSync(raw) : raw)
    .toString('utf8', 0, inferredArray.MARK_BYTES)) ?? undefined
})()
function parentIndex(rs: any[]) {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rs) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectParentRow(r, acc)
  }
  return score.finishParent(acc, ingest.finishProfile(FATHER, byChrom, bafSums, first).build.build, banner)
}

const VARIANTS: Variant[] = ['orig', 'shuffle', 'flip', 'half']
const fatherRows = rows(fatherPath)
const fathers = new Map(VARIANTS.map((v) => [v, parentIndex(transform(fatherRows, v))]))
const brief = (v: any): string => (v === undefined || v === null ? '' : typeof v === 'object'
  ? String(v.verdict ?? v.kind ?? v.cls ?? JSON.stringify(v).slice(0, 40)) : String(v))

for (const gsm of cells) {
  const cellRows = rows(file(gsm))
  for (const v of VARIANTS) {
    const pat = fathers.get(v)
    const acc = score.emptyCollected(pat as never, null, panelFor(v) as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); const builds = ingest.emptyBuildSums(); let first = ''
    for (const r of transform(cellRows, v)) {
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectRow(r, pat as never, null, acc as never)
      ingest.accumulateBuild(r as never, builds)
    }
    const profile = ingest.finishProfile(gsm, byChrom, bafSums, first, builds)
    const res: any = await score.scoreSample({ acc, profile, pat, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: gsm, log: () => {} } as never)
    appendFileSync(out, JSON.stringify({
      gsm, variant: v, father: banner ? 'reconstructed' : 'measured',
      stage: res.stage?.stage, ploidy: res.stage?.ploidy, refusedBy: res.stage?.refusedBy ?? '',
      zygosity: res.zygosity, integrity: res.integrity?.level ?? '',
      chroms: (res.chroms ?? []).filter((c: any) => /^\d+$/.test(c.chrom)).map((c: any) => [c.chrom,
        c.copyState?.state ?? '', brief(c.aneuploidy), brief(c.lostParent), brief(c.aneuploidyParent), brief(c.alleleGain)]),
      segments: (res.segments ?? []).map((s: any) => `${s.kind}:${s.chrom}`).sort(),
      findings: (res.findings ?? []).map((f: any) => `${f.cls}:${f.chrom}`).sort(),
      refusals: (res.refusals ?? []).map((r: any) => `${r.claim}:${r.chrom ?? ''}:${r.cause}`).sort(),
    }) + '\n')
  }
}
appendFileSync(out, 'DONE\n')

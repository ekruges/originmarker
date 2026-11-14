/**
 * Extra chromosomes of known parental origin, built on real arrays, to test which parent the tool
 * names. See audit/one-parent-gain-preregistration.txt.
 *
 * Five GSE148488 cells Table S3 publishes with a normal karyotype. On each of six chromosomes
 * published normal there, one at a time, an extra copy is added by re-weighting the cell's own allele
 * intensities, the way audit/gain-positive-control.ts built its control: A = (1 - BAF) 2^logR and
 * B = BAF 2^logR are the two allele intensities, one copy's worth, (A + B) / 2, is added to the allele
 * the duplicated parent carries where that parent is homozygous, and split evenly where it is not,
 * and BAF and log R are recomputed. The array's own noise stays.
 *
 * Each construction, and each cell as read, is scored by the shipped scoreSample twice: against the
 * father reconstructed from his pronuclei alone (OM_FATHER), and against both measured parents.
 *
 * Run: OM_TRIOS=<dir> OM_FATHER=<reconstruction> node --experimental-strip-types \
 *        audit/validate/constructed-gains.ts <shard> <shards> <out>
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
const CELLS = ['GSM4472410', 'GSM4472421', 'GSM4472448', 'GSM4552365', 'GSM4552369']
// Pass 1 of audit/one-parent-gain-preregistration.txt; OM_CHROMS gives the held-out set.
const CHROMS = (process.env.OM_CHROMS ?? '2,5,9,13,18,21').split(',')

const file = (gsm: string) => join(TRIOS, readdirSync(TRIOS).find((f) => f.startsWith(gsm) && /probes/.test(f))!)
const text = (p: string) => { const raw = readFileSync(p); return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8') }
function rows(path: string): any[] {
  const lines = text(path).split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = ingest.headerMap(lines[h])
  const o: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) o.push(r) }
  return o
}
const truth = new Map<string, Map<string, string>>()
for (const l of readFileSync(new URL('../zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  const p = l.split('\t'); if (!p[0]?.startsWith('GSM')) continue
  if (!truth.has(p[0])) truth.set(p[0], new Map()); truth.get(p[0])!.set(p[3], (p[4] ?? '').trim().toLowerCase())
}
const manifest = text(join(TRIOS, 'trio_manifest_full.csv')).split('\n')
const head = manifest[0].split(',')
const motherOf = new Map<string, string>()
for (const l of manifest.slice(1)) {
  const c = l.split(','); const m = (c[head.indexOf('mother_gsms')] ?? '').split(/[; ]+/)[0]
  if (m) motherOf.set(c[head.indexOf('gsm')], m)
}

const panel = keys.parsePanel(text(new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
function parentIndex(path: string) {
  const banner = inferredArray.parseInferredBanner(text(path).slice(0, inferredArray.MARK_BYTES)) ?? undefined
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  const gt = new Map<string, string>()
  for (const r of rows(path)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectParentRow(r, acc)
    gt.set(r.probesetId, r.genotype)
  }
  return { idx: score.finishParent(acc, ingest.finishProfile(path, byChrom, bafSums, first).build.build, banner), gt }
}

const measuredFather = parentIndex(file(FATHER))
const rebuiltFather = parentIndex(process.env.OM_FATHER!)

/** One extra copy of `parent`'s chromosome `chrom`, by re-weighting the allele intensities. */
function construct(cell: any[], chrom: string, parentGt: Map<string, string>): any[] {
  return cell.map((r) => {
    if (r.chrom !== chrom || r.baf === null || r.log2R === null || !Number.isFinite(r.baf) || !Number.isFinite(r.log2R)) return r
    const total = 2 ** r.log2R
    let A = (1 - r.baf) * total, B = r.baf * total
    const one = (A + B) / 2
    const g = parentGt.get(r.probesetId)
    if (g === 'AA') A += one; else if (g === 'BB') B += one; else { A += one / 2; B += one / 2 }
    return { ...r, log2R: Math.log2(A + B), baf: B / (A + B) }
  })
}

async function scoreRows(rs: any[], pat: any, mat: any, name: string) {
  const acc = score.emptyCollected(pat, mat, panel as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rs) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectRow(r, pat, mat, acc as never)
  }
  const profile = ingest.finishProfile(name, byChrom, bafSums, first)
  return await score.scoreSample({ acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: name,
    log: () => {} } as never) as any
}

const summary = (res: any, chrom: string) => {
  const c = (res.chroms ?? []).find((x: any) => x.chrom === chrom) ?? {}
  const others = (res.chroms ?? []).filter((x: any) => x.chrom !== chrom && x.alleleGain?.verdict === 'gain')
    .map((x: any) => `${x.chrom}:${x.alleleGain.parent}`)
  return {
    stage: res.stage?.stage, zygosity: res.zygosity,
    gain: c.alleleGain?.verdict ?? '', parent: c.alleleGain?.parent ?? '',
    z: c.alleleGain?.z ?? null, lrrZ: c.alleleGain?.lrrZ ?? null, aneuploidy: c.aneuploidy ?? '',
    state: c.copyState?.state ?? '', r: c.copyState?.r ?? null, otherGains: others,
    withheld: (res.refusals ?? []).filter((x: any) => x.chrom === chrom && x.claim === 'alleleGain').map((x: any) => x.gate),
  }
}

for (const gsm of CELLS.filter((_, i) => i % shards === shard)) {
  const cell = rows(file(gsm))
  const mother = parentIndex(file(motherOf.get(gsm)!))
  const chroms = CHROMS.filter((c) => truth.get(gsm)?.get(c) === 'normal')
  const jobs: [string, string, any[]][] = [['-', 'none', cell]]
  for (const c of chroms) {
    jobs.push([c, 'paternal', construct(cell, c, measuredFather.gt)])
    jobs.push([c, 'maternal', construct(cell, c, mother.gt)])
  }
  for (const [chrom, side, rs] of jobs) {
    const one = await scoreRows(rs, rebuiltFather.idx, null, gsm)
    const two = await scoreRows(rs, measuredFather.idx, mother.idx, gsm)
    for (const [config, res] of [['one', one], ['two', two]] as const) {
      const targets = chrom === '-' ? chroms : [chrom]
      for (const t of targets) {
        appendFileSync(out, JSON.stringify({ gsm, chrom: t, constructed: chrom === '-' ? 'none' : side, config,
          ...summary(res, t) }) + '\n')
      }
    }
  }
}
appendFileSync(out, 'DONE\n')

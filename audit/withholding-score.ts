/**
 * One cell, scored by the shipped scoreSample against one parent file, measured or reconstructed,
 * with the shipped founder panel; printed as one JSON line of what part C compares. The parent is
 * read as cli/om.ts loadParent reads it, so a reconstructed file's banner marks it as one.
 *
 * Run: node --experimental-strip-types audit/withholding-score.ts <cell> <parent>
 */
import { readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { lineReader } from './joy-common.ts'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const keys = await import(`${W}answerKeys.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)

const [cellPath, parentPath] = process.argv.slice(2)
const panel = keys.parsePanel(gunzipSync(readFileSync(
  new URL(`../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url))).toString('utf8'))

const text = (p: string) => { const raw = readFileSync(p); return p.endsWith('.gz') ? gunzipSync(raw) : raw }
function* rows(p: string) {
  const next = lineReader(text(p))
  let map: ReturnType<typeof ingest.headerMap> = null
  for (let l = next(); l !== null; l = next()) {
    if (!map) { if (l && !l.startsWith('#')) map = ingest.headerMap(l); continue }
    const r = ingest.parseRow(l, map)
    if (r) yield r
  }
}

const pacc = score.emptyParent()
{
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rows(parentPath)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, pacc as never)
  }
  var pat = score.finishParent(pacc as never,
    ingest.finishProfile(parentPath, byChrom as never, bafSums as never, first).build.build,
    inferredArray.parseInferredBanner(text(parentPath).toString('utf8', 0, inferredArray.MARK_BYTES)) ?? undefined)
}
const acc = score.emptyCollected(pat as never, null, panel as never)
const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
for (const r of rows(cellPath)) {
  if (!first) first = r.probesetId
  ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
  score.collectRow(r as never, pat as never, null, acc as never)
}
const profile = ingest.finishProfile(cellPath, byChrom as never, bafSums as never, first)
const res = await score.scoreSample({ acc, profile, pat, mat: null, soloRole: 'paternal', sibs: [],
  sampleName: cellPath, log: () => {} }) as Record<string, any>

const brief = (v: any): string => (v === undefined ? '' : typeof v === 'object'
  ? String(v.verdict ?? v.kind ?? JSON.stringify(v).slice(0, 40)) : String(v))
console.log(JSON.stringify({
  cell: cellPath.split('/').pop(),
  parent: parentPath.split('/').pop(),
  inferred: !!(pat as any).inferred,
  stage: res.stage?.stage, ploidy: res.stage?.ploidy, refusedBy: res.stage?.refusedBy ?? '',
  dropInCorrected: res.stage?.dropInCorrected ?? null,
  integrity: res.integrity?.level ?? '',
  normal: res.answerKeys?.normalFraction ?? null,
  chroms: (res.chroms ?? []).filter((c: any) => /^\d+$/.test(c.chrom)).map((c: any) => ({
    chrom: c.chrom, state: c.copyState?.state ?? '', r: c.copyState?.r ?? null,
    aneuploidy: brief(c.aneuploidy), lostParent: brief(c.lostParent), aneuploidyParent: brief(c.aneuploidyParent),
  })),
  copySegments: (res.segments ?? []).filter((s: any) => /^copy/.test(s.kind))
    .map((s: any) => `${s.kind}:${s.chrom}:${s.refined?.startBp ?? s.startBp ?? ''}-${s.refined?.endBp ?? s.endBp ?? ''}`),
  alleleSegments: (res.segments ?? []).filter((s: any) => !/^copy/.test(s.kind)).map((s: any) => `${s.kind}:${s.chrom}`),
  refusals: (res.refusals ?? []).map((r: any) => `${r.claim}:${r.chrom ?? ''}:${r.gate}:${brief(r.withheld)}`),
  findings: (res.findings ?? []).map((f: any) => `${f.cls}:${f.chrom}:${f.wholeChromosome ? 'whole' : `${f.startBp}-${f.endBp}`}`),
}))

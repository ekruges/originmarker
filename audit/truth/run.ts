/**
 * What the shipped scoreSample reports for each array of a list, written out whole, so every stated
 * fact in audit/local/truth can be scored against it offline.
 *
 * The list is tab-separated: array id, array path, parent path, tag, and optionally the other
 * parent's path, `maternal` where the first parent is the mother, and the paths of the other arrays
 * of the same run, separated by `;`. A parent path of `-` is no parent. Each record also carries
 * `stageFree`, the stage with no parent consulted.
 *
 * Every claim is written whether the tool printed it or withheld it. A withheld claim is in
 * `refusals` with what it would have said.
 *
 * Run: node --experimental-strip-types audit/truth/run.ts <list.tsv> <shard> <shards> <out.jsonl>
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const keys = await import(`${W}answerKeys.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)
const info = await import(`${W}informationRatio.ts`)
const stageMod = await import(`${W}stage.ts`)

const [listFile, shard, shards, out] = [process.argv[2], Number(process.argv[3]), Number(process.argv[4]), process.argv[5]]
const done = new Set<string>()
if (existsSync(out)) {
  for (const l of readFileSync(out, 'utf8').split('\n')) {
    if (l.startsWith('{')) { const r = JSON.parse(l); done.add(`${r.id}|${r.tag}`) }
  }
}
const list = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t'))
  .filter((_, i) => i % shards === shard)

const text = (p: string) => { const raw = readFileSync(p); return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8') }
function rows(path: string): any[] {
  const lines = text(path).split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = ingest.headerMap(lines[h])
  const o: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) o.push(r) }
  return o
}
const panel = keys.parsePanel(text(new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
const parents = new Map<string, any>()
function parentIndex(path: string) {
  if (path === '-') return score.noParent()
  if (parents.has(path)) return parents.get(path)
  const banner = inferredArray.parseInferredBanner(text(path).slice(0, inferredArray.MARK_BYTES)) ?? undefined
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rows(path)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectParentRow(r, acc)
  }
  const idx = score.finishParent(acc, ingest.finishProfile(path, byChrom, bafSums, first).build.build, banner)
  parents.set(path, idx)
  return idx
}

const seg = (s: any) => ({ chrom: s.chrom, kind: s.kind, start: s.startBp, end: s.endBp, absentParent: s.absentParent,
  shortBy: s.shortBy,
  score: Number.isFinite(s.score) ? +s.score.toFixed(1) : undefined, lrrShift: Number.isFinite(s.lrrShift) ? +s.lrrShift.toFixed(3) : undefined,
  lrrZ: Number.isFinite(s.lrrZ) ? +s.lrrZ.toFixed(2) : undefined })
const held = (w: any) => (w && typeof w === 'object' && 'kind' in w && 'startBp' in w ? seg(w)
  : w && typeof w === 'object' ? (w.verdict ?? w.cls ?? w.state ?? 'object') : w)

// The called autosomal genotypes of another array of the same run, as the page keeps them.
const callsOf = new Map<string, Map<string, string>>()
function calls(path: string): Map<string, string> {
  let m = callsOf.get(path)
  if (!m) {
    m = new Map()
    for (const r of rows(path)) if (/^\d+$/.test(r.chrom) && r.genotype !== 'NC') m.set(r.probesetId, r.genotype)
    if (callsOf.size > 12) callsOf.clear()
    callsOf.set(path, m)
  }
  return m
}

for (const [id, path, parentPath, tag, otherPath, role, others] of list) {
  if (done.has(`${id}|${tag}`)) continue
  try {
    const pat = parentIndex(parentPath)
    const mat = otherPath ? parentIndex(otherPath) : null
    const acc = score.emptyCollected(pat, mat, panel as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    const sameHom = keys.emptySameHom(); const homTally = info.emptyHomTally()
    for (const r of rows(path)) {
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectRow(r, pat, mat, acc as never)
      if (!/^\d+$/.test(r.chrom)) continue
      const p = panel.get(r.probesetId)
      if (p) { keys.addSameHom(sameHom, r.genotype, p, p); info.addHomozygous(homTally, r.genotype, p) }
    }
    const profile = ingest.finishProfile(id, byChrom, bafSums, first)
    // The stage with no parent consulted: the founder panel's drop-in and information ratio only.
    const free = stageMod.inferStage({ ...profile, dropIn: keys.onParentScale(keys.dropIn(sameHom), panel.calibration),
      information: info.informationRatio(homTally) }) as any
    const res = await score.scoreSample({ acc, profile, pat, mat, soloRole: role === 'maternal' ? 'maternal' : 'paternal',
      sibs: (others ?? '').split(';').filter(Boolean).map(calls),
      sampleName: id, log: () => {} } as never) as any
    appendFileSync(out, JSON.stringify({
      id, tag, parent: parentPath.split('/').pop(),
      callRate: profile.callRate, hetRate: profile.hetRate,
      stage: { stage: res.stage?.stage, refusedBy: res.stage?.refusedBy ?? '', ploidy: res.stage?.ploidy ?? '',
        dropIn: res.stage?.dropInCorrected ?? null },
      stageFree: { stage: free.stage, refusedBy: free.refusedBy ?? '', ploidy: free.ploidy ?? '',
        dropIn: free.dropInCorrected ?? null },
      information: info.informationRatio(homTally),
      parentGenome: res.parentGenome ?? null, parentGenomeOther: res.parentGenomeOther ?? null,
      zygosity: res.zygosity, originClass: res.originClass,
      sex: res.sexChromosome?.constitution ?? '',
      chroms: (res.chroms ?? []).map((c: any) => ({ chrom: c.chrom, aneuploidy: c.aneuploidy ?? '',
        lostParent: c.lostParent ?? '', aneuploidyParent: c.aneuploidyParent ?? '', gain: c.alleleGain?.verdict ?? '', gainParent: c.alleleGain?.parent ?? '',
        lrrShift: c.lrrShift ?? null, lrrZ: c.lrrZ ?? null, state: c.copyState?.state ?? '' })),
      segments: (res.segments ?? []).map(seg),
      findings: (res.findings ?? []).map((f: any) => ({ cls: f.cls, chrom: f.chrom, start: f.startBp, end: f.endBp,
        whole: !!f.wholeChromosome, originBlocked: f.originBlocked ? true : undefined, flag: f.flag })),
      origins: [...(res.oneParent ?? []).map((o: any) => ({ by: 'one-parent', where: o.where, verdict: o.verdict, parent: o.parent ?? null })),
        ...(res.dosageCalls ?? []).map((d: any) => ({ by: 'dosage', where: d.where, verdict: d.verdict }))],
      refusals: (res.refusals ?? []).map((x: any) => ({ claim: x.claim, chrom: x.chrom ?? '', gate: x.gate,
        cause: x.cause, withheld: held(x.withheld) })),
    }) + '\n')
  } catch (e) {
    appendFileSync(out, JSON.stringify({ id, tag, error: String(e).slice(0, 300) }) + '\n')
  }
}
appendFileSync(out, `DONE ${shard}\n`)

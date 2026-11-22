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
 * A file whose markers are not the proven chip's is not scored: its record carries `chip` with the
 * refusal and the causes, as the page prints them, and nothing else. See web/src/platform.ts.
 *
 * `genomeVerdicts` holds the loaded parent's genome-wide verdict and, with both parents, the
 * mother's, read as the page reads it; `changes` holds each printed change as the result card lists
 * it, with its band and the second test's verdict over it (pageSays.ts).
 *
 * With OM_SURFACES=1 each row also checks that the card, the PDF and the figure say one count, one set
 * of coordinates and one pair of numbers (pageSays.surfaceProblems), one line per row on stderr.
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
const platform = await import(`${W}platform.ts`)
const parentage = await import(`${W}parentage.ts`)
const pageSays = await import(`${W}pageSays.ts`)
const sexing = await import(`${W}sexing.ts`)
const sexMaterial = await import(`${W}sexMaterial.ts`)

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
  if (!map) throw new ChipRefused(path, { text: platform.unreadableLines(lines, panel.freqB),
    causes: platform.missingRequired(lines[h]) })
  const o: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) o.push(r) }
  checkChip(path, o, map)
  return o
}
const panel = keys.parsePanel(text(new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
class ChipRefused extends Error {
  file: string; call: any
  constructor(file: string, call: any) { super(call.text); this.file = file; this.call = call }
}
/** The same recognition the page runs, from the rows already read. Throws on another chip. */
function checkChip(path: string, rs: any[], map: any) {
  const b = ingest.emptyBuildSums()
  for (const r of rs) ingest.accumulateBuild(r, b)
  const call = platform.recogniseRows(rs, ingest.buildVerdict(b), map, panel.freqB)
  if (!call.proven) throw new ChipRefused(path, call)
}
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

// The called autosomal genotypes of another array of the same run, as the page keeps them: read by
// probe id, held as one byte per probe so a whole run's arrays fit in memory together.
const probeIndex = new Map<string, number>()
const CODES = ['', 'AA', 'AB', 'BB']
const callsOf = new Map<string, { get(id: string): string | undefined }>()
function calls(path: string): { get(id: string): string | undefined } {
  let m = callsOf.get(path)
  if (!m) {
    const rs = rows(path).filter((r) => /^\d+$/.test(r.chrom))
    for (const r of rs) if (!probeIndex.has(r.probesetId)) probeIndex.set(r.probesetId, probeIndex.size)
    const code = new Uint8Array(probeIndex.size)
    for (const r of rs) if (r.genotype !== 'NC') code[probeIndex.get(r.probesetId)!] = CODES.indexOf(r.genotype)
    m = { get: (id: string) => { const i = probeIndex.get(id); return i === undefined || i >= code.length ? undefined : CODES[code[i]] || undefined } }
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
    // The mother's genome-wide verdict, read as the page reads it, so the headline can be replayed.
    const y = sexing.sexCall((acc as any).sex)
    const maternal = (acc as any).tm && mat ? parentage.classify((acc as any).tm, mat.heterozygosity, { role: 'maternal', yBearing: y.yBearing }) : null
    if (process.env.OM_SURFACES && !res.noParent) {
      const cls = maternal ? parentage.pair(res, maternal, NaN).originClass : res.originClass
      const bad = pageSays.surfaceProblems(res, maternal ?? undefined, cls, tag.startsWith('recon'))
      process.stderr.write(`SURFACES ${tag} ${bad.length ? bad.join(' | ') : 'ok'}\n`)
    } else if (process.env.OM_SURFACES) {
      const bad = pageSays.surfaceProblems(res, undefined, res.originClass, false)
      process.stderr.write(`SURFACES ${tag} ${bad.length ? bad.join(' | ') : 'ok'}\n`)
    }
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
      // The chrY reading and what the page prints of sex: the line's reason, the sample's sex word,
      // the sperm type and the parent checks.
      y: { yBearing: y.yBearing, callRatio: y.callRatio, lrrShift: y.lrrShift, yBafShare: y.yBafShare ?? null,
        bafFloor: y.bafFloor ?? null, why: y.why ?? null },
      sexWhy: res.sexChromosome?.why ?? '', sexWord: sexMaterial.shownSex(profile.sex, res.sexMaterial, res.sexChromosome),
      spermType: res.spermType ?? null, parentSanity: res.parentSanity ?? null, matSanity: res.parentSanityOther ?? null,
      parentSex: { pat: pat.sex ?? null, mat: mat?.sex ?? null },
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
      genomeVerdicts: { loaded: res.verdict ?? null, mother: maternal?.verdict ?? null },
      // Each printed change with the band it carries and the second test's verdict over it.
      changes: pageSays.defectsForResult(res).map((d: any) => ({ chrom: d.chrom, start: d.startBp, end: d.endBp,
        kind: d.kind, origin: d.origin, step: d.step ?? null, band: d.band ?? null, second: d.second?.verdict ?? null })),
      ...(res.motherFromRun ? { motherFromRun: { ...res.motherFromRun,
        drawArrays: res.motherFromRun.drawArrays.map((i: number) => (others ?? '').split(';').filter(Boolean)[i]) } } : {}),
    }) + '\n')
  } catch (e) {
    if (e instanceof ChipRefused) {
      appendFileSync(out, JSON.stringify({ id, tag, chip: { proven: false, file: e.file.split('/').pop(),
        knownShare: e.call.knownShare ?? null, causes: e.call.causes, text: e.call.text } }) + '\n')
    } else appendFileSync(out, JSON.stringify({ id, tag, error: String(e).slice(0, 300) }) + '\n')
  }
}
appendFileSync(out, `DONE ${shard}\n`)

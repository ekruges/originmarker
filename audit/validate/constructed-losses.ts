/**
 * Losses of one parent's copy on a segment, built at a known place and size on real arrays, to
 * measure how small a loss the tool finds and how far its edges sit from the truth.
 *
 * A cell Table S3 publishes with a monosomy of a named parent holds one copy of that chromosome,
 * with its own amplification artefact. Its rows over [start, end) replace the same rows of a cell of
 * the same father whose chromosome is published normal; log R is moved by the difference of the two
 * cells' autosomal medians, so the stretch keeps its shift against its own genome. The paternal
 * haplotype of the stretch may differ from the recipient's, which a reading against the father's
 * genotype does not see.
 *
 *   direction `other`   the donor lost the MOTHER's copy: with the father loaded, the other parent's
 *   direction `own`     the donor lost the FATHER's copy: the loaded parent's own
 *   direction `control` the donor's chromosome is published normal: nothing was lost
 *
 * Recipients are the cells with at least MIN_NORMAL autosomes published normal, on the chromosomes
 * where the cell as read prints no run in either arm. The donor for a chromosome is the cell with that
 * monosomy that has the most chromosomes published normal; a control donor is a recipient read as two
 * genomes, with no drop-in, that prints no run. A stretch starts at a fixed share of the chromosome's markers, so it does not
 * open in a stretch the array does not cover. Each construction, and each recipient as read, is
 * scored by the shipped scoreSample against the measured father and against the father reconstructed
 * from his pronuclei (OM_FATHER).
 *
 * With OM_SISTERS=1 the recipients are those with another biopsy of their embryo on an array, and
 * each is scored with one array per other biopsy as the other arrays of its run: the question is
 * what a sister cell adds under 3 Mb, so only the mother's copy is built, at 3 Mb and under.
 *
 * Run: OM_TRIOS=<dir> OM_FATHER=<reconstruction> [OM_SISTERS=1] node --experimental-strip-types \
 *        audit/validate/constructed-losses.ts <shard> <shards> <out>
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
/** Autosomes published normal before a cell is a recipient. */
const MIN_NORMAL = 19
const SISTERS = process.env.OM_SISTERS === '1'
const SIZES_MB = (process.env.OM_SIZES ?? (SISTERS ? '1,1.5,2,3' : '1,1.5,2,3,5,10')).split(',').filter(Boolean).map(Number)
/** Where a stretch starts, as a share of the chromosome's markers. */
const PLACES = [0.15, 0.45, 0.75]
const CONTROL_MB = 10

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
/** gsm -> chromosome -> what Table S3 states. */
const stated = new Map<string, Map<string, string>>()
for (const l of readFileSync(new URL('../zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  const p = l.split('\t'); if (!p[0]?.startsWith('GSM')) continue
  if (!stated.has(p[0])) stated.set(p[0], new Map()); stated.get(p[0])!.set(p[3], (p[4] ?? '').trim().toLowerCase())
}
const normals = (gsm: string) => [...(stated.get(gsm)?.values() ?? [])].filter((v) => v === 'normal').length
const hasFile = (gsm: string) => readdirSync(TRIOS).some((f) => f.startsWith(gsm) && /probes/.test(f))
/** gsm -> [embryo, biopsy] from the array titles, where a title names both. */
const biopsyOf = new Map<string, [string, string]>()
{
  const lines = text(join(TRIOS, 'trio_manifest_full.csv')).split('\n')
  const head = lines[0].split(',')
  for (const l of lines.slice(1)) {
    const c = l.split(',')
    const title = c[head.indexOf('title')] ?? ''
    const e = /embryo ?([A-Za-z0-9]+)/.exec(title); const b = /biopsy ?(\d+)/.exec(title)
    if (e && b) biopsyOf.set(c[head.indexOf('gsm')], [`${c[head.indexOf('egg_donor')]}:${e[1].toUpperCase()}`, b[1]])
  }
}
/** One array for each other biopsy of the same embryo. */
function sistersOf(gsm: string): string[] {
  const me = biopsyOf.get(gsm)
  if (!me) return []
  const one = new Map<string, string>()
  for (const [g, [embryo, biopsy]] of [...biopsyOf].sort((a, b) => a[0].localeCompare(b[0]))) {
    if (embryo === me[0] && biopsy !== me[1] && hasFile(g) && !one.has(biopsy)) one.set(biopsy, g)
  }
  return [...one.values()]
}
const CELLS = [...stated.keys()].filter((g) => normals(g) >= MIN_NORMAL && hasFile(g)
  && (!SISTERS || sistersOf(g).length > 0)).sort()
/** The cell stating `what` on `chrom` with the most chromosomes published normal. */
function donorOf(chrom: string, what: string, not: string): string | undefined {
  return [...stated.keys()].filter((g) => g !== not && hasFile(g) && stated.get(g)!.get(chrom) === what)
    .sort((a, b) => normals(b) - normals(a) || a.localeCompare(b))[0]
}
/** Chromosomes with a donor of each kind, the two with the most cells stating that monosomy. */
function chromsOf(parent: 'mat' | 'pat'): string[] {
  const n = new Map<string, number>()
  for (const m of stated.values()) for (const [c, v] of m) if (v === `monosomy ${c} (${parent})`) n.set(c, (n.get(c) ?? 0) + 1)
  // Chromosome 6 carries the cut site of the experiment and is left out of both.
  return [...n].filter(([c]) => c !== '6').sort((a, b) => b[1] - a[1] || Number(a[0]) - Number(b[0])).slice(0, 2).map(([c]) => c)
}

const panel = keys.parsePanel(text(new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
function parentIndex(path: string) {
  const banner = inferredArray.parseInferredBanner(text(path).slice(0, inferredArray.MARK_BYTES)) ?? undefined
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rows(path)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectParentRow(r, acc)
  }
  return score.finishParent(acc, ingest.finishProfile(path, byChrom, bafSums, first).build.build, banner)
}
const fathers = { measured: parentIndex(file(FATHER)), rebuilt: parentIndex(process.env.OM_FATHER!) }

const median = (xs: number[]) => { const q = [...xs].sort((a, b) => a - b); return q.length ? q[q.length >> 1] : NaN }
const centre = (rs: any[]) => median(rs.filter((r) => /^\d+$/.test(r.chrom) && r.log2R !== null && Number.isFinite(r.log2R)).map((r) => r.log2R))
const donorCache = new Map<string, { byId: Map<string, any>; centre: number }>()
function donor(gsm: string) {
  let d = donorCache.get(gsm)
  if (!d) {
    const rs = rows(file(gsm))
    d = { byId: new Map(rs.map((r) => [r.probesetId, r])), centre: centre(rs) }
    donorCache.clear()   // one donor's rows in memory at a time
    donorCache.set(gsm, d)
  }
  return d
}

/** `cell` with its rows on `chrom` in [start, end) replaced by the donor's. */
function splice(cell: any[], cellCentre: number, d: { byId: Map<string, any>; centre: number }, chrom: string, start: number, end: number) {
  let n = 0; let first = Infinity; let last = -Infinity
  const rs = cell.map((r) => {
    if (r.chrom !== chrom || r.pos < start || r.pos >= end) return r
    const x = d.byId.get(r.probesetId)
    if (!x) return r
    n += 1
    if (r.pos < first) first = r.pos
    if (r.pos > last) last = r.pos
    return { ...r, genotype: x.genotype, baf: x.baf,
      log2R: x.log2R === null || !Number.isFinite(x.log2R) ? x.log2R : x.log2R - d.centre + cellCentre }
  })
  return { rs, n, first, last }
}

/** The called autosomal genotypes of an array, as the page keeps them for the other arrays of a run. */
const calls = (gsm: string) => new Map<string, string>(rows(file(gsm))
  .filter((r) => /^\d+$/.test(r.chrom) && r.genotype !== 'NC').map((r) => [r.probesetId, r.genotype]))

async function scoreRows(rs: any[], pat: any, name: string, sibs: Map<string, string>[] = []) {
  const acc = score.emptyCollected(pat, null, panel as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); const builds = ingest.emptyBuildSums(); let first = ''
  for (const r of rs) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectRow(r, pat, null, acc as never)
    ingest.accumulateBuild(r as never, builds)
  }
  const profile = ingest.finishProfile(name, byChrom, bafSums, first, builds)
  return await score.scoreSample({ acc, profile, pat, mat: null, soloRole: 'paternal', sibs, sampleName: name,
    log: () => {} } as never) as any
}

/** Every parental-absence run of a result, printed or withheld, with the parent it names. */
function absences(res: any): { chrom: string; start: number; end: number; parent: string; score: number; printed: boolean; lrrShift: number | null; lrrZ: number | null; shortBy: string | null; edges: number[] | null }[] {
  const of = (s: any, printed: boolean) => ({ chrom: s.chrom, start: s.startBp, end: s.endBp, parent: s.absentParent ?? 'paternal',
    score: s.score, printed, lrrShift: Number.isFinite(s.lrrShift) ? s.lrrShift : null, lrrZ: Number.isFinite(s.lrrZ) ? s.lrrZ : null,
    shortBy: s.shortBy ?? null,
    // The interval around each edge, where the run carries one: start low, start high, end low, end high.
    edges: s.refined ? [s.refined.startLoBp, s.refined.startHiBp, s.refined.endLoBp, s.refined.endHiBp] : null })
  return [
    ...(res.segments ?? []).filter((s: any) => s.kind === 'parental-absence').map((s: any) => of(s, true)),
    ...(res.refusals ?? []).filter((x: any) => x.claim === 'segment' && x.withheld?.kind === 'parental-absence')
      .map((x: any) => of(x.withheld, false)),
  ]
}
const key = (a: { chrom: string; start: number; end: number; parent: string }) => `${a.chrom}:${a.start}-${a.end}:${a.parent}`

/**
 * Where a control stretch comes from: recipients read as two genomes with no drop-in that print no
 * run against either father. An array the tool withholds everything on prints no run either, and
 * is not a normal stretch to splice.
 */
const cleanCells: string[] = []
for (const gsm of CELLS) {
  if (cleanCells.length >= 2) break
  const rs = rows(file(gsm))
  let clean = true
  for (const pat of Object.values(fathers)) {
    const res = await scoreRows(rs, pat, gsm)
    clean &&= res.zygosity === 'diploid' && res.stage?.dropInCorrected === undefined
      && !absences(res).some((a) => a.printed)
  }
  if (clean) cleanCells.push(gsm)
}

const jobs: { direction: 'other' | 'own' | 'control'; chrom: string; lost: string }[] = [
  ...chromsOf('mat').map((c) => ({ direction: 'other' as const, chrom: c, lost: 'maternal' })),
  ...(SISTERS ? [] : chromsOf('pat').map((c) => ({ direction: 'own' as const, chrom: c, lost: 'paternal' }))),
]
for (const gsm of CELLS.filter((_, i) => i % shards === shard)) {
  const cell = rows(file(gsm))
  const cellCentre = centre(cell)
  const sibs = SISTERS ? sistersOf(gsm).map(calls) : []
  const base: Record<string, Set<string>> = {}
  const busy = new Set<string>()
  for (const [arm, pat] of Object.entries(fathers)) {
    const res = await scoreRows(cell, pat, gsm, sibs)
    base[arm] = new Set(absences(res).map(key))
    for (const a of absences(res)) busy.add(a.chrom)
    appendFileSync(out, JSON.stringify({ gsm, arm, direction: 'as read', stage: res.stage?.stage, zygosity: res.zygosity,
      dropIn: res.stage?.dropInCorrected ?? null, absences: absences(res).length, sisters: sibs.length,
      // Runs under the span floor on the cell as read: printed on a second line, or withheld.
      short: absences(res).filter((a) => a.end - a.start < 3e6 && a.parent === 'maternal')
        .map((a) => `${a.chrom}:${(a.start / 1e6).toFixed(1)}-${(a.end / 1e6).toFixed(1)}:${a.printed ? a.shortBy : 'withheld'}`) }) + '\n')
  }
  for (const job of jobs) {
    if (stated.get(gsm)?.get(job.chrom) !== 'normal' || busy.has(job.chrom)) continue
    const places = cell.filter((r) => r.chrom === job.chrom).map((r) => r.pos).sort((a, b) => a - b)
    const kinds: { direction: string; from: string | undefined; sizes: number[] }[] = [
      { direction: job.direction, from: donorOf(job.chrom, `monosomy ${job.chrom} (${job.lost === 'maternal' ? 'mat' : 'pat'})`, gsm), sizes: SIZES_MB },
      { direction: 'control', from: cleanCells.find((g) => g !== gsm), sizes: [CONTROL_MB] },
    ]
    for (const kind of kinds) {
      if (!kind.from) continue
      const d = donor(kind.from)
      for (const place of PLACES) {
        for (const mb of kind.sizes) {
          const start = places[Math.floor(place * places.length)]; const end = start + mb * 1e6
          const { rs, n, first, last } = splice(cell, cellCentre, d, job.chrom, start, end)
          for (const [arm, pat] of Object.entries(fathers)) {
            const res = await scoreRows(rs, pat, gsm, sibs)
            const all = absences(res)
            // The run of the lost parent that covers most of the stretch, if any.
            const hits = all.filter((a) => a.chrom === job.chrom && a.parent === job.lost && a.start < end && a.end > start)
              .sort((a, b) => (Math.min(b.end, end) - Math.max(b.start, start)) - (Math.min(a.end, end) - Math.max(a.start, start)))
            const hit = hits[0]
            appendFileSync(out, JSON.stringify({
              gsm, arm, direction: kind.direction, lost: job.lost, chrom: job.chrom, donor: kind.from, mb, start, end, rows: n,
              first, last,
              found: hit ? (hit.printed ? 'printed' : 'withheld') : 'missed',
              runStart: hit?.start ?? null, runEnd: hit?.end ?? null, score: hit?.score ?? null, lrrShift: hit?.lrrShift ?? null,
              lrrZ: hit?.lrrZ ?? null, shortBy: hit?.shortBy ?? null, edges: hit?.edges ?? null,
              // Runs that the recipient as read does not have and that are not the one built.
              extra: all.filter((a) => a !== hit && !base[arm].has(key(a))).map((a) => `${a.chrom}:${(a.start / 1e6).toFixed(1)}-${(a.end / 1e6).toFixed(1)}:${a.parent}:${a.printed ? 'printed' : 'withheld'}`),
            }) + '\n')
          }
        }
      }
    }
  }
}
appendFileSync(out, 'DONE\n')

/**
 * PARENT OF ORIGIN FOR EVERY REGION IN ONE EXPERIMENT, FROM ONE GENOTYPED PARENT.
 *
 * This is the run PLAN.txt calls for and the one that had never been done: the method existed and
 * had only ever been exercised on a constructed CEPH trio. Here it meets the laboratory arrays.
 *
 *   usage: node --experimental-strip-types audit/asymmetry/one_parent_cohort.ts <dir> [out.json]
 *
 * EVERYTHING GOES THROUGH THE SHIPPED MODULES. Every artefact recorded in AUDIT-4.1 came from
 * re-deriving an analysis outside the app, where the guards do not exist. Six one-directional
 * results were produced that way and all six were wrong. Nothing is re-implemented here.
 *
 * THE ORDER MATTERS, because each step is only meaningful if the one before it held.
 *
 *   1  QUALITY, from the shipped stage inference. A failed amplification imitates haploid material
 *      and, before the gate added in 4.8.1, was handed a confident dropout. PLAN.txt names four
 *      arrays in this experiment reading 0.53-0.56 heterozygous at 0.45-0.59 call rate: no genome
 *      does that. They are excluded by the gate rather than by name, so the exclusion is a
 *      measurement rather than a list I curated.
 *
 *   2  THE REFERENCE MUST BE BULK. Measured against a single-cell reference the obligate-het
 *      statistic does not separate ploidy at all. A diploid on this platform reads 0.15-0.17
 *      heterozygous, NOT the 0.30 of a dense common-SNP panel, so a threshold carried from another
 *      array calls every diploid here haploid.
 *
 *   3  THE RELATIONSHIP MUST BE ESTABLISHED, NOT ASSUMED. Opposite-homozygote rate alone is near
 *      zero for a parent AND for the same genome twice, which is how a previous screen returned
 *      77% wrong. Biparental confirmation via hetCall is required on top of it.
 *
 *   4  REGIONS FROM THE INTENSITY CHANNEL ONLY. Amplification dropout removes genotype calls
 *      without removing DNA, so a region derived from heterozygosity is a dropout report.
 *
 *   5  ORIGIN, from callOneParentOrigin, parameterised by the dropout the stage inferred.
 *
 * MEMORY IS THE REASON FOR THE PASS STRUCTURE. Each array is 825,657 markers, and holding the
 * genotype map, the position map and the intensity channel for every array at once reached 17 GB
 * on a 16 GB machine in earlier work here. So: pass one profiles every array and keeps only the
 * profile, pass two loads the single reference, and pass three streams each candidate through
 * linkage, segmentation and the origin call before discarding it. Each file is read twice and at
 * most two arrays are resident.
 *
 * WHAT THIS CANNOT SETTLE, and it is the whole reason PLAN.txt exists. If maternal and paternal
 * events are detected with unequal sensitivity, true parity reads as asymmetry: on one cell per
 * embryo a 7.7% maternal shortfall alone manufactures an apparent paternal fraction of 0.52. This
 * run produces per-region calls on ONE experiment. It is the unit that lets the method be checked
 * against real material; it is not the cohort a claim would rest on.
 */
import { readFileSync, writeFileSync, readdirSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../../web/src/', import.meta.url).pathname
const { headerMap, parseRow, emptyBafSums, accumulateBaf, accumulate, finishProfile } =
  await import(`${W}ingest.ts`)
const { isAutosome } = await import(`${W}parentage.ts`)
const { hetCall, addOneParent, emptyHet } = await import(`${W}obligateHet.ts`)
const { scanCopyNumber, externalNull, segmentCoords } = await import(`${W}segments.ts`)
const { callOneParentOrigin } = await import(`${W}oneParentOrigin.ts`)
const { inferStage, locus } = await import(`${W}stage.ts`)

type AB = 'AA' | 'AB' | 'BB' | 'NC'
type Profile = ReturnType<typeof finishProfile>
type Stage = ReturnType<typeof inferStage>

const DIR = process.argv[2]
const OUT = process.argv[3] ?? 'audit/asymmetry/one-parent-cohort.json'
/**
 * The reference array, optionally named.
 *
 * WHICH SAMPLE IS THE PARENT IS METADATA, NOT A MEASUREMENT. Call rate separates a bulk somatic
 * sample from an amplified one, but not from a trophectoderm biopsy, which also amplifies well.
 * And heterozygosity cannot be used for it here: this experiment's cumulus sample reads 0.145
 * against 0.168 for the bulk gDNA on this platform, which the shipped stage classifier therefore
 * calls single-cell. That gap is ancestry rather than material, and it is exactly the confound
 * the stage audit quantified at +0.155 for an East Asian sample against this European anchor.
 * So the laboratory's own identification is accepted where it is given, and auto-detection only
 * proposes candidates rather than settling it silently.
 */
const REF = process.argv[4]
if (!DIR) throw new Error('usage: one_parent_cohort.ts <array-dir> [out.json] [reference-id]')

/** A reference must be bulk. These are the platform's own figures, not a general array's. */
const BULK_CALL_MIN = 0.90
const BULK_HET_MIN = 0.135
const BULK_HET_MAX = 0.190

/** Opposite-homozygote rate under which a relationship is possible. Necessary, NOT sufficient. */
const OPP_HOM_MAX = 0.020

interface CnMarker { chrom: string, pos: number, called: boolean, log2R: number | null }

interface Loaded {
  id: string
  gt: Map<string, AB>
  /** Position per marker, so a segment's markers are gathered without a second read of the file. */
  pos: Map<string, { chrom: string, pos: number }>
  /** The copy-number channel: EVERY marker, called or not, with its intensity ratio. A region
   *  that is gone stops producing calls, which a genotype-derived indicator cannot see. */
  cnByChrom: Map<string, CnMarker[]>
  profile: ReturnType<typeof finishProfile>
}

/**
 * Read one array the way the app reads it.
 *
 * The accumulators are the shipped ones. Recomputing a call rate or a heterozygous rate here, even
 * correctly, would mean the stage gate downstream is parameterised by a number the app never
 * produced, and the whole point of routing through the modules is that it cannot be.
 */
function load(path: string, id: string, full: boolean): Loaded | null {
  const lines = readFileSync(path, 'utf8').split('\n')
  let h = -1
  for (let i = 0; i < Math.min(60, lines.length); i += 1) {
    if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  }
  if (h < 0) return null
  const map = headerMap(lines[h])
  if (!map) return null

  const gt = new Map<string, AB>()
  const pos = new Map<string, { chrom: string, pos: number }>()
  const cnByChrom = new Map<string, CnMarker[]>()
  const byChrom = new Map<string, never>()
  const baf = emptyBafSums()
  let firstId = ''

  for (let i = h + 1; i < lines.length; i += 1) {
    const row = parseRow(lines[i], map)
    if (!row) continue
    if (!firstId) firstId = row.probesetId
    accumulate(row as never, byChrom as never)
    accumulateBaf(row as never, baf as never)
    if (!full || !isAutosome(row.chrom)) continue
    if (row.genotype !== 'NC') {
      gt.set(row.probesetId, row.genotype as AB)
      pos.set(row.probesetId, { chrom: row.chrom, pos: row.pos })
    }
    const cn = cnByChrom.get(row.chrom) ?? []
    cn.push({ chrom: row.chrom, pos: row.pos, called: row.genotype !== 'NC', log2R: row.log2R })
    cnByChrom.set(row.chrom, cn)
  }
  return { id, gt, pos, cnByChrom, profile: finishProfile(id, byChrom as never, baf as never, firstId) }
}

/** Opposite homozygotes: a parent and child cannot be AA and BB at the same marker. */
function oppositeHomRate(a: Map<string, AB>, b: Map<string, AB>): { rate: number, n: number } {
  let n = 0
  let opp = 0
  for (const [probe, ga] of a) {
    if (ga !== 'AA' && ga !== 'BB') continue
    const gb = b.get(probe)
    if (gb !== 'AA' && gb !== 'BB') continue
    n += 1
    if (ga !== gb) opp += 1
  }
  return { rate: n ? opp / n : NaN, n }
}

const files = readdirSync(DIR).filter((f: string) => f.endsWith('.probes')).sort()
if (!files.length) throw new Error(`no .probes files in ${DIR}`)
const idOf = (f: string) => f.replace(/_\d+\.CEL\.probes$/, '').replace(/\.probes$/, '')
process.stderr.write(`${files.length} arrays in ${DIR}\n\n`)

// --- pass 1: profile and stage every array, keeping only the profile ---------------------------
const staged: { file: string, id: string, profile: Profile, stage: Stage }[] = []
for (const f of files) {
  const L = load(join(DIR, f), idOf(f), false)
  if (!L) { process.stderr.write(`  ${idOf(f)}: unreadable\n`); continue }
  staged.push({ file: f, id: L.id, profile: L.profile, stage: inferStage(L.profile) })
}

const failed = staged.filter((s) => s.stage.stage === 'failed')
for (const s of failed) {
  process.stderr.write(`EXCLUDED ${s.id}: call ${(s.profile.callRate * 100).toFixed(1)}% `
    + `het ${(s.profile.hetRate * 100).toFixed(1)}%, amplification failure\n`)
}
const usable = staged.filter((s) => s.stage.stage !== 'failed')

// --- pass 2: the reference, found by measurement rather than by name ---------------------------
const bulkLike = usable.filter((s) => (
  s.profile.callRate >= BULK_CALL_MIN
  && s.profile.hetRate >= BULK_HET_MIN
  && s.profile.hetRate <= BULK_HET_MAX
)).sort((a, b) => b.profile.callRate - a.profile.callRate)
const named = REF ? bulkLike.filter((s) => s.id === REF) : []
if (REF && !named.length) {
  const any = usable.find((s) => s.id === REF)
  throw new Error(any
    ? `reference ${REF} is not bulk-like: call ${any.profile.callRate.toFixed(3)} `
      + `het ${any.profile.hetRate.toFixed(3)}. Measured against a single-cell reference the `
      + 'obligate-het statistic does not separate ploidy at all, so this is refused rather than run'
    : `reference ${REF} is not in ${DIR}`)
}
const refCandidates = named.length ? named : bulkLike.slice(0, 1)
process.stderr.write(`\nbulk-like candidates: ${bulkLike.length
  ? bulkLike.map((r) => r.id).join(', ') : 'NONE'}\n`)
for (const r of bulkLike) {
  process.stderr.write(`  ${r.id}: call ${r.profile.callRate.toFixed(3)} `
    + `het ${r.profile.hetRate.toFixed(3)} -> ${r.stage.stage}\n`)
}
process.stderr.write(REF ? `\nreference given: ${REF}\n`
  : `\nno reference named; using the highest call rate, ${refCandidates[0]?.id ?? 'none'}\n`)
if (!refCandidates.length) {
  writeFileSync(OUT, JSON.stringify({ dir: DIR, refs: [], excluded: failed.map((f) => f.id),
    note: 'no bulk-like array in this directory; nothing can be anchored', results: [] }, null, 1))
  process.stderr.write('\nno bulk reference here. Nothing is anchored, and nothing is guessed.\n')
  process.exit(0)
}

const results: Record<string, unknown>[] = []
for (const refMeta of refCandidates) {
  const ref = load(join(DIR, refMeta.file), refMeta.id, true)
  if (!ref) continue
  process.stderr.write(`\n=== reference ${ref.id} (${ref.gt.size} autosomal calls) ===\n`)

  // --- pass 3: one candidate at a time, all the way through -----------------------------------
  for (const cand of usable) {
    if (cand.id === refMeta.id) continue
    const c = load(join(DIR, cand.file), cand.id, true)
    if (!c) continue

    const opp = oppositeHomRate(ref.gt, c.gt)
    if (!(opp.rate <= OPP_HOM_MAX) || opp.n < 10_000) {
      process.stderr.write(`  ${c.id}: opposite-hom ${opp.rate.toFixed(4)} over ${opp.n}, `
        + 'not a child of this reference\n')
      continue
    }
    // Opposite-homozygote rate alone cannot separate a child from the SAME GENOME twice: both sit
    // near zero on this platform, which is how a screen without this second step came back 77%
    // wrong in earlier work here. At markers where the reference is HOMOZYGOUS, a child reads
    // heterozygous at the rate the OTHER parent supplies the alternative allele, while a second
    // array of the reference's own genome cannot read heterozygous at all. So the one-parent
    // ploidy call is the duplicate discriminator: 'diploid' means a second contribution exists.
    const h = emptyHet()
    for (const [probe, pg] of ref.gt) {
      const cg = c.gt.get(probe)
      if (cg) addOneParent(pg as never, cg as never, h as never)
    }
    const link = hetCall(h as never, 1) as {
      ploidy: string, fraction: number, informative: number, why: string
    }
    // 'biparental' means a second parental contribution is present, which a second array of the
    // reference's own genome cannot show. 'uniparental' is a duplicate or a haploid product of
    // this same parent, and 'uncalled' is the provisional band between them, which with one
    // parent is wide on purpose. Only the first is a child.
    if (link.ploidy !== 'biparental') {
      process.stderr.write(`  ${c.id}: SET ASIDE, opp ${opp.rate.toFixed(4)} but ${link.why}\n`)
      continue
    }

    // Regions from the INTENSITY channel, per chromosome, against a null that EXCLUDES the
    // chromosome under test. A self-derived null lets one large event set its own baseline.
    const noCall = new Map<string, [number, number]>()
    for (const [ch, ms] of c.cnByChrom) {
      noCall.set(ch, [ms.length, ms.filter((m) => !m.called).length])
    }
    const lrrAll = [...c.cnByChrom.values()].flat()
      .map((m) => m.log2R).filter((x): x is number => x !== null).sort((a, b) => a - b)
    const genomeLrr = lrrAll.length ? lrrAll[lrrAll.length >> 1] : 0
    const segs = [...c.cnByChrom]
      .flatMap(([ch, ms]) => scanCopyNumber(ms as never, externalNull(noCall, ch), genomeLrr))

    process.stderr.write(`  ${c.id}: child of ${ref.id} (opp ${opp.rate.toFixed(4)}, `
      + `second contribution at ${(link.fraction * 100).toFixed(1)}% of `
      + `${link.informative} informative), ${cand.stage.stage}, ${segs.length} region(s)\n`)
    if (!segs.length) {
      results.push({ ref: ref.id, sample: c.id, stage: cand.stage.stage, regions: 0 })
      continue
    }

    for (const sg of segs as { chrom: string, kind: string }[]) {
      const co = segmentCoords(sg as never) as { start: number, end: number }
      const pairs: [string, string][] = []
      for (const [probe, p] of c.pos) {
        if (p.chrom !== sg.chrom || p.pos < co.start || p.pos > co.end) continue
        const pg = ref.gt.get(probe)
        const cg = c.gt.get(probe)
        if (pg && cg) pairs.push([pg, cg])
      }
      const ado = Number.isFinite(cand.stage.dropout) ? cand.stage.dropout : 0.308
      const call = callOneParentOrigin(pairs as never, ado) as {
        verdict: string, posterior: number, markers: number, exclusive: number, why: string
      }
      const row = {
        ref: ref.id,
        sample: c.id,
        stage: cand.stage.stage,
        dropout: Number.isFinite(cand.stage.dropout) ? cand.stage.dropout : null,
        kind: sg.kind,
        locus: locus(sg.chrom, co.start, co.end),
        verdict: call.verdict,
        posterior: call.posterior,
        markers: call.markers,
        exclusive: call.exclusive,
        why: call.why,
      }
      results.push(row)
      process.stderr.write(`      ${row.locus} ${sg.kind}: ${call.verdict} `
        + `(posterior ${call.posterior.toFixed(3)}, ${call.markers} informative, `
        + `${call.exclusive} exclusive)\n`)
    }
  }
}

writeFileSync(OUT, JSON.stringify({
  dir: DIR,
  refs: refCandidates.map((r) => ({ id: r.id, callRate: r.profile.callRate, hetRate: r.profile.hetRate })),
  excluded: failed.map((f) => ({ id: f.id, callRate: f.profile.callRate, hetRate: f.profile.hetRate })),
  results,
}, null, 1))

const verdicts = results.filter((r) => r.verdict && r.verdict !== 'refused')
process.stderr.write(`\n${results.length} region row(s), ${verdicts.length} with a verdict -> ${OUT}\n`)
const tally = new Map<string, number>()
for (const r of results) {
  if (!r.verdict) continue
  tally.set(r.verdict as string, (tally.get(r.verdict as string) ?? 0) + 1)
}
for (const [v, n] of [...tally].sort((a, b) => b[1] - a[1])) {
  process.stderr.write(`  ${v}: ${n}\n`)
}

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
 *   4  EVENTS, OF BOTH KINDS. Segmental change comes from the intensity channel, never from
 *      heterozygosity: amplification dropout removes genotype calls without removing DNA, so a
 *      het-derived region is a dropout report. WHOLE-CHROMOSOME aneuploidy is a separate detector
 *      and must not be omitted, which the first version of this harness did. The copy-number scan
 *      excludes whole chromosomes by construction, so a trisomy 21 produced no segment and the
 *      array came back with nothing to assign. In embryos that is where most of the abnormality
 *      is, and it is the thing the question is about.
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
import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join, basename } from 'node:path'

const W = new URL('../../web/src/', import.meta.url).pathname
const { headerMap, parseRow, emptyBafSums, accumulateBaf, accumulate, finishProfile } =
  await import(`${W}ingest.ts`)
const { isAutosome, emptyTally, tallyRow, classify } = await import(`${W}parentage.ts`)
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
if (!DIR) throw new Error('usage: one_parent_cohort.ts <dir> [out.json] [ref-id-or-path]')

/** A reference given as a PATH comes from outside the directory, which is the usual case: the
 *  genotyped father is a public bulk array and the embryos are laboratory files. */
const REF_IS_PATH = !!REF && (REF.includes('/') || REF.endsWith('.probes'))

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
  /** chrom:pos -> probeset id, so the per-chromosome tally can be built from the copy-number
   *  channel without carrying a fourth map of the file. */
  probeAt: Map<string, string>
  profile: ReturnType<typeof finishProfile>
}

/**
 * Read one array the way the app reads it.
 *
 * The accumulators are the shipped ones. Recomputing a call rate or a heterozygous rate here, even
 * correctly, would mean the stage gate downstream is parameterised by a number the app never
 * produced, and the whole point of routing through the modules is that it cannot be.
 */
/**
 * Stride for the LINKAGE SCREEN, not for anything reported.
 *
 * Most arrays in a corpus are not children of a given reference, and a full read to discover that
 * costs the same as a full read to analyse one that is. The screen is an opposite-homozygote rate
 * and a heterozygous fraction, both of which are proportions: over 103,207 markers their standard
 * error is under 0.001, which is two orders below the 0.02 and 0.065-0.085 boundaries they are
 * compared against. Every array that passes the screen is then re-read in full, so nothing that
 * reaches a result was measured on a subsample. This is the same stride the shipped example data
 * uses, for the same reason.
 */
const SCREEN_STRIDE = 8

function load(path: string, id: string, full: boolean, stride = 1): Loaded | null {
  // Public arrays arrive gzipped and are kept that way: uncompressed the series is five times
  // the size and this machine does not have it to spare.
  const raw = readFileSync(path)
  const lines = (path.endsWith('.gz') ? gunzipSync(raw).toString('utf8') : raw.toString('utf8'))
    .split('\n')
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
  const probeAt = new Map<string, string>()
  const byChrom = new Map<string, never>()
  const baf = emptyBafSums()
  let firstId = ''

  for (let i = h + 1; i < lines.length; i += stride) {
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
    probeAt.set(`${row.chrom}:${row.pos}`, row.probesetId)
  }
  return { id, gt, pos, cnByChrom, probeAt, profile: finishProfile(id, byChrom as never, baf as never, firstId) }
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

/** Recursive, because the corpus keeps each experiment in its own subdirectory and flattening it
 *  on disk once cost 148 basename collisions. Paths are kept relative to DIR so ids stay short. */
function probesUnder(dir: string, rel = ''): string[] {
  const out: string[] = []
  for (const e of readdirSync(join(dir, rel)).sort()) {
    const r = rel ? join(rel, e) : e
    if (statSync(join(dir, r)).isDirectory()) out.push(...probesUnder(dir, r))
    else if (e.endsWith('.probes') || e.endsWith('.txt.gz')) out.push(r)
  }
  return out
}
const files = probesUnder(DIR)
if (!files.length) throw new Error(`no .probes files in ${DIR}`)
const idOf = (f: string) => basename(f)
  .replace(/_\d+\.CEL\.probes$/, '').replace(/\.probes$/, '')
  .replace(/\.CEL\.txt\.gz$/, '').replace(/^(GSM\d+)_.*$/, '$1')
process.stderr.write(`${files.length} arrays in ${DIR}\n\n`)

/**
 * Everything downstream of "these two arrays are a parent and a child", for one candidate.
 *
 * Kept as one function because the two reference modes differ only in how the reference was
 * found. Splitting the pipeline across them is how a guard gets applied on one path and not the
 * other, which is the failure this whole harness exists to avoid.
 */
function scoreCandidate(ref: Loaded, c: Loaded, stage: Stage): Record<string, unknown>[] {
  const opp = oppositeHomRate(ref.gt, c.gt)
  if (!(opp.rate <= OPP_HOM_MAX) || opp.n < 10_000) {
    process.stderr.write(`  ${c.id}: opposite-hom ${opp.rate.toFixed(4)} over ${opp.n}, `
      + 'not a child of this reference\n')
    return []
  }
  // Opposite-homozygote rate alone cannot separate a child from the SAME GENOME twice: both sit
  // near zero on this platform, which is how a screen without this second step came back 77%
  // wrong in earlier work here. At markers where the reference is HOMOZYGOUS, a child reads
  // heterozygous at the rate the OTHER parent supplies the alternative allele, while a second
  // array of the reference's own genome cannot read heterozygous at all. 'biparental' is
  // therefore the duplicate discriminator; 'uncalled' is the provisional band, wide on purpose
  // because one parent carries less information than two.
  const h = emptyHet()
  for (const [probe, pg] of ref.gt) {
    const cg = c.gt.get(probe)
    if (cg) addOneParent(pg as never, cg as never, h as never)
  }
  const link = hetCall(h as never, 1) as {
    ploidy: string, fraction: number, informative: number, why: string
  }
  if (link.ploidy !== 'biparental') {
    process.stderr.write(`  ${c.id}: SET ASIDE, opp ${opp.rate.toFixed(4)} but ${link.why}\n`)
    return []
  }

  // WHOLE-CHROMOSOME ANEUPLOIDY, which the segmental scan cannot see because it excludes whole
  // chromosomes by construction. classify supplies it from the call-rate collapse plus the
  // direction of the intensity shift, which is the shipped detector rather than a local one.
  const t = emptyTally()
  for (const [ch, ms] of c.cnByChrom) {
    for (const m of ms) {
      const probe = c.probeAt.get(`${ch}:${m.pos}`)
      tallyRow((probe ? ref.gt.get(probe) ?? 'NC' : 'NC') as never, {
        probesetId: probe ?? '', chrom: ch, pos: m.pos, log2R: m.log2R, baf: null,
        genotype: (probe ? c.gt.get(probe) ?? 'NC' : 'NC'), copyNumber: null,
      } as never, t as never)
    }
  }
  const cls = classify(t as never, ref.profile.hetRate, { role: 'paternal' }) as {
    chroms: { chrom: string, aneuploidy?: 'loss' | 'gain' }[]
  }
  const aneuploid = cls.chroms.filter((x) => x.aneuploidy)

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

  process.stderr.write(`  ${c.id}: CHILD of ${ref.id} (opp ${opp.rate.toFixed(4)}, second `
    + `contribution ${(link.fraction * 100).toFixed(1)}% of ${link.informative}), `
    + `${stage.stage}, ${segs.length} segment(s), ${aneuploid.length} whole chromosome(s)\n`)
  const out: Record<string, unknown>[] = []

  // Whole chromosomes first: they are the larger event and the one an embryo is most likely to
  // carry. Every marker on the chromosome is informative for the origin call, which is why these
  // are the best powered rows in the output rather than the worst.
  for (const a of aneuploid) {
    const pairs: [string, string][] = []
    let lo = Infinity
    let hi = 0
    for (const [probe, p] of c.pos) {
      if (p.chrom !== a.chrom) continue
      if (p.pos < lo) lo = p.pos
      if (p.pos > hi) hi = p.pos
      const pg = ref.gt.get(probe)
      const cg = c.gt.get(probe)
      if (pg && cg) pairs.push([pg, cg])
    }
    const ado = Number.isFinite(stage.dropout) ? stage.dropout : 0.308
    const call = callOneParentOrigin(pairs as never, ado) as {
      verdict: string, posterior: number, markers: number, exclusive: number, why: string
    }
    out.push({
      ref: ref.id, sample: c.id, stage: stage.stage, kind: `whole-chromosome ${a.aneuploidy}`,
      locus: Number.isFinite(lo) ? locus(a.chrom, lo, hi) : `chr${a.chrom}`,
      verdict: call.verdict, posterior: call.posterior, markers: call.markers,
      exclusive: call.exclusive, why: call.why,
    })
    process.stderr.write(`      chr${a.chrom} whole-chromosome ${a.aneuploidy}: ${call.verdict} `
      + `(posterior ${call.posterior.toFixed(3)}, ${call.markers} informative, `
      + `${call.exclusive} exclusive)\n`)
  }

  if (!segs.length && !aneuploid.length) {
    return [{ ref: ref.id, sample: c.id, stage: stage.stage, regions: 0 }]
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
    const ado = Number.isFinite(stage.dropout) ? stage.dropout : 0.308
    const call = callOneParentOrigin(pairs as never, ado) as {
      verdict: string, posterior: number, markers: number, exclusive: number, why: string
    }
    const row = {
      ref: ref.id,
      sample: c.id,
      stage: stage.stage,
      dropout: Number.isFinite(stage.dropout) ? stage.dropout : null,
      kind: sg.kind,
      locus: locus(sg.chrom, co.start, co.end),
      verdict: call.verdict,
      posterior: call.posterior,
      markers: call.markers,
      exclusive: call.exclusive,
      why: call.why,
    }
    out.push(row)
    process.stderr.write(`      ${row.locus} ${sg.kind}: ${call.verdict} `
      + `(posterior ${call.posterior.toFixed(3)}, ${call.markers} informative, `
      + `${call.exclusive} exclusive)\n`)
  }
  return out
}

const results: Record<string, unknown>[] = []
const excluded: Record<string, unknown>[] = []
const bulkLike: { id: string, callRate: number, hetRate: number, stage: string }[] = []
let refUsed = ''

if (REF_IS_PATH) {
  // ONE PASS PER ARRAY. The reference is external and loaded once, so profiling, staging,
  // linkage, segmentation and the origin call all come from a single read of each candidate and
  // nothing but the reference stays resident.
  const ref = load(REF, basename(REF).replace(/\.(probes|CEL\.txt|CEL\.txt\.gz)$/, '').replace(/^(GSM\d+)_.*$/, '$1'), true)
  if (!ref) throw new Error(`reference ${REF} could not be read`)
  const rs = inferStage(ref.profile)
  refUsed = ref.id
  process.stderr.write(`reference ${ref.id}: call ${ref.profile.callRate.toFixed(4)} `
    + `het ${ref.profile.hetRate.toFixed(4)} -> ${rs.stage}, ${ref.gt.size} autosomal calls\n\n`)
  if (rs.stage !== 'bulk') {
    throw new Error(`reference ${ref.id} stages as ${rs.stage}, not bulk. Measured against a `
      + 'single-cell reference the obligate-het statistic does not separate ploidy at all, so '
      + 'this is refused rather than run')
  }
  let screened = 0
  for (const f of files) {
    const id = idOf(f)
    const scr = load(join(DIR, f), id, true, SCREEN_STRIDE)
    screened += 1
    if (!scr) { process.stderr.write(`  ${id}: unreadable\n`); continue }
    const sst = inferStage(scr.profile)
    if (sst.stage === 'failed') {
      process.stderr.write(`  EXCLUDED ${id}: call ${(scr.profile.callRate * 100).toFixed(1)}% `
        + `het ${(scr.profile.hetRate * 100).toFixed(1)}%, ${sst.why}\n`)
      excluded.push({ id, callRate: scr.profile.callRate, hetRate: scr.profile.hetRate })
      continue
    }
    const opp = oppositeHomRate(ref.gt, scr.gt)
    if (!(opp.rate <= OPP_HOM_MAX) || opp.n < 10_000) continue
    // Passed the screen, so re-read in full. Nothing that reaches a result is measured on a
    // subsample: the screen only decides which arrays are worth reading properly.
    process.stderr.write(`  ${id}: screen opp ${opp.rate.toFixed(4)}, re-reading in full\n`)
    const c = load(join(DIR, f), id, true)
    if (!c) continue
    results.push(...scoreCandidate(ref, c, inferStage(c.profile)))
  }
  process.stderr.write(`\nscreened ${screened} arrays at stride ${SCREEN_STRIDE}\n`)
} else {
  // --- pass 1: profile and stage every array, keeping only the profile -------------------------
  const staged: { file: string, id: string, profile: Profile, stage: Stage }[] = []
  for (const f of files) {
    const L = load(join(DIR, f), idOf(f), false)
    if (!L) { process.stderr.write(`  ${idOf(f)}: unreadable\n`); continue }
    staged.push({ file: f, id: L.id, profile: L.profile, stage: inferStage(L.profile) })
  }
  for (const s of staged.filter((x) => x.stage.stage === 'failed')) {
    process.stderr.write(`EXCLUDED ${s.id}: call ${(s.profile.callRate * 100).toFixed(1)}% `
      + `het ${(s.profile.hetRate * 100).toFixed(1)}%, ${s.stage.why}\n`)
    excluded.push({ id: s.id, callRate: s.profile.callRate, hetRate: s.profile.hetRate })
  }
  const usable = staged.filter((s) => s.stage.stage !== 'failed')

  // --- pass 2: the reference. WHICH SAMPLE IS THE PARENT IS METADATA, NOT A MEASUREMENT. Call
  // rate separates a bulk somatic sample from an amplified one but not from a trophectoderm
  // biopsy, which also amplifies well, and heterozygosity cannot be used for it: one cumulus
  // sample here reads 0.145 against 0.168 for bulk gDNA on this platform, which the shipped
  // classifier therefore calls single-cell. That gap is ancestry rather than material. So a
  // named reference is accepted, and auto-detection only proposes.
  const cands = usable.filter((s) => (
    s.profile.callRate >= BULK_CALL_MIN
    && s.profile.hetRate >= BULK_HET_MIN
    && s.profile.hetRate <= BULK_HET_MAX
  )).sort((a, b) => b.profile.callRate - a.profile.callRate)
  for (const r of cands) {
    bulkLike.push({ id: r.id, callRate: r.profile.callRate, hetRate: r.profile.hetRate, stage: r.stage.stage })
    process.stderr.write(`  bulk-like: ${r.id} call ${r.profile.callRate.toFixed(3)} `
      + `het ${r.profile.hetRate.toFixed(3)} -> ${r.stage.stage}\n`)
  }
  const named = REF ? cands.filter((s) => s.id === REF) : []
  if (REF && !named.length) throw new Error(`${REF} is not a bulk-like array in ${DIR}`)
  const chosen = named.length ? named[0] : cands[0]
  if (!chosen) {
    process.stderr.write('\nno bulk reference here. Nothing is anchored, and nothing is guessed.\n')
  } else {
    const ref = load(join(DIR, chosen.file), chosen.id, true)
    if (ref) {
      refUsed = ref.id
      process.stderr.write(`\n=== reference ${ref.id} (${ref.gt.size} autosomal calls) ===\n`)
      for (const cand of usable) {
        if (cand.id === chosen.id) continue
        const c = load(join(DIR, cand.file), cand.id, true)
        if (c) results.push(...scoreCandidate(ref, c, cand.stage))
      }
    }
  }
}

writeFileSync(OUT, JSON.stringify({ dir: DIR, reference: refUsed, bulkLike, excluded, results }, null, 1))
const tally = new Map<string, number>()
for (const r of results) {
  if (r.verdict) tally.set(r.verdict as string, (tally.get(r.verdict as string) ?? 0) + 1)
}
process.stderr.write(`\n${results.length} region row(s) -> ${OUT}\n`)
for (const [v, n] of [...tally].sort((a, b) => b[1] - a[1])) {
  process.stderr.write(`  ${v}: ${n}\n`)
}

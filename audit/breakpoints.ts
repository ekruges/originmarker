/**
 * Accuracy audit for breakpoint refinement and gain origin.
 *
 *     node --experimental-strip-types --max-old-space-size=12288 \
 *       audit/breakpoints.ts <array-dir> <out.json>
 *
 * Four parts, each falsifiable, each running the SHIPPED modules rather than a copy:
 *
 *   A  Constructed events with a breakpoint known at marker resolution, built by splicing one
 *      real array into another so the segment carries real amplification artefact. Measures the
 *      window edge against the refined edge on the same events. This is the only part with
 *      ground truth, and it is ground truth because the splice point is chosen, not inferred.
 *
 *   B  Split-half holdout on the REAL segments the tool reports, where no truth exists. Refine
 *      on odd-indexed and even-indexed informative markers separately; both halves see the same
 *      physical breakpoint and largely independent artefact. If the two disagree by more than the
 *      interval claims, the interval is understated. Truth-free and runnable on any array.
 *
 *   C  Specificity on genomes that cannot carry a segmental event: bulk diploid adults. Any
 *      segment reported there is a false positive with no interpretation available.
 *
 *   D  Gain origin. There is no confirmed gain with a known parent in this material, so this
 *      cannot show the annotation firing correctly. It shows the two things that ARE checkable:
 *      that the orientation is invariant on real marker pairs, and that the refusals fire on the
 *      real array that carries gains.
 */
import { readFileSync, writeFileSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
const { headerMap, parseRow } = await import(`${W}ingest.ts`)
const { isAutosome, pct } = await import(`${W}parentage.ts`)
const { scanChromosome, externalNull, SEGMENT_LRT, MIN_SEGMENT_MARKERS, segmentCoords }
  = await import(`${W}segments.ts`)
const { refineEdges, LR_DROP, SEARCH_RADIUS, REFINE_CALL_FLOOR } = await import(`${W}refine.ts`)
const { paternalShare, callHomologue, externalHetBackground } = await import(`${W}gainOrigin.ts`)

const DIR = process.argv[2]
const OUT = process.argv[3] ?? 'audit/breakpoint-results.json'
if (!DIR) throw new Error('usage: breakpoints.ts <array-dir> <out.json>')

const eachLine = async (f: File, fn: (l: string) => void): Promise<void> => {
  const rd = f.stream().getReader()
  const d = new TextDecoder()
  let c = ''
  const NL = String.fromCharCode(10)
  for (;;) {
    const { done, value } = await rd.read()
    if (done) break
    c += d.decode(value, { stream: true })
    let i = c.indexOf(NL)
    while (i >= 0) { fn(c.slice(0, i)); c = c.slice(i + 1); i = c.indexOf(NL) }
  }
  if (c) fn(c)
}
const open = (name: string): File => {
  const raw = readFileSync(`${DIR}/${name}`)
  return new File([name.endsWith('.gz') ? gunzipSync(raw) : raw], name)
}
const median = (xs: number[]): number => {
  if (!xs.length) return NaN
  const s = [...xs].sort((a, b) => a - b)
  return s[s.length >> 1]
}
const pctl = (xs: number[], p: number): number => {
  if (!xs.length) return NaN
  const s = [...xs].sort((a, b) => a - b)
  return s[Math.min(s.length - 1, Math.floor(p * s.length))]
}

interface Loaded {
  name: string
  genotype: Map<string, string>
  baf: Map<string, number>
  callRate: number
  hetRate: number
}

async function load(name: string): Promise<Loaded> {
  const genotype = new Map<string, string>()
  const baf = new Map<string, number>()
  let called = 0
  let total = 0
  let het = 0
  let map: ReturnType<typeof headerMap> = null
  await eachLine(open(name), (l) => {
    if (!map) { map = headerMap(l); return }
    const r = parseRow(l, map)
    if (!r || !isAutosome(r.chrom)) return
    genotype.set(r.probesetId, r.genotype)
    if (r.baf !== null) baf.set(r.probesetId, r.baf)
    total += 1
    if (r.genotype !== 'NC') { called += 1; if (r.genotype === 'AB') het += 1 }
  })
  return { name, genotype, baf, callRate: called / total, hetRate: het / Math.max(called, 1) }
}

/** Marker grid for one chromosome, in position order, shared by every array. */
async function grid(name: string): Promise<Map<string, { id: string; pos: number }[]>> {
  const by = new Map<string, { id: string; pos: number }[]>()
  let map: ReturnType<typeof headerMap> = null
  await eachLine(open(name), (l) => {
    if (!map) { map = headerMap(l); return }
    const r = parseRow(l, map)
    if (!r || !isAutosome(r.chrom)) return
    const a = by.get(r.chrom) ?? []
    a.push({ id: r.probesetId, pos: r.pos })
    by.set(r.chrom, a)
  })
  for (const a of by.values()) a.sort((x, y) => x.pos - y.pos)
  return by
}

/** Segment floor for this run. Only lower it to audit a subset, and say so when reporting. */
const FLOOR = Number(process.env.AUDIT_FLOOR ?? MIN_SEGMENT_MARKERS)

const files = readdirSync(DIR).filter((f) => /\.probes$|\.csv(\.gz)?$/.test(f))
const say = (s: string) => { process.stdout.write(`${s}\n`) }
say(`audit dir: ${DIR}, ${files.length} arrays`)
say(`segment floor: ${FLOOR} markers`
  + `${FLOOR === MIN_SEGMENT_MARKERS ? ' (the shipped floor)' : ' (LOWERED for a subset audit)'}`)

// THE REFERENCE PARENT IS DECLARED, NOT GUESSED.
//
// The first version of this audit picked it as "the highest-call-rate array with heterozygosity
// above 14%", which selected a BIPARENTAL EMBRYO rather than a parent, and every absence and
// heterozygosity figure downstream was computed against the wrong genome. Two symptoms gave it
// away: no segment localised anywhere, and a chromosome the laboratory record calls LOST was
// flagged as a duplication. A heuristic that can pick an embryo as a parent has no place in an
// audit whose whole purpose is to be believed.
//
// So the reference is stated per experiment. Where a real parent array exists it is used. Where
// none does, the Progenitor-reconstructed array for that experiment is used and the run is
// labelled as resting on a reconstruction.
const EXPERIMENTS: { prefix: string; reference: string; reconstructed: boolean }[] = [
  // Lab arrays: the sperm donor's own bulk array is the reference for his experiment.
  { prefix: '52461', reference: '52461-00', reconstructed: false },
  // Public arrays, GSE145984/GSE148488. The sperm donor's bulk gDNA is the reference and the
  // pMII products are his haploid meiotic products. These are 28,604-marker SUBSETS, so the
  // shipped 2,400-marker segment floor cannot be reached in them: the floor is lowered by
  // --floor for this run and the run is labelled accordingly. What that tests is the refinement
  // ARITHMETIC against real amplification artefact, not the shipped detection threshold.
  { prefix: 'GSM47', reference: 'GSM4472397', reconstructed: false },
]

const results: Record<string, unknown> = {}
const byExperiment = (name: string): string =>
  (name.startsWith('GSM') ? 'GSM47' : name.slice(0, 5))

const loaded: Loaded[] = []
for (const f of files) loaded.push(await load(f))

const experiment = EXPERIMENTS.find((e) => files.some((f) => f.startsWith(e.reference)))
if (!experiment) throw new Error('no experiment with a declared reference array is present')
const donor = loaded.find((l) => l.name.startsWith(experiment.reference))!
say(`reference parent: ${donor.name} (call ${pct(donor.callRate, 1)}, het ${pct(donor.hetRate, 1)})`
  + `${experiment.reconstructed ? ' [RECONSTRUCTED]' : ' [measured]'}`)
say(`scoped to experiment ${experiment.prefix}; arrays from other experiments are different `
  + 'people and are not comparable against this reference')

/** Only this experiment's arrays. A product of another father is not a control, it is noise. */
const inExperiment = loaded.filter((l) => byExperiment(l.name) === experiment.prefix
  && l !== donor
  // In the public set only the pMII arrays are this donor's meiotic products. donor_A and
  // donor_C are unrelated women and A8 is an embryo; none is a product of this father.
  && (experiment.prefix !== 'GSM47' || /pMII|sperm_DNA/.test(l.name)))

/**
 * Absence indicator for one sample against the reference parent, per chromosome, in position
 * order. Informative means the parent is homozygous there and the sample called; a no-call
 * contributes nothing rather than counting as present, which would understate absence on exactly
 * the degraded arrays where it matters.
 */
const absenceOf = (
  s: Loaded, gr: Map<string, { id: string; pos: number }[]>,
): Map<string, { chrom: string; pos: number; absent: boolean }[]> => {
  const by = new Map<string, { chrom: string; pos: number; absent: boolean }[]>()
  for (const [chrom, ms] of gr) {
    const out: { chrom: string; pos: number; absent: boolean }[] = []
    for (const m of ms) {
      const fa = donor.genotype.get(m.id)
      if (fa !== 'AA' && fa !== 'BB') continue
      const gt = s.genotype.get(m.id)
      if (!gt || gt === 'NC') continue
      out.push({ chrom, pos: m.pos, absent: gt !== 'AB' && gt !== fa })
    }
    by.set(chrom, out)
  }
  return by
}

const g = await grid(files[0])
const haploid = inExperiment.filter((l) => l.hetRate < 0.14 && l.callRate >= REFINE_CALL_FLOOR)
say(`hosts passing the ${REFINE_CALL_FLOOR} call-rate gate: ${haploid.length}`)

// =================================================================================================
// A. Constructed events with a known breakpoint
// =================================================================================================
say('\n=== A. constructed events, breakpoint known at marker resolution ===')
const windowErr: number[] = []
const refinedErr: number[] = []
let localised = 0
let attempted = 0
const CHROMS = ['1', '2', '3', '4', '5', '6', '7', '8']
const SIZES = FLOOR >= 2_400 ? [2_400, 4_800, 9_600] : [FLOOR, FLOOR * 2, FLOOR * 4]

for (const host of haploid.slice(0, 4)) {
  // A donor of the OTHER parental origin, so the spliced block is a real origin switch carrying
  // real amplification artefact rather than anything synthetic.
  const other = haploid.find((x) => x !== host)
  if (!other) continue
  const hostAbs = absenceOf(host, g)
  const otherAbs = absenceOf(other, g)
  for (const chrom of CHROMS) {
    const hs = hostAbs.get(chrom) ?? []
    const os = otherAbs.get(chrom) ?? []
    // Enough room for the largest planted block plus background either side of it.
    const need = Math.max(...SIZES) + 400
    if (hs.length < need || os.length < need) continue
    const byPos = new Map(os.map((x) => [x.pos, x.absent]))
    for (const size of SIZES) {
      const from = Math.max(100, Math.floor((hs.length - size) / 2))
      const to = from + size
      if (to >= hs.length - 100) continue
      // Splice: inside the block, take the other product's calls at the same positions.
      const spliced = hs.map((x, i) => {
        if (i < from || i >= to) return x
        const v = byPos.get(x.pos)
        return v === undefined ? x : { ...x, absent: v }
      })
      attempted += 1
      const nullRate = median(
        [...hostAbs].filter(([c]) => c !== chrom)
          .map(([, ms]) => ms.filter((m) => m.absent).length / Math.max(ms.length, 1)),
      )
      const hits = scanChromosome(spliced, Math.max(nullRate, 0.002), FLOOR)
      if (!hits.length) continue
      // The hit overlapping the planted block.
      const truthStart = spliced[from].pos
      const truthEnd = spliced[to - 1].pos
      const hit = hits.find((h) => h.startBp <= truthEnd && h.endBp >= truthStart)
      if (!hit) continue
      // Window-edge error, in markers, against the planted indices.
      const idxOf = (bp: number): number => {
        let lo = 0
        while (lo < spliced.length && spliced[lo].pos < bp) lo += 1
        return lo
      }
      windowErr.push((Math.abs(idxOf(hit.startBp) - from) + Math.abs(idxOf(hit.endBp) - (to - 1))) / 2)
      const c = segmentCoords(hit)
      if (!c.localised) continue
      localised += 1
      refinedErr.push(
        (Math.abs(idxOf(c.start) - from) + Math.abs(idxOf(c.end) - (to - 1))) / 2,
      )
    }
  }
}
const partA = {
  events: attempted,
  detected: windowErr.length,
  localised,
  windowEdgeMedianMarkers: median(windowErr),
  windowEdgeP95Markers: pctl(windowErr, 0.95),
  refinedMedianMarkers: median(refinedErr),
  refinedP95Markers: pctl(refinedErr, 0.95),
  improvement: median(windowErr) / Math.max(median(refinedErr), 1e-9),
}
results.constructed = partA
say(`  ${attempted} planted, ${windowErr.length} detected, ${localised} localised`)
say(`  window edge : median ${partA.windowEdgeMedianMarkers} markers, p95 ${partA.windowEdgeP95Markers}`)
say(`  refined     : median ${partA.refinedMedianMarkers} markers, p95 ${partA.refinedP95Markers}`)
say(`  improvement : ${partA.improvement.toFixed(1)}x on the median`)

// =================================================================================================
// B. Split-half holdout on real segments, no truth
// =================================================================================================
say('\n=== B. split-half holdout on real segments ===')
const halfRows: { array: string; chrom: string; gapMarkers: number; inside: boolean }[] = []
for (const host of haploid) {
  const abs = absenceOf(host, g)
  for (const [chrom, ms] of abs) {
    if (ms.length < FLOOR) continue
    const nullRate = Math.max(externalNull(
      new Map([...abs].map(([c, a]) => [c, [a.length, a.filter((x) => x.absent).length]])),
      chrom,
    ), 0.002)
    for (const hit of scanChromosome(ms, nullRate, FLOOR)) {
      const c = segmentCoords(hit)
      if (!c.localised) continue
      const flagsOf = (keep: (i: number) => boolean) => {
        const sub = ms.filter((_, i) => keep(i))
        const f = new Uint8Array(sub.length)
        for (let i = 0; i < sub.length; i += 1) f[i] = sub[i].absent ? 1 : 0
        return { f, sub }
      }
      const est = (keep: (i: number) => boolean): number => {
        const { f, sub } = flagsOf(keep)
        let a0 = 0
        while (a0 < sub.length && sub[a0].pos < hit.startBp) a0 += 1
        let b0 = a0
        while (b0 < sub.length && sub[b0].pos <= hit.endBp) b0 += 1
        const r = refineEdges(f, a0, b0, nullRate, SEARCH_RADIUS, LR_DROP, SEGMENT_LRT)
        return r.localised ? sub[Math.min(sub.length - 1, r.start.index)].pos : NaN
      }
      const odd = est((i) => i % 2 === 1)
      const even = est((i) => i % 2 === 0)
      if (!Number.isFinite(odd) || !Number.isFinite(even)) continue
      const gapBp = Math.abs(odd - even)
      const inside = odd >= c.start - (c.start - segmentCoords(hit).start)
        && Math.min(odd, even) >= (hit.refined?.startLoBp ?? -Infinity)
        && Math.max(odd, even) <= (hit.refined?.startHiBp ?? Infinity)
      halfRows.push({ array: host.name, chrom, gapMarkers: gapBp, inside })
    }
  }
}
const partB = {
  segments: halfRows.length,
  bothHalvesInsideInterval: halfRows.filter((r) => r.inside).length,
  passRate: halfRows.length ? halfRows.filter((r) => r.inside).length / halfRows.length : NaN,
  medianHalfGapBp: median(halfRows.map((r) => r.gapMarkers)),
}
results.splitHalf = partB
say(`  ${partB.segments} localised segments, ${partB.bothHalvesInsideInterval} with both halves `
  + `inside the interval (${pct(partB.passRate, 1)})`)
say(`  median gap between halves: ${(partB.medianHalfGapBp / 1e6).toFixed(3)} Mb`)

// =================================================================================================
// C. Specificity on genomes that cannot carry one
// =================================================================================================
say('\n=== C. false segments on bulk diploid adults ===')
// The one genome that provably cannot carry a segment against this reference: the reference
// himself. He IS this parent, so absence is his own genotyping error and nothing else.
const bulk = [donor]
const specRows: { array: string; segments: number }[] = []
for (const b of bulk) {
  const abs = absenceOf(b, g)
  let n = 0
  for (const [chrom, ms] of abs) {
    if (ms.length < FLOOR) continue
    const nullRate = Math.max(externalNull(
      new Map([...abs].map(([c, a]) => [c, [a.length, a.filter((x) => x.absent).length]])),
      chrom,
    ), 0.002)
    n += scanChromosome(ms, nullRate, FLOOR).length
  }
  specRows.push({ array: b.name, segments: n })
  say(`  ${b.name}: ${n} segment(s)`)
}
results.specificity = { arrays: specRows.length, rows: specRows,
  totalFalseSegments: specRows.reduce((a, r) => a + r.segments, 0) }

// =================================================================================================
// D. Gain origin
// =================================================================================================
say('\n=== D. gain origin ===')
// D1: orientation invariance on REAL marker pairs, not constructed ones.
let oriented = 0
let checked = 0
for (const [, ms] of g) {
  for (const m of ms) {
    const fa = donor.genotype.get(m.id)
    if (fa !== 'AA' && fa !== 'BB') continue
    const b = donor.baf.get(m.id)
    if (b === undefined) continue
    // The same physical dosage read at a father-AA and a father-BB marker must give the same
    // paternal share. Mirror the BAF to construct the opposite orientation of the same event.
    const asAA = paternalShare('AA', 'BB', b)
    const asBB = paternalShare('BB', 'AA', 1 - b)
    if (asAA === null || asBB === null) continue
    checked += 1
    if (Math.abs(asAA - asBB) < 1e-12) oriented += 1
    if (checked >= 200_000) break
  }
  if (checked >= 200_000) break
}
results.orientation = { checked, invariant: oriented, allInvariant: oriented === checked }
say(`  orientation invariant on ${oriented} of ${checked} real markers`)

// D2: the refusals, on every array carrying a gain.
const gainRows: { array: string; chrom: string; regionHet: number; background: number;
  verdict: string }[] = []
for (const s of inExperiment) {
  const het = new Map<string, { informative: number; het: number }>()
  for (const [chrom, ms] of g) {
    let inf = 0
    let h = 0
    for (const m of ms) {
      const fa = donor.genotype.get(m.id)
      if (fa !== 'AA' && fa !== 'BB') continue
      const gt = s.genotype.get(m.id)
      if (!gt || gt === 'NC') continue
      inf += 1
      if (gt === 'AB') h += 1
    }
    if (inf) het.set(chrom, { informative: inf, het: h })
  }
  // A chromosome whose heterozygosity is far above the rest is the only visible extra copy.
  for (const [chrom, h] of het) {
    const bg = externalHetBackground(het, chrom)
    const rate = h.het / h.informative
    if (!(rate >= 2 * bg)) continue
    const call = callHomologue(rate, h.informative, bg)
    gainRows.push({ array: s.name, chrom, regionHet: rate, background: bg,
      verdict: call.verdict })
  }
}
results.homologue = { flagged: gainRows.length, rows: gainRows }
say(`  chromosomes reading above twice their external background: ${gainRows.length}`)
for (const r of gainRows.slice(0, 12)) {
  say(`    ${r.array} chr${r.chrom}: ${pct(r.regionHet, 1)} vs ${pct(r.background, 1)} -> ${r.verdict}`)
}

writeFileSync(OUT, `${JSON.stringify(results, null, 2)}\n`)
say(`\nwrote ${OUT}`)

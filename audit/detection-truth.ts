/**
 * Whether chromosomal change detection finds what is there and nothing else, against truth the
 * tool did not produce. Detection only: the call-rate and intensity channels never read the
 * parental array, so whether a change exists, and where, is independent of attribution.
 *
 * Arm 1, specificity: the adult gamete donors' bulk genomic DNA. A living adult who produced
 * viable gametes carries no autosomal whole-chromosome loss, so every whole-chromosome autosomal
 * call on those arrays is a false positive, with no injection, construction or noise model.
 *
 * Arm 2, sensitivity with a position. Zuccaro et al. cut EYS, on chromosome 6, with Cas9, so
 * chromosome 6 is where a change is expected and the other 21 autosomes are a specificity arm on
 * the same arrays: a detector firing at random shows the same per-chromosome rate on 6 as on 13.
 * The locus is used for position only, never as a filter. Segments are reported wherever they are
 * found and their distance from EYS measured afterwards.
 *
 * Arm 3, structure. Every emitted event is checked for a start past its own end, coordinates off
 * the end of its chromosome, a non-finite number in a field a reader will read, zero markers behind
 * a call, the same locus emitted twice, and a per-array event count no genome could carry. Any such
 * event, or any array that throws, fails the run.
 *
 * Run: OM_TRIOS=<dir> [OM_LAB=<dir>] node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/detection-truth.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
const LAB = process.env.OM_LAB
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('detection-truth: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}
const CONTAMINATED = 'GSM4472408'

/**
 * EYS, the gene Zuccaro et al. cut, on chromosome 6.
 *
 * Both assemblies are carried because the reader infers the build from marker positions and this
 * series is not guaranteed to be one of them. The span is only ever used to measure how far a
 * reported segment sits from the cut, never to decide whether to report it.
 */
const EYS = {
  chrom: '6',
  hg19: { start: 64_429_876, end: 66_417_118 },
  hg38: { start: 63_719_980, end: 65_707_213 },
}

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')], pronucleus: f[col('pronucleus')],
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

function rowsOfPath(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  if (h < 0) throw new Error('no header')
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

/**
 * A FIXED ARRAY IN THE PARENTAL SLOT. Both channels under test here are parent-independent, so
 * using one array for every sample means the parental slot cannot be a variable in the result.
 */
const REF = join(DIR, 'GSM4472397.probes.gz')
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOfPath(REF)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Ev {
  kind: 'whole' | 'segment'
  chrom: string
  cls: string
  startBp: number
  endBp: number
  markers: number
}
interface Scored {
  name: string
  group: string
  stage: string
  callRate: number
  build: string
  events: Ev[]
  /** Highest position the ARRAY ITSELF carries per chromosome, so a coordinate can be checked
   *  against what was measured rather than against an external table. */
  maxPos: Map<string, number>
  junk: string[]
  error?: string
}

const isAutosome = (c: string) => /^\d+$/.test(c)

async function scoreOne(path: string, name: string, group: string): Promise<Scored> {
  const maxPos = new Map<string, number>()
  const junk: string[] = []
  try {
    const acc = score.emptyCollected(refParent as never, null)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOfPath(path)()) {
      if (!first) first = r.probesetId
      if (Number.isFinite(r.pos)) {
        maxPos.set(r.chrom, Math.max(maxPos.get(r.chrom) ?? 0, r.pos))
      }
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(name, byChrom as never, bafSums as never, first)
    const out = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: name, log: () => {},
    }) as Record<string, unknown>

    // `informative`, not `markers`: a whole-chromosome row counts informative markers and has no
    // field called `markers`. Reading the wrong name defaulted every one of them to zero and
    // reported 35 perfectly good calls as structurally invalid.
    const chroms = (out.chroms ?? []) as {
      chrom: string; aneuploidy?: string; informative?: number
    }[]
    const segs = (out.segments ?? []) as {
      chrom: string; kind: string; markers: number; startBp: number; endBp: number
      refined?: { startBp: number; endBp: number }
    }[]
    const events: Ev[] = [
      ...chroms.filter((c) => c.aneuploidy).map((c) => ({
        kind: 'whole' as const, chrom: c.chrom, cls: String(c.aneuploidy),
        startBp: 0, endBp: maxPos.get(c.chrom) ?? 0, markers: c.informative ?? 0,
      })),
      ...segs.map((s) => ({
        kind: 'segment' as const, chrom: s.chrom, cls: s.kind,
        startBp: s.refined?.startBp ?? s.startBp,
        endBp: s.refined?.endBp ?? s.endBp,
        markers: s.markers,
      })),
    ]

    // ---------------------------------------------------------------- ARM 3, the junk checks
    const seen = new Set<string>()
    for (const e of events) {
      const at = `${e.kind} chr${e.chrom} ${e.cls}`
      if (!Number.isFinite(e.startBp) || !Number.isFinite(e.endBp)) {
        junk.push(`${at}: a coordinate is not a finite number`)
      } else if (e.kind === 'segment' && !(e.startBp < e.endBp)) {
        junk.push(`${at}: start ${e.startBp} is not before end ${e.endBp}`)
      }
      const top = maxPos.get(e.chrom)
      if (top === undefined) {
        junk.push(`${at}: reported on a chromosome this array carries no markers for`)
      } else if (Number.isFinite(e.endBp) && e.endBp > top * 1.05) {
        junk.push(`${at}: ends at ${e.endBp} past the last marker on it, ${top}`)
      }
      if (!Number.isFinite(e.markers) || e.markers <= 0) {
        junk.push(`${at}: ${e.markers} markers behind the call`)
      }
      const key = `${e.kind}:${e.chrom}:${e.startBp}:${e.endBp}`
      if (seen.has(key)) junk.push(`${at}: the same locus is emitted twice`)
      seen.add(key)
    }
    // A GENOME CANNOT LOSE MOST OF ITSELF AND STILL BE READ AS A GENOME. Not a threshold on the
    // biology: a call set this large on an array the tool did not refuse means the refusal, not
    // the detector, is what failed.
    const wholeAuto = events.filter((e) => e.kind === 'whole' && isAutosome(e.chrom)).length
    const stage = (out.stage as { stage?: string } | undefined)?.stage ?? 'none'
    if (wholeAuto > 15 && stage !== 'failed') {
      junk.push(`${wholeAuto} of 22 autosomes called whole-chromosome events on an array the tool `
        + `did not refuse (stage ${stage})`)
    }

    return {
      name, group, stage, callRate: profile.callRate,
      build: String((profile as { build?: { build?: string } }).build?.build ?? '?'),
      events, maxPos, junk,
    }
  } catch (e) {
    return {
      name, group, stage: 'error', callRate: NaN, build: '?', events: [], maxPos, junk,
      error: String((e as Error).message ?? e).slice(0, 90),
    }
  }
}

const results: Scored[] = []

// ================================================================ ARM 1: KNOWN-NORMAL BULK ADULTS
const donors = recs.filter((r) => r.role === 'parent' && has(r.gsm) && r.gsm !== CONTAMINATED)
console.log(`ARM 1, specificity: ${donors.length} bulk adult donor arrays`)
for (const d of donors) results.push(await scoreOne(join(DIR, `${d.gsm}.probes.gz`), d.gsm, 'donor'))

// ================================================================ ARM 2: THE EDITED EMBRYOS
const CAP = Number(process.env.OM_CAP ?? 45)
const embryos = recs.filter((r) => r.role !== 'parent' && !r.pronucleus && has(r.gsm)).slice(0, CAP)
console.log(`ARM 2, on-target: ${embryos.length} embryo arrays from the EYS editing experiment`)
for (const e of embryos) {
  results.push(await scoreOne(join(DIR, `${e.gsm}.probes.gz`), e.gsm, 'embryo'))
}

// ================================================================ ARM 3 also covers the lab corpus
if (LAB && existsSync(LAB)) {
  const walk = (dir: string, out: string[] = []): string[] => {
    for (const e of readdirSync(dir)) {
      const p = join(dir, e)
      let st
      try { st = statSync(p) } catch { continue }
      if (st.isDirectory()) walk(p, out)
      else if (e.endsWith('.probes') || e.endsWith('.probes.gz')) out.push(p)
    }
    return out
  }
  const labFiles = walk(LAB).slice(0, Number(process.env.OM_LAB_CAP ?? 60))
  console.log(`ARM 3 also over ${labFiles.length} laboratory arrays`)
  for (const p of labFiles) results.push(await scoreOne(p, p.slice(LAB.length + 1), 'lab'))
}
console.log('')

// ---------------------------------------------------------------- ARM 1 REPORT
console.log('=== ARM 1: OVERDETECTION ON KNOWN-NORMAL BULK ADULT DNA')
console.log('A living adult gamete donor carries no autosomal whole-chromosome loss, so every call')
console.log('below is a false positive. Nothing is constructed and no noise model is assumed.')
const donorRes = results.filter((r) => r.group === 'donor')
let fpWhole = 0
let fpSeg = 0
for (const r of donorRes) {
  const whole = r.events.filter((e) => e.kind === 'whole' && isAutosome(e.chrom))
  const seg = r.events.filter((e) => e.kind === 'segment' && isAutosome(e.chrom))
  fpWhole += whole.length
  fpSeg += seg.length
  console.log(`  ${r.name}  stage ${r.stage.padEnd(14)} call ${r.callRate.toFixed(3)}  `
    + `whole-chromosome ${whole.length}  segments ${seg.length}`
    + `${whole.length ? `  <-- ${whole.map((e) => `chr${e.chrom}:${e.cls}`).join(' ')}` : ''}`)
}
const donorChromObs = donorRes.length * 22
console.log(`  TOTAL over ${donorRes.length} arrays and ${donorChromObs} autosome observations: `
  + `${fpWhole} whole-chromosome false positives, ${fpSeg} segments`)
console.log(`  per-autosome false positive rate: ${(fpWhole / Math.max(donorChromObs, 1)).toFixed(5)}`)

// ---------------------------------------------------------------- ARM 2 REPORT
console.log('')
console.log('=== ARM 2: ON TARGET. Cas9 cut EYS, which is on chromosome 6')
const emb = results.filter((r) => r.group === 'embryo' && !r.error && r.stage !== 'failed')
const withChrom = (c: string) => emb.filter((r) => r.events.some((e) => e.chrom === c)).length
const chr6 = withChrom('6')
const others = Array.from({ length: 22 }, (_, i) => String(i + 1)).filter((c) => c !== '6')
const otherCounts = others.map(withChrom)
const otherMean = otherCounts.reduce((a, b) => a + b, 0) / otherCounts.length
console.log(`  usable embryo arrays: ${emb.length}`)
console.log(`  carrying an event on chr6 (the cut chromosome): ${chr6}/${emb.length} `
  + `= ${(chr6 / Math.max(emb.length, 1)).toFixed(3)}`)
console.log(`  mean over the other 21 autosomes:               ${otherMean.toFixed(2)}/${emb.length} `
  + `= ${(otherMean / Math.max(emb.length, 1)).toFixed(3)}`)
console.log(`  ENRICHMENT: ${(chr6 / Math.max(otherMean, 1e-9)).toFixed(1)}x`)
const ranked = [['6', chr6] as [string, number], ...others.map((c, i) => [c, otherCounts[i]] as [string, number])]
  .sort((a, b) => b[1] - a[1])
console.log(`  chromosomes by how many arrays carry an event, highest first:`)
console.log(`    ${ranked.slice(0, 6).map(([c, n]) => `chr${c}=${n}`).join('  ')}`)
console.log(`  chr6 rank: ${ranked.findIndex(([c]) => c === '6') + 1} of 22`)

// Position, which is a stronger claim than the chromosome.
// BOTH ASSEMBLIES, rather than picking one from a build call that can come back unresolved. EYS
// sits 0.7 Mb apart between them, so reporting both removes the dependency instead of hiding it.
console.log('')
console.log('  POSITION, against EYS. Both assemblies are shown because the span differs by 0.7 Mb')
console.log('  between them and nothing here should rest on which one was inferred.')
const chr6segs = emb.flatMap((r) => r.events
  .filter((e) => e.chrom === '6' && e.kind === 'segment')
  .map((e) => ({ name: r.name, e })))
for (const [name, span] of [['GRCh37', EYS.hg19], ['GRCh38', EYS.hg38]] as const) {
  const covers = chr6segs.filter(({ e }) => e.startBp <= span.end && e.endBp >= span.start)
  console.log(`    ${name} chr6:${(span.start / 1e6).toFixed(1)}-${(span.end / 1e6).toFixed(1)}Mb  `
    + `chr6 segments ${chr6segs.length}, spanning the cut site ${covers.length}`)
}
const chr6whole = emb.filter((r) => r.events.some((e) => e.chrom === '6' && e.kind === 'whole')).length
console.log(`    chr6 WHOLE-CHROMOSOME events: ${chr6whole} arrays (a whole chromosome contains the`)
console.log('    cut site by definition, so it is counted separately rather than as a hit)')

// ---------------------------------------------------------------- ARM 3 REPORT
console.log('')
console.log('=== ARM 3: JUNK. Structural validity of every event emitted anywhere above')
const allJunk = results.flatMap((r) => r.junk.map((j) => `${r.group}/${r.name}: ${j}`))
const errored = results.filter((r) => r.error)
const totalEvents = results.reduce((a, r) => a + r.events.length, 0)
console.log(`  arrays scored: ${results.length}, events emitted: ${totalEvents}`)
console.log(`  arrays that threw: ${errored.length}`)
for (const e of errored.slice(0, 10)) console.log(`    ${e.name}: ${e.error}`)
console.log(`  structurally invalid events: ${allJunk.length}`)
for (const j of allJunk.slice(0, 25)) console.log(`    ${j}`)
if (allJunk.length > 25) console.log(`    ... and ${allJunk.length - 25} more`)

console.log('')
const clean = allJunk.length === 0 && errored.length === 0
console.log(clean
  ? 'detection-truth: every emitted event is structurally valid and no array threw.'
  : `detection-truth: ${allJunk.length} invalid events and ${errored.length} throws. `
    + 'Either is a failure whatever the statistics say.')
process.exit(clean ? 0 : 1)

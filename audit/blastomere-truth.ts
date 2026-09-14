/**
 * The pipeline on single blastomeres, where agreement within one embryo is and is not required.
 *
 * GSE186407 and GSE290961 are single blastomeres from human embryos injected at the 2PN stage,
 * arrayed one cell at a time, several cells per embryo.
 *
 * Truth, from the filename rather than this tool:
 *   the target locus     GSE186407 cut chr16, and every file says so on its face
 *   the embryo           blastomeres carrying the same embryo tag are cells of ONE genome.
 *                        GSE290961 carries no embryo identifier, so it contributes only the
 *                        structural arm.
 *
 * Properties of the PERSON must agree across one embryo's cells: which parental genome is present,
 * and the sex unless a called sex-chromosome event explains the difference. Properties of a CELL
 * need not: these embryos are frequently mosaic (Zuccaro and others), so copy number is reported
 * as a spread with no pass mark.
 *
 * Failure: an array that throws, a structurally invalid event, or a sex difference within an
 * embryo with no called sex-chromosome event behind it.
 *
 * Run: OM_BLASTO=<dir> [OM_TARGET=16] node --experimental-strip-types \
 *        --max-old-space-size=3072 audit/blastomere-truth.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const sexing = await import(`${W}sexing.ts`)

const DIR = process.env.OM_BLASTO
if (!DIR || !existsSync(DIR)) {
  console.log('blastomere-truth: OM_BLASTO is not a readable directory. NOT RUN.')
  process.exit(0)
}
/** The chromosome the experiment cut, when the series states one. */
const TARGET = process.env.OM_TARGET || ''

function rowsOf(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

/** The embryo tag a file carries, lower-cased because the series mixes Z10 and z10. */
function embryoOf(name: string): string | null {
  const m = /_([zZ]\d+)_/.exec(name) || /_(\d{5}-\d+)_/.exec(name)
  return m ? m[1].toLowerCase() : null
}

const files = readdirSync(DIR).filter((f) => f.endsWith('.probes.txt.gz') || f.endsWith('.probes'))
  .sort()
const CAP = Number(process.env.OM_CAP ?? 400)
const use = files.slice(0, CAP)
console.log(`blastomere arrays: ${use.length} of ${files.length}`)
if (TARGET) console.log(`the experiment cut chromosome ${TARGET}; the other 21 are the control`)

/** A fixed array in the parental slot. The channels scored here never consult it. */
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR, use[0]))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Row {
  name: string; embryo: string | null; stage: string; callRate: number
  zygosity: string; y: boolean | null; eventChroms: string[]; onTarget: boolean
  sexEvent: boolean
  junk: string[]; error?: string
}
const out: Row[] = []
const isAutosome = (c: string) => /^\d+$/.test(c)

for (const f of use) {
  const path = join(DIR, f)
  try {
    const acc = score.emptyCollected(refParent as never, null)
    const sx = sexing.emptySex()
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    const maxPos = new Map<string, number>()
    for (const r of rowsOf(path)()) {
      if (!first) first = r.probesetId
      if (Number.isFinite(r.pos)) maxPos.set(r.chrom, Math.max(maxPos.get(r.chrom) ?? 0, r.pos))
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      sexing.accumulateSex(r as never, sx as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(f, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: f, log: () => {},
    }) as Record<string, unknown>
    const chroms = (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
    const segs = (res.segments ?? []) as {
      chrom: string; markers: number; startBp: number; endBp: number
      refined?: { startBp: number; endBp: number }
    }[]
    const junk: string[] = []
    for (const s of segs) {
      const st = s.refined?.startBp ?? s.startBp
      const en = s.refined?.endBp ?? s.endBp
      const top = maxPos.get(s.chrom)
      if (!Number.isFinite(st) || !Number.isFinite(en)) junk.push(`chr${s.chrom} NaN coordinate`)
      else if (!(st < en)) junk.push(`chr${s.chrom} start not before end`)
      if (top !== undefined && Number.isFinite(en) && en > top * 1.05) {
        junk.push(`chr${s.chrom} ends past its last marker`)
      }
      if (!Number.isFinite(s.markers) || s.markers <= 0) junk.push(`chr${s.chrom} zero markers`)
    }
    const eventChroms = [...new Set([
      ...chroms.filter((c) => c.aneuploidy && isAutosome(c.chrom)).map((c) => c.chrom),
      ...segs.filter((s) => isAutosome(s.chrom)).map((s) => s.chrom),
    ])]
    out.push({
      name: f, embryo: embryoOf(f),
      stage: (res.stage as { stage?: string } | undefined)?.stage ?? 'none',
      callRate: profile.callRate,
      zygosity: String(res.zygosity ?? '?'),
      y: (sexing.sexCall(sx as never) as { yBearing: boolean | null }).yBearing,
      eventChroms,
      onTarget: TARGET ? eventChroms.includes(TARGET) : false,
      sexEvent: chroms.some((c) => c.aneuploidy && !isAutosome(c.chrom))
        || segs.some((sg) => !isAutosome(sg.chrom)),
      junk,
    })
  } catch (e) {
    out.push({
      name: f, embryo: embryoOf(f), stage: 'error', callRate: NaN, zygosity: '?',
      y: null, eventChroms: [], onTarget: false, sexEvent: false, junk: [],
      error: String((e as Error).message ?? e).slice(0, 80),
    })
  }
  if (out.length % 25 === 0) console.log(`  scored ${out.length}/${use.length}`)
}

const ok = out.filter((r) => !r.error)
const usable = ok.filter((r) => r.stage !== 'failed')
console.log('')
console.log(`scored ${out.length}, threw ${out.length - ok.length}, `
  + `resolved to a stage ${usable.length}`)

// ---------------------------------------------------------------- on target
if (TARGET) {
  console.log('')
  console.log(`=== ON TARGET: the experiment cut chromosome ${TARGET}`)
  const withT = usable.filter((r) => r.eventChroms.includes(TARGET)).length
  const others = Array.from({ length: 22 }, (_, i) => String(i + 1)).filter((c) => c !== TARGET)
  const otherCounts = others.map((c) => usable.filter((r) => r.eventChroms.includes(c)).length)
  const mean = otherCounts.reduce((a, b) => a + b, 0) / otherCounts.length
  console.log(`  arrays with an event on chr${TARGET}: ${withT}/${usable.length} `
    + `= ${(withT / Math.max(usable.length, 1)).toFixed(3)}`)
  console.log(`  mean over the other 21 autosomes:    ${mean.toFixed(2)}/${usable.length} `
    + `= ${(mean / Math.max(usable.length, 1)).toFixed(3)}`)
  console.log(`  ENRICHMENT: ${(withT / Math.max(mean, 1e-9)).toFixed(1)}x`)
  const ranked = [[TARGET, withT] as [string, number],
    ...others.map((c, i) => [c, otherCounts[i]] as [string, number])]
    .sort((a, b) => b[1] - a[1])
  console.log(`  ranked: ${ranked.slice(0, 6).map(([c, n]) => `chr${c}=${n}`).join('  ')}`)
  console.log(`  chr${TARGET} rank: ${ranked.findIndex(([c]) => c === TARGET) + 1} of 22`)
}

// ---------------------------------------------------------------- 1:1 within one embryo
console.log('')
console.log('=== WITHIN ONE EMBRYO. Cells of one person, so identity must be 1:1.')
const groups = new Map<string, Row[]>()
for (const r of usable) if (r.embryo) groups.set(r.embryo, [...(groups.get(r.embryo) ?? []), r])
const multi = [...groups.entries()].filter(([, v]) => v.length >= 2)
console.log(`  embryos with 2 or more usable blastomeres: ${multi.length}`)
if (!multi.length) {
  // NOT A PASS. GSE290961 states only "tissue: Preimplantation embryo" per sample: no embryo
  // identifier anywhere in its metadata, and its file prefixes are case numbers carrying up to 57
  // cells, far more than one day-3 embryo. Printing 0/0 agree here read as a clean sweep of a test
  // that never ran. This series can only contribute the structural arm below.
  console.log('  NOT ASSESSED: this series states no embryo identifier, so cells of one genome')
  console.log('  cannot be identified and nothing below is scored. The structural arm still is.')
}

let sexAgree = 0; let sexSplit = 0; let sexExplained = 0
let zygAgree = 0; let zygSplit = 0
const cnSpread: number[] = []
for (const [emb, cells] of multi) {
  // A SEX DISAGREEMENT IS NOT AUTOMATICALLY A DEFECT, and treating it as one was wrong in the
  // first version of this file. Loss of the Y chromosome is a copy-number event like any other,
  // and these embryos are mosaic: on embryo z10 one blastomere called 0 of 812 Y probes while its
  // six siblings called about 90 percent, which is a real Y loss correctly reported. A split only
  // counts against the tool when nothing on the sex chromosomes explains it.
  const sexes = new Set(cells.map((c) => c.y).filter((v) => v !== null))
  const zygs = new Set(cells.map((c) => c.zygosity).filter((z) => z !== 'unknown'))
  const explained = cells.some((c) => c.sexEvent)
  if (sexes.size <= 1) sexAgree += 1
  else if (explained) sexExplained += 1
  else sexSplit += 1
  if (zygs.size <= 1) zygAgree += 1; else zygSplit += 1
  const sets = cells.map((c) => c.eventChroms.slice().sort().join(','))
  cnSpread.push(new Set(sets).size)
  if (sexes.size > 1) {
    console.log(`    ${emb}: sex differs across ${cells.length} cells `
      + `-> ${cells.map((c) => (c.y === null ? '?' : c.y ? 'Y' : '-')).join(' ')}`
      + `${explained ? '   explained by a sex-chromosome event' : '   UNEXPLAINED, a defect'}`)
  }
}
console.log('')
console.log(`  sex, one value per embryo      : ${sexAgree}/${multi.length} agree, ${sexSplit} split`)
console.log(`  zygosity, one value per embryo : ${zygAgree}/${multi.length} agree, ${zygSplit} split`)
console.log('  Sex is a property of the person and MUST be 1:1. A split is a defect.')
console.log('')
const distinct = cnSpread.reduce((a, b) => a + b, 0) / Math.max(cnSpread.length, 1)
console.log(`  distinct copy-number call sets per embryo, mean: ${distinct.toFixed(2)}`)
console.log('  REPORTED, NOT SCORED. These embryos are mosaic by the source study\'s own finding,')
console.log('  so blastomeres of one embryo legitimately differ here and 1:1 is not the standard.')

// ---------------------------------------------------------------- junk
console.log('')
const allJunk = out.flatMap((r) => r.junk.map((j) => `${r.name}: ${j}`))
console.log(`=== STRUCTURE: ${out.length} arrays, ${out.length - ok.length} threw, `
  + `${allJunk.length} structurally invalid events`)
for (const j of allJunk.slice(0, 10)) console.log(`  ${j}`)
for (const e of out.filter((r) => r.error).slice(0, 6)) console.log(`  THREW ${e.name}: ${e.error}`)
console.log('')
console.log((out.length - ok.length) === 0 && allJunk.length === 0 && sexSplit === 0
  ? 'blastomere-truth: nothing threw, nothing malformed, and every sex difference within an '
    + 'embryo has a called sex-chromosome event behind it.'
  : `blastomere-truth: ${out.length - ok.length} threw, ${allJunk.length} malformed, `
    + `${sexSplit} embryos with a split sex call.`)

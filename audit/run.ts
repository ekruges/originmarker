/**
 * The accuracy audit: every case, end to end, through the code the browser runs.
 *
 * Run: node --experimental-strip-types audit/run.ts <raw-dir> <out-dir>
 *
 * This drives `ingest`, `parentage` and `syngamyPdf` rather than the Python, because those are
 * what a user's browser executes. The only things replaced are the file plumbing (a gzipped file
 * off disk rather than a dropped one) and favicon.svg, read rather than fetched. A parent is held
 * as one call per marker and each sample is streamed against it, exactly as the page does.
 *
 * Ground truth comes from the bench, not from this tool. Zuccaro et al. separated the two
 * pronuclei of a 2PN zygote by micromanipulation and arrayed each alone, so a "paternal nucleus"
 * sample contains one parental genome by construction and a "maternal nucleus" the other.
 */
import { createReadStream, readFileSync, mkdirSync, writeFileSync } from 'node:fs'
import { createGunzip } from 'node:zlib'
import { createInterface } from 'node:readline'
import { createHash } from 'node:crypto'
import { basename, join } from 'node:path'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, gates, headerMap, parseRow,
  type ChromStats, type Gate, type ProbeRow, type SampleProfile,
} from '../web/src/ingest.ts'
import type { AB } from '../web/src/informativity.ts'
import {
  agreement, classify, emptyTally, isAutosome, pair, tallyRow,
  type PairResult, type ParentageResult,
} from '../web/src/parentage.ts'
import { buildReportPdf, reportId, type ReportFile } from '../web/src/syngamyPdf.ts'

const RAW = process.argv[2] ?? 'raw'
const OUT = process.argv[3] ?? 'out'
const MARK_SVG = readFileSync(new URL('../web/public/favicon.svg', import.meta.url), 'utf8')

/** Stream one file's rows, gunzipping as it goes and hashing the uncompressed bytes. */
async function eachRow(path: string, fn: (r: ProbeRow) => void): Promise<{
  profile: SampleProfile; gates: Gate[]; sha256: string; bytes: number; markers: number
}> {
  const byChrom = new Map<string, ChromStats>()
  const baf = emptyBafSums()
  const builds = emptyBuildSums()
  const hash = createHash('sha256')
  let map: ReturnType<typeof headerMap> = null
  let firstId = ''
  let n = 0
  let bytes = 0

  const gz = createReadStream(path).pipe(createGunzip())
  gz.on('data', (c: Buffer) => { hash.update(c); bytes += c.length })
  for await (const line of createInterface({ input: gz, crlfDelay: Infinity })) {
    if (!map) { map = headerMap(line); continue }
    const r = parseRow(line, map)
    if (!r) continue
    if (!firstId) firstId = r.probesetId
    accumulate(r, byChrom)
    accumulateBaf(r, baf)
    accumulateBuild(r, builds)
    fn(r)
    n += 1
  }
  if (!map) throw new Error(`${path}: no recognisable header`)
  const profile: SampleProfile = {
    ...finishProfile(basename(path), byChrom, baf, firstId, builds), build: buildVerdict(builds),
  }
  return { profile, gates: gates(profile), sha256: hash.digest('hex'), bytes, markers: n }
}

export interface Parent {
  file: string
  gt: Map<string, AB>
  heterozygosity: number
  profile: SampleProfile
  gates: Gate[]
  sha256: string
  bytes: number
  markers: number
}

/** A parent, held as one call per marker. */
export async function readParent(path: string): Promise<Parent> {
  const gt = new Map<string, AB>()
  let called = 0
  let het = 0
  const meta = await eachRow(path, (r) => {
    gt.set(r.probesetId, r.genotype)
    if (r.genotype !== 'NC' && isAutosome(r.chrom)) {
      called += 1
      if (r.genotype === 'AB') het += 1
    }
  })
  return { file: basename(path), gt, heterozygosity: called ? het / called : NaN, ...meta }
}

export interface Scored {
  file: string
  profile: SampleProfile
  gates: Gate[]
  result: ParentageResult
  maternal?: ParentageResult
  paired?: PairResult
  sha256: string
  bytes: number
  markers: number
}

/** One sample against one or two parents, streamed and discarded. */
export async function scoreSample(
  path: string, pat: Parent, mat?: Parent, agree = NaN,
): Promise<Scored> {
  const tp = emptyTally()
  const tm = mat ? emptyTally() : null
  const meta = await eachRow(path, (r) => {
    tallyRow(pat.gt.get(r.probesetId) ?? 'NC', r, tp)
    if (tm && mat) tallyRow(mat.gt.get(r.probesetId) ?? 'NC', r, tm)
  })
  const result = classify(tp, pat.heterozygosity)
  const maternal = tm && mat ? classify(tm, mat.heterozygosity, { role: 'maternal' }) : undefined
  return {
    file: basename(path),
    result,
    maternal,
    paired: maternal ? pair(result, maternal, agree) : undefined,
    ...meta,
  }
}

/** A report for one case, written to disk. Same builder the Report (PDF) button calls. */
export async function writePdf(
  name: string, pat: Parent, sample: Scored, mat: Parent | undefined, generatedAt: string,
  tool: string,
): Promise<{ path: string; bytes: number; id: string }> {
  const asFile = (p: Parent, role: 'donor' | 'oocyte'): ReportFile => ({
    name: p.file, size: p.bytes, sha256: p.sha256, role, markers: p.markers,
    profile: p.profile, gates: p.gates,
  })
  const files: ReportFile[] = [asFile(pat, 'donor')]
  if (mat) files.push(asFile(mat, 'oocyte'))
  files.push({
    name: sample.file, size: sample.bytes, sha256: sample.sha256, role: 'sample',
    markers: sample.markers, profile: sample.profile, gates: sample.gates,
    result: sample.result, maternal: sample.maternal, paired: sample.paired,
  })
  const id = reportId(files.map((f) => f.sha256).join('') + generatedAt)
  const blob = await buildReportPdf({
    files,
    donorHeterozygosity: pat.heterozygosity,
    startedAt: null,
    generatedAt,
    tool,
    reportId: id,
    fromExamples: false,
    markSvg: MARK_SVG,
  })
  const path = join(OUT, `${name}.pdf`)
  writeFileSync(path, Buffer.from(await blob.arrayBuffer()))
  return { path, bytes: blob.size, id }
}

export { eachRow, MARK_SVG, RAW, OUT }
if (import.meta.url === `file://${process.argv[1]}`) mkdirSync(OUT, { recursive: true })

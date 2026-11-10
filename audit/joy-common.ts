/** Shared loader for the three checks Joy asked for. Keeps each harness to its own question. */
import { readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
export const ingest = await import(`${W}ingest.ts`)
export const segs = await import(`${W}segments.ts`)
export const stageMod = await import(`${W}stage.ts`)

export type Marker = { chrom: string; pos: number; called: boolean; log2R: number | null }
export type Loaded = {
  id: string
  markers: Marker[]
  hetRate: number
  callRate: number
  genomeLrr: number
  nullNoCall: number
  profile: unknown
  /** Per chromosome, the called markers in order with whether each is heterozygous. */
  hetByChrom: Map<string, { pos: number; het: boolean }[]>
}

/** What a harness needs built. Everything is built by default, so no existing caller changes. */
export type LoadNeed = {
  /**
   * Per-chromosome heterozygosity in position order. Costs a push per called marker plus a sort
   * per chromosome, and the segment harnesses never read it.
   */
  hetByChrom?: boolean
}

/**
 * One line at a time off the buffer, with `String.split('\n')` semantics exactly.
 *
 * WHY NOT `split`. A 200 MB array file split on newlines holds the fully decoded text AND all
 * 825,656 line strings at once, roughly 400 MB of peak that the parse never needs more than one
 * line of. That peak is why each shard wanted a 3.5 GB heap, and why 16 GB of RAM meant six
 * parallel shards instead of twelve: it spends wall clock, not memory.
 *
 * The semantics are pinned in joy-common.check.ts rather than trusted, because an off-by-one here
 * drops the last marker of every array silently, and every recorded baseline in this repository
 * is read through this function.
 */
export const lineReader = (buf: Buffer): (() => string | null) => {
  let at = 0
  return () => {
    if (at > buf.length) return null
    // A file ending in a newline has one empty field after it, which is what `split` returns.
    if (at === buf.length) { at += 1; return '' }
    let end = buf.indexOf(10, at)
    if (end < 0) end = buf.length
    const s = buf.toString('utf8', at, end)
    at = end + 1
    return s
  }
}

export function loadArray(path: string, id: string, need: LoadNeed = {}): Loaded {
  const raw = readFileSync(path)
  const buf = path.endsWith('.gz') ? gunzipSync(raw) : raw
  const nextLine = lineReader(buf)
  // Only the first 60 lines are held, which is where the header is looked for.
  const head: string[] = []
  for (let i = 0; i < 60; i += 1) { const l = nextLine(); if (l === null) break; head.push(l) }
  let h = -1
  for (let i = 0; i < head.length; i += 1) if (head[i] && !head[i].startsWith('#')) { h = i; break }
  if (h < 0) throw new Error(`${path}: no header in the first 60 lines`)
  const map = ingest.headerMap(head[h])
  if (!map) throw new Error(`${path}: header not recognised`)

  const markers: Marker[] = []
  const hetByChrom = new Map<string, { pos: number; het: boolean }[]>()
  const byChrom = new Map()
  const bafSums = ingest.emptyBafSums()
  let first = ''
  let called = 0; let het = 0; let rows = 0
  const lrrs: number[] = []
  // The rows already read while looking for the header, then the rest of the file streamed.
  const buffered = head.slice(h + 1)
  let bi = 0
  const nextRow = (): string | null => (bi < buffered.length ? buffered[bi++] : nextLine())
  for (let line = nextRow(); line !== null; line = nextRow()) {
    const r = ingest.parseRow(line, map)
    if (!r) continue
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    rows += 1
    const isCalled = r.genotype !== 'NC'
    markers.push({ chrom: r.chrom, pos: r.pos, called: isCalled, log2R: r.log2R ?? null })
    if (r.log2R !== null && r.log2R !== undefined && Number.isFinite(r.log2R)) lrrs.push(r.log2R)
    if (!isCalled) continue
    called += 1
    const isHet = r.genotype === 'AB'
    if (isHet) het += 1
    if (need.hetByChrom === false) continue
    if (!/^\d+$/.test(r.chrom)) continue
    let a = hetByChrom.get(r.chrom)
    if (!a) hetByChrom.set(r.chrom, a = [])
    a.push({ pos: r.pos, het: isHet })
  }
  for (const a of hetByChrom.values()) a.sort((x, y) => x.pos - y.pos)
  lrrs.sort((a, b) => a - b)
  // The pipeline's own profile, passed to inferStage exactly as the application passes it. An
  // earlier version handed it `.build`, which is the array build and carries no rates, so every
  // sample staged as `unknown`.
  const profile = ingest.finishProfile(id, byChrom as never, bafSums as never, first) as {
    callRate: number; hetRate: number
  }
  return {
    id,
    markers,
    hetRate: profile.hetRate,
    callRate: profile.callRate,
    genomeLrr: lrrs.length ? lrrs[lrrs.length >> 1] : 0,
    nullNoCall: rows ? (rows - called) / rows : NaN,
    profile,
    hetByChrom,
  }
}

/** Reciprocal overlap of two intervals, the share each covers of the other, lower of the two. */
export const reciprocal = (
  aS: number, aE: number, bS: number, bE: number,
): number => {
  const inter = Math.min(aE, bE) - Math.max(aS, bS)
  if (inter <= 0) return 0
  return Math.min(inter / Math.max(1, aE - aS), inter / Math.max(1, bE - bS))
}

/** Array files, across the naming conventions in the corpora. Excludes series matrices. */
export const ARRAY_FILE = /\.(probes|CEL)(\.txt)?(\.gz)?$/i
export const isArray = (f: string): boolean => ARRAY_FILE.test(f) && !/series_matrix/i.test(f)

/**
 * Which embryo and which cell an array came from.
 *
 * GEO writes both into the filename and the trailing number is the SCAN, so the same cell appears
 * more than once: `4-2_58` and `4-2_74` are two arrays of one cell, not two sisters. Comparing
 * those as sisters would compare a cell with itself.
 */
export const cellOf = (f: string): { embryo: string; cell: string } => {
  const tok = (/^(?:GSM\d+_)?(.+?)(?:_\d+)?\.(?:probes|CEL)/i.exec(f)?.[1]) ?? f.replace(/\..*$/, '')
  const dash = tok.lastIndexOf('-')
  if (dash > 0) return { embryo: tok.slice(0, dash), cell: tok.slice(dash + 1) }
  const m = /^([A-Za-z]+)(\d+)$/.exec(tok)
  if (m) return { embryo: m[1], cell: m[2] }
  return { embryo: tok, cell: '' }
}

/**
 * The copy-number scan, invoked the way the application invokes it.
 *
 * `scoreSample` runs it PER CHROMOSOME, with a null that is the median absence rate over the
 * sample's OTHER autosomes, not a genome-wide figure:
 *
 *   [...cnByChrom].flatMap(([c, ms]) => scanCopyNumber(ms, externalNull(noCall, c), genomeLrr))
 *
 * Passing every marker in one call instead sorts positions from all chromosomes together, so the
 * windows straddle chromosome boundaries and the scan returns nothing. An earlier version of this
 * harness did exactly that and would have reported zero recall as if it were a measurement.
 */
export const scanSegments = (a: Loaded) => {
  const byChrom = new Map<string, Marker[]>()
  // PUSH, NEVER SPREAD. Rebuilding the accumulated array once per marker is quadratic: at 825,656
  // markers over 22 autosomes it is about 1.5e10 element copies and as many throwaway allocations,
  // measured at 48 seconds against 0.041 for the same grouping done by push. That one line was the
  // bulk of every benchmark run's wall clock, well past the scan it was feeding.
  for (const m of a.markers) {
    if (m.chrom.charCodeAt(0) < 48 || m.chrom.charCodeAt(0) > 57) continue
    if (!/^\d+$/.test(m.chrom)) continue
    let arr = byChrom.get(m.chrom)
    if (!arr) byChrom.set(m.chrom, arr = [])
    arr.push(m)
  }
  const noCall = new Map<string, [number, number]>()
  for (const [c, ms] of byChrom) noCall.set(c, [ms.length, ms.filter((m) => !m.called).length])
  const out: ReturnType<typeof segs.scanCopyNumber> = []
  // The intensity cut is a property of the whole autosome set, so it is taken once per sample.
  const profile = segs.highLrrProfile(byChrom)
  for (const [c, ms] of byChrom) {
    // The no-call arm first, so the gain arm can recognise what it has already raised.
    const noCallArm = segs.scanCopyNumber(
      ms, segs.externalNull(noCall, c) as number, a.genomeLrr)
    out.push(...noCallArm)
    out.push(...segs.scanCopyGain(
      ms, segs.externalNull(profile.byChrom, c) as number, a.genomeLrr, profile, noCallArm))
  }
  return out
}

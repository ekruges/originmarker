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

export function loadArray(path: string, id: string): Loaded {
  const raw = readFileSync(path)
  const text = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  if (h < 0) throw new Error(`${path}: no header in the first 60 lines`)
  const map = ingest.headerMap(lines[h])
  if (!map) throw new Error(`${path}: header not recognised`)

  const markers: Marker[] = []
  const hetByChrom = new Map<string, { pos: number; het: boolean }[]>()
  const byChrom = new Map()
  const bafSums = ingest.emptyBafSums()
  let first = ''
  let called = 0; let het = 0; let rows = 0
  const lrrs: number[] = []
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
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
    if (!/^\d+$/.test(r.chrom)) continue
    const a = hetByChrom.get(r.chrom) ?? []
    a.push({ pos: r.pos, het: isHet })
    hetByChrom.set(r.chrom, a)
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
  for (const m of a.markers) {
    if (!/^\d+$/.test(m.chrom)) continue
    byChrom.set(m.chrom, [...(byChrom.get(m.chrom) ?? []), m])
  }
  const noCall = new Map<string, [number, number]>()
  for (const [c, ms] of byChrom) noCall.set(c, [ms.length, ms.filter((m) => !m.called).length])
  const out: ReturnType<typeof segs.scanCopyNumber> = []
  for (const [c, ms] of byChrom) {
    out.push(...segs.scanCopyNumber(ms, segs.externalNull(noCall, c) as number, a.genomeLrr))
  }
  return out
}

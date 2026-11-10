/**
 * Whether an array's genotype calls carry any information about any genome.
 *
 * WHY THIS EXISTS. 5.47.0 splits a refusal into "the array failed" and "this chromosome offered
 * nothing", by call rate. Neither branch can say that the calls themselves are noise, and the
 * statistic that looked like the candidate for it does not work: on amplified material hetBafSd
 * measures how much amplification the sample went through, not how good the array is, and arrays
 * with no genotype signal at all read in the middle of the usable amplified range rather than
 * above it. Any upper threshold on it refuses good single cells and passes every failed array.
 *
 * THE STATISTIC. The information ratio R. Take the markers where the sample and a reference array
 * are both homozygous. Among those, measure how often they are homozygous for OPPOSITE alleles.
 * Divide that by the rate expected if the sample's calls said nothing about the reference, which
 * is s0(1 - r0) + (1 - s0)r0 from the two arrays' own AA fractions over the same marker set.
 *
 *   R near 1    the calls carry no information about this reference, or any genome
 *   R well below 1    the calls are informative
 *
 * The reference does not have to be a parent. Two unrelated people still share population
 * structure, so an informative array sits well below 1 against a stranger, while an array whose
 * calls are noise sits at 1 against everyone. That is what makes this a parent-free read on the
 * sample itself.
 *
 * A DNA MIXTURE READS LOW, NOT HIGH, by construction: a homozygous call in a mixture needs every
 * contributor to carry that allele. So R separates a failed array from a contaminated one, which
 * the heterozygosity ceiling on its own cannot.
 *
 * The null here is the count form. The consult response also validates a circular-shift form and
 * reports the two agreeing to within 0.0096 in R over 135 arrays, so only the count form is
 * implemented.
 *
 * Run: OM_DIR=<dir,dir,...> OM_REFS=<path,path,...> \
 *        node --experimental-strip-types audit/noise-information-ratio.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { join, basename } from 'node:path'
import { gunzipSync } from 'node:zlib'
import { ingest, isArray, cellOf } from './joy-common.ts'

const DIRS = (process.env.OM_DIR ?? '').split(',').filter((d) => d && existsSync(d))
const REFS = (process.env.OM_REFS ?? '').split(',').filter((r) => r && existsSync(r))
if (!DIRS.length) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
if (!REFS.length) { console.log('OM_REFS not readable. NOT RUN.'); process.exit(0) }

/** Bands fixed before the run, from the consult response. Reasons only; no verdict moves. */
const INFORMATIVE = 0.30
const NO_INFORMATION = 0.60
/** A reference whose A/B coding is flipped against the sample reads near 2, not 1. Reported. */
const FLIP_GUARD = 1.5

const read = (p: string): string[] => {
  const raw = readFileSync(p)
  return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8').split('\n')
}
const headerAt = (lines: string[]): number => {
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) return i
  return -1
}

/** Homozygous calls of one reference, keyed by probeset id. 0 is AA, 1 is BB. */
const loadRef = (p: string): Map<string, 0 | 1> => {
  const lines = read(p)
  const h = headerAt(lines)
  const map = h < 0 ? null : ingest.headerMap(lines[h])
  const out = new Map<string, 0 | 1>()
  if (!map) return out
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r) continue
    if (r.genotype === 'AA') out.set(r.probesetId, 0)
    else if (r.genotype === 'BB') out.set(r.probesetId, 1)
  }
  return out
}

const refs = REFS.map((p) => ({ name: basename(p).replace(/\..*$/, ''), prefix: cellOf(basename(p)).embryo,
  hom: loadRef(p) }))
console.error(`${refs.length} references loaded: ${refs.map((r) => r.hom.size).join(', ')} `
  + 'homozygous markers each')

/** The four counts the ratio needs, kept per reference and per autosome. */
type Tally = { both: number; opp: number; sAA: number; rAA: number }
const blank = (): Tally => ({ both: 0, opp: 0, sAA: 0, rAA: 0 })
const ratio = (t: Tally): number => {
  if (t.both < 500) return NaN
  const s0 = t.sAA / t.both
  const r0 = t.rAA / t.both
  const expected = s0 * (1 - r0) + (1 - s0) * r0
  return expected > 0 ? (t.opp / t.both) / expected : NaN
}
const band = (r: number): string =>
  !Number.isFinite(r) ? 'too-few-markers'
    : r >= FLIP_GUARD ? 'check-allele-coding'
      : r >= NO_INFORMATION ? 'no-information'
        : r <= INFORMATIVE ? 'informative' : 'unresolved'

const files: string[] = []
for (const d of DIRS) for (const f of readdirSync(d)) if (isArray(f)) files.push(join(d, f))
files.sort()
console.error(`${files.length} arrays to score`)

const perChrom: string[] = []
console.log(['id', 'markers', 'call_rate', 'het_rate', 'het_baf_sd', 'both_hom', 'R', 'R_min',
  'band', 'band_min', 'refs_used'].join('\t'))

for (const path of files) {
  const f = basename(path)
  const id = (/^(?:GSM\d+_)?([0-9A-Za-z]+-[0-9]+)/.exec(f)?.[1]) ?? f.replace(/\..*$/, '')
  const prefix = cellOf(f).embryo
  // A reference from the sample's own barcode prefix may be a relative, which lowers R for a
  // reason that is not array quality. Scored against strangers only.
  const use = refs.filter((r) => r.prefix !== prefix)
  if (!use.length) { console.error(`${id}: every reference shares its prefix, skipped`); continue }

  let lines: string[]
  try { lines = read(path) } catch { console.error(`${id}: unreadable, skipped`); continue }
  const h = headerAt(lines)
  const map = h < 0 ? null : ingest.headerMap(lines[h])
  if (!map) { console.error(`${id}: header not recognised, skipped`); continue }

  const byChrom = new Map()
  const bafSums = ingest.emptyBafSums()
  const whole = use.map(blank)
  const byAuto = new Map<string, Tally[]>()
  let first = ''; let rows = 0

  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    rows += 1
    const sHom = r.genotype === 'AA' ? 0 : r.genotype === 'BB' ? 1 : -1
    if (sHom < 0) continue
    const auto = /^\d+$/.test(r.chrom)
    if (auto && !byAuto.has(r.chrom)) byAuto.set(r.chrom, use.map(blank))
    const chromTal = auto ? byAuto.get(r.chrom)! : null
    for (let j = 0; j < use.length; j += 1) {
      const rHom = use[j].hom.get(r.probesetId)
      if (rHom === undefined) continue
      for (const t of chromTal ? [whole[j], chromTal[j]] : [whole[j]]) {
        t.both += 1
        if (sHom === 0) t.sAA += 1
        if (rHom === 0) t.rAA += 1
        if (sHom !== rHom) t.opp += 1
      }
    }
  }

  const profile = ingest.finishProfile(id, byChrom as never, bafSums as never, first) as
    { callRate: number; hetRate: number; hetBafSd: number }
  // The worst reference decides. An array only counts as informative if it is informative against
  // every stranger it was compared with, which is the conservative direction for that claim.
  const rs = whole.map(ratio).filter(Number.isFinite)
  const R = rs.length ? Math.max(...rs) : NaN
  // The references available in this corpus are single cells rather than bulk DNA, and a
  // reference carrying its own call errors pushes R up for every sample scored against it. The
  // best reference is reported beside the worst so that effect is visible rather than assumed.
  const Rmin = rs.length ? Math.min(...rs) : NaN
  const both = Math.max(...whole.map((t) => t.both))
  console.log([id, rows, profile.callRate.toFixed(4), profile.hetRate.toFixed(4),
    Number.isFinite(profile.hetBafSd) ? profile.hetBafSd.toFixed(4) : 'NaN',
    both, Number.isFinite(R) ? R.toFixed(4) : 'NaN',
    Number.isFinite(Rmin) ? Rmin.toFixed(4) : 'NaN',
    band(R), band(Rmin), use.length].join('\t'))

  for (const [chrom, tals] of [...byAuto].sort((a, b) => Number(a[0]) - Number(b[0]))) {
    const cr = tals.map(ratio).filter(Number.isFinite)
    perChrom.push([id, chrom, Math.max(...tals.map((t) => t.both)),
      cr.length ? Math.max(...cr).toFixed(4) : 'NaN',
      band(cr.length ? Math.max(...cr) : NaN)].join('\t'))
  }
}

if (process.env.OM_CHROM_OUT) {
  const { writeFileSync } = await import('node:fs')
  writeFileSync(process.env.OM_CHROM_OUT,
    ['id\tchrom\tboth_hom\tR\tband', ...perChrom].join('\n') + '\n')
  console.error(`per-autosome rows written: ${perChrom.length}`)
}

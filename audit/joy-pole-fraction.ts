/**
 * Homozygous-pole fraction inside each hand-called region, against the same array's genome figure.
 *
 * WHY. The gains benchmark reports 11 of 26 matched. Ten of those eleven sit in cells carrying no
 * hand-called loss at all, at log2R deviations of +1.61 to +2.83, which is above anything a third
 * copy produces: log2(3/2) is 0.585. Amplification distortion reads that way and a trisomy does
 * not, so the recall figure may be agreement on distorted material rather than detection of the
 * biology. Distorted material also collapses the homozygous BAF poles, and that is measurable here
 * without re-running detection.
 *
 * THE STATISTIC. The fraction of markers with min(BAF, 1 - BAF) below POLE, which is the share of
 * markers sitting at a homozygous pole rather than smeared toward the middle. Intact amplified
 * material reads around 0.83 genome-wide; collapsed arrays read 0.10 to 0.37.
 *
 * Detection is not re-run. Join the output on (cell, chrom, start, end) against the baseline
 * per-region table to split found from missed.
 *
 * Run: OM_DIR=<dir of .probes> OM_REGIONS=<regions-dieter.tsv> \
 *        node --experimental-strip-types audit/joy-pole-fraction.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { gunzipSync } from 'node:zlib'
import { ingest, isArray } from './joy-common.ts'

const DIR = process.env.OM_DIR
const REG = process.env.OM_REGIONS
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
if (!REG || !existsSync(REG)) { console.log('OM_REGIONS not readable. NOT RUN.'); process.exit(0) }

/** Distance from the nearer homozygous pole below which a marker counts as at the pole. */
const POLE = Number(process.env.OM_POLE ?? 0.1)

type Hand = { cell: string; chrom: string; start: number; end: number
  kind: 'gain' | 'loss'; dev: number }
const hand: Hand[] = []
for (const l of readFileSync(REG, 'utf8').split('\n')) {
  if (!l || l.startsWith('cell')) continue
  const c = l.split('\t')
  if (c.length < 6) continue
  const k = c[5].trim().toLowerCase()
  if (k !== 'gain' && k !== 'loss') continue
  hand.push({ cell: c[0].trim(), chrom: c[1].trim(), start: Number(c[2]), end: Number(c[3]),
    kind: k, dev: Number(c[4]) })
}

const byCell = new Map<string, Hand[]>()
for (const h of hand) byCell.set(h.cell, [...(byCell.get(h.cell) ?? []), h])

console.log(['cell', 'chrom', 'kind', 'mb', 'log2R_dev', 'markers_in', 'pole_in', 'pole_genome',
  'ratio'].join('\t'))

for (const f of readdirSync(DIR).filter(isArray).sort()) {
  // The leading accession-and-position, which is how the sheet names them. Same rule as
  // joy-dieter.ts, so the two harnesses key on identical cell ids and their outputs join.
  const cell = (/^(?:GSM\d+_)?([0-9]+-[0-9]+)/.exec(f)?.[1]) ?? f.replace(/\..*$/, '')
  const regions = byCell.get(cell)
  if (!regions) continue

  const raw = readFileSync(join(DIR, f))
  const text = (f.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  if (h < 0) continue
  const map = ingest.headerMap(lines[h])
  if (!map) continue

  // One pass, accumulating the genome figure and every region's figure together, because these
  // files are large enough that a second pass per region is the whole runtime.
  let gPole = 0; let gTot = 0
  const inPole = regions.map(() => 0)
  const inTot = regions.map(() => 0)
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || r.baf === null || !Number.isFinite(r.baf)) continue
    const atPole = Math.min(r.baf, 1 - r.baf) < POLE ? 1 : 0
    gTot += 1; gPole += atPole
    for (let j = 0; j < regions.length; j += 1) {
      const g = regions[j]
      if (r.chrom !== g.chrom || r.pos < g.start || r.pos > g.end) continue
      inTot[j] += 1; inPole[j] += atPole
    }
  }
  const genome = gTot ? gPole / gTot : NaN

  for (let j = 0; j < regions.length; j += 1) {
    const g = regions[j]
    const frac = inTot[j] ? inPole[j] / inTot[j] : NaN
    console.log([cell, g.chrom, g.kind, ((g.end - g.start) / 1e6).toFixed(1), g.dev.toFixed(4),
      inTot[j], frac.toFixed(4), genome.toFixed(4), (frac / genome).toFixed(4)].join('\t'))
  }
}

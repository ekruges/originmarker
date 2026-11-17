/**
 * Which chip a file comes from, and whether the tool is proven on it.
 *
 * The tool is built and measured on one chip, the UK Biobank Axiom array. Every threshold in the
 * scoring code, the founder panel and the answer keys are that chip's. On a file from another chip
 * those numbers read the wrong scale and print confident output (a "haploid" stage, a damaged-sample
 * banner) that says nothing true about the sample. So a file is scored only when its markers are
 * the proven chip's, and otherwise refused with every cause that can be named.
 *
 * RECOGNITION IS BY MARKERS, NOT BY NAME. The share of the file's autosomal markers whose id the
 * shipped founder panel holds. The panel carries every autosomal marker of the Axiom chip, so an
 * Axiom file reads all of them known and another chip reads none: see
 * audit/truth/unknown-chip-2026-11-16.txt for the measured shares and why the bar sits where it
 * does. A ratio of the rows present, so an export missing rows reads the same share as a whole one.
 *
 * This is the place a per-chip profile object belongs when chips other than Axiom are proven.
 */
import type { BuildCall, ColumnMap } from './ingest.ts'

/** The chip the tool is proven on. */
export const PROVEN_CHIP = 'UK Biobank Axiom'

/**
 * Whether the proven chip genotypes its chromosome Y markers. It does: a man's array calls nearly
 * all of its 812 and a woman's none, so on this chip a file calling none of them has no Y. A chip
 * that never genotypes chrY (the CytoScan 750K calls none on a man or a woman) would be read by
 * intensity instead; see sexChromosome.ts yPresent.
 */
export const PROVEN_CHIP_GENOTYPES_Y = true

/**
 * Least share of a file's autosomal markers the founder panel must hold for the file to be the
 * proven chip. Measured: 1.0000 on every Axiom file read (362 arrays and three rebuilt fathers),
 * 0.0000 on both outside chips (68 arrays). Half way between, so that neither a few
 * markers renamed on Axiom nor a few shared ids on another chip moves a file across.
 */
export const KNOWN_SHARE_MIN = 0.5

/** The names each column is looked for under. ingest.ts resolves headers from these lists, so a
 *  refusal names exactly what was tried. The intensity and allele-fraction names also match after
 *  an export's per-sample prefix (`GSM123.Log R Ratio`). */
export const COLUMN_NAMES: Record<string, readonly string[]> = {
  probeset_id: ['probeset_id', 'snp name', 'name'],
  chr: ['chr'],
  position: ['position'],
  genotype: ['genotype'],
  log2r: ['log2r', 'normalized_intensity', 'lrr', 'log_r_ratio', 'logrratio', 'log r ratio'],
  baf: ['baf', 'b_allele_freq', 'b allele freq', 'ballelefreq', 'b_allele_frequency', 'b allele frequency'],
}

/** Plain names for the reader, in refusals and log lines. */
const COLUMN_LABEL: Record<string, string> = {
  probeset_id: 'marker id', chr: 'chromosome', position: 'position', genotype: 'genotype call',
  log2r: 'intensity (log2 ratio)', baf: 'B allele fraction',
}

export const REQUIRED_COLUMNS = ['probeset_id', 'chr', 'position', 'genotype'] as const
export const SIGNAL_COLUMNS = ['log2r', 'baf'] as const

/** A column the reader looked for and did not find, with the names it tried. */
export const missingColumn = (key: string): string =>
  `no ${COLUMN_LABEL[key] ?? key} column (looked for ${COLUMN_NAMES[key].join(', ')})`

/** The required columns a header line lacks, for a file whose header the reader could not resolve. */
export function missingRequired(headerLine: string): string[] {
  const best = ['\t', ',']
    .map((d) => headerLine.trim().split(d).map(norm))
    .sort((a, b) => b.length - a.length)[0]
  return REQUIRED_COLUMNS.filter((k) => !COLUMN_NAMES[k].some((n) => best.includes(n))).map(missingColumn)
}

/** The intensity and allele-fraction columns a resolved header lacks. Named, never dropped. */
export const missingSignals = (map: ColumnMap): string[] =>
  SIGNAL_COLUMNS.filter((k) => typeof map[k] !== 'number').map(missingColumn)

/** Lower-cased and trimmed, as headers are compared. */
export const norm = (c: string): string => c.trim().toLowerCase()

/** A header cell matches a name exactly, or as the suffix after an export's `Sample.` prefix. */
export const matches = (cell: string, name: string): boolean =>
  cell === name || cell.endsWith(`.${name}`)

export interface ChipTally { markers: number; known: number }
export const emptyChipTally = (): ChipTally => ({ markers: 0, known: 0 })
type Known = ReadonlySet<string> | ReadonlyMap<string, unknown>

/** Autosomes only: the panel holds the founders' autosomal markers, and none on X, Y or MT. */
export function tallyChip(
  t: ChipTally, id: string, chrom: string, known: Known,
): void {
  if (!/^([1-9]|1\d|2[0-2])$/.test(chrom)) return
  t.markers += 1
  if (known.has(id)) t.known += 1
}

export type ChipCall =
  | { proven: true; chip: string; knownShare: number }
  | { proven: false; knownShare: number; causes: string[]; text: string }

const pct = (x: number) => `${(100 * x).toFixed(1)}%`
const int = (x: number) => x.toLocaleString('en-US')

/**
 * The verdict for one file. `build` and `signals` only add causes to a refusal: on the proven chip
 * an unrecognised build or a missing column is named by the reader and the gates, not refused here.
 */
export function recogniseChip(t: ChipTally, build: BuildCall, missing: string[] = []): ChipCall {
  const knownShare = t.markers ? t.known / t.markers : 0
  if (t.markers > 0 && knownShare >= KNOWN_SHARE_MIN) return { proven: true, chip: PROVEN_CHIP, knownShare }
  const causes = [
    t.markers
      ? `the founder panel holds ${int(t.known)} of its ${int(t.markers)} autosomal markers (${pct(knownShare)}); `
        + `a ${PROVEN_CHIP} file has all of its autosomal markers in the panel`
      : 'it has no markers on chromosomes 1 to 22',
  ]
  if (!build.build) {
    causes.push(build.tested > 0 && build.illegal.GRCh37 > 0 && build.illegal.GRCh38 > 0
      ? `genome build not recognised: ${int(build.illegal.GRCh37)} of ${int(build.tested)} markers fall `
        + `outside GRCh37 and ${int(build.illegal.GRCh38)} outside GRCh38`
      : `genome build not recognised: ${build.note}`)
  }
  causes.push(...missing)
  return { proven: false, knownShare, causes, text: refusalText(causes) }
}

/** The sentence a reader sees, for a chip refusal or an unreadable header. */
export const refusalText = (causes: string[]): string =>
  `Not scored: this file is not from the ${PROVEN_CHIP} chip, the only chip this tool is proven on, `
  + `so none of its numbers apply. ${causes.map((c) => c[0].toUpperCase() + c.slice(1)).join('. ')}.`

/** The sentence for a file whose header the reader could not resolve at all. */
export const unreadableText = (missing: string[]): string =>
  `Not scored: the reader could not find the columns it needs, so the chip cannot be checked. `
  + `${missing.map((c) => c[0].toUpperCase() + c.slice(1)).join('. ')}.`

/** The proven-chip verdict for rows already parsed, as every reader that holds a whole file uses it. */
export function recogniseRows(
  rows: Iterable<{ probesetId: string; chrom: string }>, build: BuildCall, map: ColumnMap, known: Known,
): ChipCall {
  const t = emptyChipTally()
  for (const r of rows) tallyChip(t, r.probesetId, r.chrom, known)
  return recogniseChip(t, build, missingSignals(map))
}

/**
 * Where a header the reader could not resolve still names the marker id (an Illumina Final
 * Report has `SNP Name` but no chromosome, position or single genotype column), so the chip can
 * be read from the ids alone.
 */
export function idColumn(headerLine: string): { index: number; delim: string } | null {
  for (const delim of ['\t', ',']) {
    const cols = headerLine.trim().split(delim).map(norm)
    const index = cols.findIndex((c) => COLUMN_NAMES.probeset_id.includes(c))
    if (index >= 0 && cols.length > 1) return { index, delim }
  }
  return null
}
export const idAt = (line: string, col: { index: number; delim: string }): string =>
  (line.split(col.delim)[col.index] ?? '').trim()

/** Share of every row's id the panel holds, with no chromosome to restrict to. Measured on Axiom
 *  exports: 0.9737, the X and Y rows being the rest. */
export const AXIOM_ALL_ROWS_SHARE = '97%'

/**
 * The refusal for a file whose header lacks a required column. Read as another chip when the ids
 * under a marker column show it, with the panel share and the missing columns as causes; otherwise
 * the columns are named and the chip is said to be unchecked.
 */
export function unreadableRefusal(missing: string[], ids: ChipTally | null): string {
  if (ids && ids.markers > 0 && ids.known / ids.markers < KNOWN_SHARE_MIN) {
    return refusalText([
      `the founder panel holds ${int(ids.known)} of its ${int(ids.markers)} markers `
        + `(${pct(ids.known / ids.markers)}); a ${PROVEN_CHIP} file has ${AXIOM_ALL_ROWS_SHARE} of its `
        + 'markers in the panel',
      ...missing,
    ])
  }
  return unreadableText(missing)
}

/** The same, for a reader holding the whole file as lines. */
export function unreadableLines(lines: readonly string[], known: Known | null | undefined): string {
  let header = -1
  let col: { index: number; delim: string } | null = null
  for (let i = 0; i < Math.min(lines.length, 200); i += 1) {
    const l = lines[i]
    if (!l || !l.trim() || l.startsWith('#')) continue
    if (header < 0) header = i
    col = idColumn(l)
    if (col) { header = i; break }
  }
  const missing = header < 0 ? REQUIRED_COLUMNS.map(missingColumn) : missingRequired(lines[header])
  if (!col || !known) return unreadableRefusal(missing, null)
  const t = emptyChipTally()
  for (let i = header + 1; i < lines.length; i += 1) {
    const id = lines[i] ? idAt(lines[i], col) : ''
    if (!id) continue
    t.markers += 1
    if (known.has(id)) t.known += 1
  }
  return unreadableRefusal(missing, t)
}

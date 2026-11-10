/**
 * The segment scan's specificity against a PUBLISHED karyotype, not against a local table.
 *
 * WHY THIS ARM EXISTS. The gains benchmark scores recall and false calls against
 * audit/local/features/regions-dieter.tsv, whose provenance is unsettled: all 93 rows are an exact
 * subset of this project's own output, the deviation column agrees to four decimals, every
 * boundary lands on a 1 Mb grid, and the whole table reproduces by filtering that output. So a
 * false call measured there may be a gap in the table rather than an error in the tool, and no
 * amount of scanning settles which.
 *
 * Zuccaro 2020 Table S3 karyotyped every arrayed cell by eye and published the result, one call
 * per autosome. 477 of those rows read `normal`. A segmental call landing on one of them disagrees
 * with a number somebody else published and defended in review, which is the specificity statement
 * worth making.
 *
 * WHAT A DISAGREEMENT IS AND IS NOT. Table S3 is a WHOLE-CHROMOSOME reading. It cannot confirm or
 * deny a segmental event, so a segment on a chromosome called normal is a disagreement and not a
 * proven error: a hand reading of whole-chromosome dosage can miss a 10 Mb duplication. It is
 * reported as a disagreement, with its span and intensity, exactly as the Dieter arm is.
 *
 * THE BAR, FIXED BEFORE THIS RUNS. The Dieter arm allows 1 to 2 false segments over 72 arrays,
 * about 1,580 autosome rows. 477 rows is 30% of that, so the proportionate bar is 1 disagreement
 * on a normal chromosome. At 2 or more the intensity nomination does not ship on this evidence.
 * Calls on chromosomes the table reads monosomic or nullisomic are reported separately and are
 * not disagreements: those are chromosomes where something IS wrong.
 *
 * Run: OM_DIR=<dir of GSM*.probes.gz> OM_TRUTH=audit/zuccaro-table-s3.tsv \
 *        node --experimental-strip-types audit/joy-public-specificity.ts
 */
import { readFileSync, existsSync, readdirSync, writeFileSync } from 'node:fs'
import { join } from 'node:path'
import { loadArray, scanSegments, isArray } from './joy-common.ts'

const DIR = process.env.OM_DIR
const TRUTH = process.env.OM_TRUTH ?? 'audit/zuccaro-table-s3.tsv'
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
if (!existsSync(TRUTH)) { console.log('OM_TRUTH not readable. NOT RUN.'); process.exit(0) }

/** cell -> chromosome -> the table's own word for it. */
const truth = new Map<string, Map<string, string>>()
const karyotype = new Map<string, string>()
for (const l of readFileSync(TRUTH, 'utf8').split('\n')) {
  if (!l || l.startsWith('#')) continue
  const c = l.split('\t')
  if (c.length < 5 || c[0] === 'gsm' || !c[0]) continue
  const cell = c[0].trim()
  if (!truth.has(cell)) truth.set(cell, new Map())
  truth.get(cell)!.set(c[3].trim().replace(/^chr/i, ''), c[4].trim().toLowerCase())
  karyotype.set(cell, c[2].trim())
}
let normalRows = 0
for (const m of truth.values()) for (const v of m.values()) if (v === 'normal') normalRows += 1
console.log(`${truth.size} cells with a published karyotype, ${normalRows} autosomes read normal\n`)

const SHARDS = Number(process.env.OM_SHARDS ?? 1)
const SHARD = Number(process.env.OM_SHARD ?? 0)
const files = readdirSync(DIR).filter(isArray).sort().filter((_, i) => i % SHARDS === SHARD)

type Row = { cell: string; chrom: string; kind: string; mb: number; verdict: string; kary: string }
const rows: Row[] = []
let scanned = 0
let withheld = 0
for (const f of files) {
  const cell = /^(GSM\d+)/.exec(f)?.[1]
  if (!cell || !truth.has(cell)) continue
  const t = truth.get(cell)!
  let segs
  try { segs = scanSegments(loadArray(join(DIR, f), cell)) } catch (e) {
    console.log(`  ${cell}: ${(e as Error).message}`); continue
  }
  scanned += 1
  // The shipped tool withholds the intensity scan on a genome with no normal majority, and the
  // harness now mirrors that, so an empty result on a chaotic cell is a refusal and not a pass.
  if (!segs.length) { withheld += 1; continue }
  for (const s of segs) {
    if (s.kind === 'parental-absence') continue
    rows.push({ cell, chrom: s.chrom, kind: s.kind,
      mb: (s.endBp - s.startBp) / 1e6,
      verdict: t.get(s.chrom) ?? 'not in table',
      kary: karyotype.get(cell) ?? '' })
  }
}
console.log(`${scanned} of the published cells scanned, ${withheld} returned no segment\n`)

const onNormal = rows.filter((r) => r.verdict === 'normal')
const onAbnormal = rows.filter((r) => r.verdict !== 'normal' && r.verdict !== 'not in table')
console.log('DISAGREEMENTS WITH THE PUBLISHED KARYOTYPE')
console.log(`  segmental calls on a chromosome published NORMAL: ${onNormal.length}   (bar: 1)`)
for (const r of onNormal) {
  console.log(`    ${r.cell} chr${r.chrom} ${r.kind} ${r.mb.toFixed(1)}Mb   published: ${r.kary}`)
}
console.log(`\n  calls on a chromosome published monosomic or nullisomic: ${onAbnormal.length}`)
const byVerdict = new Map<string, number>()
for (const r of onAbnormal) byVerdict.set(r.verdict, (byVerdict.get(r.verdict) ?? 0) + 1)
for (const [v, n] of [...byVerdict].sort((a, b) => b[1] - a[1])) console.log(`    ${n}  ${v}`)

if (process.env.OM_TSV) {
  const head = ['cell', 'chrom', 'kind', 'mb', 'published_chrom_call', 'published_karyotype']
  writeFileSync(process.env.OM_TSV, [head.join('\t'),
    ...rows.map((r) => [r.cell, r.chrom, r.kind, r.mb.toFixed(1), r.verdict, r.kary].join('\t')),
  ].join('\n') + '\n')
  console.log(`\nper-call rows written to ${process.env.OM_TSV}`)
}

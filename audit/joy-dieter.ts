/**
 * Joy's ask: benchmark the shipped segment detector against Dieter Egli's manual calls.
 *
 * His calls are intensity-based regions, 26 gains and 67 losses over 21 of 93 arrays, in
 * audit/local/features/regions-dieter.tsv. The detector run here is `scanCopyNumber`, which is the
 * same function the web tool uses and needs no parent, so the comparison is like for like: both
 * sides read log R and the no-call channel on the same arrays.
 *
 * The 72 arrays he called nothing on are the false-call arm, and they matter as much as the recall.
 *
 * Run: OM_DIR=<dir of .probes> OM_REGIONS=<regions-dieter.tsv> \
 *        node --experimental-strip-types audit/joy-dieter.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { loadArray, segs, reciprocal, isArray, scanSegments } from './joy-common.ts'

const DIR = process.env.OM_DIR
const REG = process.env.OM_REGIONS
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
if (!REG || !existsSync(REG)) { console.log('OM_REGIONS not readable. NOT RUN.'); process.exit(0) }
const OVERLAP = Number(process.env.OM_OVERLAP ?? 0.5)

type Hand = { cell: string; chrom: string; start: number; end: number; kind: 'gain' | 'loss' }
const hand: Hand[] = []
for (const l of readFileSync(REG, 'utf8').split('\n')) {
  if (!l || l.startsWith('cell')) continue
  const c = l.split('\t')
  if (c.length < 6) continue
  const kind = /gain/i.test(c[5]) ? 'gain' as const : /loss/i.test(c[5]) ? 'loss' as const : undefined
  if (!kind) continue
  hand.push({ cell: c[0].trim(), chrom: c[1].trim(), start: Number(c[2]), end: Number(c[3]), kind })
}
const handCells = new Set(hand.map((h) => h.cell))
console.log(`${hand.length} hand-called regions over ${handCells.size} cells: `
  + `${hand.filter((h) => h.kind === 'gain').length} gains, `
  + `${hand.filter((h) => h.kind === 'loss').length} losses`)
console.log(`reciprocal overlap required: ${OVERLAP}\n`)

const files = readdirSync(DIR).filter(isArray).sort()
console.log(`${files.length} arrays in ${DIR}\n`)

type Row = { cell: string; called: ReturnType<typeof segs.scanCopyNumber> }
const rows: Row[] = []
let n = 0
for (const f of files) {
  // The cell id is the leading accession-and-position, which is how his sheet names them.
  const cell = (/^(?:GSM\d+_)?([0-9]+-[0-9]+)/.exec(f)?.[1]) ?? f.replace(/\..*$/, '')
  try {
    const a = loadArray(join(DIR, f), cell)
    const called = scanSegments(a)
    rows.push({ cell, called })
    n += 1
    if (n % 10 === 0) console.log(`  ... ${n} of ${files.length} arrays scanned`)
  } catch (e) {
    console.log(`  ${cell}: ${(e as Error).message}`)
  }
}
console.log(`\n${rows.length} arrays scanned\n`)

// RECALL, against his calls.
let hit = 0; const missed: string[] = []
const byKind = new Map<string, { hit: number; total: number }>([
  ['gain', { hit: 0, total: 0 }], ['loss', { hit: 0, total: 0 }]])
for (const h of hand) {
  const k = byKind.get(h.kind)!
  k.total += 1
  const row = rows.find((r) => r.cell === h.cell)
  const want = h.kind === 'gain' ? 'copy-gain' : 'copy-loss'
  const found = row?.called.some((s) => {
    const co = segs.segmentCoords(s)
    return s.chrom === h.chrom && s.kind === want
      && reciprocal(h.start, h.end, co.start, co.end) >= OVERLAP
  })
  if (found) { hit += 1; k.hit += 1 } else {
    missed.push(`${h.cell} chr${h.chrom} ${(h.start / 1e6).toFixed(1)}-`
      + `${(h.end / 1e6).toFixed(1)}Mb ${h.kind}${row ? '' : '  (no array for this cell)'}`)
  }
}
console.log('RECALL AGAINST DIETER EGLI\'S MANUAL CALLS')
for (const [kind, k] of byKind) {
  console.log(`  ${kind.padEnd(6)} ${k.hit} of ${k.total}`
    + `  (${(100 * k.hit / Math.max(1, k.total)).toFixed(0)}%)`)
}
console.log(`  overall ${hit} of ${hand.length} `
  + `(${(100 * hit / Math.max(1, hand.length)).toFixed(0)}%)`)
console.log('\n  missed:')
for (const m of missed.slice(0, 40)) console.log(`    ${m}`)
if (missed.length > 40) console.log(`    and ${missed.length - 40} more`)

// FALSE CALLS, on the arrays he called nothing on.
const quiet = rows.filter((r) => !handCells.has(r.cell))
let extra = 0
const extraDetail: string[] = []
for (const r of quiet) {
  for (const s of r.called) {
    extra += 1
    const co = segs.segmentCoords(s)
    extraDetail.push(`${r.cell} chr${s.chrom} ${(co.start / 1e6).toFixed(1)}-`
      + `${(co.end / 1e6).toFixed(1)}Mb ${s.kind}`)
  }
}
console.log(`\nFALSE-CALL ARM: ${quiet.length} arrays he called nothing on`)
console.log(`  segments the detector raised there: ${extra}`)
console.log(`  arrays with at least one: ${new Set(extraDetail.map((d) => d.split(' ')[0])).size}`)
for (const d of extraDetail.slice(0, 25)) console.log(`    ${d}`)
if (extraDetail.length > 25) console.log(`    and ${extraDetail.length - 25} more`)

// And what it raised on his cells that he did NOT call, which is a different question from recall.
let onHisCells = 0
for (const r of rows.filter((x) => handCells.has(x.cell))) {
  for (const s of r.called) {
    const co = segs.segmentCoords(s)
    const matched = hand.some((h) => h.cell === r.cell && h.chrom === s.chrom
      && reciprocal(h.start, h.end, co.start, co.end) >= OVERLAP)
    if (!matched) onHisCells += 1
  }
}
console.log(`\nON HIS 21 CALLED CELLS, segments not matching any of his: ${onHisCells}`)
console.log('\nHis calls are a comparison, not a truth: both readings come from the same arrays and')
console.log('neither was confirmed by dissection. A disagreement is reported as a disagreement.')

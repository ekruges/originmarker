/**
 * Joy's ask: do sister cells of one embryo show RECIPROCAL events?
 *
 * Vanneste 2009 reported that segmental deletions and duplications in cleavage-stage embryos are
 * reciprocal between sister blastomeres: one cell gains the piece the other lost. That is the
 * signature of a mis-segregation AFTER the two lineages split, because the material has to go
 * somewhere. A shared event in the same direction points the other way, to something that happened
 * before the split or in the gamete.
 *
 * This groups the arrays by embryo, runs the shipped copy-number scan on each cell, and sorts the
 * overlaps into reciprocal and shared.
 *
 * Run: OM_DIR=<dir> node --experimental-strip-types audit/joy-sisters.ts
 */
import { existsSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { loadArray, segs, reciprocal, isArray, cellOf, scanSegments } from './joy-common.ts'

const DIR = process.env.OM_DIR
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
const OVERLAP = Number(process.env.OM_OVERLAP ?? 0.5)

const files = readdirSync(DIR).filter(isArray).sort()
// Group by CELL first, so two scans of one cell are replicates and not sisters.
const byCell = new Map<string, string[]>()
const cellsOfEmbryo = new Map<string, Set<string>>()
for (const f of files) {
  const { embryo, cell } = cellOf(f)
  if (!cell) continue
  const key = `${embryo}:${cell}`
  byCell.set(key, [...(byCell.get(key) ?? []), f])
  cellsOfEmbryo.set(embryo, (cellsOfEmbryo.get(embryo) ?? new Set()).add(key))
}
const multi = [...cellsOfEmbryo].filter(([, cs]) => cs.size >= 2)
  .sort((a, b) => b[1].size - a[1].size)
console.log(`${files.length} arrays, ${byCell.size} distinct cells, ${cellsOfEmbryo.size} embryos`)
console.log(`${multi.length} embryo(s) with 2 or more distinct cells: `
  + multi.map(([e, cs]) => `${e} x${cs.size}`).join(', ') + '\n')
if (!multi.length) { console.log('no embryo has two distinct cells here. NOT RUN.'); process.exit(0) }
// NOT every group is an embryo. In GSE148488 `pMII` and `mMII` are paternal and maternal
// pronuclei, each from a different oocyte, so two of them are not sister cells and a reciprocal
// pair between them would mean nothing. Those groups are reported separately rather than dropped,
// since the grouping comes from the filenames and not from a sheet.
const notSisters = (e: string) => /MII|ES$/i.test(e)
const flagged = multi.filter(([e]) => notSisters(e)).map(([e]) => e)
if (flagged.length) {
  console.log(`NOT sister cells, reported separately: ${flagged.join(', ')}. These are pronuclei or`)
  console.log('a cell line, so a reciprocal pair among them carries no meaning.\n')
}

type Seg = { chrom: string; start: number; end: number; kind: string }
let reciprocalPairs = 0
let sharedPairs = 0
const recipDetail: string[] = []
const sharedDetail: string[] = []

const scanOf = new Map<string, Seg[]>()
const scan = (f: string): Seg[] | undefined => {
  if (scanOf.has(f)) return scanOf.get(f)
  try {
    const a = loadArray(join(DIR, f), f)
    if (!(a.callRate >= 0.5)) {
      console.log(`  ${f}: call rate ${a.callRate.toFixed(2)}, skipped`)
      return undefined
    }
    const out = scanSegments(a).map((s) => {
      const co = segs.segmentCoords(s)
      return { chrom: s.chrom, start: co.start, end: co.end, kind: s.kind as string }
    })
    scanOf.set(f, out)
    return out
  } catch (e) {
    console.log(`  ${f}: ${(e as Error).message}`)
    return undefined
  }
}

for (const [embryo, cellKeys] of multi) {
  const keys = [...cellKeys].sort()
  // One array per cell for the sister comparison, the first that scans.
  const pick = new Map<string, Seg[]>()
  for (const k of keys) {
    for (const f of byCell.get(k)!) {
      const s = scan(f)
      if (s) { pick.set(k, s); break }
    }
  }
  const ids = [...pick.keys()]
  if (ids.length < 2) continue
  const total = ids.reduce((n, i) => n + pick.get(i)!.length, 0)
  const before = { r: reciprocalPairs, s: sharedPairs }
  for (let i = 0; i < ids.length; i += 1) {
    for (let j = i + 1; j < ids.length; j += 1) {
      for (const a of pick.get(ids[i])!) {
        for (const b of pick.get(ids[j])!) {
          if (a.chrom !== b.chrom) continue
          const ov = reciprocal(a.start, a.end, b.start, b.end)
          if (ov < OVERLAP) continue
          const opposite = (a.kind === 'copy-gain' && b.kind === 'copy-loss')
            || (a.kind === 'copy-loss' && b.kind === 'copy-gain')
          const line = `  ${embryo} ${ids[i]} vs ${ids[j]} chr${a.chrom} `
            + `${(a.start / 1e6).toFixed(1)}-${(a.end / 1e6).toFixed(1)}Mb: ${a.kind} / ${b.kind}, `
            + `overlap ${ov.toFixed(2)}`
          const tag = notSisters(embryo) ? ' [NOT SISTERS]' : ''
          if (opposite) { reciprocalPairs += 1; recipDetail.push(line + tag) } else {
            sharedPairs += 1; sharedDetail.push(line + tag)
          }
        }
      }
    }
  }
  // Printed per embryo rather than only at the end, so a run that is cut short still leaves
  // numbers behind. An earlier run was killed by its timeout and left none.
  console.log(`embryo ${embryo}: ${ids.length} cells, ${total} segments, `
    + `+${reciprocalPairs - before.r} reciprocal, +${sharedPairs - before.s} shared `
    + `(running ${reciprocalPairs} / ${sharedPairs})`)
}

// FREE CONTROL. Two scans of the SAME cell say how reproducible the scan itself is, which sets the
// floor for believing anything about two different cells.
let repCells = 0; let repAgree = 0; let repTotal = 0
for (const [key, fs] of byCell) {
  if (fs.length < 2) continue
  const a = scan(fs[0]); const b = scan(fs[1])
  if (!a || !b) continue
  repCells += 1
  for (const x of a) {
    repTotal += 1
    if (b.some((y) => y.chrom === x.chrom && y.kind === x.kind
      && reciprocal(x.start, x.end, y.start, y.end) >= OVERLAP)) repAgree += 1
  }
}
if (repCells) {
  console.log(`\nREPRODUCIBILITY, two scans of the same cell across ${repCells} cells`)
  console.log(`  ${repAgree} of ${repTotal} segments from the first scan reappear in the second `
    + `(${(100 * repAgree / Math.max(1, repTotal)).toFixed(0)}%)`)
  console.log('  Anything below this is noise, so it is the floor for the sister comparison.')
}

console.log('\nRECIPROCAL, one cell gained what the other lost')
console.log(`  ${reciprocalPairs} pair(s)`)
for (const d of recipDetail.slice(0, 20)) console.log(d)
if (recipDetail.length > 20) console.log(`  and ${recipDetail.length - 20} more`)

console.log('\nSHARED, both cells the same direction')
console.log(`  ${sharedPairs} pair(s)`)
for (const d of sharedDetail.slice(0, 20)) console.log(d)
if (sharedDetail.length > 20) console.log(`  and ${sharedDetail.length - 20} more`)

console.log('\nA reciprocal pair means the mis-segregation happened after those two lineages split,')
console.log('because the piece one cell lost is the piece the other gained. A shared pair in the')
console.log('same direction points earlier, to before the split or to the gamete.')
console.log(reciprocalPairs > 0
  ? `\nReciprocal events ARE present in this material: ${reciprocalPairs} against ${sharedPairs} shared.`
  : `\nNo reciprocal pair found. ${sharedPairs} shared pair(s). Either the scan is not resolving\n`
    + 'the segmental events, or this material does not carry them at a detectable size.')

/**
 * IS A DISAGREEING ROW A WRONG CALL OR A WRONG PAIRING?
 *
 * A sample reading absent against the parent assigned to it has two explanations and they are
 * distinguishable. If the pairing is wrong, some OTHER array in the directory is its parent, and
 * the opposite-homozygote rate against that one lands in the product range. If the pairing is
 * right, no array does, and the sample is simply not that parent's product.
 *
 * So this scores one named sample against EVERY array present and prints the rate for each, with
 * the relationship verdict the shipped `relate` returns. The parental array's own standing is
 * checked the same way: an array that is truly the mother of a group shows the product rate
 * against that group's maternal members, and that is what licenses reading absence against her.
 *
 * `relatednessAssessable` is reported rather than assumed: on amplified material the
 * opposite-homozygote gate does not separate a parent from a stranger, so the verdict column is
 * advisory and the RATE is the thing to read.
 *
 * Run: OM_LAB=<dir> OM_FOCUS=<barcode,barcode,...> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/spot-check-pairing.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const rel = await import(`${W}relatedness.ts`)

const DIR = process.env.OM_LAB
const FOCUS = (process.env.OM_FOCUS ?? '').split(',').map((x) => x.trim()).filter(Boolean)
if (!DIR || !existsSync(DIR) || !FOCUS.length) {
  console.log('spot-check-pairing: OM_LAB and OM_FOCUS must both be set. NOT RUN.')
  process.exit(0)
}

const fileOf = new Map<string, string>()
const walk = (d: string) => {
  for (const e of readdirSync(d)) {
    const p = join(d, e)
    if (statSync(p).isDirectory()) { walk(p); continue }
    const m = /(\d{4,6}-\d{2})/.exec(e)
    if (m && /\.probes(\.txt)?(\.gz)?$/.test(e) && !fileOf.has(m[1])) fileOf.set(m[1], p)
  }
}
walk(DIR)

function load(bc: string) {
  const path = fileOf.get(bc)!
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])!
  const gt = new Map<string, string>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  let called = 0; let het = 0; let n = 0
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    n += 1
    if (r.genotype === 'NC') continue
    called += 1
    if (r.genotype === 'AB') het += 1
    gt.set(r.probesetId, r.genotype)
    pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  return { gt, pos, callRate: called / n, het: het / called }
}

const all = [...fileOf.keys()].sort()
const pct = (x: number, d = 2) => (Number.isFinite(x) ? `${(x * 100).toFixed(d)}%` : '-')

for (const focus of FOCUS) {
  if (!fileOf.has(focus)) { console.log(`${focus}: no array found\n`); continue }
  const a = load(focus)
  console.log(`=== ${focus}  call ${pct(a.callRate, 1)}  genotype het ${pct(a.het, 1)} ===`)
  console.log(`${'against'.padEnd(10)} ${'opp hom'.padStart(8)} ${'markers'.padStart(10)} `
    + `${'win spread'.padStart(11)} relationship`)
  const scored: { bc: string; rate: number; n: number; v: string; spread: number }[] = []
  for (const bc of all) {
    if (bc === focus) continue
    const b = load(bc)
    const r = rel.relate({ gt: b.gt, pos: b.pos } as never, { gt: a.gt, pos: a.pos } as never,
      rel.OPPOSITE_HOM_MAX) as {
        opp: { rate: number; n: number }; win: { spread: number }; relationship: string
      }
    scored.push({ bc, rate: r.opp.rate, n: r.opp.n, v: r.relationship, spread: r.win.spread })
  }
  for (const s of scored.sort((x, y) => x.rate - y.rate)) {
    console.log(`${s.bc.padEnd(10)} ${pct(s.rate).padStart(8)} ${s.n.toLocaleString().padStart(10)} `
      + `${pct(s.spread, 1).padStart(11)} ${s.v}`)
  }
  console.log()
}

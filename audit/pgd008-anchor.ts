/**
 * The Mendelian anchor, measured rather than assumed.
 *
 * Every formulation of the reconstruction rests on one claim: at a marker where the known parent is
 * homozygous, a blastomere reads heterozygous far more often when the missing parent is homozygous
 * for the other allele than when she is heterozygous. This measures the three rates directly, using
 * the held-out array for DIAGNOSIS and not for building, which is why it is a separate file.
 *
 * Run: OM_DIR=<dir> OM_FAMILY=PGD008 node --experimental-strip-types audit/pgd008-anchor.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)

const DIR = process.env.OM_DIR
const FAM = process.env.OM_FAMILY ?? 'PGD008'
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: unreadable header`)
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r }
  }
}
const gtOf = (f: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(join(DIR, f))()) {
    if (!/^\d+$/.test(r.chrom)) continue
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const files = readdirSync(DIR).filter((f) => f.endsWith('.probes.gz'))
const fa = gtOf(files.find((f) => new RegExp(`Father_${FAM}\\.`).test(f))!)
const mo = gtOf(files.find((f) => new RegExp(`Mother_${FAM}\\.`).test(f))!)
const bl = files.filter((f) => new RegExp(`E\\d+_Bl\\d+_${FAM}C\\d`).test(f))
console.log(`${FAM}: known parent ${fa.size.toLocaleString()}, held-out ${mo.size.toLocaleString()}, `
  + `${bl.length} blastomere arrays\n`)

// First: is she the mother at all? The project's own statistic, on her against the known parent.
let oppN = 0; let opp = 0
for (const [probe, g] of fa) {
  if (g === 'AB') continue
  const t = mo.get(probe)
  if (!t || t === 'AB') continue
  oppN += 1
  if (t !== g) opp += 1
}
console.log(`opposite-homozygote rate between the two parents: ${(100 * opp / oppN).toFixed(2)}% of `
  + `${oppN.toLocaleString()} markers. Two unrelated adults should read high here.\n`)

// The three classes, and what a blastomere reads in each.
type A = { n: number; het: number; nc: number }
const cls = new Map<string, A>([['wanted', { n: 0, het: 0, nc: 0 }],
  ['her het', { n: 0, het: 0, nc: 0 }], ['same as known', { n: 0, het: 0, nc: 0 }]])
const classOf = new Map<string, string>()
for (const [probe, g] of fa) {
  if (g === 'AB') continue
  const t = mo.get(probe)
  if (!t) continue
  classOf.set(probe, t === 'AB' ? 'her het' : t === g ? 'same as known' : 'wanted')
}
for (const b of bl) {
  for (const r of rowsOf(join(DIR, b))()) {
    const k = classOf.get(r.probesetId)
    if (!k) continue
    const a = cls.get(k)!
    if (r.genotype === 'NC') { a.nc += 1; continue }
    a.n += 1
    if (r.genotype === 'AB') a.het += 1
  }
}
console.log(`${'the missing parent reads'.padEnd(16)} ${'markers'.padStart(9)} ${'called'.padStart(10)} `
  + `${'blastomere reads AB'.padStart(20)}`)
for (const [k, a] of cls) {
  const marks = [...classOf.values()].filter((x) => x === k).length
  console.log(`${k.padEnd(16)} ${marks.toLocaleString().padStart(9)} `
    + `${(a.n / Math.max(1, a.n + a.nc)).toFixed(3).padStart(10)} `
    + `${(a.het / Math.max(1, a.n)).toFixed(4).padStart(20)}`)
}
const w = cls.get('wanted')!; const h = cls.get('her het')!
const rw = w.het / Math.max(1, w.n); const rh = h.het / Math.max(1, h.n)
console.log(`\nthe ratio the whole method rests on: ${(rw / Math.max(1e-9, rh)).toFixed(2)}`)
console.log(rw > 1.6 * rh
  ? 'the anchor holds: a blastomere reads heterozygous far more often when she is homozygous for\n'
    + 'the other allele, so the vote carries information.'
  : 'THE ANCHOR DOES NOT HOLD on this material: a blastomere reads heterozygous at nearly the same\n'
    + 'rate whichever of the two she is, so no amount of voting can separate them. Mendel says the\n'
    + 'first rate should be about twice the second.')

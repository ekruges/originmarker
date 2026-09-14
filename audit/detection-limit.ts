/**
 * The whole-chromosome detection limit: at what effect size calls start, what they cost in false
 * calls, and what is missed below that.
 *
 * Positives and negatives are both known from outside the tool, per sample and per chromosome:
 *
 *   POSITIVE   chromosome 21 of a karyotype-confirmed trisomy 21 lymphoblast, GSE19247
 *   NEGATIVE   every autosome of a euploid reference line (Coriell 1463 and 1423) or of an adult
 *              gamete donor, who carries no autosomal whole-chromosome loss
 *
 * The other chromosomes of a trisomy 21 array are in neither set. Targeted-cut blastomeres are
 * excluded: "the experiment aimed here" is a prior about a group, not a stated answer for one
 * chromosome.
 *
 * Swept together: the two gates a whole-chromosome call must both clear, the z against the array's
 * own chromosome spread and the magnitude floor in log2 units, both from the centre parentage.ts
 * uses. It prints the curve and the cost at each point, including the shipped one, and picks no
 * threshold, so it neither passes nor fails.
 *
 * Run: OM_TRIOS=<dir> OM_GSE=<converted dir> \
 *        node --experimental-strip-types --max-old-space-size=3072 audit/detection-limit.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const nullMod = await import(`${W}intensityNull.ts`)
const parentage = await import(`${W}parentage.ts`)

const Z_SHIPPED = (nullMod as { Z_CHROMOSOME: number }).Z_CHROMOSOME
const FLOOR_SHIPPED = (parentage as { COPY_SHIFT_FLOOR: number }).COPY_SHIFT_FLOOR

/**
 * One array reduced to what the two gates actually read: a per-chromosome shift and z, computed
 * exactly as parentage.ts computes them, from the same centre and the same scale.
 */
function chromStats(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) return null
  const by = new Map<string, number[]>()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r || r.log2R === null || !/^\d+$/.test(r.chrom)) continue
    const a = by.get(r.chrom) ?? []
    a.push(r.log2R)
    by.set(r.chrom, a)
  }
  const med = (xs: number[]) => {
    const q = [...xs].sort((a, b) => a - b)
    return q.length ? q[q.length >> 1] : NaN
  }
  const chromMed = new Map<string, number>()
  for (const [c, xs] of by) if (xs.length >= 200) chromMed.set(c, med(xs))
  if (chromMed.size < 8) return null
  const scale = nullMod.nullScale([...chromMed.values()]) as { centre: number; scale: number }
  const out = new Map<string, { shift: number; z: number }>()
  for (const [c, m] of chromMed) {
    const z = nullMod.calibratedZ(m, scale as never) as number | undefined
    if (z === undefined) continue
    out.set(c, { shift: m - scale.centre, z })
  }
  return out
}

interface Obs { z: number; shift: number; positive: boolean; group: string }
const obs: Obs[] = []

const add = (path: string, group: string, positiveChrom: string | null) => {
  const s = chromStats(path)
  if (!s) return 0
  for (const [c, v] of s) {
    // On an array carrying a KNOWN event, only that chromosome has a stated answer. The rest are
    // unknown rather than clean, and admitting them as negatives is how a real second aneuploidy
    // becomes a false positive in the arithmetic.
    if (positiveChrom && c !== positiveChrom) continue
    obs.push({ z: Math.abs(v.z), shift: Math.abs(v.shift), positive: c === positiveChrom, group })
  }
  return 1
}

// ---------------------------------------------------------------- the sets
//
// ONLY WHERE THE ANSWER IS KNOWN PER SAMPLE AND PER CHROMOSOME. The first version of this file
// admitted targeted blastomeres to both sets and both were wrong. Chromosome 16 was counted as a
// positive on every array of a chr16 experiment, when a targeted cut succeeds in roughly a quarter
// of cells, so most of those "positives" carried nothing; and the other 21 chromosomes of those
// same embryos were counted as negatives, when widespread aneuploidy in edited embryos is the
// source study's own finding. The curve that came out read 11.9% sensitivity at 16.4 false calls
// per thousand, against 0 in 242 measured on genuinely known-normal donors. Both figures were
// artefacts of the sets.
//
// What survives is material whose answer is stated for the sample itself:
//   POSITIVE   chromosome 21 on a karyotype-confirmed trisomy 21 array. One per array, known.
//   NEGATIVE   every autosome of a euploid reference line or an adult gamete donor.
// Targeted blastomeres are excluded from BOTH sets and reported separately, because "the
// experiment aimed here" is a prior about the group and not a truth about the chromosome.
const GSE = process.env.OM_GSE
if (GSE && existsSync(join(GSE, 'truth.json'))) {
  interface T { gsm: string; src: string }
  const truth: T[] = JSON.parse(readFileSync(join(GSE, 'truth.json'), 'utf8'))
  let a = 0; let b = 0
  for (const t of truth.filter((x) => /family 1990/i.test(x.src)).slice(0, 60)) {
    if (existsSync(join(GSE, `${t.gsm}.probes`))) a += add(join(GSE, `${t.gsm}.probes`), 'trisomy21', '21')
  }
  for (const t of truth.filter((x) => /1463|1423/.test(x.src)).slice(0, 60)) {
    if (existsSync(join(GSE, `${t.gsm}.probes`))) b += add(join(GSE, `${t.gsm}.probes`), 'euploid', null)
  }
  console.log(`trisomy 21 arrays ${a} (chr21 positive, other 21 autosomes not scored either way)`)
  console.log(`euploid reference arrays ${b} (every autosome negative)`)
}

const DIR = process.env.OM_TRIOS
if (DIR && existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
  const head = lines[0].split(',')
  const gi = head.indexOf('gsm'); const ri = head.indexOf('role')
  let n = 0
  for (const l of lines.slice(1)) {
    const f = l.split(',')
    if (f[ri] !== 'parent') continue
    const p = join(DIR, `${f[gi]}.probes.gz`)
    if (existsSync(p)) n += add(p, 'adult donor', null)
  }
  console.log(`adult donor arrays ${n} (every autosome negative)`)
}

const pos = obs.filter((o) => o.positive)
const neg = obs.filter((o) => !o.positive)
console.log('')
console.log(`chromosome observations: ${obs.length.toLocaleString()}  `
  + `known positive ${pos.length}, known negative ${neg.length.toLocaleString()}`)
console.log('')

// ---------------------------------------------------------------- the sweep
const zs = [2, 3, 4, 5, Z_SHIPPED, 7, 9, 12]
const floors = [0.10, 0.20, 0.30, FLOOR_SHIPPED, 0.50, 0.60]
console.log('=== BOTH GATES SWEPT TOGETHER. A call needs |z| over the row AND |shift| over the column.')
console.log('   Each cell: sensitivity on known positives / false calls per 1,000 known-negative chromosomes')
let hdr = '   z\\floor  '
for (const fl of floors) hdr += `${fl === FLOOR_SHIPPED ? `[${fl.toFixed(2)}]` : ` ${fl.toFixed(2)} `}`.padStart(14)
console.log(hdr)
for (const z of zs) {
  let row = `  ${(z === Z_SHIPPED ? `[${z}]` : `${z}`).padEnd(9)}`
  for (const fl of floors) {
    const tp = pos.filter((o) => o.z > z && o.shift >= fl).length
    const fp = neg.filter((o) => o.z > z && o.shift >= fl).length
    const sens = pos.length ? tp / pos.length : NaN
    const per1k = neg.length ? (1000 * fp) / neg.length : NaN
    row += `${(sens * 100).toFixed(0)}% / ${per1k.toFixed(1)}`.padStart(14)
  }
  console.log(row)
}
console.log('')
console.log(`  [square brackets] mark the shipped operating point: z > ${Z_SHIPPED}, |shift| >= ${FLOOR_SHIPPED}`)

// ---------------------------------------------------------------- where the shipped point sits
const tpS = pos.filter((o) => o.z > Z_SHIPPED && o.shift >= FLOOR_SHIPPED).length
const fpS = neg.filter((o) => o.z > Z_SHIPPED && o.shift >= FLOOR_SHIPPED).length
console.log('')
console.log('=== THE SHIPPED OPERATING POINT')
console.log(`  sensitivity   ${tpS}/${pos.length} = ${(100 * tpS / Math.max(pos.length, 1)).toFixed(1)}%`)
console.log(`  false calls   ${fpS} over ${neg.length.toLocaleString()} known-negative chromosomes `
  + `= ${(1000 * fpS / Math.max(neg.length, 1)).toFixed(2)} per 1,000`)

// ---------------------------------------------------------------- the limit itself
console.log('')
console.log('=== THE DETECTION LIMIT: what effect size does a positive need to be called?')
const q = (xs: number[], p: number) => {
  const s = [...xs].sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
const called = pos.filter((o) => o.z > Z_SHIPPED && o.shift >= FLOOR_SHIPPED)
const missed = pos.filter((o) => !(o.z > Z_SHIPPED && o.shift >= FLOOR_SHIPPED))
console.log(`  positives CALLED  n=${called.length}  |shift| p50 ${q(called.map((o) => o.shift), 0.5).toFixed(3)}  `
  + `|z| p50 ${q(called.map((o) => o.z), 0.5).toFixed(1)}`)
console.log(`  positives MISSED  n=${missed.length}  |shift| p50 ${q(missed.map((o) => o.shift), 0.5).toFixed(3)}  `
  + `|z| p50 ${q(missed.map((o) => o.z), 0.5).toFixed(1)}`)
console.log('')
console.log('  A missed positive whose shift sits near the floor is a threshold cost and is')
console.log('  recoverable. One whose shift is near zero was never in the data and no threshold')
console.log('  recovers it. Those are different failures and a paper has to separate them.')
const recoverable = missed.filter((o) => o.shift >= 0.15).length
console.log(`  missed WITH a real shift (>= 0.15): ${recoverable}/${missed.length}`)
console.log(`  missed with essentially no signal:  ${missed.length - recoverable}/${missed.length}`)

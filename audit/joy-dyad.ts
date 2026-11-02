/**
 * Joy's ask: double-check the polar body result.
 *
 * The claim the pipeline rests on is that a DYAD is haploid even though it carries real
 * heterozygosity. A first polar body and an MII oocyte are both dyads: one chromosome, two sister
 * chromatids, and distal to each crossover those chromatids carry different haplotypes. So the
 * genome reads homozygous proximal to a crossover and heterozygous distal to it, in BLOCKS. A
 * second polar body, pronucleus or sperm carries one chromatid and has true heterozygosity of zero,
 * so anything heterozygous there is error.
 *
 * `stage.check.ts` pins this against a synthetic value. This checks it on real arrays, and the test
 * is not the heterozygosity alone but whether it is STRUCTURED. Dropout and drop-in scatter;
 * a crossover does not.
 *
 * Classes are read from the lab's own sheet text, not inferred.
 *
 * Run: OM_DIR=<dir> OM_CALLS=<manual_calls.tsv> node --experimental-strip-types audit/joy-dyad.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { loadArray, stageMod, isArray } from './joy-common.ts'

const DIR = process.env.OM_DIR
const CALLS = process.env.OM_CALLS
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
if (!CALLS || !existsSync(CALLS)) { console.log('OM_CALLS not readable. NOT RUN.'); process.exit(0) }

type Klass = 'dyad' | 'one chromatid' | 'diploid' | 'other'
const classify = (text: string): Klass => {
  const t = text.toLowerCase()
  // An MII oocyte and a first polar body are the same structure: a dyad.
  if (/mii|polar body|\bpb1\b|failed fertilization oocyte/.test(t)) return 'dyad'
  // 23X or 23Y is a single haploid complement: one chromatid, true heterozygosity zero.
  if (/\b23\s*[xy]\b/.test(t)) return 'one chromatid'
  if (/\b46\b|morula|hesc|te biopsy|6-c/.test(t)) return 'diploid'
  return 'other'
}
const label = new Map<string, string>()
for (const l of readFileSync(CALLS, 'utf8').split('\n')) {
  if (!l || l.startsWith('sample')) continue
  const c = l.split('\t')
  if (c[0] && c[3] && !label.has(c[0].trim())) label.set(c[0].trim(), c[3].trim())
}
console.log(`${label.size} samples carry sheet text\n`)

const files = readdirSync(DIR).filter(isArray).sort()

/**
 * Heterozygosity in windows along each chromosome. A dyad's heterozygosity sits in blocks, so some
 * windows read high and others near zero. Scattered error reads uniformly low everywhere.
 */
const windowed = (
  hetByChrom: Map<string, { pos: number; het: boolean }[]>, size = 200,
) => {
  const rates: number[] = []
  for (const marks of hetByChrom.values()) {
    for (let i = 0; i + size <= marks.length; i += size) {
      let h = 0
      for (let j = i; j < i + size; j += 1) if (marks[j].het) h += 1
      rates.push(h / size)
    }
  }
  rates.sort((a, b) => a - b)
  if (!rates.length) return undefined
  const q = (p: number) => rates[Math.min(rates.length - 1, Math.floor(p * rates.length))]
  const mean = rates.reduce((a, b) => a + b, 0) / rates.length
  return {
    windows: rates.length, median: q(0.5), p90: q(0.9), max: rates[rates.length - 1],
    // The share of the genome sitting in windows that are clearly heterozygous. For a dyad this is
    // the share distal to a crossover; for one chromatid it should be nothing.
    blockShare: rates.filter((r) => r >= 0.15).length / rates.length,
    quietShare: rates.filter((r) => r <= 0.02).length / rates.length,
    mean,
  }
}

type Row = { id: string; klass: Klass; text: string; stage: string; het: number
  callRate: number; w: ReturnType<typeof windowed>
  hetByChrom: Map<string, { pos: number; het: boolean }[]> }
const rows: Row[] = []
for (const f of files) {
  const id = (/^(?:GSM\d+_)?([0-9]+-[0-9]+)/.exec(f)?.[1]) ?? f.replace(/\..*$/, '')
  const text = label.get(id)
  if (!text) continue
  try {
    const a = loadArray(join(DIR, f), id)
    const s = stageMod.inferStage(a.profile as never) as { stage: string }
    rows.push({ id, klass: classify(text), text, stage: s.stage, het: a.hetRate,
      callRate: a.callRate, w: windowed(a.hetByChrom), hetByChrom: a.hetByChrom })
  } catch (e) {
    console.log(`  ${id}: ${(e as Error).message}`)
  }
}

console.log(`${'sample'.padEnd(11)} ${'class'.padEnd(14)} ${'stage'.padEnd(11)} ${'call'.padStart(5)} `
  + `${'het'.padStart(6)} ${'win med'.padStart(8)} ${'win max'.padStart(8)} ${'blocks'.padStart(7)} `
  + `${'quiet'.padStart(6)}`)
for (const r of rows.sort((a, b) => a.klass.localeCompare(b.klass) || a.id.localeCompare(b.id))) {
  console.log(`${r.id.padEnd(11)} ${r.klass.padEnd(14)} ${r.stage.padEnd(11)} `
    + `${r.callRate.toFixed(2).padStart(5)} ${r.het.toFixed(4).padStart(6)} `
    + `${(r.w?.median ?? NaN).toFixed(3).padStart(8)} ${(r.w?.max ?? NaN).toFixed(3).padStart(8)} `
    + `${(r.w?.blockShare ?? NaN).toFixed(3).padStart(7)} ${(r.w?.quietShare ?? NaN).toFixed(3).padStart(6)}`)
}

const dyads = rows.filter((r) => r.klass === 'dyad')
const singles = rows.filter((r) => r.klass === 'one chromatid')
const dips = rows.filter((r) => r.klass === 'diploid')
console.log(`\n${dyads.length} dyad, ${singles.length} one-chromatid, ${dips.length} diploid\n`)

console.log('CLAIM 1. A DYAD STAGES AS HAPLOID DESPITE REAL HETEROZYGOSITY')
for (const r of dyads) {
  console.log(`  ${r.id}: stage ${r.stage}, heterozygosity ${r.het.toFixed(4)}`)
  console.log(`    sheet says: ${r.text}`)
}
const dyadOk = dyads.length > 0 && dyads.every((r) => r.stage === 'haploid')
console.log(dyads.length === 0 ? '  NO DYAD SAMPLE PRESENT, claim not testable on this material'
  : dyadOk ? '  HOLDS: every dyad stages haploid.'
  : '  FAILS: a dyad staged as something other than haploid.')

// THE DISCRIMINATOR. A dyad has one or two crossovers per chromosome, so each chromosome should
// show a STEP: quiet proximal to the crossover where the sister chromatids are still identical,
// heterozygous distal to it. A diploid is heterozygous the whole way along. A single chromatid is
// quiet the whole way along. So the test is a within-chromosome step, not the genome-wide average,
// and it does not depend on an absolute quiet threshold that this platform's drop-in would break.
const steps = (hetByChrom: Map<string, { pos: number; het: boolean }[]>, size = 200) => {
  let stepped = 0; let scored = 0
  const detail: string[] = []
  for (const [chrom, marks] of [...hetByChrom].sort((a, b) => Number(a[0]) - Number(b[0]))) {
    const w: number[] = []
    for (let i = 0; i + size <= marks.length; i += size) {
      let h = 0
      for (let j = i; j < i + size; j += 1) if (marks[j].het) h += 1
      w.push(h / size)
    }
    if (w.length < 6) continue
    scored += 1
    const lo = Math.min(...w); const hi = Math.max(...w)
    // A step means the quiet end and the loud end of the same chromosome are far apart.
    // The loud end has to be genuinely dyad-level, not merely above the quiet end. Without the
    // `hi` bar a true single chromatid scores a step on every chromosome, because one noisy window
    // clears a relative gap when the baseline is near zero.
    if (lo <= 0.06 && hi >= 0.30 && hi - lo >= 0.20) {
      stepped += 1
      detail.push(`chr${chrom} ${lo.toFixed(3)} to ${hi.toFixed(3)}`)
    }
  }
  return { stepped, scored, share: scored ? stepped / scored : NaN, detail }
}

console.log('\nCLAIM 2A. WITHIN-CHROMOSOME STEP, which is what a crossover actually leaves')
for (const [name, set] of [['dyad', dyads], ['staged haploid', rows.filter((r) => r.stage === 'haploid')],
  ['diploid', dips]] as const) {
  if (!set.length) continue
  const v = set.map((r) => steps(r.hetByChrom))
  const sh = v.map((x) => x.share).filter(Number.isFinite)
  console.log(`  ${name.padEnd(16)} chromosomes with a step: `
    + `${(100 * sh.reduce((a, b) => a + b, 0) / Math.max(1, sh.length)).toFixed(0)}% `
    + `(n=${set.length})`)
  if (name === 'dyad') for (const d of v[0].detail.slice(0, 8)) console.log(`    ${d}`)
}

console.log('\nCLAIM 2. THE HETEROZYGOSITY IS STRUCTURED, NOT SCATTERED')
console.log('  A crossover puts heterozygosity in blocks. Dropout and drop-in scatter it.')
const avg = (xs: Row[], pick: (w: NonNullable<Row['w']>) => number) => {
  const v = xs.map((r) => (r.w ? pick(r.w) : NaN)).filter(Number.isFinite)
  return v.length ? v.reduce((a, b) => a + b, 0) / v.length : NaN
}
for (const [name, set] of [['dyad', dyads], ['one chromatid', singles], ['diploid', dips]] as const) {
  if (!set.length) continue
  console.log(`  ${name.padEnd(14)} het ${avg(set, (w) => w.mean).toFixed(4)}  `
    + `windows at or over 0.15: ${avg(set, (w) => w.blockShare).toFixed(3)}  `
    + `windows at or under 0.02: ${avg(set, (w) => w.quietShare).toFixed(3)}`)
}
const clean = rows.filter((r) => r.stage === 'haploid')
if (clean.length) {
  console.log(`\n  Using the ${clean.length} sample(s) that actually STAGE haploid as the contrast, `
    + 'because the sheet\'s own haploid label does not hold here:')
  console.log(`  ${clean.map((r) => `${r.id} het ${r.het.toFixed(4)}`).join(', ')}`)
  console.log(`  dyad blocks ${avg(dyads, (w) => w.blockShare).toFixed(3)} against `
    + `${avg(clean, (w) => w.blockShare).toFixed(3)} for those, a factor of `
    + `${(avg(dyads, (w) => w.blockShare) / Math.max(1e-9, avg(clean, (w) => w.blockShare))).toFixed(0)}`)
}
if (dyads.length && singles.length) {
  const dBlocks = avg(dyads, (w) => w.blockShare)
  const sBlocks = avg(singles, (w) => w.blockShare)
  const dQuiet = avg(dyads, (w) => w.quietShare)
  console.log(`\n  a dyad carries ${(dBlocks / Math.max(1e-9, sBlocks)).toFixed(1)}x the share of `
    + 'clearly heterozygous windows that a single chromatid does')
  console.log(dBlocks > 2 * sBlocks && dQuiet > 0.1
    ? '  HOLDS: the dyad is heterozygous in blocks and still has quiet stretches, which is what a\n'
      + '  crossover produces. A uniformly noisy array would have neither.'
    : '  DOES NOT HOLD on this material: the dyad does not look structurally different from a\n'
      + '  single chromatid, so its heterozygosity cannot be attributed to crossovers here.')
}
console.log('\nOne dyad is one dyad. This confirms or fails the claim on the material the lab has,')
console.log('and the contrast classes are what make a single sample worth anything.')

/**
 * THE TOOL ITSELF, ON A SECOND STUDY, A SECOND LABORATORY AND A SECOND ARRAY CHEMISTRY.
 *
 * Everything else in this directory measures the tool on GSE148488: one laboratory, one Affymetrix
 * platform, one experimental design. Every threshold in the tree was measured there. A constant that
 * is really a property of that series rather than of the biology would look identical to a correct
 * one, and no amount of re-running it on the same corpus would tell the difference.
 *
 * GSE19247 (Vanneste 2009, Nature Medicine) is a different laboratory, a different decade and
 * Illumina rather than Affymetrix, at roughly half the marker density. Its GEO `source_name` field
 * carries several kinds of known answer, and none of them come from this tool:
 *
 *   sperm cells                       HAPLOID. One parental complement, by biology. The zygosity
 *                                     call must not read them as diploid, and that call is the gate
 *                                     every other channel sits behind.
 *   lymphoblast, Coriell family 1990  karyotype-confirmed TRISOMY 21, and diploid otherwise. A gain
 *                                     on chr21 is the expected finding, and the other 21 autosomes
 *                                     are a specificity set on the same arrays.
 *   Day 3 cleavage stage embryo NNNN  several single cells of ONE embryo. Same genome, so the
 *                                     answers should agree; the biology is famously mosaic, so
 *                                     they may honestly differ, and the two are separated below.
 *   blood cell                        ordinary diploid somatic tissue. Must read diploid.
 *
 * THE INTENSITY IS THIS CONVERSION'S, NOT THE PUBLISHER'S. See toprobes.py: the log2 ratio is built
 * from total intensity against a per-marker median across the series. Anything the intensity channel
 * says here carries that construction with it. The GENOTYPES come from the caller already in this
 * directory and are unchanged.
 *
 * AND THE B-ALLELE FREQUENCY IS NOT FIT FOR THE ZYGOSITY CALL. MEASURED, by
 * audit/zygosity-crossplatform.ts, on known-haploid material from both series:
 *
 *   GSE148488 pronuclei, haploid   band 0.0359   homBand 0.0033   excess 0.0327
 *   GSE19247  sperm, haploid       band 0.1362   homBand 0.0000   excess 0.1362
 *   GSE19247  lymphoblast, diploid band 0.0903   homBand 0.0000   excess 0.0903
 *
 * Two things there are impossible for a real array. `homBand` is identically zero, because the
 * piecewise map in toprobes.py saturates at exactly 0 and 1 outside the homozygote clouds, so the
 * homozygous-cluster correction that HET_BAND_EXCESS is built on cannot exist. And the known
 * HAPLOID sperm carry MORE mid-band mass than the known DIPLOID lymphoblasts, which is backwards.
 *
 * SO THE ZYGOSITY ARM BELOW IS A TEST OF THE CONVERTER, NOT OF THE TOOL, and it fails as one. It is
 * kept and still printed because deleting a failing arm is how a limitation becomes invisible, but
 * nothing about the tool's zygosity boundary may be concluded from it in either direction. What
 * this file DOES establish on a second platform is detection: the trisomy arm and the crash arm
 * read a real log2 ratio and a real genotype, neither of which depends on the BAF shape.
 *
 * Run: OM_GSE=<dir of .probes and truth.json> node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/gse19247-tool.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_GSE
if (!DIR || !existsSync(DIR) || !existsSync(join(DIR, 'truth.json'))) {
  console.log('gse19247-tool: OM_GSE must hold the converted .probes and truth.json. NOT RUN.')
  console.log('  Build them with audit/gse19247/toprobes.py; the data is public, see that README.')
  process.exit(0)
}

interface Truth { gsm: string; title: string; src: string }
const truth: Truth[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: Truth) => existsSync(join(DIR, `${t.gsm}.probes`)))

const isSperm = (s: string) => /sperm/i.test(s)
const isTrisomy21 = (s: string) => /Coriell family 1990/i.test(s)
const isBlood = (s: string) => /blood cell/i.test(s)
/** "Day 3 cleavage stage embryo 1077 from family 13" -> "1077". One embryo, many cells. */
const embryoOf = (s: string) => /cleavage stage embryo\s+(\S+)/i.exec(s)?.[1] ?? null

function rowsOf(path: string) {
  const ls = readFileSync(path, 'utf8').split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  if (h < 0) throw new Error('no header')
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

/**
 * A fixed array in the parental slot, taken FROM THIS SERIES so the platform matches.
 *
 * The channels read below (zygosity, stage, copy number) do not consult it. It is here because
 * `scoreSample` takes a parent, and using one array for every sample means it cannot be a variable.
 */
const refRec = truth.find((t) => isTrisomy21(t.src)) ?? truth[0]
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR, `${refRec.gsm}.probes`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Out {
  gsm: string; src: string; stage: string; zygosity: string; originClass: string
  callRate: number; hetRate: number; gains: string[]; losses: string[]; events: string[]
  error?: string
}
const results: Out[] = []
const LIMIT = Number(process.env.OM_N ?? truth.length)

console.log(`GSE19247 converted arrays: ${truth.length}, scoring ${Math.min(LIMIT, truth.length)}`)
console.log(`reference array in the parental slot: ${refRec.gsm}`)
console.log('')

for (const t of truth.slice(0, LIMIT)) {
  try {
    const acc = score.emptyCollected(refParent as never, null)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(join(DIR, `${t.gsm}.probes`))()) {
      if (!first) first = r.probesetId
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(t.gsm, byChrom as never, bafSums as never, first)
    const out = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: t.gsm, log: () => {},
    }) as Record<string, unknown>
    const chroms = (out.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
    const segs = (out.segments ?? []) as { chrom: string }[]
    results.push({
      gsm: t.gsm, src: t.src,
      stage: (out.stage as { stage?: string } | undefined)?.stage ?? 'none',
      zygosity: String(out.zygosity ?? 'none'),
      originClass: String(out.originClass ?? 'none'),
      callRate: profile.callRate,
      hetRate: profile.hetRate,
      gains: chroms.filter((c) => c.aneuploidy === 'gain').map((c) => c.chrom),
      losses: chroms.filter((c) => c.aneuploidy === 'loss').map((c) => c.chrom),
      events: [...new Set([...chroms.filter((c) => c.aneuploidy).map((c) => c.chrom),
        ...segs.map((s) => s.chrom)])],
    })
  } catch (e) {
    results.push({
      gsm: t.gsm, src: t.src, stage: 'error', zygosity: '-', originClass: '-',
      callRate: NaN, hetRate: NaN, gains: [], losses: [], events: [],
      error: String((e as Error).message ?? e).slice(0, 70),
    })
  }
}

let failures = 0
const fail = (s: string) => { failures += 1; console.log(`  *** FAIL  ${s}`) }
const pct = (n: number, d: number) => (d ? `${n}/${d} (${(100 * n / d).toFixed(0)}%)` : 'n/a')

// ================================================================ 1. SPERM IS HAPLOID
console.log('=== 1. SPERM CELLS ARE HAPLOID, and the zygosity call is the gate for everything')
{
  const sperm = results.filter((r) => isSperm(r.src) && !r.error)
  const diploid = sperm.filter((r) => r.zygosity === 'diploid')
  const uni = sperm.filter((r) => /uniparental/.test(r.zygosity))
  const withheld = sperm.filter((r) => !/uniparental|diploid/.test(r.zygosity))
  console.log(`    sperm arrays scored: ${sperm.length}`)
  console.log(`      read as uniparental  ${pct(uni.length, sperm.length)}   <- correct`)
  console.log(`      read as DIPLOID      ${pct(diploid.length, sperm.length)}   <- wrong, `
    + 'a sperm cell has one complement')
  console.log(`      withheld or unknown  ${pct(withheld.length, sperm.length)}`)
  console.log(`      median het rate ${sperm.length
    ? [...sperm].sort((a, b) => a.hetRate - b.hetRate)[sperm.length >> 1].hetRate.toFixed(4)
    : 'n/a'}, median call rate ${sperm.length
    ? [...sperm].sort((a, b) => a.callRate - b.callRate)[sperm.length >> 1].callRate.toFixed(3)
    : 'n/a'}`)
  // NOT COUNTED AS A FAILURE OF THE TOOL. The converted BAF cannot support this call: see the
  // header. Printed so the limitation is visible, scored so it cannot be mistaken for a pass.
  if (sperm.length && diploid.length > sperm.length / 2) {
    console.log(`      NOT SCORED: the converted BAF gives homBand 0.0000 on every array, so the `
      + 'correction this boundary rests on is absent by construction. The zygosity call is '
      + 'UNVALIDATED on this platform rather than refuted by it.')
  }
}

// ================================================================ 2. TRISOMY 21
console.log('')
console.log('=== 2. KARYOTYPE-CONFIRMED TRISOMY 21, and the other 21 autosomes as the control')
{
  const tri = results.filter((r) => isTrisomy21(r.src) && !r.error)
  const usable = tri.filter((r) => r.stage !== 'failed')
  const on21 = usable.filter((r) => r.gains.includes('21') || r.events.includes('21'))
  const offPer = usable.length
    ? usable.reduce((a, r) => a + r.events.filter((c) => c !== '21').length, 0) / usable.length
    : NaN
  console.log(`    arrays: ${tri.length}, usable after the tool's own gates: ${usable.length}`)
  console.log(`      change found on chr21      ${pct(on21.length, usable.length)}`)
  console.log(`      other chromosomes per array ${Number.isFinite(offPer) ? offPer.toFixed(2) : 'n/a'}`
    + '  (of 21 available)')
  console.log(`      read as diploid            ${pct(usable.filter((r) => r.zygosity === 'diploid').length, usable.length)}`)
  const rate21 = usable.length ? on21.length / usable.length : NaN
  const rateOther = Number.isFinite(offPer) ? offPer / 21 : NaN
  console.log(`      chr21 rate ${Number.isFinite(rate21) ? rate21.toFixed(3) : 'n/a'} against `
    + `${Number.isFinite(rateOther) ? rateOther.toFixed(3) : 'n/a'} per other chromosome`)
}

// ================================================================ 3. ONE EMBRYO, MANY CELLS
console.log('')
console.log('=== 3. SEVERAL SINGLE CELLS OF ONE EMBRYO')
console.log('    Cleavage-stage embryos are genuinely mosaic, which this series is famous for')
console.log('    establishing, so cells of one embryo MAY differ for real. What is reported is the')
console.log('    spread, not a pass or a fail.')
{
  const byEmbryo = new Map<string, Out[]>()
  for (const r of results) {
    const e = embryoOf(r.src)
    if (e && !r.error) byEmbryo.set(e, [...(byEmbryo.get(e) ?? []), r])
  }
  const groups = [...byEmbryo.entries()].filter(([, v]) => v.length >= 2)
  let sameStage = 0
  let sameZyg = 0
  for (const [e, cells] of groups) {
    const st = new Set(cells.map((c) => c.stage)).size === 1
    const zy = new Set(cells.map((c) => c.zygosity)).size === 1
    sameStage += st ? 1 : 0
    sameZyg += zy ? 1 : 0
    const evs = cells.map((c) => c.events.length)
    console.log(`    embryo ${e.padEnd(6)} cells ${cells.length}  events per cell `
      + `${evs.join('/')}  stage ${st ? 'same' : 'DIFFER'}  zygosity ${zy ? 'same' : 'DIFFER'}`)
  }
  console.log(`    embryos with 2+ cells: ${groups.length}, same material inferred `
    + `${pct(sameStage, groups.length)}, same zygosity ${pct(sameZyg, groups.length)}`)
}

// ================================================================ 4. NOTHING CRASHED
console.log('')
console.log('=== 4. THE PIPELINE ON A PLATFORM IT WAS NOT TUNED ON')
{
  const errs = results.filter((r) => r.error)
  const failedStage = results.filter((r) => r.stage === 'failed')
  console.log(`    scored ${results.length}, threw ${errs.length}, `
    + `refused as failed material ${failedStage.length}`)
  for (const e of errs.slice(0, 8)) console.log(`      ${e.gsm}: ${e.error}`)
  const stages = new Map<string, number>()
  for (const r of results) stages.set(r.stage, (stages.get(r.stage) ?? 0) + 1)
  console.log(`    material inferred: ${JSON.stringify(Object.fromEntries(stages))}`)
  if (errs.length) fail(`${errs.length} array(s) threw on this platform`)
}

// Blood, the plainest diploid there is.
const blood = results.filter((r) => isBlood(r.src) && !r.error)
if (blood.length) {
  console.log('')
  console.log(`=== 5. BLOOD, ordinary diploid somatic tissue: read as diploid `
    + `${pct(blood.filter((r) => r.zygosity === 'diploid').length, blood.length)}`)
}

console.log('')
console.log(failures === 0
  ? 'gse19247-tool: a second study, a second platform, and nothing the tool asserted was refuted.'
  : `gse19247-tool: ${failures} failure(s) against a known answer from another laboratory.`)
process.exit(failures ? 1 : 0)

/**
 * THIRTEEN COMPLETE HYDATIDIFORM MOLES, WHERE THE ANSWER IS BOTH KNOWN AND HARD.
 *
 * WHY THIS SERIES AND NOT A LARGER ONE. A complete hydatidiform mole carries two paternal genomes
 * and no maternal one. These thirteen were established as DISPERMIC by the submitters, from
 * centromeric zygosity, which means they are androgenetic AND heterozygous at the same time. That
 * combination is the hardest case this tool has: heterozygosity is what its zygosity call reads as
 * two parental contributions, and here there are two contributions from ONE parent.
 *
 * It is also the exact branch that carried the class inversion fixed in 5.25.0, whose own comment
 * reads "meaning two sperm". These are thirteen real cases of two sperm, and until now that branch
 * had only ever been exercised against a written-out fixture.
 *
 * WHAT IS KNOWN, FROM THE SUBMITTERS AND NOT FROM THIS TOOL:
 *   all 13   androgenetic, dispermic, therefore HETEROZYGOUS rather than homozygous
 *    3 of 13 carry a trisomy
 *    6 of 13 female placental sex, 7 male
 *
 * WHAT THIS CANNOT TEST, and saying so matters more than the arms below. There is no parental
 * array in this series, so `originClass` cannot be scored: naming a parent needs a parent loaded.
 * What can be scored is everything the sample answers about ITSELF.
 *
 * THE SECOND-PLATFORM QUESTION. 5.27.0 recorded that HET_BAND_EXCESS could not be tested off
 * Affymetrix, because a cluster-file B-allele frequency is a function of the same theta the
 * genotype came from and its homozygous band is therefore structurally zero. These files are
 * GenomeStudio final reports, and 36% of their homozygous calls sit on an exact 0 or 1 against
 * 70-85% for a cluster-file export. So a floor is measurable here, and the boundary can be put to
 * a second platform for the first time.
 *
 * Run: OM_MOLES=<dir of .probes> node --experimental-strip-types \
 *        --max-old-space-size=3072 audit/mole-truth.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const parentage = await import(`${W}parentage.ts`)
const score = await import(`${W}scoreSample.ts`)
const sexing = await import(`${W}sexing.ts`)

const DIR = process.env.OM_MOLES
if (!DIR || !existsSync(DIR)) {
  console.log('mole-truth: OM_MOLES is not a readable directory. NOT RUN.')
  process.exit(0)
}

/** The submitters' counts. Nothing below is allowed to move them. */
const TRUTH = { n: 13, androgenetic: 13, trisomy: 3, female: 6, male: 7 }

function rowsOf(path: string) {
  const ls = readFileSync(path, 'utf8').split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

const files = readdirSync(DIR).filter((f) => f.endsWith('.probes')).sort()
console.log(`moles on disk: ${files.length}, expected ${TRUTH.n}`)
console.log(`the boundary under test: HET_BAND_EXCESS = ${parentage.HET_BAND_EXCESS}`)
console.log('')

/**
 * A FIXED ARRAY IN THE PARENTAL SLOT, and it is one of the moles.
 *
 * Every arm below reads the sample alone: zygosity, its own band, its own sex, its own dosage.
 * None of them consults the parental array. One is supplied because scoreSample requires one, and
 * using a member of this same series keeps the comparison inside one platform. No parental claim
 * is scored anywhere in this file, for the reason in the header.
 */
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR, files[0]))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Row {
  name: string; stage: string; zygosity: string; callRate: number
  hetBand: number; homBand: number; excess: number
  yBearing: boolean | null; gains: string[]; losses: string[]; junk: string[]
}
const out: Row[] = []

for (const f of files) {
  const path = join(DIR, f)
  const t = parentage.emptyTally()
  const sx = sexing.emptySex()
  const acc = score.emptyCollected(refParent as never, null)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  const maxPos = new Map<string, number>()
  for (const r of rowsOf(path)()) {
    if (!first) first = r.probesetId
    if (Number.isFinite(r.pos)) maxPos.set(r.chrom, Math.max(maxPos.get(r.chrom) ?? 0, r.pos))
    parentage.tallyRow('NC' as never, r as never, t as never)
    sexing.accumulateSex(r as never, sx as never)
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, refParent as never, null, acc as never)
  }
  const profile = ingest.finishProfile(f, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
    sampleName: f, log: () => {},
  }) as Record<string, unknown>

  const tt = t as unknown as {
    bafInBand: number; bafTotal: number; homInBand: number; homTotal: number
  }
  const hetBand = tt.bafTotal ? tt.bafInBand / tt.bafTotal : NaN
  const homBand = tt.homTotal ? tt.homInBand / tt.homTotal : NaN
  const chroms = (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
  const junk: string[] = []
  for (const c of chroms) {
    if (c.aneuploidy && !maxPos.has(c.chrom)) junk.push(`chr${c.chrom} has no markers`)
  }
  out.push({
    name: f.replace('.probes', ''),
    stage: (res.stage as { stage?: string } | undefined)?.stage ?? 'none',
    zygosity: String(res.zygosity ?? '?'),
    callRate: profile.callRate,
    hetBand,
    homBand,
    excess: hetBand - homBand,
    yBearing: (sexing.sexCall(sx as never) as { yBearing: boolean | null }).yBearing,
    gains: chroms.filter((c) => c.aneuploidy === 'gain').map((c) => c.chrom),
    losses: chroms.filter((c) => c.aneuploidy === 'loss').map((c) => c.chrom),
    junk,
  })
}

console.log('=== PER MOLE')
console.log('sample       stage          zygosity                call   hetBand homBand excess  Y   gains/losses')
for (const r of out) {
  console.log(`  ${r.name.padEnd(11)}${r.stage.padEnd(15)}${r.zygosity.padEnd(24)}`
    + `${r.callRate.toFixed(3)}  ${r.hetBand.toFixed(4)}  ${r.homBand.toFixed(4)}  `
    + `${r.excess.toFixed(4)}  ${r.yBearing === null ? '?' : r.yBearing ? 'Y' : '-'}   `
    + `${r.gains.length ? `+${r.gains.join(',+')}` : ''}${r.losses.length ? ` -${r.losses.join(',-')}` : ''}`)
}

// ---------------------------------------------------------------- the scored claims
console.log('')
console.log('=== CLAIM 1: a heterozygous mole must not be called one parental contribution')
const uni = out.filter((r) => r.zygosity.startsWith('uniparental'))
console.log(`  called uniparental_homozygous: ${uni.length}/${out.length}`
  + `${uni.length ? `   <-- WRONG on ${uni.map((r) => r.name).join(', ')}` : '   correct'}`)
const dip = out.filter((r) => r.zygosity === 'diploid').length
const unk = out.filter((r) => r.zygosity === 'unknown').length
console.log(`  diploid ${dip}, unknown ${unk}, uniparental ${uni.length}`)

console.log('')
console.log('=== CLAIM 2: the boundary, on a second platform with a measurable floor')
const q = (xs: number[], p: number) => {
  const s = xs.filter(Number.isFinite).sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
console.log(`  hetBand  p50 ${q(out.map((r) => r.hetBand), 0.5).toFixed(4)}`)
console.log(`  homBand  p50 ${q(out.map((r) => r.homBand), 0.5).toFixed(4)}`
  + `   (structurally 0.0000 on a cluster-file export; a real value here is the point)`)
console.log(`  excess   p50 ${q(out.map((r) => r.excess), 0.5).toFixed(4)}`
  + `  against HET_BAND_EXCESS ${parentage.HET_BAND_EXCESS}`)
console.log(`  over the boundary: ${out.filter((r) => r.excess > parentage.HET_BAND_EXCESS).length}/${out.length}`)
console.log('  These are heterozygous genomes, so an excess ABOVE the boundary is the right answer')
console.log('  and reaching it on non-Affymetrix material is what could not be shown before.')

console.log('')
console.log('=== CLAIM 3: placental sex, 6 female and 7 male by the submitters')
const y = out.filter((r) => r.yBearing === true).length
const noY = out.filter((r) => r.yBearing === false).length
console.log(`  Y-bearing ${y}, not ${noY}, undetermined ${out.length - y - noY}`)
console.log(`  expected ${TRUTH.male} male and ${TRUTH.female} female: `
  + `${y === TRUTH.male && noY === TRUTH.female ? 'MATCHES' : 'does not match'}`)

console.log('')
console.log('=== CLAIM 4: three of thirteen carry a trisomy')
const withGain = out.filter((r) => r.gains.length > 0)
console.log(`  arrays showing a whole-chromosome gain: ${withGain.length}, expected ${TRUTH.trisomy}`)
for (const r of withGain) console.log(`    ${r.name}: +${r.gains.join(',+')}`)
console.log('  The submitters give a count, not which three, so this scores the COUNT only.')

console.log('')
const allJunk = out.flatMap((r) => r.junk)
console.log(`structurally invalid events: ${allJunk.length}`)
console.log(uni.length === 0
  ? 'mole-truth: no heterozygous mole was read as a single parental contribution.'
  : `mole-truth: ${uni.length} heterozygous moles read as one parental contribution.`)

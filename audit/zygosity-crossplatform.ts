/**
 * DOES THE ZYGOSITY BOUNDARY TRANSFER TO A SECOND PLATFORM, OR IS IT A PROPERTY OF ONE SERIES?
 *
 * THE OBSERVATION THAT FORCED THIS. Run on GSE19247, the tool read 20 of 23 SPERM CELLS as diploid.
 * A sperm cell carries one parental complement; there is no argument about the truth here. On
 * GSE148488 the same call gets every haploid pronucleus right. So either the boundary does not
 * transfer, or the conversion that put GSE19247 into the tool's format does not produce the BAF
 * shape the boundary was measured on. Those are opposite conclusions and the difference matters:
 * one is a limit of the tool that belongs in any paper, the other is a bug in a harness.
 *
 * WHAT SEPARATES THEM. `HET_BAND_EXCESS` compares the heterozygous BAF band against the same
 * array's OWN homozygous clusters, and the whole reason it is written that way is that the LEVEL of
 * the band is amplification quality rather than ploidy.
 *
 * THE ANSWER TURNED OUT TO BE NEITHER OF THOSE TWO, and this header records it so the run is read
 * for what it found rather than for what it set out to ask. The converted BAF is correctly scaled:
 * `audit/baf-shape.ts` shows AB calls sitting at 0.4886 to 0.5036 on every material, which is what
 * a B-allele frequency means. What the converted file cannot supply is the FLOOR, because a
 * cluster-file BAF is a function of the same theta the genotype came from, so a marker called
 * homozygous can never read mid-band. Affymetrix has a floor of 0.0020 to 0.0034 only because it
 * genotypes and computes allelic ratio with separate algorithms that disagree at that rate.
 *
 * So the boundary is not disproved and it is not confirmed: it is untestable on this material, and
 * the tool now says so rather than falling through to a boundary calibrated for clean unamplified
 * arrays. That fall-through was a real bug and it read 23 of 23 known-haploid sperm as diploid.
 * The CALLED column below is what makes that visible; the bands alone never did.
 *
 * BOTH SETS ARE KNOWN HAPLOID BY BIOLOGY, NOT BY THIS TOOL:
 *   GSE148488 pronuclei    one parental complement, from the dissection
 *   GSE19247 sperm cells   one parental complement, from being sperm
 *
 * And both sets get a diploid control from their own series, so a platform difference in the
 * SEPARATION is visible rather than only a difference in level.
 *
 * Run: OM_TRIOS=<dir> OM_GSE=<converted dir> node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/zygosity-crossplatform.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const parentage = await import(`${W}parentage.ts`)

const TRIOS = process.env.OM_TRIOS
const GSE = process.env.OM_GSE
if (!TRIOS || !existsSync(TRIOS)) {
  console.log('zygosity-crossplatform: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}

function rowsOfText(t: string) {
  const ls = t.split('\n')
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
const readAny = (p: string) => (p.endsWith('.gz')
  ? gunzipSync(readFileSync(p)).toString('utf8') : readFileSync(p, 'utf8'))

/**
 * The numbers the zygosity call rests on AND the call itself, for one array.
 *
 * The verdict matters as much as the bands: an array can report a plausible band and still have
 * reached its answer down a branch that does not apply to it, which is exactly what this series
 * found. `classify` is the shipped entry point, so the call here is the call an operator gets.
 */
function bands(path: string) {
  const t = parentage.emptyTally()
  for (const r of rowsOfText(readAny(path))()) parentage.tallyRow('NC' as never, r as never, t as never)
  const tt = t as unknown as {
    bafInBand: number; bafTotal: number; homInBand: number; homTotal: number; homPinned: number
  }
  const hetBand = tt.bafTotal ? tt.bafInBand / tt.bafTotal : NaN
  const homBand = tt.homTotal ? tt.homInBand / tt.homTotal : NaN
  // A fixed 0.17 stands in for the loaded parent's heterozygosity. Only the genotype fallback
  // reads it, and using ONE value for every array means it cannot be a variable here.
  const r = parentage.classify(t as never, 0.17, { role: 'paternal' }) as {
    zygosity: string; homPinnedRate: number; limits: string[]
  }
  return {
    hetBand, homBand, excess: hetBand - homBand,
    pinned: r.homPinnedRate, zygosity: r.zygosity,
    clampNoted: r.limits.some((l) => /exactly 0 or 1/.test(l)),
  }
}

const q = (xs: number[], p: number) => {
  const s = xs.filter(Number.isFinite).sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
type Shape = ReturnType<typeof bands>
/** `expect` is what the biology says, so the CALLED column is scored rather than just printed. */
const line = (label: string, xs: Shape[], expect?: string) => {
  if (!xs.length) { console.log(`  ${label.padEnd(30)} no arrays`); return }
  const e = xs.map((x) => x.excess)
  // REFUSED IS NOT WRONG, and collapsing the two would make a tool that declines to answer look
  // identical to one that answers incorrectly. `unknown` is what classify returns when the call
  // rate is under the floor, which is a refusal by design.
  const right = expect ? xs.filter((x) => x.zygosity === expect).length : NaN
  const refused = xs.filter((x) => x.zygosity === 'unknown').length
  const wrong = expect ? xs.length - right - refused : NaN
  console.log(`  ${label.padEnd(30)} n=${String(xs.length).padStart(3)}  `
    + `band ${q(xs.map((x) => x.hetBand), 0.5).toFixed(4)}  `
    + `homBand ${q(xs.map((x) => x.homBand), 0.5).toFixed(4)}  `
    + `pinned ${(q(xs.map((x) => x.pinned), 0.5) * 100).toFixed(1)}%  `
    + `EXCESS p50 ${q(e, 0.5).toFixed(4)}  `
    + `clamp noted ${xs.filter((x) => x.clampNoted).length}/${xs.length}  `
    + (expect
      ? `${expect.slice(0, 12)}: ${right}/${xs.length} right, ${wrong} WRONG, ${refused} refused`
        + `${wrong > 0 ? '   <-- WRONG CALLS' : ''}`
      : ''))
}

console.log(`the boundary under test: HET_BAND_EXCESS = ${parentage.HET_BAND_EXCESS}`)
console.log('an excess ABOVE it is called two parental contributions, below it one')
console.log('')

// ---------------------------------------------------------------- GSE148488, where it was measured
const lines = readFileSync(join(TRIOS, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], material: f[col('material')], role: f[col('role')],
    pronucleus: f[col('pronucleus')], complete: f[col('complete_trio')] === 'True',
  }
})
const hasT = (g: string) => !!g && existsSync(join(TRIOS, `${g}.probes.gz`))

console.log('=== GSE148488, Affymetrix, the series every threshold was measured on')
{
  const pn = recs.filter((r) => r.pronucleus && hasT(r.gsm))
    .map((r) => bands(join(TRIOS, `${r.gsm}.probes.gz`)))
  line('pronuclei, KNOWN HAPLOID', pn, 'uniparental_homozygous')
  const kids = recs.filter((r) => r.complete && r.role !== 'parent' && !r.pronucleus && hasT(r.gsm))
    .slice(0, 25).map((r) => bands(join(TRIOS, `${r.gsm}.probes.gz`)))
  line('children, known two-parent', kids, 'diploid')
  const par = recs.filter((r) => r.role === 'parent' && hasT(r.gsm))
    .map((r) => bands(join(TRIOS, `${r.gsm}.probes.gz`)))
  line('adult donors, bulk diploid', par, 'diploid')
}

// ---------------------------------------------------------------- GSE19247, the second platform
if (GSE && existsSync(join(GSE, 'truth.json'))) {
  console.log('')
  console.log('=== GSE19247, Illumina, converted by audit/gse19247/toprobes.py')
  interface T { gsm: string; src: string }
  const truth: T[] = JSON.parse(readFileSync(join(GSE, 'truth.json'), 'utf8'))
    .filter((t: T) => existsSync(join(GSE, `${t.gsm}.probes`)))
  const grab = (f: (s: string) => boolean, cap = 40) => truth.filter((t) => f(t.src))
    .slice(0, cap).map((t) => bands(join(GSE, `${t.gsm}.probes`)))
  line('sperm cells, KNOWN HAPLOID', grab((s) => /sperm/i.test(s)), 'uniparental_homozygous')
  // SINGLE CELLS, not bulk: the series titles them 1 to 83. A diploid single cell after whole
  // genome amplification is the hardest case in either series, and it is scored as diploid because
  // that is what it is.
  line('lymphoblasts, SINGLE diploid', grab((s) => /lymphoblast/i.test(s), 30), 'diploid')
  line('blood, BULK diploid', grab((s) => /blood cell/i.test(s)), 'diploid')
} else {
  console.log('')
  console.log('=== GSE19247 NOT INCLUDED: set OM_GSE to the converted directory.')
}

console.log('')
console.log('HOW TO READ IT. If the KNOWN HAPLOID rows sit below the boundary on both platforms,')
console.log('the boundary transfers and anything that read a sperm cell as diploid is the')
console.log('conversion. If the haploid excess is above the boundary on one platform only, the')
console.log('number is calibrated to a series rather than to biology, and that belongs in any')
console.log('claim made from it.')

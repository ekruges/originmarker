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
 * the band is amplification quality rather than ploidy. If the excess on known-haploid material
 * reads the same on both platforms, the boundary transfers and the sperm result is a conversion
 * artefact. If the excess is genuinely larger on one platform for genomes that are equally haploid,
 * the boundary is calibrated to a series rather than to biology.
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

/** The three numbers the zygosity call actually rests on, for one array. */
function bands(path: string) {
  const t = parentage.emptyTally()
  for (const r of rowsOfText(readAny(path))()) parentage.tallyRow('NC' as never, r as never, t as never)
  const tt = t as unknown as {
    bafInBand: number; bafTotal: number; homInBand: number; homTotal: number
  }
  const hetBand = tt.bafTotal ? tt.bafInBand / tt.bafTotal : NaN
  const homBand = tt.homTotal ? tt.homInBand / tt.homTotal : NaN
  return { hetBand, homBand, excess: hetBand - homBand }
}

const q = (xs: number[], p: number) => {
  const s = xs.filter(Number.isFinite).sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
const line = (label: string, xs: { hetBand: number; homBand: number; excess: number }[]) => {
  const e = xs.map((x) => x.excess)
  console.log(`  ${label.padEnd(34)} n=${String(xs.length).padStart(3)}  `
    + `band ${q(xs.map((x) => x.hetBand), 0.5).toFixed(4)}  `
    + `homBand ${q(xs.map((x) => x.homBand), 0.5).toFixed(4)}  `
    + `EXCESS min ${q(e, 0).toFixed(4)} p50 ${q(e, 0.5).toFixed(4)} max ${q(e, 1).toFixed(4)}  `
    + `over ${parentage.HET_BAND_EXCESS}: ${e.filter((x) => x > parentage.HET_BAND_EXCESS).length}/${e.length}`)
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
  line('pronuclei, KNOWN HAPLOID', pn)
  const kids = recs.filter((r) => r.complete && r.role !== 'parent' && !r.pronucleus && hasT(r.gsm))
    .slice(0, 25).map((r) => bands(join(TRIOS, `${r.gsm}.probes.gz`)))
  line('children, known two-parent', kids)
  const par = recs.filter((r) => r.role === 'parent' && hasT(r.gsm))
    .map((r) => bands(join(TRIOS, `${r.gsm}.probes.gz`)))
  line('adult donors, bulk diploid', par)
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
  line('sperm cells, KNOWN HAPLOID', grab((s) => /sperm/i.test(s)))
  line('lymphoblasts, known diploid', grab((s) => /lymphoblast/i.test(s), 30))
  line('blood, known diploid', grab((s) => /blood cell/i.test(s)))
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

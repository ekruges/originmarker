/**
 * Whether a converted B-allele frequency is a BAF at all, checked against its definition.
 *
 * A marker called AB sits in the heterozygous cluster, so by definition its BAF centres on 0.5.
 * That expectation comes from what BAF means, not from either series or any threshold in this
 * tool, so it can judge a conversion without consulting the zygosity call the conversion feeds.
 * Failure is AB calls that do not centre on 0.5: the conversion is then broken, and nothing
 * measured downstream of it says anything about the tool.
 *
 * Saturation is reported beside it. `HET_BAND_EXCESS` subtracts the mid-band mass at homozygous
 * calls, the array's own noise floor, from the mid-band mass overall. A conversion that pins
 * homozygotes to exactly 0 and 1 makes that floor zero, so the correction cannot fire and the call
 * falls through to a flat threshold without any sign of it in the zygosity output.
 *
 * Known material, and where each fact comes from:
 *   GSE148488 adult donors     bulk gDNA of a consenting adult, so diploid, and the series every
 *                              threshold in this tool was measured on
 *   GSE148488 pronuclei        one parental complement, from the micromanipulation
 *   GSE19247  lymphoblasts     a diploid cell line, and the source of that series' trisomy truth
 *   GSE19247  sperm cells      one parental complement, from being sperm
 *
 * A diploid genome carries roughly a quarter to a third of its markers heterozygous, from
 * population genetics rather than any tool output. An AB fraction far below that points at the
 * caller, not the BAF scale.
 *
 * Run: OM_TRIOS=<dir> OM_GSE=<converted dir> node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/baf-shape.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const parentage = await import(`${W}parentage.ts`)

const TRIOS = process.env.OM_TRIOS
const GSE = process.env.OM_GSE

function rowsOf(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
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

const q = (xs: number[], p: number) => {
  const s = xs.filter(Number.isFinite).sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}

/** The shape of one array's BAF, per genotype class. Autosomes only, as the tool's own tally is. */
function shape(path: string) {
  let called = 0; let ab = 0; let hom = 0
  let satAA = 0; let satBB = 0; let aa = 0; let bb = 0
  let homBandN = 0; let bandN = 0; let total = 0
  const abBaf: number[] = []
  for (const r of rowsOf(path)()) {
    if (!parentage.isAutosome(r.chrom)) continue
    if (r.genotype !== 'NC') called += 1
    if (r.baf === null) continue
    total += 1
    const inBand = r.baf >= 0.35 && r.baf <= 0.65
    if (inBand) bandN += 1
    if (r.genotype === 'AB') { ab += 1; if (abBaf.length < 60_000) abBaf.push(r.baf) }
    if (r.genotype === 'AA' || r.genotype === 'BB') {
      hom += 1
      if (inBand) homBandN += 1
      if (r.genotype === 'AA') { aa += 1; if (r.baf === 0) satAA += 1 }
      else { bb += 1; if (r.baf === 1) satBB += 1 }
    }
  }
  return {
    pAB: called ? ab / called : NaN,
    medAB: q(abBaf, 0.5),
    iqrAB: q(abBaf, 0.75) - q(abBaf, 0.25),
    // The share of homozygous calls pinned to the very end of the scale. A real array does not
    // put half its homozygotes on an exact integer.
    sat: hom ? (satAA + satBB) / hom : NaN,
    hetBand: total ? bandN / total : NaN,
    homBand: hom ? homBandN / hom : NaN,
  }
}

type Shape = ReturnType<typeof shape>
const line = (label: string, xs: Shape[]) => {
  if (!xs.length) { console.log(`  ${label.padEnd(30)} no arrays`); return }
  const c = (f: (s: Shape) => number) => q(xs.map(f), 0.5)
  console.log(`  ${label.padEnd(30)} n=${String(xs.length).padStart(3)}  `
    + `AB frac ${c((s) => s.pAB).toFixed(4)}  `
    + `BAF at AB: med ${c((s) => s.medAB).toFixed(4)} iqr ${c((s) => s.iqrAB).toFixed(4)}  `
    + `hom pinned to 0/1 ${(c((s) => s.sat) * 100).toFixed(1)}%  `
    + `hetBand ${c((s) => s.hetBand).toFixed(4)}  homBand ${c((s) => s.homBand).toFixed(4)}`)
}

console.log('THE EXPECTATION, AND IT IS DEFINITIONAL: a marker called AB sits in the heterozygous')
console.log('cluster, so its BAF must centre on 0.5. A diploid genome carries roughly a quarter to a')
console.log('third of markers heterozygous. Neither number comes from this tool.')
console.log('')

if (TRIOS && existsSync(join(TRIOS, 'trio_manifest_full.csv'))) {
  const lines = readFileSync(join(TRIOS, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
  const head = lines[0].split(',')
  const col = (n: string) => head.indexOf(n)
  const recs = lines.slice(1).map((l) => {
    const f = l.split(',')
    return {
      gsm: f[col('gsm')], role: f[col('role')], pronucleus: f[col('pronucleus')],
      complete: f[col('complete_trio')] === 'True',
    }
  })
  const has = (g: string) => !!g && existsSync(join(TRIOS, `${g}.probes.gz`))
  const grab = (f: (r: typeof recs[0]) => boolean, cap = 40) => recs.filter((r) => f(r) && has(r.gsm))
    .slice(0, cap).map((r) => shape(join(TRIOS, `${r.gsm}.probes.gz`)))
  console.log('=== GSE148488, Affymetrix, PUBLISHER-NORMALISED, the series the thresholds came from')
  line('adult donors, DIPLOID', grab((r) => r.role === 'parent'))
  line('children, diploid', grab((r) => r.complete && r.role !== 'parent' && !r.pronucleus, 25))
  line('pronuclei, HAPLOID', grab((r) => !!r.pronucleus))
} else {
  console.log('=== GSE148488 NOT INCLUDED: OM_TRIOS unreadable.')
}

if (GSE && existsSync(join(GSE, 'truth.json'))) {
  console.log('')
  console.log('=== GSE19247, Illumina, CONVERTED BY audit/gse19247/toprobes.py')
  interface T { gsm: string; src: string }
  const truth: T[] = JSON.parse(readFileSync(join(GSE, 'truth.json'), 'utf8'))
    .filter((t: T) => existsSync(join(GSE, `${t.gsm}.probes`)))
  const grab = (f: (s: string) => boolean, cap = 30) => truth.filter((t) => f(t.src))
    .slice(0, cap).map((t) => shape(join(GSE, `${t.gsm}.probes`)))
  line('lymphoblasts, DIPLOID', grab((s) => /lymphoblast/i.test(s)))
  line('blood, diploid', grab((s) => /blood cell/i.test(s)))
  line('sperm cells, HAPLOID', grab((s) => /sperm/i.test(s)))
} else {
  console.log('')
  console.log('=== GSE19247 NOT INCLUDED: set OM_GSE to the converted directory.')
}

console.log('')
console.log('HOW TO READ IT. On the publisher-normalised series the AB median fixes the scale that')
console.log('the whole comparison is against. If the converted arrays put their AB calls anywhere')
console.log('other than 0.5, or pin most homozygotes to an exact 0 or 1, then the conversion does')
console.log('not produce a BAF and every zygosity number measured on it is about the converter.')
console.log('A diploid reading an AB fraction far off a quarter to a third localises it further:')
console.log('low means the CALLER is failing, not just the BAF scale it writes out.')

/**
 * Complete hydatidiform moles on a genotype-only platform.
 *
 * A complete hydatidiform mole is androgenetic: no maternal contribution. A monospermic mole is one
 * sperm genome duplicated, so it is homozygous at every marker and is always 46,XX, because the
 * 46,YY duplicate of a Y-bearing sperm does not survive. A dispermic mole is androgenetic and
 * heterozygous, and may be XY.
 *
 * Scored, each against the sample's own calls or against biology rather than this tool:
 *   zygosity   a mole homozygous genome-wide must be read as one parental contribution, and a
 *              heterozygous one must not be
 *   sex        a homozygous mole carrying a Y is biologically impossible
 *   findings   on a genome with one parental contribution the copy-neutral LOH and uniparental
 *              disomy classes have nothing to report, so a whole-chromosome finding is spurious
 *
 * With OM_MONOSPERMIC=1 the series is taken at its word that every mole carries a genome derived from
 * a single sperm, so every mole is expected homozygous genome-wide and 46,XX, and one that is not
 * contradicts either the series or the call. GSE12713 states this for its whole collection.
 *
 * The parental slot holds a diploid array of the same platform, because the genotype-only zygosity
 * route reads the sample's heterozygosity as a share of the loaded parent's.
 *
 * Run: OM_MOLES=<dir of .probes and truth.json> OM_REF=<a diploid .probes, same platform> \
 *        [OM_MONOSPERMIC=1] \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/mole-genotype.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const sexing = await import(`${W}sexing.ts`)

const DIR = process.env.OM_MOLES
const REF = process.env.OM_REF
if (!DIR || !REF || !existsSync(join(DIR, 'truth.json')) || !existsSync(REF)) {
  console.log('mole-genotype: OM_MOLES (with truth.json) and OM_REF are required. NOT RUN.')
  process.exit(0)
}

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
const isAutosome = (c: string) => /^\d+$/.test(c) && +c >= 1 && +c <= 22

const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(REF)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

type T = { gsm: string; title: string }
const truth: T[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(DIR, `${t.gsm}.probes`)))

interface Row { gsm: string; title: string; het: number; zyg: string; y: boolean | null
  whole: string[]; error?: string }
const out: Row[] = []
for (const t of truth) {
  try {
    let het = 0; let called = 0
    const sx = sexing.emptySex()
    const acc = score.emptyCollected(refParent as never, null)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(join(DIR, `${t.gsm}.probes`))()) {
      if (!first) first = r.probesetId
      if (r.genotype !== 'NC' && isAutosome(r.chrom)) { called += 1; if (r.genotype === 'AB') het += 1 }
      sexing.accumulateSex(r as never, sx as never)
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(t.gsm, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: t.gsm, log: () => {},
    }) as Record<string, unknown>
    const whole = [
      ...((res.findings ?? []) as { cls: string; chrom: string; wholeChromosome?: boolean }[])
        .filter((f) => f.wholeChromosome).map((f) => `${f.cls}:chr${f.chrom}`),
      ...((res.chroms ?? []) as { chrom: string; aneuploidy?: string }[])
        .filter((c) => c.aneuploidy).map((c) => `aneuploidy-${c.aneuploidy}:chr${c.chrom}`),
    ]
    out.push({ gsm: t.gsm, title: t.title, het: called ? het / called : NaN,
      zyg: String(res.zygosity ?? '?'),
      y: (sexing.sexCall(sx as never) as { yBearing: boolean | null }).yBearing, whole })
  } catch (e) {
    out.push({ gsm: t.gsm, title: t.title, het: NaN, zyg: '?', y: null, whole: [],
      error: String((e as Error).message ?? e).slice(0, 80) })
  }
}

const ok = out.filter((r) => !r.error)
const hs = ok.map((r) => r.het).filter(Number.isFinite).sort((a, b) => a - b)
const q = (p: number) => hs[Math.min(hs.length - 1, Math.floor(p * hs.length))]
console.log(`moles: ${out.length}, threw: ${out.length - ok.length}`)
console.log(`genome-wide autosomal heterozygosity: min ${q(0).toFixed(4)}  p10 ${q(0.1).toFixed(4)}  `
  + `p50 ${q(0.5).toFixed(4)}  p90 ${q(0.9).toFixed(4)}  max ${q(1).toFixed(4)}`)

/** Bands read off the sample's own calls: a monospermic mole is genotyping error only. */
const band = (h: number) => (h < 0.02 ? 'homozygous' : h > 0.10 ? 'heterozygous' : 'between')
console.log('')
console.log('=== ZYGOSITY: the tool\'s call against the sample\'s own heterozygosity')
const tab = new Map<string, number>()
for (const r of ok) tab.set(`${band(r.het)} -> ${r.zyg}`, (tab.get(`${band(r.het)} -> ${r.zyg}`) ?? 0) + 1)
for (const [k, n] of [...tab].sort()) console.log(`  ${String(n).padStart(4)}  ${k}`)
const homo = ok.filter((r) => band(r.het) === 'homozygous')
const wrongHomo = homo.filter((r) => r.zyg === 'diploid')
const wrongHet = ok.filter((r) => band(r.het) === 'heterozygous' && r.zyg.startsWith('uniparental'))
console.log(`  homozygous moles read as two parental contributions: ${wrongHomo.length}/${homo.length}`)
console.log(`  heterozygous moles read as one: ${wrongHet.length}/${ok.filter((r) => band(r.het) === 'heterozygous').length}`)
if (process.env.OM_MONOSPERMIC === '1') {
  const notHomo = ok.filter((r) => band(r.het) !== 'homozygous')
  console.log(`  the series states every mole is monospermic; moles that are not homozygous by their own `
    + `calls: ${notHomo.length}/${ok.length}`)
  for (const r of notHomo.slice(0, 8)) console.log(`    ${r.gsm} ${r.title} het ${r.het.toFixed(4)} ${r.zyg}`)
  const called = ok.filter((r) => r.zyg.startsWith('uniparental')).length
  console.log(`  moles the tool reads as one parental contribution: ${called}/${ok.length}`)
}

console.log('')
console.log('=== SEX: a homozygous complete mole is 46,XX, so a Y call on one is impossible')
const yHomo = homo.filter((r) => r.y === true)
console.log(`  homozygous moles: ${homo.length}, Y-bearing ${yHomo.length}, no Y `
  + `${homo.filter((r) => r.y === false).length}, undetermined ${homo.filter((r) => r.y === null).length}`)
for (const r of yHomo.slice(0, 10)) console.log(`    IMPOSSIBLE ${r.gsm} ${r.title} het ${r.het.toFixed(4)}`)

console.log('')
console.log('=== FINDINGS: whole-chromosome emissions on moles homozygous genome-wide')
const spur = homo.filter((r) => r.whole.length > 0)
console.log(`  homozygous moles carrying any whole-chromosome finding: ${spur.length}/${homo.length}`)
const byZyg = new Map<string, number>()
for (const r of spur) byZyg.set(r.zyg, (byZyg.get(r.zyg) ?? 0) + 1)
console.log(`  zygosity of those moles: ${[...byZyg].map(([z, n]) => `${z} ${n}`).join(', ') || '-'}`)
for (const r of spur.slice(0, 8)) console.log(`    ${r.gsm} ${r.zyg}  ${r.whole.slice(0, 6).join(' ')}`)
for (const r of out.filter((x) => x.error).slice(0, 6)) console.log(`  THREW ${r.gsm}: ${r.error}`)

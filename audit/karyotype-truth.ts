/**
 * Blastomere and blastocyst arrays scored against karyotypes their submitters published.
 *
 * GSE20975 gives a 24-chromosome karyotype per cleavage-stage blastomere. GSE18932 gives one per
 * blastocyst sample, trophectoderm and inner cell mass, several samples per embryo. Both were called
 * by the submitters' own microarray method, so agreement here is agreement with a published method
 * on the same data, not with an independent assay. GSE18932's cleavage-stage FISH field describes a
 * different cell at an earlier stage and is not used.
 *
 * These series publish genotype calls only, which fixes what can be asked:
 *   sex        read from Y probe calls, scored against the karyotype's sex chromosomes
 *   monosomy   one copy leaves no heterozygous calls, so the listed chromosome should be the
 *              array's least heterozygous autosome
 *   trisomy    a third copy does not change an AB call, so gains are not scored
 * Positives are only chromosomes a karyotype lists. Negatives are only arrays whose karyotype is
 * normal: the unlisted chromosomes of an abnormal array are unknown, not clean.
 *
 * Run: OM_KARYO=<dir with .probes and truth.json> OM_FIELD=<characteristic holding the karyotype> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/karyotype-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const sexing = await import(`${W}sexing.ts`)

const DIR = process.env.OM_KARYO
const FIELD = process.env.OM_FIELD
if (!DIR || !FIELD || !existsSync(join(DIR, 'truth.json'))) {
  console.log('karyotype-truth: OM_KARYO and OM_FIELD are required. NOT RUN.')
  process.exit(0)
}

type T = { gsm: string; title: string } & Record<string, string>
const truth: T[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(DIR, `${t.gsm}.probes`)))

interface Karyotype { sex: string | null; losses: string[]; gains: string[]; normal: boolean }
/** A karyotype string as the submitters wrote it; null where it is not a karyotype at all. */
function parseKaryotype(raw: string | undefined): Karyotype | null {
  if (!raw) return null
  const s = raw.replace(/\s+/g, '')
  if (!/^\d{2}/.test(s)) return null
  const sex = /^\d{2},?([XYO]+)/.exec(s)?.[1] ?? null
  const losses = [...s.matchAll(/-(\d{1,2}|X|Y)/g)].map((m) => m[1])
  const gains = [...s.matchAll(/\+(\d{1,2}|X|Y)/g)].map((m) => m[1])
  return { sex, losses, gains, normal: !losses.length && !gains.length && (sex === 'XX' || sex === 'XY') }
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

function wilson(k: number, n: number): [number, number] {
  if (!n) return [NaN, NaN]
  const z = 1.96; const p = k / n; const d = 1 + (z * z) / n
  const c = p + (z * z) / (2 * n)
  const h = z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))
  return [Math.max(0, (c - h) / d), Math.min(1, (c + h) / d)]
}
const pct = (x: number) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : '-')
const isAutosome = (c: string) => /^\d+$/.test(c) && +c >= 1 && +c <= 22

const parsed = truth.map((t) => ({ t, k: parseKaryotype(t[FIELD]) }))
const firstNormal = parsed.find((p) => p.k?.normal)
if (!firstNormal) {
  console.log('no normal karyotype in this series to fill the parental slot. NOT RUN.')
  process.exit(0)
}
/** A normal array of the same series in the parental slot; the arms below read the sample alone. */
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR, `${firstNormal.t.gsm}.probes`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Row {
  gsm: string; title: string; k: Karyotype | null; y: boolean | null; zyg: string; stage: string
  het: Map<string, number>; whole: Map<string, string[]>; junk: number; error?: string
}
const out: Row[] = []

for (const { t, k } of parsed) {
  try {
    const sx = sexing.emptySex()
    const counts = new Map<string, [number, number]>()
    const acc = score.emptyCollected(refParent as never, null)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(join(DIR, `${t.gsm}.probes`))()) {
      if (!first) first = r.probesetId
      sexing.accumulateSex(r as never, sx as never)
      if (r.genotype !== 'NC' && isAutosome(r.chrom)) {
        const c = counts.get(r.chrom) ?? [0, 0]
        c[1] += 1
        if (r.genotype === 'AB') c[0] += 1
        counts.set(r.chrom, c)
      }
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(t.gsm, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: t.gsm, log: () => {},
    }) as Record<string, unknown>
    const het = new Map<string, number>()
    for (const [c, [h, n]] of counts) if (n >= 200) het.set(c, h / n)
    const whole = new Map<string, string[]>()
    const add = (c: string, what: string) => whole.set(c, [...(whole.get(c) ?? []), what])
    for (const c of (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]) {
      if (c.aneuploidy) add(c.chrom, `aneuploidy-${c.aneuploidy}`)
    }
    for (const f of (res.findings ?? []) as { cls: string; chrom: string; wholeChromosome?: boolean }[]) {
      if (f.wholeChromosome) add(f.chrom, String(f.cls))
    }
    let junk = 0
    for (const sg of (res.segments ?? []) as { startBp: number; endBp: number; markers: number }[]) {
      if (!(sg.startBp < sg.endBp) || !(sg.markers > 0)) junk += 1
    }
    out.push({
      gsm: t.gsm, title: t.title, k,
      y: (sexing.sexCall(sx as never) as { yBearing: boolean | null }).yBearing,
      zyg: String(res.zygosity ?? '?'), stage: (res.stage as { stage?: string } | undefined)?.stage ?? '?',
      het, whole, junk,
    })
  } catch (e) {
    out.push({ gsm: t.gsm, title: t.title, k, y: null, zyg: '?', stage: '?', het: new Map(), whole: new Map(),
      junk: 0, error: String((e as Error).message ?? e).slice(0, 80) })
  }
}

const ok = out.filter((r) => !r.error)
const withK = ok.filter((r) => r.k)
console.log(`arrays: ${out.length}, threw: ${out.length - ok.length}, karyotype parsed: ${withK.length}, `
  + `normal: ${withK.filter((r) => r.k!.normal).length}, `
  + `listing a monosomy: ${withK.filter((r) => r.k!.losses.some(isAutosome)).length}`)

// ---------------------------------------------------------------- sex
console.log('')
console.log('=== SEX: Y probe calls against the karyotype\'s sex chromosomes')
const sexed = withK.filter((r) => r.k!.sex && r.y !== null)
const agree = sexed.filter((r) => r.y === r.k!.sex!.includes('Y')).length
const [slo, shi] = wilson(agree, sexed.length)
console.log(`  agree ${agree}/${sexed.length} = ${pct(agree / Math.max(sexed.length, 1))}  `
  + `95% interval ${pct(slo)} to ${pct(shi)}`)
console.log(`  no Y call possible (too few Y probes): ${withK.filter((r) => r.y === null).length}`)
for (const r of sexed.filter((x) => x.y !== x.k!.sex!.includes('Y')).slice(0, 12)) {
  console.log(`    ${r.gsm}  karyotype ${r.k!.sex}  tool ${r.y ? 'Y-bearing' : 'no Y'}  ${r.title.slice(0, 40)}`)
}

// ---------------------------------------------------------------- monosomy in the data
console.log('')
console.log('=== MONOSOMY IN THE DATA: heterozygosity of a listed monosomic chromosome, relative to the')
console.log('=== same array\'s median autosome. One copy leaves only genotyping error heterozygous.')
const med = (xs: number[]) => { const s = [...xs].sort((a, b) => a - b); return s.length ? s[s.length >> 1] : NaN }
const q = (xs: number[], p: number) => { const s = [...xs].sort((a, b) => a - b); return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN }
const posRatio: number[] = []; const posRank: number[] = []
for (const r of withK) {
  const m = med([...r.het.values()])
  for (const c of r.k!.losses.filter(isAutosome)) {
    const h = r.het.get(c)
    if (h === undefined || !(m > 0)) continue
    posRatio.push(h / m)
    posRank.push(1 + [...r.het.values()].filter((v) => v < h).length)
  }
}
const negRatio: number[] = []
for (const r of withK.filter((x) => x.k!.normal)) {
  const m = med([...r.het.values()])
  if (m > 0) for (const h of r.het.values()) negRatio.push(h / m)
}
console.log(`  listed monosomic chromosomes: ${posRatio.length}   het ratio p50 ${q(posRatio, 0.5).toFixed(3)}  `
  + `max ${q(posRatio, 1).toFixed(3)}`)
console.log(`  autosomes of normal arrays:   ${negRatio.length}   het ratio p01 ${q(negRatio, 0.01).toFixed(3)}  `
  + `min ${q(negRatio, 0).toFixed(3)}`)
console.log(`  listed chromosome ranked least heterozygous of its array: `
  + `${posRank.filter((x) => x === 1).length}/${posRank.length}`)

// ---------------------------------------------------------------- what the tool emits
console.log('')
console.log('=== WHAT THE TOOL EMITS, whole-chromosome, on each set')
console.log('=== A refused array (stage failed) reports nothing and says so, which is scored apart from a miss.')
const refused = (r: Row) => r.stage === 'failed'
let hit = 0; let posN = 0; let posRefused = 0
for (const r of withK) {
  for (const c of r.k!.losses.filter(isAutosome)) {
    if (refused(r)) { posRefused += 1; continue }
    posN += 1
    if (r.whole.has(c)) hit += 1
  }
}
console.log(`  listed monosomic chromosomes on refused arrays: ${posRefused}`)
console.log(`  listed monosomic chromosome carries a whole-chromosome emission, accepted arrays: ${hit}/${posN}`)
const normalAll = withK.filter((r) => r.k!.normal)
const normal = normalAll.filter((r) => !refused(r))
console.log(`  normal-karyotype arrays refused: ${normalAll.length - normal.length}/${normalAll.length}`)
const falseArrays = normal.filter((r) => r.whole.size > 0)
const [flo, fhi] = wilson(falseArrays.length, normal.length)
console.log(`  accepted normal-karyotype arrays with ANY whole-chromosome emission: ${falseArrays.length}/${normal.length}  `
  + `95% interval ${pct(flo)} to ${pct(fhi)}`)
const byCls = new Map<string, number>()
const byZyg = new Map<string, number>()
for (const r of falseArrays) {
  byZyg.set(r.zyg, (byZyg.get(r.zyg) ?? 0) + 1)
  for (const v of r.whole.values()) for (const cls of v) byCls.set(cls, (byCls.get(cls) ?? 0) + 1)
}
for (const [c, n] of [...byCls].sort((a, b) => b[1] - a[1])) console.log(`      ${n} x ${c}`)
console.log(`      zygosity of those arrays: ${[...byZyg].map(([z, n]) => `${z} ${n}`).join(', ') || '-'}`)

// ---------------------------------------------------------------- same embryo
const embryo = (title: string) => /Embryo\s*(\d+)/i.exec(title)?.[1] ?? null
const groups = new Map<string, Row[]>()
for (const r of ok) { const e = embryo(r.title); if (e) groups.set(e, [...(groups.get(e) ?? []), r]) }
const multi = [...groups.values()].filter((g) => g.length >= 2)
if (multi.length) {
  const split = multi.filter((g) => new Set(g.map((r) => r.y).filter((v) => v !== null)).size > 1)
  console.log('')
  console.log(`=== SAME EMBRYO: ${multi.length} embryos with 2+ samples, sex split on ${split.length}`)
  for (const g of split.slice(0, 8)) {
    console.log(`    ${g.map((r) => `${r.gsm}:${r.y === null ? '?' : r.y ? 'Y' : '-'}`).join(' ')}`)
  }
}

console.log('')
console.log(`structurally invalid segments: ${ok.reduce((a, r) => a + r.junk, 0)}`)
for (const r of out.filter((x) => x.error).slice(0, 6)) console.log(`  THREW ${r.gsm}: ${r.error}`)

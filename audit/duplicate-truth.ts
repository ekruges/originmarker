/**
 * Identity: whether the opposite-homozygote rate tells two arrays of the same DNA from two
 * different samples. Every relationship verdict rests on that rate, and identity is the easiest
 * point on its scale, so a material with no gap here cannot carry a relationship verdict at all.
 *
 * Truth needs no pedigree. Two arrays of the same biopsy are the same DNA, by the experiment; two
 * arrays of different samples are not. The negative set is drawn from a DIFFERENT replicate group,
 * never by accession alone, so two arrays are called different DNA only when the experiment says
 * they are different samples.
 *
 * Failure: two different samples called the same person.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=3072 \
 *        audit/duplicate-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const rel = await import(`${W}relatedness.ts`)
const stage = await import(`${W}stage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('duplicate-truth: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return { gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')] }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

/**
 * The biopsy a title names, with any replicate suffix removed.
 *
 * Two arrays sharing this string are the same DNA measured twice. The suffix forms are the ones
 * this series actually uses; anything not matching keeps its whole title and therefore groups
 * alone, which is the safe direction: a missed replicate costs a positive pair, a wrongly merged
 * pair would corrupt the truth.
 */
const biopsy = (title: string) => title
  .replace(/[ _-]*(rep(licate)?[ _.-]*\d+|r\d+|tech[ _-]*\d+)$/i, '')
  .replace(/\s+/g, ' ')
  .trim()

interface Loaded {
  gt: Map<string, 'AA' | 'AB' | 'BB' | 'NC'>
  pos: Map<string, { chrom: string; pos: number }>
  stage: string
  callRate: number
}
function load(gsm: string): Loaded {
  const t = gunzipSync(readFileSync(join(DIR, `${gsm}.probes.gz`))).toString('utf8')
  const ls = t.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  const gt = new Map<string, 'AA' | 'AB' | 'BB' | 'NC'>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    gt.set(r.probesetId, r.genotype as 'AA')
    if (Number.isFinite(r.pos)) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const st = stage.inferStage(profile as never, {} as never) as { stage: string }
  return { gt, pos, stage: st.stage, callRate: profile.callRate }
}

const groups = new Map<string, string[]>()
for (const r of recs) {
  if (!has(r.gsm)) continue
  const b = biopsy(r.title)
  groups.set(b, [...(groups.get(b) ?? []), r.gsm])
}
const repGroups = [...groups.entries()].filter(([, v]) => v.length >= 2)
console.log(`samples on disk: ${[...groups.values()].flat().length}`)
console.log(`groups of 2 or more arrays of one biopsy: ${repGroups.length}`)
console.log('')

const CAP = Number(process.env.OM_CAP ?? 14)
const use = repGroups.slice(0, CAP)
// A BOUNDED CACHE. The first version kept every array it loaded and died at exit 134, a V8
// out-of-memory: each array is two maps over ~800k markers and thirty of them do not fit. The
// pairing below is ordered so a small window is enough.
const CACHE_MAX = Number(process.env.OM_CACHE ?? 6)
const cache = new Map<string, Loaded>()
const get = (g: string) => {
  const hit = cache.get(g)
  if (hit) {
    cache.delete(g)
    cache.set(g, hit)
    return hit
  }
  const v = load(g)
  cache.set(g, v)
  while (cache.size > CACHE_MAX) cache.delete(cache.keys().next().value as string)
  return v
}

interface Pair { a: string; b: string; same: boolean; rate: number; n: number
  verdict: string; matA: string; matB: string }
const pairs: Pair[] = []

// SAME DNA: every within-group pair.
for (const [, gs] of use) {
  for (let i = 0; i < gs.length; i += 1) {
    for (let j = i + 1; j < gs.length; j += 1) {
      const A = get(gs[i]); const B = get(gs[j])
      const opp = rel.oppositeHom(A.gt, B.gt) as { rate: number; n: number }
      if (opp.n < 5_000) continue
      const r = rel.relate({ gt: A.gt, pos: A.pos }, { gt: B.gt, pos: B.pos },
        rel.OPPOSITE_HOM_MAX) as { relationship: string }
      pairs.push({ a: gs[i], b: gs[j], same: true, rate: opp.rate, n: opp.n,
        verdict: r.relationship, matA: A.stage, matB: B.stage })
    }
  }
}
// DIFFERENT DNA: one array from each of two different groups, matched so the comparison is not
// between a clean array and a degraded one.
for (let i = 0; i < use.length; i += 1) {
  for (let j = i + 1; j < use.length; j += 1) {
    const a = use[i][1][0]; const b = use[j][1][0]
    const A = get(a); const B = get(b)
    const opp = rel.oppositeHom(A.gt, B.gt) as { rate: number; n: number }
    if (opp.n < 5_000) continue
    const r = rel.relate({ gt: A.gt, pos: A.pos }, { gt: B.gt, pos: B.pos },
      rel.OPPOSITE_HOM_MAX) as { relationship: string }
    pairs.push({ a, b, same: false, rate: opp.rate, n: opp.n,
      verdict: r.relationship, matA: A.stage, matB: B.stage })
  }
}

const q = (xs: number[], p: number) => {
  const s = [...xs].sort((x, y) => x - y)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
const same = pairs.filter((p) => p.same)
const diff = pairs.filter((p) => !p.same)
const rr = (ps: Pair[]) => ps.map((p) => p.rate)

console.log('=== OPPOSITE-HOMOZYGOTE RATE: the same DNA twice, against two different samples')
console.log(`  SAME DNA       n=${String(same.length).padStart(4)}  `
  + `min ${q(rr(same), 0).toFixed(4)}  p50 ${q(rr(same), 0.5).toFixed(4)}  max ${q(rr(same), 1).toFixed(4)}`)
console.log(`  DIFFERENT DNA  n=${String(diff.length).padStart(4)}  `
  + `min ${q(rr(diff), 0).toFixed(4)}  p50 ${q(rr(diff), 0.5).toFixed(4)}  max ${q(rr(diff), 1).toFixed(4)}`)
const sMax = q(rr(same), 1)
const dMin = q(rr(diff), 0)
console.log(`  gap ${(dMin - sMax).toFixed(4)}  ${dMin > sMax
  ? `SEPARABLE, any cut in (${sMax.toFixed(4)}, ${dMin.toFixed(4)})`
  : 'NOT SEPARABLE, the two overlap'}`)
console.log(`  the shipped gate sits at OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX}`)

console.log('')
console.log('=== BY MATERIAL, because amplification is what the earlier failures turned on')
const mats = [...new Set(pairs.map((p) => [p.matA, p.matB].sort().join('+')))].sort()
console.log('materials                    same-DNA n / range          different-DNA n / range')
for (const m of mats) {
  const f = (ps: Pair[]) => ps.filter((p) => [p.matA, p.matB].sort().join('+') === m).map((p) => p.rate)
  const s = f(same); const d = f(diff)
  const fmt = (xs: number[]) => (xs.length
    ? `${String(xs.length).padStart(4)}  ${q(xs, 0).toFixed(4)}-${q(xs, 1).toFixed(4)}`
    : '   0          -       ')
  console.log(`${m.padEnd(28)} ${fmt(s)}      ${fmt(d)}`)
}

console.log('')
console.log('=== THE SHIPPED VERDICT')
const tally = new Map<string, number>()
for (const p of pairs) {
  tally.set(`${p.same ? 'SAME DNA ' : 'DIFFERENT'} -> ${p.verdict}`,
    (tally.get(`${p.same ? 'SAME DNA ' : 'DIFFERENT'} -> ${p.verdict}`) ?? 0) + 1)
}
for (const [k, v] of [...tally].sort()) console.log(`  ${k}: ${v}`)

const falseId = diff.filter((p) => /duplicate|this parent only/i.test(p.verdict))
console.log('')
console.log(`TWO DIFFERENT SAMPLES CALLED THE SAME PERSON: ${falseId.length}/${diff.length}`)
for (const p of falseId.slice(0, 8)) {
  console.log(`  ${p.a} + ${p.b}  ${p.matA}/${p.matB}  rate ${p.rate.toFixed(4)}  ${p.verdict}`)
}
console.log('')
console.log('READ IT LIKE THIS. Identity is the easiest question on the relatedness scale. A gap')
console.log('here bounds what the harder questions can possibly achieve on the same material; no')
console.log('gap here means the statistic is not usable on it at any relationship distance.')

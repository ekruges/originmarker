/**
 * Relatedness on GSE19247, where the submitters label every sample with its family.
 *
 * Truth: the submitters' family labels, on material from another laboratory and chemistry that this
 * tool was never calibrated against. What is known decides what is scored:
 *
 *   ACROSS families      unrelated. Any relative verdict is a false positive.
 *   WITHIN one family    related, but the degree is not published, and several arrays of one
 *                        family may be single cells of the SAME person rather than relatives.
 *                        So a within-family pair is scored only as "should not read unrelated",
 *                        never as a specific relationship.
 *
 * Reading a within-family pair as a specific degree would be inventing truth the metadata does not
 * carry. The across-family arm is where the falsifiable claim lives.
 *
 * Run: OM_GSE=<converted dir with truth.json> node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/gse19247-relatedness.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const rel = await import(`${W}relatedness.ts`)

const GSE = process.env.OM_GSE
if (!GSE || !existsSync(join(GSE, 'truth.json'))) {
  console.log('gse19247-relatedness: OM_GSE must hold .probes and truth.json. NOT RUN.')
  process.exit(0)
}

interface T { gsm: string; title: string; src: string }
const truth: T[] = JSON.parse(readFileSync(join(GSE, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(GSE, `${t.gsm}.probes`)))

/** The family a sample belongs to, from the submitters' own source_name. */
function family(src: string): string | null {
  const m = /family\s+(\d+)/i.exec(src)
  return m ? m[1] : null
}
/** What the material is, because a haploid product behaves differently and is reported apart. */
function material(src: string): string {
  if (/sperm/i.test(src)) return 'sperm'
  if (/blood/i.test(src)) return 'blood'
  if (/lymphoblast/i.test(src)) return 'lymphoblast'
  if (/cleavage/i.test(src)) return 'embryo'
  return 'other'
}

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Loaded { gt: Map<string, AB>; pos: Map<string, { chrom: string; pos: number }> }
function load(gsm: string): Loaded {
  const ls = readFileSync(join(GSE, `${gsm}.probes`), 'utf8').split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  const gt = new Map<string, AB>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    gt.set(r.probesetId, r.genotype as AB)
    if (Number.isFinite(r.pos)) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  return { gt, pos }
}

const pool = truth.map((t) => ({ ...t, fam: family(t.src), mat: material(t.src) }))
  .filter((t) => t.fam)
const byFam = new Map<string, typeof pool>()
for (const p of pool) byFam.set(p.fam!, [...(byFam.get(p.fam!) ?? []), p])

console.log(`samples carrying a family label: ${pool.length} across ${byFam.size} families`)
for (const [f, ps] of byFam) {
  const mats = [...new Set(ps.map((p) => p.mat))].join(', ')
  console.log(`  family ${f}: ${ps.length} arrays (${mats})`)
}
console.log('')
console.log(`the gate under test: OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX}`)
console.log('')

const PER = Number(process.env.OM_PER ?? 4)
/** A few arrays per family, so the pairing does not run to thousands of loads. */
const picked = [...byFam.entries()].flatMap(([f, ps]) => {
  const seen = new Map<string, number>()
  return ps.filter((p) => {
    const n = seen.get(p.mat) ?? 0
    if (n >= PER) return false
    seen.set(p.mat, n + 1)
    return true
  }).map((p) => ({ ...p, fam: f }))
})
console.log(`arrays loaded for pairing: ${picked.length} (up to ${PER} per material per family)`)

const cache = new Map<string, Loaded>()
const get = (g: string) => {
  const hit = cache.get(g)
  if (hit) return hit
  const v = load(g)
  cache.set(g, v)
  return v
}

interface Pair {
  a: string; b: string; sameFamily: boolean
  mats: string; rate: number; n: number; verdict: string
}
const pairs: Pair[] = []
for (let i = 0; i < picked.length; i += 1) {
  for (let j = i + 1; j < picked.length; j += 1) {
    const A = picked[i]; const B = picked[j]
    const la = get(A.gsm); const lb = get(B.gsm)
    const opp = rel.oppositeHom(la.gt, lb.gt) as { rate: number; n: number }
    if (opp.n < 5_000) continue
    const r = rel.relate({ gt: la.gt, pos: la.pos }, { gt: lb.gt, pos: lb.pos },
      rel.OPPOSITE_HOM_MAX) as { relationship: string }
    pairs.push({
      a: A.gsm, b: B.gsm, sameFamily: A.fam === B.fam,
      mats: [A.mat, B.mat].sort().join('+'),
      rate: opp.rate, n: opp.n, verdict: r.relationship,
    })
  }
}

const q = (xs: number[], p: number) => {
  const s = [...xs].sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
const same = pairs.filter((p) => p.sameFamily)
const diff = pairs.filter((p) => !p.sameFamily)

console.log('')
console.log('=== OPPOSITE-HOMOZYGOTE RATE, by whether the metadata puts them in one family')
console.log(`  SAME family     n=${String(same.length).padStart(4)}  `
  + `min ${q(same.map((p) => p.rate), 0).toFixed(4)}  p50 ${q(same.map((p) => p.rate), 0.5).toFixed(4)}  `
  + `max ${q(same.map((p) => p.rate), 1).toFixed(4)}`)
console.log(`  DIFFERENT family n=${String(diff.length).padStart(4)}  `
  + `min ${q(diff.map((p) => p.rate), 0).toFixed(4)}  p50 ${q(diff.map((p) => p.rate), 0.5).toFixed(4)}  `
  + `max ${q(diff.map((p) => p.rate), 1).toFixed(4)}`)
const sMax = q(same.map((p) => p.rate), 1)
const dMin = q(diff.map((p) => p.rate), 0)
console.log(`  gap: ${(dMin - sMax).toFixed(4)}  ${dMin > sMax
  ? `SEPARABLE, any cut in (${sMax.toFixed(4)}, ${dMin.toFixed(4)})`
  : 'NOT SEPARABLE, they overlap'}`)

console.log('')
console.log('=== THE SAME SPLIT, BY MATERIAL PAIRING')
const mats = [...new Set(pairs.map((p) => p.mats))].sort()
console.log('materials              same-family n / range            different-family n / range')
for (const m of mats) {
  const s = pairs.filter((p) => p.mats === m && p.sameFamily).map((p) => p.rate)
  const d = pairs.filter((p) => p.mats === m && !p.sameFamily).map((p) => p.rate)
  const fmt = (xs: number[]) => (xs.length
    ? `${String(xs.length).padStart(4)}  ${q(xs, 0).toFixed(4)}-${q(xs, 1).toFixed(4)}` : '   0        -      ')
  console.log(`${m.padEnd(22)} ${fmt(s)}        ${fmt(d)}`)
}

console.log('')
console.log('=== THE SHIPPED VERDICT, against what the metadata says')
const tally = new Map<string, number>()
for (const p of pairs) {
  const k = `${p.sameFamily ? 'SAME family ' : 'DIFFERENT   '} -> ${p.verdict}`
  tally.set(k, (tally.get(k) ?? 0) + 1)
}
for (const [k, v] of [...tally].sort()) console.log(`  ${k}: ${v}`)

// The falsifiable claim: a pair from two different families must never read first-degree.
const wrong = diff.filter((p) => /parent and child|first-degree|sibling|this parent only/i.test(p.verdict))
console.log('')
console.log(`FALSE RELATIVES: different-family pairs called any kind of relative: ${wrong.length}/${diff.length}`)
for (const p of wrong.slice(0, 10)) {
  console.log(`  ${p.a} + ${p.b}  ${p.mats}  rate ${p.rate.toFixed(4)}  ${p.verdict}`)
}
console.log('')
console.log('READ IT LIKE THIS. Different-family pairs are unrelated by the submitters\' own')
console.log('labelling, so any relative verdict on one is a false positive with nothing to argue')
console.log('about. Same-family pairs are related but the degree is not published and some may be')
console.log('one person twice, so they are not scored for a specific relationship.')

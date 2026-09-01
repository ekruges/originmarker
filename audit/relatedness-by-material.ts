/**
 * Can the opposite-homozygote rate tell a true parent from a stranger ON AMPLIFIED MATERIAL?
 *
 * THE REASON THIS EXISTS. `om link` separates a first-degree relative from an unrelated adult at an
 * opposite-homozygote rate of 0.020, and that number is sound on the material it was measured on:
 * bulk genomic DNA, where two homozygous calls disagreeing means the two genomes really do disagree.
 * Wiring the same gate into the scorer withheld every named parent on eight of eight CORRECT runs.
 * The first real single cell it saw, a blastomere against its genetically confirmed father, read
 * 0.0276 over 524,821 markers and came back `unrelated`.
 *
 * THE MECHANISM IS NOT SUBTLE. Amplification drops alleles. A heterozygous marker that loses one
 * allele is CALLED homozygous, and if it kept the allele the parent does not carry, the pair reads
 * as opposite homozygotes. So the statistic is inflated by exactly the thing single-cell material
 * has most of, and it is inflated in the direction that makes a true parent look like a stranger.
 *
 * WHAT THIS MEASURES. The rate for every usable trio of GSE148488, against three arrays: the child's
 * genetically confirmed father, its confirmed mother, and an unrelated adult. Split by material,
 * because dropout is a property of how much template there was. If the true-parent and unrelated
 * populations separate within a material, a per-material threshold is possible and the numbers here
 * are what it should be set from. If they overlap, this statistic cannot gate anything on that
 * material and the scorer must not pretend otherwise.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 \
 *        audit/relatedness-by-material.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const rel = await import(`${W}relatedness.ts`)
const stageMod = await import(`${W}stage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('relatedness-by-material: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

/** Excluded as a parental genotype: het 0.3455 at an 82% call rate. */
const CONTAMINATED = 'GSM4472408'

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  const first = (s: string) => (s ?? '').split(';').filter((x) => x && x !== CONTAMINATED)[0] ?? ''
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')],
    father: first(f[col('father_gsms')]), mother: first(f[col('mother_gsms')]),
    complete: f[col('complete_trio')] === 'True',
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

/**
 * ARRAYS THE TOOL ITSELF REFUSES ARE EXCLUDED, AND THAT IS NOT A CONVENIENCE.
 *
 * The first version of this table pooled every array, and trophectoderm came out unseparable with a
 * median opposite-homozygote rate of 0.53 against the true father. Chasing that number found 23
 * arrays reading 55 to 62 percent heterozygosity. No diploid genome is 61 percent heterozygous; the
 * panel's expectation is about 17. They are mixed or contaminated samples, the same condition that
 * excluded egg donor D at 34.5 percent, and `stage.inferStage` already calls every one of them
 * `failed`. In a real run they never reach an origin channel at all: the failed-array guard empties
 * the result before any of this.
 *
 * So pooling them measured a population the tool does not serve, and reported the answer as though
 * it applied to the population it does. Both are printed below.
 */
const stageOf = new Map<string, string>()
function accepted(gsm: string): boolean {
  const hit = stageOf.get(gsm)
  if (hit !== undefined) return hit !== 'failed'
  const text = gunzipSync(readFileSync(join(DIR!, `${gsm}.probes.gz`))).toString('utf8')
  const ls = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
  }
  const st = stageMod.inferStage(ingest.finishProfile(gsm, byChrom as never, bafSums as never, first),
    {}) as { stage: string }
  stageOf.set(gsm, st.stage)
  return st.stage !== 'failed'
}

type AB = 'AA' | 'AB' | 'BB' | 'NC'
const cache = new Map<string, Map<string, AB>>()
function genotypes(gsm: string, keep = false): Map<string, AB> {
  const hit = cache.get(gsm)
  if (hit) return hit
  const text = gunzipSync(readFileSync(join(DIR, `${gsm}.probes.gz`))).toString('utf8')
  const ls = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  const gt = new Map<string, AB>()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (r) gt.set(r.probesetId, r.genotype as AB)
  }
  // Parents are reused across every trio and worth holding; a child is read once and dropped.
  if (keep) cache.set(gsm, gt)
  return gt
}

const trios = recs.filter((r) => r.complete && r.role !== 'parent' && has(r.gsm)
  && has(r.father) && has(r.mother) && r.mother !== CONTAMINATED)
/**
 * WHO EACH DONOR ARRAY IS, collapsed across replicates.
 *
 * THE FIRST VERSION OF THIS FILE PICKED THE "UNRELATED" ARRAY BY GSM, taking the first donor that
 * was not this trio's father or mother. The sperm donor has four arrays, so for a trio fathered by
 * rep 1 that selected rep 2: the same person. The unrelated arm was measuring the true father
 * against himself and reported, inevitably, that a true parent and a stranger are
 * indistinguishable. Identity comes from the GEO title, which names the person.
 */
const person = (title: string) => (/sperm/i.test(title) ? 'sperm'
  : (/egg donor ([A-Z])/i.exec(title)?.[1] ?? title))
const donors = recs.filter((r) => r.role === 'parent' && has(r.gsm) && r.gsm !== CONTAMINATED)
  .map((r) => ({ id: r.gsm, who: person(r.title) }))
const whoOf = new Map(donors.map((d) => [d.id, d.who]))

console.log(`trios ${trios.length}, donor arrays ${donors.length}, `
  + `distinct people ${new Set(donors.map((d) => d.who)).size}`)
console.log(`the gate under test: OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX}`)
console.log('')

interface Row { material: string; kind: string; rate: number; n: number; usable: boolean }
const rows: Row[] = []

for (const t of trios) {
  const usable = accepted(t.gsm)
  const child = genotypes(t.gsm)
  // An unrelated adult: a donor who is neither of this child's parents. Same panel, same platform,
  // same laboratory, so the only difference from the true parents is the relationship.
  const kin = new Set([whoOf.get(t.father), whoOf.get(t.mother)])
  const stranger = donors.find((d) => !kin.has(d.who))
  if (!stranger) continue
  const strangerId = stranger.id
  // A SECOND ARRAY OF THE TRUE FATHER, where he has one. This is the arm whose absence let the
  // first version of this file go wrong: if "unrelated" and "father's replicate" read the same,
  // the unrelated arm is not measuring unrelatedness.
  const replicate = donors.find((d) => d.who === whoOf.get(t.father) && d.id !== t.father)?.id
  for (const [kind, pid] of ([
    ['true father', t.father], ['true mother', t.mother],
    ...(replicate ? [["father's replicate", replicate]] : []),
    ['unrelated', strangerId],
  ] as [string, string][])) {
    const o = rel.oppositeHom(genotypes(pid, true), child)
    rows.push({ material: t.material, kind, rate: o.rate, n: o.n, usable })
  }
}

const q = (xs: number[], p: number) => {
  if (!xs.length) return NaN
  const s = [...xs].sort((a, b) => a - b)
  return s[Math.min(s.length - 1, Math.floor(p * s.length))]
}

const failed = [...new Set(rows.filter((r) => !r.usable).map((r) => r.material))]
console.log(`arrays the tool refuses outright (stage "failed"), excluded below: `
  + `${rows.filter((r) => !r.usable && r.kind === 'true father').length} of `
  + `${rows.filter((r) => r.kind === 'true father').length}`
  + (failed.length ? `, all on ${failed.join(', ')}` : ''))
console.log('')
console.log('material         relationship   n    min      p50      max      over the 0.020 gate')
const materials = [...new Set(rows.filter((r) => r.usable).map((r) => r.material))].sort()
for (const m of materials) {
  for (const kind of ['true father', 'true mother', "father's replicate", 'unrelated']) {
    const xs = rows.filter((r) => r.usable && r.material === m && r.kind === kind)
      .map((r) => r.rate).filter((x) => Number.isFinite(x))
    if (!xs.length) continue
    const over = xs.filter((x) => x > rel.OPPOSITE_HOM_MAX).length
    console.log(`${m.padEnd(16)} ${kind.padEnd(13)} ${String(xs.length).padStart(3)}  `
      + `${q(xs, 0).toFixed(4)}   ${q(xs, 0.5).toFixed(4)}   ${q(xs, 1).toFixed(4)}   `
      + `${over}/${xs.length}`)
  }
}

// THE QUESTION THIS FILE EXISTS TO ANSWER. Within each material, is there any threshold at all that
// keeps every true parent and rejects every stranger? If the highest true-parent rate sits above the
// lowest unrelated rate, no single number can do it and the statistic cannot gate on that material.
console.log('')
console.log('CAN A THRESHOLD EXIST, per material?')
for (const m of materials) {
  // The replicate arm is a control on the harness, not a parent under test, so it is excluded
  // from the separability question.
  const t = rows.filter((r) => r.usable && r.material === m
    && (r.kind === 'true father' || r.kind === 'true mother'))
    .map((r) => r.rate).filter(Number.isFinite)
  const u = rows.filter((r) => r.usable && r.material === m && r.kind === 'unrelated')
    .map((r) => r.rate).filter(Number.isFinite)
  if (!t.length || !u.length) continue
  const hiTrue = Math.max(...t)
  const loUnrel = Math.min(...u)
  const gap = loUnrel - hiTrue
  console.log(`  ${m.padEnd(16)} worst true parent ${hiTrue.toFixed(4)}, best stranger `
    + `${loUnrel.toFixed(4)}  ->  ${gap > 0
      ? `SEPARABLE, any threshold in (${hiTrue.toFixed(4)}, ${loUnrel.toFixed(4)})`
      : `NOT SEPARABLE, they overlap by ${(-gap).toFixed(4)}`}`)
}

console.log('')
console.log('relatedness-by-material: a threshold measured on bulk DNA is not a threshold on a '
  + 'single cell, and this is the table that says so.')

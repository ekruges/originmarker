/**
 * DR EGLI'S QUESTION, PUT TO THE SHIPPED CODE ON REAL ARRAYS.
 *
 * "Whether breakpoints and losses affect paternal or maternal genomes differentially, or equally."
 *
 * THE CORPUS ANSWERS IT IN ONE PLACE AND ONE PLACE ONLY: the pronuclei. Each carries exactly one
 * parental complement, so every change found in one belongs to that parent by construction, with no
 * per-event attribution needed and no detection floor to clear. Seven are maternal and eight are
 * paternal. Every other material in this set is biparental, where attributing an event to a parent
 * is a per-event call with its own error rate, and stacking that error on top of the comparison is
 * how a difference in the ASSAY becomes a difference in the biology.
 *
 * THE COMPARISON IS `parentalBalance`, WHICH SHIPS. Nothing is re-implemented here. It carries a
 * cluster-robust permutation test, a minimum group size, and MAX_POWER_SKEW: a refusal to compare
 * two groups whose detection power differs by more than 1.4-fold, on the grounds that the better
 * measured group would show more events whatever the biology. That guard is the reason this harness
 * can be trusted to say "cannot tell" rather than inventing an answer.
 *
 * WHAT THIS PRINTS AND WHY EACH PART IS THERE.
 *   the denominator      how many genomes and how many events are actually in each group, because a
 *                        comparison of five events against six is not a finding whatever its p value
 *   the power check      the skew between the groups, and whether the shipped guard accepted them
 *   the answer           differential or not, with the interval, from the shipped test
 *   without the X        a male genome carries a maternal X and no paternal one, so leaving the sex
 *                        chromosomes in manufactures a maternal excess out of ordinary sex
 *   by egg donor         the corpus has ONE sperm donor and several egg donors, so donor identity is
 *                        confounded with parent; this splits it where the counts allow
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 \
 *        audit/parental-balance-real.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const balance = await import(`${W}parentalBalance.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('parental-balance-real: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

const CONTAMINATED = 'GSM4472408'
const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')], pronucleus: f[col('pronucleus')],
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

const textOf = (g: string) => gunzipSync(readFileSync(join(DIR, `${g}.probes.gz`))).toString('utf8')
function rowsFromText(t: string) {
  const ls = t.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}
const parentCache = new Map<string, ReturnType<typeof buildParent>>()
function buildParent(gsm: string) {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(gsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
}
const parentOf = (g: string) => {
  const hit = parentCache.get(g)
  if (hit) return hit
  const p = buildParent(g)
  parentCache.set(g, p)
  return p
}
const gtCache = new Map<string, Map<string, string>>()
function genotypes(g: string) {
  const hit = gtCache.get(g)
  if (hit) return hit
  const m = new Map<string, string>()
  for (const r of rowsFromText(textOf(g))()) m.set(r.probesetId, r.genotype)
  gtCache.set(g, m)
  return m
}

const donors = recs.filter((r) => r.role === 'parent' && has(r.gsm) && r.gsm !== CONTAMINATED)
const who = (t: string) => (/sperm/i.test(t) ? 'sperm' : (/egg donor ([A-Z])/i.exec(t)?.[1] ?? t))

function resolveDonor(childGsm: string): { id: string; who: string } | null {
  const c = genotypes(childGsm)
  const scored = donors.map((d) => {
    const p = genotypes(d.gsm)
    let opp = 0
    let n = 0
    for (const [probe, cg] of c) {
      if (cg !== 'AA' && cg !== 'BB') continue
      const pg = p.get(probe)
      if (pg !== 'AA' && pg !== 'BB') continue
      n += 1
      if (pg !== cg) opp += 1
    }
    return { id: d.gsm, who: who(d.title), rate: n ? opp / n : Infinity }
  }).sort((a, b) => a.rate - b.rate)
  const best = scored[0]
  const next = scored.find((x) => x.who !== best.who)
  if (!next || !(next.rate > best.rate * 2)) return null
  return { id: best.id, who: best.who }
}

const isAutosome = (c: string) => /^\d+$/.test(c.replace(':PAR', ''))

interface Built {
  sample: balance.BalanceSample
  donor: string
  sexEvents: number
}
const built: Built[] = []

const pronuclei = recs.filter((r) => r.pronucleus && has(r.gsm))
console.log(`pronuclei available: ${pronuclei.length}`)
console.log('')

for (const p of pronuclei) {
  const res = resolveDonor(p.gsm)
  if (!res) { console.log(`  ${p.gsm}: donor unresolved by linkage, excluded`); continue }
  // THE PARENT THIS GENOME CARRIES. A maternal pronucleus carries the egg donor's complement.
  const carries: 'maternal' | 'paternal' = p.pronucleus === 'maternal' ? 'maternal' : 'paternal'
  // Score it against the array it actually carries, under that array's true role, which is the run
  // an operator would perform once they knew what the sample was.
  const pat = parentOf(res.id)
  const role: 'paternal' | 'maternal' = res.who === 'sperm' ? 'paternal' : 'maternal'
  const acc = score.emptyCollected(pat as never, null)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(p.gsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, null, acc as never)
  }
  const profile = ingest.finishProfile(p.gsm, byChrom as never, bafSums as never, first)
  const out = await score.scoreSample({
    acc, profile, pat, mat: null, soloRole: role, sibs: [], sampleName: p.gsm, log: () => {},
  }) as Record<string, unknown>

  const segs = (out.segments ?? []) as { chrom: string; kind: string; startBp?: number
    endBp?: number; refined?: { startBp: number; endBp: number } }[]
  const chroms = (out.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
  const events = [
    ...segs.map((s) => ({
      cls: s.kind, chrom: s.chrom,
      startBp: s.refined?.startBp ?? s.startBp ?? 0,
      endBp: s.refined?.endBp ?? s.endBp ?? 0,
    })),
    ...chroms.filter((c) => c.aneuploidy).map((c) => ({
      cls: `whole-${c.aneuploidy}`, chrom: c.chrom, startBp: 0, endBp: 250e6,
    })),
  ]
  built.push({
    donor: res.who,
    sexEvents: events.filter((e) => !isAutosome(e.chrom)).length,
    sample: {
      name: p.gsm,
      parent: carries,
      originClass: String(out.originClass ?? 'unclear'),
      informative: Number(out.informative ?? NaN),
      explainable: Number(out.explainable ?? NaN),
      material: String((out.stage as { stage?: string } | undefined)?.stage ?? p.material),
      events,
    },
  })
  console.log(`  ${p.gsm} carries ${carries.padEnd(9)} donor ${res.who.padEnd(6)} `
    + `class ${String(out.originClass).padEnd(16)} events ${events.length} `
    + `(${events.filter((e) => !isAutosome(e.chrom)).length} on sex chromosomes)`)
}

// ------------------------------------------------------------------ THE DENOMINATOR
console.log('')
console.log('=== THE DENOMINATOR, before any test')
const byParent = (p: string) => built.filter((b) => b.sample.parent === p)
for (const p of ['paternal', 'maternal'] as const) {
  const g = byParent(p)
  const ev = g.reduce((a, b) => a + b.sample.events.length, 0)
  console.log(`  genomes carrying ${p.padEnd(9)} ${g.length}, events in them ${ev}, `
    + `carrying at least one ${g.filter((b) => b.sample.events.length > 0).length}`)
}
console.log(`  the shipped minimum per group is ${balance.MIN_PER_GROUP}, and the size at which it`)
console.log(`  will report rather than only refuse is ${balance.REPORTING_PER_GROUP}`)

const report = (label: string, samples: balance.BalanceSample[]) => {
  console.log('')
  console.log(`=== ${label}`)
  if (samples.length < 2) { console.log('  too few samples to compare'); return }
  const r = balance.parentalBalance(samples, { permutations: 20_000, seed: 7 })
  const j = r as unknown as Record<string, unknown>
  for (const k of Object.keys(j)) {
    const v = j[k]
    if (v === null || v === undefined) continue
    const s = typeof v === 'object' ? JSON.stringify(v) : String(v)
    console.log(`  ${k}: ${s.length > 400 ? `${s.slice(0, 400)}...` : s}`)
  }
}

report('THE ANSWER, all events, every pronucleus', built.map((b) => b.sample))

// ------------------------------------------------------------------ WITHOUT THE SEX CHROMOSOMES
//
// A male genome carries a maternal X and no paternal one, so leaving them in manufactures a
// maternal excess out of ordinary sex determination rather than out of damage.
report('AUTOSOMES ONLY, which is the comparison that means anything',
  built.map((b) => ({
    ...b.sample,
    events: b.sample.events.filter((e) => isAutosome(e.chrom)),
  })))

// ------------------------------------------------------------------ DONOR CONFOUND
//
// One sperm donor, several egg donors. Any quality difference between that one array and the egg
// donor arrays reads as a parental effect. This splits it as far as the counts allow.
console.log('')
console.log('=== THE DONOR CONFOUND')
const donorCounts = new Map<string, number>()
for (const b of built) donorCounts.set(b.donor, (donorCounts.get(b.donor) ?? 0) + 1)
console.log(`  genomes by donor: ${JSON.stringify(Object.fromEntries(donorCounts))}`)
console.log('  Every paternal genome comes from ONE person. If that person\'s arrays differ in')
console.log('  quality from the egg donors\', the difference is not parental origin. The power')
console.log('  numbers in the test above are what decide whether that has happened.')

console.log('')
console.log('parental-balance-real: the pronuclei are the only place this corpus answers the '
  + 'question without stacking a per-event attribution error on top of it.')

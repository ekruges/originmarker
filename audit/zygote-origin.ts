/**
 * THE FOUNDATIONAL FUNCTION, END TO END: does a zygote-stage genome get the right parent?
 *
 * WHY THIS EXISTS SEPARATELY FROM trio-calibration.ts. That harness measures the per-chromosome
 * Mendelian channel by calling it directly, and reports 0.8539 with zero wrong parents. But a real
 * zygote run does NOT use that channel: `scoreSample` switches it off on a uniparental sample,
 * because a genome carrying one parental contribution needs one statement about the genome rather
 * than twenty-two rows saying the same thing. What an operator actually sees for a pronucleus is the
 * GENOME-LEVEL call, `originClass` plus `zygosity`, and every event inherits its parent from that.
 *
 * So the number everyone quotes for the original use case was measured on a channel the original use
 * case does not run. This closes that gap by putting every pronucleus through the whole pipeline,
 * exactly as the browser and the command line do, and comparing the genome-level answer against the
 * dissection.
 *
 * GROUND TRUTH IS THE DISSECTION, NEVER THIS TOOL. GSE148488 labels each pronucleus maternal or
 * paternal in its GEO metadata, from which pronucleus was aspirated. A maternal pronucleus carries
 * the egg donor's genome and no paternal contribution, so the correct answer is gynogenetic; a
 * paternal one is androgenetic. The parental array is resolved independently by linkage, because
 * the manifest assigns the pronuclei no parents.
 *
 * BOTH ARMS, AND HERE THEY SHOULD AGREE. Each pronucleus is scored with its OWN donor loaded and
 * with the other side's donor loaded. Unlike a per-event parental call, the genome-level class is a
 * property of the GENOME and not of which array was loaded: loading the mother says "her genome is
 * all that is here" and loading the father says "he contributed nothing", and both mean
 * gynogenetic. So agreement between the arms is the correct outcome and disagreement is the
 * warning sign. Before the role fix the arms disagreed on 12 of 14, precisely because one of them
 * was inverted.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 \
 *        audit/zygote-origin.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('zygote-origin: OM_TRIOS is not set to a readable directory. NOT RUN.')
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

/** Genotypes only, for the linkage that resolves which donor a pronucleus came from. */
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
const who = (title: string) => (/sperm/i.test(title) ? 'sperm'
  : (/egg donor ([A-Z])/i.exec(title)?.[1] ?? title))

/**
 * Which donor a pronucleus came from, by linkage rather than from the sheet.
 *
 * The opposite-homozygote rate is Mendelianly forbidden between a true parent and its offspring, so
 * the real donor sits near zero and everyone else an order of magnitude above. Replicates of one
 * donor are collapsed first, or the runner-up is always the same person's second array.
 */
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
  const nextPerson = scored.find((x) => x.who !== best.who)
  if (!nextPerson || !(nextPerson.rate > best.rate * 2)) return null
  return { id: best.id, who: best.who }
}

const pronuclei = recs.filter((r) => r.pronucleus && has(r.gsm))
console.log(`pronuclei: ${pronuclei.length}, donor arrays: ${donors.length}`)
console.log('')
console.log('THE QUESTION: does the whole pipeline give a zygote-stage genome the right parent?')
console.log('Truth is the dissection: a MATERNAL pronucleus carries only the egg donor, so it is')
console.log('gynogenetic; a PATERNAL one is androgenetic. Nothing below comes from this tool.')
console.log('')

let correct = 0
let wrong = 0
let silent = 0
let unresolved = 0
let armCollapsed = 0
const wrongDetail: string[] = []

console.log('array        side      donor resolved   loaded    originClass        zygosity          verdict')
for (const p of pronuclei) {
  const res = resolveDonor(p.gsm)
  if (!res) { unresolved += 1; console.log(`${p.gsm}  ${p.pronucleus.padEnd(9)} NOT RESOLVED`); continue }
  // A maternal pronucleus carries the egg donor. Load that donor under its true role.
  const ownRole: 'paternal' | 'maternal' = res.who === 'sperm' ? 'paternal' : 'maternal'
  // And the other side: a donor who is a different person, under the opposite role.
  const otherId = donors.find((d) => who(d.title) !== res.who)!.gsm
  const otherRole: 'paternal' | 'maternal' = ownRole === 'paternal' ? 'maternal' : 'paternal'

  const runs: Record<string, { originClass: string; zygosity: string }> = {}
  for (const [arm, pid, role] of [
    ['own', res.id, ownRole], ['other', otherId, otherRole],
  ] as [string, string, 'paternal' | 'maternal'][]) {
    const pat = parentOf(pid)
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
    runs[arm] = {
      originClass: String(out.originClass ?? 'none'),
      zygosity: String(out.zygosity ?? 'none'),
    }
  }

  // THE ANSWER. A maternal pronucleus should read gynogenetic, a paternal one androgenetic.
  const want = p.pronucleus === 'maternal' ? 'gynogenetic' : 'androgenetic'
  const got = runs.own.originClass
  const verdict = got === want ? 'CORRECT'
    : (got === 'gynogenetic' || got === 'androgenetic') ? 'WRONG PARENT' : 'silent'
  if (verdict === 'CORRECT') correct += 1
  else if (verdict === 'WRONG PARENT') { wrong += 1; wrongDetail.push(`${p.gsm}: wanted ${want}, got ${got}`) }
  else silent += 1
  // The two arms should AGREE: see the header. Disagreement means the class moved with the loaded
  // array, which is what the role inversion did.
  if (runs.other.originClass !== runs.own.originClass) armCollapsed += 1

  console.log(`${p.gsm}  ${p.pronucleus.padEnd(9)} ${res.who.padEnd(6)} ${res.id}  `
    + `${ownRole.padEnd(9)} ${got.padEnd(18)} ${runs.own.zygosity.padEnd(17)} ${verdict}`
    + (runs.other.originClass === runs.own.originClass
      ? `   [other arm gave the SAME: ${runs.other.originClass}]` : ''))
}

const n = correct + wrong + silent
console.log('')
console.log(`resolved and scored: ${n}, unresolved by linkage: ${unresolved}`)
console.log(`  CORRECT parent at genome level  ${correct}/${n}`
  + (n ? ` (${(100 * correct / n).toFixed(0)}%)` : ''))
console.log(`  WRONG parent                    ${wrong}/${n}`)
console.log(`  silent, no uniparental call     ${silent}/${n}`)
console.log(`  the two arms DISAGREED on: ${armCollapsed}/${n}   `
  + '(agreement is correct here: the class belongs to the genome, not to the loaded array)')
for (const d of wrongDetail) console.log(`  ${d}`)
console.log('')
console.log(wrong === 0
  ? 'zygote-origin: the foundational call names no wrong parent through the whole pipeline.'
  : `zygote-origin: ${wrong} WRONG PARENT at genome level. This is the original use case.`)
process.exit(wrong ? 1 : 0)

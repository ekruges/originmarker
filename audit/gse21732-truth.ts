/**
 * Embryos of a balanced translocation carrier, scored against the karyotype published per sample.
 *
 * GSE21732 is a PGD case report: 48 embryos, both parents on their own arrays, and one ISCN
 * karyotype per sample in the series metadata, agreed by real-time PCR, SNP microarray and FISH.
 * The mother carries 46,XX,t(2;20)(q21;p12.2) and the father is 46,XY, so every derivative
 * chromosome in an embryo came from her. This is the only set in the corpus where a segmental
 * imbalance carries a stated parent of origin, and the only one whose answers were produced by a
 * published method rather than by a bench dissection.
 *
 * Scored from the ISCN string alone:
 *   whole chromosome  a +N or -N token, against the tool's own whole-chromosome call and direction
 *   derivative        a der(2) or der(20) token, against any imbalance reported on chr2 or chr20,
 *                     and against the parent named for it, which must be maternal
 *   balanced          a t(2;20) token with no derivative carries no imbalance, which is the
 *                     discriminating case this study exists for
 *   normal            46,XX or 46,XY alone is a negative array end to end
 *
 * A refusal is scored apart from a miss: an array whose stage inference failed reports nothing and
 * says so.
 *
 * Run: OM_T21732=<dir of .probes and truth.json> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/gse21732-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const sexing = await import(`${W}sexing.ts`)

const DIR = process.env.OM_T21732 as string
if (!DIR || !existsSync(join(DIR, 'truth.json'))) {
  console.log('gse21732-truth: OM_T21732 with truth.json is required. NOT RUN.')
  process.exit(0)
}
type T = { gsm: string; title: string } & Record<string, string>
const truth: T[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(DIR, `${t.gsm}.probes`)))
const FIELD = 'blastocyst snp microarray result'

interface Karyotype {
  sex: string | null
  gains: string[]
  losses: string[]
  derivatives: string[]
  balanced: boolean
  readable: boolean
}
/** The submitters' ISCN, parsed as written. Commas separate tokens; the t() notation carries none. */
function parse(raw: string | undefined): Karyotype {
  const empty = { sex: null, gains: [], losses: [], derivatives: [], balanced: false, readable: false }
  if (!raw || !/^\d{2}/.test(raw.trim())) return empty
  const k: Karyotype = { sex: null, gains: [], losses: [], derivatives: [], balanced: false, readable: true }
  for (const tok of raw.split(',').map((x) => x.trim()).filter(Boolean)) {
    if (/^\d+$/.test(tok)) continue
    if (/^X[XY]?$|^XXY$|^XYY$/.test(tok)) { k.sex = tok; continue }
    if (/der\(/.test(tok)) {
      for (const c of [...tok.matchAll(/der\((\d{1,2})\)/g)].map((m) => m[1])) k.derivatives.push(c)
      continue
    }
    if (/^t\(/.test(tok)) { k.balanced = true; continue }
    const m = /^([+-])(\d{1,2}|X|Y)$/.exec(tok)
    if (m) (m[1] === '+' ? k.gains : k.losses).push(m[2])
  }
  return k
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
function parentOf(gsm: string) {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR, `${gsm}.probes`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
}

const titled = (s: string) => truth.find((t) => t.title.toLowerCase().startsWith(s))
const mum = titled('maternal')
const dad = titled('paternal')
if (!mum || !dad) { console.log('both parental arrays are required. NOT RUN.'); process.exit(0) }
console.log(`parents: mother ${mum.gsm} "${mum[FIELD]}", father ${dad.gsm} "${dad[FIELD]}"`)
const pat = parentOf(dad.gsm)
const mat = parentOf(mum.gsm)

const wilson = (k: number, n: number): [number, number] => {
  if (!n) return [NaN, NaN]
  const z = 1.96; const p = k / n; const d = 1 + (z * z) / n
  const c = p + (z * z) / (2 * n)
  const h = z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))
  return [Math.max(0, (c - h) / d), Math.min(1, (c + h) / d)]
}
const pct = (x: number) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : '-')

const tally = {
  wholeListed: 0, wholeCalled: 0, wholeDirection: 0, wholeRefused: 0,
  derListed: 0, derCalled: 0, derMaternal: 0, derNamed: 0, derRefused: 0,
  normalN: 0, normalClean: 0, normalRefused: 0,
  balancedN: 0, balancedClean: 0, balancedRefused: 0,
  sexOk: 0, sexN: 0, threw: 0,
}
console.log('')
console.log('sample                       karyotype                                      tool')
for (const t of truth) {
  if (t.gsm === mum.gsm || t.gsm === dad.gsm) continue
  const k = parse(t[FIELD])
  if (!k.readable) continue
  try {
    const sx = sexing.emptySex()
    const acc = score.emptyCollected(pat as never, mat as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(join(DIR, `${t.gsm}.probes`))()) {
      if (!first) first = r.probesetId
      sexing.accumulateSex(r as never, sx as never)
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, pat as never, mat as never, acc as never)
    }
    const profile = ingest.finishProfile(t.gsm, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: t.gsm, log: () => {},
    }) as Record<string, any>
    const refused = res.stage?.stage === 'failed'
    const chroms = (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
    const segs = (res.segments ?? []) as { chrom: string; kind: string }[]
    const rows = (res.oneParent ?? []) as { where: string; parent?: string | null }[]
    const called = new Map(chroms.filter((c) => c.aneuploidy).map((c) => [c.chrom, c.aneuploidy as string]))
    const onChrom = (c: string) => called.has(c) || segs.some((s) => s.chrom === c)

    for (const [list, want] of [[k.gains, 'gain'], [k.losses, 'loss']] as const) {
      for (const c of list.filter((x) => /^\d+$/.test(x))) {
        if (refused) { tally.wholeRefused += 1; continue }
        tally.wholeListed += 1
        if (called.has(c)) {
          tally.wholeCalled += 1
          if (called.get(c) === want) tally.wholeDirection += 1
        }
      }
    }
    for (const c of k.derivatives) {
      if (refused) { tally.derRefused += 1; continue }
      tally.derListed += 1
      if (onChrom(c)) tally.derCalled += 1
      const named = rows.filter((r) => r.where?.includes(`chr${c}`) && r.parent)
      if (named.length) {
        tally.derNamed += 1
        if (named.every((r) => r.parent === 'maternal')) tally.derMaternal += 1
      }
    }
    const clean = called.size === 0 && segs.length === 0
    if (!k.gains.length && !k.losses.length && !k.derivatives.length) {
      if (k.balanced) {
        tally.balancedN += 1
        if (refused) tally.balancedRefused += 1
        else if (clean) tally.balancedClean += 1
      } else {
        tally.normalN += 1
        if (refused) tally.normalRefused += 1
        else if (clean) tally.normalClean += 1
      }
    }
    const y = (sexing.sexCall(sx as never) as { yBearing: boolean | null }).yBearing
    if (k.sex && y !== null) { tally.sexN += 1; if (y === k.sex.includes('Y')) tally.sexOk += 1 }

    const toolStr = refused ? 'REFUSED (stage failed)'
      : `${[...called].map(([c, a]) => `chr${c} ${a}`).join(', ') || 'no whole-chromosome call'}`
        + `${segs.length ? `; ${segs.length} segment(s) on ${[...new Set(segs.map((s) => `chr${s.chrom}`))].join(',')}` : ''}`
    console.log(`  ${t.gsm} ${t.title.slice(0, 22).padEnd(23)} ${String(t[FIELD]).slice(0, 44).padEnd(45)} ${toolStr}`)
  } catch (e) {
    tally.threw += 1
    console.log(`  THREW ${t.gsm}: ${String((e as Error).message ?? e).slice(0, 70)}`)
  }
}

const band = (k: number, n: number) => `[${pct(wilson(k, n)[0])}, ${pct(wilson(k, n)[1])}]`
console.log('')
console.log('=== WHOLE CHROMOSOMES the karyotype lists')
console.log(`  called on accepted arrays: ${tally.wholeCalled}/${tally.wholeListed} ${band(tally.wholeCalled, tally.wholeListed)}`)
console.log(`  direction agrees as well:  ${tally.wholeDirection}/${tally.wholeListed}`)
console.log(`  on arrays the tool refused: ${tally.wholeRefused}`)
console.log('')
console.log('=== DERIVATIVE CHROMOSOMES, which came from the mother by the parents\' own karyotypes')
console.log(`  imbalance reported on the derivative's chromosome: ${tally.derCalled}/${tally.derListed} ${band(tally.derCalled, tally.derListed)}`)
console.log(`  a parent named there: ${tally.derNamed}/${tally.derListed}, of which MATERNAL: ${tally.derMaternal}`)
console.log(`  on arrays the tool refused: ${tally.derRefused}`)
console.log('')
console.log('=== NEGATIVES')
console.log(`  normal karyotype, no call anywhere: ${tally.normalClean}/${tally.normalN - tally.normalRefused} (refused ${tally.normalRefused})`)
console.log(`  balanced carrier, no call anywhere: ${tally.balancedClean}/${tally.balancedN - tally.balancedRefused} (refused ${tally.balancedRefused})`)
console.log('')
console.log(`sex agrees with the karyotype: ${tally.sexOk}/${tally.sexN}   arrays that threw: ${tally.threw}`)

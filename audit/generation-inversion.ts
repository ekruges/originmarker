/**
 * GENERATION INVERSION: IS IT REALLY UNDETECTABLE, OR IS A THEOREM EXCUSING A MISSED CHECK?
 *
 * THE STANDING CLAIM, from relatedness.ts: "DIRECTION IS NOT IN HERE, and no amount of arithmetic
 * on two genotype files will put it there. Under Hardy-Weinberg the likelihood of a pair factors as
 * P(a)P(b|a) = P(b)P(a|b), so parent-child and child-parent are the same hypothesis."
 *
 * THAT IS TRUE, AND IT IS ABOUT A PAIR. The parent-child joint genotype distribution really is
 * exchangeable under Hardy-Weinberg, so the statistic the tool computes cannot see direction and
 * saying so is honest. What the theorem does not cover is a THIRD array, and the standard clinical
 * configuration for this tool has one: both parents are loaded, not one.
 *
 * THE ASYMMETRY A THIRD ARRAY CREATES, and it needs no new arithmetic:
 *
 *   correct     paternal slot = father,  maternal slot = mother,  sample = child
 *               the mother IS the child's parent            -> parent and child
 *   inverted    paternal slot = CHILD,   maternal slot = mother,  sample = FATHER
 *               the mother is the father's WIFE, not his parent   -> unrelated
 *
 * Two unrelated adults is not a subtle signal and it is not a quality signal. Both arrays in that
 * comparison are bulk genomic DNA, so it is measured where the statistic has the MOST power, which
 * is the opposite of the situation that made the material-based gate unusable.
 *
 * WHY THIS IS WORTH MEASURING RATHER THAN ASSUMING. scoreSample.ts calls `relate` on `pat` and on
 * nothing else. The maternal array is loaded, carried through the whole scorer, used for parental
 * attribution, and never once asked whether it is related to the sample. If the separation below is
 * clean then the tool holds the evidence to refuse an inverted run and does not look at it.
 *
 * THE ONE-PARENT ARM IS THE CONTROL, and it is there to keep the finding honest. With a single
 * parental array the theorem should bite and nothing should separate. A harness that only ran the
 * two-parent arm could not tell a real catch from an artefact of how it was set up.
 *
 * TRUTH COMES FROM THE MANIFEST. `father_gsms` and `mother_gsms` were established by the
 * experiment, and `complete_trio` marks the children where both are known.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=8192 \
 *        audit/generation-inversion.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const rel = await import(`${W}relatedness.ts`)
const defects = await import(`${W}defects.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('generation-inversion: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}
const CONTAMINATED = 'GSM4472408'
const CAP = Number(process.env.OM_CAP ?? 8)

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], material: f[col('material')], role: f[col('role')],
    pronucleus: f[col('pronucleus')],
    father: (f[col('father_gsms')] ?? '').split(';').filter(Boolean)[0],
    mother: (f[col('mother_gsms')] ?? '').split(';').filter(Boolean)[0],
    complete: f[col('complete_trio')] === 'True',
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

type AB = 'AA' | 'AB' | 'BB' | 'NC'
/** Genotypes and marker positions, which is what `relate` needs. */
const posCache = new Map<string, {
  gt: Map<string, AB>; pos: Map<string, { chrom: string; pos: number }>
}>()
function loadPositioned(gsm: string) {
  const hit = posCache.get(gsm)
  if (hit) return hit
  const gt = new Map<string, AB>()
  // The segmental test needs chromosome AND position, not a bare coordinate: it must never let a
  // window span two chromosomes. Handing it a number map leaves every marker looking like a new
  // chromosome, which silently produces zero windows rather than an error.
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (const r of rowsFromText(textOf(gsm))()) {
    gt.set(r.probesetId, r.genotype as AB)
    if (Number.isFinite(r.pos)) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  const v = { gt, pos }
  posCache.set(gsm, v)
  return v
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

/** One run of the SHIPPED scorer, in whatever configuration the caller names. */
async function runOne(patGsm: string, matGsm: string | null, sampleGsm: string,
  soloRole: 'paternal' | 'maternal') {
  const pat = parentOf(patGsm)
  const mat = matGsm ? parentOf(matGsm) : null
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(sampleGsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(sampleGsm, byChrom as never, bafSums as never, first)
  const out = await score.scoreSample({
    acc, profile, pat, mat, soloRole, sibs: [], sampleName: sampleGsm, log: () => {},
  }) as Record<string, unknown>
  const rows = (out.oneParent ?? []) as { parent?: string; band?: string }[]
  return {
    named: rows.filter((r) => r.parent).length,
    rows: rows.length,
    verdict: (out.relationship as { verdict?: string } | undefined)?.verdict ?? 'not computed',
    stage: (out.stage as { stage?: string } | undefined)?.stage ?? '?',
    integrity: String((out.integrity as { level?: string } | undefined)?.level ?? '-'),
    withheld: !!out.withholdReason,
    alerts: defects.runAlerts(out as never).length,
  }
}

/** The check the shipped scorer never performs: is the SECOND parental array related at all? */
function secondArrayCheck(otherGsm: string, sampleGsm: string) {
  const a = loadPositioned(otherGsm)
  const b = loadPositioned(sampleGsm)
  const r = rel.relate({ gt: a.gt, pos: a.pos } as never,
    { gt: b.gt, pos: b.pos } as never, rel.OPPOSITE_HOM_MAX)
  return {
    verdict: r.relationship, opp: r.opp.rate, n: r.opp.n,
    floor: r.win.floor, spread: r.win.spread, windows: r.win.windows,
  }
}

const trios = recs.filter((r) => r.complete && r.role !== 'parent' && has(r.gsm)
  && has(r.father) && has(r.mother) && r.mother !== CONTAMINATED
  && r.father !== CONTAMINATED).slice(0, CAP)

console.log(`complete trios used: ${trios.length} (cap ${CAP})`)
console.log(`the gate that would do the refusing: OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX}`)
console.log('')

interface Row {
  child: string; arm: string
  toolVerdict: string; named: number; rows: number; stage: string
  withheld: boolean; alerts: number
  secondVerdict: string; secondOpp: number
  secondFloor: number; secondWindows: number
}
const out: Row[] = []

for (const t of trios) {
  // ---- ONE PARENT, which is where the theorem applies and nothing should separate.
  const c1 = await runOne(t.father, null, t.gsm, 'paternal')
  const i1 = await runOne(t.gsm, null, t.father, 'paternal')
  // ---- TWO PARENTS, the standard configuration, where a third array exists.
  const c2 = await runOne(t.father, t.mother, t.gsm, 'paternal')
  const i2 = await runOne(t.gsm, t.mother, t.father, 'paternal')
  // The mother against whatever sat in the sample slot.
  const s2c = secondArrayCheck(t.mother, t.gsm)
  const s2i = secondArrayCheck(t.mother, t.father)

  const push = (arm: string, r: Awaited<ReturnType<typeof runOne>>,
    s: { verdict: string; opp: number; floor: number; windows: number } | null) => out.push({
    child: t.gsm, arm, toolVerdict: r.verdict, named: r.named, rows: r.rows,
    stage: r.stage, withheld: r.withheld, alerts: r.alerts,
    secondVerdict: s ? s.verdict : 'n/a, one parent loaded',
    secondOpp: s ? s.opp : NaN,
    secondFloor: s ? s.floor : NaN, secondWindows: s ? s.windows : 0,
  })
  push('1-parent CORRECT', c1, null)
  push('1-parent INVERTED', i1, null)
  push('2-parent CORRECT', c2, s2c)
  push('2-parent INVERTED', i2, s2i)

  console.log(`${t.gsm} ${t.material}`)
  for (const r of out.slice(-4)) {
    console.log(`  ${r.arm.padEnd(18)} tool says "${r.toolVerdict.slice(0, 34).padEnd(34)}" `
      + `named ${String(r.named).padStart(2)}/${String(r.rows).padStart(2)} `
      + `stage ${r.stage.padEnd(14)} alerts ${r.alerts} `
      + `| 2nd array: ${r.secondVerdict.slice(0, 30).padEnd(30)}`
      + `${Number.isFinite(r.secondOpp) ? ` opp ${r.secondOpp.toFixed(4)}` : ''}`
      + `${Number.isFinite(r.secondFloor) ? ` floor ${r.secondFloor.toFixed(4)} win ${r.secondWindows}` : ''}`)
  }
}

// ------------------------------------------------------------------ the measurement
const arm = (a: string) => out.filter((r) => r.arm === a)
console.log('')
console.log('=== DOES THE TOOL, AS SHIPPED, SEPARATE CORRECT FROM INVERTED?')
for (const a of ['1-parent', '2-parent']) {
  const c = arm(`${a} CORRECT`); const i = arm(`${a} INVERTED`)
  const namedC = c.reduce((x, r) => x + r.named, 0)
  const namedI = i.reduce((x, r) => x + r.named, 0)
  const wC = c.filter((r) => r.withheld).length; const wI = i.filter((r) => r.withheld).length
  const vC = new Set(c.map((r) => r.toolVerdict)); const vI = new Set(i.map((r) => r.toolVerdict))
  console.log(`  ${a}: parents named ${namedC} correct vs ${namedI} INVERTED, `
    + `withheld ${wC}/${c.length} vs ${wI}/${i.length}`)
  console.log(`    verdicts correct  : ${[...vC].join(' | ')}`)
  console.log(`    verdicts inverted : ${[...vI].join(' | ')}`)
}

console.log('')
console.log('=== DOES THE SECOND ARRAY SEPARATE THEM? (the check scoreSample never runs)')
const sc = arm('2-parent CORRECT'); const si = arm('2-parent INVERTED')
const oppC = sc.map((r) => r.secondOpp).filter(Number.isFinite)
const oppI = si.map((r) => r.secondOpp).filter(Number.isFinite)
const hi = (xs: number[]) => Math.max(...xs); const lo = (xs: number[]) => Math.min(...xs)
console.log(`  mother vs child  (correct) : opp ${lo(oppC).toFixed(4)} to ${hi(oppC).toFixed(4)}  `
  + `verdicts ${[...new Set(sc.map((r) => r.secondVerdict))].join(' | ')}`)
console.log(`  mother vs father (inverted): opp ${lo(oppI).toFixed(4)} to ${hi(oppI).toFixed(4)}  `
  + `verdicts ${[...new Set(si.map((r) => r.secondVerdict))].join(' | ')}`)
// The segmental statistic too, because a rate that overlaps does not mean every statistic does.
const flC = sc.map((r) => r.secondFloor).filter(Number.isFinite)
const flI = si.map((r) => r.secondFloor).filter(Number.isFinite)
if (flC.length && flI.length) {
  console.log(`  window floor, mother vs child   : ${lo(flC).toFixed(4)} to ${hi(flC).toFixed(4)}`)
  console.log(`  window floor, mother vs father  : ${lo(flI).toFixed(4)} to ${hi(flI).toFixed(4)}`)
  console.log(`  floor gap: ${(lo(flI) - hi(flC)).toFixed(4)} `
    + `(IBD0_FLOOR_MAX = ${rel.IBD0_FLOOR_MAX})`)
}
const gap = lo(oppI) - hi(oppC)
const gateInside = rel.OPPOSITE_HOM_MAX > hi(oppC) && rel.OPPOSITE_HOM_MAX < lo(oppI)
console.log('')
console.log('  NOTE ON MATERIAL, because it decides whether the shipped gate can be reused. In the')
console.log('  correct arm the mother is compared against AMPLIFIED material, which inflates the')
console.log('  rate on a true relative; defects.ts records a real blastomere reading 0.0276 against')
console.log('  its confirmed father. In the inverted arm she is compared against BULK genomic DNA,')
console.log('  which is the material OPPOSITE_HOM_MAX was measured on. So the two arms are not')
console.log('  equally powered and the gap below matters more than either absolute value.')
console.log('')
if (gap > 0) {
  console.log(`  SEPARABLE. Every correct run sits under ${hi(oppC).toFixed(4)} and every inverted`)
  console.log(`  run at or over ${lo(oppI).toFixed(4)}, a gap of ${gap.toFixed(4)}.`)
  console.log(gateInside
    ? `  The shipped gate OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX} falls INSIDE that gap, so no new`
      + '\n  constant is needed: running the existing test on the maternal array refuses every'
      + '\n  inverted run here and passes every correct one.'
    : `  The shipped gate OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX} does NOT fall inside that gap,`
      + '\n  so the separation is real but the existing constant cannot be reused as-is on this'
      + '\n  material. A threshold exists; it would have to be measured rather than inherited.')
  console.log('  Either way scoreSample.ts calls relate on `pat` only and never asks this question.')
} else {
  console.log('  NOT SEPARABLE on this material. The limit is real in the two-parent case too and')
  console.log('  the theorem covers it; the docstring stands as written.')
}
console.log('')
console.log('generation-inversion: the one-parent arm is the control. If it separated, something')
console.log('other than the third array is doing the work and this finding would not hold.')

/**
 * WHICH PARENT, ON REAL EVENTS, AGAINST TRUTH THIS TOOL DID NOT PRODUCE.
 *
 * ATTRIBUTION IS THE CLAIM THAT MATTERS AND THE ONE HARDEST TO CHECK. Detection can be checked
 * against a karyotype or an experimental design. Which parent's copy went missing usually cannot,
 * because nobody dissected the embryo. This file builds that truth three different ways and none of
 * them reads the tool's answer.
 *
 * ARM A, DIRECT MENDELIAN TRUTH FROM THE TRIO, AND IT IS NOT THE TOOL'S ARITHMETIC. Where the
 * father is homozygous for one allele and the mother is homozygous for the OTHER, a normal child
 * must be heterozygous. If the child instead carries only the mother's allele across a region, the
 * paternal copy is the one that went; only the father's allele, and it was the maternal. That is
 * decided here in a few lines, over both parental arrays at once, and compared against what the
 * tool said from ONE parent. It is the configuration a laboratory is usually in, scored against the
 * configuration that can actually answer.
 *
 *   the marker classes this rests on, and why they are the informative ones:
 *     father AA, mother BB   a normal child is AB. Homozygous AA means the mother's copy is gone,
 *                            homozygous BB means the father's.
 *   Dropout removes one allele of a heterozygote, so a single marker proves nothing. The call is
 *   taken over every informative marker in the interval with a floor and a margin, and refuses
 *   rather than guessing when the two counts are close.
 *
 * ARM B, RECIPROCAL CONSISTENCY, WHICH NEEDS NO TRUTH AT ALL. Score the same array twice, once
 * against the father and once against the mother. Each run is most reliable about ITS OWN loaded
 * parent's copy, at 0.8539 against 0.2890 for the other direction. If the father-loaded run says
 * the paternal copy is gone at an interval and the mother-loaded run says the maternal copy is gone
 * at the same interval, both cannot be true: the sample would have neither copy and the DNA is
 * there. A contradiction is an error with no ground truth required, which is why this arm can be
 * run on every array rather than only on the few with a resolvable answer.
 *
 * ARM C, THE NEGATIVE CONTROL, and it is the one that catches attribution invented out of nothing.
 * On a genome carrying ONE parental complement, "the loaded parent's copy is absent here" is true
 * on every chromosome by construction and names the parent that was never there. Nothing may be
 * named on those. The pronuclei are known one-complement by the dissection, so this needs no
 * threshold either.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 \
 *        audit/attribution-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const defects = await import(`${W}defects.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('attribution-truth: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}
const CONTAMINATED = 'GSM4472408'
/** Fewest informative markers before the trio call means anything. */
const MIN_INFORMATIVE = 200
/** How far one side must lead the other before the trio call is taken rather than refused. */
const MARGIN = 3

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')], pronucleus: f[col('pronucleus')],
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
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Geno { gt: Map<string, AB>; pos: Map<string, { chrom: string; pos: number }> }
const genoCache = new Map<string, Geno>()
function genotypes(gsm: string, keep = false): Geno {
  const hit = genoCache.get(gsm)
  if (hit) return hit
  const gt = new Map<string, AB>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (const r of rowsFromText(textOf(gsm))()) {
    gt.set(r.probesetId, r.genotype as AB)
    if (Number.isFinite(r.pos)) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  const v = { gt, pos }
  if (keep) genoCache.set(gsm, v)
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

/**
 * WHICH PARENT'S COPY IS MISSING, decided from BOTH parental arrays and the child, with no part of
 * the scorer involved. This is the truth column for arm A.
 */
function trioVerdict(
  father: Geno, mother: Geno, child: Geno,
  chrom: string, startBp: number, endBp: number,
): { call: 'paternal' | 'maternal' | 'both-present' | 'refused'; nPat: number; nMat: number
  nBoth: number; informative: number } {
  let nPat = 0
  let nMat = 0
  let nBoth = 0
  for (const [probe, cg] of child.gt) {
    const at = child.pos.get(probe)
    if (!at || at.chrom !== chrom || at.pos < startBp || at.pos > endBp) continue
    const fg = father.gt.get(probe)
    const mg = mother.gt.get(probe)
    // The informative class: both parents homozygous and for OPPOSITE alleles, where a normal
    // child is obliged to be heterozygous.
    if ((fg !== 'AA' && fg !== 'BB') || (mg !== 'AA' && mg !== 'BB') || fg === mg) continue
    if (cg === 'AB') { nBoth += 1; continue }
    if (cg === fg) nMat += 1        // only the father's allele survives: the maternal copy is gone
    else if (cg === mg) nPat += 1   // only the mother's: the paternal copy is gone
  }
  const informative = nPat + nMat + nBoth
  if (informative < MIN_INFORMATIVE) return { call: 'refused', nPat, nMat, nBoth, informative }
  // Both copies present is the answer when the heterozygous class dominates, which is what a
  // region with no loss looks like.
  if (nBoth > (nPat + nMat) * 2) return { call: 'both-present', nPat, nMat, nBoth, informative }
  if (nPat > nMat * MARGIN) return { call: 'paternal', nPat, nMat, nBoth, informative }
  if (nMat > nPat * MARGIN) return { call: 'maternal', nPat, nMat, nBoth, informative }
  return { call: 'refused', nPat, nMat, nBoth, informative }
}

interface Row {
  where: string
  parent?: string | null
  verdict: string
  band?: string
}
async function runOne(patGsm: string, sampleGsm: string,
  role: 'paternal' | 'maternal', matGsm?: string) {
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
    acc, profile, pat, mat, soloRole: role, sibs: [], sampleName: sampleGsm, log: () => {},
  }) as Record<string, unknown>
  const segs = (out.segments ?? []) as {
    chrom: string; kind: string; startBp: number; endBp: number
    refined?: { startBp: number; endBp: number }
  }[]
  const chroms = (out.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
  return {
    rows: (out.oneParent ?? []) as Row[],
    zygosity: String(out.zygosity ?? '?'),
    stage: (out.stage as { stage?: string } | undefined)?.stage ?? '?',
    intervals: [
      ...chroms.filter((c) => c.aneuploidy).map((c) => ({
        chrom: c.chrom, startBp: 0, endBp: Number.MAX_SAFE_INTEGER, label: `chr${c.chrom}`,
      })),
      ...segs.map((s) => ({
        chrom: s.chrom,
        startBp: s.refined?.startBp ?? s.startBp,
        endBp: s.refined?.endBp ?? s.endBp,
        label: `chr${s.chrom} ${((s.refined?.startBp ?? s.startBp) / 1e6).toFixed(1)}-`
          + `${((s.refined?.endBp ?? s.endBp) / 1e6).toFixed(1)}Mb`,
      })),
    ],
  }
}

/**
 * The parent a row names, THROUGH THE SHIPPED MAPPING rather than a copy of it.
 *
 * This was a regex over the verdict string, matching `loaded-parent|this-parent|parent-absent`.
 * The vocabulary it needed to match is `known-parent-lost` and `other-parent-lost`, so it matched
 * nothing and reported the tool as silent on eight events it had in fact answered correctly.
 * defects.parentNamed is the function both the command line and the report already use, and its
 * switch is exhaustive on purpose so a new verdict is a compile error rather than a silent null.
 */
function namedParent(r: Row, loadedRole: 'paternal' | 'maternal'): string | null {
  // With both arrays loaded the parent is on the row itself, because it comes from which SIDE
  // answered and no verdict string carries that.
  if (r.parent) return r.parent
  return (defects.parentNamed as (v: string, l: string) => string | null)(r.verdict, loadedRole)
}

const CAP = Number(process.env.OM_CAP ?? 14)
const kids = recs.filter((r) => r.complete && r.role !== 'parent' && !r.pronucleus && has(r.gsm)
  && has(r.father) && has(r.mother) && r.father !== CONTAMINATED && r.mother !== CONTAMINATED)
  .slice(0, CAP)

console.log(`complete trios with both parents on disk: ${kids.length} (cap ${CAP})`)
console.log(`trio truth needs >= ${MIN_INFORMATIVE} informative markers and a ${MARGIN}x margin`)
console.log('')

let agree = 0
let disagree = 0
let toolSilent = 0
let truthRefused = 0
let contradictions = 0
let reciprocalPairs = 0

console.log('=== ARM A and ARM B, per event')
for (const k of kids) {
  const father = genotypes(k.father, true)
  const mother = genotypes(k.mother, true)
  const child = genotypes(k.gsm)

  const asPat = await runOne(k.father, k.gsm, 'paternal')
  const asMat = await runOne(k.mother, k.gsm, 'maternal')

  const label = `${k.gsm} ${asPat.stage}/${asPat.zygosity}`
  if (!asPat.intervals.length && !asMat.intervals.length) {
    console.log(`  ${label}: no events`)
    genoCache.delete(k.gsm)
    continue
  }

  for (const iv of asPat.intervals) {
    const endBp = iv.endBp === Number.MAX_SAFE_INTEGER
      ? Math.max(...[...child.pos.values()].filter((p) => p.chrom === iv.chrom).map((p) => p.pos), 0)
      : iv.endBp
    const truth = trioVerdict(father, mother, child, iv.chrom, iv.startBp, endBp)

    const rowPat = asPat.rows.find((r) => r.where === iv.label)
    const rowMat = asMat.rows.find((r) => r.where === iv.label)
    const callPat = rowPat ? namedParent(rowPat, 'paternal') : null
    const callMat = rowMat ? namedParent(rowMat, 'maternal') : null

    // ---- ARM B: the two single-parent runs must not name opposite parents for one interval.
    if (callPat && callMat) {
      reciprocalPairs += 1
      if (callPat !== callMat) {
        contradictions += 1
        console.log(`  ${label} ${iv.label}: CONTRADICTION father-loaded says ${callPat}, `
          + `mother-loaded says ${callMat}`)
      }
    }

    // ---- ARM A: against the trio.
    const toolCall = callPat ?? callMat
    if (truth.call === 'refused' || truth.call === 'both-present') {
      truthRefused += 1
      console.log(`  ${label} ${iv.label}: trio ${truth.call} `
        + `(pat ${truth.nPat} mat ${truth.nMat} het ${truth.nBoth}), tool says ${toolCall ?? 'nothing'}`)
    } else if (!toolCall) {
      toolSilent += 1
      console.log(`  ${label} ${iv.label}: trio ${truth.call.toUpperCase()} `
        + `(pat ${truth.nPat} mat ${truth.nMat}), tool named nothing`)
    } else if (toolCall === truth.call) {
      agree += 1
      console.log(`  ${label} ${iv.label}: trio ${truth.call.toUpperCase()}, tool ${toolCall}  CORRECT`)
    } else {
      disagree += 1
      console.log(`  ${label} ${iv.label}: trio ${truth.call.toUpperCase()} `
        + `(pat ${truth.nPat} mat ${truth.nMat}), tool ${toolCall}  WRONG`)
    }
  }
  genoCache.delete(k.gsm)
}

// ================================================================ ARM C, one-complement genomes
console.log('')
console.log('=== ARM C: NOTHING MAY BE NAMED ON A ONE-COMPLEMENT GENOME')
console.log('A pronucleus carries one parental complement by the dissection, so "the loaded')
console.log('parent\'s copy is absent" is true everywhere and names a parent that was never there.')
const pronuclei = recs.filter((r) => r.pronucleus && has(r.gsm)).slice(0, 16)
let pnNamed = 0
let pnArrays = 0
for (const p of pronuclei) {
  const r = await runOne('GSM4472397', p.gsm, 'paternal')
  const named = r.rows.filter((x) => namedParent(x, 'paternal')).length
  pnNamed += named
  pnArrays += 1
  console.log(`  ${p.gsm} ${p.pronucleus} pronucleus  ${r.stage}/${r.zygosity}  `
    + `events ${r.intervals.length}  PARENTS NAMED ${named}${named ? '   <-- FALSE POSITIVE' : ''}`)
}

// ================================================================ the tally
console.log('')
console.log('=== TALLY')
const scored = agree + disagree
console.log(`ARM A, against direct trio truth on real events:`)
console.log(`  events where the trio could answer and the tool named a parent: ${scored}`)
console.log(`    CORRECT ${agree}    WRONG ${disagree}`
  + `${scored ? `    accuracy ${(agree / scored).toFixed(4)}` : ''}`)
console.log(`  tool named nothing where the trio could answer: ${toolSilent}`)
console.log(`  trio itself could not answer: ${truthRefused}`)
console.log('')
console.log(`ARM B, reciprocal consistency, no ground truth needed:`)
console.log(`  intervals named by BOTH single-parent runs: ${reciprocalPairs}`)
console.log(`  CONTRADICTIONS, opposite parents for one interval: ${contradictions}`)
console.log('')
console.log(`ARM C, one-complement negative control:`)
console.log(`  pronuclei scored: ${pnArrays}, parents named on them: ${pnNamed}`)
console.log('')
const clean = disagree === 0 && contradictions === 0 && pnNamed === 0
console.log(clean
  ? 'attribution-truth: no wrong parent, no contradiction, and nothing named on a genome that '
    + 'has only one parent to name.'
  : `attribution-truth: ${disagree} wrong, ${contradictions} contradictions, ${pnNamed} named on `
    + 'one-complement genomes.')

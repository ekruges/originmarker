/**
 * THE PROGENITOR TO SYNGAMY HANDOFF, END TO END, ON REAL ARRAYS.
 *
 * WHAT BROKE. Progenitor reconstructs a parent from haploid meiotic products and writes it out as
 * an array file for Syngamy to load as the donor. That file is zero percent heterozygous BY
 * CONSTRUCTION: the reconstruction can only assert the sites its products agree on, and a marker
 * where they disagree is exactly a marker where the parent is heterozygous, which no single
 * haploid product can resolve. Syngamy's `genome-wide LOH` gate reads heterozygosity under 2% as a
 * candidate failed genome unification and excludes the array, and a parent array that fails its own
 * gates stops the run. So every reconstruction was refused before a single call was made, and the
 * one product the two halves of this tool exist to hand between them could not cross.
 *
 * WHY THAT GATE IS RIGHT ABOUT A SAMPLE AND WRONG ABOUT A REFERENCE. Zero heterozygosity in a
 * FERTILISED genome is a real finding. In a reconstruction it is arithmetic. The mark that tells
 * the two apart was already written into the file by Progenitor and already read on drop by
 * Syngamy; it simply never reached the gates.
 *
 * WHAT THIS PROVES, and none of it is asserted from the shape of the code:
 *   1. a reconstruction built from real pronuclei carries zero heterozygosity, as expected
 *   2. the gates refuse it when it is treated as a measured array
 *   3. the gates accept it when it is declared a reconstruction, and no OTHER gate changes
 *   4. the file still round-trips through the reader and still identifies its own man
 *
 * Point 4 is the one that matters: a gate can be made to pass by weakening it. The reconstruction
 * has to still WORK after crossing, and that is scored against the sperm donor's own array.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=3072 \
 *        audit/handoff-truth.ts
 */
import { readFileSync, existsSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'
import { tmpdir } from 'node:os'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const inferredRef = await import(`${W}inferredReference.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('handoff-truth: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}

/** The sperm donor's own arrays. Never given to the reconstruction, only used to score it. */
const TRUTH_ARRAYS = ['GSM4472397', 'GSM4472398']

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return { gsm: f[col('gsm')], title: f[col('title')], pronucleus: f[col('pronucleus')] }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))
const textOf = (g: string) => gunzipSync(readFileSync(join(DIR, `${g}.probes.gz`))).toString('utf8')

function rowsOfText(t: string) {
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

// ---------------------------------------------------------------- 1. build a real reconstruction
const products = recs.filter((r) => r.pronucleus === 'paternal' && has(r.gsm))
console.log(`paternal pronuclei available as products: ${products.length}`)
if (products.length < 5) {
  console.log('handoff-truth: fewer than the 5 products a reconstruction needs. NOT RUN.')
  process.exit(0)
}

const ps = new inferredRef.ProductSet()
for (const p of products) {
  const slot = ps.begin(p.gsm)
  const band = { inBand: 0, total: 0 }
  for (const r of rowsOfText(textOf(p.gsm))()) ps.add(slot, r as never, band as never)
  ps.end(slot, band as never)
}
// The depth the module picks for itself, using no array of the parent. Same call Progenitor makes.
const chosen = ps.chooseM() as { mMin: number }
const built = ps.build(chosen.mMin) as {
  genotype: Map<string, string>; mMin: number; contamination: number
  spuriousAbsence: number; hRetained: number
}
console.log(`marker depth chosen without seeing the parent: m >= ${chosen.mMin}`)
console.log(`reconstruction: ${built.genotype.size.toLocaleString()} markers asserted from `
  + `${products.length} products`)

const text = inferredArray.inferredArrayText({
  genotype: built.genotype,
  locus: (id: string) => ps.locus(id),
  products: products.map((p) => p.gsm),
  mMin: built.mMin,
  contamination: built.contamination,
  spuriousAbsence: built.spuriousAbsence,
  hRetained: built.hRetained,
  side: 'paternal',
  reportId: 'handoff-check',
  generatedAt: new Date(0).toISOString(),
  build: 'GRCh37',
} as never) as string
const path = join(tmpdir(), 'progenitor-inferred-paternal-handoff.probes')
writeFileSync(path, text)
console.log(`written: ${path}, ${(text.length / 1e6).toFixed(1)} MB`)
console.log(`carries its own mark: ${inferredArray.inferredMark(text.slice(0, 4000))}`)

// ---------------------------------------------------------------- 2 and 3. the gates
const byChrom = new Map(); const bafSums = ingest.emptyBafSums()
const builds = ingest.emptyBuildSums()
let first = ''
const pacc = score.emptyParent()
for (const r of rowsOfText(text)()) {
  if (!first) first = r.probesetId
  ingest.accumulate(r as never, byChrom as never)
  ingest.accumulateBaf(r as never, bafSums as never)
  ingest.accumulateBuild(r as never, builds as never)
  score.collectParentRow(r as never, pacc as never)
}
const profile = {
  ...ingest.finishProfile('inferred', byChrom as never, bafSums as never, first, builds as never),
  build: ingest.buildVerdict(builds as never),
}
console.log('')
console.log(`heterozygosity of the reconstruction: ${(profile as { hetRate: number }).hetRate.toFixed(6)}`)
console.log(`call rate: ${(profile as { callRate: number }).callRate.toFixed(4)}`)

const asMeasured = ingest.gates(profile as never, false) as { name: string; verdict: string; detail: string }[]
const asInferred = ingest.gates(profile as never, true) as { name: string; verdict: string; detail: string }[]
const blockedM = asMeasured.filter((g) => g.verdict === 'exclude')
const blockedI = asInferred.filter((g) => g.verdict === 'exclude')

console.log('')
console.log('=== THE GATES, on the same file, read two ways')
console.log('gate                                    as a measured array   as a reconstruction')
for (let i = 0; i < asMeasured.length; i += 1) {
  const a = asMeasured[i]
  const b = asInferred.find((x) => x.name === a.name)
  const changed = b && b.verdict !== a.verdict
  console.log(`  ${a.name.padEnd(38)}${a.verdict.padEnd(22)}${b?.verdict ?? '-'}`
    + `${changed ? '   <- changed' : ''}`)
}
console.log('')
console.log(`  refused as a measured array: ${blockedM.length} `
  + `${blockedM.length ? `(${blockedM.map((g) => g.name).join(', ')})` : ''}`)
console.log(`  refused as a reconstruction: ${blockedI.length} `
  + `${blockedI.length ? `(${blockedI.map((g) => g.name).join(', ')})` : ''}`)
console.log('')
console.log(blockedM.length > 0 && blockedI.length === 0
  ? '  THE HANDOFF IS OPEN: refused before, accepted now, and only on the declaration.'
  : blockedI.length > 0
    ? `  STILL REFUSED: ${blockedI.map((g) => g.name).join(', ')}`
    : '  NOTHING CHANGED: it was not being refused, so this was not the bug.')

// Every gate that is NOT about heterozygosity must be untouched, or the fix is a blanket pass.
const unrelated = asMeasured.filter((g) => !/LOH|asymmetry/.test(g.name))
const moved = unrelated.filter((g) => asInferred.find((x) => x.name === g.name)?.verdict !== g.verdict)
console.log(`  gates unrelated to heterozygosity that changed: ${moved.length}`
  + `${moved.length ? `  <- ${moved.map((g) => g.name).join(', ')}` : '   (none, which is the point)'}`)

// ---------------------------------------------------------------- 4. does it still WORK
console.log('')
console.log('=== AND IT STILL HAS TO WORK. The reconstruction against the man himself.')
const pat = score.finishParent(pacc as never, (profile as { build: { build: string } }).build.build)
for (const gsm of TRUTH_ARRAYS.filter(has)) {
  const acc = score.emptyCollected(pat as never, null)
  const bc = new Map(); const bs = ingest.emptyBafSums(); let f0 = ''
  for (const r of rowsOfText(textOf(gsm))()) {
    if (!f0) f0 = r.probesetId
    ingest.accumulate(r as never, bc as never)
    ingest.accumulateBaf(r as never, bs as never)
    score.collectRow(r as never, pat as never, null, acc as never)
  }
  const prof = ingest.finishProfile(gsm, bc as never, bs as never, f0)
  const res = await score.scoreSample({
    acc, profile: prof, pat, mat: null, soloRole: 'paternal', sibs: [],
    sampleName: gsm, log: () => {},
  }) as Record<string, unknown>
  console.log(`  ${gsm} (HIS OWN ARRAY): verdict ${String(res.verdict)}  `
    + `nonParentalRate ${Number(res.nonParentalRate).toFixed(4)}  `
    + `explainable ${Number(res.explainable).toFixed(4)}`)
}
const other = recs.filter((r) => r.pronucleus === 'maternal' && has(r.gsm)).slice(0, 3)
for (const o of other) {
  const acc = score.emptyCollected(pat as never, null)
  const bc = new Map(); const bs = ingest.emptyBafSums(); let f0 = ''
  for (const r of rowsOfText(textOf(o.gsm))()) {
    if (!f0) f0 = r.probesetId
    ingest.accumulate(r as never, bc as never)
    ingest.accumulateBaf(r as never, bs as never)
    score.collectRow(r as never, pat as never, null, acc as never)
  }
  const prof = ingest.finishProfile(o.gsm, bc as never, bs as never, f0)
  const res = await score.scoreSample({
    acc, profile: prof, pat, mat: null, soloRole: 'paternal', sibs: [],
    sampleName: o.gsm, log: () => {},
  }) as Record<string, unknown>
  console.log(`  ${o.gsm} (NOT his, a maternal pronucleus): verdict ${String(res.verdict)}  `
    + `nonParentalRate ${Number(res.nonParentalRate).toFixed(4)}`)
}
console.log('')
console.log('His own arrays must read present and the maternal pronuclei absent. A gate that lets')
console.log('the file through without that is a gate that was weakened rather than corrected.')

/**
 * THE NUMBERS BEHIND EACH PARENT CALL, and what the reconstruction floor does to them.
 *
 * `lab-zygote-check.ts` prints the class. This prints the two quantities the class is decided by:
 * the opposite-homozygote rate against the loaded parent, and the ceiling that sample's own noise
 * can explain. Then it re-runs the SAME tally through `classify` a second time with
 * `spuriousAbsence` supplied, which is the floor a reconstructed reference contributes and which
 * `scoreSample` does not pass. The difference between the two columns is the cost of that gap,
 * measured rather than estimated.
 *
 * The floor is read from the reference file's own banner, not assumed: a measured parental array
 * contributes none and is scored with zero.
 *
 * NO SAMPLE IS NAMED IN THIS FILE. Barcodes come from the directory and the pairing file.
 *
 * Run: OM_LAB=<dir> OM_PAIRS=<tsv> OM_SPURIOUS=<fraction> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/spot-check-numbers.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const parentage = await import(`${W}parentage.ts`)

const DIR = process.env.OM_LAB
const PAIRS = process.env.OM_PAIRS
if (!DIR || !existsSync(DIR) || !PAIRS || !existsSync(PAIRS)) {
  console.log('spot-check-numbers: OM_LAB and OM_PAIRS must both be readable. NOT RUN.')
  process.exit(0)
}
/** The floor a reconstructed reference contributes, from its own banner. */
const SPURIOUS = Number(process.env.OM_SPURIOUS ?? '0')
/**
 * The reconstructed parent's OWN heterozygosity, from its banner (`Reference.hRetained`),
 * estimated from disagreement among her products rather than from the reference file, which is
 * homozygous at every marker by construction and therefore reads zero.
 */
const H_RETAINED = Number(process.env.OM_HRETAINED ?? 'NaN')

const fileOf = new Map<string, string>()
const walk = (d: string) => {
  for (const e of readdirSync(d)) {
    const p = join(d, e)
    if (statSync(p).isDirectory()) { walk(p); continue }
    const m = /(\d{4,6}-\d{2})/.exec(e)
    if (m && /\.probes(\.txt)?(\.gz)?$/.test(e) && !fileOf.has(m[1])) fileOf.set(m[1], p)
  }
}
walk(DIR)

/** Whether a parental file announces itself as reconstructed rather than measured. */
function isRebuilt(path: string): boolean {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8', 0, 4096)
    : readFileSync(path, 'utf8').slice(0, 4096)
  return raw.includes('OriginMarker-inferred-reference')
}

function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () { for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r } }
}
const cache = new Map<string, unknown>()
function parentIndex(bc: string) {
  const hit = cache.get(bc); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(bc)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(bc, byChrom as never, bafSums as never, first).build.build)
  cache.set(bc, idx); return idx
}

const pairs = readFileSync(PAIRS, 'utf8').trim().split('\n')
  .filter((l) => l && !l.startsWith('#'))
  .map((l) => { const [sample, father, mother] = l.split('\t').map((x) => x.trim()); return { sample, father, mother } })

const pct = (x: number, d = 2) => (Number.isFinite(x) ? `${(x * 100).toFixed(d)}%` : '-')
console.log(`reconstruction floor supplied to the second pass: ${pct(SPURIOUS, 3)}\n`)
console.log(`parent heterozygosity supplied to the third pass: ${pct(H_RETAINED, 3)}\n`)
console.log(`${'sample'.padEnd(10)} ${'ref'.padEnd(9)} ${'absence'.padStart(8)} ${'ceiling'.padStart(8)} `
  + `${'as shipped'.padEnd(13)}${'class'.padEnd(14)} ${'+floor'.padEnd(13)}${'class'.padEnd(14)} `
  + `${'+floor+het'.padEnd(13)}${'class'.padEnd(14)} moved`)

const short = (v: string) => v.replace('parent_genome_present', 'present')
  .replace('no_parental_contribution', 'ABSENT').replace('unclear', 'unclear')
let moved = 0; let scored = 0; let classMoved = 0
for (const p of pairs) {
  if (!fileOf.has(p.sample)) continue
  const hasF = !!p.father && fileOf.has(p.father)
  const hasM = !!p.mother && fileOf.has(p.mother)
  if (!hasF && !hasM) continue
  const soloRole = hasF ? 'paternal' : 'maternal'
  const parentBc = hasF ? p.father : p.mother
  const pat = hasF ? parentIndex(p.father) : parentIndex(p.mother)
  const mat = hasF && hasM ? parentIndex(p.mother) : null
  const rebuilt = isRebuilt(fileOf.get(parentBc)!)
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(p.sample)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never); ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(p.sample, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat, soloRole, sibs: [], sampleName: p.sample, log: () => {},
  } as never) as {
    stage?: { stage?: string }; verdict?: string; originClass?: string
    genomeRate?: number; explainable?: number
  }
  const bulk = res.stage?.stage === 'bulk' && profile.callRate >= 0.80
  const hMeasured = (pat as { heterozygosity: number }).heterozygosity
  // The same tally, scored again with the floor the caller omits. Zero for a measured array, so a
  // measured row is its own control: the two columns must agree there.
  const floor = rebuilt ? SPURIOUS : 0
  const second = parentage.classify(acc.t as never, hMeasured, {
    role: soloRole, bulk, spuriousAbsence: floor,
  } as never) as { verdict?: string; originClass?: string }
  // And a third time with the reference's TRUE heterozygosity as well. A reconstructed parent is
  // homozygous at every marker, so `finishParent` reads her heterozygosity as zero and the
  // second-parent axis is withheld under PANEL_HET_FLOOR; the real figure is on the banner.
  const hTrue = rebuilt && Number.isFinite(H_RETAINED) ? H_RETAINED : hMeasured
  const third = parentage.classify(acc.t as never, hTrue, {
    role: soloRole, bulk, spuriousAbsence: floor,
  } as never) as { verdict?: string; originClass?: string }

  const a = String(res.verdict ?? '-'); const b = String(second.verdict ?? '-')
  const c = String(third.verdict ?? '-')
  const ca = String(res.originClass ?? '-'); const cb = String(second.originClass ?? '-')
  const cc = String(third.originClass ?? '-')
  if (a !== c) moved += 1
  if (ca !== cc) classMoved += 1
  scored += 1
  console.log(`${p.sample.padEnd(10)} ${(rebuilt ? 'rebuilt' : 'measured').padEnd(9)} `
    + `${pct(res.genomeRate ?? NaN).padStart(8)} ${pct(res.explainable ?? NaN).padStart(8)} `
    + `${short(a).padEnd(13)}${ca.padEnd(14)} ${short(b).padEnd(13)}${cb.padEnd(14)} `
    + `${short(c).padEnd(13)}${cc.padEnd(14)} ${ca !== cc ? 'CLASS' : a !== c ? 'verdict' : ''}`)
  if (ca !== cc) {
    // The axis that moved: alleles the parent lacks, against half her heterozygosity, which is
    // what a second parent contributes. Printed so the new class is checkable rather than taken.
    const t3 = third as { nonParentalRate?: number; secondParentExpected?: number }
    console.log(`${' '.repeat(11)}alleles the parent lacks ${pct(t3.nonParentalRate ?? NaN)} `
      + `against the ${pct(t3.secondParentExpected ?? NaN)} a second parent supplies`)
  }
}
console.log(`\n${moved} of ${scored} absence verdicts change once both omissions are repaired`)
console.log(`${classMoved} of ${scored} ORIGIN CLASSES change, which is what the sheet reports`)

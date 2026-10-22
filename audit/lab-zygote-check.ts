/**
 * WHAT THE SHIPPED PIPELINE SAYS ABOUT THE LAB'S OWN PRONUCLEI, run end to end.
 *
 * The question this answers is narrow and it matters: are the parent calls that were reported on
 * these sheets the TOOL'S output, or an analysis made beside it? Those are different claims and
 * only one of them is a contribution of the tool.
 *
 * Every array is put through `scoreSample` exactly as the browser and the command line do, with the
 * parents named at run time, and the genome-level answer is printed with the stage and zygosity
 * that gate it. Nothing here interprets, classifies or thresholds anything itself.
 *
 * NO SAMPLE IS NAMED IN THIS FILE. Barcodes come from the directory and the pairing file, both
 * supplied at run time and neither tracked. See scripts/release-check.sh.
 *
 * Run: OM_LAB=<dir of .probes files> OM_PAIRS=<tsv: sample, father, mother> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/lab-zygote-check.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)

const DIR = process.env.OM_LAB
const PAIRS = process.env.OM_PAIRS
if (!DIR || !existsSync(DIR) || !PAIRS || !existsSync(PAIRS)) {
  console.log('lab-zygote-check: OM_LAB and OM_PAIRS must both be readable. NOT RUN.')
  process.exit(0)
}

/** Every .probes file under the directory, keyed by the barcode its filename starts with. */
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

function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () { for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r } }
}
/** The head of a parental file, where a reconstruction declares its own parameters. */
function headOf(path: string): string {
  const raw = readFileSync(path)
  return (path.endsWith('.gz') ? gunzipSync(raw) : raw)
    .toString('utf8', 0, inferredArray.MARK_BYTES)
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
  // Reconstructed references carry their heterozygosity and absence floor in the banner; the
  // genotypes cannot supply either. Measured arrays return null here and are unaffected.
  const banner = inferredArray.parseInferredBanner(headOf(fileOf.get(bc)!)) ?? undefined
  const idx = score.finishParent(
    acc as never,
    ingest.finishProfile(bc, byChrom as never, bafSums as never, first).build.build,
    banner as never,
  )
  cache.set(bc, idx); return idx
}

const pairs = readFileSync(PAIRS, 'utf8').trim().split('\n')
  .filter((l) => l && !l.startsWith('#'))
  .map((l) => { const [sample, father, mother] = l.split('\t').map((x) => x.trim()); return { sample, father, mother } })

console.log(`${pairs.length} samples to score; ${fileOf.size} arrays found under the directory\n`)
console.log(`${'sample'.padEnd(11)} ${'father'.padEnd(10)} ${'mother'.padEnd(10)} ${'stage'.padEnd(14)} `
  + `${'zygosity'.padEnd(13)} ${'origin class'.padEnd(16)} reported?`)

let reported = 0; let noParent = 0
for (const p of pairs) {
  if (!fileOf.has(p.sample)) { console.log(`${p.sample.padEnd(11)} no array found`); continue }
  // The first parental slot is whichever parent is present, and soloRole names which that is.
  // With the mother alone the slot holds her and the role is maternal; the pipeline asks the one
  // question that array can answer either way.
  const hasF = !!p.father && fileOf.has(p.father)
  const hasM = !!p.mother && fileOf.has(p.mother)
  if (!hasF && !hasM) {
    console.log(`${p.sample.padEnd(11)} ${(p.father || '-').padEnd(10)} ${(p.mother || '-').padEnd(10)} `
      + `${'-'.padEnd(14)} ${'-'.padEnd(13)} ${'NO PARENT ARRAY'.padEnd(16)} NO`)
    noParent += 1
    continue
  }
  const soloRole = hasF ? 'paternal' : 'maternal'
  const pat = hasF ? parentIndex(p.father) : parentIndex(p.mother)
  const mat = hasF && hasM ? parentIndex(p.mother) : null
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
    stage?: { stage?: string }; zygosity?: string; originClass?: string; origin?: { class?: string }
  }
  const stage = String(res.stage?.stage ?? '-')
  const zyg = String(res.zygosity ?? '-')
  const cls = String(res.originClass ?? res.origin?.class ?? '-')
  const says = cls !== '-' && cls !== 'unknown' && cls !== 'unclear'
  if (says) reported += 1
  console.log(`${p.sample.padEnd(11)} ${(p.father || '-').padEnd(10)} ${(p.mother || '-').padEnd(10)} `
    + `${stage.padEnd(14)} ${zyg.padEnd(13)} ${cls.padEnd(16)} ${says ? 'yes' : 'NO'}`)
}
console.log(`\nthe shipped pipeline names a genome-level class on ${reported} of ${pairs.length}`)
console.log(`${noParent} of those have no parental array at all, so the pipeline cannot be asked`)

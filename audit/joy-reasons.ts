/**
 * What a reader sees when a call is withheld, on real arrays.
 *
 * Two things shipped together. A sample now reports PLOIDY beside the material word, so a diploid
 * cell says diploid rather than only 'trophectoderm'. And a chromosome that gets no contribution
 * call now carries the reason, where before it was simply absent and nothing said why.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/joy-reasons.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const segs = await import(`${W}segments.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('contested-flag-tests: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

const fileOf = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const m = /^(GSM\d+)/.exec(f)
  if (m && (f.endsWith('.gz') || f.endsWith('.probes') || f.endsWith('.txt'))) fileOf.set(m[1], join(DIR, f))
}
function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r }
  }
}
const cache = new Map<string, unknown>()
function parentIndex(gsm: string) {
  const hit = cache.get(gsm); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  cache.set(gsm, idx)
  return idx
}











const stageMod = await import(`${W}stage.ts`)

const manifestPath = join(DIR, 'trio_manifest_full.csv')
if (!existsSync(manifestPath)) { console.log('no manifest. NOT RUN.'); process.exit(0) }
const manifest = readFileSync(manifestPath, 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''

async function scoreSolo(gsm: string, father: string) {
  const pat = parentIndex(father)
  const acc = score.emptyCollected(pat as never, undefined as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, undefined as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat: undefined, soloRole: 'paternal', sibs: [], sampleName: gsm,
    log: () => {},
  } as never) as { chroms: { chrom: string; contribution?: { ploidy: string }
    contributionWithheld?: string }[] }
  return { res, profile }
}

let samples = 0
let called = 0
let withheld = 0
const kinds = new Map<string, number>()
const shown: string[] = []

console.log('PLOIDY, now said beside the material word\n')
console.log(`${'sample'.padEnd(13)} ${'stage'.padEnd(15)} ${'ploidy'.padEnd(15)} `
  + `${'het'.padStart(7)}  chromosomes called / withheld`)
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const gsm = at(r, 'gsm')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  if (!gsm || !father || !fileOf.has(gsm)) continue
  if (samples >= 12) break
  let out
  try { out = await scoreSolo(gsm, father) } catch (e) {
    console.log(`  ${gsm}: ${(e as Error).message}`); continue
  }
  samples += 1
  const st = stageMod.inferStage(out.profile as never) as
    { stage: string; ploidy: string }
  let cN = 0; let wN = 0
  for (const c of out.res.chroms) {
    if (c.contribution) { cN += 1; called += 1 } else if (c.contributionWithheld) {
      wN += 1; withheld += 1
      const key = /informative markers/.test(c.contributionWithheld) ? 'too few markers'
        : /genotyped at/.test(c.contributionWithheld) ? 'call rate too low' : 'other'
      kinds.set(key, (kinds.get(key) ?? 0) + 1)
      if (shown.length < 6) shown.push(`  ${gsm} chr${c.chrom}: ${c.contributionWithheld}`)
    }
  }
  const prof = out.profile as { hetRate: number }
  console.log(`${gsm.padEnd(13)} ${st.stage.padEnd(15)} ${st.ploidy.padEnd(15)} `
    + `${prof.hetRate.toFixed(4).padStart(7)}  ${String(cN).padStart(3)} / ${wN}`)
}

console.log(`\n\nWITHHELD CONTRIBUTIONS, with the reason that used to be absent`)
console.log(`  chromosomes with a call:      ${called}`)
console.log(`  chromosomes withheld:         ${withheld}`)
for (const [k, n] of [...kinds].sort((a, b) => b[1] - a[1])) {
  console.log(`    ${k.padEnd(20)} ${n}`)
}
console.log(`  withheld with NO reason:      `
  + `${withheld === [...kinds.values()].reduce((a, b) => a + b, 0) ? 0 : 'SOME, which is the bug'}`)
console.log('\n  examples as a reader sees them:')
for (const s of shown) console.log(s)

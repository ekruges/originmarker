/**
 * 5.46.0 against real arrays, across every corpus on the box.
 *
 * Four things are checked, and only the first is a pass/fail invariant:
 *
 *   A. NO SILENT WITHHOLD. Every chromosome without a contribution carries a reason. This is the
 *      whole point of the release, so a single exception is a failure.
 *   B. THE REASON IS TRUE. Each stated figure is re-derived from the array rather than trusted,
 *      because a reason that is wrong is worse than none.
 *   C. PLOIDY AGAINST WHAT IS KNOWN. The lab's own sheet text where there is one.
 *   D. NOTHING RECLASSIFIED. 5.46.0's stage against 5.45.6's on the same profiles, which is the
 *      claim the changelog makes.
 *
 * Run: OM_TRIOS=<dir> [OM_CALLS=<manual_calls.tsv>] node --experimental-strip-types \
 *        audit/euploidy-validation.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const stageMod = await import(`${W}stage.ts`)
const obl = await import(`${W}obligateHet.ts`)
const oldStage = await import('./stage-5456.ts')

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('OM_TRIOS not readable. NOT RUN.'); process.exit(0) }
const LIMIT = Number(process.env.OM_LIMIT ?? 400)

const fileOf = new Map<string, string>()
// Some corpora file their arrays in nested folders by experiment date, so the walk has to descend.
// A flat readdir found nothing in two of them and reported a vacuous pass over zero samples.
const walk = (dir: string, depth = 0): void => {
  if (depth > 5) return
  for (const f of readdirSync(dir)) {
    const full = join(dir, f)
    let st
    try { st = statSync(full) } catch { continue }
    if (st.isDirectory()) { walk(full, depth + 1); continue }
    // `._name` is a macOS AppleDouble sidecar that rsync carried in. It is not gzip and not an
    // array, and because it sorts first it was being picked as the reference for every pair.
    if (f.startsWith('._') || /series_matrix/i.test(f)) continue
    if (!/\.(probes|CEL)(\.txt)?(\.gz)?$/i.test(f)) continue
    const id = (/^(GSM\d+)/.exec(f)?.[1]) ?? f.replace(/\.(probes|CEL).*$/i, '')
    if (!id || fileOf.has(id)) continue
    fileOf.set(id, full)
  }
}
walk(DIR)
console.log(`${fileOf.size} array file(s) found under ${DIR}`)
function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header not recognised')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r }
  }
}
const cache = new Map<string, unknown>()
function parentIndex(id: string) {
  const hit = cache.get(id); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(id)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never,
    ingest.finishProfile(id, byChrom as never, bafSums as never, first).build.build)
  cache.set(id, idx)
  return idx
}

// Sheet text, where the lab wrote any.
const sheet = new Map<string, string>()
if (process.env.OM_CALLS && existsSync(process.env.OM_CALLS)) {
  for (const l of readFileSync(process.env.OM_CALLS, 'utf8').split('\n')) {
    const c = l.split('\t')
    if (c[0] && c[3] && !sheet.has(c[0].trim())) sheet.set(c[0].trim(), c[3].trim())
  }
}

type Pair = { gsm: string; father?: string }
const pairs: Pair[] = []
const manifestPath = join(DIR, 'trio_manifest_full.csv')
if (existsSync(manifestPath)) {
  const m = readFileSync(manifestPath, 'utf8').trim().split('\n')
  const cols = m[0].split(',')
  const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
  for (const line of m.slice(1)) {
    const r = line.split(',')
    const g = at(r, 'gsm')
    const f = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
    if (g && fileOf.has(g)) pairs.push({ gsm: g, father: f })
  }
} else {
  // No manifest: use the lowest-numbered array of each series as the reference.
  const ids = [...fileOf.keys()].sort()
  const ref = ids[0]
  for (const g of ids.slice(1)) pairs.push({ gsm: g, father: ref })
}

let silent = 0; let withheld = 0; let called = 0
let reasonWrong = 0
let unmeasuredN = 0
let reasonChecked = 0
const kinds = new Map<string, number>()
const ploidyBy = new Map<string, number>()
const stageMoved: string[] = []
const mismatches: string[] = []
const wrongDetail: string[] = []
let n = 0

for (const p of pairs.slice(0, LIMIT)) {
  if (!p.father) continue
  let res; let profile
  try {
    const pat = parentIndex(p.father)
    const acc = score.emptyCollected(pat as never, undefined as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    // Per-chromosome marker and call counts, so the stated reason can be re-derived.
    const counts = new Map<string, { markers: number; called: number }>()
    for (const r of rowsOf(fileOf.get(p.gsm)!)()) {
      if (!first) first = r.probesetId
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, pat as never, undefined as never, acc as never)
      const c = counts.get(r.chrom) ?? { markers: 0, called: 0 }
      c.markers += 1
      if (r.genotype !== 'NC') c.called += 1
      counts.set(r.chrom, c)
    }
    profile = ingest.finishProfile(p.gsm, byChrom as never, bafSums as never, first)
    res = await score.scoreSample({
      acc, profile, pat, mat: undefined, soloRole: 'paternal', sibs: [],
      sampleName: p.gsm, log: () => {},
    } as never) as { chroms: { chrom: string; contribution?: { informative: number }
      contributionWithheld?: string }[]
      unmeasured: { chrom: string; why: string }[] }
    n += 1

    // --- A and B
    for (const c of res.chroms) {
      if (c.contribution) { called += 1; continue }
      withheld += 1
      if (!c.contributionWithheld) {
        silent += 1
        mismatches.push(`${p.gsm} chr${c.chrom}: withheld with NO reason`)
        continue
      }
      const why = c.contributionWithheld
      const kind = /autosomes only/.test(why) ? 'sex chromosome, channel does not apply'
        : /confined to this chromosome/.test(why) ? 'call rate low, this chromosome only'
        : /not good enough quality/.test(why) ? 'call rate low, the whole array'
        : /No array can do better/.test(why) ? 'marker-poor chromosome'
        : /this array did not read them/.test(why) ? 'array did not read them'
        : /informative markers/.test(why) ? 'too few markers (offer unknown)'
        : /genotyped at/.test(why) ? 'call rate too low' : 'other'
      kinds.set(kind, (kinds.get(kind) ?? 0) + 1)
      // The two new reasons must each name two counts, and the one that blames the parents must
      // name an offer below the floor while the one that blames the array names an offer above it.
      const pair = /only ([\d,]+) sites/.exec(why) ?? /offer ([\d,]+) informative/.exec(why)
      if (pair) {
        const off = Number(pair[1].replace(/,/g, ''))
        const blamesParents = /No array can do better/.test(why)
        if (blamesParents !== (off < 200)) {
          reasonWrong += 1
          wrongDetail.push(`${p.gsm} chr${c.chrom}: offer ${off} but blames `
            + `${blamesParents ? 'the parents' : 'the array'}`)
        }
        reasonChecked += 1
      }
      // B: re-derive the number the reason states.
      const stated = /genotyped at (\d+)% of its markers/.exec(why)
      if (stated) {
        reasonChecked += 1
        const cc = counts.get(c.chrom)
        const actual = cc && cc.markers ? Math.round(100 * cc.called / cc.markers) : -1
        if (Math.abs(actual - Number(stated[1])) > 1) {
          reasonWrong += 1
          wrongDetail.push(`${p.gsm} chr${c.chrom}: says ${stated[1]}%, array says ${actual}%`)
        }
      }
    }

    // A CHROMOSOME WITH NO ROW is a silent refusal too, and the one that mattered most: it is what
    // an array that failed on one chromosome produces. Counted here with the same reasons.
    for (const u of res.unmeasured) {
      unmeasuredN += 1
      const k = /No array can do better/.test(u.why) ? 'not measured: marker-poor'
        : /confined to this chromosome/.test(u.why) ? 'not measured: this chromosome only'
        : /not good enough quality/.test(u.why) ? 'not measured: the whole array'
        : /did not read them/.test(u.why) ? 'not measured: array did not read them'
        : /genotyped at/.test(u.why) ? 'not measured: call rate too low' : 'not measured: other'
      kinds.set(k, (kinds.get(k) ?? 0) + 1)
      if (!u.why) { silent += 1; mismatches.push(`${p.gsm} chr${u.chrom}: unmeasured, NO reason`) }
    }

    // --- C and D
    const now = stageMod.inferStage(profile as never) as { stage: string; ploidy: string }
    const before = oldStage.inferStage(profile as never) as { stage: string }
    ploidyBy.set(now.ploidy, (ploidyBy.get(now.ploidy) ?? 0) + 1)
    if (now.stage !== before.stage) {
      stageMoved.push(`${p.gsm}: ${before.stage} -> ${now.stage}`)
    }
    const txt = sheet.get(p.gsm.replace(/^GSM\d+_/, '')) ?? sheet.get(p.gsm)
    if (txt) {
      const saysHaploid = /\b23\s*[xy]\b|mii|polar body/i.test(txt)
      const saysDiploid = /\b46\b|morula|hesc|te biopsy/i.test(txt)
      if (saysHaploid && now.ploidy !== 'haploid') {
        mismatches.push(`${p.gsm}: sheet "${txt.slice(0, 40)}" but ploidy ${now.ploidy}`)
      }
      if (saysDiploid && now.ploidy !== 'diploid') {
        mismatches.push(`${p.gsm}: sheet "${txt.slice(0, 40)}" but ploidy ${now.ploidy}`)
      }
    }
  } catch (e) {
    console.log(`  ${p.gsm}: ${(e as Error).message}`)
  }
}

console.log(`\n${n} samples scored from ${DIR}`)
if (n === 0) {
  console.log(`NOT RUN: ${fileOf.size} files found, ${pairs.length} pair(s) built. A pass over `
    + 'zero samples is not a result.')
  process.exit(0)
}
console.log()
console.log('A. NO SILENT WITHHOLD')
console.log(`   chromosomes with a call:   ${called}`)
console.log(`   chromosomes withheld:      ${withheld}`)
console.log(`   chromosomes with no row:   ${unmeasuredN}`)
console.log(`   withheld with NO reason:   ${silent}`)
for (const [k, v] of [...kinds].sort((a, b) => b[1] - a[1])) console.log(`     ${k.padEnd(20)} ${v}`)
console.log(silent === 0 ? '   PASS: every withheld chromosome states why.'
  : `   FAIL: ${silent} chromosome(s) withheld in silence.`)

console.log('\nB. THE STATED REASON RE-DERIVED FROM THE ARRAY')
console.log(`   reasons carrying a figure:  ${reasonChecked}`)
console.log(`   figures that do not match:  ${reasonWrong}`)
for (const d of wrongDetail.slice(0, 8)) console.log(`     ${d}`)
console.log(reasonWrong === 0 ? '   PASS: every stated call rate matches the array.'
  : `   FAIL: ${reasonWrong} reason(s) state a figure the array does not support.`)

console.log('\nC. PLOIDY')
for (const [k, v] of [...ploidyBy].sort((a, b) => b[1] - a[1])) console.log(`   ${k.padEnd(18)} ${v}`)
console.log(`   disagreements with the sheet: ${mismatches.length}`)
for (const d of mismatches.slice(0, 10)) console.log(`     ${d}`)

console.log('\nD. NOTHING RECLASSIFIED, 5.46.0 against 5.45.6')
console.log(`   samples whose stage moved: ${stageMoved.length} of ${n}`)
for (const d of stageMoved.slice(0, 10)) console.log(`     ${d}`)
console.log(stageMoved.length === 0 ? '   PASS: every array reaches the stage it reached before.'
  : '   FAIL: the stage axis moved, which this release was not supposed to touch.')

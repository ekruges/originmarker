/**
 * Joy's table: every unclear call, by sample and chromosome, with the reason.
 *
 * "Unclear" is four different things and she gets all four, labelled:
 *   WITHHELD      no contribution call at all, because the markers or the call rate were too thin
 *   NOT MEASURED  the chromosome never reached the measurement list, which is what an array that
 *                 failed on one chromosome produces. These left the table silently before
 *   UNCALLED      a call was attempted and the heterozygosity sat between the two boundaries
 *   MIXED         between the boundaries with one parent dominant: a mosaic, a partial loss, or
 *                 noise
 *
 * Writes a tab-separated file she can open, plus a summary. The reasons are the ones the pipeline
 * itself produces, not text invented for the table, so what she reads here is what the tool said.
 *
 * Run: OM_TRIOS=<dir> OM_OUT=<file.tsv> node --experimental-strip-types audit/joy-unclear-table.ts
 */
import { readFileSync, writeFileSync, existsSync, readdirSync, statSync } from 'node:fs'
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


const OUT = process.env.OM_OUT ?? '/tmp/unclear.tsv'

type Row = { sample: string; chrom: string; kind: string; reason: string
  informative: string; callRate: string; stage: string; ploidy: string }
const rows: Row[] = []
let samples = 0
let clear = 0

for (const p of pairs.slice(0, Number(process.env.OM_LIMIT ?? 400))) {
  if (!p.father) continue
  try {
    const pat = parentIndex(p.father)
    const acc = score.emptyCollected(pat as never, undefined as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
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
    const profile = ingest.finishProfile(p.gsm, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat, mat: undefined, soloRole: 'paternal', sibs: [],
      sampleName: p.gsm, log: () => {},
    } as never) as { chroms: { chrom: string; mixed?: true
      contribution?: { ploidy: string; why: string; informative: number }
      contributionWithheld?: string }[]
      unmeasured: { chrom: string; why: string }[] }
    samples += 1
    const st = stageMod.inferStage(profile as never) as { stage: string; ploidy: string }
    for (const c of res.chroms) {
      const cc = counts.get(c.chrom)
      const cr = cc && cc.markers ? (cc.called / cc.markers) : NaN
      const base = {
        sample: p.gsm, chrom: c.chrom, stage: st.stage, ploidy: st.ploidy,
        callRate: Number.isFinite(cr) ? cr.toFixed(3) : '',
        informative: c.contribution ? String(c.contribution.informative) : '',
      }
      if (c.contributionWithheld) {
        rows.push({ ...base, kind: 'WITHHELD', reason: c.contributionWithheld })
      } else if (c.contribution && c.contribution.ploidy === 'uncalled') {
        rows.push({ ...base, kind: 'UNCALLED', reason: c.contribution.why })
      } else if (c.mixed) {
        rows.push({ ...base, kind: 'MIXED', reason: c.contribution?.why ?? '' })
      } else {
        clear += 1
      }
    }
    // AUTOSOMES WITH NO ROW AT ALL, which until now left the table silently. A chromosome the
    // array failed to genotype is excluded from the measurement list, correctly, and was therefore
    // excluded from this table too: the one result worth reporting was the one that did not appear.
    for (const u of res.unmeasured) {
      const cc = counts.get(u.chrom)
      const cr = cc && cc.markers ? (cc.called / cc.markers) : NaN
      rows.push({
        sample: p.gsm, chrom: u.chrom, stage: st.stage, ploidy: st.ploidy,
        callRate: Number.isFinite(cr) ? cr.toFixed(3) : '',
        informative: '', kind: 'NOT MEASURED', reason: u.why,
      })
    }
  } catch (e) {
    console.log(`  ${p.gsm}: ${(e as Error).message}`)
  }
}

const head = ['sample', 'chromosome', 'unclear_because', 'reason', 'informative_markers',
  'chromosome_call_rate', 'stage', 'ploidy']
const tsv = [head.join('\t'), ...rows.map((r) => [r.sample, `chr${r.chrom}`, r.kind, r.reason,
  r.informative, r.callRate, r.stage, r.ploidy].join('\t'))].join('\n')
writeFileSync(OUT, tsv + '\n')

console.log(`${samples} samples, ${clear} chromosomes called cleanly, ${rows.length} unclear\n`)
const by = new Map<string, number>()
for (const r of rows) by.set(r.kind, (by.get(r.kind) ?? 0) + 1)
for (const [k, v] of [...by].sort((a, b) => b[1] - a[1])) {
  console.log(`  ${k.padEnd(13)} ${String(v).padStart(5)}`)
}
console.log('\nthe distinct reasons, which is the list worth reading:')
const reasons = new Map<string, number>()
for (const r of rows) {
  const gist = r.reason.replace(/[\d.,]+/g, 'N').slice(0, 90)
  reasons.set(gist, (reasons.get(gist) ?? 0) + 1)
}
for (const [k, v] of [...reasons].sort((a, b) => b[1] - a[1]).slice(0, 10)) {
  console.log(`  ${String(v).padStart(5)}  ${k}`)
}
console.log(`\nwritten to ${OUT}`)

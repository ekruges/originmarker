/**
 * The five tests of audit/contested-flag-preregistration.txt, in one pass over the public cells.
 *
 * A: why a cell that is nullisomic almost everywhere produces no segment.
 * B: whether the share annotation reaches the segment, or the 0.1 Mb key drops it.
 * C: reciprocal overlap against bare intersection, both r values per pair.
 * D: the denominator for an absence scan, 22 deletion-bearing events rather than 27.
 * E: the run channel's own accuracy against the hand label, beside the share's.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/contested-flag-tests.ts
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
async function scoreWith(gsm: string, father: string, mother: string) {
  const pat = parentIndex(father); const mat = parentIndex(mother)
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
  } as never) as never
  return { res, profile } as { res: Res; profile: { callRate: number } }
}

type Ann = { where: string; origin: string; called: boolean }
type Seg = {
  chrom: string; startBp: number; endBp: number; kind: string; spanBp: number
  absentParent?: 'paternal' | 'maternal'; contested?: boolean
  refined?: unknown
  alleleGain?: { shift: number; informative: number }
}
type Res = {
  segments?: Seg[]; gains?: Ann[]; losses?: Ann[]
  stage?: { stage: string }
  chroms?: { chrom: string; lostParent?: string; verdict: string }[]
}

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const parentsOf = new Map<string, { father: string; mother: string }>()
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  const mother = at(r, 'mother_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  if (father && mother) parentsOf.set(at(r, 'gsm'), { father, mother })
}
const hand = new Map<string, string>()
for (const l of readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8').split('\n')) {
  if (!l || l.startsWith('#') || l.startsWith('gsm')) continue
  const [gsm, , , chrom, manual] = l.split('\t')
  hand.set(`${gsm}|${chrom}`, manual)
}
const cells = [...new Set([...hand.keys()].map((k) => k.split('|')[0]))]
  .filter((g) => parentsOf.has(g) && fileOf.has(g)).sort()
const parentOf = (h: string) => (/\(\s*pat/i.test(h) ? 'paternal' : /\(\s*mat/i.test(h) ? 'maternal' : '')
const isDup = (h: string) => /dup/i.test(h) && !/del/i.test(h)
const events = [...hand.entries()].filter(([, v]) => /del|dup/i.test(v))
  .map(([k, v]) => ({ gsm: k.split('|')[0], chrom: k.split('|')[1], hand: v }))
  .filter((w) => cells.includes(w.gsm) && parentOf(w.hand))

// The key defects.ts joins the share annotation to the segment on.
const whereOf = (sg: Seg) => {
  const co = segs.segmentCoords(sg as never) as { start: number; end: number }
  return `chr${sg.chrom} ${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb`
}
const coordsOf = (sg: Seg) => segs.segmentCoords(sg as never) as { start: number; end: number }
const rOf = (a: [number, number], b: [number, number]) => {
  const inter = Math.min(a[1], b[1]) - Math.max(a[0], b[0])
  return inter <= 0 ? 0 : inter / (Math.max(a[1], b[1]) - Math.min(a[0], b[0]))
}

const perCell: { gsm: string; stage: string; callRate: number; lost: string[];
  notMeasured: number; segs: Seg[]; anns: Ann[] }[] = []
for (const gsm of cells) {
  const p = parentsOf.get(gsm)!
  const { res, profile } = await scoreWith(gsm, p.father, p.mother)
  perCell.push({
    gsm,
    stage: res.stage?.stage ?? 'none',
    callRate: profile.callRate,
    lost: (res.chroms ?? []).filter((c) => c.lostParent).map((c) => c.chrom),
    notMeasured: (res.chroms ?? []).filter((c) => c.verdict === 'not_measured').length,
    segs: res.segments ?? [],
    anns: [...(res.gains ?? []), ...(res.losses ?? [])],
  })
  process.stdout.write('.')
}
console.log('\n')

// --- A -------------------------------------------------------------------------------------
console.log('TEST A. THE GATES, ON EVERY CELL THE REVIEW AND THE RUN DISAGREE ABOUT')
const review = ['GSM4552418', 'GSM4472409', 'GSM4552422', 'GSM4552378', 'GSM4552380', 'GSM4552382']
const contestedCells = new Set(perCell.filter((c) => c.segs.some((s) => s.contested)).map((c) => c.gsm))
console.log(`${'cell'.padEnd(12)} ${'stage'.padEnd(14)} ${'callRate'.padStart(8)} `
  + `${'segs'.padStart(5)} ${'absence'.padStart(8)} ${'contested'.padStart(9)} `
  + `${'lostByAlleles'.padStart(13)} ${'notMeas'.padStart(7)}`)
for (const c of perCell.filter((x) => review.includes(x.gsm) || contestedCells.has(x.gsm))) {
  const abs = c.segs.filter((s) => s.kind === 'parental-absence')
  console.log(`${c.gsm.padEnd(12)} ${c.stage.padEnd(14)} ${c.callRate.toFixed(2).padStart(8)} `
    + `${String(c.segs.length).padStart(5)} ${String(abs.length).padStart(8)} `
    + `${String(abs.filter((s) => s.contested).length).padStart(9)} `
    + `${String(c.lost.length).padStart(13)} ${String(c.notMeasured).padStart(7)}`)
}

// --- B -------------------------------------------------------------------------------------
console.log('\nTEST B. DOES THE SHARE ANNOTATION REACH THE SEGMENT, OR DOES THE 0.1 Mb KEY DROP IT')
let onKey = 0; let sameChromDifferentKey = 0; let noneOnChrom = 0
for (const c of perCell) {
  const byKey = new Set(c.anns.map((a) => a.where))
  for (const sg of c.segs.filter((s) => s.kind === 'parental-absence')) {
    const k = whereOf(sg)
    if (byKey.has(k)) { onKey += 1; continue }
    if (c.anns.some((a) => a.where.startsWith(`chr${sg.chrom} `))) sameChromDifferentKey += 1
    else noneOnChrom += 1
  }
}
console.log(`  share annotation on the segment's own key:          ${onKey}`)
console.log(`  annotation on the SAME CHROMOSOME, different key:   ${sameChromDifferentKey}`)
console.log(`  no annotation on that chromosome at all:            ${noneOnChrom}`)

// --- C -------------------------------------------------------------------------------------
console.log('\nTEST C. RECIPROCAL OVERLAP, ON WINDOW EDGES AND ON THE EDGES THE TEST USES')
let pairs = 0; let pairsOverThreshold = 0
for (const c of perCell) {
  const abs = c.segs.filter((s) => s.kind === 'parental-absence' && s.absentParent)
  for (let i = 0; i < abs.length; i += 1) {
    for (let j = i + 1; j < abs.length; j += 1) {
      const a = abs[i]; const b = abs[j]
      if (a.chrom !== b.chrom || a.absentParent === b.absentParent) continue
      if (a.startBp > b.endBp || b.startBp > a.endBp) continue
      pairs += 1
      const rw = rOf([a.startBp, a.endBp], [b.startBp, b.endBp])
      const ca = coordsOf(a); const cb = coordsOf(b)
      const rc = rOf([ca.start, ca.end], [cb.start, cb.end])
      if (rc >= segs.RECIPROCAL_OVERLAP) pairsOverThreshold += 1
      console.log(`  ${c.gsm} chr${a.chrom.padStart(2)}  `
        + `${(a.startBp / 1e6).toFixed(1)}-${(a.endBp / 1e6).toFixed(1)} ${a.absentParent}  vs  `
        + `${(b.startBp / 1e6).toFixed(1)}-${(b.endBp / 1e6).toFixed(1)} ${b.absentParent}  `
        + `r(window) ${rw.toFixed(3)}  r(coords) ${rc.toFixed(3)}`
        + `${rc >= segs.RECIPROCAL_OVERLAP ? '' : '   SEPARATE'}`)
    }
  }
}
const allAbs = perCell.flatMap((c) => c.segs.filter((s) => s.kind === 'parental-absence'))
console.log(`  intersecting opposite-parent pairs:   ${pairs}`)
console.log(`  of those, r >= ${segs.RECIPROCAL_OVERLAP} so ONE region:   ${pairsOverThreshold}`)
console.log(`  parental-absence segments:           ${allAbs.length}`)
console.log(`  marked contested under the new rule: ${allAbs.filter((s) => s.contested).length}`)

// --- D and E -------------------------------------------------------------------------------
console.log('\nTEST D AND E. THE DENOMINATOR, AND BOTH CHANNELS AGAINST THE HAND LABEL')
const segOf = new Map<string, Seg[]>()
for (const c of perCell) {
  for (const sg of c.segs) {
    const k = `${c.gsm}|${sg.chrom}`
    segOf.set(k, [...(segOf.get(k) ?? []), sg])
  }
}
const del = events.filter((e) => !isDup(e.hand))
const dup = events.filter((e) => isDup(e.hand))
const found = (e: { gsm: string; chrom: string }) => segOf.has(`${e.gsm}|${e.chrom}`)
const tally = (list: typeof events, role: string) =>
  list.filter((e) => parentOf(e.hand) === role)
console.log(`  hand-called segmental events:                  ${events.length}`)
console.log(`  DELETION-BEARING (an absence scan can see):    ${del.length}`)
console.log(`  duplication-only (it structurally cannot):     ${dup.length}`)
console.log(`  detected, of the deletion-bearing:             ${del.filter(found).length} of ${del.length}`)
console.log(`    paternal:                                    ${tally(del, 'paternal').filter(found).length} of ${tally(del, 'paternal').length}`)
console.log(`    maternal:                                    ${tally(del, 'maternal').filter(found).length} of ${tally(del, 'maternal').length}`)
console.log(`  detected, of the duplication-only:             ${dup.filter(found).length} of ${dup.length}`)
for (const e of dup) {
  const got = segOf.get(`${e.gsm}|${e.chrom}`)
  console.log(`    ${e.gsm} chr${e.chrom.padStart(2)}  ${e.hand.padEnd(26)} `
    + `kinds: ${got ? got.map((s) => s.kind).join(',') : 'none'}`)
}
let runRight = 0; let runWrong = 0; let shareRight = 0; let shareWrong = 0; let runSilent = 0
console.log(`\n  ${'cell'.padEnd(12)} ${'chr'.padStart(4)} ${'hand'.padEnd(10)} ${'run'.padEnd(10)} ${'share'.padEnd(10)}`)
for (const e of del) {
  const got = segOf.get(`${e.gsm}|${e.chrom}`)
  if (!got) continue
  const pr = parentOf(e.hand)
  const abs = got.filter((s) => s.kind === 'parental-absence' && !s.contested)
  const run = abs.map((s) => s.absentParent).find(Boolean)
  const shift = abs.map((s) => s.alleleGain?.shift).find((x) => x !== undefined)
  const share = shift === undefined ? '' : (shift < 0 ? 'paternal' : 'maternal')
  if (!run) runSilent += 1
  else if (run === pr) runRight += 1
  else runWrong += 1
  if (share) { if (share === pr) shareRight += 1; else shareWrong += 1 }
  console.log(`  ${e.gsm.padEnd(12)} ${e.chrom.padStart(4)} ${pr.padEnd(10)} `
    + `${(run ?? '-').padEnd(10)} ${(share || '-').padEnd(10)}`
    + `${run && run !== pr ? '  RUN WRONG' : ''}`
    + `${share && share !== pr ? '  SHARE WRONG' : ''}`)
}
console.log(`\n  the RUN channel:    right ${runRight}, wrong ${runWrong}, no parental-absence segment ${runSilent}`)
console.log(`  the SHARE channel:  right ${shareRight}, wrong ${shareWrong}`)

/**
 * Which parent a SEGMENT's change fell on, on the cells whose karyotype was read by hand.
 *
 * The whole-chromosome scan dilutes an arm-sized event over the full chromosome, so an arm that
 * moves its own region plainly can read near zero there. This runs the shipped pipeline and prints
 * every called segment with what the alleles say inside its window, beside the hand karyotype.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/segment-origin.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('segment-origin: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

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
  return score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
  } as never) as never
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

const want = [...hand.entries()].filter(([, v]) => /del|dup/i.test(v))
  .map(([k, v]) => ({ gsm: k.split('|')[0], chrom: k.split('|')[1], hand: v }))
  .filter((w) => cells.includes(w.gsm))
const parentOf = (h: string) => (/\(\s*pat/i.test(h) ? 'paternal' : /\(\s*mat/i.test(h) ? 'maternal' : '')

let total = 0; let named = 0
const rows: string[] = []
const seen = new Map<string, { kinds: string[]; gain?: string; shift?: number; info?: number }>()
// EVERY parental-absence segment, flat, because the two channels that name its parent are
// independent and the question is whether they agree. `absentParent` is which of the two absence
// runs produced the call, so it is the parent by construction. The allele share inside the window
// is a separate statistic on the same region. Agreement is checkable; disagreement is a finding.
const abs: { gsm: string; chrom: string; startBp: number; endBp: number; span: number;
  absentParent?: string; contested?: boolean; shareParent?: string; shift?: number;
  info?: number }[] = []
for (const gsm of cells) {
  const p = parentsOf.get(gsm)!
  const res = await scoreWith(gsm, p.father, p.mother) as {
    segments?: { chrom: string; startBp: number; endBp: number; kind: string; spanBp: number;
      absentParent?: string; contested?: boolean;
      alleleGain?: { verdict: string; parent?: string; informative: number; shift: number; z: number } }[]
  }
  for (const sg of res.segments ?? []) {
    total += 1
    if (sg.kind === 'parental-absence') {
      abs.push({
        gsm, chrom: sg.chrom, startBp: sg.startBp, endBp: sg.endBp, span: sg.spanBp,
        absentParent: sg.absentParent, contested: sg.contested,
        shareParent: sg.alleleGain?.shift === undefined ? undefined
          : (sg.alleleGain.shift < 0 ? 'paternal' : 'maternal'),
        shift: sg.alleleGain?.shift,
        info: sg.alleleGain?.informative,
      })
    }
    const key = `${gsm}|${sg.chrom}`
    const e = seen.get(key) ?? { kinds: [] }
    e.kinds.push(sg.kind)
    if (sg.alleleGain) {
      named += 1
      e.gain = sg.alleleGain.verdict + (sg.alleleGain.parent ? ' ' + sg.alleleGain.parent : '')
      e.shift = sg.alleleGain.shift
      e.info = sg.alleleGain.informative
    }
    seen.set(key, e)
  }
  process.stdout.write('.')
}
console.log('\n')
console.log('EVERY HAND-CALLED SEGMENTAL EVENT WITH A NAMED PARENT')
console.log(`${'cell'.padEnd(12)} ${'chr'.padStart(4)} ${'hand says'.padEnd(10)} ${'segment called'.padEnd(22)} ${'alleles'.padEnd(22)} share`)
const tally = { pat: [0, 0], mat: [0, 0] }
for (const w of want.sort((a, b) => a.gsm.localeCompare(b.gsm))) {
  const pr = parentOf(w.hand)
  if (!pr) continue
  const e = seen.get(`${w.gsm}|${w.chrom}`)
  const slot = pr === 'paternal' ? tally.pat : tally.mat
  slot[1] += 1
  if (e) slot[0] += 1
  const dir = e?.shift === undefined ? '' : (e.shift < 0 ? 'paternal depleted' : 'maternal depleted')
  const ok = dir ? ((dir.startsWith(pr)) ? ' RIGHT' : ' WRONG') : ''
  console.log(`${w.gsm.padEnd(12)} ${w.chrom.padStart(4)} ${pr.padEnd(10)} `
    + `${(e ? e.kinds.join(',') : 'none').padEnd(22)} ${((e?.gain ?? '-') + ok).padEnd(22)} `
    + `${e?.shift === undefined ? '' : e.shift.toFixed(3)}`)
}
console.log(`\nsegment called where the hand says a PATERNAL segmental event: ${tally.pat[0]} of ${tally.pat[1]}`)
console.log(`segment called where the hand says a MATERNAL segmental event: ${tally.mat[0]} of ${tally.mat[1]}`)
console.log(`\n${total} segments called across ${cells.length} cells; the alleles reach ${named}`)

// ---------------------------------------------------------------------------------------------
console.log('\nPARENTAL-ABSENCE SEGMENTS: WHICH RUN FOUND IT, AND WHAT THE SHARE SAYS')
const withParent = abs.filter((a) => a.absentParent)
const both = abs.filter((a) => a.absentParent && a.shareParent)
const agree = both.filter((a) => a.absentParent === a.shareParent)
console.log(`  parental-absence segments:                 ${abs.length}`)
console.log(`  carrying which run found them:             ${withParent.length}`)
console.log(`    of those, paternal run:                  ${withParent.filter((a) => a.absentParent === 'paternal').length}`)
console.log(`    of those, maternal run:                  ${withParent.filter((a) => a.absentParent === 'maternal').length}`)
console.log(`  with BOTH a run and a window share:        ${both.length}`)
console.log(`    the two channels agree:                  ${agree.length} of ${both.length}`)
for (const a of both.filter((x) => x.absentParent !== x.shareParent)) {
  console.log(`    DISAGREE ${a.gsm} chr${a.chrom} ${(a.span / 1e6).toFixed(1)}Mb  `
    + `run says ${a.absentParent}, share says ${a.shareParent} (shift ${a.shift?.toFixed(3)}, n ${a.info})`)
}

// THE SAME REGION CALLED BY BOTH RUNS is not a biological event. One lost copy takes one parent's
// alleles with it; both parents' alleles depleted over the same window is dropout or a copy loss,
// and the copy-number scan handles a real copy loss on its own channel. Counted rather than
// assumed, because nothing yet says it happens.
let overlaps = 0
for (let i = 0; i < abs.length; i += 1) {
  for (let j = i + 1; j < abs.length; j += 1) {
    const a = abs[i]; const b = abs[j]
    if (a.gsm !== b.gsm || a.chrom !== b.chrom) continue
    if (a.absentParent === b.absentParent) continue
    if (a.startBp > b.endBp || b.startBp > a.endBp) continue
    overlaps += 1
    console.log(`    BOTH RUNS ${a.gsm} chr${a.chrom} `
      + `${(a.startBp / 1e6).toFixed(1)}-${(a.endBp / 1e6).toFixed(1)}Mb ${a.absentParent} and `
      + `${(b.startBp / 1e6).toFixed(1)}-${(b.endBp / 1e6).toFixed(1)}Mb ${b.absentParent}`)
  }
}
console.log(`  one window called by BOTH runs:            ${overlaps}`)

// WHAT THE WITHHOLDING COSTS. A flag that suppresses a parent is only honest if the number of
// correct calls it suppresses is known. Counted against the hand karyotype, not asserted.
const flagged = abs.filter((a) => a.contested)
const flaggedCells = new Set(flagged.map((a) => a.gsm))
const handCells = new Set(want.filter((w) => parentOf(w.hand)).map((w) => w.gsm))
const detected = new Set(want.filter((w) => parentOf(w.hand) && seen.has(`${w.gsm}|${w.chrom}`))
  .map((w) => `${w.gsm}|${w.chrom}`))
const lost = [...detected].filter((k) => flagged.some((a) => `${a.gsm}|${a.chrom}` === k))
console.log(`  segments marked contested:                 ${flagged.length} of ${abs.length}`)
console.log(`    on how many cells:                       ${flaggedCells.size} of ${cells.length}`)
console.log(`    of those cells, carrying a hand-called segmental event: `
  + `${[...flaggedCells].filter((g) => handCells.has(g)).length}`)
console.log(`    hand-called events that LOSE their parent to the flag:  ${lost.length}`)
for (const k of lost) console.log(`      ${k}`)

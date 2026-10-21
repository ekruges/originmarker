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
for (const gsm of cells) {
  const p = parentsOf.get(gsm)!
  const res = await scoreWith(gsm, p.father, p.mother) as {
    segments?: { chrom: string; startBp: number; endBp: number; kind: string; spanBp: number;
      alleleGain?: { verdict: string; parent?: string; informative: number; shift: number; z: number } }[]
  }
  for (const sg of res.segments ?? []) {
    total += 1
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

/**
 * Whether each pair of a run's arrays reads as ONE genome, by both candidate measures.
 *
 * Item 3.1 needs a same-genome rule that works with NO parent loaded, because the cells it is for
 * have no parent on disk. Two measures can do it and this prints both for every pair, so the choice
 * rests on their separation rather than on taste:
 *
 *   agree  share of shared called markers that are identical (parentage.agreement). This is the
 *          page's existing definition of a unit of one embryo (abnormalities.SAME_EMBRYO_AGREEMENT).
 *   opp    share of shared HOMOZYGOUS calls that are opposite (relatedness.oppositeHom). The rule as
 *          first written in the plan, at 0.02.
 *
 * Each array is read once and kept as one byte per marker, which is why both measures cost one pass.
 * The byte loop is checked against parentage.agreement itself on a sample of pairs, so the two cannot
 * drift: there is one definition and this is an index over it.
 *
 * Run: node --experimental-strip-types audit/truth/sister-pairs.ts <list.tsv> <out.tsv> [sample]
 * list.tsv is id<TAB>path. Every id must be listed open in units.tsv by the caller.
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const { agreement } = await import(`${W}parentage.ts`)

const [listFile, out, sampleArg] = process.argv.slice(2)
const SAMPLE = Number(sampleArg ?? 3)
const text = (p: string) => {
  const raw = readFileSync(p)
  return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
}
/** 0 no call, 1 AA, 2 AB, 3 BB. Homozygous is 1 or 3, and opposite is 1 against 3. */
const CODES = ['', 'AA', 'AB', 'BB']
const index = new Map<string, number>()
let sealed = false
let added = 0

function load(path: string): Uint8Array {
  const lines = text(path).split('\n')
  let h = 0
  while (h < lines.length && (!lines[h] || lines[h].startsWith('#'))) h += 1
  const map = ingest.headerMap(lines[h])
  const pairs: [number, number][] = []
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || !/^\d+$/.test(r.chrom)) continue
    let k = index.get(r.probesetId)
    if (k === undefined) {
      if (sealed) { added += 1; continue }
      index.set(r.probesetId, k = index.size)
    }
    pairs.push([k, CODES.indexOf(r.genotype)])
  }
  sealed = true
  const code = new Uint8Array(index.size)
  for (const [k, g] of pairs) if (g > 0) code[k] = g
  return code
}

const list = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t'))
const calls = new Map<string, Uint8Array>()
for (const [id, path] of list) calls.set(id, load(path))
if (added) console.error(`WARNING ${added} probesets after the first array were dropped: not on it`)

/** Both measures in one pass over the two byte arrays. */
function measure(a: Uint8Array, b: Uint8Array) {
  let shared = 0; let same = 0; let hom = 0; let opp = 0
  const n = Math.min(a.length, b.length)
  for (let i = 0; i < n; i += 1) {
    const x = a[i]; const y = b[i]
    if (x === 0 || y === 0) continue
    shared += 1
    if (x === y) same += 1
    if (x !== 2 && y !== 2) { hom += 1; if (x !== y) opp += 1 }
  }
  return { shared, same, hom, opp, agree: shared > 1000 ? same / shared : NaN, oppRate: hom ? opp / hom : NaN }
}

const ids = [...calls.keys()]
const rows = ['a\tb\tshared\tagree\thom\topp_rate']
const sampled: [string, string][] = []
for (let i = 0; i < ids.length; i += 1) {
  for (let j = i + 1; j < ids.length; j += 1) {
    const m = measure(calls.get(ids[i])!, calls.get(ids[j])!)
    rows.push(`${ids[i]}\t${ids[j]}\t${m.shared}\t${m.agree.toFixed(5)}\t${m.hom}\t${m.oppRate.toFixed(5)}`)
    if (sampled.length < SAMPLE) sampled.push([ids[i], ids[j]])
  }
}
writeFileSync(out, `${rows.join('\n')}\n`)

// THE BYTE LOOP AGAINST THE SHIPPED FUNCTION. A disagreement here means this file is measuring
// something the tool does not, which would make every number above unattributable.
const gtOf = (path: string): Map<string, string> => {
  const lines = text(path).split('\n')
  let h = 0
  while (h < lines.length && (!lines[h] || lines[h].startsWith('#'))) h += 1
  const map = ingest.headerMap(lines[h])
  const g = new Map<string, string>()
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (r && /^\d+$/.test(r.chrom) && CODES.includes(r.genotype)) g.set(r.probesetId, r.genotype)
  }
  return g
}
const pathOf = new Map(list.map(([id, p]) => [id, p]))
for (const [a, b] of sampled) {
  const mine = measure(calls.get(a)!, calls.get(b)!).agree
  const theirs = agreement(gtOf(pathOf.get(a)!) as never, gtOf(pathOf.get(b)!) as never)
  if (!(Math.abs(mine - theirs) < 1e-12)) {
    throw new Error(`${a} ${b}: this file reads ${mine}, parentage.agreement reads ${theirs}`)
  }
}
console.log(`${ids.length} arrays, ${rows.length - 1} pairs, ${sampled.length} checked against parentage.agreement`)

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
 * COUNTS PER CHROMOSOME, NOT PER GENOME, because the aggregate has to be taken over the chromosomes
 * NEITHER cell reports a change on. Two cells of one embryo that lost opposite parental copies of a
 * region read as opposite homozygotes across that whole region (relatedness.ts, four of 122 same-embryo
 * blastomere pairs at 0.053 to 0.056), which is FIRST order and is exactly the reciprocal event this
 * rule is hunting. Aggregating over the whole genome therefore pushes a true sister pair above any
 * gate for carrying the thing being looked for, while a pair of ordinary siblings stays under it. So
 * this file emits one row per pair per chromosome and the aggregate is taken offline, where the scored
 * calls say which chromosomes to leave out.
 *
 * Each array is read once and kept as one byte per marker. The byte loop is checked against
 * parentage.agreement on a sample of pairs, so the two cannot drift: there is one definition of the
 * measure and this is an index over it.
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
/** The chromosome of each index, in the first array's order. */
const chromOf: string[] = []
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
      chromOf[k] = r.chrom
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

/** Index lists per chromosome, so a pair is counted chromosome by chromosome. */
const byChrom = new Map<string, number[]>()
for (let k = 0; k < chromOf.length; k += 1) {
  const c = chromOf[k]
  if (!c) continue
  const g = byChrom.get(c)
  if (g) g.push(k)
  else byChrom.set(c, [k])
}

/**
 * The four joint states of a marker both cells call, over one chromosome.
 *
 *   hh  both heterozygous
 *   hx  one heterozygous, one homozygous. Under ONE genome this is a single dropout, FIRST order.
 *   ss  both homozygous and the same
 *   oo  both homozygous and opposite. Under ONE genome and no copy change this needs a dropout in
 *       both cells at the same heterozygous marker in opposite directions, SECOND order.
 */
function measureOn(a: Uint8Array, b: Uint8Array, ks: readonly number[]) {
  let shared = 0; let same = 0; let hom = 0; let opp = 0; let hh = 0; let hx = 0
  for (const k of ks) {
    const x = a[k]; const y = b[k]
    if (x === undefined || y === undefined || x === 0 || y === 0) continue
    shared += 1
    if (x === y) same += 1
    const hetX = x === 2; const hetY = y === 2
    if (hetX && hetY) hh += 1
    else if (hetX !== hetY) hx += 1
    else { hom += 1; if (x !== y) opp += 1 }
  }
  return { shared, same, hom, opp, hh, hx }
}

const ids = [...calls.keys()]
const chroms = [...byChrom.keys()].filter((c) => /^\d+$/.test(c)).sort((x, y) => Number(x) - Number(y))
const rows = ['a\tb\tchrom\tshared\tsame\thom\topp\thh\thx']
const sampled: [string, string][] = []
for (let i = 0; i < ids.length; i += 1) {
  for (let j = i + 1; j < ids.length; j += 1) {
    const a = calls.get(ids[i])!; const b = calls.get(ids[j])!
    for (const c of chroms) {
      const m = measureOn(a, b, byChrom.get(c)!)
      if (!m.shared) continue
      rows.push(`${ids[i]}\t${ids[j]}\t${c}\t${m.shared}\t${m.same}\t${m.hom}\t${m.opp}\t${m.hh}\t${m.hx}`)
    }
    if (sampled.length < SAMPLE) sampled.push([ids[i], ids[j]])
  }
}
writeFileSync(out, `${rows.join('\n')}\n`)

// THE BYTE LOOP AGAINST THE SHIPPED FUNCTION. A disagreement here means this file is measuring
// something the tool does not, which would make every number above unattributable. The shipped
// function is genome-wide over every chromosome it is given, so the check sums this file's rows the
// same way, over all chromosomes rather than the autosomes alone.
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
const every = [...byChrom.keys()].filter((c) => /^\d+$/.test(c))
for (const [a, b] of sampled) {
  let shared = 0; let same = 0
  for (const c of every) {
    const m = measureOn(calls.get(a)!, calls.get(b)!, byChrom.get(c)!)
    shared += m.shared; same += m.same
  }
  const mine = shared > 1000 ? same / shared : NaN
  const theirs = agreement(gtOf(pathOf.get(a)!) as never, gtOf(pathOf.get(b)!) as never)
  if (!(Math.abs(mine - theirs) < 1e-12)) {
    throw new Error(`${a} ${b}: this file reads ${mine}, parentage.agreement reads ${theirs}`)
  }
}
console.log(`${ids.length} arrays, ${rows.length - 1} pair-chromosome rows, ${sampled.length} checked against parentage.agreement`)

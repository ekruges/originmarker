/**
 * Which arrays carry the genome of which genotyped adult.
 *
 * For each array of a list, against each reference array: over autosomal markers both call
 * homozygous, the share called opposite. A child of the reference, or one of its haploid products,
 * sits near zero plus its own dropout; an unrelated genome sits at the population rate. One line
 * per array, with its call rate and heterozygosity, so bulk material not yet used as a reference
 * can be found in the output.
 *
 * Run: node --experimental-strip-types audit/truth/kin.ts <refs.tsv> <list.tsv> <shard> <shards> <out.jsonl>
 * refs.tsv: name<TAB>path. list.tsv: id<TAB>path.
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { ingest, lineReader } from '../joy-common.ts'

const [refFile, listFile, shard, shards, out] = [process.argv[2], process.argv[3], Number(process.argv[4]), Number(process.argv[5]), process.argv[6]]
const open = (p: string) => { const raw = readFileSync(p); return lineReader(p.endsWith('.gz') ? gunzipSync(raw) : raw) }
function* rows(path: string) {
  const next = open(path); let header: any = null
  for (let line = next(); line !== null; line = next()) {
    if (!header) { if (line && !line.startsWith('#')) header = ingest.headerMap(line); continue }
    const r = ingest.parseRow(line, header); if (r) yield r
  }
}
// Each reference as one byte per marker in the order the first reference lists them: 1 AA, 2 BB, 0 otherwise.
const index = new Map<string, number>()
const refs: { name: string; hom: Uint8Array }[] = []
for (const l of readFileSync(refFile, 'utf8').split('\n').filter(Boolean)) {
  const [name, path] = l.split('\t')
  const hom = new Uint8Array(1_000_000)
  for (const r of rows(path)) {
    if (!/^\d+$/.test(r.chrom)) continue
    let k = index.get(r.probesetId)
    if (k === undefined) { k = index.size; index.set(r.probesetId, k) }
    hom[k] = r.genotype === 'AA' ? 1 : r.genotype === 'BB' ? 2 : 0
  }
  refs.push({ name, hom })
}
const done = new Set<string>()
if (existsSync(out)) for (const l of readFileSync(out, 'utf8').split('\n')) if (l.startsWith('{')) done.add(JSON.parse(l).id)
const list = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t')).filter((_, i) => i % shards === shard)
for (const [id, path] of list) {
  if (done.has(id)) continue
  try {
    const both = new Array(refs.length).fill(0); const opp = new Array(refs.length).fill(0)
    let n = 0; let called = 0; let het = 0
    for (const r of rows(path)) {
      if (!/^\d+$/.test(r.chrom)) continue
      n++
      if (r.genotype === 'NC') continue
      called++
      if (r.genotype === 'AB') { het++; continue }
      const k = index.get(r.probesetId); if (k === undefined) continue
      const g = r.genotype === 'AA' ? 1 : 2
      for (let j = 0; j < refs.length; j++) {
        const h = refs[j].hom[k]
        if (h) { both[j]++; if (h !== g) opp[j]++ }
      }
    }
    appendFileSync(out, JSON.stringify({ id, callRate: +(called / n).toFixed(4), hetRate: +(het / Math.max(1, called)).toFixed(4),
      refs: Object.fromEntries(refs.map((r, j) => [r.name, [both[j], opp[j]]])) }) + '\n')
  } catch (e) { appendFileSync(out, JSON.stringify({ id, error: String(e).slice(0, 200) }) + '\n') }
}
appendFileSync(out, `DONE ${shard}\n`)

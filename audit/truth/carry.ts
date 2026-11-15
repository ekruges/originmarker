/**
 * Does an array carry an adult's genome, chromosome by chromosome?
 *
 * Read only where the adult is homozygous for an allele the founder panel says is RARE (frequency at
 * most RARE). There a genome that inherited from the adult must hold that allele, and a genome that
 * did not almost never does. So, among the array's called genotypes at those markers, the share
 * called homozygous for the OTHER allele is near zero for one of the adult's haploid products, held
 * down to the dropout rate for a child, and near one for anyone else, whatever their ancestry.
 * The same tally at markers where the adult is homozygous for the COMMON allele is the control: a
 * stranger agrees there by chance.
 *
 * Run: node --experimental-strip-types audit/truth/carry.ts <refs.tsv> <list.tsv> <shard> <shards> <out.jsonl>
 * Each line: id, call and het rate, and per reference per autosome
 *   [rare: same-hom, het, opposite-hom, no-call,  common: same-hom, het, opposite-hom, no-call]
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { ingest, lineReader } from '../joy-common.ts'

const R = new URL('../../', import.meta.url).pathname
const keys = await import(R + 'web/src/answerKeys.ts')
const panel = keys.parsePanel(gunzipSync(readFileSync(R + 'web/public/' + keys.FOUNDER_PANEL_PATH)).toString('utf8'))
export const RARE = 0.2
// The lowest founder frequency of the adult's allele that is tallied: 0 keeps alleles no founder holds.
const RARE_MIN = Number(process.env.OM_RARE_MIN ?? 0)
const [refFile, listFile, shard, shards, out] = [process.argv[2], process.argv[3], Number(process.argv[4]), Number(process.argv[5]), process.argv[6]]
function* rows(path: string) {
  const raw = readFileSync(path); const next = lineReader(path.endsWith('.gz') ? gunzipSync(raw) : raw); let header: any = null
  for (let line = next(); line !== null; line = next()) {
    if (!header) { if (line && !line.startsWith('#')) header = ingest.headerMap(line); continue }
    const r = ingest.parseRow(line, header); if (r) yield r
  }
}
// Per reference, one byte per marker: 0 not usable, 1 the adult is AA and A is rare, 2 BB and B rare,
// 3 AA and A common, 4 BB and B common.
const index = new Map<string, number>()
const refs: { name: string; cls: Uint8Array }[] = []
for (const l of readFileSync(refFile, 'utf8').split('\n').filter(Boolean)) {
  const [name, path] = l.split('\t')
  const cls = new Uint8Array(1_000_000)
  for (const r of rows(path)) {
    if (!/^\d+$/.test(r.chrom) || (r.genotype !== 'AA' && r.genotype !== 'BB')) continue
    const fB = panel.freqB?.get(r.probesetId); if (fB === undefined) continue
    const fOwn = r.genotype === 'BB' ? fB : 1 - fB
    const c = fOwn <= RARE && fOwn >= RARE_MIN ? (r.genotype === 'AA' ? 1 : 2)
      : fOwn >= 1 - RARE && fOwn <= 1 - RARE_MIN ? (r.genotype === 'AA' ? 3 : 4) : 0
    if (!c) continue
    let k = index.get(r.probesetId)
    if (k === undefined) { k = index.size; index.set(r.probesetId, k) }
    cls[k] = c
  }
  refs.push({ name, cls })
}
const done = new Set<string>()
if (existsSync(out)) for (const l of readFileSync(out, 'utf8').split('\n')) if (l.startsWith('{')) done.add(JSON.parse(l).id)
const list = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t')).filter((_, i) => i % shards === shard)
for (const [id, path] of list) {
  if (done.has(id)) continue
  try {
    const t = refs.map(() => Array.from({ length: 23 }, () => new Array(8).fill(0)))
    let n = 0; let called = 0; let het = 0
    for (const r of rows(path)) {
      if (!/^\d+$/.test(r.chrom)) continue
      n++
      if (r.genotype !== 'NC') called++
      if (r.genotype === 'AB') het++
      const k = index.get(r.probesetId); if (k === undefined) continue
      const ch = Number(r.chrom)
      for (let j = 0; j < refs.length; j++) {
        const c = refs[j].cls[k]; if (!c) continue
        const adult = c === 1 || c === 3 ? 'AA' : 'BB'
        const base = c <= 2 ? 0 : 4
        const slot = r.genotype === 'NC' ? 3 : r.genotype === 'AB' ? 1 : r.genotype === adult ? 0 : 2
        t[j][ch][base + slot]++
      }
    }
    appendFileSync(out, JSON.stringify({ id, callRate: +(called / n).toFixed(4), hetRate: +(het / Math.max(1, called)).toFixed(4),
      refs: Object.fromEntries(refs.map((r, j) => [r.name, t[j].slice(1)])) }) + '\n')
  } catch (e) { appendFileSync(out, JSON.stringify({ id, error: String(e).slice(0, 200) }) + '\n') }
}
appendFileSync(out, `DONE ${shard}\n`)

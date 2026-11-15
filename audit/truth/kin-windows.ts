/**
 * Along the genome, how often an array's homozygous calls oppose one reference's, per 5 Mb window.
 *
 * A child of the reference is low everywhere. A grandchild, niece or nephew is low only where it
 * inherited the reference's haplotype, in long stretches that every cell of one embryo shares. An
 * unrelated genome is at the population rate everywhere. The whole-genome rate cannot tell the
 * second from a noisy child; the stretches can.
 *
 * Run: node --experimental-strip-types audit/truth/kin-windows.ts <ref path> <list.tsv> <out.jsonl>
 */
import { readFileSync, appendFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { ingest, lineReader } from '../joy-common.ts'

const [refPath, listFile, out] = process.argv.slice(2)
function* rows(path: string) {
  const raw = readFileSync(path); const next = lineReader(path.endsWith('.gz') ? gunzipSync(raw) : raw); let header: any = null
  for (let line = next(); line !== null; line = next()) {
    if (!header) { if (line && !line.startsWith('#')) header = ingest.headerMap(line); continue }
    const r = ingest.parseRow(line, header); if (r) yield r
  }
}
const ref = new Map<string, number>()
for (const r of rows(refPath)) if (/^\d+$/.test(r.chrom) && (r.genotype === 'AA' || r.genotype === 'BB')) ref.set(r.probesetId, r.genotype === 'AA' ? 1 : 2)
for (const [id, path] of readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t'))) {
  const w = new Map<string, [number, number]>()
  for (const r of rows(path)) {
    if (r.genotype !== 'AA' && r.genotype !== 'BB') continue
    const h = ref.get(r.probesetId); if (!h) continue
    const k = `${r.chrom}:${Math.floor(r.pos / 5e6)}`
    let c = w.get(k); if (!c) w.set(k, c = [0, 0])
    c[0]++; if (h !== (r.genotype === 'AA' ? 1 : 2)) c[1]++
  }
  appendFileSync(out, JSON.stringify({ id, windows: [...w].map(([k, [n, o]]) => [k, n, o]) }) + '\n')
}
appendFileSync(out, 'DONE\n')

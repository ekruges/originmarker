/**
 * Panel-monomorphic markers: where every founder in a panel of unrelated people is called the same
 * homozygote. A sample's heterozygous calls there are almost all drop-in, with no parent needed.
 *
 * WHY. Drop-in measured where both PARENTS share a homozygote needs both parents arrayed, and the
 * laboratory's material has at most one, usually reconstructed from the children. A panel of
 * unrelated founders gives the same kind of answer key from the population instead of the family:
 * where everyone in the panel carries one allele on both chromosomes, a sample is very unlikely to
 * carry the other. The residual is real heterozygosity at variants the panel happens to lack, and it
 * shrinks with every founder added; it is measured, not assumed, by scoring a held-out bulk array.
 *
 * Founders must be bulk DNA of distinct, unrelated people. Replicates of one person add nothing and
 * a relative of the sample would bias the set toward the sample's own alleles.
 *
 * Run: node --experimental-strip-types audit/founder-panel.ts <out.tsv> <founder>...
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { lineReader } from './joy-common.ts'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)

const [out, ...founders] = process.argv.slice(2)
if (!out || founders.length < 2) {
  console.log('usage: founder-panel.ts <out.tsv> <founder>... (at least two)')
  process.exit(1)
}

/** probeset -> [chrom, 'A' or 'B'] for this founder's homozygous autosomal calls. */
function homozygotes(path: string): { het: number; called: number; hom: Map<string, [string, 'A' | 'B']> } {
  const raw = readFileSync(path)
  const next = lineReader(path.endsWith('.gz') ? gunzipSync(raw) : raw)
  const head: string[] = []
  for (let i = 0; i < 60; i += 1) { const l = next(); if (l === null) break; head.push(l) }
  const h = head.findIndex((l) => l && !l.startsWith('#'))
  const map = ingest.headerMap(head[h])
  if (!map) throw new Error(`${path}: header not recognised`)
  const hom = new Map<string, [string, 'A' | 'B']>()
  let het = 0; let called = 0
  const take = (l: string) => {
    const r = ingest.parseRow(l, map)
    if (!r || !/^\d+$/.test(r.chrom) || r.genotype === 'NC') return
    called += 1
    if (r.genotype === 'AB') { het += 1; return }
    hom.set(r.probesetId, [r.chrom, r.genotype === 'AA' ? 'A' : 'B'])
  }
  for (const l of head.slice(h + 1)) take(l)
  for (let l = next(); l !== null; l = next()) take(l)
  return { het, called, hom }
}

let panel: Map<string, [string, 'A' | 'B']> | null = null
for (const f of founders) {
  const { het, called, hom } = homozygotes(f)
  console.log(`${f.split('/').pop()}: ${called.toLocaleString()} autosomal calls, `
    + `heterozygous ${(het / called).toFixed(4)}`)
  if (!panel) { panel = hom; continue }
  for (const [id, [, allele]] of panel) {
    const mine = hom.get(id)
    if (!mine || mine[1] !== allele) panel.delete(id)
  }
}
const lines = ['probeset\tchrom\tallele']
for (const [id, [chrom, allele]] of panel!) lines.push(`${id}\t${chrom}\t${allele}`)
writeFileSync(out, lines.join('\n') + '\n')
console.log(`${founders.length} founders, ${panel!.size.toLocaleString()} panel-monomorphic markers -> ${out}`)

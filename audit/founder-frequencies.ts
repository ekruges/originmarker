/**
 * Founder allele frequencies from a public genotype matrix: the parent-free reference the one-parent
 * allele state reads against, and the panel the answer keys read drop-in against.
 *
 * WHY THIS MATERIAL. GSE186441 is constitutional blood DNA of 187 adults arrayed on the same Axiom UK
 * Biobank platform as the laboratory's material, published with every genotype call. It is unrelated to
 * GSE148488, so every validation against GSE148488's families is independent of the panel. A frequency
 * per marker serves two uses: where it is 0 or 1 the marker is panel-monomorphic, the drop-in answer
 * key; and how common the allele a parent lacks is, is what separates a missing copy from a present one
 * at markers where that parent is homozygous.
 *
 * QUALITY. Arrays are kept at a call rate of 0.95 or more and a heterozygosity a single diploid reads,
 * and one of each exact duplicate pair is dropped: a person counted twice halves the evidence that a
 * marker is monomorphic and biases every frequency toward their genotype.
 *
 * Run: node --experimental-strip-types audit/founder-frequencies.ts <matrix.txt.gz> <reference.probes> <out.tsv>
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { lineReader } from './joy-common.ts'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)

const [matrix, reference, out] = process.argv.slice(2)
if (!matrix || !reference || !out) { console.log('usage: founder-frequencies.ts <matrix.gz> <reference.probes> <out.tsv>'); process.exit(1) }

/** probeset -> autosome, from an array on the same platform. */
const chromOf = new Map<string, string>()
{
  const raw = readFileSync(reference)
  const next = lineReader(reference.endsWith('.gz') ? gunzipSync(raw) : raw)
  const head: string[] = []
  for (let i = 0; i < 60; i += 1) { const l = next(); if (l === null) break; head.push(l) }
  const h = head.findIndex((l) => l && !l.startsWith('#'))
  const map = ingest.headerMap(head[h])
  const take = (l: string) => { const r = ingest.parseRow(l, map); if (r && /^\d+$/.test(r.chrom)) chromOf.set(r.probesetId, r.chrom) }
  for (const l of head.slice(h + 1)) take(l)
  for (let l = next(); l !== null; l = next()) take(l)
}
console.log(`${chromOf.size.toLocaleString()} autosomal probesets on the platform`)

const code = (g: string): number => (g === 'AA' ? 0 : g === 'AB' ? 1 : g === 'BB' ? 2 : -1)

// PASS 1: per array, call rate and heterozygosity, and a sparse sample of genotypes for duplicates.
let names: string[] = []
let calls: number[] = []; let hets: number[] = []; let rows = 0
const SAMPLE_EVERY = 40
const sampled: Int8Array[] = []
{
  const next = lineReader(gunzipSync(readFileSync(matrix)))
  for (let l = next(); l !== null; l = next()) {
    if (!l || l.startsWith('#')) continue
    const c = l.split('\t')
    if (c[0] === 'probeset_id') {
      names = c.slice(1); calls = names.map(() => 0); hets = names.map(() => 0)
      continue
    }
    if (!chromOf.has(c[0])) continue
    rows += 1
    const keep = rows % SAMPLE_EVERY === 0
    const s = keep ? new Int8Array(names.length) : null
    for (let i = 0; i < names.length; i += 1) {
      const g = code(c[i + 1])
      if (s) s[i] = g
      if (g < 0) continue
      calls[i] += 1
      if (g === 1) hets[i] += 1
    }
    if (s) sampled.push(s)
  }
}
const callRate = calls.map((n) => n / rows)
const het = calls.map((n, i) => (n ? hets[i] / n : NaN))
let kept = names.map((_, i) => callRate[i] >= 0.95 && het[i] >= 0.12 && het[i] <= 0.20)
console.log(`${names.length} arrays, ${rows.toLocaleString()} autosomal markers; quality: `
  + `${kept.filter(Boolean).length} kept (call rate 0.95 or more, heterozygosity 0.12 to 0.20)`)

// Exact duplicates: opposite homozygotes at sampled markers where both are homozygous.
const dropped: string[] = []
for (let i = 0; i < names.length; i += 1) {
  if (!kept[i]) continue
  for (let j = i + 1; j < names.length; j += 1) {
    if (!kept[j]) continue
    let both = 0; let opp = 0
    for (const s of sampled) {
      const a = s[i]; const b = s[j]
      if ((a === 0 || a === 2) && (b === 0 || b === 2)) { both += 1; if (a !== b) opp += 1 }
    }
    if (both > 1000 && opp / both < 0.005) {
      kept[j] = false
      dropped.push(`${names[j]} (duplicate of ${names[i]}, opposite homozygotes ${(opp / both).toFixed(4)})`)
    }
  }
}
for (const d of dropped) console.log(`  dropped ${d}`)
const founders = kept.filter(Boolean).length
console.log(`${founders} founders`)

// PASS 2: allele counts per marker over the founders.
const lines = [
  `# OriginMarker founder frequencies: GSE186441 (public, GEO), ${founders} unrelated adults, constitutional`,
  '# blood DNA on Axiom UK Biobank, after quality and duplicate filtering. freqB is the B allele frequency',
  '# among founders called at that marker; n is how many were called. Built by audit/founder-frequencies.ts.',
  'probeset\tfreqB\tn',
]
{
  const next = lineReader(gunzipSync(readFileSync(matrix)))
  for (let l = next(); l !== null; l = next()) {
    if (!l || l.startsWith('#') || l.startsWith('probeset_id')) continue
    const c = l.split('\t')
    if (!chromOf.has(c[0])) continue
    let n = 0; let b = 0
    for (let i = 0; i < names.length; i += 1) {
      if (!kept[i]) continue
      const g = code(c[i + 1])
      if (g < 0) continue
      n += 1; b += g
    }
    if (n) lines.push(`${c[0]}\t${(b / (2 * n)).toFixed(3)}\t${n}`)
  }
}
writeFileSync(out, lines.join('\n') + '\n')
console.log(`${(lines.length - 4).toLocaleString()} markers written to ${out}`)

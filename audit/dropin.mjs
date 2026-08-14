// Drop-in, measured on same-individual and clonal replicate arrays.
//
// Drop-in is a FALSE HETEROZYGOUS call at a truly homozygous marker. It is the one term that
// biases a sibling-referenced parental call in a single direction, because a falsely het marker
// always retains the reference allele and so always reads as loss of the other copy. Every
// threshold downstream is a function of it, so it is measured here rather than assumed.
//
// The material is ideal and was produced as a by-product: the linkage runs set aside 217 arrays
// that are close relatives but NOT children, i.e. the same individual, clonal cells, or haploid
// products. Where two arrays are the same genome, a marker the reference calls homozygous is
// homozygous in truth, and any AB call in the replicate is drop-in.
//
// usage: node audit/dropin.mjs <linkage json> <data dir>
import { readFileSync, readdirSync, statSync } from 'node:fs'

const [JSONF, DIR] = process.argv.slice(2)
const walk = (d) => readdirSync(d, { withFileTypes: true })
  .flatMap((e) => (e.isDirectory() ? walk(`${d}/${e.name}`) : [`${d}/${e.name}`]))
const files = walk(DIR).filter((f) => f.endsWith('.probes') && statSync(f).size > 1e6)
const nameOf = (p) => p.split('/').pop().replace(/\.CEL\.probes$/i, '')

const read = (p) => {
  const L = readFileSync(p, 'utf8').split('\n')
  const h = L[0].split('\t')
  const gi = h.indexOf('genotype'); const pi = h.indexOf('probeset_id')
  const ci = h.indexOf('chr')
  const m = new Map()
  for (let i = 1; i < L.length; i += 1) {
    const q = L[i].split('\t')
    if (q.length <= gi) continue
    const c = q[ci]
    if (c === 'X' || c === 'Y' || c === 'MT') continue     // autosomes only
    m.set(q[pi], q[gi])
  }
  return m
}

const d = JSON.parse(readFileSync(JSONF, 'utf8'))
const pairs = (d.skipped ?? []).filter((s) => s.fraction < 0.065)  // one-genome, not a child
console.log(`${pairs.length} same-genome pairs available in ${JSONF.split('/').pop()}`)

const out = []
for (const s of pairs.slice(0, 40)) {
  const rf = files.find((f) => nameOf(f) === s.ref)
  const cf = files.find((f) => nameOf(f) === s.cell)
  if (!rf || !cf) continue
  const R = read(rf); const C = read(cf)
  let hom = 0; let het = 0
  for (const [k, gr] of R) {
    if (gr !== '0' && gr !== '2') continue      // reference homozygous: truth is homozygous
    const gc = C.get(k)
    if (!gc || gc === '-1') continue
    hom += 1
    if (gc === '1') het += 1                    // AB where truth is hom: drop-in
  }
  if (hom > 20000) out.push({ cell: s.cell, ref: s.ref, hom, dropIn: het / hom })
}
out.sort((a, b) => a.dropIn - b.dropIn)
for (const o of out) {
  console.log(`  ${o.cell} vs ${o.ref}  ${o.hom} hom markers  drop-in ${o.dropIn.toFixed(4)}`)
}
if (out.length) {
  const v = out.map((o) => o.dropIn).sort((a, b) => a - b)
  const med = v[v.length >> 1]
  console.log(`\nn=${v.length}  min ${v[0].toFixed(4)}  median ${med.toFixed(4)}  `
    + `max ${v[v.length - 1].toFixed(4)}`)
  console.log(`the consult assumed 0.0854 median from 7 arrays; this is ${v.length} pairs`)
}

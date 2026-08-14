/**
 * How many haploid products each donor group has, which decides where a parent can be
 * RECONSTRUCTED rather than arrayed.
 *
 * Reconstruction from haploid products is not the reconstruction PLAN.txt disqualifies. That one
 * is sibling and cross-donor sharing, which fails for the mother because siblings cannot reveal
 * where she is heterozygous: it admits those markers as homozygous and they then read as maternal
 * absence, turning a true log2 ratio of -0.40 into +2.18. Haploid products do not have that
 * failure, because one product reading A and another reading B PROVES the parent heterozygous.
 * That is why PLAN.txt names polar bodies as the correct substitute for a maternal array, and why
 * Progenitor scored 35 correct, 6 refused, 0 incorrect against a father it was never shown.
 *
 * MIN_PRODUCTS is 5. Below that the method inverts and true offspring read as decisively absent,
 * so a group with four is not a group with nearly enough.
 *
 *   usage: node --experimental-strip-types audit/asymmetry/haploid_census.ts <dir> [out.json]
 */
import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs'
import { join, basename } from 'node:path'

const W = new URL('../../web/src/', import.meta.url).pathname
const { headerMap, parseRow, emptyBafSums, accumulateBaf, accumulate, finishProfile } =
  await import(`${W}ingest.ts`)
const { inferStage } = await import(`${W}stage.ts`)
const { MIN_PRODUCTS } = await import(`${W}inferredReference.ts`)

const DIR = process.argv[2]
const OUT = process.argv[3] ?? 'audit/asymmetry/haploid-census.json'
const STRIDE = 8

function probesUnder(dir: string, rel = ''): string[] {
  const out: string[] = []
  for (const e of readdirSync(join(dir, rel)).sort()) {
    const r = rel ? join(rel, e) : e
    if (statSync(join(dir, r)).isDirectory()) out.push(...probesUnder(dir, r))
    else if (e.endsWith('.probes')) out.push(r)
  }
  return out
}

const files = probesUnder(DIR)
process.stderr.write(`${files.length} arrays\n`)
const rows: { id: string, group: string, stage: string, call: number, het: number }[] = []

for (const f of files) {
  const id = basename(f).replace(/_\d+\.CEL\.probes$/, '').replace(/\.probes$/, '')
  const lines = readFileSync(join(DIR, f), 'utf8').split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  if (h < 0) continue
  const map = headerMap(lines[h])
  if (!map) continue
  const byChrom = new Map()
  const baf = emptyBafSums()
  let first = ''
  for (let i = h + 1; i < lines.length; i += STRIDE) {
    const r = parseRow(lines[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    accumulate(r as never, byChrom as never)
    accumulateBaf(r as never, baf as never)
  }
  const p = finishProfile(id, byChrom as never, baf as never, first)
  const s = inferStage(p)
  // The group is the identifier stem the laboratory numbers by. Arrays of one experiment share it.
  const group = id.replace(/[-_].*$/, '')
  rows.push({ id, group, stage: s.stage, call: p.callRate, het: p.hetRate })
}

const byGroup = new Map<string, typeof rows>()
for (const r of rows) byGroup.set(r.group, [...(byGroup.get(r.group) ?? []), r])

const usable: string[] = []
process.stderr.write(`\ngroup      arrays  haploid  failed   verdict\n`)
for (const [g, rs] of [...byGroup].sort((a, b) => a[0].localeCompare(b[0]))) {
  const hap = rs.filter((r) => r.stage === 'haploid')
  const bad = rs.filter((r) => r.stage === 'failed')
  const ok = hap.length >= MIN_PRODUCTS
  if (ok) usable.push(g)
  process.stderr.write(`${g.padEnd(10)} ${String(rs.length).padStart(6)}  `
    + `${String(hap.length).padStart(7)}  ${String(bad.length).padStart(6)}   `
    + `${ok ? 'RECONSTRUCTABLE' : ''}\n`)
}
writeFileSync(OUT, JSON.stringify({ dir: DIR, minProducts: MIN_PRODUCTS, usable, rows }, null, 1))
process.stderr.write(`\n${usable.length} group(s) with ${MIN_PRODUCTS}+ haploid products -> ${OUT}\n`)

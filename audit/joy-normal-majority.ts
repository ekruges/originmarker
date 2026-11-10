/**
 * Does this genome have a NORMAL MAJORITY? Calibrated against published karyotypes.
 *
 * WHY. The copy-number scan normalises against the sample's own per-chromosome median, so it needs
 * most chromosomes to be unaffected: where they are not, the median is itself an affected
 * chromosome and the null is wrong. scoreSample currently withholds the scan when heterozygosity
 * passes a diploid ceiling, which is a proxy for that and a poor one in both directions.
 *
 * TOO STRICT. Three cells refused by the ceiling carry an information ratio of 0.323, 0.410 and
 * 0.429, inside the range of cells that PASS it (0.305, 0.330, 0.334, 0.356, 0.363). Their calls
 * carry as much genotype information as arrays the tool accepts, and they hold 11 of 26 hand-called
 * gains. The ceiling's own caveat says a mixture, a contaminated reaction and a failed array cannot
 * be told apart; 5.48.0's information ratio tells a failed array (R near 1) from a mixture (R low),
 * so that caveat is now out of date.
 *
 * TOO LOOSE. Two cells Zuccaro 2020 hand-called chaotic aneuploid are NOT refused, and the scan
 * calls 28 segments on chromosomes that paper publishes as normal.
 *
 * THE STATISTIC. Per autosome, the median log R minus the genome-wide median log R. The concordant
 * fraction is the share of autosomes whose deviation is within a band of that median. A genome with
 * a normal majority has most chromosomes in one tight cluster, and that cluster IS the normal set.
 * A chaotic genome has them scattered, so no majority exists and the median is not a reference.
 *
 * This run only MEASURES it against Zuccaro Table S3, whose karyotypes are published and reviewed,
 * so the band and the majority threshold are chosen from the separation rather than guessed. It
 * ships nothing.
 *
 * Run: OM_DIR=<dir of arrays> OM_TRUTH=audit/zuccaro-table-s3.tsv \
 *        node --experimental-strip-types audit/joy-normal-majority.ts
 */
import { readFileSync, existsSync, readdirSync, writeFileSync } from 'node:fs'
import { join } from 'node:path'
import { loadArray, isArray, type Marker } from './joy-common.ts'

const DIR = process.env.OM_DIR
const TRUTH = process.env.OM_TRUTH
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }

/** cell -> the published whole-cell karyotype, where there is one. */
const kary = new Map<string, string>()
const abnormalChroms = new Map<string, Set<string>>()
if (TRUTH && existsSync(TRUTH)) {
  for (const l of readFileSync(TRUTH, 'utf8').split('\n')) {
    if (!l || l.startsWith('#')) continue
    const c = l.split('\t')
    if (c.length < 5 || c[0] === 'gsm' || !c[0]) continue
    kary.set(c[0].trim(), c[2].trim())
    if (!abnormalChroms.has(c[0].trim())) abnormalChroms.set(c[0].trim(), new Set())
    if (c[4].trim().toLowerCase() !== 'normal') {
      abnormalChroms.get(c[0].trim())!.add(c[3].trim().replace(/^chr/i, ''))
    }
  }
}

const SHARDS = Number(process.env.OM_SHARDS ?? 1)
const SHARD = Number(process.env.OM_SHARD ?? 0)
const files = readdirSync(DIR).filter(isArray).sort().filter((_, i) => i % SHARDS === SHARD)

const BANDS = [0.15, 0.25, 0.35, 0.50]
const out: string[] = []
for (const f of files) {
  const cell = /^(GSM\d+)/.exec(f)?.[1] ?? (/^(?:GSM\d+_)?([0-9]+-[0-9]+)/.exec(f)?.[1]) ?? f
  let a
  try { a = loadArray(join(DIR, f), cell, { hetByChrom: false }) } catch (e) {
    console.log(`  ${cell}: ${(e as Error).message}`); continue
  }
  const byChrom = new Map<string, Marker[]>()
  for (const m of a.markers) {
    if (m.chrom.charCodeAt(0) < 48 || m.chrom.charCodeAt(0) > 57) continue
    if (!/^\d+$/.test(m.chrom)) continue
    let arr = byChrom.get(m.chrom)
    if (!arr) byChrom.set(m.chrom, arr = [])
    arr.push(m)
  }
  const med = (xs: number[]): number => {
    if (!xs.length) return NaN
    const s = [...xs].sort((x, y) => x - y)
    return s[s.length >> 1]
  }
  const devs: { chrom: string; dev: number }[] = []
  for (const [c, ms] of byChrom) {
    const v = ms.map((m) => m.log2R).filter((x): x is number => x !== null)
    if (v.length < 2_000) continue
    devs.push({ chrom: c, dev: med(v) - a.genomeLrr })
  }
  const fracs = BANDS.map((b) => devs.filter((d) => Math.abs(d.dev) <= b).length / devs.length)
  const truth = kary.get(cell) ?? ''
  const nAbnormal = abnormalChroms.get(cell)?.size ?? -1
  out.push([cell, a.callRate.toFixed(3), a.hetRate.toFixed(3), devs.length,
    ...fracs.map((x) => x.toFixed(3)),
    nAbnormal, truth].join('\t'))
}
const head = ['cell', 'call_rate', 'het_rate', 'autosomes_measured',
  ...BANDS.map((b) => `concordant_${b}`), 'published_abnormal_chroms', 'published_karyotype']
const text = [head.join('\t'), ...out].join('\n') + '\n'
if (process.env.OM_TSV) { writeFileSync(process.env.OM_TSV, text); console.log(`wrote ${process.env.OM_TSV}`) }
else console.log(text)

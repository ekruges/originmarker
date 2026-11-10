/**
 * Why each hand-called gain is missed, measured rather than inferred from his deviation column.
 *
 * joy-dieter.ts attributes a miss to SEGMENT_LRR_SHIFT when the deviation HE reports is under the
 * constant. That is not a diagnosis: his number is measured over his boundaries, and it says
 * nothing about whether the scan nominated a window there at all. Lowering the gain gate to
 * log2(3/2) recovered none of the four it blamed on the gate, and all four report a best overlap
 * of 0.00, meaning no segment was raised anywhere near them.
 *
 * The scan has two stages and only the second one is a gate. A window is NOMINATED on the no-call
 * indicator by scanChromosome against SEGMENT_LRT, and only then does the intensity gate judge it.
 * This separates them: for every hand region it reports the best no-call score any window
 * overlapping it can reach with the threshold removed, beside the intensity deviation the detector
 * itself measures there. A score far under SEGMENT_LRT means the nomination channel is blind and no
 * gate setting can matter.
 *
 * Run: OM_DIR=<dir of .probes> OM_REGIONS=<regions-dieter.tsv> \
 *        node --experimental-strip-types audit/joy-gain-miss-cause.ts
 */
import { readFileSync, writeFileSync, existsSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { loadArray, segs, isArray, type Marker } from './joy-common.ts'

const DIR = process.env.OM_DIR
const REG = process.env.OM_REGIONS
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }
if (!REG || !existsSync(REG)) { console.log('OM_REGIONS not readable. NOT RUN.'); process.exit(0) }

type Hand = { cell: string; chrom: string; start: number; end: number; kind: string; dev: number }
const hand: Hand[] = []
for (const l of readFileSync(REG, 'utf8').split('\n')) {
  if (!l || l.startsWith('cell')) continue
  const c = l.split('\t')
  if (c.length < 6) continue
  const kind = /gain/i.test(c[5]) ? 'gain' : /loss/i.test(c[5]) ? 'loss' : ''
  if (!kind) continue
  hand.push({ cell: c[0].trim(), chrom: c[1].trim(), start: Number(c[2]), end: Number(c[3]),
    kind, dev: Number(c[4]) })
}

const SHARDS = Number(process.env.OM_SHARDS ?? 1)
const SHARD = Number(process.env.OM_SHARD ?? 0)
const files = readdirSync(DIR).filter(isArray).sort().filter((_, i) => i % SHARDS === SHARD)

const out: string[] = []
for (const f of files) {
  const cell = (/^(?:GSM\d+_)?([0-9]+-[0-9]+)/.exec(f)?.[1]) ?? f.replace(/\..*$/, '')
  const mine = hand.filter((h) => h.cell === cell)
  if (!mine.length) continue
  const a = loadArray(join(DIR, f), cell)
  const byChrom = new Map<string, Marker[]>()
  for (const m of a.markers) {
    if (!/^\d+$/.test(m.chrom)) continue
    byChrom.set(m.chrom, [...(byChrom.get(m.chrom) ?? []), m])
  }
  const noCall = new Map<string, [number, number]>()
  for (const [c, ms] of byChrom) noCall.set(c, [ms.length, ms.filter((m) => !m.called).length])

  for (const h of mine) {
    const ms = byChrom.get(h.chrom)
    if (!ms) { out.push([cell, h.chrom, h.kind, h.dev, 'NO MARKERS'].join('\t')); continue }
    const nul = segs.externalNull(noCall, h.chrom) as number
    const inside = ms.filter((m) => m.pos >= h.start && m.pos <= h.end)
    const ncIn = inside.length ? inside.filter((m) => !m.called).length / inside.length : NaN
    const ncChr = ms.filter((m) => !m.called).length / ms.length
    const lrrs = inside.map((m) => m.log2R).filter((x): x is number => x !== null).sort((x, y) => x - y)
    const delta = lrrs.length ? lrrs[lrrs.length >> 1] - a.genomeLrr : NaN

    // THE NOMINATION STAGE WITH THE THRESHOLD REMOVED. Same function, same null, same marker
    // floor; only the score bar is dropped to zero, so every window it can build is returned and
    // the best one overlapping his region shows what score the channel can actually reach there.
    const probe = segs.scanChromosome(
      ms.map((m) => ({ chrom: m.chrom, pos: m.pos, absent: !m.called })), nul,
      segs.MIN_SEGMENT_MARKERS, 0,
    ) as { startBp: number; endBp: number; score: number }[]
    let best = 0
    let bestSpan = ''
    for (const p of probe) {
      const ov = Math.max(0, Math.min(p.endBp, h.end) - Math.max(p.startBp, h.start))
      if (ov <= 0) continue
      if (p.score > best) {
        best = p.score
        bestSpan = `${(p.startBp / 1e6).toFixed(1)}-${(p.endBp / 1e6).toFixed(1)}`
      }
    }
    out.push([cell, h.chrom, h.kind, h.dev.toFixed(3), ((h.end - h.start) / 1e6).toFixed(1),
      inside.length, ncIn.toFixed(3), ncChr.toFixed(3), nul.toFixed(4), delta.toFixed(3),
      best.toFixed(0), bestSpan].join('\t'))
  }
}
const head = ['cell', 'chrom', 'kind', 'his_dev', 'mb', 'markers_in', 'nocall_in', 'nocall_chrom',
  'null', 'our_delta', 'best_nocall_score', 'best_window_mb'].join('\t')
const text = [head, ...out].join('\n') + '\n'
if (process.env.OM_TSV) { writeFileSync(process.env.OM_TSV, text); console.log(`wrote ${process.env.OM_TSV}`) }
else console.log(text)

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
 * IT ALSO PROBES THE REPLACEMENT INDICATOR, before any of it is shipped. The scanner takes a
 * binary per-marker flag against a rate null, so nominating on intensity needs a different flag
 * and no new machinery: a marker counts when its log R sits above the quantile of this array's own
 * autosomal distribution, which makes the null self-calibrating per array and introduces no
 * absolute intensity threshold, since a cohort offset does not transfer between laboratories.
 * Three quantiles are reported so the cut is chosen from the measurement rather than guessed.
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
  const a = loadArray(join(DIR, f), cell, { hetByChrom: false })
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
    // This array's own autosomal log R distribution, so the cut is relative to the sample.
    const allLrr = [...byChrom.values()].flat()
      .map((m) => m.log2R).filter((x): x is number => x !== null).sort((x, y) => x - y)
    const lrrQuantile = (pq: number): number =>
      allLrr[Math.min(allLrr.length - 1, Math.floor(pq * allLrr.length))]
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
    // THE REPLACEMENT INDICATOR, probed at three cuts. Same scanner, same marker floor, same
    // null machinery; only the flag and its rate change. Threshold dropped to zero again, so the
    // score reported is what the channel can reach rather than whether it passed.
    const probeHigh: string[] = []
    for (const Q of [0.70, 0.80, 0.90]) {
      const cut = lrrQuantile(Q)
      const highByChrom = new Map<string, [number, number]>()
      for (const [c, cms] of byChrom) {
        const v = cms.filter((m) => m.log2R !== null)
        highByChrom.set(c, [v.length, v.filter((m) => (m.log2R as number) > cut).length])
      }
      const hiNull = segs.externalNull(highByChrom, h.chrom) as number
      const vIn = inside.filter((m) => m.log2R !== null)
      const hiIn = vIn.length ? vIn.filter((m) => (m.log2R as number) > cut).length / vIn.length : NaN
      const hits = segs.scanChromosome(
        ms.map((m) => ({ chrom: m.chrom, pos: m.pos,
          absent: m.log2R !== null && (m.log2R as number) > cut })),
        hiNull, segs.MIN_SEGMENT_MARKERS, 0,
      ) as { startBp: number; endBp: number; score: number }[]
      let hiBest = 0
      for (const t of hits) {
        if (Math.min(t.endBp, h.end) - Math.max(t.startBp, h.start) <= 0) continue
        if (t.score > hiBest) hiBest = t.score
      }
      probeHigh.push(hiIn.toFixed(3), hiNull.toFixed(4), hiBest.toFixed(0))
    }

    out.push([cell, h.chrom, h.kind, h.dev.toFixed(3), ((h.end - h.start) / 1e6).toFixed(1),
      inside.length, ncIn.toFixed(3), ncChr.toFixed(3), nul.toFixed(4), delta.toFixed(3),
      best.toFixed(0), bestSpan, ...probeHigh].join('\t'))
  }
}
const head = ['cell', 'chrom', 'kind', 'his_dev', 'mb', 'markers_in', 'nocall_in', 'nocall_chrom',
  'null', 'our_delta', 'best_nocall_score', 'best_window_mb',
  'hi70_in', 'hi70_null', 'hi70_score', 'hi80_in', 'hi80_null', 'hi80_score',
  'hi90_in', 'hi90_null', 'hi90_score'].join('\t')
const text = [head, ...out].join('\n') + '\n'
if (process.env.OM_TSV) { writeFileSync(process.env.OM_TSV, text); console.log(`wrote ${process.env.OM_TSV}`) }
else console.log(text)

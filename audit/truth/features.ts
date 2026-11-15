/**
 * The raw signals of each array of a list, per chromosome and per 5 Mb window, with no parent
 * consulted: how many markers, how many called, how many heterozygous, the median log R, and where
 * the B allele frequencies sit. Written once so a statistic can be tried against the stated facts
 * offline, without the arrays.
 *
 * Per chromosome the B allele frequency is a 20-bin histogram for all markers, for heterozygous
 * calls and for no-calls. Per window it is five bands: [0,.15) [.15,.42) [.42,.58] (.58,.85] (.85,1],
 * for all markers and for heterozygous calls. Log R is a histogram too: 60 bins of 0.1 from -3 per
 * chromosome, 16 bins of 0.25 from -2.5 per window, the ends open. Where the file carries the array
 * software's own copy-number state per marker, its counts are kept per chromosome and per window.
 *
 * Run: node --experimental-strip-types audit/truth/features.ts <list.tsv> <shard> <shards> <out.jsonl>
 * The list is tab-separated and its first two fields are the array id and its path.
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { ingest, stageMod, lineReader } from '../joy-common.ts'

const R = new URL('../../', import.meta.url).pathname
const keys = await import(R + 'web/src/answerKeys.ts')
const info = await import(R + 'web/src/informationRatio.ts')
const panel = keys.parsePanel(gunzipSync(readFileSync(R + 'web/public/' + keys.FOUNDER_PANEL_PATH)).toString('utf8'))
const [listFile, shard, shards, out] = [process.argv[2], Number(process.argv[3]), Number(process.argv[4]), process.argv[5]]
const done = new Set<string>()
if (existsSync(out)) for (const l of readFileSync(out, 'utf8').split('\n')) if (l.startsWith('{')) done.add(JSON.parse(l).id)
const list = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t')).filter((_, i) => i % shards === shard)

export const WINDOW_BP = 5_000_000
export const band = (b: number): number => (b < 0.15 ? 0 : b < 0.42 ? 1 : b <= 0.58 ? 2 : b <= 0.85 ? 3 : 4)
const median = (a: number[]): number | null => { if (!a.length) return null; a.sort((x, y) => x - y); return a[a.length >> 1] }
const round = (x: number | null) => (x === null ? null : Math.round(x * 1e4) / 1e4)

for (const [id, path] of list) {
  if (done.has(id)) continue
  try {
    const raw = readFileSync(path)
    const next = lineReader(path.endsWith('.gz') ? gunzipSync(raw) : raw)
    let header: any = null; let cnCol = -1
    const byChromP = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    const sameHom = keys.emptySameHom(); const homTally = info.emptyHomTally()
    type C = { n: number; called: number; het: number; lrr: number[]; all: number[]; hetB: number[]; nc: number[]; lh: number[]; cn: number[] }
    type Wn = { n: number; called: number; het: number; lrr: number[]; all: number[]; hetB: number[]; lh: number[]; cn: number[] }
    const chroms = new Map<string, C>(); const windows = new Map<string, Wn>()
    for (let line = next(); line !== null; line = next()) {
      if (!header) {
        if (line && !line.startsWith('#')) { header = ingest.headerMap(line); cnCol = line.split('\t').indexOf('copy_number') }
        continue
      }
      const r = ingest.parseRow(line, header); if (!r) continue
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChromP); ingest.accumulateBaf(r, bafSums)
      if (!/^(\d+|X)$/.test(r.chrom)) continue
      if (r.chrom !== 'X') {
        const p = panel.get(r.probesetId)
        if (p) { keys.addSameHom(sameHom, r.genotype, p, p); info.addHomozygous(homTally, r.genotype, p) }
      }
      let c = chroms.get(r.chrom)
      if (!c) chroms.set(r.chrom, c = { n: 0, called: 0, het: 0, lrr: [], all: new Array(20).fill(0), hetB: new Array(20).fill(0), nc: new Array(20).fill(0), lh: new Array(60).fill(0), cn: [0, 0, 0, 0, 0] })
      const wk = `${r.chrom}:${Math.floor(r.pos / WINDOW_BP)}`
      let w = windows.get(wk)
      if (!w) windows.set(wk, w = { n: 0, called: 0, het: 0, lrr: [], all: [0, 0, 0, 0, 0], hetB: [0, 0, 0, 0, 0], lh: new Array(16).fill(0), cn: [0, 0, 0, 0, 0] })
      c.n++; w.n++
      const calledNow = r.genotype !== 'NC'; const het = r.genotype === 'AB'
      if (calledNow) { c.called++; w.called++ }
      if (het) { c.het++; w.het++ }
      if (r.log2R !== null && r.log2R !== undefined && Number.isFinite(r.log2R)) {
        c.lrr.push(r.log2R); w.lrr.push(r.log2R)
        c.lh[Math.min(59, Math.max(0, Math.floor((r.log2R + 3) * 10)))]++
        w.lh[Math.min(15, Math.max(0, Math.floor((r.log2R + 2.5) * 4)))]++
      }
      if (cnCol >= 0) {
        const state = Number(line.split('\t')[cnCol])
        if (Number.isFinite(state)) { const k = Math.min(4, Math.max(0, Math.round(state))); c.cn[k]++; w.cn[k]++ }
      }
      if (r.baf !== null && r.baf !== undefined && Number.isFinite(r.baf)) {
        const bin = Math.min(19, Math.max(0, Math.floor(r.baf * 20))); const b = band(r.baf)
        c.all[bin]++; w.all[b]++
        if (het) { c.hetB[bin]++; w.hetB[b]++ }
        if (!calledNow) c.nc[bin]++
      }
    }
    const prof = ingest.finishProfile(id, byChromP, bafSums, first)
    const dropIn = keys.dropIn(sameHom); const ratio = info.informationRatio(homTally)
    const st = stageMod.inferStage({ ...prof, dropIn: keys.onParentScale(dropIn, panel.calibration), information: ratio }) as any
    appendFileSync(out, JSON.stringify({
      id, callRate: round(prof.callRate), hetRate: round(prof.hetRate), dropIn: round(dropIn), information: round(ratio),
      stage: st.stage, refusedBy: st.refusedBy ?? '', ploidy: st.ploidy ?? '',
      chroms: Object.fromEntries([...chroms].map(([k, c]) => [k, { n: c.n, called: c.called, het: c.het, lrr: round(median(c.lrr)),
        all: c.all, het_baf: c.hetB, nc_baf: c.nc, lrr_hist: c.lh, cn: cnCol >= 0 ? c.cn : null }])),
      windows: [...windows].map(([k, w]) => [k, w.n, w.called, w.het, round(median(w.lrr)), ...w.all, ...w.hetB]),
      window_lrr: [...windows].map(([k, w]) => [k, ...w.lh]),
      window_cn: cnCol >= 0 ? [...windows].map(([k, w]) => [k, ...w.cn]) : null,
    }) + '\n')
  } catch (e) { appendFileSync(out, JSON.stringify({ id, error: String(e).slice(0, 200) }) + '\n') }
}
appendFileSync(out, `DONE ${shard}\n`)

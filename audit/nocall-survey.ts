/**
 * Every array: its call rate, and what its no-calls physically are. Signal (log R) and allele balance
 * (BAF) at the no-calls against the called markers, per array and per autosome, with the stage the
 * shipped code gives it without a parent. One line per array, one per array-chromosome. The survey
 * behind audit/call-floor-preregistration.txt.
 *
 * Run: node --experimental-strip-types audit/nocall-survey.ts <file of array paths> <shard> <shards> <out>
 */
import { readFileSync, writeFileSync, appendFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
const R = new URL('../', import.meta.url).pathname
const ingest = await import(R + 'web/src/ingest.ts')
const keys = await import(R + 'web/src/answerKeys.ts')
const info = await import(R + 'web/src/informationRatio.ts')
const stage = await import(R + 'web/src/stage.ts')
const panel = keys.parsePanel(gunzipSync(readFileSync(R + 'web/public/' + keys.FOUNDER_PANEL_PATH)).toString('utf8'))
const [listFile, shard, shards, out] = process.argv.slice(2)
const files = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).filter((_, i) => i % Number(shards) === Number(shard))
const q = (a: number[], p: number) => { if (!a.length) return NaN; const s = [...a].sort((x, y) => x - y); return s[Math.min(s.length - 1, Math.floor(p * s.length))] }
for (const path of files) {
  try {
    const raw = readFileSync(path); const text = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
    const lines = text.split('\n'); let h = 0; while (h < lines.length && (!lines[h] || lines[h].startsWith('#'))) h++
    const map = ingest.headerMap(lines[h]); if (!map) { appendFileSync(out, `ERR\t${path}\tno header\n`); continue }
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    const sh = keys.emptySameHom(); const it = info.emptyHomTally()
    const ncL: number[] = [], cL: number[] = []; let ncN = 0, ncBafMid = 0, ncBafEdge = 0, ncBafNaN = 0, ncLNaN = 0
    const per = new Map<string, { n: number; nc: number; ncL: number[]; cL: number[]; mid: number }>()
    for (let i = h + 1; i < lines.length; i++) {
      const r = ingest.parseRow(lines[i], map); if (!r) continue
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums)
      if (!/^\d+$/.test(r.chrom)) continue
      const p = panel.get(r.probesetId); if (p) { keys.addSameHom(sh, r.genotype, p, p); info.addHomozygous(it, r.genotype, p) }
      let c = per.get(r.chrom); if (!c) per.set(r.chrom, c = { n: 0, nc: 0, ncL: [], cL: [], mid: 0 })
      c.n++
      const L = r.log2R
      if (r.genotype === 'NC') {
        ncN++; c.nc++
        if (L === null || !Number.isFinite(L)) ncLNaN++; else { if (ncL.length < 200000) ncL.push(L); if (c.ncL.length < 20000) c.ncL.push(L) }
        if (r.baf === null || !Number.isFinite(r.baf)) ncBafNaN++
        else if (r.baf > 0.2 && r.baf < 0.8) { ncBafMid++; c.mid++ } else if (r.baf < 0.05 || r.baf > 0.95) ncBafEdge++
      } else if (L !== null && Number.isFinite(L)) { if (cL.length < 200000) cL.push(L); if (c.cL.length < 20000) c.cL.push(L) }
    }
    const prof = ingest.finishProfile(path, byChrom, bafSums, first)
    const d = keys.dropIn(sh); const Rr = info.informationRatio(it)
    const st = stage.inferStage({ ...prof, dropIn: keys.onParentScale(d, panel.calibration), information: Rr })
    const row = [path.split('/').slice(-2).join('/'), prof.callRate.toFixed(4), prof.hetRate.toFixed(4), Number(prof.hetBand).toFixed(4),
      d.toFixed(4), Number.isFinite(Rr) ? Rr.toFixed(3) : 'nan', st.stage, st.refusedBy ?? '', st.ploidy,
      ncN, q(cL, 0.5).toFixed(3), q(ncL, 0.1).toFixed(3), q(ncL, 0.5).toFixed(3), q(ncL, 0.9).toFixed(3),
      (ncBafMid / Math.max(1, ncN)).toFixed(4), (ncBafEdge / Math.max(1, ncN)).toFixed(4), (ncBafNaN / Math.max(1, ncN)).toFixed(4), (ncLNaN / Math.max(1, ncN)).toFixed(4)]
    appendFileSync(out, 'ARR\t' + row.join('\t') + '\n')
    for (const [c, v] of per) appendFileSync(out, ['CHR', row[0], c, v.n, (1 - v.nc / v.n).toFixed(4), q(v.cL, 0.5).toFixed(3), q(v.ncL, 0.5).toFixed(3), (v.mid / Math.max(1, v.nc)).toFixed(4)].join('\t') + '\n')
  } catch (e) { appendFileSync(out, `ERR\t${path}\t${String(e).slice(0, 120)}\n`) }
}
appendFileSync(out, 'DONE\n')

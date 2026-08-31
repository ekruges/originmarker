/**
 * The intensity statistic, old against new, on chromosomes that carry nothing.
 *
 * THE TRUE-NEGATIVE SET IS FREE AND EXACT. In a pronucleus every autosome is at ONE copy,
 * uniformly, so no chromosome can differ from any other. Every chromosome flagged is a false
 * positive by construction, with no injection, no modelling and no assumption about noise. That is
 * what makes this the right test for an error model: any harness that CONSTRUCTS an event has to
 * assume the very noise it is trying to measure.
 *
 * MEASURED, 264 event-free chromosomes over 12 pronuclei of GSE148488:
 *
 *   statistic                                  flagged   rate
 *   old, sd/sqrt(n), |z| > 2.576                 234     0.8864
 *   new, per-array MAD null, |z| > 5.49            4     0.0152
 *
 * The old statistic had a median |z| of 17.4 and a maximum of 539 on chromosomes carrying nothing,
 * and 138 of its 234 false positives, 0.590, had POSITIVE sign. A positive shift resolves the
 * copy-number class as a gain; a gain inverts the sign map that loss and copy-neutral share; so the
 * parent came out backwards on about half of a very large number of spurious calls. That is the
 * inversion mechanism, measured on real arrays with a known answer rather than argued.
 *
 * Run: OM_SRC=<web/src/> OM_TRIOS=<dir of GSE148488> node --experimental-strip-types <this file>
 */
// THE TRUE-NEGATIVE TEST. In a pronucleus every autosome is at ONE copy, uniformly, so no
// chromosome can differ from any other. Every flagged chromosome is a false positive by
// construction. Old statistic against new, on the same real arrays.
import { readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
const W = process.env.OM_SRC!
const ingest = await import(`${W}ingest.ts`)
const nul = await import(`${W}intensityNull.ts`)
const D = process.env.OM_TRIOS!

const PN = ['GSM4774673','GSM4774674','GSM4774675','GSM4774676','GSM4774677','GSM4774678',
  'GSM4774679','GSM4774680','GSM4774681','GSM4774682','GSM4774683','GSM4774685']
const CH = Array.from({length:22},(_,i)=>String(i+1))

let oldFlag=0, newFlag=0, total=0
const oldZs:number[]=[], newZs:number[]=[]
let oldPos=0
for (const id of PN) {
  const text = gunzipSync(readFileSync(`${D}/${id}.probes.gz`)).toString('utf8')
  const lines = text.split('\n'); let h=-1
  for (let i=0;i<60;i++) if (lines[i] && !lines[i].startsWith('#')) { h=i; break }
  const map = ingest.headerMap(lines[h])
  const byChrom = new Map<string, number[]>()
  for (let i=h+1;i<lines.length;i++) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || !/^\d+$/.test(r.chrom) || r.log2R === null || !Number.isFinite(r.log2R)) continue
    const a = byChrom.get(r.chrom); if (a) a.push(r.log2R); else byChrom.set(r.chrom, [r.log2R])
  }
  const meds = CH.map(c => nul.median(byChrom.get(c) ?? [])).filter(Number.isFinite)
  const scale = nul.nullScale(meds)
  for (const c of CH) {
    const inL = byChrom.get(c) ?? []
    if (inL.length < 200) continue
    const outL = CH.filter(x=>x!==c).flatMap(x=>byChrom.get(x) ?? [])
    if (outL.length < 1000) continue
    total++
    // OLD: difference of means over sd/sqrt(n), iid.
    const mean=(xs:number[])=>xs.reduce((a,x)=>a+x,0)/xs.length
    const mu=mean(outL)
    const sd=Math.sqrt(outL.reduce((a,x)=>a+(x-mu)**2,0)/(outL.length-1))
    const zOld=(mean(inL)-mu)/(sd/Math.sqrt(inL.length))
    // NEW: region median against the array's own chromosome units, MAD scale.
    const zNew=nul.calibratedZ(nul.median(inL), scale)
    oldZs.push(Math.abs(zOld)); if (zNew!==undefined) newZs.push(Math.abs(zNew))
    if (Math.abs(zOld) > 2.576) { oldFlag++; if (zOld > 0) oldPos++ }
    if (nul.intensityDetects(zNew, true)) newFlag++
  }
}
const pct=(a:number[],p:number)=>{const q=[...a].sort((x,y)=>x-y);return q[Math.floor(q.length*p)]}
console.log(`event-free chromosomes tested: ${total} over ${PN.length} pronuclei`)
console.log(`OLD iid z > 2.576 : ${oldFlag} flagged = ${(oldFlag/total).toFixed(4)}`)
console.log(`     median |z| ${pct(oldZs,0.5).toFixed(1)}  max ${Math.max(...oldZs).toFixed(0)}`)
console.log(`     of those flagged, ${oldPos} had POSITIVE sign = ${(oldPos/Math.max(1,oldFlag)).toFixed(3)} (read as GAINS)`)
console.log(`NEW calibrated z > ${nul.Z_CHROMOSOME}: ${newFlag} flagged = ${(newFlag/total).toFixed(4)}`)
console.log(`     median |z| ${pct(newZs,0.5).toFixed(2)}  99th pct ${pct(newZs,0.99).toFixed(2)}`)

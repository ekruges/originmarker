// Does the RAW posterior err under-confident or over-confident on REAL array noise?
//
// This decides whether shipping without a calibration map is safe. Under-confident means we
// understate, which is conservative. Over-confident means we overstate, which is the failure the
// whole rework exists to remove. The consult asserts under-confident; this measures it here.
//
// Real noise, synthetic event. The observed shift is a true displacement from the pooled-dosage
// algebra PLUS an actual per-chromosome null shift measured on a real array of this corpus, so the
// dispersion, the drift and the amplification artefacts are the arrays' own.
import { readFileSync } from 'node:fs'
const W = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/src/'
const { originPosterior, shiftMean } = await import(`${W}originPosterior.ts`)
const SP = '/private/tmp/claude-501/-Users-ezrakruger-claudecodegeneralworkspace/d993bfa4-9332-4450-a429-0951e6eeb4b2/scratchpad/'

const rows = readFileSync(SP + 'null_shifts.csv', 'utf8').trim().split('\n').slice(1)
const byArray = new Map<string, number[]>()
for (const l of rows) {
  const p = l.split(','); const v = +p[3]
  if (!Number.isFinite(v)) continue
  if (!byArray.has(p[0])) byArray.set(p[0], [])
  byArray.get(p[0])!.push(v)
}
const arrays = [...byArray].filter(([, v]) => v.length >= 15)
const sd = (xs: number[]) => { const m = xs.reduce((a,b)=>a+b,0)/xs.length
  return Math.sqrt(xs.reduce((a,b)=>a+(b-m)*(b-m),0)/(xs.length-1)) }

// Deterministic pseudo-random, so the run reproduces exactly.
let seed = 12345
const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff

const CLASSES = ['loss','gain','cnn-loh'] as const
const results: {conf:number, correct:boolean, cls:string, f:number, tier:string}[] = []
for (const [name, nulls] of arrays) {
  const s = sd(nulls)
  const tier = s < 0.004 ? 'clean' : s < 0.012 ? 'mid' : 'noisy'
  for (let i = 0; i < 4000; i++) {
    const cls = CLASSES[Math.floor(rnd()*3)]
    const f = 0.02 + rnd()*0.68
    const affected = rnd() < 0.5 ? 'loaded' : 'other'
    const noise = nulls[Math.floor(rnd()*nulls.length)]
    const observed = shiftMean(cls, f, affected) + noise
    const p = originPosterior({ shift: observed, shiftSd: s, material: 'trophectoderm', markers: 800 })
    if (!Number.isFinite(p.confidence)) continue
    results.push({ conf: p.confidence, correct: p.parent === affected, cls, f, tier })
  }
}
console.log(`${results.length} injections over ${arrays.length} real arrays\n`)

const report = (label: string, rs: typeof results) => {
  if (rs.length < 50) return
  const bins = [[0.5,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,0.97],[0.97,0.999],[0.999,1.0001]]
  let ece = 0, signed = 0, tot = 0
  console.log(`  ${label}  (n=${rs.length})`)
  console.log('    stated        n     stated    actual    gap')
  for (const [lo,hi] of bins) {
    const b = rs.filter(r=>r.conf>=lo && r.conf<hi)
    if (b.length < 30) continue
    const stated = b.reduce((a,r)=>a+r.conf,0)/b.length
    const actual = b.filter(r=>r.correct).length/b.length
    ece += b.length*Math.abs(stated-actual); signed += b.length*(stated-actual); tot += b.length
    console.log(`    ${lo.toFixed(3)}-${hi.toFixed(3)} ${String(b.length).padStart(7)}   ${stated.toFixed(4)}    ${actual.toFixed(4)}   ${(stated-actual>=0?'+':'')}${(stated-actual).toFixed(4)}`)
  }
  console.log(`    ECE ${(ece/tot).toFixed(4)}   signed ${(signed/tot>=0?'+':'')}${(signed/tot).toFixed(4)}  -> ${signed/tot>0?'OVER-confident':'UNDER-confident'}\n`)
}
report('all', results)
for (const t of ['clean','mid','noisy']) report(`noise tier: ${t}`, results.filter(r=>r.tier===t))
for (const c of CLASSES) report(`class: ${c}`, results.filter(r=>r.cls===c))

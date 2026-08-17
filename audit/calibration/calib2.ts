// Which OBSERVABLE predicts the confidently-inverted gain rows? If one does, it is a gate. If none
// does, the honest response is a structural cap rather than a threshold.
import { readFileSync } from 'node:fs'
const W = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/src/'
const { originPosterior, shiftMean } = await import(`${W}originPosterior.ts`)
const SP = '/private/tmp/claude-501/-Users-ezrakruger-claudecodegeneralworkspace/d993bfa4-9332-4450-a429-0951e6eeb4b2/scratchpad/'
const rows = readFileSync(SP + 'null_shifts.csv','utf8').trim().split('\n').slice(1)
const byArray = new Map<string, number[]>()
for (const l of rows) { const p=l.split(','); const v=+p[3]
  if (Number.isFinite(v)) { if(!byArray.has(p[0]))byArray.set(p[0],[]); byArray.get(p[0])!.push(v) } }
const arrays=[...byArray].filter(([,v])=>v.length>=15)
const sd=(xs:number[])=>{const m=xs.reduce((a,b)=>a+b,0)/xs.length
  return Math.sqrt(xs.reduce((a,b)=>a+(b-m)*(b-m),0)/(xs.length-1))}
let seed=12345; const rnd=()=>(seed=(seed*1103515245+12345)&0x7fffffff)/0x7fffffff
const CLASSES=['loss','gain','cnn-loh'] as const
type R={conf:number,correct:boolean,cls:string,resolved:string,hasIntensity:boolean}
const res:R[]=[]
for(const [,nulls] of arrays){ const s=sd(nulls)
  for(let i=0;i<4000;i++){
    const cls=CLASSES[Math.floor(rnd()*3)]; const f=0.02+rnd()*0.68
    const affected=rnd()<0.5?'loaded':'other'
    const observed=shiftMean(cls,f,affected)+nulls[Math.floor(rnd()*nulls.length)]
    // Half the injections get an intensity channel, so the effect of HAVING one is measurable.
    const withI = i%2===0
    const logRSd = 0.02
    const logR = withI ? (cls==='cnn-loh'?0:cls==='gain'?Math.log2((2+f)/2):Math.log2((2-f)/2))
      + (rnd()-0.5)*logRSd*2 : undefined
    const p=originPosterior({shift:observed,shiftSd:s,logR,logRSd:withI?logRSd:undefined,
      material:'trophectoderm',markers:800})
    if(!Number.isFinite(p.confidence))continue
    res.push({conf:p.confidence,correct:p.parent===affected,cls,resolved:p.classResolved,hasIntensity:withI})
  }}
const acc=(rs:R[])=>rs.length?rs.filter(r=>r.correct).length/rs.length:NaN
const show=(label:string,rs:R[])=>{
  const top=rs.filter(r=>r.conf>=0.90)
  console.log(`  ${label.padEnd(38)} n=${String(rs.length).padStart(6)}  top-band n=${String(top.length).padStart(6)}  acc(top)=${top.length?acc(top).toFixed(4):'   n/a'}`)
}
console.log(`${res.length} injections\n`)
console.log('SPLIT BY WHETHER THE INTENSITY CHANNEL WAS PRESENT')
show('no intensity, all classes', res.filter(r=>!r.hasIntensity))
show('no intensity, TRUE class = gain', res.filter(r=>!r.hasIntensity&&r.cls==='gain'))
show('with intensity, all classes', res.filter(r=>r.hasIntensity))
show('with intensity, TRUE class = gain', res.filter(r=>r.hasIntensity&&r.cls==='gain'))
console.log('\nSPLIT BY WHAT THE POSTERIOR ITSELF REPORTS ABOUT THE CLASS')
for(const rv of ['unresolved','loss','gain','cnn-loh'])
  show(`classResolved = ${rv}`, res.filter(r=>r.resolved===rv))
console.log('\nTHE CANDIDATE GATE: class unresolved, on amplified material')
show('resolved, top band', res.filter(r=>r.resolved!=='unresolved'))
show('UNRESOLVED, top band', res.filter(r=>r.resolved==='unresolved'))

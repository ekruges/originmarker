// A gain's displacement cannot exceed f/(4+2f), which at the f=0.70 ceiling is 0.1296. Above that
// the gain hypothesis is excluded by the algebra and the call is safe. Below it a gain is feasible,
// and with no intensity channel nothing can exclude one. Is "gain feasible AND blind" a usable
// gate, or does it cost more than it saves?
import { readFileSync } from 'node:fs'
const W='/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/src/'
const {originPosterior,shiftMean,shiftMagnitude}=await import(`${W}originPosterior.ts`)
const SP='/private/tmp/claude-501/-Users-ezrakruger-claudecodegeneralworkspace/d993bfa4-9332-4450-a429-0951e6eeb4b2/scratchpad/'
const rows=readFileSync(SP+'null_shifts.csv','utf8').trim().split('\n').slice(1)
const byArray=new Map<string,number[]>()
for(const l of rows){const p=l.split(',');const v=+p[3]
  if(Number.isFinite(v)){if(!byArray.has(p[0]))byArray.set(p[0],[]);byArray.get(p[0])!.push(v)}}
const arrays=[...byArray].filter(([,v])=>v.length>=15)
const sd=(xs:number[])=>{const m=xs.reduce((a,b)=>a+b,0)/xs.length
  return Math.sqrt(xs.reduce((a,b)=>a+(b-m)*(b-m),0)/(xs.length-1))}
let seed=999; const rnd=()=>(seed=(seed*1103515245+12345)&0x7fffffff)/0x7fffffff
const CLASSES=['loss','gain','cnn-loh'] as const
const GAIN_MAX = shiftMagnitude('gain', 0.70)
console.log(`gain ceiling at f=0.70: |shift| <= ${GAIN_MAX.toFixed(4)}. Above it no gain is feasible.\n`)
let bandA=0, bandAwrong=0, gated=0, gatedWrong=0
for(const [,nulls] of arrays){ const s=sd(nulls)
  for(let i=0;i<6000;i++){
    const cls=CLASSES[Math.floor(rnd()*3)]; const f=0.02+rnd()*0.68
    const affected=rnd()<0.5?'loaded':'other'
    const observed=shiftMean(cls,f,affected)+nulls[Math.floor(rnd()*nulls.length)]
    const p=originPosterior({shift:observed,shiftSd:s,material:'trophectoderm',markers:800})
    if(!Number.isFinite(p.confidence)||p.band!=='A')continue
    const correct=p.parent===affected
    bandA++; if(!correct)bandAwrong++
    if(Math.abs(observed)<=GAIN_MAX){ gated++; if(!correct)gatedWrong++ }
  }}
console.log('BLIND (no intensity), band A only:')
console.log(`  band A rows                    ${bandA}`)
console.log(`  of them wrong                  ${bandAwrong}  (${(100*bandAwrong/bandA).toFixed(3)}%)`)
console.log(`\n  gate = |shift| <= ${GAIN_MAX.toFixed(4)} (a gain is feasible):`)
console.log(`  rows the gate would demote     ${gated}  (${(100*gated/bandA).toFixed(1)}% of band A)`)
console.log(`  wrong rows it catches          ${gatedWrong} of ${bandAwrong}  (${bandAwrong?(100*gatedWrong/bandAwrong).toFixed(1):'0'}%)`)
console.log(`  correct rows it destroys       ${gated-gatedWrong}`)
console.log(`\n  band A accuracy before gate    ${((bandA-bandAwrong)/bandA).toFixed(5)}`)
const after = bandA-gated
console.log(`  band A accuracy after gate     ${after?((after-(bandAwrong-gatedWrong))/after).toFixed(5):'n/a'}`)

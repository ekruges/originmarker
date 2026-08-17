// Stage every censused array with the tool's OWN staging function, so the material each injection
// is attributed to is the material the tool would assign at runtime rather than a label I chose.
import { readFileSync, writeFileSync } from 'node:fs'
const W = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/src/'
const stageMod = await import(`${W}stage.ts`)
const dosage = await import(`${W}dosageOrigin.ts`)
const SP = '/private/tmp/claude-501/-Users-ezrakruger-claudecodegeneralworkspace/d993bfa4-9332-4450-a429-0951e6eeb4b2/scratchpad/'

const lines = readFileSync(SP + 'census.csv', 'utf8').trim().split('\n')
const head = lines[0].split(',')
const rows = lines.slice(1).map((l) => {
  const p = l.split(','); const o: Record<string, string> = {}
  head.forEach((h, i) => { o[h] = p[i] }); return o
}).filter((r) => r.callRate)

const byMaterial = new Map<string, {name:string,callRate:number,hetRate:number,hetBafSd:number}[]>()
for (const r of rows) {
  const callRate = +r.callRate, hetRate = +r.hetRate, hetBafSd = +r.hetBafSd
  const st = stageMod.inferStage({ hetRate, callRate })
  const material = dosage.materialOf(st.stage)
  if (!byMaterial.has(material)) byMaterial.set(material, [])
  byMaterial.get(material)!.push({ name: r.name, callRate, hetRate, hetBafSd })
}
console.log('material         n     hetBafSd median   callRate median')
const picked: Record<string, string[]> = {}
for (const [m, list] of [...byMaterial].sort((a,b)=>b[1].length-a[1].length)) {
  const med = (xs: number[]) => xs.slice().sort((a,b)=>a-b)[Math.floor(xs.length/2)]
  console.log(`${m.padEnd(16)} ${String(list.length).padEnd(5)} ${med(list.map(x=>x.hetBafSd)).toFixed(4)}          ${med(list.map(x=>x.callRate)).toFixed(4)}`)
  // Take arrays nearest that material's MEDIAN noise, so the calibration is fitted on typical
  // members rather than on the cleanest ones, which would make every floor look optimistic.
  const t = med(list.map(x=>x.hetBafSd))
  picked[m] = list.slice().sort((a,b)=>Math.abs(a.hetBafSd-t)-Math.abs(b.hetBafSd-t))
    .filter(x => x.callRate >= 0.40).slice(0, 4).map(x=>x.name)
}
writeFileSync(SP + 'picked.json', JSON.stringify(picked, null, 2))
console.log('\npicked for injection:')
for (const [m, ns] of Object.entries(picked)) console.log(`  ${m}: ${ns.length}`)

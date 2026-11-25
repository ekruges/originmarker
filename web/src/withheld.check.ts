// Every word a withheld row prints has an entry. Run: node --experimental-strip-types src/withheld.check.ts
// A measured key, a gate, a claim or a cause with no words must fail here, never print its raw name.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { CAUSE_SAYS, CLAIM_SAYS, MEASURE_SAYS, measuredText, withheldRows } from './withheld.ts'
import { PHYSICAL_CAUSES } from './refusals.ts'

const src = readFileSync(new URL('./scoreSample.ts', import.meta.url), 'utf8')
// The keys of every object a refusal records as `measured`, and of the two records spread into them.
const keys = new Set<string>()
for (const m of src.matchAll(/(?:measured: \{|const (?:stageMeasured|dropInMeasured) = \{)([^}]*)\}/g)) {
  for (const k of m[1].matchAll(/(?:^|[,{\s])([a-zA-Z]+):/g)) keys.add(k[1])
}
assert.ok(keys.size >= 10, `found the measured keys: ${[...keys].join(', ')}`)
for (const k of keys) assert.ok(k in MEASURE_SAYS, `measured key with no label: ${k}`)
for (const c of [...PHYSICAL_CAUSES, 'proxy']) assert.ok(c in CAUSE_SAYS, `cause with no words: ${c}`)
for (const m of src.matchAll(/claim: '(\w+)'/g)) assert.ok(m[1] in CLAIM_SAYS, `claim with no words: ${m[1]}`)
assert.throws(() => measuredText({ someNewKey: 1 }), /no words for the measurement someNewKey/)

// The gates, by the rows they make: each literal gate of the scoring and each stage refusal.
const stage = readFileSync(new URL('./stage.ts', import.meta.url), 'utf8')
const refusedBy = /refusedBy\?: ([^\n]+)/.exec(stage)![1].match(/'([\w-]+)'/g)!.map((x) => x.slice(1, -1))
const gates = [...new Set([...src.matchAll(/gate: '([\w.-]+)'/g)].map((m) => m[1]))].filter((g) => !g.startsWith('integrity.'))
for (const g of [...gates, ...refusedBy.map((b) => `stage.${b}`), 'stage.unknown']) {
  assert.doesNotThrow(() => withheldRows({ refusals: [{ claim: 'segment', gate: g, cause: 'proxy', measured: {} }] }, {}), `gate ${g}`)
}

// The words carry no name from the code and no count of changes.
const row = withheldRows({
  refusals: [{ claim: 'sample', gate: 'integrity.disturbed', cause: 'proxy', soft: true,
    measured: { callRate: 0.792, deviantFraction: 0.0455, affectedAutosomes: 1, totalAutosomes: 22 } },
  { claim: 'findings', gate: 'normal-majority', cause: 'no-intensity-reference', measured: { dropIn: NaN, heterozygositySpread: 0.257, normalFraction: 0.682 } }],
  integrity: { why: 'Changes were found and the rest of the genome is undisturbed enough to measure them against.' },
}, {})
const all = row.map((w) => `${w.what} ${w.why} ${w.measured}`).join(' ')
assert.ok(!/\b[a-z]+[A-Z]\w*\b/.test(all.replace(/\bchr[XY]\b/g, '')), all)
assert.ok(all.includes('call rate 79.2%') && all.includes('autosomes in the commonest copy state 68.2%') && all.includes('drop-in not measured'), all)
assert.ok(!/\d+ changes?\b/.test(all), 'a sample cause counts no change')
console.log('withheld.check: ok')

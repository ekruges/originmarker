// Every word a withheld row prints has an entry. Run: node --experimental-strip-types src/withheld.check.ts
// A measured key, a gate, a claim or a cause with no words must fail here, never print its raw name.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { CAUSE_SAYS, CLAIM_SAYS, MEASURE_SAYS, measuredText, sampleCause, withheldRows } from './withheld.ts'
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
const measures = { callRate: 0.792, stage: 'single-cell', affectedAutosomes: 1, totalAutosomes: 22 }
const sampleRefusal = { claim: 'sample' as const, gate: 'integrity.disturbed', cause: 'proxy' as const, soft: true,
  measured: { callRate: 0.792, deviantFraction: 0.0455, affectedAutosomes: 1, totalAutosomes: 22 } }
const finding = { claim: 'findings' as const, gate: 'normal-majority', cause: 'no-intensity-reference' as const,
  measured: { dropIn: NaN, heterozygositySpread: 0.257, normalFraction: 0.682 } }
// A sample whose changes were measured against an undisturbed remainder withholds nothing.
const fine = withheldRows({ refusals: [sampleRefusal, finding],
  integrity: { level: 'disturbed', headline: '1 change on 1 of 22 autosomes', limits: '', measures } }, {})
assert.equal(fine.length, 1, 'no entry for a whole-sample call that counts nothing against the sample')
assert.equal(sampleCause({ level: 'disturbed', headline: 'x', limits: '', measures }), null)
// Where it counts against the sample, the cause says what was found.
const nm = withheldRows({ refusals: [sampleRefusal, finding],
  integrity: { level: 'disturbed', headline: 'NO NORMAL MAJORITY: copy number is withheld', limits: 'Do not read copy number.', measures } }, {})
assert.equal(nm.length, 2)
assert.ok(nm[0].why.startsWith('Copy number is not read: fewer than 80.0% of the autosomes share one copy state'), nm[0].why)
const bad = (m: Partial<typeof measures>) => sampleCause({ level: 'uninterpretable', headline: 'DAMAGED SAMPLE', limits: 'Nothing.', measures: { ...measures, ...m } })!
assert.ok(bad({ stage: 'failed' }).includes('did not resolve to any material'))
assert.ok(bad({ callRate: 0.5 }).includes('call rate on chromosomes 1 to 22 is under the 70.0%'))
assert.ok(bad({ affectedAutosomes: 9 }).includes('9 of 22 autosomes read as changed'))
const all = nm.map((w) => `${w.what} ${w.why} ${w.measured}`).join(' ')
assert.ok(!/\b[a-z]+[A-Z]\w*\b/.test(all.replace(/\bchr[XY]\b/g, '')), all)
assert.ok(all.includes('call rate on chromosomes 1 to 22 79.2%') && all.includes('autosomes in the commonest copy state 68.2%') && all.includes('drop-in not measured'), all)
assert.ok(!/\d+ changes?\b/.test(all), 'a sample cause counts no change')
console.log('withheld.check: ok')

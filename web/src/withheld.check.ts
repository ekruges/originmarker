// Every word a withheld row prints has an entry. Run: node --experimental-strip-types src/withheld.check.ts
// A measured key, a gate, a claim or a cause with no words must fail here, never print its raw name.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { CAUSE_SAYS, CLAIM_SAYS, MEASURE_SAYS, measuredText, sampleCause, withheldRows } from './withheld.ts'
import { PHYSICAL_CAUSES } from './refusals.ts'
import { assessIntegrity } from './integrity.ts'

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
  measured: { autosomalCallRate: 0.792, deviantFraction: 0.0455, affectedAutosomes: 1, totalAutosomes: 22 } }
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
assert.ok(nm[0].why.startsWith('Copy number is not read: no normal level is shown, as fewer than 80.0% of the autosomes share one copy state'), nm[0].why)
const bad = (m: Partial<typeof measures>) => sampleCause({ level: 'uninterpretable', headline: 'DAMAGED SAMPLE', limits: 'Nothing.', measures: { ...measures, ...m } })!
assert.ok(bad({ stage: 'failed' }).includes('did not resolve to any material'))
assert.ok(bad({ callRate: 0.5 }).includes('call rate on chromosomes 1 to 22 is under the 70.0%'))
assert.ok(bad({ affectedAutosomes: 9 }).includes('9 of 22 autosomes read as changed'))

// THE WHOLE-SAMPLE NOTE NAMES ITS PHYSICAL CAUSE. Every branch of the integrity call, under every
// stage refusal: the call carries a physical cause, or `none` where it withholds nothing, and the
// printed row says that cause's words. Only the stand-ins named here may carry `proxy`: a stage
// refused on a test that names no cause (the ceiling, a declared failure, an unread stage) and a
// call rate under the bar, whose uncalled markers carry intensity. A stand-in row says what was
// measured, never that DNA is missing and never the proxy words. A new branch or stage refusal with
// no physical cause fails here.
{
  const ok = { stage: 'single-cell' as const, callRate: 0.95, affectedAutosomes: 0, totalAutosomes: 22, events: 0 }
  const cases: [string, Parameters<typeof assessIntegrity>[0], string][] = [
    ['intact', ok, 'none'],
    ['disturbed', { ...ok, affectedAutosomes: 1, events: 1 }, 'none'],
    ['no normal majority', { ...ok, allelesStand: true, normalFraction: 0.5 }, 'no-intensity-reference'],
    ['drop-in', { ...ok, allelesStand: true, normalFraction: 0.9 }, 'no-intensity-reference'],
    ['stage: no genotype signal', { ...ok, stage: 'failed', stageCause: 'no-genotype-signal' }, 'no-genotype-signal'],
    ['stage: second genome', { ...ok, stage: 'failed', stageCause: 'second-genome' }, 'second-genome'],
    ['affected autosomes', { ...ok, affectedAutosomes: 9, events: 9 }, 'no-intensity-reference'],
    ['stand-in: stage on a test that names no cause', { ...ok, stage: 'failed', stageCause: 'proxy' }, 'proxy'],
    ['stand-in: call rate under the bar', { ...ok, callRate: 0.5 }, 'proxy'],
  ]
  for (const [name, input, cause] of cases) {
    const call = assessIntegrity(input)
    assert.equal(call.cause, cause, `${name}: cause ${call.cause}`)
    const row = sampleCause(call)
    if (cause === 'none') { assert.equal(row, null, `${name}: withholds nothing, so no row`); continue }
    assert.ok(row, `${name}: a row`)
    assert.ok(!row!.includes(CAUSE_SAYS.proxy), `${name}: no proxy words: ${row}`)
    if (cause === 'proxy') {
      assert.ok(!/DNA|template/.test(row! + call.why + call.headline), `${name}: a stand-in claims no missing DNA: ${row}`)
      continue
    }
    const words = CAUSE_SAYS[cause as keyof typeof CAUSE_SAYS]
    assert.ok(row!.includes(words), `${name}: the row names its cause: ${row}`)
    assert.ok(!row!.includes(CAUSE_SAYS.proxy) && !row!.includes('did not resolve'), `${name}: no stand-in words: ${row}`)
    assert.ok(call.why.toLowerCase().includes(words.toLowerCase()) || call.level === 'disturbed' || call.cause === 'no-intensity-reference',
      `${name}: the note says its cause: ${call.why}`)
  }
  // The stage refusals that name no cause are exactly the three of scoreSample's stageCause.
  const src2 = readFileSync(new URL('./scoreSample.ts', import.meta.url), 'utf8')
  const named = /const stageCause[^=]*= ([^\n]+\n[^\n]+)/.exec(src2)![1]
  assert.deepEqual([...named.matchAll(/'([a-z-]+)'/g)].map((x) => x[1]).filter((x) => x !== 'proxy').sort(),
    ['no-genotype-signal', 'second-genome'], 'a stage refusal that names a physical failure is passed as that failure')
}
const all = nm.map((w) => `${w.what} ${w.why} ${w.measured}`).join(' ')
assert.ok(!/\b[a-z]+[A-Z]\w*\b/.test(all.replace(/\bchr[XY]\b/g, '')), all)
assert.ok(all.includes('call rate on chromosomes 1 to 22 79.2%') && all.includes('autosomes in the commonest copy state 68.2%') && !all.includes('drop-in'), all)
assert.ok(!all.includes('not measured'), 'a value that was not measured is not listed')
assert.ok(!/\d+ changes?\b/.test(all), 'a sample cause counts no change')
// A CALL RATE CARRIES ITS SCOPE. Every call-rate label names what it was counted over, and in the
// scoring the rate over all markers (profile.callRate) is recorded only as `callRate` and the rate on
// chromosomes 1 to 22 (the integrity call's) only as `autosomalCallRate`.
for (const [k, v] of Object.entries(MEASURE_SAYS)) {
  if (/callrate/i.test(k)) assert.ok(/over all markers|on chromosomes 1 to 22/.test(v.label), `call rate with no scope: ${k}`)
}
assert.equal(MEASURE_SAYS.callRate.label, 'call rate over all markers')
assert.equal(MEASURE_SAYS.autosomalCallRate.label, 'call rate on chromosomes 1 to 22')
let rates = 0
for (const o of src.matchAll(/(?:measured: \{|const (?:stageMeasured|dropInMeasured) = \{)([^}]*)\}/g)) {
  for (const m of o[1].matchAll(/\b(callRate|autosomalCallRate): ([\w.]+)/g)) {
    rates += 1
    assert.equal(m[1], m[2] === 'profile.callRate' ? 'callRate' : 'autosomalCallRate', `${m[1]} recorded from ${m[2]}`)
  }
}
assert.ok(rates >= 2, 'both call rates are recorded where this check looks')
// ONE FORMAT FOR THE CALL-RATE BAR, on the heading and on the row.
{
  const call = assessIntegrity({ stage: 'single-cell', callRate: 0.5, affectedAutosomes: 0, totalAutosomes: 22, events: 0 })
  const bar = /below the ([\d.]+%)/.exec(call.headline)![1]
  assert.ok(sampleCause(call)!.includes(`under the ${bar}`), `${call.headline} | ${sampleCause(call)}`)
}
assert.ok(!/^triploidy$/.test(CLAIM_SAYS.triploidy), 'a claim is said in words, not by its name')
console.log('withheld.check: ok')

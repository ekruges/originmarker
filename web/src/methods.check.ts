// The methods paragraph says only what the run did. Run: node --experimental-strip-types src/methods.check.ts
import assert from 'node:assert/strict'
import { barsApplied } from './bars.ts'
import { CHECKED, methodsParagraph, type MethodsRun } from './methods.ts'

const run = (o: Partial<MethodsRun> = {}): MethodsRun => ({
  tool: 'OriginMarker 9.9.9 (Check)', chip: 'UK Biobank Axiom', builds: ['GRCh37'], samples: 1, father: 'none',
  mother: false, grouped: 0, sexRows: 0, markers: [], applied: barsApplied({ parents: 0, run: false, comparison: false }), ...o,
})
{
  const m = methodsParagraph(run())
  assert.ok(m.text.includes('No parental array was loaded.') && !/father|mother|parent's/.test(m.text), 'no parent, no parental sentence')
  assert.ok(!/units of one embryo/.test(m.text), 'one cell, no origin sentence')
  assert.ok(!/entered the change list/.test(m.text), 'no sex row, no sex-row sentence')
  assert.deepEqual(m.cites, [CHECKED.peiffer2006])
  assert.ok(m.text.includes('genome build GRCh37'))
}
{
  const m = methodsParagraph(run({ sexRows: 1 }))
  assert.ok(m.text.includes('In one sample that reading entered the change list as a row'), 'a sex row is named where there is one')
  assert.ok(methodsParagraph(run({ samples: 3, sexRows: 2 })).text.includes('In 2 samples that reading'))
}
{
  const applied = barsApplied({ parents: 2, run: false, comparison: false })
  const m = methodsParagraph(run({ father: 'measured', mother: true, applied }))
  assert.ok(m.cites.includes(CHECKED.verdyck2022), 'the two-parent share is cited with both parents only')
  assert.throws(() => methodsParagraph(run({ father: 'measured', applied: barsApplied({ parents: 0, run: false, comparison: false }) })),
    'a sentence cannot quote a bar the run did not apply')
}
{
  const one = methodsParagraph(run({ markers: [{ file: 825656, autosomal: 790123 }] })).text
  assert.ok(one.includes('The file held 825,656 markers, 790,123 of them called on chromosomes 1 to 22.'), one)
  const two = methodsParagraph(run({ samples: 2, markers: [{ file: 825656, autosomal: 600000 }, { file: 825656, autosomal: 790123 }] })).text
  assert.ok(two.includes('The files held 825,656 markers, 600,000 to 790,123 of them called on chromosomes 1 to 22, by file.'), two)
  assert.ok(!/markers, .* called/.test(methodsParagraph(run()).text), 'no marker sentence without the numbers')
}
console.log('methods.check: ok')

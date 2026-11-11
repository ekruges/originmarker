/**
 * The answer key counts only markers whose right answer is known, and drop-in is what it says.
 */
import assert from 'node:assert/strict'
import { emptySameHom, addSameHom, dropIn, parsePanel } from './answerKeys.ts'
import { hetCall, emptyHet } from './obligateHet.ts'

{
  const t = emptySameHom()
  addSameHom(t, 'AA', 'AA', 'AA')   // right
  addSameHom(t, 'AB', 'AA', 'AA')   // drop-in
  addSameHom(t, 'BB', 'AA', 'AA')   // error
  addSameHom(t, 'BB', 'BB', 'BB')   // right
  addSameHom(t, 'AA', 'BB', 'BB')   // error
  addSameHom(t, 'NC', 'AA', 'AA')   // not called: not a reading
  addSameHom(t, 'AB', 'AA', 'BB')   // opposite homozygotes: informative, not an answer key
  addSameHom(t, 'AB', 'AB', 'AB')   // heterozygous parents: anything is possible
  addSameHom(t, 'AA', 'AA', 'NC')   // one parent uncalled: unknown answer
  assert.deepEqual(t, { n: 5, het: 1, wrong: 2 })
  assert.equal(dropIn(t), 1 / 5)
  assert.ok(Number.isNaN(dropIn(emptySameHom())), 'no markers is no measurement, not zero')
}

// THE CASE THIS EXISTS FOR. A one-copy chromosome on a cell with 25% drop-in: of 1,346 informative
// markers 36.5% read heterozygous, 62.1% the surviving parent alone and 1.4% the lost one. Read as
// is, it sits between the boundaries and names no one. With the drop-in removed it is 15.3%
// heterozygous and the survivor 82.8%: one parent's genome, and which one.
{
  const n = 1346
  const t = { ...emptyHet(), informative: n, het: Math.round(0.365 * n),
    onlyA: Math.round(0.014 * n), onlyB: Math.round(0.621 * n) }
  const raw = hetCall(t, 2)
  assert.equal(raw.ploidy, 'uncalled', 'read as is, drop-in leaves it between the boundaries')
  const fixed = hetCall(t, 2, 0.25)
  assert.equal(fixed.ploidy, 'uniparental')
  assert.equal(fixed.survivor, 'b', 'the second parent survives, so the first parent\'s copy is gone')
  assert.ok(Math.abs(fixed.fraction - raw.fraction) < 1e-12, 'what is read is kept as read')
  assert.ok(/after removing this cell's 25\.0% drop-in/.test(fixed.why), `and the reason says so: ${fixed.why}`)
}

// ON A CLEAN CELL THE CORRECTION IS NEGLIGIBLE, so the calls this rule was validated on stand.
{
  const n = 1000
  const t = { ...emptyHet(), informative: n, het: 120, onlyA: 10, onlyB: 870 }
  assert.equal(hetCall(t, 2, 0.002).ploidy, hetCall(t, 2).ploidy)
  assert.equal(hetCall(t, 2, 0.002).survivor, hetCall(t, 2).survivor)
}

// A TWO-COPY CHROMOSOME ON A DROP-IN CELL STAYS TWO COPIES: 70% heterozygous at 25% drop-in is 60%.
{
  const n = 1000
  const t = { ...emptyHet(), informative: n, het: 700, onlyA: 150, onlyB: 150 }
  assert.equal(hetCall(t, 2, 0.25).ploidy, 'biparental')
}

// THE PANEL IS A PARENT PAIR THAT ALWAYS AGREES. Feeding a founder panel's shared homozygote in as
// both parents is the parent-free drop-in, counted by the same tally, so it cannot drift from the
// parent-based reading through a second implementation.
{
  const viaParents = emptySameHom()
  const viaPanel = emptySameHom()
  for (const [s, g] of [['AA', 'AA'], ['AB', 'AA'], ['BB', 'AA'], ['BB', 'BB'], ['NC', 'BB']] as const) {
    addSameHom(viaParents, s, g, g)
    addSameHom(viaPanel, s, g, g)
  }
  assert.deepEqual(viaPanel, viaParents)
  assert.deepEqual(viaPanel, { n: 4, het: 1, wrong: 1 })
}

// THE SHIPPED PANEL'S FORMAT, read by one parser in the browser and on the command line. Extra
// columns are ignored, a line without a recognised allele is skipped, and the header is not a marker.
{
  const p = parsePanel('probeset\tchrom\tallele\nAX-1\t1\tA\nAX-2\t5\tB\nAX-3\t7\tN\n\nAX-4\tB\n')
  assert.deepEqual([...p], [['AX-1', 'AA'], ['AX-2', 'BB'], ['AX-4', 'BB']])
}

console.log('  ok   answerKeys.check.ts')

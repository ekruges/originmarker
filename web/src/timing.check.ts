// Self-check for the timing call. Run: node src/timing.check.ts
//
// Two units of one embryo separate a change carried by the gamete from one that arose after
// fertilisation, and the reciprocal case is the decisive one: a gain in one daughter cell against a
// loss in the other is a division error, whatever the rest of the embryo looks like.
import assert from 'node:assert/strict'
import { callMechanism, eventsOf, timeGroup, type DirectedEvent } from './timing.ts'

const whole = (chrom: string, direction: 'gain' | 'loss'): DirectedEvent =>
  ({ chrom, startBp: 0, endBp: Number.MAX_SAFE_INTEGER, direction })

// --- 1. RECIPROCAL IS A DIVISION ERROR -------------------------------------------------------
{
  const a = [whole('7', 'gain'), whole('8', 'gain')]
  const b = [whole('7', 'loss'), whole('8', 'loss')]
  const call = callMechanism(a[0], [a, b])
  assert.equal(call.mechanism, 'reciprocal')
  assert.ok(call.why.includes('after fertilisation'))
  // Both sides read the same way: the loss is as reciprocal as the gain.
  assert.equal(callMechanism(b[0], [a, b]).mechanism, 'reciprocal')
}

// --- 2. UNIFORM IS THE GAMETE, PARTIAL IS NOT ------------------------------------------------
{
  const a = [whole('21', 'gain')]
  const b = [whole('21', 'gain')]
  const c = [whole('21', 'gain')]
  assert.equal(callMechanism(a[0], [a, b, c]).mechanism, 'meiotic',
    'the same gain in every unit is what an egg or sperm delivers')
  assert.equal(callMechanism(a[0], [a, b, []]).mechanism, 'post-zygotic',
    'a unit without it means it arose in a lineage')
  assert.equal(callMechanism(a[0], [a]).mechanism, 'unresolved',
    'one unit cannot answer this at any quality')
}

// --- 3. A RECIPROCAL PAIR OUTRANKS UNIFORMITY -------------------------------------------------
//
// Three units where two mirror each other and the third carries the gain as well. Counting units
// alone would call that meiotic, which is exactly the wrong answer.
{
  const a = [whole('1', 'gain')]
  const b = [whole('1', 'loss')]
  const c = [whole('1', 'gain')]
  assert.equal(callMechanism(a[0], [a, b, c]).mechanism, 'reciprocal')
}

// --- 4. EVENTS COME OFF A SCORED SAMPLE, WHOLE CHROMOSOMES INCLUDED ---------------------------
{
  const evs = eventsOf({
    chroms: [{ chrom: '5', aneuploidy: 'loss' }, { chrom: '6' }],
    segments: [{ chrom: '9', kind: 'copy-gain', startBp: 1e6, endBp: 9e6 } as never],
    findings: [{ chrom: '7', startBp: 0, endBp: 159e6 }, { chrom: 'genome', startBp: 0, endBp: 0 }],
  })
  assert.deepEqual(evs.map((e) => `${e.chrom}:${e.direction ?? 'none'}`), ['5:loss', '9:gain', '7:none'],
    'whole chromosomes, segments with their direction, findings without one, and never the genome')
}

// --- 5. A WHOLE GROUP AT ONCE, IN THE ORDER GIVEN ---------------------------------------------
{
  const timed = timeGroup([
    { chroms: [{ chrom: '7', aneuploidy: 'gain' }] },
    { chroms: [{ chrom: '7', aneuploidy: 'loss' }] },
  ])
  assert.equal(timed.length, 2)
  assert.equal(timed[0][0].mechanism, 'reciprocal')
  assert.equal(timed[1][0].mechanism, 'reciprocal')
}

console.log('timing.check.ts OK: reciprocal outranks uniformity, and one unit still answers nothing')

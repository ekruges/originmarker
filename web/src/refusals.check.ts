/**
 * A withheld claim is always recorded, always with its cause and its numbers.
 *
 * The census that scores refusals against published truth reads `result.refusals`. A field cleared
 * without a record is a refusal the census cannot see, which is how a gate goes on throwing away
 * confirmed calls unmeasured. So the two lists are pinned to each other at the source, and every
 * record is required to carry a cause from the closed list and a non-empty measurement.
 */
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { PHYSICAL_CAUSES, isPhysical } from './refusals.ts'

const src = readFileSync(new URL('./scoreSample.ts', import.meta.url), 'utf8')

// 1. THE ERASE-ALL BLOCK: every field it clears is a field it records.
{
  const start = src.indexOf('  if (uninterpretable) {\n    log(\'WARN\', \'no chromosomal changes')
  assert.ok(start > 0, 'the stage refusal block must exist where this check looks for it')
  const end = src.indexOf('\n  }\n', start)
  const block = src.slice(start, end)
  const cleared = new Set([...block.matchAll(/c\.([A-Za-z]+) = undefined/g)].map((m) => m[1]))
  const listed = /for \(const claim of \[([^\]]+)\] as const\)/.exec(block)
  assert.ok(listed, 'the block must record its claims from an explicit list')
  const recorded = new Set([...listed![1].matchAll(/'([A-Za-z]+)'/g)].map((m) => m[1]))
  assert.deepEqual([...cleared].sort(), [...recorded].sort(),
    'a field cleared on a refused array must be a field recorded before it is cleared. Cleared: '
    + `${[...cleared].join(', ')}. Recorded: ${[...recorded].join(', ')}.`)
  // Matched as an assignment to a field, since the guard `c[claim] === undefined` also contains
  // the text "= undefined" and sits before the push.
  assert.ok(block.indexOf('refusals!.push') < block.search(/c\.[A-Za-z]+ = undefined/),
    'the records are written BEFORE the fields are cleared, or the withheld value is already gone')
}

// 2. EVERY RECORD CARRIES A CAUSE AND ITS NUMBERS.
{
  const pushes = [...src.matchAll(/refusals!\.push\(\{([\s\S]*?)\}\)/g)].map((m) => m[1])
  assert.ok(pushes.length >= 6,
    `six refusal points are recorded: stage fields, relatedness, copy scan, absence segments, `
    + `findings, integrity. Found ${pushes.length}.`)
  for (const p of pushes) {
    assert.ok(/cause:/.test(p), `a refusal record without a cause: {${p.trim().slice(0, 80)}}`)
    assert.ok(/measured:/.test(p), `a refusal record without its measurement: {${p.trim().slice(0, 80)}}`)
    const cause = /cause: '([a-z-]+)'/.exec(p)?.[1]
    assert.ok(cause === 'proxy' || (cause && isPhysical(cause)),
      `cause '${cause}' is neither physical nor 'proxy'. The closed list: ${PHYSICAL_CAUSES.join(', ')}`)
  }
}

// 3. THE LIST IS CLOSED. Adding a cause is a decision about what counts as a physical failure,
// and it should show up here as one, not slip in through a string.
assert.deepEqual([...PHYSICAL_CAUSES], ['no-template', 'no-genotype-signal', 'second-genome',
  'no-intensity-reference', 'parents-uninformative', 'parent-not-loaded'])

console.log('  ok   refusals.check.ts')

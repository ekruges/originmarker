// The registry of bars against the scoring code. Run: node --experimental-strip-types src/bars.check.ts
//
// The scoring code is every module scoreSample.ts imports, followed import by import. A named
// numeric constant is a top-level `const NAME =` (exported or not), NAME in capitals, whose value
// is a number, arithmetic on numbers, Math functions and names in capitals, or the name of another
// such constant. Each must be in bars.BARS or bars.NOT_BARS, and every entry must still exist.
import assert from 'node:assert/strict'
import { existsSync, readFileSync } from 'node:fs'
import { BARS, NOT_BARS, barValue, barsApplied, barIn } from './bars.ts'

const dir = new URL('./', import.meta.url).pathname
const closure = new Set<string>()
const stack = ['scoreSample.ts']
while (stack.length) {
  const f = stack.pop()!
  if (closure.has(f) || !existsSync(dir + f)) continue
  closure.add(f)
  for (const m of readFileSync(dir + f, 'utf8').matchAll(/from '\.\/([\w.]+?)(?:\.ts)?'/g)) stack.push(`${m[1]}.ts`)
}

const DECL = /^(?:export )?const ([A-Z][A-Z0-9_]*)\s*(?::\s*[A-Za-z.]+\s*)?=\s*(.*?)\s*(?:\/\/.*)?$/
const VALUE = /^(?:[\d_]+(?:\.\d+)?(?:e[-+]?\d+)?|[A-Z][A-Z0-9_]*|Math\.\w+|[-+*/()\s])+$/
const raw: { file: string; name: string; value: string }[] = []
for (const file of closure) {
  for (const line of readFileSync(dir + file, 'utf8').split('\n')) {
    const m = DECL.exec(line)
    if (m && m[2] && VALUE.test(m[2])) raw.push({ file, name: m[1], value: m[2] })
  }
}
// A bare name counts only where it names another constant found here: an alias of a number.
const numeric = new Set(raw.filter((c) => /[\d(*/+-]/.test(c.value) || /Math\./.test(c.value)).map((c) => c.name))
const found = raw.filter((c) => numeric.has(c.name) || numeric.has(c.value))
const key = (c: { file: string; name: string }) => `${c.file}:${c.name}`

// ---------------------------------------------------------------- every named constant is placed
{
  const placed = new Set([...BARS.map(key), ...Object.keys(NOT_BARS)])
  const missing = found.filter((c) => !placed.has(key(c))).map(key)
  assert.deepEqual(missing, [], `named constants of the scoring code in neither BARS nor NOT_BARS: ${missing.join(', ')}`)
  const live = new Set(found.map(key))
  const stale = [...placed].filter((k) => !live.has(k))
  assert.deepEqual(stale, [], `registry entries with no constant in the scoring code: ${stale.join(', ')}`)
  const both = BARS.map(key).filter((k) => k in NOT_BARS)
  assert.deepEqual(both, [], `in both lists: ${both.join(', ')}`)
  assert.equal(new Set(BARS.map(key)).size, BARS.length, 'no bar is listed twice')
  const named = BARS.map((b) => b.name)
  assert.deepEqual(named.filter((n, i) => named.indexOf(n) !== i), [], 'no two bars share a printed name')
  console.log(`${found.length} named constants: ${BARS.length} bars, ${Object.keys(NOT_BARS).length} not bars`)
}

// ---------------------------------------------------------------- no value is written in the registry
{
  const src = readFileSync(dir + 'bars.ts', 'utf8')
  assert.ok(!/value:\s*-?[\d.]/.test(src), 'a value in bars.ts is a literal')
  for (const b of BARS) assert.ok(Number.isFinite(b.value), `${b.name} reads as a number`)
  assert.ok(BARS.every((b) => b.says && !/\d/.test(b.says)), 'each line says what the bar does, with no number in it')
}

// ---------------------------------------------------------------- the same name in two files
{
  const byName = new Map<string, { file: string; value: string }[]>()
  for (const c of found) byName.set(c.name, [...byName.get(c.name) ?? [], c])
  const twice = [...byName].filter(([, v]) => v.length > 1)
  for (const [name, v] of twice) console.log(`  named twice: ${name} ${v.map((x) => `${x.file} ${x.value}`).join(' | ')}`)
}

// ---------------------------------------------------------------- what applies where
{
  const none = barsApplied({ parents: 0, run: false, comparison: false })
  assert.ok(none.every((b) => !b.need), 'with no parent, no run and no comparison only bars every sample uses apply')
  assert.ok(!none.some((b) => b.need?.parent), 'no parental bar with no parent')
  const one = barsApplied({ parents: 1, run: false, comparison: false })
  const two = barsApplied({ parents: 2, run: false, comparison: false })
  assert.ok(one.some((b) => b.need?.parent === 'one') && !one.some((b) => b.need?.parent === 'two'), 'one parent: the one-parent bars, not the two-parent ones')
  assert.ok(two.some((b) => b.need?.parent === 'two') && !two.some((b) => b.need?.parent === 'one'), 'both parents: the two-parent bars, not the one-parent ones')
  assert.ok(!one.some((b) => b.need?.run), 'no run bar with no other array in the run')
  assert.ok(barsApplied({ parents: 1, run: true, comparison: false }).some((b) => b.need?.run && b.need.parent), 'a run with a parent reads the run bars')
  assert.throws(() => barIn(none, 'ABSENCE_MARGIN'), 'a bar not applied cannot be quoted')
  assert.equal(barValue({ value: 3_000_000, unit: 'bp' }), '3 Mb')
  assert.equal(barValue({ value: 3, unit: 'x' }), '3x')
  assert.equal(barValue({ value: 2400, unit: 'markers' }), '2400 markers')
}

// ---------------------------------------------------------------- the report prints the registry, not a typed table
{
  const pdf = readFileSync(dir + 'syngamyPdf.ts', 'utf8')
  assert.ok(!pdf.includes('Constants and where they come from'), 'the typed table of constants is gone from the report')
  assert.ok(/barsApplied\(/.test(pdf), 'the report lists the bars from the registry')
}
console.log('bars.check: ok')

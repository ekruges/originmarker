/**
 * The line reader must be `String.split('\n')`, exactly.
 *
 * It replaced split to stop holding 825,656 line strings at once. Every recorded baseline in this
 * repository is read through loadArray, so a reader that drops the final line, or invents one,
 * changes every number quietly rather than failing. The cases below are the ones that do that: a
 * file with no trailing newline, one with a trailing newline, an empty file, and CRLF, whose stray
 * carriage return split also leaves in place and which must therefore survive unchanged.
 */
import assert from 'node:assert/strict'
import { lineReader } from './joy-common.ts'

const read = (text: string): string[] => {
  const next = lineReader(Buffer.from(text, 'utf8'))
  const out: string[] = []
  for (let l = next(); l !== null; l = next()) out.push(l)
  return out
}

for (const text of [
  'a\nb',           // no trailing newline: the last line must still be returned
  'a\nb\n',         // trailing newline: split yields a final empty field and so must this
  '',               // empty file
  '\n',             // one empty line plus the field after it
  '\n\n\n',         // consecutive blank lines are not collapsed
  'solo',           // a single line and no newline at all
  'a\r\nb\r\n',     // CRLF: the carriage returns stay on the line, as split leaves them
  '#c\n#c\nid\tx\n1\t2\n',
]) {
  assert.deepEqual(read(text), text.split('\n'),
    `lineReader disagrees with split on ${JSON.stringify(text)}`)
}

// A line longer than any plausible buffer chunk, to be sure nothing is chunk-aligned.
const long = 'x'.repeat(200_000)
assert.deepEqual(read(`${long}\n${long}`), [long, long])

// And the byte offsets are real offsets, not character offsets: a multi-byte character before a
// newline would misplace every later line if the reader counted characters.
assert.deepEqual(read('héllo\nwörld\n'), ['héllo', 'wörld', ''])

console.log('  ok   joy-common.check.ts')

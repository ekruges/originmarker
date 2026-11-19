/**
 * No tracked file under audit/truth/ names an array.
 *
 * WHY THIS EXISTS. The results files are written to be read by people who do not have the corpus,
 * and the standard they are held to is in the header of audit/truth/page_says.py: cells are named
 * by number and the key from number to array lives untracked beside the run. On 19 November
 * audit/truth/gains-called-2026-11-18.txt was found naming cells by the last four digits of their
 * accession ("cell 2369"), which a reader turns back into the accession by guessing the prefix,
 * and two of those tails are ambiguous between two arrays of this corpus. The file also carried
 * two different numbering schemes with neither one named. Both were fixed by hand; this is what
 * stops the next one, because the next one will be written by someone who did not read that fix.
 *
 * WHAT COUNTS AS NAMING AN ARRAY, and what does not:
 *   - a GSM accession IS an array, and is refused anywhere in a tracked file here;
 *   - a GSE accession is a SERIES, which is a published dataset rather than a sample, and every
 *     one of these files cites the series its numbers come from. Six tracked files here do it
 *     today and the shipped founder panel is named after one. It is counted and allowed, and the
 *     count is printed so that a reader sees it was considered rather than missed;
 *   - a cell label carrying four or five digits is the accession's tail wearing a number's
 *     clothes: the two schemes in use run to 31 and to 22, so no honest label needs four digits;
 *   - a laboratory barcode (five digits, a dash, two digits) is an array of the laboratory's own
 *     series, which is the one class the plan names outright;
 *   - a unit id of the published series ("series:group:index") is an array of that series under
 *     the answer key's own name for it, which is the same disclosure by another spelling.
 *
 * TRACKED, NOT PRESENT. The list comes from `git ls-files`, so a file that is merely sitting in
 * the directory is not scanned and a file that is committed cannot be hidden from it. If git
 * cannot be run the check FAILS rather than passing on an empty list.
 *
 * Run: node --experimental-strip-types audit/no-array-names.check.ts
 */
import assert from 'node:assert/strict'
import { execFileSync } from 'node:child_process'
import { readFileSync } from 'node:fs'

const DIR = 'audit/truth'

/** Each rule: what it looks for, and what to do instead. */
const REFUSED: { what: string; find: RegExp; instead: string }[] = [
  {
    // SIX DIGITS, NOT THREE. An accession of this corpus is GSM plus seven digits; the shorter
    // form appears once in these files as an illustration of a column prefix ("GSM123.Log R
    // Ratio"), which names nobody. A check that fires on that teaches a reader to skim it, which
    // is the one thing a check must never do.
    what: 'a GSM accession, which names one array',
    find: /GSM\d{6,}/g,
    instead: 'name the cell by its number and keep the key untracked beside the run',
  },
  {
    what: 'a cell label carrying four or five digits, which is an accession\'s tail',
    find: /\bcells? \d{4,5}\b/gi,
    instead: 'the schemes in use run to 31 and to 22; use one of those numbers and name the scheme',
  },
  {
    what: 'a laboratory barcode',
    find: /\b\d{5}-\d{2}\b/g,
    instead: 'name the cell by its number',
  },
  {
    what: 'a published unit id, which is the answer key\'s own name for one array',
    find: /\b[a-z]{3,}:[a-z]?\d+:\d+\b/g,
    instead: 'name the cell by its number',
  },
]

let files: string[]
try {
  files = execFileSync('git', ['ls-files', DIR], { encoding: 'utf8' })
    .split('\n').filter(Boolean)
} catch (e) {
  throw new Error(`this check cannot read the tracked list, so it refuses to pass: ${String(e)}`)
}
assert.ok(files.length > 10, `expected the tracked files of ${DIR}, got ${files.length}`)

const problems: string[] = []
let series = 0
for (const f of files) {
  const text = readFileSync(f, 'utf8')
  series += (text.match(/GSE\d{3,}/g) ?? []).length
  const lines = text.split('\n')
  for (const rule of REFUSED) {
    for (const [i, line] of lines.entries()) {
      for (const hit of line.match(rule.find) ?? []) {
        problems.push(`${f}:${i + 1}  ${hit}  — ${rule.what}. Instead: ${rule.instead}`)
      }
    }
  }
}

assert.deepEqual(problems, [],
  `a tracked file under ${DIR} names an array:\n  ${problems.join('\n  ')}\n`
  + 'The standard is in the header of audit/truth/page_says.py: cells are named by number and the '
  + 'key from number to array lives untracked beside the run.')

console.log(`no-array-names.check: ${files.length} tracked files under ${DIR} name no array `
  + `(${series} series accessions cited, which is allowed: a series is a dataset, not a sample)`)

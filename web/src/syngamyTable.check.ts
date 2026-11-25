// The table a reader joins against their own metadata, and the columns that stop it being
// miscounted.
import assert from 'node:assert/strict'
import { syngamyTable, COLUMNS } from './syngamyTable.ts'
import type { TableSample } from './syngamyTable.ts'

const ch = (over: Partial<import('./syngamyTable.ts').TableChange> = {}) => ({
  chrom: '4', startBp: 1_000_000, endBp: 3_000_000, kind: 'cnn-loh', ...over,
})

// --- 1. INHERITED ROWS ARE ONE DETERMINATION, AND THE TABLE SAYS SO ---------------------------
// The failure this column exists to stop: a reader sums the rows and reports that many independent
// parental determinations, when every one of them rests on a single genome-level call.
{
  const s: TableSample = {
    name: 'A', originClass: 'gynogenetic', zygosity: 'uniparental_homozygous', role: 'paternal',
    changes: Array.from({ length: 12 }, (_, i) =>
      ch({ startBp: i * 5e6, endBp: i * 5e6 + 2e6, origin: 'maternal', band: 'inherited',
        inheritedMargin: 12.9 })),
  }
  const rows = syngamyTable([s]).trim().split('\n')
  assert.equal(rows.length, 13, 'header plus one row per change')
  const col = COLUMNS.indexOf('independent_determination')
  const flags = rows.slice(1).map((r) => r.split('\t')[col])
  assert.equal(flags.filter((f) => f === '1').length, 1,
    `12 inherited rows must sum to ONE determination, got ${flags.filter((f) => f === '1').length}`)
  assert.equal(flags[0], '1', 'and it is carried by the first of them')
  console.log('  12 inherited rows sum to 1 independent determination')
}

// --- 2. THE BASIS IS A COLUMN, NOT A FOOTNOTE -------------------------------------------------
{
  const s: TableSample = {
    name: 'B', role: 'paternal',
    changes: [
      ch({ origin: 'maternal', band: 'A', confidence: 0.991 }),
      ch({ origin: 'maternal', band: 'inherited', inheritedMargin: 8.6 }),
      ch({ origin: 'paternal', band: 'F' }),
      ch({ origin: 'unclear' }),
    ],
  }
  const col = COLUMNS.indexOf('parent_basis')
  const got = syngamyTable([s]).trim().split('\n').slice(1).map((r) => r.split('\t')[col])
  assert.deepEqual(got,
    ['measured-on-interval', 'inherited-from-genome-call', 'direction-only', 'none'],
    'a measured parent, an inherited one and a bare direction are three different claims')
  console.log(`  parent_basis distinguishes ${new Set(got).size} kinds of claim`)
}

// --- 3. SUPPRESSED DIGITS STAY SUPPRESSED IN THE EXPORT ---------------------------------------
// The interface stops printing a number for bands C and D because a low-specificity intermediate
// category is quoted downstream regardless of its label. A table that carried the number anyway
// would be exactly the channel through which it is quoted.
{
  const s: TableSample = {
    name: 'C', role: 'paternal',
    changes: [
      ch({ origin: 'maternal', band: 'A', confidence: 0.991 }),
      ch({ origin: 'maternal', band: 'B', confidence: 0.930 }),
      ch({ origin: 'maternal', band: 'C', confidence: 0.812 }),
      ch({ origin: 'maternal', band: 'D', confidence: 0.618 }),
    ],
  }
  const col = COLUMNS.indexOf('confidence')
  const got = syngamyTable([s]).trim().split('\n').slice(1).map((r) => r.split('\t')[col])
  assert.deepEqual(got, ['0.991', '0.93', '', ''],
    'bands A and B carry their number; C and D carry none, in the export as on screen')
  console.log('  confidence exported for bands A and B only')
}

// --- 4. THE FORMAT SURVIVES THE EVIDENCE TEXT -------------------------------------------------
// Evidence sentences are long and contain punctuation. A stray tab or newline would silently shift
// every column to its right.
{
  const s: TableSample = {
    name: 'D', role: 'paternal',
    changes: [ch({ origin: 'maternal', band: 'inherited',
      why: 'line one\nline two\tand a tab\r\nand a return' })],
  }
  const rows = syngamyTable([s]).trim().split('\n')
  assert.equal(rows.length, 2, 'an embedded newline must not become a new row')
  assert.equal(rows[1].split('\t').length, COLUMNS.length,
    'and an embedded tab must not become a new column')
  console.log(`  ${COLUMNS.length} columns hold even with tabs and newlines in the evidence`)
}

// --- 5. EVERY COLUMN IS FILLED FROM SOMEWHERE --------------------------------------------------
{
  const s: TableSample = {
    name: 'E', headline: 'ANDROGENETIC', originClass: 'androgenetic', zygosity: 'uniparental_homozygous', role: 'paternal',
    stage: 'esc-single', informative: 580_133,
    changes: [ch({ origin: 'paternal', band: 'inherited', inheritedMargin: 2.8, why: 'because', word: 'consistent',
      when: 'after fertilisation', fragile: 'Break at 61.2 Mb in fragile site FRA3B' })],
  }
  const [head, row] = syngamyTable([s]).trim().split('\n')
  assert.equal(head.split('\t').length, COLUMNS.length)
  const cells = row.split('\t')
  const blank = COLUMNS.filter((_, i) => cells[i] === '')
  assert.deepEqual([...blank], ['confidence'],
    `only the suppressed confidence should be blank on a fully populated row, got ${blank}`)
  console.log('  a fully populated row fills every column but the deliberately suppressed one')
}

// --- 6. WHEN IT AROSE, AS THE PAGE SAYS IT, PER STRETCH WHERE THE STRETCHES DIFFER ------------------
{
  const { originText } = await import('./timing.ts')
  const col = COLUMNS.indexOf('when_it_arose')
  const one = originText({ chrom: '17', kind: 'copy-loss', startBp: 6157, endBp: 52.6e6, mechanism: 'consistent-with-gamete' })
  const parts = originText({ chrom: '17', kind: 'copy-loss', locus: 'chr17', startBp: 0, endBp: 0, mechanismStretches: [
    { fromBp: 0, toBp: 52.6e6, mechanism: 'consistent-with-gamete', why: '' },
    { fromBp: 52.6e6, toBp: 78.4e6, mechanism: 'mitotic-reciprocal', why: '' },
    { fromBp: 78.4e6, toBp: Number.MAX_SAFE_INTEGER, mechanism: 'none', why: '' },
  ] })
  const got = syngamyTable([{ name: 'F', changes: [ch({ when: one }), ch({ when: parts }), ch()] }])
    .trim().split('\n').slice(1).map((r) => r.split('\t')[col])
  assert.deepEqual(got, ['consistent with the gamete (two cells read)',
    'chr17 0.0 to 52.6 Mb, part of the whole-chromosome loss: consistent with the gamete (two cells read); '
      + 'chr17 52.6 to 78.4 Mb, part of the whole-chromosome loss: after fertilisation; '
      + 'chr17 78.4 Mb to the end, part of the whole-chromosome loss: not called', ''])
  console.log('  when_it_arose carries the page\'s words, stretch by stretch where they differ')
}

console.log('syngamyTable: joinable, and it cannot be miscounted as independent determinations')

// --- 5. ONE SAMPLE, AS THE CARD READS IT ------------------------------------------------------
// With no parent loaded the result's role is a placeholder, so the loaded-parent column is empty;
// on a sample not reportable every row reads "not reported", as on the card.
{
  const { tableSampleOf } = await import('./syngamyTable.ts')
  const { changeList } = await import('./pageSays.ts')
  const d = { chrom: '16', startBp: 85_629, endBp: 36_600_000, kind: 'copy-loss', locus: 'chr16:85629-36600000',
    origin: 'unclear', why: 'x', step: 'intensity' } as never
  const r = (noParent: boolean) => ({ noParent, role: 'maternal', originClass: 'unclear', zygosity: 'diploid' }) as never
  const says = (reportable: boolean) => ({ setup: 'none', head: { word: 'SINGLE CELL' }, defects: [d],
    list: changeList([d], 'none', reportable) }) as never
  const col = (s: string, c: string) => s.trim().split('\n').slice(1).map((l) => l.split('\t')[COLUMNS.indexOf(c as never)])
  assert.deepEqual(col(syngamyTable([tableSampleOf('A', r(true), says(true))]), 'loaded_parent'), [''])
  assert.deepEqual(col(syngamyTable([tableSampleOf('A', r(false), says(true))]), 'loaded_parent'), ['maternal'])
  assert.deepEqual(col(syngamyTable([tableSampleOf('A', r(true), says(true))]), 'word'), ['position only'])
  assert.deepEqual(col(syngamyTable([tableSampleOf('A', r(true), says(false))]), 'word'), ['not reported'])
  // And its class, parent and evidence state only what was seen.
  const nr = syngamyTable([tableSampleOf('A', r(true), says(false))])
  assert.deepEqual(col(nr, 'class'), ['lower than the rest of this array'])
  assert.deepEqual(col(nr, 'parent'), [''])
  assert.deepEqual(col(nr, 'evidence'), ['Not reported: this sample is not reportable.'])
  const armSays = { setup: 'none', head: { word: 'SINGLE CELL' }, defects: [{ ...(d as object), arm: 'p' }],
    list: changeList([{ ...(d as object), arm: 'p' } as never], 'none', false) } as never
  const armRow = col(syngamyTable([tableSampleOf('A', r(true), armSays)]), 'evidence')[0]
  assert.ok(armRow.includes('the 16p arm: its edge and its extent inside the arm are not placed') && !/\bcall/.test(armRow), armRow)
  console.log('  no parent loaded leaves loaded_parent empty; a sample not reportable reads "not reported"')
}

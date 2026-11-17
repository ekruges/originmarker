// The headline, the words on the changes, and one count, one set of coordinates and one pair of
// numbers per change on every surface. Plan item 3.14, P1 to P4.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import {
  changeList, changesText, COUNTED, figureTitle, headlineOf, kindOf, MEASURED, pairText, pdfCells, rowLine,
  saidOf, whereText, WORD, type HeadInput,
} from './pageSays.ts'
import type { Defect } from './defects.ts'

const root = new URL('../../audit/truth/', import.meta.url)

// ---------------------------------------------------------------- P1, the 31 published cells
{
  const cells = JSON.parse(readFileSync(new URL('page-says-cells.json', root), 'utf8')) as
    (HeadInput & { accession: string; setup: string; expected: string })[]
  assert.equal(cells.length, 93, '31 cells in three setups')
  const words: Record<string, Record<string, number>> = {}
  for (const c of cells) {
    const h = headlineOf(c, '')
    assert.equal(h.word, c.expected, `${c.accession} ${c.setup}: ${h.word}, expected ${c.expected}`)
    if (h.word === 'UNCLEAR') assert.ok(h.reason && h.line.length > 20, `an unclear card says why: ${c.accession}`)
    ;(words[c.setup] ??= {})[h.word] = (words[c.setup][h.word] ?? 0) + 1
  }
  assert.equal(words.two['BOTH PARENTS'], 24)
  assert.equal(words.father['FATHER + SECOND PARENT'], 26)
  assert.equal(words.recon['FATHER + SECOND PARENT'], 28)
  // No cell gets a headline the paper contradicts: all 31 are children of both donors.
  for (const c of cells) assert.ok(!['GYNOGENETIC', 'ANDROGENETIC', 'NEITHER PARENT'].includes(headlineOf(c, '').word), c.accession)
  console.log('  P1: 24 both parents, 26 father + second parent, 28 with the rebuilt father; the rest unclear with a reason')
}

// ---------------------------------------------------------------- the headline rule on its own
{
  const base: HeadInput = { failed: false, zygosity: 'diploid', cls: 'unclear', loaded: [
    { who: 'father', check: 'carries', present: 17, read: 22, absent: false, ratio: 2.3 },
    { who: 'mother', check: 'carries', present: 19, read: 22, absent: false }] }
  assert.equal(headlineOf(base, '').word, 'BOTH PARENTS')
  assert.equal(headlineOf({ ...base, loaded: [base.loaded[0]] }, '').word, 'FATHER + SECOND PARENT')
  assert.equal(headlineOf({ ...base, failed: true }, '').reason, 'failed')
  assert.equal(headlineOf({ ...base, zygosity: 'unknown' }, '').reason, 'copy-number')
  assert.equal(headlineOf({ ...base, loaded: [base.loaded[0], { ...base.loaded[1], check: 'unshown', present: 1 }] }, '').reason, 'not-shown')
  const dis = headlineOf({ ...base, loaded: [{ ...base.loaded[0], absent: true, ratio: 3.4 }] }, '')
  assert.equal(dis.reason, 'disagree')
  assert.ok(dis.line.includes('17 of 22') && dis.line.includes('3.4x'), dis.line)
  // A genome read as one copy keeps the class the ratio gives: a pronucleus is still a pronucleus.
  assert.equal(headlineOf({ ...base, zygosity: 'uniparental_homozygous', cls: 'androgenetic' }, '').word, 'ANDROGENETIC')
  assert.ok(!/—/.test(JSON.stringify([headlineOf(base, ''), dis])), 'no em dashes')
}

// ---------------------------------------------------------------- P2 to P4 on one cell
// The published three-event cell with both parents: 6q in two adjacent pieces found two ways, 16q
// and 18p found two ways, a 3q stretch seen by the newer scan alone, a chr18 imbalance.
const piece = (chrom: string, s: number, e: number, origin: 'paternal' | 'maternal', over: Partial<Defect> = {}): Defect => ({
  chrom, startBp: s, endBp: e, kind: 'parental-absence', locus: `chr${chrom}:${s.toLocaleString('en-US')}-${e.toLocaleString('en-US')}`,
  interval: '+/-100 kb / +/-200 kb', origin, why: 'x', step: 'run',
  second: { verdict: 'parent-lost', parent: origin }, rate: 0.105, nullRate: 0.023, ...over,
})
const cellA: Defect[] = [
  { chrom: '18', startBp: 0, endBp: 78_009_148, kind: 'whole-chromosome', locus: 'chr18', origin: 'unclear', why: 'x', step: 'imbalance' },
  piece('6', 107_034_701, 163_755_870, 'paternal'),
  piece('16', 73_603_649, 90_146_288, 'maternal', { rate: 0.114, nullRate: 0.017 }),
  piece('18', 65_863, 14_717_921, 'maternal', { rate: 0.158, nullRate: 0.017 }),
  piece('6', 71_075_343, 107_033_570, 'paternal', { rate: 0.106 }),
  piece('3', 146_367_804, 154_890_206, 'maternal', { step: 'allele-scan', second: undefined, rate: 0.063, nullRate: 0.019, lrr: -0.181 }),
]
{
  const l = changeList(cellA, 'two')
  assert.equal(l.counted.length, 3, 'one loss printed as two pieces is one change')
  assert.deepEqual(l.counted.map((g) => g.word), ['found two ways', 'found two ways', 'found two ways'])
  assert.deepEqual(l.once.map((g) => g.chrom), ['3'])
  assert.deepEqual(l.seen.map((g) => g.chrom), ['18'])
  assert.equal(changesText(l), '3 changes found; 1 more stretch seen once, not confirmed.')
  const head = headlineOf({ failed: false, zygosity: 'diploid', cls: 'unclear', loaded: [
    { who: 'father', check: 'carries', present: 17, read: 22, absent: false },
    { who: 'mother', check: 'carries', present: 19, read: 22, absent: false }] }, changesText(l))
  assert.equal(head.line, "Both parents' genomes are in this cell. 3 changes found; 1 more stretch seen once, not confirmed.")
  assert.equal(figureTitle(head, l, { failed: false, zygosity: 'diploid', cls: 'x', loaded: [] }),
    "Both parents' copies are present on every chromosome, except three stretches.")
  // One set of coordinates: the card line, the PDF cells and the figure all print the pieces' own.
  const six = l.counted.find((g) => g.chrom === '6')!
  assert.equal(whereText(six), 'chr6 71.1-107.0, 107.0-163.8 Mb')
  assert.ok(rowLine(six).startsWith(whereText(six)))
  const [where, , word, detail] = pdfCells(six, 'two')
  assert.equal(where, whereText(six))
  assert.equal(word, six.word)
  for (const p of six.pieces) {
    assert.ok(detail.includes(p.locus), `the PDF prints the refined locus: ${p.locus}`)
    // One pair of numbers: the card's pair, the PDF's and the figure line's are the same digits.
    const pair = pairText(p).match(/[\d.]+%/g)!
    assert.ok(pair.every((x) => detail.includes(x)), `the PDF quotes the card's pair ${pair}`)
    assert.ok(pair.every((x) => saidOf(p)!.includes(x)), `the figure quotes the card's pair ${pair}`)
  }
  assert.ok(detail.includes('19 of 26'), 'the measured count is in the PDF')
  // P3: no named, printed call carries a band word.
  for (const g of [...l.counted, ...l.once]) assert.ok(!/no call|weak|no parent|\b[A-F]\b/.test(g.word), g.word)
}

// ---------------------------------------------------------------- P3: the table is the results file
{
  const text = readFileSync(new URL('page-says-2026-11-17.txt', root), 'utf8')
  const lines = [...text.matchAll(/^ {2}MEASURED (\S+) (\S+) (\d+) (\d+)$/gm)]
  assert.ok(lines.length >= 15, 'the results file carries the measured table')
  for (const [, kind, setup, right, scored] of lines) {
    const got = (MEASURED as Record<string, Record<string, readonly [number, number]>>)[kind]?.[setup]
    assert.deepEqual(got, [Number(right), Number(scored)], `${kind} ${setup}`)
  }
  for (const k of Object.keys(WORD) as (keyof typeof WORD)[]) {
    assert.ok(['found two ways', 'consistent', 'one reading', 'also seen', 'position only'].includes(WORD[k]))
  }
  assert.ok(COUNTED.has('found two ways') && COUNTED.has('consistent') && !COUNTED.has('one reading'))
  assert.equal(kindOf({ ...cellA[1], second: { verdict: 'refused', parent: null } }, 'two'), 'run-only')
  console.log(`  P2 to P4: one change per chromosome and parent, one coordinate set, one pair; ${lines.length} measured counts match the results file`)
}

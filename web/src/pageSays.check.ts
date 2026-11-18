// The headline, the words on the changes, and one count, one set of coordinates and one pair of
// numbers per change on every surface. Plan item 3.14, P1 to P4.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import {
  changeList, changesText, COUNTED, defectsForResult, disturbedText, figureTitle, headlineOf, kindOf, MEASURED, pairText, pdfCells, rowLine,
  saidOf, spanMb, whereText, whoseText, WORD, type HeadInput,
} from './pageSays.ts'
import type { Defect } from './defects.ts'

const root = new URL('../../audit/truth/', import.meta.url)

// ---------------------------------------------------------------- P1, the 31 published cells
{
  const cells = JSON.parse(readFileSync(new URL('page-says-cells.json', root), 'utf8')) as
    (HeadInput & { cell: string; setup: string; expected: string })[]
  assert.equal(cells.length, 93, '31 cells in three setups')
  const words: Record<string, Record<string, number>> = {}
  for (const c of cells) {
    const h = headlineOf(c, '')
    assert.equal(h.word, c.expected, `${c.cell} ${c.setup}: ${h.word}, expected ${c.expected}`)
    if (h.word === 'UNCLEAR') assert.ok(h.reason && h.line.length > 20, `an unclear card says why: ${c.cell}`)
    ;(words[c.setup] ??= {})[h.word] = (words[c.setup][h.word] ?? 0) + 1
  }
  assert.equal(words.two['BOTH PARENTS'], 24)
  assert.equal(words.father['FATHER + SECOND PARENT'], 26)
  assert.equal(words.recon['FATHER + SECOND PARENT'], 28)
  // No cell gets a headline the paper contradicts: all 31 are children of both donors.
  for (const c of cells) assert.ok(!['GYNOGENETIC', 'ANDROGENETIC', 'NEITHER PARENT'].includes(headlineOf(c, '').word), c.cell)
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

// ---------------------------------------------------------------- a disturbed genome is a sentence
// The complex finding (abnormalities.detectComplex) is a property of the whole genome: no row, not
// counted, said once.
{
  const complex: Defect = {
    chrom: 'genome', startBp: 0, endBp: 0, kind: 'complex', locus: 'whole genome', origin: 'unclear', step: 'finding',
    why: '9 of 22 autosomes deviate from this array\'s own centre, a share of 0.41 over the 0.3 at which self-reference fails',
  }
  for (const setup of ['two', 'measured', 'none'] as const) {
    const l = changeList([...cellA, complex], setup)
    assert.ok(![...l.counted, ...l.once, ...l.seen].some((g) => g.chrom === 'genome'), `no genome row (${setup})`)
    assert.equal(l.counted.length + l.once.length, changeList(cellA, setup).counted.length + changeList(cellA, setup).once.length)
    assert.equal(disturbedText(l), '9 of 22 chromosomes carry a change: this cell\'s genome is widely disturbed.')
    assert.ok(!/letter|\u2014/.test(disturbedText(l)))
  }
  assert.equal(disturbedText(changeList(cellA, 'two')), '')
}

// ---------------------------------------------------------------- a second copy lost on a lost chromosome
// partLost.ts prints a copy-loss piece beside a whole-chromosome loss only a further copy down, so
// the piece reads as a second copy lost, with or without a parent; an allele-scan stretch does not.
{
  const wholeLoss = (origin: 'paternal' | 'unclear'): Defect => ({
    chrom: '6', startBp: 0, endBp: 170_000_000, kind: 'copy-loss', locus: 'chr6', origin, why: 'x',
    step: origin === 'unclear' ? 'whole-intensity' : 'whole-alleles',
  })
  const copyPiece: Defect = {
    chrom: '6', startBp: 120_000_000, endBp: 170_400_000, kind: 'copy-loss', locus: 'chr6:120,000,000-170,400,000',
    origin: 'unclear', why: 'x', step: 'intensity', rate: 0.4, nullRate: 0.1,
  }
  for (const [setup, origin] of [['none', 'unclear'], ['measured', 'paternal']] as const) {
    const l = changeList([wholeLoss(origin), copyPiece], setup)
    const all = [...l.counted, ...l.once]
    assert.equal(whoseText(all.find((g) => !g.whole)!), 'second copy lost too', setup)
    assert.notEqual(whoseText(all.find((g) => g.whole)!), 'second copy lost too')
  }
  const alleleStretch: Defect = { ...copyPiece, kind: 'parental-absence', step: 'allele-scan', origin: 'maternal' }
  const l = changeList([wholeLoss('paternal'), alleleStretch], 'measured')
  assert.equal(whoseText([...l.counted, ...l.once].find((g) => !g.whole)!), "mother's copy missing")
  // Without a whole-chromosome loss the piece keeps its own wording.
  assert.equal(whoseText(changeList([copyPiece], 'none').counted[0]), 'copy lost, parent not named')
}

// ---------------------------------------------------------------- one stretch, two scans
// A loss from the copy scan and from an absence scan over the same stretch is one row with both
// readings; a different stretch on the same chromosome stays its own row.
{
  const at = (s: number, e: number, over: Partial<Defect>): Defect => ({
    chrom: '15', startBp: s, endBp: e, kind: 'copy-loss', locus: `chr15:${s}-${e}`, origin: 'unclear', why: 'x', ...over,
  })
  const copy15 = at(71_500_000, 100_000_000, { step: 'intensity', rate: 0.3, nullRate: 0.1 })
  const abs15 = at(76_600_000, 98_900_000, { kind: 'parental-absence', step: 'window', rate: 0.2, nullRate: 0.02 })
  const copy11 = { ...at(18_600_000, 26_300_000, { step: 'intensity' }), chrom: '11' }
  const abs11 = { ...at(17_400_000, 130_300_000, { kind: 'parental-absence', step: 'window' }), chrom: '11' }
  const l = changeList([copy15, abs15, copy11, abs11], 'measured')
  const fifteen = l.once.filter((g) => g.chrom === '15')
  assert.equal(fifteen.length, 1, 'one stretch read by two scans is one row')
  assert.equal(fifteen[0].pieces.length, 2, 'both readings kept')
  assert.equal(fifteen[0].word, 'one reading')
  assert.equal(whoseText(fifteen[0]), 'copy lost, parent not named')
  assert.equal(spanMb(fifteen[0]), '28.5', 'overlapping readings counted once')
  assert.equal(whereText(fifteen[0]), 'chr15 71.5-100.0, 76.6-98.9 Mb')
  assert.equal(l.once.filter((g) => g.chrom === '11').length, 2, 'different stretches stay two rows')
  assert.ok(changesText(l).includes('3 more stretches'), changesText(l))
  const named = changeList([copy15, { ...abs15, origin: 'paternal', step: 'allele-scan' }], 'measured')
  assert.equal(named.once.length, 1)
  assert.equal(whoseText(named.once[0]), "father's copy missing", 'the row names the parent the call names')
  assert.ok(!named.counted.length, 'never found two ways')
}

// ---------------------------------------------------------------- one parent, no two-parent count
// The opposite-homozygote count needs both parents' arrays; a one-parent row quotes the copy-state
// score instead, with the father measured, the father rebuilt or the mother alone.
{
  const chrom = (c: string, over: object = {}) => ({
    chrom: c, aneuploidy: undefined, lostParent: undefined, callFraction: 1, lrrShift: 0,
    copyState: { state: 'disomic', base: 0.2, r: 0.45, n1: 500, n2: 500 }, ...over,
  })
  const result = (role: 'paternal' | 'maternal', twoParents: boolean) => ({
    role, twoParents, segments: [], gains: [], losses: [], oneParent: [], dosageCalls: [], findings: [],
    chroms: [chrom('6', {
      aneuploidy: 'loss', lostParent: role, callFraction: 0.75, lrrShift: -0.87,
      contribution: { informative: 26930, fraction: 0.314, onlyA: 0, onlyB: 0 },
      copyState: { state: 'parent-lost', base: 0.26, r: 0.97, n1: 4376, n2: 3238 },
    }), ...['1', '2', '3'].map((c) => chrom(c))],
  }) as never
  for (const role of ['paternal', 'maternal'] as const) {
    const one = defectsForResult(result(role, false)).find((d) => d.chrom === '6')!
    assert.ok(!one.why.includes('opposite homozygotes'), `a one-parent row carries a two-parent count: ${one.why}`)
    const who = role === 'paternal' ? 'father' : 'mother'
    assert.ok(one.why.includes(`Read with the ${who} alone`) && one.why.includes('0.97') && one.why.includes('0.45'), one.why)
    const two = defectsForResult(result(role, true)).find((d) => d.chrom === '6')!
    assert.ok(two.why.includes('Of 26,930 markers where the parents are opposite homozygotes'), two.why)
  }
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

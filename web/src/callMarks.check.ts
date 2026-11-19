// What the figure marks against what the card counts: one number a change, nothing marked twice,
// nothing marked that the card does not list, and the list saying which readings do not count.
import assert from 'node:assert'
import { readFileSync } from 'node:fs'
import { changeList, defectsForResult, type ChangeList } from './pageSays.ts'
import { directionColour, isCounted } from './plotSeries.ts'
import { direction, type Defect } from './defects.ts'

/** The figure's rule, as CellFigure applies it: every mark numbered in chromosome order. */
function numbering(marks: readonly { chrom: string; startBp: number }[]) {
  const ordered = [...marks].sort((a, b) => +a.chrom - +b.chrom || a.startBp - b.startBp)
  const numbers = new Map<string, number>()
  ordered.forEach((m, i) => numbers.set(`${m.chrom}:${m.startBp}`, i + 1))
  return { ordered, numbers }
}

// ---------------------------------------------------------------- one number a mark, in order
{
  const marks = [{ chrom: '11', startBp: 0 }, { chrom: '2', startBp: 5e6 }, { chrom: '2', startBp: 0 }]
  const { ordered, numbers } = numbering(marks)
  assert.deepEqual(ordered.map((m) => `${m.chrom}:${m.startBp}`), ['2:0', '2:5000000', '11:0'],
    'chromosome order, then position, so 2 comes before 11')
  assert.deepEqual([...numbers.values()], [1, 2, 3], 'every mark has its own number')
  assert.equal(new Set(numbers.values()).size, marks.length, 'no number is used twice')
  assert.equal(numbers.get('2:5000000'), 2, 'a chromosome carrying two changes gets two numbers')
}

// ---------------------------------------------------------------- the card's count against the marks
/** A cell as the scoring pass leaves it: three changes the rule counts, one read once, one not counted. */
const defect = (chrom: string, kind: Defect['kind'], step: string, extra: Partial<Defect> = {}): Defect => ({
  chrom, kind, locus: `chr${chrom}`, startBp: 0, endBp: 1e8, why: '', origin: 'paternal',
  step, ...extra,
} as Defect)

{
  const defects: Defect[] = [
    defect('2', 'copy-loss', 'whole'),
    defect('5', 'copy-loss', 'whole'),
    defect('9', 'parental-absence', 'second-test'),
    defect('4', 'copy-loss', 'intensity'),            // one reading
    defect('16', 'imbalance' as Defect['kind'], 'imbalance'),
  ]
  const list: ChangeList = changeList(defects, 'two')
  const all = [...list.counted, ...list.once, ...list.seen]
  const marks = all.map((g) => ({ chrom: g.chrom, startBp: g.startBp }))
  const { numbers } = numbering(marks)
  assert.equal(numbers.size, all.length, 'every row the card lists is marked exactly once')
  assert.ok(list.counted.length <= all.length, 'the headline counts a subset of what is marked')
  // The figure must be able to say, on each line, which readings the headline does not count.
  for (const g of list.once) assert.equal(g.word, 'one reading')
  for (const g of list.seen) assert.equal(g.word, 'also seen')
  // No chromosome is marked by two rows of the same kind: a change stated twice is a fault.
  const keys = marks.map((m) => `${m.chrom}:${m.startBp}`)
  assert.equal(new Set(keys).size, keys.length, 'no change is marked twice')
}

// ---------------------------------------------------------------- nothing marked that is not listed
{
  const defects: Defect[] = [defect('7', 'copy-loss', 'whole'), defect('7', 'copy-loss', 'whole')]
  const list = changeList(defects, 'two')
  const rows = [...list.counted, ...list.once, ...list.seen]
  assert.equal(rows.length, 1, 'two readings of one change are one row, so the figure marks it once')
  assert.equal(rows[0].pieces.length, 2, 'both readings stay in the row')
}

// ---------------------------------------------------------------- every row lands on exactly one list
// The figure marks what the card lists, so a row that fell off every list, or onto two, would make
// the number of marks disagree with the headline for a reason no reader could see.
{
  const defects: Defect[] = [
    defect('2', 'copy-loss', 'whole'), defect('4', 'copy-loss', 'intensity'),
    defect('6', 'parental-absence', 'second-test'), defect('9', 'copy-gain', 'whole'),
    defect('13', 'copy-loss', 'run'), defect('21', 'parental-absence', 'allele-scan'),
  ]
  const l = changeList(defects, 'two')
  const all = [...l.counted, ...l.once, ...l.seen]
  const keys = all.map((g) => `${g.chrom}:${g.startBp}:${g.call}`)
  assert.equal(new Set(keys).size, keys.length, 'no row is on two lists')
  const pieces = all.flatMap((g) => g.pieces)
  assert.equal(new Set(pieces).size, defects.length, 'every change reaches exactly one row')
  assert.ok(l.counted.length < all.length, 'this cell has readings the headline does not count')
}

// ---------------------------------------------------------------- colour is the direction, always
// One channel, one meaning: colour says which way the copies went, weight says whether the card
// counts the reading. A gain the card does not count is still drawn as a gain, because that is the
// reading a reader most needs to understand; its badge is an outline and its chip says why.
{
  const css = readFileSync(new URL('./index.css', import.meta.url), 'utf8')
  const token = (name: string) => (css.match(new RegExp(`--${name}:\\s*([^;]+);`)) ?? [])[1]?.trim()
  const loss = token('om-defect')
  const gain = token('om-gain')
  assert.ok(loss && gain && loss !== gain, 'a loss and a gain have their own colours')
  const hue = (hex: string) => {
    const [r, g, b] = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255)
    const max = Math.max(r, g, b); const min = Math.min(r, g, b)
    if (max === min) return 0
    const d = max - min
    const h = max === r ? ((g - b) / d) % 6 : max === g ? (b - r) / d + 2 : (r - g) / d + 4
    return (h * 60 + 360) % 360
  }
  const light = (hex: string) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16))
    .reduce((s, v, i) => s + v * [0.2126, 0.7152, 0.0722][i], 0)
  const apart = Math.abs(hue(loss!) - hue(gain!))
  assert.ok(Math.min(apart, 360 - apart) > 90, `a loss and a gain are far apart in hue (${apart.toFixed(0)} degrees)`)
  assert.ok(Math.abs(light(loss!) - light(gain!)) > 10,
    'and differ in lightness too, so they survive a reader who sees no colour')
  // The gain colour is the new one, so it must clear both colours the parent tracks use. The loss
  // colour is the page's own defect red, which shares a hue family with the mother's brown: it is
  // told apart by lightness instead, and it never appears inside a panel where a track is drawn.
  for (const track of ['om-blue', 'om-higher']) {
    const d = Math.abs(hue(token(track)!) - hue(gain!))
    assert.ok(Math.min(d, 360 - d) > 25, `the gain colour stays clear of --${track}`)
  }
  assert.ok(Math.abs(light(loss!) - light(token('om-higher')!)) > 25,
    "the loss colour is told from the mother's track by lightness")
  const blue = Math.abs(hue(token('om-blue')!) - hue(loss!))
  assert.ok(Math.min(blue, 360 - blue) > 25, "the loss colour stays clear of the father's track")
}

// ---------------------------------------------------------------- and weight carries the counting
{
  assert.notEqual(directionColour('gain'), directionColour('loss'), 'the two directions never share a colour')
  // A class that is neither a loss nor a gain is neither colour: an imbalance the card cannot call
  // either way must not be painted as a copy gone.
  assert.equal(direction(false, 'whole-chromosome'), 'none', 'an imbalance has no direction')
  assert.equal(direction(true, 'trisomy'), 'gain', 'three copies is a gain')
  assert.equal(direction(false, 'copy-loss'), 'loss', 'a lost copy is a loss')
  assert.ok(![directionColour('gain'), directionColour('loss')].includes(directionColour('none')),
    'and it is drawn in neither direction colour')
  // The colour rule reads the direction and nothing else, so no status can grey a direction out.
  // That matters most on the readings the card does not count: those are the lines where a reader
  // needs to see whether a copy was doubled or lost.
  const marks: [{ status?: string; withheld?: string }, boolean][] = [
    [{}, true], [{ status: 'counted' }, true], [{ status: 'one reading' }, false],
    [{ status: 'also seen' }, false], [{ withheld: 'too little evidence' }, false],
  ]
  for (const [m, counts] of marks) {
    assert.equal(directionColour('gain'), 'var(--om-gain)', 'an extra copy is teal')
    assert.equal(directionColour('loss'), 'var(--om-defect)', 'and a copy gone is red')
    assert.equal(isCounted(m), counts,
      'and whether the card counts it is a separate question, answered on its own')
  }
}


console.log('callMarks: ok')

// Self-check for rule B1, when a change arose across the units of one embryo. Run: node src/timing.check.ts
import assert from 'node:assert/strict'
import { callMechanism, ORIGIN_WORD, readingOf, timeGroup, type UnitChange, type UnitReading } from './timing.ts'

const ch = (chrom: string, kind: string, origin = 'unclear', startBp = 0, endBp = 100e6, whole = true): UnitChange =>
  ({ chrom, startBp, endBp, kind, origin, whole })
const unit = (changes: UnitChange[], o: Partial<UnitReading> = {}): UnitReading =>
  ({ changes, measured: ['1', '7', '13', '21'], reportable: true, ...o })
const call = (units: UnitReading[], u = 0, c = 0) => callMechanism(u, units[u].changes[c], units)
const normal = unit([])

// --- 1. ONE UNIT: NO ORIGIN, NO WORD ---------------------------------------------------------
{
  const m = call([unit([ch('21', 'copy-gain')])])
  assert.equal(m.mechanism, 'none')
  assert.equal(ORIGIN_WORD[m.mechanism], '')
}

// --- 2. EVERY UNIT CARRIES IT: THREE OR MORE FROM THE GAMETE, TWO ONLY CONSISTENT -------------
{
  const g = () => unit([ch('21', 'copy-gain')])
  assert.equal(call([g(), g(), g()]).mechanism, 'meiotic')
  assert.equal(ORIGIN_WORD.meiotic, 'from the gamete')
  const two = call([g(), g()])
  assert.equal(two.mechanism, 'consistent-with-gamete', 'two cells can share an error of the first division')
  assert.equal(ORIGIN_WORD[two.mechanism], 'consistent with the gamete (two cells read)')
  assert.ok(two.why.startsWith('2 of 2 units carry this gain, 0 read chromosome 21 normal'), two.why)
  // Kinds of one direction count together: a whole copy-loss and a parental absence are both losses.
  assert.equal(call([unit([ch('13', 'copy-loss')]), unit([ch('13', 'parental-absence')]),
    unit([ch('13', 'monosomy')])]).mechanism, 'meiotic')
}

// --- 3. THE DIRECTION MUST MATCH -------------------------------------------------------------
{
  // A loss beside a direction-less finding on the same chromosome is not carried by both.
  const units = [unit([ch('7', 'copy-loss')]), unit([ch('7', 'cnn-loh')]), unit([ch('7', 'copy-loss')])]
  assert.equal(call(units).mechanism, 'none')
  assert.ok(call(units).why.startsWith('2 of 3 units'))
  // A finding carries no origin of its own.
  assert.equal(call(units, 1).mechanism, 'none')
}

// --- 4. THE PARENT MUST MATCH WHERE BOTH NAME ONE ----------------------------------------------
{
  const m = (o: string) => unit([ch('13', 'copy-loss', o)])
  assert.equal(call([m('maternal'), m('paternal'), m('maternal')]).mechanism, 'none',
    'a mother\'s copy lost in one cell and a father\'s in another are two events')
  assert.equal(call([m('maternal'), m('unclear'), m('maternal')]).mechanism, 'meiotic',
    'a cell that names no parent does not contradict one that does')
  // Different parents with another cell reading normal: not in every cell, after fertilisation.
  assert.equal(call([m('maternal'), m('paternal'), normal]).mechanism, 'mitotic-not-all')
}

// --- 5. RECIPROCAL: ANOTHER UNIT GAINS WHAT THIS ONE LOSES -------------------------------------
{
  const units = [unit([ch('1', 'copy-gain')]), unit([ch('1', 'copy-loss')]), unit([ch('1', 'copy-gain')])]
  assert.equal(call(units).mechanism, 'mitotic-reciprocal', 'reciprocal outranks two carriers of three')
  assert.equal(call(units, 1).mechanism, 'mitotic-reciprocal')
  assert.equal(ORIGIN_WORD['mitotic-reciprocal'], 'after fertilisation')
  assert.ok(call(units).why.includes('opposite change'))
  // Over an overlapping interval only.
  const apart = [unit([ch('1', 'copy-gain', 'unclear', 1e6, 10e6, false)]),
    unit([ch('1', 'copy-loss', 'unclear', 50e6, 60e6, false)])]
  assert.equal(call(apart).mechanism, 'none')
  // A gain and a loss in ONE cell are not a pair of daughter cells.
  const own = [unit([ch('1', 'copy-gain'), ch('1', 'parental-absence')]), unit([ch('1', 'copy-gain')])]
  assert.equal(call(own).mechanism, 'consistent-with-gamete')
}

// --- 6. "ABSENT" NEEDS A REPORTABLE CELL THAT READS THE CHROMOSOME NORMAL ----------------------
{
  const g = unit([ch('21', 'copy-gain')])
  assert.equal(call([g, normal, normal]).mechanism, 'mitotic-not-all')
  assert.ok(call([g, normal, normal]).why.startsWith('1 of 3 units carry this gain, 2 read chromosome 21 normal'))
  assert.equal(call([g, unit([], { reportable: false })]).mechanism, 'none', 'an unreportable cell reads nothing normal')
  assert.equal(call([g, unit([], { measured: ['1'] })]).mechanism, 'none', 'an unmeasured chromosome is not normal')
  assert.equal(call([g, unit([ch('21', 'cnn-loh', 'unclear', 5e6, 9e6, false)])]).mechanism, 'none',
    'a piece of any kind on it is not normal')
}

// --- 7. PIECES OVERLAP; A WHOLE ROW SPANS ITS CHROMOSOME ---------------------------------------
{
  const p = (a: number, b: number) => ch('7', 'copy-loss', 'unclear', a, b, false)
  assert.equal(call([unit([p(5e6, 25e6)]), unit([p(4.8e6, 25.2e6)]), unit([p(6e6, 30e6)])]).mechanism, 'meiotic')
  assert.equal(call([unit([p(5e6, 25e6)]), unit([ch('7', 'copy-loss', 'unclear', 0, 0, true)])]).mechanism,
    'consistent-with-gamete', 'a whole row with no recorded extent still covers the piece')
}

// --- 8. A SCORED SAMPLE'S READING, AND A WHOLE GROUP IN ORDER ----------------------------------
{
  const r = readingOf([
    { chrom: '5', startBp: 0, endBp: 180e6, kind: 'copy-loss', origin: 'maternal', locus: 'chr5' },
    { chrom: '9', startBp: 1e6, endBp: 9e6, kind: 'copy-gain', origin: 'unclear', locus: 'chr9:1-9' },
    { chrom: 'genome', startBp: 0, endBp: 0, kind: 'complex', origin: 'unclear', locus: 'genome' },
  ], { chroms: [{ chrom: '5' }, { chrom: '9' }], stage: { stage: 'blastomere' }, integrity: { level: 'intact' } })
  assert.deepEqual(r.changes.map((c) => `${c.chrom}:${c.whole}`), ['5:true', '9:false'])
  assert.equal(r.reportable, true)
  assert.equal(readingOf([], { stage: { stage: 'failed' } }).reportable, false)
  assert.equal(readingOf([], { integrity: { level: 'disturbed' } }).reportable, false)
  const timed = timeGroup([unit([ch('7', 'copy-gain')]), unit([ch('7', 'copy-loss')])])
  assert.equal(timed.length, 2)
  assert.equal(timed[0][0].mechanism, 'mitotic-reciprocal')
  assert.equal(timed[1][0].kind, 'copy-loss')
}

console.log('timing.check.ts OK: rule B1 with direction, parent and measured normal; one cell prints no origin')

// --- THE REPORT'S TIMING TABLE PRINTS NO COORDINATES FOR A WHOLE ROW -----------------------------
// A whole row's end is the chromosome's length only where a parent was loaded; with none it is 0.
{
  const { timingWhere } = await import('./syngamyPdf.ts')
  const [[w, s]] = timeGroup([unit([ch('16', 'monosomy', 'unclear', 0, 0), ch('7', 'copy-loss', 'unclear', 5e6, 9e6, false)]), normal])
  assert.equal(timingWhere(w), 'chr16, whole chromosome')
  assert.equal(timingWhere(s), 'chr7 5.0-9.0Mb')
  console.log('timing.check.ts: the report prints a whole row as the whole chromosome')
}

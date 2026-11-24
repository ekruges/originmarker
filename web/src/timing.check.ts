// Self-check for rule B1, when a change arose across the units of one embryo. Run: node src/timing.check.ts
import assert from 'node:assert/strict'
import {
  callMechanism, ORIGIN_WORD, readingOf, SAME_BREAK_BP, split, stretchesOf, stretchText, stretchWhere, stretchWord, timeGroup, TO_END,
  type UnitChange, type UnitReading,
} from './timing.ts'

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
  // Over the same stretch only: apart, each cell reads the other's stretch normal.
  const apart = [unit([ch('1', 'copy-gain', 'unclear', 1e6, 10e6, false)]),
    unit([ch('1', 'copy-loss', 'unclear', 50e6, 60e6, false)])]
  assert.equal(call(apart).mechanism, 'mitotic-not-all')
  assert.equal(call(apart, 1).mechanism, 'mitotic-not-all')
  // A gain and a loss in ONE cell are not a pair of daughter cells.
  const own = [unit([ch('1', 'copy-gain'), ch('1', 'parental-absence')]), unit([], { reportable: false })]
  assert.equal(call(own).mechanism, 'none')
  assert.equal(call(own, 0, 1).mechanism, 'none')
  // But where a second cell carries the loss, the gain in the first faces it: reciprocal, and read
  // once for the embryo, so the first cell's loss and the second's read the same.
  const faced = [unit([ch('1', 'copy-gain'), ch('1', 'parental-absence')]), unit([ch('1', 'copy-loss')])]
  assert.equal(call(faced, 0, 1).mechanism, 'mitotic-reciprocal')
  assert.equal(call(faced, 1).mechanism, 'mitotic-reciprocal')
}

// --- 6. "ABSENT" NEEDS A REPORTABLE CELL THAT READS THE CHROMOSOME NORMAL ----------------------
{
  const g = unit([ch('21', 'copy-gain')])
  assert.equal(call([g, normal, normal]).mechanism, 'mitotic-not-all')
  assert.ok(call([g, normal, normal]).why.startsWith('1 of 3 units carry this gain, 2 read chromosome 21 normal'))
  assert.equal(call([g, unit([], { reportable: false })]).mechanism, 'none', 'an unreportable cell reads nothing normal')
  assert.equal(call([g, unit([], { measured: ['1'] })]).mechanism, 'none', 'an unmeasured chromosome is not normal')
  const f = [g, unit([ch('21', 'cnn-loh', 'unclear', 5e6, 9e6, false)])]
  assert.equal(call(f).mechanism, 'none', 'a piece of any kind is not normal where it lies')
  assert.deepEqual(stretchesOf(f[0].changes[0], f).map((s) => [s.fromBp, s.toBp, s.mechanism]),
    [[0, 5e6, 'mitotic-not-all'], [5e6, 9e6, 'none'], [9e6, TO_END, 'mitotic-not-all']])
}

// --- 7. PIECES OVERLAP; A WHOLE ROW SPANS ITS CHROMOSOME ---------------------------------------
{
  const p = (a: number, b: number) => ch('7', 'copy-loss', 'unclear', a, b, false)
  assert.equal(call([unit([p(5e6, 25e6)]), unit([p(4.8e6, 25.2e6)]), unit([p(5.5e6, 25.6e6)])]).mechanism, 'meiotic',
    'edges under 1 Mb apart are one breakpoint')
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

// --- 9. ONE STRETCH, ONE WORD IN EVERY CELL (D5) ------------------------------------------------
{
  // The S2 rows 96 and 97 shape: one cell loses chromosome 17 whole, its sister loses 0 to 52.6 Mb
  // and gains 52.6 to 78.4 Mb, both limited (so neither reads a stretch normal). Edges a few kb from
  // the chromosome start and from each other are one breakpoint.
  const whole = unit([ch('17', 'copy-loss', 'unclear', 0, 0, true)], { reportable: false })
  const parts = unit([ch('17', 'copy-loss', 'unclear', 6157, 52_598_275, false),
    ch('17', 'copy-gain', 'unclear', 52_618_539, 78_359_420, false)], { reportable: false })
  const [[w], [l, g]] = timeGroup([whole, parts])
  assert.equal(w.mechanism, 'none', 'a whole call over stretches that read differently states no one word')
  assert.equal(w.why, 'its stretches read differently')
  assert.deepEqual(w.stretches.map((s) => [stretchText(s), stretchWord(s)]), [
    ['0.0 to 52.6 Mb', 'consistent with the gamete (two cells read)'],
    ['52.6 to 78.4 Mb', 'after fertilisation'],
    ['78.4 Mb to the end', 'not called'],
  ])
  assert.equal(l.mechanism, 'consistent-with-gamete')
  assert.equal(l.stretches.length, 1)
  assert.equal(g.mechanism, 'mitotic-reciprocal')
  assert.ok(split(w.stretches) && !split(l.stretches))
  // Every stretch prints one word whichever cell it is read from.
  const word = (t: typeof w, at: number) => t.stretches.find((s) => s.fromBp <= at && at < s.toBp)?.mechanism
  for (const at of [1e6, 30e6, 60e6]) {
    const here = [w, l, g].map((t) => word(t, at)).filter(Boolean)
    assert.equal(new Set(here.map((m) => ORIGIN_WORD[m!])).size, 1, `one word at ${at}`)
  }
}

// --- 10. NEAR-EQUAL BREAKPOINTS ARE ONE, AND THE JOIN IS NOT WIDER THAN IT SAYS ---------------------
{
  const p = (a: number, b: number) => ch('7', 'copy-loss', 'unclear', a, b, false)
  const near = [unit([p(10e6, 40e6)]), unit([p(10e6 + SAME_BREAK_BP - 1, 40e6)])]
  assert.deepEqual(stretchesOf(near[0].changes[0], near).map((s) => s.mechanism), ['consistent-with-gamete'])
  const far = [unit([p(10e6, 40e6)]), unit([p(10e6 + SAME_BREAK_BP, 40e6)])]
  assert.deepEqual(stretchesOf(far[0].changes[0], far).map((s) => [s.fromBp, s.mechanism]),
    [[10e6, 'mitotic-not-all'], [11e6, 'consistent-with-gamete']])
  // A call shorter than the join is read over its own extent.
  const tiny = [unit([p(20e6, 20.4e6)]), unit([p(20.1e6, 20.5e6)])]
  assert.deepEqual(stretchesOf(tiny[0].changes[0], tiny).map((s) => [s.fromBp, s.toBp, s.mechanism]),
    [[20e6, 20.4e6, 'consistent-with-gamete']])
}

// --- 11. PARENTS ARE READ ONCE FOR THE EMBRYO --------------------------------------------------
{
  // A cell naming no parent no longer reads "from the gamete" while two others name different ones.
  const m = (o: string) => unit([ch('13', 'copy-loss', o)])
  const units = [m('maternal'), m('unclear'), m('paternal')]
  assert.deepEqual(units.map((_, i) => call(units, i).mechanism), ['none', 'none', 'none'])
  assert.ok(call(units, 1).why.includes('they name different parents'))
}

// --- 12. ORIGIN IS READ AMONG REPORTABLE CELLS ONLY (rule R) ------------------------------------
{
  // A sample not reportable gets no word, and is neither a carrier nor a cell without it.
  const g = () => unit([ch('21', 'copy-gain')])
  const hidden = unit([ch('21', 'copy-gain')], { shown: false })
  const [[a], [b], [c]] = timeGroup([g(), g(), hidden])
  assert.equal(c.mechanism, 'none')
  assert.ok(c.why.includes('not reportable'), c.why)
  assert.equal(a.mechanism, 'consistent-with-gamete', 'two shown units read as two, not three')
  assert.ok(a.why.startsWith('2 of 2 units'), a.why)
  assert.equal(b.mechanism, 'consistent-with-gamete')
  // Left with one shown unit, the embryo reads as one cell.
  const [[one], rest] = timeGroup([g(), { ...normal, shown: false }])
  assert.equal(one.mechanism, 'none', 'a cell not reportable cannot read the chromosome normal')
  assert.equal(rest.length, 0, 'it has no change to time')
  assert.equal(readingOf([], {}, false).shown, false)
  assert.equal(readingOf([], {}).shown, true)
}

console.log('timing.check.ts OK: rule B1 read per stretch of the embryo, with direction, parent and measured normal; one cell prints no origin')

// --- THE REPORT'S TIMING TABLE PRINTS NO COORDINATES FOR A WHOLE ROW -----------------------------
// A whole row's end is the chromosome's length only where a parent was loaded; with none it is 0.
{
  const { timingWhere } = await import('./syngamyPdf.ts')
  const [[w, s]] = timeGroup([unit([ch('16', 'monosomy', 'unclear', 0, 0), ch('7', 'copy-loss', 'unclear', 5e6, 9e6, false)]), normal])
  assert.equal(timingWhere(w), 'chr16, whole chromosome')
  assert.equal(timingWhere(s), 'chr7 5.0-9.0Mb')
  console.log('timing.check.ts: the report prints a whole row as the whole chromosome')
}

// --- THE OPENED ROW AND THE REPORT GIVE EACH STRETCH ITS OWN LINE -----------------------------------
{
  const { timingRows } = await import('./syngamyPdf.ts')
  const { whenLines } = await import('./pageSays.ts')
  const whole = unit([ch('17', 'copy-loss', 'unclear', 0, 0, true)], { reportable: false })
  const parts = unit([ch('17', 'copy-loss', 'unclear', 6157, 52_598_275, false),
    ch('17', 'copy-gain', 'unclear', 52_618_539, 78_359_420, false)], { reportable: false })
  const [[w], [l]] = timeGroup([whole, parts])
  const rows = timingRows([w, l])
  // The stretch first, the call named as what it belongs to: never "whole chromosome, 0.0 to 52.6 Mb".
  assert.deepEqual(rows.map((r) => r.where), [
    'chr17 0.0 to 52.6 Mb, part of the whole-chromosome loss', 'chr17 52.6 to 78.4 Mb, part of the whole-chromosome loss',
    'chr17 78.4 Mb to the end, part of the whole-chromosome loss', 'chr17 0.0-52.6Mb'])
  assert.equal(stretchWhere({ chrom: '9', kind: 'parental-absence', startBp: 7.15e6, endBp: 139.4e6, whole: false },
    { fromBp: 16.7e6, toBp: 133.9e6 }), 'chr9 16.7 to 133.9 Mb, part of the loss over 7.2 to 139.4 Mb')
  assert.deepEqual(rows.map((r) => r.mechanism), ['consistent-with-gamete', 'mitotic-reciprocal', 'none', 'consistent-with-gamete'])
  const lines = whenLines({ pieces: [{ chrom: '17', kind: 'copy-loss', locus: 'chr17', startBp: 0, endBp: 0,
    mechanismStretches: w.stretches }] } as never)
  assert.equal(lines.length, 3)
  assert.ok(lines[0].startsWith('When it arose, chr17 0.0 to 52.6 Mb, part of the whole-chromosome loss: '
    + 'consistent with the gamete (two cells read). 2 of 2 units'), lines[0])
  assert.ok(lines[2].startsWith('When it arose, chr17 78.4 Mb to the end, part of the whole-chromosome loss: not called.'), lines[2])
  assert.deepEqual(whenLines({ pieces: [{ mechanism: l.mechanism, mechanismWhy: l.why }] } as never),
    [`When it arose: consistent with the gamete (two cells read). ${l.why}.`])
  console.log('timing.check.ts: a call read by stretch gets one line per stretch in the opened row and the report')
}

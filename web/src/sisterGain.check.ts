/**
 * A sister's loss is read as a gain only where this cell's own intensity is up over that interval.
 */
import assert from 'node:assert/strict'
import {
  lrrWindows, readInterval, sisterGains, oneGenome, sisterUnitOf, WHOLE, WINDOW_BP,
  MIN_WINDOW_MARKERS, SISTER_SHIFT, type LrrWindow, type SisterUnit,
} from './sisterGain.ts'
import type { AB } from './informativity.ts'

// A deterministic stream, so a failure is the code changing and never the draw.
let seed = 99
const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 - 0.5 }

/** A cell: 22 chromosomes of 20 windows, noise of the given size, and a shift on some intervals. */
function cell(noise: number, shifts: { chrom: string; from: number; to: number; by: number }[] = []): LrrWindow[] {
  const out: LrrWindow[] = []
  for (let c = 1; c <= 22; c += 1) {
    for (let i = 0; i < 20; i += 1) {
      const s = shifts.find((x) => x.chrom === String(c) && i >= x.from && i < x.to)
      out.push([String(c), i, rnd() * noise + (s ? s.by : 0)])
    }
  }
  return out
}
const all = Array.from({ length: 22 }, (_, i) => String(i + 1))
const without = (...cs: string[]) => all.filter((c) => !cs.includes(c))

/** Homozygous calls, alternating, so two of these differ only where one is flipped. */
const genome = (n: number): Map<string, AB> => {
  const g = new Map<string, AB>()
  for (let i = 0; i < n; i += 1) g.set(`p${i}`, i % 2 ? 'BB' : 'AA')
  return g
}
/** The same genome with every `every`-th call flipped, so it agrees with `genome` at 1 - 1/every. */
const flipped = (n: number, every: number): Map<string, AB> => {
  const g = genome(n)
  for (let i = 0; i < n; i += every) g.set(`p${i}`, g.get(`p${i}`) === 'AA' ? 'BB' : 'AA')
  return g
}
/** Probe k on chromosome (k mod 22) + 1, so a chromosome can be left out of the comparison. */
const CHROM_OF = new Map(Array.from({ length: 12_000 },
  (_, i) => [`p${i}`, { chrom: String((i % 22) + 1) }] as const))
const N = 12_000
const ONE = genome(N)
const ONE_NOISY = flipped(N, 20)    // agrees 0.95: one genome read twice through dropout
const OTHER = flipped(N, 4)         // agrees 0.75: a different genome
/** Under the floor `agreement` refuses below, so the pair is unmeasurable however identical. */
const THIN = genome(1000)

// One genome read twice counts; two genomes do not; and an unmeasurable pair does not, however
// identical, because the fault this rule stops is reading a gain off a cell of another genome.
{
  const u = (gt: Map<string, AB>, changed: string[] = []): SisterUnit => ({
    windows: [], gt, pos: CHROM_OF, reference: [], losses: [], changed: new Set(changed),
  })
  assert.ok(oneGenome(u(ONE), u(ONE_NOISY)))
  assert.ok(!oneGenome(u(ONE), u(OTHER)))
  assert.ok(!oneGenome(u(THIN), u(THIN)), 'identical but unmeasurable is not a sister')

  // THE RECIPROCAL EVENT MUST NOT COST THE PAIR ITS SISTERHOOD. A real pair of sisters sits just
  // above the threshold on dropout alone, so chromosomes reading opposite end to end, as two cells
  // that lost opposite parental copies there do, push it under. Compared whole it fails; compared off
  // the chromosomes each cell reports a change on it passes.
  const RECIP = ['7', '9', '11', '13', '15']
  const recip = new Map(ONE_NOISY)
  for (const [probe, c] of CHROM_OF) {
    if (RECIP.includes(c.chrom)) recip.set(probe, ONE.get(probe) === 'AA' ? 'BB' : 'AA')
  }
  assert.ok(!oneGenome(u(ONE), u(recip)), 'chromosomes read opposite sink the whole-genome comparison')
  assert.ok(oneGenome(u(ONE, RECIP), u(recip, RECIP)), 'left out, the pair is one genome again')
  // Only a chromosome a cell actually reports is left out: unreported ones still count against it.
  assert.ok(!oneGenome(u(ONE, ['9']), u(recip, ['9'])), 'leaving out one of the three is not enough')
}

// One derivation of a unit, so the page and the harness cannot drift. A reported change takes its
// chromosome out of the reference and out of what is asked; only losses are offered to a sister.
{
  const u = sisterUnitOf({
    chroms: [{ chrom: '5', aneuploidy: 'loss' }, { chrom: '9', aneuploidy: 'gain' },
      { chrom: '7', aneuploidy: null }, { chrom: 'X', aneuploidy: 'loss' }],
    segments: [{ chrom: '16', kind: 'copy-loss', start: 60e6, end: 90e6 },
      { chrom: '2', kind: 'parental-absence', start: 0, end: 1e6 }],
  }, [], ONE, CHROM_OF)
  assert.deepEqual(u.losses.map((l) => [l.chrom, l.whole]), [['5', true], ['16', false]],
    'a gain is not a loss, a sex chromosome is not asked, and an absence is not a copy change')
  assert.deepEqual([...u.changed].sort(), ['16', '5', '9'])
  assert.ok(!u.reference.includes('5') && !u.reference.includes('9') && u.reference.includes('7'))
}

// Windows are medians of the markers in each 5 Mb, and a thin window is dropped.
{
  const ms = [
    ...Array.from({ length: MIN_WINDOW_MARKERS }, (_, i) => ({ pos: i * 1000, log2R: 0.2 })),
    ...Array.from({ length: MIN_WINDOW_MARKERS - 1 }, (_, i) => ({ pos: WINDOW_BP + i * 1000, log2R: 0.9 })),
    { pos: 3 * WINDOW_BP, log2R: null },
  ]
  assert.deepEqual(lrrWindows(new Map([['7', ms]])), [['7', 0, 0.2]])
}

// A whole chromosome up by a third copy is read; the same chromosome at the cell's own level is not.
{
  const gained = cell(0.1, [{ chrom: '5', from: 0, to: 20, by: 0.4 }])
  const r = readInterval(gained, without('5'), '5', 0, Number.MAX_SAFE_INTEGER)
  assert.ok(r && r.shift > 0.35 && r.z > 3 && r.windows === 20)
  const flat = readInterval(cell(0.1), all, '5', 0, Number.MAX_SAFE_INTEGER)
  assert.ok(flat && Math.abs(flat.shift) < 0.1 && flat.z < 3)
  // A segment: the last 40 Mb of chromosome 16.
  const seg = cell(0.1, [{ chrom: '16', from: 12, to: 20, by: 0.4 }])
  const inside = readInterval(seg, without('16'), '16', 12 * WINDOW_BP, 20 * WINDOW_BP)
  const outside = readInterval(seg, without('16'), '16', 0, 12 * WINDOW_BP)
  assert.ok(inside && inside.shift > 0.35 && inside.z > 3 && inside.windows === 8)
  assert.ok(outside && outside.shift < SISTER_SHIFT)
  // Too little to read against, or too short an interval, is not read.
  assert.equal(readInterval(seg, ['1', '2', '3'], '16', 0, Number.MAX_SAFE_INTEGER), null)
  assert.equal(readInterval(seg, without('16'), '16', 0, WINDOW_BP), null)
}

// Three cells of one embryo. The first lost chromosome 5 and the end of 16; the second holds what it
// lost; the third holds neither. Only the second reports gains, and only over what the first lost.
{
  const unit = (windows: LrrWindow[], over: Partial<SisterUnit> = {}): SisterUnit => ({
    windows, gt: ONE, pos: CHROM_OF, reference: all, losses: [], changed: new Set(), ...over,
  })
  const loser = unit(cell(0.1), {
    reference: without('5', '16'),
    losses: [WHOLE('5'), { chrom: '16', startBp: 12 * WINDOW_BP, endBp: 20 * WINDOW_BP, whole: false }],
    changed: new Set(['5', '16']),
  })
  const gainer = unit(cell(0.1, [
    { chrom: '5', from: 0, to: 20, by: 0.4 }, { chrom: '16', from: 12, to: 20, by: 0.4 },
    { chrom: '9', from: 0, to: 20, by: 0.4 },
  ]), { reference: without('5', '16', '9') })
  const bystander = unit(cell(0.1))
  const [a, b, c] = sisterGains([loser, gainer, bystander])
  assert.deepEqual(a, [], 'the cell that lost them is not asked about its own losses')
  assert.deepEqual(b.map((g) => [g.chrom, g.whole, g.from]), [['5', true, 0], ['16', false, 0]])
  assert.deepEqual(c, [])
  // Chromosome 9 is up in the gainer and no sister lost it: nothing here reports it.
  assert.ok(!b.some((g) => g.chrom === '9'))
  // A chromosome the cell already reports a change on is not asked again.
  const [, again] = sisterGains([loser, { ...gainer, changed: new Set(['5']) }])
  assert.deepEqual(again.map((g) => g.chrom), ['16'])
}

// THE FAULT THIS RULE EXISTS TO STOP. The same two cells, the same intensity, the same reported loss:
// a gain is read only where the pair reads as one genome. A different genome in the other slot, or a
// pair too thin to measure, reports nothing.
{
  const loser = (gt: Map<string, AB>): SisterUnit => ({
    windows: cell(0.1), gt, pos: CHROM_OF, reference: without('5'), losses: [WHOLE('5')],
    changed: new Set(['5']),
  })
  const gainer: SisterUnit = {
    windows: cell(0.1, [{ chrom: '5', from: 0, to: 20, by: 0.4 }]),
    gt: ONE, pos: CHROM_OF, reference: without('5'), losses: [], changed: new Set(),
  }
  const read = (gt: Map<string, AB>) => sisterGains([loser(gt), gainer])[1].map((g) => g.chrom)
  assert.deepEqual(read(ONE_NOISY), ['5'], 'one genome read twice: the gain is read')
  assert.deepEqual(read(OTHER), [], 'a different genome in the sister slot reports nothing')
  assert.deepEqual(read(THIN), [], 'too few shared calls reports nothing')
}

// A CELL WITH NO SISTER IN THE RUN IS UNTOUCHED (gate S5). One cell alone, and a cell whose only
// companion is another genome, both report nothing: this reading is additive and never subtracts.
{
  const lone: SisterUnit = {
    windows: cell(0.1, [{ chrom: '5', from: 0, to: 20, by: 0.4 }]),
    gt: ONE, pos: CHROM_OF, reference: without('5'), losses: [WHOLE('5')], changed: new Set(['5']),
  }
  assert.deepEqual(sisterGains([lone]), [[]], 'one cell alone has no sister and reports nothing')
  const stranger: SisterUnit = { ...lone, gt: OTHER, losses: [WHOLE('5')] }
  const [a, b] = sisterGains([lone, stranger])
  assert.deepEqual([a, b], [[], []], 'a companion of another genome is not a sister')
}

console.log('sisterGain.check ok')

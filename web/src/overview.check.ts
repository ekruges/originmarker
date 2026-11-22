// Self-check for the run's overview table. Run: node src/overview.check.ts
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { normaliseTrack } from './comparison.ts'
import { armsFrom, extentOf, overviewSampleOf, overviewTable, type OverviewSample } from './overview.ts'

const arms = armsFrom(normaliseTrack(JSON.parse(readFileSync(new URL('../public/hg19_features.json', import.meta.url), 'utf8')))!.gaps)
const one = arms.get('1')!
assert.ok(one.cenStart > 120e6 && one.cenEnd < 126e6 && one.end > 249e6, 'chromosome 1 from the shipped gap table')

// --- 1. WHOLE, ARM, PIECE ---------------------------------------------------------------------
{
  const at = (startBp: number, endBp: number) => extentOf({ chrom: '1', whole: false, startBp, endBp }, arms)
  assert.equal(extentOf({ chrom: '1', whole: true, startBp: 0, endBp: 0 }, arms), 'whole')
  assert.equal(at(1e6, one.cenStart - 2e6), 'arm', 'end to centromere')
  assert.equal(at(one.cenEnd + 3e6, one.end - 1e6), 'arm', 'centromere to the other end')
  assert.equal(at(1e6, 60e6), 'piece')
  assert.equal(at(one.cenEnd + 1e6, one.cenEnd + 4e6), 'piece', 'a piece beside the centromere is not an arm')
  assert.equal(extentOf({ chrom: '1', whole: false, startBp: 1e6, endBp: 120e6 }, new Map()), 'piece')
}

// --- 2. CARD ROWS IN, ONE ROW ONE CHANGE, NOT-REPORTABLE ADDS NOTHING -------------------------
{
  const piece = (kind: string, mechanism?: string) => ({ kind, startBp: 0, endBp: 249e6, mechanism }) as never
  const list = (rows: unknown[]) => ({ counted: rows, once: [], seen: [] }) as never
  const row = (pieces: unknown[], gain = false, origin = 'maternal') =>
    ({ chrom: '1', gain, whole: true, origin, startBp: 0, endBp: 249e6, pieces })
  // One loss printed as a copy loss and a parental absence is one card row and one change.
  const s = overviewSampleOf(list([row([piece('copy-loss', 'meiotic'), piece('parental-absence', 'meiotic')])]),
    { stage: { stage: 'blastomere' }, refusals: [] }, 0)
  assert.equal(s.rows.length, 1)
  assert.deepEqual(s.rows[0].when, ['from the gamete'])
  const bad = overviewSampleOf(list([row([piece('copy-gain')], true)]),
    { stage: { stage: 'blastomere' }, refusals: [{ claim: 'sample', gate: 'integrity.uninterpretable' }] })
  const t = new Map(overviewTable([s, bad] as OverviewSample[], arms))
  assert.equal(t.get('samples'), 2)
  assert.equal(t.get('not reportable (integrity uninterpretable)'), 1)
  assert.equal(t.get('changes (as the result cards count them)'), 1, 'a sample not reportable contributes no change')
  assert.equal(t.get('gains'), 0)
  assert.equal(t.get('origin: from the gamete'), 1)
  assert.equal(t.get('embryos of two or more units'), 1)
  assert.equal(t.get('samples in no embryo group'), 1)
  // A row whose pieces read differently is counted on its own line.
  const mixed = overviewSampleOf(list([row([piece('copy-loss', 'meiotic'), piece('copy-loss', 'mitotic-not-all')])]),
    { stage: { stage: 'bulk' } })
  assert.equal(new Map(overviewTable([mixed], arms)).get('origin: pieces read differently'), 1)
}

// --- 3. NO EM DASH IN ANY LABEL ----------------------------------------------------------------
for (const [label] of overviewTable([], arms)) assert.ok(!label.includes('\u2014'), label)

console.log('overview.check.ts OK: card rows counted once each, arms from the shipped gap table, unreportable samples add nothing')

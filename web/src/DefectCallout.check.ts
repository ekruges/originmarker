// Self-check for the defect display. Run: node src/DefectCallout.check.ts
//
// The component is markup and is not tested here. What IS tested is defectsFrom, which decides
// what the reader is told: which origin is shown, and what is said when there is none. Getting
// that wrong shows a confident parent for a region that never had one, which is the single most
// damaging thing this UI could do.
import assert from 'node:assert/strict'
import { defectsFrom, headline, bandColour, sentences } from './defects.ts'
import type { Segment } from './segments.ts'
import type { GainAnnotation } from './parentage.ts'

const seg = (chrom: string, startBp: number, endBp: number, kind: Segment['kind']): Segment =>
  ({ chrom, startBp, endBp, kind, score: 100, markers: 500, absent: 400, rate: 0.8 } as Segment)

// --- 1. an annotated region carries its origin through -------------------------------------------
{
  const s = seg('4', 10_000_000, 20_000_000, 'copy-loss')
  const ann: GainAnnotation[] = [{
    where: 'chr4 10.0-20.0Mb', kind: 'segment', origin: 'maternal',
    why: 'maternal alleles under-represented', called: true,
  }]
  const [d] = defectsFrom([s], [], ann)
  assert.equal(d.origin, 'maternal')
  assert.equal(d.basis, 'two-parent', 'an annotated region says what the call rests on')
  assert.ok(d.why.includes('under-represented'))
}

// --- 2. an UNannotated region says so, and does not borrow a neighbour's origin -------------------
//
// The failure this guards against: a region with no annotation picking up the previous one's
// origin through a loose match, and being displayed as a confident parental call.
{
  const s = seg('9', 30_000_000, 35_000_000, 'copy-loss')
  const elsewhere: GainAnnotation[] = [{
    where: 'chr4 10.0-20.0Mb', kind: 'segment', origin: 'paternal',
    why: 'paternal', called: true,
  }]
  const [d] = defectsFrom([s], elsewhere, [])
  assert.equal(d.origin, 'unclear', 'a region with no annotation must not borrow another one')
  assert.equal(d.basis, undefined, 'and must not claim a basis it does not have')
  assert.ok(d.why.includes('not determined'), 'and must say the parent was not determined')
  // The coordinates must survive: position needs no parent and is reported regardless.
  assert.equal(d.startBp, 30_000_000)
  assert.equal(d.endBp, 35_000_000)
}

// --- 3. origin strings that merely START with a parent name are read correctly --------------------
//
// The uniparental path writes "this parent, other homologue (meiotic)", which names no parent for
// display purposes and must not be parsed as one.
{
  const s = seg('1', 0, 5_000_000, 'copy-gain')
  const ann: GainAnnotation[] = [{
    where: 'chr1 0.0-5.0Mb', kind: 'segment',
    origin: 'this parent, other homologue (meiotic)', why: 'homologue', called: true,
  }]
  const [d] = defectsFrom([s], ann, [])
  assert.equal(d.origin, 'unclear', 'a homologue call names no parent and must not display as one')
}

// --- 4. every segment produces exactly one row, in order ------------------------------------------
{
  const segs = [
    seg('1', 0, 1_000_000, 'copy-loss'),
    seg('2', 0, 2_000_000, 'copy-gain'),
    seg('3', 0, 3_000_000, 'parental-absence'),
  ]
  const out = defectsFrom(segs, [], [])
  assert.equal(out.length, 3, 'no segment may be dropped from the display')
  assert.deepEqual(out.map((d) => d.chrom), ['1', '2', '3'])
  assert.deepEqual(out.map((d) => d.kind), ['copy-loss', 'copy-gain', 'parental-absence'])
  assert.ok(out.every((d) => d.origin === 'unclear'), 'no annotations means no named parents')
}

// --- THE CONFIDENCE TRAVELS WITH THE PARENT --------------------------------------------------------
//
// A named parent with no number reads as certain. On amplified material most are not: a top-band
// call measures 0.995 and the weakest band 0.62, and those must not look alike at a glance.
{
  const base = { chrom: '7', startBp: 1_000_000, endBp: 9_000_000, locus: 'chr7:1,000,000-9,000,000',
    kind: 'copy-loss' as const, origin: 'paternal' as const, why: 'x' }
  const strong = headline({ ...base, confidence: 0.9961, band: 'A' })
  assert.ok(strong.includes('0.996'), `the number must be in the headline: ${strong}`)
  assert.ok(strong.includes('very confident'), `and its band word: ${strong}`)

  const weak = headline({ ...base, confidence: 0.6184, band: 'D' })
  assert.ok(weak.includes('0.618'))
  assert.ok(weak.includes('not for reporting'),
    'the weakest band must say so in the headline, not only in a chip')

  // Dimmed for C and D, inked for A and B, so weak and strong differ before anything is read.
  assert.equal(bandColour({ ...base, band: 'A' }), 'var(--om-defect)')
  assert.equal(bandColour({ ...base, band: 'D' }), 'var(--om-text-dim)')
  assert.equal(bandColour({ ...base, origin: 'unclear', band: 'A' }), 'var(--om-text-dim)')

  // No confidence at all must still render, since the genotype channel supplies none.
  assert.ok(headline(base).includes('PATERNAL'))
  // And no em dash anywhere in user-facing text.
  assert.ok(!headline({ ...base, origin: 'unclear' }).includes('\u2014'),
    'no em dashes in user-facing strings')
}

console.log('DefectCallout.check.ts: all assertions passed, including the no-borrowed-origin case')

// --- 5. a one-parent call fills the origin, and maps to the right parent -------------------------
//
// The path that matters for a sperm-only run. 'known-parent-lost' means the LOADED parent's copy
// is gone, so with the sperm donor loaded it reads paternal; 'other-parent-lost' reads maternal.
// Getting that mapping backwards would invert every call in the display while looking healthy.
{
  const s = seg('7', 1_000_000, 9_000_000, 'copy-loss')
  const one = [{
    where: 'chr7 1.0-9.0Mb', verdict: 'known-parent-lost', posterior: 0.999,
    markers: 4200, exclusive: 610, why: 'known-parent-lost at posterior 0.999',
  }]
  const [d] = defectsFrom([s], [], [], one, 'paternal')
  assert.equal(d.origin, 'paternal', 'the loaded parent losing its copy is that parent')
  assert.equal(d.basis, 'one-parent')
  assert.equal(d.informative, 4200)
  assert.equal(d.posterior, 0.999)

  const [flip] = defectsFrom([s], [], [], [{ ...one[0], verdict: 'other-parent-lost' }], 'paternal')
  assert.equal(flip.origin, 'maternal', 'the other parent losing its copy is the other parent')

  // With the MOTHER loaded the same verdicts must map the other way round.
  const [asMum] = defectsFrom([s], [], [], one, 'maternal')
  assert.equal(asMum.origin, 'maternal')

  // A refusal must not fill anything in.
  const [ref] = defectsFrom([s], [], [], [{ ...one[0], verdict: 'refused', posterior: NaN }],
    'paternal')
  assert.equal(ref.origin, 'unclear', 'a refused one-parent call names nobody')
  assert.equal(ref.basis, undefined)

  // Two-parent evidence must win where both exist.
  const ann: GainAnnotation[] = [{
    where: 'chr7 1.0-9.0Mb', kind: 'segment', origin: 'maternal', why: 'dosage', called: true,
  }]
  const [both] = defectsFrom([s], [], ann, one, 'paternal')
  assert.equal(both.origin, 'maternal', 'two-parent dosage outranks the single-parent call')
  assert.equal(both.basis, 'two-parent')
}

console.log('DefectCallout.check.ts: one-parent mapping pinned in both parent roles')

// --- THE FACE MUST NOT CONTRADICT ITSELF -------------------------------------------------------
// A band F row NAMES a parent from evidence that reaches no measured band. Printing "no call"
// beside that parent says two opposite things in one line, which a real run produced.
{
  const bandWord = (band: string): string => (band === 'C' ? 'weak'
    : band === 'D' ? 'not evaluable' : band === 'F' ? 'direction only' : '')
  for (const band of ['C', 'D', 'F']) {
    assert.ok(bandWord(band).length > 0,
      `band ${band} suppresses its digits, so it must say what it is instead of falling through `
      + 'to the wording reserved for a row that named nothing')
    assert.notEqual(bandWord(band), 'no call',
      `band ${band} names a parent, so it must never read "no call"`)
  }
  console.log('  bands C, D and F each say what they are rather than reading as no call')
}

// --- 6. sentence splitting must not cut decimals in half -----------------------------------------
//
// Reason text is full of measured numbers. A splitter that breaks at every period turns "0.78x"
// into "0." and "78x": the number is destroyed on rejoin, and the shared/unique grouping compares
// fragments, so a half-sentence leads the shared block.
{
  const s = 'An intact chromosome calls at 0.78x to 1.16x of its median and never leaves '
    + '-0.79 to +0.42 log2. This sample carries ONE parental genome at 5.0x.'
  const out = sentences(s)
  assert.equal(out.length, 2, `two sentences, got ${out.length}: ${JSON.stringify(out)}`)
  assert.ok(out[0].startsWith('An intact'), 'no leading fragment')
  assert.equal(out.join(' '), s, 'rejoining must reproduce the input exactly')
  for (const n of ['0.78x', '1.16x', '-0.79', '+0.42', '5.0x']) {
    assert.ok(out.join(' ').includes(n), `${n} must survive intact`)
  }
  // A trailing sentence with no terminator is still one sentence.
  assert.deepEqual(sentences('Only this'), ['Only this'])
  assert.deepEqual(sentences(''), [])
}

console.log('DefectCallout.check.ts: decimals survive sentence splitting')

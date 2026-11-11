/**
 * Self-check for the information ratio.
 *
 * The three cases that matter are a sample that agrees with the reference, one whose calls are
 * independent of it, and a mixture. The first two fix the scale at both ends; the third is the
 * claim the shipped reason makes in words and would be the costliest to get wrong, because it is
 * the one that tells an operator their sample is not contaminated.
 */
import assert from 'node:assert'
import { readFileSync } from 'node:fs'
import {
  emptyHomTally, addHomozygous, informationRatio, informationBand, whyNoInformation,
  INFORMATIVE_MAX, NO_INFORMATION_MIN, MIN_SHARED_HOMOZYGOTES, FLIPPED_CODING_MIN,
} from './informationRatio.ts'
import type { AB } from './informativity.ts'

/** A deterministic generator, so a failure is reproducible rather than occasional. */
const rng = (seed: number) => () => {
  seed = (seed * 1103515245 + 12345) & 0x7fffffff
  return seed / 0x7fffffff
}

const N = 20_000
/** Reference genotypes at a realistic minor-allele spectrum, homozygous only. */
const refOf = (rand: () => number): AB[] =>
  Array.from({ length: N }, () => (rand() < 0.62 ? 'AA' : 'BB'))

// A SAMPLE THAT IS THE SAME PERSON never reads an opposite homozygote, so R is 0.
{
  const rand = rng(1); const ref = refOf(rand)
  const t = emptyHomTally()
  for (const g of ref) addHomozygous(t, g, g)
  assert.equal(informationRatio(t), 0, 'identical calls carry maximal information')
  assert.equal(informationBand(0), 'informative')
}

// CALLS THAT SAY NOTHING sit at 1, which is what fixes the top of the scale. Drawn independently
// of the reference, at the same AA fraction, so composition cannot explain the result.
{
  const rand = rng(2); const ref = refOf(rand)
  const t = emptyHomTally()
  for (const g of ref) addHomozygous(t, rand() < 0.62 ? 'AA' : 'BB', g)
  const r = informationRatio(t)
  assert.ok(Math.abs(r - 1) < 0.04, `independent calls read about 1, got ${r}`)
  assert.equal(informationBand(r), 'no-information')
  assert.ok(whyNoInformation(r, 'the father')?.includes('mostly noise'))
}

// A MIXTURE READS LOW, which is the claim the shipped reason makes in words. Two contributors: a
// marker is called homozygous only where both carry that allele, so the surviving homozygotes
// agree with the reference MORE than chance.
{
  const rand = rng(3); const ref = refOf(rand)
  const t = emptyHomTally()
  for (const g of ref) {
    // The reference is one of the two contributors; the other is an unrelated person.
    const other: AB = rand() < 0.62 ? 'AA' : 'BB'
    // Where the two contributors disagree the mixture reads heterozygous and is skipped.
    addHomozygous(t, g === other ? g : 'AB', g)
  }
  const r = informationRatio(t)
  assert.ok(r < INFORMATIVE_MAX, `a mixture must read low, not at 1. Got ${r}`)
  assert.equal(whyNoInformation(r, 'the father'), undefined,
    'and must never be given the no-information reason')
}

// AN INVERTED A/B CODING reads near 2 and is called out as that rather than as a failure, because
// the two call for completely different work.
{
  const rand = rng(4); const ref = refOf(rand)
  const t = emptyHomTally()
  for (const g of ref) addHomozygous(t, g === 'AA' ? 'BB' : 'AA', g)
  const r = informationRatio(t)
  assert.ok(r >= FLIPPED_CODING_MIN, `a flipped coding reads at or above ${FLIPPED_CODING_MIN}, got ${r}`)
  assert.equal(informationBand(r), 'flipped-coding')
  assert.ok(whyNoInformation(r, 'the father')?.includes('inverted A/B coding'),
    'and the reason says so instead of claiming the array failed')
}

// HETEROZYGOUS AND UNCALLED MARKERS CARRY NOTHING HERE and must not enter any of the four counts.
{
  const t = emptyHomTally()
  addHomozygous(t, 'AB', 'AA'); addHomozygous(t, 'AA', 'AB')
  addHomozygous(t, null, 'AA'); addHomozygous(t, 'AA', null)
  assert.deepEqual(t, emptyHomTally(), 'only shared homozygotes are counted')
}

// TOO FEW SHARED MARKERS IS NOT A VERDICT. A thin tally returns NaN rather than a ratio read off
// a handful of calls, and NaN never produces a reason.
{
  const t = emptyHomTally()
  for (let i = 0; i < MIN_SHARED_HOMOZYGOTES - 1; i += 1) addHomozygous(t, 'AA', 'BB')
  assert.ok(Number.isNaN(informationRatio(t)), 'under the floor there is no estimate')
  assert.equal(informationBand(NaN), 'too-few-markers')
  assert.equal(whyNoInformation(NaN, 'the father'), undefined)
}

// THE BANDS ARE CLOSED AT THE VALUES THEY NAME, so a run and the shipped path agree at the edge.
{
  assert.equal(informationBand(INFORMATIVE_MAX), 'informative')
  assert.equal(informationBand(INFORMATIVE_MAX + 1e-9), 'unresolved')
  assert.equal(informationBand(NO_INFORMATION_MIN - 1e-9), 'unresolved')
  assert.equal(informationBand(NO_INFORMATION_MIN), 'no-information')
  assert.equal(whyNoInformation(NO_INFORMATION_MIN - 1e-9, 'the father'), undefined,
    'just under the cut there is no reason')
  assert.ok(whyNoInformation(NO_INFORMATION_MIN, 'the father'), 'and at it there is')
}

// THE REASON NEVER OFFERS A RERUN, which replicate arrays of one amplification product contradict.
{
  const text = whyNoInformation(1.0, 'the father')!
  assert.ok(!/rerun/i.test(text), 'the reason must not send a reader to repeat the run')
  const flat = text.replace(/\s+/g, ' ')
  assert.ok(/unlikely to help/.test(flat) && /not a mixture of two people/.test(flat),
    'it says what the failure is not, and what will not fix it')
}

// REASONS ONLY, PINNED ON THE SHIPPED SOURCE. The branch may add a sentence to an array's limits
// and must never touch a verdict, or the new wording would silently change which chromosomes get
// called. Read off scoreSample.ts rather than asserted in prose, so the claim cannot drift from
// the code it describes.
{
  const src = readFileSync(new URL('./scoreSample.ts', import.meta.url), 'utf8')
  const start = src.indexOf('// WHETHER THE CALLS THEMSELVES SAY ANYTHING')
  assert.ok(start > 0, 'the information-ratio branch is still in scoreSample.ts')
  const end = src.indexOf('// IS THE ARRAY IN THE PARENTAL SLOT', start)
  assert.ok(end > start, 'and its end is still where this check expects it')
  // Comments are not behaviour, and this block's own comment names the things it must not do.
  const block = src.slice(start, end).split('\n')
    .filter((l) => !l.trim().startsWith('//')).join('\n')

  const touches = [...block.matchAll(/result\.([A-Za-z]+)/g)].map((m) => m[1])
  assert.deepEqual([...new Set(touches)], ['limits'],
    `the branch may only write result.limits. It also touches: ${[...new Set(touches)].join(', ')}`)
  for (const forbidden of ['aneuploidy', 'lostParent', 'contribution', 'alleleGain', 'segments',
    'findings', 'stage']) {
    assert.ok(!block.includes(forbidden),
      `a reason branch must not mention ${forbidden}: that is a verdict, not a reason`)
  }
  assert.ok(!/\breturn\b/.test(block), 'and it must not short-circuit the rest of the scoring')
}

console.log('informationRatio.check.ts: the scale, a mixture, a flipped coding, the bands, '
  + 'and that the branch changes reasons only')

// R IS THE SHARE OF RANDOM CALLS. A sample whose homozygous calls are a fraction f drawn at random,
// and otherwise the reference's own, reads R = f: the bar at a half is "more noise than genome".
{
  for (const f of [0.1, 0.3, 0.5, 0.6, 0.9]) {
    const rand = rng(Math.round(f * 1000) + 7); const ref = refOf(rand)
    const t = emptyHomTally()
    for (const g of ref) {
      const s: AB = rand() < f ? (rand() < 0.62 ? 'AA' : 'BB') : g
      addHomozygous(t, s, g)
    }
    const r = informationRatio(t)
    assert.ok(Math.abs(r - f) < 0.03, `a sample ${f} random reads R ${r.toFixed(3)}`)
  }
  assert.equal(informationBand(0.49), 'unresolved', 'under half random, the genome still outweighs the noise')
  assert.equal(informationBand(0.5), 'no-information')
}

// The channel that answers the case a real run could not: 18 samples, 1,272 detected regions, and
// not one parental origin called, because every one was a segment on single-cell material where no
// detection floor exists at any mosaic fraction.
import assert from 'node:assert/strict'
import { uniparentalOrigin, MARGIN_HALF_EVIDENCE } from './uniparentalOrigin.ts'
import { SYSTEMATIC_ERROR_BOUND } from './oneParentOrigin.ts'
import type { Zygosity } from './parentage.ts'

// THE VALUE THE APPLICATION ACTUALLY PRODUCES. Written out rather than paraphrased: the first
// version of this file used 'uniparental', which is not one of the three shipped values, so the
// channel returned null on every real sample and this file passed anyway.
const ZYGOSITY: Zygosity = 'uniparental_homozygous'
const base = { zygosity: ZYGOSITY, role: 'paternal' as const, explainable: 0.0070 }

// --- 1. THE REAL SAMPLE ------------------------------------------------------------------------
// 52461-01_76: gynogenetic, absence 9.04% against a ceiling of 0.70%, which is 12.9x. The loaded
// array is the sperm donor, so a maternal genome is the OTHER parent.
{
  const c = uniparentalOrigin({ ...base, originClass: 'gynogenetic', genomeRate: 0.0904 })
  assert.ok(c, 'a gynogenetic sample at 12.9x its ceiling must produce a call')
  assert.equal(c.parent, 'maternal', 'a gynogenetic genome is maternal')
  assert.equal(c.verdict, 'other-parent', 'with the sperm donor loaded, maternal is the other one')
  assert.ok(c.foldOverCeiling > 12 && c.foldOverCeiling < 13, `12.9x, got ${c.foldOverCeiling}`)
  console.log(`  gynogenetic at ${c.foldOverCeiling.toFixed(1)}x -> ${c.parent},`
    + ` confidence ${c.confidence.toFixed(3)} (band ${c.band})`)
}

// The mirror, and the role flip. Load an OOCYTE instead and the same maternal genome becomes the
// loaded parent. Getting this backwards names the wrong parent on every row.
{
  const pat = uniparentalOrigin({ ...base, originClass: 'androgenetic', genomeRate: 0.0904 })
  assert.equal(pat?.parent, 'paternal')
  assert.equal(pat?.verdict, 'loaded-parent', 'sperm loaded, androgenetic genome is the loaded one')
  const viaOocyte = uniparentalOrigin({
    ...base, role: 'maternal', originClass: 'gynogenetic', genomeRate: 0.0904 })
  assert.equal(viaOocyte?.parent, 'maternal')
  assert.equal(viaOocyte?.verdict, 'loaded-parent',
    'with the oocyte loaded, a gynogenetic genome is the LOADED parent, not the other one')
  console.log('  the loaded/other split follows which array was supplied, both ways')
}

// --- 2. WHERE IT MUST STAY SILENT --------------------------------------------------------------
// This channel exists because a uniparental genome has only one candidate. Anything else has two,
// and answering from zygosity there would be asserting a parent with no basis at all.
{
  for (const originClass of ['biparental', 'unclear']) {
    assert.equal(uniparentalOrigin({ ...base, originClass, genomeRate: 0.0904 }), null,
      `${originClass} has two candidate parents and must fall through to the dosage channel`)
  }
  for (const z of ['unknown', 'diploid'] satisfies Zygosity[]) {
    assert.equal(
      uniparentalOrigin({ ...base, zygosity: z, originClass: 'gynogenetic', genomeRate: 0.0904 }),
      null, `zygosity ${z} gives no basis for the inheritance`)
  }
  // 52461-12_86 from the same run: labelled unclear, absence 11.63% against a ceiling of 22.43%,
  // which is BELOW its own ceiling. Inheriting from a call that was never made would be inventing
  // the evidence.
  assert.equal(
    uniparentalOrigin({ ...base, originClass: 'gynogenetic', genomeRate: 0.1163, explainable: 0.2243 }),
    null, 'below its own ceiling the genome-level call is not made, so nothing is inherited')
  console.log('  silent on biparental, unclear, unresolved zygosity, and below-ceiling absence')
}

// --- 3. THE CEILING ----------------------------------------------------------------------------
// The channel has no injection series, so it must not be able to claim what a calibrated one can.
// Bounded by the same systematic error bound as everything else: 13 validation units, zero events.
{
  const cap = 1 - SYSTEMATIC_ERROR_BOUND
  // An absurd margin, far beyond anything measurable, still cannot buy band A or B.
  const extreme = uniparentalOrigin({ ...base, originClass: 'gynogenetic', genomeRate: 100 })
  assert.ok(extreme, 'an extreme margin still produces a call')
  assert.ok(extreme.confidence < cap + 1e-9,
    `confidence ${extreme.confidence} must stay under the systematic bound ${cap}`)
  assert.ok(extreme.band === 'C' || extreme.band === 'D',
    `an inherited answer must not reach band B, got ${extreme.band}`)
  // And it rises with the evidence rather than being one flat number.
  const weak = uniparentalOrigin({ ...base, originClass: 'gynogenetic', genomeRate: 0.0105 })
  const strong = uniparentalOrigin({ ...base, originClass: 'gynogenetic', genomeRate: 0.0904 })
  assert.ok(weak && strong && strong.confidence > weak.confidence,
    'a more decisive genome-level call must carry more confidence, not the same')
  console.log(`  bounded at ${cap.toFixed(3)}: 1.5x -> ${weak.confidence.toFixed(3)},`
    + ` 12.9x -> ${strong.confidence.toFixed(3)}, extreme -> ${extreme.confidence.toFixed(3)}`)
  assert.ok(MARGIN_HALF_EVIDENCE > 0, 'the half-evidence is a stated constant, not a magic number')
}

console.log('uniparentalOrigin: names a parent where dosage cannot, and stays silent where it must')

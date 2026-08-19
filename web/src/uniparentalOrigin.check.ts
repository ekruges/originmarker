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
const base = { zygosity: ZYGOSITY, role: 'paternal' as const, explainable: 0.0070,
  hetBand: 0.030 }

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


// --- 5. THE GATE MUST NOT BE ONE-SIDED ---------------------------------------------------------
//
// THE DEFECT THIS EXISTS TO PREVENT, measured on this project's own audit set. The first version
// tested absence of the loaded parent for every sample. Absence is large only when the loaded
// parent is MISSING, so with the sperm donor loaded the gate fired on 14 of 14 gynogenetic samples
// and 0 of 8 androgenetic ones: margins 4.77 to 17.44 against 0.09 to 0.62. Half the samples were
// dropped without a word, and any aggregation fed from these rows would have reported a maternal
// excess whatever the biology, because the paternal group could not contribute a single called
// event.
{
  // Androgenetic with the sperm donor loaded: the genome IS the loaded parent, so its absence is
  // low by construction. Real values from the audit set.
  for (const [absenceFold, hetBand] of [[0.62, 0.0285], [0.47, 0.0325], [0.09, 0.0300]]) {
    const c = uniparentalOrigin({
      ...base, originClass: 'androgenetic', role: 'paternal',
      genomeRate: 0.0070 * absenceFold, hetBand,
    })
    assert.ok(c, `androgenetic at absence fold ${absenceFold} must still be called: the loaded `
      + 'parent being present is not a reason to say nothing')
    assert.equal(c.parent, 'paternal')
    assert.equal(c.verdict, 'loaded-parent')
  }
  // And the mirror still works off absence, unchanged.
  const gyno = uniparentalOrigin({
    ...base, originClass: 'gynogenetic', role: 'paternal', genomeRate: 0.0904, hetBand: 0.030 })
  assert.equal(gyno?.parent, 'maternal')
  assert.ok(gyno.foldOverCeiling > 12, 'the absence route keeps its own margin')
  console.log('  both classes are called: absence carries the missing-parent case,'
    + ' zygosity carries the present-parent case')
}

// --- 5b. AND THE PRESENT-PARENT BRANCH MUST REFUSE A GENOME THAT LOOKS DIPLOID ------------------
// Three androgenetic arrays in the audit set sit at hetBand 0.15 to 0.16, ABOVE the 0.08 a second
// parental contribution produces. The absence route was admitting them on evidence that says the
// opposite. Zygosity is the measurement that has to carry this branch, so it is the one that
// refuses.
{
  for (const hetBand of [0.1527, 0.1542, 0.1603]) {
    assert.equal(
      uniparentalOrigin({ ...base, originClass: 'androgenetic', role: 'paternal',
        genomeRate: 0.0070 * 0.10, hetBand }),
      null,
      `hetBand ${hetBand} is above the diploid threshold, so this genome is not established as `
      + 'carrying one parental contribution and nothing may be inherited from the class')
  }
  console.log('  a genome above the diploid heterozygosity threshold is refused, not inherited from')
}

console.log('uniparentalOrigin: names a parent where dosage cannot, and stays silent where it must')

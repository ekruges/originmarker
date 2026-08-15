// Self-check for the dosage channel. Run: node src/dosageOrigin.check.ts
//
// The thing being guarded against is specific and has happened here twice. Two earlier attempts
// formed a paternal allele share and centred it on the sample's own median. For a one-parent
// informative set the majority of markers sit at the ceiling for every allele frequency under one
// half, so the median IS the ceiling, upward headroom is zero by construction, and every call came
// back the same way. Six one-directional results were produced that way and all six were wrong.
//
//   SO THE TEST THAT MATTERS IS SYMMETRY. The same evidence pointing the other way must produce
//   the mirror-image answer with the same confidence. If one direction is easier to reach than the
//   other, the statistic has a built-in preference and its output is about the statistic.
//
//   AND THE SECOND ONE IS THAT IT MUST NOT BE A CONSTANT. A caller that says the same thing on
//   every input is symmetric, consistent, and useless.
import assert from 'node:assert/strict'
import {
  callDosageOrigin, band, BAND_EXTREME, MIN_MARKERS_DOSAGE, DOSAGE_NOISE, NAMES_A_PARENT,
} from './dosageOrigin.ts'
import type { AB } from './informativity.ts'

const Q = 0.30
const N = 4_000

/** Build a region. `f` returns the sample's B-allele frequency for marker i. */
const region = (parent: AB, n: number, f: (i: number) => number): [AB, number][] =>
  Array.from({ length: n }, (_, i) => [parent, f(i)] as [AB, number])

/** The dosage a marker shows under each truth, given the parent's own allele is A. */
const bothPresent = (i: number) => (i % 10 < Q * 10 ? 0.50 : 0.02)
const knownLost = (i: number) => (i % 10 < Q * 10 ? 0.98 : 0.02)
const otherLost = () => 0.02

// --- 0. THE CLASS IS WITHHELD WHILE THE NULL IS NOT SELF-REFERENCED -----------------------------
//
// A methods review measured this statistic's null on real material: the one-parent centroid null
// median is -0.031 on trophectoderm and -0.023 on blastomeres, i.e. off-centre under NO event, in
// the direction of the parent that was not genotyped, and on TE that offset equals the shift a
// mosaic fraction of 0.117 would produce. Below f = 0.3, half to all detections name the wrong
// parent, and the wrong-sign null tail reaches z = -15.9. So the event is reported and the class
// is not, which is the MoChA precedent: 29% of its events were left unassigned because power to
// detect an imbalance exceeded power to resolve which imbalance it was.
{
  assert.equal(NAMES_A_PARENT, false,
    'flipping this back on requires a self-referenced null, not just an edit here')
  for (const f of [knownLost, otherLost]) {
    const c = callDosageOrigin(region('AA', N, f), Q)
    assert.equal(c.verdict, 'refused', 'an imbalance must not be given a parent')
    assert.ok(c.why.includes('does not name a parent'))
    assert.ok(c.posterior > 0.5, 'and the evidence for the imbalance is still reported')
  }
  // Two copies present is not a parental claim, so it is still stated.
  assert.equal(callDosageOrigin(region('AA', N, bothPresent), Q).verdict, 'both-present')
}

// --- 1. the underlying discrimination still works, which is what makes the refusal a CHOICE ------
//
// Asserted on the posterior rather than the verdict: the likelihood must still separate the three
// states, so that turning naming back on after self-referencing is a threshold change and not a
// rebuild. A statistic that had stopped discriminating would pass the section above for the wrong
// reason.
{
  const known = callDosageOrigin(region('AA', N, knownLost), Q)
  const other = callDosageOrigin(region('AA', N, otherLost), Q)
  assert.ok(known.excluded > N / 4, 'loss of the loaded copy still shows the impossible extreme')
  assert.equal(other.excluded, 0, 'loss of the other copy still shows none of it')
  assert.ok(known.posterior > 0.95 && other.posterior > 0.95)
}

// --- 2. SYMMETRY, which is the whole point -------------------------------------------------------
//
// A BB parent is the same problem with the dosage axis reversed. Flipping every reading must flip
// nothing about the answer, because the orientation is handled in one place. If these disagree,
// one direction has been given headroom the other lacks, which is exactly the failure that
// produced six wrong results here before.
{
  for (const [name, f] of [['both', bothPresent], ['known', knownLost], ['other', otherLost]] as const) {
    const aa = callDosageOrigin(region('AA', N, f as (i: number) => number), Q)
    const bb = callDosageOrigin(region('BB', N, (i) => 1 - (f as (i: number) => number)(i)), Q)
    assert.equal(aa.verdict, bb.verdict, `${name}: the mirrored region must give the same verdict`)
    // Symmetry between the two DIRECTIONS holds and is asserted below, but note it does not
    // protect against a null that is off-centre to begin with, which is why naming is withheld.
    assert.ok(Math.abs(aa.posterior - bb.posterior) < 1e-9,
      `${name}: and the same confidence, got ${aa.posterior} against ${bb.posterior}`)
    assert.equal(aa.excluded, bb.excluded, `${name}: the excluded count must mirror exactly`)
    assert.equal(aa.middle, bb.middle)
  }
}

// --- 3. the two loss directions are equally reachable ---------------------------------------------
//
// Not just mirror-consistent within one parent genotype, but equally easy to reach from a neutral
// starting point. Take the same number of decisive markers and point them each way: the posteriors
// must match. A caller that finds one direction easier will report that direction more often on
// real data and the result will be about the caller.
{
  const decisive = 900
  const towardKnown = callDosageOrigin(
    region('AA', N, (i) => (i < decisive ? 0.98 : 0.02)), Q,
  )
  const towardBoth = callDosageOrigin(
    region('AA', N, (i) => (i < decisive ? 0.50 : 0.02)), Q,
  )
  assert.equal(towardBoth.verdict, 'both-present')
  assert.equal(towardKnown.verdict, 'refused', 'named class withheld, evidence still reported')
  assert.ok(Math.abs(towardKnown.posterior - towardBoth.posterior) < 1e-6,
    'the same weight of evidence must buy the same confidence in either direction, got '
    + `${towardKnown.posterior} against ${towardBoth.posterior}`)
}

// --- 4. not a constant ---------------------------------------------------------------------------
{
  // The two decisive bands together, since neither separates all three states alone: both-present
  // and other-parent-lost share an excluded count of zero and differ only in the middle band.
  const sig = [bothPresent, knownLost, otherLost].map((f) => {
    const c = callDosageOrigin(region('AA', N, f), Q)
    return `${c.excluded}/${c.middle}`
  })
  assert.ok(new Set(sig).size === 3, `the evidence must vary with the truth, saw ${sig.join(' ')}`)
}

// --- 5. THE POINT OF THE CHANNEL: it survives a collapsed call rate -------------------------------
//
// The genotype caller cannot answer a whole-chromosome loss because the event is detected by its
// genotypes failing. Dosage is read whether or not a genotype is emitted, so the same region with
// every genotype removed is unchanged here. Modelled by passing the dosages alone: there is no
// genotype in this function's input at all, which is the design.
{
  const collapsed = callDosageOrigin(region('AA', N, knownLost), Q)
  assert.ok(collapsed.markers === N, 'every dosage is used, genotype or not')
  assert.ok(collapsed.excluded > 0, 'and the evidence survives a collapsed call rate')
}

// --- 6. the impossible bands are what carry the call ----------------------------------------------
{
  // A middle-band dosage cannot come from one copy, so a region full of them rules out both losses.
  const het = callDosageOrigin(region('AA', N, () => 0.50), Q)
  assert.equal(het.verdict, 'both-present')
  assert.equal(het.middle, N)
  assert.equal(het.excluded, 0)

  // An excluded-extreme dosage cannot occur while the loaded parent's copy is present.
  const gone = callDosageOrigin(region('AA', N, () => 0.98), Q)
  assert.equal(gone.excluded, N)
  assert.ok(gone.posterior > 0.99, 'the imbalance is certain even though the class is withheld')
}

// --- 6b. UNRESOLVED DOSAGE IS REFUSED, NOT SCORED -------------------------------------------------
//
// Found on a real blastomere. Chromosome 1 scattered 37.2% of its dosages into no band at all
// against 7.1% over the rest of the same genome, and because a middle-band reading is impossible
// with one copy, the noise that landed there read as both-copies-present at posterior 1.0000.
{
  const noisy = region('AA', N, (i) => 0.20 + ((i * 7) % 60) / 100)
  const scored = callDosageOrigin(noisy, Q)
  assert.notEqual(scored.verdict, 'refused', 'with no background there is nothing to compare to')

  const guarded = callDosageOrigin(noisy, Q, DOSAGE_NOISE, { background: { between: 0.07 } })
  assert.equal(guarded.verdict, 'refused', 'against the sample\'s own genome it must refuse')
  assert.ok(guarded.why.includes('no band at all'))
  assert.ok(guarded.between > 0)

  // A clean region is unaffected by the guard.
  const clean = callDosageOrigin(region('AA', N, knownLost), Q, DOSAGE_NOISE,
    { background: { between: 0.07 } })
  assert.ok(clean.excluded > 0 && clean.between < N / 4, 'a clean region is not caught by the guard')

  // And a sample whose whole genome is noisy is not punished for it: the comparison is relative.
  const tolerant = callDosageOrigin(noisy, Q, DOSAGE_NOISE, { background: { between: 0.30 } })
  assert.notEqual(tolerant.verdict, 'refused')
}

// --- 7. refusals ---------------------------------------------------------------------------------
{
  const few = callDosageOrigin(region('AA', MIN_MARKERS_DOSAGE - 1, knownLost), Q)
  assert.equal(few.verdict, 'refused')
  assert.ok(few.why.includes('noisier'))

  // A parent heterozygous everywhere carries nothing: every hypothesis predicts the same dosage.
  assert.equal(callDosageOrigin(region('AB', N, knownLost), Q).markers, 0)
  assert.equal(callDosageOrigin(region('AB', N, knownLost), Q).verdict, 'refused')

  // Missing dosages are dropped rather than counted as evidence.
  const half = callDosageOrigin(
    Array.from({ length: N }, (_, i) => ['AA', i % 2 ? null : 0.98] as [AB, number | null]), Q,
  )
  assert.equal(half.markers, N / 2)
  assert.ok(half.excluded > 0)
}

// --- 8. noise costs power rather than buying a wrong answer ---------------------------------------
//
// Poor amplification drags an extreme toward the middle. That moves evidence out of the decisive
// bands into the band every hypothesis allows, so it must lower confidence rather than change the
// answer. The unsafe direction is a middle dragged to an extreme, which is why the extreme bands
// are narrow.
{
  const clean = callDosageOrigin(region('AA', N, knownLost), Q)
  const smeared = callDosageOrigin(
    region('AA', N, (i) => (i % 10 < Q * 10 ? (i % 3 === 0 ? 0.72 : 0.98) : 0.02)), Q,
  )
  assert.ok(smeared.posterior <= clean.posterior + 1e-12, 'smearing must not increase confidence')
  assert.ok(smeared.excluded < clean.excluded, 'the decisive count falls, which is the mechanism')
}

// --- 9. the band boundaries are oriented, not absolute --------------------------------------------
{
  assert.equal(band('AA', 0.02), 'own')
  assert.equal(band('BB', 0.98), 'own', 'the parent\'s OWN allele is low only after orientation')
  assert.equal(band('AA', 0.98), 'excluded')
  assert.equal(band('BB', 0.02), 'excluded')
  assert.equal(band('AA', 0.50), 'middle')
  assert.equal(band('BB', 0.50), 'middle', 'the middle band is its own mirror')
  assert.equal(band('AA', 0.25), 'between', 'between the bands is not evidence for anything')
  assert.equal(band('AA', NaN), 'between')
  assert.ok(BAND_EXTREME < 0.35 && DOSAGE_NOISE > 0 && DOSAGE_NOISE < 0.1)
}

console.log('dosageOrigin.check.ts: all assertions passed, including the two loss directions being '
  + 'equally reachable and noise costing power rather than buying an answer')

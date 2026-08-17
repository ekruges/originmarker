// Self-check for the untransmitted-haplotype channel. Run: node src/untransmitted.check.ts
//
// Two properties carry this module and both are easy to lose in a refactor.
//
//   THE MARKER SET IS THE POINT. The gain over the obligate-het channel is not a cleverer
//   statistic, it is that every marker here is informative by construction. A change that lets an
//   ambiguous marker through would keep the code working and quietly give back the advantage.
//
//   BPH IS POSITIVE EVIDENCE AND SPH IS NOT. The two single-homologue bands are also populated by
//   an ordinary euploid genome, so occupancy there means nothing on its own: power 0.134 against
//   euploid. Only the both-homologues band is unreachable without both parental homologues. A
//   version that reported SPH positively would be asserting something its own data cannot support.
import assert from 'node:assert/strict'
import {
  untransmittedPairs, impossibleRate, orientUntransmitted, callMechanism,
  BAND_BPH, BAND_SPH_LOW, MIN_MARKERS_MECHANISM, F80_UNTRANSMITTED, SNR_GAIN_RANGE,
} from './untransmitted.ts'
import type { AB } from './informativity.ts'

type Row = [AB, AB, number | null]

// --- 1. only determined transmissions are kept ----------------------------------------------------
{
  const rows: Row[] = [
    ['AB', 'AA', 0.02],   // parent gave A, withheld B
    ['AB', 'BB', 0.98],   // parent gave B, withheld A
    ['AB', 'AB', 0.50],   // ambiguous: either allele could have come from either parent
    ['AA', 'AA', 0.02],   // parent homozygous: the other channel's territory, not this one
    ['BB', 'AB', 0.50],
    ['AB', 'NC', 0.40],   // no call
    ['AB', 'AA', null],   // no dosage
  ]
  const r = untransmittedPairs(rows)
  assert.equal(r.pairs.length, 2, 'only the two determined transmissions are informative here')
  assert.equal(r.considered, 5, 'every parent-heterozygous marker is counted as considered')
  assert.equal(r.ambiguous, 1, 'and the heterozygous sample is counted as ambiguous, not dropped silently')
  assert.equal(r.pairs[0].untransmitted, 'B', 'a sample reading AA means B was withheld')
  assert.equal(r.pairs[1].untransmitted, 'A')

  // Disjointness from the obligate-het channel, which is what makes this additional evidence
  // rather than a re-reading of the same markers.
  assert.equal(untransmittedPairs([['AA', 'AA', 0.1], ['BB', 'BB', 0.9]]).pairs.length, 0,
    'a parent-homozygous marker must never enter this channel')
}

// --- 2. orientation puts the untransmitted allele high, in both directions ------------------------
//
// The untransmitted allele reaches the sample only through the OTHER parent, so its dosage rising
// is evidence the loaded parent's contribution is short. Both directions must be handled by the
// same expression or the channel acquires the directional null the dosage channel had to remove.
{
  assert.ok(Math.abs(orientUntransmitted({ untransmitted: 'B', baf: 0.9 }) - 0.9) < 1e-12)
  assert.ok(Math.abs(orientUntransmitted({ untransmitted: 'A', baf: 0.1 }) - 0.9) < 1e-12)
  // A mirrored pair orients identically, which is the symmetry the check exists for.
  assert.equal(
    orientUntransmitted({ untransmitted: 'B', baf: 0.25 }),
    orientUntransmitted({ untransmitted: 'A', baf: 0.75 }),
  )
}

// --- 3. the impossible rate is a real per-array error rate ----------------------------------------
{
  const clean = [
    { untransmitted: 'B' as const, baf: 0.03 },
    { untransmitted: 'A' as const, baf: 0.97 },
  ]
  assert.equal(impossibleRate(clean), 0, 'readings agreeing with their own genotype are possible')

  // A sample calling AA whose dosage sits at the B extreme contradicts itself.
  const contradictory = [
    { untransmitted: 'B' as const, baf: 0.95 },
    { untransmitted: 'B' as const, baf: 0.03 },
  ]
  assert.equal(impossibleRate(contradictory), 0.5)
  assert.ok(Number.isNaN(impossibleRate([])))
}

// --- 4. BPH IS CALLED POSITIVELY, not-BPH ONLY BY EXCLUSION -------------------------------------------
{
  const at = (share: number, n: number) => Array.from({ length: n }, () => ({
    // orientUntransmitted returns baf when the untransmitted allele is B, so this places mass
    // directly at the requested share.
    untransmitted: 'B' as const, baf: share,
  }))

  const bph = callMechanism(at(BAND_BPH, 600))
  assert.equal(bph.mechanism, 'BPH')
  assert.ok(bph.atBph > 0.9)
  assert.ok(bph.why.includes('euploid genome does not populate'))

  const sph = callMechanism([...at(BAND_SPH_LOW, 300), ...at(1.0, 300)])
  assert.equal(sph.mechanism, 'not-BPH')
  assert.ok(sph.atBph < 0.25, 'the both-homologues band is empty under a duplicated homologue')
  assert.ok(sph.why.includes('EXCLUSION'),
    'SPH must be reported by exclusion, since its bands are also populated by a euploid genome')
  assert.ok(sph.why.includes('0.134'), 'and the measured power against euploid must travel with it')
  // The class must state what it POOLS. Calling it "SPH" invited a reader to hear "mitotic", and
  // meiosis II without recombination is indistinguishable from mitotic on genotype alone.
  assert.ok(sph.why.includes('POOLS'), 'not-BPH must say it pools MII-without-recombination with mitotic')
  assert.ok(sph.why.includes('distal to the centromere'),
    'and must say why: a recombinant MII would have shown both-homologue tracts distally')

  // A euploid genome puts mass at the same two bands SPH uses, which is exactly why SPH cannot be
  // positive evidence. The caller must not mistake one for the other.
  const euploid = callMechanism([...at(BAND_SPH_LOW, 300), ...at(1.0, 300)])
  assert.equal(euploid.mechanism, sph.mechanism,
    'a euploid genome and an SPH trisomy are indistinguishable on this statistic alone, which is '
    + 'why copy number three must already be established before it is asked')
}

// --- 4b. THE MECHANISM IS NOT ASKED UNLESS COPY NUMBER THREE IS ESTABLISHED -----------------------
//
// Found by running the channel on a real euploid chromosome of a confirmed parent-child pair: it
// returned SPH, by exclusion, confidently and wrongly. A euploid genome populates both
// single-homologue bands, so this statistic answers SPH for every normal chromosome unless the
// question is gated on the copy number first.
{
  const at = (share: number, n: number) => Array.from({ length: n }, () => ({
    untransmitted: 'B' as const, baf: share,
  }))
  const euploidLike = [...at(BAND_SPH_LOW, 300), ...at(1.0, 300)]
  assert.equal(callMechanism(euploidLike).mechanism, 'not-BPH',
    'ungated it answers not-BPH, which is the trap')
  const gated = callMechanism(euploidLike, { copyNumberThree: false })
  assert.equal(gated.mechanism, 'unresolved', 'gated on a euploid chromosome it must not answer')
  assert.ok(gated.why.includes('only exists once'))
  // And the gate does not block a genuine trisomy.
  assert.equal(callMechanism(at(BAND_BPH, 600), { copyNumberThree: true }).mechanism, 'BPH')
}

// --- 5. it is an occupancy statistic, so it needs COUNT rather than precision ----------------------
{
  const noisy = Array.from({ length: 600 }, (_, i) => ({
    untransmitted: 'B' as const,
    // Wide per-marker scatter around the both-homologues band: a per-marker three-band assignment
    // would misassign 13-22% of these, and the occupancy is unmoved.
    baf: BAND_BPH + 0.10 * Math.sin(i * 2.399963),
  }))
  assert.equal(callMechanism(noisy).mechanism, 'BPH',
    'per-marker dispersion must not defeat an occupancy statistic')

  const few = callMechanism(Array.from({ length: MIN_MARKERS_MECHANISM - 1 },
    () => ({ untransmitted: 'B' as const, baf: BAND_BPH })))
  assert.equal(few.mechanism, 'unresolved')
  assert.ok(few.why.includes('occupancy'))
}

// --- 6. the inherited floors say what this channel is FOR ------------------------------------------
//
// The blastomere entry is the reason the module exists: the obligate-het channel has no floor
// there at all, in either configuration, and this one has a number. Still too high to call, so it
// converts a material impossibility into a measurable quantity rather than into an answer.
{
  assert.ok(Number.isFinite(F80_UNTRANSMITTED.twoParents.blastomere),
    'this is the only channel that gives a blastomere a defined floor')
  assert.ok(F80_UNTRANSMITTED.twoParents.blastomere > 0.35,
    'and it is honest that the floor is still above the callable bound')
  assert.ok(F80_UNTRANSMITTED.twoParents.trophectoderm < 0.35,
    'while a trophectoderm biopsy with two parents is callable through it')
  assert.ok(SNR_GAIN_RANGE[0] > 1.3 && SNR_GAIN_RANGE[1] < 2.0)
}

console.log('untransmitted.check.ts: all assertions passed, including only determined '
  + 'transmissions entering the channel, and not-BPH reported by exclusion while stating what it '
  + 'pools')

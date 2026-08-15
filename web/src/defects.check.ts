// Self-check for the defect display mapping. Run: node src/defects.check.ts
//
// This file exists because of what the module it tests can get wrong. Showing a confident parent
// for a region that never had one, or naming the WRONG parent, is the most damaging thing this UI
// can do: the reader cannot tell a mislabelled call from a correct one by looking at it.
//
//   THE LOADED PARENT IS NOT ALWAYS THE FATHER. A single-parent verdict says whose copy is
//   missing relative to whoever was loaded. Reading 'known-parent-lost' as paternal is only
//   correct when the loaded parent WAS the father, and the display was doing exactly that.
//
//   TWO PARENTS OUTRANK ONE. Where dosage from both parents already named a side, a single-parent
//   call must not overwrite it; it fills a gap rather than competing.
import assert from 'node:assert/strict'
import { defectsFrom } from './defects.ts'
import type { Segment } from './segments.ts'
import type { GainAnnotation } from './parentage.ts'
import type { StageCall } from './stage.ts'

const seg = (chrom: string, startBp: number, endBp: number): Segment => ({
  chrom,
  startBp,
  endBp,
  kind: 'copy-loss',
  markers: 900,
  rate: 0.42,
  refined: false,
} as unknown as Segment)

// The key both sides join on. Written the way the pipeline writes it, so a change to that format
// breaks this test rather than silently unjoining the two halves of the display.
const whereOf = (s: Segment, start: number, end: number) =>
  `chr${s.chrom} ${(start / 1e6).toFixed(1)}-${(end / 1e6).toFixed(1)}Mb`

const one = (where: string, verdict: string) => ({
  where, verdict, posterior: 0.994, markers: 1_420, exclusive: 5_130, why: 'test',
})

// --- 1. a single-parent verdict is read against the parent that was actually loaded ------------
{
  const s = seg('21', 14_300_000, 46_700_000)
  const w = whereOf(s, 14_300_000, 46_700_000)

  const father = defectsFrom([s], [], [], [one(w, 'known-parent-lost')], 'paternal')
  assert.equal(father[0].origin, 'paternal', 'the loaded father lost his copy')

  // The same verdict with the MOTHER loaded is the opposite answer. This parameter was always
  // here and always correct; the caller passed the literal 'paternal' to it, because the result
  // did not carry its own role and so there was nothing truthful to pass. The result carries it
  // now, and this asserts the direction the caller depends on.
  const mother = defectsFrom([s], [], [], [one(w, 'known-parent-lost')], 'maternal')
  assert.equal(mother[0].origin, 'maternal', 'the loaded mother lost hers; it is not paternal')

  // And 'other-parent-lost' names the parent who was NOT loaded, in both directions.
  assert.equal(defectsFrom([s], [], [], [one(w, 'other-parent-lost')], 'paternal')[0].origin,
    'maternal')
  assert.equal(defectsFrom([s], [], [], [one(w, 'other-parent-lost')], 'maternal')[0].origin,
    'paternal')

  // A verdict that named nobody must not name anybody.
  for (const v of ['both-present', 'refused']) {
    const d = defectsFrom([s], [], [], [one(w, v)], 'paternal')[0]
    assert.equal(d.origin, 'unclear', `${v} is not a parent`)
    assert.equal(d.basis, undefined, 'and carries no basis, since nothing was established')
  }
}

// --- 2. two-parent dosage is not overwritten by a single-parent call ---------------------------
{
  const s = seg('16', 1_000_000, 9_000_000)
  const w = whereOf(s, 1_000_000, 9_000_000)
  const ann = { where: w, origin: 'maternal copy missing', why: 'from dosage', called: true }

  const d = defectsFrom([s], [], [ann as unknown as GainAnnotation],
    [one(w, 'known-parent-lost')], 'paternal')[0]
  assert.equal(d.origin, 'maternal', 'both parents outrank one; the single-parent call must yield')
  assert.equal(d.basis, 'two-parent')
}

// --- 3. a run with no parental genotype still reports the position ------------------------------
//
// The coordinates come from intensity and owe nothing to a parent, so withholding them because
// origin is unknown would hide a real finding behind a missing field.
{
  const s = seg('6', 39_302, 294_904)
  const d = defectsFrom([s], [], [], [], 'paternal')[0]
  assert.equal(d.origin, 'unclear')
  assert.equal(d.locus, 'chr6:39,302-294,904', 'position is unconditional and conventionally written')
  assert.ok(d.why.includes('not determined'))
  assert.equal(d.exclusive, undefined, 'and no evidence count is invented for it')
}

// --- 4. the Mendelian evidence reaches the display ----------------------------------------------
//
// The exclusive-marker count is the one number on a single-parent run that a reader can check
// without trusting the model: dropout removes alleles and cannot invent one. It was being
// computed and discarded.
{
  const s = seg('21', 14_300_000, 46_700_000)
  const w = whereOf(s, 14_300_000, 46_700_000)
  const d = defectsFrom([s], [], [], [one(w, 'known-parent-lost')], 'paternal')[0]
  assert.equal(d.exclusive, 5_130)
  assert.equal(d.informative, 1_420)
  assert.equal(d.posterior, 0.994)
  assert.equal(d.basis, 'one-parent')
}

// --- 4b. dosage fills what genotypes cannot reach, and never overrides them ----------------------
//
// The two channels answer the same question from different measurements and are not
// interchangeable. A whole chromosome is detected by its genotype call rate collapsing, so on
// those events the genotype channel has no evidence left and dosage is the only one that can
// speak. Where genotypes DID answer, dosage must not change it.
{
  const s = seg('21', 10_873_592, 48_088_571)
  const w = whereOf(s, 10_873_592, 48_088_571)
  const dose = [{ where: 'chr21', verdict: 'known-parent-lost', shift: 0.118, z: 6.4,
    impliedF: 0.382, window: 9_100, why: 'from dosage' }]

  // Genotypes silent: dosage answers, and says so.
  const filled = defectsFrom([s], [], [], [], 'paternal', undefined, dose)[0]
  assert.equal(filled.origin, 'paternal')
  assert.equal(filled.basis, 'dosage')
  assert.equal(filled.channel, 'dosage')
  assert.equal(filled.z, 6.4)
  assert.equal(filled.impliedF, 0.382)
  assert.equal(filled.shift, 0.118)

  // Genotypes answered: dosage must not touch it.
  const geno = defectsFrom([s], [], [], [one(w, 'other-parent-lost')], 'paternal', undefined, dose)[0]
  assert.equal(geno.origin, 'maternal', 'the genotype answer stands')
  assert.equal(geno.basis, 'one-parent')
  assert.equal(geno.channel, 'genotype')

  // The loaded parent decides the direction here too.
  assert.equal(defectsFrom([s], [], [], [], 'maternal', undefined, dose)[0].origin, 'maternal')

  // Every verdict short of naming a parent names nobody. 'imbalance-unassigned' is the important
  // one: an event IS present and the class is deliberately withheld, so it must not leak a parent.
  for (const v of ['imbalance-unassigned', 'no-imbalance', 'not-evaluable', 'array-excluded']) {
    const r = defectsFrom([s], [], [], [], 'paternal', undefined,
      [{ ...dose[0], verdict: v }])[0]
    assert.equal(r.origin, 'unclear', `${v} must not name a parent`)
    assert.equal(r.basis, undefined)
  }
}

// --- 5. the material is carried with the call, including when it is not usable ------------------
{
  const s = seg('21', 14_300_000, 46_700_000)
  const stage = {
    stage: 'blastomere', dropout: 0.308, basis: 'heterozygosity-shortfall',
    templates: '2', markerFloor: 200, caveat: 'c', why: 'w',
  } as StageCall
  const d = defectsFrom([s], [], [], [], 'paternal', stage)[0]
  assert.equal(d.stage, 'blastomere')
  assert.equal(d.dropout, 0.308)
  assert.equal(d.dropoutBasis, 'heterozygosity-shortfall')

  // A failed array has no dropout figure. It must arrive as absent rather than as NaN, which
  // renders as "NaN" in a chip and reads to a reader as a measurement that came out strange.
  const dead = { ...stage, stage: 'failed', dropout: NaN, basis: 'none' } as StageCall
  const f = defectsFrom([s], [], [], [], 'paternal', dead)[0]
  assert.equal(f.stage, 'failed')
  assert.equal(f.dropout, undefined)
  assert.equal(f.dropoutBasis, undefined)
}

console.log('defects.check.ts: all assertions passed, including a single-parent verdict read '
  + 'against the parent actually loaded rather than assumed paternal')

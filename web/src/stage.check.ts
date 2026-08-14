// Self-check for stage inference. Run: node src/stage.check.ts
//
// Two properties matter more than the boundaries themselves.
//
//   QUALITY BEFORE PLOIDY. A failed amplification drives heterozygosity toward zero, so without a
//   call-rate gate it is classified as haploid and handed a confident parameter set instead of
//   being rejected. The two cannot be separated by heterozygosity, so quality decides first.
//
//   HAPLOID IS A DIFFERENT AXIS, BUT NOT A UNIFORM ONE. A PB2, pronucleus or sperm carries one
//   chromatid and has true heterozygosity of zero. A FIRST polar body carries a dyad, and distal to
//   each crossover its sister chromatids differ, so roughly 44% of its genome is genuinely
//   heterozygous. An earlier version of this module asserted all haploid heterozygosity was error.
//
//   EVERY SHORTFALL ESTIMATE CARRIES ITS CONFOUNDS. Consanguinity, LOH, UPD and ancestry all
//   depress heterozygosity and are absorbed into the number, so the number must not travel alone.
import assert from 'node:assert/strict'
import {
  inferStage, locus, dropoutFromReplicates, BULK_HETEROZYGOSITY, HAPLOID_MAX_HET, QC_CALL_FLOOR,
} from './stage.ts'

const p = (hetRate: number, callRate = 0.95) => ({ hetRate, callRate })

// --- 1. each stage is reached by its measured heterozygosity -------------------------------------
{
  assert.equal(inferStage(p(0.168)).stage, 'bulk')
  assert.equal(inferStage(p(0.150)).stage, 'trophectoderm')
  assert.equal(inferStage(p(0.135)).stage, 'single-cell')
  assert.equal(inferStage(p(0.116)).stage, 'blastomere')
}

// --- 2. haploid is separated from the diploid stages ---------------------------------------------
//
// Measured haploid products run 0.002 to 0.10, which spans both a single chromatid at the drop-in
// floor and a first polar body with real heterozygosity. All must come back haploid rather than
// blastomere, and none may carry a large dropout parameter.
{
  for (const h of [0.002, 0.021, 0.055, 0.099]) {
    const s = inferStage(p(h))
    assert.equal(s.stage, 'haploid', `${h} heterozygous must read haploid, got ${s.stage}`)
    assert.ok(s.dropout < 0.05,
      'a genome with no heterozygotes has no heterozygote to drop, so dropout must stay small')
    assert.ok(s.why.includes('one genome'))
  }
  // And the boundary sits in the gap rather than on either population.
  assert.ok(HAPLOID_MAX_HET > 0.10 && HAPLOID_MAX_HET < 0.116)
}

// --- 3. dropout rises monotonically as material gets scarcer --------------------------------------
{
  const stages = [0.168, 0.150, 0.135, 0.116].map((h) => inferStage(p(h)).dropout)
  for (let i = 1; i < stages.length; i += 1) {
    assert.ok(stages[i] > stages[i - 1],
      `dropout must rise as template count falls, got ${stages.join(' ')}`)
  }
}

// --- 4. QUALITY IS GATED BEFORE PLOIDY ------------------------------------------------------------
//
// The failure an audit found in the first version: near-total dropout drives heterozygosity toward
// zero, so a failed amplification was classified as HAPLOID and handed a confident parameter set
// rather than being rejected. The two are indistinguishable by heterozygosity, so call rate must
// decide first, and the gate must fire even when the heterozygosity looks haploid.
{
  const dead = inferStage(p(0.004, 0.22))
  assert.equal(dead.stage, 'failed', 'a failed amplification must not be called haploid')
  assert.ok(!Number.isFinite(dead.dropout), 'and must not carry a usable dropout parameter')
  assert.equal(dead.basis, 'none')

  // The gate is on call rate, so a plausible heterozygosity does not rescue a dead array.
  assert.equal(inferStage(p(0.30, 0.31)).stage, 'failed')
  // And a good call rate with the same heterozygosity IS classified.
  assert.notEqual(inferStage(p(0.004, 0.95)).stage, 'failed')
  assert.equal(inferStage(p(NaN)).stage, 'unknown')
  assert.ok(QC_CALL_FLOOR > 0.2 && QC_CALL_FLOOR < 0.7)
}

// --- 4b. FIRST polar bodies are not homozygous, and are not treated as error ----------------------
//
// A PB1 carries a dyad, and distal to each crossover its two sister chromatids differ, so about 44%
// of its genome is genuinely heterozygous at an expected h near 0.074. The first version asserted
// all haploid heterozygosity was error, which mis-parameterises exactly this material.
{
  const pb1 = inferStage(p(0.074))
  assert.equal(pb1.stage, 'haploid')
  assert.ok(pb1.caveat.includes('FIRST polar body') || pb1.templates.includes('chromatids'),
    'a sample sitting where a PB1 sits must say so rather than calling its heterozygosity error')

  const pb2 = inferStage(p(0.008))
  assert.equal(pb2.stage, 'haploid')
  assert.equal(pb2.templates, '1 chromatid', 'a single chromatid is a different object from a dyad')
}

// --- 4c. every shortfall-based call carries its confounds -----------------------------------------
//
// The estimate absorbs consanguinity, LOH, UPD and ancestry, and an East Asian sample against this
// European anchor shifts by enough to change the stage. If that is not stated with the number, the
// number reads as a measurement.
{
  for (const h of [0.168, 0.150, 0.135, 0.116]) {
    const c = inferStage(p(h))
    assert.equal(c.basis, 'heterozygosity-shortfall')
    assert.ok(c.caveat.includes('ancestry'), 'ancestry bias must travel with the estimate')
    assert.ok(c.caveat.includes('replicate'), 'and the better estimator must be named')
  }
}

// --- 4d. the replicate estimator inverts correctly, and is a lower bound --------------------------
{
  for (const d of [0.05, 0.199, 0.308]) {
    const phi = (2 * d) / (1 + d)
    assert.ok(Math.abs(dropoutFromReplicates(phi) - d) < 1e-9,
      `phi = 2d/(1+d) must invert to d, failed at ${d}`)
  }
  // Correlated failure lowers the discordant fraction, so the recovered value must fall too:
  // the estimator is a lower bound rather than unbiased.
  const d = 0.199
  const phiIdeal = (2 * d) / (1 + d)
  const pBoth = d * d + 0.25 * d * (1 - d)
  const pOne = 2 * (d - pBoth)
  const phiCorr = pOne / (pOne + (1 - 2 * d + pBoth))
  assert.ok(dropoutFromReplicates(phiCorr) < dropoutFromReplicates(phiIdeal),
    'shared failures must make the estimate an underestimate, not an overestimate')
}

// --- 5. an unusually clean array of a lossy stage is not credited with bulk amplification ---------
//
// The floor exists so that one good blastomere does not get a bulk dropout parameter and with it a
// confidence its template count cannot support.
{
  const s = inferStage(p(0.120))
  assert.equal(s.stage, 'blastomere')
  assert.ok(s.dropout > 0.15, `a blastomere must not be credited with bulk dropout, got ${s.dropout}`)
}

// --- 6. the marker floor varies by stage, which was the point ------------------------------------
{
  assert.equal(inferStage(p(0.168)).markerFloor, 100, 'bulk needs fewer markers')
  assert.equal(inferStage(p(0.116)).markerFloor, 200, 'a blastomere needs more')
  assert.ok(inferStage(p(0.150)).markerFloor < inferStage(p(0.116)).markerFloor)
}

// --- 7. template counts are reported, since they are the reason for the rest ----------------------
{
  assert.ok(inferStage(p(0.168)).templates.includes('10^6'))
  assert.equal(inferStage(p(0.116)).templates, '2')
  assert.equal(inferStage(p(0.021)).templates, '1 chromatid',
    'a PB2, pronucleus or sperm is one chromatid, not loosely one genome')
}

// --- 8. locus formatting is the conventional one --------------------------------------------------
{
  assert.equal(locus('6', 39302, 294904), 'chr6:39,302-294,904')
  assert.equal(locus('X', 1_000_000, 2_500_000), 'chrX:1,000,000-2,500,000')
  // Fractional base positions cannot exist; they must not reach the string.
  assert.equal(locus('1', 100.6, 200.2), 'chr1:101-200')
}

// --- 9. the bulk figure is the one every shortfall is measured against ----------------------------
{
  assert.equal(BULK_HETEROZYGOSITY, 0.168)
  const s = inferStage(p(BULK_HETEROZYGOSITY))
  assert.ok(s.dropout < 0.03, 'a sample at the bulk rate has essentially no dropout')
}

console.log('stage.check.ts: all assertions passed, including quality gated before ploidy, '
  + 'first polar bodies not treated as error, and confounds travelling with every estimate')

// Self-check for stage inference. Run: node src/stage.check.ts
//
// Two properties matter more than the boundaries themselves.
//
//   HAPLOID IS A DIFFERENT AXIS. A polar body carries one genome and is homozygous by
//   construction, so on heterozygosity alone it looks like a catastrophically dropped-out diploid.
//   Reading it as one would attach a nonsensical dropout parameter to it, so it must be separated
//   before any diploid stage is considered.
//
//   A FAILED ARRAY MUST NOT GET A CONFIDENT STAGE. An array that failed outright can show any
//   heterozygosity at all, and calling a stage from it attaches a precise dropout to noise.
import assert from 'node:assert/strict'
import { inferStage, locus, BULK_HETEROZYGOSITY, HAPLOID_MAX_HET } from './stage.ts'

const p = (hetRate: number, callRate = 0.95) => ({ hetRate, callRate })

// --- 1. each stage is reached by its measured heterozygosity -------------------------------------
{
  assert.equal(inferStage(p(0.168)).stage, 'bulk')
  assert.equal(inferStage(p(0.150)).stage, 'trophectoderm')
  assert.equal(inferStage(p(0.135)).stage, 'single-cell')
  assert.equal(inferStage(p(0.116)).stage, 'blastomere')
}

// --- 2. haploid is separated FIRST, and not read as a lossy diploid -------------------------------
//
// Measured haploid products run 0.002 to 0.10. Every one must come back haploid rather than
// blastomere, and must not carry a large dropout parameter.
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

// --- 4. a failed array gets the conservative answer, not a confident one --------------------------
{
  const s = inferStage(p(0.30, 0.31))
  assert.equal(s.stage, 'unknown', 'a 31% call rate cannot support a stage call')
  assert.equal(s.dropout, 0.308, 'and must assume the worst rather than the best')
  assert.ok(s.why.includes('too low'))
  assert.equal(inferStage(p(NaN)).stage, 'unknown')
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
  assert.equal(inferStage(p(0.021)).templates, '1 genome')
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

console.log('stage.check.ts: all assertions passed, including haploid separated before diploid '
  + 'stages and a failed array refusing a confident call')

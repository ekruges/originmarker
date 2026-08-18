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
  BAND_DIPLOID_CERTAIN, MAX_DIPLOID_HET,
  BULK_MAX_BAF_SD, stageFacts, amplificationWord,
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

// --- 4a. THE LADDER IS CLOSED AT THE TOP AS WELL AS THE BOTTOM ----------------------------------
//
// Found on real laboratory arrays. The ladder was open above: any value at or over the bulk
// boundary matched bulk, so an array reading 52% heterozygous at a 55% call rate was classified as
// bulk genomic DNA and given a dropout of 0.008, the most confident parameter set here. That is
// the same failure as calling a dead array haploid, at the other end and with worse consequences,
// since every downstream likelihood then treats a mixture as pristine material.
{
  for (const h of [0.30, 0.52, 0.56, 0.90]) {
    const s = inferStage(p(h, 0.55))
    assert.equal(s.stage, 'failed', `${h} heterozygous is not a genome, got ${s.stage}`)
    assert.ok(!Number.isFinite(s.dropout))
    assert.ok(s.why.includes('ceiling'))
  }
  // The ceiling must clear what a real diploid can reach: its own rate plus drop-in at the
  // highest rate measured here, 0.168 + (1 - 0.168) * 0.0525 = 0.212.
  const maxReal = BULK_HETEROZYGOSITY + (1 - BULK_HETEROZYGOSITY) * 0.0525
  assert.ok(MAX_DIPLOID_HET > maxReal,
    `a diploid with maximal drop-in reaches ${maxReal.toFixed(3)} and must not be rejected`)
  assert.equal(inferStage(p(maxReal, 0.95)).stage, 'bulk')
  // And a call rate above the floor does not rescue an impossible heterozygosity.
  assert.equal(inferStage(p(0.52, 0.99)).stage, 'failed')
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

// --- 4e. the BAF band vetoes a haploid call, and only vetoes ------------------------------------
//
// Intensity sees a heterozygous locus whether or not the caller dropped it, so a band at the
// diploid line means two genomes however few heterozygotes survived. The converse does NOT hold on
// amplified material: 91 of 120 measured non-haploid arrays fall inside the haploid band range,
// so a low band must not be allowed to confirm one genome.
{
  const h = 0.04
  assert.equal(inferStage({ hetRate: h, callRate: 0.95 }).stage, 'haploid',
    'without band evidence this is the existing behaviour')

  const vetoed = inferStage({ hetRate: h, callRate: 0.95, hetBand: 0.22 })
  assert.notEqual(vetoed.stage, 'haploid', 'a diploid-range BAF band must veto the haploid call')
  assert.ok(vetoed.dropout > 0.5, 'and the sample it produces instead is barely powered')
  assert.ok(vetoed.why.includes('one template'))

  // A low band is uninformative on WGA material, so it changes nothing in the other direction.
  assert.equal(inferStage({ hetRate: h, callRate: 0.95, hetBand: 0.03 }).stage, 'haploid')
  assert.equal(inferStage({ hetRate: h, callRate: 0.95, hetBand: NaN }).stage, 'haploid')

  // The veto never reclassifies material that was already diploid by heterozygosity.
  assert.equal(inferStage({ hetRate: 0.150, callRate: 0.95, hetBand: 0.30 }).stage, 'trophectoderm')
  // Quality still outranks it: a dead array is failed, not rescued into a stage by its band.
  assert.equal(inferStage({ hetRate: h, callRate: 0.22, hetBand: 0.30 }).stage, 'failed')
  assert.ok(BAND_DIPLOID_CERTAIN >= 0.15)
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

// --- THE STAGE LADDER NEEDS BOTH AXES, NOT JUST HETEROZYGOSITY -----------------------------------
//
// Heterozygosity says how many heterozygotes SURVIVED, which tracks template count. BAF spread at
// heterozygous calls says how far the survivors SCATTERED, which tracks amplification. Every
// constant keyed to stage encodes the second: drift alone spans 0.0011 to 0.0193 across these
// classes. Measured over 877 arrays, the ladder on heterozygosity alone put 176 amplified arrays
// on the BULK rung, median scatter 0.232 against unamplified DNA's 0.088, and handed them the
// tightest drift constant in the module. Adding the second axis leaves 18 on that rung at 0.0772.
{
  // Bulk-level heterozygosity. Without the spread, this is bulk, as it always was.
  const noSd = inferStage({ hetRate: 0.168, callRate: 0.97 })
  assert.equal(noSd.stage, 'bulk', 'absent a spread the ladder behaves exactly as before')

  // Same heterozygosity, unamplified scatter: still bulk.
  const clean = inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.088 })
  assert.equal(clean.stage, 'bulk', 'genuine genomic DNA keeps the rung it earned')

  // Same heterozygosity, AMPLIFIED scatter: demoted, not rejected.
  const amp = inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.232 })
  assert.notEqual(amp.stage, 'bulk',
    'an array scattering like amplified material must not be handed the bulk parameter set on the '
    + 'strength of having retained its heterozygotes')
  assert.equal(amp.stage, 'trophectoderm',
    'and it is DEMOTED by heterozygosity rather than discarded: this one kept the most template, '
    + 'so it lands on the highest amplified rung')

  // The boundary itself, and that it is the one already carried for excluding an array outright.
  assert.equal(BULK_MAX_BAF_SD, 0.11)
  assert.equal(inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.109 }).stage, 'bulk')
  assert.notEqual(inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.111 }).stage, 'bulk')

  // It is a dial, and it moves the answer.
  assert.equal(inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.232 },
    { bulkMaxBafSd: 0.5 }).stage, 'bulk', 'a loosened dial restores the old behaviour')

  // THE EXPLANATION MUST BE THE REASON THE RUNG WAS REACHED, not a reason that argued against it.
  // Before this, a demoted sample was told "17.2% heterozygous ... which is where trophectoderm
  // material sits", crediting heterozygosity for a rung heterozygosity had said bulk about. That is
  // a plausible wrong answer, which is the thing this project treats as worse than an error.
  assert.ok(amp.why.includes('would have read as bulk'),
    `a demoted sample must say what it would otherwise have been: ${amp.why}`)
  assert.ok(amp.why.includes('spreads by'), 'and the measurement that demoted it')
  assert.ok(amp.why.includes('placed at trophectoderm by its heterozygosity instead'),
    'and how the new rung was chosen')
  assert.ok(!clean.why.includes('would have read as'),
    'a sample that was never demoted must not carry demotion language')

  // THE SECOND AXIS TOUCHES ONLY THE TOP RUNG. A blastomere is already placed by heterozygosity
  // and amplified scatter is expected there, so nothing below bulk may move.
  for (const h of [0.150, 0.130, 0.110]) {
    const before = inferStage({ hetRate: h, callRate: 0.9 }).stage
    const after = inferStage({ hetRate: h, callRate: 0.9, hetBafSd: 0.25 }).stage
    assert.equal(before, after, `the ${before} rung must not move: it never claimed to be unamplified`)
  }
}

// --- THE RUN PANEL SHOWS BOTH AXES, AND CANNOT DISAGREE WITH THE CALL ----------------------------
//
// The stage sets the drift floor, the variance inflation, the detection floors and the dropout that
// parameterises every directional call below it, across a seventeen-fold range. A reader who cannot
// see which stage was inferred cannot judge any number underneath it, so it is the first thing in
// a run rather than a line in a log.
{
  const amp = inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.232 })
  const f = stageFacts(amp, { hetRate: 0.168, hetBafSd: 0.232 })
  assert.equal(f.stage, 'trophectoderm')
  assert.equal(f.hetRate, 0.168, 'the survived axis')
  assert.equal(f.hetBafSd, 0.232, 'and the scattered axis, separately')
  assert.equal(f.amplified, true)
  assert.equal(f.demotedFromBulk, true, 'the panel must say when heterozygosity was overruled')

  const clean = inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.088 })
  const g = stageFacts(clean, { hetRate: 0.168, hetBafSd: 0.088 })
  assert.equal(g.stage, 'bulk')
  assert.equal(g.amplified, false)
  assert.equal(g.demotedFromBulk, false)
  assert.ok(amplificationWord(g).includes('unamplified'))
  assert.ok(amplificationWord(f).includes('amplified'))

  // DEMOTION IS READ FROM THE CALL'S OWN REASON, not recomputed, so the panel and the call cannot
  // drift apart and tell a reader two different stories about the same sample.
  assert.equal(f.demotedFromBulk, amp.why.includes('would have read as bulk'))

  // A file with no allele fractions must say the axis is absent rather than imply unamplified,
  // since an absent measurement and a low one lead to opposite conclusions.
  const blind = stageFacts(inferStage({ hetRate: 0.168, callRate: 0.97 }), { hetRate: 0.168 })
  assert.equal(blind.amplified, null)
  assert.ok(Number.isNaN(blind.hetBafSd))
  assert.ok(amplificationWord(blind).includes('not measured'))
}

console.log('stage.check.ts: all assertions passed, including quality gated before ploidy, '
  + 'first polar bodies not treated as error, and confounds travelling with every estimate')

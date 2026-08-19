// Self-check for the tunable thresholds. Run: node cli/om.check.ts
//
// The command line exists to move thresholds the web tool does not expose. That creates one
// danger worth testing for directly.
//
//   OMITTING EVERY DIAL MUST BE THE SHIPPED CONFIGURATION. If the optional arguments changed any
//   default, the browser tool and the command line would quietly disagree, and every validation
//   figure in the audit would apply to neither. The web tool passes nothing, so "nothing" has to
//   mean exactly what the constants say.
//
//   A DIAL THAT DOES NOTHING IS WORSE THAN NO DIAL. The first version of the CLI advertised
//   --max-region-het and silently ignored it, because that threshold lived inside the module and
//   was not a parameter. The flag parsed, the run completed, and the number came out unchanged.
//   Each one is asserted to actually move the answer.
import assert from 'node:assert/strict'
import {
  inferStage, QC_CALL_FLOOR, MAX_DIPLOID_HET, HAPLOID_MAX_HET, BULK_HETEROZYGOSITY,
  BAND_DIPLOID_CERTAIN,
} from '../web/src/stage.ts'
import {
  callOneParentOrigin, MIN_MARKERS, MAX_REGION_HET, CALL_POSTERIOR, DEFAULT_Q, GENOTYPE_ERROR,
  DROP_IN,
} from '../web/src/oneParentOrigin.ts'
import type { AB } from '../web/src/informativity.ts'

const p = (hetRate: number, callRate = 0.95, hetBand?: number) => ({ hetRate, callRate, hetBand })

// --- 1. no options is the shipped configuration, on every branch --------------------------------
{
  const explicit = {
    callFloor: QC_CALL_FLOOR,
    maxDiploidHet: MAX_DIPLOID_HET,
    haploidMaxHet: HAPLOID_MAX_HET,
    bulkHeterozygosity: BULK_HETEROZYGOSITY,
    bandDiploidCertain: BAND_DIPLOID_CERTAIN,
  }
  // One case per branch: not-a-genome, failed, haploid, PB1, band veto, and each diploid rung.
  for (const c of [
    p(0.52, 0.55), p(0.004, 0.22), p(0.004), p(0.074), p(0.04, 0.95, 0.22),
    p(0.168), p(0.150), p(0.135), p(0.116), p(NaN),
  ]) {
    assert.deepEqual(inferStage(c), inferStage(c, explicit),
      `passing the shipped constants explicitly must equal passing nothing, at het ${c.hetRate}`)
    assert.deepEqual(inferStage(c), inferStage(c, {}))
  }
}

const pairs = (n: number, f: (i: number) => [AB, AB]): [AB, AB][] =>
  Array.from({ length: n }, (_, i) => f(i))

{
  const cases: [AB, AB][][] = [
    pairs(2_000, (i) => ['AA', i % 10 < 3 ? 'AB' : 'AA']),
    pairs(2_000, (i) => ['AA', i % 10 < 3 ? 'BB' : 'AA']),
    pairs(2_000, () => ['AA', 'AA']),
    pairs(20, () => ['AA', 'AB']),
  ]
  for (const c of cases) {
    assert.deepEqual(
      callOneParentOrigin(c, 0.20),
      callOneParentOrigin(c, 0.20, DEFAULT_Q, GENOTYPE_ERROR, DROP_IN, {
        minMarkers: MIN_MARKERS, maxRegionHet: MAX_REGION_HET, callPosterior: CALL_POSTERIOR,
      }),
      'the shipped constants passed explicitly must equal passing nothing',
    )
    assert.deepEqual(callOneParentOrigin(c, 0.20), callOneParentOrigin(c, 0.20, DEFAULT_Q,
      GENOTYPE_ERROR, DROP_IN, {}))
  }
}

// --- 2. every dial actually moves the answer ----------------------------------------------------
{
  // A dead array is failed by default; lowering the floor lets it through.
  assert.equal(inferStage(p(0.004, 0.22)).stage, 'failed')
  assert.notEqual(inferStage(p(0.004, 0.22), { callFloor: 0.10 }).stage, 'failed')

  // A mixture is failed by default; raising the ceiling admits it as a genome.
  assert.equal(inferStage(p(0.52, 0.55)).stage, 'failed')
  assert.equal(inferStage(p(0.52, 0.55), { maxDiploidHet: 0.70, callFloor: 0.10 }).stage, 'bulk')

  // The haploid boundary moves in both directions.
  assert.equal(inferStage(p(0.099)).stage, 'haploid')
  assert.notEqual(inferStage(p(0.099), { haploidMaxHet: 0.05 }).stage, 'haploid')

  // The anchor changes the dropout the shortfall implies.
  assert.ok(inferStage(p(0.135), { bulkHeterozygosity: 0.30 }).dropout
    > inferStage(p(0.135)).dropout, 'a higher anchor implies more dropout for the same sample')

  // The band veto can be switched off by raising it out of reach.
  assert.notEqual(inferStage(p(0.04, 0.95, 0.22)).stage, 'haploid')
  assert.equal(inferStage(p(0.04, 0.95, 0.22), { bandDiploidCertain: 0.99 }).stage, 'haploid')
}

{
  // The marker floor.
  const few = pairs(20, () => ['AA', 'BB'])
  assert.equal(callOneParentOrigin(few, 0.20).verdict, 'refused')
  assert.notEqual(callOneParentOrigin(few, 0.20, DEFAULT_Q, GENOTYPE_ERROR, DROP_IN,
    { minMarkers: 10 }).verdict, 'refused')

  // The region heterozygosity ceiling. This is the one that silently did nothing.
  const junk = pairs(2_000, (i) => ['AA', i % 10 < 6 ? 'AB' : 'AA'])
  assert.equal(callOneParentOrigin(junk, 0.30).verdict, 'refused')
  assert.notEqual(callOneParentOrigin(junk, 0.30, DEFAULT_Q, GENOTYPE_ERROR, DROP_IN,
    { maxRegionHet: 0.95 }).verdict, 'refused',
  'raising the ceiling must actually admit the region, or the flag is decoration')

  // The posterior a call must reach. Asserted against the posterior actually achieved rather
  // than against 1.0: with two thousand clean markers the posterior IS 1.0 in floating point,
  // so 1.0 is a reachable bar and testing against it proves nothing.
  const clear = pairs(2_000, (i) => ['AA', i % 10 < 3 ? 'BB' : 'AA'])
  const got = callOneParentOrigin(clear, 0.20)
  assert.notEqual(got.verdict, 'refused')
  // THE GATE AND THE DISPLAYED NUMBER ARE ON DIFFERENT SCALES, DELIBERATELY, and this test has to
  // use the gate's own. `callPosterior` asks whether the MODEL is sure enough to name a hypothesis,
  // and the model is decisive over hundreds of Mendelian markers: it reaches 1.0 in floating point.
  // The reported confidence asks something else, how likely the CALL is right once model risk is
  // included, and is capped at 1 - 0.206 by what the validation supports. Testing the gate against
  // the reported number would compare a raw likelihood against a bounded one and refuse everything.
  assert.equal(callOneParentOrigin(clear, 0.20, DEFAULT_Q, GENOTYPE_ERROR, DROP_IN,
    { callPosterior: 1 + 1e-9 }).verdict, 'refused',
  'a bar above any reachable likelihood must refuse the call')
  assert.notEqual(callOneParentOrigin(clear, 0.20, DEFAULT_Q, GENOTYPE_ERROR, DROP_IN,
    { callPosterior: 0.99 }).verdict, 'refused',
  'and a bar the likelihood clears must not')
  assert.ok(got.posterior < 1 - 0.2,
    `the DISPLAYED number is bounded by the validation, got ${got.posterior}. If this ever `
    + 'approaches 1 the cap has been lost and the number is claiming certainty nobody validated')
}

// --- 3. the dials cannot be used to manufacture a call from nothing ------------------------------
//
// Loosening a threshold may admit a region the defaults refuse. It must never turn an absence of
// evidence into a named parent: with no informative markers there is nothing for any setting to
// work on, and the refusal has to survive every dial being opened at once.
{
  const none = callOneParentOrigin([], 0.20, DEFAULT_Q, GENOTYPE_ERROR, DROP_IN,
    { minMarkers: 0, maxRegionHet: 1, callPosterior: 0 })
  assert.equal(none.markers, 0)
  assert.notEqual(none.verdict, 'known-parent-lost')
  assert.notEqual(none.verdict, 'other-parent-lost')
}

console.log('om.check.ts: all assertions passed, including no-options equalling the shipped '
  + 'configuration on every branch, and every dial actually moving the answer')

// Self-check for the browser parentage layer. Run: node src/parentage.check.ts
//
// The arithmetic here is duplicated from `origin.py` deliberately, so it is checked against the
// same real measurements the Python side was validated on rather than against itself. Every
// number in section 1 came off a real file during that validation.
import assert from 'node:assert/strict'
import { readFileSync, createReadStream } from 'node:fs'
import { createGunzip } from 'node:zlib'
import { createInterface } from 'node:readline'
import { headerMap, parseRow, type ProbeRow } from './ingest.ts'
import {
  ABSENCE_MARGIN, BAF_EXTREME_FLOOR, CALL_COLLAPSE, COPY_SHIFT_FLOOR, MOSAIC_Z, absenceExplainable,
  agreement,
  classify,
  emptyTally, COPY_SHIFT_FLOOR_BULK, BULK_SPREAD_MAX,
  HET_BAND_DIPLOID,
  isAutosome, pair, pct, secondParentSignal, tallyRow, type Tally,
  BAF_CLAMPED_FLOOR, GT_FALLBACK_HAPLOID_MAX, GT_FALLBACK_DIPLOID_MIN,
} from './parentage.ts'
import { calibratedZ, intensityDetects, nullScale } from './intensityNull.ts'

import type { AB } from './informativity.ts'

const row = (
  chrom: string, pos: number, genotype: AB, baf: number | null = null,
  log2R: number | null = null,
): ProbeRow => ({
  probesetId: `m${pos}`, chrom, pos, log2R, baf, copyNumber: null,
  genotype, bestProbeset: true,
})

// --- 1. the noise bound, against pairs of known relationship ----------------------------------
//
// The bound is the product and neither term alone: a homozygous genome is immune however much
// drops out, which is why an androgenote holds at 0.16% absence at a 13.6% no-call rate while a
// diploid embryo at 46.9% reaches 9.69% against its own father.
assert.equal(absenceExplainable(0.5, 0), 0, 'homozygous: immune to dropout')
assert.equal(absenceExplainable(0, 0.5), 0, 'no dropout: nothing to inflate')
for (const [nc, het, observed] of [
  [0.004, 0.318, 0.0005], // PennCNV offspring
  [0.136, 0.013, 0.0016], // 52461 androgenote 02
  [0.137, 0.020, 0.0022], // 52461 androgenote 04
  [0.469, 0.222, 0.0969], // Zuccaro A8, a TRUE pair scoring above unrelated ones
] as const) {
  assert.ok(absenceExplainable(nc, het) >= observed,
    `bound ${absenceExplainable(nc, het)} must cover observed ${observed}`)
}
assert.ok(Number.isNaN(absenceExplainable(NaN, 0.2)))

// --- 2. the second-parent signal is derived, not fitted ---------------------------------------
assert.ok(Math.abs(secondParentSignal(0.170) - 0.085) < 1e-9)
assert.ok(Math.abs(secondParentSignal(0.322) - 0.161) < 1e-9)
assert.ok(Number.isNaN(secondParentSignal(NaN)))

// --- 3. the three classes, on synthetic genomes with the real shapes ---------------------------
function build(opts: {
  absent: number; het: number; nonParental: number; n?: number; chrom?: string
}): { t: Tally } {
  const t = emptyTally()
  const n = opts.n ?? 4000
  for (let i = 0; i < n; i += 1) {
    const s = (i * 7919) % 1000
    const parent: AB = 'AA'
    let gt: AB = 'AA'
    if (s < opts.absent * 1000) gt = 'BB'
    else if (s >= 1000 - opts.nonParental * 1000) gt = 'AB'
    tallyRow(parent, row(opts.chrom ?? '1', 1000 + i * 1000, gt,
      s < opts.het * 1000 ? 0.5 : 0.0), t)
  }
  return { t }
}

const andro = classify(build({ absent: 0.002, het: 0.02, nonParental: 0.04 }).t, 0.17)
assert.equal(andro.verdict, 'parent_genome_present')
assert.equal(andro.zygosity, 'uniparental_homozygous')
assert.equal(andro.originClass, 'androgenetic')

const gyno = classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
assert.equal(gyno.verdict, 'no_parental_contribution')
assert.equal(gyno.originClass, 'gynogenetic')

// A DIPLOID GENOME'S MID-BAND READINGS SIT AT ITS HETEROZYGOUS CALLS. That is what an AB call is,
// and it is what separates a second parental contribution from a smeared homozygous cluster: the
// band is read against the array's own homozygous clusters, not against a flat level. The knobs
// above place the band on homozygous calls, which is a smeared array rather than a diploid one, so
// this case is built where a real array puts it. See HET_BAND_EXCESS.
const bipT = ((): Tally => {
  const t = emptyTally()
  for (let i = 0; i < 4000; i += 1) {
    const s = (i * 7919) % 1000
    const gt: AB = s < 2 ? 'BB' : s < 300 ? 'AB' : 'AA'
    // 0.4% of the homozygous calls read mid-band anyway, which is the amplification smear every
    // real array carries and the baseline the boundary is measured against.
    tallyRow('AA', row('1', 1000 + i * 1000, gt, gt === 'AB' || s % 250 === 7 ? 0.5 : 0.0), t)
  }
  return t
})()
const bip = classify(bipT, 0.17)
assert.ok(bip.homBand > 0 && bip.homBand < 0.01, `a real smear, got ${bip.homBand}`)
assert.equal(bip.zygosity, 'diploid')
assert.equal(bip.originClass, 'biparental')

// THE SAME BAND ON A SMEARED ARRAY IS NOT A SECOND CONTRIBUTION. Same 30% of markers reading
// mid-band, but at HOMOZYGOUS calls, which is amplification quality rather than ploidy. A flat
// threshold on the band cannot tell these two apart and called both diploid; that read 316 of 726
// lab arrays, and 74 of 81 blastomeres, as carrying one parental contribution.
const smeared = classify(build({ absent: 0.002, het: 0.30, nonParental: 0.02 }).t, 0.17)
assert.ok(smeared.hetBand > HET_BAND_DIPLOID, 'the band alone would have said diploid')
assert.equal(smeared.zygosity, 'uniparental_homozygous', 'but it sits at homozygous calls')

// A homozygous genome never consults the narrow axis: one allele per locus, and if it is the
// parent's there is no room for a second complement.
const loud = classify(build({ absent: 0.002, het: 0.02, nonParental: 0.30 }).t, 0.17)
assert.ok(loud.nonParentalRate > loud.secondParentExpected, 'the narrow axis would say biparental')
assert.equal(loud.originClass, 'androgenetic', 'but zygosity settles it')

// --- 4. the uncalled band exists and says what would settle it --------------------------------
const mid = classify(build({ absent: 0.012, het: 0.30, nonParental: 0.30 }).t, 0.17)
assert.equal(mid.verdict, 'unclear')
assert.equal(mid.originClass, 'unclear')
assert.ok(mid.limits.some((l) => l.includes('array would measure dropout directly')))

// --- 5. a male offspring loses the paternal X by biology, not by defect ------------------------
{
  const t = build({ absent: 0.002, het: 0.02, nonParental: 0.04 }).t
  for (let i = 0; i < 4000; i += 1) {
    const s = (i * 7919) % 1000
    tallyRow('AA', row('X', 1000 + i * 1000, s < 315 ? 'BB' : 'AA', 0.0), t)
  }
  // The chrY this sample is asserted to carry. Without it the fixture only said "the X is
  // missing", and the exemption used to read that as a Y all by itself.
  for (let i = 0; i < 600; i += 1) tallyRow('AA', row('Y', 1000 + i * 1000, 'AA'), t)
  const pat = classify(t, 0.17, { role: 'paternal' })
  const x = pat.chroms.find((c) => c.chrom === 'X')
  assert.equal(x?.verdict, 'expected_absent', 'the father sent a Y instead')
  assert.equal(pat.spermType, 'Y_bearing')

  // A mother transmits an X to a child of either sex, so the same exemption must not apply.
  const mat = classify(t, 0.17, { role: 'maternal' })
  assert.equal(mat.chroms.find((c) => c.chrom === 'X')?.verdict, 'absent')
}

// --- 6. the X-bearing case a chrY test cannot reach -------------------------------------------
{
  const t = build({ absent: 0.002, het: 0.02, nonParental: 0.04 }).t
  for (let i = 0; i < 4000; i += 1) tallyRow('AA', row('X', 1000 + i * 1000, 'AA', 0.0), t)
  const r = classify(t, 0.17)
  assert.equal(r.spermType, 'X_bearing')
  assert.ok(r.notes.some((n) => n.includes('a chrY test cannot call')))
}

// --- 7. housekeeping ---------------------------------------------------------------------------
assert.ok(isAutosome('1') && isAutosome('22') && !isAutosome('X') && !isAutosome('23'))
assert.equal(pct(0.0016), '0.16%')
assert.equal(pct(NaN), 'n/a')
assert.equal(ABSENCE_MARGIN, 3)

// No-calls are excluded rather than counted as evidence in either direction.
{
  const t = emptyTally()
  for (let i = 0; i < 500; i += 1) tallyRow('AA', row('1', i * 1000, 'NC'), t)
  assert.equal(t.byChrom.size, 0, 'a no-call is not an informative marker')
  assert.equal(t.called, 0)
  assert.ok(t.markers === 500)
}


// --- 8. both parents, mirroring origin.both_parents --------------------------------------------
{
  const present = () => classify(build({ absent: 0.002, het: 0.30, nonParental: 0.30 }).t, 0.17)
  const absent = () => classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
  const middle = () => classify(build({ absent: 0.012, het: 0.30, nonParental: 0.30 }).t, 0.17)

  assert.equal(pair(present(), present(), NaN).originClass, 'biparental')
  assert.equal(pair(present(), absent(), NaN).originClass, 'androgenetic')
  assert.equal(pair(absent(), present(), NaN).originClass, 'gynogenetic')

  // The outcome a father-only run cannot reach: his absence alone reads as gynogenetic.
  const neither = pair(absent(), absent(), NaN)
  assert.equal(neither.originClass, 'neither_parent')
  assert.equal(classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
    .originClass, 'gynogenetic', 'which is exactly what one parent would have called it')
  assert.ok(neither.notes.some((n) => n.includes('check the pairing')))

  // An unresolved parent is never resolved by the other one.
  assert.equal(pair(middle(), present(), NaN).originClass, 'unclear')
  assert.equal(pair(present(), middle(), NaN).originClass, 'unclear')

  assert.ok(pair(present(), present(), 0.998).notes.some((n) => n.includes('relabelled file')))
  assert.ok(!pair(present(), present(), 0.62).notes.some((n) => n.includes('relabelled file')))
  assert.ok(!pair(present(), present(), NaN).notes.some((n) => n.includes('relabelled file')),
    'too few shared markers must not raise a duplicate warning')
}

// --- 9. parent agreement -----------------------------------------------------------------------
{
  const gt = (n: number, f: (i: number) => AB): Map<string, AB> =>
    new Map(Array.from({ length: n }, (_, i) => [`m${i}`, f(i)]))
  const a = gt(2000, (i) => (i % 3 === 0 ? 'AB' : 'AA'))
  assert.equal(agreement(a, a), 1, 'a file against itself')
  const third = gt(2000, (i) => (i % 3 === 0 ? 'AB' : 'BB'))
  assert.equal(agreement(a, third), 667 / 2000, 'agrees only where both say AB')
  assert.ok(Number.isNaN(agreement(a, gt(50, () => 'AA'))), 'too few shared markers')
  // A no-call on either side is not a disagreement, it is no observation.
  assert.equal(agreement(a, new Map([...a].map(([k, v], i) => [k, i % 4 ? v : 'NC' as AB]))), 1)
}


// --- 10. below the call-rate floor, zygosity is withheld ---------------------------------------
// The gate ingest.gates already excludes on, applied to the axis that rests on it. Three isolated
// paternal pronuclei at 53.8-59.1% call rate each showed an 18-27% heterozygous band, which a
// haploid genome cannot produce, and each was called biparental on the strength of it.
{
  const noisy = (): Tally => {
    const t = emptyTally()
    for (let i = 0; i < 4000; i += 1) {
      const s = (i * 7919) % 1000
      // Half the markers fail to call, and the band is full of spurious heterozygotes.
      tallyRow('AA', row('1', 1000 + i * 1000, s < 500 ? 'NC' : 'AA', s < 250 ? 0.5 : 0.0), t)
    }
    return t
  }
  const r = classify(noisy(), 0.17)
  assert.ok(r.hetBand > HET_BAND_DIPLOID, 'the band alone would have said diploid')
  assert.equal(r.zygosity, 'unknown', 'but the call rate says it is artefact')
  assert.equal(r.originClass, 'unclear', 'so the class cannot be asserted')
  assert.ok(r.limits.some((l) => l.includes('below 60%')), 'and the reason is stated')
  // Absence is Mendelian and survives: a genome with no paternal contribution still says so.
  const gone = emptyTally()
  for (let i = 0; i < 4000; i += 1) {
    const s = (i * 7919) % 1000
    tallyRow('AA', row('1', 1000 + i * 1000, s < 500 ? 'NC' : 'BB', s < 250 ? 0.5 : 0.0), gone)
  }
  assert.equal(classify(gone, 0.17).verdict, 'no_parental_contribution')
}


// --- 11. the chrX exemption needs a chrY MEASUREMENT, not an absent X --------------------------
// Inferring "he sent a Y" from a missing paternal X is circular. On a sample with no chrY called
// it exempted a real paternal X loss as ordinary sex determination and asserted a Y the file does
// not show, which is a silent false negative on exactly the event this tool exists to detect.
{
  const withY = (y: boolean): ReturnType<typeof classify> => {
    const t = emptyTally()
    for (let i = 0; i < 20_000; i += 1) {
      const s = (i * 7919) % 1000
      tallyRow('AA', row(String(1 + (i % 22)), i * 1000, s < 3 ? 'BB' : 'AA', s < 300 ? 0.5 : 0), t)
    }
    for (let i = 0; i < 4000; i += 1) tallyRow('AA', row('X', i * 1000, 'BB', 0.5), t)
    for (let i = 0; i < 600; i += 1) tallyRow('AA', row('Y', i * 1000, y ? 'AA' : 'NC'), t)
    return classify(t, 0.17)
  }
  const male = withY(true)
  assert.equal(male.verdict, 'parent_genome_present', 'the autosomes carry his genome')
  assert.equal(male.chroms.find((c) => c.chrom === 'X')?.verdict, 'expected_absent')
  assert.equal(male.spermType, 'Y_bearing')

  const noY = withY(false)
  const x = noY.chroms.find((c) => c.chrom === 'X')
  assert.equal(x?.verdict, 'absent', 'with no chrY this is a loss, not sex determination')
  assert.equal(noY.spermType, 'unknown', 'and no sperm type may be asserted from it')
  assert.ok(x?.note?.includes('no chrY was called'), 'and the reason is stated')
}

// --- 12. an unresolved pair still names the parent it resolved ---------------------------------
{
  const present = classify(build({ absent: 0.002, het: 0.30, nonParental: 0.30 }).t, 0.17)
  const absent = classify(build({ absent: 0.068, het: 0.03, nonParental: 0.12 }).t, 0.17)
  const middle = classify(build({ absent: 0.012, het: 0.30, nonParental: 0.30 }).t, 0.17)

  const patOpen = pair(middle, present, NaN)
  assert.equal(patOpen.originClass, 'unclear')
  assert.ok(patOpen.notes.some((n) => n.includes("oocyte donor's contribution is confirmed present")))

  const matOpen = pair(absent, middle, NaN)
  assert.ok(matOpen.notes.some((n) => n.includes("sperm donor's contribution is confirmed absent")))

  assert.ok(pair(middle, middle, NaN).notes.some((n) => n.includes('Neither parent is resolved')))
}

// --- 13. the two implementations agree, which nothing checked until now ------------------------
// tests/fixtures/parentage_cross.json is written by tools/gen_parentage_fixture.py from origin.py
// on the shipped example arrays. An independent audit found the two sides differing on every pair
// it tried, including fifteen chromosomes on one sample, because each was pinned only against
// itself. Regenerate the fixture when the CLI's numbers legitimately change; never edit it to
// match the browser.
{
  const fx = JSON.parse(
    readFileSync(new URL('../../tests/fixtures/parentage_cross.json', import.meta.url), 'utf8'),
  ) as {
    donor_heterozygosity: number
    cases: { sample: string; verdict: string; origin_class: string; zygosity: string
      sperm_type: string; genome_rate: number; explainable: number; nonpaternal_rate: number
      second_parent_expected: number; het_band: number; no_call_rate: number
      dispersion: number; min_chrom_rate: number }[]
  }
  assert.ok(fx.cases.length >= 4, 'the fixture lost cases')

  const donorGt = await readExample('GSM4472397_sperm_DNA_71.subset.csv.gz')
  for (const c of fx.cases) {
    const t = emptyTally()
    await eachExampleRow(`${c.sample}.subset.csv.gz`, (r) => {
      tallyRow(donorGt.get(r.probesetId) ?? 'NC', r, t)
    })
    const got = classify(t, fx.donor_heterozygosity)
    const near = (a: number, b: number, what: string) => assert.ok(
      Math.abs(a - b) < 1e-6, `${c.sample} ${what}: browser ${a}, CLI ${b}`,
    )
    assert.equal(got.originClass, c.origin_class, `${c.sample} class`)
    assert.equal(got.zygosity, c.zygosity, `${c.sample} zygosity`)
    assert.equal(got.spermType, c.sperm_type, `${c.sample} sperm type`)
    assert.equal(got.verdict, c.verdict, `${c.sample} verdict`)
    near(got.genomeRate, c.genome_rate, 'absence')
    near(got.explainable, c.explainable, 'ceiling')
    near(got.nonParentalRate, c.nonpaternal_rate, 'alleles the donor lacks')
    near(got.secondParentExpected, c.second_parent_expected, 'second-parent expectation')
    near(got.hetBand, c.het_band, 'BAF band')
    near(got.dispersion, c.dispersion, 'dispersion')
    near(got.minChromRate, c.min_chrom_rate, 'cleanest chromosome')
    near(got.noCallRate, c.no_call_rate, 'no-call rate')
  }
}

// --- 14. the pseudoautosomal region is scored apart from the rest of chrX ----------------------
//
// The PAR is on both the X and the Y, so a Y-bearing sperm delivers it while the rest of the X
// is legitimately gone. Pooled into one bucket that positive control was invisible.
{
  const t = emptyTally()
  t.build = 'GRCh37'
  // PAR1 is 60001-2699520 on GRCh37: present, as it must be from either sperm type.
  for (let i = 0; i < 600; i += 1) tallyRow('AA', row('X', 70_000 + i * 4000, 'AA', 0.0), t)
  // The rest of chrX, absent, as it is from a Y-bearing sperm.
  for (let i = 0; i < 600; i += 1) tallyRow('AA', row('X', 20_000_000 + i * 1000, 'BB', 0.0), t)
  for (let i = 0; i < 600; i += 1) tallyRow('AA', row('Y', 1000 + i * 1000, 'AA'), t)
  for (let i = 0; i < 4000; i += 1) {
    tallyRow('AA', row('1', i * 1000, 'AA', i % 100 < 16 ? 0.5 : 0.0), t)
  }
  const r = classify(t, 0.17, { role: 'paternal' })
  const par = r.chroms.find((c) => c.chrom === 'X:PAR')
  const rest = r.chroms.find((c) => c.chrom === 'X')
  assert.equal(par?.verdict, 'present')
  assert.equal(rest?.verdict, 'expected_absent')
  assert.equal(r.spermType, 'Y_bearing')
  // The PAR is not an autosome and must not move the genome-wide figure either way.
  assert.ok(r.genomeRate < 0.001)
}

// --- 15. uniform absence everywhere reads as a mixture, patchy absence does not ----------------
//
// Measured: two arrays of one man 1.11, a degraded true offspring 0.76, an unrelated adult 0.11,
// two genomes blended 0.10. The rule fires only inside the uncalled band, where it decides
// nothing on its own and only says which of two shapes the sample has.
{
  const uniform = emptyTally()
  const patchy = emptyTally()
  for (let c = 1; c <= 22; c += 1) {
    for (let i = 0; i < 400; i += 1) {
      const s = (i * 7919) % 1000
      // Same rate on every chromosome, against the same total confined to four of them. The
      // band is held at the 16% a real diploid reads: 0.5 everywhere would be 100%, which the
      // plausibility gate refuses and rightly so.
      const baf = i % 100 < 16 ? 0.5 : 0.0
      tallyRow('AA', row(String(c), 1000 + i * 1000, s < 25 ? 'BB' : 'AA', baf), uniform)
      tallyRow('AA', row(String(c), 1000 + i * 1000,
        c <= 4 && s < 138 ? 'BB' : 'AA', baf), patchy)
    }
  }
  // No-calls the parent does not carry, so they raise the dropout ceiling without entering the
  // absence count. 500 of them against 8,800 informative markers gives a 5.4% no-call rate,
  // which at a 16% band puts the ceiling at 1.36% and leaves 2.5% absence inside the uncalled
  // band at 1.8x. An earlier version reached the same place by setting every BAF to 0.5, which
  // is a 100% band and not a thing a genome does.
  for (let i = 0; i < 500; i += 1) tallyRow('NC', row('1', 9e8 + i, 'NC'), uniform)
  for (let i = 0; i < 500; i += 1) tallyRow('NC', row('1', 9e8 + i, 'NC'), patchy)
  const u = classify(uniform, 0.17)
  const q = classify(patchy, 0.17)
  assert.equal(u.verdict, 'unclear')
  assert.equal(q.verdict, 'unclear')
  assert.ok(u.dispersion < q.dispersion)
  const mixed = (r: typeof u) => r.limits.some((l) => l.includes('mixed sample'))
  assert.ok(mixed(u), 'uniform absence should read as a mixture')
  assert.ok(!mixed(q), 'absence confined to four chromosomes is a loss, not a mixture')
}

// --- 16. a reconstructed reference raises the ceiling by its own known false absence -----------
//
// Measured on the five public pronuclei: the leave-one-out references carry 2.3 to 3.6%
// contamination, so a true product picks up over 1% absence through no fault of its own, against
// a ceiling of about 1%. Three of the five true products read unclear until the reference's own
// floor was admitted into the ceiling. The floor is `Reference.spuriousAbsence`, which is
// contamination/2 because a haploid product carries the other allele half the time.
{
  // 1.1% absence on a clean sample: the shape a true product has against a dirty reference.
  const t = build({ absent: 0.011, het: 0.03, nonParental: 0.05 }).t
  assert.equal(classify(t, NaN).verdict, 'unclear',
    'a measured reference contributes nothing, so this sits above the sample\'s own noise')
  assert.equal(classify(t, NaN, { spuriousAbsence: 0.0165 }).verdict, 'parent_genome_present',
    'against a reference whose own contamination explains it, the same sample is present')
  assert.ok(classify(t, NaN, { spuriousAbsence: 0.0165 }).explainable
    > classify(t, NaN).explainable + 0.016, 'the floor is added, not merged into the noise term')

  // The term must not rescue a genuinely unrelated sample. Unrelated haploids read 9.9 to 16.1%
  // opposite homozygotes, far above 3x even the widened ceiling.
  assert.equal(classify(build({ absent: 0.099, het: 0.03, nonParental: 0.12 }).t, NaN,
    { spuriousAbsence: 0.018 }).verdict, 'no_parental_contribution')

  // Omitted or unmeasurable, it changes nothing: the diploid path passes no term and is
  // calibrated without one.
  assert.equal(classify(t, NaN, { spuriousAbsence: NaN }).explainable, classify(t, NaN).explainable)
  assert.equal(classify(t, NaN, {}).explainable, classify(t, NaN).explainable)
}

/** Stream one shipped example, gunzipping, exactly as the browser parses a dropped file. */
async function eachExampleRow(name: string, fn: (r: ProbeRow) => void): Promise<void> {
  const path = new URL(`../public/examples/${name}`, import.meta.url)
  const gz = createReadStream(path).pipe(createGunzip())
  let map: ReturnType<typeof headerMap> = null
  for await (const line of createInterface({ input: gz, crlfDelay: Infinity })) {
    if (!map) { map = headerMap(line); continue }
    const r = parseRow(line, map)
    if (r) fn(r)
  }
}

async function readExample(name: string): Promise<Map<string, AB>> {
  const gt = new Map<string, AB>()
  await eachExampleRow(name, (r) => { gt.set(r.probesetId, r.genotype) })
  return gt
}

console.log('parentage.check.ts OK')

// --- 17. a chromosome whose calls are not measuring it gets no verdict -------------------------
//
// Array mis-clustering makes a chromosome's genotypes systematically wrong. Absence then rises to
// something indistinguishable from a real chromosome-scale loss, and every genotype-derived
// measure agrees with it, because they are all reading the same broken calls. Only the allelic
// ratio shows the cause.
//
// Measured on the five paternal pronuclei of GSE148488, each a meiotic product of the sperm donor
// and therefore present on every autosome: real chromosomes run 0.752 to 0.976 B-allele
// frequencies at an extreme, while GSM4774681 chr1 runs 0.130 and was reported ABSENT at 18.36%
// on a sample whose genome-wide verdict is "parent genome present".
{
  const build = (extreme: number, absent: number) => {
    const t = emptyTally()
    for (let i = 0; i < 4_000; i += 1) {
      const s = (i * 7919) % 1000
      // The allelic ratio, which is what mis-clustering corrupts, set independently of the
      // genotype so the two signals cannot be confused for one another.
      const baf = s < extreme * 1000 ? (s % 2 ? 0.02 : 0.98) : 0.5
      tallyRow('AA', row('7', 1000 + i * 1000, s < absent * 1000 ? 'BB' : 'AA', baf), t)
    }
    return classify(t, 0.17).chroms.find((c) => c.chrom === '7')!
  }

  // Clean clustering, heavy absence: a real loss, and it is still called.
  const real = build(0.90, 0.18)
  assert.ok(real.bafExtreme > BAF_EXTREME_FLOOR)
  assert.equal(real.verdict, 'absent', 'a correctly clustered chromosome is still called')

  // The same absence, on a chromosome whose calls are not measuring it.
  const broken = build(0.13, 0.18)
  assert.ok(broken.bafExtreme < BAF_EXTREME_FLOOR)
  assert.equal(broken.verdict, 'not_measured',
    'mis-clustering must withhold the verdict, not produce a loss')
  assert.ok(broken.note!.includes('not measuring this chromosome'),
    'and the reason has to travel with it')
  // The rate is still reported: it is an observation. What is withheld is the verdict.
  assert.ok(Math.abs(broken.rate - real.rate) < 1e-9, 'the measured rate is unchanged')

  // The floor sits in an empty gap, so neither edge of the real range can cross it.
  assert.equal(BAF_EXTREME_FLOOR, 0.40)
  assert.ok(BAF_EXTREME_FLOOR > 0.130 * 2, 'clear of the mis-clustered chromosome')
  assert.ok(BAF_EXTREME_FLOOR < 0.752 / 1.8, 'and clear of the worst correctly clustered one')
}

// --- 18. a mosaic is seen in the allelic ratio, where the genotype structurally cannot ----------
//
// A genotype call is a threshold on that ratio. The AB-against-BB boundary sits near 0.917 and a
// chromosome at mosaic fraction f has a true ratio of 1/(2-f), which does not cross it until
// f = 0.909, so below that the CALL never changes and more markers buy nothing. The ratio itself
// is not thresholded. Measured on four bulk diploid arrays, 88 chromosome observations with no
// mosaic anywhere, z runs -1.65 to 5.18; a mosaic titrated onto one of their chromosomes reads
// 9.7/4.6/22.7/12.5 at f=0.30 and 21.3/11.7/48.5/28.6 at f=0.50.
{
  // A diploid genome: heterozygous everywhere, ratios at 0.5 with real spread.
  const build = (mosaicChrom: string | null, f: number) => {
    const t = emptyTally()
    for (let c = 1; c <= 12; c += 1) {
      for (let i = 0; i < 1_000; i += 1) {
        const s = (i * 7919) % 1000
        // Spread about 0.5, and a small per-chromosome offset. Without the offset every
        // chromosome has an identical mean, the standard deviation of the others is exactly zero
        // and the contrast is undefined: the statistic correctly returns NaN, and the fixture
        // rather than the code was wrong. Real chromosomes differ.
        const noise = ((s % 21) - 10) / 200 + ((c * 37) % 11) / 4_000
        const shift = String(c) === mosaicChrom ? ((i % 2 ? 1 : -1) * (1 / (2 - f) - 0.5)) : 0
        const baf = Math.min(1, Math.max(0, 0.5 + noise + shift))
        tallyRow('AA', row(String(c), 1000 + i * 1000, 'AB', baf), t)
      }
    }
    return classify(t, 0.17)
  }

  const clean = build(null, 0)
  assert.equal(clean.zygosity, 'diploid', 'the fixture has to be diploid for this to run at all')
  const worst = Math.max(...clean.chroms.map((c) => c.mosaicZ).filter(Number.isFinite))
  assert.ok(worst < MOSAIC_Z, `a genome with no mosaic must not report one: worst z ${worst}`)
  assert.ok(!clean.limits.some((l) => l.includes('MIXTURE OF CELL LINEAGES')))

  const half = build('7', 0.5)
  const hit = half.chroms.find((c) => c.chrom === '7')!
  assert.ok(hit.mosaicZ >= MOSAIC_Z, `a 50% mosaic must be seen: z ${hit.mosaicZ}`)
  assert.ok(half.limits.some((l) => l.includes('MIXTURE OF CELL LINEAGES') && l.includes('chr7')),
    'and it must be reported, naming the chromosome')
  // No fraction is claimed anywhere, because inverting this statistic is biased low by half.
  assert.ok(half.limits.some((l) => l.includes('No fraction is reported')))

  // Stronger mosaic, stronger signal. A statistic that saturated would hide the severe case.
  assert.ok(build('7', 0.9).chroms.find((c) => c.chrom === '7')!.mosaicZ > hit.mosaicZ)

  // A UNIPARENTAL genome is homozygous by construction, so it has no heterozygous sites for a
  // mixture to shift and any deviation there is artefact. The axis is withheld, not guessed.
  const hap = emptyTally()
  for (let c = 1; c <= 12; c += 1) {
    for (let i = 0; i < 1_000; i += 1) {
      tallyRow('AA', row(String(c), 1000 + i * 1000, 'AA', i % 2 ? 0.02 : 0.98), hap)
    }
  }
  const uni = classify(hap, 0.17)
  assert.notEqual(uni.zygosity, 'diploid')
  assert.ok(uni.chroms.every((c) => Number.isNaN(c.mosaicZ)),
    'a uniparental genome has no heterozygous sites, so no mosaic contrast exists')
}

// --- 19. a chromosome that is GONE is a loss, not a chromosome that failed to measure ----------
//
// This pins a correction. Three chromosomes were previously diagnosed as array mis-clustering and
// withheld, on the reasoning that a paternal pronucleus is its father's on every autosome by
// construction. That was wrong: the series exists because chromosomes ARE lost in these embryos.
// Both a loss and a mis-clustering depress the allelic ratio; the CALL RATE separates them,
// because a mis-clustered chromosome still calls and an absent one cannot.
//
// Measured over 1,012 chromosome observations from 46 arrays: eleven sit at 0.20x to 0.41x the
// genome's median call rate and the other 1,001 at 0.78x to 1.16x. Every one of the eleven also
// carries |log2R| of 1.59 to 2.04 where the rest never leave -0.79 to +0.42.
{
  // `gone` collapses the call rate and the intensity, as an absent chromosome does. `broken`
  // keeps both and only corrupts the allelic ratio, as mis-clustering does.
  const build = (mode: 'clean' | 'gone' | 'broken', shift = -1.9) => {
    const t = emptyTally()
    for (let c = 1; c <= 12; c += 1) {
      const target = c === 7
      for (let i = 0; i < 3_000; i += 1) {
        const s = (i * 7919) % 1000
        const bad = target && mode !== 'clean'
        // A lost chromosome yields no DNA: four fifths of its probes say nothing at all.
        const silent = target && mode === 'gone' && s < 800
        const baf = bad ? 0.5 : (s % 2 ? 0.02 : 0.98)
        tallyRow('AA', row(String(c), 1000 + i * 1000, silent ? 'NC' : 'AA', baf, target && mode === 'gone' ? shift : 0), t)
      }
    }
    return classify(t, 0.17).chroms.find((x) => x.chrom === '7')!
  }

  const gone = build('gone')
  assert.equal(gone.aneuploidy, 'loss', `a collapsed chromosome is a LOSS: ${JSON.stringify(gone)}`)
  assert.equal(gone.verdict, 'absent', 'and it is reported, not withheld')
  assert.ok(gone.callFraction < CALL_COLLAPSE)
  assert.ok(gone.note!.includes('LOST'))

  // The same collapse with the intensity going the other way is a gain, not a loss.
  assert.equal(build('gone', 1.9).aneuploidy, 'gain')

  // Mis-clustering still withholds: the ratio is corrupt but the chromosome is being genotyped.
  const broken = build('broken')
  assert.equal(broken.aneuploidy, undefined, 'a chromosome that still calls has not been lost')
  assert.equal(broken.verdict, 'not_measured', 'so the ratio may still diagnose mis-clustering')

  const clean = build('clean')
  assert.equal(clean.aneuploidy, undefined)
  assert.ok(clean.callFraction > 0.9, 'an intact chromosome calls with the rest of the genome')

  // The threshold sits in the empty gap between the two measured populations.
  assert.equal(CALL_COLLAPSE, 0.60)
  assert.ok(CALL_COLLAPSE > 0.41 && CALL_COLLAPSE < 0.78)
}

// --- 20. a chromosome that still GENOTYPES, and which way it went -----------------------------
//
// The call rate is one entry channel and it sees one thing: a chromosome absent from the array.
// A chromosome at one copy still genotypes, and on bulk DNA it genotypes perfectly; three copies
// genotype better than two. With the call rate as the only gate, 0 of 271 constructed
// whole-chromosome events entered at all and 0 of 235 constructed gains were ever called a gain,
// because the gain test was a magnitude gate at 1.0 and a trisomy is log2(3/2) = 0.585.
//
// So intensity is the second channel, read against the array's own null (intensityNull.ts). Its
// false-positive rate is measured on chromosomes that are event-free by construction:
// audit/intensity-null-truenegatives.ts, 4 of 264 = 0.0152.
{
  // Twenty-two autosomes, every call intact, with per-chromosome intensity noise so the array has
  // a null of its own. The offsets are a fixed pattern rather than a generator: a check whose
  // fixture is random is a check that fails on somebody else's afternoon.
  const buildAll = (shift: number) => {
    const t = emptyTally()
    for (let c = 1; c <= 22; c += 1) {
      const noise = ((c * 37) % 11 - 5) / 100  // -0.05 to +0.05 between chromosomes
      for (let i = 0; i < 2_000; i += 1) {
        const s = (i * 7919) % 1000
        tallyRow('AA', row(String(c), 1000 + i * 1000, 'AA', s % 2 ? 0.02 : 0.98,
          noise + (c === 7 ? shift : 0) + (s % 7 - 3) / 200), t)
      }
    }
    return classify(t, 0.17).chroms
  }
  const build = (shift: number) => buildAll(shift).find((x) => x.chrom === '7')!

  // What the CALIBRATED z ALONE would say, built from the same units the detector builds. Every
  // reported lrrShift is its chromosome's median minus the genome's, so the set differs from the
  // set of chromosome medians by one constant; a median-centred scale is unmoved by that, and the
  // z is identical. This exists so the floor's assertion below demonstrates the gap rather than
  // asserting it.
  const zClears = (shift: number) => {
    const cs = buildAll(shift)
    const n = nullScale(cs.map((x) => x.lrrShift))
    return intensityDetects(calibratedZ(cs.find((x) => x.chrom === '7')!.lrrShift, n), true)
  }

  // Three copies: +log2(3/2). The call rate never moves.
  const tri = build(Math.log2(3 / 2))
  assert.equal(tri.aneuploidy, 'gain',
    `a trisomy is a GAIN, not a loss and not silence: ${JSON.stringify(tri)}`)
  assert.ok(tri.callFraction > 0.9, 'and it genotypes with the rest of the genome')
  assert.ok(tri.note!.includes('GAINED'))

  // One copy: -1.0, the call rate intact. This is the case the old gate could not see at all.
  const mono = build(-1.0)
  assert.equal(mono.aneuploidy, 'loss', `a genotyping monosomy is a LOSS: ${JSON.stringify(mono)}`)
  assert.ok(mono.callFraction > 0.9)
  // And it must NOT short-circuit the parental verdict: the allelic ratio is still measuring this
  // chromosome, so whose copy survived is read from it rather than asserted from the copy number.
  assert.notEqual(mono.verdict, 'absent',
    'a chromosome that still calls has a measurable allelic ratio; do not assert absence from '
    + 'copy number alone')

  // A GAIN NAMES NO PARENT. The absence statistic that names one on a loss reads nothing on a
  // gain: a third copy makes the sample miss FEWER alleles whichever parent supplied it, so both
  // answers land in the same low-absence branch. This fixture is homozygous throughout, which is
  // that branch, and the old code returned 'other' from it.
  assert.equal(tri.aneuploidyParent, undefined,
    'whose extra copy it is comes from the allelic ratio, which is not measured here; naming a '
    + 'parent from absence on a gain names one at random')
  assert.ok(/NOT RESOLVED/.test(tri.note!), 'and the note has to say so, not imply a parent')
  assert.ok(!/the copy that went/.test(tri.note!), 'loss language must not appear on a gain')

  // An intact chromosome, on the same fixture, stays silent.
  assert.equal(build(0).aneuploidy, undefined)
  // And a shift under the calibrated threshold stays silent too, which is what keeps the
  // false-positive rate where it was measured.
  assert.equal(build(0.05).aneuploidy, undefined,
    'a shift inside the array\'s own spread is not an event')

  // THE FLOOR, and this is the one the corpus paid for. A relative threshold has no physical
  // scale: 0.30 log2 on this fixture clears the calibrated z several times over, and no
  // whole-chromosome copy state lives there. Over all 137 arrays of GSE148488 the z alone added
  // four events at |log2R| 0.20 to 0.33, and four arrays of ONE trophectoderm biopsy read +0.174,
  // +0.201, +0.163 and +0.168 on chr19 with exactly one called a gain. See COPY_SHIFT_FLOOR.
  const belowFloor = build(0.30)
  assert.ok(zClears(0.30), 'the calibrated z alone would admit 0.30, which is the whole point')
  assert.equal(belowFloor.aneuploidy, undefined,
    'a shift no copy state can produce is not an event, however far it sits from this array\'s '
    + 'own noise')
  assert.ok(0.33 < COPY_SHIFT_FLOOR && COPY_SHIFT_FLOOR < 0.559,
    'the floor must sit above the largest false shift measured on real arrays and below the '
    + 'smallest true one, which is a trisomy read on trophectoderm')
  // And the smallest event that physically exists still gets through, on all three materials.
  for (const observed of [0.585, 0.559, 0.615]) {
    assert.equal(build(observed).aneuploidy, 'gain',
      `a trisomy reads ${observed} log2 on real material and must still be detected`)
  }
}

// --- 21. THE GENOME-LEVEL CLASS MUST FOLLOW WHICH PARENT WAS LOADED ----------------------------
//
// EVERY CHECK IN THIS FILE LOADED A PATERNAL PARENT, AND THAT IS HOW THIS SURVIVED. `classify`
// hardcoded `androgenetic` for a genome carrying the loaded parent's alleles and `gynogenetic` for
// one that did not, which is right only when the loaded parent is the father. The whole tool
// supports loading the mother alone and says so in its own comments.
//
// Measured on the seven maternal pronuclei of GSE148488, each resolved to egg donor B by linkage
// and loaded under its true maternal role: all seven came back ANDROGENETIC. For a zygote-stage
// genome that word IS the answer, so the foundational call was inverted on every maternal run.
{
  // A genome carrying ONE parental complement, and it is the loaded parent's: every call is the
  // parent's own homozygous genotype, so nothing the other parent could have supplied is present.
  const uniparental = () => {
    const t = emptyTally()
    for (let c = 1; c <= 22; c += 1) {
      for (let i = 0; i < 2_000; i += 1) {
        const g = (c + i) % 2 ? 'AA' : 'BB'
        tallyRow(g, row(String(c), 1000 + i * 1000, g, g === 'AA' ? 0.02 : 0.98,
          ((c * 37) % 11 - 5) / 100), t)
      }
    }
    return t
  }
  const t = uniparental()
  const asPaternal = classify(t, 0.17, { role: 'paternal' })
  const asMaternal = classify(t, 0.17, { role: 'maternal' })

  assert.equal(asPaternal.zygosity, 'uniparental_homozygous',
    `the fixture must reach the uniparental branch, got ${asPaternal.zygosity}`)
  assert.equal(asMaternal.zygosity, asPaternal.zygosity,
    'and zygosity does not depend on which parent was loaded')

  // THE POINT. Same genome, same data, one loaded array. Which class it is depends entirely on
  // whose array that was.
  assert.equal(asPaternal.originClass, 'androgenetic',
    'the father\'s genome alone is androgenetic')
  assert.equal(asMaternal.originClass, 'gynogenetic',
    `the MOTHER'S genome alone is GYNOGENETIC, not androgenetic. Got `
    + `${asMaternal.originClass}. This is the inversion that reported all seven maternal `
    + 'pronuclei of GSE148488 as androgenetic.')
  assert.notEqual(asPaternal.originClass, asMaternal.originClass,
    'a class that does not move when the loaded parent changes is not reading parentage')

  // AND THE MIRROR: a genome carrying NONE of the loaded parent belongs to the other one.
  const absent = () => {
    const t2 = emptyTally()
    for (let c = 1; c <= 22; c += 1) {
      for (let i = 0; i < 2_000; i += 1) {
        // The parent is homozygous for one allele and the sample is homozygous for the other, so
        // nothing the parent could have transmitted is present anywhere.
        const pg = (c + i) % 2 ? 'AA' : 'BB'
        const sg = pg === 'AA' ? 'BB' : 'AA'
        tallyRow(pg, row(String(c), 1000 + i * 1000, sg, sg === 'AA' ? 0.02 : 0.98,
          ((c * 37) % 11 - 5) / 100), t2)
      }
    }
    return t2
  }
  const a = absent()
  const noPat = classify(a, 0.17, { role: 'paternal' })
  const noMat = classify(a, 0.17, { role: 'maternal' })
  assert.equal(noPat.originClass, 'gynogenetic',
    `no paternal contribution means gynogenetic, got ${noPat.originClass}`)
  assert.equal(noMat.originClass, 'androgenetic',
    `no MATERNAL contribution means ANDROGENETIC, got ${noMat.originClass}`)

  // AND NO NOTE MAY NAME THE WRONG PARENT EITHER. The notes were keyed on the class name, so they
  // inverted with it.
  assert.ok(!asMaternal.notes.some((n) => /two sperm/.test(n)),
    `a maternal run must not explain a duplicated maternal genome as "two sperm": `
    + JSON.stringify(asMaternal.notes))
}


// --- 22. A ZERO NOISE FLOOR IS NOT THE SAME FACT AS A CLEAN ARRAY ----------------------------
//
// THE RUN THAT FOUND THIS. GSE19247 converted from its raw intensities reads 23 of 23 known-haploid
// sperm cells as diploid, while the same call gets every haploid pronucleus in GSE148488 right.
// The cause is not the boundary.
//
// WHY THE FLOOR CAN BE UNMEASURABLE. `HET_BAND_EXCESS` subtracts the mid-band mass at HOMOZYGOUS
// calls, and that quantity only exists if the genotype and the B-allele frequency were derived
// INDEPENDENTLY. Affymetrix calls genotypes with one algorithm and computes the allelic ratio with
// another, so the two disagree at a small rate and that disagreement is the floor: 0.0020 to 0.0034
// on the arrays every threshold here was measured on. A cluster-file BAF is a function of the same
// theta the genotype came from, so a marker called homozygous CANNOT read mid-band and the floor is
// identically 0.0000. On the converted series it is exactly that on all three materials, alongside
// 70 to 85 percent of homozygous calls sitting on an exact 0 or 1.
//
// The branch read that zero as "this array has no amplification smear to correct for" and applied
// the UNAMPLIFIED boundary to single sperm cells, the most heavily amplified material in the series.
// Absence of a floor and inability to measure one are different facts. This pins the difference.
{
  /**
   * A HAPLOID genome, every marker homozygous, with an amplification smear.
   *
   * The smear is carried by markers the caller FAILED, which is where it really sits: a mid-scale
   * reading on a one-parent genome is an artefact, and a caller that derives both quantities from
   * it produces a no-call rather than a homozygous call at 0.5. Those markers still carry a BAF,
   * which is why the band sees them and `HET_BAND_EXCESS` was written to correct for them.
   */
  const haploid = (clamped: boolean) => {
    const t = emptyTally()
    for (let c = 1; c <= 22; c += 1) {
      for (let i = 0; i < 2_000; i += 1) {
        const homAllele: AB = (c + i) % 2 ? 'AA' : 'BB'
        // One marker in eight reads mid-band, in both arms, so the BAND IS IDENTICAL and nothing
        // about its level distinguishes them. 0.125 is over HET_BAND_DIPLOID, so a flat threshold
        // calls this haploid genome diploid outright.
        const smear = i % 8 === 0
        if (clamped) {
          // A CALLER THAT DERIVES BOTH QUANTITIES FROM ONE THETA. A mid-scale reading cannot come
          // back as a homozygous call, so the smear lands on no-calls and the homozygous calls sit
          // on exact endpoints. The floor has nowhere to be measured.
          if (smear) tallyRow(homAllele, row(String(c), 1000 + i * 1000, 'NC', 0.5, 0), t)
          else {
            tallyRow(homAllele,
              row(String(c), 1000 + i * 1000, homAllele, homAllele === 'AA' ? 0 : 1, 0), t)
          }
        } else if (smear) {
          // INDEPENDENTLY DERIVED. The genotype algorithm calls this marker homozygous and the
          // allelic ratio disagrees with it, which is exactly the disagreement the floor measures.
          tallyRow(homAllele, row(String(c), 1000 + i * 1000, homAllele, 0.5, 0), t)
        } else {
          tallyRow(homAllele,
            row(String(c), 1000 + i * 1000, homAllele, homAllele === 'AA' ? 0.01 : 0.99, 0), t)
        }
      }
    }
    return t
  }

  const clamped = classify(haploid(true), 0.17, { role: 'paternal' })
  const continuous = classify(haploid(false), 0.17, { role: 'paternal' })

  // The fixture must actually reproduce the situation, or the rest proves nothing.
  assert.ok(clamped.hetBand > HET_BAND_DIPLOID,
    `the band must be over the flat boundary, got ${clamped.hetBand}`)
  assert.equal(clamped.homBand, 0, `the clamped arm must have no measurable floor, `
    + `got ${clamped.homBand}`)
  assert.ok(continuous.homBand > 0, `the continuous arm must have one, got ${continuous.homBand}`)
  assert.ok(clamped.homPinnedRate >= BAF_CLAMPED_FLOOR,
    `a clamped file must be recognised as clamped, got ${clamped.homPinnedRate}`)
  assert.ok(continuous.homPinnedRate < BAF_CLAMPED_FLOOR,
    `a continuous file must not be, got ${continuous.homPinnedRate}`)

  // THE REGRESSION. This genome is haploid by construction. Before the fix the clamped arm took
  // the flat unamplified boundary, on a band of 0.125 against a 0.08 threshold, and called it
  // diploid.
  assert.equal(clamped.zygosity, 'uniparental_homozygous',
    'a clamped BAF on a HAPLOID genome must not be called diploid by the unamplified boundary, '
    + `got ${clamped.zygosity} at hetBand ${clamped.hetBand} homBand ${clamped.homBand}`)
  // And the operator is told why the weaker measure was used, rather than the file silently
  // getting a boundary that does not apply to it.
  assert.ok(clamped.limits.some((l) => /exactly 0 or 1/.test(l)),
    `the clamped run must say why: ${JSON.stringify(clamped.limits)}`)

  // THE CONTROL. Where the floor is measurable the correction still runs, and it still reads this
  // genome as one parent's: the smear appears at homozygous calls too, so the EXCESS is small.
  assert.equal(continuous.zygosity, 'uniparental_homozygous',
    `the correction must still work on a continuous file, got ${continuous.zygosity} `
    + `at excess ${continuous.hetBand - continuous.homBand}`)
  assert.ok(!continuous.limits.some((l) => /exactly 0 or 1/.test(l)),
    'a continuous file must not be told its BAF is clamped')
}

// --- 23. THE GENOTYPE FALLBACK MUST REFUSE THE MIDDLE, NOT GUESS IT ---------------------------
//
// THE RUN THAT FOUND THIS, and it was a regression introduced by section 22's own fix. Sending a
// clamped-BAF file to the genotype fallback took the converted second platform from 23 wrong calls
// to 38: it read 15 of 23 known-haploid sperm right, and 29 of 30 known-DIPLOID single lymphoblasts
// wrong. Both classes sit at the same place on this measure, 0.32 and 0.41 of the loaded parent's
// heterozygosity, because what put them there is dropout rather than ploidy. A single split at the
// midpoint cannot do anything but trade one error for the other.
//
// Bulk material with the same clamped BAF sits at 0.971 and is answered correctly. So the fallback
// keeps the case it can measure and refuses the case it cannot.
{
  const parentHet = 0.17
  /** `abEvery` controls the sample's heterozygosity, which is the only thing under test here. */
  const mid = (abEvery: number) => {
    const t = emptyTally()
    for (let c = 1; c <= 22; c += 1) {
      for (let i = 0; i < 2_000; i += 1) {
        const homAllele: AB = (c + i) % 2 ? 'AA' : 'BB'
        // No B-allele frequency at all, so the fallback is the only branch available.
        const gt: AB = i % abEvery === 0 ? 'AB' : homAllele
        tallyRow(homAllele, row(String(c), 1000 + i * 1000, gt, null, 0), t)
      }
    }
    return t
  }

  // One AB in fifteen is 6.7% heterozygous, which is 0.39 of the parent's: squarely in the gap
  // where the two amplified classes overlapped.
  const ambiguous = classify(mid(15), parentHet, { role: 'paternal' })
  // With no BAF at all, `hetFraction` IS the genotype heterozygosity the fallback reads.
  const share = ambiguous.hetFraction / parentHet
  assert.ok(share > GT_FALLBACK_HAPLOID_MAX && share < GT_FALLBACK_DIPLOID_MIN,
    `the fixture must land in the refused gap, got ${share}`)
  assert.equal(ambiguous.zygosity, 'unknown',
    `heterozygosity at ${share.toFixed(3)} of the parent's separates nothing, so it must be `
    + `refused rather than called, got ${ambiguous.zygosity}`)
  assert.ok(ambiguous.limits.some((l) => /about amplification rather than ploidy/.test(l)),
    `the refusal must say why: ${JSON.stringify(ambiguous.limits)}`)

  // BOTH ENDS STILL ANSWER. One AB in 200 is 0.5%, well under a quarter of the parent's.
  const clearlyHaploid = classify(mid(200), parentHet, { role: 'paternal' })
  assert.equal(clearlyHaploid.zygosity, 'uniparental_homozygous',
    `almost no heterozygosity is one parent's genome, got ${clearlyHaploid.zygosity}`)
  // One AB in five is 20%, over three quarters of the parent's 17%.
  const clearlyDiploid = classify(mid(5), parentHet, { role: 'paternal' })
  assert.equal(clearlyDiploid.zygosity, 'diploid',
    `heterozygosity above the parent's own is two contributions, got ${clearlyDiploid.zygosity}`)
}

// --- THE MAGNITUDE FLOOR FOLLOWS THE MATERIAL -----------------------------------------------------
//
// COPY_SHIFT_FLOOR is 0.40 because four arrays of ONE biopsy disagree with each other by 0.33, and
// that disagreement is an amplification artefact. Bulk DNA does not carry it: over 60 arrays of a
// series stating each sample free of clinical variation, the largest shift on any autosome is
// 0.070, while a stated trisomy on the same chemistry shifts by 0.299 to 0.357. Held at 0.40, the
// gate refuses every one of those positives.
{
  const tally = (shifted: string, shift: number, wobbleStep = 0.002): Tally => {
    const t = emptyTally()
    for (let c = 1; c <= 22; c += 1) {
      const chrom = String(c)
      for (let i = 0; i < 300; i += 1) {
        // A quiet array: chromosomes sit a few thousandths apart, which is what the measured
        // spread between one bulk array's chromosomes looks like, and the null needs a spread.
        const wobble = ((c % 7) - 3) * wobbleStep
        const l2 = (chrom === shifted ? shift : 0) + wobble + ((i % 5) - 2) * 0.0005
        tallyRow('AA', row(chrom, 1000 + i * 1000, i % 4 === 0 ? 'AB' : 'AA', 0.5, l2), t)
      }
    }
    return t
  }
  const gainOf = (r: ReturnType<typeof classify>, chrom: string) =>
    r.chroms.find((x) => x.chrom === chrom)?.aneuploidy

  // BULK AND QUIET: the stated trisomies of the bulk series, called.
  assert.equal(gainOf(classify(tally('21', 0.30), 0.17, { bulk: true }), '21'), 'gain',
    'a stated trisomy shifts by 0.30 on bulk, and a quiet bulk array must call it')
  assert.equal(gainOf(classify(tally('21', 0.070), 0.17, { bulk: true }), '21'), undefined,
    'the worst drift measured on a normal bulk array must not be called')

  // NOT READ AS BULK: the drift floor stands however quiet the array is, because an amplified
  // blastomere can be as quiet as a bulk array.
  assert.equal(gainOf(classify(tally('21', 0.30), 0.17, { bulk: false }), '21'), undefined)
  assert.equal(gainOf(classify(tally('21', 0.30), 0.17), '21'), undefined,
    'an unstated material is treated as amplified, which is the conservative side')

  // READ AS BULK BUT NOISY: single cells amplified by MDA land on the bulk rung through drop-in,
  // and their own spread is what gives them away.
  assert.equal(gainOf(classify(tally('21', 0.30, 0.013), 0.17, { bulk: true }), '21'), undefined,
    'an array with a single cell\'s spread keeps the drift floor even when read as bulk')
  assert.ok(COPY_SHIFT_FLOOR_BULK >= 2 * 0.070 && COPY_SHIFT_FLOOR_BULK < 0.299)
  assert.ok(BULK_SPREAD_MAX > 0.0249 && BULK_SPREAD_MAX < 0.0302,
    'the boundary sits between the largest bulk spread and the smallest amplified one')
}

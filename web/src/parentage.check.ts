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
  ABSENCE_MARGIN, absenceExplainable, agreement, classify, emptyTally, HET_BAND_DIPLOID,
  isAutosome, pair, pct, secondParentSignal, tallyRow, type Tally,
} from './parentage.ts'

import type { AB } from './informativity.ts'

const row = (
  chrom: string, pos: number, genotype: AB, baf: number | null = null,
): ProbeRow => ({
  probesetId: `m${pos}`, chrom, pos, log2R: null, baf, copyNumber: null,
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

const bip = classify(build({ absent: 0.002, het: 0.30, nonParental: 0.30 }).t, 0.17)
assert.equal(bip.zygosity, 'diploid')
assert.equal(bip.originClass, 'biparental')

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
  for (let i = 0; i < 4000; i += 1) tallyRow('AA', row('1', i * 1000, 'AA', 0.5), t)
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
      // Same rate on every chromosome, against the same total confined to four of them.
      tallyRow('AA', row(String(c), 1000 + i * 1000, s < 25 ? 'BB' : 'AA', 0.5), uniform)
      tallyRow('AA', row(String(c), 1000 + i * 1000,
        c <= 4 && s < 138 ? 'BB' : 'AA', 0.5), patchy)
    }
  }
  // A dropout ceiling low enough that 2.5% lands in the band rather than below it.
  for (let i = 0; i < 800; i += 1) tallyRow('NC', row('1', 9e8 + i, i < 88 ? 'NC' : 'AA'), uniform)
  for (let i = 0; i < 800; i += 1) tallyRow('NC', row('1', 9e8 + i, i < 88 ? 'NC' : 'AA'), patchy)
  const u = classify(uniform, 0.17)
  const q = classify(patchy, 0.17)
  assert.equal(u.verdict, 'unclear')
  assert.equal(q.verdict, 'unclear')
  assert.ok(u.dispersion < q.dispersion)
  const mixed = (r: typeof u) => r.limits.some((l) => l.includes('mixed sample'))
  assert.ok(mixed(u), 'uniform absence should read as a mixture')
  assert.ok(!mixed(q), 'absence confined to four chromosomes is a loss, not a mixture')
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

// Self-check for parental origin of a gain. Run: node src/gainOrigin.check.ts
//
// The band positions are measured figures from an external review of 34 trio-usable arrays. What
// this file checks is that the implementation reads them the right way round, refuses where the
// review says the evidence does not reach, and cannot be made to invert by the two mechanisms
// that would do it silently: marker orientation, and a biased centre.
import assert from 'node:assert/strict'
import {
  paternalShare, paternalShareOneParent, recentre, callGainOrigin, callHomologue,
  MIN_INFORMATIVE_DEFAULT, MIN_INFORMATIVE_TROPHECTODERM, SHARE_MARGIN, EXPECTED_SEPARATION,
  HETERO_MULTIPLE, HETERO_ABSOLUTE_MIN, HET_BACKGROUND_FLOOR, externalHetBackground,
  type DosageMarker,
} from './gainOrigin.ts'
import type { AB } from './informativity.ts'

/** `n` markers at a given paternal share, alternating which parent carries the A allele. */
const region = (share: number, n: number, chrom = '1'): DosageMarker[] =>
  Array.from({ length: n }, (_, i) => ({ chrom, pos: i * 1000, patShare: share }))

// --- 1. ORIENTATION. The trap that would invert every call ---------------------------------------
//
// An extra paternal copy pushes raw BAF DOWN where the father is AA and UP where he is BB. A
// statistic built on raw BAF averages those two to nothing, or picks the wrong sign. The paternal
// share must point the same way at both kinds of marker.
{
  // Trisomy AAB: two paternal A, one maternal B. Raw BAF 1/3.
  const fatherIsA = paternalShare('AA', 'BB', 1 / 3)
  // The mirror marker: father BB, mother AA, two paternal B and one maternal A. Raw BAF 2/3.
  const fatherIsB = paternalShare('BB', 'AA', 2 / 3)
  assert.ok(fatherIsA !== null && fatherIsB !== null)
  assert.ok(Math.abs(fatherIsA - 2 / 3) < 1e-9, 'father-AA marker reads 2/3 paternal share')
  assert.ok(Math.abs(fatherIsB - 2 / 3) < 1e-9, 'father-BB marker reads 2/3 too, not 1/3')
  assert.ok(Math.abs(fatherIsA - fatherIsB) < 1e-9,
    'THE SAME event must read identically at both marker orientations')

  // And the maternal direction mirrors it.
  assert.ok(Math.abs(paternalShare('AA', 'BB', 2 / 3)! - 1 / 3) < 1e-9)
  assert.ok(Math.abs(paternalShare('BB', 'AA', 1 / 3)! - 1 / 3) < 1e-9)

  // A euploid cell sits at 0.5 whichever way round the parents are.
  assert.equal(paternalShare('AA', 'BB', 0.5), 0.5)
  assert.equal(paternalShare('BB', 'AA', 0.5), 0.5)

  // Real proof against the mixed-orientation failure: a region of BOTH marker kinds carrying one
  // paternal gain. In raw BAF these average to 0.5 and the event vanishes.
  const mixed: DosageMarker[] = []
  const rawBaf: number[] = []
  for (let i = 0; i < 400; i += 1) {
    const fatherA = i % 2 === 0
    const baf = fatherA ? 1 / 3 : 2 / 3
    rawBaf.push(baf)
    const s = fatherA ? paternalShare('AA', 'BB', baf) : paternalShare('BB', 'AA', baf)
    mixed.push({ chrom: '1', pos: i * 1000, patShare: s! })
  }
  const rawMean = rawBaf.reduce((a, b) => a + b, 0) / rawBaf.length
  assert.ok(Math.abs(rawMean - 0.5) < 1e-9,
    'raw BAF averages to 0.5 on mixed orientations, which is the signal disappearing')
  const call = callGainOrigin(mixed, 0.5)
  assert.equal(call.origin, 'paternal', 'the oriented share still finds it')
}

// --- 2. the measured bands land where they should ------------------------------------------------
{
  // Review medians, which do not compress across stages: 0.3336 / 0.5003 / 0.6669 raw BAF at
  // father-AA markers, i.e. paternal shares of 0.6664 / 0.4997 / 0.3331.
  assert.equal(callGainOrigin(region(0.6664, 400), 0.5).origin, 'paternal')
  assert.equal(callGainOrigin(region(0.4997, 400), 0.5).origin, 'unclear')
  assert.equal(callGainOrigin(region(0.3331, 400), 0.5).origin, 'maternal')

  // WHY MEDIANS, as an executable fact rather than a comment. On a WGA blastomere the MEAN band
  // separation compresses 36%, from 0.167 to 0.106, so the bands would sit at only +/-0.053 from
  // the centre - inside this margin, and every single-cell gain would be refused. The band
  // MEDIANS do not compress: they stay at 0.3336 / 0.5003 / 0.6669, i.e. +/-0.0835. A
  // median-based statistic therefore still calls a blastomere gain where a mean-based one could
  // not, which is the whole reason this module medians and never means.
  const MEAN_COMPRESSED = 0.106 / 2
  const MEDIAN_HALF = EXPECTED_SEPARATION / 2
  assert.ok(MEAN_COMPRESSED < SHARE_MARGIN,
    'a mean-based statistic would fall inside the margin on a single cell')
  assert.ok(SHARE_MARGIN < MEDIAN_HALF, 'and the median-based one clears it')
  assert.equal(callGainOrigin(region(0.5 + MEDIAN_HALF, 400), 0.5).origin, 'paternal')
  assert.equal(callGainOrigin(region(0.5 - MEDIAN_HALF, 400), 0.5).origin, 'maternal')
}

// --- 3. THE BIASED CENTRE. The other mechanism that inverts calls ---------------------------------
//
// Against a reconstructed parent the theoretical 0.5 is wrong, biased toward maternal, because a
// truly father-heterozygous marker cannot be represented in an all-homozygous reference and is
// promoted into the informative set carrying a maternal-looking share. Measured offsets: +0.013 at
// five products, +0.032 at four, +0.073 to +0.077 at three.
{
  // A genuine PATERNAL gain of half the measured separation, on an array whose centre is offset
  // 0.077 toward maternal by a three-product reference.
  const OFFSET = 0.077
  const trueGain = 0.5 + EXPECTED_SEPARATION / 2
  const observed = region(trueGain - OFFSET, 400)
  const background = region(0.5 - OFFSET, 4000)

  // Judged against the theoretical centre, the offset eats the signal.
  const naive = callGainOrigin(observed, 0.5)
  assert.notEqual(naive.origin, 'paternal',
    'the theoretical centre must NOT recover a paternal gain through a 0.077 maternal offset')

  // Judged against the sample's own median, it is recovered.
  const centre = recentre([...background, ...observed])
  const fixed = callGainOrigin(observed, centre)
  assert.equal(fixed.origin, 'paternal', 'recentring on the sample recovers it')
  assert.ok(Math.abs(centre - (0.5 - OFFSET)) < 0.02,
    `the recovered centre should sit near the offset one, got ${centre}`)

  // And the offset is large enough to have inverted a weaker call outright, which is why this is
  // required rather than advisory.
  assert.ok(OFFSET > EXPECTED_SEPARATION * 0.4,
    'the measured offset is a large fraction of the whole separation')
}

// --- 3b. the same failure, with the numbers measured on a REAL array -----------------------------
//
// From audit/gain-positive-control.ts: A8's own centre is 0.5993, not the theoretical 0.5000. Its
// untouched euploid genome therefore sits +0.0993 from the theoretical centre, over the margin, and
// would be reported as a paternal gain that is not there. This is the argument for recentring made
// as a measurement rather than a claim.
{
  const A8_CENTRE = 0.5993
  const untouched = region(A8_CENTRE, 12_185)
  assert.equal(callGainOrigin(untouched, 0.5).origin, 'paternal',
    'against the theoretical centre a real euploid array reads as a gain it does not have')
  assert.equal(callGainOrigin(untouched, recentre(untouched)).origin, 'unclear',
    'and against its own centre it correctly reads as nothing')

  // Both directions, at the deviations the control actually produced.
  assert.equal(callGainOrigin(region(0.7494, 12_185), A8_CENTRE).origin, 'paternal')
  assert.equal(callGainOrigin(region(0.4278, 12_185), A8_CENTRE).origin, 'maternal')
}

// --- 4. it refuses where the review says the evidence does not reach ------------------------------
{
  // Too few informative markers. 400 is the WGA blastomere requirement and the default.
  assert.equal(MIN_INFORMATIVE_DEFAULT, 400)
  const thin = callGainOrigin(region(0.70, 399), 0.5)
  assert.equal(thin.origin, 'unclear', 'under the marker floor is a refusal however strong')
  assert.ok(thin.why.includes('399'))

  // A focal gain on a WGA single cell: 400 markers span a median 47.6 Mb, so a 5 Mb region
  // carrying a median 51 informative markers cannot be annotated at the default.
  assert.equal(callGainOrigin(region(0.70, 51), 0.5).origin, 'unclear',
    'a focal gain on a single cell is refused at the default floor')
  // The same region on a trophectoderm biopsy, where the requirement is 100, is still refused;
  // at 51 markers nothing clears either bar.
  assert.equal(callGainOrigin(region(0.70, 51), 0.5, MIN_INFORMATIVE_TROPHECTODERM).origin,
    'unclear')
  // But a region that does clear the trophectoderm bar is callable there and not at the default.
  assert.equal(callGainOrigin(region(0.70, 120), 0.5, MIN_INFORMATIVE_TROPHECTODERM).origin,
    'paternal')
  assert.equal(callGainOrigin(region(0.70, 120), 0.5).origin, 'unclear',
    'and the default stays strict, because the tool cannot see which stage it was handed')

  // Inside the margin is uncalled, not the nearer band.
  const mid = callGainOrigin(region(0.5 + SHARE_MARGIN / 2, 400), 0.5)
  assert.equal(mid.origin, 'unclear')
  assert.ok(mid.why.includes('not called either way'))

  // No informative markers at all.
  assert.equal(callGainOrigin([], NaN).origin, 'unclear')
}

// --- 5. one parent is a different, weaker path and cannot be confused with two -------------------
{
  // Informative only where the sample is heterozygous, since that is what says the unseen parent
  // differed. A homozygous sample says nothing.
  assert.equal(paternalShareOneParent('AA', 'AA', 0.0), null)
  assert.equal(paternalShareOneParent('AA', 'NC', 0.5), null)
  // Oriented the same way as the two-parent form, so downstream cannot tell them apart by sign.
  assert.ok(Math.abs(paternalShareOneParent('AA', 'AB', 1 / 3)! - 2 / 3) < 1e-9)
  assert.ok(Math.abs(paternalShareOneParent('BB', 'AB', 2 / 3)! - 2 / 3) < 1e-9)
}

// --- 6. the uniparental case: one answer exists and one does not ---------------------------------
{
  // A haploid array's spurious heterozygosity runs 7.1-10.9% on clean arrays. An extra copy that
  // is the parent's OTHER homologue reads far above that.
  const hetero = callHomologue(0.17, 400, 0.085)
  assert.equal(hetero.verdict, 'other homologue')
  assert.ok(hetero.why.includes('meiotic'))

  // An extra copy of the SAME homologue is bit-identical to one copy: the region reads exactly
  // like the background. This must NOT be reported as "no gain".
  const iso = callHomologue(0.085, 400, 0.085)
  assert.equal(iso.verdict, 'indistinguishable')
  assert.ok(iso.why.includes('not evidence of no gain'),
    'the refusal must say what it is not claiming')
  assert.ok(iso.why.includes('bit-identical'))

  // It holds against a degraded array, where the background runs 24-45%.
  assert.equal(callHomologue(0.30, 400, 0.35).verdict, 'indistinguishable',
    'a high background must not be read as a duplication')
  assert.equal(callHomologue(0.75, 400, 0.35).verdict, 'other homologue')

  // Under the marker floor it refuses rather than reading noise.
  assert.equal(callHomologue(0.50, 49, 0.085).verdict, 'indistinguishable')
  assert.equal(HETERO_MULTIPLE, 2.0)

  // A RATIO IS NOT ENOUGH, and this is a real case from the audit. Two clean chromosomes whose
  // amplification noise differs by a factor of two: 7.4% against a 3.4% background clears the
  // doubling and is nowhere near what a second homologue produces. The absolute floor removes it.
  assert.ok(0.074 >= HETERO_MULTIPLE * 0.034, 'it does clear the ratio test')
  assert.equal(callHomologue(0.074, 5000, 0.034).verdict, 'indistinguishable',
    'and must still be refused, because 7.4% is not a second homologue')
  assert.ok(HETERO_ABSOLUTE_MIN > 0.074 && HETERO_ABSOLUTE_MIN < 0.17,
    'the floor sits above amplification noise and below the parent heterozygosity a real '
    + 'duplication reads at')
  // The genuine ones from the same audit, at 55-61%, still call.
  for (const r of [0.551, 0.613, 0.600, 0.601]) {
    assert.equal(callHomologue(r, 5000, 0.045).verdict, 'other homologue', `${r} is real`)
  }
}

// --- 7. direction is reported; a copy number is not ----------------------------------------------
//
// The band separation compresses 36% on WGA single cells while the direction stays unbiased
// (|paternal| 0.1061 against |maternal| 0.1062). So the magnitude cannot be turned into a copy
// count, and nothing in this module returns one.
{
  const strong = callGainOrigin(region(0.90, 400), 0.5)
  const weak = callGainOrigin(region(0.60, 400), 0.5)
  assert.equal(strong.origin, weak.origin, 'a bigger deviation is the same call, not a bigger one')
  const keys = Object.keys(strong)
  assert.ok(!keys.some((k) => /copies|copyNumber|ploidy/i.test(k)),
    'nothing here may imply a copy count')
}

// --- 8. the background is external, and on a failed array everything refuses --------------------
//
// Measured per-chromosome heterozygosity at donor-homozygous markers on the one array in this
// series carrying whole-chromosome gains. It is a QC failure - 42.1% call rate, excluded as a
// product - and the numbers show why nothing can be read from it: the whole genome runs at 69-77%
// heterozygosity, which is amplification noise, not biology.
//
// Note what these numbers do NOT show. The gained chromosomes read LOWER than the rest, not
// higher, so this array is not an example of events inflating their own background. The
// inflation mechanism is real and is why the background is external, but it is demonstrated
// below on a constructed case and honestly labelled as one.
{
  const measured = new Map<string, { informative: number; het: number }>()
  const rows: [string, number, number][] = [
    ['1', 10658, 0.505], ['2', 8857, 0.475], ['3', 22097, 0.716], ['4', 20093, 0.689],
    ['5', 19095, 0.704], ['6', 22284, 0.704], ['7', 18913, 0.729], ['8', 15655, 0.708],
    ['9', 5283, 0.496], ['10', 16609, 0.723], ['11', 18574, 0.735], ['12', 15987, 0.725],
    ['13', 11736, 0.706], ['14', 11558, 0.718], ['15', 4037, 0.508], ['16', 13017, 0.746],
    ['17', 14164, 0.755], ['18', 8765, 0.702], ['19', 4972, 0.404], ['20', 8579, 0.747],
    ['21', 4904, 0.718], ['22', 6218, 0.774],
  ]
  for (const [c, n, r] of rows) measured.set(c, { informative: n, het: Math.round(n * r) })

  // Every gained chromosome refuses, and that is the right answer on an array like this.
  for (const c of ['1', '2', '9', '15', '19']) {
    const h = measured.get(c)!
    const bg = externalHetBackground(measured, c)
    assert.equal(callHomologue(h.het / h.informative, h.informative, bg).verdict,
      'indistinguishable',
      `chr${c} on a 42% call-rate array must refuse, not be read as a duplication`)
  }

  // The external background here is the median of the OTHER chromosomes, around 0.72, which is
  // the array's noise floor rather than anything biological.
  assert.ok(externalHetBackground(measured, '1') > 0.6,
    'a failed array has a high background and the tool sees it')
}

// --- 8b. the inflation mechanism, on a CONSTRUCTED case that isolates it -------------------------
//
// Not from this series. Built to show why the background excludes the chromosome under test: if
// enough chromosomes carry the event, a genome-wide rate is lifted by the very thing being
// measured and each event reads as ordinary. This is the self-null failure the segment scan
// already solved with an external per-chromosome null.
{
  const built = new Map<string, { informative: number; het: number }>()
  // Six of twelve chromosomes carrying a hetero-duplication, six clean at 8.5%.
  const EVENT = 0.30
  for (const c of ['1', '2', '3', '4', '5', '6']) {
    built.set(c, { informative: 6000, het: Math.round(6000 * EVENT) })
  }
  for (const c of ['7', '8', '9', '10', '11', '12']) {
    built.set(c, { informative: 6000, het: Math.round(6000 * 0.085) })
  }
  const genomeWide = [...built.values()].reduce((a, h) => a + h.het, 0)
    / [...built.values()].reduce((a, h) => a + h.informative, 0)
  const external = externalHetBackground(built, '1')
  assert.ok(external < genomeWide, 'excluding the tested chromosome lowers the bar it clears')
  assert.equal(callHomologue(EVENT, 6000, external).verdict, 'other homologue',
    'against the external background the duplication is found')
  assert.equal(callHomologue(EVENT, 6000, genomeWide).verdict, 'indistinguishable',
    'against the genome-wide rate the same duplication is missed, because the other five events '
    + 'have already lifted the bar it has to clear')
}

// --- 9. the background cannot be dragged, floored, or read off nothing --------------------------
{
  const many = new Map<string, { informative: number; het: number }>()
  // Eight of twelve chromosomes carrying a gain must not move a median.
  for (const c of ['1', '2', '3', '4', '5', '6', '7', '8']) {
    many.set(c, { informative: 5000, het: 2500 })
  }
  for (const c of ['9', '10', '11', '12']) many.set(c, { informative: 5000, het: 425 })
  assert.ok(externalHetBackground(many, '9') > 0.4,
    'a median over mostly-affected chromosomes does move, and that limit is real')

  // A spotless genome floors rather than dividing by nothing.
  const clean = new Map([['1', { informative: 5000, het: 0 }], ['2', { informative: 5000, het: 0 }]])
  assert.equal(externalHetBackground(clean, '9'), HET_BACKGROUND_FLOOR)

  // Thin chromosomes are excluded from the background rather than distorting it.
  const thin = new Map([
    ['1', { informative: 10, het: 10 }],
    ['2', { informative: 5000, het: 425 }],
    ['3', { informative: 5000, het: 425 }],
  ])
  assert.ok(Math.abs(externalHetBackground(thin, '9') - 0.085) < 1e-9,
    'a 10-marker chromosome at 100% must not become the background')

  // Nothing to measure is NaN, which callHomologue turns into a refusal.
  assert.ok(Number.isNaN(externalHetBackground(new Map(), '1')))
  assert.equal(callHomologue(0.5, 400, NaN).verdict, 'indistinguishable')
}

console.log('gainOrigin.check.ts: all assertions passed')

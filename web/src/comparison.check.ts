// Self-check for the optional feature comparison. Run: node src/comparison.check.ts
//
// Three properties carry this module.
//
//   THE CAVEAT MUST TRAVEL. This analysis is positional and says nothing about parental origin,
//   and that is its most tempting misreading: the fragile compartment is established on both
//   parental genomes from the first cell cycle. A result that could be shown without the caveat
//   attached is a result that will be.
//
//   TOO FEW REGIONS IS NOT A NEGATIVE RESULT. Reporting "no coincidence" from three regions reads
//   as evidence of independence when it is evidence of nothing.
//
//   THE THRESHOLD IS CORRECTED FOR HOW MANY FEATURES WERE TESTED. Testing every track against one
//   set of regions is otherwise a way of finding whichever one happens to fit.
import assert from 'node:assert/strict'
import {
  compare, meaningFor, methodsText, enrichmentBars, regionGrid, foldChart, NOT_ABOUT_PARENTS,
  comparisonRegions, listOf, normaliseTrack, nullHistograms, pooledRegions, pooledMarkers,
  regionFlags, relatedCount, RELATED_MIN_FOLD,
  ALPHA,
} from './comparison.ts'
import type { FeatureTrack, Region } from './features.ts'
import { readFileSync } from 'node:fs'

// A track whose fragile sites sit exactly where the test regions do, and whose long genes do not.
const track: FeatureTrack = {
  build: 'hg19',
  fragile: Array.from({ length: 12 }, (_, i) => ({ chrom: '1', startBp: i * 10e6, endBp: i * 10e6 + 9e6, name: `FRA1${i}` })),
  longGenes: [{ chrom: '1', startBp: 195e6, endBp: 199e6, name: 'FAR' }],
  gaps: [],
  geneDensityPerMb: {},
  lateReplicationValleysES: [],
  lateReplicationValleysConstitutive: [],
}
const regions: Region[] = Array.from({ length: 8 }, (_, i) => ({
  chrom: '1', startBp: i * 10e6 + 1e6, endBp: i * 10e6 + 4e6,
}))
const markers = new Map<string, number[]>([['1',
  Array.from({ length: 4000 }, (_, i) => i * 50_000)]])

// --- 1. THE CAVEAT TRAVELS ON EVERY RESULT --------------------------------------------------------
{
  const c = compare(track, regions, markers, { permutations: 200 })
  assert.equal(c.caveat, NOT_ABOUT_PARENTS)
  assert.ok(c.caveat.includes('never be read as support for a parental call'))
  assert.ok(c.caveat.includes('both parental genomes'),
    'the REASON must travel too, not just the prohibition: a reader who knows the fragile '
    + 'compartment is biparental does not need to be told twice')

  // An empty run still carries it, because that is exactly when a reader improvises.
  assert.equal(compare(track, [], markers, { permutations: 200 }).caveat, NOT_ABOUT_PARENTS)
}

// --- 2. TOO FEW REGIONS IS NOT A NEGATIVE RESULT ---------------------------------------------------
{
  const few = compare(track, regions.slice(0, 3), markers, { permutations: 200 })
  assert.equal(few.verdict, 'underpowered')
  assert.ok(few.headline.includes('No conclusion either way'),
    `an underpowered run must not read as independence: ${few.headline}`)
  assert.ok(!few.headline.includes('do not coincide'))

  const none = compare(track, [], markers, { permutations: 200 })
  assert.equal(none.verdict, 'underpowered')
  assert.equal(none.features.length, 0)
}

// --- 3. THE THRESHOLD IS CORRECTED FOR THE NUMBER OF FEATURES TESTED -------------------------------
{
  const c = compare(track, regions, markers, { permutations: 400 })
  assert.ok(c.features.length >= 2, 'the track carries several feature sets')
  // The methods text must state the corrected bound and the count it was divided by, so a reader
  // can check the correction by counting the rows in front of them.
  assert.ok(c.methods.includes('0.05 divided by'))
  assert.ok(c.methods.includes(`${c.features.length} feature set`))
  assert.ok(c.methods.includes('SAME chromosome'), 'the null must be described, not just named')
  assert.ok(c.methods.includes(String(400)) || c.methods.includes('400'),
    'and the permutation count must be the one actually used')

  // Generated from the run: a different run must produce different methods.
  assert.notEqual(methodsText(8, 400, 5, 0.01), methodsText(8, 2000, 5, 0.01))
}

// --- 4. THE VERDICT FOLLOWS THE FEATURES, AND SAYS WHAT TO DO --------------------------------------
{
  const c = compare(track, regions, markers, { permutations: 400 })
  assert.ok(['feature-coincident', 'independent'].includes(c.verdict))
  if (c.verdict === 'feature-coincident') {
    assert.ok(c.headline.includes('mundane explanation'),
      'a coincidence must be phrased as an alternative explanation rather than a refutation')
    assert.ok(c.features.some((f) => f.significant))
  } else {
    assert.ok(c.headline.includes('leaves their detection statistics standing'))
    assert.ok(!c.features.some((f) => f.significant))
  }
  // Every feature carries what a coincidence with it would MEAN, or the table is unreadable.
  for (const f of c.features) {
    assert.ok(f.means.length > 20, `${f.feature} must say what coinciding with it implies`)
    assert.ok(f.label.length > 0)
  }
  assert.ok(meaningFor('common fragile site').means.includes('replication stress'))
  assert.ok(meaningFor('not-a-real-track').label === 'not-a-real-track',
    'an unknown track must still be presentable rather than throwing')

  // THE MEANING TABLE MUST BE KEYED ON WHAT scoreAll ACTUALLY EMITS. It was keyed on the track's
  // FIELD names, which looked right and matched nothing: every feature fell through to the generic
  // fallback and the table lost the column that makes it readable. Nothing failed, because a
  // fallback is not an error. This is the assertion that would have caught it.
  for (const f of c.features) {
    assert.notEqual(f.label, f.feature,
      `"${f.feature}" has no entry in FEATURE_MEANING, so it fell through to the raw name. The `
      + 'keys are the names scoreAll emits, not the track object\'s field names')
    assert.ok(f.means.length > 40, `"${f.feature}" fell through to the generic meaning`)
  }
}

// --- 5. THE CHARTS SHARE ONE AXIS, WHICH IS WHY THEY ARE READABLE ----------------------------------
{
  const c = compare(track, regions, markers, { permutations: 400 })
  const bars = enrichmentBars(c)
  assert.equal(bars.length, c.features.length)
  const axes = new Set(bars.map((b) => b.axisMax))
  assert.equal(axes.size, 1,
    'ONE axis across every row. Per-row axes make every feature fill its own width and look '
    + 'equally enriched, which is the easiest way to draw a chart that says the opposite of its data')
  assert.ok([...axes][0] >= Math.max(...bars.map((b) => Math.max(b.observed, b.hi))))

  // The grid is the plain answer to "overlap or lack thereof", checkable against the enrichment.
  const names = ['chr1 1.0-4.0Mb', 'chr1 11.0-14.0Mb']
  const grid = regionGrid(c, names)
  assert.equal(grid.length, names.length * c.features.length)
  assert.ok(grid.every((g) => typeof g.touched === 'boolean'))

  // THE GRID MUST AGREE WITH THE ENRICHMENT BESIDE IT. It read region names out of `hits`, which
  // holds FEATURE names, so it matched nothing and drew an empty matrix while the enrichment on
  // the same data reported coincidence. Two panels of one figure contradicting each other is worse
  // than either panel alone.
  const full = compare(track, regions, markers, { permutations: 300 })
  for (const f of full.features) {
    assert.equal(f.regionHits.length, regions.length,
      'every region needs an entry, indexed as it was passed in')
    const share = f.regionHits.filter(Boolean).length / regions.length
    assert.ok(Math.abs(share - f.observed) < 1e-9,
      `panel D says ${share} of regions touch ${f.feature} while panel A says ${f.observed}. `
      + 'These are the same quantity and must be computed from the same overlap')
  }

  const fold = foldChart(c)
  assert.ok(fold.every((f) => Number.isFinite(f.fold)), 'an undefined fold must be dropped, not drawn')
  for (let i = 1; i < fold.length; i += 1) assert.ok(fold[i - 1].fold >= fold[i].fold, 'sorted')
}

// --- 6. THE NUMBERS ARE THE ENRICHMENT'S OWN, IN ITS OWN UNITS -------------------------------------
{
  const c = compare(track, regions, markers, { permutations: 400 })
  for (const f of c.features) {
    assert.ok(f.observed >= 0 && f.observed <= 1, 'observed is a FRACTION of regions, as shipped')
    assert.equal(f.observedCount, Math.round(f.observed * c.regions),
      'and the count beside it must be that fraction, not a second measurement')
    assert.ok(f.nullLo <= f.nullHi, 'the null band must not be inverted')
  }
  assert.equal(ALPHA, 0.05)
}

// --- 7. THE SHIPPED TRACK MUST ACTUALLY MATCH SOMETHING -------------------------------------------
//
// THIS IS THE TEST THAT WAS MISSING, and its absence hid a live defect for the life of the feature.
// hg19_features.json stores each interval as ["1", 61300000, 84900000, "FRA1B"], while every
// consumer reads f.chrom and f.startBp. On a tuple those are undefined, so the overlap test
// compared undefined against a chromosome name and returned false every time. The enrichment
// reported ZERO overlap with every feature on every real run, whatever the regions overlapped, and
// attached a p value to it.
//
// Nothing failed, because every existing test built intervals as objects: the suite exercised a
// shape the application never sees. So this one reads the file the application actually loads.
{
  const raw = JSON.parse(readFileSync(
    new URL('../public/hg19_features.json', import.meta.url), 'utf8'))
  const real = normaliseTrack(raw)
  assert.ok(real, 'the shipped track must normalise')
  assert.ok(real!.fragile.length > 10, 'and carry its fragile sites')

  for (const f of real!.fragile.slice(0, 5)) {
    assert.equal(typeof f.chrom, 'string', 'a decoded interval must have a chromosome')
    assert.ok(Number.isFinite(f.startBp) && Number.isFinite(f.endBp),
      'and finite coordinates, which a tuple read as an object does not')
    assert.ok(f.endBp > f.startBp)
  }

  // The end-to-end property: regions placed ON real fragile sites must be found, and the SAME
  // regions moved away must not be. Either half alone can pass with a broken track.
  const markers = new Map<string, number[]>()
  for (const c of ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']) {
    markers.set(c, Array.from({ length: 6000 }, (_, i) => i * 40_000))
  }
  const on = real!.fragile.filter((f) => markers.has(f.chrom)).slice(0, 12)
    .map((f) => ({ chrom: f.chrom, startBp: f.startBp + 5e5, endBp: f.startBp + 3e6 }))
  assert.ok(on.length >= 8, 'enough real fragile sites to test against')

  const hit = compare(real!, on, markers, { permutations: 400 })
  const fra = hit.features.find((f) => f.feature === 'common fragile site')!
  assert.ok(fra.observed > 0.9,
    `regions placed ON fragile sites must overlap them, got ${fra.observed}. Zero here is the `
    + 'defect this test exists for: a track that matches nothing reports no overlap with a p value')
  assert.ok(fra.significant, 'and the enrichment must clear the corrected threshold')

  const away = on.map((r) => ({ chrom: r.chrom, startBp: 2e6, endBp: 4.5e6 }))
  const miss = compare(real!, away, markers, { permutations: 400 })
  const fraAway = miss.features.find((f) => f.feature === 'common fragile site')!
  assert.ok(fraAway.observed < 0.2,
    'and the same regions moved away must NOT overlap them, or the test is passing on a statistic '
    + 'that says yes to everything')
}

// --- 8. FINDINGS ARE CHROMOSOMAL CHANGES TOO ------------------------------------------------------
//
// Taking only the segments left the panel reporting "0 regions" on a run whose changes were all
// taxonomy findings, while the defect list beside it said 23. A copy-neutral event and an isodisomy
// are changes exactly as a deletion is, and the question here is identical for all of them.
{
  const coords = (sg: unknown) => {
    const x = sg as { startBp: number; endBp: number }
    return { start: x.startBp, end: x.endBp }
  }
  const run = {
    segments: [{ chrom: '4', startBp: 1e6, endBp: 9e6 }] as never,
    findings: [
      { chrom: '7', startBp: 5e6, endBp: 25e6 },
      { chrom: '9', startBp: 2e6, endBp: 14e6 },
      // Genome-wide: a triploidy has no interval, so it cannot be compared against a place.
      { chrom: 'genome', startBp: 0, endBp: 0 },
    ],
  }
  const rs = comparisonRegions(run, coords)
  assert.equal(rs.length, 3, 'one segment plus two interval findings, and NOT the genome-wide one')
  assert.ok(rs.every((r) => r.region.chrom !== 'genome'),
    'a whole-genome finding has no place to compare, so it must be dropped rather than compared '
    + 'against coordinates it does not have')
  assert.ok(rs.some((r) => r.region.chrom === '7' && r.region.startBp === 5e6),
    'a taxonomy finding must be included')
  assert.ok(rs.some((r) => r.region.chrom === '4'), 'and so must a segment')
  assert.ok(rs.every((r) => /^chr\d+ [\d.]+-[\d.]+Mb$/.test(r.name)), `names: ${rs.map((r) => r.name)}`)

  // A run with findings and no segments must still have something to compare, which is the case
  // that was returning zero.
  const onlyFindings = comparisonRegions({ findings: run.findings }, coords)
  assert.equal(onlyFindings.length, 2)
}

// --- 9. THE HEADLINE READS AS A SENTENCE ----------------------------------------------------------
//
// Five coincidences joined with "and" between every pair reads as one long chain, which is what a
// real run on the bundled examples produced.
{
  assert.equal(listOf([]), '')
  assert.equal(listOf(['a']), 'a')
  assert.equal(listOf(['a', 'b']), 'a and b')
  assert.equal(listOf(['a', 'b', 'c']), 'a, b and c')
  assert.equal(listOf(['a', 'b', 'c', 'd', 'e']), 'a, b, c, d and e')
}

// --- 10. THE PERMUTATION NULL IS DRAWABLE, AND SAYS WHEN THE OBSERVATION IS OUTSIDE IT -----------
//
// A p value is a summary of this distribution, and summaries of it mislead both ways: a ratio near
// one reaches a small p when the null is tight, and a large ratio reaches nothing when it is broad.
// Drawing the null is what lets a reader judge rather than take the p on trust.
{
  const c = compare(track, regions, markers, { permutations: 400 })
  const hs = nullHistograms(c)
  assert.ok(hs.length, 'the enrichment already computes the distribution, so it must be drawable')
  for (const h of hs) {
    assert.ok(h.bins.length > 1, 'a histogram needs bins')
    assert.ok(h.observedAt >= 0 && h.observedAt <= 1, 'the marker must stay on the axis')
    assert.ok(h.hi >= h.lo)
    assert.ok(h.observedBin === -1 || (h.observedBin >= 0 && h.observedBin < h.bins.length))
  }

  // AN OBSERVATION THE PERMUTATION NEVER REACHED is the strongest result here, so it is reported
  // as outside rather than clamped into the end bin and drawn as if it were merely extreme.
  const outside = nullHistograms({
    ...c,
    features: [{ ...c.features[0], observed: 99, nullHist: { lo: 0, hi: 1, counts: [1, 2, 3] } }],
  })[0]
  assert.equal(outside.observedBin, -1, 'outside the null must be flagged, not clamped')
  assert.equal(outside.observedAt, 1, 'and still drawable, pinned at the edge')
}

// --- 11. THE RUN IS POOLED, BECAUSE ONE CHIP IS NOT A DATASET -------------------------------------
{
  const coords = (sg: unknown) => {
    const x = sg as { startBp: number; endBp: number }
    return { start: x.startBp, end: x.endBp }
  }
  const entries = [
    { file: { name: 'A_sample.csv.gz' },
      result: { segments: [{ chrom: '1', startBp: 1e6, endBp: 9e6 }] as never,
        markerPositions: new Map([['1', [1, 2, 3]]]) } },
    { file: { name: 'B_sample.csv.gz' },
      result: { findings: [{ chrom: '2', startBp: 5e6, endBp: 20e6 }],
        markerPositions: new Map([['1', [2, 3, 4]], ['2', [9, 8]]]) } },
    { file: { name: 'C_unrun.csv.gz' } },
  ]
  const pooled = pooledRegions(entries as never, coords)
  assert.equal(pooled.length, 2, 'both samples contribute; the unrun one does not')
  assert.ok(pooled[0].name.startsWith('A_sample'), `names must carry the sample: ${pooled[0].name}`)
  assert.ok(pooled[1].name.startsWith('B_sample'))

  // Marker positions are UNIONED, not concatenated: two arrays of one platform share almost every
  // marker, and counting each twice would inflate the density the null is matched on.
  const m = pooledMarkers(entries as never)
  assert.deepEqual(m.get('1'), [1, 2, 3, 4], 'shared positions must appear once, sorted')
  assert.deepEqual(m.get('2'), [8, 9])
}

// --- 12. THE STAR IS NARROW, AND THAT IS THE WHOLE POINT ------------------------------------------
//
// A region is flagged only where it overlaps a class that CLEARED THE CORRECTED THRESHOLD, not one
// it merely touches. Late-replication valleys cover most of the genome, so flagging on contact
// would star nearly every region and mean nothing. The flag says the region sits in a class whose
// coincidence with this set of regions is more than the matched null produces.
{
  const c = compare(track, regions, markers, { permutations: 400 })
  const flags = regionFlags(c)
  assert.equal(flags.length, regions.length, 'one flag per region, in order')

  for (const f of flags) {
    // Everything in `related` must be significant AND overlapped. Neither alone is enough.
    for (const r of f.related) {
      const feat = c.features.find((x) => x.label === r)!
      assert.ok(feat.significant, `${r} is starred but did not clear the threshold`)
      assert.ok(feat.regionHits[f.index], `${r} is starred on a region it does not overlap`)
    }
    // `overlaps` is the wider set and must contain `related`.
    for (const r of f.related) {
      assert.ok(f.overlaps.includes(r), 'a related class must also be an overlapped one')
    }
    assert.ok(f.overlaps.length >= f.related.length)
  }

  // TOUCHING IS NOT ENOUGH. If any feature is overlapped but not significant, no region may be
  // starred on its account, which is what stops the star meaning "touched something".
  const insig = c.features.find((f) => f.testable && !f.significant && f.regionHits.some(Boolean))
  if (insig) {
    assert.ok(flags.every((f) => !f.related.includes(insig.label)),
      `${insig.label} is overlapped but not significant, so it must never star a region`)
  }

  // AND SIGNIFICANCE IS NOT ENOUGH EITHER, which shipping it that way proved: 27 of 27 regions
  // were starred on a real run. Significance is a statement about the SET, and it can hold for a
  // class so broad that the null already expects three regions in four to overlap it. Then a single
  // region overlapping it is the default outcome, not a finding. An effect-size floor is what makes
  // the star mean what a reader takes it to mean.
  for (const f of flags) {
    for (const r of f.related) {
      const feat = c.features.find((x) => x.label === r)!
      assert.ok(feat.fold >= RELATED_MIN_FOLD,
        `${r} stars a region at ${feat.fold.toFixed(2)}x, under the ${RELATED_MIN_FOLD}x floor. `
        + 'A class the null already expects to be hit says nothing about one region')
    }
  }
  const broad = c.features.find((f) => f.testable && f.significant && f.fold < RELATED_MIN_FOLD)
  if (broad) {
    assert.ok(flags.every((f) => !f.related.includes(broad.label)),
      `${broad.label} is significant but only ${broad.fold.toFixed(2)}x, so it must not star`)
  }
  // The floor must be a real bar, not a formality.
  assert.ok(RELATED_MIN_FOLD >= 2, 'at least twice the chance rate')

  assert.equal(relatedCount(c), flags.filter((f) => f.related.length).length)
  assert.ok(relatedCount(c) <= regions.length)
}

console.log('comparison.check.ts: all assertions passed, including the parental caveat travelling '
  + 'on every result, too few regions refusing to read as independence, and one shared axis across '
  + 'the enrichment chart')

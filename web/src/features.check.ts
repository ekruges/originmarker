// Self-check for the positional feature test. Run: node src/features.check.ts
//
// The thing worth testing here is not whether an overlap is computed correctly, which is two
// comparisons, but whether the NULL is the marker-matched one. A uniform null is the failure this
// module exists to avoid, and it is invisible in the output: it produces a confident p from a
// biased comparison. So the cases below construct a genome where marker density and feature
// position are deliberately confounded, and require the module to report no enrichment.
import assert from 'node:assert/strict'
import {
  enrichment, enrichmentBy, scoreAll, markersIn, geneDensity, MIN_REGIONS,
  buildIndex, countIn, coveredBp, overlapsIndexed,
  type FeatureTrack, type FeatureInterval, type Region,
} from './features.ts'

const feat = (chrom: string, startBp: number, endBp: number, name: string): FeatureInterval =>
  ({ chrom, startBp, endBp, name })

// --- 1. markersIn counts, and counts the right half-open interval -------------------------------
{
  const m = [10, 20, 30, 40, 50]
  assert.equal(markersIn(m, 20, 40), 3, 'inclusive of both ends')
  assert.equal(markersIn(m, 21, 39), 1)
  assert.equal(markersIn(m, 0, 100), 5)
  assert.equal(markersIn(m, 60, 70), 0)
}

// --- 2. the case this module exists for: density confounded with the feature --------------------
//
// Markers exist ONLY in the first 10 Mb, and the feature also sits there. Every real region is in
// that window because nowhere else is callable, so a uniform null would find a large enrichment.
// Matched on marker count, the null lands in the same window too, and the answer is no enrichment.
{
  const markers: number[] = []
  for (let p = 0; p < 10_000_000; p += 1000) markers.push(p)
  const byChrom = new Map([['1', markers]])
  const features = [feat('1', 0, 10_000_000, 'DENSE')]
  const regions: Region[] = Array.from({ length: 20 }, (_, i) => ({
    chrom: '1', startBp: i * 400_000, endBp: i * 400_000 + 300_000,
  }))

  const e = enrichment('dense', features, regions, byChrom, 300)
  assert.equal(e.observed, 1, 'every region touches the feature, by construction')
  assert.ok(e.nullMean > 0.95,
    `a matched null must also land in the callable window, got ${e.nullMean}`)
  assert.ok(e.ratio < 1.05, `no enrichment may be reported here, got ratio ${e.ratio}`)
  assert.ok(e.p > 0.05, `and it must not be significant, got p ${e.p}`)
}

// --- 3. a real enrichment is still found ---------------------------------------------------------
//
// Markers uniform across the whole chromosome, feature in a small part of it, and every region
// placed inside the feature. The matched null now has somewhere else to go, so this must fire.
{
  const markers: number[] = []
  for (let p = 0; p < 100_000_000; p += 1000) markers.push(p)
  const byChrom = new Map([['1', markers]])
  const features = [feat('1', 0, 5_000_000, 'HOT')]
  const regions: Region[] = Array.from({ length: 20 }, (_, i) => ({
    chrom: '1', startBp: i * 200_000, endBp: i * 200_000 + 100_000,
  }))

  const e = enrichment('hot', features, regions, byChrom, 300)
  assert.equal(e.observed, 1)
  assert.ok(e.nullMean < 0.2, `null must be low when the feature is 5% of the genome, got ${e.nullMean}`)
  assert.ok(e.ratio > 4, `a real enrichment must be reported, got ${e.ratio}`)
  assert.ok(e.p < 0.05, `and be significant, got ${e.p}`)
  assert.deepEqual(e.hits, ['HOT'], 'the feature hit is named, not just counted')
}

// --- 4. depletion is reportable, not silently one-sided ------------------------------------------
{
  const markers: number[] = []
  for (let p = 0; p < 100_000_000; p += 1000) markers.push(p)
  const byChrom = new Map([['1', markers]])
  // feature covers most of the chromosome; regions all avoid it
  const features = [feat('1', 20_000_000, 100_000_000, 'MOST')]
  const regions: Region[] = Array.from({ length: 20 }, (_, i) => ({
    chrom: '1', startBp: i * 500_000, endBp: i * 500_000 + 200_000,
  }))
  const e = enrichment('most', features, regions, byChrom, 300)
  assert.equal(e.observed, 0, 'regions avoid the feature by construction')
  assert.ok(e.nullMean > 0.5, 'but a matched null lands in it most of the time')
  assert.ok(e.p < 0.05, 'so depletion must be significant rather than ignored')
}

// --- 5. it refuses on too few regions rather than reporting a rate -------------------------------
{
  const markers = Array.from({ length: 1000 }, (_, i) => i * 1000)
  const byChrom = new Map([['1', markers]])
  const few: Region[] = [{ chrom: '1', startBp: 0, endBp: 10_000 }]
  const e = enrichment('x', [feat('1', 0, 1_000_000, 'F')], few, byChrom, 100)
  assert.ok(Number.isNaN(e.p), 'no p from one region')
  assert.ok(e.why.includes(String(MIN_REGIONS)), 'and it says what it needed')
}

// --- 6. the same inputs give the same p ----------------------------------------------------------
{
  const markers = Array.from({ length: 5000 }, (_, i) => i * 10_000)
  const byChrom = new Map([['1', markers]])
  const features = [feat('1', 0, 10_000_000, 'F')]
  const regions: Region[] = Array.from({ length: 10 }, (_, i) => ({
    chrom: '1', startBp: i * 3_000_000, endBp: i * 3_000_000 + 500_000,
  }))
  const a = enrichment('f', features, regions, byChrom, 200)
  const b = enrichment('f', features, regions, byChrom, 200)
  assert.equal(a.p, b.p, 'a reported p must be reproducible from the same inputs')
}

// --- 7. gene density reads the binned track ------------------------------------------------------
{
  const track = {
    build: 'hg19', fragile: [], longGenes: [], gaps: [],
    lateReplicationValleysES: [], lateReplicationValleysConstitutive: [],
    geneDensityPerMb: { '1': { '0': 10, '1': 20 } },
  } as unknown as FeatureTrack
  assert.equal(geneDensity(track, { chrom: '1', startBp: 0, endBp: 1_500_000 }), 15)
  assert.ok(Number.isNaN(geneDensity(track, { chrom: '9', startBp: 0, endBp: 1000 })))
  assert.equal(scoreAll(track, [], new Map()).length, 5, 'five tracks scored')
}

// --- 8. THE NULL IS NOT FORKED --------------------------------------------------------------
//
// This is the assertion the widened track set rests on. enrichmentBy() adds density and coverage
// statistics and an interval index for speed; if any of that perturbed the permutation - a
// different draw order, a re-seeded generator, a differently rounded marker count - then every
// number it reports would be incomparable with the validated figures, and the difference would be
// INVISIBLE in the output. So mode 'binary' must reproduce enrichment() EXACTLY, to the last bit
// of the p value, on inputs where the answer is non-trivial.
{
  const markers: number[] = []
  for (let p = 0; p < 100_000_000; p += 997) markers.push(p)
  const byChrom = new Map([['1', markers], ['2', markers.map((m) => m + 13)]])
  const features = [
    feat('1', 3_000_000, 9_000_000, 'A'), feat('1', 40_000_000, 41_000_000, 'B'),
    feat('2', 0, 2_000_000, 'C'), feat('1', 8_500_000, 8_600_000, 'D'),
  ]
  const regions: Region[] = Array.from({ length: 30 }, (_, i) => ({
    chrom: i % 3 === 0 ? '2' : '1',
    startBp: i * 1_700_000,
    endBp: i * 1_700_000 + 900_000,
  }))
  const a = enrichment('t', features, regions, byChrom, 500)
  const b = enrichmentBy('t', features, regions, byChrom, 'binary', 500)
  assert.equal(b.observed, a.observed, 'binary mode must reproduce the observed rate exactly')
  assert.equal(b.nullMean, a.nullMean, 'and the SAME null mean: a different null is a fork')
  assert.equal(b.p, a.p, 'and therefore the identical p')
  assert.deepEqual(b.nullQuantiles, a.nullQuantiles, 'and the identical null distribution')
}

// --- 9. the interval index answers the same questions as a linear scan ------------------------
//
// The index exists only because the widened set carries ~500,000 G4 intervals and a linear scan
// per permutation is 10^11 comparisons. A fast wrong answer is worse than a slow one, so it is
// checked against brute force on overlapping, nested and abutting intervals.
{
  const fs = [
    feat('1', 100, 200, 'a'), feat('1', 150, 400, 'b'), feat('1', 160, 170, 'nested'),
    feat('1', 400, 500, 'abut'), feat('1', 900, 1000, 'far'), feat('2', 0, 50, 'other'),
  ]
  const idx = buildIndex(fs)
  const brute = (c: string, s: number, e: number) =>
    fs.filter((f) => f.chrom === c && f.startBp < e && s < f.endBp)
  for (const [s, e] of [[0, 50], [90, 110], [150, 165], [199, 201], [395, 405], [500, 900],
                        [0, 2000], [1000, 1100], [999, 1000]] as [number, number][]) {
    const want = brute('1', s, e)
    assert.equal(countIn(idx, '1', s, e), want.length, `count on [${s},${e})`)
    assert.equal(overlapsIndexed(idx, '1', s, e), want.length > 0, `overlap on [${s},${e})`)
    // coverage counts each base once however many intervals cover it
    const covered = new Set<number>()
    for (const f of want) {
      for (let p = Math.max(s, f.startBp); p < Math.min(e, f.endBp); p += 1) covered.add(p)
    }
    assert.equal(coveredBp(idx, '1', s, e), covered.size, `coverage on [${s},${e})`)
  }
  assert.equal(countIn(idx, '9', 0, 1000), 0, 'a chromosome with no features is not an error')
}

// --- 10. density and coverage detect what binary overlap CANNOT ------------------------------
//
// The failure this whole mode system exists to prevent. A track of small intervals spread evenly
// over the chromosome is touched by every large region AND by every matched null interval, so
// binary overlap reports 1.00 and hides a real difference. Here the regions genuinely carry
// TWICE the interval density of the rest of the chromosome; binary must miss it and density must
// find it.
{
  const markers: number[] = []
  for (let p = 0; p < 60_000_000; p += 500) markers.push(p)
  const byChrom = new Map([['1', markers]])
  const features: FeatureInterval[] = []
  // baseline: one interval every 100 kb across the chromosome
  for (let p = 0; p < 60_000_000; p += 100_000) features.push(feat('1', p, p + 500, 'bg'))
  // the observed regions sit in the first 12 Mb, which carries an extra interval per 100 kb
  for (let p = 0; p < 12_000_000; p += 100_000) features.push(feat('1', p + 50_000, p + 50_500, 'hot'))
  const regions: Region[] = Array.from({ length: 12 }, (_, i) => ({
    chrom: '1', startBp: i * 1_000_000, endBp: i * 1_000_000 + 800_000,
  }))

  const bin = enrichmentBy('dense', features, regions, byChrom, 'binary', 400)
  assert.equal(bin.observed, 1, 'every region contains an interval, by construction')
  assert.ok(bin.ratio > 0.99 && bin.ratio < 1.01,
    `binary overlap must be blind here, got ratio ${bin.ratio}`)

  const den = enrichmentBy('dense', features, regions, byChrom, 'density', 400)
  assert.ok(den.ratio > 1.6 && den.ratio < 2.4,
    `density must recover the two-fold difference binary missed, got ${den.ratio}`)
  assert.ok(den.p < 0.05, `and call it, got p ${den.p}`)
}

// --- 11. an extra track is scored through scoreAll under its declared mode --------------------
{
  const track = {
    build: 'hg19', fragile: [], longGenes: [], gaps: [],
    lateReplicationValleysES: [], lateReplicationValleysConstitutive: [],
    geneDensityPerMb: {},
    extra: { imprintedDomain: [feat('1', 0, 1_000_000, 'H19')], g4seqK: [] },
    modes: { imprintedDomain: 'binary' as const, g4seqK: 'coverage' as const },
  } as unknown as FeatureTrack
  const out = scoreAll(track, [], new Map())
  assert.equal(out.length, 7, 'five built-in tracks plus the two extra ones')
  assert.deepEqual(out.slice(5).map((e) => e.feature), ['g4seqK', 'imprintedDomain'],
    'extra tracks are scored in a stable (sorted) order')
}

console.log('features.check.ts: all assertions passed, including the confounded-density case')

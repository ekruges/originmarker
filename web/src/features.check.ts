// Self-check for the positional feature test. Run: node src/features.check.ts
//
// The thing worth testing here is not whether an overlap is computed correctly, which is two
// comparisons, but whether the NULL is the marker-matched one. A uniform null is the failure this
// module exists to avoid, and it is invisible in the output: it produces a confident p from a
// biased comparison. So the cases below construct a genome where marker density and feature
// position are deliberately confounded, and require the module to report no enrichment.
import assert from 'node:assert/strict'
import {
  enrichment, scoreAll, markersIn, geneDensity, MIN_REGIONS,
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

console.log('features.check.ts: all assertions passed, including the confounded-density case')

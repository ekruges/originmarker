// Self-check for the segment scan. Run: node src/segments.check.ts
//
// Every figure quoted here was measured on GSE148488, on real genomes rather than simulated ones:
// the null on five paternal pronuclei of the sperm donor, which are present on every autosome, and
// the positive class by splicing a maternal pronucleus of the same series into one of them so the
// segment carries real amplification artefact.
import assert from 'node:assert/strict'
import {
  MIN_SEGMENT_MARKERS, NULL_FLOOR, SEGMENT_LRT, SEGMENT_LRR_SHIFT, externalNull,
  scanChromosome, scanCopyNumber,
  type MarkerAbsence,
} from './segments.ts'

/** Markers on one chromosome, with absence planted over a given index range. */
const chrom = (n: number, from = -1, to = -1, rate = 0.11, bg = 0.005): MarkerAbsence[] =>
  Array.from({ length: n }, (_, i) => {
    const inSeg = i >= from && i <= to
    // Deterministic, and spread rather than contiguous, so nothing here depends on run structure.
    const hit = ((i * 7919) % 1000) < (inSeg ? rate : bg) * 1000
    return { chrom: '2', pos: 1_000_000 + i * 5_000, absent: hit }
  })

// --- 1. the null excludes the chromosome under test ---------------------------------------------
//
// A whole-chromosome event otherwise raises the bar it has to clear, which is how a self-null
// scan comes back negative on the largest events in the sample.
{
  const by = new Map<string, [number, number]>([
    ['1', [50_000, 250]], ['2', [50_000, 25_000]], ['3', [50_000, 250]], ['4', [50_000, 250]],
  ])
  assert.ok(Math.abs(externalNull(by, '2') - 0.005) < 1e-9,
    'chr2 is half absent and must not appear in its own null')

  // Exclusion bites once the event chromosome would sit at or past the median. With two of three
  // affected, scanning one of them against a null that still contains it more than doubles the
  // bar it has to clear, which is how a self-null scan comes back negative on its largest events.
  const two = new Map<string, [number, number]>([
    ['1', [50_000, 250]], ['2', [50_000, 25_000]], ['3', [50_000, 25_000]],
  ])
  assert.ok(externalNull(two, '2') < externalNull(two, '9') / 1.9,
    'excluding the chromosome under test must lower the bar it is judged against')

  // Several events cannot move a median the way they move a mean.
  const many = new Map<string, [number, number]>([
    ['1', [50_000, 25_000]], ['2', [50_000, 25_000]], ['3', [50_000, 250]],
    ['4', [50_000, 250]], ['5', [50_000, 250]], ['6', [50_000, 250]], ['7', [50_000, 250]],
  ])
  assert.ok(Math.abs(externalNull(many, '9') - 0.005) < 1e-9,
    'two whole-chromosome events must not drag the median')

  // A clean genome would otherwise divide by nearly zero and make every window catastrophic.
  const clean = new Map<string, [number, number]>([['1', [50_000, 0]], ['2', [50_000, 0]]])
  assert.equal(externalNull(clean, '9'), NULL_FLOOR)
  assert.equal(externalNull(new Map(), '9'), NULL_FLOOR, 'and nothing to measure floors too')
}

// --- 2. a real segment is found, and a clean chromosome yields nothing --------------------------
{
  const clean = scanChromosome(chrom(40_000), 0.005)
  assert.deepEqual(clean, [], 'a chromosome at the null rate carries no segment')

  const hit = scanChromosome(chrom(40_000, 5_000, 14_999), 0.005)
  assert.equal(hit.length, 1, 'one planted event is one segment')
  assert.ok(hit[0].score > SEGMENT_LRT)
  assert.ok(hit[0].rate > 0.08, 'and it reports the local rate, not the genome rate')
  assert.equal(hit[0].nullRate, 0.005, 'with the null it was scored against')
}

// --- 3. localisation: the best window, not the union of everything overlapping it ---------------
//
// This is a fixed defect, not a hypothetical. Merging every overlapping hit reported a 12 Mb loss
// as spanning 231 Mb, because a whole-chromosome window carrying a diluted version of the same
// event overlaps the tight one. The scan must describe the event, not the chromosome.
{
  const n = 40_000
  const hit = scanChromosome(chrom(n, 5_000, 14_999), 0.005)
  const seg = hit[0]
  const whole = (n - 1) * 5_000
  assert.ok(seg.spanBp < whole / 2,
    `a 10,000-marker event must not be reported as most of the chromosome: ${seg.spanBp} of ${whole}`)
  // The planted block sits at markers 5,000 to 14,999, which is 25 Mb to 75 Mb here.
  assert.ok(seg.startBp >= 1_000_000 + 3_000 * 5_000, 'starts near the event, not at the telomere')
  assert.ok(seg.endBp <= 1_000_000 + 17_000 * 5_000, 'and ends near it')
}

// --- 4. the marker floor is a refusal, not a weak answer ----------------------------------------
{
  // Below the floor there is no scan at all, however strong the signal.
  const tiny = chrom(2_000, 0, 1_999, 0.5)
  assert.deepEqual(scanChromosome(tiny, 0.005), [],
    `a chromosome under ${MIN_SEGMENT_MARKERS} informative markers is not scanned`)

  // At 1,200 markers a real spliced event scored 152 to 197 against a null maximum of 139, which
  // is 1.09x and not a detection. At 2,400 the weakest scored 431, which is 3.1x. The floor sits
  // where separation becomes real.
  assert.equal(MIN_SEGMENT_MARKERS, 2_400)
  assert.ok(SEGMENT_LRT > 139, 'above the measured null maximum')
  assert.ok(SEGMENT_LRT < 431, 'and below the weakest real detection at the floor')
}

// --- 5. a window CLEANER than the genome is not a finding ---------------------------------------
//
// The statistic is two-sided by construction and a fraction of the genome being unusually clean
// says nothing about the parent being missing there.
{
  const quiet = scanChromosome(chrom(40_000, 5_000, 14_999, 0.0, 0.11), 0.11)
  assert.deepEqual(quiet, [], 'a below-null window must never be reported as absence')
}

// --- 6. more absence scores higher, monotonically ------------------------------------------------
{
  const score = (rate: number): number =>
    scanChromosome(chrom(40_000, 5_000, 14_999, rate), 0.005, MIN_SEGMENT_MARKERS, 0)
      .reduce((a, s) => Math.max(a, s.score), 0)
  const at = [0.03, 0.06, 0.11, 0.2].map(score)
  for (let i = 1; i < at.length; i += 1) {
    assert.ok(at[i] > at[i - 1], `a stronger event must score higher: ${at}`)
  }
}

console.log('segments.check.ts: all assertions passed')

// --- 7. copy number is a different indicator from parental absence -----------------------------
//
// `scanChromosome` asks where a PARENT's alleles are missing, which a region can do with its DNA
// entirely present. `scanCopyNumber` asks where the DNA itself is gone, read from the array no
// longer calling. Both events are real and they are not the same event.
//
// Measured across 46 arrays: without the intensity requirement the copy scan returns 31 regions
// including a 44.9 Mb one on a bulk diploid adult, who cannot have it. With it, 17 survive and
// none fall on any of the six bulk diploid arrays. Rejected regions carry -0.72 to +0.38 log2;
// kept ones carry -0.97 to -1.94.
{
  const cn = (n: number, from: number, to: number, noCall: number, lrr: number, bg = 0.10) =>
    Array.from({ length: n }, (_, i) => {
      const inSeg = i >= from && i <= to
      const s = (i * 7919) % 1000
      return {
        chrom: '4',
        pos: 1_000_000 + i * 5_000,
        called: !(s < (inSeg ? noCall : bg) * 1000),
        log2R: inSeg ? lrr : 0,
      }
    })

  // A real deletion: the array stops calling and the intensity agrees.
  const lost = scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, -1.9), 0.10, 0)
  assert.equal(lost.length, 1, `a real deletion is one region: ${JSON.stringify(lost)}`)
  assert.equal(lost[0].kind, 'copy-loss')
  assert.ok(lost[0].spanBp > 20e6 && lost[0].spanBp < 120e6, 'and it is localised')

  // The same call-rate collapse with no intensity behind it is an array failing, not a deletion.
  // This is the case that produced a 44.9 Mb false loss on a bulk diploid adult.
  assert.deepEqual(scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, -0.3), 0.10, 0), [],
    'a region the array merely failed on is not a copy-number change')

  // Extra copies are the same machinery with the shift reversed. IMPLEMENTED AND UNVALIDATED: no
  // segmental gain occurs in the 46 arrays, so this has never fired on a true positive.
  const gained = scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, 1.9), 0.10, 0)
  assert.equal(gained.length, 1)
  assert.equal(gained[0].kind, 'copy-gain')

  // A clean chromosome yields nothing in either direction.
  assert.deepEqual(scanCopyNumber(cn(40_000, -1, -1, 0, 0), 0.10, 0), [])

  // The two channels are independent: a region whose DNA is present but whose PARENT is absent
  // is found by the allelic scan and not by this one.
  const parental = scanChromosome(
    Array.from({ length: 40_000 }, (_, i) => ({
      chrom: '4', pos: 1_000_000 + i * 5_000,
      absent: ((i * 7919) % 1000) < (i >= 5_000 && i < 15_000 ? 110 : 5),
    })), 0.005)
  assert.equal(parental.length, 1)
  assert.equal(parental[0].kind, 'parental-absence')
  assert.equal(SEGMENT_LRR_SHIFT, 1.0)
}

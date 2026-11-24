// Self-check for the segment scan. Run: node src/segments.check.ts
//
// Every figure quoted here was measured on GSE148488, on real genomes rather than simulated ones:
// the null on five paternal pronuclei of the sperm donor, and the positive class by splicing a
// maternal pronucleus of the same series into one of them so the segment carries real
// amplification artefact.
//
// CORRECTION, 3.7.2. That null set was described here as "present on every autosome" and it is
// not. GSM4774681 (pMII-2) carries a real chr1 event: an external review measured absence 0.190
// against 0.0105 genome-wide with log2R -1.79 and a multiscale LRT of 5,411 against a
// leave-one-out threshold of 134. Checked against the subset this repo ships, scoring pMII-2
// against a reference built from the other four products: chr1 is its worst chromosome at 9.78%
// absence against 1.04% genome-wide, 9.4x, which corroborates the direction on 460 informative
// markers. The same check flags pMII-6 chr20 at 12x, so the number of affected arrays is not
// established either.
//
// What follows from that, and it is not cosmetic: a real event sat inside the set the null was
// calibrated on, counted as noise. Every false-positive rate derived from that set is therefore
// OPTIMISTIC by an unknown margin. The thresholds below are not changed here, because changing
// them without recalibrating on a verified-clean set would be trading a known bias for an
// unknown one. They are flagged instead, and the recalibration is tracked as an open item.
import assert from 'node:assert/strict'
import { readFileSync, readdirSync } from 'node:fs'
import {
  MIN_SEGMENT_MARKERS, MIN_SEGMENT_SPAN_BP, NULL_FLOOR, SEGMENT_LRT, SEGMENT_LRR_SHIFT,
  SEGMENT_GAIN_SHIFT,
  RECIPROCAL_OVERLAP,
  canDetectLoss,
  externalNull,
  scanChromosome, scanCopyNumber, scanCopyGain, highLrrProfile, segmentCoords, MIN_GAIN_SPAN_BP,
  type MarkerAbsence,
} from './segments.ts'
import { OTHER_MIN_SPAN_BP } from './otherParentLoss.ts'
import { COPY_SHIFT_GAIN_FLOOR } from './parentage.ts'

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

// --- A SCAN THAT CANNOT SEE A REAL EVENT MUST NOT PASS FOR A CLEAN CHROMOSOME ------------------
//
// The null is the median absence rate over the sample's other autosomes, which holds only while
// most of them are unaffected. Past that the bar rises toward the rate a true event produces and
// no window can clear the threshold, which looks exactly like a chromosome carrying nothing.
// Measured on GSE148488: a true single-parent loss reads 0.1153 to 0.1244 against the parent whose
// copy went, and the three cells whose hand-called terminal deletion this scan misses read chr16
// maternal nulls of 0.0140, 0.0280 and 0.0282.
{
  const EVENT = 0.1244
  assert.ok(canDetectLoss(NULL_FLOOR, EVENT), 'a clean cell can see a real loss')
  assert.ok(!canDetectLoss(0.0280, EVENT), 'null 0.0280 is past the blind boundary')
  assert.ok(!canDetectLoss(0.0282, EVENT), 'null 0.0282 is past the blind boundary')
  // The third cell reads a lower null and is NOT explained this way: its event signal is weak.
  assert.ok(canDetectLoss(0.0140, EVENT), 'null 0.0140 is still inside the detectable range')
  assert.ok(canDetectLoss(0.0250, EVENT) && !canDetectLoss(0.0256, EVENT),
    'the boundary sits at a null of about 0.0253')
  // A larger window buys range back, which is the only lever that moves no threshold.
  assert.ok(canDetectLoss(0.0280, EVENT, 4_800), 'twice the window clears the same null')
  assert.ok(!canDetectLoss(0, EVENT), 'a zero null is a refusal, not infinite power')
  assert.ok(!canDetectLoss(0.20, EVENT), 'a null above the event rate can never detect it')
}

// --- TWO SEPARATE EVENTS ON ONE CHROMOSOME ARE NOT ONE CONTESTED REGION ------------------------
//
// The first version of the contested test was a bare intersection, which at the near-chromosome
// window lengths this corpus produces is satisfied by any two events on one chromosome. One arm
// losing the father's copy and the other losing the mother's is two findings; collapsing them
// into one region with no parent loses both.
{
  assert.equal(RECIPROCAL_OVERLAP, 0.5)
  const r = (a: [number, number], b: [number, number]) => {
    const inter = Math.min(a[1], b[1]) - Math.max(a[0], b[0])
    return inter <= 0 ? 0 : inter / (Math.max(a[1], b[1]) - Math.min(a[0], b[0]))
  }
  assert.ok(r([0, 50e6], [50e6, 100e6]) < RECIPROCAL_OVERLAP, 'abutting arms are two events')
  assert.ok(r([0, 80.9e6], [0, 76.6e6]) >= RECIPROCAL_OVERLAP, 'near-identical windows are one')
  // The two measured pairs that fall below on window edges, which is what the rule is for.
  assert.ok(r([40.7e6, 100.8e6], [20.0e6, 72.6e6]) < RECIPROCAL_OVERLAP)
  assert.ok(r([0.2e6, 135.0e6], [4.6e6, 35.4e6]) < RECIPROCAL_OVERLAP)
  assert.ok(r([0.2e6, 135.0e6], [43.0e6, 137.9e6]) >= RECIPROCAL_OVERLAP)
}

// --- EVERY SURFACE THAT LABELS A PARENTAL ABSENCE MUST READ WHOSE IT WAS -----------------------
//
// `absentParent` is which of the two absence runs produced the segment, so it is the parent the
// label has to name. Three surfaces print that label and all three hardcoded or assumed the
// loaded role, which was correct only while one run existed. This asserts the seam rather than
// the three instances: a fourth surface added later fails here instead of shipping a label that
// names the wrong parent on every segment the second run finds.
{
  const dir = new URL('./', import.meta.url)
  const missing = readdirSync(dir)
    .filter((f) => /\.(ts|tsx)$/.test(f) && !f.endsWith('.check.ts'))
    .filter((f) => {
      const src = readFileSync(new URL(f, dir), 'utf8')
      // A segment label is built from a template ("${...} alleles absent"); a genome-wide
      // measurement named "Paternal alleles absent" is not one.
      return /\$\{[^}]*\} alleles absent/.test(src) && !src.includes('absentParent')
    })
  assert.deepEqual(missing, [],
    `these label a parental-absence segment without reading whose alleles are missing: ${missing}`)
  // The label every surface prints is pageSays.segmentKindText, which reads it.
  assert.ok(/absentParent \?\? role\} alleles absent/.test(readFileSync(new URL('pageSays.ts', dir), 'utf8')))
}

console.log('segments.check.ts: all assertions passed, including that every surface labelling a '
  + 'parental absence reads which run found it')

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

  // Extra copies are the same machinery with the shift reversed. Validated on a constructed
  // positive class, since no segmental gain occurs in the 46 arrays: a block of a
  // whole-chromosome-gain array spliced into a clean array of the same series scores 3,162 at
  // the marker floor and 25,007 at 19,200 markers, copy-gain at all four sizes, while the same
  // splice from a EUPLOID array returns nothing at any of them. See scanCopyNumber's docstring.
  const gained = scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, 1.9), 0.10, 0)
  assert.equal(gained.length, 1)
  assert.equal(gained[0].kind, 'copy-gain')

  // The direction is read from the SIGN, so the two must never collapse onto each other.
  assert.equal(scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, -1.9), 0.10, 0)[0].kind, 'copy-loss')
  // And a shift that does not clear the bar is neither, in either direction.
  // THE ARMS ARE NOT SYMMETRIC, and the previous version of this check asserted that they were.
  // A hemizygous loss is log2(1/2) = 1.0 and a third copy is log2(3/2) = 0.585, so one gate for
  // both refused real gains: four of the hand-called regions carry 0.815 to 0.998, one of them
  // 145 Mb. A deviation of +0.9 is a gain and must be called.
  const nearGain = scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, 0.9), 0.10, 0)
  assert.equal(nearGain.length, 1, '+0.9 is over the gain gate and under the loss gate')
  assert.equal(nearGain[0].kind, 'copy-gain')
  // The same magnitude on the loss arm is still refused, because that gate is untouched. These two
  // lines together are the whole of the change: one direction moved and the other did not.
  assert.deepEqual(scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, -0.9), 0.10, 0), [])
  // And the gain arm is a gate, not an opening. Below log2(3/2) it still refuses.
  assert.deepEqual(scanCopyNumber(cn(40_000, 5_000, 14_999, 0.85, 0.4), 0.10, 0), [],
    '+0.4 is under the gain gate')

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
  // A window too short to be a segmental event is refused on BOTH arms, whatever its intensity.
  // At this platform's mean density the marker floor spans about 8.7 Mb, so a window holding it
  // inside 1 Mb is a probe cluster. cn() lays markers 5 kb apart, so 2,400 of them span 12 Mb;
  // the same count packed into 150 bp steps spans 0.36 Mb and must be refused.
  // 20 bp apart, so all 40,000 markers span 0.8 Mb and NO window the scan can build reaches the
  // floor. The identical event at the 5 kb spacing cn() uses is called above, so what this pins is
  // the span requirement and not that the event is unfindable.
  const dense = Array.from({ length: 40_000 }, (_, i) => ({
    chrom: '9', pos: 1_000_000 + i * 20,
    called: !(i >= 5_000 && i <= 14_999 && ((i * 7919) % 1000) < 850),
    log2R: i >= 5_000 && i <= 14_999 ? 1.9 : 0,
  }))
  assert.ok((dense[dense.length - 1].pos - dense[0].pos) < MIN_SEGMENT_SPAN_BP,
    'the dense fixture must be shorter than the floor for this to test anything')
  assert.deepEqual(scanCopyNumber(dense, 0.10, 0), [],
    'a window inside a probe cluster is not a segmental call at any intensity')

  // THE INVARIANT, over every fixture above: nothing this scan reports claims less than the floor.
  // Stated as a property because the span that matters is the REPORTED one, and a window of
  // several megabases can refine down to a tenth of one.
  for (const fixture of [
    cn(40_000, 5_000, 14_999, 0.85, -1.9), cn(40_000, 5_000, 14_999, 0.85, 1.9),
    cn(40_000, 5_000, 14_999, 0.85, 0.9), cn(40_000, 5_000, 14_999, 0.85, 0.4),
    cn(40_000, -1, -1, 0, 0), dense,
  ]) {
    for (const sg of scanCopyNumber(fixture, 0.10, 0)) {
      assert.ok(segmentCoords(sg).spanBp >= MIN_SEGMENT_SPAN_BP,
        `reported a ${segmentCoords(sg).spanBp} bp segment, under the ${MIN_SEGMENT_SPAN_BP} floor`)
    }
  }

  assert.equal(SEGMENT_LRR_SHIFT, 1.0)
  // The gain gate is DERIVED from the loss gate, so the two cannot drift apart by editing one.
  assert.ok(Math.abs(SEGMENT_GAIN_SHIFT - 0.5849625007211562) < 1e-12,
    'the gain gate is the loss gate times log2(3/2)')
  assert.ok(SEGMENT_GAIN_SHIFT < SEGMENT_LRR_SHIFT,
    'a third copy is a smaller deviation than a lost copy, so its gate is lower')
}

// --- THE INTENSITY-NOMINATED GAIN ARM: ITS FLOOR AND ITS OFF-POLE TEST ------------------------
//
// Two changes of 18 November, both in audit/truth/gains-called-2026-11-18.txt:
//   - the magnitude bar on this arm is the whole-chromosome gain floor, not the full 0.585 a third
//     copy carries in theory. Real gain pieces of this corpus read 0.23 to 0.50, so the old bar sat
//     above the entire observed population and stopped 9 pieces whose region test already cleared.
//   - a raised stretch must carry markers OFF the homozygous poles, which one homologue at twice
//     the intensity does not, and which an isodisomic third copy does not either: hence the bound
//     at ONE_EXTRA_COPY_MAX, above which the test applies and below which it does not.
{
  const N = 12_000
  const FROM = 4_000
  const TO = 7_000
  /** One chromosome: `delta` over the window, `off` whether the window's BAF sits off the poles. */
  const markers = (delta: number, off: boolean, chromName = '2') =>
    Array.from({ length: N }, (_, i) => {
      const inSeg = i >= FROM && i <= TO
      const jitter = ((i * 7919) % 11 - 5) / 200
      return {
        chrom: chromName,
        pos: 1_000_000 + i * 5_000,
        log2R: (inSeg ? delta : 0) + jitter,
        // Off the poles only inside the window when asked; the rest of the genome always carries
        // its own off-pole share, so the reference below is a real one.
        baf: inSeg && !off ? (i % 2 ? 0.02 : 0.98) : (i % 3 === 0 ? 0.5 : (i % 2 ? 0.02 : 0.98)),
      }
    })
  /** The reference share the caller passes: a third of the genome's markers sit off the poles. */
  const REFERENCE_OFF_POLE = 1 / 3
  const run = (delta: number, off: boolean, reference = REFERENCE_OFF_POLE) => {
    const ms = markers(delta, off)
    const byChrom = new Map<string, typeof ms>([['2', ms]])
    // Two quiet chromosomes so the nomination null is not this chromosome's own rate.
    for (const c of ['3', '4']) byChrom.set(c, markers(0, true, c))
    const profile = highLrrProfile(byChrom)
    return scanCopyGain(ms, externalNull(profile.byChrom, '2'), 0, profile, [], reference)
  }

  // THE MAGNITUDE BAR IS UNCHANGED, and this pins what that costs as well as what it buys. A piece
  // at 0.35 is inside the population of real third copies measured on this corpus (0.23 to 0.50)
  // and is still refused; lowering the bar to the whole-chromosome gain floor was run once and
  // printed 21 gains on 12 of the 72 arrays stated free of variation, against 9 on 2 today, in the
  // dense-probe windows this scan has always had to hold off.
  assert.equal(run(0.35, true).length, 0,
    'a gain piece under the full magnitude of a third copy is not printed')
  assert.equal(run(SEGMENT_GAIN_SHIFT + 0.05, true).length, 1, 'and one over it is')
  assert.ok(COPY_SHIFT_GAIN_FLOOR < SEGMENT_GAIN_SHIFT,
    'the whole-chromosome arm asks less, because a whole chromosome is not a window that probe '
    + 'density can manufacture')

  // The off-pole test, above the one-extra-copy bound.
  assert.equal(run(1.1, true).length, 1, 'a doubled stretch that still carries two alleles stands')
  assert.equal(run(1.1, false).length, 0,
    'and one at the poles is withheld: one homologue at twice the intensity reads the same step')
  // Below the bound it does not apply: an isodisomic third copy is at the poles and is real.
  assert.equal(run(SEGMENT_GAIN_SHIFT + 0.05, false).length, 1,
    'one extra copy of the same homologue carries one allele and must still be reported')
  // No reference to read against: the veto does not act, because it is not the evidence.
  assert.equal(run(1.1, false, NaN).length, 1,
    'an unreadable off-pole reference withholds nothing')
}

// --- A GAIN PIECE SPANS WHAT A GAIN THAT WAS EVER RIGHT SPANS ---------------------------------
//
// Measured over 636 records of 391 open arrays carrying 43 printed gain pieces: of the 6 under
// 3 Mb none is right by the paper and every one sits at 2.9x to 9.4x the genome's probe density,
// while no piece that IS right is under 10 Mb and all of them sit at 0.7x to 1.7x.
{
  assert.equal(MIN_GAIN_SPAN_BP, OTHER_MIN_SPAN_BP,
    'the gain floor is the figure the likelihood loss scans already use, not a new one; they are '
    + 'written out separately only because importing one from the other is a module cycle')
  assert.ok(MIN_GAIN_SPAN_BP > MIN_SEGMENT_SPAN_BP,
    'and it is above the span floor the LOSS arm keeps, which this does not touch')

  const N = 12_000
  const markers = (from: number, to: number, delta: number, chromName = '2') =>
    Array.from({ length: N }, (_, i) => {
      const inSeg = i >= from && i <= to
      const jitter = ((i * 7919) % 11 - 5) / 200
      return { chrom: chromName, pos: 1_000_000 + i * 5_000, log2R: (inSeg ? delta : 0) + jitter,
        baf: i % 3 === 0 ? 0.5 : (i % 2 ? 0.02 : 0.98) }
    })
  const run = (from: number, to: number) => {
    const ms = markers(from, to, 0.9)
    const byChrom = new Map<string, typeof ms>([['2', ms]])
    for (const c of ['3', '4']) byChrom.set(c, markers(0, -1, 0, c))
    const profile = highLrrProfile(byChrom)
    return scanCopyGain(ms, externalNull(profile.byChrom, '2'), 0, profile, [], 1 / 3)
  }
  // 2,401 markers over 12 Mb: it clears the marker floor and the span floor.
  assert.equal(run(4_000, 7_000).length, 1, 'a 15 Mb gain piece is printed')
  // The same marker count packed into 2.4 Mb, which is what a probe cluster looks like.
  const dense = (() => {
    const ms = Array.from({ length: N }, (_, i) => {
      const inSeg = i >= 4_000 && i <= 7_000
      // 800 bp apart inside the window, 5 kb outside: a 9x cluster, as the HLA region reads.
      const pos = i <= 4_000 ? 1_000_000 + i * 5_000
        : i <= 7_000 ? 1_000_000 + 4_000 * 5_000 + (i - 4_000) * 800
          : 1_000_000 + 4_000 * 5_000 + 3_000 * 800 + (i - 7_000) * 5_000
      return { chrom: '2', pos, log2R: (inSeg ? 0.9 : 0) + ((i * 7919) % 11 - 5) / 200,
        baf: i % 3 === 0 ? 0.5 : (i % 2 ? 0.02 : 0.98) }
    })
    const byChrom = new Map<string, typeof ms>([['2', ms]])
    for (const c of ['3', '4']) byChrom.set(c, markers(0, -1, 0, c))
    const profile = highLrrProfile(byChrom)
    return scanCopyGain(ms, externalNull(profile.byChrom, '2'), 0, profile, [], 1 / 3)
  })()
  assert.equal(dense.length, 0,
    'and the same markers packed into 2.4 Mb are not: every gain piece under 3 Mb on the open '
    + 'corpus is a probe cluster and none of them is right')
}

console.log('segments.check (the gain span floor): ok')
console.log('segments.check (intensity gain arm): ok')

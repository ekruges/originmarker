// Self-check for the phase-2 run-length discriminator. Run: node src/runlength.check.ts
//
// Section 1 is the acceptance test: all 12 cases of golden_test_vectors.csv against the expected
// summary in golden_test_expected_outputs.csv. Those files were generated independently of this
// code, so agreement is a real check rather than a restatement.
//
// Sections 2 to 5 pin the properties the fixture cannot exercise: numerical stability where the
// Erdos-Renyi form cancels, r_min falling out of q rather than being hardcoded, the max-not-sum
// composition of q, and the below-resolution verdict that stops a small event reading as absent.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import {
  KOTHIYAL_FLOOR, MAX_CLUSTERING, MIN_EXPECTED_RUNS, analyseRuns, measureClustering, parseLocus,
  qHat, rMin, runLengthP, scoreMarker, type Clustering,
} from './runlength.ts'
import type { AB, Marker } from './informativity.ts'

const load = (name: string) => {
  const text = readFileSync(new URL(`../../tests/fixtures/${name}`, import.meta.url), 'utf8')
  const rows = text.split('\n').filter((l) => l.trim() !== '')
  const head = rows[0].split(',')
  return rows.slice(1).map((l) => {
    // No quoted commas in the observation columns; the note column is last and unused here.
    const cells = l.split(',')
    return Object.fromEntries(head.map((h, i) => [h, cells[i] ?? ''])) as Record<string, string>
  })
}

// --- 1. the 12 golden cases -----------------------------------------------------------
const vectors = load('golden_test_vectors.csv')
const expected = load('golden_test_expected_outputs.csv')
assert.equal(expected.length, 12, 'twelve acceptance cases')

const CHROM = '11'
const VARIANT = 1_090_000        // between-marker, matching the fixture's "variant at marker 10"

/** A sample that demonstrated independence: enough expected runs to measure, and no excess. */
const INDEPENDENT: Clustering = {
  n: 500_000, q: 0.005, observed: 12, expected: 12, ratio: 1, independent: true,
}

let cases = 0
for (const exp of expected) {
  const rows = vectors.filter((v) => v.case === exp.case)
  assert.ok(rows.length === 20, `${exp.case}: 20 markers`)

  const markers: Marker[] = rows.map((r) => ({ rsid: r.marker, chrom: CHROM, pos: Number(r.pos) }))
  const father = new Map<string, AB>(rows.map((r) => [r.marker, r.father_gt as AB]))
  const mother = new Map<string, AB>(rows.map((r) => [r.marker, r.mother_gt as AB]))
  const embryo = new Map<string, AB>(rows.map((r) => [r.marker, r.embryo_gt as AB]))

  // The fixture's statistics are over the whole panel, so read the widest window.
  //
  // Independence is asserted here rather than measured, because what this fixture pins is the
  // ARITHMETIC of the run-length tail, shared with the Python. Whether a real sample has earned
  // that tail is a separate property and section 6 pins it.
  const res = analyseRuns(markers, father, mother, embryo, CHROM, VARIANT,
    { clustering: INDEPENDENT })
    .find((w) => w.window === 'whole_chromosome')!

  const where = exp.case
  assert.equal(res.nL3, Number(exp.n_L3), `n_L3: ${where}`)
  assert.equal(res.zSum, Number(exp.z_sum), `z_sum: ${where}`)
  assert.equal(res.longestRun, Number(exp.longest_run), `longest_run: ${where}`)
  assert.equal(res.longestRunMaternal, Number(exp.longest_run_maternal),
    `longest_run_maternal: ${where}`)
  assert.equal(res.rMin, Number(exp.r_min), `r_min: ${where}`)

  // run_p to the fixture's two significant figures, EXCEPT where the fixture recorded the
  // catastrophically-cancelled value: see section 2.
  const expP = Number(exp.run_p)
  if (expP > 0) {
    assert.ok(Math.abs(res.runP - expP) <= 0.02 * expP,
      `run_p: ${where} mine ${res.runP.toExponential(3)} vs fixture ${exp.run_p}`)
  }
  cases++
}
assert.equal(cases, 12)

// The pair that carries the design lesson: C3 and C5 have IDENTICAL run statistics, because this
// statistic cannot tell copy-loss from copy-neutral LOH. Only LRR separates them, in phase 3.
{
  const stats = (name: string) => {
    const rows = vectors.filter((v) => v.case === name)
    const markers: Marker[] = rows.map((r) => ({ rsid: r.marker, chrom: CHROM, pos: Number(r.pos) }))
    const mk = (col: string) => new Map<string, AB>(rows.map((r) => [r.marker, r[col] as AB]))
    return analyseRuns(markers, mk('father_gt'), mk('mother_gt'), mk('embryo_gt'), CHROM, VARIANT)
      .find((w) => w.window === 'whole_chromosome')!
  }
  const c3 = stats('C3_H3c_segmental_loss'), c5 = stats('C5_H3e_copyneutral_LOH')
  assert.equal(c3.longestRun, c5.longestRun)
  assert.equal(c3.runP, c5.runP)
  assert.equal(c3.verdict, c5.verdict)
  assert.match(c3.note, /does NOT distinguish copy-loss from copy-neutral/)
}

// --- 2. numerical stability where the direct form fails --------------------------------
// 1 - (1 - q^r)^m cancels catastrophically: at q = 0.0063 and r = 10 the true value is 9.8e-23,
// but `1 - x` rounds to 1 in float64 and the subtraction returns 0. The fixture's run_p of
// 0.00e+00 for C4 is that artefact, not a bug in its answer key.
{
  const naive = (r: number, n: number, q: number) => 1 - (1 - q ** r) ** (n - r + 1)
  assert.equal(naive(10, 10, KOTHIYAL_FLOOR), 0, 'the direct form really does return zero')
  const stable = runLengthP(10, 10, KOTHIYAL_FLOOR)
  assert.ok(stable > 0, 'the stable form does not')
  assert.ok(Math.abs(stable - 9.849e-23) < 1e-25, `got ${stable.toExponential(4)}`)
  // And where the direct form is fine the two must agree, or the rewrite changed the maths.
  for (const [r, n, q] of [[2, 100, 0.05], [4, 10, 0.0063], [3, 50, 0.10]] as const) {
    assert.ok(Math.abs(runLengthP(r, n, q) - naive(r, n, q)) < 1e-12, `r=${r} n=${n} q=${q}`)
  }
}
assert.equal(runLengthP(0, 10, 0.1), 1)
assert.equal(runLengthP(11, 10, 0.1), 0, 'a run longer than the marker count is impossible')
assert.equal(runLengthP(3, 10, 0), 0)

// --- 3. r_min is computed, never hardcoded --------------------------------------------
// The spec states these, and reproducing them independently is what makes the threshold
// auditable rather than magic.
assert.equal(rMin(100, 0.01, 0.05, 3), 2, 'spec 5.3: q=0.01, n=100, 3 windows')
assert.equal(rMin(100, 0.05, 0.05, 3), 3, 'spec 5.3: q=0.05')
assert.equal(rMin(100, 0.15, 0.05, 3), 5, 'spec 5.3: q=0.15')
assert.equal(rMin(10, KOTHIYAL_FLOOR, 0.05, 3), 2, 'spec 5.3: Kothiyal floor, n=10')
// A noisier sample needs a longer run, monotonically. That relationship is the whole reason the
// threshold cannot be a constant.
let prev = 0
for (const q of [0.0063, 0.01, 0.05, 0.10, 0.15, 0.25]) {
  const r = rMin(100, q, 0.05, 3)!
  assert.ok(r >= prev, `r_min must not fall as q rises: q=${q} gave ${r} after ${prev}`)
  prev = r
}
// And each extra marker in a run buys roughly twenty-fold, which is why runs are decisive.
const ratio = runLengthP(4, 100, 0.05) / runLengthP(5, 100, 0.05)
assert.ok(ratio > 10 && ratio < 40, `per-marker gain was ${ratio.toFixed(1)}x`)

// --- 4. q composes as max, not sum ----------------------------------------------------
// Kothiyal is a genome-wide trio floor that already contains error classes the empirical
// estimator also captures, so summing double-counts. It is a LOWER BOUND on q.
assert.deepEqual(qHat(null), { q: KOTHIYAL_FLOOR, source: 'kothiyal_floor' })
assert.deepEqual(qHat(0.001), { q: KOTHIYAL_FLOOR, source: 'kothiyal_floor' })
assert.deepEqual(qHat(0.05), { q: 0.05, source: 'empirical' })
assert.equal(qHat(0.0063).source, 'kothiyal_floor', 'at the floor exactly, the floor is used')
assert.equal(qHat(NaN).q, KOTHIYAL_FLOOR, 'a failed fit is not a q of NaN')
// Not the sum: 0.05 + 0.0063 would be 0.0563 and would inflate every r_min.
assert.equal(qHat(0.05).q, 0.05)

// --- 5. scoring, and the two ways a marker can say nothing ----------------------------
// s_pat exists only where the father is HOMOZYGOUS: a heterozygous father is compatible with
// either observed allele, so he cannot testify to presence or absence at all.
assert.deepEqual(scoreMarker('AA', 'BB', 'AB'), { sPat: 1, sMat: 1 })
assert.deepEqual(scoreMarker('AA', 'BB', 'BB'), { sPat: 0, sMat: 1 }, 'paternal allele absent')
assert.deepEqual(scoreMarker('AA', 'BB', 'AA'), { sPat: 1, sMat: 0 }, 'the mirror')
assert.deepEqual(scoreMarker('BB', 'AA', 'AB'), { sPat: 1, sMat: 1 }, 'symmetric under A/B')
assert.deepEqual(scoreMarker('BB', 'AA', 'AA'), { sPat: 0, sMat: 1 })
assert.equal(scoreMarker('AB', 'AA', 'AB').sPat, null, 'father het: s_pat undefined')
assert.equal(scoreMarker('AA', 'AB', 'BB').sMat, null, 'mother het: s_mat undefined')
assert.deepEqual(scoreMarker('AA', 'BB', 'NC'), { sPat: null, sMat: null }, 'excluded, not scored')
assert.deepEqual(scoreMarker('NC', 'BB', 'AB'), { sPat: null, sMat: null }, 'father NC')
assert.equal(scoreMarker('AA', null, 'BB').sMat, null, 'no mother: no maternal side at all')
// The missing-mother mode still scores the paternal side, which is why it is worth running.
assert.equal(scoreMarker('AA', null, 'BB').sPat, 0, 'paternal absence provable without her')

// --- 6. a window with no L3-capable markers is not a negative result ------------------
{
  const markers: Marker[] = Array.from({ length: 20 }, (_, i) =>
    ({ rsid: `h${i}`, chrom: CHROM, pos: 1_000_000 + i * 10_000 }))
  const allHet = new Map<string, AB>(markers.map((m) => [m.rsid, 'AB' as AB]))
  const mo = new Map<string, AB>(markers.map((m) => [m.rsid, 'AA' as AB]))
  const res = analyseRuns(markers, allHet, mo, allHet, CHROM, VARIANT)
    .find((w) => w.window === 'whole_chromosome')!
  assert.equal(res.nL3, 0)
  assert.equal(res.verdict, 'undefined_father_heterozygous')
  assert.match(res.note, /not a negative result/)
  assert.equal(res.rMin, null, 'no threshold exists where the statistic is undefined')
}

// --- 7. below resolution is not the same as absent ------------------------------------
// The likeliest way to misuse the statistic. Here r_min x median L3 spacing is 40 kb, so a 20 kb
// event cannot produce a significant run and "absent" would be a false negative.
{
  const rows = vectors.filter((v) => v.case === 'C1_clean_H1')
  const markers: Marker[] = rows.map((r) => ({ rsid: r.marker, chrom: CHROM, pos: Number(r.pos) }))
  const mk = (col: string) => new Map<string, AB>(rows.map((r) => [r.marker, r[col] as AB]))
  const opts = { eventSizeOfInterestBp: 20_000 }
  const res = analyseRuns(markers, mk('father_gt'), mk('mother_gt'), mk('embryo_gt'),
    CHROM, VARIANT, opts).find((w) => w.window === 'whole_chromosome')!
  assert.equal(res.medianSpacingL3Bp, 20_000, 'L3 markers are every OTHER marker, so 20 kb apart')
  assert.equal(res.resolutionFloorBp, 40_000, 'r_min 2 x 20 kb')
  assert.equal(res.verdict, 'below_resolution')
  assert.match(res.note, /never as absent/)
  // An event above the floor gets the ordinary no-run verdict instead.
  const big = analyseRuns(markers, mk('father_gt'), mk('mother_gt'), mk('embryo_gt'),
    CHROM, VARIANT, { eventSizeOfInterestBp: 500_000 }).find((w) => w.window === 'whole_chromosome')!
  assert.equal(big.verdict, 'no_significant_run')
}

// --- 8. markers need not arrive sorted -----------------------------------------------
// Runs are physical, so an unsorted input would silently fragment them and understate every run.
{
  const rows = vectors.filter((v) => v.case === 'C3_H3c_segmental_loss')
  const mk = (col: string) => new Map<string, AB>(rows.map((r) => [r.marker, r[col] as AB]))
  const inOrder: Marker[] = rows.map((r) => ({ rsid: r.marker, chrom: CHROM, pos: Number(r.pos) }))
  const shuffled = [...inOrder].reverse()
  const run = (ms: Marker[]) => analyseRuns(ms, mk('father_gt'), mk('mother_gt'), mk('embryo_gt'),
    CHROM, VARIANT).find((w) => w.window === 'whole_chromosome')!.longestRun
  assert.equal(run(shuffled), run(inOrder), 'reversed input must give the same run')
  assert.equal(run(inOrder), 4)
  // Off-chromosome markers must not join a run either.
  const offChrom = inOrder.map((m, i) => (i % 2 ? { ...m, chrom: '7' } : m))
  assert.ok(run(offChrom) <= run(inOrder), 'chr7 markers cannot extend a chr11 run')
}

console.log(`runlength.check.ts: all assertions passed (${cases} golden cases)`)

// --- 8. presence needs the mother, absence does not -------------------------------------------
// The half of scoreMarker that `origin.score_paternal` exists to get right. At a paternal
// deletion the embryo is hemizygous for the maternal allele, so wherever a heterozygous mother
// transmitted the allele the father also carries, the call still contains his allele. Reading
// that as presence breaks the run at roughly half of all mother-heterozygous markers.
{
  // Absence is Mendelian and needs nothing from her.
  for (const mo of ['AA', 'AB', 'BB', 'NC', null] as const) {
    assert.equal(scoreMarker('AA', mo, 'BB').sPat, 0, `absence must hold with mother ${mo}`)
  }
  // Presence needs her homozygous for the other allele.
  assert.equal(scoreMarker('AA', 'BB', 'AB').sPat, 1, 'she could not have supplied the A')
  assert.equal(scoreMarker('AA', 'AB', 'AA').sPat, null, 'she could have supplied it: silent')
  assert.equal(scoreMarker('AA', 'AA', 'AA').sPat, null)
  assert.equal(scoreMarker('AA', null, 'AA').sPat, null, 'no mother: presence unobservable')
  assert.equal(scoreMarker('AA', 'NC', 'AA').sPat, null)
  // A heterozygous father testifies to nothing in either direction.
  assert.equal(scoreMarker('AB', 'BB', 'AA').sPat, null)

  // The run that motivated the rule: a paternal deletion under a heterozygous mother. Every
  // marker is absence or silence, and none may read as presence.
  const f = new Map<string, AB>()
  const m = new Map<string, AB>()
  const e = new Map<string, AB>()
  const markers = []
  for (let i = 0; i < 40; i += 1) {
    const id = `m${i}`
    f.set(id, 'AA')
    // She transmits A at every other marker, which is where the old rule saw "present".
    m.set(id, 'AB')
    e.set(id, i % 2 ? 'AA' : 'BB')
    markers.push({ rsid: id, chrom: '1', pos: 1_000_000 + i * 1000 })
  }
  const scores = markers.map((k) => scoreMarker(f.get(k.rsid)!, m.get(k.rsid)!, e.get(k.rsid)!))
  assert.ok(scores.every((s) => s.sPat !== 1), 'no marker may read as paternal presence here')
  assert.equal(scores.filter((s) => s.sPat === 0).length, 20, 'the BB half is provable absence')
}

// --- 9. r_min never drops to 1, and a locus is a human chromosome -----------------------------
// Contiguity is the entire discriminator, so a run of one has nothing to offer. Below n = 3 the
// uncorrected threshold came back as 1, which made a single absent marker a significant run.
for (const n of [1, 2, 3, 5, 50, 500]) {
  const r = rMin(n, 0.0063)
  assert.ok(r === null || r >= 2, `r_min(${n}) = ${r}`)
}
// A higher fitted q demands a longer run. The Kothiyal floor is a LOWER bound on the spurious
// rate, so leaving q there on amplified material shrinks r_min and over-calls.
assert.ok(rMin(200, 0.30)! > rMin(200, 0.05)!, 'more dropout, longer run required')
assert.ok(rMin(200, 0.05)! > rMin(200, KOTHIYAL_FLOOR)!)
assert.equal(qHat(0.002).source, 'kothiyal_floor', 'below the floor falls back to it')
assert.equal(qHat(0.30).source, 'empirical')

for (const [text, want] of [
  ['chr7:117559590', { chrom: '7', pos: 117559590 }],
  ['7:117559590', { chrom: '7', pos: 117559590 }],
  ['chr7 117,559,590', { chrom: '7', pos: 117559590 }],
  ['chrX:2,700,151', { chrom: 'X', pos: 2700151 }],
  ['22:1', { chrom: '22', pos: 1 }],
  // Not chromosomes. These used to parse and then match no marker at all.
  ['chr0:100', null], ['chr23:1', null], ['chr45:100', null], ['chr99:1', null],
  ['rs334', null], ['', null], ['chr7:', null], ['7', null],
] as const) {
  assert.deepEqual(parseLocus(text), want, `parseLocus(${JSON.stringify(text)})`)
}

// --- 6. the p-value is withheld unless the sample demonstrates independence -------------------
//
// The tail in `runLengthP` assumes a dropout at marker i says nothing about marker i+1. Measured
// on the public series that is false on every material class: WGA single cells run 6.1x to 10.3x
// more maximal runs of two than independence predicts, and even bulk gDNA runs 2.0x. So a run of
// artefact reaches a length the model calls impossible and the p-value reports significance for
// it. The run is still reported; the significance is not.
{
  // Alternating singletons: absences at a real rate, never two adjacent, so no runs at all.
  const spread = measureClustering([
    Array.from({ length: 20_000 }, (_, i) => i % 7 === 0)])
  assert.equal(spread.observed, 0, 'singletons make no runs of two')
  assert.ok(spread.expected >= MIN_EXPECTED_RUNS, 'and enough expected to measure against')
  assert.equal(spread.independent, true, 'so this sample has demonstrated independence')

  // The same rate, in pairs. Identical q, every absence contiguous.
  const paired = measureClustering([
    Array.from({ length: 20_000 }, (_, i) => i % 14 === 0 || i % 14 === 1)])
  assert.ok(paired.ratio > MAX_CLUSTERING, `pairs must read clustered, got ${paired.ratio}`)
  assert.equal(paired.independent, false)

  // Too few expected runs to settle anything is not a pass: absence of evidence is not evidence.
  const sparse = measureClustering([Array.from({ length: 300 }, (_, i) => i % 150 === 0)])
  assert.ok(sparse.expected < MIN_EXPECTED_RUNS)
  assert.equal(sparse.independent, false, 'an unmeasurable sample cannot demonstrate independence')

  // A run long enough to be significant, under each of the three conditions.
  const markers2: Marker[] = Array.from({ length: 60 }, (_, i) => ({
    rsid: `x${i}`, chrom: CHROM, pos: 1_000_000 + i * 10_000,
  }))
  const father2 = new Map<string, AB>(markers2.map((m) => [m.rsid, 'AA' as AB]))
  const embryo2 = new Map<string, AB>(markers2.map((m, i) =>
    [m.rsid, (i >= 20 && i < 32 ? 'BB' : 'AA') as AB]))
  const at = (c: Clustering | null) => analyseRuns(markers2, father2, undefined, embryo2,
    CHROM, 1_450_000, c ? { clustering: c } : {})
    .find((w) => w.window === 'whole_chromosome')!

  const shown = at(INDEPENDENT)
  assert.equal(shown.verdict, 'significant_run')
  assert.ok(shown.runP > 0 && shown.runP < 1e-6, 'a demonstrated sample keeps its p-value')

  for (const [label, c] of [
    ['clustered', { n: 500_000, q: 0.005, observed: 120, expected: 12, ratio: 10, independent: false }],
    ['not supplied', null],
  ] as const) {
    const r = at(c as Clustering | null)
    assert.equal(r.verdict, 'independence_not_demonstrated', label)
    assert.ok(Number.isNaN(r.runP), `${label}: the p-value must be withheld, not weakened`)
    assert.equal(r.longestRun, shown.longestRun, `${label}: the run itself is still reported`)
    assert.ok(!r.note.includes('p = '), `${label}: and no p-value survives in the prose`)
  }
}

// Self-check for the phase-0 ingestion layer. Run: node src/ingest.check.ts
//
// Fixtures are synthetic on purpose. The real files this module reads are patient genotypes
// from a consented protocol, and a test that depends on one would put that data in the repo.
//
// The externally-anchored test is section 6: the critical no-call counts are reproduced
// independently here and must agree with the values derived in the ingestion research. Two
// independent implementations of the same exact binomial agreeing is the only check available,
// since no published method exists for this scan.
import assert from 'node:assert/strict'
import {
  accumulate, accumulateBaf, accumulateBuild, assessRelatedness, binomUpperTail, buildVerdict, callSex, criticalNocalls, detectBuild, detectProduct, emptyBafSums, emptyBuildSums, estimateDropout, finishProfile, gates, headerMap, idPrefixOf, parseRow, profileText, scanNocallClusters, verifyCoding, type ChromStats, type ProbeRow, type SampleProfile,
} from './ingest.ts'
import type { AB } from './informativity.ts'

const HEADER = 'probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset'

// --- 1. columns are found by NAME, not by position ------------------------------------
// A reordered export is a plausible variation, and reading column 6 blindly would silently
// swap copy_number for genotype - two numeric columns that would not fail any type check.
const m = headerMap(HEADER)!
assert.ok(m)
assert.equal(m['genotype'], 6)
const reordered = headerMap('genotype\tchr\tposition\tprobeset_id\tbaf')!
assert.equal(reordered['genotype'], 0)
assert.equal(reordered['probeset_id'], 3)
// A file missing a required column is not partially usable.
assert.equal(headerMap('chr\tposition\tbaf'), null, 'no probeset_id or genotype: refuse')
// The four required columns alone are enough; log2R, baf and copy_number are optional, and
// their absence costs specific capabilities rather than the whole file.
assert.ok(headerMap('probeset_id\tchr\tposition\tgenotype'), 'required four suffice')
assert.ok(headerMap('PROBESET_ID\tCHR\tPOSITION\tGENOTYPE'), 'header match is case-insensitive')

// --- 2. the numeric genotype coding ---------------------------------------------------
const row = (s: string) => parseRow(s, m)!
assert.equal(row('AX-1\t1\t1000\t0.1\t0.01\t2.0\t0\t1').genotype, 'AA')
assert.equal(row('AX-1\t1\t1000\t0.1\t0.50\t2.0\t1\t1').genotype, 'AB')
assert.equal(row('AX-1\t1\t1000\t0.1\t0.99\t2.0\t2\t1').genotype, 'BB')
assert.equal(row('AX-1\t1\t1000\t0.1\t\t2.0\t-1\t1').genotype, 'NC')
// A blank field is missing data, not zero. The real files leave baf empty at ~4% of markers,
// and reading that as BAF 0.0 would place them on the homozygous band.
assert.equal(row('AX-1\t1\t1000\t0.1\t\t2.0\t1\t1').baf, null)
assert.equal(row('AX-1\t1\t1000\t\t0.5\t\t1\t1').log2R, null)
assert.equal(row('AX-1\t1\t1000\t0.0\t0.5\t2.0\t1\t1').log2R, 0, 'a real zero is not missing')
assert.equal(parseRow('\t1\t1000\t0.1\t0.5\t2\t1\t1', m), null, 'no probeset id: unusable row')
assert.equal(parseRow('AX-1\t1\t\t0.1\t0.5\t2\t1\t1', m), null, 'no position: unusable row')

// --- 3. sex, from the sample's own chrX/autosome ratio --------------------------------
// Dropout-robust because it scales numerator and denominator together. What must NOT happen is
// the middle band being resolved: that is where a chrX abnormality, a pooled file and
// contamination all land, and each is a finding rather than a rounding error.
const chroms = (autoHet: number, xHet: number): Map<string, ChromStats> => {
  const b = new Map<string, ChromStats>()
  for (let i = 1; i <= 22; i++) {
    b.set(String(i), { markers: 1000, called: 1000, het: Math.round(1000 * autoHet), nocall: 0 })
  }
  b.set('X', { markers: 1000, called: 1000, het: Math.round(1000 * xHet), nocall: 0 })
  return b
}
assert.equal(callSex(chroms(0.30, 0.006)).sex, 'male', 'PAR-only het is a male')
assert.equal(callSex(chroms(0.30, 0.24)).sex, 'female', 'ratio 0.8 is a female')
assert.equal(callSex(chroms(0.30, 0.21)).sex, 'female', 'ratio 0.7 is a female')
assert.equal(callSex(chroms(0.30, 0.14)).sex, 'ambiguous', 'ratio 0.46 is NOT a female')
assert.equal(callSex(chroms(0.30, 0.09)).sex, 'ambiguous', 'ratio 0.3 is neither')
// The ratio survives dropout: halve every het rate and the call is unchanged.
assert.equal(callSex(chroms(0.15, 0.003)).sex, 'male')
assert.equal(callSex(chroms(0.11, 0.088)).sex, 'female', 'heavy dropout, ratio 0.8, still female')
assert.equal(callSex(chroms(0.30, 0)).ratio, 0)
assert.equal(callSex(new Map()).sex, 'ambiguous', 'no data is not a sex call')

// --- 4. product, and the methylation trap ---------------------------------------------
assert.equal(detectProduct(825_656, 'AX-'), 'UK Biobank Axiom Array')
assert.equal(detectProduct(823_046, 'AX-'), 'UK Biobank Axiom Array', 'bi-allelic-only export')
assert.match(detectProduct(1_400_000, 'AX-'), /product unidentified/)
// A UCSC track named snpArrayIllumina850k is a METHYLATION array. File and track naming cannot
// be trusted, so cg/ch prefixes are rejected outright rather than treated as genotypes.
assert.match(detectProduct(850_000, 'cg'), /^REJECT/)
assert.match(detectProduct(850_000, 'ch'), /^REJECT/)
assert.equal(idPrefixOf('AX-13216142'), 'AX-')
assert.equal(idPrefixOf('rs429358'), 'rs')
assert.equal(idPrefixOf('cg00000029'), 'cg')

// --- 5. the coding check: detectable, unlike the nucleotide convention ----------------
// Which NUCLEOTIDE is called A is unrecoverable from a file with no nucleotides. Which NUMERAL
// is the B-homozygote is a different question, and BAF answers it.
const sums = (b0: number, b2: number) =>
  ({ hom0: b0 * 100, n0: 100, hom2: b2 * 100, n2: 100, het: 50, nHet: 100 })
assert.equal(verifyCoding(sums(0.04, 0.97)).verdict, 'ok')
assert.equal(verifyCoding(sums(0.97, 0.04)).verdict, 'inverted', 'code 2 carrying no B allele')
assert.equal(verifyCoding(sums(0.45, 0.55)).verdict, 'untestable', 'unseparated means: BAF unreliable')
assert.equal(verifyCoding(emptyBafSums()).verdict, 'untestable', 'no BAF column')
// An inverted file is a hard exclude, because the two channels contradict each other and one is
// wrong. (A file inverted UNIFORMLY relative to nothing is harmless - the informativity logic
// is invariant under a consistent A/B relabel - which is exactly why this checks BAF and not
// the genotypes alone.)
const invertedProfile = finishProfile('x', chroms(0.30, 0.24), sums(0.97, 0.04), 'AX-1')
assert.equal(gates(invertedProfile).find((g) => g.name === 'numeric genotype coding')!.verdict, 'exclude')

// --- 6. dropout, fitted from the heterozygote deficit --------------------------------
// At father-het x mother-hom markers Mendel gives exactly half AA and half AB regardless of
// allele frequency, so d = 1 - 2h needs no reference panel and no assumed population het rate.
const dropoutCase = (n: number, hetFraction: number) => {
  const fa = new Map<string, AB>(), mo = new Map<string, AB>(), em = new Map<string, AB>()
  for (let i = 0; i < n; i++) {
    fa.set(`r${i}`, 'AB'); mo.set(`r${i}`, 'AA')
    em.set(`r${i}`, i < Math.round(n * hetFraction) ? 'AB' : 'AA')
  }
  return estimateDropout(fa, mo, em)!
}
assert.equal(dropoutCase(10_000, 0.50).d, 0, 'no deficit is no dropout')
assert.ok(Math.abs(dropoutCase(10_000, 0.25).d - 0.50) < 1e-9, 'half the heterozygotes gone')
assert.ok(Math.abs(dropoutCase(10_000, 0.185).d - 0.63) < 1e-9)
// The precision is load-bearing: a 2 percentage point dropout mismatch between two samples
// inverts genome-wide IBS0, so the estimator must resolve better than that before any
// count-based relatedness statistic may be used at all.
assert.ok(dropoutCase(10_000, 0.185).se < 0.02, `se was ${dropoutCase(10_000, 0.185).se}`)
assert.ok(dropoutCase(1_000, 0.185).se > 0.02, 'and at 1k markers it does NOT resolve that')
// Only father-het x mother-hom markers count. Everything else carries no known expectation.
{
  const fa = new Map<string, AB>([['a', 'AB'], ['b', 'AA'], ['c', 'AB']])
  const mo = new Map<string, AB>([['a', 'AA'], ['b', 'AA'], ['c', 'AB']])
  const em = new Map<string, AB>([['a', 'AB'], ['b', 'AB'], ['c', 'AB']])
  assert.equal(estimateDropout(fa, mo, em)!.n, 1, 'only marker a qualifies')
}
// A no-call is missing data, not a homozygote, so it leaves the denominator rather than
// inflating the dropout estimate.
{
  const fa = new Map<string, AB>([['a', 'AB'], ['b', 'AB']])
  const mo = new Map<string, AB>([['a', 'AA'], ['b', 'AA']])
  assert.equal(estimateDropout(fa, mo, new Map([['a', 'AB'], ['b', 'NC']]))!.n, 1)
}
assert.equal(estimateDropout(new Map([['a', 'AA']]), new Map([['a', 'AA']]), new Map([['a', 'AA']])), null)

// --- 7. the local no-call scan --------------------------------------------------------
// Independent reproduction of the critical counts derived in the ingestion research. Two
// separate implementations of the same exact binomial agreeing is the available check, since no
// published method exists for this statistic.
const ALPHA = 0.05 / 825_656
for (const [win, expected] of [[20, 12], [50, 20], [100, 30], [200, 47], [500, 90], [1000, 155]] as const) {
  assert.equal(criticalNocalls(win, 0.10, ALPHA), expected, `critical count at window ${win}`)
  assert.ok(binomUpperTail(expected, win, 0.10) < ALPHA, `tail below alpha at window ${win}`)
  assert.ok(binomUpperTail(expected - 1, win, 0.10) > ALPHA, `and ${expected} is the SMALLEST such`)
}
// The tail must not underflow to a wrong answer at large n, which is why it is computed in log
// space and summed from the top.
assert.ok(binomUpperTail(500, 1000, 0.10) > 0)
assert.equal(binomUpperTail(0, 200, 0.1), 1)
assert.equal(binomUpperTail(201, 200, 0.1), 0)

// An absolute ceiling is vacuous on amplified material; a local excess against the sample's own
// baseline is not. A sample at 10% genome-wide with a planted 60% run must flag the run only.
{
  const n = 5_000, flags: (0 | 1)[] = []
  for (let i = 0; i < n; i++) flags.push(i % 10 === 0 ? 1 : 0)             // 10% baseline
  for (let i = 2_000; i < 2_300; i++) flags[i] = i % 10 < 6 ? 1 : 0        // local 60%
  const hits = scanNocallClusters(flags, 0.10, 200, 0.05)
  assert.ok(hits.length >= 1, 'the planted cluster is found')
  assert.ok(hits.some((h) => h.startIndex <= 2_300 && h.endIndex >= 2_000), 'and it is in the right place')
  assert.ok(hits.every((h) => h.rate > 0.10), 'every hit exceeds the baseline it was tested against')
  // Overlapping windows over one event report one finding, not one per offset.
  assert.ok(hits.length <= 3, `one event should not fan out: got ${hits.length}`)
}
// A uniformly noisy sample has no LOCAL excess and must produce nothing, which is the whole
// point: the old absolute 5% ceiling flagged every window of exactly this sample.
{
  const flags: (0 | 1)[] = []
  for (let i = 0; i < 5_000; i++) flags.push(i % 8 === 0 ? 1 : 0)          // 12.5% throughout
  assert.deepEqual(scanNocallClusters(flags, 0.125, 200, 0.05), [],
    'uniform 12.5% against a 12.5% baseline is not a finding')
}
assert.deepEqual(scanNocallClusters([1, 1, 1], 0.1, 200), [], 'fewer markers than the window')

// --- 8. gates refuse what they must, and only report what has no threshold ------------
const prof = (callRate: number, hetRate: number) => {
  const b = new Map<string, ChromStats>()
  const called = Math.round(1000 * callRate)
  for (let i = 1; i <= 22; i++) {
    b.set(String(i), { markers: 1000, called, het: Math.round(called * hetRate), nocall: 1000 - called })
  }
  b.set('X', { markers: 1000, called, het: Math.round(called * hetRate * 0.8), nocall: 1000 - called })
  return finishProfile('s', b, sums(0.04, 0.97), 'AX-1')
}
const verdict = (p: ReturnType<typeof prof>, name: string) =>
  gates(p).find((g) => g.name === name)!.verdict

assert.equal(verdict(prof(0.55, 0.30), 'call rate'), 'exclude', 'below 60%: the one published cut')
assert.equal(verdict(prof(0.68, 0.30), 'call rate'), 'marginal')
assert.equal(verdict(prof(0.907, 0.122), 'call rate'), 'usable', 'the real cohort is inside the band')
assert.equal(verdict(prof(0.99, 0.30), 'call rate'), 'usable')
// The vendor gate is named in the message rather than applied, because it was established on
// unamplified bulk DNA and applying it rejects the entire target cohort.
assert.match(gates(prof(0.907, 0.122)).find((g) => g.name === 'call rate')!.detail, /vendor gate of 97%/)

// The het-to-hom asymmetry is an approximation on amplified material, not a theorem: erroneous
// heterozygous calls become common below 60% call rate. It must be gated, never assumed.
assert.equal(verdict(prof(0.55, 0.30), 'het-to-hom asymmetry valid'), 'exclude')
assert.equal(verdict(prof(0.907, 0.122), 'het-to-hom asymmetry valid'), 'usable')

// Genome-wide LOH is qualitative in the source and stays qualitative here.
assert.equal(verdict(prof(0.90, 0.005), 'genome-wide LOH'), 'exclude')
assert.equal(verdict(prof(0.90, 0.122), 'genome-wide LOH'), 'report_only')
assert.match(gates(prof(0.90, 0.005)).find((g) => g.name === 'genome-wide LOH')!.detail,
  /qualitative and embryo-clustered/)

// --- 9. end to end on a synthetic file -----------------------------------------------
{
  const lines = [HEADER]
  for (let i = 0; i < 2_000; i++) {
    const gt = i % 10 === 0 ? '-1' : i % 3 === 0 ? '1' : i % 3 === 1 ? '0' : '2'
    const baf = gt === '1' ? '0.51' : gt === '0' ? '0.03' : gt === '2' ? '0.97' : ''
    lines.push(`AX-${i}\t${(i % 22) + 1}\t${1000 + i * 500}\t0.01\t${baf}\t2.0\t${gt}\t1`)
  }
  const p = profileText('synthetic', lines.join('\n'))!
  assert.equal(p.markers, 2_000)
  assert.equal(p.nocallRate, 0.10)
  assert.equal(p.coding.verdict, 'ok')
  assert.equal(p.idPrefix, 'AX-')
  assert.match(p.product, /product unidentified/, '2000 markers is not the UKBB array')
  // No chrX in this fixture, so sex is not callable and must say so rather than default.
  assert.equal(p.sex, 'ambiguous')
  assert.equal(p.chrXHetRatio, null)
}
assert.equal(profileText('empty', 'nonsense\nlines\n'), null, 'no locatable header: refuse the file')

// --- 10. genome build from positions alone --------------------------------------------
// No rsIDs, no manifest, no chain file. A marker cannot sit past the end of its chromosome or
// inside an assembly N-gap, so the correct build scores exactly zero illegal placements.
const at = (chrom: string, pos: number) => ({ chrom, pos })
const spread = (chrom: string, pos: number, n: number) =>
  Array.from({ length: n }, (_, i) => at(chrom, pos + i * 7))

// chr1 runs to 249,250,621 in GRCh37 but only 248,956,422 in GRCh38, so this stretch exists in
// one assembly and is off the end of the other. Positions chosen empirically to be clear of
// N-gaps in the build that contains them: picking a round number lands in an hg19 gap and reads
// as illegal under BOTH, which is a fixture bug rather than a detector bug.
{
  const b = detectBuild(spread('1', 249_081_422, 2_000))
  assert.equal(b.build, 'GRCh37')
  assert.equal(b.illegal.GRCh37, 0)
  assert.equal(b.illegal.GRCh38, 2_000, 'every one is past the GRCh38 end of chr1')
}
// And the other direction, so a detector biased toward one build cannot pass.
{
  const b = detectBuild(spread('3', 198_072_430, 2_000))
  assert.equal(b.build, 'GRCh38')
  assert.equal(b.illegal.GRCh38, 0)
  assert.equal(b.illegal.GRCh37, 2_000)
}
// 'chr1' and '1' are one chromosome; the bare-integer form carries no build information, since
// Ensembl uses it for both assemblies and the naming reflects the exporting software.
assert.equal(detectBuild(spread('chr1', 249_081_422, 2_000)).build, 'GRCh37')
assert.equal(detectBuild(spread('CHR1', 249_081_422, 2_000)).build, 'GRCh37')

// The test is ASYMMETRIC and the implementation must respect it. A nonzero count excludes a
// build at ANY marker count, because the correct-build rate is exactly zero rather than small.
// A zero count only supports a build once enough markers were examined to have expected one.
{
  const few = detectBuild(spread('2', 50_000_000, 20))
  assert.equal(few.build, null, 'both builds clean at 20 markers: undetermined, not a coin flip')
  assert.match(few.note, /too few|both builds are clean/)
  // But 20 markers past every chromosome end still excludes both, which IS informative.
  const bad = detectBuild(spread('1', 260_000_000, 20))   // past BOTH chromosome ends
  assert.equal(bad.build, null)
  assert.ok(bad.illegal.GRCh37 > 0 && bad.illegal.GRCh38 > 0)
  assert.match(bad.note, /both builds place markers illegally/)
}
// Assembly N-gaps, not just chromosome ends: chr1 opens with a 10 kb telomere gap in both.
{
  const g = detectBuild([at('1', 5_000)])
  assert.ok(g.illegal.GRCh37 > 0 && g.illegal.GRCh38 > 0, 'inside the telomere gap in both')
}
assert.equal(detectBuild([at('MT', 500)]).tested, 0, 'unrecognised chromosome is not tested')
assert.equal(detectBuild([]).build, null)

// --- 11. relatedness: verification only, and only what the statistic supports ----------
// Simulated WITH LINKAGE, because the dispersion the statistic reads is created by IBD blocks.
// Independent transmission per marker would average every window to the same value and the
// sibling signature would vanish - which is a property of the simulation, not of the method.
{
  let seed = 7
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  const NCHR = 22, PER = 4_000, XO = 1 / 2_000, WIN = 200
  const ids: string[] = []
  for (let c = 0; c < NCHR; c++) for (let i = 0; i < PER; i++) ids.push(`c${c}_${i}`)
  const hap = () => ids.map(() => (rnd() < 0.5 ? 'A' : 'B'))
  const transmit = (h1: string[], h2: string[]) => {
    const o: string[] = []
    let k = 0
    for (let c = 0; c < NCHR; c++) {
      let cur = rnd() < 0.5 ? 0 : 1
      for (let i = 0; i < PER; i++, k++) { if (rnd() < XO) cur = 1 - cur; o.push(cur ? h2[k] : h1[k]) }
    }
    return o
  }
  const geno = (x: string[], y: string[]) => new Map<string, AB>(
    ids.map((id, i) => [id, (x[i] === y[i] ? (x[i] === 'A' ? 'AA' : 'BB') : 'AB') as AB]),
  )
  const drop = (m: Map<string, AB>, d: number) => new Map<string, AB>(
    [...m].map(([k, v]) => [k, v === 'AB' && rnd() < d ? (rnd() < 0.5 ? 'AA' : 'BB') : v]),
  )
  const P1 = [hap(), hap()], P2 = [hap(), hap()], P3 = [hap(), hap()]
  const kid = (A: string[][], B: string[][]) => [transmit(A[0], A[1]), transmit(B[0], B[1])]
  const c1 = kid(P1, P2), c2 = kid(P1, P2)
  const parent = geno(P1[0], P1[1]), child = geno(c1[0], c1[1])
  const sib = geno(c2[0], c2[1]), stranger = geno(P3[0], P3[1])

  // The separation that matters holds at every dropout level: a sibling pair's dispersion stays
  // well above a parent-offspring pair's, while the genome-wide rate stops separating anything.
  for (const d of [0, 0.30, 0.63]) {
    const po = assessRelatedness(drop(parent, d), drop(child, d), ids, WIN)
    const ss = assessRelatedness(drop(sib, d), drop(child, d), ids, WIN)
    assert.equal(ss.relationship, 'sibling_or_more_distant', `siblings at dropout ${d}`)
    assert.notEqual(po.relationship, 'sibling_or_more_distant', `parent-offspring at dropout ${d}`)
    assert.ok(ss.windowedP0Sd! > po.windowedP0Sd! * 2,
      `dispersion separation at dropout ${d}: ${ss.windowedP0Sd} vs ${po.windowedP0Sd}`)
  }

  // What the statistic CANNOT do, asserted so nobody later mistakes it for something it isn't:
  // parent-offspring and unrelated both have flat dispersion, one because an allele is shared
  // everywhere and the other because nothing is shared anywhere. Under dropout the genome-wide
  // rate no longer separates them either, and a family trio cannot supply the cohort of
  // mostly-unrelated pairs that READ's normalisation needs.
  const poFlat = assessRelatedness(drop(parent, 0.63), drop(child, 0.63), ids, WIN)
  const unFlat = assessRelatedness(drop(parent, 0.63), drop(stranger, 0.63), ids, WIN)
  assert.equal(poFlat.relationship, 'flat_dispersion_parent_offspring_or_unrelated')
  assert.equal(unFlat.relationship, 'flat_dispersion_parent_offspring_or_unrelated')
  assert.match(poFlat.note, /consistent with parent-offspring AND with unrelated/)
  // At zero dropout the rate does separate them, which is why it is reported rather than hidden.
  assert.ok(assessRelatedness(parent, child, ids, WIN).ibs0Rate < 0.001)
  assert.ok(assessRelatedness(parent, stranger, ids, WIN).ibs0Rate > 0.05)

  // A duplicate must be caught before any role logic, or it reads as a first-degree relative.
  assert.equal(assessRelatedness(parent, parent, ids, WIN).relationship, 'same_individual')

  // Floors: below the degree-call floor the answer is indeterminate, not a guess.
  assert.equal(assessRelatedness(parent, child, ids.slice(0, 500), WIN).relationship, 'indeterminate')
  const small = assessRelatedness(parent, child, ids.slice(0, 5_000), WIN)
  assert.equal(small.relationship, 'indeterminate')
  assert.match(small.note, /below the 10000 floor/)
}

// --- 12. values measured on real public data, pinned as regressions -------------------
// Sources: Zuccaro 2020 GSE148488 and Turocy 2026 GSE186407, both UK Biobank Axiom, 825,656
// markers. The files are not in the repo - they are 40 MB each and are patient-derived even
// though public - so the MEASUREMENTS are pinned instead. Each of these caught a real defect.

// Known-sex bulk samples. The female expectation is ~1.0 on this array, NOT the 0.7-0.8
// predicted from the X's smaller effective population size: that prediction does not survive
// contact with this array's ascertainment.
const sexRatio = (ratio: number) => callSex(chroms(0.30, 0.30 * ratio)).sex
assert.equal(sexRatio(0.113), 'male', 'Zuccaro sperm donor rep 1, known male')
assert.equal(sexRatio(0.115), 'male', 'rep 2, same individual, same call')
assert.equal(sexRatio(1.060), 'female', 'Zuccaro egg donor A, known female')
assert.equal(sexRatio(1.144), 'female', 'egg donor C, known female')
// And the user's own sample still lands between the bands, which real data confirms is a real
// gap rather than a threshold artefact: measured males are 0.11 and females 1.06-1.14.
assert.equal(sexRatio(0.462), 'ambiguous')

// Technical replicates of the same bulk DNA concord at 95.8% on this array, not 99%. A 0.99
// duplicate threshold MISSED a genuine duplicate; the next class down (a real parent-offspring
// pair on the same data) concords at 54.9%, so 0.90 separates them with room to spare.
{
  const n = 20_000
  const ids = Array.from({ length: n }, (_, i) => `m${i}`)
  const a = new Map<string, AB>(ids.map((id, i) => [id, (i % 3 === 0 ? 'AB' : i % 3 === 1 ? 'AA' : 'BB') as AB]))
  // 95.8% concordance, and disagreements pushed to het->hom so IBS0 stays low, as dropout does.
  const b = new Map<string, AB>([...a].map(([k, v], i) => [k, i % 24 === 0 && v === 'AB' ? 'AA' : v]))
  const r = assessRelatedness(a, b, ids, 500)
  assert.ok(r.concordance > 0.90 && r.concordance < 0.99, `concordance was ${r.concordance}`)
  assert.equal(r.relationship, 'same_individual', 'real replicates must be caught')
}

// A WGA'd blastomere at 67% call rate and 56.0% heterozygosity passed every other gate as merely
// "marginal", and produced a sex call contradicting the other blastomere of the same embryo. A
// diploid genome cannot be 56% heterozygous, so an upper bound catches it where call rate did not.
{
  const bad = prof(0.670, 0.560)
  assert.equal(verdict(bad, 'heterozygosity plausible'), 'exclude')
  assert.match(gates(bad).find((g) => g.name === 'heterozygosity plausible')!.detail,
    /not attainable by a diploid genome/)
  // The call-rate gate alone called this merely marginal, which is why the new gate was needed.
  assert.equal(verdict(bad, 'call rate'), 'marginal')
  // And the samples that should pass, do.
  for (const h of [0.163, 0.153, 0.312]) {
    assert.equal(verdict(prof(0.90, h), 'heterozygosity plausible'), 'usable', `het ${h}`)
  }
}

// The degree call is REFUSED under asymmetric chemistry, because on real data it is not merely
// noisy but INVERTED: a bulk parent against a WGA'd embryo (dropout 0.00 vs 0.33) gave
// parent-offspring dispersion 0.073-0.098 while real full siblings gave 0.041. Simulation with
// equal dropout in both samples cannot see this, which is why it shipped wrong.
{
  const n = 20_000
  const ids = Array.from({ length: n }, (_, i) => `p${i}`)
  const a = new Map<string, AB>(ids.map((id, i) => [id, (i % 2 ? 'AB' : 'AA') as AB]))
  const b = new Map<string, AB>(ids.map((id, i) => [id, (i % 5 ? 'AB' : 'BB') as AB]))
  const guarded = assessRelatedness(a, b, ids, 500, { a: 0.00, b: 0.33 })
  assert.equal(guarded.relationship, 'indeterminate')
  assert.match(guarded.note, /ordering inverts/)
  // Without the dropout arguments the refusal cannot fire, so it must not be silent about that.
  assert.notEqual(assessRelatedness(a, b, ids, 500).relationship, 'indeterminate')
  // Comparable chemistry is still allowed through.
  assert.notEqual(assessRelatedness(a, b, ids, 500, { a: 0.30, b: 0.32 }).relationship,
    'indeterminate')
}

// Both real column conventions parse, from two papers by the same lab on the same array.
{
  const turocy = headerMap('probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset')!
  assert.equal(turocy.__delim, '\t')
  assert.equal(turocy['log2r'], 3)
  // Zuccaro: COMMA separated, different order, and `normalized_intensity` for the same channel.
  const zuccaro = headerMap('probeset_id,copy_number,chr,position,probe_classification,baf,genotype,normalized_intensity')!
  assert.equal(zuccaro.__delim, ',')
  assert.equal(zuccaro['log2r'], 7, 'normalized_intensity must fold onto log2r')
  assert.equal(zuccaro['copy_number'], 1)
  // Matching one spelling would have silently dropped the whole intensity channel, and phase 3
  // would then report copy-neutral LOH untestable for a reason unrelated to the sample.
  const row = parseRow('AX-13216142,2.0,1,86028,Other,0.009244147,0,0.2614052', zuccaro)!
  assert.equal(row.log2R, 0.2614052)
  assert.equal(row.genotype, 'AA')
  assert.equal(row.copyNumber, 2)
}

console.log('ingest.check.ts: all assertions passed')

// --- genotype dialects ------------------------------------------------------------------------
// AB space is the commonest spelling there is and was not read at all until 2.2.3: an ordinary
// file came back 100% no-call and was excluded on its call rate, with nothing saying the file was
// fine and only the dialect unknown. Nucleotide space is a different case and must stay refused:
// which nucleotide is allele A is a per-marker convention that needs pooling across every sample
// in the run, and assigning it per file can give two files opposite conventions at one marker.
{
  const m = headerMap('probeset_id,chr,position,genotype,baf')!
  const gt = (v: string) => parseRow(`m1,1,1000,${v},0.5`, m)!
  for (const [raw, want] of [
    ['0', 'AA'], ['1', 'AB'], ['2', 'BB'], ['-1', 'NC'],
    ['AA', 'AA'], ['AB', 'AB'], ['BA', 'AB'], ['BB', 'BB'],
    ['NC', 'NC'], ['--', 'NC'], ['---', 'NC'], ['?', 'NC'], ['.', 'NC'], ['', 'NC'],
  ] as const) {
    assert.equal(gt(raw).genotype, want, `genotype ${JSON.stringify(raw)}`)
    assert.equal(gt(raw).unreadableGenotype, false, `${JSON.stringify(raw)} is readable`)
  }
  // Nucleotide pairs read as no-call, but are FLAGGED so the reader is told which dialect it is.
  for (const raw of ['AG', 'GG', 'CT', 'TT']) {
    assert.equal(gt(raw).genotype, 'NC')
    assert.ok(gt(raw).unreadableGenotype, `${raw} must be flagged, not silently dropped`)
  }

  const profileOf = (v: (i: number) => string): SampleProfile => {
    const byChrom = new Map<string, ChromStats>()
    const baf = emptyBafSums()
    const builds = emptyBuildSums()
    for (let i = 0; i < 3000; i += 1) {
      const r = parseRow(`m${i},1,${i * 1000},${v(i)},0.5`, m)!
      accumulate(r, byChrom); accumulateBaf(r, baf); accumulateBuild(r, builds)
    }
    return { ...finishProfile('x', byChrom, baf, 'm0', builds), build: buildVerdict(builds) }
  }
  const ab = profileOf((i) => (i % 3 === 0 ? 'AB' : 'AA'))
  assert.equal(ab.callRate, 1, 'an AB-space file is fully called')
  assert.ok(!gates(ab).some((g) => g.name === 'genotype format'), 'and raises no format gate')

  const nuc = profileOf((i) => (i % 3 === 0 ? 'AG' : 'GG'))
  const fmt = gates(nuc).find((g) => g.name === 'genotype format')
  assert.equal(fmt?.verdict, 'exclude', 'a nucleotide file is excluded, not silently empty')
  assert.ok(fmt?.detail.includes('nucleotide pairs'), 'and the dialect is named')
}

// Self-check for the phase-1 informativity engine. Run: node src/informativity.check.ts
//
// The load-bearing test is section 1: every one of the 48 rows of informativity_table.csv,
// which was enumerated independently of this code. That independence is the point. A
// Mendelian classifier is easy to get consistently wrong - the same sign error applied
// everywhere passes any test written from the same understanding - so the assertion here is
// against someone else's enumeration, not against my own restatement of it.
//
// The other sections pin the failures that would quietly corrupt a real run: a father no-call
// being imputed rather than excluded, a haplotype identity being asserted with no phase
// source, and a strand-ambiguous SNP being assigned to AB space by guessing.
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import {
  abFromAlleles, capability, classifyMarker, isStrandAmbiguous, karyomappingClass, keySnp,
  summarizeWindows, type AB, type Cohort, type H3Value, type Marker, type PhaseDeclaration,
} from './informativity.ts'
import { finishProfile, type ChromStats } from './ingest.ts'
import type { SampleProfile } from './ingest.ts'

// The table records hap_M/hap_W, which presupposes a phase label. Its convention is that the
// A allele sits on the mutant homologue; supplying that is what makes those columns
// reproducible, and withholding it is what section 3 tests.
const TABLE_PHASE: PhaseDeclaration = { mutantAllele: 'A', route: 'long_read_father' }

// --- 1. the 48-row table --------------------------------------------------------------
// Fields are compared exactly where the table's vocabulary is an enum, and semantically
// where it is prose (h3_diagnostic_value), because matching prose word-for-word would test
// the wording rather than the logic.

/** Minimal RFC4180-ish split: the table quotes fields containing commas. */
function splitCsvLine(line: string): string[] {
  const out: string[] = []
  let cur = '', inQ = false
  for (let i = 0; i < line.length; i++) {
    const c = line[i]
    if (inQ) {
      if (c === '"') {
        if (line[i + 1] === '"') { cur += '"'; i++ } else inQ = false
      } else cur += c
    } else if (c === '"') inQ = true
    else if (c === ',') { out.push(cur); cur = '' }
    else cur += c
  }
  out.push(cur)
  return out
}

/** The table's prose h3 field to the enum this module reports. */
function h3FromProse(s: string): H3Value {
  const t = s.trim()
  if (t === '' || t === 'none') return 'none'
  if (t.includes('anti-H3')) return 'anti_h3'            // before POSITIVE: "POSITIVE anti-H3"
  if (t.startsWith('MIRROR')) return 'mirror'
  if (t.startsWith('POSITIVE H3')) return 'positive_h3'
  if (t.startsWith('weak')) return 'weak_diploidy_only'
  if (t.startsWith('ambiguous')) return 'ambiguous'
  throw new Error(`unclassified h3 prose: ${t}`)
}

const raw = readFileSync(new URL('../../tests/fixtures/informativity_table.csv', import.meta.url), 'utf8')
const lines = raw.split('\n').filter((l) => l.trim() !== '')
const header = splitCsvLine(lines[0])
const rows = lines.slice(1).map((l) => {
  const cells = splitCsvLine(l)
  return Object.fromEntries(header.map((h, i) => [h, cells[i] ?? ''])) as Record<string, string>
})

assert.equal(rows.length, 48, 'the table enumerates 3 father x 4 mother x 4 embryo states')

const asAB = (s: string): AB | null => (s === 'missing' ? null : (s as AB))
const tri = (s: string): boolean | null => (s === 'yes' ? true : s === 'no' ? false : null)

let checked = 0
for (const r of rows) {
  const father = r.father_gt as AB
  const mother = asAB(r.mother_gt)
  const embryo = r.embryo_gt as AB
  const where = `father ${r.father_gt} x mother ${r.mother_gt} x embryo ${r.embryo_gt}`
  const v = classifyMarker(father, mother, embryo, TABLE_PHASE)

  assert.equal(v.mendelianConsistent, tri(r.mendelian_consistent), `mendelian: ${where}`)
  assert.equal(v.violationImplicates ?? '', r.violation_implicates, `implicates: ${where}`)
  assert.equal(v.paternalAlleleDeducible ?? 'none', r.paternal_allele_deducible, `deducible: ${where}`)
  assert.equal(v.paternalHaplotypeCall ?? 'none', r.paternal_haplotype_call, `haplotype: ${where}`)
  assert.equal(v.informativeHaplotypeOrigin, tri(r.informative_haplotype_origin), `L2: ${where}`)
  assert.equal(v.paternalPresencePower, r.paternal_presence_power, `pat presence: ${where}`)
  assert.equal(v.maternalPresencePower, r.maternal_presence_power, `mat presence: ${where}`)
  assert.equal(v.adoVulnerability, r.ado_vulnerability, `ADO: ${where}`)
  assert.equal(v.karyomappingClass, r.karyomapping_class, `class: ${where}`)
  assert.equal(v.h3, h3FromProse(r.h3_diagnostic_value), `h3: ${where}`)
  checked++
}
assert.equal(checked, 48)

// --- 2. the two structural claims the table encodes -----------------------------------
// These are asserted separately because they are the reason the engine exists, and a future
// edit could keep all 48 rows passing while breaking the property they were chosen to show.

// Exactly two combinations do both jobs at once: father het x mother homozygous x embryo het.
const both = rows.filter((r) => {
  const v = classifyMarker(r.father_gt as AB, asAB(r.mother_gt), r.embryo_gt as AB, TABLE_PHASE)
  return v.informativeHaplotypeOrigin && v.provesPaternalPresence
})
assert.equal(both.length, 2, 'exactly two fully-informative combinations (spec 4)')
for (const r of both) {
  assert.equal(r.father_gt, 'AB')
  assert.ok(r.mother_gt === 'AA' || r.mother_gt === 'BB')
  assert.equal(r.embryo_gt, 'AB')
}

// Dropout converts het to hom and never the reverse, so a het embryo call is robust to every
// H3 mechanism and a hom call never is. Any drift here silently promotes maskable markers.
for (const r of rows) {
  const v = classifyMarker(r.father_gt as AB, asAB(r.mother_gt), r.embryo_gt as AB, TABLE_PHASE)
  if (r.embryo_gt === 'AB') assert.equal(v.adoVulnerability, 'no_het_call_is_robust', r.embryo_gt)
  else if (r.embryo_gt === 'NC') assert.equal(v.adoVulnerability, 'na')
  else assert.equal(v.adoVulnerability, 'yes_hom_call_maskable')
}

// key/non-key is computable from the PARENTS alone, before any embryo data exist. Same
// parents must give the same partition whatever the embryo turned out to be.
for (const father of ['AA', 'AB', 'BB'] as AB[]) {
  for (const mother of ['AA', 'AB', 'BB', null] as (AB | null)[]) {
    const partitions = new Set(
      (['AA', 'AB', 'BB', 'NC'] as AB[]).map((e) => classifyMarker(father, mother, e).keySnp),
    )
    assert.equal(partitions.size, 1, `key/non-key must not depend on the embryo: ${father} x ${mother}`)
    assert.equal([...partitions][0], keySnp(father, mother))
  }
}
// And it is exactly the fully-informative parental class.
assert.ok(keySnp('AB', 'AA') && keySnp('AB', 'BB'))
assert.ok(!keySnp('AB', 'AB'), 'both parents het: one dropout mimics the opposite phase')
assert.ok(!keySnp('AA', 'BB'), 'father homozygous: his homologues are indistinguishable')
assert.ok(!keySnp('AB', null), 'no mother: no mother-impossible allele exists')

// --- 3. haplotype identity is refused without a declared phase source -----------------
// Fixture case C11: observations identical, only the declaration differs. The refusal must be
// driven by the declaration and not by the genotypes.
const withPhase = classifyMarker('AB', 'AA', 'AB', TABLE_PHASE)
const without = classifyMarker('AB', 'AA', 'AB')
assert.equal(withPhase.paternalHaplotypeCall, 'hap_W')
assert.equal(without.paternalHaplotypeCall, 'refused_no_phase_source')
// Everything not requiring the label is unchanged by withholding it.
assert.equal(without.paternalAlleleDeducible, withPhase.paternalAlleleDeducible)
assert.equal(without.provesPaternalPresence, withPhase.provesPaternalPresence)
assert.equal(without.informativeHaplotypeOrigin, withPhase.informativeHaplotypeOrigin)
// The mutant-allele convention is honoured, not assumed: flipping it flips the call.
assert.equal(
  classifyMarker('AB', 'AA', 'AB', { mutantAllele: 'B', route: 'single_sperm' }).paternalHaplotypeCall,
  'hap_M',
)

// --- 4. father no-call is excluded, never imputed --------------------------------------
// Imputing the father from the embryos would reintroduce the circularity the whole design
// avoids: the embryo panel is the thing under test.
for (const embryo of ['AA', 'AB', 'BB', 'NC'] as AB[]) {
  const v = classifyMarker('NC', 'AA', embryo, TABLE_PHASE)
  assert.equal(v.informativeHaplotypeOrigin, false, 'father NC is never informative')
  assert.equal(v.provesPaternalPresence, false)
  assert.equal(v.paternalAlleleDeducible, null)
  assert.equal(v.h3, 'none', 'father NC yields no H3 evidence in either direction')
  assert.equal(v.mendelianConsistent, null, 'unknowable, not violated')
}

// --- 5. the missing-mother mode, and the one class that survives it -------------------
// Without the mother, L3 presence proof is gone at every marker.
for (const father of ['AA', 'AB', 'BB'] as AB[]) {
  for (const embryo of ['AA', 'AB', 'BB'] as AB[]) {
    assert.equal(
      classifyMarker(father, null, embryo, TABLE_PHASE).provesPaternalPresence, false,
      'no mother: nothing proves paternal PRESENCE',
    )
  }
}
// But paternal ABSENCE is still provable where the father is homozygous, which is what makes
// the degraded mode worth running (fixture case C8).
for (const [father, embryo] of [['AA', 'BB'], ['BB', 'AA']] as [AB, AB][]) {
  const v = classifyMarker(father, null, embryo, TABLE_PHASE)
  assert.equal(v.h3, 'positive_h3', `${father} x - x ${embryo} proves paternal absence`)
  assert.equal(v.paternalPresencePower, 'obligate_violation_paternal_allele_ABSENT')
  assert.match(v.violationImplicates ?? '', /provable without the mother/)
}
// A heterozygous father gives nothing in this mode: any embryo genotype is explainable.
for (const embryo of ['AA', 'AB', 'BB'] as AB[]) {
  assert.equal(classifyMarker('AB', null, embryo, TABLE_PHASE).mendelianConsistent, true)
}

// --- 5b. a maternal no-call is not a homozygous mother ---------------------------------
// The table enumerates mother AA/AB/BB/missing and has no maternal-NC row, so this is a rule
// the table does not supply. Treating NC as homozygous would promote the marker to
// fully-informative on the strength of a measurement that FAILED - the one error this module
// exists to prevent - because a homozygous mother is exactly what creates L3 power.
assert.equal(karyomappingClass('AB', 'NC'), 'father het, mother unavailable')
assert.notEqual(karyomappingClass('AB', 'NC'), karyomappingClass('AB', 'AA'))
assert.equal(keySnp('AB', 'NC'), false, 'a failed maternal call cannot make a key SNP')
for (const embryo of ['AA', 'AB', 'BB'] as AB[]) {
  const nc = classifyMarker('AB', 'NC', embryo, TABLE_PHASE)
  const absent = classifyMarker('AB', null, embryo, TABLE_PHASE)
  assert.equal(nc.provesPaternalPresence, false, 'no maternal constraint, so no presence proof')
  // Locally the two are the same inference, and that is the point: NC yields no more than
  // absence does. The distinction between them is reported as a tally, not as a class.
  assert.equal(nc.karyomappingClass, absent.karyomappingClass)
  assert.equal(nc.paternalPresencePower, absent.paternalPresencePower)
}
// The distinct fact is still counted per window; asserted in section 7 where the panel exists.

// --- 6. the mirror class is never pooled with H3 ---------------------------------------
// Maternal absence indicates maternal dropout, loss of the maternal homologue, or PATERNAL
// isodisomy. Pooling it with H3 would report a maternal event as evidence about the father.
const mirror = classifyMarker('AA', 'BB', 'AA', TABLE_PHASE)
assert.equal(mirror.h3, 'mirror')
assert.equal(mirror.paternalPresencePower, 'proves_paternal_presence')
assert.equal(mirror.maternalPresencePower, 'maternal_allele_ABSENT')
assert.notEqual(mirror.h3, 'positive_h3')

// --- 7. windows, spacing, and the father-no-call ceiling -------------------------------
const VARIANT = 1_000_000
// 2.5 kb spacing so the whole panel sits inside the 25 kb local window: that leaves enough
// markers per flank to test the father-no-call ceiling WITHOUT also tripping the
// two-per-flank rule, which is the only way to show the two conditions are independent.
const mk = (i: number): Marker => ({ rsid: `rs${i}`, chrom: '11', pos: VARIANT - 25_000 + i * 2_500 })
const markers = Array.from({ length: 21 }, (_, i) => mk(i))   // index 10 sits on the variant

const father = new Map<string, AB>(markers.map((m) => [m.rsid, 'AB']))
const mother = new Map<string, AB>(markers.map((m) => [m.rsid, 'AA']))
const cohort: Cohort = {
  father, mother, variantChrom: '11', variantPos: VARIANT,
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
}
const wins = summarizeWindows(markers, cohort)
const local = wins.find((w) => w.name === 'local')!
// The edited site is excluded from the marker set (spec 4.4), so 21 markers give 20.
assert.equal(local.markers, 20, 'the variant position itself is never used as a marker')
assert.equal(local.keyLower, 10)
assert.equal(local.keyUpper, 10)
// The gap spanning the excluded variant is double, and the median is unmoved by it.
assert.equal(local.medianSpacingBp, 2_500)
assert.equal(local.lowConfidence, false)

// One-sided coverage is a candidate breakpoint, not a call: it must not read as confident.
const oneSided: Cohort = {
  ...cohort,
  father: new Map<string, AB>(markers.map((m) => [m.rsid, m.pos < VARIANT ? 'AB' : 'AA'])),
}
const lopsided = summarizeWindows(markers, oneSided).find((w) => w.name === 'local')!
assert.equal(lopsided.keyUpper, 0)
assert.equal(lopsided.fatherNocallRate, 0, 'no no-calls here: the flank rule alone must fire')
assert.equal(lopsided.lowConfidence, true, 'no key markers on one flank is low-confidence')

// Father no-call above 5% makes the window low-confidence whatever else holds: a non-random
// no-call cluster is itself evidence of a structural event in the FATHER's genome, which
// would invalidate the phase. Both flanks stay well covered here, so the ceiling is what
// fires and not the flank rule.
const noisy = new Map<string, AB>(father)
noisy.set('rs3', 'NC')
noisy.set('rs4', 'NC')
const noisyWin = summarizeWindows(markers, { ...cohort, father: noisy })
  .find((w) => w.name === 'local')!
assert.equal(noisyWin.fatherNocall, 2)
assert.ok(noisyWin.fatherNocallRate > 0.05, `rate was ${noisyWin.fatherNocallRate}`)
assert.ok(noisyWin.keyLower >= 2 && noisyWin.keyUpper >= 2, 'flanks still covered')
assert.equal(noisyWin.lowConfidence, true)

// A maternal no-call is counted as its own fact (see 5b) and drops that marker from the key
// set, because a failed maternal call cannot establish the mother-impossible allele.
const moNc = new Map<string, AB>(mother)
moNc.set('rs3', 'NC')
const moNcWin = summarizeWindows(markers, { ...cohort, mother: moNc })
  .find((w) => w.name === 'local')!
assert.equal(moNcWin.motherNocall, 1)
assert.equal(moNcWin.keySnps, 19, 'the maternal-NC marker drops out of the key set')
assert.equal(moNcWin.fatherNocall, 0, 'a maternal failure is not a paternal one')
// With no mother at all there is nothing to fail, so the tally is zero rather than everything:
// "never sampled" and "sampled and failed" must not report as the same number.
const noMotherWin = summarizeWindows(markers, { ...cohort, mother: undefined })
  .find((w) => w.name === 'local')!
assert.equal(noMotherWin.motherNocall, 0, 'absent is not the same fact as failed')
assert.equal(noMotherWin.keySnps, 0, 'no mother: no key SNPs anywhere')

// --- 8. cohort composition decides which claims exist ---------------------------------
// The spec's data-availability matrix as a computation, so a missing sample produces a stated
// refusal instead of a quiet degradation.
const full = capability({
  father, mother, variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  fatherAmplification: 'unamplified_bulk', motherAmplification: 'unamplified_bulk',
  derivation: { blastocystsAttempted: 12, linesObtained: 6, failuresGenotyped: true },
  subjects: [
    { id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() },
    { id: 'B1', material: 'blastomere', amplification: 'unamplified_bulk', calls: new Map() },
    { id: 'L1', material: 'embryo_line', amplification: 'unamplified_bulk', calls: new Map() },
    { id: 'P0', material: 'parental_line_passage_matched', amplification: 'unamplified_bulk', calls: new Map() },
    { id: 'C1', material: 'control_embryo', amplification: 'unamplified_bulk', calls: new Map() },
  ],
})
assert.ok(full.haplotypeIdentity && full.paternalPresenceProof && full.cohortOutcomeRate)
assert.ok(full.cultureSubtraction && full.falseLohFloor)
assert.deepEqual(full.refusals, [], 'a complete cohort refuses nothing')

// Lines only: per-line state is reportable, the experiment's outcome RATE is not, because a
// derived line is survivorship-selected against the artefact class.
const linesOnly = capability({
  father, mother, variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  subjects: [{ id: 'L1', material: 'embryo_line', amplification: 'unamplified_bulk', calls: new Map() }],
})
assert.equal(linesOnly.cohortOutcomeRate, false)
assert.ok(linesOnly.refusals.some((r) => /outcome\s+RATE is not/.test(r)))

// No phase source: identity refused, presence untouched.
const noPhase = capability({
  father, mother, variantChrom: '11', variantPos: VARIANT,
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
})
assert.equal(noPhase.haplotypeIdentity, false)
assert.equal(noPhase.paternalPresenceProof, true)

// No mother: presence proof and MCC screening both go, together, for the same reason.
const noMother = capability({
  father, variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
})
assert.equal(noMother.paternalPresenceProof, false)
assert.equal(noMother.mccScreening, false)
assert.ok(noMother.refusals.some((r) => /ABSENCE is still provable/.test(r)))

// --- 9. AB space, and the strand hazard 2.0 inherits by asking a new question ----------
// 1.x could ignore strand because heterozygosity is strand-invariant. L3 power is an identity
// claim, so A/T and C/G sites cannot be assigned by guessing - a wrong assignment inverts
// every informativity class at that marker.
assert.ok(isStrandAmbiguous('A', 'T') && isStrandAmbiguous('C', 'G'))
assert.ok(isStrandAmbiguous('T', 'A'), 'order must not matter')
assert.ok(!isStrandAmbiguous('A', 'G') && !isStrandAmbiguous('C', 'T'))
assert.equal(abFromAlleles('AA', 'A', 'G'), 'AA')
assert.equal(abFromAlleles('AG', 'A', 'G'), 'AB')
assert.equal(abFromAlleles('GG', 'A', 'G'), 'BB')
assert.equal(abFromAlleles('GA', 'A', 'G'), 'AB', 'allele order in the file is not meaningful')
assert.equal(abFromAlleles('--', 'A', 'G'), 'NC')
assert.equal(abFromAlleles('AT', 'A', 'T'), null, 'strand-ambiguous: refuse rather than guess')
assert.equal(abFromAlleles('AC', 'A', 'G'), null, 'genotype disagrees with the declared alleles')
// Case and separators in the file are not meaningful; the manifest's case is not either.
assert.equal(abFromAlleles('ag', 'A', 'G'), 'AB')
assert.equal(abFromAlleles('A G', 'A', 'G'), 'AB')
assert.equal(abFromAlleles('AG', 'a', 'g'), 'AB')
// I and D are the manifest's own symbols for indel alleles. Stripping them turned every
// indel marker into a silent no-call, which loses real markers without saying so.
assert.equal(abFromAlleles('II', 'I', 'D'), 'AA')
assert.equal(abFromAlleles('ID', 'I', 'D'), 'AB')
assert.equal(abFromAlleles('DD', 'I', 'D'), 'BB')
// A diploid genotype is exactly two symbols. Reading the first two of three would report a
// confident answer from a string nobody parsed.
assert.equal(abFromAlleles('AGG', 'A', 'G'), null, 'three symbols is malformed, not a genotype')
assert.equal(abFromAlleles('A', 'A', 'G'), 'NC', 'a half-call is not two alleles agreeing')
assert.equal(abFromAlleles('', 'A', 'G'), 'NC')
// A manifest declaring both alleles identical is malformed: AB there would invent a
// heterozygote at a monomorphic site.
assert.equal(abFromAlleles('AA', 'A', 'A'), null)
assert.equal(abFromAlleles('AG', 'A', ''), null)

// --- 10. a window is an interval on ONE chromosome -------------------------------------
// Filtering on |pos - variantPos| alone admits a marker from another chromosome that happens
// to sit at a similar coordinate, which would then be counted as flanking and as a key SNP:
// a marker with no linkage to the locus at all, presented as evidence about it.
const mixed = [
  { rsid: 'x1', chrom: '11', pos: VARIANT - 5_000 },
  { rsid: 'x2', chrom: '7', pos: VARIANT + 5_000 },   // same coordinate, wrong chromosome
  { rsid: 'x3', chrom: 'chr11', pos: VARIANT + 6_000 },  // 'chr11' and '11' are one chromosome
]
const mixedFather = new Map<string, AB>([['x1', 'AB'], ['x2', 'AB'], ['x3', 'AB']])
const mixedMother = new Map<string, AB>([['x1', 'AA'], ['x2', 'AA'], ['x3', 'AA']])
const mixedWin = summarizeWindows(mixed, {
  ...cohort, father: mixedFather, mother: mixedMother,
}).find((w) => w.name === 'local')!
assert.equal(mixedWin.markers, 2, 'the chr7 marker is not in a chr11 window')
assert.equal(mixedWin.offChromosome, 1)
assert.equal(mixedWin.keySnps, 2, 'and it is not counted as a key SNP either')
assert.equal(mixedWin.keyUpper, 1, 'only the real chr11 marker flanks upward')

// A repeated rsID would reuse one parental genotype at two positions, inflating every count
// that follows. Keep the first, count the rest, say so.
const dup = [
  { rsid: 'd1', chrom: '11', pos: VARIANT - 1_000 },
  { rsid: 'd1', chrom: '11', pos: VARIANT + 1_000 },
]
const dupWin = summarizeWindows(dup, {
  ...cohort,
  father: new Map<string, AB>([['d1', 'AB']]),
  mother: new Map<string, AB>([['d1', 'AA']]),
}).find((w) => w.name === 'local')!
assert.equal(dupWin.markers, 1, 'one rsID is one marker')
assert.equal(dupWin.duplicateRsid, 1)

// The window bound is inclusive, which is a choice and is asserted so a change is deliberate.
const onEdge = summarizeWindows(
  [{ rsid: 'e1', chrom: '11', pos: VARIANT - 25_000 }, { rsid: 'e2', chrom: '11', pos: VARIANT + 25_001 }],
  { ...cohort, father: new Map<string, AB>([['e1', 'AB'], ['e2', 'AB']]), mother: new Map<string, AB>([['e1', 'AA'], ['e2', 'AA']]) },
).find((w) => w.name === 'local')!
assert.equal(onEdge.markers, 1, 'exactly 25 kb is inside the local window; 25 kb + 1 bp is not')

// karyomappingClass is a pure function of the parents and must agree with the table's four
// values exhaustively.
assert.equal(karyomappingClass('AB', 'AA'), 'fully informative (father het, mother homozygous)')
assert.equal(karyomappingClass('AB', null), 'father het, mother unavailable')
assert.equal(karyomappingClass('AA', 'AB'), 'non-informative (father homozygous)')
assert.match(karyomappingClass('AB', 'AB'), /^semi-informative/)

// --- 11. defects found by probing the surface the table does not cover -----------------
// Each of these produced a WRONG ANSWER, not a style complaint, and each is pinned so it
// cannot come back.

// L3 capability takes both parents. Father AA x mother AA leaves B impossible for the mother
// and unavailable to the father, so no embryo genotype there can ever prove paternal presence.
// Counting the mother alone made the commonest class on any real array look L3-capable.
const l3panel = (fa: AB, mo: AB) => summarizeWindows(markers, {
  ...cohort,
  father: new Map<string, AB>(markers.map((m) => [m.rsid, fa])),
  mother: new Map<string, AB>(markers.map((m) => [m.rsid, mo])),
}).find((w) => w.name === 'local')!
assert.equal(l3panel('AA', 'AA').informativeL3, 0, 'father cannot transmit the impossible allele')
assert.equal(l3panel('BB', 'BB').informativeL3, 0)
assert.equal(l3panel('AA', 'BB').informativeL3, 20, 'here the father CAN transmit it')
assert.equal(l3panel('AB', 'AA').informativeL3, 20)
// And the window agrees with the per-marker verdict, which is the invariant that matters.
assert.equal(classifyMarker('AA', 'AA', 'AB', TABLE_PHASE).provesPaternalPresence, false)
assert.equal(classifyMarker('AA', 'BB', 'AB', TABLE_PHASE).provesPaternalPresence, true)

// Flank support is either axis. Counting key SNPs alone made every father-homozygous panel
// permanently low-confidence - including the no-mother configuration the spec calls analysable
// and this panel, where every marker proves paternal presence.
const c8panel = l3panel('AA', 'BB')
assert.equal(c8panel.keySnps, 0, 'a homozygous father yields no key SNPs, by definition')
assert.ok(c8panel.supportLower >= 2 && c8panel.supportUpper >= 2, 'but presence evidence is evidence')
assert.equal(c8panel.lowConfidence, false, 'so this panel is analysable, not low-confidence')

// A father no-call is not a homozygous father. Reporting it as one asserts a genotype nobody
// measured - the exact mirror of the maternal bug closed in 5b.
assert.equal(karyomappingClass('NC', 'AA'), 'excluded: father not called')
assert.notEqual(karyomappingClass('NC', 'AA'), karyomappingClass('AA', 'AA'))

// The strand guard is case-folded: a manifest writing its allele columns lowercase describes
// the same site, and a case-sensitive test let exactly the sites it exists to refuse past it.
assert.ok(isStrandAmbiguous('a', 't') && isStrandAmbiguous('c', 'g') && isStrandAmbiguous('A', 't'))
assert.equal(abFromAlleles('AT', 'a', 't'), null)
assert.equal(abFromAlleles('cg', 'c', 'g'), null)

// capability counts usable calls, not map presence: an empty or all-NC map is a supplied file
// with nothing in it.
const emptyFather = capability({
  father: new Map(), mother, variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
})
assert.equal(emptyFather.haplotypeIdentity, false)
assert.equal(emptyFather.cohortOutcomeRate, false, 'nothing is claimable without the father')
assert.match(emptyFather.refusals[0], /No paternal genotypes/)
const emptyMother = capability({
  father, mother: new Map(), variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
})
assert.equal(emptyMother.paternalPresenceProof, false, 'an empty maternal map proves nothing')
assert.equal(
  capability({
    father, mother: new Map<string, AB>([['rs0', 'NC']]), variantChrom: '11', variantPos: VARIANT,
    subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
  }).mccScreening,
  false, 'nor does an all-NC one',
)

// One rule for "is there usable maternal information", applied in every branch. An earlier
// version tested `mother === null` in some and folded NC into unavailable in others, so the
// same evidential state produced two different verdicts depending on which branch ran.
for (const [fa, em] of [['AA', 'BB'], ['BB', 'AA'], ['AB', 'AB'], ['AA', 'NC']] as [AB, AB][]) {
  const viaNc = classifyMarker(fa, 'NC', em, TABLE_PHASE)
  const viaAbsent = classifyMarker(fa, null, em, TABLE_PHASE)
  assert.deepEqual(viaNc, viaAbsent, `mother NC and mother absent must agree: ${fa} x ${em}`)
}

// --- 12. the no-call rule is a LOCAL excess, not an absolute ceiling -------------------
// An absolute 5% ceiling was vacuous: real amplified material runs 8-13% genome-wide, so every
// window failed the gate while nothing was detected. The signal was always a local excess
// against the sample's OWN baseline.
const nocallWindow = (nocallCount: number, baseline: number | undefined) => {
  const fa = new Map<string, AB>(markers.map((m) => [m.rsid, 'AB' as AB]))
  // Put the no-calls on markers below the variant so the flank rule does not also fire.
  markers.filter((m) => m.pos < VARIANT).slice(0, nocallCount).forEach((m) => fa.set(m.rsid, 'NC'))
  return summarizeWindows(markers, { ...cohort, father: fa, fatherNocallBaseline: baseline })
    .find((w) => w.name === 'local')!
}
// A sample running 10% genome-wide, with a window ALSO at 10%: not a finding. This is the exact
// case the old ceiling rejected and the reason it had to go.
const atBaseline = nocallWindow(2, 0.10)
assert.ok(atBaseline.nocallExcessP !== null)
assert.ok(atBaseline.nocallExcessP! > 0.05 / 20, `p was ${atBaseline.nocallExcessP}`)
assert.ok(!atBaseline.lowConfidenceCauses.some((c) => /no-call excess/.test(c)),
  'a window at the sample baseline is not an excess')
// The same absolute rate against a LOW baseline is a finding.
const aboveBaseline = nocallWindow(8, 0.001)
assert.ok(aboveBaseline.nocallExcessP! < 0.05 / 20)
assert.ok(aboveBaseline.lowConfidenceCauses.some((c) => /no-call excess/.test(c)))
assert.equal(aboveBaseline.lowConfidence, true)
// With no baseline the tool says the baseline is missing rather than guessing a threshold.
const noBaseline = nocallWindow(4, undefined)
assert.equal(noBaseline.nocallExcessP, null)
assert.ok(noBaseline.lowConfidenceCauses.some((c) => /no genome-wide baseline/.test(c)))

// Causes are named, because the 2.4 no-call signal and the marker-support floor are different
// findings with different remedies and a bare boolean cannot tell them apart.
assert.ok(lopsided.lowConfidenceCauses.some((c) => /telomeric flank/.test(c)))
assert.ok(!lopsided.lowConfidenceCauses.some((c) => /no-call/.test(c)))
assert.deepEqual(local.lowConfidenceCauses, [], 'a good window names no causes')
assert.equal(local.lowConfidence, false)

// --- 13. phase-0 profiles gate phase-1 guarantees, and can contradict a declaration ----
const profileWith = (id: string, callRate: number, xHet: number, autoHet: number): SampleProfile => {
  const b = new Map<string, ChromStats>()
  const called = Math.round(1000 * callRate)
  for (let i = 1; i <= 22; i++) {
    b.set(String(i), { markers: 1000, called, het: Math.round(called * autoHet), nocall: 1000 - called })
  }
  b.set('X', { markers: 1000, called, het: Math.round(called * xHet), nocall: 1000 - called })
  return finishProfile(id, b, { hom0: 4, n0: 100, hom2: 97, n2: 100, het: 50, nHet: 100 }, 'AX-1')
}
const base = {
  father, mother, variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  subjects: [{ id: 'E1', material: 'embryo' as const, amplification: 'unamplified_bulk' as const, calls: new Map<string, AB>() }],
}
// Below 60% call rate the het-to-hom asymmetry is suspended, so key SNPs stop guaranteeing that
// a phase assignment survives one dropout. That must be refused, not silently relied on.
const degraded = capability({ ...base, fatherProfile: profileWith('F', 0.55, 0.24, 0.30) })
assert.ok(degraded.refusals.some((r) => /asymmetry is SUSPENDED/.test(r)))
const healthy = capability({ ...base, fatherProfile: profileWith('F', 0.91, 0.006, 0.30) })
assert.ok(!healthy.refusals.some((r) => /SUSPENDED/.test(r)), 'inside the band, no suspension')

// A declared father reading genetically female is a sample swap, and saying so is the one thing
// inference is for here. It never assigns a role; it contradicts one.
const swapped = capability({ ...base, fatherProfile: profileWith('F', 0.91, 0.24, 0.30) })
assert.ok(swapped.refusals.some((r) => /CONTRADICTION.*declared as the father/.test(r)))
const swappedMother = capability({ ...base, motherProfile: profileWith('M', 0.91, 0.006, 0.30) })
assert.ok(swappedMother.refusals.some((r) => /CONTRADICTION.*oocyte donor.*reads/.test(r)))
// An uncallable sex cannot verify OR refute, and must say that rather than pass silently.
const middling = capability({ ...base, fatherProfile: profileWith('F', 0.91, 0.14, 0.30) })
assert.ok(middling.refusals.some((r) => /not callable/.test(r)))
// A correctly declared male father produces no sex refusal at all.
assert.ok(!healthy.refusals.some((r) => /CONTRADICTION|not callable/.test(r)))

// --- 14. the two declarations that cannot be measured ----------------------------------
const cohortWith = (over: Partial<Cohort>): Cohort => ({
  father, mother, variantChrom: '11', variantPos: VARIANT, phase: TABLE_PHASE,
  fatherAmplification: 'unamplified_bulk', motherAmplification: 'unamplified_bulk',
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() }],
  ...over,
})

// Chemistry is refused rather than defaulted. Guessing bulk for an MDA sample understates
// dropout roughly sixtyfold, and guessing the other way discards a clean sample's claims.
assert.ok(capability(cohortWith({ fatherAmplification: undefined })).refusals
  .some((r) => /Amplification chemistry undeclared for father/.test(r)))
assert.ok(capability(cohortWith({
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'unknown', calls: new Map() }],
})).refusals.some((r) => /undeclared for E1/.test(r)))
assert.ok(!capability(cohortWith({})).refusals.some((r) => /Amplification/.test(r)))

// Mixed chemistry is legal and common. It is also the hazard that inverts count-based
// cross-sample comparison, so it is surfaced rather than silently tolerated.
const mixedChem = capability(cohortWith({
  subjects: [{ id: 'E1', material: 'embryo', amplification: 'mda', calls: new Map() }],
}))
assert.ok(mixedChem.refusals.some((r) => /Mixed amplification/.test(r)))
assert.ok(mixedChem.refusals.some((r) => /NOT equally reliable/.test(r)))
assert.equal(mixedChem.cohortOutcomeRate, true, 'mixed chemistry does not remove the rate, only trust')

// A cohort rate over derived lines needs the derivation attempts, or it measures survivorship.
const withLines: Partial<Cohort> = {
  subjects: [
    { id: 'E1', material: 'embryo', amplification: 'unamplified_bulk', calls: new Map() },
    { id: 'L1', material: 'embryo_line', amplification: 'unamplified_bulk', calls: new Map() },
  ],
}
const noRecords = capability(cohortWith(withLines))
assert.equal(noRecords.cohortOutcomeRate, false)
assert.ok(noRecords.refusals.some((r) => /no derivation-attempt records/.test(r)))
// Records that admit the failures were never genotyped are worse than absent: the coefficient
// is gone rather than merely missing, and the refusal has to say which.
const lostRecords = capability(cohortWith({
  ...withLines,
  derivation: { blastocystsAttempted: 12, linesObtained: 6, failuresGenotyped: false },
}))
assert.equal(lostRecords.cohortOutcomeRate, false)
assert.ok(lostRecords.refusals.some((r) => /gone rather than merely missing/.test(r)))
const good = capability(cohortWith({
  ...withLines,
  derivation: { blastocystsAttempted: 12, linesObtained: 6, failuresGenotyped: true },
}))
assert.equal(good.cohortOutcomeRate, true)
assert.ok(!good.refusals.some((r) => /derivation/.test(r)))
// Embryos only, no lines in the denominator: no records needed, nothing to correct for.
assert.equal(capability(cohortWith({})).cohortOutcomeRate, true)

console.log(`informativity.check.ts: all assertions passed (${checked} table rows)`)

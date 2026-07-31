/**
 * informativity - which markers can testify about an embryo's paternal contribution.
 *
 * Phase 1 of OriginMarker 2.0 (build spec 10.1). Pure genotype arithmetic: no fitted
 * parameter, no signal model, no network. It answers "is this marker set worth anything"
 * BEFORE any embryo data exist, which is also what decides whether the experiment can be
 * analysed at all.
 *
 * The novel part is that informativity is TWO orthogonal properties, not one (spec 4). The
 * literature's single "informative" label conflates them:
 *
 *   L2 power - does the observation say WHICH paternal homologue was transmitted?
 *              Requires the father heterozygous.
 *   L3 power - does it PROVE a paternal allele is physically present?
 *              Requires the embryo to carry an allele the mother cannot supply.
 *
 * Separating them is what makes an artefactual correction detectable, because H3 is exactly
 * the case where an L2 answer is confident and vacuous.
 *
 * Everything here works in AB space (AA/AB/BB/NC), which is what an Illumina GType column
 * already reports and what BAF is defined against. That is deliberate: 1.x's `callFromArray`
 * throws allele identity away because heterozygosity is strand-invariant, but L3 power is an
 * identity claim ("the mother cannot supply this allele"), so 2.0 crosses the line where
 * strand convention starts to matter. AB is cluster-relative and therefore strand-free;
 * nucleotide input needs the manifest and is refused for A/T and C/G sites (`abFromAlleles`).
 *
 * This module never names a haplotype "mutant" on its own. Which of the father's two
 * homologues carries the pathogenic allele is an EXTERNAL fact supplied by a phase source
 * (spec 3), and absent one the transmitted allele is reported and the haplotype identity is
 * refused (spec 2.1, fixture case C11).
 *
 * Self-check:  node src/informativity.check.ts
 */

// Explicit extension so this module runs under `node` as well as through the bundler: its
// self-check executes it directly, and node does not resolve extensionless specifiers.
import { normChrom } from './genotypes.ts'
import { binomUpperTail, type SampleProfile } from './ingest.ts'

/** A genotype in AB space. NC is informative and is never silently dropped. */
export type AB = 'AA' | 'AB' | 'BB' | 'NC'

/** What kind of material a sample is. This is not cosmetic: only embryos and blastomeres
 *  can support a cohort outcome rate, because a derived line is survivorship-selected
 *  (spec 8.9, and Zuccaro's 1/36 clones against 17/20 embryos). */
export type Material =
  | 'embryo'
  | 'blastomere'
  | 'embryo_line'
  | 'control_embryo'
  | 'control_line'
  | 'parental_line_passage_matched'

export const MEASURES_OUTCOME_RATE: ReadonlySet<Material> = new Set<Material>([
  'embryo', 'blastomere', 'control_embryo',
])

/**
 * Amplification chemistry, declared per sample because it cannot be measured from the file and
 * because the error regime it implies differs by two orders of magnitude.
 *
 * `unknown` is a legal value and triggers a refusal rather than a default. Guessing bulk when a
 * sample was MDA-amplified understates dropout roughly sixty-fold; guessing the other way throws
 * away every claim a clean sample could support. Neither error is recoverable downstream.
 */
export type Amplification =
  | 'unamplified_bulk'
  | 'mda' | 'malbac' | 'picoplex' | 'sureplex' | 'pta' | 'dop_pcr' | 'other'
  | 'unknown'

/**
 * Derivation-attempt records, required before any cohort FREQUENCY (spec 2.1, 8.9).
 *
 * Without them a rate computed over derived lines is a survivorship measurement wearing an
 * outcome rate's clothes: Zuccaro found the paternal-loss event in 1/36 ESC clones against 17/20
 * embryos, a thirty-fold depletion, and recorded that no lines came from two of the
 * apparent-wild-type blastocysts at all. If the failures were never genotyped that information
 * is simply gone, and the honest output is per-line verdicts with no cohort rate.
 */
export interface DerivationRecords {
  blastocystsAttempted: number
  linesObtained: number
  /** Were the blastocysts that failed derivation genotyped? If not, the coefficient is lost. */
  failuresGenotyped: boolean
}

/** Which of the father's homologues carries the pathogenic allele, and how that was learnt.
 *  `route` is recorded rather than checked: the spec ranks routes but the tool's job is to
 *  state which one the user declared, and to refuse when none was. */
export interface PhaseDeclaration {
  /** The allele, in AB space, sitting on the father's MUTANT homologue. */
  mutantAllele: 'A' | 'B'
  route:
    | 'long_read_father'
    | 'linked_read_father'
    | 'informative_relative'
    | 'single_sperm'
    | 'embryo_panel_external_anchor'
}

export type PresencePower =
  /** The embryo carries an allele the mother cannot supply: that allele IS paternal. */
  | 'proves_paternal_presence'
  | 'proves_maternal_presence'
  /** Two homologues are present but neither is attributable to a parent. */
  | 'diploidy_only'
  /** Mendelian-impossible with the parent's allele detected. Positive evidence of absence.
   *  The paternal side carries the longer name because it is the load-bearing one: it is the
   *  only positive-H3 reading, and the vocabulary is the table's, not invented here. */
  | 'obligate_violation_paternal_allele_ABSENT'
  | 'paternal_allele_ABSENT'
  | 'maternal_allele_ABSENT'
  /** The embryo carries an allele NEITHER parent could supply, so no parental-origin
   *  argument survives: error, contamination, or a complex event. */
  | 'both_sides_violated'
  /** Consistent with the allele being there, and equally with it having dropped out. */
  | 'consistent_with_absence'
  | 'consistent'
  | 'none'
  /** No maternal genotype, so no "the mother cannot supply this" argument exists. */
  | 'none_mother_unavailable'

export type KaryomappingClass =
  | 'fully informative (father het, mother homozygous)'
  | 'semi-informative (both parents het): informative only on homozygous embryo call'
  | 'non-informative (father homozygous)'
  | 'father het, mother unavailable'
  /** The table has no such row: a father no-call is not a homozygous father, and reporting it
   *  as one asserts a genotype nobody measured. */
  | 'excluded: father not called'

/** Semantic category of what the marker says about H3, as an enum rather than prose. */
export type H3Value =
  /** Mendelian-impossible with the paternal allele present: direct evidence of absence. */
  | 'positive_h3'
  /** The maternal allele is the one undetected. Reported separately, never pooled with H3
   *  (spec 4.3): maternal dropout, loss of the maternal homologue, or paternal isodisomy. */
  | 'mirror'
  /** The paternal allele is physically demonstrated, which is evidence AGAINST H3. */
  | 'anti_h3'
  /** Two homologues present, origin unattributable. Does not exclude maternal heterodisomy. */
  | 'weak_diploidy_only'
  /** Transmission of a mother-shared allele is indistinguishable from paternal loss. */
  | 'ambiguous'
  | 'none'

export interface MarkerVerdict {
  /** null when an NC makes consistency unknowable rather than violated. */
  mendelianConsistent: boolean | null
  violationImplicates: string | null
  /** Which paternal allele reached the embryo. 'A|B' means both are consistent. */
  paternalAlleleDeducible: 'A' | 'B' | 'A|B' | null
  /** Only ever hap_M/hap_W when a phase source was declared. */
  paternalHaplotypeCall:
    | 'hap_M' | 'hap_W' | 'ambiguous'
    | 'NA_father_homozygous'
    | 'refused_no_phase_source'
    | null
  /** L2 power. */
  informativeHaplotypeOrigin: boolean
  paternalPresencePower: PresencePower
  maternalPresencePower: PresencePower
  /** Whether a single dropout could have produced this genotype from a different truth.
   *  A het call is robust to every H3 mechanism; a hom call never is (spec 4.2). */
  adoVulnerability: 'yes_hom_call_maskable' | 'no_het_call_is_robust' | 'na'
  karyomappingClass: KaryomappingClass
  h3: H3Value
  /** L3 power: does this marker prove the paternal allele is physically present. */
  provesPaternalPresence: boolean
  /** A key SNP's phase cannot have been flipped by one dropout (Natesan 2014). Computable
   *  from the PARENTAL genotypes alone, before any embryo data exist. */
  keySnp: boolean
}

const alleles = (g: AB): string[] => (g === 'NC' ? [] : g.split(''))
const isHom = (g: AB): boolean => g === 'AA' || g === 'BB'
const carries = (g: AB, a: string): boolean => alleles(g).includes(a)

/** The allele a homozygote must transmit. */
const only = (g: AB): string => g[0]

/**
 * key/non-key, from the parents alone.
 *
 * A marker is key when the father is heterozygous and the mother is homozygous. There, the
 * embryo's informative reading is a HET call (it carries the mother-impossible allele), and
 * a het call cannot be manufactured by dropout - dropout only ever turns het into hom. So no
 * single dropout can flip the phase assignment. Where both parents are het, the informative
 * reading is a HOM call, and one dropout produces exactly the genotype the opposite phase
 * would give: non-key.
 */
export function keySnp(father: AB, mother: AB | null): boolean {
  if (father !== 'AB' || mother === null || mother === 'NC') return false
  return isHom(mother)
}

/**
 * The marker's class, from the parents alone.
 *
 * A maternal no-call is treated as maternal-unavailable, because at that marker it is: no
 * maternal constraint exists, so no "the mother cannot supply this allele" argument can be
 * made. `informativity_table.csv` enumerates mother AA/AB/BB/missing and has no maternal-NC
 * row, so this is a rule the table does not supply - the same gap the spec closes for the
 * father in 2.4, closed the same way. Calling it homozygous instead would promote the marker
 * to fully-informative on the strength of a measurement that failed, which is the one error
 * this module exists to prevent. The distinction between "she was never sampled" and "she was
 * sampled and this marker failed" is real and is reported as a separate tally, not folded into
 * the class.
 */
export function karyomappingClass(father: AB, mother: AB | null): KaryomappingClass {
  // Checked before the homozygous branch: 'not AB' is not the same as 'homozygous', and the
  // difference matters because a cluster of paternal no-calls is itself evidence of a
  // structural event in the father's genome (spec 2.4).
  if (father === 'NC') return 'excluded: father not called'
  if (father !== 'AB') return 'non-informative (father homozygous)'
  if (mother === null || mother === 'NC') return 'father het, mother unavailable'
  if (mother === 'AB') {
    return 'semi-informative (both parents het): informative only on homozygous embryo call'
  }
  return 'fully informative (father het, mother homozygous)'
}

/**
 * Classify one marker for one subject sample.
 *
 * `mother === null` is the declared degraded mode (spec 2.1 / 11 config 3), not an error.
 * It is materially weaker: without her, "the mother cannot supply this allele" is unavailable,
 * so L3 presence proof disappears at every marker. One class survives, and it is the reason
 * the degraded mode is worth running at all - a homozygous father cannot transmit a genotype
 * lacking his allele, so father-hom x embryo-opposite-hom still proves paternal ABSENCE
 * without any maternal data (spec 4.3, fixture case C8).
 */
export function classifyMarker(
  father: AB,
  mother: AB | null,
  subject: AB,
  phase?: PhaseDeclaration,
): MarkerVerdict {
  const cls = karyomappingClass(father, mother)
  const key = keySnp(father, mother)

  const base: MarkerVerdict = {
    mendelianConsistent: null,
    violationImplicates: null,
    paternalAlleleDeducible: null,
    paternalHaplotypeCall: null,
    informativeHaplotypeOrigin: false,
    paternalPresencePower: 'none',
    maternalPresencePower: 'none',
    adoVulnerability: 'na',
    karyomappingClass: cls,
    h3: 'none',
    provesPaternalPresence: false,
    keySnp: key,
  }

  // A father no-call excludes the marker outright: it is not imputed from the embryos,
  // because the embryo panel is the thing under test (spec 2.4).
  if (father === 'NC') return base

  // One rule for "is there usable maternal information here", applied everywhere. An earlier
  // version tested `mother === null` in some branches and folded 'NC' into unavailable in
  // others, so the same evidential state - no maternal genotype at this marker - produced two
  // different verdicts depending on which branch ran.
  const motherKnown = mother !== null && mother !== 'NC'

  if (subject === 'NC') {
    return { ...base, paternalPresencePower: motherKnown ? 'none' : 'none_mother_unavailable' }
  }

  const ado = alleles(subject)[0] === alleles(subject)[1]
    ? 'yes_hom_call_maskable' as const
    : 'no_het_call_is_robust' as const

  // --- what each parent could have contributed -----------------------------------------
  // A homozygote transmits one known allele. A heterozygote could transmit either. With no
  // maternal genotype there is no maternal constraint at all.
  const patCan = father === 'AB' ? ['A', 'B'] : [only(father)]
  const matCan = mother === null || mother === 'NC' ? ['A', 'B'] : mother === 'AB' ? ['A', 'B'] : [only(mother)]

  const sub = alleles(subject)
  // Every ordered (paternal, maternal) pair that could produce the observed genotype.
  const explanations: Array<[string, string]> = []
  for (const p of patCan) {
    for (const m of matCan) {
      const produced = [p, m].sort().join('')
      if (produced === sub.slice().sort().join('')) explanations.push([p, m])
    }
  }

  const consistent = explanations.length > 0

  // --- Mendelian-inconsistent readings --------------------------------------------------
  if (!consistent) {
    // Two different failures hide behind "Mendelian-inconsistent", and only one of them is
    // evidence about the paternal allele.
    //
    // If the embryo carries an allele NEITHER parent could supply, nothing about parental
    // origin survives: that is error, contamination or a complex event. Note this is not the
    // same test as "is each parent's allele present" - father AA x mother AA x embryo AB has
    // an A that either parent could have given, so both look supplied, yet the genotype is
    // still impossible because the B comes from nowhere.
    const supplyable = new Set([...patCan, ...matCan])
    const unexplained = sub.filter((a) => !supplyable.has(a))

    if (unexplained.length > 0) {
      const het = ado === 'no_het_call_is_robust'
      return {
        ...base,
        mendelianConsistent: false,
        violationImplicates: 'genotyping error, contamination, or complex event',
        // A het call still demonstrates two homologues are physically present, which is worth
        // recording even when neither is attributable. Maternal heterodisomy is not excluded.
        paternalPresencePower: het ? 'diploidy_only' : 'none',
        maternalPresencePower: 'both_sides_violated',
        adoVulnerability: ado,
        h3: het ? 'weak_diploidy_only' : 'none',
      }
    }

    // Otherwise one parent's obligate allele is missing from the embryo. Which one.
    const patSupplied = patCan.some((a) => carries(subject, a))
    if (!patSupplied) {
      // The paternal allele is not in the embryo. With a het mother she supplies at most one
      // of the two required alleles, so the second cannot be maternal either.
      const implicates = !motherKnown
        ? 'paternal allele ABSENT - provable without the mother (father homozygous cannot '
          + 'transmit the observed genotype)'
        : 'loss/dropout of paternal homologue, or maternal iso-UPD'
      return {
        ...base,
        mendelianConsistent: false,
        violationImplicates: implicates,
        paternalPresencePower: 'obligate_violation_paternal_allele_ABSENT',
        maternalPresencePower: motherKnown ? 'paternal_allele_ABSENT' : 'none',
        adoVulnerability: ado,
        h3: 'positive_h3',
      }
    }
    return {
      ...base,
      mendelianConsistent: false,
      violationImplicates: 'loss/dropout of maternal homologue, or paternal iso-UPD',
      paternalPresencePower: 'proves_paternal_presence',
      maternalPresencePower: 'maternal_allele_ABSENT',
      adoVulnerability: ado,
      h3: 'mirror',
      provesPaternalPresence: true,
    }
  }

  // --- consistent readings --------------------------------------------------------------
  const patSet = new Set(explanations.map(([p]) => p))
  const deducible = patSet.size === 1 ? ([...patSet][0] as 'A' | 'B') : ('A|B' as const)

  // L3: the embryo carries an allele the mother demonstrably cannot supply.
  const motherImpossible = motherKnown && isHom(mother)
    ? (only(mother) === 'A' ? 'B' : 'A')
    : null
  const proves = motherImpossible !== null && carries(subject, motherImpossible)

  // L2: a single paternal allele is deducible AND the father is heterozygous, so the two
  // homologues are distinguishable. A homozygous father transmits a known allele but his
  // homologues are indistinguishable, so there is nothing to identify.
  const l2 = father === 'AB' && patSet.size === 1

  let haplotype: MarkerVerdict['paternalHaplotypeCall']
  if (father !== 'AB') haplotype = 'NA_father_homozygous'
  else if (patSet.size > 1) haplotype = 'ambiguous'
  else if (!phase) haplotype = 'refused_no_phase_source'
  else haplotype = deducible === phase.mutantAllele ? 'hap_M' : 'hap_W'

  let presence: PresencePower = 'none'
  if (proves) presence = 'proves_paternal_presence'
  else if (!motherKnown) presence = 'none_mother_unavailable'
  else if (subject === 'AB') presence = 'diploidy_only'
  else if (l2 && motherKnown && isHom(mother) && carries(subject, only(mother))) {
    // The embryo is homozygous for an allele the mother also carries, so the paternal
    // contribution is inferred, not demonstrated: this reading is what paternal loss looks
    // like too.
    presence = 'consistent_with_absence'
  }

  let maternalPresence: PresencePower = 'none'
  if (motherKnown) {
    const patImpossible = father !== 'AB' ? (only(father) === 'A' ? 'B' : 'A') : null
    if (patImpossible !== null && carries(subject, patImpossible)) {
      maternalPresence = 'proves_maternal_presence'
    } else {
      maternalPresence = 'consistent'
    }
  }

  let h3: H3Value = 'none'
  if (proves) h3 = 'anti_h3'
  else if (subject === 'AB' && motherKnown) h3 = 'weak_diploidy_only'
  else if (presence === 'consistent_with_absence') h3 = 'ambiguous'

  return {
    ...base,
    mendelianConsistent: true,
    paternalAlleleDeducible: deducible,
    paternalHaplotypeCall: haplotype,
    informativeHaplotypeOrigin: l2,
    paternalPresencePower: presence,
    maternalPresencePower: maternalPresence,
    adoVulnerability: ado,
    h3,
    provesPaternalPresence: proves,
    keySnp: key,
  }
}

// --- cohort-level summary ----------------------------------------------------------------

export interface Marker {
  rsid: string
  chrom: string
  pos: number
}

export interface SubjectSample {
  id: string
  material: Material
  /** Declared, never inferred. 'unknown' is legal and refuses rather than defaulting. */
  amplification: Amplification
  /** Phase-0 profile for this sample, if it has been ingested. */
  profile?: SampleProfile
  /** Genotypes keyed by rsID, in AB space. Missing entries are treated as NC. */
  calls: Map<string, AB>
}

export interface Cohort {
  /** REQUIRED. Defines both paternal homologues and therefore all informativity. */
  father: Map<string, AB>
  /** Optional. Absent is the declared degraded mode, never entered silently. */
  mother?: Map<string, AB>
  subjects: SubjectSample[]
  phase?: PhaseDeclaration
  /** The edited site's chromosome. Required, because a window is an interval on ONE
   *  chromosome and a position alone does not say which. */
  variantChrom: string
  /** The edited position, excluded from the marker set (spec 4.4). */
  variantPos: number
  /** The FATHER's genome-wide no-call rate, from the phase-0 profile. The window test is a
   *  local excess against this, never against an absolute ceiling: on amplified material a
   *  fixed 5% limit rejects every window of a sample running 8-13% while detecting nothing. */
  fatherNocallBaseline?: number
  /** Phase-0 profiles. Roles are DECLARED, never inferred: two parents of one embryo are both
   *  exactly first-degree to it with identical kinship 0.25 and IBD1 = 1, so no autosomal
   *  statistic can order them, and the pairwise sex test INVERTS under mixed chemistry. What a
   *  profile is for is the other direction - an inferred sex can CONTRADICT a declaration, and
   *  a contradiction is a robust finding where an assignment would not be. */
  fatherProfile?: SampleProfile
  motherProfile?: SampleProfile
  /** Declared chemistry for the parental samples, same rule as the subjects. */
  fatherAmplification?: Amplification
  motherAmplification?: Amplification
  /** Required before any cohort frequency is reported. */
  derivation?: DerivationRecords
}

export interface WindowSummary {
  name: string
  halfWidthBp: number
  markers: number
  informativeL2: number
  informativeL3: number
  fullyInformative: number
  keySnps: number
  /** Per side of the variant, because a call resting on one side is a candidate crossover
   *  or breakpoint, not a segment call (spec 8.5). */
  keyLower: number
  keyUpper: number
  /** Markers able to testify on EITHER axis - key (haplotype) or L3-capable (presence).
   *  This is what the two-per-side floor counts, because presence evidence is evidence. */
  supportLower: number
  supportUpper: number
  medianSpacingBp: number | null
  fatherNocall: number
  fatherNocallRate: number
  /** Markers where the mother WAS sampled and the call failed. Locally equivalent to having
   *  no mother, but a different fact, and a cluster of them is a signal in its own right. */
  motherNocall: number
  /** Markers dropped because they sit on a different chromosome from the edited site. */
  offChromosome: number
  /** Markers dropped because an earlier marker already claimed that rsID. The rsID is the
   *  join key into the parental genotype maps, so a repeat would reuse one genotype at two
   *  positions and inflate every count that follows. */
  duplicateRsid: number
  /** Exact binomial upper tail for this window's father no-call count against the father's own
   *  genome-wide rate. null when no baseline was supplied. A LOCAL excess is the signal; a
   *  globally high rate is just the chemistry. */
  nocallExcessP: number | null
  /** Why the window is low-confidence, named rather than collapsed. The 2.4 no-call signal and
   *  the marker-support floor are different findings with different remedies, and a consumer
   *  keying off a bare boolean cannot tell them apart. */
  lowConfidenceCauses: string[]
  lowConfidence: boolean
}

export const WINDOWS: ReadonlyArray<{ name: string; halfWidthBp: number }> = [
  { name: 'local', halfWidthBp: 25_000 },
  { name: 'segmental', halfWidthBp: 10_000_000 },
  { name: 'whole_chromosome', halfWidthBp: Number.POSITIVE_INFINITY },
]

// An absolute father no-call ceiling used to live here at 0.05. It was vacuous: real amplified
// material runs 8-13% genome-wide, so every window failed the gate while nothing was detected.
// Replaced by the local-excess test in summarizeWindows. Kept only as the fallback trigger when
// no baseline is available, where it at least says the baseline is missing.

const median = (xs: number[]): number | null => {
  if (xs.length === 0) return null
  const s = [...xs].sort((a, b) => a - b)
  const m = s.length >> 1
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}

/**
 * Per-window marker adequacy, from the parents alone.
 *
 * This is the phase 1 deliverable that is useful with no embryo data at all: it says whether
 * the panel could support a conclusion, and the answer is frequently no.
 */
export function summarizeWindows(markers: Marker[], cohort: Cohort): WindowSummary[] {
  const want = normChrom(cohort.variantChrom)

  // A window is an interval on ONE chromosome. Filtering on |pos - variantPos| alone would
  // admit a marker from another chromosome that happens to sit at a similar coordinate, and
  // it would then be counted as flanking and as a key SNP - a marker with no linkage to the
  // locus at all presented as evidence about it.
  const onChrom = markers.filter((m) => normChrom(m.chrom) === want)
  const offChromosome = markers.length - onChrom.length

  // The rsID is the join key into the parental maps, so a repeated rsID reuses one genotype
  // at two positions. Keep the first and count the rest rather than silently double-count.
  const seen = new Set<string>()
  const unique = onChrom.filter((m) => (seen.has(m.rsid) ? false : (seen.add(m.rsid), true)))
  const duplicateRsid = onChrom.length - unique.length

  const usable = unique.filter((m) => m.pos !== cohort.variantPos)   // spec 4.4
  return WINDOWS.map(({ name, halfWidthBp }) => {
    const inWin = usable.filter((m) => Math.abs(m.pos - cohort.variantPos) <= halfWidthBp)
    let l2 = 0, l3 = 0, both = 0, keys = 0, keyLower = 0, keyUpper = 0, nocall = 0, moNocall = 0
    let supportLower = 0, supportUpper = 0
    for (const m of inWin) {
      const fa = cohort.father.get(m.rsid) ?? 'NC'
      if (fa === 'NC') { nocall++; continue }
      const mo = cohort.mother ? (cohort.mother.get(m.rsid) ?? 'NC') : null
      if (cohort.mother && mo === 'NC') moNocall++
      const hasL2 = fa === 'AB'
      // L3 capability takes BOTH parents. A homozygous mother creates a mother-impossible
      // allele, but that only helps if the father can actually transmit it: father AA x mother
      // AA leaves B impossible for the mother and unavailable to the father, so no embryo
      // genotype at that marker can ever prove paternal presence. Testing the mother alone
      // counted the commonest class on any real array as L3-capable.
      const impossible = mo !== null && mo !== 'NC' && isHom(mo)
        ? (only(mo) === 'A' ? 'B' : 'A')
        : null
      const hasL3 = impossible !== null && (fa === 'AB' || only(fa) === impossible)
      if (hasL2) l2++
      if (hasL3) l3++
      if (hasL2 && hasL3) both++
      const isKey = keySnp(fa, mo)
      if (isKey) keys++
      // Flank support is either axis, not key SNPs alone. A key SNP testifies about which
      // homologue came through; an L3-capable marker testifies that a paternal allele is there
      // at all. Counting only key SNPs made every father-homozygous panel permanently
      // low-confidence - including the no-mother configuration the spec calls analysable, and
      // the father-hom x mother-hom-opposite panel where every marker proves presence.
      if (isKey || hasL3) {
        if (m.pos < cohort.variantPos) supportLower++
        else supportUpper++
      }
      if (isKey) {
        if (m.pos < cohort.variantPos) keyLower++
        else keyUpper++
      }
    }
    const positions = inWin.map((m) => m.pos).sort((a, b) => a - b)
    const gaps = positions.slice(1).map((p, i) => p - positions[i])
    const rate = inWin.length ? nocall / inWin.length : 0
    // Local excess against the father's OWN baseline. Bonferroni across the markers actually
    // examined, so the correction matches the search that was performed.
    const baseline = cohort.fatherNocallBaseline
    const nocallExcessP = baseline !== undefined && baseline > 0 && baseline < 1 && inWin.length > 0
      ? binomUpperTail(nocall, inWin.length, baseline)
      : null
    const alpha = 0.05 / Math.max(1, usable.length)
    const causes: string[] = []
    if (nocallExcessP !== null && nocallExcessP < alpha) {
      causes.push(`father no-call excess: ${nocall}/${inWin.length} against a genome-wide `
        + `${(baseline! * 100).toFixed(1)}%, p=${nocallExcessP.toExponential(1)}. A non-random `
        + 'cluster is evidence of a structural event in the FATHER\'s genome, which would '
        + 'invalidate the phase.')
    }
    if (baseline === undefined && rate > 0.05) {
      causes.push(`father no-call rate ${(rate * 100).toFixed(1)}% with no genome-wide baseline `
        + 'supplied, so a local excess cannot be distinguished from the sample\'s chemistry. '
        + 'Supply the phase-0 profile.')
    }
    if (supportLower < 2) causes.push('fewer than 2 informative markers on the centromeric flank')
    if (supportUpper < 2) causes.push('fewer than 2 informative markers on the telomeric flank')
    return {
      name,
      halfWidthBp,
      markers: inWin.length,
      informativeL2: l2,
      informativeL3: l3,
      fullyInformative: both,
      keySnps: keys,
      keyLower,
      keyUpper,
      supportLower,
      supportUpper,
      medianSpacingBp: median(gaps),
      fatherNocall: nocall,
      fatherNocallRate: rate,
      motherNocall: moNocall,
      offChromosome,
      duplicateRsid,
      nocallExcessP,
      lowConfidenceCauses: causes,
      lowConfidence: causes.length > 0,
    }
  })
}

// --- what the cohort's composition permits ------------------------------------------------

export interface CohortCapability {
  /** Which paternal homologue was transmitted. Needs a declared phase source. */
  haplotypeIdentity: boolean
  /** Positive proof a paternal allele is present. Needs the mother. */
  paternalPresenceProof: boolean
  /** Maternal contamination screening. Needs the mother. */
  mccScreening: boolean
  /** The experiment's outcome RATE, as opposed to per-sample state. Needs embryos or
   *  blastomeres, plus derivation records if lines are the denominator. */
  cohortOutcomeRate: boolean
  /** Culture-acquired change separable from editing. Needs the passage-matched parental line. */
  cultureSubtraction: boolean
  /** An empirical false-LOH floor. Needs unedited controls. */
  falseLohFloor: boolean
  refusals: string[]
}

/**
 * What this cohort can and cannot support, from its composition alone.
 *
 * The spec's data-availability matrix (11) as a computation rather than a table to read, so
 * that a missing sample produces a stated refusal instead of a quiet degradation.
 */
export function capability(cohort: Cohort): CohortCapability {
  const materials = new Set(cohort.subjects.map((s) => s.material))
  // Count usable calls, not map presence. An empty or all-NC map is a supplied file with
  // nothing in it, and treating it as "the sample exists" let a cohort with no paternal
  // genotypes at all claim every capability and refuse nothing.
  const called = (m: Map<string, AB> | undefined): number =>
    m === undefined ? 0 : [...m.values()].filter((g) => g !== 'NC').length
  const fatherCalls = called(cohort.father)
  const hasMother = called(cohort.mother) > 0
  const hasPhase = cohort.phase !== undefined
  const hasOutcomeMaterial = [...materials].some((m) => MEASURES_OUTCOME_RATE.has(m))
  // Lines in the denominator make the rate a survivorship measurement unless the derivation
  // attempts are known.
  const hasLines = [...materials].some((m) => !MEASURES_OUTCOME_RATE.has(m))
  const derivationOk = cohort.derivation !== undefined && cohort.derivation.failuresGenotyped
  const rateAvailable = hasOutcomeMaterial && (!hasLines || derivationOk)
  const hasParentalLine = materials.has('parental_line_passage_matched')
  const hasControls = materials.has('control_embryo') || materials.has('control_line')

  const refusals: string[] = []
  // The father is the one required sample, so an empty paternal map is not a degraded mode -
  // nothing in this module works without it, and every other claim below is moot.
  if (fatherCalls === 0) {
    refusals.push(
      'No paternal genotypes: the father defines both homologues and therefore all '
      + 'informativity. Nothing can be computed and no claim below is available.',
    )
    return {
      haplotypeIdentity: false, paternalPresenceProof: false, mccScreening: false,
      cohortOutcomeRate: false, cultureSubtraction: false, falseLohFloor: false, refusals,
    }
  }
  if (!hasPhase) {
    refusals.push(
      'No phase source declared: which paternal homologue was transmitted is REFUSED. '
      + 'Paternal presence and copy state remain reportable.',
    )
  }
  if (!hasMother) {
    refusals.push(
      'No maternal genotype: L3 paternal-presence proof is unavailable at every marker, and '
      + 'maternal-cell-contamination screening cannot run. Paternal ABSENCE is still provable '
      + 'at father-homozygous markers.',
    )
  }
  if (!hasOutcomeMaterial) {
    refusals.push(
      'No embryos or blastomeres: per-sample state is reportable, the experiment\'s outcome '
      + 'RATE is not. Derived lines are survivorship-selected against the artefact class.',
    )
  } else if (hasLines && !derivationOk) {
    refusals.push(
      cohort.derivation === undefined
        ? 'Derived lines are in the cohort with no derivation-attempt records, so a cohort rate '
          + 'over them cannot be corrected for survivorship and is REFUSED. Per-line verdicts and '
          + 'a rate over the embryos alone remain available.'
        : 'Derivation records state the failed blastocysts were not genotyped, so the survivorship '
          + 'coefficient is unrecoverable and a cohort rate including the lines is REFUSED. That '
          + 'information is gone rather than merely missing.',
    )
  }

  // Chemistry is declared, not measured, and an undeclared one is refused rather than defaulted:
  // guessing bulk for an MDA sample understates dropout roughly sixty-fold.
  const declared: Array<[string, Amplification | undefined]> = [
    ['father', cohort.fatherAmplification],
    ...(cohort.mother ? [['oocyte donor', cohort.motherAmplification] as [string, Amplification | undefined]] : []),
    ...cohort.subjects.map((s) => [s.id, s.amplification] as [string, Amplification | undefined]),
  ]
  const undeclared = declared.filter(([, a]) => a === undefined || a === 'unknown').map(([id]) => id)
  if (undeclared.length > 0) {
    refusals.push(
      `Amplification chemistry undeclared for ${undeclared.join(', ')}. It cannot be read from the `
      + 'file and it sets the error regime, so it is refused rather than assumed. Declare '
      + '"unamplified_bulk" if that is what it was.',
    )
  }
  // Mixed chemistry within a cohort is legal and common, and it is also the specific hazard that
  // inverts count-based relatedness: with bulk parents and an amplified embryo the mother can
  // read as LESS related than the father. Per-sample fitted dropout handles the model; the
  // report has to say the samples are not equally reliable.
  const chemistries = new Set(declared.map(([, a]) => a).filter((a): a is Amplification =>
    a !== undefined && a !== 'unknown'))
  if (chemistries.size > 1) {
    refusals.push(
      `Mixed amplification within the cohort (${[...chemistries].join(', ')}). Supported, but the `
      + 'samples are NOT equally reliable: a het-to-hom conversion means far less in the noisier '
      + 'one, and any count-based cross-sample comparison can invert. Dropout is fitted per '
      + 'sample so the model accounts for it; the interpretation must too.',
    )
  }
  if (!hasParentalLine) {
    refusals.push(
      'No passage-matched unedited parental line: culture-acquired change cannot be separated '
      + 'from an editing outcome.',
    )
  }
  if (!hasControls) {
    refusals.push('No unedited controls: no empirical false-LOH floor for this platform and lab.')
  }

  // The het-to-hom asymmetry that key SNPs rest on is an approximation on amplified material,
  // not a theorem: erroneous heterozygous calls become common below 60% call rate, and gain of
  // heterozygosity has been measured directly. Below that the key/non-key partition loses its
  // guarantee, so it must be gated rather than assumed.
  const profiles = [
    cohort.fatherProfile, cohort.motherProfile, ...cohort.subjects.map((s) => s.profile),
  ].filter((p): p is SampleProfile => p !== undefined)
  const degraded = profiles.filter((p) => p.callRate < 0.60)
  if (degraded.length > 0) {
    refusals.push(
      `Call rate below 60% in ${degraded.map((p) => p.id).join(', ')}: the het-to-hom asymmetry `
      + 'is SUSPENDED there, so the key/non-key partition no longer guarantees a phase '
      + 'assignment survives a single dropout. Those samples cannot carry an L2 claim.',
    )
  }

  // Verification, not assignment. A declared father who reads genetically female, or a declared
  // oocyte donor who reads male, is a sample swap or a mislabelled file.
  if (cohort.fatherProfile?.sex === 'female') {
    refusals.push(
      `CONTRADICTION: the sample declared as the father (${cohort.fatherProfile.id}) has a `
      + `chrX/autosome heterozygosity ratio of ${cohort.fatherProfile.chrXHetRatio?.toFixed(3)}, `
      + 'inside the female band. Check for a sample swap before anything else.',
    )
  }
  if (cohort.motherProfile?.sex === 'male') {
    refusals.push(
      `CONTRADICTION: the sample declared as the oocyte donor (${cohort.motherProfile.id}) reads `
      + `male (chrX/autosome ratio ${cohort.motherProfile.chrXHetRatio?.toFixed(3)}). Check for a `
      + 'sample swap before anything else.',
    )
  }
  if (cohort.fatherProfile?.sex === 'ambiguous') {
    refusals.push(
      `Sex of the declared father (${cohort.fatherProfile.id}) is not callable: chrX/autosome `
      + `ratio ${cohort.fatherProfile.chrXHetRatio?.toFixed(3) ?? 'n/a'} falls between the bands, `
      + 'so the declaration cannot be verified either way. A chrX abnormality, a pooled file and '
      + 'contamination all land here.',
    )
  }

  return {
    haplotypeIdentity: hasPhase,
    paternalPresenceProof: hasMother,
    mccScreening: hasMother,
    cohortOutcomeRate: rateAvailable,
    cultureSubtraction: hasParentalLine,
    falseLohFloor: hasControls,
    refusals,
  }
}

// --- nucleotide input --------------------------------------------------------------------

/** A/T and C/G SNPs are strand-ambiguous: the complement of the pair is the pair. Identity
 *  cannot be recovered from the genotype alone, and L3 power is an identity claim. */
export const isStrandAmbiguous = (a: string, b: string): boolean => {
  // Case-folded: a manifest writing its allele columns lowercase describes the same site, and
  // a case-sensitive test would let exactly the sites this guard exists to refuse slip past it.
  const p = [a.trim().toUpperCase(), b.trim().toUpperCase()].sort().join('')
  return p === 'AT' || p === 'CG'
}

/**
 * A nucleotide genotype to AB space, given the manifest's A/B allele definitions.
 *
 * Returns null rather than guessing when the site is strand-ambiguous or the genotype does
 * not match the declared alleles. A wrong AB assignment silently inverts every informativity
 * class at that marker, so guessing here is worse than dropping the marker.
 */
export function abFromAlleles(
  genotype: string,
  aAllele: string,
  bAllele: string,
): AB | null {
  const A = aAllele.trim().toUpperCase(), B = bAllele.trim().toUpperCase()
  // A manifest declaring both alleles the same is malformed. Returning AB for it would invent
  // a heterozygote at a monomorphic site.
  if (!A || !B || A === B) return null
  if (isStrandAmbiguous(A, B)) return null

  // I and D are the manifest's own symbols for an insertion and a deletion allele, so they
  // are kept: stripping them turned every indel marker into a silent no-call.
  const g = genotype.trim().toUpperCase().replace(/[^ACGTID]/g, '')
  if (g.length === 0) return 'NC'
  // A diploid genotype is exactly two allele symbols. One is a half-call, which is not two
  // alleles agreeing; three or more is malformed, and reading only the first two would report
  // a confident answer from a string nobody understood.
  if (g.length === 1) return 'NC'
  if (g.length > 2) return null

  const n = [g[0], g[1]].filter((x) => x === A).length
  const m = [g[0], g[1]].filter((x) => x === B).length
  if (n + m !== 2) return null
  return n === 2 ? 'AA' : m === 2 ? 'BB' : 'AB'
}

/**
 * The full taxonomy of chromosomal abnormality this platform can and cannot see, and the detectors
 * for the classes that were previously missing.
 *
 * WHY A TAXONOMY AND NOT JUST MORE DETECTORS. Until now this tool found three things: a parental
 * absence, a copy loss and a copy gain. Everything else an embryologist might be looking at was
 * simply not looked for, and silence read as absence. A reader had no way to tell "this array has
 * no uniparental disomy" from "this tool has never checked for uniparental disomy". So every class
 * is enumerated here, including the ones that cannot be answered, and each carries what it rests on
 * and what would be needed to reach it. A class that is impossible on this material says so rather
 * than not appearing.
 *
 * THE ORIGIN MACHINERY IS SHARED, NOT DUPLICATED. Every finding produced here is scored for
 * parental origin through the same posterior, the same four bands and the same class-inversion veto
 * as the events that came before it, because a taxonomy whose new classes carried a different kind
 * of confidence would be worse than no taxonomy: a reader compares rows.
 *
 * WHAT IS INHERITED. The detectability of each class, the floors, and the three limits that cannot
 * be lifted come from a methods review of this platform and its literature
 * (audit/CONSULT-calibrated-origin.txt, section 3). The thresholds carry their sources inline. The
 * arithmetic and the guards are computed here and checked in abnormalities.check.ts.
 */
import type { Material } from './dosageOrigin.ts'

export type AbnormalityClass =
  | 'monosomy'
  | 'trisomy'
  | 'segmental-deletion'
  | 'segmental-duplication'
  | 'cnn-loh'
  | 'isodisomy'
  | 'heterodisomy'
  | 'segmental-upd'
  | 'triploidy'
  | 'haploidy'
  | 'complex'
  | 'gamete-de-novo'
  | 'reverse-segregation'
  | 'tandem-vs-inserted'

export interface TaxonomyEntry {
  cls: AbnormalityClass
  label: string
  /** Whether this tool attempts it at all, and if not, why not. */
  detectable: 'yes' | 'partly' | 'no'
  /** The measurement it rests on. */
  by: string
  /** What a parental origin for it rests on, or why none is available. */
  origin: string
  /** Stated where the answer is a hard limit rather than a floor that better data would clear. */
  limit?: string
}

/**
 * Every class, including the four that cannot be answered.
 *
 * The impossible entries are the point of the table rather than an appendix to it. Three of them
 * are hard limits that no array of this platform clears at any quality, and naming them stops both
 * a reader and a future contributor from chasing them.
 */
export const TAXONOMY: readonly TaxonomyEntry[] = [
  {
    cls: 'monosomy', label: 'Whole-chromosome monosomy', detectable: 'yes',
    by: 'call-rate collapse with the direction of the log2R shift',
    origin: 'allele dosage, self-referenced against the array\'s own genome',
  },
  {
    cls: 'trisomy', label: 'Whole-chromosome trisomy', detectable: 'yes',
    by: 'call-rate collapse with the direction of the log2R shift',
    origin: 'allele dosage, plus band occupancy at parent-heterozygous markers for the mechanism: '
      + 'both parental homologues reported positively, the alternative only by exclusion',
  },
  {
    cls: 'segmental-deletion', label: 'Segmental deletion', detectable: 'yes',
    by: 'per-chromosome scan against an external null, at 12 Mb and above',
    origin: 'allele dosage over the interval',
  },
  {
    cls: 'segmental-duplication', label: 'Segmental duplication', detectable: 'yes',
    by: 'per-chromosome scan against an external null, at 12 Mb and above',
    origin: 'allele dosage over the interval, with the class-inversion veto where a small gain and '
      + 'a small loss are both consistent',
  },
  {
    cls: 'cnn-loh', label: 'Copy-neutral loss of heterozygosity', detectable: 'yes',
    by: 'heterozygosity depletion against the array\'s own background, with copy number unchanged',
    origin: 'allele dosage. This is the LARGEST-signal class, deviation f/2 against f/(4-2f) for a '
      + 'loss, because copy number never changes and the caller never degrades',
  },
  {
    cls: 'isodisomy', label: 'Uniparental isodisomy, whole chromosome', detectable: 'yes',
    by: 'a long run of homozygosity spanning the chromosome, with copy number two',
    origin: 'allele dosage, and the untransmitted-haplotype channel identifies WHICH homologue',
  },
  {
    cls: 'heterodisomy', label: 'Uniparental heterodisomy, whole chromosome', detectable: 'no',
    by: 'nothing available from the embryo alone: copy number is two and heterozygosity is normal',
    origin: 'not reachable',
    limit: 'A HARD LIMIT, not a floor. Both parents plus the embryo are required, and 32% of '
      + 'confirmed uniparental disomy shows no significant homozygous stretch at all, so absence '
      + 'of a stretch is not evidence of absence of disomy.',
  },
  {
    cls: 'segmental-upd', label: 'Segmental uniparental disomy', detectable: 'partly',
    by: 'a run of homozygosity above the reporting length, with copy number unchanged',
    origin: 'allele dosage over the run',
    limit: 'Threshold-bound and overlapping the normal population: runs of this length occur '
      + 'without disomy, and consanguinity produces them across many chromosomes at once. A '
      + 'consanguinity flag accompanies the call rather than being folded into it.',
  },
  {
    cls: 'triploidy', label: 'Triploidy', detectable: 'yes',
    by: 'genome-wide allele-fraction band structure: mass at one third and two thirds with the '
      + 'half band vacated, which a diploid genome cannot produce',
    origin: 'the PLOIDY is called with no parent at all. Saying whether the extra set is maternal '
      + 'or paternal needs parental DNA',
  },
  {
    cls: 'haploidy', label: 'Haploidy', detectable: 'yes',
    by: 'genome-wide heterozygosity at or below the haploid ceiling',
    origin: 'which parent contributed the single set, at 0.9911 measured over 112 chromosomes',
  },
  {
    cls: 'complex', label: 'Complex or chaotic genome', detectable: 'yes',
    by: 'the share of autosomes deviating from the array\'s own centre',
    origin: 'NONE AT ANY CONFIDENCE. Every origin statistic here is self-referenced against the '
      + 'array\'s own genome, and a genome this disturbed has no undisturbed part to reference '
      + 'against. The event is reported and no parent is named',
  },
  {
    cls: 'gamete-de-novo', label: 'Segmental change introduced by the gamete', detectable: 'partly',
    by: 'a segmental event whose origin resolves to one parent',
    origin: 'the same as any segmental event',
    limit: 'Separating a gamete-borne change from a post-zygotic one needs a SECOND UNIT of the '
      + 'same embryo. Uniformity across units is the only channel measured to break the tie: '
      + 'meiotic 64 of 64 uniform, mitotic 6 of 7 non-uniform. Genotype alone cannot do it.',
  },
  {
    cls: 'reverse-segregation', label: 'Reverse segregation', detectable: 'no',
    by: 'nothing: in 20 of 26 observed cases the copy number was NORMAL',
    origin: 'not reachable',
    limit: 'A NAMED BLIND SPOT. Both homologues separate sister chromatids at the first division '
      + 'and non-sisters at the second, leaving pericentromeric heterozygosity that the classical '
      + 'rule reads as a first-division error when the physical event was equational. It was the '
      + 'most frequent non-canonical pattern across 23 complete meioses. No marker density fixes '
      + 'this: the rule is wrong rather than underpowered.',
  },
  {
    cls: 'tandem-vs-inserted', label: 'Tandem versus inserted duplication', detectable: 'no',
    by: 'nothing: neither the allelic nor the intensity channel carries positional information '
      + 'about where a duplicated segment landed',
    origin: 'not applicable',
    limit: 'A HARD LIMIT of the platform rather than of this tool.',
  },
]

export const taxonomyFor = (cls: AbnormalityClass): TaxonomyEntry | undefined =>
  TAXONOMY.find((t) => t.cls === cls)

// ---------------------------------------------------------------------------------------------
// Thresholds, each with the measurement behind it.

/**
 * Length at which a homozygous run is worth reporting, in megabases.
 *
 * From the clinical cytogenetics literature rather than chosen here. Below it, runs are common in
 * the normal population and reporting them would generate far more false uniparental disomy than
 * real.
 */
export const LCSH_REPORT_MB = 13.5

/**
 * Share of confirmed uniparental disomy that shows NO significant homozygous stretch.
 *
 * Carried so it can be stated wherever a negative result is shown. A run-based detector cannot see
 * heterodisomy at all, so "no run found" must never be presented as "no disomy".
 */
export const UPD_NO_STRETCH_RATE = 0.32

/**
 * The false-positive generator for anything built on heterozygosity depletion.
 *
 * Measured: 1% of windows with no event present already read 40% het-depleted, from dropout that
 * clusters rather than scattering. So a depletion threshold set at or below 0.40 would call
 * copy-neutral LOH on one window in a hundred of a clean array. The shipped threshold sits above
 * it, and the gap between them is the whole margin this detector has.
 */
export const CLUSTERED_DROPOUT_DEPLETION = 0.40
export const CLUSTERED_DROPOUT_RATE = 0.01

/** Depletion a window must reach before copy-neutral LOH is called, set clear of the above. */
export const LOH_DEPLETION = 0.65

/**
 * Share of deviant autosomes above which NO origin is called, at any confidence.
 *
 * Every origin statistic in this project is self-referenced: a chromosome is measured against the
 * rest of the array's own genome. Past this point there is no undisturbed remainder to reference
 * against, and the statistic measures the reference as much as the target. All 70 whole-chromosome
 * losses in the series this was measured on sit here, which is why they were never assignable.
 */
export const COMPLEX_DEVIANT_FRACTION = 0.30
/** Call rate below which the same applies, for the same reason. */
export const COMPLEX_CALL_RATE = 0.70

/** Allele-fraction bands a triploid genome occupies, and the one it vacates. */
export const TRIPLOID_BANDS = [1 / 3, 2 / 3] as const
export const DIPLOID_BAND = 0.5
export const BAND_HALF_WIDTH = 0.08
/** Occupancy at the two triploid bands, above which a third set is present. */
export const TRIPLOID_OCCUPANCY = 0.30
/** Occupancy at the half band, below which a diploid genome is excluded. */
export const DIPLOID_VACATED = 0.12
/** Markers below which the occupancy statistic is not attempted. */
export const MIN_MARKERS_PLOIDY = 400

/**
 * Runs of homozygosity across this many chromosomes at once, which reads as consanguinity.
 *
 * Uniparental disomy affects one chromosome. Shared parental ancestry produces long homozygous
 * runs on many, and the two are indistinguishable one chromosome at a time. So the count is what
 * separates them, and where it is reached the disomy calls are flagged rather than withdrawn: the
 * finding is still real, its interpretation is what changes.
 */
export const CONSANGUINITY_CHROMOSOMES = 4

// ---------------------------------------------------------------------------------------------
// Detectors.

export interface Finding {
  cls: AbnormalityClass
  chrom: string
  startBp: number
  endBp: number
  /** True where the finding spans the whole chromosome rather than an interval within it. */
  wholeChromosome: boolean
  /** What was measured, in the units it was measured in. */
  evidence: string
  /** Set where the class itself blocks an origin call, independently of the evidence. */
  originBlocked?: string
  /** Set where the finding is real but its interpretation is contested by something else. */
  flag?: string
}

export interface WindowStat {
  chrom: string
  startBp: number
  endBp: number
  /**
   * Whether this window spans the whole chromosome, which decides which detection floor applies.
   *
   * NOT COSMETIC. A 12 Mb-scale interval on amplified material has NO floor at any mosaic fraction
   * with one parent, while a whole chromosome on the same material has 0.399. Reporting a
   * whole-chromosome copy-neutral event as a segment sent every one of them to the segment floor
   * and came back not-evaluable, so the class was detected and then refused an origin for a reason
   * that did not apply to it.
   */
  wholeChromosome?: boolean
  /** Markers called in the window, and how many of those are heterozygous. */
  called: number
  het: number
  /** Mean window log2R, self-referenced. Used only to keep copy number out of the LOH call. */
  logR?: number
}

/**
 * Copy-neutral loss of heterozygosity: heterozygosity gone while the DNA stayed.
 *
 * BOTH HALVES ARE REQUIRED AND THE SECOND IS WHAT MAKES IT COPY-NEUTRAL. Depleted heterozygosity
 * on its own is equally the signature of a deletion, which removes one copy and every heterozygote
 * with it. The intensity channel is what separates them, so a window whose log2R has moved is not
 * called here: it belongs to the deletion detector. Where no intensity is supplied the class is
 * reported as loss-or-copy-neutral rather than guessed at.
 */
export function detectLoh(
  windows: readonly WindowStat[],
  opts: { backgroundHet?: number; depletion?: number; logRTolerance?: number } = {},
): Finding[] {
  const need = opts.depletion ?? LOH_DEPLETION
  const tol = opts.logRTolerance ?? 0.15
  const usable = windows.filter((w) => w.called > 0)
  if (usable.length < 3) return []
  // The array's own background, so a globally low-heterozygosity sample is not called end to end.
  const background = opts.backgroundHet
    ?? (usable.reduce((a, w) => a + w.het, 0) / usable.reduce((a, w) => a + w.called, 0))
  if (!(background > 0)) return []
  const out: Finding[] = []
  for (const w of usable) {
    const rate = w.het / w.called
    const depletion = 1 - rate / background
    if (!(depletion >= need)) continue
    const copyMoved = w.logR !== undefined && Number.isFinite(w.logR) && Math.abs(w.logR) > tol
    if (copyMoved) continue     // a deletion, and the deletion detector owns it
    out.push({
      cls: 'cnn-loh', chrom: w.chrom, startBp: w.startBp, endBp: w.endBp,
      wholeChromosome: w.wholeChromosome ?? false,
      evidence: `heterozygosity ${(100 * rate).toFixed(1)}% against this array's own `
        + `${(100 * background).toFixed(1)}%, a depletion of ${(100 * depletion).toFixed(0)}% over `
        + `${w.called} called markers`
        + (w.logR !== undefined
          ? `, with window log2R at ${w.logR.toFixed(3)} inside the ${tol} band that keeps copy `
            + 'number unchanged'
          : '. NO INTENSITY was supplied, so this is loss-or-copy-neutral rather than established '
            + 'as copy-neutral'),
      flag: `dropout that clusters rather than scattering already puts `
        + `${(100 * CLUSTERED_DROPOUT_RATE).toFixed(0)}% of event-free windows at `
        + `${(100 * CLUSTERED_DROPOUT_DEPLETION).toFixed(0)}% depletion, which is the false-positive `
        + 'floor this threshold sits above',
    })
  }
  return out
}

export interface RunOfHomozygosity {
  chrom: string
  startBp: number
  endBp: number
  markers: number
  /** True where the run covers essentially the whole chromosome. */
  wholeChromosome: boolean
}

/**
 * Uniparental disomy from runs of homozygosity, whole-chromosome and segmental.
 *
 * WHAT THIS CANNOT DO IS THE MORE IMPORTANT HALF. It sees ISOdisomy, where the two copies are the
 * same homologue and therefore identical. HETEROdisomy leaves copy number at two and heterozygosity
 * normal, so nothing here touches it, and 32% of confirmed uniparental disomy shows no significant
 * run at all. A negative result from this function therefore means "no run was found", never "no
 * disomy is present", and it is returned in those words.
 *
 * Consanguinity produces the same runs on many chromosomes at once. Rather than fold that into the
 * threshold, runs are counted across chromosomes and the calls are flagged when the count is
 * reached: the runs are real either way, it is what they imply that changes.
 */
export function detectUpd(
  runs: readonly RunOfHomozygosity[],
  opts: { reportMb?: number; consanguinityChromosomes?: number } = {},
): Finding[] {
  const minMb = opts.reportMb ?? LCSH_REPORT_MB
  const many = opts.consanguinityChromosomes ?? CONSANGUINITY_CHROMOSOMES
  const long = runs.filter((r) => (r.endBp - r.startBp) / 1e6 >= minMb)
  const chroms = new Set(long.map((r) => r.chrom))
  const consanguineous = chroms.size >= many
  return long.map((r) => {
    const mb = ((r.endBp - r.startBp) / 1e6).toFixed(1)
    return {
      cls: r.wholeChromosome ? 'isodisomy' : 'segmental-upd',
      chrom: r.chrom, startBp: r.startBp, endBp: r.endBp, wholeChromosome: r.wholeChromosome,
      evidence: `a homozygous run of ${mb} Mb over ${r.markers} markers, past the ${minMb} Mb at `
        + 'which such runs are worth reporting. This is ISOdisomy: the two copies are the same '
        + `homologue. Heterodisomy leaves heterozygosity normal and is invisible here, and `
        + `${(100 * UPD_NO_STRETCH_RATE).toFixed(0)}% of confirmed uniparental disomy carries no `
        + 'significant run at all, so the absence of one is not the absence of disomy',
      flag: consanguineous
        ? `runs of this length are present on ${chroms.size} chromosomes at once, which is the `
          + 'pattern of shared parental ancestry rather than of disomy. The runs are real; whether '
          + 'they mean disomy is what this puts in doubt'
        : undefined,
    }
  })
}

/**
 * Find runs of consecutive homozygous markers, which is what a disomy detector needs as input.
 *
 * TOLERATES A FEW HETEROZYGOTES INSIDE A RUN, and that tolerance is not slack. On amplified
 * material a genuine homozygous stretch still throws occasional heterozygous calls from allele
 * drop-in, measured at 4.35% of truly homozygous markers over 113 same-genome pairs. A run-finder
 * that broke on the first heterozygote would find nothing on exactly the material this tool is for.
 * The tolerance is a rate rather than a count, so a long run is not held to the same absolute
 * budget as a short one.
 */
export function runsOfHomozygosity(
  markers: readonly { chrom: string; pos: number; het: boolean }[],
  opts: { tolerance?: number; minMarkers?: number; chromEndBp?: Map<string, number> } = {},
): RunOfHomozygosity[] {
  const tol = opts.tolerance ?? 0.05
  const minMarkers = opts.minMarkers ?? 100
  const byChrom = new Map<string, { pos: number; het: boolean }[]>()
  for (const m of markers) {
    if (!byChrom.has(m.chrom)) byChrom.set(m.chrom, [])
    byChrom.get(m.chrom)!.push({ pos: m.pos, het: m.het })
  }
  const out: RunOfHomozygosity[] = []
  for (const [chrom, ms] of byChrom) {
    ms.sort((a, b) => a.pos - b.pos)
    let start = 0
    let hets = 0
    for (let i = 0; i <= ms.length; i += 1) {
      const over = i === ms.length || (hets + (ms[i].het ? 1 : 0)) / (i - start + 1) > tol
      if (over) {
        const n = i - start
        if (n >= minMarkers) {
          const span = ms[i - 1].pos - ms[start].pos
          const chromEnd = opts.chromEndBp?.get(chrom)
          out.push({
            chrom, startBp: ms[start].pos, endBp: ms[i - 1].pos, markers: n,
            // Whole-chromosome only when the run covers essentially all of what was assayed, which
            // is what separates an isodisomy from a long segmental run.
            wholeChromosome: chromEnd !== undefined
              ? span >= 0.9 * chromEnd
              : n >= 0.9 * ms.length,
          })
        }
        if (i === ms.length) break
        start = i
        hets = ms[i].het ? 1 : 0
      } else if (ms[i].het) hets += 1
    }
  }
  return out
}

/**
 * Ploidy from allele-fraction band occupancy, which needs no parent at all.
 *
 * A diploid genome puts its heterozygotes at one half. A triploid puts them at one third and two
 * thirds and VACATES the half band, and that vacancy is the discriminating half: occupancy at a
 * third alone is reachable by a noisy diploid, while a genuinely empty half band is not.
 *
 * The PARENT of the extra set is a different question and this does not answer it. Digynic and
 * diandric triploidy are indistinguishable from band structure, and separating them needs parental
 * DNA. Reported as a ploidy finding with the origin explicitly blocked rather than left blank.
 */
export function detectTriploidy(
  bafs: readonly number[],
  opts: { minMarkers?: number; halfWidth?: number; occupancy?: number; vacated?: number } = {},
): Finding | null {
  const need = opts.minMarkers ?? MIN_MARKERS_PLOIDY
  const w = opts.halfWidth ?? BAND_HALF_WIDTH
  const usable = bafs.filter((b) => Number.isFinite(b) && b > 0.1 && b < 0.9)
  if (usable.length < need) return null
  let tri = 0
  let di = 0
  for (const b of usable) {
    if (TRIPLOID_BANDS.some((c) => Math.abs(b - c) <= w)) tri += 1
    if (Math.abs(b - DIPLOID_BAND) <= w) di += 1
  }
  const atTri = tri / usable.length
  const atDi = di / usable.length
  if (!(atTri >= (opts.occupancy ?? TRIPLOID_OCCUPANCY))) return null
  if (!(atDi <= (opts.vacated ?? DIPLOID_VACATED))) return null
  return {
    cls: 'triploidy', chrom: 'genome', startBp: 0, endBp: 0, wholeChromosome: true,
    evidence: `${(100 * atTri).toFixed(1)}% of intermediate readings sit at the one-third and `
      + `two-thirds bands and only ${(100 * atDi).toFixed(1)}% at the half band, over `
      + `${usable.length} markers. A diploid genome populates the half band and cannot vacate it`,
    originBlocked: 'the ploidy is established without any parent, but whether the extra set is '
      + 'maternal or paternal cannot be read from band structure. Digynic and diandric triploidy '
      + 'are indistinguishable here and separating them needs parental DNA',
  }
}

/**
 * A genome too disturbed to reference anything against, which blocks every origin call on it.
 *
 * NOT A QUALITY GATE AND THE DIFFERENCE MATTERS TO A READER. A quality failure says this array was
 * poor. This says the array may be perfect and the GENOME is the thing that is disturbed, so the
 * self-referenced statistics every origin call here depends on have no undisturbed remainder to
 * measure against. It is reported as a finding in its own right rather than as a refusal.
 */
export function detectComplex(
  deviantAutosomes: number, totalAutosomes: number, callRate: number,
  opts: { deviantFraction?: number; callRate?: number } = {},
): Finding | null {
  const maxDev = opts.deviantFraction ?? COMPLEX_DEVIANT_FRACTION
  const minCall = opts.callRate ?? COMPLEX_CALL_RATE
  const frac = totalAutosomes > 0 ? deviantAutosomes / totalAutosomes : NaN
  const byDeviance = Number.isFinite(frac) && frac > maxDev
  const byCallRate = Number.isFinite(callRate) && callRate < minCall
  if (!byDeviance && !byCallRate) return null
  return {
    cls: 'complex', chrom: 'genome', startBp: 0, endBp: 0, wholeChromosome: true,
    evidence: byDeviance
      ? `${deviantAutosomes} of ${totalAutosomes} autosomes deviate from this array's own centre, `
        + `a share of ${frac.toFixed(2)} over the ${maxDev} at which self-reference fails`
      : `call rate ${callRate.toFixed(3)} is under the ${minCall} at which self-reference fails`,
    originBlocked: 'every origin statistic here is measured against the rest of this array\'s own '
      + 'genome, and this genome has no undisturbed part left to measure against. No parent is '
      + 'named on any event of this array, at any confidence. This is a property of the genome '
      + 'rather than a fault of the array',
  }
}

/**
 * Whether a segmental change came from the gamete or arose after fertilisation.
 *
 * ONLY ANSWERABLE WITH A SECOND UNIT OF THE SAME EMBRYO, and that is the whole content of this
 * function. Genotype cannot separate the two at any material quality: the measured discriminator is
 * uniformity across independently sampled units, meiotic 64 of 64 uniform against mitotic 6 of 7
 * non-uniform. The tool cannot ask for a second biopsy, but it can use one the user already has,
 * so this costs nothing at runtime and is the only route to the answer.
 */
export function callUniformity(
  presentInUnits: number, totalUnits: number,
): { mechanism: 'meiotic' | 'post-zygotic' | 'unresolved'; why: string } {
  if (totalUnits < 2) {
    return {
      mechanism: 'unresolved',
      why: 'only one unit of this embryo was loaded. Whether a segmental change came from the '
        + 'gamete or arose after fertilisation is not separable on genotype at any material '
        + 'quality; uniformity across independently sampled units is the only channel measured to '
        + 'break the tie. A second array of the same embryo answers it and costs nothing at runtime',
    }
  }
  if (presentInUnits === totalUnits) {
    return {
      mechanism: 'meiotic',
      why: `present in all ${totalUnits} units of this embryo. A change carried by the gamete is `
        + 'in every cell; one arising after fertilisation is in a lineage. Measured: meiotic 64 of '
        + '64 uniform, post-zygotic 6 of 7 non-uniform',
    }
  }
  return {
    mechanism: 'post-zygotic',
    why: `present in ${presentInUnits} of ${totalUnits} units of this embryo, so it is confined to `
      + 'a lineage rather than carried by the gamete. Measured: post-zygotic 6 of 7 non-uniform, '
      + 'meiotic 64 of 64 uniform',
  }
}

/**
 * The classes this tool looked for and cannot answer, so silence is never read as absence.
 *
 * Returned alongside the findings rather than buried in documentation. A report that lists what was
 * found tells a reader nothing about what was never checked, and on this platform the difference
 * between "no uniparental heterodisomy present" and "uniparental heterodisomy is unreachable from
 * an embryo alone" is the difference between a result and a misunderstanding.
 */
export const unanswerable = (twoParents: boolean, units: number): TaxonomyEntry[] =>
  TAXONOMY.filter((t) => {
    // THE CONDITIONAL WALLS ARE TESTED FIRST, and the order is the whole point. Heterodisomy is
    // marked undetectable because it is undetectable FROM THE EMBRYO ALONE; a second parental array
    // lifts it. Testing the blanket flag before the condition would report it as a permanent limit
    // to a user who had already supplied what it takes to clear it, which is the opposite of
    // actionable.
    if (t.cls === 'heterodisomy') return !twoParents
    if (t.cls === 'gamete-de-novo') return units < 2
    return t.detectable === 'no'
  })

/** Classes whose origin cannot be named even where the class itself is detected. */
export const ORIGIN_UNREACHABLE: ReadonlySet<AbnormalityClass> = new Set<AbnormalityClass>([
  'triploidy', 'complex', 'heterodisomy', 'reverse-segregation', 'tandem-vs-inserted',
])

/**
 * Material classes on which a copy-neutral event is the EASIEST rather than the hardest to assign.
 *
 * Kept here because it inverts the intuition the rest of this file might otherwise leave. A
 * copy-neutral event displaces the parent-allele share by f/2 against f/(4-2f) for a loss, the
 * largest of the three classes at every fraction, precisely because copy number never changes and
 * so the genotype caller never degrades.
 */
export const CNN_LOH_IS_EASIEST = true

export type { Material }

// ---------------------------------------------------------------------------------------------
// Units of one embryo, and what having more than one settles.

/**
 * Genotype concordance above which two arrays are the SAME individual rather than two.
 *
 * Not chosen here. It is the figure this project already measured for catching a parent labelling
 * accident: replicate arrays of one person concord at 95.8% of shared called markers, and a
 * parent-offspring pair at 54.9%. Two biopsies of one embryo are the same genome, so they sit with
 * the replicates. Siblings from the same two parents sit far below, which is what makes this
 * separable without asking the user to label anything.
 */
export const SAME_EMBRYO_AGREEMENT = 0.90

/**
 * Group arrays into embryos by genotype concordance, so nobody has to declare the grouping.
 *
 * ASKING THE USER WOULD HAVE BEEN EASIER AND WORSE. A mislabelled group produces a uniformity call
 * that is confidently wrong about mechanism, and this is measurable rather than declarable: the
 * separation between the same genome and a sibling is 41 points wide. `agree` is passed in so this
 * module does not have to own a genotype representation.
 */
export function groupUnits<T>(
  units: readonly T[],
  agree: (a: T, b: T) => number,
  threshold = SAME_EMBRYO_AGREEMENT,
): T[][] {
  const groups: T[][] = []
  for (const u of units) {
    // Join the first group whose members all read as the same genome. Requiring ALL rather than
    // ANY keeps a chain of near-misses from merging two embryos through an intermediate.
    const hit = groups.find((g) => g.every((m) => {
      const a = agree(m, u)
      return Number.isFinite(a) && a > threshold
    }))
    if (hit) hit.push(u)
    else groups.push([u])
  }
  return groups
}

export interface EventLocation { chrom: string; startBp: number; endBp: number }

/** Whether two intervals on the same chromosome overlap at all. */
export const overlaps = (a: EventLocation, b: EventLocation): boolean =>
  a.chrom === b.chrom && a.startBp <= b.endBp && b.startBp <= a.endBp

/**
 * How many units of one embryo carry an event overlapping this one.
 *
 * Overlap rather than an exact match, because two biopsies of one embryo do not place a breakpoint
 * identically: the interval is measured to about 8.7 kb where only informative markers count, and
 * the events being compared span tens of megabases. Requiring identical edges would report every
 * genuinely uniform event as non-uniform, which inverts the answer.
 */
export const unitsCarrying = (
  event: EventLocation, perUnit: readonly (readonly EventLocation[])[],
): number => perUnit.filter((evs) => evs.some((e) => overlaps(e, event))).length

/** Markers per sliding window when copy-neutral LOH is scanned below chromosome scale. */
export const LOH_SEGMENT_MARKERS = 600

/**
 * Collapse the redundancy that overlapping windows and a whole-chromosome pass produce.
 *
 * BOTH SCALES ARE SCANNED ON PURPOSE and that leaves duplicates to resolve. The whole-chromosome
 * pass is the one with a detection floor on amplified material, since a 12 Mb-scale interval has
 * none at any mosaic fraction with one parent; the sliding windows are what give a partial event
 * its own extent instead of reporting the entire chromosome. So the rule is: keep the widest
 * interval covering a position, and drop a segment that is only its chromosome restated.
 *
 * The wider interval is kept rather than the narrower because it is the one that can be scored: a
 * precise extent that comes back not-evaluable is worth less to a reader than a coarse one that
 * carries a parent and a confidence.
 */
export function mergeLoh(findings: readonly Finding[]): Finding[] {
  const byChrom = new Map<string, Finding[]>()
  for (const f of findings) {
    if (!byChrom.has(f.chrom)) byChrom.set(f.chrom, [])
    byChrom.get(f.chrom)!.push(f)
  }
  const out: Finding[] = []
  for (const [, fs] of byChrom) {
    // Widest first, so a narrower one is only kept where it reaches somewhere none of the wider
    // ones did.
    const sorted = fs.slice().sort((a, b) => (b.endBp - b.startBp) - (a.endBp - a.startBp))
    const kept: Finding[] = []
    for (const f of sorted) {
      const covered = kept.some((k) => k.startBp <= f.startBp && f.endBp <= k.endBp)
      if (!covered) kept.push(f)
    }
    // Adjacent survivors from half-overlapping windows are one event seen twice; join them.
    kept.sort((a, b) => a.startBp - b.startBp)
    for (const f of kept) {
      const prev = out[out.length - 1]
      if (prev && prev.chrom === f.chrom && !prev.wholeChromosome && !f.wholeChromosome
        && f.startBp <= prev.endBp) {
        out[out.length - 1] = { ...prev, endBp: Math.max(prev.endBp, f.endBp) }
      } else out.push(f)
    }
  }
  return out
}

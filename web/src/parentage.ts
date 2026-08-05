/**
 * parentage - did this sample inherit a paternal genome, and which parts of one.
 *
 * The browser half of `origin.parental_origin`. It runs here rather than on the server for the
 * same reason the carrier's genotypes do: these are a family's arrays, and the terms promise
 * nothing about a family is submitted or retained. There is no endpoint to send them to.
 *
 * Three measurements, and each one's reference is derived from the data in front of it rather
 * than being a constant fitted elsewhere:
 *
 *   absence     how often the parent's obligate allele is missing, against the rate THIS
 *               sample's own noise can manufacture. Dropout only produces paternal absence by
 *               turning a heterozygous call homozygous and discarding the parent's allele, so
 *               that bound is no-call rate times heterozygous fraction, and neither alone: a
 *               homozygous genome is immune however much drops out.
 *   other parent  how often the sample carries an allele the parent cannot supply, against half
 *               the parent's own heterozygosity, which is what a second parent contributes.
 *   zygosity    the fraction of B-allele frequencies in the heterozygous band, which says
 *               whether the genome is diploid at all. Dropout does not empty that band.
 *
 * The arithmetic is duplicated from the Python deliberately. NOTE, from an independent audit:
 * this docstring used to claim both sides are pinned by a shared fixture so a divergence fails a
 * test. That is not implemented. `parentage.check.ts` pins this side against measured values, and
 * `origin.py` has its own suite, but nothing compares the two, and the audit found them differing
 * on all 16 pairs it tried, including a 15-chromosome disagreement on one sample. The per-locus
 * layer IS pinned that way (`runlength.check.ts` loads the shared fixture); the parentage layer
 * is not. Do not restate the guarantee until there is a test behind it.
 */
import { type ProbeRow } from './ingest.ts'
import { type AB } from './informativity.ts'

/** Residual absence on clean data, from genotyping error alone. Measured at 0.03% and 0.05%. */
export const ABSENCE_ERROR_FLOOR = 0.005

/** How far above the explainable level an absence must sit before it is called. */
export const ABSENCE_MARGIN = 3.0

/** Diploid genomes run 15-16% in the heterozygous BAF band, uniparental ones 1.3-3.4%. */
export const HET_BAND_DIPLOID = 0.08

/**
 * Call rate below which heterozygous calls stop being trustworthy, so nothing derived from them
 * may be asserted. Natesan 2014, the only published threshold measured on amplified material,
 * and the same figure `ingest.gates` already excludes on.
 *
 * Absence is unaffected and still called: it is Mendelian, and a homozygous father must transmit
 * his allele whatever the call rate. What is withheld is ZYGOSITY, which is read from the
 * heterozygous band, and with it the androgenetic/biparental split that depends on it. Measured
 * on three isolated paternal pronuclei at 53.8%, 55.6% and 59.1% call rate: each is haploid and
 * therefore homozygous by construction, each showed an 18-27% heterozygous band that a haploid
 * genome cannot produce, and each was called biparental on the strength of it.
 */
export const CALL_RATE_FLOOR = 0.60

/**
 * Parent heterozygosity below which the marker set is not a polymorphic panel.
 *
 * Every rate here is calibrated on common-SNP arrays, where a parent runs 15 to 19 per cent
 * heterozygous. A whole-genome variant callset runs far lower because most of its sites are rare:
 * measured at 3.2% on a 1000 Genomes chr22 VCF, where the parent is homozygous reference at 95%
 * of markers against 75% on an array of the same chromosome. Those sites cannot show absence in
 * anyone, so they dilute the denominator and pull an unrelated pair toward a related one: the
 * same pair that reads 4.1% absence on the array reads 0.69% in the callset.
 *
 * The gap between the two populations is wide and empty, which is why a single figure can sit in
 * it. Nothing is adjusted on the strength of it; the run is annotated so the number is read
 * against the right calibration.
 */
export const PANEL_HET_FLOOR = 0.10

/**
 * Above this the heterozygous band is not measuring zygosity, and the noise ceiling built from
 * it bounds nothing.
 *
 * A diploid genome reads 15 to 16%. Dropout inflates that: Zuccaro A8 sits at 22.2% on a 53.0%
 * call rate and is still read correctly, present against its own father and refused against two
 * unrelated egg donors. Arrays measured at 36.9 to 43.0% are failures, and their ceilings near
 * 19% swallowed the absence of genomes they were unrelated to and reported them as present.
 *
 * The bound is on the band rather than the call rate, because the call rate does not separate
 * those two cases: A8 and the failed arrays sit within three points of each other on it.
 */
export const HET_BAND_IMPLAUSIBLE = 0.30

/**
 * Coefficient of variation of the per-chromosome absence rate below which the difference is
 * UNIFORM across the genome rather than confined to part of it.
 *
 * Measured on the shipped arrays: an unrelated adult 0.112, a 50:50 blend of two unrelated
 * genomes 0.104, a degraded but true parent-offspring pair 0.762, two arrays of one person 1.106.
 * The gap between 0.11 and 0.76 is empty, which is why a single figure can sit in it.
 *
 * Uniformity alone identifies nothing: an unrelated adult and a blend are indistinguishable by it
 * (0.112 against 0.104). It becomes informative only together with WHERE the rate sits. A blend
 * lands between the present and absent expectations and stays there on every chromosome; a real
 * genome with a lost segment is patchy, with some chromosomes clean.
 */
export const UNIFORM_CV = 0.35

export type OriginClass = 'androgenetic' | 'gynogenetic' | 'biparental' | 'unclear'

/** One line per class, so the card, the table and the PDF cannot describe a call differently. */
export const GLOSS: Record<OriginClass, string> = {
  androgenetic: 'Paternal genome only. No maternal complement.',
  gynogenetic: 'No paternal contribution. Maternal in origin.',
  biparental: 'Both parents contributed.',
  unclear: 'Not separated by the available evidence.',
}
export type Verdict = 'parent_genome_present' | 'no_parental_contribution' | 'unclear'
export type Zygosity = 'diploid' | 'uniparental_homozygous' | 'unknown'
export type SpermType = 'X_bearing' | 'Y_bearing' | 'unknown'

/**
 * Pseudoautosomal boundaries, from the GRC assembly region reports.
 *
 * PAR1 and PAR2 sit on both the X and the Y and recombine between them, so a Y-bearing sperm
 * still delivers paternal PAR alleles. That makes the PAR a positive control on the very sample
 * where the rest of chrX is legitimately absent: pooled into one chrX bucket it was invisible.
 */
const PAR: Record<string, [number, number][]> = {
  GRCh37: [[60001, 2699520], [154931044, 155260560]],
  GRCh38: [[10001, 2781479], [155701383, 156030895]],
}

/** False when the assembly is undetermined: PAR1 is nearly the same address in both builds but
 *  PAR2 is not, so guessing would put ordinary chrX markers into a positive control. */
export const inPar = (pos: number, build: string | null): boolean =>
  (PAR[build ?? ''] ?? []).some(([a, b]) => pos >= a && pos <= b)

const AUTOSOME = /^(?:[1-9]|1[0-9]|2[0-2])$/
export const isAutosome = (c: string): boolean => AUTOSOME.test(c)

/**
 * The highest absence rate a TRUE parent-offspring pair can produce by noise alone.
 *
 * Bounds every related pair measured, and tightly: 0.13% predicted against 0.05% observed,
 * 0.18% against 0.16%, 10.4% against 9.69%. It replaces comparing against an "unrelated" rate,
 * which cannot be a constant because it depends on the ancestry the two people share and ranged
 * from 6.8% to 50% across the pairs measured.
 */
export function absenceExplainable(noCallRate: number, hetFraction: number): number {
  if (!Number.isFinite(noCallRate) || !Number.isFinite(hetFraction)) return NaN
  return Math.max(0, noCallRate) * Math.max(0, hetFraction)
}

/**
 * What a second parent contributes on the allele-the-parent-lacks axis.
 *
 * Derived, not fitted. A second parent supplies such an allele only where the first is
 * homozygous, at the sum over markers of p^2q + q^2p = pq, and the parent's own heterozygosity
 * is the sum of 2pq over those same markers. So it is exactly half of it.
 */
export const secondParentSignal = (parentHeterozygosity: number): number =>
  Number.isFinite(parentHeterozygosity) ? parentHeterozygosity / 2 : NaN

/** Running totals for one sample against one parent. Counters only: memory stays flat. */
export interface Tally {
  /** Per chromosome: [informative markers, markers missing the parent's allele]. */
  byChrom: Map<string, [number, number]>
  nonParental: number
  nonParentalDen: number
  called: number
  het: number
  markers: number
  bafInBand: number
  bafTotal: number
  /** chrY markers the sample called, and chrY markers the parent's file carries. A Y-bearing
   *  sperm is measured from these, never inferred from an absent X. */
  yCalled: number
  yTotal: number
  /** Assembly, so the pseudoautosomal boundaries used are the right ones. GRCh37 if unknown. */
  build: string | null
}

export const emptyTally = (): Tally => ({
  byChrom: new Map(), nonParental: 0, nonParentalDen: 0,
  called: 0, het: 0, markers: 0, bafInBand: 0, bafTotal: 0, yCalled: 0, yTotal: 0, build: null,
})

/** One marker of the sample, against the parent's call at the same marker. */
export function tallyRow(parent: AB, row: ProbeRow, t: Tally): void {
  t.markers += 1
  if (row.chrom === 'Y' || row.chrom === '24') {
    // Whether the SAMPLE carries a Y is a property of the sample alone, so the denominator is
    // chrY probes on its array. Gating this on the parent's own call collapsed it to zero
    // whenever the parent was female, since she calls no chrY, and a Y-bearing sample then read
    // as not Y-bearing against an oocyte donor.
    t.yTotal += 1
    if (row.genotype !== 'NC') t.yCalled += 1
  }
  // Counted BEFORE the no-call check, to match the Python: a failed genotype still has an
  // intensity reading, and a dropped heterozygote sits mid-band, which is what separates a
  // homozygous genome from a diploid one that lost calls. Excluding them moved the band from
  // 1.27% to 0.92% on a real sample, and the noise ceiling with it.
  if (row.baf !== null && isAutosome(row.chrom)) {
    t.bafTotal += 1
    if (row.baf >= 0.35 && row.baf <= 0.65) t.bafInBand += 1
  }
  if (row.genotype === 'NC') return
  t.called += 1
  if (row.genotype === 'AB') t.het += 1
  if (parent === 'NC') return

  if (isAutosome(row.chrom)) {
    t.nonParentalDen += 1
    // An allele the parent does not carry at all had to come from somewhere else. This is the
    // only measure separating a parent-only genome from a biparental one: the parent's alleles
    // are present in both, so absence cannot tell them apart.
    for (const allele of row.genotype) if (!parent.includes(allele)) { t.nonParental += 1; break }
  }
  if (parent !== 'AA' && parent !== 'BB') return
  // chrX is split, because its two halves answer different questions: the pseudoautosomal region
  // is delivered by a sperm of either type, the rest only by an X-bearing one.
  const key = (row.chrom === 'X' || row.chrom === '23') && inPar(row.pos, t.build)
    ? 'X:PAR' : row.chrom
  const d = t.byChrom.get(key) ?? [0, 0]
  d[0] += 1
  if (!row.genotype.includes(parent[0])) d[1] += 1
  t.byChrom.set(key, d)
}

export interface ChromResult {
  chrom: string
  informative: number
  absent: number
  rate: number
  verdict: 'present' | 'absent' | 'expected_absent' | 'unclear'
  note?: string
}

export interface ParentageResult {
  verdict: Verdict
  originClass: OriginClass
  zygosity: Zygosity
  spermType: SpermType
  genomeRate: number
  explainable: number
  informative: number
  nonParentalRate: number
  secondParentExpected: number
  hetBand: number
  noCallRate: number
  /** Second factor of the ceiling, beside noCallRate: a dropped call only fakes absence where
   *  the genotype was heterozygous, so the ceiling is their product plus the error floor. */
  hetFraction: number
  /** Spread of the per-chromosome rate: low is uniform, high is confined to part of the genome. */
  dispersion: number
  /** The cleanest chromosome. Near zero means some of the genome is untouched. */
  minChromRate: number
  chroms: ChromResult[]
  notes: string[]
  limits: string[]
}

/** Turn the tallies into the three axes and the class they imply. */
export function classify(
  t: Tally,
  parentHeterozygosity: number,
  opts: { role?: 'paternal' | 'maternal'; spuriousAbsence?: number } = {},
): ParentageResult {
  const role = opts.role ?? 'paternal'
  const notes: string[] = []
  const limits: string[] = []

  let nTot = 0
  let aTot = 0
  for (const [c, [n, a]] of t.byChrom) if (isAutosome(c) && n >= 200) { nTot += n; aTot += a }
  const genomeRate = nTot ? aTot / nTot : NaN
  const noCallRate = t.markers ? 1 - t.called / t.markers : NaN
  const hetBand = t.bafTotal ? t.bafInBand / t.bafTotal : NaN
  const gtHet = t.called ? t.het / t.called : NaN
  const hetFraction = Number.isFinite(hetBand) ? hetBand : gtHet
  // A measured parental array contributes no absence of its own, so the ceiling is the sample's
  // noise alone. A RECONSTRUCTED reference does: at a contaminated marker the parent is really
  // heterozygous and a haploid product carries the other allele half the time, which reads as an
  // opposite homozygote through no fault of the sample. That floor is known rather than guessed
  // (`Reference.spuriousAbsence`), and leaving it out judges every true product against a ceiling
  // roughly its own size, so true products read unclear. Diploid samples need no such term and
  // pass nothing: their added absence carries a second factor of the maternal allele frequency,
  // and the diploid path is calibrated to 3.37x against 3.35x on the real array without it.
  const explainable = absenceExplainable(noCallRate, hetFraction) + ABSENCE_ERROR_FLOOR
    + (Number.isFinite(opts.spuriousAbsence ?? 0) ? opts.spuriousAbsence ?? 0 : 0)
  const nonParentalRate = t.nonParentalDen ? t.nonParental / t.nonParentalDen : NaN
  // The axis is DERIVED from the parent's own heterozygosity, so a parent showing none cannot
  // supply it. A genotype reconstructed from haploid meiotic products is homozygous everywhere
  // by construction, which drives the expectation to zero and makes any nonzero rate read as
  // biparental on no evidence. Withholding it routes to the existing unclear branch.
  const secondParentExpected =
    Number.isFinite(parentHeterozygosity) && parentHeterozygosity >= PANEL_HET_FLOOR
      ? secondParentSignal(parentHeterozygosity) : NaN

  if (Number.isFinite(parentHeterozygosity) && parentHeterozygosity < PANEL_HET_FLOOR) {
    limits.push(
      `The parent is heterozygous at ${pct(parentHeterozygosity, 1)} of called autosomal markers. `
      + 'Every rate here is calibrated on common-SNP arrays, where that figure runs 15 to 19%. A '
      + 'marker set this monomorphic is a whole-genome variant callset rather than a polymorphic '
      + 'panel: most of its sites are rare, they cannot show absence in anyone, and they dilute '
      + 'the denominator. Related and unrelated pairs move closer together, so a call here is '
      + 'weaker than the same call on an array and an unrelated pair may read as unclear. '
      + 'Restrict such a file to common variants before relying on it.',
    )
  }

  // The gate `ingest.gates` already applies, applied to the two axes that depend on it: below
  // this call rate the heterozygous band is artefact, so zygosity is withheld, and the noise
  // ceiling grows large enough that "under the ceiling" stops meaning "related".
  const callRate = t.markers ? t.called / t.markers : NaN
  const hetUsable = !Number.isFinite(callRate) || callRate >= CALL_RATE_FLOOR

  let verdict: Verdict
  if (Number.isFinite(hetBand) && hetBand > HET_BAND_IMPLAUSIBLE) {
    verdict = 'unclear'
    limits.push(
      `The heterozygous band is ${pct(hetBand, 1)}. A diploid genome reads 15 to 16% and dropout `
      + `inflates that into the high twenties, so ${pct(HET_BAND_IMPLAUSIBLE, 0)} is above `
      + 'anything a real genome produces: this array has failed rather than merely degraded. The '
      + 'band is one of the two factors in the noise ceiling, so the ceiling here bounds nothing '
      + 'and neither presence nor absence is reported. Repeat the array.',
    )
  } else if (!nTot) {
    verdict = 'unclear'
    limits.push('No markers where the parent is homozygous and the sample is called.')
  } else if (!Number.isFinite(explainable)) {
    verdict = 'unclear'
    limits.push('Sample quality could not be measured, so there is no bound on noise-driven '
      + 'absence and no call is made.')
  } else if (genomeRate <= explainable) {
    verdict = 'parent_genome_present'
  } else if (genomeRate >= explainable * ABSENCE_MARGIN) {
    verdict = 'no_parental_contribution'
  } else {
    verdict = 'unclear'
    limits.push(
      `Absence of ${pct(genomeRate)} sits between what this sample's own noise can explain `
      + `(${pct(explainable)}) and ${ABSENCE_MARGIN}x that. A noisy diploid sample can reach the `
      + 'unrelated range against its own parent, so this is left uncalled rather than guessed. '
      + 'The other parent\'s array would measure dropout directly and settle it.',
    )
  }

  let zygosity: Zygosity = 'unknown'
  if (!hetUsable) {
    limits.push(
      `Call rate ${pct(callRate, 1)} is below ${pct(CALL_RATE_FLOOR, 0)}, where erroneous `
      + 'heterozygous calls become common. Zygosity is read from the heterozygous band, so it is '
      + 'not reported here and the genome cannot be separated into uniparental or biparental. '
      + 'The presence or absence of this parent\'s contribution is unaffected: that is Mendelian '
      + 'and does not rest on heterozygous calls.',
    )
  } else if (Number.isFinite(hetBand)) {
    zygosity = hetBand > HET_BAND_DIPLOID ? 'diploid' : 'uniparental_homozygous'
  } else if (Number.isFinite(gtHet) && Number.isFinite(parentHeterozygosity)) {
    zygosity = gtHet > parentHeterozygosity / 2 ? 'diploid' : 'uniparental_homozygous'
    limits.push('No B-allele frequencies in this file, so zygosity comes from genotype '
      + 'heterozygosity rather than the BAF band. That is the weaker of the two measures.')
  }

  const present = verdict === 'parent_genome_present'
  let originClass: OriginClass
  if (zygosity === 'uniparental_homozygous') {
    // One allele per locus. If those alleles are the parent's there is no room for a second
    // complement, so this settles it without consulting the narrower axis.
    originClass = present ? 'androgenetic'
      : verdict === 'no_parental_contribution' ? 'gynogenetic' : 'unclear'
  } else if (zygosity === 'diploid' && present) {
    if (!Number.isFinite(secondParentExpected)) originClass = 'unclear'
    else if (nonParentalRate > secondParentExpected) originClass = 'biparental'
    else {
      originClass = 'androgenetic'
      notes.push(
        `Diploid but carrying ${pct(nonParentalRate)} alleles the parent lacks, below the `
        + `${pct(secondParentExpected)} a second parent would contribute. Consistent with a `
        + 'parent-only genome that is heterozygous rather than duplicated, meaning two sperm. '
        + 'This axis separates by about 1.6x against thirty-fold for absence, so treat a '
        + 'near-boundary call as provisional.',
      )
    }
  } else if (verdict === 'no_parental_contribution') {
    originClass = 'gynogenetic'
  } else {
    originClass = 'unclear'
  }

  const chroms: ChromResult[] = []
  let spermType: SpermType = 'unknown'
  const yBearing = t.yTotal > 0 && t.yCalled / t.yTotal > 0.5
  for (const [c, [n, a]] of [...t.byChrom].sort(byChromName)) {
    if (n < 200) continue
    const rate = a / n
    const sex = c === 'X' || c === '23'  // 'X:PAR' is deliberately not sex here
    // A male offspring has no paternal X at all: his father sent a Y instead. Flagging that as a
    // loss would report ordinary sex determination as an anomaly on half of all samples. A
    // MOTHER transmits an X either way, so the same exemption must not apply to her.
    // Gated on a chrY MEASUREMENT, not on the absent X itself. Inferring "he sent a Y" from a
    // missing paternal X is circular, and on a sample with no chrY at all it exempts a real
    // paternal X loss as ordinary sex determination while asserting a Y the file does not show.
    const xLost = rate > explainable * ABSENCE_MARGIN
    const expected = sex && role === 'paternal' && present && xLost && yBearing
    if (sex && present) {
      spermType = xLost ? (yBearing ? 'Y_bearing' : 'unknown') : 'X_bearing'
    }
    chroms.push({
      chrom: c,
      informative: n,
      absent: a,
      rate,
      verdict: expected ? 'expected_absent'
        : rate >= explainable * ABSENCE_MARGIN ? 'absent'
          : rate <= explainable ? 'present' : 'unclear',
      note: c === 'X:PAR' ? 'pseudoautosomal: on both the X and the Y, so a sperm of either type '
        + 'delivers it. Present here is the positive control on the rest of chrX being absent'
        : expected ? "no paternal X: this sample carries the parent's Y"
        : sex && role === 'paternal' && present && xLost
          ? 'the paternal X is absent and no chrY was called, so this is not ordinary sex '
            + 'determination. Reported as a loss.'
          : undefined,
    })
  }

  // Dispersion across chromosomes: uniform, or confined to part of the genome. Reported always,
  // consulted only where the genome-wide rate is already ambiguous.
  const rates = chroms.filter((c) => isAutosome(c.chrom)).map((c) => c.rate)
  const mean = rates.length ? rates.reduce((a, b) => a + b, 0) / rates.length : NaN
  const sd = rates.length > 1
    ? Math.sqrt(rates.reduce((a, b) => a + (b - mean) ** 2, 0) / (rates.length - 1)) : NaN
  const dispersion = Number.isFinite(sd) && mean > 0 ? sd / mean : NaN
  const minChrom = rates.length ? Math.min(...rates) : NaN

  // Neither present nor absent, and differing by the same amount on EVERY chromosome, is two
  // genomes blended rather than one genome missing pieces. A lost segment is patchy and leaves
  // the rest clean; a blend cannot be, since every chromosome carries both contributions.
  if (verdict === 'unclear' && Number.isFinite(dispersion) && dispersion < UNIFORM_CV
      && minChrom > explainable) {
    limits.push(
      `Absence is uniform across the genome: coefficient of variation ${dispersion.toFixed(2)} `
      + `between chromosomes, and the cleanest chromosome still sits at ${pct(minChrom)}, above `
      + `this sample's own ceiling of ${pct(explainable)}. A lost segment would leave the rest of `
      + 'the genome clean and does not look like this; two genomes blended in one tube do, '
      + 'because every chromosome then carries both contributions. Read this as a mixed sample '
      + 'rather than a partial loss. No reanalysis separates them.',
    )
  }

  if (originClass === 'androgenetic') {
    notes.push('Every allele traces to this parent and the genome is homozygous: a parent-only '
      + 'complement, duplicated. Nothing here required chrY.')
  }
  if (originClass === 'gynogenetic' && spermType === 'unknown') {
    notes.push('No contribution from this parent anywhere. chrY alone would not have separated '
      + 'this from an X-bearing sperm carrying a full paternal genome.')
  }
  if (spermType === 'X_bearing') {
    notes.push("No chrY, yet the parent's alleles are present on chrX as well as the autosomes: "
      + 'an X-bearing sperm. This is the case a chrY test cannot call.')
  }

  return {
    verdict, originClass, zygosity, spermType, genomeRate, explainable, informative: nTot,
    nonParentalRate, secondParentExpected, hetBand, noCallRate, hetFraction,
    dispersion, minChromRate: minChrom, chroms, notes, limits,
  }
}

const byChromName = (a: [string, unknown], b: [string, unknown]): number =>
  a[0].length - b[0].length || a[0].localeCompare(b[0])

/* --- both parents ---------------------------------------------------------------------------- */

export type PairClass = OriginClass | 'neither_parent'

export interface PairResult {
  originClass: PairClass
  /** Fraction of shared called markers where the two declared parents agree. NaN if too few. */
  agreement: number
  notes: string[]
}

/**
 * How often the two declared parents carry the same genotype, over markers both called.
 *
 * Only ever used to catch a labelling accident. Two arrays of one person agree at nearly every
 * marker, and that produces a confident "biparental" because both tests pass against the same
 * genome. The 1,000-marker floor is there because a handful of shared markers says nothing.
 */
export function agreement(a: Map<string, AB>, b: Map<string, AB>): number {
  let shared = 0
  let same = 0
  // Walk the smaller map: the intersection is the same either way and this is the cheaper scan.
  const [small, large] = a.size <= b.size ? [a, b] : [b, a]
  for (const [k, v] of small) {
    const w = large.get(k)
    if (w === undefined || v === 'NC' || w === 'NC') continue
    shared += 1
    if (v === w) same += 1
  }
  return shared > 1000 ? same / shared : NaN
}

/**
 * The class of a sample tested against both parents, mirroring `origin.both_parents`.
 *
 * With only the father, a maternal origin is inferred from his absence. With her array it is
 * measured, which is what separates "no paternal contribution" from "belongs to neither of these
 * two people" - a distinction a single-parent run cannot make, and one that a mislabelled file
 * or the wrong donor pair produces just as readily as biology does.
 */
export function pair(
  paternal: ParentageResult, maternal: ParentageResult, parentAgreement: number,
): PairResult {
  const P = paternal.verdict === 'parent_genome_present'
  const M = maternal.verdict === 'parent_genome_present'
  const unresolved = paternal.verdict === 'unclear' || maternal.verdict === 'unclear'
  const notes: string[] = []

  const originClass: PairClass = unresolved ? 'unclear'
    : P && M ? 'biparental'
      : P ? 'androgenetic'
        : M ? 'gynogenetic' : 'neither_parent'

  // 0.90, not 0.99, matching `ingest.relatedness`: real replicate arrays of one person concord
  // at 95.8% and a parent-offspring pair at 54.9%, so 0.99 sat above the thing it was meant to
  // catch and missed every genuine duplicate.
  if (parentAgreement > 0.90) {
    notes.push(`The two declared parents agree at ${pct(parentAgreement, 1)} of shared markers, `
      + 'which is a duplicate or a relabelled file rather than two people. Every conclusion here '
      + 'rests on them being different, so resolve this before reading it.')
  }
  if (originClass === 'neither_parent') {
    notes.push('Neither declared parent accounts for this genome. Before reading that as '
      + 'biology, check the pairing: a mislabelled sample or the wrong donor pair produces '
      + 'exactly this, and so does contamination.')
  }
  if (originClass === 'androgenetic') {
    notes.push('The maternal complement is absent and the oocyte donor\'s own array confirms it '
      + 'directly, rather than it being inferred from allele sharing.')
  }
  if (originClass === 'gynogenetic') {
    notes.push('The oocyte donor accounts for this genome and the sperm donor does not. Measured '
      + 'against both arrays, so this is not an absence standing in for a presence.')
  }
  if (originClass === 'unclear') {
    // Name the half that IS settled. The class has to stay unclear, but saying only "at least one
    // parent is unresolved" throws away a contribution that was measured with margin.
    const resolved = (r: ParentageResult, who: string): string =>
      `the ${who}'s contribution is confirmed ${r.verdict === 'parent_genome_present'
        ? 'present' : 'absent'}`
    const both = paternal.verdict === 'unclear' && maternal.verdict === 'unclear'
    notes.push(both
      ? 'Neither parent is resolved, so the pair cannot be classified. The per-parent rates say '
        + 'how far each sits from its own ceiling.'
      : `${paternal.verdict === 'unclear'
        ? `${resolved(maternal, 'oocyte donor')}, but the sperm donor's is unresolved`
        : `${resolved(paternal, 'sperm donor')}, but the oocyte donor's is unresolved`}`
        + ', so the pair cannot be classified. What is open is only whether the other parent '
        + 'contributed; the settled half stands.')
  }
  return { originClass, agreement: parentAgreement, notes }
}

/** Percentages are formatted here so the figure, the table and the prose cannot disagree. */
export const pct = (x: number, dp = 2): string =>
  Number.isFinite(x) ? `${(x * 100).toFixed(dp)}%` : 'n/a'

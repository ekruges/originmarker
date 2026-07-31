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
 * The arithmetic is duplicated from the Python deliberately and bounded the same way the
 * run-length statistic is: both sides are pinned by the same fixture, so a divergence fails a
 * test rather than reaching a result.
 */
import { type ProbeRow } from './ingest.ts'
import { type AB } from './informativity.ts'

/** Residual absence on clean data, from genotyping error alone. Measured at 0.03% and 0.05%. */
export const ABSENCE_ERROR_FLOOR = 0.005

/** How far above the explainable level an absence must sit before it is called. */
export const ABSENCE_MARGIN = 3.0

/** Diploid genomes run 15-16% in the heterozygous BAF band, uniparental ones 1.3-3.4%. */
export const HET_BAND_DIPLOID = 0.08

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
}

export const emptyTally = (): Tally => ({
  byChrom: new Map(), nonParental: 0, nonParentalDen: 0,
  called: 0, het: 0, markers: 0, bafInBand: 0, bafTotal: 0,
})

/** One marker of the sample, against the parent's call at the same marker. */
export function tallyRow(parent: AB, row: ProbeRow, t: Tally): void {
  t.markers += 1
  // The BAF band is counted BEFORE the no-call check, deliberately and to match the Python. A
  // marker whose genotype failed still has an intensity reading, and a dropped heterozygote sits
  // in the middle of the band, which is exactly the evidence that distinguishes a genuinely
  // homozygous genome from a diploid one that lost calls. Excluding them shifted the band from
  // 1.27% to 0.92% on a real sample and moved the noise ceiling with it.
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
  const d = t.byChrom.get(row.chrom) ?? [0, 0]
  d[0] += 1
  if (!row.genotype.includes(parent[0])) d[1] += 1
  t.byChrom.set(row.chrom, d)
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
  chroms: ChromResult[]
  notes: string[]
  limits: string[]
}

/** Turn the tallies into the three axes and the class they imply. */
export function classify(
  t: Tally,
  parentHeterozygosity: number,
  opts: { role?: 'paternal' | 'maternal' } = {},
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
  const explainable = absenceExplainable(noCallRate, hetFraction) + ABSENCE_ERROR_FLOOR
  const nonParentalRate = t.nonParentalDen ? t.nonParental / t.nonParentalDen : NaN
  const secondParentExpected = secondParentSignal(parentHeterozygosity)

  let verdict: Verdict
  if (!nTot) {
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
  if (Number.isFinite(hetBand)) {
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
  for (const [c, [n, a]] of [...t.byChrom].sort(byChromName)) {
    if (n < 200) continue
    const rate = a / n
    const sex = c === 'X' || c === '23'
    // A male offspring has no paternal X at all: his father sent a Y instead. Flagging that as a
    // loss would report ordinary sex determination as an anomaly on half of all samples. A
    // MOTHER transmits an X either way, so the same exemption must not apply to her.
    const expected = sex && role === 'paternal' && present && rate > explainable * ABSENCE_MARGIN
    if (sex && present) spermType = rate > explainable * ABSENCE_MARGIN ? 'Y_bearing' : 'X_bearing'
    chroms.push({
      chrom: c,
      informative: n,
      absent: a,
      rate,
      verdict: expected ? 'expected_absent'
        : rate >= explainable * ABSENCE_MARGIN ? 'absent'
          : rate <= explainable ? 'present' : 'unclear',
      note: expected ? "no paternal X: this sample carries the parent's Y" : undefined,
    })
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
    nonParentalRate, secondParentExpected, hetBand, noCallRate, chroms, notes, limits,
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

  if (parentAgreement > 0.99) {
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
    notes.push('At least one parent is unresolved, so the pair cannot be classified. The '
      + 'per-parent rates say which one and by how much.')
  }
  return { originClass, agreement: parentAgreement, notes }
}

/** Percentages are formatted here so the figure, the table and the prose cannot disagree. */
export const pct = (x: number, dp = 2): string =>
  Number.isFinite(x) ? `${(x * 100).toFixed(dp)}%` : 'n/a'

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
import type { Segment } from './segments.ts'
import type { HetCall } from './obligateHet.ts'
import type { Enrichment } from './features.ts'
import type { StageCall } from './stage.ts'

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
 *
 * TWO LIMITS ON WHAT THIS MAY BE USED FOR, both measured rather than assumed.
 *
 * It does NOT detect MOSAICISM, and must never be extended to. A blend of two whole genomes and a
 * mosaic are different objects: the blend contributes on every chromosome at a full fraction,
 * while a mosaic shifts one chromosome by the mosaic fraction. Measured against uniform artefact,
 * the separation ratio of this statistic is 0.58 at f=0.10, 0.69 at 0.20, 0.83 at 0.30, 0.99 at
 * 0.50 and 1.00 at 0.80. A ratio at or below 1 is no separation at any fraction. The 8.36 that
 * appears at f=1.00 is a whole extra genome, which is the blend case above and not a mosaic.
 *
 * The 1.106 figure is a same-individual REPLICATE, which shares the whole genome and is the
 * easiest control in the set. It is reported for completeness and carries no weight: the lower
 * edge of the gap this threshold sits in is the 0.762 degraded true pair, which does not.
 */
export const UNIFORM_CV = 0.35

/**
 * Fraction of a chromosome's B-allele frequencies that must sit at an extreme, outside 0.15 to
 * 0.85, before any verdict is reported for that chromosome.
 *
 * This catches array MIS-CLUSTERING, which no allelic statistic can see. When a chromosome's
 * probes cluster badly the genotype calls come back systematically wrong, absence rises to
 * something that looks exactly like a real chromosome-scale loss, and every measure built on
 * genotypes agrees with it because they are all reading the same broken calls. Only the
 * intensity-derived allelic ratio shows the cause.
 *
 * Measured across the five paternal pronuclei of GSE148488, every one a meiotic product of the
 * sperm donor and therefore present on every autosome:
 *
 *     real autosomes, all five samples, 110 chromosomes   0.752 to 0.976
 *     GSM4774681 chr1, reported ABSENT at 18.36%          0.130
 *
 * That sample's genome-wide verdict is "parent genome present" and chr1 is his son's. The gap
 * between 0.130 and 0.752 is empty by a factor of 5.8, which is why a single figure sits in it.
 *
 * A genuinely trisomic chromosome also fails this gate, since its allelic ratios cluster at a
 * third and two thirds rather than at the extremes. That is the correct outcome and not a cost:
 * gains are refused on this platform in both channels, so the alternative to withholding is a
 * wrong answer rather than a right one.
 */
export const BAF_EXTREME_FLOOR = 0.40

/**
 * How far a chromosome's allelic ratio has to sit from the rest of the genome, in standard
 * deviations, before a mosaic is reported.
 *
 * The GENOTYPE cannot see a mosaic at all below a fraction of about 0.91. A genotype call is a
 * threshold on the allelic ratio, the AB-against-BB boundary sits near 0.917, and a chromosome
 * present at fraction f has a true ratio of 1/(2-f), which does not cross that boundary until
 * f = 0.909. Below it the call does not change and the marker carries nothing, so adding markers
 * buys nothing against a structural floor. The ratio ITSELF is not thresholded and moves from the
 * first percent, which is why this reads the ratio instead.
 *
 * The contrast is INTERNAL, against the sample's own other chromosomes, and that is load-bearing
 * rather than tidy. A degraded array has a globally shifted ratio, so any absolute deviation
 * measure calls the entire genome mosaic; a contrast cannot, because the shift sits in both terms.
 *
 * Measured on four bulk diploid arrays with no mosaic anywhere, 88 chromosome observations: z runs
 * -1.65 to 5.18. Titrating a mosaic loss onto one chromosome of those same arrays, shifting the
 * observed ratio and keeping each array's own noise:
 *
 *     fraction   z across the four arrays        detected at this threshold
 *       0.15      3.4   1.1   7.9   4.0            0 of 4
 *       0.20      5.2   2.1  12.2   6.3            1 of 4
 *       0.30      9.7   4.6  22.7  12.5            3 of 4
 *       0.50     21.3  11.7  48.5  28.6            4 of 4
 *
 * So the honest floor is 0.50 for a reliable call and 0.30 for a partial one, not the 0.20 an
 * earlier estimate suggested. Array quality dominates: the same fraction reads 33.8 on one array
 * and 128.4 on another. NO FRACTION IS REPORTED. Inverting the statistic to a fraction is biased
 * low by roughly half, because the truncation that hides the mosaic from the genotype also removes
 * the most-shifted markers from the mean.
 */
export const MOSAIC_Z = 8

/**
 * Fraction of the genome's median call rate below which a chromosome is not merely noisy: it is
 * not there to be genotyped.
 *
 * A chromosome that has been lost yields no DNA, so the array reads nothing at its probes and the
 * call rate collapses. That is a GENOTYPE-level signal and needs no intensity, which matters
 * because intensity on amplified material cannot be trusted for fine copy-number work.
 *
 * Measured over 1,012 chromosome observations from 46 arrays across three experiments:
 *
 *     eleven chromosomes      0.20x to 0.41x the genome median
 *     the other 1,001         0.78x to 1.16x
 *
 * The gap is empty by a factor of 1.9. Every one of the eleven also carries a large intensity
 * shift, |log2R| of 1.59 to 2.04 against a spread of -0.79 to +0.42 on the rest, from two
 * channels that fail in unrelated ways.
 *
 * THIS CORRECTS AN EARLIER ERROR. Three of those chromosomes were previously diagnosed as array
 * mis-clustering and withheld, on the reasoning that a paternal pronucleus is its father's on
 * every autosome by construction. That reasoning was wrong: the series exists because chromosomes
 * ARE lost in these embryos, so a pronucleus can genuinely lack one. Mis-clustering and loss both
 * depress the allelic ratio, and the call rate is what separates them, because a mis-clustered
 * chromosome still calls.
 */
export const CALL_COLLAPSE = 0.60

/**
 * Which way a collapsed chromosome went. Only the SIGN is used, and only once the call rate has
 * already established that something is wrong.
 *
 * Measured on the same eleven: six sit at log2R -1.59 to -2.04 against the genome, and five at
 * +1.60 to +1.95. Nothing else in 1,012 observations leaves -0.79 to +0.42. So the sign is not a
 * threshold on a noisy quantity, it is the direction of a shift already known to be real.
 *
 * A gain is reported more cautiously than a loss and says so: fine copy-number work on amplified
 * material is refused elsewhere in this tool for good measured reasons, and what makes this
 * different is only that the effect is an order of magnitude larger than the noise those
 * measurements were about.
 */
export const LRR_SHIFT = 1.0

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
  /** Per chromosome: [B-allele frequencies at an extreme, B-allele frequencies read]. The
   *  mis-clustering check; see BAF_EXTREME_FLOOR. */
  bafByChrom: Map<string, [number, number]>
  /** Per chromosome: [sum of |BAF - 0.5| over heterozygous calls, count]. The mosaic contrast;
   *  see MOSAIC_Z. Heterozygous calls only, because a homozygous site sits at an extreme
   *  already and its deviation says nothing about a mixture. */
  hetDevByChrom: Map<string, [number, number]>
  /** Per chromosome: [markers called, markers on the array]. A chromosome that is not there
   *  cannot be genotyped, so this is what an aneuploidy shows up in; see CALL_COLLAPSE. */
  callByChrom: Map<string, [number, number]>
  /** Per chromosome, every log2 intensity ratio read. Only the SIGN of the median is used, to
   *  tell a loss from a gain once the call rate has already said something is wrong. */
  lrrByChrom: Map<string, number[]>
  /** chrY markers the sample called, and chrY markers the parent's file carries. A Y-bearing
   *  sperm is measured from these, never inferred from an absent X. */
  yCalled: number
  yTotal: number
  /** Assembly, so the pseudoautosomal boundaries used are the right ones. GRCh37 if unknown. */
  build: string | null
}

export const emptyTally = (): Tally => ({
  byChrom: new Map(), nonParental: 0, nonParentalDen: 0,
  called: 0, het: 0, markers: 0, bafInBand: 0, bafTotal: 0, bafByChrom: new Map(),
  hetDevByChrom: new Map(), callByChrom: new Map(), lrrByChrom: new Map(),
  yCalled: 0, yTotal: 0, build: null,
})

/** One marker of the sample, against the parent's call at the same marker. */
export function tallyRow(parent: AB, row: ProbeRow, t: Tally): void {
  t.markers += 1
  if (isAutosome(row.chrom)) {
    const cc = t.callByChrom.get(row.chrom) ?? [0, 0]
    cc[1] += 1
    if (row.genotype !== 'NC') cc[0] += 1
    t.callByChrom.set(row.chrom, cc)
    if (row.log2R !== null) {
      (t.lrrByChrom.get(row.chrom) ?? t.lrrByChrom.set(row.chrom, []).get(row.chrom)!)
        .push(row.log2R)
    }
  }
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
    const b = t.bafByChrom.get(row.chrom) ?? [0, 0]
    b[1] += 1
    if (row.baf < 0.15 || row.baf > 0.85) b[0] += 1
    t.bafByChrom.set(row.chrom, b)
    if (row.genotype === 'AB') {
      const d = t.hetDevByChrom.get(row.chrom) ?? [0, 0]
      d[0] += Math.abs(row.baf - 0.5)
      d[1] += 1
      t.hetDevByChrom.set(row.chrom, d)
    }
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
    /** The chromosome's genotype calls are not measuring it. See BAF_EXTREME_FLOOR. */
    | 'not_measured'
  /** Fraction of this chromosome's B-allele frequencies at an extreme. The mis-clustering check. */
  bafExtreme: number
  /** How far this chromosome's heterozygous allelic ratio sits from the sample's own others, in
   *  standard deviations. NaN unless the sample is diploid; see MOSAIC_Z. */
  mosaicZ: number
  /** Whole-chromosome copy change, from the call rate collapsing; see CALL_COLLAPSE. */
  aneuploidy?: 'loss' | 'gain'
  /**
   * Whether this chromosome carries one parental contribution or two, from heterozygosity at
   * markers where a parent is homozygous. See obligateHet.ts.
   *
   * PER CHROMOSOME rather than per array, because a cell can be biparental overall and
   * uniparental on part of its genome, which is the event a whole-array verdict would hide.
   * Absent when no parent made enough markers informative here.
   */
  contribution?: HetCall
  /** Whose copy changed, where determinable. 'this' is the parent this result is scored against;
   *  'other' means their alleles survive on what is left, so the copy that went was the other
   *  parent's. Undefined where the sample carries none of this parent's genome anywhere. */
  aneuploidyParent?: 'this' | 'other'
  /** This chromosome's call rate as a fraction of the genome's median. */
  callFraction: number
  /** Median log2 intensity ratio against the genome's. Its SIGN tells a loss from a gain. */
  lrrShift: number
  note?: string
}

/** A gain, annotated with where the extra copy came from or why that cannot be said. */
export interface GainAnnotation {
  where: string
  kind: 'whole chromosome' | 'segment'
  /** Which parent, when the cell carries both and the region has the markers for it. */
  origin: string
  /** The evidence, or the reason there is none. */
  why: string
  /** True when a direction was established. */
  called: boolean
}

export interface ParentageResult {
  /** Which parent this result is ABOUT. Carried rather than assumed: a run with one parent
   *  loaded can be either, and a display that assumes paternal names the wrong parent. */
  role: 'paternal' | 'maternal'
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
  /** Where along a chromosome the parental genome is missing, rather than whether. Populated by
   *  the caller, which is the only place marker positions are in hand; see `segments.ts`. */
  segments: Segment[]
  /** One entry per gain found, whole-chromosome or segmental. Empty when there are no gains,
   *  which is not the same as gains with no origin: those appear here saying so. */
  gains: GainAnnotation[]
  /**
   * Which parent's copy is MISSING, for every loss called. Separate from `gains` because the two
   * read the same measurement in opposite directions: an under-represented parent is the one lost,
   * and the one that GAINED when the event is a gain. See callLossOrigin.
   */
  losses: GainAnnotation[]
  /**
   * Where the called regions sit against fragile sites, long genes and centromeres. Positional
   * only: it needs no parental genotype, and says nothing about parent of origin.
   */
  placement?: Enrichment[]
  /**
   * Per-segment call from the embryo's own event-free cells, where the run holds siblings.
   * Says whether a copy is genuinely missing rather than which parent's: naming a side needs
   * phase, and naming the parent needs an anchor. See siblingOrigin.ts.
   */
  /**
   * Per-segment origin from a SINGLE loaded parent. 'known-parent-lost' means the loaded parent's
   * copy is missing; with the sperm donor loaded that is paternal. See oneParentOrigin.ts.
   */
  /** Developmental stage inferred from this array, and the dropout its template count implies. */
  stage?: StageCall
  oneParent?: {
    where: string
    verdict: string
    posterior: number
    markers: number
    exclusive: number
    why: string
  }[]
  /**
   * Whole-chromosome origin, read from ALLELE DOSAGE rather than genotypes.
   *
   * Separate from `oneParent` because the two answer the same question from different channels and
   * are not interchangeable. A whole chromosome is detected by the collapse of its genotype call
   * rate, so the genotype channel has no evidence left on exactly those events; dosage is read
   * whether or not a genotype is emitted. Keyed by `chrN` rather than by an interval, since the
   * unit here is the whole chromosome.
   */
  dosageCalls?: {
    where: string
    verdict: string
    /** Self-referenced centroid shift. Positive means the loaded parent's copy is short. */
    shift: number
    z: number
    /** Mosaic fraction the shift implies, via f = 4d/(1+2d). */
    impliedF: number
    /** Markers in the central window, the only denominator that carries the signal. */
    window: number
    markers: number
    material: string
    /** Smallest fraction this material and width could detect. NaN where none could. */
    floor: number
    why: string
  }[]
  siblingCalls?: {
    where: string
    hypothesis: string
    posterior: number
    markers: number
    phi: number
    why: string
  }[]
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

  // The mosaic contrast, per chromosome against the sample's own others. Only meaningful on a
  // DIPLOID sample: a uniparental genome is homozygous by construction, so it has no heterozygous
  // sites for a mixture to shift and any deviation there is artefact.
  const hetDev = new Map<string, number>()
  for (const [c, [sum, n]] of t.hetDevByChrom) {
    if (isAutosome(c) && n >= 500) hetDev.set(c, sum / n)
  }
  const mosaicZ = new Map<string, number>()
  if (zygosity === 'diploid' && hetDev.size >= 4) {
    for (const [c, v] of hetDev) {
      const others = [...hetDev].filter(([x]) => x !== c).map(([, y]) => y)
      const mu = others.reduce((a, b) => a + b, 0) / others.length
      const sd = Math.sqrt(others.reduce((a, b) => a + (b - mu) ** 2, 0) / (others.length - 1))
      if (sd > 0) mosaicZ.set(c, (v - mu) / sd)
    }
  }

  // Aneuploidy, from the call rate rather than from the alleles. A chromosome that is gone yields
  // no DNA and cannot be genotyped, so it collapses here while its allelic statistics only look
  // noisy. Computed before the per-chromosome verdicts because it changes what they may say.
  const callPerChrom = new Map<string, number>()
  for (const [c, [k, n]] of t.callByChrom) if (isAutosome(c) && n >= 200) callPerChrom.set(c, k / n)
  const callSorted = [...callPerChrom.values()].sort((x, y) => x - y)
  const callMedian = callSorted.length ? callSorted[callSorted.length >> 1] : NaN
  const medianOf = (xs: number[]): number => {
    if (!xs.length) return NaN
    const q = [...xs].sort((x, y) => x - y)
    return q[q.length >> 1]
  }
  const genomeLrr = medianOf([...t.lrrByChrom.values()].flat())
  const aneuploidy = new Map<string, 'loss' | 'gain'>()
  const callFrac = new Map<string, number>()
  const lrrShift = new Map<string, number>()
  for (const [c, r] of callPerChrom) {
    const frac = callMedian > 0 ? r / callMedian : NaN
    callFrac.set(c, frac)
    const shift = medianOf(t.lrrByChrom.get(c) ?? []) - genomeLrr
    lrrShift.set(c, shift)
    if (!(frac < CALL_COLLAPSE)) continue
    // Which way it went. Without a usable intensity the collapse is still reported, as a loss,
    // because that is what a chromosome absent from the array looks like from the genotypes alone.
    aneuploidy.set(c, Number.isFinite(shift) && shift >= LRR_SHIFT ? 'gain' : 'loss')
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
    // Mis-clustering first, because every genotype-derived measure below reads the same broken
    // calls and would agree with each other about a chromosome that was never measured.
    const [ext, nBaf] = t.bafByChrom.get(c) ?? t.bafByChrom.get(c.replace(':PAR', '')) ?? [0, 0]
    const extreme = nBaf ? ext / nBaf : NaN
    // A collapsed call rate says the chromosome is not there, which ALSO depresses the allelic
    // ratio. Both look identical to the ratio alone, so the ratio may only diagnose mis-clustering
    // where the chromosome is still being genotyped. Withholding a real loss was the earlier bug.
    const aneu = aneuploidy.get(c)
    const clustered = !Number.isFinite(extreme) || extreme >= BAF_EXTREME_FLOOR || aneu !== undefined

    chroms.push({
      chrom: c,
      informative: n,
      absent: a,
      bafExtreme: extreme,
      mosaicZ: mosaicZ.get(c) ?? NaN,
      aneuploidy: aneu,
      // Whose copy went. Determinable only where the parent is measurable on this chromosome at
      // all: a sample carrying none of this parent's genome anywhere cannot say anything per
      // chromosome, and a chromosome with nothing left has no surviving copy to attribute.
      aneuploidyParent: !aneu || !present ? undefined
        : rate <= explainable ? 'other'
          : rate >= explainable * ABSENCE_MARGIN ? 'this' : undefined,
      callFraction: callFrac.get(c) ?? NaN,
      lrrShift: lrrShift.get(c) ?? NaN,
      verdict: !clustered ? 'not_measured'
        : aneu === 'loss' ? 'absent'
          : expected ? 'expected_absent'
            : rate >= explainable * ABSENCE_MARGIN ? 'absent'
              : rate <= explainable ? 'present' : 'unclear',
      rate,
      note: aneu
        ? `${aneu === 'loss' ? 'LOST' : 'GAINED'}: this chromosome calls at `
          + `${pct(callFrac.get(c) ?? NaN, 0)} of the genome's median rate and its intensity sits `
          + `${(lrrShift.get(c) ?? NaN).toFixed(2)} log2 from the rest. A chromosome that is not `
          + 'there cannot be genotyped, which is what the call rate is reading; the intensity says '
          + `which way it went. Measured over 1,012 chromosomes, an intact one calls at 0.78x to `
          + '1.16x and never leaves -0.79 to +0.42 log2.'
          + (!present
            ? ` No parent is attached: this sample carries none of the ${role} genome anywhere, `
              + 'so there is nothing to attribute per chromosome.'
            : rate <= explainable
              ? ` The ${role} alleles ARE present on what remains, so the copy that went was the `
                + "other parent's."
              : rate >= explainable * ABSENCE_MARGIN
                ? ` The ${role} alleles are absent from what remains, so the copy that went was `
                  + `the ${role} one.`
                : ' Which parent lost the copy is not resolved: absence here sits between what '
                  + "this sample's own noise explains and the margin needed to call it.")
        : !clustered
        ? `${pct(extreme, 1)} of this chromosome's B-allele frequencies sit at an extreme, `
          + `against a ${pct(BAF_EXTREME_FLOOR, 0)} floor and ${pct(0.752, 0)} to `
          + `${pct(0.976, 0)} on chromosomes that clustered correctly. The genotype calls here `
          + 'are not measuring this chromosome, so no verdict is reported for it. Absence would '
          + 'read as a chromosome-scale loss and every genotype-derived measure would agree, '
          + 'because they all read the same broken calls.'
        : c === 'X:PAR' ? 'pseudoautosomal: on both the X and the Y, so a sperm of either type '
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
  const mosaicHits = chroms.filter((c) => c.mosaicZ >= MOSAIC_Z)
  if (mosaicHits.length) {
    limits.push(
      `The allelic ratio on ${mosaicHits.map((c) => `chr${c.chrom}`).join(', ')} sits `
      + `${mosaicHits.map((c) => c.mosaicZ.toFixed(1)).join(', ')} standard deviations from the `
      + "rest of this genome, measured over heterozygous sites against the sample's own other "
      + 'chromosomes. That is the signature of a MIXTURE OF CELL LINEAGES carrying different copy '
      + 'numbers there, and the genotype calls cannot show it: a call is a threshold on that '
      + 'ratio, and the threshold is not crossed until roughly 90% of cells are affected. No '
      + 'fraction is reported, because inverting this statistic to one is biased low by about '
      + 'half. Reliable detection begins near 50% of cells and is partial from 30%, so this is a '
      + 'floor on what was seen rather than a measure of how much.',
    )
  }

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
    role,
    verdict, originClass, zygosity, spermType, genomeRate, explainable, informative: nTot,
    segments: [],
    gains: [],
    losses: [],
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

/**
 * What kind of material this array is, inferred from the array itself.
 *
 * READ THE LIMITS BEFORE THE LADDER. An external audit (audit/STAGE-AUDIT.txt) found the first
 * version of this module overstated its mechanism and shipped an estimator that is biased in ways
 * invisible in its output. What follows is the corrected account.
 *
 * WHY DROPOUT HAPPENS, which is the part that holds. A heterozygous locus in a single cell is
 * represented by one molecule per allele before amplification. If that molecule fails to prime in
 * the early cycles the allele is absent from everything downstream and the marker reads homozygous.
 * Nothing about the genome changed; the reaction lost it.
 *
 * WHY THE LADDER IS NOT A COPY-NUMBER MODEL, which is where the first version was wrong. If each
 * template failed independently with probability f, dropout would be f^n for n starting molecules.
 * Calibrating f on the single-cell rung predicts 5.5e-6 for a 15-template biopsy against a measured
 * 0.050, four orders of magnitude out, and exactly 0 for bulk against a measured 0.013. The bulk and
 * biopsy figures are therefore NOT template loss: losing every one of a million copies is
 * impossible, and 0.013 sits inside this platform's own replicate error, measured at 3.31% marker
 * disagreement between technical replicates of one bulk sample, 98.4% of it heterozygote-to-
 * homozygote. The values below are CALIBRATION CONSTANTS measured on this material, not predictions
 * of a mechanism, and the low rungs are dominated by genotyping error rather than dropout.
 *
 * WHY SINGLE ES CELLS AND BLASTOMERES DIFFER IS UNKNOWN. Both begin from two molecules. A published
 * experiment varying reaction conditions across >3000 single-cell amplifications found amplicon
 * size, DNA degradation, freeze-thaw and cell number affected dropout, while CELL TYPE had little or
 * no effect (Piyamongkol 2003). An earlier version of this file attributed the difference to
 * chromatin state; that is the mechanism that study looked for and did not find, and it has been
 * removed. A same-laboratory comparison also reports first polar bodies, with one genome, showing
 * the LOWEST dropout of three cell types, so template count does not order these data either
 * (Rechitsky 1998). The parameter is empirical.
 *
 * HETEROZYGOSITY SHORTFALL IS A CORRELATE, NOT A MEASUREMENT. Consanguinity, copy-neutral loss of
 * heterozygosity, uniparental disomy and ancestry differing from the anchor population all depress
 * heterozygosity and are absorbed into the estimate: first-cousin consanguinity biases it +0.050,
 * 20% genome LOH +0.160, and an East Asian sample against a European anchor +0.155, which is enough
 * to classify bulk DNA as a single cell with nothing in the output indicating a problem. Where two
 * amplifications of one genome exist, `dropoutFromReplicates` needs no anchor and should be
 * preferred.
 *
 * HAPLOID IS A DIFFERENT AXIS, AND FIRST POLAR BODIES ARE NOT HOMOZYGOUS. A PB2, pronucleus or
 * sperm carries one chromatid and has true heterozygosity of zero. A FIRST polar body carries a
 * dyad, two sister chromatids of one homologue, and distal to every crossover those sisters carry
 * different haplotypes: roughly 44% of a PB1 genome is genuinely heterozygous, expected h about
 * 0.074 with no dropout at all. Treating its heterozygous calls as error mis-parameterises it.
 *
 * QUALITY IS GATED FIRST, ALWAYS. A failed amplification drives heterozygosity toward zero, which
 * without a gate reads as haploid and receives a confident parameter set instead of a rejection.
 *
 * THE BAF BAND VETOES HAPLOID BUT CANNOT CONFIRM IT. The audit asked for B-allele-frequency band
 * position as a required second signal, on the reasoning that intensity sees a heterozygous locus
 * whether or not the genotype caller dropped it. That holds one way only, and this project has the
 * measurement: across 120 non-haploid arrays spanning blastomeres, trophectoderm biopsies and ESC
 * lines, 91 fall INSIDE the haploid band range and 22 blastomeres sit below the diploid-exclusion
 * floor, because whole-genome amplification widens the heterozygote band and fills the homozygote
 * band until the two are one distribution. So a low band is no evidence of one genome on amplified
 * material and is not treated as any. A HIGH band is still decisive: intensity cannot manufacture
 * intermediate signal from a single template, so a sample above the diploid line carries two
 * genomes no matter how few heterozygous calls survived, and that is the one path by which a
 * heavily dropped-out diploid was still being called haploid after the quality gate.
 */
import type { SampleProfile } from './ingest.ts'
import { NO_INFORMATION_MIN } from './informationRatio.ts'
import { BULK_HETEROZYGOSITY, DROP_IN } from './oneParentOrigin.ts'
import { NORMAL_MAJORITY_MIN, type CopyState } from './oneParentState.ts'

export type Stage =
  | 'bulk'
  | 'trophectoderm'
  | 'single-cell'
  | 'blastomere'
  | 'haploid'
  | 'failed'
  | 'unknown'

/**
 * How many genomes the array carries, which is a different question from what kind of material it
 * is. `stage` answers the second: blastomere, trophectoderm, bulk. Nothing answered the first, so a
 * diploid cell was only ever described by the material word.
 */
export type Ploidy = 'haploid' | 'diploid' | 'not one genome' | 'unknown'

export interface StageCall {
  stage: Stage
  /** How many genomes, separate from the material `stage` names. */
  ploidy: Ploidy
  /**
   * Which test refused the array, on a `failed` stage only. Carried so a refusal record can name
   * the gate and its numbers without reading them back out of `why`.
   */
  refusedBy?: 'heterozygosity-ceiling' | 'declared' | 'no-genotype-signal' | 'second-genome'
  /**
   * The drop-in, measured on this array against a founder panel, that was removed before it was
   * staged. Present only on an array that read over the heterozygosity ceiling and came back under
   * it with drop-in removed: one genome, whose excess heterozygosity is calls that cannot be real.
   */
  dropInCorrected?: number
  /** Expected allele dropout, used to parameterise downstream likelihoods. */
  dropout: number
  /** How the dropout figure was arrived at, which decides how much it can be trusted. */
  basis: 'replicate-discordance' | 'heterozygosity-shortfall' | 'stage-default' | 'none'
  /** Order-of-magnitude template molecules per locus. Reported, NOT used to predict dropout. */
  templates: string
  markerFloor: number
  /** Confounds that would be absorbed into a shortfall estimate, stated rather than hidden. */
  caveat: string
  why: string
}

/** Defined once beside the dropout it is read against: oneParentOrigin.BULK_HETEROZYGOSITY. */
export { BULK_HETEROZYGOSITY }

/**
 * Chromosomes that must state a copy state against the loaded parent before the alleles may say how
 * many genomes an array carries: most of the 22.
 */
export const ALLELE_PLOIDY_MIN_STATED = 11

/**
 * How many genomes the alleles show against a loaded parent, measured or reconstructed.
 *
 * A haploid product carries one parent's genome on every chromosome: one copy of that parent's
 * (`other-lost` against it) or of the other's (`parent-lost`). A diploid carries a second genome on
 * most chromosomes, and allele dropout does not hide it, since dropout removes alleles at random and
 * leaves the other parent's in the calls. So two genomes when most stated chromosomes read disomic,
 * and one when the commonest state is one copy of one parent on NORMAL_MAJORITY_MIN of them. A chaotic
 * genome holding one copy of some chromosomes from each parent is neither. No template and undecided
 * state nothing and are left out.
 *
 * WHY THIS REPLACED A CALL-RATE FLOOR. The floor refused every array under 40% called, because a
 * failed reaction was taken to drive heterozygosity toward zero and read as one genome. Over 1,001
 * arrays it fired on none, and the arrays with no DNA read 58% to 67% heterozygous, at random, and are
 * refused as noise. The case it stood in for exists at high call rate, a chaotic diploid of GSE148488
 * staged haploid at 88% called, and the alleles are what reach it.
 */
export function allelePloidy(states?: readonly CopyState[]): 'one genome' | 'two genomes' | undefined {
  if (!states) return undefined
  const stated = states.filter((s) => s !== 'undecided' && s !== 'no-template')
  if (stated.length < ALLELE_PLOIDY_MIN_STATED) return undefined
  const n = (s: CopyState): number => stated.filter((x) => x === s).length
  if (n('disomic') / stated.length >= 0.5) return 'two genomes'
  if (Math.max(n('other-lost'), n('parent-lost')) / stated.length >= NORMAL_MAJORITY_MIN) return 'one genome'
  return undefined
}

/**
 * Heterozygosity above which the array is not one genome at all.
 *
 * The ladder is open at the bottom by design, since dropout drives heterozygosity down. It was
 * open at the TOP by oversight, and that is worse: `BOUNDS.find` matched bulk for any value at or
 * above 0.158, so an array reading 52% heterozygous at a 55% call rate was classified as bulk
 * genomic DNA and handed a dropout of 0.008, the most confident parameter set in this module.
 *
 * THE CEILING A SINGLE DIPLOID CAN REACH. This panel's diploid rate is 0.168 and the European
 * anchor is already the highest of the ancestries measured against it, so the rate itself cannot
 * rise much. Drop-in can: at the top of the range measured here, 0.0525 across 113 same-genome
 * pairs, a diploid reads 0.168 + (1 - 0.168) x 0.0525 = 0.212. This sits above that with margin.
 * Beyond it the array is two individuals mixed, a contaminated reaction, or a failure, and the
 * four arrays that prompted this read 0.53-0.56, which is over twice the maximum.
 */
/**
 * Per-chromosome heterozygosity spread, as a coefficient of variation, above which the departure
 * from a uniform genome is LOCALISED to some chromosomes rather than shared by all of them.
 *
 * WHAT THIS ACTUALLY FIRES ON, CORRECTED. The first version of this comment said it detects a
 * localised heterozygosity EXCESS. It does not, and the arithmetic says so without any
 * measurement: at the panel's 0.168 diploid rate across 22 autosomes, elevating k chromosomes by
 * 1.5x reaches a spread of at most 0.204 and by 2.0x at most 0.353, while driving k chromosomes to
 * hemizygous level reaches 0.217 at k=1, 0.538 at k=5 and 0.988 at k=11. This threshold is
 * unreachable by any elevation a copy number can produce and is crossed by COLLAPSE at five of
 * twenty-two. The two cells it classifies read 0.69 and 0.78, far above anything elevation
 * reaches. The constant's own derivation gave this away and was not read: half the chromosomes at
 * one copy is the deficit branch, because a hemizygous chromosome has one haplotype and no
 * heterozygotes to call.
 *
 * SO THE SIGN MATTERS AND THE CONFOUND IS REAL. A monosomy collapses heterozygosity on its
 * chromosome and so does a region that did not amplify. One is a genome and the other is not, and
 * a coefficient of variation is blind to the difference because it measures dispersion without
 * regard to cause. The discriminator is where the markers went: a chromosome present at one copy
 * is still genotyped at a normal rate, while a chromosome that failed to amplify loses its call
 * rate along with its heterozygotes. So the spread is measured over chromosomes that are being
 * genotyped and over no others; see `hetSpread`.
 *
 * Measured on the 31 cells of GSE148488 that carry a hand karyotype, against the authors' own
 * reading. The two cells the ceiling refuses that the hand calls chaotic aneuploid read 0.69 and
 * 0.78 on the unconditioned statistic. The seven other refused cells read 0.014 to 0.155, and the
 * two with the highest heterozygosity, 0.62 and 0.64, read the LOWEST spread at 0.014 and 0.028,
 * which is the uniform signature. The value sits inside that gap and is not fitted: any value in
 * it gives the same answer on this corpus, which is equally why the corpus cannot calibrate it
 * further. What the corpus has NOT established is the threshold's behaviour on the conditioned
 * statistic, which excludes collapsed chromosomes and can therefore only lower a cell's reading.
 * The direction is safe and the magnitude is not measured.
 */
export const HET_SPREAD_LOCALISED = 0.50

/**
 * Coefficient of variation of the heterozygous rate across autosomes.
 *
 * NaN where fewer than five autosomes carry enough calls to measure, which is a refusal to report
 * a spread rather than a spread of zero: a uniform reading and an unmeasurable one must not look
 * alike to the caller, because one of them admits a sample and the other must not.
 */
export const hetSpread = (
  byChrom: Map<string, { markers: number; called: number; het: number }>,
  minCallRate = 0.60,
  minCalled = 500,
): number => {
  // ONLY CHROMOSOMES THAT ARE BEING GENOTYPED. A chromosome whose call rate has collapsed has
  // lost its heterozygotes to dropout, not to hemizygosity, and including it lets amplification
  // failure masquerade as a genome with chromosomes at one copy. Those two want opposite
  // treatment, which is the whole reason this statistic exists. The floor is the same 0.60 the
  // per-chromosome mis-clustering gate uses for a collapsed call rate.
  const rates = [...byChrom.entries()]
    .filter(([c, v]) => /^\d+$/.test(c) && v.called >= minCalled
      && v.markers > 0 && v.called / v.markers >= minCallRate)
    .map(([, v]) => v.het / v.called)
  if (rates.length < 5) return NaN
  const mean = rates.reduce((a, b) => a + b, 0) / rates.length
  if (!(mean > 0)) return NaN
  const sd = Math.sqrt(rates.reduce((a, b) => a + (b - mean) ** 2, 0) / rates.length)
  return sd / mean
}

export const MAX_DIPLOID_HET = 0.25

/**
 * BAF spread at heterozygous calls above which a sample is NOT unamplified genomic DNA.
 *
 * THE STAGE LADDER MEASURED THE WRONG AXIS ON ITS TOP RUNG, and this is the correction. Stage was
 * assigned from heterozygosity alone, which tracks how many heterozygotes SURVIVED and therefore
 * how much template there was. It says nothing about how far the survivors SCATTERED, which is what
 * amplification does and what every constant keyed to stage actually encodes: drift, variance
 * inflation and the detection floors span a seventeen-fold range across these classes.
 *
 * Those two axes come apart. Measured across 877 arrays of this platform, the four classes the
 * ladder assigns have BAF spreads of 0.229, 0.232, 0.241 and 0.251 at the median, which is no
 * separation at all against a within-class spread of 0.10 to 0.13. Worse, 176 arrays reading
 * amplified-level scatter of 0.232 were called BULK on heterozygosity and handed a drift constant
 * of 0.0011, seventeen times tighter than a blastomere's, because they had retained enough
 * heterozygotes to clear the top rung.
 *
 * So bulk now requires BOTH: enough heterozygosity AND unamplified scatter. The threshold is the
 * one this project already carries from MoChA for excluding an array outright, 0.11, sitting
 * between unamplified genomic DNA at about 0.088 and amplified material at 0.21 to 0.30. An array
 * failing it is not rejected here, it is DEMOTED to the amplified rungs, where its heterozygosity
 * decides which one. Nothing is lost and the confident parameter set stops being handed out on the
 * strength of the wrong measurement.
 */
export const BULK_MAX_BAF_SD = 0.11

/**
 * Below this the sample carries one genome rather than a heavily dropped-out two.
 *
 * NOT SAFE IN BOTH DIRECTIONS, and the audit quantified which way it fails. A blastomere at 0.308
 * dropout sits at h = 0.1163, eleven thousandths above the line; add first-cousin consanguinity and
 * it falls to 0.1090, effectively on it; at 0.40 dropout it crosses and is called haploid. A PB1 at
 * expected h 0.074 is below the line, but its upper tail reaches 0.092 before drop-in is added.
 */
export const HAPLOID_MAX_HET = 0.105

/**
 * BAF band at or above which the sample carries two genomes whatever its genotype calls say.
 *
 * Taken from the diploid arm of the triage gate rather than invented here, so the two modules
 * cannot drift apart. Used ONLY to veto a haploid call: see the header for why the converse fails
 * on amplified material.
 */
export const BAND_DIPLOID_CERTAIN = 0.15

/**
 * Dropout from two independent amplifications of ONE genome, which needs no population anchor.
 *
 * Among markers called heterozygous in at least one replicate, the discordant fraction is
 * 2d/(1+d), so d = phi/(2-phi). Immune to consanguinity, LOH, UPD and ancestry, because those
 * change WHICH markers are heterozygous but not the chance a heterozygous one survives twice.
 *
 * TWO DOCUMENTED LIMITS. Correlated failure biases it low in proportion: recovery falls to 0.90x at
 * 10% shared failures and 0.50x at 50%, so it is a lower bound. And it has a floor set by the
 * platform's own genotyping error, which on bulk replicates returns about 0.10 where the true
 * dropout is near zero, so it must not be applied to bulk or ES-line material where that floor
 * dominates the signal.
 */
export function dropoutFromReplicates(discordantFraction: number): number {
  const phi = Math.min(0.999, Math.max(0, discordantFraction))
  return phi / (2 - phi)
}

/**
 * Bulk heterozygosity of THE PANEL IN FRONT OF THE TOOL, taken from the parental array.
 *
 * Every rung in BOUNDS and every ceiling here is a heterozygosity, and heterozygosity is a property
 * of which markers a panel carries. Measured on bulk arrays of this corpus:
 *
 *     panel                       bulk heterozygosity
 *     Axiom (GPL28377)            0.141 to 0.170
 *     Affymetrix 250K Nsp         0.201 to 0.308
 *     Illumina CytoSNP-12         0.301 to 0.316
 *
 * Applied as the Axiom figure everywhere, MAX_DIPLOID_HET sits BELOW the normal bulk range of two
 * of those three, so an undamaged bulk array on either is refused as not a stage: a newborn's own
 * gDNA at a 96.9% call rate and 30.8% heterozygous reads `failed`. The parental array is on the
 * same panel as the sample, so its heterozygosity anchors the table to the panel in use.
 *
 * SCALED UP ONLY, AND CAPPED. A parent's own dropout can only depress its heterozygosity, so a
 * ratio under one is evidence about that array rather than about the panel and the table is left
 * alone. The cap is there because one bad parental array also reads high: the richest panel
 * measured here is 1.85x the anchor.
 *
 * WHAT IT SCALES IS THE CEILING AND THE DROPOUT ANCHOR, NOT THE RUNGS. See inferStage.
 */
/**
 * Call rate an array must reach before its heterozygosity is allowed to anchor the panel.
 *
 * An array that is itself heavily dropped out says nothing about which markers a panel carries.
 * The case that forces this is an operator dropping the SAMPLE into the parental slot: the sample
 * then anchors its own ceiling, and an array at 28 to 31% heterozygous, which on its panel is a
 * failed reaction, certifies itself as ordinary. Bulk arrays across the four panels measured here
 * call 0.81 to 0.99; the arrays that wrongly inflated the anchor called 0.53 to 0.59.
 */
export const PANEL_ANCHOR_CALL_RATE = 0.80

export const PANEL_SCALE_MAX = 2.0
export const panelScale = (panelHeterozygosity?: number): number => (
  Number.isFinite(panelHeterozygosity as number) && (panelHeterozygosity as number) > 0
    ? Math.min(PANEL_SCALE_MAX, Math.max(1, (panelHeterozygosity as number) / BULK_HETEROZYGOSITY))
    : 1
)

/**
 * Whether a PARENTAL array is too heterozygous to be one person's genotype, judged against the
 * panel the other loaded parent defines.
 *
 * Measured on GSE148488: six egg donors and the sperm donor read 0.141 to 0.174 heterozygous at
 * 94 to 98 percent call rate; the seventh donor reads 0.3455 at 82 percent, which is two people or
 * a contaminated reaction, and anything read from her heterozygosity is read against alleles she
 * does not transmit. The ceiling is MAX_DIPLOID_HET scaled by the panel, and the panel is the
 * OTHER parent's heterozygosity, because a parent judged against its own figure cannot exceed it.
 * With one parent loaded there is nothing to judge against, and this is false. The application
 * reports it as a limit and withholds nothing: the array's homozygous calls still make markers
 * informative, and the allele count reads only those.
 */
export const parentOverCeiling = (het: number, panelHet: number): boolean =>
  Number.isFinite(het) && Number.isFinite(panelHet) && het > MAX_DIPLOID_HET * panelScale(panelHet)

const BOUNDS: { stage: Stage, minHet: number, dropout: number, templates: string, floor: number }[] = [
  { stage: 'bulk', minHet: 0.158, dropout: 0.013, templates: '~10^6', floor: 100 },
  { stage: 'trophectoderm', minHet: 0.145, dropout: 0.050, templates: '10-20', floor: 100 },
  { stage: 'single-cell', minHet: 0.126, dropout: 0.199, templates: '2', floor: 200 },
  { stage: 'blastomere', minHet: HAPLOID_MAX_HET, dropout: 0.308, templates: '2', floor: 200 },
]

/**
 * The shipped parameters for a NAMED stage, with no array to read them off.
 *
 * For the path where the operator declares the material rather than the tool inferring it. The
 * figures are the same ones inferStage assigns, read from the same table, so a declared stage and
 * an inferred one of the same name are parameterised identically.
 */
export function stageDefaults(stage: Stage): StageCall {
  const hit = BOUNDS.find((b) => b.stage === stage)
  if (hit) {
    return {
      stage,
      ploidy: 'diploid',
      dropout: hit.dropout,
      basis: 'stage-default',
      templates: hit.templates,
      markerFloor: hit.floor,
      caveat: '',
      why: `${stage} parameters, from the shipped table`,
    }
  }
  if (stage === 'haploid') {
    return {
      stage, ploidy: 'haploid', dropout: 0.02, basis: 'stage-default',
      templates: '1 chromatid', markerFloor: 200,
      caveat: 'a heavily dropped-out or consanguineous diploid can also read this way',
      why: 'haploid parameters, from the shipped table',
    }
  }
  // failed and unknown name the absence of a call, so they carry the most conservative dropout
  // measured on any stage rather than a figure that would flatter a likelihood.
  return {
    // 'failed' here is a declared absence rather than the measured ceiling breach inferStage
    // reports, so it claims nothing about how many genomes are present.
    stage, ploidy: 'unknown', dropout: stage === 'failed' ? NaN : 0.308, basis: 'none',
    templates: stage === 'failed' ? 'not applicable' : 'unknown', markerFloor: 200,
    caveat: '', why: `${stage} carries no stage parameters`,
    ...(stage === 'failed' ? { refusedBy: 'declared' as const } : {}),
  }
}

/**
 * Call the stage from a sample's own profile.
 *
 * QUALITY IS GATED BEFORE PLOIDY, which is the ordering the audit required. A failed amplification
 * drives heterozygosity toward zero, so without this gate it is classified as haploid and given a
 * confident parameter set rather than being rejected. The two are indistinguishable by
 * heterozygosity alone, so call rate has to decide first.
 */
export function inferStage(
  profile: Pick<SampleProfile, 'hetRate' | 'callRate'> & {
    hetBand?: number, hetBafSd?: number,
    /**
     * Parent-free answer keys, read against a founder panel (answerKeys.ts): this array's drop-in,
     * on the scale of the drop-in read between measured parents, and its information ratio. Both
     * optional, and with neither the ceiling behaves exactly as it always has.
     */
    dropIn?: number, information?: number,
    /**
     * The one-parent copy state of each autosome, read against a loaded parent, measured or
     * reconstructed (oneParentState.ts). Supplies the ploidy genotypes cannot: see allelePloidy.
     */
    copyStates?: readonly CopyState[],
  },
  /**
   * The boundaries, so a caller that knows what it is doing can move them.
   *
   * The web tool never passes this and offers no knob for it: the whole point there is that the
   * stage is inferred rather than declared. The command line passes it, because the question
   * asked there is often what happens at a different threshold. Defaults are the measured
   * values, so omitting this argument is exactly the shipped configuration.
   */
  opts: {
    maxDiploidHet?: number, haploidMaxHet?: number,
    bulkHeterozygosity?: number, bandDiploidCertain?: number,
    /** Bulk heterozygosity of this panel, which scales every boundary below. See panelScale. */
    panelHeterozygosity?: number,
    /** BAF spread above which the bulk rung is refused. See BULK_MAX_BAF_SD. */
    bulkMaxBafSd?: number,
  } = {},
): StageCall {
  let h = profile.hetRate
  const call = profile.callRate
  let dropInCorrected: number | undefined
  const band = profile.hetBand
  // An explicit threshold from a caller is used as given; only the defaults follow the panel.
  const scale = panelScale(opts.panelHeterozygosity)
  // THE CEILING AND THE DROPOUT ANCHOR FOLLOW THE PANEL; THE RUNGS DO NOT, and that asymmetry is
  // measured rather than chosen. The ceiling is the defect: it refuses undamaged bulk arrays on a
  // rich panel. The rungs decide which MATERIAL a sample is treated as, and scaling them moved one
  // bulk array of the stress corpus onto an amplified rung, where the relationship verdict is
  // withheld: the tool then failed to warn that an array had been put in its own parental slot, on
  // one run in eight. A gate that scales a ceiling must not quietly disarm a safety check.
  const maxDiploid = opts.maxDiploidHet ?? MAX_DIPLOID_HET * scale
  const haploidMax = opts.haploidMaxHet ?? HAPLOID_MAX_HET
  const bulkHet = opts.bulkHeterozygosity ?? BULK_HETEROZYGOSITY * scale
  const bandCertain = opts.bandDiploidCertain ?? BAND_DIPLOID_CERTAIN

  if (!Number.isFinite(h) || !Number.isFinite(call)) {
    return {
      stage: 'unknown', ploidy: 'unknown', dropout: 0.308, basis: 'none', templates: 'unknown', markerFloor: 200,
      caveat: 'no usable profile',
      why: 'the array does not report a heterozygous rate and a call rate, so no stage is inferred',
    }
  }
  // OVER THE CEILING, THE ANSWER KEYS SAY WHY, where they are supplied. Heterozygosity over what one
  // genome can reach has three physical causes and they call for different answers. Calls that carry
  // nothing are refused. Heterozygosity that is still too high once drop-in is removed is a second
  // genome, refused. And heterozygosity that is excess DROP-IN, calls at markers where this array
  // cannot be heterozygous, is one genome read badly, which is staged with the drop-in removed. On
  // GSE148488 seven of the nine arrays the ceiling refuses are the third kind: drop-in of 0.24 to 0.32
  // against the 0.02 typical of admitted arrays, single genomes whose allele calls agree with the
  // published karyotype at 133 of 138.
  //
  // CALLS THAT ARE MOSTLY NOISE ARE REFUSED WHATEVER THEIR HETEROZYGOSITY. An array with no DNA does
  // not stop calling: the seven of the 1,001 surveyed that are recorded as having none call 42% to 75%
  // of markers, at random. Random calls read heterozygous, so in practice they land over the ceiling,
  // but what refuses them is that they are random, and that holds under the ceiling too.
  if (Number.isFinite(profile.information)
    && (profile.information as number) >= NO_INFORMATION_MIN) {
    return {
      stage: 'failed', ploidy: 'unknown', dropout: NaN, basis: 'none',
      templates: 'not applicable', markerFloor: 200, refusedBy: 'no-genotype-signal',
      caveat: 'measured against a founder panel, so no parent is needed for it',
      why: `${(h * 100).toFixed(1)}% heterozygous, and against markers where unrelated people all `
        + 'carry the same homozygote this array\'s calls agree no better than chance would put them '
        + `(information ratio ${(profile.information as number).toFixed(2)}, at or over `
        + `${NO_INFORMATION_MIN.toFixed(2)}): about ${Math.round(Math.min(1, profile.information as number) * 100)}% `
        + 'of its homozygous calls are random, more noise than genome, so nothing read from them stands',
    }
  }
  if (h > maxDiploid && Number.isFinite(profile.dropIn)
    && (profile.dropIn as number) > 0 && (profile.dropIn as number) < 1) {
    const d = profile.dropIn as number
    const corrected = Math.max(0, (h - d) / (1 - d))
    if (corrected > maxDiploid) {
      return {
        stage: 'failed', ploidy: 'not one genome', dropout: NaN, basis: 'none',
        templates: 'not applicable', markerFloor: 200, refusedBy: 'second-genome',
        caveat: 'drop-in rules out an unrelated contaminant only at markers the panel covers, and a '
          + 'second genome sharing alleles with this one, a relative, adds fewer heterozygotes',
        why: `${(h * 100).toFixed(1)}% heterozygous, and ${(corrected * 100).toFixed(1)}% once this `
          + `array's ${(d * 100).toFixed(1)}% drop-in is removed, still over the `
          + `${(maxDiploid * 100).toFixed(0)}% one genome reaches. More than one genome is present`,
      }
    }
    // ONE GENOME. Staged at the heterozygosity it would show at the drop-in the rungs were measured
    // with, so the rungs need no second calibration.
    dropInCorrected = d
    h = corrected * (1 - DROP_IN) + DROP_IN
  }
  if (h > maxDiploid) {
    return {
      stage: 'failed', ploidy: 'not one genome', dropout: NaN, basis: 'none',
      templates: 'not applicable', markerFloor: 200, refusedBy: 'heterozygosity-ceiling',
      caveat: 'two individuals mixed, a contaminated reaction and a failed array are not '
        + 'distinguishable from each other here, only from a genome',
      why: `${(h * 100).toFixed(1)}% heterozygous exceeds the `
        + `${(maxDiploid * 100).toFixed(0)}% ceiling a single diploid can reach on this `
        + 'platform, which is its 16.8% rate plus drop-in at the highest rate measured here. No '
        + 'genome reads this way, so this is not a stage',
    }
  }

  // Intensity outranks the genotype calls in this one direction. A sample whose BAF band is at
  // the diploid line has intermediate signal at a large fraction of its markers, which one
  // template cannot produce, so it is diploid however few heterozygotes the caller kept. Without
  // this a blastomere at 0.40 dropout crosses the heterozygosity boundary and is called haploid.
  // NOT ON A DROP-IN CORRECTED ARRAY. Drop-in puts heterozygous calls at intermediate BAF, so it
  // fills the band exactly as it inflates heterozygosity: on GSE148488 two paternal pronuclei, haploid
  // by micromanipulation, read the band as diploid at drop-in of 0.30 and 0.43. A band is evidence of
  // two copies only where drop-in is not supplying it.
  const bandSaysDiploid = dropInCorrected === undefined
    && Number.isFinite(band as number) && (band as number) >= bandCertain
  const alleles = allelePloidy(profile.copyStates)

  // DROP-IN CORRECTED HETEROZYGOSITY DOES NOT SEPARATE ONE COPY FROM TWO on this material: haploid
  // pronuclei read 0.002 to 0.014 with drop-in removed and diploid progeny from 0.015. So a corrected
  // array does not reach the haploid rung on heterozygosity; only the BAF band can call its ploidy.
  if (h <= haploidMax && !bandSaysDiploid && dropInCorrected === undefined && alleles !== 'two genomes') {
    // PB2, pronuclei and sperm carry one chromatid and are homozygous throughout. A FIRST polar
    // body carries a dyad and is genuinely heterozygous distal to each crossover, expected h about
    // 0.074, so the two are distinguished by where in the band the sample sits rather than lumped.
    const looksPB1 = h >= 0.055
    return {
      stage: 'haploid',
      ploidy: 'haploid',
      // One chromatid has no heterozygote to drop, so the residual is genotyping error. A PB1's
      // heterozygosity is largely real, so no dropout is inferred from it either.
      dropout: 0.02,
      basis: 'stage-default',
      templates: looksPB1 ? '1 homologue, 2 chromatids' : '1 chromatid',
      markerFloor: 200,
      caveat: looksPB1
        ? 'consistent with a FIRST polar body, whose sister chromatids differ distal to each '
          + 'crossover, so roughly 44% of its genome is genuinely heterozygous and those calls are '
          + 'not error'
        : 'a heavily dropped-out or consanguineous diploid can also fall below this boundary',
      why: `${(h * 100).toFixed(1)}% heterozygous is below the ${(haploidMax * 100).toFixed(1)}% `
        + `a diploid reaches at any stage, so this carries one genome`
        + (looksPB1 ? ', and sits where a first polar body sits rather than at the drop-in floor' : ''),
    }
  }

  let hit = BOUNDS.find((b) => h >= (b.stage === 'blastomere' ? haploidMax : b.minHet))
    ?? BOUNDS[BOUNDS.length - 1]
  // BULK NEEDS THE SECOND AXIS. Heterozygosity alone put 176 amplified arrays on this rung and gave
  // them the tightest drift constant in the module. Where a BAF spread is available and says the
  // sample was amplified, it is demoted to the rungs below and its heterozygosity picks which.
  const sd = profile.hetBafSd
  const amplified = sd !== undefined && Number.isFinite(sd) && sd > (opts.bulkMaxBafSd ?? BULK_MAX_BAF_SD)
  let demotedFrom: string | null = null
  if (amplified && hit.stage === 'bulk') {
    demotedFrom = 'bulk'
    hit = BOUNDS.find((b) => b.stage !== 'bulk'
      && h >= (b.stage === 'blastomere' ? haploidMax : b.minHet)) ?? BOUNDS[BOUNDS.length - 1]
  }
  // ONE GENOME BY ITS ALLELES, UNDER THE CEILING TOO. Drop-in and failed calls lift a haploid's
  // heterozygosity past the haploid line without taking it over the ceiling, and lift its BAF band.
  // Where the alleles read one parent's genome, heterozygosity is read with the measured drop-in
  // removed: at the haploid level it is one genome. A cell carrying two copies of one parent's genome
  // reads diploid heterozygosity with drop-in removed, and stays diploid.
  if (alleles === 'one genome' && dropInCorrected === undefined && Number.isFinite(profile.dropIn)
    && (profile.dropIn as number) > 0 && (profile.dropIn as number) < 1) {
    const d = profile.dropIn as number
    const hc = Math.max(0, (h - d) / (1 - d))
    if (hc <= haploidMax) {
      return {
        stage: 'haploid', ploidy: 'haploid', dropout: 0.02, basis: 'stage-default',
        templates: '1 chromatid', markerFloor: 200,
        caveat: 'staged by its alleles against the loaded parent and its heterozygosity with drop-in removed',
        why: `${(h * 100).toFixed(1)}% heterozygous as read, ${(hc * 100).toFixed(1)}% with this array's `
          + `${(d * 100).toFixed(1)}% drop-in removed, under the ${(haploidMax * 100).toFixed(1)}% a `
          + 'diploid reaches; and against the loaded parent most chromosomes carry one parent\'s copy '
          + 'alone. One genome',
      }
    }
  }
  if ((bandSaysDiploid || alleles === 'two genomes') && h <= haploidMax && dropInCorrected === undefined) {
    return {
      stage: hit.stage,
      ploidy: 'diploid',
      ...(dropInCorrected !== undefined ? { dropInCorrected } : {}),
      // The shortfall formula is calibrated against a diploid's heterozygosity, and this sample
      // is diploid, so it applies. It saturates at the 0.6 ceiling here, which is the point.
      dropout: Math.min(0.6, Math.max(0.005, 1 - h / bulkHet)),
      basis: 'heterozygosity-shortfall',
      templates: hit.templates,
      markerFloor: hit.floor,
      caveat: `genotype heterozygosity alone would have called this one genome; the ${bandSaysDiploid
        ? 'BAF band' : 'alleles against the parent'} overrode that. Dropout is at or near its ceiling, so `
        + 'every downstream call on this sample is weakly powered even where it is admitted',
      why: bandSaysDiploid
        ? `${(h * 100).toFixed(1)}% heterozygous is in the haploid range, but a BAF band of `
          + `${((band as number) * 100).toFixed(1)}% is at or above the ${(bandCertain * 100)
            .toFixed(0)}% diploid line. Intermediate intensity at that many markers cannot come from `
          + 'one template, so this is a heavily dropped-out diploid rather than one genome'
        : `${(h * 100).toFixed(1)}% heterozygous is in the haploid range, but against the loaded `
          + 'parent most chromosomes carry a second genome\'s alleles, which one genome cannot. Allele '
          + 'dropout removes alleles at random and leaves the other parent\'s in the calls, so this is '
          + 'a heavily dropped-out diploid rather than one genome',
    }
  }
  const implied = Math.min(0.6, Math.max(0.005, 1 - h / bulkHet))
  // ONE GENOME BY ITS ALLELES, where drop-in has hidden it from heterozygosity and the band. A haploid
  // product reads one parent's copy alone on every chromosome against a loaded parent, and drop-in
  // cannot add the other parent's alleles: the four drop-in paternal pronuclei of GSE148488 read so.
  if (dropInCorrected !== undefined && alleles === 'one genome') {
    return {
      stage: 'haploid', ploidy: 'haploid', dropInCorrected, dropout: 0.02, basis: 'stage-default',
      templates: '1 chromatid', markerFloor: 200,
      caveat: 'staged by its alleles against the loaded parent, since drop-in hides ploidy from '
        + 'heterozygosity and the BAF band',
      why: `${(profile.hetRate * 100).toFixed(1)}% heterozygous as read, over the ceiling from `
        + `${(dropInCorrected * 100).toFixed(1)}% drop-in. Against the loaded parent most chromosomes `
        + 'carry one parent\'s copy alone, which is one genome',
    }
  }
  if (dropInCorrected !== undefined) {
    return {
      stage: hit.stage,
      // Neither heterozygosity nor the BAF band says how many copies once drop-in is this heavy:
      // drop-in corrupts both. The alleles do where a parent is loaded; otherwise unknown.
      ploidy: alleles === 'two genomes' ? 'diploid' : 'unknown',
      dropInCorrected,
      dropout: Math.max(implied, hit.dropout * 0.6),
      basis: 'heterozygosity-shortfall',
      templates: hit.templates,
      markerFloor: hit.floor,
      caveat: 'staged with drop-in removed. Copy number is not read from this array: whether it has a '
        + 'normal majority to normalise against has not been shown',
      why: `${(profile.hetRate * 100).toFixed(1)}% heterozygous as read, over the ceiling, but `
        + `${(dropInCorrected * 100).toFixed(1)}% of calls at markers where it cannot be heterozygous `
        + 'are heterozygous anyway: drop-in. With that removed it reads as one genome, placed at '
        + `${hit.stage} with ${(implied * 100).toFixed(0)}% dropout. How many copies it carries is not `
        + 'separable from its heterozygosity once drop-in this heavy is removed',
    }
  }
  return {
    stage: hit.stage,
    ploidy: 'diploid',
    dropout: Math.max(implied, hit.dropout * 0.6),
    basis: 'heterozygosity-shortfall',
    templates: hit.templates,
    markerFloor: hit.floor,
    caveat: 'consanguinity, copy-neutral LOH, UPD and ancestry differing from the anchor all '
      + 'depress heterozygosity and are absorbed into this estimate; an East Asian sample against '
      + 'this European-derived anchor is biased by about +0.155, enough to shift the stage by one '
      + 'rung. Where a same-genome replicate exists, prefer dropoutFromReplicates',
    why: `${(h * 100).toFixed(1)}% heterozygous against ${(bulkHet * 100).toFixed(1)}% `
      + `for bulk DNA on this panel implies ${(implied * 100).toFixed(1)}% dropout`
      + (demotedFrom
        ? `, which alone would have read as ${demotedFrom}. But the allele fraction at those `
          + `heterozygous calls spreads by ${(sd as number).toFixed(4)}, over the `
          + `${(opts.bulkMaxBafSd ?? BULK_MAX_BAF_SD)} that separates unamplified DNA from `
          + 'amplified material, so this sample was amplified and cannot take the bulk parameter '
          + `set. It is placed at ${hit.stage} by its heterozygosity instead`
        : `, which is where ${hit.stage} material sits`)
      + '. Template count is reported for context and does not predict this figure',
  }
}

/**
 * A region in the conventional genome-browser form, chr6:39,302-294,904.
 *
 * Thousands separators are included because these are read by people rather than parsed, and an
 * unseparated nine-digit coordinate is where a reader loses an order of magnitude.
 */
export const locus = (chrom: string, startBp: number, endBp: number): string =>
  `chr${chrom}:${Math.round(startBp).toLocaleString('en-US')}`
  + `-${Math.round(endBp).toLocaleString('en-US')}`

/**
 * The stage call reduced to what a reader needs to see, without the component knowing the rules.
 *
 * SHOWN AT THE TOP OF A RUN BECAUSE IT SETS EVERYTHING BELOW IT. The stage picks the drift floor,
 * the variance inflation, the detection floors and the dropout that parameterises every directional
 * call, and those differ by up to seventeen-fold across the classes. A reader who does not know
 * which stage was inferred cannot judge any number underneath it.
 *
 * BOTH AXES ARE RETURNED SEPARATELY, and that is the point rather than a formatting choice. They
 * are easy to confuse and they measure different things: heterozygosity counts how many
 * heterozygous markers SURVIVED, which tracks template count, while the allele-fraction spread at
 * those markers says how far the survivors SCATTERED, which tracks amplification. Reading only the
 * first is what put amplified samples on the bulk rung.
 */
export interface StageFacts {
  stage: Stage
  /** Survived: how many heterozygous markers are left, which tracks template count. */
  hetRate: number
  /** Scattered: how far they spread, which tracks amplification. NaN where no BAF was supplied. */
  hetBafSd: number
  amplified: boolean | null
  /** Set when heterozygosity alone would have said bulk and the spread overruled it. */
  demotedFromBulk: boolean
  dropout: number
  basis: string
  templates: string
  markerFloor: number
  caveat: string
  why: string
}

export function stageFacts(
  call: StageCall,
  profile: { hetRate: number; hetBafSd?: number },
  bulkMaxBafSd = BULK_MAX_BAF_SD,
): StageFacts {
  const sd = profile.hetBafSd
  const known = sd !== undefined && Number.isFinite(sd)
  return {
    stage: call.stage,
    hetRate: profile.hetRate,
    hetBafSd: known ? (sd as number) : NaN,
    amplified: known ? (sd as number) > bulkMaxBafSd : null,
    // Read from the reason rather than recomputed, so the panel can never disagree with the call.
    demotedFromBulk: call.why.includes('would have read as bulk'),
    dropout: call.dropout,
    basis: call.basis,
    templates: call.templates,
    markerFloor: call.markerFloor,
    caveat: call.caveat,
    why: call.why,
  }
}

/** Plain-language reading of the amplification axis, for a reader who does not know the scale. */
export const amplificationWord = (f: StageFacts): string => {
  if (f.amplified === null) return 'not measured, no allele fractions in this file'
  return f.amplified
    ? 'whole-genome amplified, which scatters at 0.21 to 0.30'
    : 'unamplified, consistent with genomic DNA at about 0.088'
}

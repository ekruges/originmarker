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

export type Stage =
  | 'bulk'
  | 'trophectoderm'
  | 'single-cell'
  | 'blastomere'
  | 'haploid'
  | 'failed'
  | 'unknown'

export interface StageCall {
  stage: Stage
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

/**
 * Heterozygosity of a diploid genome on this platform, measured on bulk gDNA under this project's
 * marker QC.
 *
 * A PANEL AND ANCESTRY PROPERTY, NOT A CONSTANT OF NATURE. Common-SNP panels run 0.32 to 0.44
 * expected heterozygosity in Europeans, so 0.168 implies roughly half this panel is rare content;
 * any change to marker QC moves it. Relative to a European anchor, mean expected heterozygosity is
 * 0.99 in South Asian, 0.93 in African and 0.81 in East Asian samples, and that last ratio alone
 * biases a shortfall estimate by +0.155.
 */
export const BULK_HETEROZYGOSITY = 0.168

/** Call rate below which no stage is inferred, because a failed reaction imitates haploid material. */
export const QC_CALL_FLOOR = 0.40

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

const BOUNDS: { stage: Stage, minHet: number, dropout: number, templates: string, floor: number }[] = [
  { stage: 'bulk', minHet: 0.158, dropout: 0.013, templates: '~10^6', floor: 100 },
  { stage: 'trophectoderm', minHet: 0.145, dropout: 0.050, templates: '10-20', floor: 100 },
  { stage: 'single-cell', minHet: 0.126, dropout: 0.199, templates: '2', floor: 200 },
  { stage: 'blastomere', minHet: HAPLOID_MAX_HET, dropout: 0.308, templates: '2', floor: 200 },
]

/**
 * Call the stage from a sample's own profile.
 *
 * QUALITY IS GATED BEFORE PLOIDY, which is the ordering the audit required. A failed amplification
 * drives heterozygosity toward zero, so without this gate it is classified as haploid and given a
 * confident parameter set rather than being rejected. The two are indistinguishable by
 * heterozygosity alone, so call rate has to decide first.
 */
export function inferStage(
  profile: Pick<SampleProfile, 'hetRate' | 'callRate'> & { hetBand?: number, hetBafSd?: number },
  /**
   * The boundaries, so a caller that knows what it is doing can move them.
   *
   * The web tool never passes this and offers no knob for it: the whole point there is that the
   * stage is inferred rather than declared. The command line passes it, because the question
   * asked there is often what happens at a different threshold. Defaults are the measured
   * values, so omitting this argument is exactly the shipped configuration.
   */
  opts: {
    callFloor?: number, maxDiploidHet?: number, haploidMaxHet?: number,
    bulkHeterozygosity?: number, bandDiploidCertain?: number,
    /** BAF spread above which the bulk rung is refused. See BULK_MAX_BAF_SD. */
    bulkMaxBafSd?: number,
  } = {},
): StageCall {
  const h = profile.hetRate
  const call = profile.callRate
  const band = profile.hetBand
  const callFloor = opts.callFloor ?? QC_CALL_FLOOR
  const maxDiploid = opts.maxDiploidHet ?? MAX_DIPLOID_HET
  const haploidMax = opts.haploidMaxHet ?? HAPLOID_MAX_HET
  const bulkHet = opts.bulkHeterozygosity ?? BULK_HETEROZYGOSITY
  const bandCertain = opts.bandDiploidCertain ?? BAND_DIPLOID_CERTAIN

  if (!Number.isFinite(h) || !Number.isFinite(call)) {
    return {
      stage: 'unknown', dropout: 0.308, basis: 'none', templates: 'unknown', markerFloor: 200,
      caveat: 'no usable profile',
      why: 'the array does not report a heterozygous rate and a call rate, so no stage is inferred',
    }
  }
  if (h > maxDiploid) {
    return {
      stage: 'failed', dropout: NaN, basis: 'none', templates: 'not applicable', markerFloor: 200,
      caveat: 'two individuals mixed, a contaminated reaction and a failed array are not '
        + 'distinguishable from each other here, only from a genome',
      why: `${(h * 100).toFixed(1)}% heterozygous exceeds the `
        + `${(maxDiploid * 100).toFixed(0)}% ceiling a single diploid can reach on this `
        + 'platform, which is its 16.8% rate plus drop-in at the highest rate measured here. No '
        + 'genome reads this way, so this is not a stage',
    }
  }
  if (call < callFloor) {
    return {
      stage: 'failed', dropout: NaN, basis: 'none', templates: 'not applicable', markerFloor: 200,
      caveat: 'amplification failure and haploid material are indistinguishable by heterozygosity',
      why: `call rate ${(call * 100).toFixed(1)}% is below the ${(callFloor * 100).toFixed(0)}% `
        + 'floor. This is an amplification failure rather than a stage: near-total dropout drives '
        + 'heterozygosity toward zero and would otherwise be read as one genome',
    }
  }

  // Intensity outranks the genotype calls in this one direction. A sample whose BAF band is at
  // the diploid line has intermediate signal at a large fraction of its markers, which one
  // template cannot produce, so it is diploid however few heterozygotes the caller kept. Without
  // this a blastomere at 0.40 dropout crosses the heterozygosity boundary and is called haploid.
  const bandSaysDiploid = Number.isFinite(band as number) && (band as number) >= bandCertain

  if (h <= haploidMax && !bandSaysDiploid) {
    // PB2, pronuclei and sperm carry one chromatid and are homozygous throughout. A FIRST polar
    // body carries a dyad and is genuinely heterozygous distal to each crossover, expected h about
    // 0.074, so the two are distinguished by where in the band the sample sits rather than lumped.
    const looksPB1 = h >= 0.055
    return {
      stage: 'haploid',
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
  if (bandSaysDiploid && h <= haploidMax) {
    return {
      stage: hit.stage,
      // The shortfall formula is calibrated against a diploid's heterozygosity, and this sample
      // is diploid, so it applies. It saturates at the 0.6 ceiling here, which is the point.
      dropout: Math.min(0.6, Math.max(0.005, 1 - h / bulkHet)),
      basis: 'heterozygosity-shortfall',
      templates: hit.templates,
      markerFloor: hit.floor,
      caveat: 'genotype heterozygosity alone would have called this one genome; the BAF band '
        + 'overrode that. Dropout is at or near its ceiling, so every downstream call on this '
        + 'sample is weakly powered even where it is admitted',
      why: `${(h * 100).toFixed(1)}% heterozygous is in the haploid range, but a BAF band of `
        + `${((band as number) * 100).toFixed(1)}% is at or above the ${(bandCertain * 100)
          .toFixed(0)}% diploid line. Intermediate intensity at that many markers cannot come from `
        + 'one template, so this is a heavily dropped-out diploid rather than one genome',
    }
  }
  const implied = Math.min(0.6, Math.max(0.005, 1 - h / bulkHet))
  return {
    stage: hit.stage,
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

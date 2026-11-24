/**
 * Every bar the scoring uses, by name, with its value read from the constant itself.
 *
 * A bar is a named number a measurement is compared against to decide something: a call, a
 * refusal, a grade. The report and the methods paragraph print bars from here and nowhere else, so
 * a printed value cannot drift from the value the scoring used. Every other named numeric constant
 * of the scoring code is listed in NOT_BARS with the reason it decides nothing on its own.
 * bars.check.ts fails when a named constant of the scoring code is in neither list.
 *
 * `need` is the setup in which the step holding the bar runs: a parent loaded (any, one, both),
 * other arrays in the run, or the feature comparison run. A bar with no `need` applies to every
 * sample. Bare numbers written inside a condition are not named and are out of reach here.
 */
import * as abnormalities from './abnormalities.ts'
import * as answerKeys from './answerKeys.ts'
import * as comparison from './comparison.ts'
import * as dosageOrigin from './dosageOrigin.ts'
import * as features from './features.ts'
import * as gainOrigin from './gainOrigin.ts'
import * as informationRatio from './informationRatio.ts'
import * as ingest from './ingest.ts'
import * as intensityNull from './intensityNull.ts'
import * as motherFromRun from './motherFromRun.ts'
import * as obligateHet from './obligateHet.ts'
import * as oneParentOrigin from './oneParentOrigin.ts'
import * as oneParentState from './oneParentState.ts'
import * as originPosterior from './originPosterior.ts'
import * as otherParentLoss from './otherParentLoss.ts'
import * as parentGenome from './parentGenome.ts'
import * as parentage from './parentage.ts'
import * as partLost from './partLost.ts'
import * as platform from './platform.ts'
import * as refine from './refine.ts'
import * as relatedness from './relatedness.ts'
import * as scoreSample from './scoreSample.ts'
import * as segments from './segments.ts'
import * as sexChromosome from './sexChromosome.ts'
import * as sexing from './sexing.ts'
import * as siblingOrigin from './siblingOrigin.ts'
import * as stage from './stage.ts'
import * as timing from './timing.ts'
import * as untransmitted from './untransmitted.ts'

const MODULES: Record<string, Record<string, unknown>> = {
  'abnormalities.ts': abnormalities, 'answerKeys.ts': answerKeys, 'comparison.ts': comparison,
  'dosageOrigin.ts': dosageOrigin, 'features.ts': features, 'gainOrigin.ts': gainOrigin,
  'informationRatio.ts': informationRatio, 'ingest.ts': ingest, 'intensityNull.ts': intensityNull,
  'motherFromRun.ts': motherFromRun, 'obligateHet.ts': obligateHet, 'oneParentOrigin.ts': oneParentOrigin,
  'oneParentState.ts': oneParentState, 'originPosterior.ts': originPosterior,
  'otherParentLoss.ts': otherParentLoss, 'parentGenome.ts': parentGenome, 'parentage.ts': parentage,
  'partLost.ts': partLost, 'platform.ts': platform, 'refine.ts': refine, 'relatedness.ts': relatedness,
  'scoreSample.ts': scoreSample, 'segments.ts': segments, 'sexChromosome.ts': sexChromosome,
  'sexing.ts': sexing, 'siblingOrigin.ts': siblingOrigin, 'stage.ts': stage, 'timing.ts': timing,
  'untransmitted.ts': untransmitted,
}

export interface Need { parent?: 'any' | 'one' | 'two'; run?: true; comparison?: true }
export interface Bar {
  file: string
  name: string
  /** Read from the module's own export, never written here. */
  value: number
  /** How the value prints: 'bp' prints in Mb; any other unit follows the number. */
  unit?: string
  need?: Need
  says: string
}

const P: Need = { parent: 'any' }
const ONE: Need = { parent: 'one' }
const TWO: Need = { parent: 'two' }
const RUN_P: Need = { parent: 'any', run: true }
const RUN: Need = { run: true }
const CMP: Need = { comparison: true }

type Spec = [file: string, name: string, says: string, need?: Need, unit?: string]
const SPECS: Spec[] = [
  // abnormalities.ts
  ['abnormalities.ts', 'LCSH_REPORT_MB', 'Shortest run of homozygosity reported, in Mb.', undefined, 'Mb'],
  ['abnormalities.ts', 'LOH_DEPLETION', 'Heterozygosity depletion a window must reach before copy-neutral loss of heterozygosity is called.'],
  ['abnormalities.ts', 'LOH_MIN_BACKGROUND_HET', 'Least background heterozygosity at which a depletion is read at all.'],
  ['abnormalities.ts', 'LOH_MIN_EXPECTED_HET', 'Least expected heterozygous markers in a window before its rate is used.', undefined, 'markers'],
  ['abnormalities.ts', 'COMPLEX_DEVIANT_FRACTION', 'Share of autosomes changed at or over which the sample is not reportable.'],
  ['abnormalities.ts', 'COMPLEX_CALL_RATE', 'Call rate under which the sample is not reportable.'],
  ['abnormalities.ts', 'DIPLOID_BAND', 'Centre of the allele-fraction band a diploid genome fills.'],
  ['abnormalities.ts', 'PLOIDY_BAND_HALF_WIDTH', 'Half-width of the allele-fraction bands read for ploidy.'],
  ['abnormalities.ts', 'TRIPLOID_OCCUPANCY', 'Occupancy of the two triploid bands over which a third genome is present.'],
  ['abnormalities.ts', 'DIPLOID_VACATED', 'Occupancy of the half band under which a diploid genome is excluded.'],
  ['abnormalities.ts', 'MIN_MARKERS_PLOIDY', 'Fewest markers before ploidy is read from the bands.', undefined, 'markers'],
  ['abnormalities.ts', 'CONSANGUINITY_CHROMOSOMES', 'Chromosomes with long homozygous runs at once that read as shared ancestry, not one-parent disomy.'],
  ['abnormalities.ts', 'SAME_EMBRYO_AGREEMENT', 'Genotype agreement over which two arrays of a run are units of one embryo.', RUN],
  ['abnormalities.ts', 'MIN_MARKERS_TRIPLOID_ORIGIN', 'Fewest markers before the parent of a third genome is read.', P, 'markers'],
  // answerKeys.ts
  ['answerKeys.ts', 'MIN_SAME_HOM', 'Fewest same-homozygote panel markers before drop-in is measured.', undefined, 'markers'],
  ['answerKeys.ts', 'MIN_FOUNDERS_CALLED', 'Founders called at a marker before it counts as monomorphic.'],
  // comparison.ts
  ['comparison.ts', 'ALPHA', 'Significance at which a coincidence with a feature set is called, before correction.', CMP],
  ['comparison.ts', 'REGION_JOIN_BP', 'Regions closer than this on one chromosome are joined.', CMP, 'bp'],
  ['comparison.ts', 'RELATED_MIN_FOLD', 'Least fold enrichment before an overlap says anything about one region.', CMP],
  // dosageOrigin.ts
  ['dosageOrigin.ts', 'Z_DETECT', 'Least |z| of an allele-dosage shift before an imbalance is reported.', P],
  ['dosageOrigin.ts', 'MAX_HET_BAF_SD', 'Allele-fraction spread at heterozygous calls over which the dosage channel is not run.', P],
  // features.ts
  ['features.ts', 'MIN_PERMUTATIONS', 'Fewest permutations that can reach significance at all.', CMP],
  ['features.ts', 'MIN_REGIONS', 'Fewest regions before an enrichment is reported.', CMP],
  // gainOrigin.ts
  ['gainOrigin.ts', 'MIN_INFORMATIVE_DEFAULT', 'Informative markers a window needs before a direction is read, single cells.', P, 'markers'],
  ['gainOrigin.ts', 'MIN_INFORMATIVE_TROPHECTODERM', 'The same, for a trophectoderm biopsy.', P, 'markers'],
  ['gainOrigin.ts', 'MIN_INFORMATIVE_BULK', 'The same, for unamplified DNA.', P, 'markers'],
  ['gainOrigin.ts', 'SHARE_MARGIN', 'Distance of a parent share from the cell\'s own centre before a gain\'s parent is named.', TWO],
  ['gainOrigin.ts', 'HETERO_MULTIPLE', 'Heterozygosity over the background, as a multiple, before an extra copy is the parent\'s other homologue.', P],
  ['gainOrigin.ts', 'HETERO_ABSOLUTE_MIN', 'Least heterozygosity before an extra copy is the parent\'s other homologue.', P],
  ['gainOrigin.ts', 'HET_BACKGROUND_FLOOR', 'Floor on the background heterozygosity a region is judged against.', P],
  ['gainOrigin.ts', 'ALLELE_GAIN_Z', 'Least z of a chromosome\'s parent share before the alleles call an extra whole copy.', TWO],
  ['gainOrigin.ts', 'ONE_PARENT_ORIENT_Z', 'Least z of the share before a gain already called by intensity names a parent, one parent loaded.', ONE],
  ['gainOrigin.ts', 'ALLELE_GAIN_SHIFT', 'Smallest share shift that is a gain however wide the cell\'s own scale.', P],
  ['gainOrigin.ts', 'ALLELE_GAIN_MIN_INFORMATIVE', 'Informative markers a chromosome needs before its share is read for a gain.', P, 'markers'],
  ['gainOrigin.ts', 'ALLELE_GAIN_MIN_REFERENCE', 'Other chromosomes with both parents a cell needs before it is its own reference.', P],
  ['gainOrigin.ts', 'ALLELE_GAIN_SCALE_FLOOR', 'Floor on the cell\'s own share spread.', P],
  ['gainOrigin.ts', 'ALLELE_GAIN_LRR_SCALE_FLOOR', 'Floor on the cell\'s own intensity spread.', P, 'log2'],
  ['gainOrigin.ts', 'ALLELE_GAIN_LRR_Z', 'Least intensity z, in the cell\'s own spread, before an allele gain is reported.', P],
  ['gainOrigin.ts', 'PAIR_DIRECTION_P', 'Posterior at which an extra copy against a missing copy is read as a call.', TWO],
  // informationRatio.ts
  ['informationRatio.ts', 'NO_INFORMATION_MIN', 'Share of random homozygous calls at or over which the genotype calls carry no information.'],
  ['informationRatio.ts', 'FLIPPED_CODING_MIN', 'Ratio over which the reference\'s allele coding reads flipped against the sample.'],
  ['informationRatio.ts', 'MIN_SHARED_HOMOZYGOTES', 'Fewest shared homozygous markers before the information ratio is read.', undefined, 'markers'],
  // ingest.ts
  ['ingest.ts', 'CALL_RATE_USABLE', 'Call rate at or over which a sample\'s calls are usable without remark in the quality gates.'],
  // intensityNull.ts
  ['intensityNull.ts', 'Z_CHROMOSOME', 'Least calibrated intensity z of a whole chromosome against the array\'s own null.'],
  ['intensityNull.ts', 'Z_WINDOW', 'The same for a window inside a chromosome.'],
  ['intensityNull.ts', 'STATE_GAP', 'Distance between chromosome medians over which they are not one copy state.', undefined, 'log2'],
  // motherFromRun.ts
  ['motherFromRun.ts', 'SHOWS_OTHER_MIN', 'Fitted signal at or over which a unit shows the other parent\'s alleles.', RUN_P],
  ['motherFromRun.ts', 'SAME_MOTHER_MEASURED', 'Window share over which two units have one mother, measured father.', RUN_P],
  ['motherFromRun.ts', 'SAME_MOTHER_REBUILT', 'The same with a rebuilt father.', RUN_P],
  ['motherFromRun.ts', 'SHARED_COPY_K', 'Window score that reads as a shared copy.', RUN_P],
  ['motherFromRun.ts', 'POSTERIOR_MIN', 'Posterior over which a marker enters as her opposite homozygote.', RUN_P],
  ['motherFromRun.ts', 'MIN_DRAWS', 'Fewest draws of her before the run is read for her genotype.', RUN_P],
  // obligateHet.ts
  ['obligateHet.ts', 'HAPLOID_MAX', 'Share under which a chromosome carries one parent\'s genome alone, both parents loaded.', TWO],
  ['obligateHet.ts', 'DIPLOID_MIN', 'Share at or over which both parents are present, both parents loaded.', TWO],
  ['obligateHet.ts', 'ONE_PARENT_HAPLOID_MAX', 'Share under which one genome is present, one parent loaded.', ONE],
  ['obligateHet.ts', 'ONE_PARENT_DIPLOID_MIN', 'Share at or over which two are present, one parent loaded.', ONE],
  ['obligateHet.ts', 'MIN_INFORMATIVE', 'Fewest informative markers before a chromosome\'s share means anything.', P, 'markers'],
  ['obligateHet.ts', 'SURVIVOR_MIN', 'Share reading one parent alone before that parent is the survivor and the other\'s copy is lost.', TWO],
  ['obligateHet.ts', 'MIXED_SHARE', 'Share at which a chromosome between the bounds is reported as a mixed signal.', TWO],
  ['obligateHet.ts', 'ABSENT_PARENT_OPP_MIN', 'Opposite-homozygote rate over which the loaded parent\'s own copy is the absent one.', P],
  // oneParentOrigin.ts
  ['oneParentOrigin.ts', 'DROPIN_CEILING', 'Top of the measured drop-in range a region\'s heterozygosity is sized against.', ONE],
  ['oneParentOrigin.ts', 'MAX_REGION_HET', 'Heterozygosity over which a region is not scored for a missing copy, one parent loaded.', ONE],
  ['oneParentOrigin.ts', 'CALL_POSTERIOR', 'Posterior a hypothesis must reach before the genotype channel calls a region.', P],
  ['oneParentOrigin.ts', 'REGION_MIN_MARKERS', 'Fewest informative markers before the genotype channel scores a region.', P, 'markers'],
  // oneParentState.ts
  ['oneParentState.ts', 'COMMON_ALLELE_MIN', 'Founder frequency at or over which the allele the parent lacks is common enough to read a copy by.', P],
  ['oneParentState.ts', 'MIN_SHARE_MARKERS', 'Markers needed in each set before a chromosome\'s copy state is read.', P, 'markers'],
  ['oneParentState.ts', 'NO_TEMPLATE_BASE', 'Base over which a chromosome carries no template.', P],
  ['oneParentState.ts', 'OTHER_LOST_MAX', 'Ratio under which the other parent\'s copy is gone.', P],
  ['oneParentState.ts', 'PARENT_LOST_MIN', 'Ratio at or over which the loaded parent\'s copy is gone.', P],
  ['oneParentState.ts', 'PARENT_LOST_LOW', 'Low end of where a chromosome missing the loaded parent\'s copy reads.', P],
  ['oneParentState.ts', 'DISOMIC_LOW', 'Low end of where a chromosome with both copies reads.', P],
  ['oneParentState.ts', 'NORMAL_MAJORITY_MIN', 'Share of autosomes in the commonest state at or over which copy number is read.', P],
  // originPosterior.ts
  ['originPosterior.ts', 'BAND_F_MIN', 'Probability under which an origin is graded F.', P],
  ['originPosterior.ts', 'BAND_A_MIN', 'Probability at or over which an origin is graded A.', P],
  ['originPosterior.ts', 'BAND_B_MIN', 'Probability at or over which an origin is graded B.', P],
  ['originPosterior.ts', 'BAND_C_MIN', 'Probability at or over which an origin is graded C.', P],
  ['originPosterior.ts', 'VETO_MAX_F', 'Implied fraction under which a possible gain cannot be told from a possible loss.', P],
  ['originPosterior.ts', 'VETO_DIRECTION_Z', 'Least intensity |z| before the intensity has resolved a direction.', P],
  ['originPosterior.ts', 'MAX_CLASS_PRIOR', 'Cap on any one copy-number class\'s prior.', P],
  // otherParentLoss.ts
  ['otherParentLoss.ts', 'OTHER_LRT', 'Divided score a run missing the other parent\'s copy must reach.', P],
  ['otherParentLoss.ts', 'OTHER_MIN_SPAN_BP', 'Shortest run missing the other parent\'s copy that is reported.', P, 'bp'],
  ['otherParentLoss.ts', 'SAME_GENOME_SHARE', 'Share of blocks over which another array of the run is the same genome.', RUN_P],
  ['otherParentLoss.ts', 'MIN_BLOCKS', 'Blocks needed before the clustering of a run is read.', P],
  ['otherParentLoss.ts', 'NOISE_FLOOR', 'Lowest intercept taken from a fit.', P],
  ['otherParentLoss.ts', 'MIN_FIT_MARKERS', 'Homozygous-class markers a reference chromosome needs before its line is fitted.', P, 'markers'],
  ['otherParentLoss.ts', 'MIN_HET_MARKERS', 'Heterozygous-class markers a reference chromosome needs before its share is read.', P, 'markers'],
  ['otherParentLoss.ts', 'MIN_OPP_CALLS', 'Calls at the parent\'s homozygous markers a run needs before the parent\'s own copy is there.', P],
  // parentGenome.ts
  ['parentGenome.ts', 'RARE_MAX', 'Founder frequency at or under which the parent\'s homozygous allele is rare.', P],
  ['parentGenome.ts', 'RARE_MIN', 'Founder frequency at or over which that allele is tallied.', P],
  ['parentGenome.ts', 'MIN_CARRY_MARKERS', 'Calls needed in each class before a chromosome is read for the parent\'s genome.', P],
  ['parentGenome.ts', 'CONTROL_MAX', 'Share of rare homozygotes, where the parent is the common one, over which the calls are noise.', P],
  ['parentGenome.ts', 'CARRY_Z', 'Standard errors from zero before a chromosome carries the parent\'s genome or not.', P],
  ['parentGenome.ts', 'PRESENT_MIN', 'Chromosomes that must carry the parent\'s genome before the array does.', P],
  ['parentGenome.ts', 'ENRICHMENT_MIN', 'Enrichment at or over which the array carries the parent\'s genome.', P],
  ['parentGenome.ts', 'MIN_ENRICHMENT_HETS', 'Heterozygous calls at the parent\'s rare homozygotes before the enrichment is read.', P],
  ['parentGenome.ts', 'MIN_ENRICHMENT_HET_SHARE', 'Their share of the calls there before the enrichment is read.', P],
  ['parentGenome.ts', 'ABSENT_MIN', 'Chromosomes that must read without the parent\'s genome before the array holds none.', P],
  // parentage.ts
  ['parentage.ts', 'ABSENCE_ERROR_FLOOR', 'Residual absence from genotyping error, added to the noise bound.', P],
  ['parentage.ts', 'ABSENCE_MARGIN', 'How many times the noise bound an absence must reach before it is called.', P, 'x'],
  ['parentage.ts', 'HET_BAND_DIPLOID', 'Heterozygous band a clean diploid array reaches, the reference the report prints beside it.'],
  ['parentage.ts', 'HET_BAND_LO', 'Low edge of the heterozygous allele-fraction band.'],
  ['parentage.ts', 'HET_BAND_HI', 'High edge of the heterozygous allele-fraction band.'],
  ['parentage.ts', 'CHROM_MIN_MARKERS', 'Fewest markers a chromosome needs before its own rate, call rate or intensity median is used.', undefined, 'markers'],
  ['parentage.ts', 'HET_BAND_EXCESS', 'How far the heterozygous band must sit over the homozygous clusters before two genomes are present.'],
  ['parentage.ts', 'BAF_CLAMPED_FLOOR', 'Share of homozygous calls pinned to an exact pole over which the file\'s allele fractions are clamped.'],
  ['parentage.ts', 'GT_FALLBACK_HAPLOID_MAX', 'Heterozygosity, as a share of the parent\'s, under which the genotype fallback calls one genome.', P],
  ['parentage.ts', 'GT_FALLBACK_DIPLOID_MIN', 'The share at or over which it calls two.', P],
  ['parentage.ts', 'CALL_RATE_FLOOR', 'Call rate under which heterozygous calls are not trusted.'],
  ['parentage.ts', 'PANEL_HET_FLOOR', 'Parent heterozygosity under which the marker set is not a polymorphic panel.', P],
  ['parentage.ts', 'HET_BAND_IMPLAUSIBLE', 'Heterozygous band over which it is not measuring zygosity.'],
  ['parentage.ts', 'UNIFORM_CV', 'Spread of the per-chromosome absence rate under which the difference is genome-wide.', P],
  ['parentage.ts', 'BAF_EXTREME_FLOOR', 'Share of a chromosome\'s allele fractions at an extreme needed before any verdict on it.'],
  ['parentage.ts', 'MOSAIC_Z', 'Standard deviations from the rest of the genome before a mosaic is reported.', undefined, 'SD'],
  ['parentage.ts', 'CALL_COLLAPSE', 'Share of the genome\'s median call rate under which a chromosome is not there to genotype.'],
  ['parentage.ts', 'LRR_SHIFT', 'Intensity shift that gives a collapsed chromosome its direction where no null can be formed.', undefined, 'log2'],
  ['parentage.ts', 'COPY_SHIFT_FLOOR', 'Smallest intensity shift that can be a whole-chromosome loss, amplified material.', undefined, 'log2'],
  ['parentage.ts', 'COPY_SHIFT_GAIN_FLOOR', 'The same for a whole-chromosome gain.', undefined, 'log2'],
  ['parentage.ts', 'GAIN_OFF_POLE_MIN', 'Share of the reference off-pole markers a raised stretch must keep before it is an extra copy.'],
  ['parentage.ts', 'ONE_EXTRA_COPY_MAX', 'Intensity step up to which a raised stretch is read as one extra copy.', undefined, 'log2'],
  ['parentage.ts', 'OFF_POLE_MIN_MARKERS', 'Allele fractions a stretch needs before its off-pole share is read.', undefined, 'markers'],
  ['parentage.ts', 'OFF_POLE_MIN_REFERENCE', 'Reference chromosomes needed before an off-pole share is read.'],
  ['parentage.ts', 'COPY_SHIFT_FLOOR_BULK', 'Smallest whole-chromosome intensity shift on unamplified DNA.', undefined, 'log2'],
  ['parentage.ts', 'BULK_SPREAD_MAX', 'Spread between chromosomes over which an array is not read as unamplified.'],
  ['parentage.ts', 'NO_TEMPLATE_HET', 'Heterozygous share at or over which a chromosome holds no template.'],
  ['parentage.ts', 'GAIN_GROUP_MIN', 'Chromosomes over the majority at or over which the heterozygosity test is asked.'],
  ['parentage.ts', 'GAIN_HET_RATIO_MAX', 'Heterozygosity of those chromosomes, as a share of the majority\'s, over which they are not gained.'],
  // partLost.ts
  ['partLost.ts', 'RESTATED_SHARE', 'Share of a chromosome past which its loss pieces are the whole loss.'],
  ['partLost.ts', 'REST_MIN_MARKERS', 'Fewest markers a chromosome, or its rest, needs before it is read as lost whole.', undefined, 'markers'],
  // platform.ts
  ['platform.ts', 'KNOWN_SHARE_MIN', 'Share of a file\'s autosomal markers the founder panel must hold for the file to be read as the proven chip.'],
  // refine.ts
  ['refine.ts', 'REFINE_CALL_FLOOR', 'Call rate under which a window edge is not refined.'],
  ['refine.ts', 'LR_DROP', 'Drop in the likelihood ratio that bounds an edge interval.'],
  ['refine.ts', 'MIN_INTERVAL_MARKERS', 'Narrowest edge interval, in markers.', undefined, 'markers'],
  // relatedness.ts
  ['relatedness.ts', 'OPPOSITE_HOM_MAX', 'Opposite-homozygote rate over which the loaded parent is not this sample\'s parent.', P],
  ['relatedness.ts', 'IBD0_SPREAD_MIN', 'Window spread at or over which the pair is not parent and child.', P],
  ['relatedness.ts', 'IBD0_FLOOR_MAX', 'Window floor under which the spread test is read.', P],
  ['relatedness.ts', 'IBD0_MIN_WINDOWS', 'Fewest windows before the spread test is read.', P],
  // scoreSample.ts
  ['scoreSample.ts', 'BIPARENTAL_AUTOSOMES_MIN', 'Autosomes reading both parents before any chromosome is read as missing one.', P],
  // segments.ts
  ['segments.ts', 'MIN_SEGMENT_MARKERS', 'Fewest called informative markers in a scan window.', undefined, 'markers'],
  ['segments.ts', 'NULL_FLOOR', 'Rate under which a chromosome adds nothing to a scan\'s null.'],
  ['segments.ts', 'SEGMENT_LRT', 'Log-likelihood ratio a segment must clear.'],
  ['segments.ts', 'SEGMENT_LRR_SHIFT', 'Intensity shift a lost segment must carry on top of its call-rate rise.', undefined, 'log2'],
  ['segments.ts', 'SEGMENT_GAIN_SHIFT', 'Intensity shift a gained segment must carry.', undefined, 'log2'],
  ['segments.ts', 'RECIPROCAL_OVERLAP', 'Overlap, as a share of each, at which two segments are one.'],
  ['segments.ts', 'GAIN_NOMINATION_QUANTILE', 'Intensity quantile over which a stretch is put forward as a gain.'],
  ['segments.ts', 'MIN_SEGMENT_SPAN_BP', 'Smallest span a segmental call may claim.', undefined, 'bp'],
  ['segments.ts', 'MIN_GAIN_SPAN_BP', 'Smallest span a gain piece may claim.', undefined, 'bp'],
  // sexChromosome.ts
  ['sexChromosome.ts', 'Y_CN_MIN_PROBES', 'Fewest chrY markers for intensity alone to read the Y.', undefined, 'markers'],
  ['sexChromosome.ts', 'Y_BESIDE_X', 'How far under the array\'s chrX level a chrY may sit and still be a chrY.', undefined, 'log2'],
  // sexing.ts
  ['sexing.ts', 'Y_CALL_MIN', 'Y call rate, as a share of the autosomal call rate, for a Y to be present.'],
  ['sexing.ts', 'Y_LRR_MIN', 'Y intensity against the autosomal median for a Y to be present.', undefined, 'log2'],
  ['sexing.ts', 'Y_MIN_PROBES', 'Fewest Y markers that can say there is no Y.', undefined, 'markers'],
  ['sexing.ts', 'Y_BAF_READABLE', 'Share of homozygous autosomal calls at a pole for allele fractions to read a Y.'],
  ['sexing.ts', 'Y_BAF_PRESENT', 'Share of chrY markers at a pole that is a Y.'],
  ['sexing.ts', 'Y_BAF_ABSENT', 'Share of chrY markers at a pole that is no Y.'],
  ['sexing.ts', 'POLE', 'Allele fraction under which, or over one less which, a marker sits at a pole.'],
  // siblingOrigin.ts
  ['siblingOrigin.ts', 'PHI_WALL', 'False-heterozygote share over which no sister-referenced call is made.', RUN],
  ['siblingOrigin.ts', 'SISTER_CALL_POSTERIOR', 'Posterior a sister-referenced call must reach.', RUN],
  ['siblingOrigin.ts', 'SISTER_MIN_MARKERS', 'Fewest informative markers before a sister-referenced region is scored.', RUN, 'markers'],
  ['siblingOrigin.ts', 'MIN_SIBLINGS', 'Fewest other cells of the run before a sister-referenced call.', RUN],
  // stage.ts
  ['stage.ts', 'ALLELE_PLOIDY_MIN_STATED', 'Chromosomes stating a copy state before the alleles say how many genomes are present.', P],
  ['stage.ts', 'HET_SPREAD_LOCALISED', 'Heterozygosity over which the array is not one genome.'],
  ['stage.ts', 'MAX_DIPLOID_HET', 'Heterozygosity over which the stage is not a diploid cell.'],
  ['stage.ts', 'BULK_MAX_BAF_SD', 'Allele-fraction spread at heterozygous calls over which a sample is not unamplified DNA.'],
  ['stage.ts', 'HAPLOID_MAX_HET', 'Heterozygosity under which a cell carries one genome.'],
  ['stage.ts', 'BAND_DIPLOID_CERTAIN', 'Heterozygous band at or over which two genomes are present whatever the calls say.'],
  ['stage.ts', 'PANEL_ANCHOR_CALL_RATE', 'Call rate a parent array needs before it anchors the panel\'s heterozygosity.', P],
  ['stage.ts', 'PANEL_SCALE_MAX', 'Most the panel scale may move the stage bars.', P],
  // timing.ts
  ['timing.ts', 'SAME_BREAK_BP', 'Breakpoints closer than this are one, when cells of one embryo are compared.', RUN, 'bp'],
  // untransmitted.ts
  ['untransmitted.ts', 'TRISOMY_BAND_HALF_WIDTH', 'Half-width around a trisomy band centre.', P],
  ['untransmitted.ts', 'MIN_MARKERS_MECHANISM', 'Fewest markers before the mechanism of an extra copy is read.', P, 'markers'],
]

/** Named numeric constants of the scoring code that decide nothing on their own, and why. */
export const NOT_BARS: Readonly<Record<string, string>> = {
  'abnormalities.ts:UPD_NO_STRETCH_RATE': 'a published rate printed beside a negative result',
  'abnormalities.ts:CLUSTERED_DROPOUT_DEPLETION': 'a measured artefact level that LOH_DEPLETION is set clear of',
  'abnormalities.ts:CLUSTERED_DROPOUT_RATE': 'the share of windows that artefact reaches, a measurement',
  'abnormalities.ts:LOH_SEGMENT_MARKERS': 'the size of a scan window',
  'abnormalities.ts:TRIPLOID_LOADED_TWO_THIRDS': 'where an allele sits under one hypothesis, a model value',
  'abnormalities.ts:TRIPLOID_LOADED_ONE_THIRD': 'the same, under the other hypothesis',
  'answerKeys.ts:PANEL_DROPIN_SLOPE': 'a calibration slope',
  'answerKeys.ts:PANEL_DROPIN_INTERCEPT': 'a calibration intercept',
  'dosageOrigin.ts:WINDOW_LO': 'the edge of the window of markers read, not a bar on a result',
  'dosageOrigin.ts:WINDOW_HI': 'the other edge of that window',
  'features.ts:DEFAULT_PERMUTATIONS': 'how many permutations are drawn',
  'features.ts:HIT_NAMES': 'how many names a list prints',
  'features.ts:MAX_FEATURE_SPAN': 'where a name search stops',
  'gainOrigin.ts:EXPECTED_SEPARATION': 'where the bands sit, printed for context and never compared against',
  'inferredArray.ts:MARK_BYTES': 'how much of a file is read for its banner',
  'informationRatio.ts:INFORMATIVE_MAX': 'not used to decide anything',
  'ingest.ts:DEGREE_CALL_FLOOR': 'the relatedness degree call of the germline module, not run by this one',
  'ingest.ts:DROPOUT_MISMATCH_CEILING': 'the same germline degree call',
  'ingest.ts:SIBLING_DISPERSION': 'the same germline degree call',
  'motherFromRun.ts:WINDOW_BP': 'the size of a window',
  'motherFromRun.ts:P_MAX': 'a cap on a probability, a numerical guard',
  'oneParentOrigin.ts:DEFAULT_Q': 'a model allele frequency',
  'oneParentOrigin.ts:GENOTYPE_ERROR': 'a model error rate',
  'oneParentOrigin.ts:DROP_IN': 'a model drop-in rate',
  'oneParentOrigin.ts:VALIDATION_UNITS': 'a count printed with the accuracy',
  'oneParentOrigin.ts:SYSTEMATIC_ERROR_BOUND': 'a measured error printed with the accuracy',
  'oneParentOrigin.ts:HALF_EVIDENCE': 'the shape of the confidence curve',
  'oneParentOrigin.ts:BULK_HETEROZYGOSITY': 'a measured heterozygosity used as a reference level',
  'otherParentLoss.ts:BLOCK_MARKERS': 'the size of a block',
  'otherParentLoss.ts:CHI2_MEDIAN': 'a mathematical constant',
  'parentage.ts:OFF_POLE_LO': 'the edge of the window an allele fraction is off the poles in',
  'parentage.ts:OFF_POLE_HI': 'the other edge of that window',
  'otherParentLoss.ts:EDGE_DROP': 'LR_DROP, the edge bound of the window scan, applied on this scan\'s scale',
  'refine.ts:SEARCH_RADIUS': 'how far a search looks',
  'segments.ts:GAIN_SHIFT_FRACTION': 'a factor; the bar it makes is SEGMENT_GAIN_SHIFT',
  'sexing.ts:FLOOR_STRIDE': 'how many markers are sampled',
  'sexing.ts:LRR_STRIDE': 'how many markers are sampled',
  'siblingOrigin.ts:WRONG_HOMOZYGOTE_RATE': 'a model error rate',
  'untransmitted.ts:BAND_BPH': 'a band centre under one hypothesis, a model value',
  'untransmitted.ts:BAND_SPH_LOW': 'a band centre under the other hypothesis',
  'untransmitted.ts:BAND_SPH_HIGH': 'a band centre under the other hypothesis',
}

export const BARS: readonly Bar[] = SPECS.map(([file, name, says, need, unit]) => {
  const value = MODULES[file]?.[name]
  if (typeof value !== 'number') throw new Error(`bars.ts: ${file} exports no number ${name}`)
  return { file, name, value, says, ...(need ? { need } : {}), ...(unit ? { unit } : {}) }
})

/** A bar's value as every surface prints it. */
export function barValue(b: Pick<Bar, 'value' | 'unit'>): string {
  if (b.unit === 'bp') return `${+(b.value / 1e6).toPrecision(4)} Mb`
  const v = String(+b.value.toPrecision(4))
  return b.unit === 'x' ? `${v}x` : b.unit ? `${v} ${b.unit}` : v
}

export interface Setup {
  /** Parental arrays loaded: 0, 1 or 2. */
  parents: 0 | 1 | 2
  /** Other arrays in the run, read beside this one. */
  run: boolean
  /** The feature comparison was run. */
  comparison: boolean
}

const fits = (n: Need | undefined, s: Setup): boolean => !n || (
  (!n.parent || (n.parent === 'any' ? s.parents > 0 : n.parent === 'one' ? s.parents === 1 : s.parents === 2))
  && (!n.run || s.run) && (!n.comparison || s.comparison))

/** The bars a sample in this setup was decided by, in registry order. */
export const barsApplied = (s: Setup): Bar[] => BARS.filter((b) => fits(b.need, s))

/** One bar by name, applied in this setup; throws where it was not, so no text can quote it. */
export function barIn(applied: readonly Bar[], name: string, file?: string): Bar {
  const b = applied.find((x) => x.name === name && (!file || x.file === file))
  if (!b) throw new Error(`bar ${name} was not applied to this run`)
  return b
}

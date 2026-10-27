/**
 * Ploidy from markers where a parent's contribution is obligate, not from the BAF band.
 *
 * WHY THIS EXISTS. The BAF-band gate in `productTriage.ts` separates a haploid meiotic product
 * from a bulk diploid adult, which is the only distinction it was ever measured on, and it does
 * that job: on the 46 arrays this tool has run, every pronucleus reads haploid and every bulk
 * adult array reads diploid. It does NOT generalise to post-zygotic material. An external review
 * measured it on 120 non-haploid arrays spanning cleavage blastomeres, trophectoderm biopsies and
 * ESC lines: 91 of the 120 fall inside the haploid BAF-band range, and 22 blastomeres sit below
 * the 10% diploid-exclusion floor entirely. The cause is whole-genome amplification, which both
 * widens the true heterozygote band and fills the homozygote band, so the two become one
 * distribution. Refining to a BAF core fraction fails the same way.
 *
 * So the band gate is not broken for what it does today. It is unfit for the stage extension, and
 * this module is the replacement that the extension needs.
 *
 * THE STATISTIC. At a marker where one parent is homozygous, a cell carrying that parent's genome
 * cannot be heterozygous unless it also carries a different allele from somewhere else. A haploid
 * meiotic product has nowhere else. So the heterozygous fraction AT PARENT-HOMOZYGOUS MARKERS is
 * an error rate in a haploid and a real signal in a biparental diploid. Measured by the same
 * review at obligate-heterozygous markers (both parents homozygous and opposite):
 *
 *     stage                     AB fraction
 *     paternal pronucleus            0.0734
 *     maternal pronucleus            0.0705
 *     cleavage blastomere            0.5244
 *     single ESC (WGA)               0.7818
 *     trophectoderm biopsy           0.8950
 *     ESC line (bulk)                0.9868
 *
 * Twelve haploid arrays span 0.047 to 0.100 and the lowest diploid is 0.4213: a gap of four times
 * the whole haploid range, where the BAF band has no gap at all. Boundaries are placed at 0.20 and
 * 0.45, well inside it, and 68 of 68 callable arrays classified correctly with none misclassified.
 *
 * ONE PARENT IS ENOUGH, WITH A CAVEAT. The measured figures above use both parents, which this
 * tool often does not have. With one parent the informative set is that parent's homozygous
 * markers, and the diploid signal is diluted: a biparental diploid only reads heterozygous where
 * the OTHER parent differs, which is roughly that parent's heterozygosity plus the minor-allele
 * mismatch, not every marker. The direction and the ordering survive; the boundary does not
 * transfer unchanged. `oneParent` therefore carries its own, wider band and refuses more often.
 *
 * AND THE REFERENCE PARENT MUST BE BULK. This is a hard precondition, measured externally rather
 * than argued, on GSE19247, where the truth is known by construction: single sperm are haploid,
 * Day 3 blastomeres are biparental. Scored against a SINGLE-CELL reference parent the statistic
 * does not separate them at all.
 *
 *     single sperm, haploid        n=23   median 0.0548   range 0.0021 - 0.5132
 *     Day 3 blastomere, biparental n=28   median 0.0962   range 0.0177 - 0.3004
 *
 * Four of the 23 haploid sperm read 27% to 51% heterozygous, which one genome cannot do, and they
 * carry the HIGHEST call rates in the set, so no call-rate gate removes them. The mechanism is the
 * reference: dropout turns its true heterozygous markers into apparent homozygotes, so its
 * "homozygous" set becomes most of the genome and conditioning on it stops doing anything.
 * Measured, obligate-het and plain genome-wide heterozygosity agreed to within 0.01 on all 51
 * arrays. A reference that is itself amplified single-cell material carries no information here.
 * See audit/BREAKPOINTS-AUDIT.txt section E2.
 *
 * PER CHROMOSOME, NOT PER ARRAY. This series is chromosome-loss experiments and the review found
 * 7 arrays that are genome-wide partial: 6 to 12 uniparental chromosomes inside an otherwise
 * biparental cell. A per-array verdict would call those diploid and hide exactly the event the
 * lab is looking for. 123 of 1,513 diploid chromosome-observations fall below the haploid
 * boundary and 39 of 40 with technical replicates are concordant, so they are real.
 */
import type { AB } from './informativity.ts'

/** Below this, the cell carries one parent's genome at this locus and nothing else. */
export const HAPLOID_MAX = 0.20

/** At or above this, two different parental contributions are present. */
export const DIPLOID_MIN = 0.45

/**
 * With one parent only, the diploid signal is diluted to the markers where the unseen parent
 * differs. These were previously 0.12 and 0.30, set by scaling the two-parent figures and
 * described as provisional. MEASURED NOW, and the dilution is far larger than that scaling
 * assumed, so the old values called every biparental child of a single genotyped parent
 * uniparental.
 *
 * Against a BULK reference parent, genome-wide at his 639,543 autosomal homozygous markers, on
 * nine of his children found across four experiments:
 *
 *     his one-genome products    0.0428  0.0514  0.0570  0.0587
 *     his biparental children    0.0894  0.0973  and above
 *
 * Both classes sat under the old 0.12 bound. The real separation is an order of magnitude lower
 * than assumed, and the boundaries below bracket it with the gap left uncalled.
 *
 * THE PRECONDITION IS UNCHANGED AND STILL BINDING: this is measured against a BULK parent, and
 * against a single-cell reference the same statistic does not separate at all. See the section
 * above and audit/regions/FINDINGS.txt.
 */
export const ONE_PARENT_HAPLOID_MAX = 0.065
export const ONE_PARENT_DIPLOID_MIN = 0.085

/** Fewest called informative markers before a fraction means anything. The review's per-chromosome
 *  informative minimum was 304 with a median of 918 to 1,362, so this refuses the thin tail. */
export const MIN_INFORMATIVE = 200

/**
 * Whether a chromosome is being genotyped well enough for its heterozygous fraction to mean
 * anything.
 *
 * AN INFORMATIVE COUNT CANNOT TELL A THIN CHROMOSOME FROM A SMALL ONE, which is why this is a call
 * rate and not a count. Measured on one laboratory series, 176 chromosome-observations over 8
 * samples: one chromosome carries 2,034 called markers of the 11,655 the array holds for it, a
 * call rate of 0.175, and reports 55.1 percent heterozygous against a median of 0.07 for its
 * siblings. It clears the absolute floor of 200 with 1,605 informative markers, and its share of
 * informative markers among those called is 0.789, which is the median for the series. Nothing
 * about the informative count is unusual. The chromosome was barely genotyped.
 *
 * The call rate separates it with room to spare: 0.175 against a median of 0.864, a fifth
 * percentile of 0.807, and nothing at all between 0.175 and 0.813.
 *
 * A first attempt used the informative count as a fraction of the sample's median autosome, and a
 * regression on the public corpus refused chr22 of one cell for carrying 17 percent of what the
 * median autosome carries. chr22 is simply a small chromosome. A count compared across
 * chromosomes of different sizes cannot work, and the call rate needs no such comparison.
 *
 * The threshold is `CALL_COLLAPSE`, which the per-chromosome mis-clustering gate in parentage.ts
 * already uses for the same question, so this adds no constant.
 */
export const informativeEnough = (
  informative: number,
  callRate: number,
  floor = MIN_INFORMATIVE,
  collapse = 0.60,
): boolean => informative >= floor && (!Number.isFinite(callRate) || callRate >= collapse)

/**
 * One parental contribution against two, from the markers where the loaded parent is HETEROZYGOUS.
 *
 * `oneParentOrigin.ts` excludes that class, correctly for the three-way copy question: at such a
 * marker both loss hypotheses predict the same thing, so they cannot be told apart there. For the
 * narrower question of whether a second contribution is present, two hypotheses suffice and the
 * excluded class carries far more signal:
 *
 *   parent homozygous    two contributions  het = q(1 - ado) + d,  q = 0.1217
 *   parent heterozygous  two contributions  het = 0.5(1 - ado) + d
 *   either class         one contribution   het = d
 *
 * Four times the prediction, and it does not depend on the drop-in parameter, which is measured on
 * bulk material and reads 1.66 times higher on uniparental amplified arrays.
 *
 * WHY A LIKELIHOOD RATIO AND NOT A BOUNDARY. The two rate distributions OVERLAP on real material:
 * over 1,021 chromosome-observations of GSE148488 against the authors' published karyotype, the
 * one-contribution class runs to 0.4320 and the two-contribution class down to 0.1334. A pair of
 * boundaries cannot separate them and a ratio can, because it reads the marker count as well as
 * the rate. Measured, same rule applied to both classes at the shipped posterior:
 *
 *   class                  one contribution   two contributions   overall   answered
 *   parent-homozygous        50 of 102 0.490   874 of 917 0.953    0.9068     99.8%
 *   parent-heterozygous      80 of 103 0.777   907 of 918 0.988    0.9667    100.0%
 *
 * Better on both, and the calls on chromosomes the hand karyotype reads as normal fall from 43 in
 * 918 to 11. Against this project's objective, detections less twice the false calls at k = 2, the
 * excluded class scores 965 against 838 while answering more chromosomes.
 *
 * WHAT IT DOES NOT SAY, AND WHY IT IS NOT THE VERDICT'S CALL. Whether the second contribution is
 * the OTHER parent's. A sample carrying two haplotypes of ONE parent is heterozygous wherever a
 * crossover separated them, so it reads high here while no second parent contributed. Measured on
 * real arrays of a series with a genotyped egg donor and no father, by the opposite-homozygote rate
 * against her:
 *
 *   opp-hom   this class   the shipped class   carries
 *     0.55%   biparental   uncalled            her genome only
 *     0.16%   biparental   uniparental         her genome only
 *     2.35%   biparental   biparental          two parental genomes
 *
 * 0.3408 and 0.3734 against 0.3974, so the classes are separated by under 0.03 and no boundary
 * divides them. The shipped call is right about the second row and declines the first; this one
 * asserts a paternal contribution that is not there. So the verdict keeps the parent-homozygous
 * call despite the table above, and this function answers only whether a SECOND HAPLOTYPE is
 * present. `obligateHet.check.ts` pins the overlap so the substitution cannot be made silently.
 */
export function hetCallOneParent(
  het: number,
  informative: number,
  ado: number,
  dropIn = 0.0435,
  posterior = 0.95,
  minMarkers = 50,
): HetCall {
  const fraction = informative ? het / informative : NaN
  // `onlyA` and `onlyB` are not defined on this class: where the loaded parent is heterozygous
  // there is no "that parent's allele alone" to count, because it carries both. They are NaN
  // rather than zero so a reader cannot mistake an absent quantity for a measured one, and
  // `provisional` is true because one parent defined the informative set.
  const base = { fraction, informative, onlyA: NaN, onlyB: NaN, provisional: true }
  if (informative < minMarkers) {
    return { ...base, ploidy: 'uncalled',
      why: `${informative} markers where the loaded parent is heterozygous is under the `
        + `${minMarkers} this needs` }
  }
  const a = Number.isFinite(ado) ? Math.min(0.95, Math.max(0, ado)) : 0.3
  const pTwo = Math.min(0.999, 0.5 * (1 - a) + dropIn)
  const pOne = Math.max(1e-9, dropIn)
  const lp = (p: number) => het * Math.log(p) + (informative - het) * Math.log(1 - p)
  const l2 = lp(pTwo)
  const l1 = lp(pOne)
  const m = Math.max(l1, l2)
  const post = Math.exp(l2 - m) / (Math.exp(l2 - m) + Math.exp(l1 - m))
  const pct2 = (x: number) => `${(100 * x).toFixed(1)}%`
  if (post >= posterior) {
    return { ...base, ploidy: 'biparental',
      why: `${pct2(fraction)} heterozygous over ${informative.toLocaleString()} markers where the `
        + `loaded parent is heterozygous. Two contributions predict ${pct2(pTwo)} at this `
        + `sample's dropout and one predicts ${pct2(pOne)}, which is drop-in alone` }
  }
  if (1 - post >= posterior) {
    return { ...base, ploidy: 'uniparental',
      why: `${pct2(fraction)} heterozygous over ${informative.toLocaleString()} markers where the `
        + `loaded parent is heterozygous, against ${pct2(pTwo)} expected if a second parent had `
        + 'contributed here' }
  }
  return { ...base, ploidy: 'uncalled',
    why: `${pct2(fraction)} heterozygous sits between the two predictions, ${pct2(pOne)} and `
      + `${pct2(pTwo)}, at a posterior of ${Math.max(post, 1 - post).toFixed(2)}` }
}

export type Ploidy2 = 'uniparental' | 'biparental' | 'uncalled'

export interface HetTally {
  /** Markers where the reference parent(s) made this marker informative and the sample called. */
  informative: number
  /** Of those, how many the sample called heterozygous. */
  het: number
  /**
   * Of those, how many read the FIRST parent's allele alone, and how many the second's alone.
   * Two parents only: with one parent loaded the other allele has no owner, so both stay at zero.
   */
  onlyA: number
  onlyB: number
}

export const emptyHet = (): HetTally => ({ informative: 0, het: 0, onlyA: 0, onlyB: 0 })

/**
 * Share of informative markers reading ONE parent's allele alone before that parent is named the
 * survivor, which names the other parent's copy as the one lost.
 *
 * Measured on GSE148488 against the hand karyotypes of Zuccaro 2020 Table S3, 31 arrays with both
 * parents, 595 autosomes: on the 105 chromosomes read as lost the surviving parent's share was
 * never under 0.735, and on the 490 intact chromosomes neither parent's share reached 0.627. The
 * floor sits in that gap. Dropout cannot cross it, because dropout removes a heterozygote's alleles
 * at random rather than one parent's allele at every marker of a chromosome.
 *
 * Two BULK parental arrays are a precondition, for the same reason as the boundaries above: a
 * reconstructed or single-cell reference is homozygous at markers where the parent is not, and
 * every such marker reads as the other parent's absence. audit/zuccaro-table-s3.ts reproduces the
 * figures.
 */
export const SURVIVOR_MIN = 0.70

/**
 * One parent's share at which a chromosome reading BETWEEN the boundaries is reported as a mixed
 * signal rather than passed over in silence. Measured on GSE148488 against Table S3: of the 30
 * chromosomes with 200 or more informative markers and a heterozygous share between 0.20 and
 * 0.45, 13 carry one parent's allele alone at 0.55 or more, and of those 8 are hand-called
 * monosomies, 4 hand-called normal and 1 something else. Below it the two classes mix evenly. The
 * row names no parent: it carries the counts and says it was not called.
 */
export const MIXED_SHARE = 0.55

/**
 * One marker, with one parent known.
 *
 * Informative means that parent is homozygous there. A no-call in the sample contributes nothing
 * either way rather than counting as homozygous, which would drag the fraction toward uniparental
 * on exactly the degraded arrays that need the opposite.
 */
export function addOneParent(parent: AB, sample: AB, t: HetTally): void {
  if (parent !== 'AA' && parent !== 'BB') return
  if (sample === 'NC') return
  t.informative += 1
  if (sample === 'AB') t.het += 1
}

/**
 * One marker, with both parents known. Informative means both are homozygous AND opposite, which
 * is the configuration a biparental diploid must read heterozygous at.
 */
export function addTwoParent(a: AB, b: AB, sample: AB, t: HetTally): void {
  const oppositeHom = (a === 'AA' && b === 'BB') || (a === 'BB' && b === 'AA')
  if (!oppositeHom) return
  if (sample === 'NC') return
  t.informative += 1
  if (sample === 'AB') t.het += 1
  else if (sample === a) t.onlyA += 1
  else if (sample === b) t.onlyB += 1
}

export interface HetCall {
  ploidy: Ploidy2
  fraction: number
  informative: number
  /** Share of informative markers reading the first parent's allele alone, and the second's. */
  onlyA: number
  onlyB: number
  /**
   * Which parent's allele is all that remains, on a uniparental chromosome with two parents
   * loaded and that share at or over SURVIVOR_MIN. Absent otherwise: a uniparental reading with no
   * survivor is a chromosome the arrays are not measuring, not a loss.
   */
  survivor?: 'a' | 'b'
  /** True when only one parent defined the informative set, so the boundaries are the wider,
   *  provisional pair rather than the measured ones. */
  provisional: boolean
  why: string
}

export function hetCall(t: HetTally, parents: 1 | 2): HetCall {
  const fraction = t.informative ? t.het / t.informative : NaN
  const onlyA = t.informative ? (t.onlyA ?? 0) / t.informative : NaN
  const onlyB = t.informative ? (t.onlyB ?? 0) / t.informative : NaN
  const provisional = parents === 1
  const lo = provisional ? ONE_PARENT_HAPLOID_MAX : HAPLOID_MAX
  const hi = provisional ? ONE_PARENT_DIPLOID_MIN : DIPLOID_MIN
  if (t.informative < MIN_INFORMATIVE) {
    return {
      ploidy: 'uncalled',
      fraction,
      informative: t.informative,
      onlyA,
      onlyB,
      provisional,
      why: `${t.informative} informative markers is under the ${MIN_INFORMATIVE} this needs`,
    }
  }
  if (fraction < lo) {
    const survivor = parents === 2 && onlyA >= SURVIVOR_MIN ? 'a' as const
      : parents === 2 && onlyB >= SURVIVOR_MIN ? 'b' as const : undefined
    return {
      ploidy: 'uniparental',
      fraction,
      informative: t.informative,
      onlyA,
      onlyB,
      survivor,
      provisional,
      why: `${(fraction * 100).toFixed(1)}% heterozygous where a parent is homozygous, under `
        + `${(lo * 100).toFixed(0)}%: one parent's genome and nothing else`,
    }
  }
  if (fraction >= hi) {
    return {
      ploidy: 'biparental',
      fraction,
      informative: t.informative,
      onlyA,
      onlyB,
      provisional,
      why: `${(fraction * 100).toFixed(1)}% heterozygous where a parent is homozygous, at or over `
        + `${(hi * 100).toFixed(0)}%: a second parental contribution is present`,
    }
  }
  return {
    ploidy: 'uncalled',
    fraction,
    informative: t.informative,
    onlyA,
    onlyB,
    provisional,
    why: `${(fraction * 100).toFixed(1)}% sits between the ${(lo * 100).toFixed(0)}% and `
      + `${(hi * 100).toFixed(0)}% boundaries, which is neither`,
  }
}


/**
 * The rate at which a sample that has lost ONE parent's copy reads an opposite homozygote against
 * that parent, computed from the two parental genotypes alone.
 *
 * At a marker where the lost parent is homozygous, the sample carries only the other parent's
 * transmitted allele, and that is the OPPOSITE allele exactly when the other parent transmitted
 * it. A heterozygous other parent transmits either with probability one half; a homozygous one
 * transmits what it has. So the rate is a sum over known genotype pairs with nothing assumed about
 * the allele-frequency spectrum.
 *
 * This is what a real event looks like, and it is needed to ask whether a scan could detect one at
 * all. Measured on GSE148488 it comes to 0.1244 against the mother and 0.1153 against the father
 * over roughly 650,000 markers each, which agrees with the 11.6 to 12.0 percent that dissected
 * paternal pronuclei read against the egg donors, a measurement of the same quantity on real
 * single-parent material.
 */
export function lossSignalRate(
  lost: Map<string, string>,
  other: Map<string, string>,
): number {
  let n = 0
  let opp = 0
  for (const [id, lg] of lost) {
    if (lg !== 'AA' && lg !== 'BB') continue
    const og = other.get(id)
    if (!og || og === 'NC') continue
    const otherAllele = lg === 'AA' ? 'B' : 'A'
    n += 1
    if (og === 'AB') opp += 0.5
    else if (og[0] === otherAllele) opp += 1
  }
  return n ? opp / n : NaN
}

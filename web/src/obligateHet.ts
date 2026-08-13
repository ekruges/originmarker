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

export type Ploidy2 = 'uniparental' | 'biparental' | 'uncalled'

export interface HetTally {
  /** Markers where the reference parent(s) made this marker informative and the sample called. */
  informative: number
  /** Of those, how many the sample called heterozygous. */
  het: number
}

export const emptyHet = (): HetTally => ({ informative: 0, het: 0 })

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
}

export interface HetCall {
  ploidy: Ploidy2
  fraction: number
  informative: number
  /** True when only one parent defined the informative set, so the boundaries are the wider,
   *  provisional pair rather than the measured ones. */
  provisional: boolean
  why: string
}

export function hetCall(t: HetTally, parents: 1 | 2): HetCall {
  const fraction = t.informative ? t.het / t.informative : NaN
  const provisional = parents === 1
  const lo = provisional ? ONE_PARENT_HAPLOID_MAX : HAPLOID_MAX
  const hi = provisional ? ONE_PARENT_DIPLOID_MIN : DIPLOID_MIN
  if (t.informative < MIN_INFORMATIVE) {
    return {
      ploidy: 'uncalled',
      fraction,
      informative: t.informative,
      provisional,
      why: `${t.informative} informative markers is under the ${MIN_INFORMATIVE} this needs`,
    }
  }
  if (fraction < lo) {
    return {
      ploidy: 'uniparental',
      fraction,
      informative: t.informative,
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
      provisional,
      why: `${(fraction * 100).toFixed(1)}% heterozygous where a parent is homozygous, at or over `
        + `${(hi * 100).toFixed(0)}%: a second parental contribution is present`,
    }
  }
  return {
    ploidy: 'uncalled',
    fraction,
    informative: t.informative,
    provisional,
    why: `${(fraction * 100).toFixed(1)}% sits between the ${(lo * 100).toFixed(0)}% and `
      + `${(hi * 100).toFixed(0)}% boundaries, which is neither`,
  }
}

/**
 * What kind of material this array is, inferred from the array itself.
 *
 * THE PHYSICAL CAUSE, which is what everything here is grounded in. Allele dropout is a sampling
 * failure during amplification, and its rate is set by HOW MANY TEMPLATE MOLECULES the reaction
 * started from. A diploid locus in bulk genomic DNA is present in millions of copies, so losing
 * every copy of one allele is impossible. The same locus in a single cell is present in exactly
 * two molecules, one per homologue, and a heterozygote survives only if BOTH amplify: if the one
 * copy of an allele fails to prime in the first cycles, that allele is absent from everything
 * downstream and the marker reads homozygous. Nothing about the genome changed; the reaction lost
 * it. So dropout is a function of starting template count, and that count is a property of the
 * developmental stage the sample was taken at:
 *
 *   bulk genomic DNA        > 10^6 cells      ~2 x 10^6 templates per locus    dropout ~0.013
 *   ES cell line, bulk      cultured colony   same order                      ~0.013
 *   trophectoderm biopsy    5-10 cells        10-20 templates                 ~0.050
 *   single ES cell, WGA     1 cell            2 templates                     ~0.199
 *   cleavage blastomere     1 cell            2 templates                     ~0.308
 *
 * The last two are both one cell and both start from two molecules, yet differ nearly two-fold.
 * That is not template count but chromatin: a cleavage-stage blastomere is in a rapid cell cycle
 * with a decondensed, replication-active genome and carries more single-stranded and partially
 * replicated template, which primes less reliably. A cultured ES cell is a more ordinary
 * interphase nucleus. So the ordering is template count first, chromatin state second.
 *
 * HAPLOID MATERIAL IS A DIFFERENT AXIS ENTIRELY and must not be read as heavy dropout. A polar
 * body, a pronucleus or a sperm carries ONE genome, so it is homozygous everywhere by construction
 * and its heterozygous calls are all error. Judged by heterozygosity alone it would look like a
 * catastrophically dropped-out diploid, and treating it as one would put the dropout parameter at
 * a nonsensical value. It is separated first, by how far below any diploid stage it sits.
 *
 * WHY THIS IS INFERRED RATHER THAN ASKED. The tool's contract is that a file is dropped in and
 * nothing else is required. Stage is knowable from the array, so asking would be asking the user
 * to restate something already in front of the program.
 */
import type { SampleProfile } from './ingest.ts'

export type Stage =
  | 'bulk'
  | 'trophectoderm'
  | 'single-cell'
  | 'blastomere'
  | 'haploid'
  | 'unknown'

export interface StageCall {
  stage: Stage
  /** Expected allele dropout for this stage, used to parameterise every downstream likelihood. */
  dropout: number
  /** Order-of-magnitude template molecules per locus the amplification started from. */
  templates: string
  /** Informative markers a directional call needs at this stage. */
  markerFloor: number
  why: string
}

/**
 * Heterozygosity of a diploid genome on this platform's marker set, measured on bulk gDNA.
 * Every dropout estimate below is a shortfall against it.
 */
export const BULK_HETEROZYGOSITY = 0.168

/**
 * Below this, the sample carries one genome rather than a heavily dropped-out two.
 *
 * Measured: haploid meiotic products run 0.002 to 0.10 genome-wide heterozygosity, and the lowest
 * diploid material sits far above. The boundary is placed in the gap rather than at either edge.
 */
export const HAPLOID_MAX_HET = 0.105

/**
 * Stage boundaries in heterozygosity, derived from the dropout rates above.
 *
 * A diploid reading h has lost the fraction 1 - h/H of its heterozygotes to dropout, so each
 * dropout rate implies a heterozygosity and the boundaries sit between them rather than on them.
 */
const BOUNDS: { stage: Stage, minHet: number, dropout: number, templates: string, floor: number }[] = [
  { stage: 'bulk', minHet: 0.158, dropout: 0.013, templates: '~10^6', floor: 100 },
  { stage: 'trophectoderm', minHet: 0.145, dropout: 0.050, templates: '10-20', floor: 100 },
  { stage: 'single-cell', minHet: 0.126, dropout: 0.199, templates: '2', floor: 200 },
  { stage: 'blastomere', minHet: HAPLOID_MAX_HET, dropout: 0.308, templates: '2', floor: 200 },
]

/**
 * Call the stage from a sample's own profile.
 *
 * `callRate` is used only as a guard: an array that failed outright can show any heterozygosity at
 * all, and calling a stage from it would attach a confident dropout parameter to noise.
 */
export function inferStage(profile: Pick<SampleProfile, 'hetRate' | 'callRate'>): StageCall {
  const h = profile.hetRate
  const call = profile.callRate

  if (!Number.isFinite(h) || !Number.isFinite(call) || call < 0.40) {
    return {
      stage: 'unknown',
      dropout: 0.308,
      templates: 'unknown',
      markerFloor: 200,
      why: `call rate ${Number.isFinite(call) ? (call * 100).toFixed(1) : '?'}% is too low for the `
        + 'stage to be read from this array; the most conservative dropout is assumed',
    }
  }

  if (h <= HAPLOID_MAX_HET) {
    return {
      stage: 'haploid',
      // One genome cannot be heterozygous, so there is no heterozygote to drop. The residual rate
      // is genotyping error, not dropout, and the dosage channel is what such a sample is read by.
      dropout: 0.02,
      templates: '1 genome',
      markerFloor: 200,
      why: `${(h * 100).toFixed(1)}% heterozygous is below the ${(HAPLOID_MAX_HET * 100).toFixed(1)}% `
        + 'a diploid reaches at any stage, so this carries one genome: a polar body, a pronucleus '
        + 'or a sperm. Its heterozygous calls are error rather than biology',
    }
  }

  const hit = BOUNDS.find((b) => h >= b.minHet) ?? BOUNDS[BOUNDS.length - 1]
  const implied = Math.min(0.6, Math.max(0.005, 1 - h / BULK_HETEROZYGOSITY))
  return {
    stage: hit.stage,
    // The sample's own implied rate, floored at the stage's expectation so that an unusually clean
    // array of a lossy stage is not credited with bulk-quality amplification.
    dropout: Math.max(implied, hit.dropout * 0.6),
    templates: hit.templates,
    markerFloor: hit.floor,
    why: `${(h * 100).toFixed(1)}% heterozygous against ${(BULK_HETEROZYGOSITY * 100).toFixed(1)}% `
      + `for bulk DNA implies ${(implied * 100).toFixed(1)}% allele dropout, which is `
      + `${hit.stage} material amplified from ${hit.templates} template copies per locus`,
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

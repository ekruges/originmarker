/**
 * What the defect display is told, separate from how it is drawn.
 *
 * Kept out of the .tsx for two reasons: node cannot strip JSX, so a component file cannot carry a
 * self-check, and this is the layer where getting it wrong matters. Showing a confident parent for
 * a region that never had one is the most damaging thing this UI could do, so the mapping that
 * decides it is testable on its own.
 */
import type { Segment } from './segments.ts'
import type { GainAnnotation } from './parentage.ts'
import { segmentCoords } from './segments.ts'
import { locus } from './stage.ts'
import type { StageCall } from './stage.ts'

/** What is known about one defect, gathered from whichever channels produced it. */
export interface Defect {
  chrom: string
  startBp: number
  endBp: number
  kind: 'copy-loss' | 'copy-gain' | 'parental-absence' | 'whole-chromosome'
  /** The region in genome-browser form, chr6:39,302-294,904. */
  locus: string
  /** Refined interval text where the breakpoint was localised, e.g. '+/- 61 kb'. */
  interval?: string
  /** Origin as the tool determined it, with the reason. */
  origin: 'paternal' | 'maternal' | 'unclear'
  why: string
  /** How the origin was reached, which decides how much weight it carries. */
  basis?: 'two-parent' | 'one-parent' | 'sibling' | 'homologue' | 'dosage'
  /** Which channel produced it. Genotypes where they exist, dosage where they have collapsed. */
  channel?: 'genotype' | 'dosage'
  /** Dosages at the extreme the loaded parent cannot produce. The dosage channel's own evidence,
   *  the counterpart of `exclusive` in the genotype channel. */
  excludedDosage?: number
  informative?: number
  posterior?: number
  /** Markers carrying an allele the loaded parent does not have. This is the Mendelian evidence
   *  itself rather than a summary of it: a parent who is AA has no B to give, and dropout removes
   *  alleles without inventing one, so a non-zero count cannot come from amplification loss. */
  exclusive?: number
  /** The material this was read from, carried so a chip can state it without a second lookup. */
  stage?: string
  /** Dropout the call was parameterised with, and how that figure was arrived at. */
  dropout?: number
  dropoutBasis?: string
}

/** Build the display list from a run's segments and whatever origin annotations exist for them. */
export function defectsFrom(
  segments: readonly Segment[],
  gains: readonly GainAnnotation[],
  losses: readonly GainAnnotation[],
  /** Single-parent calls, keyed by the same `where` string. The loaded parent's role decides
   *  whether 'known-parent-lost' reads paternal or maternal. */
  oneParent: readonly {
    where: string, verdict: string, posterior: number, markers: number, exclusive: number,
    why: string,
  }[] = [],
  loadedParent: 'paternal' | 'maternal' = 'paternal',
  stage?: StageCall,
  /** Whole-chromosome calls from allele dosage, keyed `chrN`. Fill a gap the genotype channel
   *  cannot fill rather than competing with it: they never cover the same event. */
  dosageCalls: readonly {
    where: string, verdict: string, posterior: number, markers: number, excluded: number,
    why: string,
  }[] = [],
): Defect[] {
  const byOne = new Map(oneParent.map((o) => [o.where, o]))
  const byDosage = new Map(dosageCalls.map((d) => [d.where, d]))
  const other = loadedParent === 'paternal' ? 'maternal' : 'paternal'
  const byWhere = new Map<string, GainAnnotation>()
  for (const a of [...gains, ...losses]) byWhere.set(a.where, a)
  return segments.map((sg) => {
    const co = segmentCoords(sg)
    const where = `chr${sg.chrom} ${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb`
    const ann = byWhere.get(where)
    let origin: Defect['origin'] = ann?.origin?.startsWith('paternal') ? 'paternal'
      : ann?.origin?.startsWith('maternal') ? 'maternal' : 'unclear'
    // A single-parent call names a parent where two-parent dosage could not, so it fills the gap
    // rather than overriding: two parents remain the stronger evidence where both were loaded.
    const one = byOne.get(where)
    if (origin === 'unclear' && one) {
      if (one.verdict === 'known-parent-lost') origin = loadedParent
      else if (one.verdict === 'other-parent-lost') origin = other
    }
    // Dosage fills what genotypes could not reach. A whole chromosome is detected by its genotype
    // call rate collapsing, so on those events the genotype channel has no evidence left and this
    // is the only channel that can speak. It never overrides a genotype answer.
    const dose = byDosage.get(`chr${sg.chrom}`)
    const usedDosage = origin === 'unclear' && !!dose
      && (dose.verdict === 'known-parent-lost' || dose.verdict === 'other-parent-lost')
    if (usedDosage && dose) {
      origin = dose.verdict === 'known-parent-lost' ? loadedParent : other
    }
    return {
      chrom: sg.chrom,
      startBp: co.start,
      endBp: co.end,
      locus: locus(sg.chrom, co.start, co.end),
      kind: sg.kind === 'copy-gain' ? 'copy-gain'
        : sg.kind === 'copy-loss' ? 'copy-loss' : 'parental-absence',
      interval: sg.refined ? co.interval : undefined,
      origin,
      why: ann?.why ?? (usedDosage ? dose?.why : undefined) ?? one?.why
        ?? 'No parental genotype was loaded for this run, so allele dosage cannot be oriented and '
          + 'the parent is not determined. The position and its interval are unaffected.',
      basis: ann ? 'two-parent'
        : usedDosage ? 'dosage'
          : (one && origin !== 'unclear') ? 'one-parent' : undefined,
      channel: usedDosage ? 'dosage' : (one ? 'genotype' : undefined),
      excludedDosage: usedDosage && dose ? dose.excluded : undefined,
      informative: one?.markers,
      posterior: one && Number.isFinite(one.posterior) ? one.posterior : undefined,
      exclusive: one?.exclusive,
      stage: stage?.stage,
      dropout: stage && Number.isFinite(stage.dropout) ? stage.dropout : undefined,
      dropoutBasis: stage?.basis === 'none' ? undefined : stage?.basis,
    } as Defect
  })
}

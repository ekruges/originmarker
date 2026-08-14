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
  basis?: 'two-parent' | 'one-parent' | 'sibling' | 'homologue'
  informative?: number
  posterior?: number
  phi?: number
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
): Defect[] {
  const byOne = new Map(oneParent.map((o) => [o.where, o]))
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
    return {
      chrom: sg.chrom,
      startBp: co.start,
      endBp: co.end,
      locus: locus(sg.chrom, co.start, co.end),
      kind: sg.kind === 'copy-gain' ? 'copy-gain'
        : sg.kind === 'copy-loss' ? 'copy-loss' : 'parental-absence',
      interval: sg.refined ? co.interval : undefined,
      origin,
      why: ann?.why ?? one?.why
        ?? 'No parental genotype was loaded for this run, so allele dosage cannot be oriented and '
          + 'the parent is not determined. The position and its interval are unaffected.',
      basis: ann ? 'two-parent' : (one && origin !== 'unclear') ? 'one-parent' : undefined,
      informative: one?.markers,
      posterior: one && Number.isFinite(one.posterior) ? one.posterior : undefined,
    } as Defect
  })
}

/**
 * Whether an array can be one haploid meiotic product, and whether it is good enough to be one.
 *
 * Plain TypeScript rather than part of the panel that renders it, because this is the gate the
 * whole feature stands on: nothing is compared and nothing is built until it passes. Living in a
 * .tsx put it out of reach of every check that runs under node, including the accuracy audit,
 * which then could not exercise the shipped rule at all.
 */
import { CALL_RATE_FLOOR, HET_BAND_DIPLOID, pct } from './parentage.ts'
import type { SampleProfile } from './ingest.ts'

/**
 * Three bands, not two. The cut between haploid and diploid is not clean: confirmed products
 * run 1.0% to 7.8% and arrays excluded as diploid measured 10.0% and 10.7%.
 *
 * Two signals, because neither alone is sufficient. A known unrelated diploid adult measures a
 * 7.17% BAF band, under the 8% gate, and passes on the band alone. It cost nothing in theory
 * and everything in practice: admitted as a product it read 10.9% against a genuine product,
 * which fragmented one father's group into two and produced three groups where there are two.
 *
 * Genotype heterozygosity catches it. A haploid genome cannot be heterozygous at all, so every
 * such call is error: twelve confirmed products ran 5.3% to 13.7% and that adult reads 15.0%.
 * The margin is 1.3 points, which is thin, so a sample that fails EITHER signal is withheld
 * rather than argued about.
 *
 * SCOPE, and it is narrower than it looks. Every figure above is a haploid meiotic product against
 * a BULK diploid adult, and that is the only distinction this gate has ever been measured on. On
 * the 46 arrays this tool has run it does that job with nothing misread. It does NOT generalise to
 * post-zygotic material: measured on 120 non-haploid arrays spanning cleavage blastomeres,
 * trophectoderm biopsies and ESC lines, 91 of the 120 fall inside the haploid BAF-band range and
 * 22 blastomeres sit below the 10% diploid-exclusion floor. Whole-genome amplification widens the
 * heterozygote band and fills the homozygote band until the two are one distribution, and a BAF
 * core fraction fails the same way. So this gate must not be used as a general ploidy branch.
 * `obligateHet.ts` is the statistic that separates there, and it needs at least one parent, which
 * is why this one still exists for triage that has none.
 */
export type Ploidy = 'haploid' | 'borderline' | 'diploid'

/** Above this, genotype heterozygosity is diploid rather than a haploid's error rate. */
export const HET_RATE_DIPLOID = 0.14

export const ploidyOf = (band: number, hetRate: number): Ploidy => {
  if (!Number.isFinite(band)) return 'borderline'
  if (band >= 0.15 || hetRate >= 0.17) return 'diploid'
  if (band > HET_BAND_DIPLOID || hetRate >= HET_RATE_DIPLOID) return 'borderline'
  return 'haploid'
}

export interface Triage {
  id: string
  name: string
  profile: SampleProfile
  ploidy: Ploidy
  usable: boolean
  why: string
}

export function triage(id: string, name: string, p: SampleProfile): Triage {
  const ploidy = ploidyOf(p.hetBand, p.hetRate)
  if (p.callRate < CALL_RATE_FLOOR) {
    return {
      id,
      name,
      profile: p,
      ploidy,
      usable: false,
      why: `call rate ${pct(p.callRate, 1)} is below the ${pct(CALL_RATE_FLOOR, 0)} floor`,
    }
  }
  if (ploidy !== 'haploid') {
    return {
      id,
      name,
      profile: p,
      ploidy,
      usable: false,
      why: ploidy === 'diploid'
        ? `BAF band ${pct(p.hetBand, 1)} and ${pct(p.hetRate, 1)} heterozygous calls read `
          + 'diploid, not one meiotic product'
        : `BAF band ${pct(p.hetBand, 1)} with ${pct(p.hetRate, 1)} heterozygous calls sits `
          + 'between the haploid and diploid ranges',
    }
  }
  return { id, name, profile: p, ploidy, usable: true, why: '' }
}

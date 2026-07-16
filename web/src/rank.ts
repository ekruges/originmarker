import type { Marker } from './api'

/** Expected heterozygosity 2pq: a population prior, never a given carrier's genotype. */
export const twoPQ = (p: number) => 2 * p * (1 - p)

/** Ancestry-matched 2pq prior, or null when that population has no MAF for this marker. */
export function ancestryHet(m: Marker, ancestry: string | null): number | null {
  if (!ancestry) return null
  const p = m.per_pop_maf?.[ancestry]
  return p == null ? null : twoPQ(p)
}

/**
 * The 2pq prior used for ranking under the selected ancestry.
 *
 * Must mirror panelbuilder._rank_key, whose fallback is the GLOBAL 2pq and never
 * het_max_pop: the PDF stamps "global 2pq (het)" as the ranking key, and a mirror that
 * drifts makes screen and export order the same panel differently.
 */
export const rankHet = (m: Marker, ancestry: string | null) => ancestryHet(m, ancestry) ?? m.het

/**
 * The 2pq prior a chart may DRAW under an axis titled with the selected ancestry.
 *
 * Not rankHet: an axis naming one population must not carry another's figure, so there is
 * no global fallback here. Null means no MAF for that population, which is not zero and
 * has no height.
 */
export const shownHet = (m: Marker, ancestry: string | null): number | null =>
  ancestry ? ancestryHet(m, ancestry) : m.het

/** Client-side mirror of panelbuilder._rank_key: (-het, -global het, |dist|). */
export function compareRank(a: Marker, b: Marker, ancestry: string | null): number {
  return (
    rankHet(b, ancestry) - rankHet(a, ancestry) ||
    b.het - a.het ||
    Math.abs(a.dist) - Math.abs(b.dist)
  )
}

export const rankMarkers = (ms: Marker[], ancestry: string | null): Marker[] =>
  [...ms].sort((a, b) => compareRank(a, b, ancestry))

// --- column presets (sort orders only; never LD) ------------------------------
export type Preset = 'ranked' | 'closest' | 'heterozygous' | 'robust'

export const PRESETS: { value: Preset; label: string }[] = [
  { value: 'ranked', label: 'Engine rank' },
  { value: 'closest', label: 'Closest' },
  { value: 'heterozygous', label: 'Most heterozygous' },
  { value: 'robust', label: 'Cross-ancestry robust' },
]

/** Lowest 2pq prior across populations: high => useful in every ancestry. */
export function minPopHet(m: Marker): number {
  const vs = Object.values(m.per_pop_maf ?? {})
  return vs.length ? Math.min(...vs.map(twoPQ)) : 0
}

export function applyPreset(ms: Marker[], preset: Preset, ancestry: string | null): Marker[] {
  switch (preset) {
    case 'closest':
      return [...ms].sort((a, b) => Math.abs(a.dist) - Math.abs(b.dist))
    case 'heterozygous':
      return [...ms].sort((a, b) => rankHet(b, ancestry) - rankHet(a, ancestry))
    case 'robust':
      return [...ms].sort((a, b) => minPopHet(b) - minPopHet(a) || Math.abs(a.dist) - Math.abs(b.dist))
    default:
      return rankMarkers(ms, ancestry)
  }
}

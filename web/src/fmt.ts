/** Display helpers. Coordinates never render without their build. */

export const int = (n: number) => n.toLocaleString('en-US')

/** Signed bp distance, always with an explicit sign so side is unambiguous. */
export const signedBp = (d: number) => `${d > 0 ? '+' : d < 0 ? '−' : ''}${int(Math.abs(d))}`

export const num = (n: number | null | undefined, digits = 3) =>
  n == null ? '-' : n.toFixed(digits)

export const sci = (n: number | null | undefined) => {
  if (n == null) return '-'
  if (n === 0) return '0'
  return n < 1e-3 ? n.toExponential(2) : n.toFixed(5)
}

/**
 * Two significant figures, for quantities read off the genetic map.
 *
 * Not fixed decimals: the map is sex-averaged and interpolated between markers, and spans
 * 5e-5 to ~0.8 cM here, so any fixed decimal count claims a resolution it lacks at one end.
 */
export const sig2 = (n: number | null | undefined) =>
  n == null ? '-' : n === 0 ? '0' : String(Number(n.toPrecision(2)))

/**
 * A missing fact renders as the word, never as a neighbouring value in scope.
 *
 * For provenance, which may not change after the fact. Spelled out rather than '-',
 * which reads as formatting rather than as an absence.
 */
export const orUnknown = (v: string | number | null | undefined) =>
  v == null || v === '' ? 'unknown' : String(v)

export const coord = (chrom: string, pos: number, build = 'GRCh38') =>
  `${build} chr${chrom.replace(/^chr/, '')}:${int(pos)}`

export const shortCoord = (chrom: string, pos: number) => `chr${chrom.replace(/^chr/, '')}:${int(pos)}`

export const strandLabel = (s: number | null | undefined) =>
  s === 1 ? '+ (plus)' : s === -1 ? '− (minus)' : 'unknown'

/** Axis tick label: 17.40 Mb / 17,397 kb / bare bp, whichever suits the span. */
export function axisLabel(pos: number, spanBp: number) {
  if (spanBp > 2e6) return `${(pos / 1e6).toFixed(2)} Mb`
  if (spanBp > 2e3) return `${(pos / 1e3).toFixed(spanBp > 2e5 ? 0 : 1)} kb`
  return int(pos)
}

export const utc = (iso: string) => {
  const d = new Date(iso)
  return isNaN(d.getTime()) ? iso : `${d.toISOString().slice(0, 16).replace('T', ' ')} UTC`
}

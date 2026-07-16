/**
 * The OriginMarker lockup: the monogram beside the wordmark.
 *
 * Merriweather Italic, self-hosted from src/assets under the SIL OFL. The display face is
 * used only here: coordinates and allele counts stay in the system sans.
 */
import { Mark } from './Mark'

export function Logo({ size = 'md', mark = true }: { size?: 'md' | 'hero'; mark?: boolean }) {
  const hero = size === 'hero'
  const fontSize = hero ? 46 : 19
  return (
    <span
      className="om-logo"
      style={{
        display: 'inline-flex',
        alignItems: 'center',
        // The wordmark is italic, so its first stroke leans away from the mark and opens
        // an optical gap wider than the metric one. Set tight and let the lean do the rest.
        gap: hero ? 8 : 4,
        fontSize,
        lineHeight: hero ? 1.15 : 1.5,
      }}
    >
      {/* Optical, not metric: the ring's round outer edge reads smaller than the wordmark's
          flat cap height at the same pixel size, so the mark is set slightly larger. */}
      {mark && <Mark size={Math.round(fontSize * 1.16)} />}
      <span>
        <span className="om-logo-origin">Origin</span>
        <span className="om-logo-marker">Marker</span>
      </span>
    </span>
  )
}

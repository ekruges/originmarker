/**
 * The OriginMarker wordmark.
 *
 * Merriweather Italic, self-hosted from src/assets under the SIL OFL. The display face is
 * used only here: coordinates and allele counts stay in the system sans.
 */

export function Logo({ size = 'md' }: { size?: 'md' | 'hero' }) {
  const hero = size === 'hero'
  return (
    <span
      className="om-logo"
      style={{
        fontSize: hero ? 46 : 19,
        lineHeight: hero ? 1.15 : 1.5,
        display: 'inline-block',
      }}
    >
      <span className="om-logo-origin">Origin</span>
      <span className="om-logo-marker">Marker</span>
    </span>
  )
}

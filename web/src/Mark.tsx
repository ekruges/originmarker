/**
 * The OriginMarker monogram: a geometric sans ring around the wordmark's own serif M.
 *
 * The M is the real Merriweather Italic Bold glyph, outlined from the font in
 * web/src/assets (SIL OFL) rather than redrawn, so the mark and the wordmark are the same
 * letterform. Geometry is shared with web/public/favicon.svg: change one, change both.
 */

const D = 'M-79 0 -68 93 92 120 411 1364 223 1387 237 1486H837L950 583L973 398L1041 588L1403 1486H1936L1924 1387L1750 1364L1715 120L1888 93L1881 0H1243L1251 95L1435 120L1483 973L1505 1370L1361 1021L921 -1L777 0L630 997L573 1368L490 987L285 120L496 93L487 0Z'
const TRANSFORM = 'translate(26.500,68.793) scale(0.025310,-0.025310)'

export function Mark({ size = 22 }: { size?: number }) {
  return (
    <svg viewBox="0 0 100 100" width={size} height={size} aria-hidden focusable="false"
         style={{ display: 'block', flex: 'none' }}>
      <circle cx="50.0" cy="50.0" r="38.0" fill="none" stroke="var(--om-blue)" strokeWidth="10" />
      <path d={D} fill="var(--om-blue-light)" transform={TRANSFORM} />
    </svg>
  )
}

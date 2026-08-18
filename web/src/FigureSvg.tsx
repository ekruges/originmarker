/**
 * Render a Figure to SVG. It knows how to draw four primitives and nothing else.
 *
 * Every decision about where a mark goes was made in comparisonFigures.ts, which is what makes this
 * and the report's renderer the same figure rather than two drawings of one idea.
 */
import type { Figure, Mark } from './figures.ts'
import { FIG } from './figures.ts'

const draw = (m: Mark, i: number) => {
  if (m.k === 'line') {
    return (
      <line
        key={i} x1={m.x1} y1={m.y1} x2={m.x2} y2={m.y2}
        stroke={m.colour} strokeWidth={m.width}
        strokeDasharray={m.dash ? '2 3' : undefined}
      />
    )
  }
  if (m.k === 'rect') {
    return <rect key={i} x={m.x} y={m.y} width={m.w} height={m.h} fill={m.colour} />
  }
  if (m.k === 'dot') {
    return (
      <circle
        key={i} cx={m.x} cy={m.y} r={m.r}
        fill={m.hollow ? 'none' : m.colour} stroke={m.hollow ? m.colour : 'none'} strokeWidth={0.9}
      />
    )
  }
  return (
    <text
      key={i} x={m.x} y={m.y} fontSize={m.size} fill={m.colour}
      textAnchor={m.anchor === 'start' ? 'start' : m.anchor === 'end' ? 'end' : 'middle'}
      fontWeight={m.bold ? 700 : 400}
      fontStyle={m.italic ? 'italic' : undefined}
    >
      {m.s}
    </text>
  )
}

export function FigureSvg({ figure }: { figure: Figure }) {
  return (
    <figure style={{ margin: '0 0 16px' }}>
      <svg
        viewBox={`0 0 ${figure.w} ${figure.h}`}
        width="100%"
        style={{ maxWidth: figure.w * 1.35, display: 'block' }}
        role="img"
        aria-label={figure.caption}
      >
        {figure.marks.map(draw)}
      </svg>
      <figcaption style={{
        fontSize: 11, color: 'var(--om-text-dim)', lineHeight: 1.45, marginTop: 2,
        maxWidth: figure.w * 1.35,
      }}
      >
        {figure.caption}
      </figcaption>
    </figure>
  )
}

export { FIG }

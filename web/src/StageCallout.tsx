/**
 * The inferred stage, at the top of a run because it sets every number below it.
 *
 * The stage picks the drift floor, the variance inflation, the detection floors and the dropout
 * that parameterises every directional call, and those differ by up to seventeen-fold across the
 * classes. A reader who cannot see which stage was inferred cannot judge anything underneath it.
 *
 * Both measurements are shown side by side on purpose. They are easy to confuse and they measure
 * different things, and reading only the first is what put amplified samples on the bulk rung.
 */
import { Text } from '@mantine/core'
import { amplificationWord, type StageFacts } from './stage.ts'

const pct = (x: number, d = 1) => (Number.isFinite(x) ? `${(100 * x).toFixed(d)}%` : '-')

const Axis = ({ label, value, meaning }: { label: string, value: string, meaning: string }) => (
  <div style={{ flex: '1 1 200px', minWidth: 190 }}>
    <div style={{ fontSize: 10, textTransform: 'uppercase', letterSpacing: 0.4,
      color: 'var(--om-text-dim)' }}
    >
      {label}
    </div>
    <div style={{ fontFamily: 'var(--om-mono)', fontSize: 15, fontWeight: 700 }}>{value}</div>
    <div style={{ fontSize: 11, color: 'var(--om-text-dim)', lineHeight: 1.4 }}>{meaning}</div>
  </div>
)

export function StageCallout({ facts }: { facts: StageFacts | null }) {
  if (!facts) return null
  const f = facts
  return (
    <div style={{
      border: '1px solid var(--om-line)', borderLeft: '4px solid var(--om-blue)',
      padding: '10px 14px', margin: '10px 0 4px',
    }}
    >
      <Text style={{ fontSize: 15, fontWeight: 700, lineHeight: 1.25 }}>
        Material: {f.stage}
        {f.demotedFromBulk ? ' (heterozygosity alone would have read as bulk)' : ''}
      </Text>

      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 16, marginTop: 8 }}>
        <Axis
          label="Heterozygous, of called"
          value={pct(f.hetRate, 2)}
          meaning="how many heterozygous markers survived, which tracks how much template there was"
        />
        <Axis
          label="Allele-fraction spread"
          value={Number.isFinite(f.hetBafSd) ? f.hetBafSd.toFixed(4) : '-'}
          meaning={`how far the survivors scattered: ${amplificationWord(f)}`}
        />
        <Axis
          label="Allele dropout"
          value={Number.isFinite(f.dropout) ? f.dropout.toFixed(3) : '-'}
          meaning={`basis: ${f.basis.replace(/-/g, ' ')}. Parameterises every directional call below`}
        />
      </div>

      <Text style={{ fontSize: 11, color: 'var(--om-text-dim)', marginTop: 8, lineHeight: 1.45 }}>
        {f.why}.
      </Text>
      <Text style={{ fontSize: 11, color: 'var(--om-text-dim)', marginTop: 4, lineHeight: 1.45 }}>
        Limits of this figure: {f.caveat}.
      </Text>
    </div>
  )
}

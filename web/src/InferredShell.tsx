import { useState, type ReactNode } from 'react'
import { Badge, Group, Paper, Text } from '@mantine/core'

/**
 * The shell the inferred-reference screens sit in.
 *
 * An earlier version laid every step out at full length down one page, which reads as a report
 * rather than a tool: a user scanning for the answer had to read the reasoning first. Here each
 * stage states its outcome in one line and keeps the evidence behind a disclosure, matching the
 * verdict-then-detail idiom the sample cards already use.
 */

export type StageState = 'ok' | 'attention' | 'blocked' | 'pending'

const TONE: Record<StageState, { bar: string; text: string }> = {
  ok: { bar: 'var(--om-blue)', text: 'var(--om-blue)' },
  attention: { bar: 'var(--om-higher)', text: 'var(--om-higher)' },
  blocked: { bar: 'var(--om-higher)', text: 'var(--om-higher)' },
  pending: { bar: 'var(--om-border-strong)', text: 'var(--om-text-dim)' },
}

/**
 * One gate. The headline is the outcome, not the question, so a reader who only skims the
 * headlines still gets the whole answer in order.
 */
export function Stage({ n, title, state, headline, badges = [], children, defaultOpen = false }: {
  n: number
  title: string
  state: StageState
  headline: ReactNode
  badges?: string[]
  children?: ReactNode
  defaultOpen?: boolean
}) {
  const [open, setOpen] = useState(defaultOpen)
  const tone = TONE[state]
  return (
    <Paper
      withBorder
      radius={0}
      style={{ borderLeft: `3px solid ${tone.bar}`, marginTop: 10, overflow: 'hidden' }}
    >
      <div style={{ padding: '10px 14px' }}>
        <Group gap={8} align="baseline" wrap="nowrap">
          <span style={{
            fontFamily: 'var(--om-mono)', fontSize: 10.5, color: 'var(--om-text-dim)',
            flex: 'none',
          }}
          >
            {String(n).padStart(2, '0')}
          </span>
          <span style={{
            fontSize: 11, textTransform: 'uppercase', letterSpacing: '0.06em',
            color: 'var(--om-text-dim)', fontWeight: 600, flex: 'none',
          }}
          >
            {title}
          </span>
          {badges.map((b) => (
            <Badge key={b} size="xs" variant="outline" color="genomeGrey" radius={2}>{b}</Badge>
          ))}
        </Group>
        <Text
          mt={4}
          style={{ fontSize: 15, fontWeight: 600, letterSpacing: '-0.01em', lineHeight: 1.35 }}
          c={state === 'ok' ? undefined : tone.text}
        >
          {headline}
        </Text>
        {children && (
          <button
            type="button"
            className="om-primer-toggle"
            aria-expanded={open}
            onClick={() => setOpen((o) => !o)}
            style={{ marginTop: 6 }}
          >
            <span className="om-primer-toggle-caret">{open ? '▾' : '▸'}</span>
            <span>{open ? 'Hide the evidence' : 'Show the evidence'}</span>
          </button>
        )}
      </div>
      {open && children && (
        <div style={{
          padding: '2px 14px 14px', borderTop: '1px solid var(--om-border)',
          background: 'var(--om-zebra)',
        }}
        >
          {children}
        </div>
      )}
    </Paper>
  )
}

/** The one-line answer, before any of the stages. */
export function Summary({ verdict, sub, stats }: {
  verdict: string
  sub: string
  stats: { label: string; value: string }[]
}) {
  return (
    <Paper withBorder radius={0} style={{ padding: '14px 16px', background: 'var(--om-head-bg)' }}>
      <Text style={{ fontSize: 22, fontWeight: 700, letterSpacing: '-0.015em', lineHeight: 1.15 }}>
        {verdict}
      </Text>
      <Text size="sm" c="dimmed" mt={2} style={{ maxWidth: 720 }}>{sub}</Text>
      <div style={{
        display: 'flex', gap: 22, flexWrap: 'wrap', marginTop: 12, paddingTop: 10,
        borderTop: '1px solid var(--om-border)',
      }}
      >
        {stats.map((s) => (
          <div key={s.label}>
            <div style={{
              fontFamily: 'var(--om-mono)', fontSize: 16, fontVariantNumeric: 'tabular-nums',
              letterSpacing: '-0.01em', lineHeight: 1.2,
            }}
            >
              {s.value}
            </div>
            <div style={{ fontSize: 10, color: 'var(--om-text-dim)', letterSpacing: '0.02em' }}>
              {s.label}
            </div>
          </div>
        ))}
      </div>
    </Paper>
  )
}

/** A short explanation that only appears when asked for, so the reasoning never crowds out the
 *  result it explains. */
export function Why({ children }: { children: ReactNode }) {
  const [open, setOpen] = useState(false)
  return (
    <div style={{ marginTop: 8 }}>
      <button
        type="button" className="om-primer-toggle" aria-expanded={open}
        onClick={() => setOpen((o) => !o)}
      >
        <span className="om-primer-toggle-caret">{open ? '▾' : '▸'}</span>
        <span>Why this is the rule</span>
      </button>
      {open && (
        <Text size="xs" c="dimmed" mt={4} style={{ maxWidth: 760, lineHeight: 1.55 }}>
          {children}
        </Text>
      )}
    </div>
  )
}

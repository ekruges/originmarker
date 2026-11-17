/**
 * The changes of one sample, as chips that open.
 *
 * A chip carries WHERE, WHOSE and one word for the kind of call: found two ways, consistent, or
 * one reading. The measured count behind the word is one click away, in the opened call, with the
 * scan's own coordinates and its own pair of numbers. No band letter: on the answer key the band
 * did not track right and wrong, the kind of call did (audit/truth/page-says-2026-11-17.txt).
 *
 * Counted changes first. Stretches seen by one reading follow under their own label and are not in
 * the count; a whole-chromosome imbalance is listed last as read, not counted.
 */
import { useState, type ReactNode } from 'react'
import { Text } from '@mantine/core'
import {
  KIND_SAYS, lrrText, locusText, measuredSentence, pairText, spanMb, whereText, whoseText,
  type ChangeList, type ChangeRow, type Setup,
} from './pageSays.ts'
import { sentences } from './defects.ts'

const strong = (g: ChangeRow): boolean => g.word === 'found two ways' || g.word === 'consistent'
  || g.word === 'position only'

function OneChange({ g, setup, plot }: { g: ChangeRow; setup: Setup; plot?: (g: ChangeRow) => ReactNode }) {
  const [open, setOpen] = useState(false)
  const ink = strong(g) ? 'var(--om-defect)' : 'var(--om-text-dim)'
  const why = [...new Set(g.pieces.flatMap((p) => sentences(p.why ?? '')))].join(' ')
  return (
    <div style={{
      border: `1px solid ${open ? ink : 'var(--om-defect)'}`, background: 'var(--om-defect-chip)',
      gridColumn: open ? '1 / -1' : undefined, minWidth: 0,
    }}
    >
      <button
        type="button" onClick={() => setOpen((v) => !v)} aria-expanded={open} title={whereText(g)}
        style={{
          display: 'flex', alignItems: 'baseline', gap: 7, width: '100%', border: 0, background: 'none',
          cursor: 'pointer', padding: '3px 7px', font: 'inherit', textAlign: 'left', lineHeight: 1.5,
        }}
      >
        <span style={{ fontFamily: 'var(--om-mono)', fontSize: 12, fontWeight: 600, whiteSpace: 'nowrap' }}>
          {g.chrom === 'genome' ? 'genome' : `chr${g.chrom}`}
        </span>
        <span style={{ fontFamily: 'var(--om-mono)', fontSize: 12, fontVariantNumeric: 'tabular-nums', whiteSpace: 'nowrap' }}>
          {g.whole ? 'whole' : `${spanMb(g)} Mb`}
        </span>
        <span style={{ fontSize: 12 }}>{whoseText(g)}</span>
        <span style={{ flex: 1 }} />
        <span style={{ fontSize: 11, fontWeight: 700, color: ink, whiteSpace: 'nowrap' }}>{g.word}</span>
      </button>
      {open ? (
        <div style={{ padding: '2px 7px 7px', borderTop: '1px solid var(--om-defect)' }}>
          <Text style={{ fontSize: 13, fontWeight: 700, lineHeight: 1.3, color: ink, margin: '6px 0 4px' }}>
            {whereText(g)}: {whoseText(g)}.
          </Text>
          <Text size="sm" style={{ maxWidth: 780, lineHeight: 1.5 }}>
            {KIND_SAYS[g.call]} {measuredSentence(g.call, setup)}
          </Text>
          <div style={{ marginTop: 5, display: 'grid', gap: 2 }}>
            {g.pieces.map((p) => (
              <Text key={`${p.startBp}`} size="xs" ff="monospace">
                {[locusText(p), pairText(p), lrrText(p)].filter(Boolean).join(' · ')}
              </Text>
            ))}
          </div>
          {plot?.(g)}
          {why ? <Text size="xs" c="dimmed" mt={6} style={{ maxWidth: 780, lineHeight: 1.5 }}>{why}</Text> : null}
        </div>
      ) : null}
    </div>
  )
}

const Grid = ({ rows, setup, plot }: { rows: ChangeRow[]; setup: Setup; plot?: (g: ChangeRow) => ReactNode }) => (
  <div style={{
    marginTop: 6, display: 'grid', gap: 4,
    gridTemplateColumns: 'repeat(auto-fill, minmax(280px, 1fr))', alignItems: 'start',
  }}
  >
    {rows.map((g) => <OneChange key={`${g.chrom}:${g.startBp}:${g.call}:${g.origin}`} g={g} setup={setup} plot={plot} />)}
  </div>
)

export function DefectCallout({ list, setup, plot }: {
  list: ChangeList
  setup: Setup
  /** The plot for a row, shown when the row is open. */
  plot?: (g: ChangeRow) => ReactNode
}) {
  if (!list.counted.length && !list.once.length && !list.seen.length) return null
  return (
    <div style={{ margin: '8px 0 2px' }}>
      {list.counted.length ? <Grid rows={list.counted} setup={setup} plot={plot} /> : null}
      {list.once.length ? (
        <>
          <Text size="xs" fw={700} c="dimmed" mt={8}>Seen once, not confirmed</Text>
          <Grid rows={list.once} setup={setup} plot={plot} />
        </>
      ) : null}
      {list.seen.length ? (
        <>
          <Text size="xs" fw={700} c="dimmed" mt={8}>Also read, not counted as a change</Text>
          <Grid rows={list.seen} setup={setup} plot={plot} />
        </>
      ) : null}
    </div>
  )
}

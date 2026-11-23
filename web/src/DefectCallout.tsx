/**
 * The changes of one sample, as one list of rows that open.
 *
 * A row carries WHERE, WHOSE and one word for the kind of call: found two ways, consistent, or
 * one reading. The measured count behind the word is one click away, in the opened call, with the
 * scan's own coordinates and its own pair of numbers. No band letter: on the answer key the band
 * did not track right and wrong, the kind of call did (audit/truth/page-says-2026-11-17.txt).
 *
 * Counted changes first, outlined. Stretches seen by one reading follow under their own label in
 * plain muted text and are not in the count; a whole-chromosome imbalance is listed last as read,
 * not counted.
 */
import { useState, type ReactNode } from 'react'
import { Text } from '@mantine/core'
import {
  disturbedText, KIND_SAYS, lrrText, locusText, measuredSentence, pairText, spanMb, whenLines, whereText, whoseText,
  type ChangeList, type ChangeRow, type Setup,
} from './pageSays.ts'
import { sentences } from './defects.ts'

const strong = (g: ChangeRow): boolean => g.word === 'found two ways' || g.word === 'consistent'
  || g.word === 'position only'

/** One set of columns for every row, so the list reads down: chromosome, size, what, word. */
const COLUMNS = '3.9em 5em minmax(0, 1fr) auto'

function OneChange({ g, setup, plot }: { g: ChangeRow; setup: Setup; plot?: (g: ChangeRow) => ReactNode }) {
  const [open, setOpen] = useState(false)
  const loud = strong(g)
  const why = [...new Set(g.pieces.flatMap((p) => sentences(p.why ?? '')))].join(' ')
  const figure = open ? plot?.(g) : null
  // An arm-level piece has no measured edge or extent inside its arm, so its size cell names the
  // arm rather than printing the piece's span, which would read as a measured length.
  const arm = g.pieces.find((p) => p.arm)?.arm
  const when = whenLines(g)
  return (
    <div style={{
      border: open ? '1px solid var(--om-border-strong)' : loud ? '1px solid var(--om-defect)' : '1px solid transparent',
      background: open ? '#fff' : undefined, minWidth: 0,
    }}
    >
      <button
        type="button" onClick={() => setOpen((v) => !v)} aria-expanded={open} title={whereText(g)}
        style={{
          display: 'grid', gridTemplateColumns: COLUMNS, columnGap: 8, alignItems: 'baseline', width: '100%',
          border: 0, background: loud ? 'var(--om-defect-chip)' : 'none', color: loud ? undefined : 'var(--om-text-dim)',
          cursor: 'pointer', padding: '3px 8px', font: 'inherit', textAlign: 'left', lineHeight: 1.5,
        }}
      >
        <span style={{ fontFamily: 'var(--om-mono)', fontSize: 12, fontWeight: loud ? 600 : 400, whiteSpace: 'nowrap' }}>
          {g.chrom === 'genome' ? 'genome' : `chr${g.chrom}`}
        </span>
        <span style={{ fontFamily: 'var(--om-mono)', fontSize: 12, fontVariantNumeric: 'tabular-nums', whiteSpace: 'nowrap' }}>
          {g.whole ? 'whole' : arm ? `${g.chrom}${arm} arm` : `${spanMb(g)} Mb`}
        </span>
        <span style={{ fontSize: 12 }}>{whoseText(g)}</span>
        <span style={{ fontSize: 11, fontWeight: loud ? 700 : 400, color: loud ? 'var(--om-defect)' : undefined, whiteSpace: 'nowrap' }}>
          {g.word}
        </span>
      </button>
      {open ? (
        <div style={{ padding: '4px 8px 8px' }}>
          {/* The figure's title states the finding; without a figure the row says it here. */}
          {figure ? null : (
            <Text style={{ fontSize: 13, fontWeight: 700, lineHeight: 1.3, margin: '6px 0 4px' }}>
              {whereText(g)}: {whoseText(g)}.
            </Text>
          )}
          <Text size="sm" mt={4} style={{ maxWidth: 780, lineHeight: 1.5 }}>
            {KIND_SAYS[g.call]} {measuredSentence(g.call, setup)}
          </Text>
          {when.map((w) => (
            <Text key={w} size="sm" mt={4} style={{ maxWidth: 780, lineHeight: 1.5 }}>{w}</Text>
          ))}
          {figure}
          <details style={{ marginTop: 6 }}>
            <summary style={{ fontSize: 12, cursor: 'pointer', color: 'var(--om-text-dim)' }}>Details</summary>
            <div style={{ marginTop: 4, display: 'grid', gap: 2 }}>
              {g.pieces.map((p) => (
                <Text key={`${p.startBp}`} size="xs" ff="monospace">
                  {[locusText(p), pairText(p), lrrText(p)].filter(Boolean).join(' · ')}
                </Text>
              ))}
            </div>
            {why ? <Text size="xs" c="dimmed" mt={6} style={{ maxWidth: 780, lineHeight: 1.5 }}>{why}</Text> : null}
          </details>
        </div>
      ) : null}
    </div>
  )
}

const Rows = ({ rows, setup, plot }: { rows: ChangeRow[]; setup: Setup; plot?: (g: ChangeRow) => ReactNode }) => (
  <div style={{ marginTop: 6, display: 'grid', gap: 3 }}>
    {rows.map((g) => <OneChange key={`${g.chrom}:${g.startBp}:${g.call}:${g.origin}`} g={g} setup={setup} plot={plot} />)}
  </div>
)

export function DefectCallout({ list, setup, plot }: {
  list: ChangeList
  setup: Setup
  /** The plot for a row, shown when the row is open. */
  plot?: (g: ChangeRow) => ReactNode
}) {
  if (!list.counted.length && !list.once.length && !list.seen.length && !list.disturbed) return null
  return (
    <div style={{ margin: '8px 0 2px' }}>
      {list.counted.length ? <Rows rows={list.counted} setup={setup} plot={plot} /> : null}
      {list.once.length ? (
        <>
          <Text size="xs" fw={700} c="dimmed" mt={12}>Seen once, not confirmed</Text>
          <Rows rows={list.once} setup={setup} plot={plot} />
        </>
      ) : null}
      {list.seen.length ? (
        <>
          <Text size="xs" fw={700} c="dimmed" mt={12}>Also read, not counted as a change</Text>
          <Rows rows={list.seen} setup={setup} plot={plot} />
        </>
      ) : null}
      {list.disturbed ? <Text size="xs" c="dimmed" mt={12}>{disturbedText(list)}</Text> : null}
    </div>
  )
}

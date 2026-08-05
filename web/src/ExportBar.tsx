import { useState } from 'react'
import { Button, Text } from '@mantine/core'
import { download } from './inferredExports'

/**
 * The exports, one per row.
 *
 * A user choosing between these should not have to guess which one their next tool reads, so
 * the format, the shape of the file and its consumer all sit on the row rather than in
 * documentation somewhere else.
 */

export interface ExportItem {
  label: string
  filename: string
  hint: string
  mime: string
  /** Async so a PDF, which has to typeset before it exists, fits the same row as a CSV. */
  build: () => string | Promise<Blob>
}

function Row({ item }: { item: ExportItem }) {
  const [busy, setBusy] = useState(false)
  const go = async (): Promise<void> => {
    setBusy(true)
    try {
      const out = await item.build()
      if (typeof out === 'string') download(out, item.filename, item.mime)
      else {
        const url = URL.createObjectURL(out)
        const a = document.createElement('a')
        a.href = url
        a.download = item.filename
        document.body.appendChild(a)
        a.click()
        a.remove()
        URL.revokeObjectURL(url)
      }
    } finally {
      setBusy(false)
    }
  }
  return (
    <div style={{
      display: 'flex', alignItems: 'center', gap: 14, padding: '10px 0',
      borderTop: '1px solid var(--om-border)',
    }}
    >
      <Button
        size="md"
        radius={2}
        loading={busy}
        onClick={go}
        style={{ minWidth: 210, fontFamily: 'var(--om-mono)', fontSize: 12, flex: 'none' }}
      >
        {item.filename}
      </Button>
      <div style={{ minWidth: 0 }}>
        <div style={{ fontSize: 12, fontWeight: 600, lineHeight: 1.35 }}>{item.label}</div>
        <div style={{ fontSize: 11, color: 'var(--om-text-dim)', lineHeight: 1.45 }}>
          {item.hint}
        </div>
      </div>
    </div>
  )
}

export function ExportBar({ items, note }: { items: ExportItem[]; note?: string }) {
  return (
    <div style={{ marginTop: 16, border: '1px solid var(--om-border)', padding: '14px 16px' }}>
      <Text style={{ fontSize: 14, fontWeight: 700, letterSpacing: '-0.01em' }}>Export</Text>
      <Text size="xs" c="dimmed" mt={2} mb={4} style={{ maxWidth: 760, lineHeight: 1.5 }}>
        Everything is written in this browser; no genotype is uploaded. Every file repeats that
        the reference was reconstructed rather than measured. In the CSVs that sits in a leading
        {' '}<code>#</code> block, which pandas and R skip by default, so the file still loads in
        one call.
      </Text>
      {items.map((it) => <Row key={it.filename} item={it} />)}
      {note && (
        <Text
          size="xs"
          c="dimmed"
          mt={12}
          style={{
            maxWidth: 760, lineHeight: 1.5, borderTop: '1px solid var(--om-border)',
            paddingTop: 10,
          }}
        >
          {note}
        </Text>
      )}
    </div>
  )
}

import { useState, type ReactNode } from 'react'
import { Group, Text } from '@mantine/core'

/**
 * The masthead and the file intake, shared by the features that have both.
 *
 * Syngamy and Progenitor are siblings rather than modes of one another, so each says its own
 * name before anything else. They also take the same kind of file in the same way, and two
 * implementations of that would drift: one gained drag-and-drop and the other did not.
 */

export function FeatureHeader({ name, tagline }: { name: string; tagline: string }) {
  return (
    <div style={{ borderBottom: '2px solid var(--om-blue)', paddingBottom: 8, marginBottom: 12 }}>
      <Group gap={10} align="baseline" wrap="nowrap">
        <span style={{
          fontFamily: 'Merriweather, Georgia, serif', fontSize: 30, fontWeight: 700,
          fontStyle: 'italic', letterSpacing: '-0.02em', color: 'var(--om-blue)', lineHeight: 1.05,
        }}
        >
          {name}
        </span>
        <span style={{
          fontSize: 11, textTransform: 'uppercase', letterSpacing: '0.08em',
          color: 'var(--om-text-dim)', fontWeight: 600,
        }}
        >
          {tagline}
        </span>
      </Group>
    </div>
  )
}

/**
 * The dashed intake area. Collapses to a thin strip once files are in, because at that point the
 * chips are the content and the invitation is not.
 */
export function DropZone({ empty, prompt, disabled, onFiles, children }: {
  empty: boolean
  prompt: string
  disabled?: boolean
  onFiles: (files: FileList) => void
  children?: ReactNode
}) {
  const [over, setOver] = useState(false)
  return (
    <div
      onDragOver={(e) => { e.preventDefault(); if (!disabled) setOver(true) }}
      onDragLeave={() => setOver(false)}
      onDrop={(e) => {
        e.preventDefault()
        setOver(false)
        if (!disabled && e.dataTransfer.files.length) onFiles(e.dataTransfer.files)
      }}
      style={{
        border: `1px dashed ${over ? 'var(--om-blue)' : 'var(--om-border-strong)'}`,
        background: over ? 'var(--om-head-bg)' : 'transparent',
        borderRadius: 2,
        padding: empty ? '20px 8px' : '6px',
        textAlign: empty ? 'center' : 'left',
      }}
    >
      {empty ? <Text size="xs" c="dimmed">{prompt}</Text> : children}
    </div>
  )
}

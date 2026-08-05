import { useLayoutEffect, useRef } from 'react'
import { Button, Group, Paper, Text } from '@mantine/core'

/**
 * The run log, shared by the features that have one.
 *
 * It follows itself while a run streams, and stops following the moment the reader scrolls up.
 * Both halves of that are load-bearing and neither is obvious, which is why there is one
 * implementation rather than one per feature.
 */

export interface LogLine<T extends string = string> { tag: T; text: string }

export function RunLog<T extends string>({ lines, colours, onDownload, title = 'Run log' }: {
  lines: LogLine<T>[]
  colours: Record<T, string>
  onDownload?: () => void
  title?: string
}) {
  const box = useRef<HTMLDivElement>(null)
  const pinned = useRef(true)
  // Set while this component is doing the scrolling, so its own assignment is not mistaken for
  // the reader scrolling away. Without it the log stops following itself: `scrollTop =
  // scrollHeight` fires onScroll, and the geometry read there lands mid-append often enough to
  // compute a false "they have scrolled up", after which new lines arrive below the fold in
  // silence.
  const selfScroll = useRef(false)

  // Layout effect, not effect: the scroll lands in the same frame the line is painted, so the
  // log never shows the previous bottom for a frame before jumping.
  useLayoutEffect(() => {
    const el = box.current
    if (!el || !pinned.current) return
    selfScroll.current = true
    el.scrollTop = el.scrollHeight
    requestAnimationFrame(() => { selfScroll.current = false })
  }, [lines])

  return (
    <Paper p="sm" mb={10}>
      <Group justify="space-between" align="center" mb={4}>
        <Text fw={600} size="xs" tt="uppercase" c="dimmed" style={{ letterSpacing: '0.04em' }}>
          {title}
        </Text>
        {onDownload && (
          <Button
            variant="default" size="compact-xs" px={6}
            aria-label="Download run log as a text file"
            onClick={onDownload}
          >
            <span aria-hidden style={{ fontSize: 13, lineHeight: 1 }}>&#x2913;</span>
          </Button>
        )}
      </Group>
      <div
        ref={box} className="om-mono"
        onScroll={() => {
          if (selfScroll.current) return
          const el = box.current
          // 24px, not 4: sub-pixel scroll positions and a line arriving mid-measure both put the
          // reader a few pixels off the exact bottom while they are plainly still at it.
          if (el) pinned.current = el.scrollHeight - el.scrollTop - el.clientHeight < 24
        }}
        style={{
          maxHeight: 220, overflowY: 'auto', padding: '4px 6px',
          background: 'var(--om-zebra)', border: '1px solid var(--om-border)',
          borderRadius: 2, fontSize: 11, lineHeight: 1.5, color: 'var(--om-text-dim)',
        }}
      >
        {lines.map((l, i) => (
          <div key={i} style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>
            <span style={{ color: colours[l.tag], display: 'inline-block', width: '8ch' }}>
              [{l.tag}]
            </span>
            {l.text}
          </div>
        ))}
      </div>
    </Paper>
  )
}

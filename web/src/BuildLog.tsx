import { useEffect, useId, useRef, useState } from 'react'
import { ActionIcon, Button, Tooltip } from '@mantine/core'
import type { LogLine, LogTag } from './api'
import { buildLogText, type ExportCtx, type LogMeta } from './logfile'

export type { LogMeta } from './logfile'

/** Lines kept and rendered. The engine caps its own buffer; this bounds the DOM. */
export const LOG_CAP = 500

// The tag is coloured, the line is not: a console of coloured prose is unreadable, and
// FETCH and CACHE are most of a build. WARN carries attention, not alarm: the build went
// on, and red belongs to the errors that stop one. SKIP is the account of what the build
// dropped on purpose, so it reads as ordinary, not as damage. One entry per LOG_TAGS
// member is what the compiler checks here, so a tag the engine gains cannot render
// uncoloured.
const TAG_COLOR: Record<LogTag, string> = {
  FETCH: 'var(--om-text-dim)',
  CACHE: 'var(--om-text-dim)',
  INFO: 'var(--om-blue)',
  WARN: 'var(--om-higher)',
  SKIP: 'var(--om-text-dim)',
  DONE: 'var(--om-blue)',
}

/** The build log, closed by default. Kept mounted after the build: that is when someone
 *  asks why it took 40 seconds. `openSignal` force-opens it: a verification run bumps it,
 *  so the pairs it appends are not streaming into a panel the reader has collapsed. */
export function BuildLog({
  lines, openSignal, meta,
}: {
  lines: LogLine[]
  openSignal?: number
  meta?: LogMeta
}) {
  const [open, setOpen] = useState(false)
  const box = useRef<HTMLDivElement>(null)
  // Follow the newest line only while the user is already at the bottom. Once they scroll
  // up to read something, a new line must not drag them back down.
  const pinned = useRef(true)
  const id = useId()

  useEffect(() => {
    const el = box.current
    if (el && pinned.current) el.scrollTop = el.scrollHeight
  }, [lines, open])

  // Force-open on a new signal, never force-CLOSE: the reader may have opened it themselves,
  // and a bump must not shut what they are reading. Skips the initial 0/undefined.
  useEffect(() => {
    if (openSignal) setOpen(true)
  }, [openSignal])

  const download = () => {
    const ctx: ExportCtx = {
      now: new Date().toISOString(),
      url: window.location.href,
      agent: navigator.userAgent,
    }
    const text = buildLogText(lines, meta ?? {}, ctx)
    const url = URL.createObjectURL(new Blob([text], { type: 'text/plain;charset=utf-8' }))
    const a = document.createElement('a')
    a.href = url
    a.download = `originmarker-buildlog-${meta?.jobId ?? 'panel'}.txt`
    document.body.appendChild(a)
    a.click()
    a.remove()
    URL.revokeObjectURL(url)
  }

  return (
    <div>
      <Button
        variant="subtle"
        color="gray"
        size="compact-xs"
        px={4}
        className="om-mono"
        onClick={() => setOpen((o) => !o)}
        aria-expanded={open}
        aria-controls={id}
      >
        <span aria-hidden style={{ marginRight: 6 }}>{open ? '▾' : '▸'}</span>
        Build log ({lines.length})
      </Button>
      {open && (
        // position:relative so the download button pins to the box's own top-right and does
        // not scroll away with the lines under it.
        <div style={{ position: 'relative', marginTop: 4 }}>
          <Tooltip label="Download this log as a .txt, with a build/debug footer" withArrow>
            <ActionIcon
              size="sm"
              variant="default"
              onClick={download}
              disabled={lines.length === 0}
              aria-label="Download build log as a text file"
              style={{ position: 'absolute', top: 4, right: 4, zIndex: 1 }}
            >
              <span aria-hidden style={{ fontSize: 13, lineHeight: 1 }}>⤓</span>
            </ActionIcon>
          </Tooltip>
          {/* A labelled region, not role="log": role="log" is a live region by default, and a
              screen reader reading every gnomAD chunk aloud as it lands is worse than silence.
              tabIndex makes it scrollable without a mouse. paddingRight clears the button. */}
          <div
            id={id}
            ref={box}
            role="region"
            aria-label="Build log"
            tabIndex={0}
            onScroll={() => {
              const el = box.current
              if (el) pinned.current = el.scrollHeight - el.scrollTop - el.clientHeight < 4
            }}
            className="om-mono"
            style={{
              maxHeight: 220,
              overflowY: 'auto',
              padding: '4px 6px',
              paddingRight: 34,
              background: 'var(--om-zebra)',
              border: '1px solid var(--om-border)',
              borderRadius: 2,
              fontSize: 11,
              lineHeight: 1.5,
              color: 'var(--om-text-dim)',
            }}
          >
            {lines.length === 0 ? (
              <div>no events yet</div>
            ) : (
              lines.map((l, i) => (
                <div key={i} style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>
                  {/* 8ch: the widest tag is [FETCH], so every line's text starts on the same
                      column whatever the tag. */}
                  <span style={{ color: TAG_COLOR[l.tag], display: 'inline-block', width: '8ch' }}>
                    [{l.tag}]
                  </span>
                  {l.text}
                </div>
              ))
            )}
          </div>
        </div>
      )}
    </div>
  )
}

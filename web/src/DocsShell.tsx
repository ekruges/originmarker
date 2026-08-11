import { useEffect, type ReactNode } from 'react'
import { Anchor, Button, Stack, Text, Title } from '@mantine/core'
import type { Health } from './api'

/**
 * The furniture every documentation page in this tool shares: the sticky numbered nav, the
 * article column, the deep-link behaviour, and the section headings.
 *
 * It exists because there are now three of these pages and the numbering, the scroll-margin and
 * the hash routing were being retyped each time. Section numbers come from the section list, so
 * inserting a section renumbers every cross-reference to it rather than leaving one behind.
 */

export const REPO_URL = 'https://github.com/ekruges/originmarker'
export const HOME_URL = 'https://ezrakruger.cc/'

/** GitHub's mark, inlined. currentColor so one CSS rule drives its resting and hover grey. */
export const GithubMark = () => (
  <svg viewBox="0 0 16 16" width="17" height="17" fill="currentColor" aria-hidden focusable="false">
    <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38
      0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01
      1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95
      0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27
      2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82
      2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0
      .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z" />
  </svg>
)

/** Person in a circle: the link home. Same 16-unit box and currentColor as the mark beside
 *  it, so one CSS rule greys both. The shoulders are clipped by the ring rather than drawn
 *  to meet it, which is what keeps the silhouette reading at 17px. */
export const AvatarMark = () => (
  <svg viewBox="0 0 16 16" width="17" height="17" aria-hidden focusable="false">
    <defs>
      <clipPath id="om-av-clip">
        <circle cx="8" cy="8" r="7.25" />
      </clipPath>
    </defs>
    <circle cx="8" cy="8" r="7.25" fill="none" stroke="currentColor" strokeWidth="1.5" />
    <g clipPath="url(#om-av-clip)" fill="currentColor">
      <circle cx="8" cy="6.1" r="2.6" />
      <ellipse cx="8" cy="14.6" rx="4.7" ry="3.6" />
    </g>
  </svg>
)

export interface DocSection {
  id: string
  label: string
}

export interface DocSibling {
  label: string
  href: string
}

/**
 * The other two documentation pages, as buttons.
 *
 * Named and nothing more. They used to carry a sentence each explaining what the other page
 * covered, which is the one thing a reader who is already inside the documentation does not need
 * spelled out, and it buried the link it was describing.
 */
export function DocsSiblingLinks({ siblings }: { siblings: DocSibling[] }) {
  if (!siblings.length) return null
  return (
    <Stack gap={8} mt={14}>
      {siblings.map((s) => (
        <Button
          key={s.href} component="a" href={s.href} size="md" radius={2} fullWidth
          styles={{ label: { whiteSpace: 'normal', lineHeight: 1.25 } }}
          style={{ fontSize: 12, fontWeight: 700, height: 'auto', padding: '10px 12px' }}
        >
          {s.label}
        </Button>
      ))}
    </Stack>
  )
}

/** The section id in a hash like `#/progenitor-docs/membership`, or '' for the bare route. */
export function docsSectionFromHash(prefix: string, hash: string): string {
  const m = new RegExp(`^#/${prefix}/([\\w:.-]+)$`).exec(hash)
  return m ? m[1] : ''
}

/** Scroll a section into view. A frame late, because on a cold deep link the section may not be
 *  mounted yet. */
export function jumpToSection(id: string): void {
  requestAnimationFrame(() => {
    document.getElementById(id)?.scrollIntoView({ block: 'start', behavior: 'smooth' })
  })
}

/**
 * The `Section` and `SecRef` a page uses, bound to its own route and its own section list.
 *
 * Returned from a call rather than exported directly because the numbering is per-page: section 7
 * of the Syngamy docs and section 7 of these are different sections.
 */
export function docsHelpers(prefix: string, sections: DocSection[]) {
  const href = (id: string): string => `#/${prefix}/${id}`
  const no = (id: string): number => sections.findIndex((s) => s.id === id) + 1

  const Section = ({ id, title, children }: {
    id: string; title: string; children: ReactNode
  }) => (
    <section id={id} style={{ scrollMarginTop: 12, marginBottom: 22 }}>
      <Title order={2} mb={6} pb={3} style={{ borderBottom: '1px solid var(--om-border)' }}>
        {no(id)} · {title}
      </Title>
      {children}
    </section>
  )

  const SecRef = ({ id }: { id: string }) => (
    <Anchor href={href(id)}>section {no(id)}</Anchor>
  )

  return { Section, SecRef, href, no }
}

export function DocsShell({ prefix, sections, title, subtitle, health, siblings, children }: {
  prefix: string
  sections: DocSection[]
  title: string
  subtitle: ReactNode
  health: Health | null
  /** The other documentation pages. Each is a button; the name is the whole of it. */
  siblings: DocSibling[]
  children: ReactNode
}) {
  useEffect(() => {
    const jump = (): void => {
      const id = docsSectionFromHash(prefix, window.location.hash)
      if (id) jumpToSection(id)
    }
    jump()
    window.addEventListener('hashchange', jump)
    return () => window.removeEventListener('hashchange', jump)
  }, [prefix])

  return (
    <div
      className="om-docs-wrap"
      style={{ display: 'flex', gap: 24, margin: '0 auto', padding: 12, alignItems: 'flex-start' }}
    >
      <nav
        className="om-docs-nav"
        aria-label={`${title} sections`}
        style={{ position: 'sticky', top: 12, flex: '0 0 200px', alignSelf: 'flex-start' }}
      >
        <ol style={{ listStyle: 'none', margin: 0, padding: 0 }}>
          {sections.map((s, i) => (
            <li key={s.id}>
              <a href={`#/${prefix}/${s.id}`}>
                <span className="om-mono" style={{ marginRight: 6 }}>{i + 1}</span>
                {s.label}
              </a>
            </li>
          ))}
        </ol>
        <Text size="xs" c="dimmed" mt={10} pl={8} className="om-mono">
          {health ? `${health.version} · ${health.release_codename}` : 'browser-only'} · in-tab
        </Text>
        <div className="om-docs-links">
          <a
            href={REPO_URL} target="_blank" rel="noreferrer"
            aria-label="Source on GitHub" title="Source on GitHub"
          >
            <GithubMark />
          </a>
          <a href={HOME_URL} aria-label="ezrakruger.cc" title="ezrakruger.cc">
            <AvatarMark />
          </a>
        </div>
        <DocsSiblingLinks siblings={siblings} />
      </nav>

      <article className="om-docs-body" style={{ flex: 1, minWidth: 0 }}>
        <Title order={1} mb={4}>{title}</Title>
        <Text size="xs" c="dimmed" mb="md">{subtitle}</Text>
        {children}

        <div
          style={{
            borderTop: '1px solid var(--om-border)',
            marginTop: 28,
            paddingTop: 10,
            fontSize: 11,
            color: 'var(--om-text-dim)',
            lineHeight: 1.6,
          }}
        >
          Developed by &amp; for the{' '}
          <a href="https://eglilab.com" target="_blank" rel="noopener noreferrer">Egli Lab</a>
          {' '}at Columbia University Irving Medical Center, and the Columbia Stem Cell
          Initiative. Research use only: candidate markers and parental-origin calls require
          validation in a qualified genetics laboratory. Not a clinical diagnostic.
        </div>
      </article>
    </div>
  )
}

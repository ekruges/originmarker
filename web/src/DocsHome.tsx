import { Anchor, Text, Title } from '@mantine/core'
import { REPO_URL, type DocSection } from './DocsShell'
import { SECTIONS as PANEL } from './DocsPage'
import { SECTIONS as SYNGAMY } from './SyngamyDocs'
import { SECTIONS as PROGENITOR } from './ProgenitorDocs'
import type { Health } from './api'

/**
 * The documentation index: one entry per manual, and each manual's contents by part.
 *
 * Built from the manuals' own section lists, so a section added, renamed or moved on a manual
 * appears here with the right number and nothing on this page has to be edited to follow it.
 * `#/docs/<section>` still opens the panel builder manual at that section.
 */

interface Manual {
  name: string
  prefix: string
  summary: string
  sections: DocSection[]
}

const MANUALS: Manual[] = [
  {
    name: 'Panel builder',
    prefix: 'docs',
    summary: 'Builds a ranked list of SNP markers around a pathogenic variant, for PGT-M linkage '
      + 'testing.',
    sections: PANEL,
  },
  {
    name: 'Syngamy',
    prefix: 'syngamy-docs',
    summary: 'Reads SNP arrays from an embryo and one or both parents. Reports whole-chromosome and '
      + 'segmental changes, and which parent each came from.',
    sections: SYNGAMY,
  },
  {
    name: 'Progenitor',
    prefix: 'progenitor-docs',
    summary: 'Reconstructs a parent’s genotype from that parent’s haploid products, such as '
      + 'sperm or pronuclei, for use as the parental array in Syngamy.',
    sections: PROGENITOR,
  },
]

/** Consecutive sections sharing a group, in order, each keeping its manual-wide number. */
function parts(sections: DocSection[]): { name: string; items: { s: DocSection; n: number }[] }[] {
  const out: { name: string; items: { s: DocSection; n: number }[] }[] = []
  sections.forEach((s, i) => {
    const name = s.group ?? ''
    const last = out[out.length - 1]
    if (last && last.name === name) last.items.push({ s, n: i + 1 })
    else out.push({ name, items: [{ s, n: i + 1 }] })
  })
  return out
}

export function DocsHome({ health }: { health: Health | null }) {
  return (
    <div style={{ maxWidth: 920, margin: '0 auto', padding: 12 }}>
      <Title order={1} mb={4}>Documentation</Title>
      <Text size="sm" c="dimmed" mb={20}>
        Reference manuals for the three parts of OriginMarker. Section numbers match the numbers in
        each manual.
      </Text>

      {MANUALS.map((m) => (
        <section key={m.prefix} style={{ marginBottom: 24 }}>
          <Title order={2} mb={2}>
            <Anchor href={`#/${m.prefix}/${m.sections[0].id}`} inherit>{m.name}</Anchor>
          </Title>
          <Text size="sm" mb={8}>{m.summary}</Text>
          <table className="om-doc-index">
            <tbody>
              {parts(m.sections).map((p) => (
                <tr key={p.name}>
                  <th scope="row">{p.name}</th>
                  <td>
                    {p.items.map(({ s, n }) => (
                      <a key={s.id} href={`#/${m.prefix}/${s.id}`}>
                        <span className="om-mono">{n}</span> {s.label}
                      </a>
                    ))}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </section>
      ))}

      <section style={{ marginBottom: 24 }}>
        <Title order={2} mb={6}>Other</Title>
        <Text size="sm" mb={2}>
          <Anchor href={`${REPO_URL}/blob/main/docs/validation.md`} target="_blank" rel="noreferrer">
            Validation record
          </Anchor>
          : what the tool has been tested against, and the results.
        </Text>
        <Text size="sm" mb={2}>
          <Anchor href="#/terms">Terms</Anchor>: data handling and conditions of use.
        </Text>
        <Text size="sm">
          <Anchor href={REPO_URL} target="_blank" rel="noreferrer">Source code</Anchor> on GitHub.
        </Text>
      </section>

      <Text size="xs" c="dimmed" pt={10} style={{ borderTop: '1px solid var(--om-border)' }}>
        <span className="om-mono">
          {health ? `${health.version} (${health.release_codename})` : 'Version unavailable'}
        </span>
        . Developed for the{' '}
        <a href="https://eglilab.com" target="_blank" rel="noopener noreferrer">Egli Lab</a>
        {' '}at Columbia University Irving Medical Center and the Columbia Stem Cell Initiative. For
        research use only. Results require confirmation by a qualified genetics laboratory. Not a
        clinical diagnostic.
      </Text>
    </div>
  )
}

import { Table, Text } from '@mantine/core'
import { HET_BAND_DIPLOID, pct } from './parentage'
import { CALL_RATE_FLOOR } from './parentage'
import { MIN_PRODUCTS, SAME_PARENT_MAX } from './inferredReference'
import type { SampleProfile } from './ingest'

/**
 * The two screens that come before any reference exists: quality and ploidy triage, then the
 * concordance matrix that decides which products share a parent.
 *
 * The order is a safety property rather than tidiness. Running concordance before the
 * call-rate gate splits one father into several, because an array below the floor reads high
 * against everyone: a true product at a 53.9% call rate reads 12.7% to 14.4% against its own
 * co-products, which the same-parent threshold calls a different man. Running it before the
 * ploidy gate merges unrelated fathers, because a diploid compared with a haploid has a
 * different expected rate entirely and bridges the groups.
 */

/**
 * Three bands, not two. The cut between haploid and diploid is not clean: confirmed products
 * run 1.0% to 7.8% and arrays excluded as diploid measured 10.0% and 10.7%.
 *
 * Two signals, because neither alone is sufficient. A known unrelated diploid adult measures a
 * 7.17% BAF band, under the 8% gate, and passes on the band alone. It cost nothing in theory
 * and everything in practice: admitted as a product it read 10.9% against a genuine product,
 * which fragmented one father's group into two and produced three groups where there are two.
 *
 * Genotype heterozygosity catches it. A haploid genome cannot be heterozygous at all, so every
 * such call is error: twelve confirmed products ran 5.3% to 13.7% and that adult reads 15.0%.
 * The margin is 1.3 points, which is thin, so a sample that fails EITHER signal is withheld
 * rather than argued about.
 */
export type Ploidy = 'haploid' | 'borderline' | 'diploid'

/** Above this, genotype heterozygosity is diploid rather than a haploid's error rate. */
export const HET_RATE_DIPLOID = 0.14

export const ploidyOf = (band: number, hetRate: number): Ploidy => {
  if (!Number.isFinite(band)) return 'borderline'
  if (band >= 0.15 || hetRate >= 0.17) return 'diploid'
  if (band > HET_BAND_DIPLOID || hetRate >= HET_RATE_DIPLOID) return 'borderline'
  return 'haploid'
}

export interface Triage {
  id: string
  name: string
  profile: SampleProfile
  ploidy: Ploidy
  usable: boolean
  why: string
}

export function triage(id: string, name: string, p: SampleProfile): Triage {
  const ploidy = ploidyOf(p.hetBand, p.hetRate)
  if (p.callRate < CALL_RATE_FLOOR) {
    return {
      id,
      name,
      profile: p,
      ploidy,
      usable: false,
      why: `call rate ${pct(p.callRate, 1)} is below the ${pct(CALL_RATE_FLOOR, 0)} floor`,
    }
  }
  if (ploidy !== 'haploid') {
    return {
      id,
      name,
      profile: p,
      ploidy,
      usable: false,
      why: ploidy === 'diploid'
        ? `BAF band ${pct(p.hetBand, 1)} and ${pct(p.hetRate, 1)} heterozygous calls read `
          + 'diploid, not one meiotic product'
        : `BAF band ${pct(p.hetBand, 1)} with ${pct(p.hetRate, 1)} heterozygous calls sits `
          + 'between the haploid and diploid ranges',
    }
  }
  return { id, name, profile: p, ploidy, usable: true, why: '' }
}

const CELL: React.CSSProperties = {
  fontFamily: 'var(--om-mono)', fontSize: 11.5, textAlign: 'right', padding: '3px 9px',
  fontVariantNumeric: 'tabular-nums',
}

export function TriageTable({ rows }: { rows: Triage[] }) {
  if (!rows.length) return null
  const usable = rows.filter((r) => r.usable)
  return (
    <div>
      <Table withTableBorder style={{ fontSize: 11.5 }}>
        <Table.Thead>
          <Table.Tr>
            <Table.Th>File</Table.Th>
            <Table.Th style={{ textAlign: 'right' }}>Call rate</Table.Th>
            <Table.Th style={{ textAlign: 'right' }}>BAF band</Table.Th>
            <Table.Th>Reading</Table.Th>
          </Table.Tr>
        </Table.Thead>
        <Table.Tbody>
          {rows.map((r) => (
            <Table.Tr
              key={r.id}
              style={{ background: r.usable ? undefined : 'var(--om-zebra)' }}
            >
              <Table.Td style={{
                fontFamily: 'var(--om-mono)', fontSize: 11.5, fontWeight: r.usable ? 600 : 400,
                color: r.usable ? undefined : 'var(--om-text-dim)',
                borderLeft: `3px solid ${r.usable ? 'var(--om-blue)' : 'transparent'}`,
              }}
              >
                {r.name}
              </Table.Td>
              <Table.Td style={CELL}>{pct(r.profile.callRate, 1)}</Table.Td>
              <Table.Td style={CELL}>{pct(r.profile.hetBand, 2)}</Table.Td>
              <Table.Td style={{ fontSize: 11.5 }}>
                {r.usable
                  ? <span style={{ color: 'var(--om-blue)', fontWeight: 600 }}>usable product</span>
                  : <span style={{ color: 'var(--om-text-dim)' }}>{r.why}</span>}
              </Table.Td>
            </Table.Tr>
          ))}
        </Table.Tbody>
      </Table>
      {usable.length < MIN_PRODUCTS && (
        <Text size="xs" mt={4} style={{ color: 'var(--om-higher)' }}>
          Below the {MIN_PRODUCTS} products this method needs, so no reference is built.
        </Text>
      )}
    </div>
  )
}

/** Shading on the rate alone. A per-pair label is deliberately not the visual: measured over 91
 *  pairs, one genuine different-father pair reads 9.88% and would show as same-parent. */
function shade(rate: number): string {
  if (!Number.isFinite(rate)) return 'transparent'
  const t = Math.max(0, Math.min(1, (rate - 0.04) / 0.13))
  return `rgba(46, 109, 164, ${(0.45 * (1 - t) + 0.03).toFixed(3)})`
}

/**
 * The rule this matrix exists to show is "same parent or not", so the stripe has to separate
 * every group from every other, not just the first from the rest. `groupByParent` returns an
 * unbounded number of groups and a two-colour rule made groups 2 and 3 identical on the one
 * display that decides which products belong to which person.
 */
const GROUP_COLOURS = [
  'var(--om-blue)', 'var(--om-higher)', 'var(--om-text-dim)', 'var(--om-blue-light)',
]
const groupColour = (g: number): string => GROUP_COLOURS[g % GROUP_COLOURS.length]

export function ConcordanceMatrix({ names, rate, groups }: {
  names: string[]
  rate: (a: number, b: number) => number
  groups: number[][]
}) {
  if (names.length < 2) return null
  const order = groups.flat()
  const groupOf = new Map<number, number>()
  groups.forEach((g, i) => g.forEach((x) => groupOf.set(x, i)))

  const stat = (g: number[]) => {
    let worst = 0
    for (const a of g) for (const b of g) if (a < b) worst = Math.max(worst, rate(a, b))
    return worst
  }
  let closestOut = Infinity
  for (const a of order) {
    for (const b of order) {
      if (a < b && groupOf.get(a) !== groupOf.get(b)) closestOut = Math.min(closestOut, rate(a, b))
    }
  }

  return (
    <div style={{ marginTop: 14 }}>
      <div style={{ overflowX: 'auto' }}>
        <table style={{ borderCollapse: 'collapse', fontSize: 10 }}>
          <thead>
            <tr>
              <th />
              {order.map((c) => (
                <th key={c} style={{ ...CELL, fontWeight: 500, color: 'var(--om-text-dim)' }}>
                  {names[c]}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {order.map((r) => (
              <tr key={r}>
                <td style={{
                  ...CELL, fontWeight: 500, textAlign: 'left', whiteSpace: 'nowrap',
                  borderLeft: `3px solid ${groupColour(groupOf.get(r) ?? 0)}`,
                }}
                >
                  {names[r]}
                </td>
                {order.map((c) => (
                  <td
                    key={c}
                    style={{
                      ...CELL,
                      background: r === c ? 'transparent' : shade(rate(r, c)),
                      color: r === c ? 'var(--om-text-dim)' : undefined,
                    }}
                  >
                    {r === c ? '.' : (rate(r, c) * 100).toFixed(1)}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {groups.length > 1 ? (
        <div style={{
          marginTop: 10, padding: '8px 10px', border: '1px solid var(--om-border)',
          background: 'var(--om-warn-bg)',
        }}
        >
          <Text size="xs">
            {groups.length} groups are internally consistent and inconsistent with each other.
            {groups.map((g, i) => (
              <span key={i}>
                {' '}Group {i + 1}: {g.map((x) => names[x]).join(', ')}, weakest pair inside{' '}
                {pct(stat(g), 2)}.
              </span>
            ))}
            {' '}The closest pair across groups is {pct(closestOut, 2)}.
          </Text>
          <Text size="xs" mt={6}>
            No reference is built until you say which group you mean. One product of another man
            among five destroyed a reference we measured: every member then failed its own
            leave-one-out check, including the two cleanest.
          </Text>
        </div>
      ) : (
        <Text size="xs" mt={8}>
          One group, weakest pair inside {pct(stat(groups[0] ?? []), 2)}. Consistent with a
          single parent.
        </Text>
      )}
      <Text size="xs" c="dimmed" mt={6}>
        Membership requires every pair inside a group to read same parent. A chain of pairwise
        links is not enough: one genuine different-father pair measured {pct(0.0988, 2)}, under
        the {pct(SAME_PARENT_MAX, 1)} cut, and grouping by chained links merges two men through
        that single edge. Requiring all pairs rejects it, because the same product reads 12.1% to
        12.4% against the rest of that group.
      </Text>
    </div>
  )
}


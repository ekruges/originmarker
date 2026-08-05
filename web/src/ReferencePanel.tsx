import { Text } from '@mantine/core'
import { pct } from './parentage'
import { MIN_ASCERTAINMENT } from './inferredReference'

/**
 * The reference itself: how deep it is, how far its marker set has drifted from the genome, and
 * whether each product it was built from still reads present when built without it.
 *
 * A reader has to be able to tell an inferred reference from a measured one at a glance and
 * judge it without leaving the page, so every quantity the call depends on is on screen rather
 * than in an export.
 */

const MONO = 'var(--om-mono)'

const num: React.CSSProperties = {
  fontFamily: MONO, fontSize: 11.5, textAlign: 'right', padding: '4px 10px',
  whiteSpace: 'nowrap', fontVariantNumeric: 'tabular-nums',
}
const head: React.CSSProperties = {
  fontSize: 10, textTransform: 'uppercase', letterSpacing: '0.04em',
  color: 'var(--om-text-dim)', fontWeight: 600, padding: '3px 8px', textAlign: 'right',
}

/** One measured quantity, label under value, so a row of them reads as a set. */
function Stat({ label, value, note }: { label: string; value: string; note?: string }) {
  return (
    <div style={{ minWidth: 118 }}>
      <div style={{
        fontFamily: MONO, fontSize: 17, lineHeight: 1.2, fontVariantNumeric: 'tabular-nums',
        letterSpacing: '-0.01em',
      }}
      >
        {value}
      </div>
      <div style={{ fontSize: 10, color: 'var(--om-text-dim)', letterSpacing: '0.02em' }}>
        {label}
      </div>
      {note && <div style={{ fontSize: 10, color: 'var(--om-text-dim)' }}>{note}</div>}
    </div>
  )
}

/**
 * The ascertainment ladder, which is how m is chosen.
 *
 * Raising m lowers contamination and narrows the marker set toward probes nearly every product
 * called. Those are enriched for sites where the parent is homozygous because the minor allele
 * is rare, and unrelated people carry the common allele there too, so past a point the
 * reference stops excluding them. The rule takes the deepest setting still holding 90%.
 */
export function AscertainmentLadder({ ratios, chosen }: {
  ratios: Record<string, number>
  chosen: number
}) {
  const entries = Object.entries(ratios).map(([m, r]) => [Number(m), r] as const)
    .sort((a, b) => a[0] - b[0])
  return (
    <div style={{ marginTop: 10 }}>
      <Text size="xs" fw={600} mb={2}>How deep the reference goes, and why</Text>
      <Text size="xs" c="dimmed" mb={6}>
        Each row requires that many products to have spoken at a marker and agreed. Deeper is
        cleaner but narrower, and a narrow marker set loses its grip on unrelated genomes. The
        rule takes the deepest row still holding {pct(MIN_ASCERTAINMENT, 0)}.
      </Text>
      <table style={{ borderCollapse: 'collapse' }}>
        <tbody>
          {entries.map(([m, r]) => {
            const pick = m === chosen
            const ok = r >= MIN_ASCERTAINMENT
            return (
              <tr key={m} style={{ opacity: ok ? 1 : 0.5 }}>
                <td style={{ ...num, fontWeight: pick ? 700 : 400, width: 52 }}>m &ge; {m}</td>
                <td style={{ padding: '3px 8px', width: 240 }}>
                  <div style={{ position: 'relative', height: 9, background: 'var(--om-border)' }}>
                    <div style={{
                      position: 'absolute', inset: 0, width: `${Math.max(0, r) * 100}%`,
                      background: ok ? 'var(--om-blue)' : 'var(--om-text-dim)',
                    }}
                    />
                    <div style={{
                      position: 'absolute', top: -2, bottom: -2,
                      left: `${MIN_ASCERTAINMENT * 100}%`, width: 1, background: 'var(--om-higher)',
                    }}
                    />
                  </div>
                </td>
                <td style={{ ...num, fontWeight: pick ? 700 : 400 }}>{pct(r, 0)}</td>
                <td style={{ fontSize: 10, padding: '3px 8px', color: 'var(--om-text-dim)' }}>
                  {pick ? 'chosen' : ok ? '' : 'too narrow'}
                </td>
              </tr>
            )
          })}
        </tbody>
      </table>
      <Text size="xs" c="dimmed" mt={4}>
        The mark is the {pct(MIN_ASCERTAINMENT, 0)} floor. An earlier rule fixed m at one below
        the product count; it is right at five products and wrong at eight, where it demanded
        seven and left the marker set at roughly 62% of the genome's heterozygosity.
      </Text>
    </div>
  )
}

export interface Scored {
  name: string
  absence: number
  ceiling?: number
  ratio: number
  verdict: string
  role?: string
}

const VERDICT: Record<string, { label: string; color: string }> = {
  parent_genome_present: { label: 'present', color: 'var(--om-blue)' },
  no_parental_contribution: { label: 'absent', color: 'var(--om-text-dim)' },
  unclear: { label: 'unclear', color: 'var(--om-higher)' },
}

function ScoreTable({ rows, caption }: { rows: Scored[]; caption: string }) {
  if (!rows.length) return null
  const hasCeiling = rows.some((r) => r.ceiling !== undefined)
  const hasRole = rows.some((r) => r.role !== undefined)
  return (
    <div style={{ marginTop: 12 }}>
      <Text size="xs" fw={600} mb={4}>{caption}</Text>
      <table style={{ borderCollapse: 'collapse', width: '100%', maxWidth: 620 }}>
        <thead>
          <tr style={{ borderBottom: '1px solid var(--om-border)' }}>
            <th style={{ ...head, textAlign: 'left' }}>Sample</th>
            {hasRole && <th style={{ ...head, textAlign: 'left' }}>Role</th>}
            <th style={head}>Absent</th>
            {hasCeiling && <th style={head}>Ceiling</th>}
            <th style={head}>Ratio</th>
            <th style={{ ...head, textAlign: 'left' }}>Reads</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((r) => {
            const v = VERDICT[r.verdict] ?? { label: r.verdict, color: 'inherit' }
            return (
              <tr key={r.name} style={{ borderBottom: '1px solid var(--om-border)' }}>
                <td style={{ ...num, textAlign: 'left', fontWeight: 500 }}>{r.name}</td>
                {hasRole && (
                  <td style={{ fontSize: 10, padding: '3px 8px', color: 'var(--om-text-dim)' }}>
                    {r.role ?? ''}
                  </td>
                )}
                <td style={num}>{pct(r.absence, 2)}</td>
                {hasCeiling && <td style={num}>{r.ceiling === undefined ? '' : pct(r.ceiling, 2)}</td>}
                <td style={{ ...num, fontWeight: 600 }}>{r.ratio.toFixed(2)}x</td>
                <td style={{ fontSize: 11, padding: '3px 8px', color: v.color, fontWeight: 500 }}>
                  {v.label}
                </td>
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}

export function ReferenceBlock({ group, mMin, markers, meanM, hRetained, contamination,
  spuriousAbsence, ratios, members, controls }: {
  group: string[]
  mMin: number
  markers: number
  meanM: number
  hRetained: number
  contamination: number
  spuriousAbsence: number
  ratios: Record<string, number>
  members: Scored[]
  controls: Scored[]
}) {
  const held = members.filter((m) => m.verdict === 'parent_genome_present').length
  return (
    <div>
      <div style={{
        display: 'flex', gap: 8, flexWrap: 'wrap', marginBottom: 6,
        border: '1px solid var(--om-border)', padding: '10px 14px',
      }}
      >
        <Stat label="products" value={String(group.length)} note={group.join(' ')} />
        <Stat label="agreeing calls required" value={`≥ ${mMin}`} note={`mean ${meanM.toFixed(2)}`} />
        <Stat label="markers" value={markers.toLocaleString()} />
        <Stat label="ascertainment" value={pct(hRetained, 1)} note="parent heterozygosity" />
        <Stat label="contamination" value={pct(contamination, 2)} note="het read as hom" />
        <Stat label="added absence" value={pct(spuriousAbsence, 2)} note="to a haploid product" />
      </div>

      <AscertainmentLadder ratios={ratios} chosen={mMin} />

      <ScoreTable
        caption={`Every product, scored against a reference built without it (${held} of ${members.length} present)`}
        rows={members}
      />
      <Text size="xs" c="dimmed" mt={2}>
        A product scored against a reference containing itself reads exactly zero absence, which
        is a bias the size of the whole signal, so each one is judged by the other
        {' '}{group.length - 1}. That reference is built from one fewer product and is dirtier
        for it, so each row&rsquo;s ceiling carries that reference&rsquo;s own added absence
        rather than the full build&rsquo;s.
      </Text>

      <ScoreTable caption="Controls, none of them a member" rows={controls} />
    </div>
  )
}

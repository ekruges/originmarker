/**
 * The optional feature comparison, as a panel a reader chooses to open.
 *
 * OPTIONAL BY DESIGN, not by omission. This asks whether a change sits where the genome breaks
 * anyway, which is a different question from whose the change is, and presenting the two together
 * invites a reader to treat the first as evidence about the second. It is not: the fragile
 * compartment is established on both parental genomes from the first cell cycle. So it stays shut
 * until someone opens it, and its result carries that sentence wherever it goes.
 *
 * The charts are deliberately plain. Three of them, one shared axis, no colour carrying meaning on
 * its own, because the reader's question is "do these overlap or not" and anything more elaborate
 * answers a question nobody asked.
 */
import { useState } from 'react'
import { Text } from '@mantine/core'
import {
  compare, enrichmentBars, foldChart, regionGrid, nullHistograms,
  type ComparisonResult,
} from './comparison.ts'
import type { FeatureTrack, Region } from './features.ts'
import { comparisonPdf } from './comparisonPdf.ts'

const INK = 'var(--om-text)'
const DIM = 'var(--om-text-dim)'
const HIT = 'var(--om-defect)'

/**
 * Observed against its own null, one row per feature, on ONE shared axis.
 *
 * The bar is the middle 95% of the matched null and the dot is what was observed. A permutation p
 * is uninterpretable without the distribution it came from, so the distribution is drawn rather
 * than summarised: a ratio near 1 can reach a small p when the null is tight, and a reader shown
 * only the p would call that a finding.
 */
function EnrichmentChart({ c }: { c: ComparisonResult }) {
  const bars = enrichmentBars(c)
  if (!bars.length) return null
  const W = 260
  const rowH = 26
  const x = (v: number) => (v / bars[0].axisMax) * W
  return (
    <svg width="100%" viewBox={`0 0 ${W + 230} ${bars.length * rowH + 26}`} role="img"
      aria-label="Observed overlap against the matched null, one row per feature"
    >
      {bars.map((b, i) => {
        const y = i * rowH + 16
        return (
          <g key={b.label}>
            <text x={0} y={y + 4} fontSize={9} fill={DIM}>{b.label}</text>
            <line x1={168} x2={168 + W} y1={y} y2={y} stroke="var(--om-line)" strokeWidth={1} />
            <line
              x1={168 + x(b.lo)} x2={168 + x(b.hi)} y1={y} y2={y}
              stroke={DIM} strokeWidth={6}
            />
            <circle cx={168 + x(b.observed)} cy={y} r={4}
              fill={b.significant ? HIT : INK}
            />
            <text x={168 + W + 6} y={y + 4} fontSize={9}
              fill={b.significant ? HIT : DIM}
            >
              {`${(100 * b.observed).toFixed(0)}% vs ${(100 * b.expected).toFixed(0)}%`}
              {b.significant ? `  p=${b.p.toExponential(1)}` : ''}
            </text>
          </g>
        )
      })}
      <text x={168} y={bars.length * rowH + 20} fontSize={8} fill={DIM}>
        bar = middle 95% of the matched null, dot = observed
      </text>
    </svg>
  )
}

/**
 * The permutation distribution itself, with the observation marked in it.
 *
 * THE MOST DEFENSIBLE FIGURE HERE, and it costs nothing because the enrichment already computed it.
 * A p value is a summary of this picture, and summaries of it can mislead in both directions: a
 * ratio near one reaches a small p when the null is tight, and a large ratio reaches nothing when
 * the null is broad. Drawing the null lets a reader make that judgement instead of taking the p on
 * trust. An observation the permutation never reached is drawn at the edge and said to be outside.
 */
function NullHistograms({ c }: { c: ComparisonResult }) {
  const hs = nullHistograms(c)
  if (!hs.length) return null
  const W = 150
  const H = 34
  return (
    <div style={{ display: 'flex', flexWrap: 'wrap', gap: 14 }}>
      {hs.map((h) => {
        const max = Math.max(1, ...h.bins)
        const bw = W / h.bins.length
        return (
          <div key={h.label} style={{ width: W + 12 }}>
            <div style={{ fontSize: 9, color: DIM, marginBottom: 2 }}>{h.label}</div>
            <svg width={W} height={H + 12} role="img"
              aria-label={`Permutation null for ${h.label}, with the observed value marked`}
            >
              {h.bins.map((v, i) => (
                <rect
                  key={i} x={i * bw} y={H - (v / max) * H}
                  width={Math.max(0.7, bw - 0.6)} height={(v / max) * H}
                  fill={i === h.observedBin ? HIT : 'var(--om-line)'}
                />
              ))}
              <line x1={h.observedAt * W} x2={h.observedAt * W} y1={0} y2={H}
                stroke={h.significant ? HIT : INK} strokeWidth={1.4}
              />
              <text x={0} y={H + 10} fontSize={7.5} fill={DIM}>
                {h.observedBin === -1 ? 'observed OUTSIDE the null' : `p=${h.p.toExponential(1)}`}
              </text>
            </svg>
          </div>
        )
      })}
    </div>
  )
}

/** Fold enrichment, which puts every feature on one scale so they can be ranked at a glance. */
function FoldChart({ c }: { c: ComparisonResult }) {
  const rows = foldChart(c)
  if (!rows.length) return null
  const max = Math.max(2, ...rows.map((r) => r.fold))
  const W = 260
  const rowH = 20
  return (
    <svg width="100%" viewBox={`0 0 ${W + 230} ${rows.length * rowH + 26}`} role="img"
      aria-label="Fold enrichment per feature against the matched null"
    >
      <line x1={168 + (1 / max) * W} x2={168 + (1 / max) * W} y1={8}
        y2={rows.length * rowH + 4} stroke="var(--om-line)" strokeDasharray="3 3"
      />
      {rows.map((r, i) => {
        const y = i * rowH + 16
        return (
          <g key={r.label}>
            <text x={0} y={y + 3} fontSize={9} fill={DIM}>{r.label}</text>
            <rect x={168} y={y - 5} width={Math.max(1, (r.fold / max) * W)} height={10}
              fill={r.fold >= 1 ? INK : DIM} opacity={0.8}
            />
            <text x={168 + W + 6} y={y + 3} fontSize={9} fill={DIM}>
              {`${r.fold.toFixed(2)}x`}
            </text>
          </g>
        )
      })}
      <text x={168} y={rows.length * rowH + 20} fontSize={8} fill={DIM}>
        dashed line = the null, 1.00x. Left of it is LESS overlap than chance
      </text>
    </svg>
  )
}

/**
 * Which region touches which feature, as a grid.
 *
 * The plainest answer to "overlap or lack thereof", and the one a reader checks the enrichment
 * against: a feature can clear its null on two regions out of twenty, and this is where that
 * becomes visible rather than hiding inside a p value.
 */
function RegionGrid({ c, names }: { c: ComparisonResult, names: string[] }) {
  if (!names.length || !c.features.length) return null
  const cells = regionGrid(c, names)
  return (
    <div style={{ overflowX: 'auto' }}>
      <table style={{ borderCollapse: 'collapse', fontSize: 10 }}>
        <thead>
          <tr>
            <th style={{ textAlign: 'left', padding: '2px 8px 2px 0', color: DIM }}>Region</th>
            {c.features.map((f) => (
              <th key={f.feature} style={{ padding: '2px 6px', color: DIM, fontWeight: 600,
                writingMode: 'vertical-rl', whiteSpace: 'nowrap', height: 86 }}
              >
                {f.label}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {names.map((n) => (
            <tr key={n}>
              <td style={{ padding: '2px 8px 2px 0', fontFamily: 'var(--om-mono)' }}>{n}</td>
              {c.features.map((f) => {
                const hit = cells.find((x) => x.region === n && x.feature === f.label)?.touched
                return (
                  <td key={f.feature} style={{ textAlign: 'center', padding: '2px 6px',
                    color: hit ? HIT : 'var(--om-line)', fontWeight: hit ? 700 : 400 }}
                  >
                    {hit ? '●' : '·'}
                  </td>
                )
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

export function ComparisonPanel({
  regions, regionNames, markerPositions, loadTrack, existing, onDone, sampleName, build,
}: {
  regions: Region[]
  regionNames: string[]
  markerPositions?: Map<string, number[]>
  loadTrack: () => Promise<FeatureTrack | null>
  existing?: ComparisonResult
  onDone: (c: ComparisonResult) => void
  sampleName: string
  build: string
}) {
  const [open, setOpen] = useState(false)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const c = existing

  const run = async () => {
    setBusy(true); setError(null)
    try {
      const track = await loadTrack()
      if (!track) { setError('the feature track could not be loaded'); return }
      if (!markerPositions?.size) {
        setError('this run kept no marker positions, so the matched null cannot be built')
        return
      }
      onDone(compare(track, regions, markerPositions, { regionNames }))
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e))
    } finally { setBusy(false) }
  }

  return (
    <div style={{ border: '1px solid var(--om-line)', margin: '10px 0 4px' }}>
      <button
        type="button"
        onClick={() => setOpen(!open)}
        style={{ width: '100%', textAlign: 'left', background: 'none', border: 'none',
          padding: '9px 14px', cursor: 'pointer', fontSize: 13, fontWeight: 700, color: INK }}
      >
        {open ? '▾' : '▸'} Feature comparison
        <span style={{ fontWeight: 400, color: DIM }}>
          {c ? `  ${c.verdict.replace(/-/g, ' ')}` : '  optional, does not run by default'}
        </span>
      </button>

      {open && (
        <div style={{ padding: '0 14px 12px' }}>
          {!c && (
            <>
              <Text size="sm" mb={8} style={{ color: DIM }}>
                Asks whether these regions sit where the genome breaks anyway: common fragile sites,
                long genes, assembly gaps and late-replication valleys. A change at one of those has
                a mundane explanation available that a change elsewhere does not. It is a prior on
                how to read the result, not a verdict on it.
              </Text>
              <button type="button" onClick={run} disabled={busy || !regions.length}
                style={{ padding: '5px 12px', fontSize: 12, cursor: 'pointer' }}
              >
                {busy ? 'Comparing...' : `Run comparison on ${regions.length} region`
                  + `${regions.length === 1 ? '' : 's'}`}
              </button>
              {!regions.length && (
                <Text size="xs" style={{ color: DIM, marginTop: 6 }}>
                  This run found no regions to compare.
                </Text>
              )}
            </>
          )}
          {error && (
            <Text size="xs" style={{ color: 'var(--om-higher)', marginTop: 6 }}>{error}</Text>
          )}

          {c && (
            <>
              <Text style={{ fontSize: 13, fontWeight: 700, lineHeight: 1.35, marginBottom: 8 }}>
                {c.headline}
              </Text>

              <Text size="xs" fw={700} mt={10} mb={2}>Overlap against a matched null</Text>
              <EnrichmentChart c={c} />

              <Text size="xs" fw={700} mt={12} mb={4}>
                The permutation null, with the observation marked
              </Text>
              <NullHistograms c={c} />

              <Text size="xs" fw={700} mt={10} mb={2}>Fold enrichment</Text>
              <FoldChart c={c} />

              <Text size="xs" fw={700} mt={10} mb={4}>Which region touches which feature</Text>
              <RegionGrid c={c} names={regionNames} />

              <Text size="xs" fw={700} mt={12} mb={2}>What each feature would mean</Text>
              {c.features.map((f) => (
                <div key={f.feature} style={{ fontSize: 11, marginBottom: 4, lineHeight: 1.45 }}>
                  <span style={{ fontWeight: 600 }}>{f.label}</span>
                  <span style={{ color: DIM }}>{`: ${f.means}`}</span>
                </div>
              ))}

              <Text size="xs" fw={700} mt={12} mb={2}>Methods</Text>
              <Text style={{ fontSize: 11, color: DIM, lineHeight: 1.5 }}>{c.methods}</Text>

              <Text style={{ fontSize: 11, color: INK, lineHeight: 1.5, marginTop: 10,
                borderLeft: '3px solid var(--om-line)', paddingLeft: 8 }}
              >
                {c.caveat}
              </Text>

              <button
                type="button"
                onClick={() => {
                  const blob = comparisonPdf(c, regionNames, sampleName, build)
                  const url = URL.createObjectURL(blob)
                  const a = document.createElement('a')
                  a.href = url
                  a.download = `${sampleName.replace(/\.[^.]+$/, '')}-feature-comparison.pdf`
                  a.click()
                  URL.revokeObjectURL(url)
                }}
                style={{ marginTop: 12, padding: '5px 12px', fontSize: 12, cursor: 'pointer' }}
              >
                Comparison report (PDF)
              </button>
              <Text size="xs" style={{ color: DIM, marginTop: 6 }}>
                This section is also folded into the main report, with its methods, once it has been
                run. It is absent from that report until then, rather than appearing empty.
              </Text>
            </>
          )}
        </div>
      )}
    </div>
  )
}

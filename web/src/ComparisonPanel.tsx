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
import { compare, type ComparisonResult } from './comparison.ts'
import { GenomeViewer } from './GenomeViewer.tsx'
import type { FeatureTrack, Region } from './features.ts'
import { comparisonPdf } from './comparisonPdf.ts'

const INK = 'var(--om-text)'
const DIM = 'var(--om-text-dim)'





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

              {/*
                THE MAP IS THE HEADLINE AND THE PANELS ARE NOT HERE. A reader meeting a result wants
                to see the genome and what happened to it; a reader checking one wants the four
                quantitative panels, and those are in the report where they can be read at length
                beside their methods. Putting all five on screen buries the one that answers "where"
                under four that answer "how much".
              */}
              <div style={{ marginTop: 12 }}>
                <GenomeViewer c={c} />
              </div>
              <Text size="xs" style={{ color: DIM, marginTop: 2 }}>
                The enrichment, the permutation nulls, the fold and the region-by-feature matrix are
                in the report, with the methods that produced them.
              </Text>

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

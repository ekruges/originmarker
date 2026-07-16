import { useCallback, useEffect, useRef, useState } from 'react'
import {
  Alert, Anchor, Badge, Button, Group, List, Loader, Paper, Progress,
  Select, Skeleton, Text,
} from '@mantine/core'
import {
  ANCESTRIES, api, ApiError, type Health, type NLResponse, type PanelResult,
  type ResolveResponse, type StructuredQuery,
} from './api'
import { int, orUnknown, utc } from './fmt'
import { SearchPanel } from './SearchPanel'
import { VariantCard } from './VariantCard'
import { LocusTrack } from './LocusTrack'
import { PanelTable } from './PanelTable'
import { DocsPage } from './DocsPage'
import { TermsPage } from './TermsPage'
import { Logo } from './Logo'

type Phase = 'idle' | 'resolving' | 'resolved' | 'building' | 'done' | 'error'

/** The current hash route. Hash, not path: the deployment has no server-side rewrite. */
function useHashRoute(): string {
  const [hash, setHash] = useState(() => window.location.hash || '#/')
  useEffect(() => {
    const on = () => setHash(window.location.hash || '#/')
    window.addEventListener('hashchange', on)
    return () => window.removeEventListener('hashchange', on)
  }, [])
  return hash
}

export default function App() {
  const [health, setHealth] = useState<Health | null>(null)
  const [phase, setPhase] = useState<Phase>('idle')
  const [error, setError] = useState<string | null>(null)
  const [query, setQuery] = useState<StructuredQuery | null>(null)
  // Held here, not in SearchPanel: that component unmounts when the landing view collapses
  // into the results view, and the model-chosen-variant caveat has to outlive it.
  const [nl, setNl] = useState<NLResponse | null>(null)
  const [resolved, setResolved] = useState<ResolveResponse | null>(null)
  const [result, setResult] = useState<PanelResult | null>(null)
  const [jobId, setJobId] = useState<string | null>(null)
  const [progress, setProgress] = useState({ stage: '', fraction: 0 })
  const [ancestry, setAncestry] = useState<string | null>(null)
  const esRef = useRef<EventSource | null>(null)
  const pollRef = useRef<ReturnType<typeof setInterval> | null>(null)

  const route = useHashRoute()
  const atDocs = route.startsWith('#/docs')
  const atTerms = route.startsWith('#/terms')
  // Keyed on having no data, not on phase === 'idle': the hero must stay mounted while a
  // resolve is in flight, or the layout tears down mid-click.
  const atHome = !atDocs && !atTerms && !resolved && !result

  const stopWatch = () => {
    esRef.current?.close()
    esRef.current = null
    if (pollRef.current) clearInterval(pollRef.current)
    pollRef.current = null
  }

  useEffect(() => {
    api.health().then(setHealth).catch(() => setHealth(null))
    return stopWatch
  }, [])

  /** Back to the landing view, with nothing left running or left over. */
  const clearAll = useCallback(() => {
    // stopWatch first: an in-flight build would otherwise land setResult() on the cleared page.
    stopWatch()
    setPhase('idle')
    setError(null)
    setQuery(null)
    setNl(null)
    setResolved(null)
    setResult(null)
    setJobId(null)
    setProgress({ stage: '', fraction: 0 })
    setAncestry(null)
    if (window.location.hash && window.location.hash !== '#/') window.location.hash = '#/'
  }, [])

  const doResolve = useCallback(async (q: StructuredQuery, parsedBy?: NLResponse) => {
    stopWatch()
    setPhase('resolving')
    setError(null)
    setResolved(null)
    setResult(null)
    setJobId(null)
    setQuery(q)
    // Set unconditionally: a typed identifier must clear the previous query's model caveat.
    setNl(parsedBy ?? null)
    setAncestry(q.ancestry ?? null)
    try {
      const r = await api.resolve(q.variant, q.build)
      setResolved(r)
      setPhase('resolved')
    } catch (e) {
      setError(e instanceof ApiError ? e.message : 'Could not resolve that variant.')
      setPhase('error')
    }
  }, [])

  const settle = useCallback((job: { status: string; result?: PanelResult; error?: string }) => {
    if (job.status === 'done' && job.result) {
      setResult(job.result)
      setAncestry(job.result.provenance.ancestry_rank ?? null)
      setPhase('done')
    } else {
      setError(job.error ?? 'The panel job finished without a result.')
      setPhase('error')
    }
  }, [])

  // `q` becomes the stored query: `query` must describe the panel actually on screen, or
  // the next build reverts to the previous parameters.
  const buildPanel = useCallback(async (q?: StructuredQuery) => {
    const use = q ?? query
    if (!use) return
    setQuery(use)
    setPhase('building')
    setError(null)
    setProgress({ stage: 'submitting job', fraction: 0.02 })

    let id: string
    try {
      id = (await api.panel(use)).job_id
    } catch (e) {
      setError(e instanceof ApiError ? e.message : 'Could not start the panel build.')
      setPhase('error')
      return
    }
    setJobId(id)

    stopWatch()
    const es = new EventSource(api.streamUrl(id))
    esRef.current = es

    es.addEventListener('progress', (e) => {
      try {
        const d = JSON.parse((e as MessageEvent).data)
        setProgress({ stage: d.stage ?? '', fraction: typeof d.fraction === 'number' ? d.fraction : 0 })
      } catch { /* a malformed frame is not fatal; the next one lands shortly */ }
    })

    es.addEventListener('done', () => {
      stopWatch()
      api.job(id).then(settle).catch((e) => {
        setError(e instanceof ApiError ? e.message : 'Could not fetch the finished panel.')
        setPhase('error')
      })
    })

    es.addEventListener('error', (e) => {
      const data = (e as MessageEvent).data
      if (!data) return // no payload => transport hiccup; EventSource retries on its own
      stopWatch()
      try {
        setError(JSON.parse(data).message ?? 'The panel build failed.')
      } catch {
        setError('The panel build failed.')
      }
      setPhase('error')
    })

    // Fallback poller: covers a proxy that buffers SSE, or a dropped stream.
    pollRef.current = setInterval(async () => {
      try {
        const job = await api.job(id)
        if (job.status !== 'running') {
          stopWatch()
          settle(job)
        }
      } catch { /* transient failure; keep polling */ }
    }, 3000)
  }, [query, settle])

  const prov = result?.provenance

  return (
    <div style={{ display: 'flex', flexDirection: 'column', minHeight: '100%' }}>
      {!atHome && (
        <header
          style={{
            borderBottom: '1px solid var(--om-border-strong)',
            background: 'var(--om-head-bg)',
            padding: '4px 12px',
          }}
        >
          <Group justify="space-between" wrap="nowrap">
            <Anchor
              href="#/"
              underline="never"
              aria-label="OriginMarker home"
              onClick={(e) => {
                e.preventDefault()
                clearAll()
              }}
              style={{ display: 'inline-flex' }}
            >
              <Logo />
            </Anchor>
            <Group gap={4}>
              {(resolved || result || error) && (
                <Button variant="default" size="xs" onClick={clearAll}>
                  Clear
                </Button>
              )}
              <Button variant="subtle" size="xs" component="a" href="#/docs">
                Documentation
              </Button>
              {!health && (
                <Badge size="xs" variant="light" color="red">API unreachable</Badge>
              )}
            </Group>
          </Group>
        </header>
      )}

      <main
        style={{
          flex: 1,
          padding: atHome ? 0 : 12,
          maxWidth: atDocs ? 1180 : atTerms ? 900 : 1500,
          width: '100%',
          margin: '0 auto',
        }}
      >
        {atDocs ? (
          <DocsPage health={health} />
        ) : atTerms ? (
          <TermsPage />
        ) : (
          <>
            {!health && (
              <Alert color="orange" mb="sm" title="Backend not reachable" role="alert">
                <Text size="xs">/api/health did not respond. Queries will fail until the API is running.</Text>
              </Alert>
            )}

            {atHome ? (
              <div className="om-hero">
                <Logo size="hero" />
                <Text size="sm" c="dimmed" mt={2} mb={22} ta="center" style={{ maxWidth: '52ch' }}>
                  Candidate flanking-SNP marker panels for determining which parental allele
                  an embryo inherited.
                </Text>
                <div style={{ width: '100%', maxWidth: 820 }}>
                  <SearchPanel
                    health={health}
                    busy={phase === 'resolving' || phase === 'building'}
                    onResolve={doResolve}
                    hero
                  />
                </div>
                <Anchor href="#/docs" size="xs" mt={18} c="dimmed">
                  How it works, data sources and references
                </Anchor>
              </div>
            ) : (
              <SearchPanel
                health={health}
                busy={phase === 'resolving' || phase === 'building'}
                onResolve={doResolve}
              />
            )}

        {error && (
          <Alert color="red" mb="sm" title="Error" role="alert" withCloseButton onClose={() => setError(null)}>
            <Text size="xs">{error}</Text>
          </Alert>
        )}

        {phase === 'resolving' && (
          <Paper mb="sm" p={8}>
            <Group gap={8}>
              <Loader size="xs" />
              <Text size="xs">Resolving the variant against ClinVar and Ensembl…</Text>
            </Group>
            <Skeleton height={8} mt={8} radius={2} />
            <Skeleton height={8} mt={4} width="70%" radius={2} />
          </Paper>
        )}

        {resolved && phase !== 'resolving' && (
          // Wrapped, not passed by reference: the click event would arrive as buildPanel's
          // query override.
          <VariantCard
            data={resolved}
            nl={nl}
            onBuild={() => buildPanel()}
            busy={phase === 'building'}
          />
        )}

        {phase === 'building' && (
          <Paper mb="sm" p={8}>
            <Group justify="space-between" mb={4}>
              <Group gap={8}>
                <Loader size="xs" />
                <Text size="xs">{progress.stage || 'starting…'}</Text>
              </Group>
              <Text size="xs" c="dimmed" className="om-mono">
                {Math.round(progress.fraction * 100)}%
              </Text>
            </Group>
            <Progress
              value={progress.fraction * 100}
              size="sm"
              radius={2}
              striped
              animated
              aria-label="Panel build progress"
            />
          </Paper>
        )}

        {phase === 'done' && result && prov && (
          <>
            <Paper mb="sm">
              <Group justify="space-between" className="om-section-title" wrap="nowrap">
                <span>Panel</span>
                <Text size="xs" c="dimmed">
                  {int(prov.candidate_n)} candidates · {result.recommended.length} shortlisted ·
                  window ±{int(prov.window_bp)} bp · MAF ≥ {prov.common_maf}
                  {typeof prov.elapsed_s === 'number' ? ` · built in ${prov.elapsed_s.toFixed(1)} s` : ''}
                </Text>
              </Group>
              <Group p={8} gap={8} wrap="wrap" align="flex-end">
                <Select
                  size="xs"
                  w={230}
                  label="Ancestry for the 2pq prior"
                  description="re-orders this page; does not re-select the shortlist"
                  placeholder="Global (best population)"
                  clearable
                  value={ancestry}
                  onChange={setAncestry}
                  data={[...ANCESTRIES]}
                  aria-label="Ancestry for the expected-heterozygosity prior"
                />
                <Group gap={6} ml="auto">
                  {jobId &&
                    (['csv', 'xlsx', 'json', 'pdf'] as const).map((ext) => (
                      <Button
                        key={ext}
                        variant="default"
                        size="xs"
                        component="a"
                        href={api.exportUrl(jobId, ext)}
                        download
                      >
                        {ext.toUpperCase()}
                      </Button>
                    ))}
                </Group>
              </Group>
              {ancestry && ancestry !== prov.ancestry_rank && (
                <Alert color="yellow" variant="light" mx={8} mb={8} title={`This page is ordered on ${ancestry}; the shortlist and the exports are not`}>
                  <Text size="xs">
                    Your browser re-orders the rows on gnomAD {ancestry} frequencies, and that is
                    all it can do. Which markers were shortlisted was decided when the panel was
                    built, on {prov.ancestry_rank ?? 'the global prior'}, and the exports carry
                    that selection in that order. Rebuild to select and rank on {ancestry}{' '}
                    throughout.
                  </Text>
                  {query && (
                    <Button
                      size="xs"
                      variant="default"
                      mt={6}
                      onClick={() => buildPanel({ ...query, ancestry })}
                    >
                      Rebuild on {ancestry}
                    </Button>
                  )}
                </Alert>
              )}
            </Paper>

            <Coverage result={result} />
            <LocusTrack result={result} ancestry={ancestry} />
            <PanelTable result={result} ancestry={ancestry} />

            <Alert color="blue" variant="light" mb="sm" title="These are candidate markers, not a usable panel yet">
              <List size="xs" spacing={2}>
                <List.Item>
                  Genotype the carrier and keep only the markers where that carrier is actually
                  heterozygous. The 2pq values here are population priors, not this carrier's genotype.
                </List.Item>
                <List.Item>
                  Phase the retained markers against the pathogenic allele using family samples.
                  <Text span fw={600}> This app cannot determine phase.</Text>
                </List.Item>
                <List.Item>
                  Parental origin follows from phasing, never from population LD with the pathogenic
                  variant.
                </List.Item>
              </List>
              <Anchor size="xs" href="#/docs/layerb" mt={4} style={{ display: 'inline-block' }}>
                Read the full protocol
              </Anchor>
            </Alert>

            <ResultProvenance result={result} />
          </>
        )}

          </>
        )}
      </main>

      <footer
        style={{
          borderTop: '1px solid var(--om-border)',
          background: 'var(--om-head-bg)',
          padding: '8px 12px',
          marginTop: 'auto',
        }}
      >
        {/* R8: the disclaimer appears on every page, including the landing view. */}
        <Group justify="space-between" wrap="nowrap" gap={12}>
          <Text size="xs" style={{ color: '#4d545c' }}>
            {prov?.disclaimer ?? health?.disclaimer ?? FALLBACK_DISCLAIMER}
          </Text>
          <Group gap={10} wrap="nowrap">
            {health?.release && (
              <Text size="xs" c="dimmed" style={{ whiteSpace: 'nowrap' }} title={health.release_gloss}>
                {health.release}
              </Text>
            )}
            <Anchor href="#/docs" size="xs" style={{ whiteSpace: 'nowrap' }}>
              Documentation
            </Anchor>
            <Text size="xs" c="dimmed" aria-hidden>|</Text>
            <Anchor href="#/terms" size="xs" style={{ whiteSpace: 'nowrap' }}>
              Terms
            </Anchor>
          </Group>
        </Group>
      </footer>
    </div>
  )
}

/**
 * Provenance for one result, rendered with that result.
 *
 * Takes `result` and nothing else, and that is the invariant: every fact must come off the
 * panel, and one the panel does not carry renders as `unknown`. Reading a live source such
 * as /api/health here would let a frozen panel report versions it was not built against.
 */
function ResultProvenance({ result }: { result: PanelResult }) {
  const p = result.provenance
  const approx = result.candidates.some((m) => m.map_approx)
  return (
    <Paper mb="sm">
      <div className="om-section-title">Provenance</div>
      <div style={{ padding: 8 }}>
        <Text size="xs" c="dimmed" className="om-mono">
          Data as of {utc(p.queried_utc)} · {orUnknown(p.sources?.gnomad)} · Ensembl{' '}
          {orUnknown(p.ensembl_release)} · {p.build}
          {p.source_responses_from_cache ? ` · ${p.source_responses_from_cache} of ${
            (p.source_responses_from_cache ?? 0) + (p.source_responses_from_network ?? 0)
          } source responses served from cache` : ''}
          {p.requested_build && p.requested_build !== p.build
            ? ` (input given as ${p.requested_build}, converted)`
            : ''}
        </Text>
        <Text size="xs" c="dimmed" className="om-mono" mt={2}>
          Genetic map: {orUnknown(p.sources?.genetic_map)}
        </Text>
        <Text size="xs" c="dimmed" mt={2}>
          Window ±{int(p.window_bp)} bp · MAF floor {p.common_maf} · {int(p.candidate_n)}{' '}
          candidates · built in {p.elapsed_s}s
        </Text>
        {approx && (
          <Text size="xs" c="orange.8" mt={4}>
            Some markers fall outside the bundled genetic map: their cM values are a uniform
            1 cM/Mb approximation, marked with a tilde in the table.
          </Text>
        )}
        <Anchor href="#/docs" size="xs" mt={4} style={{ display: 'inline-block' }}>
          Data sources and references →
        </Anchor>
      </div>
    </Paper>
  )
}

function Coverage({ result }: { result: PanelResult }) {
  const { coverage: c, recommended } = result
  // Not every flag is a coverage flag: hotspot flags fire with both sides fully covered,
  // so the title and the advice below key off the side test, not off flags.length.
  const thinSide = c.lower_core_near < 2 || c.higher_core_near < 2
  const byTier = (cen: boolean) => {
    const side = recommended.filter((m) => m.dist > 0 === cen)
    const band = (lo: number, hi: number) =>
      side.filter((m) => Math.abs(m.dist) >= lo && Math.abs(m.dist) < hi).length
    return { A: band(0, 2_000), B: band(2_000, 30_000), C: band(30_000, Infinity) }
  }

  return (
    <Paper mb="sm">
      <div className="om-section-title">Flanking coverage</div>
      <div style={{ padding: 8 }}>
        {c.flags.length > 0 && (
          <Alert
            color="red"
            mb={8}
            title={thinSide ? 'Flanking coverage is incomplete' : 'Warnings about this marker set'}
            role="alert"
          >
            <List size="xs" spacing={2}>
              {c.flags.map((f) => (
                <List.Item key={f}>{f}</List.Item>
              ))}
            </List>
            {thinSide && (
              <Text size="xs" mt={4}>
                A one-sided marker set cannot detect a recombination between the variant and the
                markers. Widen the window or lower the MAF floor.
              </Text>
            )}
          </Alert>
        )}
        {/* Labelled by coordinate, not by arm: the engine's tel_/cen_ field names do not
            track which flank runs toward the telomere. That depends on the centromere
            position, which nothing here looks up. */}
        <Group gap={24} wrap="wrap">
          <SideStat
            label="Lower coordinate"
            total={c.lower_count}
            near={c.lower_core_near}
            tiers={byTier(false)}
            color="var(--om-lower)"
            ok={c.lower_core_near >= 2}
          />
          <SideStat
            label="Higher coordinate"
            total={c.higher_count}
            near={c.higher_core_near}
            tiers={byTier(true)}
            color="var(--om-higher)"
            ok={c.higher_core_near >= 2}
          />
        </Group>
      </div>
    </Paper>
  )
}

const SideStat = ({
  label, total, near, tiers, color, ok,
}: {
  label: string
  total: number
  near: number
  tiers: { A: number; B: number; C: number }
  color: string
  ok: boolean
}) => (
  <div style={{ borderLeft: `3px solid ${color}`, paddingLeft: 8, minWidth: 230 }}>
    <Group gap={6}>
      <Text size="xs" fw={600}>{label}</Text>
      <Badge size="xs" color={ok ? 'green' : 'red'} variant="light">
        {ok ? 'covered' : 'under-covered'}
      </Badge>
    </Group>
    <Text size="xs" c="dimmed" className="om-mono">
      {total} shortlisted · {near} within 30 kb
    </Text>
    <Text size="xs" c="dimmed" className="om-mono">
      A(&lt;2kb) {tiers.A} · B(2–30kb) {tiers.B} · C(30kb+) {tiers.C}
    </Text>
  </div>
)

// R8: a byte-for-byte copy of pb.DISCLAIMER, never a paraphrase. Used only when the API is
// unreachable and cannot supply the canonical string itself.
const FALLBACK_DISCLAIMER =
  'Research use only. Candidate markers require validation and per-family phasing in a qualified genetics laboratory. Not a clinical diagnostic.'


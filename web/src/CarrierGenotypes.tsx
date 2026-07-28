import { useRef, useState } from 'react'
import { Alert, Anchor, Badge, Button, Group, Paper, Progress, Select, Text } from '@mantine/core'
import {
  checkBuild, emptySet, informativeCoverage, ingestLine,
  type GenotypeSet, type ParseOpts,
} from './genotypes'
import type { Marker } from './api'
import { int } from './fmt'

/**
 * Load the carrier's own genotypes and keep the markers they are actually heterozygous at.
 *
 * The panel ranks on expected heterozygosity, which is a population average and says
 * nothing about the person being tested; the lab protocol's first two steps are to genotype
 * that carrier and drop the markers where they turn out to be homozygous. This does that
 * step, against their real file.
 *
 * The file is read here, in this browser, and there is no endpoint to send it to. That is
 * not a preference: a carrier's WGS is the most identifying artefact in this workflow and
 * the terms promise nothing about a family is submitted or retained.
 */
export function CarrierGenotypes({
  markers, onLoad, set,
}: {
  markers: Marker[]
  onLoad: (s: GenotypeSet | null) => void
  set: GenotypeSet | null
}) {
  const [busy, setBusy] = useState(false)
  const [read, setRead] = useState(0)
  const [error, setError] = useState<string | null>(null)
  const [sample, setSample] = useState<string | null>(null)
  const fileRef = useRef<HTMLInputElement>(null)
  const lastFile = useRef<File | null>(null)

  // The panel's own span. A whole-genome file is discarded down to this as it streams, so
  // the memory cost is the window and not the file.
  const span = markers.length
    ? {
        chrom: markers[0].chrom,
        lo: Math.min(...markers.map((m) => m.pos)),
        hi: Math.max(...markers.map((m) => m.pos)),
      }
    : null

  const parse = async (file: File, wantSample?: string) => {
    if (!span) return
    setBusy(true); setError(null); setRead(0)
    lastFile.current = file
    const opts: ParseOpts = { ...span, sample: wantSample }
    const next = emptySet()
    try {
      // Streamed, not read whole: a WGS VCF is tens of gigabytes and would not fit in
      // memory, and every line outside the window is discarded as it arrives.
      const reader = file.stream().getReader()
      const decoder = new TextDecoder()
      let tail = ''
      let bytes = 0
      for (;;) {
        const { done, value } = await reader.read()
        if (done) break
        bytes += value.byteLength
        const text = tail + decoder.decode(value, { stream: true })
        const lines = text.split('\n')
        tail = lines.pop() ?? ''
        for (const line of lines) ingestLine(line.replace(/\r$/, ''), next, opts)
        setRead(Math.round((bytes / file.size) * 100))
      }
      if (tail) ingestLine(tail, next, opts)
      if (next.linesRead === 0) {
        setError('No genotype lines were recognised. Expected a VCF, or a tab-separated '
                 + 'array export with rsID, chromosome, position and genotype columns.')
        onLoad(null)
      } else {
        setSample(next.sample)
        onLoad(next)
      }
    } catch (e) {
      setError(e instanceof Error ? `Could not read the file: ${e.message}` : 'Could not read the file.')
      onLoad(null)
    } finally {
      setBusy(false)
    }
  }

  const cov = set ? informativeCoverage(markers, set) : null
  const build = set ? checkBuild(markers, set) : null
  const thinSide = cov ? cov.lower < 2 || cov.higher < 2 : false

  return (
    <Paper mb="sm">
      <Group justify="space-between" className="om-section-title" wrap="nowrap">
        <span>Carrier genotypes</span>
        <Text size="xs" c="dimmed">read in this browser, never uploaded</Text>
      </Group>

      <div style={{ padding: 8 }}>
        <Text size="xs" c="dimmed" mb={8}>
          The 2pq figures on this page are population priors. Load this carrier's own VCF or
          SNP-array export and the panel can say which markers they are actually
          heterozygous at, which is the first thing the lab protocol asks for. The file is
          parsed here and there is no endpoint that receives it.{' '}
          <Anchor href="#/docs/carrier" size="xs">What this does and does not settle</Anchor>
        </Text>

        <Group gap={8} align="center">
          <Button
            size="xs"
            variant="default"
            loading={busy}
            onClick={() => fileRef.current?.click()}
          >
            {set ? 'Load a different file' : 'Load genotypes'}
          </Button>
          <input
            ref={fileRef}
            type="file"
            accept=".vcf,.txt,.tsv,.csv,.gz"
            style={{ display: 'none' }}
            onChange={(e) => {
              const f = e.currentTarget.files?.[0]
              if (f) parse(f)
              e.currentTarget.value = ''      // same file twice must re-trigger
            }}
          />
          {set && (
            <Button size="xs" variant="subtle" color="gray" onClick={() => { onLoad(null); setSample(null) }}>
              Clear
            </Button>
          )}
          {set && set.samples.length > 1 && (
            <Select
              size="xs"
              w={200}
              aria-label="VCF sample column"
              label={undefined}
              allowDeselect={false}
              value={sample}
              data={set.samples}
              onChange={(v) => {
                setSample(v)
                if (lastFile.current && v) parse(lastFile.current, v)
              }}
            />
          )}
          {busy && <Progress value={read} size="sm" w={160} striped animated aria-label="Reading file" />}
        </Group>

        {error && (
          <Alert color="orange" p={6} mt={8} role="alert">
            <Text size="xs">{error}</Text>
          </Alert>
        )}

        {/* A .gz cannot be read as text, and saying so beats a file that parses to nothing. */}
        {set === null && !busy && !error && (
          <Text size="xs" c="dimmed" mt={6}>
            Plain text only: gunzip a .vcf.gz first. A whole-genome VCF works but is read at
            disk speed, so slicing the window first is faster:{' '}
            <code className="om-mono">
              tabix -h in.vcf.gz {span ? `${span.chrom}:${int(span.lo)}-${int(span.hi)}` : 'chr:lo-hi'} &gt; window.vcf
            </code>
          </Text>
        )}

        {build?.mismatch && (
          <Alert color="red" mt={8} role="alert" title="This file looks like a different genome assembly">
            <Text size="xs">
              {build.agreed} of {build.checked} markers matched by rsID sit at the coordinate
              this panel expects; the rest are offset (for example by{' '}
              {build.offsets.slice(0, 3).map((o) => (o > 0 ? `+${int(o)}` : int(o))).join(', ')} bp).
              This panel is GRCh38. Markers still match by rsID, but any site your file leaves
              unnamed is being matched on a coordinate that means something else here. Lift
              the file over to GRCh38, or rely only on the rsID-matched rows.
            </Text>
          </Alert>
        )}

        {set && cov && (
          <div style={{ marginTop: 10 }}>
            <Group gap={8} mb={6} wrap="wrap">
              <Badge size="sm" variant="light" color={cov.het ? 'green' : 'red'}>
                {cov.het} heterozygous
              </Badge>
              <Text size="xs" c="dimmed" className="om-mono">
                {cov.lower} lower coord · {cov.higher} higher coord · {cov.hom} homozygous ·{' '}
                {cov.nocall} no-call · {cov.absent} not in file
              </Text>
              <Text size="xs" c="dimmed">
                {set.format === 'vcf' ? `VCF${set.sample ? `, sample ${set.sample}` : ''}` : 'array export'}
                {' · '}{int(set.sitesKept)} sites in this window
              </Text>
            </Group>

            {thinSide && (
              <Alert color="orange" p={8} role="alert" title="Not enough informative markers on one side">
                <Text size="xs">
                  A side with fewer than two markers this carrier is heterozygous at cannot
                  survive a single dropout or a crossover between the marker and the variant.
                  Widen the window, lower the MAF floor, or genotype more densely there.
                </Text>
              </Alert>
            )}

            {/* Absent is not homozygous-reference, and a reader counting it as one would
                write off markers that are simply not in a variants-only file. */}
            {cov.absent > 0 && (
              <Text size="xs" c="dimmed" mt={6}>
                {cov.absent} marker{cov.absent === 1 ? ' is' : 's are'} not in the file. In a
                variants-only VCF that usually means homozygous reference, and so
                uninformative, but a region nobody called looks identical. Only an all-sites
                or gVCF file tells those apart, so they are left uncounted rather than
                assumed.
              </Text>
            )}

            <Alert color="blue" variant="light" p={8} mt={8} title="This still does not give you phase">
              <Text size="xs">
                Knowing the carrier is heterozygous says which markers are worth typing. It
                does not say which of their two chromosomes carries the pathogenic allele,
                and short-read sequencing cannot bridge tens of kilobases to tell you. That
                needs an informative relative, or reads long enough to span the distance.{' '}
                <Anchor href="#/docs/carrier" size="xs">Choosing a platform</Anchor>
              </Text>
            </Alert>
          </div>
        )}
      </div>
    </Paper>
  )
}

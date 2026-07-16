import { useEffect, useMemo, useState } from 'react'
import {
  Alert, Anchor, Badge, Group, List, Pagination, Paper, Select, Table, Text, TextInput,
  Tooltip,
} from '@mantine/core'
import { isUpper, links, posMismatch, TIER_LABEL, type Marker, type PanelResult } from './api'
import { int, num, sig2, signedBp } from './fmt'
import { ancestryHet, applyPreset, type Preset, PRESETS } from './rank'

const PER_PAGE = 50

/** "MISMATCH:17401130" -> 17401130, the position Ensembl claims (panelbuilder's format).
 *  Anything unparseable yields null and is rendered raw rather than as NaN. */
const ensemblPos = (check: string | null) => {
  const n = Number((check ?? '').split(':')[1])
  return Number.isFinite(n) ? n : null
}

/**
 * Mirror of panelbuilder.select_panel's `assessed`.
 *
 * On an approximated map stretch the uniform 1 cM/Mb fallback can never reach the hotspot
 * threshold, so `hotspot_between === false` there means "not assessed", not "no hotspot".
 */
const hotspotAssessed = (m: Marker) => m.hotspot_between != null && m.map_approx === false

type SortKey = 'pos' | 'dist' | 'maf' | 'het' | 'anc' | 'cm' | 'theta'

interface Props {
  result: PanelResult
  ancestry: string | null
}

export function PanelTable({ result, ancestry }: Props) {
  const { candidates, recommended } = result
  const [scope, setScope] = useState<'recommended' | 'all'>('recommended')
  const [preset, setPreset] = useState<Preset>('ranked')
  const [sort, setSort] = useState<{ key: SortKey; desc: boolean } | null>(null)
  const [q, setQ] = useState('')
  const [page, setPage] = useState(1)

  const recSet = useMemo(() => new Set(recommended.map((m) => m.variant_id)), [recommended])
  const rows = scope === 'recommended' ? recommended : candidates

  // Derived from the whole result, never from `rows` or the current page: this alarm must
  // not be silenceable by filtering or switching scope. The engine cross-checks positions
  // after selecting the panel, so a disputed marker is still shortlisted.
  const disputed = useMemo(() => candidates.filter(posMismatch), [candidates])

  useEffect(() => setPage(1), [scope, preset, q, ancestry, sort])

  const filtered = useMemo(() => {
    const t = q.trim().toLowerCase()
    if (!t) return rows
    return rows.filter(
      (m) =>
        m.rsid?.toLowerCase().includes(t) ||
        String(m.pos).includes(t) ||
        m.tier.toLowerCase().includes(t) ||
        m.side.toLowerCase().includes(t),
    )
  }, [rows, q])

  const sorted = useMemo(() => {
    if (!sort) return applyPreset(filtered, preset, ancestry)
    const val = (m: Marker): number => {
      switch (sort.key) {
        case 'pos': return m.pos
        case 'dist': return m.dist
        case 'maf': return m.maf
        case 'het': return m.het
        case 'anc': return ancestryHet(m, ancestry) ?? -1
        case 'cm': return m.cm ?? -1
        case 'theta': return m.recomb_fraction ?? -1
      }
    }
    return [...filtered].sort((a, b) => (sort.desc ? val(b) - val(a) : val(a) - val(b)))
  }, [filtered, sort, preset, ancestry])

  const pages = Math.max(1, Math.ceil(sorted.length / PER_PAGE))
  const pageRows = sorted.slice((page - 1) * PER_PAGE, page * PER_PAGE)

  const th = (key: SortKey, label: string, help?: string) => {
    const active = sort?.key === key
    const el = (
      <Table.Th
        className="om-sortable"
        onClick={() => setSort(active && !sort.desc ? { key, desc: true } : { key, desc: false })}
        role="columnheader"
        aria-sort={active ? (sort.desc ? 'descending' : 'ascending') : 'none'}
        tabIndex={0}
        onKeyDown={(e) =>
          e.key === 'Enter' && setSort(active && !sort.desc ? { key, desc: true } : { key, desc: false })
        }
      >
        {label}
        {active ? (sort.desc ? ' ▾' : ' ▴') : ''}
      </Table.Th>
    )
    return help ? <Tooltip key={key} label={help} withArrow multiline w={260}>{el}</Tooltip> : el
  }

  return (
    <Paper mb="sm">
      <Group justify="space-between" className="om-section-title" wrap="nowrap">
        <span>Candidate markers</span>
        <Text size="xs" c="dimmed">
          scored on a population prior: expected heterozygosity, then proximity
        </Text>
      </Group>

      {disputed.length > 0 && (
        <Alert
          color="red"
          mx={8}
          mt={8}
          role="alert"
          title={
            disputed.length === 1
              ? 'A marker in this panel is at a disputed position'
              : `${disputed.length} markers in this panel are at disputed positions`
          }
        >
          <Text size="xs" mb={6}>
            gnomAD and Ensembl report different GRCh38 coordinates for{' '}
            {disputed.length === 1 ? 'this marker' : 'these markers'}. Only one source can be
            right, and every figure on the row was computed from the gnomAD position:
            distance to the variant, cM, θ, and which flank it counts toward. The panel was
            selected before the positions were re-checked, so{' '}
            {disputed.length === 1 ? 'it is' : 'they are'} still shortlisted. Resolve the
            disagreement at dbSNP before ordering an assay, or drop{' '}
            {disputed.length === 1 ? 'the marker' : 'them'}.
          </Text>
          <List size="xs" spacing={2}>
            {disputed.map((m) => {
              const e = ensemblPos(m.ensembl_pos_check)
              return (
                <List.Item key={m.variant_id}>
                  <span className="om-mono">
                    {m.rsid}: gnomAD {int(m.pos)} · Ensembl{' '}
                    {e == null ? m.ensembl_pos_check : `${int(e)} (${int(Math.abs(e - m.pos))} bp apart)`}
                  </span>
                  {recSet.has(m.variant_id) && (
                    <Text span size="xs" fw={600}> · in the shortlist</Text>
                  )}
                </List.Item>
              )
            })}
          </List>
        </Alert>
      )}

      <Group gap={8} p={8} wrap="wrap">
        <Select
          size="xs"
          w={190}
          label={undefined}
          aria-label="Marker set"
          allowDeselect={false}
          value={scope}
          onChange={(v) => setScope(v as 'recommended' | 'all')}
          data={[
            { value: 'recommended', label: `Shortlist (${recommended.length})` },
            { value: 'all', label: `All candidates (${int(candidates.length)})` },
          ]}
        />
        <Select
          size="xs"
          w={190}
          aria-label="Column preset"
          allowDeselect={false}
          value={preset}
          onChange={(v) => { setPreset(v as Preset); setSort(null) }}
          data={PRESETS}
        />
        <TextInput
          size="xs"
          w={190}
          aria-label="Filter markers"
          placeholder="filter rsID / position / tier"
          value={q}
          onChange={(e) => setQ(e.currentTarget.value)}
        />
        {sort && (
          <Anchor size="xs" component="button" type="button" onClick={() => setSort(null)}>
            clear column sort
          </Anchor>
        )}
        <Text size="xs" c="dimmed" ml="auto">
          {int(sorted.length)} rows
        </Text>
      </Group>

      <div style={{ overflowX: 'auto' }}>
        <Table className="om-table" striped="even" highlightOnHover stickyHeader>
          <Table.Thead>
            <Table.Tr>
              <Table.Th>rsID</Table.Th>
              {th('pos', 'GRCh38 coord', 'Genomic position on GRCh38. All coordinates on this page are GRCh38.')}
              {th('dist', 'dist (bp)', 'Signed distance from the pathogenic variant. Negative = lower GRCh38 coordinate.')}
              <Table.Th>side</Table.Th>
              <Table.Th>tier</Table.Th>
              <Table.Th>alleles</Table.Th>
              {th('maf', 'MAF', 'gnomAD v4 global minor allele frequency.')}
              {th('het', '2pq global', 'Expected heterozygosity under Hardy-Weinberg: a POPULATION PRIOR, not this carrier\'s genotype. The carrier must be genotyped.')}
              {th('anc', ancestry ? `2pq ${ancestry}` : '2pq best pop', ancestry
                ? `Expected heterozygosity in ${ancestry}: a population prior, not this carrier's genotype.`
                : 'Highest expected heterozygosity across gnomAD populations. A population prior. Select an ancestry to see a matched value.')}
              {th('cm', 'cM', 'Genetic distance to the variant, interpolated from the bundled map. Shown to two significant figures: the map is population- and sex-averaged, and does not resolve finer than that.')}
              {th('theta', 'θ', 'Haldane recombination fraction over that genetic distance, to two significant figures. It carries the map\'s averaging: it is the expected rate across many meioses, not what happened in this one.')}
              <Table.Th>hotspot</Table.Th>
              <Tooltip
                label="Position re-checked against Ensembl, a source independent of the gnomAD record it came from. Only the nearest shortlisted markers are checked, so most rows are blank."
                withArrow
                multiline
                w={280}
              >
                <Table.Th>x-check</Table.Th>
              </Tooltip>
              <Table.Th>links</Table.Th>
            </Table.Tr>
          </Table.Thead>
          <Table.Tbody>
            {pageRows.map((m) => {
              const anc = ancestryHet(m, ancestry)
              return (
                <Table.Tr
                  key={m.variant_id}
                  // Inline, not a class: it has to beat the zebra striping rather than
                  // alternate with it.
                  style={posMismatch(m) ? { background: 'rgba(224, 49, 49, 0.10)' } : undefined}
                >
                  <Table.Td className="om-mono">
                    <Group gap={4} wrap="nowrap">
                      {recSet.has(m.variant_id) && scope === 'all' && (
                        <Tooltip label="shortlisted" withArrow>
                          <span style={{ color: 'var(--om-blue)' }} aria-label="shortlisted">●</span>
                        </Tooltip>
                      )}
                      {m.rsid}
                    </Group>
                  </Table.Td>
                  <Table.Td className="om-mono">{int(m.pos)}</Table.Td>
                  <Table.Td className="om-mono" style={{ textAlign: 'right' }}>{signedBp(m.dist)}</Table.Td>
                  <Table.Td>
                    <Badge
                      size="xs"
                      variant="light"
                      color={isUpper(m) ? 'orange' : 'blue'}
                    >
                      {isUpper(m) ? 'higher' : 'lower'}
                    </Badge>
                  </Table.Td>
                  <Table.Td>
                    <Text size="xs" c="dimmed">{TIER_LABEL[m.tier] ?? m.tier}</Text>
                  </Table.Td>
                  <Table.Td className="om-mono">{m.ref}&gt;{m.alt}</Table.Td>
                  <Table.Td className="om-mono">{num(m.maf)}</Table.Td>
                  <Table.Td className="om-mono">{num(m.het)}</Table.Td>
                  <Table.Td className="om-mono">
                    {anc == null ? (
                      <Tooltip label={ancestry ? `no gnomAD ${ancestry} frequency for this marker` : ''} withArrow disabled={!ancestry}>
                        <span>{ancestry ? '-' : num(m.het_max_pop)}</span>
                      </Tooltip>
                    ) : (
                      num(anc)
                    )}
                  </Table.Td>
                  <Table.Td className="om-mono">
                    {sig2(m.cm)}
                    {m.map_approx && (
                      <Tooltip label="cM approximated at 1 cM/Mb: no map data here" withArrow>
                        <span style={{ color: 'var(--om-text-dim)' }}>~</span>
                      </Tooltip>
                    )}
                  </Table.Td>
                  <Table.Td className="om-mono">{sig2(m.recomb_fraction)}</Table.Td>
                  <Table.Td>
                    {m.hotspot_between ? (
                      <Tooltip label="A recombination hotspot lies between this marker and the variant" withArrow>
                        <Badge size="xs" color="red" variant="light">hotspot</Badge>
                      </Tooltip>
                    ) : hotspotAssessed(m) ? (
                      <Tooltip label="The map records no hotspot between this marker and the variant" withArrow>
                        <Text size="xs" c="dimmed" style={{ cursor: 'help' }}>none</Text>
                      </Tooltip>
                    ) : (
                      <Tooltip
                        label="Not assessed: without real map data here, a hotspot cannot be ruled in or out"
                        withArrow
                      >
                        <Text size="xs" c="dimmed" style={{ cursor: 'help' }}>not assessed</Text>
                      </Tooltip>
                    )}
                  </Table.Td>
                  <Table.Td>
                    {m.ensembl_pos_check == null ? (
                      // null is silence, not a pass: only the nearest shortlisted markers
                      // are sent to Ensembl.
                      <Tooltip label="Not checked: this marker's position rests on the gnomAD record alone" withArrow>
                        <Text size="xs" c="dimmed" style={{ cursor: 'help' }}>not checked</Text>
                      </Tooltip>
                    ) : m.ensembl_pos_check === 'ok' ? (
                      <Tooltip label="Ensembl reports the same GRCh38 position" withArrow>
                        <Text size="xs" c="green">match</Text>
                      </Tooltip>
                    ) : (
                      <Tooltip label={`Ensembl disagrees: ${m.ensembl_pos_check}`} withArrow>
                        <Badge size="xs" color="red">mismatch</Badge>
                      </Tooltip>
                    )}
                  </Table.Td>
                  <Table.Td>
                    <Group gap={5} wrap="nowrap">
                      <Ext href={links.dbsnp(m.rsid)} label={`${m.rsid} at dbSNP`}>dbSNP</Ext>
                      <Ext href={links.gnomad(m)} label={`${m.rsid} at gnomAD`}>gnomAD</Ext>
                      <Ext href={links.ensembl(m.rsid)} label={`${m.rsid} at Ensembl`}>Ens</Ext>
                      <Ext href={links.ucsc(m.chrom, m.pos)} label={`${m.rsid} at UCSC`}>UCSC</Ext>
                    </Group>
                  </Table.Td>
                </Table.Tr>
              )
            })}
          </Table.Tbody>
        </Table>
      </div>

      {!pageRows.length && (
        <Text size="xs" c="dimmed" ta="center" py="md">
          No markers match this filter.
        </Text>
      )}

      <Group justify="space-between" p={8}>
        {pages > 1 && (
          <Pagination size="xs" total={pages} value={page} onChange={setPage} withEdges />
        )}
      </Group>
    </Paper>
  )
}

const Ext = ({ href, label, children }: { href: string; label: string; children: React.ReactNode }) => (
  <Anchor href={href} target="_blank" rel="noreferrer" size="xs" aria-label={label}>
    {children}
  </Anchor>
)

import { Fragment, useEffect, useMemo, useState } from 'react'
import {
  Alert, Anchor, Badge, Group, List, Pagination, Paper, Select, Table, Text, TextInput,
  Tooltip,
} from '@mantine/core'
import {
  flankingRule, hasPair, isUpper, links, pcrDanger, pcrUnchecked, posMismatch, primerBuild,
  starred, TIER_LABEL, type Marker, type PanelResult, type Primer, type PrimerBuild,
  type PrimerResult, type PrimerWarning,
} from './api'
import { int, num, sig2, signedBp } from './fmt'
import { PrimerChip } from './PrimerOptions'
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

/** What the star claims, for a reader who cannot see it. Short because a screen reader
 *  reads it once per starred row; the engine's full wording is the legend, in page text. */
const STAR_CLAIM = 'meets ESHRE flanking criteria'

/** The engine's legend opens with a literal ★ naming the glyph. The words are rendered
 *  verbatim; only that one character is dropped, because the real glyph is drawn in its
 *  place. A legend without the prefix loses nothing. */
const legendWords = (legend: string) => legend.replace(/^\s*★\s*/, '')

/** Shape, not colour, carries the star: it survives a greyscale print and a colour-blind
 *  reader. Yellow is decoration only, and the dark rim is what makes it legible on white. */
const Star = () => (
  <svg width={11} height={11} viewBox="0 0 20 20" aria-hidden="true" focusable="false">
    <path
      d="M10 1.6l2.6 5.2 5.8.9-4.2 4.1 1 5.7-5.2-2.7-5.2 2.7 1-5.7-4.2-4.1 5.8-.9z"
      fill="#f5b301"
      stroke="#6b5200"
      strokeWidth={1.4}
      strokeLinejoin="round"
    />
  </svg>
)

/** One oligo, as ordered. Monospace and selectable: this string gets copied onto an order
 *  form, so it is rendered whole and never truncated. */
const Oligo = ({ side, p }: { side: string; p: Primer }) => (
  <Group gap={8} wrap="nowrap" align="baseline">
    <Text size="xs" fw={700} style={{ width: 12, flex: 'none' }}>{side}</Text>
    <Text
      size="xs"
      className="om-mono"
      style={{ userSelect: 'text', wordBreak: 'break-all' }}
    >
      5'-{p.seq}-3'
    </Text>
    {/* Its own length, not seq.length: the engine states it, and a second way of counting
        is a second answer. GRCh38 on the coordinate, because every coordinate here is. */}
    <Text size="xs" c="dimmed" className="om-mono" style={{ whiteSpace: 'nowrap', flex: 'none' }}>
      {p.length} nt · Tm {p.tm.toFixed(1)} °C · GC {p.gc.toFixed(1)}% · GRCh38 {int(p.pos)}
    </Text>
  </Group>
)

/**
 * One line of the engine's own words, at the length a table can carry.
 *
 * `short` and a link, never `long`: `long` is a paragraph, it is the same paragraph on every
 * row, and the reader who wants it has the docs and the PDF. Nothing here rewrites either
 * one, and the link is the deferral the short form is written to make.
 */
const Note = ({ w, tone }: { w: PrimerWarning; tone: 'warn' | 'danger' | 'ok' }) => (
  <Group gap={6} mb={4} wrap="nowrap" align="flex-start">
    {tone !== 'ok' && (
      <Badge
        size="xs"
        variant="filled"
        color={tone === 'danger' ? 'red' : 'yellow.4'}
        c={tone === 'danger' ? undefined : '#3d2f00'}
        style={{ flex: 'none' }}
      >
        {tone === 'danger' ? 'DANGER' : 'WARNING'}
      </Badge>
    )}
    <Text size="xs" c={tone === 'ok' ? 'green.9' : undefined}>
      {noteWords(w)}{' '}
      <Anchor href={w.docs || '#/docs/primers'} size="xs">why</Anchor>
    </Text>
  </Group>
)

/**
 * A note's words, and never nothing.
 *
 * `short` is what belongs in a table. The rest is not decoration: a badge reading WARNING
 * beside an empty line is a warning that says nothing, which reads as a pair with nothing
 * wrong with it. That is the one thing this module refuses, and it is exactly what an older
 * server's shape rendered here.
 */
const noteWords = (w: PrimerWarning) =>
  w.short || w.long
  || 'This pair carries a warning this page could not read. Treat it as unchecked.'

/** The engine promises never to hand back a pair with nothing said about it. This is that
 *  contract broken, so it says so rather than rendering a clean-looking pair. */
const NO_WARNINGS: PrimerWarning = {
  code: 'no_warnings',
  short: 'This pair arrived with no warnings, which the design does not do. Treat it as unchecked.',
  long: '',
  docs: '#/docs/primers',
}

/**
 * One marker's pair, under its row, behind a line you open.
 *
 * Collapsed by default: the pair is four lines of detail per row, and a table that shows all
 * of it for every marker is one nobody reads. What collapsing must never hide is a finding,
 * so the summary line carries the verdict badge, and a dangerous pair opens itself. The
 * panel-level alert lists every dangerous pair regardless of what is open.
 *
 * Every state renders the pair underneath: a warned primer is information, a hidden one is a
 * decision taken for the reader.
 */
export const PrimerDetail = ({ d, defaultOpen }: {
  d: PrimerResult
  /** Overrides the rule below. The rule is the default; this exists for a caller that
   *  already knows what it wants open. */
  defaultOpen?: boolean
}) => {
  const danger = pcrDanger(d)
  const pair = hasPair(d)
  const warnings = d.warnings ?? []
  // Danger and a failed design are the two states worth the reader's attention before they
  // ask for it. Everything else starts shut.
  const [open, setOpen] = useState(defaultOpen ?? (!!danger || !!d.error))
  return (
    <Table.Tr>
      {/* Wider than the table on purpose: the header is the one place the column count is
          written down, and a copy of it here would drift on the next column. */}
      <Table.Td
        colSpan={99}
        p={0}
        style={{ background: 'var(--om-head-bg)', borderLeft: '3px solid var(--om-blue)' }}
      >
      {/* Pinned and bounded, because the cell is as wide as the table and the table is wider
          than its scroll box: unpinned, a warning runs off the right edge and the reader has
          to scroll sideways to find out what is wrong. `normal` undoes the table's nowrap,
          which is right for a row of figures and silently truncates a sentence. */}
      <div
        style={{
          position: 'sticky',
          left: 0,
          width: 'min(900px, calc(100vw - 48px))',
          whiteSpace: 'normal',
        }}
      >
        <button
          type="button"
          className="om-primer-toggle"
          aria-expanded={open}
          onClick={() => setOpen((o) => !o)}
        >
          <span className="om-primer-toggle-caret">{open ? '▾' : '▸'}</span>
          <span>{open ? 'Hide primer design' : 'Open primer design'}</span>
          {/* On the line, not inside the box: a finding that only exists once the reader
              opens the row is a finding hidden behind a click. */}
          {danger && (
            <Badge size="xs" color="red" variant="filled" style={{ flex: 'none' }}>DANGER</Badge>
          )}
          {!danger && d.error && (
            <Badge size="xs" color="orange" variant="light" style={{ flex: 'none' }}>no pair</Badge>
          )}
          {!danger && pair && !pcrUnchecked(d) && (
            <Badge size="xs" color="green" variant="light" style={{ flex: 'none' }}>
              one product
            </Badge>
          )}
        </button>
        {open && (
          <div style={{ padding: '2px 8px 8px' }}>
            {d.error && (
              <Alert color="orange" role="alert" p={6} mb={6} title="No primer pair for this marker">
                <Text size="xs">{d.error}</Text>
              </Alert>
            )}
            {pair && danger && (
              <Alert
                color="red"
                role="alert"
                p={6}
                mb={6}
                title="Do not order this pair without redesigning it"
              >
                <Text size="xs">
                  {danger.short}{' '}
                  <Anchor href={danger.docs || '#/docs/primers'} size="xs">why</Anchor>
                </Text>
              </Alert>
            )}
            {/* Not alerts: these are the normal state of every pair the design hands back, and
                a red banner on all of them is one nobody reads by the third row. A pass is
                still the engine's sentence, not a summary of it. */}
            {pair && !danger && (warnings.length ? warnings : [NO_WARNINGS]).map((w) => (
              <Note
                key={w.code}
                w={w}
                tone={warnings.length && !pcrUnchecked(d) ? 'ok' : 'warn'}
              />
            ))}
            {d.fwd && <Oligo side="F" p={d.fwd} />}
            {d.rev && <Oligo side="R" p={d.rev} />}
            {(d.product_size != null || d.mask_note) && (
              <Text size="xs" c="dimmed" mt={3}>
                {[
                  d.product_size != null ? `product ${int(d.product_size)} bp` : null,
                  d.mask_note ? noteWords(d.mask_note) : null,
                ].filter(Boolean).join(' · ')}
                {d.mask_note && (
                  <>
                    {' '}
                    <Anchor href={d.mask_note.docs || '#/docs/primers'} size="xs">why</Anchor>
                  </>
                )}
              </Text>
            )}
          </div>
        )}
      </div>
      </Table.Td>
    </Table.Tr>
  )
}

type SortKey = 'pos' | 'dist' | 'maf' | 'het' | 'anc' | 'cm' | 'theta'

interface Props {
  result: PanelResult
  ancestry: string | null
  /** Rebuild this panel on new primer settings. Optional: a page that cannot build renders
   *  the same settings read-only, which is what they are to a panel already built. */
  onRebuild?: (build: PrimerBuild) => void
  /** Start a UCSC check of every designed pair. Undefined where the server has no key, and
   *  then the box says so rather than offering a button that cannot work. */
  onVerify?: () => void
  verify?: { running: boolean; stage: string; error?: string }
}

export function PanelTable({ result, ancestry, onRebuild, onVerify, verify }: Props) {
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

  // As above, and for the same reason: a pair that in-silico PCR contradicts must not be
  // hideable by paging past it. Walked over the candidates rather than the primer map, so
  // the rsID and the panel's own order come with it.
  const dangerous = useMemo(
    () => candidates.flatMap((m) => {
      const why = m.primer && pcrDanger(m.primer)
      return why ? [{ m, why }] : []
    }),
    [candidates],
  )

  // The rule gates the glyph and its words together: no rule, no stars at all, whatever the
  // rows happen to carry. A star the page cannot explain is worse than no star.
  const rule = flankingRule(result.provenance)

  // No recorded settings, no primer UI at all: the server has no primer module, and an
  // empty form over an absent feature invites a rebuild that would change nothing.
  const build = primerBuild(result)

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
            gnomAD and Ensembl report different GRCh38 coordinates, so every figure on the row
            rests on a position only one of them agrees with. Resolve it at dbSNP before
            ordering an assay, or drop{' '}
            {disputed.length === 1 ? 'the marker' : 'them'}.{' '}
            <Anchor href="#/docs/conventions" size="xs">Why this happens</Anchor>
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

      {dangerous.length > 0 && (
        <Alert
          color="red"
          mx={8}
          mt={8}
          role="alert"
          title={
            dangerous.length === 1
              ? 'In-silico PCR contradicts a primer pair in this panel'
              : `In-silico PCR contradicts ${dangerous.length} primer pairs in this panel`
          }
        >
          <Text size="xs" mb={6}>
            A pair that does not amplify exactly the marker it was designed for cannot
            genotype it.{' '}
            {dangerous.length === 1 ? 'It is' : 'They are'} still shown, under{' '}
            {dangerous.length === 1 ? 'its marker' : 'their markers'}, with the finding.{' '}
            <Anchor href="#/docs/primers" size="xs">What this means</Anchor>
          </Text>
          <List size="xs" spacing={2}>
            {dangerous.map(({ m, why }) => (
              <List.Item key={m.variant_id}>
                <span className="om-mono">{m.rsid}</span>: {why.short}
              </List.Item>
            ))}
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
        {build && <PrimerChip build={build} onRebuild={onRebuild} onVerify={onVerify} verify={verify} />}
        {sort && (
          <Anchor size="xs" component="button" type="button" onClick={() => setSort(null)}>
            clear column sort
          </Anchor>
        )}
        <Text size="xs" c="dimmed" ml="auto">
          {int(sorted.length)} rows
        </Text>
      </Group>

      {rule && (
        <Group gap={6} px={8} pb={8} wrap="nowrap" align="flex-start">
          <span style={{ lineHeight: 0, paddingTop: 2 }}><Star /></span>
          {/* A link, but not link-blue: it sits under a table as a key, and colouring it
              would pull the eye off the panel. The affordance is the hover. */}
          <a href={rule.docs_href ?? '#/docs'} className="om-star-key">
            <Text size="xs" c="dimmed" component="span">{legendWords(rule.legend)}</Text>
          </a>
        </Group>
      )}

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
              const design = m.primer
              return (
                <Fragment key={m.variant_id}>
                <Table.Tr
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
                      {rule && starred(m) && (
                        <Tooltip label={rule.summary ?? rule.legend} withArrow multiline w={300}>
                          <a
                            href={rule.docs_href ?? '#/docs'}
                            // The claim, and only the claim: an <a href> is already announced as a link, so
                            // saying so here would make a screen reader read "link" twice.
                            aria-label={STAR_CLAIM}
                            style={{ lineHeight: 0, display: 'inline-flex' }}
                          >
                            <Star />
                          </a>
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
                {design && <PrimerDetail d={design} />}
                </Fragment>
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

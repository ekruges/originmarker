import { useRef, useState, type ReactNode } from 'react'
import {
  Alert, Badge, Button, Group, Paper, Progress, SegmentedControl, Table, Text,
} from '@mantine/core'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, gates, headerMap, parseRow,
  type ChromStats, type Gate, type ProbeRow, type SampleProfile,
} from './ingest'
import { type AB } from './informativity'
import {
  agreement, classify, emptyTally, GLOSS, isAutosome, pair, pct, tallyRow,
  type ChromResult, type PairResult, type ParentageResult,
} from './parentage'
import type { Health } from './api'
import { int, utc } from './fmt'
import { EXAMPLES, EXAMPLE_CITATION, EXAMPLE_MARKERS, loadExample } from './examples'
import { analyseRuns, parseLocus, type RunResult } from './runlength'
import type { Marker } from './informativity'
import { buildReportPdf, reportId, sha256, type ReportFile } from './syngamyPdf'
import { syngamyLogText } from './logfile'
import { FeatureHeader, DropZone } from './FeatureHeader'
import { RunLog } from './RunLog'

/**
 * Syngamy - whether the two gametic genomes fused, and which parts of each survived.
 *
 * Files are read in this browser and there is no endpoint to send them to: these are a family's
 * arrays, and the terms promise nothing about a family is submitted or retained. The donor is
 * held as one call per marker and each sample streams against it, so memory stays flat.
 */

type Tag = 'READ' | 'PARSE' | 'CALL' | 'WARN' | 'DONE'
interface Line { tag: Tag; text: string }

const TAG_COLOR: Record<Tag, string> = {
  READ: 'var(--om-text-dim)',
  PARSE: 'var(--om-text-dim)',
  CALL: 'var(--om-blue)',
  WARN: 'var(--om-higher)',
  DONE: 'var(--om-blue)',
}

type Role = 'donor' | 'oocyte' | 'sample'
type State = 'profiling' | 'waiting' | 'running' | 'done' | 'failed'

interface Entry {
  id: string
  file: File
  role: Role
  state: State
  profile?: SampleProfile
  gates?: Gate[]
  result?: ParentageResult
  maternal?: ParentageResult
  paired?: PairResult
  error?: string
}

/** One parent held as a call per marker, plus the heterozygosity the second-parent axis needs.
 *  ponytail: two of these is ~130 MB on 825k-marker arrays. Fine on a desktop; if a laptop
 *  starts thrashing, the fix is a packed typed array keyed on a sorted probe list, not a cap. */
interface DonorIndex { gt: Map<string, AB>; heterozygosity: number; build: string | null }

/** A finished per-locus test, lifted to the page so the report can carry it. */
export interface LocusRun {
  chrom: string
  pos: number
  eventSizeBp?: number
  bySample: { name: string; results: RunResult[] }[]
}
interface Stage { id: string; markers: number; bytes: number; total: number }

const mb = (b: number): string => `${(b / 1e6).toFixed(1)} MB`

/** Hand a blob to the browser as a download. The only way anything leaves this page. */
function grab(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url
  a.download = filename
  document.body.appendChild(a)
  a.click()
  a.remove()
  URL.revokeObjectURL(url)
}


async function eachLine(
  file: File, fn: (line: string) => void, onChunk?: (bytesRead: number) => void,
): Promise<void> {
  const reader = file.stream().getReader()
  const decoder = new TextDecoder()
  let carry = ''
  let read = 0
  for (;;) {
    const { done, value } = await reader.read()
    if (done) break
    read += value.byteLength
    carry += decoder.decode(value, { stream: true })
    let i = carry.indexOf('\n')
    while (i >= 0) {
      fn(carry.slice(0, i))
      carry = carry.slice(i + 1)
      i = carry.indexOf('\n')
    }
    onChunk?.(read)
    // Yield, so the log and the bars paint while a large file streams.
    await new Promise((r) => { setTimeout(r, 0) })
  }
  if (carry) fn(carry)
}

async function profileFile(
  file: File,
  onRow: ((r: ProbeRow) => void) | null,
  tick: (s: Partial<Stage>) => void,
  log: (tag: Tag, text: string) => void,
): Promise<{ profile: SampleProfile; gates: Gate[] }> {
  const byChrom = new Map<string, ChromStats>()
  const baf = emptyBafSums()
  const builds = emptyBuildSums()
  let map: ReturnType<typeof headerMap> = null
  let firstId = ''
  let n = 0
  let skipped = 0
  log('READ', `${file.name}  ${mb(file.size)}`)
  await eachLine(file, (line) => {
    if (!map) {
      map = headerMap(line)
      if (map) log('PARSE', `header resolved, ${Object.keys(map).length - 1} columns mapped`)
      return
    }
    const r = parseRow(line, map)
    if (!r) { skipped += 1; return }
    if (!firstId) firstId = r.probesetId
    accumulate(r, byChrom)
    accumulateBaf(r, baf)
    accumulateBuild(r, builds)
    onRow?.(r)
    n += 1
    if (n % 100_000 === 0) { tick({ markers: n }); log('PARSE', `${int(n)} markers`) }
  }, (bytes) => tick({ bytes }))
  if (!map) throw new Error('no recognisable header row')
  tick({ markers: n })
  if (skipped) log('WARN', `${int(skipped)} rows unparseable, skipped`)

  const profile: SampleProfile = {
    ...finishProfile(file.name, byChrom, baf, firstId, builds), build: buildVerdict(builds),
  }
  const g = gates(profile)
  log('PARSE', `${int(n)} markers, ${byChrom.size} chromosomes`)
  log('CALL', `call ${pct(profile.callRate, 1)}, het ${pct(profile.hetRate, 1)}, `
    + `sex ${profile.sex}, ${profile.build.build ?? 'assembly undetermined'}`)
  for (const x of g) {
    if (x.verdict === 'exclude' || x.verdict === 'marginal') {
      log('WARN', `gate ${x.name}: ${x.verdict}`)
    }
  }
  return { profile, gates: g }
}

export function SyngamyPage({ health }: { health?: Health | null }) {
  const [entries, setEntries] = useState<Entry[]>([])
  const [donor, setDonor] = useState<DonorIndex | null>(null)
  const [oocyte, setOocyte] = useState<DonorIndex | null>(null)
  const [stage, setStage] = useState<Stage | null>(null)
  const [lines, setLines] = useState<Line[]>([])
  const [busy, setBusy] = useState(false)
  const [startedAt, setStartedAt] = useState<string | null>(null)
  const [examples, setExamples] = useState(false)
  const [saving, setSaving] = useState(false)
  const [locus, setLocus] = useState<LocusRun | null>(null)
  const pick = useRef<HTMLInputElement>(null)

  const log = (tag: Tag, text: string) => setLines((p) => [...p.slice(-499), { tag, text }])
  const patch = (id: string, d: Partial<Entry>) =>
    setEntries((p) => p.map((e) => (e.id === id ? { ...e, ...d } : e)))

  const add = async (files: FileList | File[], roleFor?: (name: string) => Role | undefined) => {
    // Built OUTSIDE the state updater. Populating a local array from inside one and then looping
    // over it reads an array React has not filled yet, and every chip sits on "profiling".
    const seen = new Set(entries.map((e) => e.id))
    const fresh: Entry[] = []
    for (const f of [...files]) {
      const id = `${f.name}:${f.size}`
      if (seen.has(id)) continue
      seen.add(id)
      // The first file in is the likely donor and the rest samples. A guess, with the control
      // sitting on the chip to correct it.
      fresh.push({
        id,
        file: f,
        role: roleFor?.(f.name) ?? (entries.length + fresh.length === 0 ? 'donor' : 'sample'),
        state: 'profiling',
      })
    }
    if (!fresh.length) return
    setEntries((prev) => [...prev, ...fresh])

    // Profile immediately, so the chip carries a shape and a quality read before anything is
    // compared. Nothing here needs the donor, which is why it can run on drop.
    for (const e of fresh) {
      setStage({ id: e.id, markers: 0, bytes: 0, total: e.file.size })
      try {
        const { profile, gates: g } = await profileFile(
          e.file, null,
          (x) => setStage((p) => ({
            id: e.id,
            total: e.file.size,
            markers: x.markers ?? (p?.id === e.id ? p.markers : 0),
            bytes: x.bytes ?? (p?.id === e.id ? p.bytes : 0),
          })),
          log,
        )
        patch(e.id, { state: 'waiting', profile, gates: g })
      } catch (err) {
        const m = err instanceof Error ? err.message : String(err)
        log('WARN', `${e.file.name}: ${m}`)
        patch(e.id, { state: 'failed', error: m })
      }
    }
    setStage(null)
  }

  const loadExamples = async () => {
    setBusy(true)
    setExamples(true)
    log('READ', `fetching ${EXAMPLES.length} public example arrays, GEO GSE148488`)
    try {
      const files = await Promise.all(EXAMPLES.map((x) => loadExample(x)))
      const roles = new Map(EXAMPLES.map((x) => [x.file.replace(/\.gz$/, ''), x.role]))
      setBusy(false)
      await add(files, (n) => roles.get(n))
    } catch (e) {
      setBusy(false)
      log('WARN', `examples: ${e instanceof Error ? e.message : String(e)}`)
    }
  }

  const save = async () => {
    setSaving(true)
    try {
      const files: ReportFile[] = []
      for (const e of entries) {
        files.push({
          name: e.file.name,
          size: e.file.size,
          sha256: await sha256(e.file),
          role: e.role,
          markers: e.profile?.markers ?? 0,
          profile: e.profile,
          gates: e.gates,
          result: e.result,
          maternal: e.maternal,
          paired: e.paired,
          error: e.error,
        })
      }
      const rank = { donor: 0, oocyte: 1, sample: 2 }
      files.sort((a, b) => rank[a.role] - rank[b.role])
      const generatedAt = `${new Date().toISOString().replace('T', ' ').slice(0, 19)}Z`
      const id = reportId(files.map((f) => f.sha256).join('') + (startedAt ?? generatedAt))
      const blob = await buildReportPdf({
        files,
        donorHeterozygosity: donor?.heterozygosity ?? NaN,
        startedAt,
        generatedAt,
        tool: health ? `OriginMarker ${health.version} (${health.release_codename})`
          : 'OriginMarker',
        reportId: id,
        fromExamples: examples,
        locus,
      })
      grab(blob, `syngamy-report-${id}.pdf`)
    } finally {
      setSaving(false)
    }
  }

  const run = async () => {
    setBusy(true)
    setStartedAt(new Date().toISOString())
    let index = donor
    let maternalIndex = oocyte

    const bar = (id: string, total: number) => (s: Partial<Stage>) =>
      setStage((p) => ({
        id,
        total,
        markers: s.markers ?? (p?.id === id ? p.markers : 0),
        bytes: s.bytes ?? (p?.id === id ? p.bytes : 0),
      }))

    // Both parents are indexed the same way and differ only in what they are then compared as,
    // which is the point: a maternal absence has to be measured by the same instrument as a
    // paternal one or the two rates cannot be set against each other.
    const indexParent = async (e: Entry, what: string): Promise<DonorIndex | null> => {
      patch(e.id, { state: 'running' })
      setStage({ id: e.id, markers: 0, bytes: 0, total: e.file.size })
      try {
        const gt = new Map<string, AB>()
        let called = 0
        let het = 0
        const { profile, gates: g } = await profileFile(e.file, (r) => {
          gt.set(r.probesetId, r.genotype)
          if (r.genotype !== 'NC' && isAutosome(r.chrom)) {
            called += 1
            if (r.genotype === 'AB') het += 1
          }
        }, bar(e.id, e.file.size), log)
        const h = called ? het / called : NaN
        log('DONE', `${what} indexed, ${int(gt.size)} markers, autosomal het ${pct(h, 3)}`)
        patch(e.id, { state: 'done', profile, gates: g })
        return { gt, heterozygosity: h, build: profile.build.build }
      } catch (err) {
        const m = err instanceof Error ? err.message : String(err)
        log('WARN', `${e.file.name}: ${m}`)
        patch(e.id, { state: 'failed', error: m })
        return null
      }
    }

    const donorEntry = entries.find((e) => e.role === 'donor')
    if (donorEntry && donorEntry.state !== 'done') {
      index = await indexParent(donorEntry, 'sperm donor')
      setDonor(index)
    }
    const oocyteEntry = entries.find((e) => e.role === 'oocyte')
    if (oocyteEntry && oocyteEntry.state !== 'done') {
      maternalIndex = await indexParent(oocyteEntry, 'oocyte donor')
      setOocyte(maternalIndex)
    }

    // Computed once, not per sample: it is a property of the two parents, not of the embryo.
    const agree = index && maternalIndex ? agreement(index.gt, maternalIndex.gt) : NaN
    if (agree > 0.99) {
      log('WARN', `the two declared parents agree at ${pct(agree, 1)} of shared markers, which `
        + 'is one person in two files rather than two people')
    }

    if (index) {
      const pat = index
      const mat = maternalIndex
      for (const s of entries.filter((e) => e.role === 'sample' && e.state !== 'done')) {
        patch(s.id, { state: 'running' })
        setStage({ id: s.id, markers: 0, bytes: 0, total: s.file.size })
        try {
          const t = emptyTally()
          const tm = mat ? emptyTally() : null
          // The parent's assembly, so the pseudoautosomal boundaries are the right ones. Known
          // before the sample streams, and the same across one experiment's arrays.
          t.build = pat.build
          if (tm) tm.build = mat?.build ?? null
          const { profile, gates: g } = await profileFile(s.file, (r) => {
            tallyRow(pat.gt.get(r.probesetId) ?? 'NC', r, t)
            if (tm && mat) tallyRow(mat.gt.get(r.probesetId) ?? 'NC', r, tm)
          }, bar(s.id, s.file.size), log)
          const result = classify(t, pat.heterozygosity)
          const maternal = tm && mat
            ? classify(tm, mat.heterozygosity, { role: 'maternal' }) : undefined
          const paired = maternal ? pair(result, maternal, agree) : undefined
          log('DONE', `${s.file.name}: ${paired?.originClass ?? result.originClass}, paternal `
            + `absent ${pct(result.genomeRate)} vs ceiling ${pct(result.explainable)}`
            + (maternal ? `, maternal absent ${pct(maternal.genomeRate)} vs ceiling `
              + `${pct(maternal.explainable)}` : ''))
          patch(s.id, { state: 'done', profile, gates: g, result, maternal, paired })
        } catch (e) {
          const m = e instanceof Error ? e.message : String(e)
          log('WARN', `${s.file.name}: ${m}`)
          patch(s.id, { state: 'failed', error: m })
        }
      }
    }
    setStage(null)
    setBusy(false)
  }

  const pending = entries.some((e) => e.state === 'waiting')
  const hasDonor = entries.some((e) => e.role === 'donor')

  return (
    <>
      <FeatureHeader name="Syngamy" tagline="parent of origin from SNP arrays" />
      <Paper p="sm" mb={10}>
        <Group justify="flex-end" align="center" mb={8}>
          <Group gap={6}>
            <input
              ref={pick} type="file" multiple style={{ display: 'none' }}
              onChange={(e) => { if (e.target.files?.length) void add(e.target.files) }}
            />
            <Button variant="default" size="xs" disabled={busy}
              onClick={() => pick.current?.click()}
            >
              Add files&hellip;
            </Button>
            {entries.length === 0 && (
              <Button variant="default" size="xs" disabled={busy}
                onClick={() => { void loadExamples() }}
              >
                Examples
              </Button>
            )}
            <Button size="xs" disabled={busy || !pending || !hasDonor} onClick={() => { void run() }}>
              Run
            </Button>
            {entries.some((e) => e.result) && (
              <Button variant="default" size="xs" disabled={busy || saving}
                onClick={() => { void save() }}
              >
                {saving ? 'Building\u2026' : 'Report (PDF)'}
              </Button>
            )}
            {entries.length > 0 && (
              <Button
                variant="subtle" size="xs" disabled={busy}
                onClick={() => {
                  setEntries([]); setDonor(null); setOocyte(null); setLines([])
                  setStartedAt(null); setExamples(false); setLocus(null)
                }}
              >
                Clear
              </Button>
            )}
          </Group>
        </Group>

        <DropZone
          empty={entries.length === 0}
          prompt="Drop array exports here. Label one donor, the rest samples."
          onFiles={(f) => { void add(f) }}
        >
            <div style={{ display: 'flex', flexWrap: 'wrap', gap: 6 }}>
              {entries.map((e) => (
                <Chip
                  key={e.id} entry={e} busy={busy} stage={stage?.id === e.id ? stage : null}
                  onRole={(role) => patch(e.id, {
                    role,
                    state: e.profile ? 'waiting' : e.state,
                    result: undefined,
                    maternal: undefined,
                    paired: undefined,
                  })}
                  onRemove={() => setEntries((p) => p.filter((x) => x.id !== e.id))}
                />
              ))}
          </div>
        </DropZone>

        {entries.filter((e) => e.role === 'oocyte').length > 1 && (
          <Alert color="orange" mt={8} p="xs">
            <Text size="xs">Only one file can be the oocyte donor.</Text>
          </Alert>
        )}
        {entries.length > 0 && !hasDonor && (
          <Alert color="orange" mt={8} p="xs">
            <Text size="xs">
              One file must be labelled sperm. An oocyte donor is optional and measures the
              maternal side directly instead of inferring it.
            </Text>
          </Alert>
        )}
      </Paper>

      {entries.filter((e) => e.result).map((e) => (
        <ResultCard
          key={e.id} entry={e}
          donorName={entries.find((x) => x.role === 'donor')?.file.name ?? ''}
          oocyteName={entries.find((x) => x.role === 'oocyte')?.file.name ?? ''}
        />
      ))}

      {entries.some((e) => e.result) && (
        <LocusTest
          entries={entries} donor={donor} oocyte={oocyte} log={log}
          onResult={setLocus}
        />
      )}

      {lines.length > 0 && (
        <RunLog
          lines={lines}
          colours={TAG_COLOR}
          onDownload={() => grab(
            new Blob([syngamyLogText(lines, {
              tool: health ? `OriginMarker ${health.version} (${health.release_codename})`
                : 'OriginMarker',
              started: startedAt,
            }, {
              now: new Date().toISOString(),
              url: window.location.href,
              agent: navigator.userAgent,
            })], { type: 'text/plain;charset=utf-8' }),
            `syngamy-runlog-${new Date().toISOString().slice(0, 19).replace(/[:T]/g, '')}.txt`,
          )}
        />
      )}
      <RunInformation
        entries={entries} startedAt={startedAt} health={health} examples={examples}
      />
    </>
  )
}

/* --- chip ----------------------------------------------------------------------------------- */

const STATE_COLOUR: Record<State, string> = {
  profiling: 'gray', waiting: 'gray', running: 'genomeBlue', done: 'genomeBlue', failed: 'red',
}

function Chip({ entry, busy, stage, onRole, onRemove }: {
  entry: Entry
  busy: boolean
  stage: Stage | null
  onRole: (r: Role) => void
  onRemove: () => void
}) {
  const { file, role, state, profile, result } = entry
  const verdict = entry.paired?.originClass ?? result?.originClass
  const running = state === 'profiling' || state === 'running'

  return (
    <div
      className="om-chip"
      style={{
        border: '1px solid var(--om-border)',
        borderLeft: `3px solid ${verdict && verdict !== 'unclear' ? 'var(--om-blue)'
          : verdict === 'unclear' ? 'var(--om-higher)'
            : state === 'failed' ? '#c92a2a' : 'var(--om-border-strong)'}`,
        borderRadius: 2,
        background: '#fff',
        width: 318,
        maxWidth: '100%',
        padding: '6px 8px',
      }}
    >
      <Group gap={5} wrap="nowrap" align="center" mb={4}>
        <Text size="xs" fw={600} truncate style={{ flex: 1, minWidth: 0 }} title={file.name}>
          {file.name}
        </Text>
        <Badge
          size="xs" variant="light" style={{ flex: 'none', maxWidth: 'none' }}
          color={verdict ? (verdict === 'unclear' ? 'orange' : 'genomeBlue') : STATE_COLOUR[state]}
        >
          {verdict?.replace(/_/g, ' ') ?? state}
        </Badge>
        <button
          type="button" onClick={onRemove} disabled={busy || running} aria-label="Remove"
          style={{
            border: 0, background: 'none', cursor: 'pointer', padding: 0, lineHeight: 1,
            color: 'var(--om-text-dim)', flex: 'none', fontSize: 13,
          }}
        >
          ×
        </button>
      </Group>

      {/* On upload this is the sample's own per-chromosome heterozygosity; after a run it becomes
          paternal absence against that sample's noise ceiling. Both answer "what shape is this
          genome" at a glance, which is what the chip is for. */}
      <div style={{ height: 42, marginBottom: 4 }}>
        {running ? (
          <Progress value={stage?.total ? (stage.bytes / stage.total) * 100 : 0} size="sm" animated />
        ) : result ? (
          <ChromStrip chroms={result.chroms} ceiling={result.explainable} />
        ) : profile ? (
          <ProfileStrip profile={profile} />
        ) : null}
      </div>

      <Group gap={4} wrap="nowrap" align="center" className="om-chip-foot">
        <SegmentedControl
          size="xs" value={role} disabled={busy || running} onChange={(v) => onRole(v as Role)}
          data={[{ label: 'sperm', value: 'donor' }, { label: 'oocyte', value: 'oocyte' },
            { label: 'sample', value: 'sample' }]}
          styles={{ root: { flex: 'none' }, label: { padding: '1px 6px', fontSize: 10 } }}
        />
        <Text
          size="xs" c="dimmed" ff="monospace" ta="right"
          style={{ flex: 1, minWidth: 0, whiteSpace: 'nowrap' }}
        >
          {running && stage
            ? `${int(stage.markers)} · ${mb(stage.bytes)}/${mb(stage.total)}`
            : result
              ? `${pct(result.genomeRate)} of ${pct(result.explainable)}`
              : profile
                ? `${int(profile.markers)} mk · ${pct(profile.callRate, 0)} · ${profile.sex}`
                : mb(file.size)}
        </Text>
      </Group>

      {entry.error && <Text size="xs" c="red" mt={2}>{entry.error}</Text>}
    </div>
  )
}

/** One measured parent, as the line of numbers the call was made on. */
function Axis({ label, r }: { label: string; r: ParentageResult }) {
  const margin = r.genomeRate / r.explainable
  return (
    <Text size="xs" c="dimmed" ff="monospace">
      {label && `${label} `}absent {pct(r.genomeRate)} &middot; ceiling {pct(r.explainable)} &middot;{' '}
      {Number.isFinite(margin) ? `${margin.toFixed(1)}x` : '-'} &middot;{' '}
      {int(r.informative)} informative
      <br />
      {/* The ceiling's own two factors, so a ratio that fell is readable as the ceiling rising
          rather than the signal shrinking. They are what it is the product of. */}
      &nbsp;&nbsp;from {pct(r.noCallRate, 1)} no-call &times; {pct(r.hetFraction, 1)} het
      {Number.isFinite(r.dispersion) && (
        <>
          {' '}&middot; spread {r.dispersion.toFixed(2)} across chromosomes, cleanest{' '}
          {pct(r.minChromRate)}
        </>
      )}
    </Text>
  )
}

/** The verdict, at the size it deserves, with the detail folded underneath it. */
function ResultCard({ entry, donorName, oocyteName }: {
  entry: Entry; donorName: string; oocyteName: string
}) {
  const [open, setOpen] = useState(false)
  const r = entry.result!
  const m = entry.maternal
  // With both arrays the class comes from the pair: each parent is measured, so a maternal
  // origin is a presence rather than the shape left behind by a paternal absence.
  const originClass = entry.paired?.originClass ?? r.originClass
  const decided = originClass !== 'unclear'
  const gloss = originClass === 'neither_parent'
    ? 'Neither declared parent accounts for this genome.'
    : GLOSS[originClass as keyof typeof GLOSS]

  return (
    <Paper p={0} mb={12}>
      <div style={{
        padding: '13px 16px',
        borderBottom: open ? '1px solid var(--mantine-color-gray-3)' : undefined,
        background: decided
          ? 'var(--mantine-color-genomeBlue-0)' : 'var(--mantine-color-orange-0)',
      }}
      >
        <Text size="xs" c="dimmed" mb={2}>
          {entry.file.name} &middot; vs {donorName}{m && ` and ${oocyteName}`}
        </Text>
        <Group gap={10} align="baseline">
          <Text
            style={{ fontSize: 26, fontWeight: 700, letterSpacing: '-0.01em', lineHeight: 1.1 }}
            c={decided ? 'genomeBlue.9' : 'orange.9'}
          >
            {originClass.replace(/_/g, ' ').toUpperCase()}
          </Text>
          {r.spermType !== 'unknown' && (
            <Badge size="sm" variant="outline" color="genomeGrey">
              {r.spermType.replace('_', '-')} sperm
            </Badge>
          )}
          <Badge size="sm" variant="outline" color="genomeGrey">
            {r.zygosity.replace(/_/g, ' ')}
          </Badge>
          {m && (
            <Badge size="sm" variant="light" color="genomeBlue">both parents measured</Badge>
          )}
        </Group>
        <Text size="sm" mt={3}>{gloss}</Text>
        <div style={{ marginTop: 6 }}>
          <Axis label={m ? 'paternal' : ''} r={r} />
          {m && <Axis label="maternal" r={m} />}
        </div>
        <Group gap={12} align="center" mt={4}>
          <button
            type="button" className="om-primer-toggle" aria-expanded={open}
            onClick={() => setOpen((o) => !o)}
          >
            <span className="om-primer-toggle-caret">{open ? '▾' : '▸'}</span>
            <span>{open ? 'Hide detail' : 'Detail'}</span>
          </button>
        </Group>
      </div>
      {open && (
        <div style={{ padding: '2px 16px 14px' }}>
          {entry.paired && entry.paired.notes.length > 0 && (
            <Section title="Both parents">
              {entry.paired.notes.map((n) => <Text key={n} size="xs" mb={3}>{n}</Text>)}
              {Number.isFinite(entry.paired.agreement) && (
                <Text size="xs" c="dimmed" mt={3}>
                  The two parents agree at {pct(entry.paired.agreement, 1)} of shared markers.
                </Text>
              )}
            </Section>
          )}
          <SampleDetail
            result={r} maternal={m} profile={entry.profile!} gates={entry.gates!}
          />
        </div>
      )}
    </Paper>
  )
}
function DotTrack({ points, rules, band, height = 42, ariaLabel }: {
  points: { key: string; y: number; weight: number; over: boolean; title: string }[]
  rules: { y: number }[]
  /** Everything at or below this height is inside the explainable range. Shaded, so a point
   *  clearing it is visible as leaving a region rather than as crossing a thin line. */
  band?: number
  height?: number
  ariaLabel: string
}) {
  if (!points.length) return null
  const W = 320
  const H = height
  const pad = 5
  const span = H - pad * 2
  const maxW = Math.max(...points.map((p) => p.weight), 1)
  const step = W / points.length
  const yPix = (v: number) => pad + span * (1 - Math.max(0, Math.min(1, v)))
  const floor = yPix(0)

  return (
    <svg
      viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="none" role="img" aria-label={ariaLabel}
      style={{ width: '100%', height: H, display: 'block' }}
    >
      {points.map((p, i) => (i % 2 === 1 ? (
        <rect key={`z${p.key}`} x={i * step} y={0} width={step} height={H} fill="var(--om-zebra)" />
      ) : null))}

      {band !== undefined && (
        <rect
          x={0} y={yPix(band)} width={W} height={floor - yPix(band)}
          fill="var(--om-blue)" opacity={0.07}
        />
      )}

      {rules.map((r, i) => (
        <line
          key={i} x1={0} y1={yPix(r.y)} x2={W} y2={yPix(r.y)}
          stroke="var(--om-border-strong)" strokeWidth={0.6} strokeDasharray="3,2"
        />
      ))}

      {/* Stems. A dot alone floats; a stem ties it to the baseline and the profile of the whole
          genome becomes one silhouette. */}
      {points.map((p, i) => (
        <line
          key={`s${p.key}`}
          x1={i * step + step / 2} y1={floor} x2={i * step + step / 2} y2={yPix(p.y)}
          stroke={p.over ? 'var(--om-higher)' : 'var(--om-blue)'}
          strokeWidth={1} opacity={p.over ? 0.5 : 0.25}
        />
      ))}

      {points.map((p, i) => {
        const r = 1.6 + 2.2 * Math.sqrt(p.weight / maxW)
        return (
          <g key={p.key}>
            {p.over && (
              <circle
                cx={i * step + step / 2} cy={yPix(p.y)} r={r + 1.8}
                fill="none" stroke="var(--om-higher)" strokeWidth={0.8} opacity={0.45}
              />
            )}
            <circle
              cx={i * step + step / 2} cy={yPix(p.y)} r={r}
              fill={p.over ? 'var(--om-higher)' : 'var(--om-blue)'}
              opacity={p.over ? 1 : 0.75}
              stroke="#fff" strokeWidth={0.5}
            >
              <title>{p.title}</title>
            </circle>
          </g>
        )
      })}
      <line x1={0} y1={floor} x2={W} y2={floor} stroke="var(--om-border)" strokeWidth={0.6} />
    </svg>
  )
}

/** After a run: paternal absence as a multiple of this sample's own noise ceiling. */
function ChromStrip({ chroms, ceiling }: { chroms: ChromResult[]; ceiling: number }) {
  if (!chroms.length || !Number.isFinite(ceiling) || ceiling <= 0) return null
  const CAP = 10
  return (
    <DotTrack
      ariaLabel="per-chromosome paternal absence as a multiple of the noise ceiling"
      rules={[{ y: 1 / CAP }, { y: 3 / CAP }]}
      band={3 / CAP}
      points={chroms.map((c) => {
        const mult = c.rate / ceiling
        return {
          key: c.chrom,
          y: Math.min(1, mult / CAP),
          weight: c.informative,
          over: c.verdict !== 'present',
          title: `chr${c.chrom}: ${pct(c.rate, 2)} absent, ${mult.toFixed(1)}x ceiling, `
            + `${int(c.informative)} informative`,
        }
      })}
    />
  )
}

/** On upload: this sample's own per-chromosome heterozygosity, before anything is compared. A
 *  uniparental genome already sits flat along the floor here. */
function ProfileStrip({ profile }: { profile: SampleProfile }) {
  const chroms = [...profile.byChrom.entries()]
    .filter(([, v]) => v.called > 50)
    .sort((a, b) => a[0].length - b[0].length || a[0].localeCompare(b[0]))
  if (!chroms.length) return null
  // 2pq is bounded above by 0.5, so that is the axis top rather than a chosen one.
  const TOP = 0.5
  return (
    <DotTrack
      ariaLabel="per-chromosome heterozygosity"
      rules={[{ y: profile.hetRate / TOP }]}
      points={chroms.map(([c, v]) => ({
        key: c,
        y: (v.het / v.called) / TOP,
        weight: v.called,
        over: false,
        title: `chr${c}: het ${pct(v.het / v.called, 1)}, ${int(v.called)} called`,
      }))}
    />
  )
}

/* --- detail ---------------------------------------------------------------------------------- */

/** A table wider than the column it sits in scrolls inside its own box rather than pushing
 *  the card off the screen. Every table below is wider than a phone. */
const Scroll = ({ children }: { children: ReactNode }) => (
  <div className="om-scroll-x">{children}</div>
)

function Section({ title, children }: { title: string; children: ReactNode }) {
  return (
    <div style={{ marginTop: 10 }}>
      <Text fw={600} size="xs" tt="uppercase" c="dimmed" mb={3} style={{ letterSpacing: '0.04em' }}>
        {title}
      </Text>
      {children}
    </div>
  )
}

function Facts({ items }: { items: [string, ReactNode][] }) {
  return (
    <Scroll>
      <Table withTableBorder withColumnBorders>
        <Table.Tbody>
          {items.map(([k, v]) => (
            <Table.Tr key={k}>
              <Table.Td style={{ width: '30%' }}>{k}</Table.Td>
              <Table.Td ff="monospace">{v}</Table.Td>
            </Table.Tr>
          ))}
        </Table.Tbody>
      </Table>
    </Scroll>
  )
}

const GATE_COLOUR: Record<string, string> = {
  usable: 'green', marginal: 'orange', exclude: 'red', report_only: 'gray',
}

function Quality({ profile: p, gates: g }: { profile: SampleProfile; gates: Gate[] }) {
  return (
    <>
      <Facts items={[
        ['Markers', int(p.markers)],
        ['Called', `${int(p.called)}   ${pct(p.callRate, 1)}`],
        ['No-call', pct(p.nocallRate, 2)],
        ['Heterozygous, of called', pct(p.hetRate, 2)],
        ['chrX / autosomal het', p.chrXHetRatio === null ? '-' : p.chrXHetRatio.toFixed(3)],
        ['Sex', p.sex],
        ['Product', p.product],
        ['Assembly', p.build.build
          ? `${p.build.build}   ${int(p.build.illegal.GRCh37)} illegal GRCh37, `
            + `${int(p.build.illegal.GRCh38)} GRCh38, ${int(p.build.tested)} tested`
          : `undetermined   ${p.build.note}`],
        ['Coding', `${p.coding.verdict}   ${p.coding.note}`],
        ['Mean BAF, AA / AB / BB', [p.coding.meanBafHom0, p.coding.meanBafHet, p.coding.meanBafHom2]
          .map((x) => (x === null ? '-' : x.toFixed(3))).join('   ')],
      ]}
      />
      <Scroll>
        <Table withTableBorder withColumnBorders mt={6}>
          <Table.Thead>
            <Table.Tr>
              <Table.Th style={{ width: '24%' }}>Gate</Table.Th>
              <Table.Th ta="right" style={{ width: '10%' }}>Value</Table.Th>
              <Table.Th style={{ width: '13%' }}>Verdict</Table.Th>
              <Table.Th>Basis</Table.Th>
            </Table.Tr>
          </Table.Thead>
          <Table.Tbody>
            {g.map((x) => (
              <Table.Tr key={x.name}>
                <Table.Td>{x.name}</Table.Td>
                <Table.Td ta="right" ff="monospace">{x.value === null ? '-' : pct(x.value, 1)}</Table.Td>
                <Table.Td>
                  <Badge size="xs" variant="light" color={GATE_COLOUR[x.verdict] ?? 'gray'}>
                    {x.verdict.replace(/_/g, ' ')}
                  </Badge>
                </Table.Td>
                <Table.Td><Text size="xs" c="dimmed">{x.detail}</Text></Table.Td>
              </Table.Tr>
            ))}
          </Table.Tbody>
        </Table>
      </Scroll>
    </>
  )
}

function SampleDetail({ result: r, maternal, profile, gates: g }: {
  result: ParentageResult; maternal?: ParentageResult; profile: SampleProfile; gates: Gate[]
}) {
  const sex = r.chroms.filter((c) => c.chrom === 'X' || c.chrom === 'Y' || c.chrom === '23')
  return (
    <>
      <Section title={maternal ? 'Evidence, sperm donor' : 'Evidence'}>
        <Scroll>
          <Table withTableBorder withColumnBorders>
            <Table.Thead>
              <Table.Tr>
                <Table.Th style={{ width: '30%' }}>Measurement</Table.Th>
                <Table.Th ta="right" style={{ width: '14%' }}>Observed</Table.Th>
                <Table.Th ta="right" style={{ width: '14%' }}>Reference</Table.Th>
                <Table.Th>Reference is</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              <Row label="Paternal alleles absent" obs={pct(r.genomeRate)}
                ref_={pct(r.explainable)} basis="no-call rate x heterozygous fraction"
              />
              <Row label="Alleles the donor lacks" obs={pct(r.nonParentalRate)}
                ref_={pct(r.secondParentExpected)} basis="half the donor's heterozygosity"
              />
              <Row label="Heterozygous BAF band" obs={pct(r.hetBand)} ref_="8%"
                basis="diploid 15-16%, uniparental 1.3-3.4%"
              />
            </Table.Tbody>
          </Table>
        </Scroll>
      </Section>

      {maternal && (
        <Section title="Evidence, oocyte donor">
          <Scroll>
            <Table withTableBorder withColumnBorders>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th style={{ width: '30%' }}>Measurement</Table.Th>
                  <Table.Th ta="right" style={{ width: '14%' }}>Observed</Table.Th>
                  <Table.Th ta="right" style={{ width: '14%' }}>Reference</Table.Th>
                  <Table.Th>Reference is</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <Row label="Maternal alleles absent" obs={pct(maternal.genomeRate)}
                  ref_={pct(maternal.explainable)} basis="no-call rate x heterozygous fraction"
                />
                <Row label="Alleles the donor lacks" obs={pct(maternal.nonParentalRate)}
                  ref_={pct(maternal.secondParentExpected)} basis="half the donor's heterozygosity"
                />
              </Table.Tbody>
            </Table>
          </Scroll>
        </Section>
      )}

      {sex.length > 0 && (
        <Section title="Sex chromosomes">
          <ChromTable chroms={sex} ceiling={r.explainable} />
        </Section>
      )}

      <Section title={maternal ? 'Chromosomes, sperm donor' : 'Chromosomes'}>
        <ChromTable chroms={r.chroms} ceiling={r.explainable} />
      </Section>

      {maternal && (
        <Section title="Chromosomes, oocyte donor">
          <ChromTable chroms={maternal.chroms} ceiling={maternal.explainable} />
        </Section>
      )}

      <Section title="Sample quality"><Quality profile={profile} gates={g} /></Section>

      {[...new Set([...r.limits, ...maternal?.limits ?? []])].length > 0 && (
        <Section title="Limits">
          {[...new Set([...r.limits, ...maternal?.limits ?? []])].map((l) => (
            <Alert key={l} color="orange" p="xs" mb={4}><Text size="xs">{l}</Text></Alert>
          ))}
        </Section>
      )}
      {[...new Set([...r.notes, ...maternal?.notes ?? []])].length > 0 && (
        <Section title="Findings">
          {[...new Set([...r.notes, ...maternal?.notes ?? []])].map((n) => (
            <Text key={n} size="xs" mb={3}>{n}</Text>
          ))}
        </Section>
      )}
    </>
  )
}

function ChromTable({ chroms, ceiling }: { chroms: ChromResult[]; ceiling: number }) {
  return (
    <Scroll>
      <Table withTableBorder withColumnBorders>
        <Table.Thead>
          <Table.Tr>
            <Table.Th>Chr</Table.Th>
            <Table.Th ta="right">Absent</Table.Th>
            <Table.Th ta="right">Missing / informative</Table.Th>
            <Table.Th ta="right">x ceiling</Table.Th>
            <Table.Th>Verdict</Table.Th>
          </Table.Tr>
        </Table.Thead>
        <Table.Tbody>
          {chroms.map((c) => (
            <Table.Tr key={c.chrom}>
              <Table.Td>chr{c.chrom}</Table.Td>
              <Table.Td ta="right" ff="monospace">{pct(c.rate, 2)}</Table.Td>
              <Table.Td ta="right" ff="monospace" c="dimmed">
                {int(c.absent)} / {int(c.informative)}
              </Table.Td>
              <Table.Td ta="right" ff="monospace" c="dimmed">{(c.rate / ceiling).toFixed(1)}</Table.Td>
              <Table.Td>
                {c.verdict === 'present'
                  ? <Text span size="xs" c="dimmed">present</Text>
                  : (
                    <Badge size="xs" variant="light"
                      color={c.verdict === 'expected_absent' ? 'gray' : 'orange'}
                    >
                      {c.verdict.replace(/_/g, ' ')}
                    </Badge>
                  )}
                {c.note && <Text span size="xs" c="dimmed"> &middot; {c.note}</Text>}
              </Table.Td>
            </Table.Tr>
          ))}
        </Table.Tbody>
      </Table>
    </Scroll>
  )
}

function Row({ label, obs, ref_, basis }: {
  label: string; obs: string; ref_: string; basis: string
}) {
  return (
    <Table.Tr>
      <Table.Td>{label}</Table.Td>
      <Table.Td ta="right" ff="monospace" fw={600}>{obs}</Table.Td>
      <Table.Td ta="right" ff="monospace" c="dimmed">{ref_}</Table.Td>
      <Table.Td><Text size="xs" c="dimmed">{basis}</Text></Table.Td>
    </Table.Tr>
  )
}

/* --- per-locus deletion test ---------------------------------------------------------------- */

const VERDICT_COLOUR: Record<string, string> = {
  significant_run: 'orange',
  no_significant_run: 'gray',
  below_resolution: 'gray',
  undefined_father_heterozygous: 'gray',
}


/**
 * The per-locus test, asked of a position rather than of the genome.
 *
 * The genome-wide call says whether a paternal contribution arrived. This asks whether it is
 * specifically absent AROUND ONE SITE, from the length of the run of consecutive markers where
 * the father's obligate allele is undetected.
 *
 * It requires the oocyte donor and refuses without her, for the reason `origin.py` refuses:
 * paternal PRESENCE cannot be established at any single marker without someone who could not
 * have supplied the allele, so the informative set would hold absences only, nothing could break
 * a run, and the statistic would have no null to be significant against.
 *
 * The sample is re-read rather than retained. Its genotypes were streamed and discarded to keep
 * memory flat, and reading one chromosome again costs less than holding every marker of every
 * sample against the chance that someone asks this question.
 */
function LocusTest({ entries, donor, oocyte, log, onResult }: {
  entries: Entry[]
  donor: DonorIndex | null
  oocyte: DonorIndex | null
  log: (tag: Tag, text: string) => void
  onResult: (r: LocusRun | null) => void
}) {
  const [text, setText] = useState('')
  const [size, setSize] = useState('')
  const [rows, setRows] = useState<{ name: string; results: RunResult[] }[] | null>(null)
  const [busy, setBusy] = useState(false)
  const [err, setErr] = useState('')

  const samples = entries.filter((e) => e.result && e.state === 'done')
  const locus = parseLocus(text)

  const run = async () => {
    if (!locus || !donor || !oocyte) return
    setBusy(true)
    setErr('')
    setRows(null)
    // Dropped before the run, not only on failure: a report built mid-test must not carry the
    // previous locus beside the new one's heading.
    onResult(null)
    const eventSizeOfInterestBp = size.trim() ? Number(size.replace(/[,_ ]/g, '')) : undefined
    log('CALL', `per-locus test at chr${locus.chrom}:${int(locus.pos)}, `
      + `${samples.length} sample(s)`)
    try {
      const out: { name: string; results: RunResult[] }[] = []
      for (const e of samples) {
        const markers: Marker[] = []
        const embryo = new Map<string, AB>()
        await profileFile(e.file, (r) => {
          if (r.chrom.replace(/^chr/i, '').toUpperCase() !== locus.chrom) return
          markers.push({ rsid: r.probesetId, chrom: r.chrom, pos: r.pos })
          embryo.set(r.probesetId, r.genotype)
        }, () => {}, () => {})
        // q is the rate at which a marker shows paternal absence while the paternal genome IS
        // present, and this sample's own genome-wide absence rate measures exactly that. The
        // Kothiyal floor is a LOWER bound on it, so leaving q at the floor on an amplified sample
        // shrinks r_min and makes significance easier to reach: it over-calls. Usable only where
        // the paternal genome is present genome-wide, since where it is absent that rate is not a
        // spurious rate at all; qHat falls back to the floor for those.
        const qEmpirical = e.result?.verdict === 'parent_genome_present'
          ? e.result.genomeRate : null
        out.push({
          name: e.file.name,
          results: analyseRuns(markers, donor.gt, oocyte.gt, embryo, locus.chrom, locus.pos, {
            qEmpirical,
            ...(eventSizeOfInterestBp !== undefined && Number.isFinite(eventSizeOfInterestBp)
              ? { eventSizeOfInterestBp } : {}),
          }),
        })
        log('DONE', `${e.file.name}: ${out[out.length - 1].results
          .map((r) => `${r.window} ${r.verdict.replace(/_/g, ' ')}`).join(', ')}`)
      }
      setRows(out)
      onResult({
        chrom: locus.chrom,
        pos: locus.pos,
        eventSizeBp: eventSizeOfInterestBp !== undefined
          && Number.isFinite(eventSizeOfInterestBp) ? eventSizeOfInterestBp : undefined,
        bySample: out,
      })
    } catch (x) {
      const m = x instanceof Error ? x.message : String(x)
      setErr(m)
      log('WARN', `per-locus test: ${m}`)
    } finally {
      setBusy(false)
    }
  }

  return (
    <Paper p="sm" mb={10}>
      <Text fw={600} size="xs" tt="uppercase" c="dimmed" mb={4} style={{ letterSpacing: '0.04em' }}>
        Per-locus deletion test
      </Text>
      {!oocyte ? (
        <Alert color="orange" p="xs">
          <Text size="xs">
            This test needs the oocyte donor and is not run without her. Paternal presence cannot
            be established at any single marker without someone who could not have supplied the
            allele: the informative set would hold absences only, nothing could break a run, and
            the run-length statistic would have no null to be significant against. Label one file
            oocyte and run again. Parent of origin, sperm type and segmental loss above need no
            oocyte donor and are unaffected.
          </Text>
        </Alert>
      ) : (
        <>
          <Group gap={6} align="flex-end" wrap="wrap" mb={6}>
            <input
              value={text}
              onChange={(ev) => setText(ev.target.value)}
              onKeyDown={(ev) => { if (ev.key === 'Enter' && locus && !busy) void run() }}
              placeholder="chr7:117559590"
              aria-label="Variant position"
              className="om-mono"
              style={{
                flex: '1 1 190px', minWidth: 0, fontSize: 13, padding: '5px 8px',
                border: '1px solid var(--om-border-strong)', borderRadius: 2,
              }}
            />
            <input
              value={size}
              onChange={(ev) => setSize(ev.target.value)}
              placeholder="event bp (optional)"
              aria-label="Event size of interest, in base pairs"
              className="om-mono"
              style={{
                flex: '0 1 150px', minWidth: 0, fontSize: 13, padding: '5px 8px',
                border: '1px solid var(--om-border-strong)', borderRadius: 2,
              }}
            />
            <Button size="xs" disabled={!locus || busy} onClick={() => { void run() }}>
              {busy ? 'Reading\u2026' : 'Test'}
            </Button>
          </Group>
          <Text size="xs" c="dimmed" mb={6}>
            The position of the variant, in this array&apos;s assembly. An event size turns an
            absent run into an explicit &quot;below resolution&quot; rather than letting it read
            as no event.
          </Text>
          {text && !locus && <Text size="xs" c="orange">Not a position: try chr7:117559590.</Text>}
          {err && <Text size="xs" c="red">{err}</Text>}
          {rows?.map((r) => (
            <div key={r.name} style={{ marginTop: 10 }}>
              <Text size="xs" fw={600} mb={3}>{r.name}</Text>
              <Scroll>
                <Table withTableBorder withColumnBorders>
                  <Table.Thead>
                    <Table.Tr>
                      <Table.Th>Window</Table.Th>
                      <Table.Th ta="right">L3 markers</Table.Th>
                      <Table.Th ta="right">Longest run</Table.Th>
                      <Table.Th ta="right">r_min</Table.Th>
                      <Table.Th ta="right">p</Table.Th>
                      <Table.Th>Verdict</Table.Th>
                    </Table.Tr>
                  </Table.Thead>
                  <Table.Tbody>
                    {r.results.map((w) => (
                      <Table.Tr key={w.window}>
                        <Table.Td>{w.window.replace(/_/g, ' ')}</Table.Td>
                        <Table.Td ta="right" ff="monospace">{int(w.nL3)}</Table.Td>
                        <Table.Td ta="right" ff="monospace">{int(w.longestRun)}</Table.Td>
                        <Table.Td ta="right" ff="monospace">{w.rMin ?? '-'}</Table.Td>
                        <Table.Td ta="right" ff="monospace">
                          {w.nL3 ? w.runP.toExponential(1) : '-'}
                        </Table.Td>
                        <Table.Td>
                          <Badge
                            size="xs" variant="light"
                            color={VERDICT_COLOUR[w.verdict] ?? 'gray'}
                          >
                            {w.verdict.replace(/_/g, ' ')}
                          </Badge>
                        </Table.Td>
                      </Table.Tr>
                    ))}
                  </Table.Tbody>
                </Table>
              </Scroll>
              {r.results.map((w) => (
                <Text key={w.window} size="xs" c="dimmed" mt={3}>
                  <b>{w.window.replace(/_/g, ' ')}:</b> {w.note}
                </Text>
              ))}
            </div>
          ))}
        </>
      )}
    </Paper>
  )
}

/* --- log and run information ------------------------------------------------------------------ */


/** What a citation needs and nothing else. The reasoning behind the numbers is in the PDF. */
function RunInformation({ entries, startedAt, health, examples }: {
  entries: Entry[]
  startedAt: string | null
  health?: Health | null
  examples: boolean
}) {
  const done = entries.filter((e) => e.state === 'done')
  if (!done.length) return null
  const samples = done.filter((e) => e.result)
  const donorEntry = done.find((e) => e.role === 'donor')
  const builds = [...new Set(done.map((e) => e.profile?.build.build ?? 'undetermined'))]

  return (
    <Paper p="sm" mb={10}>
      <Text fw={600} size="xs" tt="uppercase" c="dimmed" mb={4} style={{ letterSpacing: '0.04em' }}>
        Run information
      </Text>
      <Facts items={[
        ['Started', startedAt ? utc(startedAt) : '-'],
        ['Tool', health
          ? `OriginMarker ${health.version} (${health.release_codename})` : 'OriginMarker'],
        ['Donor', donorEntry?.file.name ?? '-'],
        ['Samples', String(samples.length)],
        ['Assembly', builds.join(', ')],
        ...(examples
          ? [['Example data', `${EXAMPLE_CITATION} Every eighth marker, `
              + `${int(EXAMPLE_MARKERS)} per file, values unaltered.`] as [string, ReactNode]]
          : []),
      ]}
      />
    </Paper>
  )
}

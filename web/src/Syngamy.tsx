import { breathe } from './scan'
import {
  collectParentRow, collectRow, emptyCollected, emptyParent, finishParent, scoreSample,
  type ParentIndex,
} from './scoreSample'
import { syngamyTable } from './syngamyTable'
import { ParentalBalancePanel } from './ParentalBalancePanel'
import type { BalanceSample, PairedSample } from './parentalBalance'
import { useEffect, useRef, useState, type ReactNode } from 'react'
import { StageCallout } from './StageCallout'
import { ComparisonPanel } from './ComparisonPanel'
import {
  normaliseTrack, pooledRegions, pooledMarkers, type ComparisonResult,
} from './comparison'
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
  agreement, classify, GLOSS, isAutosome, pair, pct,
  type ChromResult, type PairResult, type ParentageResult,
} from './parentage'
import {
  segmentCoords, MIN_SEGMENT_MARKERS, SEGMENT_LRT,
  type Segment, type SegmentKind,
} from './segments'
import type { Health } from './api'
import { int, utc } from './fmt'
import { EXAMPLES, EXAMPLE_CITATION, EXAMPLE_MARKERS, loadExample } from './examples'
import { analyseRuns, measureClustering, parseLocus, type RunResult } from './runlength'
import type { GainAnnotation } from './parentage'
import { type FeatureTrack } from './features'
import type { Marker } from './informativity'
import { buildReportPdf, reportId, sha256, type ReportFile } from './syngamyPdf'
import { isInferredFile } from './inferredArray'
import { receiveHandoff, wantsHandoff } from './handoff'
import { syngamyLogText } from './logfile'
import { FeatureHeader, DropZone } from './FeatureHeader'
import { RunLog } from './RunLog'
import { DefectCallout } from './DefectCallout'
import { defectsFrom, findingToDefect, parentNamed, withMechanism } from './defects'
import type { DosageVerdict } from './dosageOrigin'
import type { OneParentVerdict } from './oneParentOrigin'
import { groupUnits, unitsCarrying, callUniformity } from './abnormalities'
import { stageFacts } from './stage'
/**
 * Syngamy - whether the two gametic genomes fused, and which parts of each survived.
 *
 * Files are read in this browser and there is no endpoint to send them to: these are a family's
 * arrays, and the terms promise nothing about a family is submitted or retained. The donor is
 * held as one call per marker and each sample streams against it, so memory stays flat.
 */

type Tag = 'READ' | 'PARSE' | 'CALL' | 'WARN' | 'DONE' | 'SCAN'
interface Line { tag: Tag; text: string }

const TAG_COLOR: Record<Tag, string> = {
  READ: 'var(--om-text-dim)',
  PARSE: 'var(--om-text-dim)',
  CALL: 'var(--om-blue)',
  WARN: 'var(--om-higher)',
  DONE: 'var(--om-blue)',
  SCAN: 'var(--om-text-dim)',
}

/**
 * An F-grade call for a finding that covers the whole genome rather than an interval.
 *
 * There is no interval to measure, so the evidence is the array-wide absence of the loaded parent,
 * which is the right evidence for a genome-scoped statement and the wrong evidence for a located
 * one. Graded F because it is a direction and not a measurement of anything in particular.
 */
const genomeGrade = (r: {
  genomeRate?: number; explainable?: number
}): { verdict: string; shift: number; z: number; impliedF: number; why: string
  confidence: number; band: string } | undefined => {
  const { genomeRate, explainable } = r
  if (!Number.isFinite(genomeRate) || !Number.isFinite(explainable)) return undefined
  return {
    verdict: (genomeRate as number) > (explainable as number) ? 'other-parent' : 'loaded-parent',
    // No interval, so there is no interval statistic. Left absent rather than filled with the
    // array-wide numbers, which would read as a measurement of this finding.
    shift: NaN,
    z: NaN,
    impliedF: NaN,
    confidence: 0.5,
    band: 'F',
    why: 'Graded F and named anyway: this finding covers the whole genome and has no interval to '
      + `measure, so the direction is the array's own, from absence of the loaded parent at `
      + `${((genomeRate as number) * 100).toFixed(2)}% against a `
      + `${((explainable as number) * 100).toFixed(2)}% ceiling. Do not report or count this row`,
  }
}

/**
 * The change list for one sample, built once so the panel and the export cannot disagree.
 *
 * Previously assembled inline in the render, which meant any machine-readable output would have
 * had to reassemble it and could have drifted from what the reader was looking at.
 */
/** The last marker position seen on a chromosome, for a whole-chromosome event's extent. */
function chromEndOf(r: ParentageResult, chrom: string): number {
  const pos = r.markerPositions?.get(chrom)
  return pos && pos.length ? pos[pos.length - 1] : 0
}

function defectsForResult(r: ParentageResult) {
  // WHOLE-CHROMOSOME ANEUPLOIDY IS A CHROMOSOMAL CHANGE, and belongs in the one list with the rest.
  // It used to have a panel of its own, so a reader comparing a whole-chromosome loss against a
  // segmental one was comparing two boxes with different vocabularies and no shared ordering.
  const aneu = (r.chroms ?? []).filter((c) => c.aneuploidy).map((c) => {
    const scored = (r.dosageCalls ?? []).find((d: { where: string }) => d.where === `chr${c.chrom}`)
    const other = r.role === 'paternal' ? 'maternal' : 'paternal'
    // THE MENDELIAN CHANNEL FIRST, measured at 0.920 on obvious whole-chromosome events against
    // the dosage channel's band D on the same material. It asks whether an allele is there at all
    // rather than whether a mean has moved, so it needs no detection floor.
    const mendel = (r.oneParent ?? []).find((o: { where: string }) => o.where === `chr${c.chrom}`)
    const mendelParent = mendel
      ? parentNamed(mendel.verdict as OneParentVerdict, r.role) : null
    const origin = mendelParent
      ?? (c.aneuploidyParent === 'this' ? r.role
        : c.aneuploidyParent === 'other' ? other
          : (scored && parentNamed(scored.verdict as DosageVerdict, r.role)) || 'unclear')
    return {
      chrom: c.chrom,
      startBp: 0,
      endBp: chromEndOf(r, c.chrom),
      kind: c.aneuploidy === 'gain' ? 'copy-gain' : 'copy-loss',
      locus: `chr${c.chrom}`,
      origin,
      band: mendelParent ? mendel?.band : scored?.band,
      confidence: mendelParent ? mendel?.posterior : scored?.confidence,
      inheritedMargin: (scored as { inheritedMargin?: number } | undefined)?.inheritedMargin,
      stage: r.stage?.stage,
      why: `calls at ${c.callFraction.toFixed(2)}x the genome rate with intensity `
        + `${c.lrrShift > 0 ? '+' : ''}${c.lrrShift.toFixed(2)} log2 from the rest. An intact `
        + 'chromosome calls at 0.78x to 1.16x of its genome median and never leaves -0.79 to +0.42 '
        + `log2, measured over 1,012 chromosomes.`
        + `${mendelParent ? ` ${mendel?.why}` : ''}${scored?.why ? ` ${scored.why}` : ''}`,
    } as never
  })

  return withMechanism([
    ...aneu,
      ...defectsFrom(r.segments, r.gains, r.losses ?? [], r.oneParent ?? [], r.role,
        r.stage, r.dosageCalls ?? []),
      // One list, not two. A copy-neutral event and a deletion are different measurements
      // of the same kind of thing and a reader compares them against each other.
      ...(r.findings ?? []).map((f) => findingToDefect(f, r.stage,
        // The call scored over THIS finding's own interval, matched by the same label the
        // scorer wrote, so a finding shows the origin that was actually measured for it
        // rather than one borrowed from its chromosome.
        // A GENOME-SCOPED FINDING HAS NO INTERVAL TO SCORE, so it never reached the
        // scorer and was the one row still leaving without a grade. Its evidence is the
        // whole array's, which is exactly what a genome-scoped finding should be judged
        // on, so it is graded F from the same array-wide absence measurement the other
        // ungradable rows fall back to.
        f.chrom === 'genome'
          ? genomeGrade(r)
          : (r.dosageCalls ?? []).find((d: { where: string }) => d.where
            === `chr${f.chrom} ${(f.startBp / 1e6).toFixed(1)}-${(f.endBp / 1e6).toFixed(1)}Mb`),
        r.role)),
    ], r.uniformity)
}

/**
 * The whole run reduced to what the parental comparison needs, one row per sample.
 *
 * ONLY MEASURED OR INHERITED PARENTS COUNT. A row graded F names no parent at all, and one whose
 * class blocks an origin never had one, so neither can contribute a side. The sample's parent is
 * the genome-level class where it has one, because that is what makes a uniparental genome the
 * clean comparison in the first place.
 */
function balanceRowsFor(entries: readonly Entry[]): BalanceSample[] {
  const out: BalanceSample[] = []
  for (const e of entries) {
    const r = e.result
    if (!r) continue
    const parent = r.originClass === 'gynogenetic' ? 'maternal' as const
      : r.originClass === 'androgenetic' ? 'paternal' as const : null
    // A sample that lost its class still declares which group it WOULD have joined, so an
    // exclusion correlated with quality is visible rather than silently repairing the balance.
    const declaredParent = parent ?? (r.role === 'paternal' ? 'maternal' as const : 'paternal' as const)
    const events = defectsForResult(r)
      .filter((d) => d.origin && d.origin !== 'unclear' && d.band !== 'F')
      .map((d) => ({ cls: d.kind, chrom: d.chrom, startBp: d.startBp, endBp: d.endBp }))
    out.push({
      name: e.file.name,
      parent,
      declaredParent,
      originClass: r.originClass,
      informative: r.informative,
      explainable: r.explainable,
      material: r.stage?.stage ?? 'unknown',
      events,
    })
  }
  return out
}

/**
 * Biparental samples reduced to a within-sample maternal/paternal count.
 *
 * ONLY MEASURED ORIGINS COUNT HERE, and that is what keeps the design honest. The whole value of
 * comparing within one array is that both sides were read on the same chemistry in the same
 * reaction, so nothing about the array can favour one over the other. An inherited origin comes
 * from a genome-level call and would put every event in a sample on the same side; a direction-only
 * one carries no side at all. Either would import the between-sample confounds this design exists
 * to remove.
 */
function pairedRowsFor(entries: readonly Entry[]): PairedSample[] {
  const out: PairedSample[] = []
  for (const e of entries) {
    const r = e.result
    if (!r || r.originClass !== 'biparental') continue
    const measured = defectsForResult(r).filter((d) => d.band === 'A' || d.band === 'B'
      || d.band === 'C' || d.band === 'D')
    const maternalEvents = measured.filter((d) => d.origin === 'maternal').length
    const paternalEvents = measured.filter((d) => d.origin === 'paternal').length
    if (!maternalEvents && !paternalEvents) continue
    out.push({ name: e.file.name, maternalEvents, paternalEvents })
  }
  return out
}

/** How many log lines are kept. The oldest are dropped, so a run cannot grow the page without end. */
const LOG_LINES = 500

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
  /** This file is a genotype Progenitor reconstructed, not an array anyone measured. */
  inferred?: boolean
}

/** One parent held as a call per marker. The shape is scoreSample's, so the two cannot diverge.
 *  ponytail: two of these is ~130 MB on 825k-marker arrays. Fine on a desktop; if a laptop
 *  starts thrashing, the fix is a packed typed array keyed on a sorted probe list, not a cap. */
type DonorIndex = ParentIndex

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
    + `BAF spread at het calls ${Number.isFinite(profile.hetBafSd)
      ? profile.hetBafSd.toFixed(4) : 'n/a'}, `
    + `sex ${profile.sex}, ${profile.build.build ?? 'assembly undetermined'}`)
  for (const x of g) {
    if (x.verdict === 'exclude' || x.verdict === 'marginal') {
      log('WARN', `gate ${x.name}: ${x.verdict}`)
    }
  }
  return { profile, gates: g }
}

/**
 * hg19 feature intervals, fetched once and shared. 13 KB gzipped, so it is fetched rather than
 * bundled, and a failure is not an error: placement is an addition to a run, never a gate on one.
 */
let featureTrack: FeatureTrack | null = null
let trackTried = false
async function loadFeatureTrack(): Promise<FeatureTrack | null> {
  if (featureTrack || trackTried) return featureTrack
  trackTried = true
  try {
    const r = await fetch(`${import.meta.env.BASE_URL}hg19_features.json`)
    // Normalised rather than cast. The shipped file stores intervals as tuples while every
    // consumer reads object fields, so a cast produced a track that matched nothing and said so
    // with a p value.
    if (r.ok) featureTrack = normaliseTrack(await r.json())
  } catch { /* offline or blocked; the run is unaffected */ }
  return featureTrack
}

export function SyngamyPage({ health }: { health?: Health | null }) {
  const [entries, setEntries] = useState<Entry[]>([])
  /** Run-wide, because the comparison is a cohort question rather than a per-sample one. */
  const [comparison, setComparison] = useState<ComparisonResult | undefined>(undefined)
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

  // BUFFERED, because a line used to cost a React render and a 500-element array copy, and a run
  // writes hundreds of them in a burst. Three hundred lines was three hundred renders and 150,000
  // object copies, which is time spent on the log rather than on the analysis it describes.
  // The buffer flushes on a microtask, so a burst inside one chunk of work lands as one update and
  // is on screen before the next chunk starts.
  const logBuffer = useRef<Line[]>([])
  const flushLog = () => {
    if (!logBuffer.current.length) return
    const add = logBuffer.current
    logBuffer.current = []
    setLines((p) => [...p, ...add].slice(-LOG_LINES))
  }
  const downloadTable = () => {
    const rows = entries.filter((e) => e.result).map((e) => ({
      name: e.file.name,
      originClass: e.result!.originClass,
      zygosity: e.result!.zygosity,
      role: e.result!.role,
      genomeRate: e.result!.genomeRate,
      explainable: e.result!.explainable,
      informative: e.result!.informative,
      stage: e.result!.stage?.stage,
      changes: defectsForResult(e.result!).map((d) => ({
        chrom: d.chrom, startBp: d.startBp, endBp: d.endBp, kind: d.kind,
        origin: d.origin, band: d.band, confidence: d.confidence,
        inheritedMargin: d.inheritedMargin, stage: d.stage, informative: d.informative,
        why: d.why,
      })),
    }))
    const blob = new Blob([syngamyTable(rows)], { type: 'text/tab-separated-values' })
    const a = document.createElement('a')
    a.href = URL.createObjectURL(blob)
    a.download = `syngamy-changes-${new Date().toISOString().slice(0, 10)}.tsv`
    a.click()
    URL.revokeObjectURL(a.href)
  }

  const log = (tag: Tag, text: string) => {
    logBuffer.current.push({ tag, text })
    if (logBuffer.current.length === 1) queueMicrotask(flushLog)
  }
  const patch = (id: string, d: Partial<Entry>) =>
    setEntries((p) => p.map((e) => (e.id === id ? { ...e, ...d } : e)))

  // Opened from Progenitor with an array in hand. The reconstruction is passed across in
  // memory rather than through storage or a file picker, and lands here as the donor.
  const [handoff, setHandoff] = useState(() => wantsHandoff(window.location.hash))
  useEffect(() => {
    const stop = receiveHandoff((file, note) => {
      log('READ', `handed over from Progenitor: ${note}`)
      setHandoff(false)
      void add([file], () => 'donor')
    })
    return stop
    // Once, on mount. The opener only answers the first request.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

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

    // A reconstructed parent is a legitimate donor here and is the reason Progenitor writes one.
    // It is not a measured array, and every artefact of this run has to say which it was, so the
    // mark is read before anything is profiled rather than trusted to the file name.
    for (const e of fresh) {
      if (await isInferredFile(e.file)) {
        patch(e.id, { inferred: true })
        log('WARN', `${e.file.name}: this is a RECONSTRUCTED genotype, not a measured array. `
          + 'Every call made against it inherits that, and the report says so throughout.')
      }
    }

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
          inferred: e.inferred,
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
        // The run-wide feature comparison, if someone ran it. Without this line the field exists,
        // the drawing code exists, and the section never prints: input.comparison is always
        // undefined and the bundling silently does not happen.
        comparison,
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
        const pacc = emptyParent()
        const { profile, gates: g } = await profileFile(
          e.file, (r) => collectParentRow(r, pacc), bar(e.id, e.file.size), log,
        )
        // THE ONE REFUSAL THAT SHOULD STOP A RUN. Every downstream channel measures against the
        // loaded parent, so a parent array that failed its own gates cannot support any of them,
        // and continuing produces a report whose every row is a refusal traceable to this line.
        const parentBlocked = (g ?? []).filter((x) => x.verdict === 'exclude')
        if (parentBlocked.length) {
          log('WARN', `${what} failed ${parentBlocked.length} of its own gates: `
            + `${parentBlocked.map((x) => x.name).join(', ')}. Every channel below measures `
            + 'against this array, so nothing downstream can be trusted and the run stops here '
            + 'rather than producing a report of refusals')
          throw new Error(`parent array failed ${parentBlocked.length} gate(s)`)
        }
        const idx = finishParent(pacc, profile.build.build)
        log('DONE', `${what} indexed, ${int(idx.gt.size)} markers, `
          + `autosomal het ${pct(idx.heterozygosity, 3)}`)
        patch(e.id, { state: 'done', profile, gates: g })
        return idx
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

    // EITHER parent alone is enough. Nothing below this line is paternal except the label it is
    // given: the tally, the segments and the Mendelian exclusion all ask whether the LOADED
    // parent's copy is present, and which parent that is only decides what the answer is called.
    // Gating the whole run on a sperm donor meant an oocyte-only run silently did nothing at all.
    const solo = index ?? maternalIndex
    const soloRole: 'paternal' | 'maternal' = index ? 'paternal' : 'maternal'
    if (solo) {
      const pat = solo
      const mat = index ? maternalIndex : undefined
      // Genotypes retained per sample so that a sample's SIBLINGS in the same run can establish
      // which markers the embryo is heterozygous at. This is what removes the parental-array
      // requirement: unaffected cells of one embryo share both parents exactly. Kept as a plain
      // map per sample rather than streamed twice, since the run already holds every file.
      const sampleGt = new Map<string, Map<string, string>>()
      /** Finished samples, kept locally because patch() writes React state we cannot read back. */
      const finished: { id: string; result: ParentageResult }[] = []
      const queue = entries.filter((e) => e.role === 'sample' && e.state !== 'done')
      let doneCount = 0
      for (const s of queue) {
        // Three things an operator can act on: which file, how many remain, and any refusal that
        // changes the final answer. Everything else belongs in the report.
        log('READ', `sample ${doneCount + 1} of ${queue.length}: ${s.file.name}`)
        doneCount += 1
        patch(s.id, { state: 'running' })
        setStage({ id: s.id, markers: 0, bytes: 0, total: s.file.size })
        try {
          // ONE IMPLEMENTATION OF THE ANSWER, in scoreSample.ts, which the command line and the
          // cross-surface check also call. What is left here is what only a browser has to do:
          // read a File, paint progress, and hold React state.
          const acc = emptyCollected(pat, mat)
          const { profile, gates: g } = await profileFile(
            s.file, (r) => collectRow(r, pat, mat, acc), bar(s.id, s.file.size), log,
          )
          sampleGt.set(s.id, acc.myGt)
          const result = await scoreSample({
            acc,
            profile,
            pat,
            mat,
            soloRole,
            // Which arrays are units of one embryo is a property of the run, not of this sample.
            sibs: [...sampleGt.entries()].filter(([id]) => id !== s.id).map(([, m]) => m),
            sampleName: s.file.name,
            log,
          })
          const maternal = acc.tm && mat
            ? classify(acc.tm, mat.heterozygosity, { role: 'maternal' }) : undefined
          const paired = maternal ? pair(result, maternal, agree) : undefined
          log('DONE', `${s.file.name}: ${paired?.originClass ?? result.originClass}, ${soloRole} `
            + `absent ${pct(result.genomeRate)} vs ceiling ${pct(result.explainable)}`
            + (maternal ? `, maternal absent ${pct(maternal.genomeRate)} vs ceiling `
              + `${pct(maternal.explainable)}` : ''))
          finished.push({ id: s.id, result })
          await breathe()
          patch(s.id, { state: 'done', profile, gates: g, result, maternal, paired })
        } catch (e) {
          const m = e instanceof Error ? e.message : String(e)
          log('WARN', `${s.file.name}: ${m}`)
          patch(s.id, { state: 'failed', error: m })
        }
      }

      // UNITS OF ONE EMBRYO, AND THE TIMING QUESTION THEY SETTLE.
      //
      // Whether a segmental change came from the gamete or arose after fertilisation is not
      // separable on genotype at any material quality. Uniformity across independently sampled
      // units of the same embryo is the only channel measured to break the tie: meiotic 64 of 64
      // uniform, post-zygotic 6 of 7 non-uniform. The tool cannot ask for a second biopsy, but a
      // user who already arrayed two units gets the answer for free.
      //
      // THE GROUPING IS MEASURED RATHER THAN DECLARED. Two biopsies of one embryo are the same
      // genome and concord at 95.8% of shared called markers, against 54.9% for a parent-offspring
      // pair, so nobody has to label anything and a labelling slip cannot produce a confidently
      // wrong mechanism.
      if (finished.length > 1) {
        const groups = groupUnits(finished, (a, b) => {
          const ga = sampleGt.get(a.id)
          const gb = sampleGt.get(b.id)
          return ga && gb ? agreement(ga as never, gb as never) : NaN
        })
        for (const group of groups) {
          if (group.length < 2) continue
          const perUnit = group.map((u) => [
            ...(u.result.segments ?? []).map((sg: { chrom: string }) => {
              const co = segmentCoords(sg as never)
              return { chrom: sg.chrom, startBp: co.start, endBp: co.end }
            }),
            ...(u.result.findings ?? []).map((f) => ({
              chrom: f.chrom, startBp: f.startBp, endBp: f.endBp,
            })),
          ])
          log('DONE', `${group.length} arrays are units of one embryo, by genotype concordance. `
            + 'Segmental changes can now be separated into gamete-borne and post-zygotic')
          for (const u of group) {
            u.result.units = group.length
            u.result.uniformity = [
              ...(u.result.segments ?? []).map((sg: { chrom: string }) => {
                const co = segmentCoords(sg as never)
                return { chrom: sg.chrom, startBp: co.start, endBp: co.end }
              }),
              ...(u.result.findings ?? []).map((f) => ({
                chrom: f.chrom, startBp: f.startBp, endBp: f.endBp,
              })),
            ].filter((e) => e.chrom !== 'genome').map((e) => {
              const carried = unitsCarrying(e, perUnit)
              const call = callUniformity(carried, group.length)
              return { ...e, mechanism: call.mechanism, why: call.why }
            })
            patch(u.id, { result: { ...u.result } })
            for (const m of u.result.uniformity) {
              log(m.mechanism === 'unresolved' ? 'WARN' : 'DONE',
                `timing chr${m.chrom}: ${m.mechanism}. ${m.why}`)
            }
          }
        }
      }
    }
    setStage(null)
    setBusy(false)
  }

  const pending = entries.some((e) => e.state === 'waiting')
  // Either parent alone is enough to run. The analysis asks whether the LOADED parent's copy is
  // present, which is the same question whichever parent that is, so requiring a sperm donor
  // turned an oocyte-only run into a disabled button with no explanation.
  const hasDonor = entries.some((e) => e.role === 'donor')
  const hasParent = hasDonor || entries.some((e) => e.role === 'oocyte')

  return (
    <>
      {handoff && (
        <Alert color="blue" p="xs" mb={8}>
          <Text size="xs">
            Waiting for the reconstructed array from Progenitor. If nothing arrives, the tab that
            opened this one was closed: save the array there and drop it in here as the donor.
          </Text>
        </Alert>
      )}
      <FeatureHeader name="Syngamy" tagline="parent of origin from SNP arrays" />

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
            <Button size="xs" disabled={busy || !pending || !hasParent} onClick={() => { void run() }}>
              Run
            </Button>
            {entries.some((e) => e.result) && (
              <Button variant="default" size="xs" disabled={busy || saving}
                onClick={() => { void save() }}
              >
                {saving ? 'Building\u2026' : 'Report (PDF)'}
              </Button>
            )}
            {/* THE JOINABLE OUTPUT, beside the record rather than behind it. A reader of this tool
                joins these calls against their own sample metadata, which a PDF cannot support. */}
            {entries.some((e) => e.result) && (
              <Button variant="default" size="xs" disabled={busy}
                onClick={() => { downloadTable() }}
                title="One row per change, tab separated"
              >
                Table (TSV)
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
        {entries.length > 0 && !hasParent && (
          <Alert color="orange" mt={8} p="xs">
            <Text size="xs">
              One file must be labelled as a parent, sperm or oocyte. Either alone is enough: the
              run asks whether that parent&rsquo;s copy is present, and names the other parent by
              exclusion where the alleles allow it. Loading both measures each side directly
              rather than inferring one.
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

      {/*
        ONE COMPARISON OVER THE WHOLE RUN, NOT ONE PER SAMPLE. A single chip contributes a handful
        of regions, and under five the matched null has no shape at all: the answer would be "no
        conclusion" on almost every card. Pooled across every sample the run becomes a dataset, and
        the question is a cohort question anyway. Do the changes THIS EXPERIMENT found sit where the
        genome breaks anyway? That is not a property of one embryo.

        The regions carry their sample name so the grid stays readable and a coincidence can be
        traced back to the array it came from.
      */}
      {/* THE QUESTION THE TOOL EXISTS FOR, across the whole run rather than one sample. */}
      {entries.some((e) => e.result) && (
        <ParentalBalancePanel
          samples={balanceRowsFor(entries)}
          paired={pairedRowsFor(entries)}
        />
      )}
      {entries.some((e) => e.result) && (
        <ComparisonPanel
          regions={pooledRegions(entries, (sg) => segmentCoords(sg as never)).map((x) => x.region)}
          regionNames={pooledRegions(entries, (sg) => segmentCoords(sg as never)).map((x) => x.name)}
          markerPositions={pooledMarkers(entries)}
          loadTrack={loadFeatureTrack}
          existing={comparison}
          onDone={setComparison}
          sampleName={`${entries.filter((e) => e.result).length} samples in this run`}
          build={entries.find((e) => e.profile)?.profile?.build.build ?? 'assembly undetermined'}
        />
      )}

      {entries.some((e) => e.result) && (
        <LocusTest
          entries={entries} donor={donor} oocyte={oocyte} log={log}
          onResult={setLocus}
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
          {r.chroms.some((c) => c.aneuploidy) && (
            <Badge size="sm" variant="filled" color="orange">
              {r.chroms.filter((c) => c.aneuploidy).map((c) =>
                `chr${c.chrom} ${c.aneuploidy}`).join(', ')}
            </Badge>
          )}
          {r.segments.length > 0 && (
            <Badge size="sm" variant="filled" color="orange">
              {r.segments.length} chromosomal change{r.segments.length === 1 ? '' : 's'}
            </Badge>
          )}
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
          <StageCallout facts={r.stage && entry.profile
            ? stageFacts(r.stage, entry.profile) : null}
          />
          {r.gains.length > 0 && (
            <div style={{ marginTop: 10 }}>
              <Text size="xs" fw={700} mb={4}>Extra copies, and where they came from</Text>
              {r.gains.map((g) => (
                <div key={g.where} style={{ fontSize: 11, marginBottom: 5, lineHeight: 1.45 }}>
                  <span style={{ fontFamily: 'var(--om-mono)', fontWeight: 600 }}>{g.where}</span>
                  {'  '}
                  <span style={{ color: g.called ? 'var(--om-blue)' : 'var(--om-higher)' }}>
                    {g.origin}
                  </span>
                  <div style={{ color: 'var(--om-text-dim)' }}>{g.why}</div>
                </div>
              ))}
            </div>
          )}
          <SegmentCallout segments={r.segments} role={r.role} />
          <DefectCallout defects={defectsForResult(r)} />
          <GainCallout gains={r.gains} />
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
        // Two different measurements, so they sit together: how many heterozygotes SURVIVED, and
        // how far the survivors SCATTERED. The second is what decides the amplification class.
        ['Allele-fraction spread at het calls', Number.isFinite(p.hetBafSd)
          ? `${p.hetBafSd.toFixed(4)}   ${p.hetBafSd <= 0.11 ? 'unamplified' : 'amplified'}`
          : '-'],
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

/**
 * Chromosomal change, said once and prominently, at the top of the detail.
 *
 * A partly lost chromosome is the finding a reader is least likely to go looking for and most
 * likely to need: the whole-chromosome verdict for it is "unclear", which reads like an absence
 * of information rather than a located event. So this states it before anything else in the
 * detail, in the same visual language the tool uses for a warning, and the table below it carries
 * the numbers rather than replacing the sentence.
 */
/** What each segment kind is, in one line, because they are different events and a reader must
 *  not read a copy loss as a parental one or the reverse. */
const kindLabel = (k: SegmentKind, role: 'paternal' | 'maternal'): string => (
  k === 'copy-loss' ? 'DNA absent'
    : k === 'copy-gain' ? 'extra copies'
      : `${role} alleles absent`)

/**
 * Where each extra copy came from, or why that cannot be said.
 *
 * A gain with no origin is still shown. Hiding it would leave the reader believing the tool had
 * nothing to say about a finding it did make, and the reason for the refusal is the useful part:
 * an extra copy of the same homologue is bit-identical to a single copy and no channel this
 * array carries can see it.
 */
function GainCallout({ gains }: { gains: GainAnnotation[] }) {
  if (!gains.length) return null
  return (
    <div style={{ marginTop: 10, border: '1px solid var(--om-border)', padding: '10px 12px' }}>
      <Text size="xs" fw={700} mb={6}>
        Extra copies, and where they came from
      </Text>
      {gains.map((g) => (
        <div key={g.where} style={{ marginBottom: 8 }}>
          <div style={{ fontSize: 12 }}>
            <span style={{ fontFamily: 'var(--om-mono)', fontWeight: 600 }}>{g.where}</span>
            <span style={{ color: 'var(--om-text-dim)' }}>{' '}{g.kind}{' \u00b7 '}</span>
            <span style={{ color: g.called ? 'var(--om-blue)' : 'var(--om-higher)', fontWeight: 600 }}>
              {g.origin}
            </span>
          </div>
          <div style={{ fontSize: 11, color: 'var(--om-text-dim)', lineHeight: 1.45 }}>{g.why}</div>
        </div>
      ))}
      <Text size="xs" c="dimmed" style={{ lineHeight: 1.45 }}>
        Not validated on a true positive: no confirmed gain with a known origin exists in the
        material behind these thresholds, so this has been shown to refuse correctly and never to
        fire correctly. See the audit.
      </Text>
    </div>
  )
}

/**
 * Where the regions landed. Positional only, and the null is stated on the page rather than left
 * to the reader, because the number means nothing without it: a region is only callable where
 * markers are, so the comparison is against intervals carrying the same number of markers.
 */

function SegmentCallout({ segments, role }: {
  segments: Segment[], role: 'paternal' | 'maternal'
}) {
  if (!segments.length) return null
  const mb = (x: number): string => `${(x / 1e6).toFixed(1)} Mb`
  const total = segments.reduce((a, sg) => a + sg.spanBp, 0)
  const chroms = [...new Set(segments.map((sg) => `chr${sg.chrom}`))]
  return (
    <div style={{
      border: '1px solid var(--om-higher)', borderLeft: '4px solid var(--om-higher)',
      background: 'var(--om-warn-bg)', padding: '11px 14px', margin: '10px 0 4px',
    }}
    >
      <Text style={{ fontSize: 15, fontWeight: 700, color: 'var(--om-higher)', lineHeight: 1.25 }}>
        Chromosomal change on {chroms.join(', ')}
      </Text>
      <Text size="sm" mt={3} style={{ maxWidth: 760, lineHeight: 1.5 }}>
        {segments.length === 1 ? 'One region' : `${segments.length} regions`} totalling{' '}
        {mb(total)}, at a rate the rest of this genome does not reach. The whole-chromosome verdict
        cannot show this: a chromosome that is partly changed reads as neither present nor absent.
        Two different things are listed and the difference matters: <b>DNA absent</b> means the
        region is not there at all, read from the array no longer calling it with the intensity
        agreeing, while <b>{role} alleles absent</b> means the DNA is present but came from the
        other parent.
      </Text>
      <div style={{ marginTop: 7, display: 'grid', gap: 3 }}>
        {segments.map((sg) => (
          <Text key={`${sg.chrom}:${sg.startBp}`} size="xs" ff="monospace">
            chr{sg.chrom}&nbsp;{int(sg.startBp)}&ndash;{int(sg.endBp)}
            {'  '}&middot;{'  '}{mb(sg.spanBp)}
            {'  '}&middot;{'  '}{kindLabel(sg.kind, role)}
            {'  '}&middot;{'  '}{pct(sg.rate, 1)} against {pct(sg.nullRate, 1)}
          </Text>
        ))}
      </div>
      <Text size="xs" c="dimmed" mt={6} style={{ maxWidth: 760, lineHeight: 1.5 }}>
        A copy-number region is called only where the array stopped genotyping AND the intensity
        agrees by at least 1.0 log2. Measured over 46 arrays that requirement takes the scan from
        31 regions, including a 44.9 Mb one on a bulk diploid adult who cannot have it, to 17 with
        none on any of the six bulk arrays. Extra copies are found by the same machinery with the
        shift reversed, validated on a constructed positive class since no segmental gain occurs in
        those 46 arrays: a block of a whole-chromosome-gain array spliced into a clean one is
        called copy-gain at all four sizes tested, and the same splice from a euploid array returns
        nothing.
      </Text>
    </div>
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

      {r.segments.length > 0 && (
        <Section title={`Segments where the ${r.role} genome is missing`}>
          <Text size="xs" c="dimmed" mb={6} style={{ maxWidth: 780, lineHeight: 1.5 }}>
            A chromosome that is partly lost reads as neither present nor absent, so the whole
            chromosome comes back unclear and the missing part goes unreported. These are the
            regions inside a chromosome where the {r.role} allele is absent at a rate the rest of
            this genome does not reach. Each is scored against the median rate of the sample&rsquo;s
            OTHER chromosomes, which one event cannot move.
          </Text>
          <Table striped withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>Chromosome</Table.Th>
                <Table.Th ta="right">Span</Table.Th>
                <Table.Th ta="right">Markers</Table.Th>
                <Table.Th ta="right">Absent</Table.Th>
                <Table.Th ta="right">Against</Table.Th>
                <Table.Th ta="right">Score</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              {r.segments.map((sg) => (
                <Table.Tr key={`${sg.chrom}:${sg.startBp}`}>
                  <Table.Td>
                    chr{sg.chrom}
                    <Text span size="xs" c="dimmed" ff="monospace">
                      {' '}{int(segmentCoords(sg).start)}&ndash;{int(segmentCoords(sg).end)}
                    </Text>
                    <Text
                      span
                      size="xs"
                      ff="monospace"
                      c={segmentCoords(sg).localised ? 'dimmed' : 'var(--om-higher)'}
                    >
                      {' '}{segmentCoords(sg).interval}
                    </Text>
                  </Table.Td>
                  <Table.Td ta="right" ff="monospace">
                    {(segmentCoords(sg).spanBp / 1e6).toFixed(1)} Mb
                  </Table.Td>
                  <Table.Td ta="right" ff="monospace" c="dimmed">{int(sg.markers)}</Table.Td>
                  <Table.Td ta="right" ff="monospace">{pct(sg.rate, 2)}</Table.Td>
                  <Table.Td ta="right" ff="monospace" c="dimmed">{pct(sg.nullRate, 2)}</Table.Td>
                  <Table.Td ta="right" ff="monospace">{sg.score.toFixed(0)}</Table.Td>
                </Table.Tr>
              ))}
            </Table.Tbody>
          </Table>
          <Text size="xs" c="dimmed" mt={4} style={{ maxWidth: 780, lineHeight: 1.5 }}>
            The span is the resolution, not the event: the smallest region that can be called here
            is {int(MIN_SEGMENT_MARKERS)} informative markers, which on this array covers anywhere
            from a few hundred kilobases to tens of megabases depending on marker density. A real
            event smaller than that is reported at the size of the window that found it. The score
            is a log-likelihood ratio against a threshold of {SEGMENT_LRT}, and one just over it is
            a weak finding rather than a strong one.
          </Text>
        </Section>
      )}

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
              <Table.Td ta="right" ff="monospace" c="dimmed">
                {/* A ratio against the ceiling is a finding, and a chromosome whose calls are not
                    measuring it has none. Printing 18.4x beside "not measured" invites the
                    reader to take the number anyway. */}
                {c.verdict === 'not_measured' ? '\u2014' : (c.rate / ceiling).toFixed(1)}
              </Table.Td>
              <Table.Td>
                {c.verdict === 'present'
                  ? <Text span size="xs" c="dimmed">present</Text>
                  : (
                    <Badge size="xs" variant="light"
                      color={c.verdict === 'expected_absent' || c.verdict === 'not_measured'
                        ? 'gray' : 'orange'}
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
        // The absence indicator on every OTHER chromosome, in marker order, which is the null the
        // run-length p-value is read against. Measured away from the locus on purpose: a real
        // event is contiguous by definition and would inflate the null it is judged against.
        const elsewhere = new Map<string, { pos: number; absent: boolean }[]>()
        await profileFile(e.file, (r) => {
          const c = r.chrom.replace(/^chr/i, '').toUpperCase()
          if (c !== locus.chrom) {
            if (!isAutosome(c)) return
            const fa = donor.gt.get(r.probesetId)
            if (fa !== 'AA' && fa !== 'BB') return
            if (r.genotype === 'NC') return
            const arr = elsewhere.get(c) ?? []
            arr.push({ pos: r.pos, absent: r.genotype !== 'AB' && r.genotype !== fa })
            elsewhere.set(c, arr)
            return
          }
          markers.push({ rsid: r.probesetId, chrom: r.chrom, pos: r.pos })
          embryo.set(r.probesetId, r.genotype)
        }, () => {}, () => {})
        const clustering = measureClustering([...elsewhere.values()].map((xs) =>
          xs.sort((a, b) => a.pos - b.pos).map((x) => x.absent)))
        log('CALL', `${e.file.name}: artefact runs ${clustering.ratio.toFixed(1)}x independence `
          + `over ${int(clustering.n)} markers off chr${locus.chrom}`
          + `${clustering.independent ? '' : ', so no run-length p-value is reported'}`)
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
            clustering,
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
      {/*
        * BOTH parents, and the gate has to say so. This test needs a marker where one parent
        * could not have supplied the allele, so it needs two people, and the run below returns
        * immediately without both. The form used to be gated on the oocyte alone, which was
        * correct until oocyte-only runs became possible: after that a run with only an oocyte
        * rendered the form, enabled the button, and did nothing when it was pressed. A control
        * that is live in a state where it cannot act is worse than one that is absent.
        */}
      {!oocyte || !donor ? (
        <Alert color="orange" p="xs">
          <Text size="xs">
            This test needs <b>both</b> parents and is not run without them
            {!donor && !oocyte ? ', and neither is loaded'
              : !donor ? ', and the sperm donor is missing' : ', and the oocyte donor is missing'}.
            Presence cannot be established at any single marker without someone who could not have
            supplied the allele: the informative set would hold absences only, nothing could break
            a run, and the run-length statistic would have no null to be significant against. Label
            one file {!donor ? 'sperm' : 'oocyte'} and run again. Parent of origin, sperm type and
            segmental loss above need one parent only and are unaffected.
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

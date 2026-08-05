import { useRef, useState } from 'react'
import { Button, Group, Paper, Progress, Text } from '@mantine/core'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, headerMap, parseRow,
  type ChromStats, type ProbeRow, type SampleProfile,
} from './ingest'
import { ProductSet, groupByParent, MIN_PRODUCTS } from './inferredReference'
import { classify, emptyTally, tallyRow, pct } from './parentage'
import { ReferenceBlock, type Scored } from './ReferencePanel'
import { LimitsPanel, ProvenanceStamp, LIMITS } from './InferredLimits'
import { ExportBar } from './ExportBar'
import { buildInferredPdf } from './inferredPdf'
import {
  concordanceLongCsv, concordanceMatrixCsv, sampleResultsCsv, runManifestJson,
  type RunProvenance, type PairRate, type SampleResult,
} from './inferredExports'
import { SAME_PARENT_MAX, DIFFERENT_PARENT_MIN } from './inferredReference'
import { reportId, sha256 } from './syngamyPdf'
import { utc } from './fmt'
import { TriageTable, ConcordanceMatrix, triage, type Triage } from './ProductPanel'
import { Stage, Summary, Why } from './InferredShell'
import { int } from './fmt'
import { FeatureHeader, DropZone } from './FeatureHeader'
import { RunLog, type LogLine } from './RunLog'
import {
  PROGENITOR_EXAMPLES, PROGENITOR_CITATION, loadProductExample,
} from './progenitorExamples'

/**
 * Progenitor: reconstruct a parent's genotype from the haploid cells it produced.
 *
 * The case this exists for is an experiment with no array of the parent. Each haploid product
 * carries one allele per locus, so where the parent is homozygous every product agrees and the
 * allele is recoverable; where the parent is heterozygous the products disagree, which is both
 * how heterozygosity is detected and the source of the only error the method has.
 *
 * Files are read in this browser and there is no endpoint to send them to. A product set of
 * eight arrays costs about 6.6 MB here, because each product is a Uint8Array of allele codes
 * over a shared probe index rather than a map per file.
 */

type Stage0 = 'idle' | 'reading' | 'ready' | 'failed'

type Tag = 'READ' | 'GATE' | 'PAIR' | 'BUILD' | 'WARN' | 'DONE'
const TAG_COLOUR: Record<Tag, string> = {
  READ: 'var(--om-text-dim)',
  GATE: 'var(--om-text-dim)',
  PAIR: 'var(--om-blue)',
  BUILD: 'var(--om-blue)',
  WARN: 'var(--om-higher)',
  DONE: 'var(--om-blue)',
}

/** Concordance and grouping run on demand rather than on every render, so the log has a run to
 *  narrate and a large product set is not recomputed while files are still arriving. */
interface Analysis {
  groups: number[][]
  usableIds: string[]
}

interface Loaded {
  id: string
  name: string
  profile: SampleProfile
  slot: number
  /** Kept so a product can be re-streamed to score it against a reference built without it.
   *  The allele codes in the ProductSet are enough to BUILD, not to score: scoring needs the
   *  B-allele frequencies too, which is what the noise ceiling is computed from. */
  file: File
  sha: string
}

const mb = (b: number): string => `${(b / 1e6).toFixed(1)} MB`

async function eachLine(
  file: File, fn: (line: string) => void, onChunk?: (bytes: number) => void,
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
    // Yield, so the progress bar paints while a 40 MB file streams.
    await new Promise((r) => { setTimeout(r, 0) })
  }
  if (carry) fn(carry)
}

export function ProgenitorPage() {
  const [products, setProducts] = useState<Loaded[]>([])
  const [set, setSet] = useState<ProductSet | null>(null)
  const [state, setState] = useState<Stage0>('idle')
  const [busyName, setBusyName] = useState('')
  const [pctDone, setPctDone] = useState(0)
  const [errors, setErrors] = useState<string[]>([])
  const pick = useRef<HTMLInputElement>(null)
  const [building, setBuilding] = useState(false)
  const [lines, setLines] = useState<LogLine<Tag>[]>([])
  const [analysis, setAnalysis] = useState<Analysis | null>(null)
  const [examples, setExamples] = useState(false)
  const log = (tag: Tag, text: string): void =>
    setLines((p) => [...p.slice(-499), { tag, text }])
  const [ref, setRef] = useState<{
    group: number
    stats: ReturnType<ProductSet['build']>
    ratios: Map<number, number>
    members: Scored[]
  } | null>(null)

  const add = async (files: FileList | File[]): Promise<void> => {
    const ps = set ?? new ProductSet()
    const next = [...products]
    const errs: string[] = []
    setState('reading')
    for (const f of [...files]) {
      const id = `${f.name}:${f.size}`
      if (next.some((p) => p.id === id)) continue
      setBusyName(`${f.name}  ${mb(f.size)}`)
      setPctDone(0)
      log('READ', `${f.name}  ${mb(f.size)}`)
      const name = f.name.replace(/\.(csv|txt|probes)(\.gz)?$/i, '')
      try {
        const slot = ps.begin(name)
        const baf = { inBand: 0, total: 0 }
        const byChrom = new Map<string, ChromStats>()
        const bafSums = emptyBafSums()
        const builds = emptyBuildSums()
        let map: ReturnType<typeof headerMap> = null
        let firstId = ''
        let n = 0
        await eachLine(f, (line) => {
          if (!map) { map = headerMap(line); return }
          const r: ProbeRow | null = parseRow(line, map)
          if (!r) return
          if (!firstId) firstId = r.probesetId
          accumulate(r, byChrom)
          accumulateBaf(r, bafSums)
          accumulateBuild(r, builds)
          ps.add(slot, r, baf)
          n += 1
        }, (bytes) => setPctDone(Math.round((bytes / f.size) * 100)))
        if (!map) throw new Error('no recognisable header row')
        ps.end(slot, baf)
        const profile: SampleProfile = {
          ...finishProfile(f.name, byChrom, bafSums, firstId, builds),
          build: buildVerdict(builds),
        }
        const t = triage(id, name, profile)
        log(t.usable ? 'GATE' : 'WARN',
          `${t.name}: call ${pct(profile.callRate, 1)}, band ${pct(profile.hetBand, 2)}`
          + `${t.usable ? ', usable product' : `, set aside: ${t.why}`}`)
        next.push({
          id, name: t.name, profile, slot, file: f, sha: await sha256(f),
        })
      } catch (err) {
        const m = err instanceof Error ? err.message : String(err)
        errs.push(`${f.name}: ${m}`)
        log('WARN', `${f.name}: ${m}`)
      }
    }
    setSet(ps)
    setProducts(next)
    setErrors(errs)
    setBusyName('')
    setState(next.length ? 'ready' : 'failed')
  }

  const EXPORT_NOTE = 'The reconstructed genotype itself is not exported. The build is '
    + 'deterministic, so the product files and their checksums are the reference and re-running '
    + 'reproduces it exactly. A file of homozygous calls belonging to an identifiable person, '
    + 'which could be re-imported as though it were a measured array, is a hazard with no '
    + 'matching gain.'

  const rows: Triage[] = products.map((p) => triage(p.id, p.name, p.profile))
  const usable = rows.filter((r) => r.usable)
  const usableSlots = usable.map((r) => products.find((p) => p.id === r.id)!.slot)
  const rate = (a: number, b: number): number =>
    (set ? set.opposite(usableSlots[a], usableSlots[b]).rate : NaN)
  const groups = analysis?.groups ?? []
  const names = usable.map((r) => r.name)
  const biggest = groups[0]?.length ?? 0

  /** Compare every pair and split into groups that share a parent. Nothing is built here: a
   *  reference is only reconstructed once a group has been chosen. */
  const run = (): void => {
    if (!set || usable.length < 2) return
    setRef(null)
    log('PAIR', `comparing ${usable.length * (usable.length - 1) / 2} pairs of products`)
    let lo = 1
    let hi = 0
    for (let a = 0; a < usable.length; a += 1) {
      for (let b = a + 1; b < usable.length; b += 1) {
        const r = rate(a, b)
        lo = Math.min(lo, r)
        hi = Math.max(hi, r)
      }
    }
    log('PAIR', `opposite-homozygote rate ${pct(lo, 2)} to ${pct(hi, 2)}`)
    const g = groupByParent(usable.length, rate)
    g.forEach((grp, i) => {
      log(g.length > 1 ? 'WARN' : 'DONE',
        `group ${i + 1}: ${grp.map((x) => usable[x].name).join(', ')}`)
    })
    if (g.length > 1) {
      log('WARN', `${g.length} parent groups. A reference built from all of these would be `
        + 'built from more than one person.')
    }
    setAnalysis({ groups: g, usableIds: usable.map((u) => u.id) })
  }

  const loadExamples = async (): Promise<void> => {
    setExamples(true)
    log('READ', `fetching ${PROGENITOR_EXAMPLES.length} public example products`)
    const files: File[] = []
    for (const ex of PROGENITOR_EXAMPLES) {
      try {
        files.push(await loadProductExample(ex))
      } catch (err) {
        log('WARN', `${ex.file}: ${err instanceof Error ? err.message : String(err)}`)
      }
    }
    if (files.length) await add(files)
  }

  // Membership is computed from the products that passed the gates, so adding or clearing a file
  // invalidates it rather than leaving a stale group on screen.
  const stale = analysis !== null
    && (analysis.usableIds.length !== usable.length
      || analysis.usableIds.some((id, i) => usable[i]?.id !== id))

  /**
   * Build the reference for a group and verify every member against a reference built without
   * it. A product scored against a reference containing itself reads exactly zero absence, a
   * bias the size of the whole signal, so leave-one-out is the only honest check.
   */
  const buildFor = async (gi: number): Promise<void> => {
    if (!set) return
    setBuilding(true)
    try {
      const chosen = groups[gi].map((i) => usable[i])
      const keep = new Set(chosen.map((c) => products.find((p) => p.id === c.id)!.name))
      const drop = set.ids.filter((id) => !keep.has(id))
      log('BUILD', `${chosen.length} products in the group, ${drop.length} of `
        + `${set.ids.length} loaded excluded`)
      const { mMin, ratios } = set.chooseM(drop)
      const full = set.build(mMin, drop)
      log('BUILD', `m >= ${mMin}, ${int(full.markers)} markers, `
        + `ascertainment ${pct(ratios.get(mMin) ?? NaN, 1)}, `
        + `contamination ${pct(full.contamination, 2)}`)
      const scored: Scored[] = []
      for (const c of chosen) {
        const load = products.find((p) => p.id === c.id)!
        // n-1 products can never reach a threshold of n, so the leave-one-out build has to sit
        // one below the chosen depth when the chosen depth is the whole group.
        const ref = set.build(Math.min(mMin, chosen.length - 1), [...drop, load.name])
        const t = emptyTally()
        let map: ReturnType<typeof headerMap> = null
        await eachLine(load.file, (line) => {
          if (!map) { map = headerMap(line); return }
          const r = parseRow(line, map)
          if (r) tallyRow(ref.genotype.get(r.probesetId) ?? 'NC', r, t)
        })
        // The leave-one-out reference is built from one fewer product, so it is dirtier than the
        // full one and carries its own floor rather than the full build's.
        const res = classify(t, ref.hRetained, { spuriousAbsence: ref.spuriousAbsence })
        log(res.verdict === 'parent_genome_present' ? 'DONE' : 'WARN',
          `${load.name}: ${pct(res.genomeRate, 2)} absent of ${pct(res.explainable, 2)} `
          + `= ${(res.genomeRate / res.explainable).toFixed(2)}x, ${res.verdict.replace(/_/g, ' ')}`)
        scored.push({
          name: load.name,
          absence: res.genomeRate,
          ceiling: res.explainable,
          ratio: res.genomeRate / res.explainable,
          verdict: res.verdict,
        })
      }
      setRef({ group: gi, stats: full, ratios, members: scored })
    } finally {
      setBuilding(false)
    }
  }

  const verdict = groups.length > 1 ? `${groups.length} parent groups, not one`
      : biggest >= MIN_PRODUCTS ? 'One parent, enough products to reconstruct'
        : `One parent, but only ${biggest} usable product${biggest === 1 ? '' : 's'}`

  const sub = groups.length > 1
      ? 'These arrays do not all come from the same parent. A single reference would be built '
        + 'from two people, so a group has to be chosen before anything is reconstructed.'
      : biggest >= MIN_PRODUCTS
        ? 'Every pair agrees at the rate two products of one parent do.'
        : `Below ${MIN_PRODUCTS} products the method inverts: true offspring read as decisively `
          + 'absent rather than uncalled. Nothing is built.'

  const chosenNames = ref ? groups[ref.group].map((i) => usable[i].name) : []
  const provenance: RunProvenance | null = ref ? {
    tool: 'OriginMarker',
    generatedAt: utc(new Date().toISOString()),
    reportId: reportId(products.map((p) => p.sha).join('') + chosenNames.join('')),
    experiment: products[0]?.profile.idPrefix ?? '',
    products: chosenNames,
    mMin: ref.stats.mMin,
    markers: ref.stats.markers,
    meanM: ref.stats.meanM,
    ascertainment: ref.ratios.get(ref.stats.mMin) ?? NaN,
    contamination: ref.stats.contamination,
    spuriousAbsence: ref.stats.spuriousAbsence,
    sameParentMax: SAME_PARENT_MAX,
    differentParentMin: DIFFERENT_PARENT_MIN,
  } : null

  const pairs: PairRate[] = usable.flatMap((a, i) => usable.slice(i + 1).map((b, j) => {
    const gi = groups.findIndex((g) => g.includes(i))
    const gj = groups.findIndex((g) => g.includes(i + 1 + j))
    return { a: a.name, b: b.name, rate: rate(i, i + 1 + j), sameGroup: gi === gj && gi >= 0 }
  }))

  const sampleRows: SampleResult[] = rows.map((r) => {
    const gi = groups.findIndex((g) => g.some((x) => usable[x].id === r.id))
    const scored = ref?.members.find((m) => m.name === r.name)
    return {
      sample: r.name,
      role: chosenNames.includes(r.name) ? 'reference product'
        : r.usable ? 'product, other group' : 'excluded',
      group: gi >= 0 ? `group ${gi + 1}` : '',
      callRate: r.profile.callRate,
      hetBand: r.profile.hetBand,
      absence: scored?.absence,
      ceiling: scored?.ceiling,
      ratio: scored?.ratio,
      verdict: scored?.verdict,
      excludedBecause: r.usable ? '' : r.why,
    }
  })

  const exportItems = provenance ? [
    {
      label: 'The report',
      filename: `progenitor-report-${provenance.reportId}.pdf`,
      mime: 'application/pdf',
      hint: 'Letter PDF in the same format as the Syngamy report. Every page states that the '
        + 'reference was inferred, so a page lifted out of context still says so.',
      build: () => buildInferredPdf({
        provenance,
        groups: groups.map((g) => g.map((i) => usable[i].name)),
        samples: sampleRows,
        pairs,
        matrixNames: names,
        matrixRate: rate,
        ascertainment: Object.fromEntries(ref!.ratios),
        members: ref!.members.map((m) => ({
          name: m.name, absence: m.absence, ceiling: m.ceiling ?? NaN, ratio: m.ratio,
          verdict: m.verdict,
        })),
        controls: [],
        limits: LIMITS,
        files: products.map((pr) => ({
          name: pr.file.name,
          role: chosenNames.includes(pr.name) ? 'product' : 'set aside',
          size: pr.file.size,
          markers: pr.profile.markers,
          sha256: pr.sha,
        })),
        fromExamples: examples,
        markSvg: 'favicon.svg',
      }),
    },
    {
      label: 'Pairwise concordance',
      filename: `progenitor-concordance-${provenance.reportId}.csv`,
      mime: 'text/csv',
      hint: 'One row per pair of products, with the rate and whether they grouped together. '
        + 'The shape a plotting or clustering script wants.',
      build: () => concordanceLongCsv(pairs, provenance),
    },
    {
      label: 'Concordance matrix',
      filename: `progenitor-matrix-${provenance.reportId}.csv`,
      mime: 'text/csv',
      hint: 'The square matrix as printed, ordered by group. Loads with index_col=0 for a '
        + 'heatmap or a supplementary table.',
      build: () => concordanceMatrixCsv(names, rate, provenance),
    },
    {
      label: 'Per-sample results',
      filename: `progenitor-samples-${provenance.reportId}.csv`,
      mime: 'text/csv',
      hint: 'One row per array including the ones excluded, each with its reason, so the table '
        + 'accounts for every file submitted.',
      build: () => sampleResultsCsv(sampleRows, provenance),
    },
    {
      label: 'Run manifest',
      filename: `progenitor-run-${provenance.reportId}.json`,
      mime: 'application/json',
      hint: 'The whole run as JSON: reference parameters, membership, every sample, every pair, '
        + 'and what was withheld. For a pipeline or a LIMS.',
      build: () => runManifestJson(provenance, groups.map((g) => g.map((i) => usable[i].name)),
        sampleRows, pairs, Object.fromEntries(ref!.ratios), LIMITS),
    },
  ] : []

  return (
    <div style={{ paddingBottom: 60 }}>
      <FeatureHeader name="Progenitor" tagline="reconstructed parental genotype" />

      {lines.length > 0 && (
        <RunLog
          lines={lines}
          colours={TAG_COLOUR}
          onDownload={() => {
            const head = [
              'Progenitor run log',
              `generated ${new Date().toISOString()}`,
              `url ${window.location.href}`,
              examples ? `examples: ${PROGENITOR_CITATION}` : '',
              '',
            ].filter(Boolean).join('\n')
            const body = lines.map((l) => `[${l.tag}] ${l.text}`).join('\n')
            const blob = new Blob([`${head}${body}\n`], { type: 'text/plain;charset=utf-8' })
            const url = URL.createObjectURL(blob)
            const a = document.createElement('a')
            a.href = url
            a.download =
              `progenitor-runlog-${new Date().toISOString().slice(0, 19)
                .replace(/[:T]/g, '')}.txt`
            document.body.appendChild(a)
            a.click()
            a.remove()
            URL.revokeObjectURL(url)
          }}
        />
      )}

      {analysis && !stale && (
      <Summary
        verdict={verdict}
        sub={sub}
        stats={[
          { label: 'arrays in', value: String(products.length) },
          { label: 'usable products', value: String(usable.length) },
          { label: 'parent groups', value: String(groups.length) },
          { label: 'largest group', value: String(biggest) },
        ]}
      />
      )}

      <Paper p="sm" mt={10}>
        <Group justify="flex-end" align="center" mb={8}>
          <Group gap={6}>
            <input
              ref={pick}
              type="file"
              multiple
              accept=".csv,.txt,.probes,.gz,text/csv,text/plain"
              style={{ display: 'none' }}
              onChange={(e) => { if (e.target.files?.length) void add(e.target.files) }}
            />
            <Button
              variant="default" size="xs" disabled={state === 'reading'}
              onClick={() => pick.current?.click()}
            >
              Add files&hellip;
            </Button>
            {products.length === 0 && (
              <Button
                variant="default" size="xs" disabled={state === 'reading'}
                onClick={() => { void loadExamples() }}
              >
                Examples
              </Button>
            )}
            <Button
              size="xs"
              disabled={state === 'reading' || usable.length < 2 || (!!analysis && !stale)}
              onClick={run}
            >
              Run
            </Button>
            {products.length > 0 && (
              <Button
                variant="subtle" size="xs" disabled={state === 'reading'}
                onClick={() => {
                  setProducts([]); setSet(null); setErrors([]); setState('idle')
                  setRef(null); setAnalysis(null); setLines([]); setExamples(false)
                }}
              >
                Clear
              </Button>
            )}
          </Group>
        </Group>

        <DropZone
          empty={products.length === 0}
          prompt={'Drop the haploid products of one parent. One file per cell: a pronucleus, '
            + 'a polar body, a single sperm. Read in this browser; nothing is uploaded.'}
          disabled={state === 'reading'}
          onFiles={(f) => { void add(f) }}
        >
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 6 }}>
            {rows.map((r) => (
              <div
                key={r.id}
                style={{
                  border: `1px solid ${r.usable ? 'var(--om-blue)' : 'var(--om-border)'}`,
                  borderRadius: 2,
                  padding: '3px 8px',
                  background: r.usable ? 'transparent' : 'var(--om-zebra)',
                  fontSize: 11,
                }}
              >
                <div style={{
                  fontFamily: 'var(--om-mono)', fontWeight: 600,
                  color: r.usable ? undefined : 'var(--om-text-dim)',
                }}
                >
                  {r.name}
                </div>
                <div style={{
                  fontSize: 10, color: r.usable ? 'var(--om-blue)' : 'var(--om-text-dim)',
                }}
                >
                  {r.usable ? `${pct(r.profile.callRate, 0)} call, band ${pct(r.profile.hetBand, 1)}`
                    : 'set aside'}
                </div>
              </div>
            ))}
          </div>
        </DropZone>

        {state === 'reading' && (
          <>
            <Text size="xs" c="dimmed" mt={6}>Reading {busyName}</Text>
            <Progress value={pctDone} size="sm" mt={4} radius={0} />
          </>
        )}
        {errors.map((e) => (
          <Text key={e} size="xs" mt={4} style={{ color: 'var(--om-higher)' }}>{e}</Text>
        ))}
      </Paper>

      {products.length > 0 && (
        <Stage
          n={1}
          title="Quality and ploidy"
          state={usable.length === products.length ? 'ok' : 'attention'}
          badges={[`${products.length - usable.length} set aside`]}
          headline={`${usable.length} of ${products.length} arrays can be a product`}
          defaultOpen={usable.length !== products.length}
        >
          <TriageTable rows={rows} />
          <Why>
            A product is one haploid meiotic cell. An array below the call-rate floor cannot be
            judged at all, and a genome that is not haploid is not one meiotic draw whatever its
            alleles say. Both gates run before anything is measured against anything else:
            reversing them lets a failed array read as a separate parent, and lets a diploid
            bridge two unrelated ones.
          </Why>
        </Stage>
      )}

      {analysis && !stale && usable.length > 1 && (
        <Stage
          n={2}
          title="Membership"
          state={groups.length > 1 ? 'attention' : 'ok'}
          badges={groups.map((g, i) => `group ${i + 1}: ${g.length}`)}
          headline={groups.length > 1
            ? `These ${usable.length} products are not all from one parent`
            : `All ${usable.length} products agree on one parent`}
          defaultOpen
        >
          <ConcordanceMatrix names={names} rate={rate} groups={groups} />
          <Why>
            Two haploid products of one parent differ only where that parent is heterozygous and
            the two drew differently. Measured across two experiments: 4.68% to 9.70% within one
            parent over 46 pairs, and 9.88% to 16.10% between parents over 45. Every pair inside
            a group has to read same parent, because one genuine cross-parent pair measured
            9.88%, under the cut, and grouping by chained links merges two parents through that
            single edge.
          </Why>
        </Stage>
      )}

      {analysis && !stale && groups.length > 0 && (
        <Stage
          n={3}
          title="Reconstruction"
          state={ref ? 'ok' : biggest >= MIN_PRODUCTS ? 'pending' : 'blocked'}
          badges={ref ? [`m ≥ ${ref.stats.mMin}`,
            `${pct(ref.stats.contamination, 2)} contamination`] : []}
          headline={ref
            ? `Built from ${ref.stats.nProducts} products, and `
              + `${ref.members.filter((m) => m.verdict === 'parent_genome_present').length} of `
              + `${ref.members.length} verify against a reference built without them`
            : biggest >= MIN_PRODUCTS
              ? 'Choose which parent to reconstruct'
              : `No group has the ${MIN_PRODUCTS} products this method needs`}
          defaultOpen
        >
          {biggest >= MIN_PRODUCTS && !ref && (
            <Group gap={8} mt={8}>
              {groups.map((g, i) => (
                <Button
                  key={i}
                  size="sm"
                  radius={2}
                  variant={g.length >= MIN_PRODUCTS ? 'filled' : 'default'}
                  disabled={g.length < MIN_PRODUCTS || building}
                  loading={building}
                  onClick={() => void buildFor(i)}
                >
                  Reconstruct group {i + 1} ({g.length} products)
                </Button>
              ))}
            </Group>
          )}
          {ref && (
            <ReferenceBlock
              group={groups[ref.group].map((i) => usable[i].name)}
              mMin={ref.stats.mMin}
              markers={ref.stats.markers}
              meanM={ref.stats.meanM}
              hRetained={ref.stats.hRetained}
              contamination={ref.stats.contamination}
              spuriousAbsence={ref.stats.spuriousAbsence}
              ratios={Object.fromEntries(ref.ratios)}
              members={ref.members}
              controls={[]}
            />
          )}
        </Stage>
      )}

      {ref && (
        <>
          <Stage
            n={4}
            title="Limits"
            state="blocked"
            badges={[`${LIMITS.length} withheld`]}
            headline={`${LIMITS.length} things this reference cannot tell you`}
          >
            <LimitsPanel />
          </Stage>
          <ProvenanceStamp p={provenance!} />
          <ExportBar items={exportItems} note={EXPORT_NOTE} />
        </>
      )}

      {products.length === 0 && (
        <Text size="xs" c="dimmed" mt={10} style={{ maxWidth: 760, lineHeight: 1.55 }}>
          Takes the same exports the rest of the tool does: Axiom <code>.probes</code>, or a
          genotype table as CSV or TSV, gzipped or not. {int(825657)} markers streams in a few
          seconds per file, and eight products cost about 6.6 MB, because each is held as one byte
          per marker rather than as a map. At least {MIN_PRODUCTS} products of one parent are
          needed before anything is reconstructed.
        </Text>
      )}
    </div>
  )
}

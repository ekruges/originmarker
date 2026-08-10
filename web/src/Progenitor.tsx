import { useRef, useState } from 'react'
import { Button, Group, Paper, Progress, Table, Text } from '@mantine/core'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, headerMap, parseRow,
  type ChromStats, type ProbeRow, type SampleProfile,
} from './ingest'
import { ProductSet, groupByParent, MIN_PRODUCTS } from './inferredReference'
import {
  emptySex, accumulateSex, sexCall, paternalGroup, type SexCall,
} from './sexing'
import { inferredArrayText, isInferredFile } from './inferredArray'
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

/** One array's parental-origin call, which is what the whole feature exists to produce. */
interface Origin {
  name: string
  /** 'paternal' | 'maternal' when the father's group was named, 'this parent' | 'the other
   *  parent' when it was not, 'unclear' when the array's own noise can manufacture its absence. */
  origin: string
  absence: number
  ceiling: number
  ratio: number
  usable: boolean
  inReference: boolean
}

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
  /** Index into `groups` of the father's group, or null when the products did not say. */
  paternal: number | null
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
  /** Whether this product carries a Y, which is the only thing that names the reconstructed
   *  parent when there is no array of that parent to compare against. */
  sex: SexCall
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

/** The answer, one row per array. Everything else in this page is how it was arrived at. */
function OriginTable({ origins, side }: { origins: Origin[]; side: string }) {
  const colour = (o: string): string =>
    (o === 'unclear' ? 'var(--om-text-dim)'
      : o === 'paternal' || o === 'this parent' ? 'var(--om-blue)' : 'var(--om-higher)')
  return (
    <Table striped withTableBorder fz={11} mt={8}>
      <Table.Thead>
        <Table.Tr>
          <Table.Th>Array</Table.Th>
          <Table.Th>Parental origin</Table.Th>
          <Table.Th ta="right">Absence</Table.Th>
          <Table.Th ta="right">Ceiling</Table.Th>
          <Table.Th ta="right">x ceiling</Table.Th>
          <Table.Th>Note</Table.Th>
        </Table.Tr>
      </Table.Thead>
      <Table.Tbody>
        {origins.map((o) => (
          <Table.Tr key={o.name}>
            <Table.Td style={{ fontFamily: 'var(--om-mono)' }}>{o.name}</Table.Td>
            <Table.Td style={{ color: colour(o.origin), fontWeight: 600 }}>{o.origin}</Table.Td>
            <Table.Td ta="right">{pct(o.absence, 2)}</Table.Td>
            <Table.Td ta="right">{pct(o.ceiling, 2)}</Table.Td>
            <Table.Td ta="right">{o.ratio.toFixed(2)}x</Table.Td>
            <Table.Td style={{ color: 'var(--om-text-dim)' }}>
              {!o.usable ? 'excluded as a product, still called'
                : o.inReference ? `in the ${side || 'reference'} set, scored without itself` : ''}
            </Table.Td>
          </Table.Tr>
        ))}
      </Table.Tbody>
    </Table>
  )
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
    origins: Origin[]
    side: string
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
      // Refused before it is read, on the mark rather than on its numbers. The ploidy gate does
      // reject the file as written, but only incidentally: the export carries no BAF column, so
      // the band is undefined and the gate lands on "borderline". That is a column happening to
      // be absent, not a reconstruction being recognised. On the genotypes themselves it is
      // homozygous at every marker, which is 0% heterozygosity and the cleanest haploid product
      // ever submitted. Add a BAF column, or convert the file, and it would pass. Building from
      // it would fold a reference into itself, agree with itself everywhere, and report a
      // contamination of nothing.
      if (await isInferredFile(f)) {
        const why = 'this is a reconstructed genotype, not a haploid cell. It cannot be a '
          + 'product: it is homozygous everywhere by construction, so it would agree with a '
          + 'reference built from it at every marker and report a contamination of nothing.'
        errs.push(`${f.name}: ${why}`)
        log('WARN', `${name}: ${why}`)
        continue
      }
      try {
        const slot = ps.begin(name)
        const baf = { inBand: 0, total: 0 }
        const sx = emptySex()
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
          accumulateSex(r, sx)
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
        const sex = sexCall(sx)
        if (sex.yBearing) {
          log('GATE', `${t.name}: carries a Y, so its group is the father's`)
        }
        next.push({
          id, name: t.name, profile, slot, file: f, sha: await sha256(f), sex,
        })
        // Published per file rather than per batch. A folder of twenty arrays is a minute of
        // reading, and showing nothing until the last one lands makes a working page look hung.
        setProducts([...next])
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

  const EXPORT_NOTE = 'The reconstructed genotype is exported as an array file, so it can be '
    + 'used as the reference in Syngamy or anywhere else an array goes. It opens with a banner '
    + 'saying it is not a measured array and carries a machine-readable mark to the same effect, '
    + 'because a file of homozygous calls that reads like a real person\'s array is the one '
    + 'hazard this feature has.'

  const rows: Triage[] = products.map((p) => triage(p.id, p.name, p.profile))
  const usable = rows.filter((r) => r.usable)
  const usableSlots = usable.map((r) => products.find((p) => p.id === r.id)!.slot)
  const rate = (a: number, b: number): number =>
    (set ? set.opposite(usableSlots[a], usableSlots[b]).rate : NaN)
  const sexOf = (t: Triage): SexCall => products.find((p) => p.id === t.id)!.sex
  const groups = analysis?.groups ?? []
  const names = usable.map((r) => r.name)
  const biggest = groups[0]?.length ?? 0

  /**
   * Split the products into groups sharing a parent, name the father's group, and reconstruct it.
   *
   * All of it, in one pass, because that is what the tool is for. An experiment arrives as a
   * folder of arrays with nothing saying which pronucleus came from where; every step between
   * that and a parental-origin call is mechanical, so none of them is a question to put to the
   * person holding the folder.
   */
  const run = async (): Promise<void> => {
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
    const pat = paternalGroup(g, usable.map((u) => sexOf(u).yBearing))
    g.forEach((grp, i) => {
      log(g.length > 1 ? 'WARN' : 'DONE',
        `group ${i + 1}: ${grp.map((x) => usable[x].name).join(', ')}`
        + (i === pat ? '  <- the father, by the Y its products carry' : ''))
    })
    if (g.length > 1) {
      log('WARN', `${g.length} parent groups. A reference built from all of these would be `
        + 'built from more than one person.')
    }
    if (pat === null) {
      // Three different reasons, and they are not interchangeable. A file with no Y probes at
      // all was never asked the question; a set that answered "no" everywhere was.
      const noProbes = usable.every((u) => sexOf(u).yBearing === null)
      const anyY = usable.some((u) => sexOf(u).yBearing)
      log('WARN', noProbes
        ? 'none of these files carries chromosome Y probes, so nothing here can say which group '
          + 'is the father\'s. The split still holds; only the naming is withheld.'
        : anyY
          ? 'more than one group carries a Y, which one sperm donor cannot do. Which side is '
            + 'which is left unnamed.'
          : 'no product carries a Y, so which group is the father cannot be established. A '
            + 'paternal group of n products is all X-bearing 2^-n of the time. The split still '
            + 'holds; only the naming is withheld.')
    }
    setAnalysis({ groups: g, usableIds: usable.map((u) => u.id), paternal: pat })
    // Reconstruct the largest group that clears the floor. Whether it is the father's or the
    // mother's does not change the call: an array either carries that parent's genome or it does
    // not, and the other side follows.
    const bi = g.findIndex((x) => x.length >= MIN_PRODUCTS)
    if (bi < 0) {
      log('WARN', `no group reaches the ${MIN_PRODUCTS} products this method needs; `
        + 'nothing is reconstructed')
      return
    }
    await buildFor(bi, g, pat)
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
  const buildFor = async (
    gi: number, gs: number[][] = groups, pat: number | null = analysis?.paternal ?? null,
  ): Promise<void> => {
    if (!set) return
    setBuilding(true)
    try {
      const chosen = gs[gi].map((i) => usable[i])
      const keep = new Set(chosen.map((c) => products.find((p) => p.id === c.id)!.name))
      const drop = set.ids.filter((id) => !keep.has(id))
      log('BUILD', `${chosen.length} products in the group, ${drop.length} of `
        + `${set.ids.length} loaded excluded`)
      const { mMin, ratios } = set.chooseM(drop)
      const full = set.build(mMin, drop)
      const side = pat === null ? '' : gi === pat ? 'paternal' : 'maternal'
      const other = side === 'paternal' ? 'maternal' : 'paternal'
      log('BUILD', `m >= ${mMin}, ${int(full.markers)} markers, `
        + `ascertainment ${pct(ratios.get(mMin) ?? NaN, 1)}, `
        + `contamination ${pct(full.contamination, 2)}`)
      const scored: Scored[] = []
      const called: Origin[] = []
      const inGroup = new Set(chosen.map((c) => c.id))
      // EVERY array is scored, not only the ones the reference was built from. The arrays that
      // did not go into it are the ones the answer is wanted for: an array that carries this
      // parent's genome came from this parent, and one that decisively does not came from the
      // other. Excluded arrays are scored too, and say on their own row that they were excluded,
      // because a fused or half-failed zygote is the case this is most often asked about and
      // silently dropping it answers nothing.
      for (const load of products) {
        const mine = inGroup.has(load.id)
        // A product scored against a reference containing itself reads exactly zero absence, a
        // bias the size of the whole signal, so a member is scored leave-one-out. n-1 products
        // cannot reach a threshold of n, so that build sits one below the chosen depth.
        const ref = mine
          ? set.build(Math.min(mMin, chosen.length - 1), [...drop, load.name]) : full
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
        const ratio = res.genomeRate / res.explainable
        const origin = res.verdict === 'unclear' ? 'unclear'
          : side === '' ? (res.verdict === 'parent_genome_present'
            ? 'this parent' : 'the other parent')
            : res.verdict === 'parent_genome_present' ? side : other
        log(res.verdict === 'parent_genome_present' ? 'DONE' : 'WARN',
          `${load.name}: ${pct(res.genomeRate, 2)} absent of ${pct(res.explainable, 2)} `
          + `= ${ratio.toFixed(2)}x, ${origin}`)
        called.push({
          name: load.name,
          origin,
          absence: res.genomeRate,
          ceiling: res.explainable,
          ratio,
          usable: rows.find((r) => r.id === load.id)?.usable ?? false,
          inReference: mine,
        })
        if (mine) {
          scored.push({
            name: load.name,
            absence: res.genomeRate,
            ceiling: res.explainable,
            ratio,
            verdict: res.verdict,
          })
        }
      }
      const decided = called.filter((c) => c.origin !== 'unclear').length
      log('DONE', `${decided} of ${called.length} arrays called`
        + (side ? `, against a ${side} reference` : ', as one side or the other'))
      setRef({ group: gi, stats: full, ratios, members: scored, origins: called, side })
    } finally {
      setBuilding(false)
    }
  }

  const decided = ref ? ref.origins.filter((o) => o.origin !== 'unclear').length : 0
  const verdict = ref
    ? `${decided} of ${ref.origins.length} arrays called`
    : groups.length > 1 ? `${groups.length} parent groups, not one`
      : biggest >= MIN_PRODUCTS ? 'One parent, enough products to reconstruct'
        : `One parent, but only ${biggest} usable product${biggest === 1 ? '' : 's'}`

  const sub = ref
    ? (ref.side
      ? `Parental origin, against a ${ref.side} genotype reconstructed from `
        + `${ref.stats.nProducts} of these arrays. No array of either parent was used.`
      : 'The arrays split into groups by shared parent, but nothing in them says which group is '
        + 'the father\'s: naming that needs a product carrying a Y, which a maternal cell cannot. '
        + 'The split holds; the two sides are named "this parent" and "the other parent" rather '
        + 'than guessed at.')
    : groups.length > 1
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
    const called = ref?.origins.find((o) => o.name === r.name)
    return {
      sample: r.name,
      origin: called?.origin ?? '',
      role: chosenNames.includes(r.name) ? 'reference product'
        : r.usable ? 'product, other group' : 'excluded',
      group: gi >= 0 ? `group ${gi + 1}` : '',
      callRate: r.profile.callRate,
      hetBand: r.profile.hetBand,
      absence: scored?.absence ?? called?.absence,
      ceiling: scored?.ceiling ?? called?.ceiling,
      ratio: scored?.ratio ?? called?.ratio,
      verdict: scored?.verdict,
      excludedBecause: r.usable ? '' : r.why,
    }
  })

  const exportItems = provenance ? [
    {
      label: 'The inferred array',
      filename: `progenitor-inferred-${ref?.side || 'parent'}-${provenance.reportId}.probes`,
      mime: 'text/plain',
      hint: 'The reconstructed genotype as an array file, in the same four columns every other '
        + 'export in this family uses. Drop it into Syngamy as the donor to call parental origin '
        + 'on arrays this run never saw. Marked on its face as inferred rather than measured.',
      build: () => inferredArrayText({
        genotype: ref!.stats.genotype,
        locus: (id) => set!.locus(id),
        products: chosenNames,
        mMin: ref!.stats.mMin,
        contamination: ref!.stats.contamination,
        spuriousAbsence: ref!.stats.spuriousAbsence,
        hRetained: ref!.stats.hRetained,
        side: ref!.side,
        reportId: provenance.reportId,
        generatedAt: provenance.generatedAt,
        build: products[0]?.profile.build.build ?? null,
      }),
    },
    {
      label: 'The report',
      filename: `progenitor-report-${provenance.reportId}.pdf`,
      mime: 'application/pdf',
      hint: 'Letter PDF in the same format as the Syngamy report. Every page states that the '
        + 'reference was inferred, so a page lifted out of context still says so.',
      build: () => buildInferredPdf({
        provenance,
        origins: ref!.origins,
        side: ref!.side,
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
        sampleRows, pairs, Object.fromEntries(ref!.ratios), LIMITS, {
          side: ref!.side,
          paternalGroup: analysis?.paternal ?? null,
          yBearing: Object.fromEntries(usable.map((u) => [u.name, sexOf(u).yBearing])),
        }),
    },
  ] : []

  return (
    <div style={{ paddingBottom: 60 }}>
      <FeatureHeader name="Progenitor" tagline="parental origin with no parent array" />

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
        stats={ref ? (() => {
          const n = (k: string): string =>
            String(ref.origins.filter((o) => o.origin === k).length)
          const [a, b] = ref.side ? ['paternal', 'maternal'] : ['this parent', 'the other parent']
          return [
            { label: 'arrays in', value: String(products.length) },
            { label: a, value: n(a) },
            { label: b, value: n(b) },
            { label: 'unclear', value: n('unclear') },
          ]
        })() : [
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
              disabled={state === 'reading' || building || usable.length < 2
                || (!!analysis && !stale)}
              loading={building}
              onClick={() => { void run() }}
            >
              Run
            </Button>
            {products.length > 0 && (
              <Button
                variant="subtle" size="xs" disabled={state === 'reading' || building}
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
          prompt={'Drop every array from the experiment. Which came from which parent is what '
            + 'this works out; they do not need sorting first, and an array of either parent is '
            + 'not needed. Read in this browser; nothing is uploaded.'}
          disabled={state === 'reading'}
          onFiles={(f) => { void add(f) }}
        >
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 6 }}>
            {rows.map((r) => (
              <div
                key={r.id}
                style={{
                  border: `1px solid ${(() => {
                    const o = ref?.origins.find((x) => x.name === r.name)?.origin
                    if (o && o !== 'unclear') {
                      return o === 'paternal' || o === 'this parent'
                        ? 'var(--om-blue)' : 'var(--om-higher)'
                    }
                    return r.usable ? 'var(--om-blue)' : 'var(--om-border)'
                  })()}`,
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
                  {ref?.origins.find((o) => o.name === r.name)?.origin
                    ?? (r.usable
                      ? `${pct(r.profile.callRate, 0)} call, `
                        + `band ${pct(r.profile.hetBand, 1)}`
                        + (sexOf(r).yBearing ? ', Y' : '')
                      : 'set aside')}
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

      {ref && (
        <Stage
          n={1}
          title="Parental origin"
          state={decided === ref.origins.length ? 'ok' : 'attention'}
          badges={ref.side ? [`${ref.side} reference`] : ['sides not named']}
          headline={`${decided} of ${ref.origins.length} arrays called, with no array of `
            + 'either parent'}
          defaultOpen
        >
          <OriginTable origins={ref.origins} side={ref.side} />
          <Why>
            An array either carries the reconstructed parent's genome or it does not, and the
            other side follows from that. The reference is built from the largest group of
            products sharing a parent, whichever parent that turns out to be; a group is named
            the father's when one of its products carries a Y, which a maternal pronucleus
            cannot. Products that went into the reference are scored against a reference built
            without them, because an array compared with something it helped build reads zero
            absence whatever it is.
          </Why>
        </Stage>
      )}

      {products.length > 0 && (
        <Stage
          n={2}
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
          n={3}
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
          n={4}
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
            n={5}
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
          Drop an experiment and press Run. The arrays are gated, split into groups sharing a
          parent, the father's group named by the Y its products carry, that parent's genotype
          reconstructed, and every array called against it, in one pass and with no array of
          either parent required. At least {MIN_PRODUCTS} products of one parent have to be
          present before anything is reconstructed. Takes the same exports the rest of the tool
          does: Axiom <code>.probes</code>, or a genotype table as CSV or TSV, gzipped or not.
          {int(825657)} markers streams in a few seconds per file.
        </Text>
      )}
    </div>
  )
}

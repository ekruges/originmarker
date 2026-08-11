import { useRef, useState } from 'react'
import { Button, Group, Paper, Progress, Table, Text } from '@mantine/core'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, headerMap, parseRow,
  type ChromStats, type ProbeRow, type SampleProfile,
} from './ingest'
import {
  ProductSet, groupByParent, kinship, MIN_ASCERTAINMENT, MIN_PRODUCTS,
} from './inferredReference'
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

type Tag = 'READ' | 'GATE' | 'PAIR' | 'GROUP' | 'BUILD' | 'CALL' | 'WARN' | 'DONE'
const TAG_COLOUR: Record<Tag, string> = {
  READ: 'var(--om-text-dim)',
  GATE: 'var(--om-text-dim)',
  PAIR: 'var(--om-text-dim)',
  GROUP: 'var(--om-blue)',
  BUILD: 'var(--om-blue)',
  CALL: 'var(--om-text-dim)',
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

/**
 * Hand the main thread back so the browser can paint.
 *
 * Comparing every pair and choosing a depth are hundreds of millions of typed-array reads with
 * no I/O in them, so nothing yields on its own and the tab freezes for the duration. A zero
 * timeout between chunks costs about a millisecond each and turns a frozen page into one with a
 * moving bar. `setTimeout` rather than a microtask: a resolved promise does not let the renderer
 * in, and the whole point is that it does.
 */
const breathe = (): Promise<void> => new Promise((r) => { setTimeout(r, 0) })

/** Chunk size for the pairwise pass. At roughly 3ms per pair this paints about every 40ms. */
const PAIR_CHUNK = 12

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

/**
 * Who ended up in which group, said in words.
 *
 * The concordance matrix below this is the evidence and is unreadable at a glance past about
 * eight products. This states the grouping itself: how many parents the files represent, which
 * arrays belong to each, which one is the father and on what basis, and which group the
 * reference was built from.
 */
function GroupRoster({ groups, names, paternal, yBearing, built }: {
  groups: number[][]
  names: string[]
  paternal: number | null
  yBearing: (boolean | null)[]
  built: number | null
}) {
  if (!groups.length) return null
  return (
    <div style={{ marginTop: 8 }}>
      <Text size="xs" fw={600} mb={6}>
        {groups.length === 1
          ? `One parent across ${names.length} products`
          : `${groups.length} parents across ${names.length} products`}
      </Text>
      {groups.map((g, i) => {
        const ys = g.filter((x) => yBearing[x] === true)
        const isDad = i === paternal
        return (
          <div
            key={i}
            style={{
              borderLeft: `3px solid ${isDad ? 'var(--om-blue)' : 'var(--om-border)'}`,
              padding: '4px 0 4px 8px',
              marginBottom: 6,
            }}
          >
            <div style={{ fontSize: 11, fontWeight: 600 }}>
              Group {i + 1}
              {': '}
              {g.length} product{g.length === 1 ? '' : 's'}
              {isDad ? ' \u00b7 the father' : paternal === null ? '' : ' \u00b7 not the father'}
              {built === i ? ' \u00b7 reference built from this group' : ''}
            </div>
            <div style={{ fontSize: 11, fontFamily: 'var(--om-mono)', marginTop: 2 }}>
              {g.map((x) => names[x] + (yBearing[x] === true ? ' (Y)' : '')).join(', ')}
            </div>
            <div style={{ fontSize: 10, color: 'var(--om-text-dim)', marginTop: 2 }}>
              {isDad
                ? `named the father because ${ys.length} of its ${g.length} products carry a `
                  + 'whole Y, which a maternal cell cannot'
                : paternal === null
                  ? 'which parent this is was not established'
                  : 'no product here carries a Y'}
              {g.length >= MIN_PRODUCTS
                ? `. Clears the ${MIN_PRODUCTS}-product floor.`
                : `. Under the ${MIN_PRODUCTS}-product floor, so it cannot be reconstructed.`}
            </div>
          </div>
        )
      })}
    </div>
  )
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
  // What the page is doing right now, and how far in. Every phase reports through this, not
  // just file reading: a run spends most of its wall clock on work that used to leave the page
  // looking hung because nothing said otherwise.
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
  /**
   * STAGE 1 of three, on its own button.
   *
   * Compare every pair and split the products into groups sharing a parent. Nothing is built
   * and no file is re-read: this is the cheap step, and it is the one whose answer decides
   * whether the expensive ones are worth running at all.
   */
  const compare = async (): Promise<void> => {
    if (!set || usable.length < 2) return
    setRef(null)
    setBuilding(true)
    try {
      const total = usable.length * (usable.length - 1) / 2
      log('PAIR', `comparing ${total} pairs across ${usable.length} products, `
        + `${int(set.markerCount)} markers each`)
      setBusyName(`comparing ${total} pairs of products`)
      setPctDone(0)
      let lo = 1
      let hi = 0
      let done = 0
      for (let a = 0; a < usable.length; a += 1) {
        for (let b = a + 1; b < usable.length; b += 1) {
          const o = set.opposite(usableSlots[a], usableSlots[b])
          lo = Math.min(lo, o.rate)
          hi = Math.max(hi, o.rate)
          log('PAIR', `${usable[a].name} vs ${usable[b].name}: ${pct(o.rate, 2)} opposite `
            + `over ${int(o.shared)} shared homozygous markers, `
            + `${kinship(o.rate).replace('_', ' ')}`)
          done += 1
          // Yield between chunks. Without this the whole pass is one synchronous block and
          // the tab is frozen for the length of it.
          if (done % PAIR_CHUNK === 0) {
            setPctDone(Math.round((done / total) * 100))
            await breathe()
          }
        }
      }
      setPctDone(100)
      log('PAIR', `opposite-homozygote rate ${pct(lo, 2)} to ${pct(hi, 2)} across all `
        + `${total} pairs`)

      setBusyName('grouping products by shared parent')
      await breathe()
      const g = groupByParent(usable.length, rate)
      const pat = paternalGroup(g, usable.map((u) => sexOf(u).yBearing))
      log('GROUP', `${g.length} parent group${g.length === 1 ? '' : 's'} from `
        + `${usable.length} products, by exact maximum clique over the pair table`)
      g.forEach((grp, i) => {
        const names = grp.map((x) => usable[x].name)
        const ys = grp.filter((x) => sexOf(usable[x]).yBearing).map((x) => usable[x].name)
        log(i === pat ? 'DONE' : 'GROUP',
          `group ${i + 1}, ${grp.length} product${grp.length === 1 ? '' : 's'}`
          + `${i === pat ? ', THE FATHER' : ''}: ${names.join(', ')}`)
        log('GROUP', `  group ${i + 1} Y-bearing: `
          + (ys.length ? `${ys.join(', ')} (${ys.length} of ${grp.length})` : 'none'))
        log('GROUP', `  group ${i + 1} ${grp.length >= MIN_PRODUCTS
          ? `clears the ${MIN_PRODUCTS}-product floor and can be reconstructed`
          : `is under the ${MIN_PRODUCTS}-product floor and cannot be reconstructed`}`)
      })
      if (g.length > 1) {
        log('WARN', `${g.length} parent groups. A reference built from all of these would be `
          + 'built from more than one person, which is why only one group is used.')
      }
      if (pat === null) {
        // Three different reasons, and they are not interchangeable. A file with no Y probes at
        // all was never asked the question; a set that answered "no" everywhere was.
        const noProbes = usable.every((u) => sexOf(u).yBearing === null)
        const anyY = usable.some((u) => sexOf(u).yBearing)
        log('WARN', noProbes
          ? 'none of these files carries chromosome Y probes, so nothing here can say which '
            + "group is the father's. The split still holds; only the naming is withheld."
          : anyY
            ? 'more than one group carries a Y, which one sperm donor cannot do. Which side is '
              + 'which is left unnamed.'
            : 'no product carries a Y, so which group is the father cannot be established. A '
              + 'paternal group of n products is all X-bearing 2^-n of the time. The split '
              + 'still holds; only the naming is withheld.')
      }
      setAnalysis({ groups: g, usableIds: usable.map((u) => u.id), paternal: pat })
      const bi = g.findIndex((x) => x.length >= MIN_PRODUCTS)
      log(bi < 0 ? 'WARN' : 'DONE', bi < 0
        ? `no group reaches the ${MIN_PRODUCTS} products this method needs; nothing can be built`
        : `ready to reconstruct group ${bi + 1}`)
    } finally {
      setBuilding(false)
      setBusyName('')
      setPctDone(0)
    }
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
  /**
   * STAGE 2 of three, on its own button.
   *
   * Choose the agreement depth and reconstruct one group's genotype. Reads no files: everything
   * it needs is already in memory as allele codes. Nothing is scored here, because scoring
   * re-reads every array and that is the expensive step this stage exists to let you decide
   * about first.
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
      const side = pat === null ? '' : gi === pat ? 'paternal' : 'maternal'
      log('BUILD', `reconstructing group ${gi + 1}, ${chosen.length} products`
        + `${side ? ` (the ${side} parent)` : ' (which parent is not established)'}`)
      log('BUILD', `products in: ${chosen.map((c) => c.name).join(', ')}`)
      log('BUILD', `${drop.length} of ${set.ids.length} loaded arrays excluded from the build`)

      // The depth ladder, one rung at a time so the bar moves and the log shows the trade the
      // threshold is making rather than only its answer.
      setBusyName('choosing the agreement depth')
      setPctDone(0)
      const ratios = new Map<number, number>()
      const base = set.hRetainedAt(2, drop)
      let mMin = 2
      for (let m = 2; m <= chosen.length; m += 1) {
        const r = set.hRetainedAt(m, drop) / base
        ratios.set(m, r)
        if (Number.isFinite(r) && r >= MIN_ASCERTAINMENT) mMin = m
        log('BUILD', `  m >= ${m}: ascertainment ${pct(r, 1)} of the genome's heterozygosity`
          + `${r >= MIN_ASCERTAINMENT ? '' : `, under the ${pct(MIN_ASCERTAINMENT, 0)} floor`}`)
        setPctDone(Math.round(((m - 1) / Math.max(chosen.length - 1, 1)) * 100))
        await breathe()
      }
      log('BUILD', `depth chosen: m >= ${mMin}, the deepest still retaining `
        + `${pct(MIN_ASCERTAINMENT, 0)} of the genome`)

      setBusyName(`building the reference over ${int(set.markerCount)} markers`)
      setPctDone(0)
      await breathe()
      const full = set.build(mMin, drop)
      setPctDone(100)
      log('BUILD', `${int(full.markers)} markers retained of ${int(set.markerCount)}, `
        + `mean ${full.meanM.toFixed(1)} products per marker`)
      log('BUILD', `parent heterozygosity over those markers ${pct(full.hRetained, 3)}`)
      log('BUILD', `contamination ${pct(full.contamination, 3)}: heterozygous sites every `
        + 'product happened to agree on, so the reference asserts them as homozygous')
      log('BUILD', `a true haploid offspring therefore reads ${pct(full.spuriousAbsence, 3)} `
        + 'absence against this reference for that reason alone')
      setRef({ group: gi, stats: full, ratios, members: [], origins: [], side })
      log('DONE', `reference built. ${products.length} arrays are ready to be called against it`)
    } finally {
      setBuilding(false)
      setBusyName('')
      setPctDone(0)
    }
  }

  /**
   * STAGE 3 of three, on its own button.
   *
   * Score every array against the reference. This is the slow one and it is slow for a reason
   * that cannot be optimised away: an array's absence rate needs its genotypes AND its B-allele
   * frequencies, and the allele codes held in memory carry neither, so every file is streamed
   * again. It is gated behind its own button and reports per file, because a run that silently
   * re-reads twenty arrays looks like a hung page.
   */
  const callOrigins = async (): Promise<void> => {
    if (!set || !ref) return
    setBuilding(true)
    try {
      const gi = ref.group
      const chosen = groups[gi].map((i) => usable[i])
      const keep = new Set(chosen.map((c) => products.find((p) => p.id === c.id)!.name))
      const drop = set.ids.filter((id) => !keep.has(id))
      const { mMin, side } = { mMin: ref.stats.mMin, side: ref.side }
      const other = side === 'paternal' ? 'maternal' : 'paternal'
      const full = ref.stats
      const scored: Scored[] = []
      const called: Origin[] = []
      const inGroup = new Set(chosen.map((c) => c.id))
      log('CALL', `scoring ${products.length} arrays against the reference. Each is streamed `
        + 'again: the absence rate needs B-allele frequencies, which the in-memory allele codes '
        + 'do not carry.')
      // EVERY array is scored, not only the ones the reference was built from. The arrays that
      // did not go into it are the ones the answer is wanted for: an array that carries this
      // parent's genome came from this parent, and one that decisively does not came from the
      // other. Excluded arrays are scored too, and say on their own row that they were excluded,
      // because a fused or half-failed zygote is the case this is most often asked about and
      // silently dropping it answers nothing.
      for (const [n, load] of products.entries()) {
        const mine = inGroup.has(load.id)
        setBusyName(`calling ${load.name}  (${n + 1} of ${products.length})`)
        setPctDone(0)
        await breathe()
        // A product scored against a reference containing itself reads exactly zero absence, a
        // bias the size of the whole signal, so a member is scored leave-one-out. n-1 products
        // cannot reach a threshold of n, so that build sits one below the chosen depth.
        const ref1 = mine
          ? set.build(Math.min(mMin, chosen.length - 1), [...drop, load.name]) : full
        if (mine) {
          log('CALL', `  ${load.name} is in the reference, so it is scored against one rebuilt `
            + `without it: m >= ${Math.min(mMin, chosen.length - 1)}, `
            + `${int(ref1.markers)} markers`)
        }
        const t = emptyTally()
        let map: ReturnType<typeof headerMap> = null
        await eachLine(load.file, (line) => {
          if (!map) { map = headerMap(line); return }
          const r = parseRow(line, map)
          if (r) tallyRow(ref1.genotype.get(r.probesetId) ?? 'NC', r, t)
        }, (bytes) => setPctDone(Math.round((bytes / load.file.size) * 100)))
        // The leave-one-out reference is built from one fewer product, so it is dirtier than the
        // full one and carries its own floor rather than the full build's.
        const res = classify(t, ref1.hRetained, { spuriousAbsence: ref1.spuriousAbsence })
        const ratio = res.genomeRate / res.explainable
        const origin = res.verdict === 'unclear' ? 'unclear'
          : side === '' ? (res.verdict === 'parent_genome_present'
            ? 'this parent' : 'the other parent')
            : res.verdict === 'parent_genome_present' ? side : other
        log(res.verdict === 'unclear' ? 'WARN' : 'DONE',
          `${load.name}: ${pct(res.genomeRate, 2)} absent against a ${pct(res.explainable, 2)} `
          + `ceiling = ${ratio.toFixed(2)}x -> ${origin}`)
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
      setRef({ ...ref, members: scored, origins: called })
    } finally {
      setBuilding(false)
      setBusyName('')
      setPctDone(0)
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
        paternalGroup: analysis?.paternal ?? null,
        builtGroup: ref!.group,
        yBearing: Object.fromEntries(usable.map((u) => [u.name, sexOf(u).yBearing])),
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
      <FeatureHeader name="Progenitor" tagline="builds a parent\u2019s SNP array from the cells it produced" />

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
              loading={building && !ref}
              onClick={() => { void compare() }}
            >
              1. Compare products
            </Button>
            <Button
              size="xs"
              disabled={state === 'reading' || building || !analysis || stale
                || biggest < MIN_PRODUCTS || !!ref}
              onClick={() => {
                const bi = groups.findIndex((g) => g.length >= MIN_PRODUCTS)
                if (bi >= 0) void buildFor(bi)
              }}
            >
              2. Build the reference
            </Button>
            <Button
              size="xs"
              disabled={state === 'reading' || building || !ref || ref.origins.length > 0}
              loading={building && !!ref && !ref.origins.length}
              onClick={() => { void callOrigins() }}
            >
              3. Call {products.length} array{products.length === 1 ? '' : 's'}
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
          prompt={'Drop every array from the experiment. These are the haploid cells; this '
            + 'builds the SNP array of the parent behind them. They do not need sorting first '
            + 'and no array of either parent is needed. Read in this browser; nothing is '
            + 'uploaded.'}
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

        {(state === 'reading' || building) && busyName && (
          <>
            <Text size="xs" c="dimmed" mt={6}>
              {state === 'reading' ? 'Reading ' : ''}{busyName}
            </Text>
            <Progress
              value={pctDone}
              size="sm"
              mt={4}
              radius={0}
              animated={pctDone === 0}
            />
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
          <GroupRoster
            groups={groups}
            names={names}
            paternal={analysis?.paternal ?? null}
            yBearing={usable.map((u) => sexOf(u).yBearing)}
            built={ref?.group ?? null}
          />
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
          {/* Only when the choice is real. With one buildable group the "2. Build" button above
              already does this, and two buttons for one action is worse than none. */}
          {biggest >= MIN_PRODUCTS && !ref
            && groups.filter((g) => g.length >= MIN_PRODUCTS).length > 1 && (
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
          Builds a SNP array for a parent nobody genotyped, out of the haploid cells that parent
          produced. Three steps, each on its own button so you can see what each one costs:{' '}
          <b>compare</b> the products and split them into groups sharing a parent,{' '}
          <b>build</b> the array for one group, then <b>call</b> every input against it for
          parental origin. At least {MIN_PRODUCTS} products of one parent are needed before
          anything is built. The array itself exports as a file you can use anywhere an array
          goes. Takes the same formats the rest of the tool does: Axiom <code>.probes</code>, or
          a genotype table as CSV or TSV, gzipped or not. {int(825657)} markers streams in a few
          seconds per file.
        </Text>
      )}
    </div>
  )
}

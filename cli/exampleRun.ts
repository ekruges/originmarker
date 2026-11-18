/**
 * The example run, off disk, assembled from the modules the page runs and in the page's order.
 *
 * Syngamy.tsx's Run reads every file labelled sperm as a candidate when there are several and
 * scores the samples against the one the search names, indexes the oocyte, then scores each sample
 * against both with the run's other samples as siblings. This does the same with the same calls,
 * so cli/examples.check.ts can pin what each example card prints. Usage, to print a run of any
 * files: node cli/exampleRun.ts sperm:<file> [sperm:<file>] oocyte:<file> sample:<file> ...
 */
import { readFileSync } from 'node:fs'
import { basename, join } from 'node:path'
import { gunzipSync } from 'node:zlib'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, headerMap, parseRow, type ChromStats, type ProbeRow, type SampleProfile,
} from '../web/src/ingest.ts'
import {
  collectParentRow, collectRow, emptyCollected, emptyParent, finishParent, noParent, scoreSample,
  type ParentIndex,
} from '../web/src/scoreSample.ts'
import { agreement, classify, pair, type ParentageResult } from '../web/src/parentage.ts'
import {
  addCandidateRow, candidateGenomes, emptyCandidateTallies, fatherSearchText, findParent,
  type ParentSearchText,
} from '../web/src/findParent.ts'
import { FOUNDER_PANEL_PATH, parsePanel } from '../web/src/answerKeys.ts'
import { runAlerts } from '../web/src/defects.ts'
import { segmentCoords } from '../web/src/segments.ts'
import { cardSays } from '../web/src/pageSays.ts'

export type Role = 'donor' | 'oocyte' | 'sample'

const PANEL = parsePanel(gunzipSync(readFileSync(
  join(new URL('../web/public/', import.meta.url).pathname, FOUNDER_PANEL_PATH))).toString('utf8'))

function read(path: string, onRow: (r: ProbeRow) => void): SampleProfile {
  const raw = readFileSync(path)
  const lines = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8').split('\n')
  let map: ReturnType<typeof headerMap> = null
  const byChrom = new Map<string, ChromStats>()
  const baf = emptyBafSums()
  const builds = emptyBuildSums()
  let first = ''
  for (const line of lines) {
    if (!map) { if (line && !line.startsWith('#')) map = headerMap(line); continue }
    const r = parseRow(line, map)
    if (!r) continue
    if (!first) first = r.probesetId
    accumulate(r, byChrom)
    accumulateBaf(r, baf)
    accumulateBuild(r, builds)
    onRow(r)
  }
  if (!map) throw new Error(`${path}: no recognisable header`)
  const name = basename(path).replace(/\.gz$/, '')
  return { ...finishProfile(name, byChrom, baf, first, builds), build: buildVerdict(builds) }
}

function indexParent(path: string): ParentIndex {
  const acc = emptyParent()
  const profile = read(path, (r) => collectParentRow(r, acc))
  return finishParent(acc, profile.build.build)
}

export interface Card {
  name: string
  /** The pair's class with both parents, else the one-parent class. */
  origin: string
  result: ParentageResult
  maternal?: ParentageResult
  alerts: string[]
}

export interface ExampleRun { father?: ParentSearchText & { named?: string }; cards: Card[] }

export async function runExamples(files: { path: string; role: Role }[]): Promise<ExampleRun> {
  const out: ExampleRun = { cards: [] }
  const sperm = files.filter((f) => f.role === 'donor')
  const samples = files.filter((f) => f.role === 'sample')
  let pat: ParentIndex | null = sperm.length === 1 ? indexParent(sperm[0].path) : null
  if (sperm.length > 1) {
    const found = sperm.map((f) => ({ f, p: indexParent(f.path) }))
    const gts = found.map((x) => x.p.gt)
    const cells = samples.map((s) => {
      const t = emptyCandidateTallies(found.length)
      read(s.path, (r) => addCandidateRow(t, gts, r, PANEL.freqB!))
      return candidateGenomes(t)
    })
    const search = findParent(cells, found.length)
    const names = found.map((x) => basename(x.f.path).replace(/\.gz$/, ''))
    const named = search.named === null ? null : found[search.named]
    pat = named?.p ?? null
    out.father = { ...fatherSearchText(search, names), named: named ? names[search.named!] : undefined }
  }
  const oocyte = files.find((f) => f.role === 'oocyte')
  const mat = oocyte ? indexParent(oocyte.path) : null
  // With no file labelled as a parent the page reads what needs none, as here.
  if (!pat && !sperm.length && !mat) pat = noParent()
  if (!pat) throw new Error('no father: this harness reads the example run, which names one')
  const agree = mat ? agreement(pat.gt, mat.gt) : NaN

  const calls = new Map<string, Map<string, string>>()
  const finished: { s: (typeof samples)[number]; result: ParentageResult; card: Card }[] = []
  const score = async (s: (typeof samples)[number]): Promise<Card> => {
    const acc = emptyCollected(pat!, mat, PANEL)
    const profile = read(s.path, (r) => collectRow(r, pat!, mat, acc))
    const sibs = [...calls].filter(([k]) => k !== s.path).map(([, m]) => m)
    const result = await scoreSample({
      acc, profile, pat: pat!, mat, soloRole: 'paternal', sibs, sampleName: profile.name, log: () => {},
    })
    calls.set(s.path, acc.myGt)
    const maternal = acc.tm && mat ? classify(acc.tm, mat.heterozygosity, { role: 'maternal' }) : undefined
    const paired = maternal ? pair(result, maternal, agree) : undefined
    return {
      name: basename(s.path).replace(/\.gz$/, ''),
      origin: paired?.originClass ?? result.originClass,
      result,
      maternal,
      alerts: runAlerts(result as never).map((a) => a.headline),
    }
  }
  for (const s of samples) {
    const card = await score(s)
    finished.push({ s, result: card.result, card })
  }
  // The page reads a sample again once the run holds every array, where its run was under the span.
  for (const f of finished.slice(0, -1)) {
    if (!f.result.refusals?.some((x) => x.gate === 'other-parent.under-span')) continue
    f.card = await score(f.s)
  }
  out.cards = finished.map((f) => f.card)
  return out
}

/** What a card says, in the terms the example lines use. */
export function summary(c: Card) {
  const r = c.result
  const says = cardSays(r, c.maternal, c.origin, false)
  const rows = (k: 'counted' | 'once' | 'seen') => says.list[k].map((g) => `chr${g.chrom} `
    + `${g.whole ? 'whole' : `${(g.startBp / 1e6).toFixed(1)}-${(g.endBp / 1e6).toFixed(1)}Mb`} `
    + `${g.gain ? 'gain' : 'loss'} ${g.origin} | ${g.word}`)
  return {
    name: c.name,
    head: says.head.word,
    decided: says.head.decided,
    line: says.head.line,
    counted: rows('counted'),
    once: rows('once'),
    seen: rows('seen'),
    origin: c.origin,
    father: r.parentGenome?.verdict,
    mother: r.parentGenomeOther?.verdict,
    lostParent: r.chroms.filter((x) => x.lostParent).map((x) => `chr${x.chrom} ${x.lostParent}`),
    copies: r.chroms.filter((x) => x.aneuploidy).map((x) => `chr${x.chrom} ${x.aneuploidy}`),
    mixed: r.chroms.filter((x) => x.mixed).map((x) => `chr${x.chrom}`),
    // The edges the card, the change list, the PDF and the plot print.
    segments: r.segments.map((g) => segmentCoords(g)).map((co, i) => `chr${r.segments[i].chrom} `
      + `${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb `
      + `${r.segments[i].kind} ${r.segments[i].absentParent ?? ''}${r.segments[i].contested ? ' contested' : ''}`.trim()),
    sex: r.sexChromosome?.constitution,
    /** The run's count of changes, the figure in the card's summary line. */
    changes: (r.integrity as { measures?: { events?: number } } | undefined)?.measures?.events ?? 0,
    stage: r.stage?.stage,
    alerts: c.alerts,
    limits: r.limits,
  }
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const files = process.argv.slice(2).map((a) => {
    const [role, path] = a.split(':')
    return { role: (role === 'sperm' ? 'donor' : role) as Role, path }
  })
  const run = await runExamples(files)
  console.log(JSON.stringify({
    father: run.father && { line: run.father.line, rows: run.father.rows },
    cards: run.cards.map(summary),
  }, null, 1))
}

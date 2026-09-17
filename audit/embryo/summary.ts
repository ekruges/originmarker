/**
 * Stage, abnormality, parent of origin and mechanism, over a whole series.
 *
 * The four questions a batch of embryo arrays is asked: how many samples at each developmental
 * stage, how many chromosome abnormalities at each, which parent each one came from, and whether
 * it was carried by the gamete or arose after fertilisation.
 *
 * Stage and embryo identity are not in the arrays. They come from the series metadata, read by
 * geo_sheet.py, and are the submitters' statements rather than anything this tool inferred.
 *
 * Parent of origin needs both parents genotyped. Where a series arrayed them the parental columns
 * are filled; where it did not they are withheld rather than guessed. A series scored against a
 * stand-in array in the parental slot can still be read for copy number, which is measured against
 * the sample's own chromosomes, but naming a parent from a stranger's genotypes means nothing.
 *
 * Mechanism needs two or more units of ONE embryo; see web/src/timing.ts. An embryo represented by
 * a single array yields 'unresolved', which is counted rather than dropped.
 *
 * Run: OM_SHEET=sheet.json OM_DIRS='GSE148488=/root/trios,GSE186407=/root/GSE186407' \
 *        [OM_TRIOS=trio_manifest_full.csv] [OM_CAP=999] [OM_OUT=summary.json] \
 *        node --experimental-strip-types --max-old-space-size=3072 audit/embryo/summary.ts
 */
import { readFileSync, existsSync, readdirSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const timing = await import(`${W}timing.ts`)

const SHEET = process.env.OM_SHEET ?? ''
const DIRS = new Map<string, string>(
  (process.env.OM_DIRS ?? '').split(',').filter(Boolean)
    .map((p) => p.split('=') as [string, string]),
)
if (!SHEET || !existsSync(SHEET) || !DIRS.size) {
  console.log('embryo summary: needs OM_SHEET and OM_DIRS. NOT RUN.')
  process.exit(0)
}

interface SheetRow {
  gsm: string; series: string; stage: string; embryo: string; role: string
  cells?: number; parentOfOrigin?: string; title?: string
}
const sheet: SheetRow[] = JSON.parse(readFileSync(SHEET, 'utf8'))

/** Every array the series directories hold, by the sample accession its filename starts with. */
const fileOf = new Map<string, string>()
for (const [series, dir] of DIRS) {
  if (!existsSync(dir)) throw new Error(`${series}: ${dir} is not readable`)
  for (const f of readdirSync(dir)) {
    const m = /^(GSM\d+)/.exec(f)
    if (m && (f.endsWith('.gz') || f.endsWith('.probes'))) fileOf.set(m[1], join(dir, f))
  }
}

/** Which arrays are this sample's parents, where the series genotyped them. */
const parentsOf = new Map<string, { father?: string; mother?: string }>()
if (process.env.OM_TRIOS && existsSync(process.env.OM_TRIOS)) {
  const lines = readFileSync(process.env.OM_TRIOS, 'utf8').trim().split('\n')
  const cols = lines[0].split(',')
  const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
  for (const line of lines.slice(1)) {
    const r = line.split(',')
    const pick = (k: string) => at(r, k).split(/[; ]+/).filter((g) => fileOf.has(g))[0]
    parentsOf.set(at(r, 'gsm'), { father: pick('father_gsms'), mother: pick('mother_gsms') })
  }
}

function rowsOf(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

const parentCache = new Map<string, unknown>()
function parentIndex(gsm: string) {
  if (parentCache.has(gsm)) return parentCache.get(gsm)
  const path = fileOf.get(gsm)
  if (!path) return null
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(path)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never,
    ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  // ponytail: no eviction. Five parents across these series, each a few hundred thousand markers.
  parentCache.set(gsm, idx)
  return idx
}

interface Event {
  chrom: string; direction?: 'gain' | 'loss'; whole: boolean
  origin?: string; mechanism?: string
}
interface Row {
  gsm: string; series: string; stage: string; embryo: string; cells?: number
  callRate: number; called: string; originAssessable: boolean
  events: Event[]; error?: string
  result?: unknown
}

const cap = Number(process.env.OM_CAP ?? 9999)
const samples = sheet.filter((s) => s.role === 'sample' && fileOf.has(s.gsm)).slice(0, cap)
console.log(`arrays to score: ${samples.length} of ${sheet.length} samples in the sheet`)
for (const [series, dir] of DIRS) {
  const n = samples.filter((s) => s.series === series).length
  const withParents = samples.filter(
    (s) => s.series === series && (parentsOf.get(s.gsm)?.father || parentsOf.get(s.gsm)?.mother),
  ).length
  console.log(`  ${series}: ${n} arrays from ${dir}, ${withParents} with a parent arrayed`)
}

/** The stand-in for a series with no parents. Copy number does not consult it; origin is withheld. */
const standIn = new Map<string, unknown>()

const rows: Row[] = []
for (const s of samples) {
  const p = parentsOf.get(s.gsm) ?? {}
  let pat = p.father ? parentIndex(p.father) : null
  let mat = p.mother ? parentIndex(p.mother) : null
  const originAssessable = !!(pat && mat)
  if (!pat && !mat) {
    if (!standIn.has(s.series)) {
      standIn.set(s.series, parentIndex(samples.find((x) => x.series === s.series)!.gsm))
    }
    pat = standIn.get(s.series)
  }
  if (!pat && mat) { pat = mat; mat = null }
  try {
    const acc = score.emptyCollected(pat as never, mat as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(fileOf.get(s.gsm)!)()) {
      if (!first) first = r.probesetId
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, pat as never, mat as never, acc as never)
    }
    const profile = ingest.finishProfile(s.gsm, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat, mat, soloRole: p.mother && !p.father ? 'maternal' : 'paternal',
      sibs: [], sampleName: s.gsm, log: () => {},
    }) as Record<string, unknown>
    // Whose copy it is, from the two lists that carry an origin. Keyed by chromosome and
    // direction, which is what the timing events are keyed by.
    const origins = new Map<string, string>()
    for (const [dir, list] of [['gain', res.gains], ['loss', res.losses]] as const) {
      for (const g of (list ?? []) as { where: string; origin: string }[]) {
        const m = /^chr(\S+)/.exec(g.where)
        if (m && !origins.has(`${m[1]}|${dir}`)) origins.set(`${m[1]}|${dir}`, g.origin)
      }
    }
    const events = timing.eventsOf(res as never).map((e: Event) => ({
      chrom: e.chrom,
      direction: e.direction,
      whole: (e as unknown as { endBp: number }).endBp === Number.MAX_SAFE_INTEGER,
      origin: originAssessable && e.direction
        ? origins.get(`${e.chrom}|${e.direction}`) ?? 'unclear' : undefined,
    }))
    rows.push({
      gsm: s.gsm, series: s.series, stage: s.stage, embryo: s.embryo, cells: s.cells,
      callRate: profile.callRate,
      called: (res.stage as { stage?: string } | undefined)?.stage ?? 'none',
      originAssessable, events, result: res,
    })
  } catch (e) {
    rows.push({
      gsm: s.gsm, series: s.series, stage: s.stage, embryo: s.embryo, callRate: NaN,
      called: 'error', originAssessable: false, events: [],
      error: String((e as Error).message ?? e).slice(0, 90),
    })
  }
  if (rows.length % 20 === 0) console.log(`  scored ${rows.length}/${samples.length}`)
}

// --------------------------------------------------------- mechanism, one embryo at a time
const byEmbryo = new Map<string, Row[]>()
for (const r of rows) {
  if (!r.error && r.embryo) {
    const k = `${r.series}|${r.embryo}`
    byEmbryo.set(k, [...(byEmbryo.get(k) ?? []), r])
  }
}
for (const [, group] of byEmbryo) {
  const timed = timing.timeGroup(group.map((g) => g.result as never))
  group.forEach((g, i) => {
    timed[i].forEach((t: { mechanism: string }, j: number) => {
      if (g.events[j]) g.events[j].mechanism = t.mechanism
    })
  })
}
for (const r of rows) for (const e of r.events) if (!e.mechanism) e.mechanism = 'unresolved'

// --------------------------------------------------------------------------- the table
const scored = rows.filter((r) => !r.error)
const usable = scored.filter((r) => r.called !== 'failed')
const stages = [...new Set(rows.map((r) => r.stage))].sort()
const pad = (s: string | number, n: number) => String(s).padEnd(n)
const num = (s: string | number, n: number) => String(s).padStart(n)

console.log('')
console.log('=== SAMPLES BY STAGE. The stage is the submitters\', not this tool\'s.')
console.log(`  ${pad('stage', 22)}${num('arrays', 7)}${num('scored', 8)}${num('usable', 8)}`)
for (const st of stages) {
  const all = rows.filter((r) => r.stage === st)
  console.log(`  ${pad(st, 22)}${num(all.length, 7)}`
    + `${num(all.filter((r) => !r.error).length, 8)}`
    + `${num(all.filter((r) => !r.error && r.called !== 'failed').length, 8)}`)
}
console.log(`  ${pad('TOTAL', 22)}${num(rows.length, 7)}${num(scored.length, 8)}`
  + `${num(usable.length, 8)}`)

console.log('')
console.log('=== ABNORMALITIES BY STAGE. Arrays that resolved to a stage are the denominator.')
console.log(`  ${pad('stage', 22)}${num('arrays', 7)}${num('affected', 9)}${num('rate', 7)}`
  + `${num('events', 8)}${num('whole', 7)}${num('segment', 8)}`)
for (const st of stages) {
  const u = usable.filter((r) => r.stage === st)
  if (!u.length) continue
  const ev = u.flatMap((r) => r.events)
  const aff = u.filter((r) => r.events.length).length
  console.log(`  ${pad(st, 22)}${num(u.length, 7)}${num(aff, 9)}`
    + `${num((aff / u.length).toFixed(2), 7)}${num(ev.length, 8)}`
    + `${num(ev.filter((e) => e.whole).length, 7)}`
    + `${num(ev.filter((e) => !e.whole).length, 8)}`)
}

console.log('')
const assessable = usable.filter((r) => r.originAssessable)
console.log('=== PARENT OF ORIGIN. Needs both parents genotyped, so only where they were.')
console.log(`  arrays with both parents on the array: ${assessable.length} of ${usable.length}`)
for (const whole of assessable.length ? [true, false] : []) {
  const ev = assessable.flatMap((r) => r.events)
    .filter((e) => e.origin && e.whole === whole)
  console.log(`  ${whole ? 'whole chromosomes' : 'segments'}:`)
  const tally = (re: RegExp) => ev.filter((e) => re.test(e.origin ?? '')).length
  console.log(`  events with an origin question: ${ev.length}`)
  console.log(`    paternal:      ${tally(/^paternal/i)}`)
  console.log(`    maternal:      ${tally(/^maternal/i)}`)
  console.log(`    not called:    ${ev.length - tally(/^paternal/i) - tally(/^maternal/i)}`)
}
for (const [series] of DIRS) {
  if (!usable.some((r) => r.series === series && r.originAssessable)) {
    console.log(`  ${series}: NOT ASSESSED. This series did not genotype the parents, so which `
      + 'parent a change came from cannot be said and is left blank.')
  }
}

console.log('')
console.log('=== MECHANISM. Reads two or more units of one embryo; see web/src/timing.ts.')
const multi = [...byEmbryo.values()].filter((g) => g.length >= 2)
console.log(`  embryos represented: ${byEmbryo.size}, of which `
  + `${multi.length} by two or more arrays`)
// WHOLE CHROMOSOMES AND SEGMENTS ARE COUNTED APART. Their false positive rates are not the same
// thing: the whole-chromosome call is measured against stated karyotypes, and a segment call on an
// amplified single cell is not. Pooling them lets segment noise decide the mechanism column,
// because an event present in one unit and absent in the next reads post-zygotic by construction.
const mech = (list: Row[], m: string, whole: boolean) =>
  list.flatMap((r) => r.events).filter((e) => e.mechanism === m && e.whole === whole).length
for (const whole of [true, false]) {
  console.log('')
  console.log(`  ${whole ? 'WHOLE CHROMOSOMES' : 'SEGMENTS'}`)
  console.log(`  ${pad('stage', 22)}${num('meiotic', 9)}${num('post-zyg', 10)}`
    + `${num('reciprocal', 12)}${num('unresolved', 12)}`)
  for (const st of stages) {
    const u = usable.filter((r) => r.stage === st)
    if (!u.some((r) => r.events.some((e) => e.whole === whole))) continue
    console.log(`  ${pad(st, 22)}${num(mech(u, 'meiotic', whole), 9)}`
      + `${num(mech(u, 'post-zygotic', whole), 10)}`
      + `${num(mech(u, 'reciprocal', whole), 12)}${num(mech(u, 'unresolved', whole), 12)}`)
  }
  console.log(`  ${pad('TOTAL', 22)}${num(mech(usable, 'meiotic', whole), 9)}`
    + `${num(mech(usable, 'post-zygotic', whole), 10)}`
    + `${num(mech(usable, 'reciprocal', whole), 12)}`
    + `${num(mech(usable, 'unresolved', whole), 12)}`)
}

console.log('')
console.log('=== CHROMOSOME 21')
const t21 = usable.filter((r) => r.events.some((e) => e.chrom === '21' && e.direction === 'gain'))
const m21 = usable.filter((r) => r.events.some((e) => e.chrom === '21' && e.direction === 'loss'))
console.log(`  arrays with a chr21 gain: ${t21.length} of ${usable.length}`)
console.log(`  arrays with a chr21 loss: ${m21.length} of ${usable.length}`)
for (const r of t21) {
  const e = r.events.find((x) => x.chrom === '21' && x.direction === 'gain')!
  console.log(`    ${r.gsm}  ${r.stage}  embryo ${r.embryo || '?'}  `
    + `${e.whole ? 'whole chromosome' : 'segment'}  ${e.mechanism}`
    + `${e.origin ? `  origin ${e.origin}` : ''}`)
}

const errs = rows.filter((r) => r.error)
if (errs.length) {
  console.log('')
  console.log(`=== THREW: ${errs.length}`)
  for (const r of errs.slice(0, 10)) console.log(`  ${r.gsm}  ${r.error}`)
}

if (process.env.OM_OUT) {
  writeFileSync(process.env.OM_OUT, JSON.stringify(
    rows.map(({ result, ...r }) => r), null, 1,
  ))
  console.log(`\nwrote ${process.env.OM_OUT}`)
}

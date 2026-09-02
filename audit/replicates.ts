/**
 * Same DNA, twice. Does the tool say the same thing?
 *
 * THIS IS THE FIRST QUESTION A REVIEWER ASKS AND THE TOOL HAD NEVER BEEN ASKED IT. Every other
 * harness here compares an answer against a truth set. This one compares the tool against ITSELF,
 * on material where there is no biological difference to find: GSE148488 carries technical
 * replicates, arrays run from one biopsy's amplified DNA, labelled `_rep1`, `_rep2` and so on. Any
 * disagreement between them is the assay and the tool, because the genome is identical.
 *
 * WHAT DISAGREEMENT WOULD MEAN. A clinic re-running a sample and getting a different chromosome, or
 * a different parent, has no way to tell which run to believe. Reproducibility is not a nicety here;
 * it is the difference between a measurement and an opinion. So this reports the rate rather than
 * asserting a threshold nobody has measured, and it reports WHERE the disagreements fall, since a
 * split over whether a marginal event exists is a different problem from a split over which parent
 * it came from.
 *
 * NOT EXPECTED TO BE PERFECT, AND THAT IS THE POINT. Amplification from a single cell is a stochastic
 * process; two arrays of one biopsy differ in which alleles dropped out. What must NOT happen is the
 * two runs naming DIFFERENT PARENTS for the same event, which is a contradiction rather than a
 * difference in sensitivity.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 audit/replicates.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('replicates: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

/** The reference every sample is scored against, so the parent is never the variable. */
const REF = 'GSM4472397'

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')],
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

/**
 * Replicate groups, from the GEO titles.
 *
 * A title ending `_repN` names one biopsy run more than once. Stripping that suffix is the whole
 * grouping rule; nothing here decides which arrays are alike by looking at the data, which would be
 * circular.
 */
const groups = new Map<string, string[]>()
for (const r of recs) {
  if (r.role === 'parent' || !has(r.gsm)) continue
  const m = /^(.*)_rep\d+$/.exec(r.title)
  if (!m) continue
  const key = m[1]
  groups.set(key, [...(groups.get(key) ?? []), r.gsm])
}
const usable = [...groups.entries()].filter(([, ids]) => ids.length >= 2)

const textOf = (gsm: string) =>
  gunzipSync(readFileSync(join(DIR, `${gsm}.probes.gz`))).toString('utf8')
function rowsFromText(t: string) {
  const ls = t.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

const parent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(REF))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile(REF, byChrom as never, bafSums as never, first).build.build)
})()

interface Answer {
  gsm: string
  stage: string
  zygosity: string
  integrity: string
  /** Whole-chromosome calls, as `chrN:loss`. */
  aneuploidy: string[]
  /** Segment calls, by chromosome only: two runs will not agree on a breakpoint to the base. */
  segments: string[]
  /** Every named parent, keyed by exact locus. */
  parents: Record<string, string>
  /**
   * The same, keyed by CHROMOSOME.
   *
   * TWO REPLICATES DO NOT AGREE ON A BREAKPOINT TO THE BASE and are not expected to: the edges are
   * refined from where markers stop supporting the call, and which markers dropped out differs per
   * array. Comparing exact interval labels therefore reports a difference of a megabase as though
   * it were a difference of opinion. The chromosome is the unit an operator acts on.
   */
  parentsByChrom: Record<string, string>
  callRate: number
}

async function answerFor(gsm: string): Promise<Answer> {
  const acc = score.emptyCollected(parent as never, null)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(gsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, parent as never, null, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const res = await score.scoreSample({
    acc, profile, pat: parent, mat: null, soloRole: 'paternal', sibs: [],
    sampleName: gsm, log: () => {},
  }) as Record<string, unknown>

  const parents: Record<string, string> = {}
  for (const o of (res.oneParent ?? []) as {
    where: string; parent?: string | null; twoParents?: boolean; verdict: string
  }[]) {
    const named = o.twoParents ? o.parent
      : o.verdict === 'known-parent-lost' ? 'paternal'
        : o.verdict === 'other-parent-lost' ? 'maternal' : null
    if (named) parents[o.where] = named
  }
  for (const d of (res.dosageCalls ?? []) as { where: string; parent?: string | null }[]) {
    if (d.parent) parents[d.where] = d.parent
  }
  const parentsByChrom: Record<string, string> = {}
  for (const [where, who] of Object.entries(parents)) {
    const c = /^chr([\dXY]+)/.exec(where)?.[1]
    if (c) parentsByChrom[`chr${c}`] = who
  }
  return {
    gsm,
    parentsByChrom,
    stage: (res.stage as { stage?: string } | undefined)?.stage ?? 'none',
    zygosity: String(res.zygosity ?? 'none'),
    integrity: (res.integrity as { level?: string } | undefined)?.level ?? 'none',
    aneuploidy: ((res.chroms ?? []) as { chrom: string; aneuploidy?: string }[])
      .filter((c) => c.aneuploidy).map((c) => `chr${c.chrom}:${c.aneuploidy}`).sort(),
    segments: [...new Set(((res.segments ?? []) as { chrom: string; kind: string }[])
      .map((s) => `chr${s.chrom}:${s.kind}`))].sort(),
    parents,
    callRate: (profile as { callRate: number }).callRate,
  }
}

console.log(`replicate groups: ${usable.length}, arrays: `
  + `${usable.reduce((a, [, ids]) => a + ids.length, 0)}, scored against ${REF}`)
console.log('')

let groupsRun = 0
let sameStage = 0
let sameZygosity = 0
let sameIntegrity = 0
let sameAneuploidy = 0
let sameSegments = 0
let sameParents = 0
let sameParentChroms = 0
/** The one that is not a difference in sensitivity but a contradiction. */
let contradictions = 0
const contradictionDetail: string[] = []

for (const [key, ids] of usable) {
  const answers: Answer[] = []
  for (const g of ids) {
    try { answers.push(await answerFor(g)) } catch (e) {
      console.log(`  ${g}: threw ${String((e as Error).message ?? e).slice(0, 60)}`)
    }
  }
  if (answers.length < 2) continue
  groupsRun += 1
  const allSame = <T>(f: (a: Answer) => T) =>
    new Set(answers.map((a) => JSON.stringify(f(a)))).size === 1
  const st = allSame((a) => a.stage)
  const zy = allSame((a) => a.zygosity)
  const ig = allSame((a) => a.integrity)
  const an = allSame((a) => a.aneuploidy)
  const sg = allSame((a) => a.segments)
  const pa = allSame((a) => a.parents)
  const pc = allSame((a) => a.parentsByChrom)
  sameStage += st ? 1 : 0
  sameZygosity += zy ? 1 : 0
  sameIntegrity += ig ? 1 : 0
  sameAneuploidy += an ? 1 : 0
  sameSegments += sg ? 1 : 0
  sameParents += pa ? 1 : 0
  sameParentChroms += pc ? 1 : 0

  // A CONTRADICTION is the same locus given two DIFFERENT parents. That is not one run seeing less
  // than another; it is the tool answering the question two ways.
  const byLocus = new Map<string, Set<string>>()
  for (const a of answers) {
    for (const [where, who] of Object.entries(a.parents)) {
      byLocus.set(where, (byLocus.get(where) ?? new Set()).add(who))
    }
  }
  for (const [where, whos] of byLocus) {
    if (whos.size > 1) {
      contradictions += 1
      contradictionDetail.push(`${key.slice(-40)} ${where}: ${[...whos].join(' vs ')}`)
    }
  }

  const short = key.length > 46 ? `...${key.slice(-43)}` : key
  const flags = [st ? '' : 'stage', zy ? '' : 'zygosity', ig ? '' : 'integrity',
    an ? '' : 'aneuploidy', sg ? '' : 'segments',
    pa ? '' : 'parent-intervals', pc ? '' : 'PARENT-CHROMOSOMES'].filter(Boolean)
  console.log(`${short.padEnd(48)} n=${answers.length} `
    + `call ${answers.map((a) => a.callRate.toFixed(2)).join('/')} `
    + `${flags.length ? `DIFFER: ${flags.join(',')}` : 'identical'}`)
  if (flags.length) {
    for (const a of answers) {
      console.log(`      ${a.gsm} ${a.stage}/${a.zygosity}/${a.integrity} `
        + `aneu[${a.aneuploidy.join(' ')}] seg[${a.segments.join(' ')}] `
        + `parents${JSON.stringify(a.parents)}`)
    }
  }
}

const pct = (n: number) => (groupsRun ? `${n}/${groupsRun} (${(100 * n / groupsRun).toFixed(0)}%)` : 'n/a')
console.log('')
console.log(`groups compared: ${groupsRun}`)
console.log(`  same material inferred     ${pct(sameStage)}`)
console.log(`  same zygosity              ${pct(sameZygosity)}`)
console.log(`  same integrity verdict     ${pct(sameIntegrity)}`)
console.log(`  same whole-chromosome set  ${pct(sameAneuploidy)}`)
console.log(`  same segment set           ${pct(sameSegments)}`)
console.log(`  same named parents, exact  ${pct(sameParents)}`)
console.log(`  same named parents, by chr ${pct(sameParentChroms)}`
  + '   <- the unit an operator acts on')
console.log('')
console.log(`CONTRADICTIONS, one locus given two different parents: ${contradictions}`)
for (const d of contradictionDetail) console.log(`  ${d}`)
console.log('')
console.log('replicates: a difference in what was DETECTED is the assay; a difference in WHICH '
  + 'PARENT is the tool answering one question two ways.')
process.exit(contradictions ? 1 : 0)

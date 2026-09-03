/**
 * DOES THE TOOL FIND THE EVENT WHERE THE EXPERIMENT PUT IT?
 *
 * THE STRONGEST TRUTH AVAILABLE, AND IT IS NOT SYNTHETIC. These experiments targeted a named locus
 * with an editing reagent. The chromosome that was cut is known before the array is read, from the
 * experimental design rather than from this tool, so every array carries a prior about where a
 * change should be. That makes the whole corpus a sensitivity set with a POSITION attached, and the
 * chromosomes that were not targeted are a specificity set on the same arrays at the same time.
 *
 * Two more groups carry a different kind of known answer: a sample karyotyped euploid, and a set of
 * experimental controls. Both should come back quiet.
 *
 * WHAT ON-TARGET ENRICHMENT MEANS HERE. A tool that reports changes at random would show the same
 * event rate on the targeted chromosome as on the other twenty-one. A tool that finds the biology
 * shows more on the target. That ratio is the measurement, and it needs no injection, no
 * construction and no assumption about what the answer should be.
 *
 * NO SAMPLE IS NAMED IN THIS FILE. The rules below match on locus and karyotype words, which are
 * biology, and the directory is supplied at run time. Per-array detail goes to the run's own log on
 * the machine holding the data and is never committed. See scripts/release-check.sh, which enforces
 * that no lab identifier is tracked.
 *
 * DETECTION ONLY, AND DELIBERATELY. Most of these arrays have no paired parental array, so parent of
 * origin cannot be asked of them. What can be asked is whether the change is seen at all and where.
 * The call-rate and intensity channels do not read the parental array, which is why a fixed
 * reference can sit in that slot without touching the answer.
 *
 * Run: OM_LAB=<dir> OM_TRIOS=<dir> node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/lab-targets.ts
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const LAB = process.env.OM_LAB
const DIR = process.env.OM_TRIOS
if (!LAB || !existsSync(LAB) || !DIR || !existsSync(DIR)) {
  console.log('lab-targets: OM_LAB and OM_TRIOS must both point at readable directories. NOT RUN.')
  process.exit(0)
}

/**
 * What each experiment aimed at, from words that name a locus rather than a sample.
 *
 * `chrom` is where a change is EXPECTED. Everything else on that array is the specificity arm.
 */
const TARGETS: { match: RegExp; chrom: string; label: string }[] = [
  { match: /chr\s*16|16\s*p|16\s*q/i, chrom: '16', label: 'chromosome 16 targeted' },
  { match: /dmd/i, chrom: 'X', label: 'DMD, on chromosome X' },
  { match: /hbb/i, chrom: '11', label: 'HBB, on chromosome 11' },
]
/** Groups whose known answer is "nothing here". */
const QUIET: { match: RegExp; label: string }[] = [
  { match: /euploid/i, label: 'karyotyped euploid' },
  { match: /control/i, label: 'experimental control' },
]

function walk(dir: string, out: string[] = []): string[] {
  for (const e of readdirSync(dir)) {
    const p = join(dir, e)
    let st
    try { st = statSync(p) } catch { continue }
    if (st.isDirectory()) walk(p, out)
    else if (e.endsWith('.probes') || e.endsWith('.probes.gz')) out.push(p)
  }
  return out
}

function rowsOf(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8')
    : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  if (h < 0) throw new Error('no header')
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) {
      const r = ingest.parseRow(ls[i], map)
      if (r) yield r
    }
  }
}

/**
 * A fixed array in the parental slot.
 *
 * The channels this harness reads, the call rate and the intensity, never consult it: both were
 * measured to be parent-independent. It is here because `scoreSample` takes a parent, not because
 * the answer depends on one, and using ONE for every array means it cannot be a variable.
 */
const REF = join(DIR, 'GSM4472397.probes.gz')
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(REF)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Seen {
  path: string
  group: string
  target: string | null
  stage: string
  callRate: number
  sexY: boolean | null
  eventChroms: string[]
  onTarget: boolean
  offTargetCount: number
  failed: boolean
  error?: string
}

const files = walk(LAB)
const perGroupCap = Number(process.env.OM_PER_GROUP ?? 20)

/** Which known-answer group a path belongs to, if any. Directory words, never a sample name. */
function groupOf(path: string): { group: string; target: string | null } | null {
  for (const q of QUIET) if (q.match.test(path)) return { group: q.label, target: null }
  for (const t of TARGETS) if (t.match.test(path)) return { group: t.label, target: t.chrom }
  return null
}

const buckets = new Map<string, string[]>()
for (const f of files) {
  const g = groupOf(f)
  if (!g) continue
  const key = g.group
  const cur = buckets.get(key) ?? []
  if (cur.length < perGroupCap) buckets.set(key, [...cur, f])
}

console.log(`lab arrays found: ${files.length}`)
console.log(`groups with a known answer: ${buckets.size}, at most ${perGroupCap} arrays each`)
for (const [g, fs] of buckets) console.log(`  ${g}: ${fs.length}`)
console.log('')

const results: Seen[] = []
for (const [group, fs] of buckets) {
  const target = groupOf(fs[0])!.target
  for (const path of fs) {
    try {
      const acc = score.emptyCollected(refParent as never, null)
      const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
      for (const r of rowsOf(path)()) {
        if (!first) first = r.probesetId
        ingest.accumulate(r as never, byChrom as never)
        ingest.accumulateBaf(r as never, bafSums as never)
        score.collectRow(r as never, refParent as never, null, acc as never)
      }
      const profile = ingest.finishProfile('x', byChrom as never, bafSums as never, first)
      const out = await score.scoreSample({
        acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
        sampleName: 'x', log: () => {},
      }) as Record<string, unknown>
      const chroms = (out.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
      const segs = (out.segments ?? []) as { chrom: string }[]
      const eventChroms = [...new Set([
        ...chroms.filter((c) => c.aneuploidy).map((c) => c.chrom),
        ...segs.map((s) => s.chrom),
      ])]
      const stage = (out.stage as { stage?: string } | undefined)?.stage ?? 'none'
      results.push({
        path, group, target, stage,
        callRate: profile.callRate,
        sexY: (out.sex as { yBearing?: boolean | null } | undefined)?.yBearing ?? null,
        eventChroms,
        onTarget: target ? eventChroms.includes(target) : false,
        offTargetCount: eventChroms.filter((c) => c !== target).length,
        failed: stage === 'failed',
      })
    } catch (e) {
      results.push({
        path, group, target, stage: 'error', callRate: NaN, sexY: null, eventChroms: [],
        onTarget: false, offTargetCount: 0, failed: false,
        error: String((e as Error).message ?? e).slice(0, 70),
      })
    }
  }
}

// ---------------------------------------------------------------- per array, to the local log only
console.log('=== PER ARRAY (this log stays on the machine holding the data)')
for (const r of results) {
  const name = r.path.slice(LAB.length + 1)
  console.log(`  ${r.stage.padEnd(14)} call ${Number.isFinite(r.callRate) ? r.callRate.toFixed(3) : ' n/a '} `
    + `${r.target ? (r.onTarget ? 'ON-TARGET ' : 'no-target ') : '          '}`
    + `off ${String(r.offTargetCount).padStart(2)}  ${r.eventChroms.length ? `chr${r.eventChroms.join(',chr')}` : '-'}`
    + `   ${name}${r.error ? `  ERROR ${r.error}` : ''}`)
}

// ---------------------------------------------------------------- the measurement
console.log('')
console.log('=== ON-TARGET ENRICHMENT, which is the whole question')
console.log('group                          n  usable  on-target  off-target chroms/array  crashed')
let anyCrash = 0
for (const [group] of buckets) {
  const g = results.filter((r) => r.group === group)
  const crashed = g.filter((r) => r.error).length
  anyCrash += crashed
  const usable = g.filter((r) => !r.failed && !r.error)
  const onT = usable.filter((r) => r.onTarget).length
  const offPer = usable.length
    ? usable.reduce((a, r) => a + r.offTargetCount, 0) / usable.length : NaN
  const tgt = g[0]?.target
  console.log(`${group.padEnd(30)} ${String(g.length).padStart(2)}  ${String(usable.length).padStart(6)}  `
    + `${tgt ? `${onT}/${usable.length}`.padStart(9) : '        -'}  `
    + `${Number.isFinite(offPer) ? offPer.toFixed(2).padStart(20) : ''.padStart(20)}  `
    + `${crashed}`)
}

console.log('')
console.log('READ IT LIKE THIS. For a targeted group, on-target counts arrays showing a change on the')
console.log('chromosome the experiment cut, and off-target is the average number of OTHER chromosomes')
console.log('carrying one on the same arrays. On-target well above the off-target rate per chromosome')
console.log('is the tool finding the biology. A quiet group should show few of either.')
console.log('')
console.log(anyCrash === 0
  ? 'lab-targets: no array crashed the pipeline.'
  : `lab-targets: ${anyCrash} array(s) threw. That is a failure whatever the biology says.`)
process.exit(anyCrash ? 1 : 0)

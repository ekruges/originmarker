/**
 * How close is each call to the threshold that produced it?
 *
 * A TOOL WHOSE ANSWERS SIT ON ITS OWN BOUNDARIES IS NOT DEFENSIBLE, however well those boundaries
 * were measured. If a whole-chromosome loss is called because a statistic cleared its threshold by
 * one percent, then re-running the same biopsy, or moving the constant by a rounding error, changes
 * the answer. Every threshold in this tree was measured on real material and has a documented
 * separation; what none of them had was a check that the SHIPPED CALLS actually live in the wide
 * part of that separation rather than crowded against the line.
 *
 * WHAT THIS REPORTS. For every call the tool makes across the corpus, the deciding statistic and
 * the threshold it was compared against, as a ratio. A ratio of 1.0 is exactly on the line. The
 * distribution of those ratios is the answer: a tool whose calls cluster at 3x its thresholds is
 * robust to any plausible revision of them, and one whose calls cluster at 1.05x is not.
 *
 * IT ASSERTS NOTHING IT HAS NOT MEASURED. There is no published figure for how much margin is
 * enough. The number this prints is the evidence for that judgement, not a pass or a fail, and the
 * one thing it does fail on is a call whose statistic is on the wrong side of its own threshold,
 * which would mean the reported reason is not the reason.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 audit/margins.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const parentage = await import(`${W}parentage.ts`)
const stageMod = await import(`${W}stage.ts`)
const nullMod = await import(`${W}intensityNull.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('margins: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

const REF = 'GSM4472397'
const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return { gsm: f[col('gsm')], material: f[col('material')], role: f[col('role')] }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

const textOf = (g: string) => gunzipSync(readFileSync(join(DIR, `${g}.probes.gz`))).toString('utf8')
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

const refGt = (() => {
  const m = new Map<string, string>()
  for (const r of rowsFromText(textOf(REF))()) m.set(r.probesetId, r.genotype)
  return m
})()

interface Margin {
  gsm: string
  material: string
  what: string
  statistic: number
  threshold: number
  /** How many times the threshold the statistic reached. 1.0 is exactly on the line. */
  ratio: number
}
const margins: Margin[] = []
let failures = 0

const arrays = recs.filter((r) => r.role !== 'parent' && has(r.gsm))
  .slice(0, Number(process.env.OM_ALL_N ?? 1000))
console.log(`arrays: ${arrays.length}, scored against ${REF}`)
console.log('')

for (const a of arrays) {
  const t = parentage.emptyTally()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  let pHet = 0
  let pCalled = 0
  for (const r of rowsFromText(textOf(a.gsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    const pg = refGt.get(r.probesetId) ?? 'NC'
    if (parentage.isAutosome(r.chrom) && pg !== 'NC') {
      pCalled += 1
      if (pg === 'AB') pHet += 1
    }
    parentage.tallyRow(pg as never, r as never, t as never)
  }
  const profile = ingest.finishProfile(a.gsm, byChrom as never, bafSums as never, first)
  const st = stageMod.inferStage(profile as never, {}) as { stage: string }
  // A refused array has no calls to have margins on, and its statistics are not measuring a genome.
  if (st.stage === 'failed') continue
  const res = parentage.classify(t as never, pCalled ? pHet / pCalled : NaN,
    { role: 'paternal' }) as never as {
      chroms: {
        chrom: string; aneuploidy?: string; callFraction: number; lrrShift: number
      }[]
      hetBand: number
      homBand?: number
    }

  // --- THE CALL RATE that separates a collapsed chromosome from an intact one.
  for (const c of res.chroms) {
    if (!c.aneuploidy) continue
    const shift = Math.abs(c.lrrShift)
    if (Number.isFinite(shift)) {
      margins.push({
        gsm: a.gsm, material: a.material, what: `chr${c.chrom} ${c.aneuploidy} intensity`,
        statistic: shift, threshold: parentage.COPY_SHIFT_FLOOR,
        ratio: shift / parentage.COPY_SHIFT_FLOOR,
      })
      // A call whose own statistic sits under the floor it is supposed to have cleared means the
      // reported reason is not the reason. That is a failure, not a thin margin, UNLESS the call
      // came in on the collapsed-call-rate channel instead, which has its own threshold.
      if (shift < parentage.COPY_SHIFT_FLOOR && !(c.callFraction < parentage.CALL_COLLAPSE)) {
        failures += 1
        console.log(`  *** ${a.gsm} chr${c.chrom}: called ${c.aneuploidy} at |lrr| `
          + `${shift.toFixed(3)}, under the ${parentage.COPY_SHIFT_FLOOR} floor, and its call `
          + `rate ${c.callFraction.toFixed(3)} did not collapse either`)
      }
    }
    if (c.callFraction < parentage.CALL_COLLAPSE) {
      margins.push({
        gsm: a.gsm, material: a.material, what: `chr${c.chrom} ${c.aneuploidy} call rate`,
        statistic: parentage.CALL_COLLAPSE - c.callFraction, threshold: parentage.CALL_COLLAPSE,
        // How far BELOW the collapse threshold, as a fraction of it. Bigger is safer.
        ratio: 1 + (parentage.CALL_COLLAPSE - c.callFraction) / parentage.CALL_COLLAPSE,
      })
    }
  }

  // --- ZYGOSITY, the band excess that decides one parental contribution from two.
  const excess = res.hetBand - (res.homBand ?? 0)
  if (Number.isFinite(excess)) {
    margins.push({
      gsm: a.gsm, material: a.material, what: 'zygosity band excess',
      statistic: excess, threshold: parentage.HET_BAND_EXCESS,
      ratio: excess / parentage.HET_BAND_EXCESS,
    })
  }

  // --- THE ARRAY'S OWN CALL RATE against the floor below which nothing is asserted.
  margins.push({
    gsm: a.gsm, material: a.material, what: 'array call rate',
    statistic: profile.callRate, threshold: parentage.CALL_RATE_FLOOR,
    ratio: profile.callRate / parentage.CALL_RATE_FLOOR,
  })
}

const q = (xs: number[], p: number) => {
  if (!xs.length) return NaN
  const s = [...xs].sort((x, y) => x - y)
  return s[Math.min(s.length - 1, Math.floor(p * s.length))]
}

console.log('what                              n     min     p10     p50     max   within 10% of line')
const kinds = [...new Set(margins.map((m) => m.what.replace(/^chr\S+ /, '')))].sort()
for (const k of kinds) {
  const xs = margins.filter((m) => m.what.replace(/^chr\S+ /, '') === k).map((m) => m.ratio)
    .filter(Number.isFinite)
  if (!xs.length) continue
  const knife = xs.filter((x) => x > 0.9 && x < 1.1).length
  console.log(`${k.padEnd(32)} ${String(xs.length).padStart(4)}  ${q(xs, 0).toFixed(2)}   `
    + `${q(xs, 0.1).toFixed(2)}   ${q(xs, 0.5).toFixed(2)}   ${q(xs, 1).toFixed(2)}   `
    + `${knife}/${xs.length}`)
}

// THE CALLS NEAREST THEIR OWN LINE, named, so a reader can go and look at them.
console.log('')
console.log('the ten calls closest to the threshold that produced them:')
const near = margins.filter((m) => Number.isFinite(m.ratio))
  .sort((x, y) => Math.abs(x.ratio - 1) - Math.abs(y.ratio - 1)).slice(0, 10)
for (const m of near) {
  console.log(`  ${m.ratio.toFixed(3)}x  ${m.gsm} ${m.material.padEnd(14)} ${m.what} `
    + `(${m.statistic.toFixed(4)} against ${m.threshold})`)
}

console.log('')
console.log(failures === 0
  ? 'margins: every call cleared the threshold its reason names.'
  : `margins: ${failures} call(s) did not clear the threshold their reason names.`)
process.exit(failures ? 1 : 0)

/**
 * Whether any statistic separates a true parent from a stranger on amplified material, where the
 * pooled opposite-homozygote rate does not.
 *
 * `OPPOSITE_HOM_MAX` holds on bulk DNA only. On amplified cells dropout turns a heterozygote into
 * whichever homozygote survived, manufacturing opposite homozygotes out of a real relative. Those
 * are scattered at random, while unrelatedness makes them dense everywhere: the same average, a
 * different shape. The hypothesis under test is that the window floor from `windowedIbs0` keeps the
 * shape. Against a true parent some 200-marker windows escape dropout and read clean, so the median
 * window rate stays low; against a stranger no window can read clean.
 *
 * Truth: `father_gsms` and `mother_gsms` in the manifest, established by the experiment.
 *
 * Three arms, per material, since a statistic that separates on bulk and not on trophectoderm is
 * not a fix and pooling the materials would hide it:
 *   TRUE        the child against its own manifest parent
 *   STRANGER    the child against a donor who is not its parent, with every replicate of the true
 *               parent excluded by identity, never by accession
 *   REPLICATE   the child against a DIFFERENT array of its true parent. It must read with TRUE;
 *               if it reads with STRANGER, the stranger arm is measuring batch and no conclusion
 *               here holds.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=8192 \
 *        audit/relationship-separability.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const rel = await import(`${W}relatedness.ts`)
const stageMod = await import(`${W}stage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('relationship-separability: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}
const CONTAMINATED = 'GSM4472408'

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')], pronucleus: f[col('pronucleus')],
    fathers: (f[col('father_gsms')] ?? '').split(';').filter(Boolean),
    mothers: (f[col('mother_gsms')] ?? '').split(';').filter(Boolean),
    complete: f[col('complete_trio')] === 'True',
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Loaded {
  gt: Map<string, AB>
  pos: Map<string, { chrom: string; pos: number }>
  stage: string
}
const cache = new Map<string, Loaded>()
function load(gsm: string, keep = false): Loaded {
  const hit = cache.get(gsm)
  if (hit) return hit
  const text = gunzipSync(readFileSync(join(DIR, `${gsm}.probes.gz`))).toString('utf8')
  const ls = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  const gt = new Map<string, AB>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    if (!first) first = r.probesetId
    gt.set(r.probesetId, r.genotype as AB)
    if (Number.isFinite(r.pos)) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
  }
  // The tool's OWN material call, so the split below is the one an operator would land in rather
  // than a label from the manifest.
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  const st = stageMod.inferStage(profile as never) as { stage?: string }
  const v: Loaded = { gt, pos, stage: st.stage ?? 'unknown' }
  if (keep) cache.set(gsm, v)
  return v
}

/** Who a donor array actually is, collapsed across its replicates, so a stranger is a stranger. */
const who = (t: string) => (/sperm/i.test(t) ? 'sperm' : (/egg donor ([A-Z])/i.exec(t)?.[1] ?? t))
const donors = recs.filter((r) => r.role === 'parent' && has(r.gsm) && r.gsm !== CONTAMINATED)
const identityOf = new Map(donors.map((d) => [d.gsm, who(d.title)]))

const children = recs.filter((r) => r.complete && r.role !== 'parent' && has(r.gsm))

interface Pair {
  arm: 'TRUE' | 'STRANGER' | 'REPLICATE'
  material: string
  opp: number; n: number
  floor: number; spread: number; windows: number
  verdict: string
}
const pairs: Pair[] = []

function measure(parentGsm: string, childGsm: string, arm: Pair['arm'], material: string) {
  const p = load(parentGsm, true)
  const c = load(childGsm)
  const r = rel.relate({ gt: p.gt, pos: p.pos }, { gt: c.gt, pos: c.pos }, rel.OPPOSITE_HOM_MAX)
  pairs.push({
    arm, material, opp: r.opp.rate, n: r.opp.n,
    floor: r.win.floor, spread: r.win.spread, windows: r.win.windows,
    verdict: r.relationship,
  })
}

let done = 0
for (const c of children) {
  const trueParents = [...c.fathers, ...c.mothers].filter((g) => has(g) && g !== CONTAMINATED)
  if (!trueParents.length) continue
  const child = load(c.gsm)
  const material = child.stage
  // Every identity that is a true parent of THIS child, so no replicate leaks into the stranger arm.
  const trueIds = new Set(trueParents.map((g) => identityOf.get(g)).filter(Boolean))

  measure(trueParents[0], c.gsm, 'TRUE', material)
  // A second array of the SAME person, when the series carries one.
  const repl = trueParents.slice(1).find((g) => identityOf.get(g) === identityOf.get(trueParents[0]))
  if (repl) measure(repl, c.gsm, 'REPLICATE', material)
  // A donor who is genuinely not a parent of this child, by identity rather than by accession.
  const stranger = donors.find((d) => !trueIds.has(identityOf.get(d.gsm)!))
  if (stranger) measure(stranger.gsm, c.gsm, 'STRANGER', material)

  done += 1
  if (done % 10 === 0) console.log(`  measured ${done}/${children.length} children`, )
  // Children are read once; only donors are worth holding.
  cache.delete(c.gsm)
}

console.log('')
console.log(`children measured: ${done}, pairs: ${pairs.length}`)
console.log(`the gate in question: OPPOSITE_HOM_MAX = ${rel.OPPOSITE_HOM_MAX}, `
  + `IBD0_FLOOR_MAX = ${rel.IBD0_FLOOR_MAX}`)
console.log('')

const materials = [...new Set(pairs.map((p) => p.material))].sort()
const lo = (xs: number[]) => (xs.length ? Math.min(...xs) : NaN)
const hi = (xs: number[]) => (xs.length ? Math.max(...xs) : NaN)

/** Does `stat` put every TRUE pair below every STRANGER pair, within one material? */
function separability(name: string, pick: (p: Pair) => number) {
  console.log('')
  console.log(`=== ${name}`)
  console.log('material          TRUE n  range                REPLICATE            '
    + 'STRANGER n  range                gap        separable')
  for (const m of materials) {
    const inM = pairs.filter((p) => p.material === m)
    const t = inM.filter((p) => p.arm === 'TRUE').map(pick).filter(Number.isFinite)
    const s = inM.filter((p) => p.arm === 'STRANGER').map(pick).filter(Number.isFinite)
    const rp = inM.filter((p) => p.arm === 'REPLICATE').map(pick).filter(Number.isFinite)
    if (!t.length || !s.length) {
      console.log(`${m.padEnd(16)} ${String(t.length).padStart(6)}  `
        + `${t.length ? `${lo(t).toFixed(4)}-${hi(t).toFixed(4)}` : 'none'.padEnd(19)}`
        + '                     '
        + `${String(s.length).padStart(10)}  ${s.length ? '' : 'none'}   no comparison`)
      continue
    }
    const gap = lo(s) - hi(t)
    console.log(`${m.padEnd(16)} ${String(t.length).padStart(6)}  `
      + `${`${lo(t).toFixed(4)}-${hi(t).toFixed(4)}`.padEnd(19)}  `
      + `${(rp.length ? `${lo(rp).toFixed(4)}-${hi(rp).toFixed(4)} (n=${rp.length})` : 'none').padEnd(19)}  `
      + `${String(s.length).padStart(10)}  ${`${lo(s).toFixed(4)}-${hi(s).toFixed(4)}`.padEnd(19)}  `
      + `${gap >= 0 ? `+${gap.toFixed(4)}` : gap.toFixed(4)}   `
      + `${gap > 0 ? `YES, any cut in (${hi(t).toFixed(4)}, ${lo(s).toFixed(4)})` : 'NO, they overlap'}`)
  }
}

separability('POOLED OPPOSITE-HOMOZYGOTE RATE, the statistic that ships as the gate',
  (p) => p.opp)
separability('WINDOW FLOOR, the median 200-marker window. Never asked this question before',
  (p) => p.floor)
separability('WINDOW SPREAD, for completeness: it is the sibling statistic, not a parentage one',
  (p) => p.spread)

// ------------------------------------------------------------------ what the shipped verdict does
console.log('')
console.log('=== WHAT THE SHIPPED VERDICT ACTUALLY SAYS, per arm and material')
for (const m of materials) {
  for (const a of ['TRUE', 'REPLICATE', 'STRANGER'] as const) {
    const g = pairs.filter((p) => p.material === m && p.arm === a)
    if (!g.length) continue
    const counts = new Map<string, number>()
    for (const p of g) counts.set(p.verdict, (counts.get(p.verdict) ?? 0) + 1)
    console.log(`  ${m.padEnd(16)} ${a.padEnd(10)} n=${String(g.length).padStart(3)}  `
      + [...counts].map(([v, n]) => `${n}x ${v.slice(0, 38)}`).join('  |  '))
  }
}

console.log('')
console.log('HOW TO READ IT. The REPLICATE column is the control: it is the same person as TRUE on a')
console.log('different array, so it must sit with TRUE. If it sits with STRANGER then the stranger')
console.log('arm is measuring batch rather than relatedness and no conclusion here is safe. Where a')
console.log('statistic shows a positive gap within a material, a threshold exists for that material;')
console.log('where it does not, the overlap is real and the tool cannot gate on it whatever the')
console.log('constant is set to.')

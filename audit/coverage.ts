/**
 * Coverage of every abnormality class `om taxonomy` names, on real arrays.
 *
 * Two quantities per class, kept apart:
 *   exercised   the class was emitted on real material
 *   confirmed   an instance carries an answer from outside this tool
 *
 * Adult gamete donors are the negative control for the taxonomy findings as well as for
 * aneuploidy. Whole-chromosome uniparental disomy is rare in the population, triploidy is not
 * survivable, a donor genome is not chaotic and a donor is diploid, so those classes on a donor are
 * false calls. Segmental runs of homozygosity and copy-neutral LOH are reported for donors without
 * being scored, because the taxonomy records that they overlap the normal population.
 *
 * Classes the taxonomy declares unreachable are listed and not counted against coverage.
 *
 * Run: OM_TRIOS=<dir> [OM_MOLES=<dir>] [OM_BLASTO=<dir>] [OM_GSE=<dir>] [OM_CAP=n] \
 *        node --experimental-strip-types --max-old-space-size=8192 audit/coverage.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const TAXONOMY: { cls: string; label: string; reachable: boolean; truth: string }[] = [
  { cls: 'monosomy', label: 'Whole-chromosome monosomy', reachable: true, truth: 'targeted cuts; lab loci' },
  { cls: 'trisomy', label: 'Whole-chromosome trisomy', reachable: true, truth: 'GSE19247 karyotype; 3 of 13 moles' },
  { cls: 'segmental-deletion', label: 'Segmental deletion', reachable: true, truth: 'EYS and chr16 cuts; trio-resolved' },
  { cls: 'segmental-duplication', label: 'Segmental duplication', reachable: true, truth: 'none' },
  { cls: 'cnloh', label: 'Copy-neutral loss of heterozygosity', reachable: true, truth: 'none' },
  { cls: 'isodisomy', label: 'Uniparental isodisomy, whole chromosome', reachable: true, truth: 'none' },
  { cls: 'heterodisomy', label: 'Uniparental heterodisomy', reachable: false, truth: 'declared unreachable' },
  { cls: 'segmental-upd', label: 'Segmental uniparental disomy', reachable: true, truth: 'none' },
  { cls: 'triploidy', label: 'Triploidy', reachable: true, truth: 'none' },
  { cls: 'haploidy', label: 'Haploidy', reachable: true, truth: 'pronuclei; sperm' },
  { cls: 'complex', label: 'Complex or chaotic genome', reachable: true, truth: 'none per sample' },
  { cls: 'gamete-segmental', label: 'Segmental change introduced by the gamete', reachable: true, truth: 'needs a second unit' },
  { cls: 'reverse-segregation', label: 'Reverse segregation', reachable: false, truth: 'declared blind spot' },
  { cls: 'tandem-vs-inserted', label: 'Tandem versus inserted duplication', reachable: false, truth: 'declared platform limit' },
]
/** Classes a healthy adult cannot carry, so any call on a donor is false. */
const DONOR_FALSE = new Set(['monosomy', 'trisomy', 'isodisomy', 'triploidy', 'complex', 'haploidy'])

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

const DIR = process.env.OM_TRIOS
const donors = new Set<string>()
if (DIR && existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
  const head = lines[0].split(',')
  const gi = head.indexOf('gsm'); const ri = head.indexOf('role')
  for (const l of lines.slice(1)) {
    const f = l.split(',')
    if (f[ri] === 'parent') donors.add(f[gi])
  }
}
const DONOR = 'adult donors (control)'
const inputs: { path: string; set: string }[] = []
if (DIR && existsSync(DIR)) {
  for (const f of readdirSync(DIR).filter((x) => x.endsWith('.probes.gz'))) {
    inputs.push({ path: join(DIR, f), set: donors.has(f.replace('.probes.gz', '')) ? DONOR : 'GSE148488' })
  }
}
for (const [env, tag] of [['OM_MOLES', 'GSE117672 moles'], ['OM_BLASTO', 'GSE186407 blastomeres'],
  ['OM_GSE', 'GSE19247']] as const) {
  const d = process.env[env]
  if (!d || !existsSync(d)) continue
  for (const f of readdirSync(d).filter((x) => x.endsWith('.probes') || x.endsWith('.probes.txt.gz'))) {
    inputs.push({ path: join(d, f), set: tag })
  }
}
const use = inputs.slice(0, Number(process.env.OM_CAP ?? 600))

const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(use[0].path)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

interface Tally { arrays: Set<string>; instances: number; bySet: Map<string, Set<string>>; byZyg: Map<string, number> }
const seen = new Map<string, Tally>()
const rawCls = new Map<string, number>()
const driven = new Map<string, number>()
/** Distinct arrays per class, so the column cannot exceed the number of arrays driven. */
function bump(cls: string, id: string, set: string, zyg: string, n = 1) {
  const t = seen.get(cls) ?? { arrays: new Set<string>(), instances: 0, bySet: new Map(), byZyg: new Map() }
  if (!t.arrays.has(id)) t.byZyg.set(zyg, (t.byZyg.get(zyg) ?? 0) + 1)
  t.arrays.add(id)
  t.instances += n
  const s = t.bySet.get(set) ?? new Set<string>()
  s.add(id)
  t.bySet.set(set, s)
  seen.set(cls, t)
}
let threw = 0

for (const { path, set } of use) {
  try {
    const acc = score.emptyCollected(refParent as never, null)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(path)()) {
      if (!first) first = r.probesetId
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(path, byChrom as never, bafSums as never, first)
    const res = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: path, log: () => {},
    }) as Record<string, unknown>
    driven.set(set, (driven.get(set) ?? 0) + 1)
    const zyg = String(res.zygosity ?? '?')

    const chroms = (res.chroms ?? []) as { aneuploidy?: string }[]
    const loss = chroms.filter((c) => c.aneuploidy === 'loss').length
    const gain = chroms.filter((c) => c.aneuploidy === 'gain').length
    if (loss) bump('monosomy', path, set, zyg, loss)
    if (gain) bump('trisomy', path, set, zyg, gain)
    const segs = (res.segments ?? []) as { kind: string }[]
    const del = segs.filter((s) => /absence|copy-loss/.test(s.kind)).length
    const dup = segs.filter((s) => /gain/.test(s.kind)).length
    if (del) bump('segmental-deletion', path, set, zyg, del)
    if (dup) bump('segmental-duplication', path, set, zyg, dup)

    for (const f of (res.findings ?? []) as { cls: string; wholeChromosome?: boolean }[]) {
      const c = String(f.cls)
      rawCls.set(`${c}${f.wholeChromosome ? ' (whole)' : ''}`, (rawCls.get(`${c}${f.wholeChromosome ? ' (whole)' : ''}`) ?? 0) + 1)
      const lc = c.toLowerCase()
      if (/triploid/.test(lc)) bump('triploidy', path, set, zyg)
      else if (/complex|chaotic/.test(lc)) bump('complex', path, set, zyg)
      else if (/isodisom|upd|disomy/.test(lc)) bump(f.wholeChromosome ? 'isodisomy' : 'segmental-upd', path, set, zyg)
      else if (/loh|heterozygosity/.test(lc)) bump('cnloh', path, set, zyg)
    }
    if (zyg.startsWith('uniparental')) bump('haploidy', path, set, zyg)
  } catch {
    threw += 1
  }
}

const total = [...driven.values()].reduce((a, b) => a + b, 0)
console.log(`arrays driven: ${total}, threw: ${threw}`)
for (const [s, n] of driven) console.log(`  ${s}: ${n}`)
console.log('')
console.log('class                                     arrays  instances  independent answer')
for (const t of TAXONOMY) {
  const v = seen.get(t.cls)
  const mark = !t.reachable ? '[X]' : v ? '   ' : '[!]'
  console.log(`${mark} ${t.label.padEnd(40)}${String(v?.arrays.size ?? 0).padStart(6)}`
    + `${String(v?.instances ?? 0).padStart(11)}  ${t.truth}`)
}

console.log('')
console.log('=== WHERE EACH EXERCISED CLASS CAME FROM, and the zygosity of the arrays emitting it')
for (const t of TAXONOMY) {
  const v = seen.get(t.cls)
  if (!v) continue
  const sets = [...v.bySet].map(([s, ids]) => `${s} ${ids.size}/${driven.get(s) ?? 0}`).join('; ')
  const zygs = [...v.byZyg].map(([z, n]) => `${z} ${n}`).join(', ')
  console.log(`  ${t.label}`)
  console.log(`      sets: ${sets}`)
  console.log(`      zygosity: ${zygs}`)
}

console.log('')
console.log(`=== NEGATIVE CONTROL: ${driven.get(DONOR) ?? 0} adult donor arrays`)
let falseCalls = 0
for (const t of TAXONOMY) {
  const ids = seen.get(t.cls)?.bySet.get(DONOR)
  if (!ids?.size) continue
  const scored = DONOR_FALSE.has(t.cls)
  if (scored) falseCalls += ids.size
  console.log(`  ${t.label.padEnd(40)} ${ids.size} donor arrays  ${scored ? 'FALSE CALLS' : 'reported, not scored'}`)
}
console.log(`  donor arrays carrying a class no healthy adult can carry: ${falseCalls}`)

console.log('')
console.log('=== RAW FINDING CLASSES AS EMITTED')
for (const [c, n] of [...rawCls].sort((a, b) => b[1] - a[1])) console.log(`  ${String(n).padStart(6)}  ${c}`)

/**
 * Coriell cell lines in GSE18932 whose karyotypes are published, as positive controls.
 *
 * Two kinds of truth sit in the series beside its embryo samples:
 *   uniparental isodisomy   GM11496 is 46,XX with isodisomy of chromosome 7 and GM15603 is 46,XY
 *                           with isodisomy of chromosome 8. Whole-chromosome isodisomy on a stated
 *                           chromosome is a class the rest of the corpus has no outside answer for.
 *   mixtures                GM00323 and GM00321 combined in stated proportions. A two-person
 *                           mixture is the signature the integrity gates exist to catch.
 *
 * Line sex comes from Cellosaurus, not from the series: GM00323 (CVCL_7280) is male and GM00321
 * (CVCL_7279) is female. The series' "cell line karyotype" field reads 46,XX on the pure GM00323 array
 * and 46,XY on the pure GM00321 array, and its mixture records disagree with each other.
 *
 * The series publishes genotype calls only. Isodisomy leaves a chromosome with no heterozygous
 * calls and copy number unchanged, so here it is visible as heterozygosity and not as dosage, and it
 * cannot be told from monosomy without intensity. The arm scores whether the stated chromosome is
 * the one flagged, not which of the two it is.
 *
 * A refusal is scored apart from a miss. An array whose stage inference failed reports nothing and
 * says so; that is not the tool reading the chromosome as normal.
 *
 * Run: OM_CELLS=<dir with .probes and truth.json> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/cellline-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_CELLS as string
if (!DIR || !existsSync(join(DIR, 'truth.json'))) {
  console.log('cellline-truth: OM_CELLS with truth.json is required. NOT RUN.')
  process.exit(0)
}

type T = { gsm: string; title: string } & Record<string, string>
const truth: T[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(DIR, `${t.gsm}.probes`)))

const LINE_SEX = new Map([['GM00323', 'XY'], ['GM00321', 'XX']])
const upd = truth.flatMap((t) => {
  const c = /UPID\s*chr(\w+)/i.exec(t.karyotype ?? '')?.[1]
  return c ? [{ t, chrom: c }] : []
})
const mix = truth.flatMap((t) => {
  const m = /mixture\s*(\d+)%(GM\d+)_(\d+)%(GM\d+)/i.exec(t.title)
  return m ? [{ t, a: m[2], pa: +m[1], b: m[4], pb: +m[3] }] : []
})
const xyShare = (m: typeof mix[number]) =>
  (LINE_SEX.get(m.a) === 'XY' ? m.pa : 0) + (LINE_SEX.get(m.b) === 'XY' ? m.pb : 0)

console.log(`isodisomy arrays: ${upd.length}   mixture arrays: ${mix.length}`)
console.log(`line sexes (Cellosaurus): ${[...LINE_SEX].map(([l, s]) => `${l} ${s}`).join(', ')}`)

function rowsOf(path: string) {
  const ls = readFileSync(path, 'utf8').split('\n')
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
const isAutosome = (c: string) => /^\d+$/.test(c) && +c >= 1 && +c <= 22
const median = (xs: number[]) => { const s = [...xs].sort((a, b) => a - b); return s.length ? s[s.length >> 1] : NaN }

const pure = mix.find((m) => m.pa === 100 || m.pb === 100)
const refGsm = pure?.t.gsm ?? upd[0]?.t.gsm ?? truth[0].gsm
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR, `${refGsm}.probes`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

async function scoreOne(gsm: string) {
  const counts = new Map<string, [number, number]>()
  const acc = score.emptyCollected(refParent as never, null)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  let het = 0; let called = 0
  for (const r of rowsOf(join(DIR, `${gsm}.probes`))()) {
    if (!first) first = r.probesetId
    if (r.genotype !== 'NC' && isAutosome(r.chrom)) {
      const c = counts.get(r.chrom) ?? [0, 0]
      c[1] += 1; called += 1
      if (r.genotype === 'AB') { c[0] += 1; het += 1 }
      counts.set(r.chrom, c)
    }
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, refParent as never, null, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first) as Record<string, any>
  const allGates = ingest.gates(profile as never) as { name: string; verdict: string }[]
  const res = await score.scoreSample({
    acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
    sampleName: gsm, log: () => {},
  }) as Record<string, any>
  const hetMap = new Map<string, number>()
  for (const [c, [h, n]] of counts) if (n >= 200) hetMap.set(c, h / n)
  const whole = new Map<string, string[]>()
  const add = (c: string, w: string) => whole.set(c, [...(whole.get(c) ?? []), w])
  for (const c of (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]) {
    if (c.aneuploidy) add(c.chrom, `aneuploidy-${c.aneuploidy}`)
  }
  const genomeScope: string[] = []
  const segmental: string[] = []
  for (const f of (res.findings ?? []) as { cls: string; chrom: string; wholeChromosome?: boolean }[]) {
    if (f.chrom === 'genome') genomeScope.push(f.cls)
    else if (f.wholeChromosome) add(f.chrom, String(f.cls))
    else segmental.push(`${f.cls}:${f.chrom}`)
  }
  return {
    het: hetMap, genomeHet: called ? het / called : NaN,
    xRatio: profile.chrXHetRatio as number | null, sex: String(profile.sex),
    zyg: String(res.zygosity ?? '?'),
    stage: res.stage?.stage ?? '?',
    integrity: res.integrity?.level ?? '?',
    whole, genomeScope, segmental,
    gates: allGates.filter((g) => g.verdict === 'exclude' || g.verdict === 'marginal').map((g) => `${g.name}:${g.verdict}`),
    sexAmbiguous: allGates.some((g) => g.name === 'sex call' && g.verdict === 'report_only'),
  }
}

console.log('')
console.log('=== UNIPARENTAL ISODISOMY on a stated chromosome')
let flagged = 0; let accepted = 0; let leastHet = 0; let elsewhere = 0; let genomeN = 0; let segN = 0
for (const u of upd) {
  try {
    const s = await scoreOne(u.t.gsm)
    const hc = s.het.get(u.chrom)
    const rank = hc === undefined ? NaN : 1 + [...s.het.values()].filter((v) => v < hc).length
    const refused = s.stage === 'failed'
    const onStated = s.whole.has(u.chrom)
    const others = [...s.whole.keys()].filter((c) => c !== u.chrom)
    if (rank === 1) leastHet += 1
    if (!refused) {
      accepted += 1
      if (onStated) flagged += 1
      elsewhere += others.length
      genomeN += s.genomeScope.length
      segN += s.segmental.length
    }
    console.log(`  ${u.t.gsm} ${u.t['corielle catalog id'] ?? ''} chr${u.chrom}: het there `
      + `${hc === undefined ? '-' : hc.toFixed(4)} against array median ${median([...s.het.values()]).toFixed(4)}, `
      + `rank ${rank} of ${s.het.size}; stage ${s.stage}${refused ? ' (REFUSED, nothing reported)' : ''}; `
      + `flagged ${onStated ? s.whole.get(u.chrom)!.join(',') : 'no'}; other chromosomes whole ${others.length ? others.map((c) => `chr${c}`).join(',') : '-'}; `
      + `genome-scope ${s.genomeScope.join(',') || '-'}; segmental ${s.segmental.join(',') || '-'}; zyg ${s.zyg}`)
  } catch (e) {
    console.log(`  THREW ${u.t.gsm}: ${String((e as Error).message ?? e).slice(0, 80)}`)
  }
}
console.log(`  stated chromosome least heterozygous in the data: ${leastHet}/${upd.length}`)
console.log(`  arrays refused by the stage inference: ${upd.length - accepted}/${upd.length}`)
console.log(`  stated chromosome flagged whole-chromosome, accepted arrays: ${flagged}/${accepted}`)
console.log(`  on accepted arrays, false: whole-chromosome on another chromosome ${elsewhere}, `
  + `genome-scope ${genomeN}, segmental ${segN}`)

console.log('')
console.log('=== MIXTURES in stated proportions, ordered by the share of the XY line')
const tally = { pure: { n: 0, refused: 0, sex: 0, gate: 0, any: 0 }, mixed: { n: 0, refused: 0, sex: 0, gate: 0, any: 0 } }
for (const m of [...mix].sort((x, y) => xyShare(x) - xyShare(y))) {
  try {
    const s = await scoreOne(m.t.gsm)
    const t = m.pa === 100 || m.pb === 100 ? tally.pure : tally.mixed
    const refused = s.stage === 'failed'
    t.n += 1
    if (refused) t.refused += 1
    if (s.sexAmbiguous) t.sex += 1
    if (s.gates.length) t.gate += 1
    if (refused || s.sexAmbiguous || s.gates.length) t.any += 1
    console.log(`  ${m.t.title.padEnd(34)} XY ${String(xyShare(m)).padStart(3)}%  het ${s.genomeHet.toFixed(4)}  `
      + `xRatio ${s.xRatio?.toFixed(3) ?? '-'} ${s.sex.padEnd(9)} stage ${s.stage.padEnd(6)} `
      + `gates ${s.gates.join(' ') || '-'}`)
  } catch (e) {
    console.log(`  THREW ${m.t.gsm}: ${String((e as Error).message ?? e).slice(0, 80)}`)
  }
}
for (const [k, t] of Object.entries(tally)) {
  console.log(`  ${k.padEnd(5)} n ${t.n}: stage refused ${t.refused}, sex call ambiguous ${t.sex}, `
    + `exclude or marginal gate ${t.gate}, any of the three ${t.any}`)
}

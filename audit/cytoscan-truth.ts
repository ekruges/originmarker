/**
 * Bulk arrays whose abnormality the submitters state per sample, on Affymetrix CytoScan 750K.
 *
 * GSE207887 is miscarriage material, chorionic villus and fetal skin, each sample carrying a
 * `chromosomal abnormality` field: whole-chromosome trisomies, 45,X, triploidy, and a set stated
 * free of clinical variation. GSE163799 is fetal material whose `diagnosis` field names a
 * microduplication or microdeletion with its cytogenetic band.
 *
 * Two things here exist nowhere else in the corpus. The DNA is BULK, so every sensitivity figure
 * measured on amplified single cells gets a companion measured without amplification. And triploidy
 * and segmental duplication are classes the taxonomy names that no real array had ever produced.
 *
 * Positions are the submitters' own, annotated against GRCh37, and the band coordinates come from
 * the UCSC cytoBand table rather than from this tool. No parental array exists in either series, so
 * nothing about parent of origin is scored: a normal sample fills the parental slot only because
 * the reader needs one, and the channels scored here read the sample alone.
 *
 * Run: OM_CYTO=<dir with .probes and truth.json> OM_BANDS=<cytoBand.txt.gz> \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/cytoscan-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_CYTO as string
const BANDS = process.env.OM_BANDS as string
if (!DIR || !existsSync(join(DIR, 'truth.json'))) {
  console.log('cytoscan-truth: OM_CYTO with truth.json is required. NOT RUN.')
  process.exit(0)
}
type T = { gsm: string; title: string; label: string }
const truth: T[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(DIR, `${t.gsm}.probes`)))

/** chrom -> band -> [start, end], from the UCSC table. */
const bandTable = new Map<string, Map<string, [number, number]>>()
if (BANDS && existsSync(BANDS)) {
  for (const line of gunzipSync(readFileSync(BANDS)).toString('utf8').split('\n')) {
    const p = line.split('\t')
    if (p.length < 4) continue
    const c = p[0].replace('chr', '')
    if (!bandTable.has(c)) bandTable.set(c, new Map())
    bandTable.get(c)!.set(p[3], [Number(p[1]), Number(p[2])])
  }
}
/** "16p13.13p13.12" -> the span from the first band's start to the last band's end. */
function bandSpan(chrom: string, label: string): [number, number] | null {
  const table = bandTable.get(chrom)
  if (!table) return null
  const parts = [...label.matchAll(/[pq][\d.]+/g)].map((m) => m[0])
  const hits = parts.map((b) => table.get(b)).filter(Boolean) as [number, number][]
  if (!hits.length) return null
  return [Math.min(...hits.map((h) => h[0])), Math.max(...hits.map((h) => h[1]))]
}

type Want =
  | { kind: 'gain' | 'loss'; chrom: string }
  /** A stated pathogenic CNV with neither a chromosome nor a size: presence only. */
  | { kind: 'cnv-present' }
  | { kind: 'triploidy' }
  | { kind: 'loh'; chrom: string }
  | { kind: 'seg-gain' | 'seg-loss'; chrom: string; span: [number, number] | null; band: string }
  | { kind: 'negative' }
  | { kind: 'unscored' }

function want(label: string): Want {
  const l = (label || '').trim()
  let m = /^(\d{1,2}) Trisomy$/.exec(l)
  if (m) return { kind: 'gain', chrom: m[1] }
  if (/^45,X$/.test(l)) return { kind: 'loss', chrom: 'X' }
  if (/^Triploidy$/i.test(l)) return { kind: 'triploidy' }
  if (/^No clinical Variations$/i.test(l)) return { kind: 'negative' }
  m = /Heterozygosity is absent on chromosome (\d{1,2}|X|Y)/i.exec(l)
  if (m) return { kind: 'loh', chrom: m[1] }
  // The series states that these arrays carry a pathogenic imbalance and nothing about where it
  // is or how big, so the only thing scoreable is whether the tool reports an imbalance at all.
  if (/^Chromosomal deletions\/duplications$/i.test(l)) return { kind: 'cnv-present' }
  m = /micro(duplication|deletion) of ([\dXY]{1,2})([pq][\d.pq]*)/i.exec(l)
    ?? /^([\dXY]{1,2})()([pq][\d.pq]*) micro(duplication)/i.exec(l)
  if (m) {
    const dup = /duplication/i.test(l)
    const chrom = m[2] || m[1]
    const band = m[3] ?? ''
    return { kind: dup ? 'seg-gain' : 'seg-loss', chrom, span: bandSpan(chrom, band), band }
  }
  return { kind: 'unscored' }
}

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
const wilson = (k: number, n: number): [number, number] => {
  if (!n) return [NaN, NaN]
  const z = 1.96; const p = k / n; const d = 1 + (z * z) / n
  const c = p + (z * z) / (2 * n)
  const h = z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))
  return [Math.max(0, (c - h) / d), Math.min(1, (c + h) / d)]
}
const pct = (x: number) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : '-')
const band = (k: number, n: number) => `${k}/${n} ${pct(k / Math.max(n, 1))} [${pct(wilson(k, n)[0])}, ${pct(wilson(k, n)[1])}]`

// The parental slot takes a sample stated free of variation. A series without one, such as a
// cohort assembled around a single finding, names any array of the same platform through OM_CYTO_REF.
const negative = truth.find((t) => want(t.label).kind === 'negative')
const refPath = process.env.OM_CYTO_REF
  ?? (negative ? join(DIR, `${negative.gsm}.probes`) : '')
if (!refPath || !existsSync(refPath)) {
  console.log('no array to fill the parental slot: set OM_CYTO_REF. NOT RUN.')
  process.exit(0)
}
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(refPath)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()
console.log(`parental slot: ${refPath.split('/').pop()}, heterozygosity `
  + `${(refParent.heterozygosity as number).toFixed(4)}`)

interface Score { n: number; hit: number; refused: number }
const byClass = new Map<string, Score>()
const bump = (k: string, hit: boolean, refused: boolean) => {
  const s = byClass.get(k) ?? { n: 0, hit: 0, refused: 0 }
  if (refused) s.refused += 1
  else { s.n += 1; if (hit) s.hit += 1 }
  byClass.set(k, s)
}
let negN = 0; let negClean = 0; let negRefused = 0; let threw = 0
const negFalse = new Map<string, number>()
/** One row per array: what the sex chromosomes read against what the series states. */
const sexRows: { gsm: string; label: string; kind: string; copies: number | null;
  constitution: string }[] = []
console.log('')
for (const t of truth) {
  const w = want(t.label)
  if (w.kind === 'unscored') continue
  try {
    const acc = score.emptyCollected(refParent as never, null)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    for (const r of rowsOf(join(DIR, `${t.gsm}.probes`))()) {
      if (!first) first = r.probesetId
      ingest.accumulate(r as never, byChrom as never)
      ingest.accumulateBaf(r as never, bafSums as never)
      score.collectRow(r as never, refParent as never, null, acc as never)
    }
    const profile = ingest.finishProfile(t.gsm, byChrom as never, bafSums as never, first) as Record<string, any>
    const res = await score.scoreSample({
      acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
      sampleName: t.gsm, log: () => {},
    }) as Record<string, any>
    const refused = res.stage?.stage === 'failed'
    const chroms = (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
    const segs = (res.segments ?? []) as { chrom: string; kind: string; startBp: number; endBp: number }[]
    const finds = (res.findings ?? []) as { cls: string; chrom: string }[]
    const whole = new Map(chroms.filter((c) => c.aneuploidy).map((c) => [c.chrom, c.aneuploidy as string]))
    const overlaps = (s: { startBp: number; endBp: number }, span: [number, number] | null) =>
      !span || (s.startBp <= span[1] && s.endBp >= span[0])

    // What the sex chromosomes were found to be. A stated 45,X is a LOSS of a chromosome that the
    // whole-chromosome test deliberately does not run on, so it is scored where it is decided.
    const sx = res.sexChromosome as
      { xCopies: number | null; constitution: string; why: string } | undefined
    sexRows.push({ gsm: t.gsm, label: String(t.label), kind: w.kind,
      copies: sx?.xCopies ?? null, constitution: sx?.constitution ?? 'none' })

    let hit = false
    if (w.kind === 'cnv-present') {
      hit = whole.size > 0 || segs.length > 0 || sx?.constitution === '45,X'
    } else if (w.kind === 'loss' && w.chrom === 'X') hit = sx?.constitution === '45,X'
    else if (w.kind === 'gain' || w.kind === 'loss') hit = whole.get(w.chrom) === w.kind
    else if (w.kind === 'triploidy') hit = finds.some((f) => f.cls === 'triploidy')
    else if (w.kind === 'loh') hit = finds.some((f) => f.chrom === w.chrom && /loh|isodisomy|upd/i.test(f.cls))
    else if (w.kind === 'seg-gain') {
      hit = whole.get(w.chrom) === 'gain'
        || segs.some((s) => s.chrom === w.chrom && /gain/i.test(s.kind) && overlaps(s, w.span))
    } else if (w.kind === 'seg-loss') {
      hit = whole.get(w.chrom) === 'loss'
        || segs.some((s) => s.chrom === w.chrom && /absence|copy-loss|loss/i.test(s.kind) && overlaps(s, w.span))
    }

    if (w.kind === 'negative') {
      if (refused) negRefused += 1
      else {
        negN += 1
        const bad = [...whole.keys()].map((c) => `whole chr${c}`)
          .concat(finds.filter((f) => f.cls === 'triploidy').map(() => 'triploidy'))
          .concat(segs.length ? [`${segs.length} segment(s)`] : [])
        if (!bad.length) negClean += 1
        for (const b of bad) negFalse.set(b.replace(/\d+/g, 'N'), (negFalse.get(b.replace(/\d+/g, 'N')) ?? 0) + 1)
      }
    } else {
      bump(w.kind === 'seg-gain' || w.kind === 'seg-loss' ? `${w.kind} ${w.band}`
        : `${w.kind}${'chrom' in w ? ` chr${w.chrom}` : ''}`,
        hit, refused)
      const shape = refused ? 'REFUSED'
        : `sex ${String(profile.sex).padEnd(9)} ${[...whole].map(([c, a]) => `chr${c} ${a}`).join(',') || '-'}`
          + `${finds.length ? ` | ${finds.map((f) => f.cls).slice(0, 3).join(',')}` : ''}`
          + `${segs.length ? ` | ${segs.length} seg` : ''}`
      console.log(`  ${t.gsm} ${String(t.label).slice(0, 34).padEnd(35)} ${hit ? 'HIT ' : '    '} ${shape}`)
    }
  } catch (e) {
    threw += 1
    console.log(`  THREW ${t.gsm}: ${String((e as Error).message ?? e).slice(0, 70)}`)
  }
}

console.log('')
console.log('=== WHAT THE SUBMITTERS STATE, AND WHETHER THE TOOL SAYS IT, ON BULK DNA')
const groups = new Map<string, Score>()
for (const [k, s] of byClass) {
  const g = k.startsWith('seg-') ? k.split(' ')[0] : k.split(' ')[0]
  const cur = groups.get(g) ?? { n: 0, hit: 0, refused: 0 }
  cur.n += s.n; cur.hit += s.hit; cur.refused += s.refused
  groups.set(g, cur)
}
for (const [k, s] of [...byClass].sort()) {
  console.log(`  ${k.padEnd(28)} ${band(s.hit, s.n)}${s.refused ? `   refused ${s.refused}` : ''}`)
}
console.log('')
for (const [k, s] of [...groups].sort()) {
  console.log(`  ALL ${k.padEnd(24)} ${band(s.hit, s.n)}${s.refused ? `   refused ${s.refused}` : ''}`)
}
console.log('')
console.log('=== NEGATIVES, stated free of clinical variation')
console.log(`  clean end to end: ${band(negClean, negN)}   refused ${negRefused}`)
for (const [k, n] of [...negFalse].sort((a, b) => b[1] - a[1])) console.log(`    ${n} x ${k}`)
console.log('')
console.log(`arrays that threw: ${threw}`)

console.log('')
console.log('=== SEX CHROMOSOMES. One X or two, and what the tool says that is.')
const sexBy = new Map<string, { n: number; copies: Map<number | null, number>;
  calls: Map<string, number> }>()
for (const r of sexRows) {
  const k = /^45,X$/.test(r.label) ? '45,X (stated)'
    : /No clinical Variations/i.test(r.label) ? 'no clinical variation'
      : /Triploidy/i.test(r.label) ? 'triploidy (stated)' : 'other stated abnormality'
  const g = sexBy.get(k) ?? { n: 0, copies: new Map(), calls: new Map() }
  g.n += 1
  g.copies.set(r.copies, (g.copies.get(r.copies) ?? 0) + 1)
  g.calls.set(r.constitution, (g.calls.get(r.constitution) ?? 0) + 1)
  sexBy.set(k, g)
}
for (const [k, g] of [...sexBy].sort()) {
  const copies = [...g.copies].sort().map(([c, n]) => `${c ?? 'not measured'}:${n}`).join('  ')
  const calls = [...g.calls].sort((a, b) => b[1] - a[1]).map(([c, n]) => `${c} ${n}`).join(', ')
  console.log(`  ${k.padEnd(26)} n ${String(g.n).padStart(3)}   X copies ${copies}`)
  console.log(`  ${''.padEnd(26)}         called   ${calls}`)
}

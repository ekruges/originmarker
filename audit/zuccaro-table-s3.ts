/**
 * Whole-chromosome parent of origin on single cells, against a published method's own hand calls.
 *
 * Zuccaro 2020 (Cell, doi 10.1016/j.cell.2020.10.025) karyotyped every arrayed cell by eye from
 * BAF and log R and recorded, per chromosome, whether it is normal, which parent's copy a monosomy
 * lost, or whether the chromosome is absent altogether. Those calls are audit/zuccaro-table-s3.tsv.
 * This runs the shipped pipeline on the same arrays with both parents loaded and scores the
 * allele-count detector (ChromResult.lostParent) and the Mendelian channel against them.
 *
 * TWO ARMS, as in trio-calibration.ts. The true parents, and then a mother who is a different egg
 * donor. A measurement of parent of origin must collapse in the second arm; the father cannot be
 * swapped because the series has one sperm donor.
 *
 * WHAT A DISAGREEMENT IS NOT. The hand calls are the comparison, not the truth: both readings come
 * from the same amplified array and neither was confirmed by dissection. A chromosome the alleles
 * read as missing one parent at normal intensity is reported here exactly as the page reports it,
 * a parental absence, and scored against a hand call of monosomy on the parent alone.
 *
 * Egg donor D (GSM4472408) is heterozygous at 0.3455 on an 82 percent call rate and is excluded as
 * a parental genotype by trio-calibration.ts. It is the mother of 18 of these 31 cells, so results
 * are reported with and without her.
 *
 * Run: OM_TRIOS=<dir with GSM*.probes.gz and trio_manifest_full.csv> [OM_OUT=out.json] \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/zuccaro-table-s3.ts
 */
import { readFileSync, existsSync, readdirSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const defects = await import(`${W}defects.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('zuccaro-table-s3: OM_TRIOS is not a readable directory. NOT RUN.')
  process.exit(0)
}
const EXCLUDED_MOTHER = 'GSM4472408'

const fileOf = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const m = /^(GSM\d+)/.exec(f)
  if (m && (f.endsWith('.gz') || f.endsWith('.probes') || f.endsWith('.txt'))) fileOf.set(m[1], join(DIR, f))
}

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
const parentsOf = new Map<string, { father: string; mother: string }>()
const mothers = new Set<string>()
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const father = at(r, 'father_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  const mother = at(r, 'mother_gsms').split(/[; ]+/).filter((g) => fileOf.has(g))[0]
  if (father && mother) parentsOf.set(at(r, 'gsm'), { father, mother })
  if (mother) mothers.add(mother)
}

interface Truth { gsm: string; row: string; chrom: string; manual: string }
const truth: Truth[] = readFileSync(new URL('./zuccaro-table-s3.tsv', import.meta.url), 'utf8')
  .split('\n').filter((l) => l && !l.startsWith('#') && !l.startsWith('gsm'))
  .map((l) => { const [gsm, row, , chrom, manual] = l.split('\t'); return { gsm, row, chrom, manual } })
const byGsm = new Map<string, Truth[]>()
for (const t of truth) byGsm.set(t.gsm, [...(byGsm.get(t.gsm) ?? []), t])

function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r }
  }
}
const parentCache = new Map<string, unknown>()
function parentIndex(gsm: string) {
  if (parentCache.has(gsm)) return parentCache.get(gsm)
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  parentCache.set(gsm, idx)
  return idx
}
async function scoreWith(gsm: string, father: string, mother: string) {
  const pat = parentIndex(father); const mat = parentIndex(mother)
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  return score.scoreSample({ acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {} }) as Promise<Record<string, unknown>>
}

/** The parent the hand call says is lost, or the class when it is not a parent-named monosomy. */
const manualClass = (m: string): { kind: 'intact' | 'monosomy' | 'other'; parent?: 'paternal' | 'maternal' } => {
  if (m === 'normal') return { kind: 'intact' }
  const mono = /^monosomy \d+ \((pat|mat)\)/.exec(m)
  if (mono) return { kind: 'monosomy', parent: mono[1] === 'pat' ? 'paternal' : 'maternal' }
  return { kind: 'other' }
}

interface Scored {
  gsm: string; mother: string; arm: 'true' | 'wrong'; stage: string; zygosity: string
  chrom: string; manual: string; lostParent?: string; aneuploidy?: string; mendel?: string | null
  onlyFather?: number; onlyMother?: number; ab?: number; informative?: number
  verdict: string
}
const scored: Scored[] = []
const samples = [...byGsm.keys()].filter((g) => fileOf.has(g) && parentsOf.has(g))
console.log(`arrays with hand calls and both parents on disk: ${samples.length} of ${byGsm.size}`)
for (const gsm of samples) {
  const p = parentsOf.get(gsm)!
  const wrongMother = [...mothers].find((m) => m !== p.mother && m !== EXCLUDED_MOTHER)!
  for (const [arm, mother] of [['true', p.mother], ['wrong', wrongMother]] as const) {
    const res = await scoreWith(gsm, p.father, mother)
    const chroms = (res.chroms ?? []) as { chrom: string; lostParent?: string; aneuploidy?: string; contribution?: { fraction: number; onlyA: number; onlyB: number; informative: number } }[]
    const one = (res.oneParent ?? []) as { where: string; verdict: string; parent?: 'paternal' | 'maternal' | null; twoParents?: boolean }[]
    for (const t of byGsm.get(gsm)!) {
      const c = chroms.find((x) => x.chrom === t.chrom)
      const row = one.find((o) => o.where === `chr${t.chrom}`)
      const mendel = row ? defects.mendelParent(row, 'paternal', c?.aneuploidy === 'gain') : null
      const m = manualClass(t.manual)
      let verdict: string
      if (m.kind === 'other') verdict = 'not scored'
      else if (m.kind === 'intact') verdict = c?.lostParent ? 'EXTRA' : 'agree intact'
      else if (!c?.lostParent) verdict = 'MISSED'
      else verdict = c.lostParent === m.parent ? 'agree parent' : 'WRONG parent'
      const k = c?.contribution
      scored.push({
        gsm, mother, arm, stage: String((res.stage as { stage?: string })?.stage ?? ''),
        zygosity: String(res.zygosity ?? ''), chrom: t.chrom, manual: t.manual,
        lostParent: c?.lostParent, aneuploidy: c?.aneuploidy, mendel,
        onlyFather: k?.onlyA, onlyMother: k?.onlyB, ab: k?.fraction, informative: k?.informative, verdict,
      })
    }
  }
  process.stdout.write('.')
}
console.log()

const tally = (rows: Scored[]) => {
  const out: Record<string, number> = {}
  for (const r of rows) out[r.verdict] = (out[r.verdict] ?? 0) + 1
  return out
}
const fmt = (o: Record<string, number>) => Object.entries(o).sort().map(([k, v]) => `${k} ${v}`).join(', ')
for (const arm of ['true', 'wrong'] as const) {
  const rows = scored.filter((r) => r.arm === arm)
  console.log(`\n${arm === 'true' ? 'TRUE PARENTS' : 'WRONG MOTHER'} (${new Set(rows.map((r) => r.gsm)).size} arrays, ${rows.length} chromosomes)`)
  console.log(`  all mothers:            ${fmt(tally(rows))}`)
  console.log(`  without egg donor D:    ${fmt(tally(rows.filter((r) => r.mother !== EXCLUDED_MOTHER && parentsOf.get(r.gsm)!.mother !== EXCLUDED_MOTHER)))}`)
  const named = rows.filter((r) => r.lostParent).length
  console.log(`  chromosomes with a parent named from the alleles: ${named}`)
  if (arm === 'true') {
    const agreeMendel = rows.filter((r) => r.verdict === 'agree parent' && r.mendel === r.lostParent).length
    const refusedMendel = rows.filter((r) => r.verdict === 'agree parent' && !r.mendel).length
    console.log(`  of the agreed parents, the Mendelian channel agrees on ${agreeMendel} and declines on ${refusedMendel}`)
  }
}
console.log('\nDISAGREEMENTS, true parents:')
for (const r of scored.filter((r) => r.arm === 'true' && /EXTRA|MISSED|WRONG/.test(r.verdict))) {
  console.log(`  ${r.verdict.padEnd(13)} ${r.gsm} chr${r.chrom.padStart(2)}  hand: ${r.manual.padEnd(20)} alleles: ${r.lostParent ?? '-'}`
    + `  AB ${((r.ab ?? NaN) * 100).toFixed(1)}%  father-only ${((r.onlyFather ?? NaN) * 100).toFixed(1)}%  mother-only ${((r.onlyMother ?? NaN) * 100).toFixed(1)}%`
    + `  n ${r.informative ?? '-'}  intensity ${r.aneuploidy ?? 'none'}  mendel ${r.mendel ?? 'declined'}  stage ${r.stage}`)
}
const stages = tally(scored.filter((r) => r.arm === 'true').map((r) => ({ ...r, verdict: `${r.stage}/${r.zygosity}` })))
console.log(`\nstage/zygosity of the scored arrays (per chromosome): ${fmt(stages)}`)
if (process.env.OM_OUT) writeFileSync(process.env.OM_OUT, JSON.stringify(scored, null, 1))

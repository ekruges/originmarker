/**
 * Detection on GSE19247 (Vanneste et al.), a public series from another laboratory on another
 * manufacturer's chemistry, where the answer is a karyotype. Both halves of the question are in
 * the one series:
 *
 *   Coriell family 1990    lymphoblasts of a karyotype-confirmed TRISOMY 21 individual
 *   Coriell family 1463    lymphoblasts of a euploid reference pedigree
 *   Coriell family 1423    lymphoblasts of a euploid reference pedigree
 *   blood, family 231      ordinary diploid somatic tissue
 *   sperm, family 231      one parental complement, so every chromosome is at one copy at once
 *
 * Truth: the GEO `source_name` field and the cell lines' published karyotypes. Nothing in it comes
 * from this tool, from the Egli series, or from any threshold measured here.
 *
 * The euploid lines are the specificity set: same material, laboratory, chips and preparation as
 * the trisomies, so every aneuploidy call on them is a false positive. Chromosome 21 on the trisomy
 * line is read only beside that set, since calling 21 everywhere would score perfectly on the
 * trisomies alone. A throw or a structurally invalid event is a failure whatever the rates say.
 *
 * Sperm are excluded from both scored arms and reported apart: a haploid cell is at one copy on
 * every chromosome, so no chromosome differs from the rest. Cleavage-stage embryos are reported
 * and not scored: their cells are aneuploid at a high, per-cell unknown rate, so no array has its
 * own answer.
 *
 * Run: OM_GSE=<converted dir with truth.json> node --experimental-strip-types \
 *        --max-old-space-size=6144 audit/gse19247-detection.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const GSE = process.env.OM_GSE
if (!GSE || !existsSync(join(GSE, 'truth.json'))) {
  console.log('gse19247-detection: OM_GSE must hold .probes files and truth.json. NOT RUN.')
  process.exit(0)
}

/** What each group's karyotype is, and whether an aneuploidy call on it is right or wrong. */
const GROUPS: {
  match: RegExp; label: string; expect: string | null; scored: boolean
}[] = [
  { match: /family 1990/i, label: 'TRISOMY 21, karyotype-confirmed', expect: '21', scored: true },
  { match: /family 1463/i, label: 'euploid Coriell 1463', expect: null, scored: true },
  { match: /family 1423/i, label: 'euploid Coriell 1423', expect: null, scored: true },
  { match: /blood cell/i, label: 'blood, diploid somatic', expect: null, scored: true },
  { match: /sperm/i, label: 'sperm, HAPLOID (not scored)', expect: null, scored: false },
  { match: /cleavage stage/i, label: 'cleavage embryos (not scored)', expect: null, scored: false },
]

interface T { gsm: string; title: string; src: string }
const truth: T[] = JSON.parse(readFileSync(join(GSE, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(GSE, `${t.gsm}.probes`)))

function rowsOf(path: string) {
  const ls = readFileSync(path, 'utf8').split('\n')
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
 * A FIXED ARRAY IN THE PARENTAL SLOT, taken from this same series so nothing crosses studies.
 * Detection reads the call rate and the intensity, neither of which consults the parental array.
 */
const refGsm = truth.find((t) => /blood cell/i.test(t.src))?.gsm ?? truth[0].gsm
const refParent = (() => {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(GSE, `${refGsm}.probes`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never,
    ingest.finishProfile('ref', byChrom as never, bafSums as never, first).build.build)
})()

const isAutosome = (c: string) => /^\d+$/.test(c)
interface Res {
  gsm: string; group: string; stage: string; callRate: number
  aneuChroms: string[]; segChroms: string[]; junk: string[]; error?: string
}
const out: Res[] = []
const CAP = Number(process.env.OM_PER_GROUP ?? 60)

const bucket = new Map<string, T[]>()
for (const t of truth) {
  const g = GROUPS.find((x) => x.match.test(t.src))
  if (!g) continue
  const cur = bucket.get(g.label) ?? []
  if (cur.length < CAP) bucket.set(g.label, [...cur, t])
}
console.log(`converted arrays available: ${truth.length}`)
for (const [g, ts] of bucket) console.log(`  ${g}: ${ts.length} (cap ${CAP})`)
console.log(`parental slot filled with ${refGsm} from this same series`)
console.log('')

for (const [group, ts] of bucket) {
  for (const t of ts) {
    const junk: string[] = []
    try {
      const maxPos = new Map<string, number>()
      const acc = score.emptyCollected(refParent as never, null)
      const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
      for (const r of rowsOf(join(GSE, `${t.gsm}.probes`))()) {
        if (!first) first = r.probesetId
        if (Number.isFinite(r.pos)) {
          maxPos.set(r.chrom, Math.max(maxPos.get(r.chrom) ?? 0, r.pos))
        }
        ingest.accumulate(r as never, byChrom as never)
        ingest.accumulateBaf(r as never, bafSums as never)
        score.collectRow(r as never, refParent as never, null, acc as never)
      }
      const profile = ingest.finishProfile(t.gsm, byChrom as never, bafSums as never, first)
      const res = await score.scoreSample({
        acc, profile, pat: refParent, mat: null, soloRole: 'paternal', sibs: [],
        sampleName: t.gsm, log: () => {},
      }) as Record<string, unknown>
      const chroms = (res.chroms ?? []) as { chrom: string; aneuploidy?: string }[]
      const segs = (res.segments ?? []) as {
        chrom: string; startBp: number; endBp: number; markers: number
        refined?: { startBp: number; endBp: number }
      }[]
      // Junk: the same structural checks, on a platform none of them were written against.
      for (const s of segs) {
        const st = s.refined?.startBp ?? s.startBp
        const en = s.refined?.endBp ?? s.endBp
        const top = maxPos.get(s.chrom)
        if (!Number.isFinite(st) || !Number.isFinite(en)) junk.push(`chr${s.chrom}: coordinate NaN`)
        else if (!(st < en)) junk.push(`chr${s.chrom}: start ${st} not before end ${en}`)
        if (top === undefined) junk.push(`chr${s.chrom}: no markers on that chromosome`)
        else if (Number.isFinite(en) && en > top * 1.05) {
          junk.push(`chr${s.chrom}: ends ${en} past last marker ${top}`)
        }
        if (!Number.isFinite(s.markers) || s.markers <= 0) junk.push(`chr${s.chrom}: 0 markers`)
      }
      out.push({
        gsm: t.gsm, group,
        stage: (res.stage as { stage?: string } | undefined)?.stage ?? 'none',
        callRate: profile.callRate,
        aneuChroms: chroms.filter((c) => c.aneuploidy && isAutosome(c.chrom)).map((c) => c.chrom),
        segChroms: segs.filter((s) => isAutosome(s.chrom)).map((s) => s.chrom),
        junk,
      })
    } catch (e) {
      out.push({
        gsm: t.gsm, group, stage: 'error', callRate: NaN, aneuChroms: [], segChroms: [], junk,
        error: String((e as Error).message ?? e).slice(0, 80),
      })
    }
  }
  const done = out.filter((r) => r.group === group).length
  console.log(`  scored ${done} of ${group}`)
}

console.log('')
console.log('=== SENSITIVITY AND SPECIFICITY, on karyotypes this tool did not produce')
console.log('group                              n  usable  chr21 called  OTHER autosomes called  crashed')
let falsePositives = 0
let tpr = { hit: 0, n: 0 }
for (const g of GROUPS) {
  const rows = out.filter((r) => r.group === g.label)
  if (!rows.length) continue
  const crashed = rows.filter((r) => r.error).length
  const usable = rows.filter((r) => !r.error && r.stage !== 'failed')
  const c21 = usable.filter((r) => r.aneuChroms.includes('21')).length
  const otherPer = usable.length
    ? usable.reduce((a, r) => a + r.aneuChroms.filter((c) => c !== '21').length, 0) / usable.length
    : NaN
  console.log(`${g.label.padEnd(34)} ${String(rows.length).padStart(2)}  ${String(usable.length).padStart(6)}  `
    + `${`${c21}/${usable.length}`.padStart(12)}  ${otherPer.toFixed(3).padStart(22)}  ${crashed}`)
  if (!g.scored) continue
  if (g.expect === '21') { tpr.hit += c21; tpr.n += usable.length } else {
    // A euploid line has no aneuploidy at all, so EVERY autosomal call on it is a false positive.
    falsePositives += usable.reduce((a, r) => a + r.aneuChroms.length, 0)
  }
}

const euploid = out.filter((r) => GROUPS.some((g) => g.label === r.group && g.scored && !g.expect))
  .filter((r) => !r.error && r.stage !== 'failed')
console.log('')
console.log(`TRUE POSITIVE RATE, chr21 on karyotype-confirmed trisomy 21: ${tpr.hit}/${tpr.n}`
  + `${tpr.n ? ` = ${(tpr.hit / tpr.n).toFixed(4)}` : ''}`)
console.log(`FALSE POSITIVES on euploid material: ${falsePositives} autosomal calls over `
  + `${euploid.length} arrays and ${euploid.length * 22} autosome observations`)
console.log(`  per-autosome false positive rate: `
  + `${(falsePositives / Math.max(euploid.length * 22, 1)).toFixed(5)}`)

console.log('')
console.log('=== JUNK, on a platform none of these checks were written against')
const allJunk = out.flatMap((r) => r.junk.map((j) => `${r.gsm}: ${j}`))
const errs = out.filter((r) => r.error)
console.log(`  arrays scored: ${out.length}, threw: ${errs.length}, invalid events: ${allJunk.length}`)
for (const e of errs.slice(0, 8)) console.log(`    ${e.gsm}: ${e.error}`)
for (const j of allJunk.slice(0, 15)) console.log(`    ${j}`)

console.log('')
const clean = errs.length === 0 && allJunk.length === 0
console.log(clean
  ? 'gse19247-detection: no array threw and every emitted event is structurally valid.'
  : `gse19247-detection: ${errs.length} throws, ${allJunk.length} invalid events.`)

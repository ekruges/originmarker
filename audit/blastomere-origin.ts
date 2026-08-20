/**
 * Can a parent be named on a BIPARENTAL sample from allele absence, and how often is it right.
 *
 * WHY THIS IS A DIFFERENT QUESTION FROM EVERY MEASUREMENT BEFORE IT. A blastomere carries both
 * parental genomes. At a marker where the loaded parent is HOMOZYGOUS, a sample carrying both
 * copies reads heterozygous exactly when the other parent transmitted the alternate allele. Remove
 * the loaded parent's copy and those markers read homozygous for the OTHER parent's allele; remove
 * the other parent's copy and they read homozygous for the LOADED parent's. The two are opposite
 * and the difference is Mendelian rather than statistical: dropout removes an allele, it cannot
 * invent one.
 *
 * That is why this does not need the dosage channel's detection floors. Those floors exist because
 * a mosaic fraction shifts a mean; this asks whether an allele is there at all.
 *
 * THE POSITIVE CONTROL. Truth is constructed rather than inferred: a real biparental array has one
 * parent's contribution removed across a chromosome, by replacing the sample's genotype at
 * parent-homozygous markers with the genotype that survives when that parent's copy is gone. Both
 * directions are built from the same array, so nothing but the removed parent differs between them.
 */
import { readFileSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { homedir } from 'node:os'
import { callOneParentOrigin } from '../web/src/oneParentOrigin.ts'
import { inferStage } from '../web/src/stage.ts'

const DIR = process.argv[2] ?? join(homedir(), 'Downloads', 'probes')
const AUTO = new Set(Array.from({ length: 22 }, (_, i) => String(i + 1)))

type Row = { chrom: string; parent: 'AA' | 'BB'; gt: string }

/** The loaded parent's homozygous calls, which are the only informative markers here. */
function readParent(path: string): Map<string, 'AA' | 'BB'> {
  const out = new Map<string, 'AA' | 'BB'>()
  for (const line of readFileSync(path, 'utf8').split('\n')) {
    const p = line.split('\t')
    if (p.length < 7 || !AUTO.has(p[1])) continue
    if (p[6] === '0') out.set(p[0], 'AA')
    else if (p[6] === '2') out.set(p[0], 'BB')
  }
  return out
}

function readSample(path: string, parent: Map<string, 'AA' | 'BB'>) {
  const rows: Row[] = []
  let called = 0
  let het = 0
  let total = 0
  for (const line of readFileSync(path, 'utf8').split('\n')) {
    const p = line.split('\t')
    if (p.length < 7 || !AUTO.has(p[1])) continue
    total += 1
    if (p[6] !== '-1') { called += 1; if (p[6] === '1') het += 1 }
    const pg = parent.get(p[0])
    if (!pg || p[6] === '-1') continue
    rows.push({ chrom: p[1], parent: pg, gt: p[6] === '0' ? 'AA' : p[6] === '1' ? 'AB' : 'BB' })
  }
  return { rows, callRate: total ? called / total : 0, het: called ? het / called : 0 }
}

/**
 * Remove one parent's copy across a chromosome.
 *
 * At a parent-homozygous marker the sample carries the parent's allele plus whatever the other
 * parent sent. Removing the LOADED parent leaves only the other parent's, so a marker that read
 * heterozygous now reads homozygous for the opposite allele, which is an allele the loaded parent
 * does not have. Removing the OTHER parent leaves the loaded parent's allele, so the same marker
 * reads homozygous for it. Markers already homozygous are unchanged either way: they carry no
 * information about which copy went.
 */
function removeCopy(rows: readonly Row[], chrom: string, which: 'loaded' | 'other'): Row[] {
  return rows.map((r) => {
    if (r.chrom !== chrom || r.gt !== 'AB') return r
    const loadedAllele = r.parent === 'AA' ? 'AA' : 'BB'
    const otherAllele = r.parent === 'AA' ? 'BB' : 'AA'
    return { ...r, gt: which === 'loaded' ? otherAllele : loadedAllele }
  })
}

const files = readdirSync(DIR).filter((f) => f.endsWith('.probes')).sort()
// The parent array is the bulk diploid one, chosen by measurement rather than by name.
let parentFile = files[0]
let best = -1
for (const f of files) {
  const txt = readFileSync(join(DIR, f), 'utf8')
  let c = 0; let h = 0; let t = 0
  for (const line of txt.split('\n').slice(1, 200_000)) {
    const p = line.split('\t')
    if (p.length < 7) continue
    t += 1
    if (p[6] === '-1') continue
    c += 1
    if (p[6] === '1') h += 1
  }
  const cr = t ? c / t : 0
  const hr = c ? h / c : 0
  const score = cr > 0.95 && hr > 0.14 ? cr * hr : -1
  if (score > best) { best = score; parentFile = f }
}
const parent = readParent(join(DIR, parentFile))
process.stderr.write(`parent: ${parentFile}, ${parent.size} homozygous sites\n`)

let right = 0
let wrong = 0
let refused = 0
const perArray: { name: string; right: number; n: number; usable: boolean; stage: string }[] = []
console.log('sample            het    callRate  chr  removed   ->  called            band')
for (const f of files) {
  if (f === parentFile) continue
  const { rows, callRate, het } = readSample(join(DIR, f), parent)
  // BIPARENTAL ONLY. A uniparental genome has no second copy to remove, and this experiment would
  // be constructing a state the material cannot be in.
  if (het < 0.105) continue
  // AND THE ARRAY MUST BE READING A GENOME. The stage inference already decides this, and an array
  // it rejects is not weak evidence about a parent, it is not evidence about anything. Measured
  // both ways below rather than assumed.
  const st = inferStage({ hetRate: het, callRate })
  const usable = st.stage !== 'failed' && st.stage !== 'unknown'
  // The array's own dropout, which the caller needs to weigh a homozygous read. Estimated from
  // how far its heterozygosity falls below the parent's, since a lost allele turns a heterozygote
  // homozygous and nothing else here does.
  const ado = Math.min(0.5, Math.max(0.02, 1 - het / 0.17))
  let r = 0
  let n = 0
  for (const chrom of ['1', '2', '7', '13', '19']) {
    for (const which of ['loaded', 'other'] as const) {
      const mutated = removeCopy(rows, chrom, which)
      const onChrom = mutated.filter((x) => x.chrom === chrom)
      if (onChrom.length < 200) continue
      // `ado` is POSITIONAL and is a number. Passing an options object here made it NaN and
      // refused all 160 calls, which is a harness fault and not a caller one.
      const call = callOneParentOrigin(
        onChrom.map((x) => [x.parent, x.gt]) as never,
        ado,
      )
      const expect = which === 'loaded' ? 'known-parent-lost' : 'other-parent-lost'
      const got = call.verdict
      n += 1
      if (got === expect) { right += 1; r += 1 }
      else if (got === 'refused' || got === 'both-present') refused += 1
      else wrong += 1
      if (n <= 2) {
        console.log(`${f.split('.')[0].padEnd(16)} ${het.toFixed(3)}  ${callRate.toFixed(3)}    `
          + `${chrom.padStart(2)}  ${which.padEnd(7)}  ->  ${String(got).padEnd(18)}${call.band ?? ''}`)
      }
    }
  }
  if (n) perArray.push({ name: f, right: r, n, usable, stage: st.stage })
}
console.log()
console.log(`correct ${right}, wrong ${wrong}, refused ${refused}`)
const report = (label: string, rows: typeof perArray) => {
  const accs = rows.filter((a) => a.n >= 4).map((a) => a.right / a.n)
  if (accs.length < 2) { console.log(`${label}: too few arrays`); return }
  const m = accs.reduce((a, x) => a + x, 0) / accs.length
  const sd = Math.sqrt(accs.reduce((a, x) => a + (x - m) ** 2, 0) / (accs.length - 1))
  const tot = rows.reduce((a, x) => a + x.n, 0)
  const cor = rows.reduce((a, x) => a + x.right, 0)
  console.log(`${label.padEnd(34)} ${cor}/${tot} calls, per-array `
    + `${m.toFixed(4)} +/- ${(1.96 * sd / Math.sqrt(accs.length)).toFixed(4)} over ${accs.length} arrays`)
}
report('ALL biparental arrays', perArray)
report('arrays that resolve to a stage', perArray.filter((a) => a.usable))
report('arrays their own inference rejects', perArray.filter((a) => !a.usable))

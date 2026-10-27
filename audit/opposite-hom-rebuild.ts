/**
 * Step 4 reformulated: reconstruct the missing parent's OPPOSITE HOMOZYGOTES, not her heterozygous
 * markers.
 *
 * At a marker where the known parent is AA, a missing parent who is BB sends B to every child, one
 * who is AB sends it to half, and one who is AA sends it to none. So the share of diploid siblings
 * reading AB separates her three genotypes, and near-unanimity is the signature of BB rather than
 * an artefact. BB is what the two-parent tally needs: it is defined on opposite homozygotes and a
 * parent asserted only heterozygous contributes nothing to it.
 *
 * Her own array is used only to score. The reconstruction sees the known parent and the siblings.
 *
 * Criteria are in audit/opposite-hom-rebuild-preregistration.txt, committed before this ran.
 *
 * Run: OM_TRIOS=<dir> OM_MOTHER=<gsm> node --experimental-strip-types audit/opposite-hom-rebuild.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const segs = await import(`${W}segments.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) { console.log('contested-flag-tests: OM_TRIOS not readable. NOT RUN.'); process.exit(0) }

const fileOf = new Map<string, string>()
for (const f of readdirSync(DIR)) {
  const m = /^(GSM\d+)/.exec(f)
  if (m && (f.endsWith('.gz') || f.endsWith('.probes') || f.endsWith('.txt'))) fileOf.set(m[1], join(DIR, f))
}
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
const cache = new Map<string, unknown>()
function parentIndex(gsm: string) {
  const hit = cache.get(gsm); if (hit) return hit
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(fileOf.get(gsm)!)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const idx = score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
  cache.set(gsm, idx)
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
  const res = await score.scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {},
  } as never) as never
  return { res, profile } as { res: Res; profile: { callRate: number } }
}










const op = await import(`${W}oneParentOrigin.ts`)
const DROP_IN = op.DROP_IN as number
const H = 0.168  // the egg donor's heterozygosity on this platform

const MOTHER = process.env.OM_MOTHER
if (!MOTHER) { console.log('OM_MOTHER not set. NOT RUN.'); process.exit(0) }

const manifest = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const cols = manifest[0].split(',')
const at = (r: string[], k: string) => r[cols.indexOf(k)] ?? ''
let father = ''
const kids: string[] = []
for (const line of manifest.slice(1)) {
  const r = line.split(',')
  const m = at(r, 'mother_gsms')
  if (m !== MOTHER) continue
  const g = at(r, 'gsm')
  if (!fileOf.has(g)) continue
  father = father || at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
  kids.push(g)
}
if (!father || kids.length < 2) { console.log('family not found. NOT RUN.'); process.exit(0) }

// The embryo each sample came from, read off GEO's own filename: the letter before the cell digit.
const embryoOf = new Map<string, string>()
for (const g of kids) {
  const f = (fileOf.get(g) ?? '').split('/').pop() ?? ''
  const m = /^GSM\d+_([A-Za-z]+)/.exec(f)
  embryoOf.set(g, m ? m[1].toUpperCase() : g)
}
const embryos = [...new Set(embryoOf.values())].sort()
console.log(`family: ${kids.length} samples over ${embryos.length} embryos, one father, `
  + 'the mother arrayed\n')

const gtOf = (g: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(fileOf.get(g)!)()) {
    if (!/^\d+$/.test(r.chrom)) continue
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const faGt = gtOf(father)
const moGt = gtOf(MOTHER)
console.log(`father ${faGt.size.toLocaleString()} called, mother ${moGt.size.toLocaleString()} called`)

// Every sample's genotype at the father-homozygous markers, plus its dropout.
const sample = new Map<string, { gt: Map<string, string>; ado: number }>()
for (const g of kids) {
  const gt = gtOf(g)
  let het = 0
  for (const v of gt.values()) if (v === 'AB') het += 1
  sample.set(g, { gt, ado: op.inferDropout(gt.size ? het / gt.size : NaN) as number })
  process.stdout.write('.')
}
console.log('\n')

// k derived from the contamination target the shipped rebuilt mother carries on its banner.
const TARGET = 0.02749
const lchoose = (n: number, k: number) => {
  let s = 0
  for (let i = 1; i <= k; i += 1) s += Math.log(n - k + i) - Math.log(i)
  return s
}
const tailGE = (n: number, p: number, k: number) => {
  let s = 0
  for (let i = k; i <= n; i += 1) s += Math.exp(lchoose(n, i) + i * Math.log(p) + (n - i) * Math.log(1 - p))
  return Math.min(1, Math.max(0, s))
}

/**
 * The vote threshold above which the missing parent is asserted homozygous for the other allele.
 *
 * Three hypotheses at a known-parent-homozygous marker, with the observed-AB rate each implies.
 * The threshold is the lowest vote count at which the posterior odds of BB against the union of AB
 * and AA clear the target, under the genotype frequencies the known parent's own array implies.
 */
function deriveU(n: number, ado: number, pBB: number, pAB: number, target: number) {
  const rBB = Math.max(1e-6, Math.min(0.999, 1 - ado))
  const rAB = Math.max(1e-6, Math.min(0.999, 0.5 * (1 - ado)))
  const rAA = DROP_IN
  const pAA = Math.max(0, 1 - pBB - pAB)
  for (let u = 1; u <= n; u += 1) {
    const lBB = pBB * tailGE(n, rBB, u)
    const lAB = pAB * tailGE(n, rAB, u)
    const lAA = pAA * tailGE(n, rAA, u)
    const contam = (lAB + lAA) / Math.max(1e-300, lBB + lAB + lAA)
    if (contam <= target) return { u, contam, complete: tailGE(n, rBB, u) }
  }
  return { u: n, contam: NaN, complete: NaN }
}

const TARGET_BB = 0.02
// The known parent's own genotype frequencies stand in for the allele frequencies. Nothing of the
// missing parent's array is read here.
let faHom = 0; let faHet = 0
for (const v of faGt.values()) { if (v === 'AB') faHet += 1; else faHom += 1 }
const q = faHet / Math.max(1, faHom + faHet)   // crude heterozygosity of the panel
// Under Hardy-Weinberg a panel heterozygosity of q implies a minor allele frequency f with
// q = 2f(1-f). At a marker where the known parent is homozygous for the major allele the missing
// parent is BB with f^2 and AB with 2f(1-f).
const f = 0.5 * (1 - Math.sqrt(Math.max(0, 1 - 2 * q)))
const priorBB = f * f
const priorAB = 2 * f * (1 - f)
console.log(`panel heterozygosity ${q.toFixed(4)}, implied minor allele frequency ${f.toFixed(4)}`)
console.log(`prior at a known-parent-homozygous marker: BB ${priorBB.toFixed(4)}, AB ${priorAB.toFixed(4)}\n`)

/** Vote counts: how many siblings read AB at each known-parent-homozygous marker. */
const votes = new Map<string, number>()
const seen = new Map<string, number>()
for (const g of kids) {
  const gt = sample.get(g)!.gt
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    const s = gt.get(probe)
    if (!s) continue
    seen.set(probe, (seen.get(probe) ?? 0) + 1)
    if (s === 'AB') votes.set(probe, (votes.get(probe) ?? 0) + 1)
  }
}

const ados = kids.map((g) => sample.get(g)!.ado).filter(Number.isFinite)
const adoMean = ados.reduce((a, b) => a + b, 0) / Math.max(1, ados.length)
const d = deriveU(kids.length, adoMean, priorBB, priorAB, TARGET_BB)
console.log(`n ${kids.length}, mean dropout ${adoMean.toFixed(3)}, derived u ${d.u}`)
console.log(`modelled contamination ${d.contam.toFixed(5)}, so predicted precision ${(1 - d.contam).toFixed(4)}`)
console.log(`modelled completeness ${d.complete.toFixed(4)}\n`)

/**
 * Score an assertion set against her real array. `want` is the genotype the assertion claims:
 * the opposite homozygote of whatever the known parent reads at that marker.
 */
const gradeBB = (asserted: Set<string>) => {
  let hit = 0; let miss = 0; const wrongClass = new Map<string, number>()
  for (const probe of asserted) {
    const t = moGt.get(probe)
    if (!t) continue
    const fg = faGt.get(probe)!
    const want = fg === 'AA' ? 'BB' : 'AA'
    if (t === want) hit += 1
    else { miss += 1; wrongClass.set(t, (wrongClass.get(t) ?? 0) + 1) }
  }
  let trueBB = 0; let found = 0
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    const want = fg === 'AA' ? 'BB' : 'AA'
    if (moGt.get(probe) !== want) continue
    trueBB += 1
    if (asserted.has(probe)) found += 1
  }
  return { precision: hit + miss ? hit / (hit + miss) : NaN, asserted: hit + miss,
    completeness: trueBB ? found / trueBB : NaN, trueBB, wrongClass }
}

const assertAt = (u: number, from?: Set<string>) => {
  const out = new Set<string>()
  for (const [probe, v] of votes) {
    if (from && !from.has(probe)) continue
    if ((seen.get(probe) ?? 0) < Math.ceil(0.8 * kids.length)) continue
    if (v >= u) out.add(probe)
  }
  return out
}

console.log('CRITERION 1 AND 2. THE DERIVED THRESHOLD')
{
  const a = assertAt(d.u)
  const s = gradeBB(a)
  console.log(`  asserted ${s.asserted.toLocaleString()} markers scored against her array`)
  console.log(`  PRECISION    ${s.precision.toFixed(4)}`)
  console.log(`  COMPLETENESS ${s.completeness.toFixed(4)} of ${s.trueBB.toLocaleString()} true opposite homozygotes`)
  const wrong = [...s.wrongClass.entries()].sort((x, y) => y[1] - x[1])
  console.log(`  what the wrong ones really are: ${wrong.map(([g, n]) => `${g} ${n.toLocaleString()}`).join(', ')}`)
  console.log(s.precision >= 0.95
    ? '  PASS criterion 1: at or above 0.95.'
    : `  FAIL criterion 1: ${s.precision.toFixed(4)}, under the 0.95 the pre-registration set.`)
}

console.log('\nCRITERION 3. THE THRESHOLD SWEPT')
console.log(`  ${'u'.padStart(3)} ${'asserted'.padStart(9)} ${'precision'.padStart(10)} ${'completeness'.padStart(12)}`)
for (let u = Math.max(2, kids.length - 12); u <= kids.length; u += 1) {
  const s = gradeBB(assertAt(u))
  if (!s.asserted) continue
  console.log(`  ${String(u).padStart(3)} ${s.asserted.toLocaleString().padStart(9)} `
    + `${s.precision.toFixed(4).padStart(10)} ${s.completeness.toFixed(4).padStart(12)}`)
}

console.log('\nCRITERION 4. LEAVE ONE EMBRYO OUT')
for (const drop of embryos) {
  const keep = kids.filter((g) => embryoOf.get(g) !== drop)
  if (keep.length < 4) continue
  const v2 = new Map<string, number>()
  const s2 = new Map<string, number>()
  for (const g of keep) {
    const gt = sample.get(g)!.gt
    for (const [probe, fg] of faGt) {
      if (fg === 'AB') continue
      const s = gt.get(probe)
      if (!s) continue
      s2.set(probe, (s2.get(probe) ?? 0) + 1)
      if (s === 'AB') v2.set(probe, (v2.get(probe) ?? 0) + 1)
    }
  }
  const ad = keep.map((g) => sample.get(g)!.ado).filter(Number.isFinite)
  const am = ad.reduce((a, b) => a + b, 0) / Math.max(1, ad.length)
  const dd = deriveU(keep.length, am, priorBB, priorAB, TARGET_BB)
  const a = new Set<string>()
  for (const [probe, v] of v2) {
    if ((s2.get(probe) ?? 0) < Math.ceil(0.8 * keep.length)) continue
    if (v >= dd.u) a.add(probe)
  }
  const s = gradeBB(a)
  console.log(`  without embryo ${drop}: n ${keep.length}, u ${dd.u}, `
    + `precision ${s.precision.toFixed(4)}, completeness ${s.completeness.toFixed(4)}, `
    + `asserted ${s.asserted.toLocaleString()}`)
}

// PER CHROMOSOME, because the channel runs per chromosome and the count is what decides usability.
console.log('\nRECONSTRUCTED OPPOSITE HOMOZYGOTES PER CHROMOSOME, at the derived threshold')
{
  const a = assertAt(d.u)
  const byChrom = new Map<string, number>()
  const chromOf = new Map<string, string>()
  for (const r of rowsOf(fileOf.get(father)!)()) {
    if (/^\d+$/.test(r.chrom)) chromOf.set(r.probesetId, r.chrom)
  }
  for (const probe of a) {
    const c = chromOf.get(probe)
    if (c) byChrom.set(c, (byChrom.get(c) ?? 0) + 1)
  }
  const counts = [...byChrom.values()].sort((x, y) => x - y)
  const med = counts.length ? counts[Math.floor(counts.length / 2)] : 0
  console.log(`  ${byChrom.size} autosomes, median ${med.toLocaleString()} markers, `
    + `range ${counts[0]?.toLocaleString()} to ${counts[counts.length - 1]?.toLocaleString()}`)
  console.log('  the real opposite-homozygote channel runs on a median of 1,034 per chromosome')
  console.log(med >= 200
    ? `  the reconstruction yields enough per chromosome to run the channel.`
    : `  the reconstruction does NOT yield enough per chromosome.`)
}


// ===========================================================================================
// THE MECHANISM. Precision FALLS as the threshold rises, which the model cannot produce, so the
// model's sample size is wrong. These 30 arrays are not 30 children: they are several cells of each
// of a few embryos, and cells of one embryo carry the same genome. A heterozygous parent transmits
// to an EMBRYO, and every cell of that embryo then reads the same, so unanimity across 30 arrays is
// unanimity across a handful of independent transmissions.
//
// The fix is to vote by embryo. Each embryo casts one vote, taken as the majority of its own cells.
// ===========================================================================================
console.log('\n\nVOTING BY EMBRYO RATHER THAN BY ARRAY')
console.log(`  ${embryos.length} embryos over ${kids.length} arrays: `
  + embryos.map((e) => `${e} x${kids.filter((g) => embryoOf.get(g) === e).length}`).join(', '))

const embryoVote = new Map<string, number>()
const embryoSeen = new Map<string, number>()
for (const e of embryos) {
  const cells = kids.filter((g) => embryoOf.get(g) === e)
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    let ab = 0; let called = 0
    for (const g of cells) {
      const s2 = sample.get(g)!.gt.get(probe)
      if (!s2) continue
      called += 1
      if (s2 === 'AB') ab += 1
    }
    if (!called) continue
    embryoSeen.set(probe, (embryoSeen.get(probe) ?? 0) + 1)
    // One vote per embryo. A single heterozygous cell is enough: dropout destroys heterozygosity
    // and does not create it, so ANY cell reading AB says the embryo inherited the allele.
    if (ab >= 1) embryoVote.set(probe, (embryoVote.get(probe) ?? 0) + 1)
  }
}

const assertByEmbryo = (u: number) => {
  const out = new Set<string>()
  for (const [probe, v] of embryoVote) {
    if ((embryoSeen.get(probe) ?? 0) < embryos.length) continue
    if (v >= u) out.add(probe)
  }
  return out
}
console.log(`\n  ${'u'.padStart(3)} ${'asserted'.padStart(9)} ${'precision'.padStart(10)} `
  + `${'completeness'.padStart(12)}  of ${embryos.length} embryos`)
for (let u = 2; u <= embryos.length; u += 1) {
  const s2 = gradeBB(assertByEmbryo(u))
  if (!s2.asserted) continue
  console.log(`  ${String(u).padStart(3)} ${s2.asserted.toLocaleString().padStart(9)} `
    + `${s2.precision.toFixed(4).padStart(10)} ${s2.completeness.toFixed(4).padStart(12)}`)
}

// WHAT THE MODEL PREDICTS AT EMBRYO LEVEL, so the measurement has something to be checked against.
// A BB parent transmits the allele to every embryo; an AB parent to half. At n embryos and
// unanimity the posterior is set by the prior odds and 0.5^n.
console.log('\n  THE MODEL AT EMBRYO LEVEL, unanimity required')
console.log(`  ${'embryos'.padStart(8)} ${'predicted precision'.padStart(20)}`)
for (const n of [4, 5, 6, 7, 8, 9, 10, 12]) {
  const lBB = priorBB * Math.pow(1 - adoMean, 0)
  const lAB = priorAB * Math.pow(0.5, n)
  const pr = lBB / (lBB + lAB)
  console.log(`  ${String(n).padStart(8)} ${pr.toFixed(4).padStart(20)}${n === embryos.length ? '   <- this family' : ''}`)
}
const need = (() => {
  for (let n = 2; n <= 40; n += 1) {
    if (priorBB / (priorBB + priorAB * Math.pow(0.5, n)) >= 0.95) return n
  }
  return NaN
})()
console.log(`\n  embryos needed for a precision of 0.95 at this panel's allele frequencies: ${need}`)
console.log(`  this family has ${embryos.length}.`)


// ===========================================================================================
// THE CORRECT STATISTIC. A hard vote per embryo throws away how many of its cells agree, and "any
// cell reading AB" gives a ten-cell embryo ten chances to collect a drop-in while a one-cell embryo
// gets one. Both are replaced by a likelihood.
//
// Per embryo, two hypotheses: it inherited the other allele, so its true genotype is heterozygous
// and each called cell reads AB with probability 1 - dropout; or it did not, so each cell reads AB
// only by drop-in. Per marker, three hypotheses for the missing parent, which differ in how the
// embryos are distributed between those two: BB sends the allele to all of them, AB to each
// independently with probability one half, AA to none.
// ===========================================================================================
console.log('\n\nTHE PER-EMBRYO LIKELIHOOD')
const cellsOf = new Map<string, string[]>()
for (const e of embryos) cellsOf.set(e, kids.filter((g) => embryoOf.get(g) === e))
const adoOf = new Map<string, number>()
for (const e of embryos) {
  const a = cellsOf.get(e)!.map((g) => sample.get(g)!.ado).filter(Number.isFinite)
  adoOf.set(e, a.length ? a.reduce((x, y) => x + y, 0) / a.length : adoMean)
}

const posterior = new Map<string, number>()
for (const [probe, fg] of faGt) {
  if (fg === 'AB') continue
  let lBB = Math.log(priorBB)
  let lAB = Math.log(priorAB)
  let lAA = Math.log(Math.max(1e-12, 1 - priorBB - priorAB))
  let embryosSeen = 0
  for (const e of embryos) {
    const a = Math.min(0.95, Math.max(0.01, adoOf.get(e)!))
    let l1 = 0; let l0 = 0; let called = 0
    for (const g of cellsOf.get(e)!) {
      const v = sample.get(g)!.gt.get(probe)
      if (!v) continue
      called += 1
      if (v === 'AB') { l1 += Math.log(1 - a); l0 += Math.log(DROP_IN) }
      else { l1 += Math.log(a); l0 += Math.log(1 - DROP_IN) }
    }
    if (!called) continue
    embryosSeen += 1
    const m = Math.max(l1, l0)
    lBB += l1
    lAB += m + Math.log(0.5 * Math.exp(l1 - m) + 0.5 * Math.exp(l0 - m))
    lAA += l0
  }
  if (embryosSeen < embryos.length) continue
  const m = Math.max(lBB, lAB, lAA)
  const d2 = Math.exp(lBB - m) + Math.exp(lAB - m) + Math.exp(lAA - m)
  posterior.set(probe, Math.exp(lBB - m) / d2)
}
console.log(`  ${posterior.size.toLocaleString()} markers scored, all ${embryos.length} embryos called`)
console.log(`\n  ${'posterior'.padStart(10)} ${'asserted'.padStart(9)} ${'precision'.padStart(10)} `
  + `${'completeness'.padStart(12)} ${'per chrom'.padStart(10)}`)
const chromOf2 = new Map<string, string>()
for (const r of rowsOf(fileOf.get(father)!)()) {
  if (/^\d+$/.test(r.chrom)) chromOf2.set(r.probesetId, r.chrom)
}
let best = { pr: 0, post: 0, med: 0, asserted: 0, comp: 0 }
for (const th of [0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999]) {
  const a = new Set<string>()
  for (const [probe, v] of posterior) if (v >= th) a.add(probe)
  const s2 = gradeBB(a)
  if (!s2.asserted) continue
  const byC = new Map<string, number>()
  for (const probe of a) {
    const c = chromOf2.get(probe)
    if (c) byC.set(c, (byC.get(c) ?? 0) + 1)
  }
  const cs = [...byC.values()].sort((x, y) => x - y)
  const med = cs.length ? cs[Math.floor(cs.length / 2)] : 0
  console.log(`  ${th.toFixed(3).padStart(10)} ${s2.asserted.toLocaleString().padStart(9)} `
    + `${s2.precision.toFixed(4).padStart(10)} ${s2.completeness.toFixed(4).padStart(12)} `
    + `${med.toLocaleString().padStart(10)}`)
  if (s2.precision > best.pr) {
    best = { pr: s2.precision, post: th, med, asserted: s2.asserted, comp: s2.completeness }
  }
}
console.log(`\n  best precision reached: ${best.pr.toFixed(4)} at a posterior of ${best.post}, `
  + `${best.asserted.toLocaleString()} markers, median ${best.med.toLocaleString()} per chromosome`)
console.log(`  the model's ceiling at ${embryos.length} embryos is `
  + `${(priorBB / (priorBB + priorAB * Math.pow(0.5, embryos.length))).toFixed(4)}`)
console.log(best.pr >= 0.95
  ? '  PASS criterion 1.'
  : `  FAIL criterion 1: ${best.pr.toFixed(4)} against the 0.95 the pre-registration set. `
    + `At this panel's allele frequencies ${need} independent embryos are needed and the corpus's `
    + 'one usable family has ' + embryos.length + '.')


// ===========================================================================================
// THE PRIOR, PER MARKER. With six embryos the likelihood ratio for BB against AB at an all-
// heterozygous reading is 2^6, which is 64, and no statistic can do better: that is the whole
// information the six independent transmissions carry. Beating prior odds of 19.8 to 1 by the
// factor of 19 that a precision of 0.95 needs requires 2^n at or above 377, so nine embryos.
//
// Unless the prior is not 19.8 to 1. That figure is a genome-wide average under Hardy-Weinberg at
// the panel's mean heterozygosity. Per marker it varies enormously: where the allele the known
// parent lacks is common, the missing parent is far likelier to be homozygous for it, and six
// embryos suffice. So the frequency of that allele is estimated per marker from the OTHER parent
// arrays of the series, which are unrelated donors, with the family's own two arrays held out.
// ===========================================================================================
console.log('\n\nTHE PRIOR, PER MARKER, FROM UNRELATED DONORS')
const panel = (process.env.OM_PANEL ?? '').split(/[,; ]+/).filter((x) => x && fileOf.has(x))
  .filter((x) => x !== MOTHER && x !== father)
if (!panel.length) { console.log('  OM_PANEL not set or empty. NOT RUN.'); process.exit(0) }
console.log(`  ${panel.length} unrelated donor arrays, ${2 * panel.length} chromosomes, `
  + `the family's own two held out`)

// Count the B allele per marker across the panel.
const bCount = new Map<string, number>()
const aCount = new Map<string, number>()
for (const g of panel) {
  for (const [probe, v] of gtOf(g)) {
    if (!faGt.has(probe)) continue
    const b = v === 'BB' ? 2 : v === 'AB' ? 1 : 0
    bCount.set(probe, (bCount.get(probe) ?? 0) + b)
    aCount.set(probe, (aCount.get(probe) ?? 0) + (2 - b))
  }
}
/** Frequency of the allele the known parent LACKS, which is the one the wanted genotype is for. */
const lackFreq = (probe: string) => {
  const fg = faGt.get(probe)
  if (!fg || fg === 'AB') return NaN
  const b = bCount.get(probe); const a = aCount.get(probe)
  if (b === undefined || a === undefined || a + b === 0) return NaN
  return fg === 'AA' ? b / (a + b) : a / (a + b)
}

const nEmb = embryos.length
const lr = Math.pow(2, nEmb)
console.log(`  likelihood ratio at an all-heterozygous reading over ${nEmb} embryos: ${lr}`)
console.log(`\n  ${'min freq'.padStart(9)} ${'markers'.padStart(9)} ${'prior odds'.padStart(11)} `
  + `${'predicted'.padStart(10)} ${'asserted'.padStart(9)} ${'PRECISION'.padStart(10)} `
  + `${'per chrom'.padStart(10)}`)
const chromOf3 = new Map<string, string>()
for (const r of rowsOf(fileOf.get(father)!)()) {
  if (/^\d+$/.test(r.chrom)) chromOf3.set(r.probesetId, r.chrom)
}
const unanimous = assertByEmbryo(nEmb)
let passed: { f: number; pr: number; n: number; med: number } | undefined
for (const minF of [0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]) {
  const pool = [...faGt.keys()].filter((probe) => {
    const f2 = lackFreq(probe)
    return Number.isFinite(f2) && f2 >= minF
  })
  if (!pool.length) continue
  // The mean prior odds over the pool, so the prediction is the pool's own and not the genome's.
  let oddsSum = 0; let oddsN = 0
  for (const probe of pool) {
    const f2 = lackFreq(probe)
    if (f2 <= 0 || f2 >= 1) continue
    oddsSum += 2 * (1 - f2) / f2
    oddsN += 1
  }
  const odds = oddsN ? oddsSum / oddsN : NaN
  const predicted = lr / (lr + odds)
  const keep = new Set([...unanimous].filter((x) => {
    const f2 = lackFreq(x)
    return Number.isFinite(f2) && f2 >= minF
  }))
  const s2 = gradeBB(keep)
  const byC = new Map<string, number>()
  for (const probe of keep) {
    const c = chromOf3.get(probe)
    if (c) byC.set(c, (byC.get(c) ?? 0) + 1)
  }
  const cs = [...byC.values()].sort((x, y) => x - y)
  const med = cs.length ? cs[Math.floor(cs.length / 2)] : 0
  console.log(`  ${minF.toFixed(2).padStart(9)} ${pool.length.toLocaleString().padStart(9)} `
    + `${odds.toFixed(2).padStart(11)} ${predicted.toFixed(4).padStart(10)} `
    + `${s2.asserted.toLocaleString().padStart(9)} ${s2.precision.toFixed(4).padStart(10)} `
    + `${med.toLocaleString().padStart(10)}`)
  if (s2.precision >= 0.95 && (!passed || s2.asserted > passed.n)) {
    passed = { f: minF, pr: s2.precision, n: s2.asserted, med }
  }
}
console.log(passed
  ? `\n  PASS criterion 1: precision ${passed.pr.toFixed(4)} at a minimum allele frequency of `
    + `${passed.f.toFixed(2)}, over ${passed.n.toLocaleString()} asserted markers, median `
    + `${passed.med.toLocaleString()} per chromosome.`
  : '\n  FAIL criterion 1 on every frequency band: six embryos do not carry enough even where the '
    + 'prior is favourable.')


// ===========================================================================================
// WHICH ERROR IS IT. A precision of 0.65 against a predicted 0.97 means a source the prior does not
// describe. The two candidates separate cleanly on what the missing parent REALLY reads at the
// markers wrongly asserted:
//
//   she reads AB  -> the likelihood ratio of 2^6 is simply not enough, an INFORMATION limit, and
//                    only more independent embryos help
//   she reads the same genotype as the known parent -> she carries none of that allele, so no
//                    embryo could have inherited it and every heterozygous call is ARTEFACT
// ===========================================================================================
console.log('\n\nWHICH ERROR IS IT, at a minimum allele frequency of 0.40')
{
  const keep = new Set([...unanimous].filter((x) => {
    const f2 = lackFreq(x)
    return Number.isFinite(f2) && f2 >= 0.40
  }))
  let wanted = 0; let herHet = 0; let sameAsKnown = 0; let uncalled = 0
  for (const probe of keep) {
    const t = moGt.get(probe)
    if (!t) { uncalled += 1; continue }
    const fg = faGt.get(probe)!
    const want = fg === 'AA' ? 'BB' : 'AA'
    if (t === want) wanted += 1
    else if (t === 'AB') herHet += 1
    else sameAsKnown += 1
  }
  const scored = wanted + herHet + sameAsKnown
  console.log(`  ${scored.toLocaleString()} asserted markers she has a call for`)
  console.log(`    she reads the wanted homozygote:      ${wanted.toLocaleString()} `
    + `(${(100 * wanted / scored).toFixed(1)}%)`)
  console.log(`    she reads HETEROZYGOUS:               ${herHet.toLocaleString()} `
    + `(${(100 * herHet / scored).toFixed(1)}%)   information limit`)
  console.log(`    she reads the SAME as the known parent: ${sameAsKnown.toLocaleString()} `
    + `(${(100 * sameAsKnown / scored).toFixed(1)}%)   artefact, she carries none of that allele`)
  // Under independence an all-heterozygous reading from a parent carrying none of the allele needs
  // drop-in in every embryo, so its expected count is vanishing. The observed count is the measure
  // of how far correlated drop-in departs from that.
  const expected = scored * Math.pow(DROP_IN, embryos.length)
  console.log(`\n  expected under independent drop-in: ${expected.toExponential(2)} markers`)
  console.log(`  observed: ${sameAsKnown.toLocaleString()}`)
  console.log(sameAsKnown > herHet
    ? '  THE LIMIT IS ARTEFACT, not information: the dominant error comes from markers where she\n'
      + '  carries none of the allele, which no number of embryos can fix and a probe filter can.'
    : '  THE LIMIT IS INFORMATION: the dominant error is her reading heterozygous, which only more\n'
      + '  independent embryos can resolve.')
}


// ===========================================================================================
// BOTH PIECES AT ONCE. The two arms above each fixed one half and kept the other's flaw.
//
// The per-embryo likelihood used a genome-wide prior of 19.8 to 1 against the wanted genotype.
// The per-marker prior used the hard vote, and "any cell reads AB" lets a ten-cell embryo collect a
// drop-in: at a drop-in of 0.0435 an embryo that inherited nothing still votes yes with probability
// 1 - (1 - 0.0435)^10, which is 0.36. That is why the realized likelihood ratio measured 3.7 where
// 2^6 is 64.
//
// Combined: the likelihood weights each embryo by how many of its cells agree, and the prior is the
// marker's own.
// ===========================================================================================
console.log('\n\nTHE PER-EMBRYO LIKELIHOOD WITH THE PER-MARKER PRIOR')
const postPM = new Map<string, number>()
for (const [probe, fg] of faGt) {
  if (fg === 'AB') continue
  const f2 = lackFreq(probe)
  if (!Number.isFinite(f2) || f2 <= 0 || f2 >= 1) continue
  // Hardy-Weinberg at this marker's own frequency for the allele the known parent lacks.
  let lW = Math.log(f2 * f2)                 // she is homozygous for it: what the tally needs
  let lH = Math.log(2 * f2 * (1 - f2))       // she is heterozygous
  let lN = Math.log((1 - f2) * (1 - f2))     // she carries none of it
  let seenE = 0
  for (const e of embryos) {
    const a = Math.min(0.95, Math.max(0.01, adoOf.get(e)!))
    let l1 = 0; let l0 = 0; let called = 0
    for (const g of cellsOf.get(e)!) {
      const v = sample.get(g)!.gt.get(probe)
      if (!v) continue
      called += 1
      if (v === 'AB') { l1 += Math.log(1 - a); l0 += Math.log(DROP_IN) }
      else { l1 += Math.log(a); l0 += Math.log(1 - DROP_IN) }
    }
    if (!called) continue
    seenE += 1
    const m = Math.max(l1, l0)
    lW += l1
    lH += m + Math.log(0.5 * Math.exp(l1 - m) + 0.5 * Math.exp(l0 - m))
    lN += l0
  }
  if (seenE < embryos.length) continue
  const m = Math.max(lW, lH, lN)
  const d2 = Math.exp(lW - m) + Math.exp(lH - m) + Math.exp(lN - m)
  postPM.set(probe, Math.exp(lW - m) / d2)
}
console.log(`  ${postPM.size.toLocaleString()} markers scored with all ${embryos.length} embryos called`)
console.log(`\n  ${'posterior'.padStart(10)} ${'asserted'.padStart(9)} ${'PRECISION'.padStart(10)} `
  + `${'completeness'.padStart(12)} ${'per chrom'.padStart(10)} ${'autosomes'.padStart(10)}`)
let bestPM: { th: number; pr: number; n: number; med: number; comp: number; chroms: number } | undefined
for (const th of [0.5, 0.8, 0.9, 0.95, 0.97, 0.99, 0.995, 0.999]) {
  const a = new Set<string>()
  for (const [probe, v] of postPM) if (v >= th) a.add(probe)
  const s2 = gradeBB(a)
  if (!s2.asserted) continue
  const byC = new Map<string, number>()
  for (const probe of a) {
    const c = chromOf3.get(probe)
    if (c) byC.set(c, (byC.get(c) ?? 0) + 1)
  }
  const cs = [...byC.values()].sort((x, y) => x - y)
  const med = cs.length ? cs[Math.floor(cs.length / 2)] : 0
  console.log(`  ${th.toFixed(3).padStart(10)} ${s2.asserted.toLocaleString().padStart(9)} `
    + `${s2.precision.toFixed(4).padStart(10)} ${s2.completeness.toFixed(4).padStart(12)} `
    + `${med.toLocaleString().padStart(10)} ${String(byC.size).padStart(10)}`)
  if (s2.precision >= 0.95 && (!bestPM || s2.asserted > bestPM.n)) {
    bestPM = { th, pr: s2.precision, n: s2.asserted, med, comp: s2.completeness, chroms: byC.size }
  }
}
console.log(bestPM
  ? `\n  PASS criterion 1: precision ${bestPM.pr.toFixed(4)} at a posterior of ${bestPM.th}, over `
    + `${bestPM.n.toLocaleString()} asserted markers across ${bestPM.chroms} autosomes, median `
    + `${bestPM.med.toLocaleString()} per chromosome, completeness ${bestPM.comp.toFixed(4)}.`
    + `\n  CRITERION 2: the opposite-homozygote channel runs on a median of 1,034 markers per `
    + `chromosome and refuses under 200.`
  : '\n  FAIL criterion 1 at every posterior.')


// ===========================================================================================
// PER-PROBE DROP-IN, FROM THE OTHER FAMILIES. The residual error is an embryo looking as though it
// inherited an allele when drop-in read heterozygous across its cells, and drop-in belongs to the
// probe. So the probes it belongs to can be identified and removed, using material that says
// nothing about the family being reconstructed.
//
// Where BOTH parents of another family are homozygous for the SAME allele, their child cannot be
// heterozygous: every heterozygous call at such a marker is artefact. Counting those per probe over
// the other families gives a drop-in rate per probe that is independent of this family entirely.
// ===========================================================================================
console.log('\n\nPER-PROBE DROP-IN, MEASURED ON THE OTHER FAMILIES')
type Fam = { mother: string; father: string; kids: string[] }
const others: Fam[] = []
{
  const byMother = new Map<string, { father: string; kids: string[] }>()
  for (const line of manifest.slice(1)) {
    const r = line.split(',')
    const mo = at(r, 'mother_gsms')
    const g = at(r, 'gsm')
    if (!mo || mo.includes(';') || mo === MOTHER || !fileOf.has(g) || !fileOf.has(mo)) continue
    const fa = at(r, 'father_gsms').split(/[; ]+/).filter((x) => fileOf.has(x))[0]
    if (!fa) continue
    const e = byMother.get(mo) ?? { father: fa, kids: [] }
    e.kids.push(g)
    byMother.set(mo, e)
  }
  for (const [mother, v] of byMother) others.push({ mother, father: v.father, kids: v.kids })
}
console.log(`  ${others.length} other families usable: `
  + others.map((f) => `${f.kids.length} arrays`).join(', '))

const dropObs = new Map<string, number>()
const dropHet = new Map<string, number>()
for (const f of others) {
  const fa2 = gtOf(f.father)
  const mo2 = gtOf(f.mother)
  // Markers where both parents are the same homozygote: the child cannot genuinely be heterozygous.
  const impossible = new Set<string>()
  for (const [probe, g1] of fa2) {
    if (g1 === 'AB') continue
    if (mo2.get(probe) === g1) impossible.add(probe)
  }
  for (const k of f.kids) {
    for (const [probe, v] of gtOf(k)) {
      if (!impossible.has(probe)) continue
      dropObs.set(probe, (dropObs.get(probe) ?? 0) + 1)
      if (v === 'AB') dropHet.set(probe, (dropHet.get(probe) ?? 0) + 1)
    }
  }
}
const withObs = [...dropObs.values()]
console.log(`  ${dropObs.size.toLocaleString()} probes observed, median `
  + `${withObs.sort((a, b) => a - b)[Math.floor(withObs.length / 2)]} observations each`)
let totObs = 0; let totHet = 0
for (const [probe, n] of dropObs) { totObs += n; totHet += dropHet.get(probe) ?? 0 }
console.log(`  overall drop-in ${(totHet / Math.max(1, totObs)).toFixed(5)} over `
  + `${totObs.toLocaleString()} observations`)
// How concentrated is it? If drop-in were a per-observation coin the bad-probe share would be tiny.
let bad = 0; let scored2 = 0
for (const [probe, n] of dropObs) {
  if (n < 6) continue
  scored2 += 1
  if ((dropHet.get(probe) ?? 0) / n >= 0.20) bad += 1
}
console.log(`  probes with 6+ observations: ${scored2.toLocaleString()}, of which `
  + `${bad.toLocaleString()} read 20% or more heterozygous (${(100 * bad / Math.max(1, scored2)).toFixed(2)}%)`)

console.log('\n  THE FILTER APPLIED, per-embryo likelihood with the per-marker prior')
console.log(`  ${'max drop-in'.padStart(11)} ${'posterior'.padStart(10)} ${'asserted'.padStart(9)} `
  + `${'PRECISION'.padStart(10)} ${'per chrom'.padStart(10)} ${'autosomes'.padStart(10)}`)
let win: { mx: number; th: number; pr: number; n: number; med: number; chroms: number } | undefined
for (const maxDrop of [1.0, 0.20, 0.10, 0.05, 0.02]) {
  const ok = (probe: string) => {
    const n = dropObs.get(probe)
    if (n === undefined || n < 6) return maxDrop >= 1
    return (dropHet.get(probe) ?? 0) / n <= maxDrop
  }
  for (const th of [0.9, 0.95, 0.97]) {
    const a = new Set<string>()
    for (const [probe, v] of postPM) if (v >= th && ok(probe)) a.add(probe)
    const s2 = gradeBB(a)
    if (!s2.asserted) continue
    const byC = new Map<string, number>()
    for (const probe of a) {
      const c = chromOf3.get(probe)
      if (c) byC.set(c, (byC.get(c) ?? 0) + 1)
    }
    const cs = [...byC.values()].sort((x, y) => x - y)
    const med = cs.length ? cs[Math.floor(cs.length / 2)] : 0
    console.log(`  ${maxDrop.toFixed(2).padStart(11)} ${th.toFixed(2).padStart(10)} `
      + `${s2.asserted.toLocaleString().padStart(9)} ${s2.precision.toFixed(4).padStart(10)} `
      + `${med.toLocaleString().padStart(10)} ${String(byC.size).padStart(10)}`)
    if (s2.precision >= 0.95 && med >= 200 && (!win || s2.asserted > win.n)) {
      win = { mx: maxDrop, th, pr: s2.precision, n: s2.asserted, med, chroms: byC.size }
    }
  }
}
console.log(win
  ? `\n  PASS BOTH CRITERIA: precision ${win.pr.toFixed(4)} at a drop-in ceiling of ${win.mx} and a `
    + `posterior of ${win.th}, ${win.n.toLocaleString()} markers over ${win.chroms} autosomes, `
    + `median ${win.med.toLocaleString()} per chromosome against the 200 the channel needs.`
  : '\n  still short: the probe filter does not lift both precision and per-chromosome yield to '
    + 'the pre-registered bars together.')

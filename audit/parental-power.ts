/**
 * IS THE TOOL EQUALLY ABLE TO SEE A PATERNAL LOSS AND A MATERNAL ONE?
 *
 * THIS IS THE QUESTION THAT DECIDES WHETHER THE BIOLOGICAL QUESTION CAN BE ASKED AT ALL. The lab
 * wants to know whether breakpoints and losses affect the paternal and maternal genomes differently.
 * A detector with unequal power on the two sides answers that question by construction, whatever the
 * biology does. Nobody has checked which way this one leans.
 *
 * THERE IS A MECHANICAL REASON TO EXPECT A LEAN. The Mendelian channel works at markers where the
 * LOADED parent is homozygous, because only there is the transmitted allele known. How many such
 * markers a parent supplies is one minus that parent's own heterozygosity. The donors here are not
 * equally heterozygous, so the two sides do not get the same number of chances. Worse, the corpus
 * has ONE sperm donor and several egg donors, so "paternal" is one person's array quality and
 * "maternal" is a mixture, and any quality difference between them becomes a parental effect that is
 * not biology.
 *
 * TWO ARMS, AND THE FIRST ONE IS REAL.
 *
 *   1. THE PRONUCLEI, SPLIT BY WHICH SIDE IS MISSING. A pronucleus is a natural, complete loss of
 *      one parental complement, and the dissection says which. Seven are maternal and eight are
 *      paternal, so the recovery rate can be computed separately for each. Nothing is constructed
 *      and nothing is injected. If those two rates differ, that difference is the tool's own
 *      asymmetry, measured on real biology.
 *
 *   2. A SYMMETRIC INJECTION ON REAL BIPARENTAL ARRAYS, for the sample size the pronuclei cannot
 *      give. The same chromosome is emptied of the paternal copy in one run and of the maternal copy
 *      in the other, on the same array, using each parent's own homozygous markers so the two
 *      constructions are mirror images. This is CONSTRUCTED and is reported as such: an earlier
 *      harness in this tree built a child's post-loss genotype out of the loaded parent's alleles
 *      and scored a signal that existed whoever was loaded.
 *
 *      THE GUARD AGAINST THAT IS THE UNMODIFIED CHILD, and the first version of this file got the
 *      guard wrong too. It scored against "a donor who is not the father", and since this corpus
 *      has exactly ONE sperm donor that fell back to the father himself, so the control was the
 *      true arm again and passed 12 of 12 while proving nothing. The control that actually
 *      controls is the same array with nothing removed, scored against the same two true parents
 *      on the same chromosome. It has to name nobody.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 \
 *        audit/parental-power.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const one = await import(`${W}oneParentOrigin.ts`)
const both = await import(`${W}bothParentsOrigin.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('parental-power: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

const CONTAMINATED = 'GSM4472408'
const ADO = 0.199
const CHROMS = Array.from({ length: 22 }, (_, i) => String(i + 1))

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  const first = (s: string) => (s ?? '').split(';').filter((x) => x && x !== CONTAMINATED)[0] ?? ''
  return {
    gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
    role: f[col('role')], pronucleus: f[col('pronucleus')],
    father: first(f[col('father_gsms')]), mother: first(f[col('mother_gsms')]),
    complete: f[col('complete_trio')] === 'True',
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Arr { gt: Map<string, AB>; chrom: Map<string, string> }
const cache = new Map<string, Arr>()
function load(gsm: string, keep = false): Arr {
  const hit = cache.get(gsm)
  if (hit) return hit
  const text = gunzipSync(readFileSync(join(DIR, `${gsm}.probes.gz`))).toString('utf8')
  const ls = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  const gt = new Map<string, AB>()
  const chrom = new Map<string, string>()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    gt.set(r.probesetId, r.genotype as AB)
    chrom.set(r.probesetId, r.chrom)
  }
  const a = { gt, chrom }
  if (keep) cache.set(gsm, a)
  return a
}

const donors = recs.filter((r) => r.role === 'parent' && has(r.gsm) && r.gsm !== CONTAMINATED)
const who = (t: string) => (/sperm/i.test(t) ? 'sperm' : (/egg donor ([A-Z])/i.exec(t)?.[1] ?? t))

/** Pairs for one chromosome: [parent genotype, child genotype] at markers both called. */
function pairsOn(parent: Arr, child: Arr, chrom: string): [AB, AB][] {
  const out: [AB, AB][] = []
  for (const [probe, cg] of child.gt) {
    if (child.chrom.get(probe) !== chrom) continue
    const pg = parent.gt.get(probe)
    if (pg) out.push([pg, cg])
  }
  return out
}

// ============================================================ 0. HOW MANY CHANCES EACH SIDE GETS
//
// The count of markers where a parent is homozygous and the child is called. That is the
// denominator of the whole channel, and if the two sides do not get the same one they do not have
// the same power however good the statistics are.
console.log('=== 0. INFORMATIVE MARKERS PER PARENT, which is the power each side actually has')
console.log('    A marker counts for a parent only where THAT parent is homozygous, so this is')
console.log('    driven by each donor\'s own heterozygosity and is not the same on the two sides.')
{
  const trios = recs.filter((r) => r.complete && r.role !== 'parent' && has(r.gsm)
    && has(r.father) && has(r.mother) && r.mother !== CONTAMINATED)
    .slice(0, Number(process.env.OM_N ?? 1000))
  const skews: number[] = []
  const perDonor = new Map<string, number[]>()
  for (const t of trios) {
    const child = load(t.gsm)
    const f = load(t.father, true)
    const m = load(t.mother, true)
    let fi = 0
    let mi = 0
    for (const [probe, cg] of child.gt) {
      if (cg === 'NC') continue
      const fg = f.gt.get(probe)
      const mg = m.gt.get(probe)
      if (fg === 'AA' || fg === 'BB') fi += 1
      if (mg === 'AA' || mg === 'BB') mi += 1
    }
    const skew = mi ? fi / mi : NaN
    if (Number.isFinite(skew)) {
      skews.push(skew)
      const key = who(recs.find((r) => r.gsm === t.mother)?.title ?? '?')
      perDonor.set(key, [...(perDonor.get(key) ?? []), skew])
    }
  }
  skews.sort((a, b) => a - b)
  const q = (p: number) => skews[Math.min(skews.length - 1, Math.floor(p * skews.length))]
  console.log(`    trios ${skews.length}: paternal informative / maternal informative`)
  console.log(`      min ${q(0).toFixed(3)}  p50 ${q(0.5).toFixed(3)}  max ${q(1).toFixed(3)}`)
  for (const [d, xs] of [...perDonor].sort()) {
    const s = [...xs].sort((a, b) => a - b)
    console.log(`      against egg donor ${d}: n=${xs.length} median `
      + `${s[s.length >> 1].toFixed(3)}`)
  }
  console.log('    A ratio of 1.000 is equal power. Away from 1 in either direction, the side with')
  console.log('    FEWER informative markers has less chance of being seen, whatever the biology.')
}

// ============================================================ 1. THE PRONUCLEI, SPLIT BY SIDE
//
// Real, complete, natural losses of one parental complement, with the side known from dissection.
// This is the measurement that matters: no construction, no injection, no assumption.
console.log('')
console.log('=== 1. REAL LOSSES: the pronuclei, recovery split by WHICH parent is missing')
console.log('    A maternal pronucleus is missing the PATERNAL complement and the reverse. If the')
console.log('    two recovery rates differ, that gap is the tool\'s own asymmetry on real biology.')
{
  function resolveDonor(childGsm: string): { id: string; who: string } | null {
    const c = load(childGsm)
    const scored = donors.map((d) => {
      const p = load(d.gsm, true)
      let opp = 0
      let n = 0
      for (const [probe, cg] of c.gt) {
        if (cg !== 'AA' && cg !== 'BB') continue
        const pg = p.gt.get(probe)
        if (pg !== 'AA' && pg !== 'BB') continue
        n += 1
        if (pg !== cg) opp += 1
      }
      return { id: d.gsm, who: who(d.title), rate: n ? opp / n : Infinity }
    }).sort((a, b) => a.rate - b.rate)
    const best = scored[0]
    const next = scored.find((x) => x.who !== best.who)
    if (!next || !(next.rate > best.rate * 2)) return null
    return { id: best.id, who: best.who }
  }

  const pronuclei = recs.filter((r) => r.pronucleus && has(r.gsm))
  const bySide: Record<string, number[]> = { 'paternal missing': [], 'maternal missing': [] }
  let unresolved = 0
  for (const p of pronuclei) {
    const res = resolveDonor(p.gsm)
    if (!res) { unresolved += 1; continue }
    // A MATERNAL pronucleus carries the egg donor, so the PATERNAL complement is what is gone.
    const missing = p.pronucleus === 'maternal' ? 'paternal missing' : 'maternal missing'
    // The array whose complement is absent: load it, and its own copy should read as lost.
    const absentId = p.pronucleus === 'maternal'
      ? donors.find((d) => who(d.title) === 'sperm')!.gsm
      : donors.find((d) => who(d.title) !== 'sperm')!.gsm
    const child = load(p.gsm)
    const absent = load(absentId, true)
    let right = 0
    let n = 0
    for (const c of CHROMS) {
      const pr = pairsOn(absent, child, c)
      if (pr.length < 200) continue
      const call = one.callOneParentOrigin(pr as never, ADO) as { verdict: string }
      n += 1
      if (call.verdict === 'known-parent-lost') right += 1
    }
    if (n) bySide[missing].push(right / n)
  }
  const stat = (xs: number[]) => {
    if (!xs.length) return { mean: NaN, lo: NaN, hi: NaN }
    const mean = xs.reduce((a, x) => a + x, 0) / xs.length
    if (xs.length < 2) return { mean, lo: NaN, hi: NaN }
    const sd = Math.sqrt(xs.reduce((a, x) => a + (x - mean) ** 2, 0) / (xs.length - 1))
    const se = sd / Math.sqrt(xs.length)
    return { mean, lo: Math.max(0, mean - 1.96 * se), hi: Math.min(1, mean + 1.96 * se) }
  }
  const pat = stat(bySide['paternal missing'])
  const mat = stat(bySide['maternal missing'])
  console.log(`    PATERNAL complement missing  recovered ${pat.mean.toFixed(4)} `
    + `[${pat.lo.toFixed(4)}, ${pat.hi.toFixed(4)}]  over ${bySide['paternal missing'].length} arrays`)
  console.log(`    MATERNAL complement missing  recovered ${mat.mean.toFixed(4)} `
    + `[${mat.lo.toFixed(4)}, ${mat.hi.toFixed(4)}]  over ${bySide['maternal missing'].length} arrays`)
  console.log(`    unresolved by linkage: ${unresolved}`)
  const overlap = Number.isFinite(pat.lo) && Number.isFinite(mat.lo)
    && pat.lo <= mat.hi && mat.lo <= pat.hi
  console.log(`    intervals ${overlap ? 'OVERLAP, so no asymmetry is demonstrated at this n'
    : 'DO NOT OVERLAP: the tool sees one side better than the other'}`)
}

// ============================================================ 2. SYMMETRIC INJECTION, higher n
//
// Constructed, and labelled constructed. One chromosome of a real biparental array is emptied of
// the paternal copy in one run and the maternal copy in the other, built the same way on both
// sides so the only difference is which parent is removed.
console.log('')
console.log('=== 2. CONSTRUCTED, MIRROR-IMAGE LOSSES on real biparental arrays')
console.log('    Same array, same chromosome, same construction, only the side differs. The')
console.log('    wrong-parent arm has to collapse or the signal is the construction, not the loss.')
{
  const trios = recs.filter((r) => r.complete && r.role !== 'parent' && has(r.gsm)
    && has(r.father) && has(r.mother) && r.mother !== CONTAMINATED)
    .slice(0, Number(process.env.OM_INJECT_N ?? 12))
  const TARGET = '7'
  let patNamed = 0
  let matNamed = 0
  let patWrong = 0
  let matWrong = 0
  let runs = 0
  let strangerNamed = 0

  for (const t of trios) {
    const child = load(t.gsm)
    const f = load(t.father, true)
    const m = load(t.mother, true)
    // Build the two mirrored children. Only markers where BOTH parents are homozygous are used, so
    // the surviving genotype is known exactly and the two constructions are the same shape.
    const patLost = new Map<string, AB>(child.gt)
    const matLost = new Map<string, AB>(child.gt)
    let built = 0
    for (const [probe, cg] of child.gt) {
      if (child.chrom.get(probe) !== TARGET || cg === 'NC') continue
      const fg = f.gt.get(probe)
      const mg = m.gt.get(probe)
      if ((fg !== 'AA' && fg !== 'BB') || (mg !== 'AA' && mg !== 'BB')) continue
      // The paternal copy is gone: only the mother's transmitted allele survives.
      patLost.set(probe, mg)
      // And the mirror.
      matLost.set(probe, fg)
      built += 1
    }
    if (built < 200) continue
    runs += 1
    const armFor = (gt: Map<string, AB>, patId: string, matId: string) => {
      const kid = { gt, chrom: child.chrom }
      const pp = pairsOn(load(patId, true), kid, TARGET)
      const mp = pairsOn(load(matId, true), kid, TARGET)
      if (pp.length < 200 || mp.length < 200) return null
      return both.callBothParentsOrigin(
        one.callOneParentOrigin(pp as never, ADO), one.callOneParentOrigin(mp as never, ADO),
      ) as { parent: string | null }
    }
    const a = armFor(patLost, t.father, t.mother)
    if (a?.parent === 'paternal') patNamed += 1
    else if (a?.parent) patWrong += 1
    const b = armFor(matLost, t.father, t.mother)
    if (b?.parent === 'maternal') matNamed += 1
    else if (b?.parent) matWrong += 1
    // THE NEGATIVE CONTROL, AND THE FIRST VERSION OF IT WAS NOT ONE. It picked "a donor who is not
    // the father", and since this corpus has exactly one sperm donor that fell back to the father
    // himself, so the control was the true arm again and reported 12 of 12. The construction is
    // what has to be controlled for, so the control is the UNMODIFIED child on the same chromosome
    // against the same two true parents. It must name nobody: chr7 is intact in these embryos, so
    // anything named there is the harness and not the loss.
    const c = armFor(child.gt, t.father, t.mother)
    if (c?.parent) strangerNamed += 1
  }
  console.log(`    arrays built: ${runs}, chromosome ${TARGET}`)
  console.log(`      paternal copy removed -> named PATERNAL  ${patNamed}/${runs}   `
    + `wrong parent ${patWrong}`)
  console.log(`      maternal copy removed -> named MATERNAL  ${matNamed}/${runs}   `
    + `wrong parent ${matWrong}`)
  console.log(`      NEGATIVE CONTROL, unmodified chr7: named anything ${strangerNamed}/${runs}`)
  console.log(`    ${patNamed === matNamed
    ? 'The two sides recovered EQUALLY on this construction.'
    : `The two sides differ by ${Math.abs(patNamed - matNamed)} of ${runs}.`}`)
}

console.log('')
console.log('parental-power: an unequal detector answers a question about inequality by '
  + 'construction. These are the numbers that say whether this one is equal.')

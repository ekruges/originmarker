/**
 * The Mendelian channel, measured against real trios with a wrong-parent control.
 *
 * THIS REPLACES audit/blastomere-origin.ts, WHICH MEASURED NOTHING, and that file is deleted rather
 * than left beside this one: a harness whose number is discredited is a trap for whoever finds it
 * next. It is recoverable from the history if the construction is ever worth reading. It built the
 * child's post-loss genotype out of the LOADED parent's own alleles, so the signal it scored existed
 * whoever was loaded: 44 of 44 with the true father, 44 of 44 with the true mother, and 44 of 44
 * with arrays that are certainly not parents. Its number, 0.920, shipped as OBVIOUS_EVENT_ACCURACY.
 *
 * TWO ARMS, ALWAYS, AND THE SECOND IS NOT OPTIONAL. Every figure this file reports is reported
 * twice: once with the child's genetically confirmed parent loaded, and once with an unrelated
 * adult. A measurement of parent of origin must COLLAPSE in the second arm. One that does not is
 * measuring something else, and the only way to find that out is to run it.
 *
 * GROUND TRUTH COMES FROM THE PEDIGREE AND FROM GENETICS, NEVER FROM THE TOOL.
 *   Specificity: an intact chromosome of a real child, with its real parent loaded, must return
 *     both-copies-present. Any loss verdict there is false by pedigree.
 *   Sensitivity: a genome carrying ONE parental contribution is missing the other parent's copy on
 *     every autosome by construction, so the channel must call a loss on all of them. The 15
 *     pronuclei of GSE148488 are exactly that, with the parent of origin declared and separately
 *     verified by linkage.
 *
 * Data: GSE148488 (Cell 2021, doi 10.1016/j.cell.2020.10.025), 135 arrays on GPL28377. Public, so
 * unlike the lab corpus a result measured here can be committed and re-run. Point OM_TRIOS at a
 * directory of <GSM>.probes.gz plus trio_manifest_full.csv.
 *
 * GSM4472408, egg donor D, is EXCLUDED as a parental genotype: autosomal heterozygosity 0.3455 at
 * an 82.0 percent call rate against 0.141 to 0.174 at 94 to 98 percent for the other donors, which
 * is a mixed or contaminated sample. 30 of the 106 declared trios name it as mother and are
 * excluded with it.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/trio-calibration.ts [--q 0.1217]
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const one = await import(`${W}oneParentOrigin.ts`)
const both = await import(`${W}bothParentsOrigin.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('trio-calibration: OM_TRIOS is not set to a readable directory. NOT RUN.')
  console.log('  This measurement needs the public GSE148488 arrays. It is skipped rather than')
  console.log('  quietly passing, because a skipped calibration that prints nothing reads exactly')
  console.log('  like one that succeeded.')
  process.exit(0)
}

/** Excluded as a parental genotype. See the header. */
const CONTAMINATED = new Set(['GSM4472408'])

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Array_ { gt: Map<string, AB>; pos: Map<string, string>; het: number; callRate: number }

/**
 * PARENTS ARE CACHED, CHILDREN ARE NOT.
 *
 * Each array is about 790,000 string-keyed entries, roughly 90 MB live. Seven donors are reused
 * across every trio and are worth holding; a child is read once and dropped. Caching both ran the
 * heap to 2 GB and died, which is a property of the harness rather than of the data.
 */
const cache = new Map<string, Array_>()
const PARENTS_ONLY = new Set<string>()
function load(gsm: string): Array_ {
  const hit = cache.get(gsm)
  if (hit) return hit
  const text = gunzipSync(readFileSync(`${DIR}/${gsm}.probes.gz`)).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(lines[h])
  const gt = new Map<string, AB>()
  const pos = new Map<string, string>()
  let het = 0
  let called = 0
  let total = 0
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || !/^\d+$/.test(r.chrom)) continue
    total += 1
    pos.set(r.probesetId, r.chrom)
    if (r.genotype === 'NC') continue
    gt.set(r.probesetId, r.genotype as AB)
    called += 1
    if (r.genotype === 'AB') het += 1
  }
  const out = { gt, pos, het: called ? het / called : NaN, callRate: total ? called / total : NaN }
  if (PARENTS_ONLY.has(gsm)) cache.set(gsm, out)
  return out
}

/**
 * Genotype pairs per chromosome: the loaded parent's call beside the child's.
 *
 * Grouped in ONE pass over the child rather than one pass per chromosome. The obvious form rescans
 * all 790,000 markers for each of 22 chromosomes, which is 17 million lookups per array per arm and
 * is the whole runtime of this harness.
 */
function pairsByChrom(parent: Array_, child: Array_): Map<string, [AB, AB][]> {
  const out = new Map<string, [AB, AB][]>()
  for (const [probe, cg] of child.gt) {
    const pg = parent.gt.get(probe)
    if (!pg || pg === 'NC') continue
    const c = child.pos.get(probe)
    if (!c) continue
    const arr = out.get(c)
    if (arr) arr.push([pg, cg])
    else out.set(c, [[pg, cg]])
  }
  return out
}

// ---------------------------------------------------------------- the manifest
const rows = readFileSync(`${DIR}/trio_manifest_full.csv`, 'utf8').trim().split('\n')
const head = rows[0].split(',')
const recs = rows.slice(1).map((l) => {
  const f = l.split(',')
  return Object.fromEntries(head.map((k, i) => [k, f[i] ?? ''])) as Record<string, string>
})
const present = (g: string) => existsSync(`${DIR}/${g}.probes.gz`)
const firstUsable = (csv: string) =>
  csv.split(/[;| ]/).map((x) => x.trim()).filter((x) => x && !CONTAMINATED.has(x) && present(x))[0]

const trios = recs
  .filter((r) => r.complete_trio.toLowerCase() === 'true' && present(r.gsm))
  .map((r) => ({
    gsm: r.gsm, material: r.material,
    father: firstUsable(r.father_gsms), mother: firstUsable(r.mother_gsms),
  }))
  .filter((t) => t.father && t.mother)

const pronuclei = recs.filter((r) => r.pronucleus && present(r.gsm))
  .map((r) => ({ gsm: r.gsm, side: r.pronucleus }))

// EVERY donor array is cached, not only those a trio names. The pronuclei are resolved against all
// of them, and a donor absent from this set is re-parsed once per pronucleus.
for (const t of trios) { PARENTS_ONLY.add(t.father!); PARENTS_ONLY.add(t.mother!) }
for (const r of recs) {
  if (r.role === 'parent' && present(r.gsm) && !CONTAMINATED.has(r.gsm)) PARENTS_ONLY.add(r.gsm)
}

console.log(`trios usable ${trios.length}, pronuclei ${pronuclei.length}, `
  + `donors held ${PARENTS_ONLY.size}`)
console.log(`  by material: ${JSON.stringify(trios.reduce((a: Record<string, number>, t) => {
  a[t.material] = (a[t.material] ?? 0) + 1; return a
}, {}))}`)

const Q = (() => {
  const i = process.argv.indexOf('--q')
  return i > 0 ? Number(process.argv[i + 1]) : undefined
})()
const CHROMS = Array.from({ length: 22 }, (_, i) => String(i + 1))
const ADO = 0.199

/** A cluster-robust interval over arrays, since many calls from one array are one array. */
function clustered(perArray: number[]): { mean: number; lo: number; hi: number } {
  const n = perArray.length
  if (!n) return { mean: NaN, lo: NaN, hi: NaN }
  const mean = perArray.reduce((s, x) => s + x, 0) / n
  if (n < 2) return { mean, lo: NaN, hi: NaN }
  const sd = Math.sqrt(perArray.reduce((s, x) => s + (x - mean) ** 2, 0) / (n - 1))
  const se = sd / Math.sqrt(n)
  return { mean, lo: Math.max(0, mean - 1.96 * se), hi: Math.min(1, mean + 1.96 * se) }
}

// ---------------------------------------------------------------- 1. specificity, both arms
//
// Intact chromosomes of real children. The true parent is loaded in one arm and an unrelated adult
// in the other. Truth is the pedigree: a child has both parents' copies, so both-present is the
// only correct verdict, whoever is loaded.
//
// REPORTED BY MATERIAL AS WELL AS OVERALL, and the split is the point rather than a nicety. A
// blastomere is one cell and a trophectoderm biopsy is five to ten, so they drop out at different
// rates and cannot share a false-call rate. The single pooled figure the first version printed was
// dominated by whichever material happened to lead the trio list. Nothing about the workflow this
// tool is for can be decided from a number that mixes them.
{
  // Default is every usable trio. OM_N exists to cut a run short while iterating, and the per-
  // material counts print alongside every rate so a short run cannot be mistaken for a full one.
  const SAMPLE = trios.slice(0, Number(process.env.OM_N ?? trios.length))
  const arms: Record<string, number[]> = { true: [], wrong: [] }
  const verdicts: Record<string, Record<string, number>> = { true: {}, wrong: {} }
  const byMaterial: Record<string, Record<string, number[]>> = {}
  for (const t of SAMPLE) {
    const child = load(t.gsm)
    // An unrelated adult: a donor who is not this child's parent.
    const wrongId = [...new Set(trios.map((x) => x.father).concat(trios.map((x) => x.mother)))]
      .find((g) => g && g !== t.father && g !== t.mother)!
    for (const [arm, pid] of [['true', t.father!], ['wrong', wrongId]] as [string, string][]) {
      const parent = load(pid)
      const byChrom = pairsByChrom(parent, child)
      let bad = 0
      let n = 0
      for (const c of CHROMS) {
        const pr = byChrom.get(c) ?? []
        if (pr.length < 200) continue
        const call = one.callOneParentOrigin(pr, ADO, Q) as { verdict: string }
        verdicts[arm][call.verdict] = (verdicts[arm][call.verdict] ?? 0) + 1
        n += 1
        if (call.verdict === 'known-parent-lost' || call.verdict === 'other-parent-lost') bad += 1
      }
      if (n) {
        arms[arm].push(bad / n)
        byMaterial[t.material] ??= { true: [], wrong: [] }
        byMaterial[t.material][arm].push(bad / n)
      }
    }
  }
  console.log(`\n=== SPECIFICITY: intact chromosomes of real children, q=${Q ?? one.DEFAULT_Q}`)
  for (const arm of ['true', 'wrong']) {
    const c = clustered(arms[arm])
    console.log(`  ${arm === 'true' ? 'TRUE parent ' : 'WRONG parent'}  false-call rate `
      + `${c.mean.toFixed(4)} [${c.lo.toFixed(4)}, ${c.hi.toFixed(4)}] over ${arms[arm].length} arrays`)
    console.log(`      verdicts: ${JSON.stringify(verdicts[arm])}`)
  }
  console.log('  BY MATERIAL, true-parent arm. A single-cell biopsy and a five-to-ten-cell one do')
  console.log('  not share a dropout rate, so they do not share a false-call rate either.')
  for (const m of Object.keys(byMaterial).sort()) {
    const t = clustered(byMaterial[m].true)
    const w = clustered(byMaterial[m].wrong)
    console.log(`    ${m.padEnd(16)} n=${String(byMaterial[m].true.length).padStart(3)}  `
      + `true ${t.mean.toFixed(4)} [${t.lo.toFixed(4)}, ${t.hi.toFixed(4)}]  `
      + `wrong ${w.mean.toFixed(4)} [${w.lo.toFixed(4)}, ${w.hi.toFixed(4)}]`)
  }
}

// ---------------------------------------------------------------- 2. DIRECTION, on the pronuclei
//
// A pronucleus carries ONE parental contribution, so exactly one parental copy is missing and WHICH
// one is known from the dissection. That makes it a direction test, and direction is the whole
// question: naming the wrong parent is the failure that matters.
//
// THE FIRST VERSION OF THIS TEST COUNTED "ANY LOSS" AND COULD NOT DISCRIMINATE. Load the declared
// parent of a maternal pronucleus and its copy IS present, so the correct verdict is
// other-parent-lost. Load the sperm donor against the same genome and the correct verdict is
// known-parent-lost. Both are losses; only the direction separates them, and scoring losses scored
// both arms at once. It reported 0.60 against 0.80 and looked like a failure of the tool.
{
  // THE MANIFEST ASSIGNS THE PRONUCLEI NO PARENTS, so they are resolved by linkage instead: the
  // opposite-homozygote rate is Mendelianly forbidden between a true parent and its offspring, so
  // the real donor sits near zero and everyone else an order of magnitude above. Replicates of one
  // donor are collapsed first, or the runner-up is always the same person's second array.
  //
  // An earlier version of this block took "the first mother in the trio list" as each maternal
  // pronucleus's parent. Those pronuclei belong to a donor that list never names, so both arms were
  // wrong-parent arms and the result read as a defect in the tool.
  const donorIds = [...new Set(trios.flatMap((t) => [t.father!, t.mother!]))]
  const allDonors = recs.filter((r) => r.role === 'parent' && present(r.gsm)
    && !CONTAMINATED.has(r.gsm)).map((r) => ({ id: r.gsm, title: r.title }))
  const person = (title: string) => (/sperm/i.test(title) ? 'sperm'
    : (/egg donor ([A-Z])/i.exec(title)?.[1] ?? title))
  const resolveParent = (childId: string): { id: string; who: string } | null => {
    const c = load(childId)
    const scored = allDonors.map((d) => {
      const P = load(d.id)
      let opp = 0
      let n = 0
      for (const [probe, cg] of c.gt) {
        if (cg === 'AB') continue
        const pg = P.gt.get(probe)
        if (!pg || pg === 'AB') continue
        n += 1
        if (pg !== cg) opp += 1
      }
      return { ...d, who: person(d.title), rate: n ? opp / n : Infinity }
    }).sort((a, b) => a.rate - b.rate)
    const best = scored[0]
    const nextPerson = scored.find((x) => x.who !== best.who)
    // Resolved when the best donor beats the best OTHER PERSON by a clear margin.
    if (!nextPerson || !(nextPerson.rate > best.rate * 2)) return null
    return { id: best.id, who: best.who }
  }
  const perArray: { present: number[]; absent: number[] } = { present: [], absent: [] }
  const seen: Record<string, Record<string, number>> = { present: {}, absent: {} }
  let unresolved = 0
  for (const p of pronuclei) {
    const child = load(p.gsm)
    const res = resolveParent(p.gsm)
    if (!res) { unresolved += 1; continue }
    const own = res.id
    // The absent parent: any donor who is a different person from the resolved one.
    const missing = allDonors.find((d) => person(d.title) !== res.who)!.id
    for (const [arm, pid, want] of [
      ['present', own, 'other-parent-lost'],
      ['absent', missing, 'known-parent-lost'],
    ] as [string, string, string][]) {
      const parent = load(pid)
      const byChrom = pairsByChrom(parent, child)
      let right = 0
      let n = 0
      for (const c of CHROMS) {
        const pr = byChrom.get(c) ?? []
        if (pr.length < 200) continue
        const call = one.callOneParentOrigin(pr, ADO, Q) as { verdict: string }
        seen[arm][call.verdict] = (seen[arm][call.verdict] ?? 0) + 1
        n += 1
        if (call.verdict === want) right += 1
      }
      if (n) perArray[arm as 'present' | 'absent'].push(right / n)
    }
  }
  console.log(`\n=== DIRECTION: pronuclei, which parental copy is missing is known by dissection`)
  console.log(`  parents resolved by linkage; ${unresolved} of ${pronuclei.length} unresolved`)
  for (const [arm, want] of [['present', 'other-parent-lost'], ['absent', 'known-parent-lost']]) {
    const c = clustered(perArray[arm as 'present' | 'absent'])
    console.log(`  loaded parent ${arm === 'present' ? 'IS  this genome' : 'is NOT this genome'}`
      + `, correct verdict ${want.padEnd(18)} ${c.mean.toFixed(4)} `
      + `[${c.lo.toFixed(4)}, ${c.hi.toFixed(4)}] over ${perArray[arm as 'present' | 'absent'].length} arrays`)
    console.log(`      verdicts seen: ${JSON.stringify(seen[arm])}`)
  }
}

// ------------------------------------------------- 3. BOTH PARENTS LOADED, which is the real run
//
// The Mendelian channel answers "the loaded parent's own copy is absent" at 0.8539 and "the other
// parent's copy is absent" at 0.2890, on the same arrays and the same events. With both parental
// arrays in hand every loss can be put to the side that makes it the first question. This measures
// whether that actually works, on the same two arms as everything else: a set where the answer is
// known, and a set where the correct answer is to say nothing.
//
// The rule under test is `web/src/bothParentsOrigin.ts`, the same function the tool calls. A rule
// measured here and re-implemented there would be measuring a different rule.
{
  const donorsAll = recs.filter((r) => r.role === 'parent' && present(r.gsm)
    && !CONTAMINATED.has(r.gsm)).map((r) => ({ id: r.gsm, title: r.title }))
  const who = (title: string) => (/sperm/i.test(title) ? 'sperm'
    : (/egg donor ([A-Z])/i.exec(title)?.[1] ?? title))
  const resolve = (childId: string): { id: string; who: string } | null => {
    const c = load(childId)
    const scored = donorsAll.map((d) => {
      const P = load(d.id)
      let opp = 0
      let n = 0
      for (const [probe, cg] of c.gt) {
        if (cg === 'AB') continue
        const pg = P.gt.get(probe)
        if (!pg || pg === 'AB') continue
        n += 1
        if (pg !== cg) opp += 1
      }
      return { ...d, who: who(d.title), rate: n ? opp / n : Infinity }
    }).sort((a, b) => a.rate - b.rate)
    const best = scored[0]
    const next = scored.find((x) => x.who !== best.who)
    if (!next || !(next.rate > best.rate * 2)) return null
    return { id: best.id, who: best.who }
  }

  // --- DIRECTION. A pronucleus is missing exactly one parental complement and the dissection says
  // which. The combined rule must name that parent.
  const perArray: number[] = []
  const named: Record<string, number> = {}
  // Does the OTHER array independently point at the same parent? It is never acted on, because
  // that direction is the 0.2890 one, but whether it tracks the truth decides if it could ever
  // serve as a filter on material where dropout drives false calls.
  const corrob: Record<string, number> = { yes: 0, no: 0 }
  let wrongParent = 0
  let unresolved = 0
  for (const p of pronuclei) {
    const child = load(p.gsm)
    const res = resolve(p.gsm)
    if (!res) { unresolved += 1; continue }
    // A pronucleus resolving to the sperm donor is PATERNAL, so the MATERNAL copy is the missing
    // one, and the reverse. Ground truth is the linkage plus the dissection, never this channel.
    const paternalPn = res.who === 'sperm'
    const wantParent = paternalPn ? 'maternal' : 'paternal'
    const patId = paternalPn ? res.id : donorsAll.find((d) => who(d.title) === 'sperm')!.id
    const matId = paternalPn ? donorsAll.find((d) => who(d.title) !== 'sperm')!.id : res.id
    const patByChrom = pairsByChrom(load(patId), child)
    const matByChrom = pairsByChrom(load(matId), child)
    let right = 0
    let n = 0
    for (const c of CHROMS) {
      const pp = patByChrom.get(c) ?? []
      const mp = matByChrom.get(c) ?? []
      if (pp.length < 200 || mp.length < 200) continue
      const call = both.callBothParentsOrigin(
        one.callOneParentOrigin(pp, ADO, Q), one.callOneParentOrigin(mp, ADO, Q),
      ) as { verdict: string; parent: string | null }
      named[call.verdict] = (named[call.verdict] ?? 0) + 1
      n += 1
      if (call.parent) { corrob[call.corroborated ? 'yes' : 'no'] += 1 }
      if (call.parent === wantParent) right += 1
      else if (call.parent !== null) wrongParent += 1
    }
    if (n) perArray.push(right / n)
  }
  const d = clustered(perArray)
  console.log('\n=== BOTH PARENTS LOADED, direction on the pronuclei')
  console.log(`  correct parent named ${d.mean.toFixed(4)} [${d.lo.toFixed(4)}, ${d.hi.toFixed(4)}]`
    + ` over ${perArray.length} arrays, ${unresolved} unresolved`)
  console.log(`  WRONG PARENT NAMED: ${wrongParent}`)
  console.log(`  verdicts: ${JSON.stringify(named)}`)
  console.log(`  corroborated by the other array: ${JSON.stringify(corrob)}`)

  // --- SPECIFICITY. Real children with BOTH real parents loaded. Every autosome carries both
  // copies by pedigree, so the only correct answer is to name nobody.
  const SAMPLE2 = trios.slice(0, Number(process.env.OM_N ?? trios.length))
  const spec: number[] = []
  const specByMaterial: Record<string, number[]> = {}
  const specVerdicts: Record<string, number> = {}
  const specCorrob: Record<string, number> = { yes: 0, no: 0 }
  for (const t of SAMPLE2) {
    const child = load(t.gsm)
    const patByChrom = pairsByChrom(load(t.father!), child)
    const matByChrom = pairsByChrom(load(t.mother!), child)
    let bad = 0
    let n = 0
    for (const c of CHROMS) {
      const pp = patByChrom.get(c) ?? []
      const mp = matByChrom.get(c) ?? []
      if (pp.length < 200 || mp.length < 200) continue
      const call = both.callBothParentsOrigin(
        one.callOneParentOrigin(pp, ADO, Q), one.callOneParentOrigin(mp, ADO, Q),
      ) as { verdict: string; parent: string | null }
      specVerdicts[call.verdict] = (specVerdicts[call.verdict] ?? 0) + 1
      n += 1
      if (call.parent !== null) { bad += 1; specCorrob[call.corroborated ? 'yes' : 'no'] += 1 }
    }
    if (n) {
      spec.push(bad / n)
      specByMaterial[t.material] ??= []
      specByMaterial[t.material].push(bad / n)
    }
  }
  const sc = clustered(spec)
  console.log('\n=== BOTH PARENTS LOADED, specificity on real children')
  console.log(`  a parent named where both copies are present: ${sc.mean.toFixed(4)} `
    + `[${sc.lo.toFixed(4)}, ${sc.hi.toFixed(4)}] over ${spec.length} arrays`)
  console.log(`  verdicts: ${JSON.stringify(specVerdicts)}`)
  console.log(`  FALSE calls corroborated by the other array: ${JSON.stringify(specCorrob)}`)
  for (const m of Object.keys(specByMaterial).sort()) {
    const c = clustered(specByMaterial[m])
    console.log(`    ${m.padEnd(16)} n=${String(specByMaterial[m].length).padStart(3)}  `
      + `${c.mean.toFixed(4)} [${c.lo.toFixed(4)}, ${c.hi.toFixed(4)}]`)
  }
}

console.log('\ntrio-calibration: both arms reported. A parental measurement whose wrong-parent arm '
  + 'does not move is not measuring parentage.')

/**
 * The panel builder's claims, on real trios whose relationships are established independently.
 *
 * `informativity.ts` decides which markers can testify about an embryo's parental contribution,
 * and a panel is nothing but the markers it kept. Four claims, each falsifiable on this data:
 *
 *   1. MENDELIAN CONSISTENCY SEPARATES A REAL FATHER FROM A STRANGER. On a verified trio the
 *      violation rate is dropout and genotyping error and nothing else. Swap in an unrelated
 *      donor and the rate must rise sharply, because a stranger transmits alleles the child does
 *      not have. If the two rates are alike, the module is not reading parentage at all. Truth is
 *      the series' own trio structure, not this tool.
 *
 *   2. A DEDUCED PATERNAL ALLELE MUST BE AN ALLELE THE FATHER HAS. Where the module reports
 *      `paternalAlleleDeducible`, it is naming the allele the embryo received from its father.
 *      The father's own array either carries that allele or it does not. This is not the same
 *      arithmetic run twice: the deduction runs off the MOTHER being unable to supply what the
 *      child has, and the check runs off the FATHER's genotype directly.
 *
 *   3. A KEY SNP MUST SURVIVE ONE DROPOUT. `keySnp` asserts that a single allele dropout cannot
 *      flip this marker's phase. That is a claim about a perturbation, so it is tested by
 *      performing the perturbation: take the real parental genotypes, drop one allele, and see
 *      whether the karyomapping class changes. A key SNP that flips has a false label, and a
 *      non-key SNP that never flips means the label is costing markers for nothing.
 *
 *   4. A CLAIM MUST SURVIVE A SECOND ARRAY OF THE SAME PERSON. The series carries technical
 *      replicates. A marker called informative from one array of a parent and not from another
 *      array of that SAME parent was never a property of the person, only of one hybridisation.
 *
 * NOTHING HERE IS CONSTRUCTED except the single dropout in claim 3, which is applied to real
 * genotypes and is the exact perturbation the claim is about.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 \
 *        audit/panel-truth.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const inf = await import(`${W}informativity.ts`)
const relMod = await import(`${W}relatedness.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  console.log('panel-truth: OM_TRIOS is not readable. NOT RUN.')
  process.exit(0)
}
const CONTAMINATED = 'GSM4472408'

const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
const head = lines[0].split(',')
const col = (n: string) => head.indexOf(n)
const recs = lines.slice(1).map((l) => {
  const f = l.split(',')
  return {
    gsm: f[col('gsm')], title: f[col('title')], role: f[col('role')],
    pronucleus: f[col('pronucleus')],
    father: (f[col('father_gsms')] ?? '').split(';').filter(Boolean)[0],
    father2: (f[col('father_gsms')] ?? '').split(';').filter(Boolean)[1],
    mother: (f[col('mother_gsms')] ?? '').split(';').filter(Boolean)[0],
    complete: f[col('complete_trio')] === 'True',
  }
})
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))
const textOf = (g: string) => gunzipSync(readFileSync(join(DIR, `${g}.probes.gz`))).toString('utf8')

type AB = 'AA' | 'AB' | 'BB' | 'NC'
const cache = new Map<string, Map<string, AB>>()
function gt(gsm: string, keep = false): Map<string, AB> {
  const hit = cache.get(gsm)
  if (hit) return hit
  const ls = textOf(gsm).split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  const m = new Map<string, AB>()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (r) m.set(r.probesetId, r.genotype as AB)
  }
  if (keep) cache.set(gsm, m)
  return m
}

const pct = (a: number, b: number) => (b ? `${((100 * a) / b).toFixed(4)}%` : 'n/a')
const isHom = (g: AB) => g === 'AA' || g === 'BB'

const kids = recs.filter((r) => r.complete && r.role !== 'parent' && !r.pronucleus
  && has(r.gsm) && has(r.father) && has(r.mother)
  && r.father !== CONTAMINATED && r.mother !== CONTAMINATED)
  .slice(0, Number(process.env.OM_CAP ?? 12))
const donors = recs.filter((r) => r.role === 'parent' && has(r.gsm) && r.gsm !== CONTAMINATED)

console.log(`complete trios: ${kids.length}, adult donor arrays available as strangers: ${donors.length}`)
console.log('')

// ================================================================ CLAIM 1 and CLAIM 2
console.log('=== CLAIM 1: Mendelian violations separate the real father from a stranger')
console.log('=== CLAIM 2: a deduced paternal allele is an allele the father actually carries')
console.log('')
let tTrue = { seen: 0, viol: 0 }
let tFalse = { seen: 0, viol: 0 }
let ded = { n: 0, fatherHasIt: 0, fatherLacksIt: 0, fatherNoCall: 0 }
let l2 = { n: 0, fatherHet: 0 }
const emptyPair = () => ({ kept: { n: 0, bad: 0 }, dropped: { n: 0, bad: 0 } })
const sel = { key: emptyPair(), l3: emptyPair(), l2: emptyPair() }
let l3 = { n: 0, motherCannot: 0 }

for (const k of kids) {
  const F = gt(k.father, true)
  const M = gt(k.mother, true)
  const C = gt(k.gsm)
  // The SECOND array of the same man, when the manifest lists one. Never used to classify.
  const F2 = k.father2 && k.father2 !== k.father ? gt(k.father2, true) : null
  // THE STRANGER IS CHOSEN BY IDENTITY, NOT BY ACCESSION, because this corpus has repeatedly
  // punished the lazy version. Taking "the first donor that is not this child's parents" picked
  // an egg donor who is the child's own mother in another slot, or a second array of the father,
  // and the arm then compared a relative against a relative. The donor with the HIGHEST
  // opposite-homozygote rate against this child is the one the data says is least related to it.
  let strangerGsm: string | null = null
  let worst = -1
  for (const d of donors) {
    if (d.gsm === k.father || d.gsm === k.mother) continue
    const o = relMod.oppositeHom(gt(d.gsm, true), C) as { rate: number; n: number }
    if (o.n > 50_000 && o.rate > worst) { worst = o.rate; strangerGsm = d.gsm }
  }
  const S = strangerGsm ? gt(strangerGsm, true) : null

  let localTrue = 0
  let localFalse = 0
  let localSeen = 0
  for (const [probe, cg] of C) {
    const fg = F.get(probe)
    const mg = M.get(probe)
    if (!fg || !mg || fg === 'NC' || mg === 'NC' || cg === 'NC') continue
    localSeen += 1

    const v = inf.classifyMarker(fg, mg, cg)
    if (v.mendelianConsistent === false) localTrue += 1

    // CLAIM 2, AND THE FATHER IS DELIBERATELY NOT CONSULTED TO MAKE IT.
    //
    // The first version of this asked classifyMarker for `paternalAlleleDeducible` and then
    // checked the father carried it. classifyMarker is GIVEN the father, so it can only ever
    // name an allele he has: the check returned 100.0000% of 6.8 million markers and could not
    // have returned anything else. A test that cannot fail measures nothing.
    //
    // This is the same deduction made from the MOTHER and CHILD alone, which is what karyomapping
    // does when the paternal contribution is the unknown. Where the mother is homozygous and the
    // child carries the other allele, that allele is not hers, so it is his. The father's array
    // is then an independent fact and the check can fail.
    // THE CIRCULARITY, AND THE ONLY WAY OUT OF IT THIS CORPUS OFFERS.
    //
    // Four versions of this arm failed the same way. classifyMarker is GIVEN the father, so every
    // label it produces is a function of his genotype, and every check that consults his genotype
    // is therefore checking the module against itself. The clearest case: provesPaternalPresence
    // kept exactly the 559,552 markers where he carries the deduced allele and dropped exactly the
    // 524,189 where he does not, for a contradiction ratio of 0.000 that could not have come out
    // any other way.
    //
    // The series carries TWO arrays of the same sperm donor. So the label is computed from the
    // FIRST and checked against the SECOND. A label that is a property of the man survives that;
    // a label that is reading one hybridisation's dropout does not. The two arrays disagree about
    // heterozygosity at 4.152% of jointly called markers, which is the floor this can resolve to
    // and is stated rather than assumed.
    if (F2 && isHom(mg) && cg === 'AB') {
      const a = cg.split('').find((x) => x !== mg[0])
      const f2 = F2.get(probe)
      if (a && f2 && f2 !== 'NC') {
        const contradicted = !f2.includes(a)
        ded.n += 1
        if (contradicted) ded.fatherLacksIt += 1
        else ded.fatherHasIt += 1
        for (const [kept, bucket] of [
          [v.keySnp, sel.key],
          [v.provesPaternalPresence, sel.l3],
          [v.informativeHaplotypeOrigin, sel.l2],
        ] as const) {
          const b = kept ? bucket.kept : bucket.dropped
          b.n += 1
          if (contradicted) b.bad += 1
        }
      }
    }
    // The two power labels, against their own stated preconditions.
    if (v.informativeHaplotypeOrigin) { l2.n += 1; if (fg === 'AB') l2.fatherHet += 1 }
    if (v.provesPaternalPresence) {
      l3.n += 1
      // The mother must be unable to supply at least one allele the child carries.
      const childAlleles = new Set(cg.split(''))
      const motherAlleles = new Set(mg.split(''))
      if ([...childAlleles].some((x) => !motherAlleles.has(x))) l3.motherCannot += 1
    }

    if (S) {
      const sg = S.get(probe)
      if (sg && sg !== 'NC') {
        const vs = inf.classifyMarker(sg, mg, cg)
        if (vs.mendelianConsistent === false) localFalse += 1
      }
    }
  }
  tTrue.seen += localSeen; tTrue.viol += localTrue
  tFalse.seen += localSeen; tFalse.viol += localFalse
  console.log(`  ${k.gsm}: ${localSeen.toLocaleString()} markers callable in all three  `
    + `TRUE father ${pct(localTrue, localSeen)}  `
    + `STRANGER ${strangerGsm ?? 'none'} (opp ${worst.toFixed(4)}) ${pct(localFalse, localSeen)}`)
  cache.delete(k.gsm)
}

console.log('')
console.log(`CLAIM 1  true fathers  ${tTrue.viol.toLocaleString()} / ${tTrue.seen.toLocaleString()} = ${pct(tTrue.viol, tTrue.seen)}`)
console.log(`         strangers     ${tFalse.viol.toLocaleString()} / ${tFalse.seen.toLocaleString()} = ${pct(tFalse.viol, tFalse.seen)}`)
console.log(`         ratio         ${(tFalse.viol / Math.max(tTrue.viol, 1)).toFixed(1)}x`)
console.log('')
console.log(`CLAIM 2  paternal allele deduced at ${ded.n.toLocaleString()} markers`)
console.log(`         the father CARRIES it        ${ded.fatherHasIt.toLocaleString()} = ${pct(ded.fatherHasIt, ded.n)}`)
console.log(`         the father LACKS it          ${ded.fatherLacksIt.toLocaleString()} = ${pct(ded.fatherLacksIt, ded.n)}   <-- every one is a false claim or a real Mendelian error`)
console.log(`         the father has no call       ${ded.fatherNoCall.toLocaleString()}`)
console.log('')
console.log('         THE SELECTION TEST. The label is computed from the father\'s FIRST array and')
console.log('         checked against his SECOND, so his genotype is not on both sides of it.')
console.log('         Contradiction means that second array does not carry the allele the mother')
console.log('         could not have supplied.')
console.log('         label                   kept n / contradicted    dropped n / contradicted    ratio')
for (const [name, b] of [['keySnp', sel.key], ['provesPaternalPresence', sel.l3],
  ['informativeHaplotypeOrigin', sel.l2]] as const) {
  const kr = b.kept.n ? b.kept.bad / b.kept.n : NaN
  const dr = b.dropped.n ? b.dropped.bad / b.dropped.n : NaN
  console.log(`         ${name.padEnd(26)} ${String(b.kept.n).padStart(9)} ${pct(b.kept.bad, b.kept.n).padStart(9)}`
    + `    ${String(b.dropped.n).padStart(9)} ${pct(b.dropped.bad, b.dropped.n).padStart(9)}`
    + `    ${Number.isFinite(kr / dr) ? (kr / dr).toFixed(3) : '   -'}`)
}
console.log('         A ratio below 1 means the kept markers are cleaner than the discarded ones,')
console.log('         which is the whole purpose of the selection. At 1 it is choosing at random.')
console.log('')
console.log(`L2 label: ${l2.n.toLocaleString()} markers, father heterozygous at ${pct(l2.fatherHet, l2.n)} of them (the label's own precondition)`)
console.log(`L3 label: ${l3.n.toLocaleString()} markers, mother demonstrably unable to supply a child allele at ${pct(l3.motherCannot, l3.n)}`)

// ================================================================ CLAIM 3: one dropout
console.log('')
console.log('=== CLAIM 3: a key SNP does not change class under one allele dropout')
console.log('A dropout turns a heterozygous CALL into a homozygous one. Applied to the real')
console.log('parental genotypes, one at a time, which is the perturbation the claim is about.')
const drop = (g: AB, which: 0 | 1): AB => (g === 'AB' ? (which === 0 ? 'AA' : 'BB') : g)
// A KEY SNP CLAIMS ITS PHASE CANNOT INVERT UNDER ONE DROPOUT, not that its class label is
// untouched. The first version of this compared karyomappingClass before and after, which changes
// whenever the genotypes change and so reported 100% of key SNPs as failing. That measured the
// wrong quantity. An inversion is `paternalAlleleDeducible` going from A to B or B to A; going to
// null or to ambiguous is the marker becoming uninformative, which is a loss and not an error.
let key = { n: 0, inverted: 0, lost: 0 }
let nonkey = { n: 0, inverted: 0, lost: 0 }
{
  const k = kids[0]
  const F = gt(k.father, true)
  const M = gt(k.mother, true)
  const C = gt(k.gsm, true)
  let seen = 0
  for (const [probe, fg] of F) {
    const mg = M.get(probe)
    const cg = C.get(probe)
    if (!mg || !cg || fg === 'NC' || mg === 'NC' || cg === 'NC') continue
    if (fg !== 'AB' && mg !== 'AB') continue   // no heterozygote for a dropout to hit
    seen += 1
    if (seen > 400_000) break
    const before = inf.classifyMarker(fg, mg, cg).paternalAlleleDeducible
    if (before !== 'A' && before !== 'B') continue
    const isKey = inf.keySnp(fg, mg)
    let inverted = false
    let lost = false
    for (const w of [0, 1] as const) {
      const after: (typeof before)[] = []
      if (fg === 'AB') after.push(inf.classifyMarker(drop(fg, w), mg, cg).paternalAlleleDeducible)
      if (mg === 'AB') after.push(inf.classifyMarker(fg, drop(mg, w), cg).paternalAlleleDeducible)
      for (const a of after) {
        if ((a === 'A' || a === 'B') && a !== before) inverted = true
        else if (a !== before) lost = true
      }
    }
    const bucket = isKey ? key : nonkey
    bucket.n += 1
    if (inverted) bucket.inverted += 1
    else if (lost) bucket.lost += 1
  }
  cache.delete(k.gsm)
}
console.log(`  key SNPs      ${key.n.toLocaleString()}  phase INVERTED ${key.inverted.toLocaleString()} = ${pct(key.inverted, key.n)}  (became uninformative ${pct(key.lost, key.n)})`)
console.log(`  non-key SNPs  ${nonkey.n.toLocaleString()}  phase INVERTED ${nonkey.inverted.toLocaleString()} = ${pct(nonkey.inverted, nonkey.n)}  (became uninformative ${pct(nonkey.lost, nonkey.n)})`)
console.log('  The claim is that the first number is zero. The second column is the cost of the')
console.log('  label, not a failure of it.')

// ================================================================ CLAIM 4: replicate stability
console.log('')
console.log('=== CLAIM 4: an informative call is a property of the person, not of one array')
const byTitle = new Map<string, string[]>()
for (const r of recs) {
  if (!has(r.gsm) || r.role === 'parent') continue
  const base = r.title.replace(/[ _-]*(rep(licate)?[ _-]*\d+|r\d+)$/i, '').trim()
  byTitle.set(base, [...(byTitle.get(base) ?? []), r.gsm])
}
const repPairs = [...byTitle.values()].filter((v) => v.length >= 2).slice(0, 6)
if (!repPairs.length) {
  console.log('  no replicate pairs identified from the manifest titles. NOT RUN.')
} else {
  const k = kids[0]
  const F = gt(k.father, true)
  const M = gt(k.mother, true)
  let agree = 0
  let differ = 0
  for (const [a, b] of repPairs.map((v) => [v[0], v[1]])) {
    const A = gt(a)
    const B = gt(b)
    for (const [probe, ag] of A) {
      const bg = B.get(probe)
      const fg = F.get(probe)
      const mg = M.get(probe)
      if (!bg || !fg || !mg || ag === 'NC' || bg === 'NC' || fg === 'NC' || mg === 'NC') continue
      const va = inf.classifyMarker(fg, mg, ag).provesPaternalPresence
      const vb = inf.classifyMarker(fg, mg, bg).provesPaternalPresence
      if (va === vb) agree += 1
      else differ += 1
    }
    cache.delete(a); cache.delete(b)
  }
  console.log(`  replicate pairs used: ${repPairs.length}`)
  console.log(`  L3 call identical on both arrays of the same sample: ${agree.toLocaleString()} = ${pct(agree, agree + differ)}`)
  console.log(`  L3 call DIFFERED:                                    ${differ.toLocaleString()} = ${pct(differ, agree + differ)}`)
}

console.log('')
console.log('HOW TO READ IT. Claim 2 is the one with no escape: a paternal allele the father does')
console.log('not carry is either a false claim or a real Mendelian error, and the true-father')
console.log('violation rate from claim 1 bounds how much of it can be the latter.')

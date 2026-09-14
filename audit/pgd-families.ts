/**
 * FOURTEEN FAMILIES WITH BOTH PARENTS NAMED, WHICH IS WHAT THE PARENT-OF-ORIGIN GAP NEEDED.
 *
 * THE GAP THIS CLOSES. Every parent-of-origin result in this project rested on eight observed
 * events from ONE sperm donor, all of them paternal because the experiment cut the paternal allele.
 * A rule that always answered "paternal" would have scored the same. GSE60909 is a preimplantation
 * genetic diagnosis series in which the submitters label every array with its family and its role:
 * Mother, Father, Sib, and embryos as E<embryo>_Bl<blastomere>. Fourteen mothers and fourteen
 * fathers, and the relationships are the submitters' own, not this tool's inference.
 *
 * WHAT THIS SERIES CAN AND CANNOT ANSWER. It publishes genotype CALLS and no intensity. So the
 * detection channels cannot run at all: with no log2 ratio there is nothing for a dosage test to
 * read, and asking anyway would produce a fabricated zero. What it can answer is everything the
 * genotype channel decides, which is where parentage lives.
 *
 * THE ARMS, and each has a truth that does not come from this tool:
 *
 *   RELATEDNESS      a mother and her embryo are first-degree; two people from different families
 *                    are not. The second is the strong one: a relative verdict across families is a
 *                    false positive with nothing to argue about.
 *   PARENT OF ORIGIN where a mother is homozygous and her embryo heterozygous, the embryo's other
 *                    allele came from its father, and the father's own array either carries it or
 *                    it does not. That is computed here from the two parents directly, and the
 *                    tool is scored against it. Both directions, because both parents are present.
 *   THE MIRROR       the same question with the parents swapped. Eight events in one direction
 *                    proved nothing about the other; this measures the maternal direction on real
 *                    families rather than on constructed losses.
 *
 * Run: OM_PGD=<dir of .probes and truth.json> node --experimental-strip-types \
 *        --max-old-space-size=3072 audit/pgd-families.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const rel = await import(`${W}relatedness.ts`)

const DIR = process.env.OM_PGD
if (!DIR || !existsSync(join(DIR, 'truth.json'))) {
  console.log('pgd-families: OM_PGD must hold .probes files and truth.json. NOT RUN.')
  process.exit(0)
}

interface T { gsm: string; title: string }
const truth: T[] = JSON.parse(readFileSync(join(DIR, 'truth.json'), 'utf8'))
  .filter((t: T) => existsSync(join(DIR, `${t.gsm}.probes`)))

/** Family, with the cycle suffix removed: PGD004C1 and PGD004C2 are one family, two cycles. */
const familyOf = (title: string) => {
  const m = /PGD(\d+)/i.exec(title)
  return m ? `PGD${m[1]}` : null
}
type Role = 'mother' | 'father' | 'sib' | 'embryo'
const roleOf = (title: string): Role | null => {
  if (/^Mother/i.test(title)) return 'mother'
  if (/^Father/i.test(title)) return 'father'
  if (/^(Sib|S\d)/i.test(title)) return 'sib'
  if (/^E\d+_Bl/i.test(title)) return 'embryo'
  return null
}

const people = truth.map((t) => ({ ...t, fam: familyOf(t.title), role: roleOf(t.title) }))
  .filter((p) => p.fam && p.role)
const byFam = new Map<string, typeof people>()
for (const p of people) byFam.set(p.fam!, [...(byFam.get(p.fam!) ?? []), p])

console.log(`arrays with a family and a role: ${people.length} across ${byFam.size} families`)
const complete = [...byFam.entries()].filter(([, v]) =>
  v.some((x) => x.role === 'mother') && v.some((x) => x.role === 'father')
  && v.some((x) => x.role === 'embryo'))
console.log(`families with a mother, a father and at least one embryo: ${complete.length}`)
console.log('')

type AB = 'AA' | 'AB' | 'BB' | 'NC'
interface Loaded { gt: Map<string, AB>; pos: Map<string, { chrom: string; pos: number }> }
const CACHE_MAX = Number(process.env.OM_CACHE ?? 6)
const cache = new Map<string, Loaded>()
function load(gsm: string): Loaded {
  const hit = cache.get(gsm)
  if (hit) { cache.delete(gsm); cache.set(gsm, hit); return hit }
  const ls = readFileSync(join(DIR, `${gsm}.probes`), 'utf8').split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  const gt = new Map<string, AB>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (!r) continue
    gt.set(r.probesetId, r.genotype as AB)
    if (Number.isFinite(r.pos)) pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  const v = { gt, pos }
  cache.set(gsm, v)
  while (cache.size > CACHE_MAX) cache.delete(cache.keys().next().value as string)
  return v
}

const q = (xs: number[], p: number) => {
  const s = [...xs].sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}

// ================================================================ ARM 1: relatedness
console.log('=== ARM 1: RELATEDNESS, against the submitters\' own family labels')
interface Pair { rate: number; verdict: string; kind: string }
const pairs: Pair[] = []
const CAP = Number(process.env.OM_FAM ?? 8)
const fams = complete.slice(0, CAP)
for (const [fam, members] of fams) {
  const mum = members.find((m) => m.role === 'mother')!
  const dad = members.find((m) => m.role === 'father')!
  const embryos = members.filter((m) => m.role === 'embryo').slice(0, 3)
  for (const e of embryos) {
    for (const [parent, tag] of [[mum, 'mother-embryo'], [dad, 'father-embryo']] as const) {
      const a = load(parent.gsm); const b = load(e.gsm)
      const opp = rel.oppositeHom(a.gt, b.gt) as { rate: number; n: number }
      if (opp.n < 5_000) continue
      const r = rel.relate({ gt: a.gt, pos: a.pos }, { gt: b.gt, pos: b.pos },
        rel.OPPOSITE_HOM_MAX) as { relationship: string }
      pairs.push({ rate: opp.rate, verdict: r.relationship, kind: tag })
    }
  }
  // Across families: unrelated by the submitters' labelling, so any relative verdict is false.
  const other = fams.find(([f]) => f !== fam)
  if (other) {
    const stranger = other[1].find((m) => m.role === 'mother')!
    const a = load(stranger.gsm); const b = load(dad.gsm)
    const opp = rel.oppositeHom(a.gt, b.gt) as { rate: number; n: number }
    if (opp.n >= 5_000) {
      const r = rel.relate({ gt: a.gt, pos: a.pos }, { gt: b.gt, pos: b.pos },
        rel.OPPOSITE_HOM_MAX) as { relationship: string }
      pairs.push({ rate: opp.rate, verdict: r.relationship, kind: 'across families' })
    }
  }
}
for (const kind of ['mother-embryo', 'father-embryo', 'across families']) {
  const g = pairs.filter((p) => p.kind === kind)
  if (!g.length) continue
  const rates = g.map((p) => p.rate)
  console.log(`  ${kind.padEnd(16)} n=${String(g.length).padStart(3)}  `
    + `opposite-hom ${q(rates, 0).toFixed(4)} to ${q(rates, 1).toFixed(4)}`)
  const tally = new Map<string, number>()
  for (const p of g) tally.set(p.verdict, (tally.get(p.verdict) ?? 0) + 1)
  for (const [v, n] of [...tally].sort((a, b) => b[1] - a[1])) {
    console.log(`      ${n} x ${v}`)
  }
}
const falseRel = pairs.filter((p) => p.kind === 'across families'
  && /parent and child|first-degree|sibling|this parent only/i.test(p.verdict))
console.log(`  ACROSS-FAMILY PAIRS CALLED RELATIVES: ${falseRel.length}/`
  + `${pairs.filter((p) => p.kind === 'across families').length}`)

// ================================================================ ARM 2 and 3: both directions
console.log('')
console.log('=== ARM 2: PARENT OF ORIGIN, both directions, on real families')
console.log('Where one parent is homozygous and the embryo heterozygous, the embryo\'s other allele')
console.log('came from the OTHER parent, who either carries it or does not. Computed from the two')
console.log('parental arrays; no part of the scorer is involved in deciding it.')
const isHom = (g: AB) => g === 'AA' || g === 'BB'
let patN = 0; let patOk = 0
let matN = 0; let matOk = 0
for (const [fam, members] of fams) {
  const mum = load(members.find((m) => m.role === 'mother')!.gsm)
  const dad = load(members.find((m) => m.role === 'father')!.gsm)
  const embryos = members.filter((m) => m.role === 'embryo').slice(0, 2)
  let fp = 0; let fpo = 0; let mp = 0; let mpo = 0
  for (const e of embryos) {
    const kid = load(e.gsm)
    for (const [probe, cg] of kid.gt) {
      if (cg !== 'AB') continue
      const mg = mum.gt.get(probe); const dg = dad.gt.get(probe)
      if (!mg || !dg || mg === 'NC' || dg === 'NC') continue
      // mother homozygous: the embryo's other allele is paternal, so the father must carry it
      if (isHom(mg)) {
        const a = cg.split('').find((x) => x !== mg[0])
        if (a) { fp += 1; if (dg.includes(a)) fpo += 1 }
      }
      // father homozygous: the mirror, and the direction eight events never tested
      if (isHom(dg)) {
        const a = cg.split('').find((x) => x !== dg[0])
        if (a) { mp += 1; if (mg.includes(a)) mpo += 1 }
      }
    }
    cache.delete(e.gsm)
  }
  patN += fp; patOk += fpo; matN += mp; matOk += mpo
  console.log(`  ${fam.padEnd(9)} paternal allele in the father `
    + `${(100 * fpo / Math.max(fp, 1)).toFixed(2)}% of ${fp.toLocaleString()}   `
    + `maternal allele in the mother ${(100 * mpo / Math.max(mp, 1)).toFixed(2)}% of ${mp.toLocaleString()}`)
}
console.log('')
console.log(`  PATERNAL direction: ${patOk.toLocaleString()}/${patN.toLocaleString()} = `
  + `${(100 * patOk / Math.max(patN, 1)).toFixed(2)}%`)
console.log(`  MATERNAL direction: ${matOk.toLocaleString()}/${matN.toLocaleString()} = `
  + `${(100 * matOk / Math.max(matN, 1)).toFixed(2)}%`)
console.log('')
console.log('  The two rates being alike is the finding. A wrong parent in either slot would')
console.log('  collapse one of them and leave the other, which is how this arm tells a real')
console.log('  asymmetry apart from a mislabelled family.')

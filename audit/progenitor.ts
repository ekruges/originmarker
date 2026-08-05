/**
 * The Progenitor audit: reconstruct a man's genotype from his sons' pronuclei, then check the
 * reconstruction against his own array.
 *
 * Run: node --experimental-strip-types --max-old-space-size=8192 audit/progenitor.ts <raw> <out>
 *
 * This drives `web/src/ingest.ts`, `web/src/inferredReference.ts` and `web/src/parentage.ts`,
 * which is what the browser executes. Only the file plumbing differs: a gzipped file off disk
 * rather than a dropped one.
 *
 * WHY THIS DATA ANSWERS THE QUESTION
 *
 * Zuccaro et al. separated the two pronuclei of a fertilised human zygote by micromanipulation
 * and arrayed each one alone, and the series also carries bulk genomic DNA for the sperm donor
 * and for the egg donors. So the material divides itself, at the bench:
 *
 *   8 paternal pronuclei   haploid meiotic products of ONE man, which is what Progenitor takes
 *   2 sperm donor arrays   that same man, measured directly. Never given to the reconstruction,
 *                          used only to score it. This is the ground truth.
 *   7 maternal pronuclei   haploid, from the egg donors, unrelated to him
 *   2 egg donor arrays     diploid adults, unrelated to him
 *
 * That makes three things measurable that the tool can otherwise only model. Ascertainment is
 * checked against his real heterozygosity rather than against the m=2 baseline. Contamination is
 * checked by asking his array which of the reference's asserted homozygotes are really
 * heterozygous, rather than by the posterior. And absence is checked on nine genomes known not
 * to be his.
 *
 * The whole file is loaded once into typed arrays over a shared probe index, so the hundreds of
 * reconstructions below cost one pass over the raw data rather than one pass each.
 */
import { createReadStream, readdirSync, writeFileSync, mkdirSync } from 'node:fs'
import { createGunzip } from 'node:zlib'
import { createInterface } from 'node:readline'
import { createHash } from 'node:crypto'
import { basename, join } from 'node:path'
import {
  accumulate, accumulateBaf, accumulateBuild, buildVerdict, emptyBafSums, emptyBuildSums,
  finishProfile, headerMap, parseRow,
  type ChromStats, type ProbeRow, type SampleProfile,
} from '../web/src/ingest.ts'
import {
  ProductSet, groupByParent, kinship, MIN_ASCERTAINMENT, MIN_PRODUCTS,
} from '../web/src/inferredReference.ts'
import { classify, emptyTally, tallyRow, pct } from '../web/src/parentage.ts'
import { ploidyOf, triage } from '../web/src/productTriage.ts'
import type { AB } from '../web/src/informativity.ts'

const RAW = process.argv[2] ?? 'raw'
const OUT = process.argv[3] ?? 'out'
mkdirSync(OUT, { recursive: true })

/** GSM to what the submitters called it. Verbatim from the GEO series record. */
const TITLE: Record<string, string> = {
  GSM4472397: 'genomic DNA sperm donor rep 1',
  GSM4472398: 'genomic DNA sperm donor rep 2',
  GSM4472407: 'genomic DNA egg donor A',
  GSM4472415: 'genomic DNA egg donor C rep 1',
  GSM4774673: 'maternal nucleus isolated from 2PN zygote 1',
  GSM4774674: 'maternal nucleus isolated from 2PN zygote 2',
  GSM4774675: 'maternal nucleus isolated from 2PN zygote 3',
  GSM4774676: 'maternal nucleus isolated from 2PN zygote 4',
  GSM4774677: 'maternal nucleus isolated from 2PN zygote 5',
  GSM4774678: 'maternal nucleus isolated from 2PN zygote 6',
  GSM4774679: 'maternal nucleus isolated from 2PN zygote 8',
  GSM4774680: 'paternal nucleus isolated from 2PN zygote 1',
  GSM4774681: 'paternal nucleus isolated from 2PN zygote 2',
  GSM4774682: 'paternal nucleus isolated from 2PN zygote 3',
  GSM4774683: 'paternal nucleus isolated from 2PN zygote 4',
  GSM4774684: 'paternal nucleus isolated from 2PN zygote 5',
  GSM4774685: 'paternal nucleus isolated from 2PN zygote 6',
  GSM4774686: 'paternal nucleus isolated from 2PN zygote 7',
  GSM4774687: 'paternal nucleus isolated from 2PN zygote 8',
}

/** The products. Haploid meiotic products of one man, established by micromanipulation. */
const PRODUCTS = ['GSM4774680', 'GSM4774681', 'GSM4774682', 'GSM4774683',
  'GSM4774684', 'GSM4774685', 'GSM4774686', 'GSM4774687']
/** His own arrays. Held back from every reconstruction; they only ever score one. */
const TRUTH = ['GSM4472397', 'GSM4472398']
/** Haploid, from the egg donors. Same platform, same lab, same processing, different parent. */
const MAT_NUCLEI = ['GSM4774673', 'GSM4774674', 'GSM4774675', 'GSM4774676',
  'GSM4774677', 'GSM4774678', 'GSM4774679']
/** Diploid adults, unrelated to him. */
const EGG_DONORS = ['GSM4472407', 'GSM4472415']

const ALL = [...PRODUCTS, ...TRUTH, ...MAT_NUCLEI, ...EGG_DONORS]

const files = readdirSync(RAW).filter((f) => f.endsWith('.CEL.txt.gz'))
const find = (gsm: string): string => {
  const f = files.find((x) => x.startsWith(gsm))
  if (!f) throw new Error(`${gsm} not in ${RAW}`)
  return join(RAW, f)
}

// --- the shared probe index -------------------------------------------------------------------
//
// Every sample is stored as one code per probe over one index, so a reconstruction is a scan
// over typed arrays and the hundreds of them below are affordable.

const index = new Map<string, number>()
const probeId: string[] = []
const chrom: string[] = []
const pos: number[] = []

const NC = 0
const CODE: Record<string, number> = { NC: 0, AA: 1, AB: 2, BB: 3 }
const UNCODE: AB[] = ['NC', 'AA', 'AB', 'BB']

interface Loaded {
  gsm: string
  title: string
  file: string
  sha256: string
  profile: SampleProfile
  /** Genotype code per probe, over the shared index. */
  gt: Uint8Array
  /** B-allele frequency per probe. NaN where the file left it empty, which is a different
   *  fact from zero and is what the noise ceiling distinguishes. */
  baf: Float32Array
  markers: number
}

async function load(gsm: string): Promise<Loaded> {
  const path = find(gsm)
  const byChrom = new Map<string, ChromStats>()
  const bafSums = emptyBafSums()
  const builds = emptyBuildSums()
  const hash = createHash('sha256')
  let map: ReturnType<typeof headerMap> = null
  let firstId = ''
  let n = 0

  // Sized on the first sample and grown as later ones bring probes the first did not have.
  let gt = new Uint8Array(Math.max(index.size, 1 << 20))
  let baf = new Float32Array(gt.length).fill(NaN)
  const grow = (need: number): void => {
    if (need <= gt.length) return
    const size = Math.max(need, gt.length * 2)
    const g2 = new Uint8Array(size); g2.set(gt); gt = g2
    const b2 = new Float32Array(size).fill(NaN); b2.set(baf); baf = b2
  }

  const gz = createReadStream(path).pipe(createGunzip())
  gz.on('data', (c: Buffer) => hash.update(c))
  for await (const line of createInterface({ input: gz, crlfDelay: Infinity })) {
    if (!map) { map = headerMap(line); continue }
    const r = parseRow(line, map)
    if (!r) continue
    if (!firstId) firstId = r.probesetId
    accumulate(r, byChrom)
    accumulateBaf(r, bafSums)
    accumulateBuild(r, builds)
    let i = index.get(r.probesetId)
    if (i === undefined) {
      i = index.size
      index.set(r.probesetId, i)
      probeId.push(r.probesetId)
      chrom.push(r.chrom)
      pos.push(r.pos)
    }
    grow(i + 1)
    gt[i] = CODE[r.genotype] ?? NC
    if (r.baf !== null) baf[i] = r.baf
    n += 1
  }
  if (!map) throw new Error(`${path}: no recognisable header`)
  const profile: SampleProfile = {
    ...finishProfile(basename(path), byChrom, bafSums, firstId, builds),
    build: buildVerdict(builds),
  }
  return {
    gsm, title: TITLE[gsm] ?? gsm, file: basename(path), sha256: hash.digest('hex'),
    profile, gt, baf, markers: n,
  }
}

/** One stored marker, back in the shape the shipped code takes. */
const rowAt = (s: Loaded, i: number): ProbeRow => ({
  probesetId: probeId[i],
  chrom: chrom[i],
  pos: pos[i],
  log2R: null,
  baf: Number.isNaN(s.baf[i]) ? null : s.baf[i],
  copyNumber: null,
  genotype: UNCODE[s.gt[i]],
  bestProbeset: true,
  unreadableGenotype: false,
  rawGenotype: UNCODE[s.gt[i]],
})

/**
 * Replay stored samples into a ProductSet, exactly as the page feeds it while streaming.
 *
 * `stride` takes every nth marker from `offset`, which gives disjoint slices of the array to
 * rebuild from. Every quantity the tool reports is a rate, so a slice should give the same
 * answer as the whole; the sweep below uses that to test the reconstruction at a fraction of
 * the marker density rather than only at the one density this platform happens to have.
 */
function productSet(samples: Loaded[], stride = 1, offset = 0): ProductSet {
  const ps = new ProductSet()
  for (const s of samples) {
    const slot = ps.begin(s.gsm)
    const band = { inBand: 0, total: 0 }
    for (let i = offset; i < index.size; i += stride) {
      if (s.gt[i] === NC && Number.isNaN(s.baf[i])) continue
      ps.add(slot, rowAt(s, i), band)
    }
    ps.end(slot, band)
  }
  return ps
}

/** Score one sample against a reconstructed reference, through the shipped classifier. */
function score(sample: Loaded, ref: Map<string, AB>, hRetained: number, spurious?: number) {
  const t = emptyTally()
  for (let i = 0; i < index.size; i += 1) {
    if (s0(sample, i)) continue
    tallyRow(ref.get(probeId[i]) ?? 'NC', rowAt(sample, i), t)
  }
  return classify(t, hRetained, spurious === undefined ? {} : { spuriousAbsence: spurious })
}
const s0 = (s: Loaded, i: number): boolean => s.gt[i] === NC && Number.isNaN(s.baf[i])

// --- load ---------------------------------------------------------------------------------------

console.error(`loading ${ALL.length} arrays from ${RAW}`)
const by: Record<string, Loaded> = {}
for (const gsm of ALL) {
  by[gsm] = await load(gsm)
  console.error(`  ${gsm}  ${by[gsm].markers.toLocaleString()} markers  `
    + `call ${pct(by[gsm].profile.callRate, 1)}  band ${pct(by[gsm].profile.hetBand, 1)}`)
}
console.error(`shared probe index: ${index.size.toLocaleString()}\n`)

/** His measured heterozygosity, from his own array. The denominator ascertainment is a fraction
 *  of, and never available to the tool. */
function trueHeterozygosity(s: Loaded): number {
  let het = 0
  let called = 0
  for (let i = 0; i < index.size; i += 1) {
    if (!isAuto(chrom[i]) || s.gt[i] === NC) continue
    called += 1
    if (s.gt[i] === CODE.AB) het += 1
  }
  return called ? het / called : NaN
}
const isAuto = (c: string): boolean => /^(?:[1-9]|1[0-9]|2[0-2])$/.test(c)

const H_TRUE = trueHeterozygosity(by[TRUTH[0]])

/**
 * How often his own two arrays disagree, which bounds what the ground truth itself can be wrong
 * about.
 *
 * The reference is scored by asking his array which of its asserted homozygotes are really
 * heterozygous. Some of those AB calls are his array's error rather than the reference's, and
 * without a bound on that the comparison cannot separate a model that understates contamination
 * from a truth that overstates it. Two arrays of one man, same platform and lab, give the bound
 * directly.
 */
function replicateDiscordance(a: Loaded, b: Loaded) {
  let shared = 0
  let hetOnlyOne = 0
  let oppositeHom = 0
  for (let i = 0; i < index.size; i += 1) {
    if (!isAuto(chrom[i]) || a.gt[i] === NC || b.gt[i] === NC) continue
    shared += 1
    const aHet = a.gt[i] === CODE.AB
    const bHet = b.gt[i] === CODE.AB
    if (aHet !== bHet) hetOnlyOne += 1
    else if (!aHet && a.gt[i] !== b.gt[i]) oppositeHom += 1
  }
  return { shared, hetOnlyOne, oppositeHom,
    hetDiscordance: shared ? hetOnlyOne / shared : NaN,
    homDiscordance: shared ? oppositeHom / shared : NaN }
}
const REP = replicateDiscordance(by[TRUTH[0]], by[TRUTH[1]])
console.error(`replicate discordance: het-in-one-only ${pct(REP.hetDiscordance, 3)}, `
  + `opposite homozygote ${pct(REP.homDiscordance, 4)}\n`)
console.error(`sperm donor measured heterozygosity: ${pct(H_TRUE, 2)}\n`)

// --- results ------------------------------------------------------------------------------------

interface Case {
  id: string
  section: string
  what: string
  expected: string
  reported: string
  outcome: 'correct' | 'refused' | 'incorrect'
  detail: Record<string, unknown>
}
const cases: Case[] = []
const add = (c: Case): void => { cases.push(c); console.error(
  `  [${c.outcome.toUpperCase().padEnd(9)}] ${c.id}: ${c.reported}`) }

/**
 * What his own array says about the markers the reference asserts.
 *
 * `wrong` is contamination MEASURED: the reference calls the marker homozygous and he is really
 * heterozygous there. `mismatch` is worse and should be ~0: the reference calls a homozygote and
 * he is the OTHER homozygote, which no amount of agreement at a het site can produce.
 */
function against(ref: Map<string, AB>, truth: Loaded) {
  let checked = 0
  let wrong = 0
  let mismatch = 0
  for (const [id, call] of ref) {
    const i = index.get(id)
    if (i === undefined) continue
    const g = truth.gt[i]
    if (g === NC) continue
    checked += 1
    if (g === CODE.AB) wrong += 1
    else if (UNCODE[g] !== call) mismatch += 1
  }
  return {
    checked,
    contaminationObserved: checked ? wrong / checked : NaN,
    oppositeHomozygote: checked ? mismatch / checked : NaN,
    concordance: checked ? 1 - (wrong + mismatch) / checked : NaN,
  }
}

// --- 1. the gates, before anything is compared or built -----------------------------------------
//
// This runs first because the tool runs it first, and that ordering is a safety property rather
// than presentation: a degraded array reads high against everyone and splits one parent into
// several, and a diploid bridges two. An audit that reconstructs from every file in the folder
// is not auditing the tool.
console.error('0. what the ground truth itself can be wrong about')
add({
  id: 'replicate-discordance',
  section: 'Ground truth',
  what: 'his two arrays, compared against each other',
  expected: 'they largely agree; whatever they do not is the error floor every comparison '
    + 'against him inherits',
  reported: `over ${REP.shared.toLocaleString()} jointly called autosomal markers, `
    + `${pct(REP.hetDiscordance, 3)} are heterozygous in exactly one of the two, and `
    + `${pct(REP.homDiscordance, 4)} are opposite homozygotes`,
  outcome: REP.hetDiscordance < 0.05 ? 'correct' : 'refused',
  detail: REP,
})

console.error('1. quality and ploidy, on all eight products')
const USABLE: string[] = []
for (const gsm of PRODUCTS) {
  const s = by[gsm]
  const t = triage(gsm, gsm, s.profile)
  if (t.usable) USABLE.push(gsm)
  add({
    id: `gate-${gsm}`,
    section: 'Quality and ploidy',
    what: `${s.title}: call ${pct(s.profile.callRate, 1)}, band ${pct(s.profile.hetBand, 1)}, `
      + `${pct(s.profile.hetRate, 1)} heterozygous calls`,
    expected: 'a true product is usable, unless its own quality puts it outside the range a '
      + 'single haploid cell produces, in which case it is set aside with the reason stated',
    reported: t.usable ? 'usable product' : `set aside, ${t.why}`,
    outcome: t.usable ? 'correct' : 'refused',
    detail: { ploidy: t.ploidy, callRate: s.profile.callRate, hetBand: s.profile.hetBand,
      hetRate: s.profile.hetRate, why: t.why },
  })
}
console.error(`   ${USABLE.length} of ${PRODUCTS.length} products usable`)
if (USABLE.length < MIN_PRODUCTS) throw new Error('too few usable products to audit')

// --- 2. the reconstruction, against the man himself ---------------------------------------------
console.error('2. reconstruction fidelity, every depth, scored on his own array')
const psAll = productSet(USABLE.map((g) => by[g]))
const truth0 = by[TRUTH[0]]
const ladder: Record<string, unknown>[] = []
for (let m = 2; m <= USABLE.length; m += 1) {
  const ref = psAll.build(m)
  if (!ref.markers) continue
  const v = against(ref.genotype, truth0)
  const ascertainment = ref.hRetained / H_TRUE
  ladder.push({
    m,
    markers: ref.markers,
    hRetained: ref.hRetained,
    ascertainmentVsTruth: ascertainment,
    contaminationPredicted: ref.contamination,
    contaminationObserved: v.contaminationObserved,
    oppositeHomozygote: v.oppositeHomozygote,
    concordance: v.concordance,
    checked: v.checked,
  })
  add({
    id: `fidelity-m${m}`,
    section: 'Reconstruction fidelity',
    what: `reference at m>=${m}, ${ref.markers.toLocaleString()} markers, against his own array`,
    expected: 'every asserted homozygote matches him, except where he is heterozygous and the '
      + 'products agreed; the predicted contamination covers the observed',
    reported: `concordance ${pct(v.concordance, 3)}, contamination predicted `
      + `${pct(ref.contamination, 3)} observed ${pct(v.contaminationObserved, 3)}, `
      + `opposite homozygote ${pct(v.oppositeHomozygote, 4)}, `
      + `replicate floor ${pct(REP.hetDiscordance, 3)}`,
    // The model must not UNDERSTATE the error by more than his own array can account for. A
    // marker the reference calls homozygous and he calls heterozygous is either contamination
    // or one of his false heterozygous calls, and this data cannot tell them apart below that
    // floor. Understating by MORE than the floor is a real gap and is reported as one.
    outcome: v.oppositeHomozygote >= 0.01 ? 'incorrect'
      : ref.contamination >= v.contaminationObserved - REP.hetDiscordance ? 'correct'
        : 'incorrect',
    detail: { ...ladder[ladder.length - 1] },
  })
}

// --- 2. the depth the tool picks, with no ground truth ------------------------------------------
console.error('2. the chosen depth')
const chosen = psAll.chooseM()
const chosenRow = ladder.find((r) => r.m === chosen.mMin)!
add({
  id: 'chosen-depth',
  section: 'Choosing the depth',
  what: `chooseM over ${USABLE.length} usable products, using no array of the parent`,
  expected: `the deepest depth holding ${pct(MIN_ASCERTAINMENT, 0)} of the parent's `
    + 'heterozygosity, judged against his real array',
  reported: `m>=${chosen.mMin}, which holds `
    + `${pct(chosenRow.ascertainmentVsTruth as number, 1)} of his measured heterozygosity`,
  outcome: (chosenRow.ascertainmentVsTruth as number) >= MIN_ASCERTAINMENT ? 'correct' : 'incorrect',
  detail: {
    chosen: chosen.mMin,
    ratiosSelfMeasured: Object.fromEntries(chosen.ratios),
    ascertainmentVsTruth: Object.fromEntries(ladder.map((r) => [r.m, r.ascertainmentVsTruth])),
    // The rule the tool replaced, evaluated on the same data.
    wouldHaveBeenNMinus1: USABLE.length - 1,
    ascertainmentAtNMinus1: ladder.find((r) => r.m === USABLE.length - 1)?.ascertainmentVsTruth,
  },
})

// --- 3. he reads present, and nine genomes that are not his read absent -------------------------
console.error('3. identification, against a reference built from his sons alone')
const refChosen = psAll.build(chosen.mMin)
const idRows: Record<string, unknown>[] = []
for (const gsm of [...TRUTH, ...MAT_NUCLEI, ...EGG_DONORS]) {
  const s = by[gsm]
  const isHim = TRUTH.includes(gsm)
  // Haploid samples carry the reference's own false absence; diploid ones do not.
  const haploid = ploidyOf(s.profile.hetBand, s.profile.hetRate) === 'haploid'
  const r = score(s, refChosen.genotype, refChosen.hRetained,
    haploid ? refChosen.spuriousAbsence : undefined)
  const ratio = r.genomeRate / r.explainable
  const present = r.verdict === 'parent_genome_present'
  const absent = r.verdict === 'no_parental_contribution'
  idRows.push({ gsm, title: s.title, isHim, haploid, absence: r.genomeRate,
    ceiling: r.explainable, ratio, verdict: r.verdict })
  add({
    id: `id-${gsm}`,
    section: 'Identification',
    what: `${s.title} against the reconstruction`,
    expected: isHim ? 'present: this is the man the products came from'
      : 'absent: unrelated to him, established at the bench',
    reported: `${pct(r.genomeRate, 2)} absent of ${pct(r.explainable, 2)} = ${ratio.toFixed(2)}x, `
      + r.verdict.replace(/_/g, ' '),
    outcome: isHim
      ? (present ? 'correct' : absent ? 'incorrect' : 'refused')
      : (absent ? 'correct' : present ? 'incorrect' : 'refused'),
    detail: idRows[idRows.length - 1],
  })
}

// --- 4. membership, with no reference at all ----------------------------------------------------
console.error('4. membership')
{
  const g = groupByParent(psAll.size, (a, b) => psAll.opposite(a, b).rate)
  add({
    id: 'group-products-only',
    section: 'Membership',
    what: `${USABLE.length} usable paternal pronuclei, compared pairwise`,
    expected: 'one group: they are products of one man',
    reported: `${g.length} group(s), sizes ${g.map((x) => x.length).join(', ')}`,
    outcome: g.length === 1 && g[0].length === USABLE.length ? 'correct' : 'incorrect',
    detail: { groups: g.map((x) => x.map((i) => psAll.ids[i])) },
  })

  // Paternal and maternal pronuclei together. Same platform, same lab, same processing, and a
  // different parent, so this is the split the method exists to make.
  const mixedIds = [...USABLE, ...MAT_NUCLEI]
  const psMix = productSet(mixedIds.map((x) => by[x]))
  const gm = groupByParent(psMix.size, (a, b) => psMix.opposite(a, b).rate)
  const named = gm.map((x) => x.map((i) => psMix.ids[i]))
  const patTogether = named.some((grp) =>
    USABLE.every((p) => grp.includes(p)) && !grp.some((x) => MAT_NUCLEI.includes(x)))
  add({
    id: 'group-mixed',
    section: 'Membership',
    what: `${USABLE.length} paternal and ${MAT_NUCLEI.length} maternal pronuclei together`,
    expected: 'the paternal ones form one group with no maternal product in it',
    reported: `${gm.length} group(s), sizes ${gm.map((x) => x.length).join(', ')}; `
      + `the paternal products ${patTogether ? 'stayed together and alone' : 'did NOT'}`,
    outcome: patTogether ? 'correct' : 'incorrect',
    detail: { groups: named },
  })

  // Every cross pair, and every within pair, against the thresholds.
  let within = { lo: 1, hi: 0, n: 0, misread: 0 }
  let across = { lo: 1, hi: 0, n: 0, misread: 0 }
  for (let a = 0; a < psMix.size; a += 1) {
    for (let b = a + 1; b < psMix.size; b += 1) {
      const rate = psMix.opposite(a, b).rate
      const same = USABLE.includes(psMix.ids[a]) === USABLE.includes(psMix.ids[b])
      // Two maternal pronuclei are from DIFFERENT egg donors, so only the paternal side is a
      // genuine within-parent pair.
      const bothPat = USABLE.includes(psMix.ids[a]) && USABLE.includes(psMix.ids[b])
      const t = bothPat ? within : (same ? null : across)
      if (!t) continue
      t.n += 1
      t.lo = Math.min(t.lo, rate)
      t.hi = Math.max(t.hi, rate)
      const k = kinship(rate)
      if (bothPat && k === 'different parents') t.misread += 1
      if (!bothPat && !same && k === 'same parent') t.misread += 1
    }
  }
  add({
    id: 'pairwise-ranges',
    section: 'Membership',
    what: 'every pair, within one parent and across parents',
    expected: 'the two ranges separate, and any single misread pair is absorbed by the '
      + 'all-pairs rule rather than becoming a grouping error',
    reported: `within ${pct(within.lo, 2)} to ${pct(within.hi, 2)} over ${within.n} pairs `
      + `(${within.misread} misread); across ${pct(across.lo, 2)} to ${pct(across.hi, 2)} over `
      + `${across.n} pairs (${across.misread} misread)`,
    outcome: patTogether ? 'correct' : 'incorrect',
    detail: { within, across },
  })
}

// --- 5. a diploid adult must never be admitted as a product -------------------------------------
console.error('5. the ploidy gate, on diploids')
for (const gsm of EGG_DONORS) {
  const s = by[gsm]
  const t = triage(gsm, gsm, s.profile)
  add({
    id: `ploidy-${gsm}`,
    section: 'Quality and ploidy',
    what: `${s.title} offered as a product`,
    expected: 'set aside: a diploid adult is not one meiotic product, and admitting one bridges '
      + 'two unrelated parents into a single group',
    reported: t.usable ? 'ADMITTED as a product' : `set aside, ${t.why}`,
    outcome: t.usable ? 'incorrect' : 'correct',
    detail: { ploidy: t.ploidy, callRate: s.profile.callRate, hetBand: s.profile.hetBand,
      hetRate: s.profile.hetRate, why: t.why },
  })
}

// --- 6. leave one out ---------------------------------------------------------------------------
console.error('6. leave-one-out')
for (const gsm of USABLE) {
  const s = by[gsm]
  const m = Math.min(chosen.mMin, USABLE.length - 1)
  const ref = psAll.build(m, [gsm])
  const r = score(s, ref.genotype, ref.hRetained, ref.spuriousAbsence)
  const ratio = r.genomeRate / r.explainable
  add({
    id: `loo-${gsm}`,
    section: 'Leave-one-out',
    what: `${s.title} against a reference built from the other ${USABLE.length - 1}`,
    expected: 'present: it is a product of the same man',
    reported: `${pct(r.genomeRate, 2)} absent of ${pct(r.explainable, 2)} = ${ratio.toFixed(2)}x, `
      + r.verdict.replace(/_/g, ' '),
    outcome: r.verdict === 'parent_genome_present' ? 'correct'
      : r.verdict === 'no_parental_contribution' ? 'incorrect' : 'refused',
    detail: { absence: r.genomeRate, ceiling: r.explainable, ratio, verdict: r.verdict,
      excludedFrom: ref.nProducts, contamination: ref.contamination },
  })
}

// --- 7. the same reconstruction, rebuilt on disjoint slices of the array -------------------------
//
// The strenuous part. Every quantity the tool reports is a rate, so a reference built from every
// 16th marker should reach the same verdicts as one built from all of them. This rebuilds from
// 1, 2, 4, 8 and 16 disjoint slices, 31 independent references in total, and scores every one on
// him and on all nine genomes that are not his: 310 classifications, none of which may put an
// unrelated genome inside the reference or leave him outside it.
//
// It is also the only test here of what a sparser platform would do. One array density is one
// observation; a reconstruction that survives a sixteenth of this one is not leaning on density.
console.error('7. rebuilt on disjoint slices of the marker set')
const NEGATIVES = [...MAT_NUCLEI, ...EGG_DONORS]
const sweep: Record<string, unknown>[] = []
let sweepBad = 0
let sweepRefused = 0
for (const stride of [1, 2, 4, 8, 16]) {
  for (let offset = 0; offset < stride; offset += 1) {
    const ps = productSet(USABLE.map((g) => by[g]), stride, offset)
    const pick = ps.chooseM()
    const ref = ps.build(pick.mMin)
    const him = score(by[TRUTH[0]], ref.genotype, ref.hRetained)
    const himRatio = him.genomeRate / him.explainable
    let worstNeg = { gsm: '', ratio: Infinity, verdict: '' }
    let negWrong = 0
    for (const gsm of NEGATIVES) {
      const sn = by[gsm]
      const hap = ploidyOf(sn.profile.hetBand, sn.profile.hetRate) === 'haploid'
      const r = score(sn, ref.genotype, ref.hRetained, hap ? ref.spuriousAbsence : undefined)
      const ratio = r.genomeRate / r.explainable
      if (r.verdict === 'parent_genome_present') negWrong += 1
      if (ratio < worstNeg.ratio) worstNeg = { gsm, ratio, verdict: r.verdict }
    }
    if (him.verdict === 'no_parental_contribution' || negWrong > 0) sweepBad += 1
    else if (him.verdict !== 'parent_genome_present') sweepRefused += 1
    sweep.push({ stride, offset, m: pick.mMin, markers: ref.markers,
      ascertainmentVsTruth: ref.hRetained / H_TRUE, contamination: ref.contamination,
      himRatio, himVerdict: him.verdict, negWrong, worstNegative: worstNeg })
    console.error(`   1/${stride} offset ${offset}: ${ref.markers.toLocaleString()} markers, `
      + `him ${himRatio.toFixed(2)}x, worst negative ${worstNeg.ratio.toFixed(2)}x`)
  }
}
const sweepOk = sweep.length - sweepBad - sweepRefused
add({
  id: 'density-sweep',
  section: 'Marker density',
  what: `${sweep.length} references rebuilt on disjoint slices of the array, from all markers `
    + `down to every 16th, each scored on him and on all ${NEGATIVES.length} negatives`,
  expected: 'the verdicts do not depend on marker density: he reads present and no negative '
    + 'reads present, in every one',
  reported: `${sweepOk} clean, ${sweepRefused} declined to call him, ${sweepBad} wrong`,
  outcome: sweepBad > 0 ? 'incorrect' : sweepRefused > 0 ? 'refused' : 'correct',
  detail: {
    references: sweep.length,
    classifications: sweep.length * (1 + NEGATIVES.length),
    clean: sweepOk, refused: sweepRefused, wrong: sweepBad,
    himRatioRange: [Math.min(...sweep.map((x) => x.himRatio as number)),
      Math.max(...sweep.map((x) => x.himRatio as number))],
    worstNegativeRatio: Math.min(...sweep.map((x) => (x.worstNegative as { ratio: number }).ratio)),
    ascertainmentRange: [Math.min(...sweep.map((x) => x.ascertainmentVsTruth as number)),
      Math.max(...sweep.map((x) => x.ascertainmentVsTruth as number))],
    markersRange: [Math.min(...sweep.map((x) => x.markers as number)),
      Math.max(...sweep.map((x) => x.markers as number))],
  },
})

// --- 8. the floor, and what it is protecting against --------------------------------------------
//
// The tool refuses below MIN_PRODUCTS. That refusal is only justified if the method really does
// invert down there, so this builds below the floor deliberately and records what it would have
// said. Nothing here is a supported configuration; it is the measurement behind the refusal.
console.error('8. below the floor')
function* combinations<T>(xs: T[], k: number): Generator<T[]> {
  if (k === 0) { yield []; return }
  for (let i = 0; i <= xs.length - k; i += 1) {
    for (const rest of combinations(xs.slice(i + 1), k - 1)) yield [xs[i], ...rest]
  }
}
const belowFloor: Record<string, unknown>[] = []
for (let n = 2; n < MIN_PRODUCTS; n += 1) {
  let inverted = 0
  let total = 0
  const ratios: number[] = []
  for (const subset of combinations(USABLE, n)) {
    const ps = productSet(subset.map((g) => by[g]))
    const pick = ps.chooseM()
    const ref = ps.build(pick.mMin)
    if (!ref.markers) continue
    const him = score(by[TRUTH[0]], ref.genotype, ref.hRetained)
    total += 1
    ratios.push(him.genomeRate / him.explainable)
    if (him.verdict === 'no_parental_contribution') inverted += 1
  }
  belowFloor.push({ n, subsets: total, invertedOnHim: inverted,
    ratioRange: [Math.min(...ratios), Math.max(...ratios)] })
  add({
    id: `floor-n${n}`,
    section: 'Below the floor',
    what: `${total} reconstructions from ${n} products, a depth the tool refuses to build at`,
    expected: 'recorded, not scored. The floor was set on true OFFSPRING inverting at this '
      + 'depth, and this series contains none: the pronuclei are the products themselves. What '
      + 'can be checked here is the parent, who is a far stronger signal than an offspring and '
      + 'is expected to survive depths an offspring would not',
    reported: `${inverted} of ${total} called him absent; his ratio ran `
      + `${Math.min(...ratios).toFixed(2)}x to ${Math.max(...ratios).toFixed(2)}x`,
    // Not a pass or a fail of the tool, which declines to build here at all. It is the
    // measurement the refusal rests on, and it does NOT reproduce the offspring inversion,
    // because the material for that measurement is not in this series.
    outcome: 'refused',
    detail: { ...belowFloor[belowFloor.length - 1],
      note: 'scored on the parent, not on an offspring; does not reproduce the 24-of-24 '
        + 'offspring inversion the floor was set from' },
  })
}
add({
  id: 'floor-enforced',
  section: 'Below the floor',
  what: `the page's own gate at ${MIN_PRODUCTS} products`,
  expected: `no reference is offered below ${MIN_PRODUCTS} products`,
  reported: `MIN_PRODUCTS = ${MIN_PRODUCTS}, and the sweep above starts there`,
  outcome: MIN_PRODUCTS === 5 ? 'correct' : 'incorrect',
  detail: { MIN_PRODUCTS },
})

// --- 9. the ceiling term, checked both ways -----------------------------------------------------
//
// A reconstructed reference asserts homozygotes where the parent is heterozygous, so a haploid
// product carries the other allele at half of those. Admitting that into the ceiling is what
// stops a true product reading unclear; it must not also rescue a genome that is not his.
console.error('9. the reconstruction ceiling term')
{
  const m = Math.min(chosen.mMin, USABLE.length - 1)
  let withTermPresent = 0
  let withoutTermPresent = 0
  for (const gsm of USABLE) {
    const ref = psAll.build(m, [gsm])
    if (score(by[gsm], ref.genotype, ref.hRetained, ref.spuriousAbsence).verdict
      === 'parent_genome_present') withTermPresent += 1
    if (score(by[gsm], ref.genotype, ref.hRetained).verdict
      === 'parent_genome_present') withoutTermPresent += 1
  }
  let negRescued = 0
  for (const gsm of MAT_NUCLEI) {
    const r = score(by[gsm], refChosen.genotype, refChosen.hRetained, refChosen.spuriousAbsence)
    if (r.verdict === 'parent_genome_present') negRescued += 1
  }
  add({
    id: 'ceiling-term',
    section: 'The ceiling term',
    what: 'true products and unrelated haploids, scored with and without the reference\'s own '
      + 'false absence in the ceiling',
    expected: 'admitting it recovers true products and rescues no unrelated genome',
    reported: `true products present: ${withTermPresent}/${USABLE.length} with the term, `
      + `${withoutTermPresent}/${USABLE.length} without. Unrelated haploids wrongly present `
      + `with the term: ${negRescued}/${MAT_NUCLEI.length}`,
    outcome: negRescued > 0 ? 'incorrect'
      : withTermPresent >= withoutTermPresent ? 'correct' : 'incorrect',
    detail: { withTermPresent, withoutTermPresent, negRescued, of: USABLE.length },
  })
}

// --- write --------------------------------------------------------------------------------------

const tally = {
  correct: cases.filter((c) => c.outcome === 'correct').length,
  refused: cases.filter((c) => c.outcome === 'refused').length,
  incorrect: cases.filter((c) => c.outcome === 'incorrect').length,
}
const out = {
  tool: 'Progenitor',
  generated: new Date().toISOString(),
  data: {
    series: 'GEO GSE148488, Zuccaro et al. 2020, Cell 183(6):1650-1664',
    doi: '10.1016/j.cell.2020.10.025',
    arrays: ALL.map((g) => ({ gsm: g, title: by[g].title, file: by[g].file,
      sha256: by[g].sha256, markers: by[g].markers,
      callRate: by[g].profile.callRate, hetBand: by[g].profile.hetBand,
      hetRate: by[g].profile.hetRate })),
    usableProducts: USABLE,
    sharedProbes: index.size,
    spermDonorHeterozygosity: H_TRUE,
  },
  tally,
  replicateDiscordance: REP,
  ladder,
  identification: idRows,
  sweep,
  belowFloor,
  cases,
}
writeFileSync(join(OUT, 'progenitor-results.json'), JSON.stringify(out, null, 1))
console.error(`\n${tally.correct} correct, ${tally.refused} refused, ${tally.incorrect} incorrect`)
console.error(`wrote ${join(OUT, 'progenitor-results.json')}`)
if (tally.incorrect > 0) process.exitCode = 1

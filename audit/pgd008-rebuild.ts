/**
 * Step 4 on material that carries enough independent embryos.
 *
 * The reconstruction's first bound is that a parent transmits to an EMBRYO, so the information is
 * the count of independent embryos and not of arrays. GSE148488 has one usable family with six.
 * GSE60909 (Zamani Esteki 2015, PMID 25983246) arrays both parents of a couple alongside
 * blastomeres of fourteen embryos over two IVF cycles, mostly one blastomere each, which is
 * fourteen near-independent transmissions.
 *
 * The missing parent's own array is held out and used only to score.
 *
 * Criteria are in audit/opposite-hom-rebuild-preregistration.txt, committed before any of this ran.
 *
 * Run: OM_DIR=<dir of split .probes.gz> OM_FAMILY=PGD008 \
 *        node --experimental-strip-types audit/pgd008-rebuild.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const op = await import(`${W}oneParentOrigin.ts`)

const DIR = process.env.OM_DIR
const FAM = process.env.OM_FAMILY ?? 'PGD008'
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }

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
const gtOf = (file: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(join(DIR, file))()) {
    if (!/^\d+$/.test(r.chrom)) continue
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const files = readdirSync(DIR).filter((f) => f.endsWith('.probes.gz'))
const find = (pat: RegExp) => files.find((f) => pat.test(f))
const fatherFile = find(new RegExp(`Father_${FAM}\\.`))
const motherFile = find(new RegExp(`Mother_${FAM}\\.`))
if (!fatherFile || !motherFile) { console.log(`family ${FAM} incomplete. NOT RUN.`); process.exit(0) }

// Blastomeres of this couple, grouped by CYCLE AND EMBRYO: the numbering restarts each cycle, so
// C1 E06 and C2 E06 are different embryos.
const cells = files.filter((f) => new RegExp(`_E\\d+_\\S*_${FAM}C\\d`).test(f)
  || new RegExp(`E\\d+_Bl\\d+_${FAM}C\\d`).test(f))
const embryoOf = new Map<string, string>()
for (const f of cells) {
  const m = new RegExp(`(E\\d+)_Bl\\d+_(${FAM}C\\d)`).exec(f)
  if (m) embryoOf.set(f, `${m[2]}:${m[1]}`)
}
const embryos = [...new Set(embryoOf.values())].sort()
console.log(`family ${FAM}: ${cells.length} blastomere arrays over ${embryos.length} embryos`)
console.log(`  ${embryos.map((e) => `${e.replace(`${FAM}`, '')} x${[...embryoOf.values()].filter((x) => x === e).length}`).join(', ')}\n`)
if (embryos.length < 2) { console.log('not enough embryos. NOT RUN.'); process.exit(0) }

const faGt = gtOf(fatherFile)
const moGt = gtOf(motherFile)
console.log(`known parent ${faGt.size.toLocaleString()} called, held-out parent ${moGt.size.toLocaleString()} called`)

// Each blastomere's calls and its own dropout, inferred from its heterozygosity as the pipeline does.
const cellGt = new Map<string, { gt: Map<string, string>; ado: number; callRate: number }>()
for (const f of cells) {
  const gt = gtOf(f)
  let het = 0
  for (const v of gt.values()) if (v === 'AB') het += 1
  const hetRate = gt.size ? het / gt.size : NaN
  cellGt.set(f, { gt, ado: op.inferDropout(hetRate) as number, callRate: gt.size / 298_563 })
}
const ados = [...cellGt.values()].map((x) => x.ado).filter(Number.isFinite)
console.log(`  mean blastomere dropout ${(ados.reduce((a, b) => a + b, 0) / ados.length).toFixed(3)}, `
  + `call rates ${Math.min(...[...cellGt.values()].map((x) => x.callRate)).toFixed(2)} to `
  + `${Math.max(...[...cellGt.values()].map((x) => x.callRate)).toFixed(2)}`)

// The per-marker prior, from the OTHER couples' parent arrays. Nothing of this family is read.
const panel = files.filter((f) => /(Mother|Father)_PGD/.test(f) && !f.includes(FAM))
const bCount = new Map<string, number>()
const nCount = new Map<string, number>()
for (const f of panel) {
  for (const [probe, v] of gtOf(f)) {
    if (!faGt.has(probe)) continue
    bCount.set(probe, (bCount.get(probe) ?? 0) + (v === 'BB' ? 2 : v === 'AB' ? 1 : 0))
    nCount.set(probe, (nCount.get(probe) ?? 0) + 2)
  }
}
console.log(`  frequency panel: ${panel.length} unrelated donor arrays, ${2 * panel.length} chromosomes`)
/** Frequency of the allele the known parent LACKS, which the wanted genotype is homozygous for. */
const lackFreq = (probe: string) => {
  const fg = faGt.get(probe)
  if (!fg || fg === 'AB') return NaN
  const b = bCount.get(probe); const n = nCount.get(probe)
  if (b === undefined || !n) return NaN
  return fg === 'AA' ? b / n : (n - b) / n
}

/**
 * Drop-in measured on THIS material, from the other families of the series.
 *
 * 0.0435 was measured on bulk Axiom arrays. These are MDA-amplified single blastomeres on an
 * Illumina panel at a mean dropout of 0.474, and a model carrying the wrong drop-in is overconfident
 * in exactly the way the posterior was. Where both parents of ANOTHER family are homozygous for the
 * same allele their blastomere cannot be heterozygous, so every heterozygous call there is drop-in.
 * Nothing of the family being reconstructed is read.
 */
const measureDropIn = () => {
  const fams = [...new Set(files.map((f) => /_(?:E\d+_Bl\d+|Mother|Father|Sib\w*|Foetus|Grand\w+|Maternal\w+)_(PGD\d+)/.exec(f)?.[1]).filter(Boolean))]
  let obs = 0; let het = 0; let used = 0
  for (const fam of fams as string[]) {
    if (fam === FAM) continue
    const fa = files.find((f) => new RegExp(`Father_${fam}\\.`).test(f))
    const mo = files.find((f) => new RegExp(`Mother_${fam}\\.`).test(f))
    const bl = files.filter((f) => new RegExp(`E\\d+_Bl\\d+_${fam}(C\\d)?\\.`).test(f))
    if (!fa || !mo || !bl.length) continue
    const fg = gtOf(fa); const mg = gtOf(mo)
    const impossible = new Set<string>()
    for (const [probe, g] of fg) if (g !== 'AB' && mg.get(probe) === g) impossible.add(probe)
    if (impossible.size < 10_000) continue
    used += 1
    for (const b of bl) {
      for (const [probe, v] of gtOf(b)) {
        if (!impossible.has(probe)) continue
        obs += 1
        if (v === 'AB') het += 1
      }
    }
  }
  return { rate: obs ? het / obs : NaN, obs, families: used }
}
/**
 * The three rates the likelihood needs, measured on the OTHER families instead of modelled.
 *
 * The model assumed a blastomere truly heterozygous reads AB with probability 1 - dropout, with
 * dropout inferred from its own heterozygosity. On MDA single blastomeres that is wrong in three
 * places at once, and all three are measurable where both parents are arrayed: how often a
 * blastomere reads AB when the missing parent is homozygous for the allele the known parent lacks,
 * when she is heterozygous, and when she carries none of it. The last is drop-in.
 */
const calibrate = () => {
  const fams = [...new Set(files.map((f) => /_(?:E\d+_Bl\d+|Mother|Father|Sib\w*|Foetus|Grand\w+|Maternal\w+)_(PGD\d+)/.exec(f)?.[1]).filter(Boolean))]
  const acc = { w: [0, 0], h: [0, 0], n: [0, 0] }
  const marks = { w: 0, h: 0, n: 0 }
  let used = 0
  for (const fam of fams as string[]) {
    if (fam === FAM) continue
    const faF = files.find((f) => new RegExp(`Father_${fam}\\.`).test(f))
    const moF = files.find((f) => new RegExp(`Mother_${fam}\\.`).test(f))
    const bl = files.filter((f) => new RegExp(`E\\d+_Bl\\d+_${fam}(C\\d)?\\.`).test(f))
    if (!faF || !moF || !bl.length) continue
    const fg = gtOf(faF); const mg = gtOf(moF)
    const cls = new Map<string, 'w' | 'h' | 'n'>()
    for (const [probe, g] of fg) {
      if (g === 'AB') continue
      const t = mg.get(probe)
      if (!t) continue
      const k = t === 'AB' ? 'h' : t === g ? 'n' : 'w'
      cls.set(probe, k); marks[k] += 1
    }
    if (cls.size < 20_000) continue
    used += 1
    for (const b of bl) {
      for (const [probe, v] of gtOf(b)) {
        const k = cls.get(probe)
        if (!k) continue
        acc[k][0] += 1
        if (v === 'AB') acc[k][1] += 1
      }
    }
  }
  const r = (k: 'w' | 'h' | 'n') => acc[k][0] ? acc[k][1] / acc[k][0] : NaN
  const tot = marks.w + marks.h + marks.n
  return { pW: r('w'), pH: r('h'), pN: r('n'), families: used,
    obs: acc.w[0] + acc.h[0] + acc.n[0],
    priorW: marks.w / tot, priorH: marks.h / tot, priorN: marks.n / tot }
}
/**
 * PER-PROBE drop-in, which is the term the independence assumption gets wrong.
 *
 * The calibrated rates above are averages over probes. The likelihood then treats each embryo's
 * reading as independent given the hypothesis, which is false when the tendency to read
 * heterozygous belongs to the PROBE: such a probe reads AB in every blastomere of every embryo and
 * produces unanimity whatever the missing parent carries. That is the residual the earlier rounds
 * could not remove, and 11 clean trios of this series can measure it.
 *
 * Shrunk toward the global rate with a pseudo-count, because a probe seen fifteen times cannot
 * support an unshrunk estimate.
 */
const perProbeDropIn = (global: number, pseudo = 12) => {
  const fams = [...new Set(files.map((f) => /_(?:E\d+_Bl\d+|Mother|Father|Sib\w*|Foetus|Grand\w+|Maternal\w+)_(PGD\d+)/.exec(f)?.[1]).filter(Boolean))]
  const obs = new Map<string, number>()
  const het = new Map<string, number>()
  for (const fam of fams as string[]) {
    if (fam === FAM) continue
    const faF = files.find((f) => new RegExp(`Father_${fam}\\.`).test(f))
    const moF = files.find((f) => new RegExp(`Mother_${fam}\\.`).test(f))
    const bl = files.filter((f) => new RegExp(`E\\d+_Bl\\d+_${fam}(C\\d)?\\.`).test(f))
    if (!faF || !moF || !bl.length) continue
    const fg = gtOf(faF); const mg = gtOf(moF)
    const none = new Set<string>()
    for (const [probe, g] of fg) if (g !== 'AB' && mg.get(probe) === g) none.add(probe)
    if (none.size < 20_000) continue
    for (const b of bl) {
      for (const [probe, v] of gtOf(b)) {
        if (!none.has(probe)) continue
        obs.set(probe, (obs.get(probe) ?? 0) + 1)
        if (v === 'AB') het.set(probe, (het.get(probe) ?? 0) + 1)
      }
    }
  }
  const rate = (probe: string) => {
    const n = obs.get(probe) ?? 0
    return ((het.get(probe) ?? 0) + pseudo * global) / (n + pseudo)
  }
  // How concentrated is drop-in? If it were a per-observation coin, few probes would stand out.
  let seen = 0; let hot = 0
  for (const [probe, n] of obs) {
    if (n < 8) continue
    seen += 1
    if ((het.get(probe) ?? 0) / n >= 0.15) hot += 1
  }
  return { rate, obs, seen, hot }
}
const cal = calibrate()
console.log(`  calibrated on ${cal.families} other families, ${cal.obs.toLocaleString()} observations:`)
console.log(`    blastomere reads AB when she is homozygous for the lacking allele: ${cal.pW.toFixed(4)}`)
console.log(`    when she is heterozygous:                                         ${cal.pH.toFixed(4)}`)
console.log(`    when she carries none of it (drop-in):                            ${cal.pN.toFixed(4)}`)
console.log(`    likelihood ratio per heterozygous reading: ${(cal.pW / cal.pH).toFixed(2)}, `
  + `per homozygous reading ${((1 - cal.pW) / (1 - cal.pH)).toFixed(2)}`)
console.log(`    class shares among known-parent-homozygous markers: wanted ${cal.priorW.toFixed(4)}, `
  + `het ${cal.priorH.toFixed(4)}, none ${cal.priorN.toFixed(4)}`)
const measured = measureDropIn()
const pp = perProbeDropIn(cal.pN)
console.log(`  per-probe drop-in: ${pp.seen.toLocaleString()} probes with 8+ observations, `
  + `${pp.hot.toLocaleString()} reading 15% or more heterozygous `
  + `(${(100 * pp.hot / Math.max(1, pp.seen)).toFixed(2)}%)`)
console.log(`  drop-in measured on ${measured.families} other families of this series: `
  + `${measured.rate.toFixed(5)} over ${measured.obs.toLocaleString()} observations`)
console.log(`  (0.0435 was measured on bulk Axiom arrays and is what the earlier runs assumed)`)
const DROP_IN = Number.isFinite(measured.rate) && measured.rate > 0 ? measured.rate : 0.0435
const cellsOfEmbryo = new Map<string, string[]>()
for (const e of embryos) cellsOfEmbryo.set(e, cells.filter((f) => embryoOf.get(f) === e))

const MIN_EMBRYOS = Number(process.env.OM_MIN_EMBRYOS ?? 6)
const post = new Map<string, number>()
const spoke = new Map<string, number>()
for (const [probe, fg] of faGt) {
  if (fg === 'AB') continue
  const f2 = lackFreq(probe)
  if (!Number.isFinite(f2) || f2 <= 0 || f2 >= 1) continue
  // Hardy-Weinberg at this marker's own frequency, rescaled so the three sum to one, with the
  // measured genome-wide class shares as the check that the frequencies are not far out.
  let lW = Math.log(f2 * f2)
  let lH = Math.log(2 * f2 * (1 - f2))
  let lN = Math.log((1 - f2) * (1 - f2))
  let seen = 0
  // The three rates are per CELL and measured, so each hypothesis is scored directly against them
  // rather than through a dropout parameter. An embryo's cells are pooled within the embryo, which
  // is where the clonality lives; the mixture over embryos is what the heterozygous hypothesis is.
  for (const e of embryos) {
    let lIn = 0; let lOut = 0; let called = 0
    for (const f of cellsOfEmbryo.get(e)!) {
      const v = cellGt.get(f)!.gt.get(probe)
      if (!v) continue
      called += 1
      // "In" means this embryo inherited the lacking allele, so its true genotype is heterozygous
      // and it reads AB at the measured wanted-class rate. "Out" means it did not.
      // The out-rate is this PROBE's, which is what breaks the false unanimity. The in-rate stays
      // global: the wanted class carries about two observations per probe, too few to estimate.
      const pn = Math.min(0.9, Math.max(1e-4, pp.rate(probe)))
      if (v === 'AB') { lIn += Math.log(cal.pW); lOut += Math.log(pn) }
      else { lIn += Math.log(1 - cal.pW); lOut += Math.log(1 - pn) }
    }
    if (!called) continue
    seen += 1
    const m = Math.max(lIn, lOut)
    lW += lIn
    lH += m + Math.log(0.5 * Math.exp(lIn - m) + 0.5 * Math.exp(lOut - m))
    lN += lOut
  }
  // An embryo that called nothing contributes nothing to any of the three hypotheses, so the
  // likelihood ratio is 2^seen and the posterior already reflects how many spoke. The gate is only
  // a floor, and the posterior threshold does the rest.
  if (seen < MIN_EMBRYOS) continue
  spoke.set(probe, seen)
  const m = Math.max(lW, lH, lN)
  const d = Math.exp(lW - m) + Math.exp(lH - m) + Math.exp(lN - m)
  post.set(probe, Math.exp(lW - m) / d)
}
{
  const hist = new Map<number, number>()
  for (const v of spoke.values()) hist.set(v, (hist.get(v) ?? 0) + 1)
  const ks = [...hist.keys()].sort((a, b) => b - a)
  console.log(`\n${post.size.toLocaleString()} markers scored with at least ${MIN_EMBRYOS} of `
    + `${embryos.length} embryos calling`)
  console.log(`  embryos calling: ${ks.slice(0, 9).map((k) => `${k}:${hist.get(k)!.toLocaleString()}`).join('  ')}`)
}

const chromOf = new Map<string, string>()
for (const r of rowsOf(join(DIR, fatherFile))()) {
  if (/^\d+$/.test(r.chrom)) chromOf.set(r.probesetId, r.chrom)
}
const grade = (asserted: Set<string>) => {
  let hit = 0; let herHet = 0; let sameAsKnown = 0
  for (const probe of asserted) {
    const t = moGt.get(probe)
    if (!t) continue
    const want = faGt.get(probe) === 'AA' ? 'BB' : 'AA'
    if (t === want) hit += 1
    else if (t === 'AB') herHet += 1
    else sameAsKnown += 1
  }
  const scored = hit + herHet + sameAsKnown
  let trueW = 0; let found = 0
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    const want = fg === 'AA' ? 'BB' : 'AA'
    if (moGt.get(probe) !== want) continue
    trueW += 1
    if (asserted.has(probe)) found += 1
  }
  const byC = new Map<string, number>()
  for (const probe of asserted) {
    const c = chromOf.get(probe)
    if (c) byC.set(c, (byC.get(c) ?? 0) + 1)
  }
  const cs = [...byC.values()].sort((a, b) => a - b)
  return { precision: scored ? hit / scored : NaN, scored, herHet, sameAsKnown,
    completeness: trueW ? found / trueW : NaN, trueW, chroms: byC.size,
    med: cs.length ? cs[Math.floor(cs.length / 2)] : 0, min: cs[0] ?? 0 }
}

console.log(`\n  ${'posterior'.padStart(10)} ${'asserted'.padStart(9)} ${'PRECISION'.padStart(10)} `
  + `${'completeness'.padStart(12)} ${'per chrom'.padStart(10)} ${'autosomes'.padStart(10)}`)
let win: { th: number; pr: number; n: number; med: number; comp: number; chroms: number } | undefined
for (const th of [0.5, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999]) {
  const a = new Set<string>()
  for (const [probe, v] of post) if (v >= th) a.add(probe)
  const s = grade(a)
  if (!s.scored) continue
  console.log(`  ${th.toFixed(4).padStart(10)} ${s.scored.toLocaleString().padStart(9)} `
    + `${s.precision.toFixed(4).padStart(10)} ${s.completeness.toFixed(4).padStart(12)} `
    + `${s.med.toLocaleString().padStart(10)} ${String(s.chroms).padStart(10)}`)
  if (s.precision >= 0.95 && s.med >= 200 && (!win || s.scored > win.n)) {
    win = { th, pr: s.precision, n: s.scored, med: s.med, comp: s.completeness, chroms: s.chroms }
  }
}

// The error breakdown, which says which of the two bounds is still binding.
{
  const a = new Set<string>()
  for (const [probe, v] of post) if (v >= 0.99) a.add(probe)
  const s = grade(a)
  console.log(`\n  at a posterior of 0.99, of ${s.scored.toLocaleString()} asserted markers:`)
  console.log(`    she reads the wanted homozygote:        ${(s.scored - s.herHet - s.sameAsKnown).toLocaleString()}`)
  console.log(`    she reads HETEROZYGOUS:                 ${s.herHet.toLocaleString()}   information limit`)
  console.log(`    she reads the SAME as the known parent: ${s.sameAsKnown.toLocaleString()}   artefact`)
}

// WHAT MATERIAL WOULD PASS. Every rate is now measured, so the requirement follows from them
// rather than from a guess. With a likelihood ratio of pW/pN per heterozygous reading and
// (1-pW)/(1-pN) per homozygous one, a marker clears a precision target when enough of the embryos
// that CALLED it read heterozygous; how often that happens is binomial in the wanted class.
{
  // The hypothesis that COMPETES is that she is heterozygous, which is a mixture: each embryo is
  // independently in or out with probability one half. So the per-embryo ratio is against that
  // mixture and not against "she carries none of it", which is the far weaker rival. Comparing
  // against the wrong rival predicted 580 markers per chromosome where the measurement gives 10.
  const pAbHet = 0.5 * cal.pW + 0.5 * cal.pN
  const lrIn = Math.log(cal.pW / pAbHet)
  const lrOut = Math.log((1 - cal.pW) / (1 - pAbHet))
  const priorOdds = cal.priorH / cal.priorW
  const need = Math.log(0.95 / 0.05 * priorOdds)
  const lchoose = (n: number, k: number) => {
    let x = 0
    for (let i = 1; i <= k; i += 1) x += Math.log(n - k + i) - Math.log(i)
    return x
  }
  const tail = (n: number, pr: number, k: number) => {
    let x = 0
    for (let i = k; i <= n; i += 1) x += Math.exp(lchoose(n, i) + i * Math.log(pr) + (n - i) * Math.log(1 - pr))
    return Math.min(1, Math.max(0, x))
  }
  console.log(`\n  WHAT MATERIAL WOULD CLEAR 0.95, from the measured rates`)
  console.log(`  prior odds against the wanted genotype ${priorOdds.toFixed(2)}, so the likelihood `
    + `ratio must reach ${Math.exp(need).toFixed(1)}`)
  console.log(`  ${'embryos calling'.padStart(15)} ${'heterozygous needed'.padStart(20)} `
    + `${'share of wanted markers'.padStart(24)} ${'markers per chromosome'.padStart(23)}`)
  // The wanted class holds this many markers per autosome when every marker is usable.
  const wantedPerChrom = (cal.priorW * faGt.size) / 22
  console.log(`  per-embryo ratio against the heterozygous rival: `
    + `${Math.exp(lrIn).toFixed(3)} per heterozygous reading, ${Math.exp(lrOut).toFixed(3)} per `
    + 'homozygous one')
  for (const n of [8, 12, 20, 30, 40, 50, 60]) {
    let k = 0
    while (k <= n && k * lrIn + (n - k) * lrOut < need) k += 1
    if (k > n) { console.log(`  ${String(n).padStart(15)} ${'unreachable'.padStart(20)}`); continue }
    const share = tail(n, cal.pW, k)
    console.log(`  ${String(n).padStart(15)} ${`${k} of ${n}`.padStart(20)} `
      + `${share.toFixed(4).padStart(24)} ${(share * wantedPerChrom).toFixed(0).padStart(23)}`)
  }
  console.log(`  the channel refuses below 200 markers per chromosome.`)
  const meanCall = [...cellGt.values()].reduce((a, c) => a + c.callRate, 0) / cellGt.size
  console.log(`  this family's mean blastomere call rate is ${meanCall.toFixed(2)}, so having N `
    + `embryos call at one marker takes about ${(1 / meanCall).toFixed(1)}N embryos arrayed.`)
}

console.log(win
  ? `\n  PASS BOTH CRITERIA: precision ${win.pr.toFixed(4)} at a posterior of ${win.th}, `
    + `${win.n.toLocaleString()} markers over ${win.chroms} autosomes, median ${win.med.toLocaleString()} `
    + `per chromosome against the 200 the channel refuses below, completeness ${win.comp.toFixed(4)}.`
  : '\n  still short of one or both bars.')

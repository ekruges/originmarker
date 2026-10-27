/**
 * The reconstruction with LINKAGE, which is the term every earlier configuration threw away.
 *
 * Ten configurations treated each marker on its own. At a marker where the known parent is
 * homozygous, the rival hypothesis is that the missing parent is heterozygous, and marker by marker
 * that rival is a mixture over which embryos inherited which of her two haplotypes. The mixture is
 * what costs the power: the per-embryo likelihood ratio against it is 1.891, where the ratio
 * against a parent carrying none of the allele is 17.2.
 *
 * But she has only TWO haplotypes, and an embryo that inherited one carries it across a whole
 * region. So the split of embryos between her haplotypes is a property of the REGION, estimable
 * from every marker in it at once. With the split known the rival stops being a mixture: under
 * "she is heterozygous" the embryos on one side read heterozygous and the embryos on the other read
 * drop-in, and the discriminating evidence is heterozygous calls among embryos on the side that
 * should not have them. That evidence carries the 17.2.
 *
 * Her own array is held out and used only to score. Rates come from audit/pgd008-rebuild.ts's
 * calibration on the other families of the series.
 *
 * Run: OM_DIR=<dir> OM_FAMILY=PGD008 node --experimental-strip-types audit/pgd008-phased.ts
 */
import { readFileSync, existsSync, readdirSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)

const DIR = process.env.OM_DIR
const FAM = process.env.OM_FAMILY ?? 'PGD008'
const WINDOW = Number(process.env.OM_WINDOW ?? 400)
if (!DIR || !existsSync(DIR)) { console.log('OM_DIR not readable. NOT RUN.'); process.exit(0) }

function rowsOf(path: string) {
  const raw = path.endsWith('.gz') ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error(`${path}: unreadable header`)
  return function* () {
    for (let i = h + 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r }
  }
}
const gtOf = (f: string) => {
  const m = new Map<string, string>()
  for (const r of rowsOf(join(DIR, f))()) {
    if (!/^\d+$/.test(r.chrom)) continue
    if (r.genotype === 'AA' || r.genotype === 'BB' || r.genotype === 'AB') m.set(r.probesetId, r.genotype)
  }
  return m
}
const files = readdirSync(DIR).filter((f) => f.endsWith('.probes.gz'))
const fatherFile = files.find((f) => new RegExp(`Father_${FAM}\\.`).test(f))!
const motherFile = files.find((f) => new RegExp(`Mother_${FAM}\\.`).test(f))!
const cells = files.filter((f) => new RegExp(`E\\d+_Bl\\d+_${FAM}C\\d`).test(f))
const embryoOf = new Map<string, string>()
for (const f of cells) {
  const m = new RegExp(`(E\\d+)_Bl\\d+_(${FAM}C\\d)`).exec(f)
  if (m) embryoOf.set(f, `${m[2]}:${m[1]}`)
}
const embryos = [...new Set(embryoOf.values())].sort()
const faGt = gtOf(fatherFile)
const moGt = gtOf(motherFile)
const cellGt = new Map<string, Map<string, string>>()
for (const f of cells) cellGt.set(f, gtOf(f))
console.log(`${FAM}: ${embryos.length} embryos, ${cells.length} arrays, window ${WINDOW} markers\n`)

// Rates measured on the other eleven families by audit/pgd008-rebuild.ts.
const P_IN = 0.3604      // an embryo that inherited the lacking allele reads heterozygous
const P_OUT = 0.0209     // one that did not, by drop-in

// Marker order along each chromosome, restricted to known-parent-homozygous markers.
type Mk = { probe: string; chrom: string; pos: number }
const marks: Mk[] = []
for (const r of rowsOf(join(DIR, fatherFile))()) {
  if (!/^\d+$/.test(r.chrom)) continue
  const g = faGt.get(r.probesetId)
  if (!g || g === 'AB') continue
  marks.push({ probe: r.probesetId, chrom: r.chrom, pos: r.pos })
}
marks.sort((a, b) => (a.chrom === b.chrom ? a.pos - b.pos : Number(a.chrom) - Number(b.chrom)))
console.log(`${marks.length.toLocaleString()} known-parent-homozygous markers in order`)

// The frequency panel, for the per-marker prior. Unrelated donors only.
const panel = files.filter((f) => /(Mother|Father)_PGD/.test(f) && !f.includes(FAM))
const bC = new Map<string, number>(); const nC = new Map<string, number>()
for (const f of panel) {
  for (const [probe, v] of gtOf(f)) {
    if (!faGt.has(probe)) continue
    bC.set(probe, (bC.get(probe) ?? 0) + (v === 'BB' ? 2 : v === 'AB' ? 1 : 0))
    nC.set(probe, (nC.get(probe) ?? 0) + 2)
  }
}
const lackFreq = (probe: string) => {
  const fg = faGt.get(probe)
  if (!fg || fg === 'AB') return NaN
  const b = bC.get(probe); const n = nC.get(probe)
  if (b === undefined || !n) return NaN
  return fg === 'AA' ? b / n : (n - b) / n
}

/** An embryo reads heterozygous at a marker: 1 yes, 0 no, -1 not called. */
const readAt = (e: string, probe: string): number => {
  let anyCalled = false
  for (const f of cells) {
    if (embryoOf.get(f) !== e) continue
    const v = cellGt.get(f)!.get(probe)
    if (!v) continue
    anyCalled = true
    if (v === 'AB') return 1
  }
  return anyCalled ? 0 : -1
}

const post = new Map<string, number>()
let windows = 0; let split2 = 0
const diag = { hi: 0, lo: 0, n: 0 }
for (let start = 0; start < marks.length; start += WINDOW) {
  const win = marks.slice(start, start + WINDOW).filter((m, i, a) => a[0].chrom === m.chrom)
  if (win.length < 50) continue
  windows += 1
  const obs = win.map((m) => embryos.map((e) => readAt(e, m.probe)))

  // THE SPLIT. At a marker where she is heterozygous, roughly half the embryos inherited the
  // lacking allele, so markers with an intermediate heterozygous share carry her phase. Embryos on
  // the same haplotype read heterozygous together at those markers.
  const phaseRows = obs.filter((row) => {
    const called = row.filter((x) => x >= 0).length
    const het = row.filter((x) => x === 1).length
    return called >= Math.max(4, 0.4 * embryos.length) && het >= 1 && het <= called - 1
  })
  if (phaseRows.length < 10) continue
  const n = embryos.length
  // THE AFFINITY THAT CARRIES PHASE is co-occurrence of HETEROZYGOUS calls, not overall agreement.
  // Two embryos on the same maternal haplotype both read heterozygous with probability
  // P_IN squared, which is 0.1298. On opposite haplotypes one of them is reading drop-in, so
  // P_IN times P_OUT, which is 0.0075. That is a 17-fold separation. Overall agreement runs the
  // other way and is useless here: at a P_IN of 0.36 the "both homozygous" case dominates, so
  // opposite-haplotype pairs agree MORE often (0.64) than same-haplotype pairs (0.54), and seeding
  // on least agreement picks the wrong pair.
  const coAb = Array.from({ length: n }, () => new Array(n).fill(0))
  const bothCalled = Array.from({ length: n }, () => new Array(n).fill(0))
  for (const row of phaseRows) {
    for (let i = 0; i < n; i += 1) {
      if (row[i] < 0) continue
      for (let j = i + 1; j < n; j += 1) {
        if (row[j] < 0) continue
        bothCalled[i][j] += 1
        if (row[i] === 1 && row[j] === 1) coAb[i][j] += 1
      }
    }
  }
  const aff = (i: number, j: number) => {
    const a = i < j ? i : j; const b = i < j ? j : i
    return bothCalled[a][b] >= 10 ? coAb[a][b] / bothCalled[a][b] : NaN
  }
  // Seed group 0 with the pair of highest co-heterozygosity, which share a haplotype.
  let hi = -1; let hj = -1; let best = -1
  for (let i = 0; i < n; i += 1) {
    for (let j = i + 1; j < n; j += 1) {
      const a = aff(i, j)
      if (Number.isFinite(a) && a > best) { best = a; hi = i; hj = j }
    }
  }
  if (hi < 0 || best < 0.06) continue
  // Seed group 1 with whichever embryo shares heterozygosity least with that pair.
  let oi = -1; let low = 2
  for (let k = 0; k < n; k += 1) {
    if (k === hi || k === hj) continue
    const a0 = aff(k, hi); const a1 = aff(k, hj)
    const v = [a0, a1].filter(Number.isFinite) as number[]
    if (!v.length) continue
    const mean = v.reduce((x, y) => x + y, 0) / v.length
    if (mean < low) { low = mean; oi = k }
  }
  if (oi < 0) continue
  split2 += 1
  // ITERATIVE REASSIGNMENT, because fixed seeds bias the result. With two seeds on one side and
  // one on the other, every embryo compares a mean over two against a mean over one, and the
  // two-member side wins ties. Back-solving the measured rates showed the cost: the richer side
  // came out 36 percent contaminated with embryos that belong on the other, while the poorer side
  // was 1.3 percent. A misassigned embryo that did inherit the allele reads heterozygous at P_IN on
  // the side that should be reading drop-in, which is indistinguishable from the parent being
  // homozygous for it, so split error MANUFACTURES the signal this is testing for.
  const side = new Array<number>(n).fill(-1)
  side[hi] = 0; side[hj] = 0; side[oi] = 1
  const meanAff = (k: number, group: number) => {
    const v: number[] = []
    for (let m2 = 0; m2 < n; m2 += 1) {
      if (m2 === k || side[m2] !== group) continue
      const a = aff(k, m2)
      if (Number.isFinite(a)) v.push(a)
    }
    return v.length ? v.reduce((x, y) => x + y, 0) / v.length : NaN
  }
  for (let k = 0; k < n; k += 1) {
    if (side[k] >= 0) continue
    const m0 = meanAff(k, 0); const m1 = meanAff(k, 1)
    if (!Number.isFinite(m0) && !Number.isFinite(m1)) continue
    side[k] = (Number.isFinite(m0) ? m0 : -1) >= (Number.isFinite(m1) ? m1 : -1) ? 0 : 1
  }
  // Now refine every embryo, seeds included, until nothing moves. Each embryo is compared on the
  // same footing, so neither side is favoured by how it was seeded.
  for (let pass = 0; pass < 12; pass += 1) {
    let moved = 0
    for (let k = 0; k < n; k += 1) {
      if (side[k] < 0) continue
      const m0 = meanAff(k, 0); const m1 = meanAff(k, 1)
      if (!Number.isFinite(m0) || !Number.isFinite(m1)) continue
      const want = m0 >= m1 ? 0 : 1
      // Never empty a side.
      if (want !== side[k] && side.filter((x) => x === side[k]).length > 1) {
        side[k] = want; moved += 1
      }
    }
    if (!moved) break
  }

  // SPLIT QUALITY, against the held-out array. Diagnosis only: at a marker where she is genuinely
  // heterozygous, the side that inherited the lacking allele should read near P_IN and the other
  // near P_OUT. A split that is no better than chance leaves both near their average.
  for (let wi = 0; wi < win.length; wi += 1) {
    const t = moGt.get(win[wi].probe)
    if (t !== 'AB') continue
    const row = obs[wi]
    let a0 = 0; let n0 = 0; let a1 = 0; let n1 = 0
    for (let k = 0; k < n; k += 1) {
      if (row[k] < 0 || side[k] < 0) continue
      if (side[k] === 0) { n0 += 1; if (row[k] === 1) a0 += 1 } else { n1 += 1; if (row[k] === 1) a1 += 1 }
    }
    if (!n0 || !n1) continue
    // Record the HIGHER and LOWER side, since which label carries the allele flips per marker.
    const hiR = Math.max(a0 / n0, a1 / n1); const loR = Math.min(a0 / n0, a1 / n1)
    diag.hi += hiR; diag.lo += loR; diag.n += 1
  }

  // WITH THE SPLIT KNOWN the rival is no longer a mixture. Three hypotheses per marker:
  //   she is homozygous for the lacking allele -> EVERY embryo inherited it, both sides read P_IN
  //   she is heterozygous                      -> one side reads P_IN, the other P_OUT
  //   she carries none of it                   -> both sides read P_OUT
  const lp = (x: number, p: number) => (x === 1 ? Math.log(p) : Math.log(1 - p))
  for (let wi = 0; wi < win.length; wi += 1) {
    const probe = win[wi].probe
    const row = obs[wi]
    const f2 = lackFreq(probe)
    if (!Number.isFinite(f2) || f2 <= 0 || f2 >= 1) continue
    let called = 0
    for (let k = 0; k < n; k += 1) if (row[k] >= 0 && side[k] >= 0) called += 1
    if (called < 6) continue
    let lW = Math.log(f2 * f2)
    let lN = Math.log((1 - f2) * (1 - f2))
    // Both phase orientations of the heterozygous hypothesis, since which side carries the allele
    // is not fixed by the split alone.
    let lH0 = Math.log(f2 * (1 - f2))
    let lH1 = Math.log(f2 * (1 - f2))
    for (let k = 0; k < n; k += 1) {
      if (row[k] < 0 || side[k] < 0) continue
      lW += lp(row[k], P_IN)
      lN += lp(row[k], P_OUT)
      lH0 += lp(row[k], side[k] === 0 ? P_IN : P_OUT)
      lH1 += lp(row[k], side[k] === 1 ? P_IN : P_OUT)
    }
    const m = Math.max(lW, lN, lH0, lH1)
    const d = Math.exp(lW - m) + Math.exp(lN - m) + Math.exp(lH0 - m) + Math.exp(lH1 - m)
    post.set(probe, Math.exp(lW - m) / d)
  }
}
console.log(`${windows} windows, ${split2} with a usable haplotype split`)
console.log(`  co-heterozygosity expected ${(P_IN * P_IN).toFixed(4)} within a haplotype and `
  + `${(P_IN * P_OUT).toFixed(4)} across, a ${(P_IN / P_OUT).toFixed(0)}-fold separation`)
console.log(`${post.size.toLocaleString()} markers scored`)
if (diag.n) {
  const hi = diag.hi / diag.n; const lo = diag.lo / diag.n
  console.log(`\nSPLIT QUALITY at her ${diag.n.toLocaleString()} genuinely heterozygous markers:`)
  console.log(`  the richer side reads ${hi.toFixed(4)} heterozygous, the poorer side ${lo.toFixed(4)}`)
  console.log(`  a correct split would read about ${P_IN.toFixed(3)} and ${P_OUT.toFixed(3)}; `
    + `a split no better than chance would read both near `
    + `${(0.5 * (P_IN + P_OUT)).toFixed(3)} before the max-min spread`)
  console.log(`  separation achieved: ${(hi / Math.max(1e-9, lo)).toFixed(2)}-fold against the `
    + `${(P_IN / P_OUT).toFixed(0)}-fold a correct split gives`)
}
console.log()

const chromOf = new Map<string, string>()
for (const m of marks) chromOf.set(m.probe, m.chrom)
const grade = (a: Set<string>) => {
  let hit = 0; let herHet = 0; let same = 0
  for (const probe of a) {
    const t = moGt.get(probe)
    if (!t) continue
    const want = faGt.get(probe) === 'AA' ? 'BB' : 'AA'
    if (t === want) hit += 1; else if (t === 'AB') herHet += 1; else same += 1
  }
  const scored = hit + herHet + same
  let trueW = 0; let found = 0
  for (const [probe, fg] of faGt) {
    if (fg === 'AB') continue
    if (moGt.get(probe) !== (fg === 'AA' ? 'BB' : 'AA')) continue
    trueW += 1
    if (a.has(probe)) found += 1
  }
  const byC = new Map<string, number>()
  for (const probe of a) { const c = chromOf.get(probe); if (c) byC.set(c, (byC.get(c) ?? 0) + 1) }
  const cs = [...byC.values()].sort((x, y) => x - y)
  return { precision: scored ? hit / scored : NaN, scored, herHet, same,
    completeness: trueW ? found / trueW : NaN, chroms: byC.size,
    med: cs.length ? cs[Math.floor(cs.length / 2)] : 0 }
}
console.log(`  ${'posterior'.padStart(10)} ${'asserted'.padStart(9)} ${'PRECISION'.padStart(10)} `
  + `${'completeness'.padStart(12)} ${'per chrom'.padStart(10)} ${'autosomes'.padStart(10)}`)
let win2: { th: number; pr: number; n: number; med: number; comp: number; chroms: number } | undefined
for (const th of [0.5, 0.9, 0.95, 0.99, 0.995, 0.999]) {
  const a = new Set<string>()
  for (const [probe, v] of post) if (v >= th) a.add(probe)
  const s = grade(a)
  if (!s.scored) continue
  console.log(`  ${th.toFixed(4).padStart(10)} ${s.scored.toLocaleString().padStart(9)} `
    + `${s.precision.toFixed(4).padStart(10)} ${s.completeness.toFixed(4).padStart(12)} `
    + `${s.med.toLocaleString().padStart(10)} ${String(s.chroms).padStart(10)}`)
  if (s.precision >= 0.95 && s.med >= 200 && (!win2 || s.scored > win2.n)) {
    win2 = { th, pr: s.precision, n: s.scored, med: s.med, comp: s.completeness, chroms: s.chroms }
  }
}
console.log(win2
  ? `\n  PASS BOTH CRITERIA: precision ${win2.pr.toFixed(4)} at a posterior of ${win2.th}, `
    + `${win2.n.toLocaleString()} markers over ${win2.chroms} autosomes, median ${win2.med.toLocaleString()} `
    + `per chromosome against the 200 the channel refuses below, completeness ${win2.comp.toFixed(4)}.`
  : '\n  still short of one or both bars.')

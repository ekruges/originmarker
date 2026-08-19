/**
 * Measure the accuracy of EVERY grade the tool can emit, on real arrays, against known truth.
 *
 * WHY THIS HAD TO BE REBUILT. `audit/bands_measured.csv` measured band D as everything below 0.75,
 * because that was the whole tail when it was made. A grade F was later carved out of that tail at
 * 0.55, and `BAND_ACCURACY` went on reporting the pooled figure for D. So the shipped D accuracy
 * was measured over a range that now includes F: the real D must be better than the number beside
 * it and F must be worse, and neither was true. A grade that cannot state its own accuracy is not
 * a grade, it is a label.
 *
 * HOW TRUTH IS ESTABLISHED. An injection: take a real array, choose a region, and displace the
 * allele dosage there by exactly the amount a known parent's copy being affected would produce,
 * under a known class at a known mosaic fraction. `shiftMean` is the tool's own forward model, so
 * the injection and the caller agree on what the world would look like; what is being measured is
 * whether the caller recovers the parent it was given, not whether the model is right. Noise comes
 * from the array itself rather than from a distribution, which is the point of using real files.
 *
 * CLUSTERED BY ARRAY. Rows from one array are not independent: the same amplification, the same
 * dropout structure, the same donor. A naive interval over pooled rows was measured 3.5 to 5.8x too
 * narrow on this material, so every interval here is over the per-array accuracies.
 */
import { readFileSync, writeFileSync, readdirSync } from 'node:fs'
import { join } from 'node:path'
import { homedir } from 'node:os'
import { shiftMean, bandOf, type EventClass } from '../web/src/originPosterior.ts'
import { callDosageOrigin, WINDOW_LO, WINDOW_HI, type Material } from '../web/src/dosageOrigin.ts'

const DIR = process.argv[2] ?? join(homedir(), 'Downloads', 'probes')
// A SEPARATE FILE from bands_measured.csv, which holds the original series. This one's injections
// use the caller's own forward model, so its A, B and C figures carry no model misspecification and
// are optimistic; what it establishes is the comparison WITHIN itself, where D and F differ by 0.2
// and F sits below chance. Overwriting the original with it would silently replace a
// harder-measured table with an easier one.
const OUT = join(import.meta.dirname, 'bands_measured_f_series.csv')

/** Deterministic, so a reported accuracy can be reproduced from the seed alone. */
const rng = (seed: number) => () => {
  seed = (seed * 1664525 + 1013904223) >>> 0
  return seed / 4294967296
}

/**
 * THE MARKER SET THE CHANNEL ACTUALLY READS, which the first version of this harness got backwards.
 *
 * The centroid is taken at sites where the LOADED PARENT is homozygous, over the sample's B-allele
 * frequency, restricted to [0.20, 0.80]. At a site where the parent is AA, a sample carrying both
 * parental copies reads near 0.5 exactly when the other parent transmitted a B, so the window
 * selects the sites where the sample is HETEROZYGOUS. Selecting the sample's own homozygous calls
 * instead, as the first version did, picks the complementary set: every one of them falls outside
 * the window and nothing informative survives. It produced band A accuracies of 0.09 to 0.38, worse
 * than guessing, which is the signature of a harness fault rather than a caller fault.
 */
type Site = { parent: 'AA' | 'BB'; baf: number }

/** The loaded parent's genotype per marker, keyed by probeset. */
function readParent(path: string): Map<string, 'AA' | 'BB'> {
  const out = new Map<string, 'AA' | 'BB'>()
  const txt = readFileSync(path, 'utf8')
  let i = txt.indexOf('\n') + 1
  while (i < txt.length) {
    const j = txt.indexOf('\n', i)
    if (j < 0) break
    const line = txt.slice(i, j); i = j + 1
    const p = line.split('\t')
    if (p.length < 7) continue
    const c = p[1]
    if (c === 'X' || c === 'Y' || c === 'MT') continue
    if (p[6] === '0') out.set(p[0], 'AA')
    else if (p[6] === '2') out.set(p[0], 'BB')
  }
  return out
}

/** The sample's dosage at those same sites, with the sample's own call rate recovered. */
function readSample(path: string, parent: Map<string, 'AA' | 'BB'>) {
  const sites: Site[] = []
  let total = 0
  let called = 0
  let het = 0
  const txt = readFileSync(path, 'utf8')
  let i = txt.indexOf('\n') + 1
  while (i < txt.length) {
    const j = txt.indexOf('\n', i)
    if (j < 0) break
    const line = txt.slice(i, j); i = j + 1
    const p = line.split('\t')
    if (p.length < 7) continue
    const c = p[1]
    if (c === 'X' || c === 'Y' || c === 'MT') continue
    total += 1
    if (p[6] !== '-1') { called += 1; if (p[6] === '1') het += 1 }
    const pg = parent.get(p[0])
    if (!pg) continue
    const baf = Number(p[4])
    if (!Number.isFinite(baf)) continue
    sites.push({ parent: pg, baf })
  }
  return { sites, callRate: total ? called / total : 0, het: called ? het / called : 0 }
}

const CLASSES: EventClass[] = ['loss', 'gain', 'cnn-loh']
const FRACTIONS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60]
const REGION_MARKERS = 600
const PER_ARRAY = 3000

/** Which material a file is, from its own quality. The bands are measured per material. */
function materialOf(callRate: number, het: number): Material {
  if (callRate > 0.95 && het > 0.14) return 'bulk'
  if (het > 0.12) return 'trophectoderm'
  if (callRate > 0.86) return 'esc-single'
  return 'blastomere'
}

const files = readdirSync(DIR).filter((f) => f.endsWith('.probes')).sort()
// The loaded parent. This series is measured against a real parent array rather than a stand-in,
// because the informative marker set is defined BY that parent's homozygous sites.
// The parent array is named on the command line, or is the one whose call rate and heterozygosity
// mark it as bulk diploid. It is never hard-coded: sample identifiers do not belong in this repo.
const PARENT_FILE = process.argv[3] ?? (() => {
  let best = files[0]
  let bestScore = -1
  for (const f of files) {
    const txt = readFileSync(join(DIR, f), 'utf8')
    let called = 0
    let het = 0
    let total = 0
    for (const line of txt.split('\n').slice(1, 200_000)) {
      const p3 = line.split('\t')
      if (p3.length < 7) continue
      total += 1
      if (p3[6] === '-1') continue
      called += 1
      if (p3[6] === '1') het += 1
    }
    const cr = total ? called / total : 0
    const h = called ? het / called : 0
    // Bulk diploid: nearly everything called, heterozygosity at the population rate.
    const score = cr > 0.95 && h > 0.14 ? cr * h : -1
    if (score > bestScore) { bestScore = score; best = f }
  }
  return best
})()
const parent = readParent(join(DIR, PARENT_FILE))
process.stderr.write(`parent ${PARENT_FILE}: ${parent.size} homozygous autosomal sites\n`)

const perArray = new Map<string, Map<string, { n: number; right: number }>>()
let seedBase = 20260819

for (const f of files) {
  if (f === PARENT_FILE) continue
  const { sites, callRate, het } = readSample(join(DIR, f), parent)
  // Only sites the caller would keep: parent homozygous AND the sample reading in the central band.
  const informative = sites.filter((s2) => {
    const o = s2.parent === 'BB' ? 1 - s2.baf : s2.baf
    return o >= WINDOW_LO && o <= WINDOW_HI
  })
  if (informative.length < REGION_MARKERS * 4) {
    process.stderr.write(`${f}  skipped, only ${informative.length} informative sites\n`)
    continue
  }
  const material = materialOf(callRate, het)
  const next = rng(seedBase += 1)
  const byBand = new Map<string, { n: number; right: number }>()
  const background: [string, number][] = informative.map((s2) => [s2.parent, s2.baf])

  for (let t = 0; t < PER_ARRAY; t += 1) {
    const cls = CLASSES[Math.floor(next() * CLASSES.length)]
    const fr = FRACTIONS[Math.floor(next() * FRACTIONS.length)]
    const affected = next() < 0.5 ? 'loaded' : 'other'
    // The displacement that parent's copy being affected produces, in ORIENTED space, applied back
    // through the orientation so the caller sees an ordinary array.
    const delta = shiftMean(cls, fr, affected)
    const start = Math.floor(next() * (informative.length - REGION_MARKERS))
    const region: [string, number][] = []
    for (let k = 0; k < REGION_MARKERS; k += 1) {
      const s2 = informative[start + k]
      const o = (s2.parent === 'BB' ? 1 - s2.baf : s2.baf) + delta
      region.push([s2.parent, s2.parent === 'BB' ? 1 - o : o])
    }
    const call = callDosageOrigin(region as never, background as never, material,
      { wholeChromosome: true, state: cls, parents: 1 })
    const p2 = call.posterior
    if (!p2 || p2.parent === 'withheld') continue
    const band = p2.band ?? bandOf(p2.confidence)
    const trueVerdict = affected === 'loaded' ? 'loaded-parent' : 'other-parent'
    const cur = byBand.get(`${material}|${band}`) ?? { n: 0, right: 0 }
    cur.n += 1
    if (call.verdict === trueVerdict) cur.right += 1
    byBand.set(`${material}|${band}`, cur)
  }
  perArray.set(f, byBand)
  const tot = [...byBand.values()].reduce((a2, x) => a2 + x.n, 0)
  process.stderr.write(`${f}  ${material}  ${informative.length} informative  ${tot} calls\n`)
}

// Cluster by array: the accuracy of each array, then the interval over those.
const keys = new Set<string>()
for (const m of perArray.values()) for (const k of m.keys()) keys.add(k)
const rows: string[] = ['material,band,n,n_arrays,acc,lo,hi']
for (const key of [...keys].sort()) {
  const [material, band] = key.split('|')
  const accs: number[] = []
  let n = 0
  let right = 0
  for (const m of perArray.values()) {
    const v = m.get(key)
    if (!v || v.n < 30) continue
    accs.push(v.right / v.n)
    n += v.n
    right += v.right
  }
  if (accs.length < 2) continue
  const mean = accs.reduce((a, x) => a + x, 0) / accs.length
  const sd = Math.sqrt(accs.reduce((a, x) => a + (x - mean) ** 2, 0) / (accs.length - 1))
  const se = sd / Math.sqrt(accs.length)
  rows.push([material, band, n, accs.length, (right / n).toFixed(4),
    Math.max(0, mean - 1.96 * se).toFixed(4), Math.min(1, mean + 1.96 * se).toFixed(4)].join(','))
}
writeFileSync(OUT, `${rows.join('\n')}\n`)
process.stderr.write(`\nwrote ${OUT}\n`)
console.log(rows.join('\n'))

/**
 * The zygosity boundary, measured against arrays whose ploidy is known from the pedigree.
 *
 * WHAT IS BEING DECIDED. `parentage.classify` reads whether a genome carries one parental
 * contribution or two from the fraction of B-allele frequencies in the heterozygous band. A flat
 * boundary on that fraction was calibrated on unamplified DNA, and on amplified material it fails
 * in BOTH directions at once, so a fix tested one way looks finished while half the problem
 * remains. This harness reports both directions on every run.
 *
 * GROUND TRUTH IS THE PEDIGREE, NEVER THE TOOL.
 *   Two contributions: a child of a COMPLETE trio, on amplified material. Both donors are
 *     identified, so the embryo has a maternal and a paternal complement by construction.
 *   One contribution: the 15 pronuclei of GSE148488. A pronucleus is at one copy on every
 *     autosome, which makes it an exact one-parent genome with nothing injected.
 *
 * THE WRONG-PARENT ARM. Zygosity is a property of the sample alone, so loading an unrelated adult
 * must not move it. That is asserted here rather than assumed: the same children are run with
 * their own father and with a donor who is not their parent, and the two must agree exactly. A
 * statistic that moved would be reading the parent rather than the genome.
 *
 * Data: GSE148488 (Cell 2021, doi 10.1016/j.cell.2020.10.025). Public, so the numbers this prints
 * can be committed and re-run. Point OM_TRIOS at a directory of <GSM>.probes.gz plus
 * trio_manifest_full.csv.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types audit/zygosity-boundary.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const P = await import(`${W}parentage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('zygosity-boundary: OM_TRIOS is not set to a readable directory. NOT RUN.')
  console.log('  This measurement needs the public GSE148488 arrays. It is skipped rather than')
  console.log('  quietly passing, because a skipped calibration that prints nothing reads exactly')
  console.log('  like one that succeeded.')
  process.exit(0)
}

/** Excluded as a parental genotype: heterozygosity 0.3455 at 82% call rate, a mixed sample. */
const CONTAMINATED = new Set(['GSM4472408'])

/** Material that went through whole-genome amplification, which is what the boundary must hold
 *  on. A cultured line and bulk genomic DNA did not and are not the failing case. */
const AMPLIFIED = new Set(['blastomere', 'trophectoderm', 'esc-single', 'zygote'])

type AB = 'AA' | 'AB' | 'BB' | 'NC'

const rows = readFileSync(`${DIR}/trio_manifest_full.csv`, 'utf8').trim().split('\n')
const head = rows[0].split(',')
const recs = rows.slice(1).map((l) => {
  const f = l.split(',')
  return Object.fromEntries(head.map((k, i) => [k, f[i] ?? ''])) as Record<string, string>
})
const present = (g: string): boolean => existsSync(`${DIR}/${g}.probes.gz`)
const firstUsable = (csv: string): string | undefined =>
  csv.split(/[;| ]/).map((x) => x.trim()).filter((x) => x && !CONTAMINATED.has(x) && present(x))[0]

/** One array, tallied exactly as the shipped code tallies it. */
function measure(gsm: string, parentGt?: Map<string, AB>): {
  band: number; homBand: number; callRate: number; zygosity: string
} {
  const text = gunzipSync(readFileSync(`${DIR}/${gsm}.probes.gz`)).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(lines[h])
  const t = P.emptyTally()
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (r) P.tallyRow(parentGt?.get(r.probesetId) ?? 'NC', r, t)
  }
  return {
    band: t.bafTotal ? t.bafInBand / t.bafTotal : NaN,
    homBand: t.homTotal ? t.homInBand / t.homTotal : NaN,
    callRate: t.markers ? t.called / t.markers : NaN,
    // The shipped answer, gates included, so the table shows what the tool actually says.
    zygosity: P.classify(t, 0.17).zygosity,
  }
}

function genotypes(gsm: string): Map<string, AB> {
  const text = gunzipSync(readFileSync(`${DIR}/${gsm}.probes.gz`)).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(lines[h])
  const gt = new Map<string, AB>()
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (r && r.genotype !== 'NC') gt.set(r.probesetId, r.genotype as AB)
  }
  return gt
}

// ---------------------------------------------------------------- the two truth sets
const N = Number(process.env.OM_N ?? 25)
const allChildren = recs
  .filter((r) => r.complete_trio.toLowerCase() === 'true' && AMPLIFIED.has(r.material)
    && present(r.gsm) && firstUsable(r.father_gsms) && firstUsable(r.mother_gsms))
  .sort((a, b) => a.gsm.localeCompare(b.gsm))
/** The headline truth set: the first N by accession, so the figure quoted in parentage.ts is
 *  reproducible without depending on how many arrays the directory happens to hold. */
const children = allChildren.slice(0, N)
const pronuclei = recs.filter((r) => r.pronucleus && present(r.gsm))
  .sort((a, b) => a.gsm.localeCompare(b.gsm))

console.log(`two-contribution set ${children.length} amplified children of complete trios`)
console.log(`one-contribution set ${pronuclei.length} pronuclei`)
console.log('\ngsm         truth        material        call    band    B_hom   excess  shipped')

interface Row {
  gsm: string; two: boolean; material: string; call: number
  band: number; homBand: number; excess: number; shipped: string
}
const table: Row[] = []
const seen = new Set<string>()
for (const [set, two] of [[children, true], [pronuclei, false], [allChildren, true]] as
     [typeof recs, boolean][]) {
  for (const r of set) {
    if (seen.has(r.gsm)) continue
    seen.add(r.gsm)
    const m = measure(r.gsm)
    const excess = m.band - m.homBand
    table.push({ gsm: r.gsm, two, material: r.material, call: m.callRate,
      band: m.band, homBand: m.homBand, excess, shipped: m.zygosity })
    console.log(`${r.gsm}  ${two ? 'two   ' : 'one   '}      ${r.material.padEnd(14)}  `
      + `${m.callRate.toFixed(3)}  ${m.band.toFixed(4)}  ${m.homBand.toFixed(4)}  `
      + `${excess.toFixed(4)}  ${m.zygosity}`)
  }
}

// ---------------------------------------------------------------- the two boundaries
const headline = table.filter((r) => !r.two || children.some((c) => c.gsm === r.gsm))
/**
 * Both arms of the boundary, on both truth sets.
 *
 * GATED is what the tool actually answers: below CALL_RATE_FLOOR zygosity is withheld, so a
 * refusal counts against sensitivity and is not a wrong answer on the other side. RAW is the
 * statistic alone, which is the only way to see what the boundary itself does rather than what
 * the gate in front of it does.
 */
const report = (name: string, diploid: (r: Row) => boolean, set: Row[], gate: boolean): void => {
  const call = (r: Row): boolean => (gate && r.call < P.CALL_RATE_FLOOR ? false : diploid(r))
  const two = set.filter((r) => r.two)
  const one = set.filter((r) => !r.two)
  const hit = two.filter(call).length
  const wrong = one.filter(call).length
  console.log(`  ${name.padEnd(22)} ${gate ? 'gated' : 'raw  '}  ${String(two.length).padStart(3)}`
    + ` two / ${one.length} one   sens ${(hit / two.length).toFixed(3)} (${hit}/${two.length})`
    + `   spec ${((one.length - wrong) / one.length).toFixed(3)} (${one.length - wrong}/`
    + `${one.length})   one-parent genomes called two-parent ${wrong}`)
}

for (const [label, set] of [['HEADLINE SET', headline], ['EVERY AMPLIFIED CHILD', table]] as
     [string, Row[]][]) {
  console.log(`\n=== BOUNDARIES, ${label}`)
  for (const gate of [true, false]) {
    report(`flat > ${P.HET_BAND_DIPLOID}`, (r) => r.band > P.HET_BAND_DIPLOID, set, gate)
    report(`band - B_hom > ${P.HET_BAND_EXCESS}`, (r) => r.excess > P.HET_BAND_EXCESS, set, gate)
  }
}
console.log('\n=== THE BOUNDARY SWEEP, headline set, gated')
for (const th of [0.02, 0.03, 0.04, 0.045, 0.05, 0.055, 0.06, 0.08]) {
  report(`band - B_hom > ${th}`, (r) => r.excess > th, headline, true)
}

const sorted = (f: (r: Row) => boolean): string =>
  headline.filter(f).map((r) => r.excess).sort((a, b) => a - b).map((x) => x.toFixed(4)).join(' ')
console.log(`\nexcess, one-parent genomes: ${sorted((r) => !r.two)}`)
console.log(`excess, two-parent genomes: ${sorted((r) => r.two)}`)

// ---------------------------------------------------------------- what changes, by material
//
// The number that decides whether a single-cell workflow works at all. A genome read as carrying
// ONE parental contribution switches off the Mendelian origin channel and both taxonomy detectors
// before any floor is consulted, so every child below that is read that way is a refusal.
console.log('\n=== EVERY AMPLIFIED CHILD OF A COMPLETE TRIO, all known to carry two contributions')
const byMaterial = new Map<string, Row[]>()
for (const r of table.filter((x) => x.two)) {
  const a = byMaterial.get(r.material) ?? []
  a.push(r)
  byMaterial.set(r.material, a)
}
console.log('material          n   read one-parent: flat / band-B_hom   withheld   changed to two')
for (const [mat, rs] of [...byMaterial].sort()) {
  const live = rs.filter((r) => r.call >= P.CALL_RATE_FLOOR)
  const flat = live.filter((r) => !(r.band > P.HET_BAND_DIPLOID)).length
  const exc = live.filter((r) => !(r.excess > P.HET_BAND_EXCESS)).length
  const gated = rs.length - live.length
  const flip = live.filter((r) => !(r.band > P.HET_BAND_DIPLOID)
    && r.excess > P.HET_BAND_EXCESS).length
  console.log(`${mat.padEnd(15)} ${String(rs.length).padStart(3)}   ${String(flat).padStart(14)}`
    + ` / ${String(exc).padStart(9)}   ${String(gated).padStart(8)}   ${String(flip).padStart(15)}`)
}

// ---------------------------------------------------------------- the wrong-parent arm
//
// Zygosity is read from the sample alone. If loading an unrelated adult moved it, the statistic
// would be measuring the comparison rather than the genome, and every number above would be about
// the harness. Three children, each run against its own father and against a donor who is not a
// parent of it.
console.log('\n=== WRONG-PARENT ARM (band and B_hom must not move)')
const donors = [...new Set(recs.filter((r) => r.role === 'parent' && present(r.gsm)
  && !CONTAMINATED.has(r.gsm)).map((r) => r.gsm))]
for (const c of children.slice(0, 3)) {
  const father = firstUsable(c.father_gsms)!
  const stranger = donors.find((g) => g !== father && g !== firstUsable(c.mother_gsms))!
  const t = measure(c.gsm, genotypes(father))
  const w = measure(c.gsm, genotypes(stranger))
  const same = t.band === w.band && t.homBand === w.homBand && t.zygosity === w.zygosity
  console.log(`${c.gsm}  true parent band ${t.band.toFixed(4)} B_hom ${t.homBand.toFixed(4)} `
    + `${t.zygosity}   wrong parent band ${w.band.toFixed(4)} B_hom ${w.homBand.toFixed(4)} `
    + `${w.zygosity}   ${same ? 'IDENTICAL' : 'MOVED - the statistic reads the parent'}`)
}

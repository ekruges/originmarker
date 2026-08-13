/**
 * Parental linkage, and the parental direction of every region, THROUGH THE SHIPPED CODE.
 *
 * WHY THIS FILE EXISTS RATHER THAN A SCRIPT. Every artefact recorded in AUDIT-4.1 came from
 * re-deriving this analysis outside the app, where none of the shipped guards exist: regions taken
 * from heterozygosity instead of intensity, self-derived nulls instead of external ones, a marker
 * floor picked rather than measured. Six one-directional results were produced that way and all
 * six were wrong. This imports the same modules the app runs, so a guard cannot be omitted by
 * forgetting it.
 *
 *   usage: node --experimental-strip-types audit/linkage.ts <array-dir> [out.json]
 *
 * WHAT IT DOES, IN THE ORDER THE ANSWERS DEPEND ON EACH OTHER.
 *
 *   1  Linkage. Every array is scored against each available BULK reference by
 *      opposite-homozygote rate, and a parent is assigned only where that rate is under the
 *      platform's own parent-child bound. Nothing downstream runs on an assumed relationship: on
 *      this platform unrelated adults sit near 0.07 rather than the 0.30 of a dense panel, so a
 *      threshold carried from elsewhere admits strangers, and did.
 *   2  Regions, from the INTENSITY channel only, via scanCopyNumber and the per-chromosome
 *      external null. Heterozygosity is never used to find a region. Amplification dropout removes
 *      genotype calls without removing DNA, so a het-derived region is a dropout report.
 *   3  Direction, via callGainOrigin and callLossOrigin, which read the same allele share in
 *      OPPOSITE directions and are therefore never interchangeable.
 *
 * Only bulk references are accepted. Measured against a single-cell reference the obligate-het
 * statistic does not separate ploidy at all, and a reconstructed mother inverts the very ratio
 * this analysis reports.
 */
import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
const { headerMap, parseRow } = await import(`${W}ingest.ts`)
const { isAutosome, pct } = await import(`${W}parentage.ts`)
const { scanCopyNumber, externalNull, segmentCoords } = await import(`${W}segments.ts`)
const {
  paternalShare, callGainOrigin, callLossOrigin, recentre, MIN_INFORMATIVE_DEFAULT,
} = await import(`${W}gainOrigin.ts`)

const DIR = process.argv[2]
const OUT = process.argv[3] ?? 'audit/linkage-results.json'
if (!DIR) throw new Error('usage: linkage.ts <array-dir> [out.json]')

/**
 * Parent-child bound for THIS platform, not carried from another.
 * Measured: parent-child 0.0002-0.02, unrelated adults 0.0727 on GPL28377, which is enriched for
 * low-frequency variants and so compresses every relatedness figure toward zero.
 */
const KIN_MAX = 0.02

type Row = { probesetId: string, chrom: string, pos: number, genotype: string, baf: number | null,
  log2R: number | null }

const readArray = (path: string): Row[] => {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const lines = raw.split('\n')
  const map = headerMap(lines[0])
  const out: Row[] = []
  for (let i = 1; i < lines.length; i += 1) {
    if (!lines[i]) continue
    const r = parseRow(lines[i], map)
    if (r) out.push(r as Row)
  }
  return out
}

const gtOf = (rows: Row[]): Map<string, string> => {
  const m = new Map<string, string>()
  for (const r of rows) m.set(r.probesetId, r.genotype)
  return m
}

/**
 * Opposite-homozygote rate. The one relatedness statistic that needs no allele frequencies.
 *
 * Screened on every STRIDE-th marker rather than all 825,657. Relatedness is a genome-wide
 * property and its standard error at 50,000 informative markers is under 0.002, which is two
 * orders of magnitude below the gap being resolved (parent-child under 0.02 against unrelated
 * 0.07 on this platform). At full density the screen is references x cells full-genome map
 * lookups, which does not finish; at this stride it does, and resolves the same question.
 */
const STRIDE = 16
const oppositeHom = (a: Map<string, string>, b: Map<string, string>): [number, number] => {
  let opp = 0
  let tot = 0
  let seen = 0
  for (const [k, ga] of a) {
    seen += 1
    if (seen % STRIDE) continue
    const gb = b.get(k)
    if (!gb) continue
    if ((ga === 'AA' || ga === 'BB') && (gb === 'AA' || gb === 'BB')) {
      tot += 1
      if (ga !== gb) opp += 1
    }
  }
  return [tot ? opp / tot : NaN, tot]
}

// Recursive: the laboratory exports sit in nested per-experiment folders, and a flat listing
// silently finds nothing rather than failing, which reads as "no data" instead of "wrong path".
const walk = (dir: string): string[] => readdirSync(dir, { withFileTypes: true })
  .flatMap((e) => (e.isDirectory() ? walk(`${dir}/${e.name}`) : [`${dir}/${e.name}`]))
const files = walk(DIR)
  .filter((f) => (f.endsWith('.probes') || f.endsWith('.txt') || f.endsWith('.gz')))
  .filter((f) => statSync(f).size > 1_000_000)
console.log(`${files.length} arrays in ${DIR}`)

// --- 1. references, and they must be bulk -------------------------------------------------------
//
// A reference is admitted on its own numbers rather than on its filename: bulk gDNA on this
// platform runs about 95% called and 15-17% heterozygous, and a single cell does not.
const refs: { name: string, gt: Map<string, string>, rows: Row[] }[] = []
for (const f of files) {
  const rows = readArray(f)
  if (!rows.length) continue
  const called = rows.filter((r) => r.genotype !== 'NC')
  const callRate = called.length / rows.length
  const het = called.filter((r) => r.genotype === 'AB').length / Math.max(called.length, 1)
  if (callRate >= 0.93 && het >= 0.12 && het <= 0.22) {
    refs.push({ name: f.split('/').pop()!.split('_')[0], gt: gtOf(rows), rows })
    console.log(`  reference candidate ${f.split('/').pop()}  call ${pct(callRate, 1)} `
      + `het ${pct(het, 1)}`)
  }
}
if (!refs.length) {
  console.log('\nNo bulk reference in this directory. Nothing downstream is computable, and a '
    + 'single-cell reference is not a substitute: measured, it does not separate ploidy at all.')
  writeFileSync(OUT, JSON.stringify({ refs: [], linked: [], regions: [] }, null, 2))
  process.exit(0)
}

// --- 2. linkage ---------------------------------------------------------------------------------
const linked: { cell: string, parent: string, rate: number, markers: number }[] = []
const cells: { name: string, rows: Row[], parent: string }[] = []
for (const f of files) {
  const name = f.split('/').pop()!.split('_')[0]
  if (refs.some((r) => r.name === name)) continue
  const rows = readArray(f)
  if (rows.length < 100_000) continue
  const g = gtOf(rows)
  for (const ref of refs) {
    const [rate, n] = oppositeHom(ref.gt, g)
    if (n > 3_000 && rate <= KIN_MAX) {
      linked.push({ cell: name, parent: ref.name, rate, markers: n })
      cells.push({ name, rows, parent: ref.name })
      break
    }
  }
}
console.log(`\n${linked.length} arrays linked to a bulk parent, of ${files.length - refs.length}`)
for (const l of linked.slice(0, 20)) {
  console.log(`  ${l.cell} -> ${l.parent}  ${l.rate.toFixed(4)} over ${l.markers} markers`)
}

// --- 3. regions from the intensity channel, then direction ---------------------------------------
const out: unknown[] = []
// A second reference, when one exists, is what makes a parental direction reportable at all.
const mother = refs.find((r) => r.name !== refs[0].name) ?? null
for (const c of cells) {
  const ref = refs.find((r) => r.name === c.parent)!
  const cn = new Map<string, { chrom: string, pos: number, called: boolean, log2R: number | null }[]>()
  const noCall = new Map<string, [number, number]>()
  const dosage = new Map<string, { chrom: string, pos: number, patShare: number }[]>()
  for (const r of c.rows) {
    if (!isAutosome(r.chrom)) continue
    const arr = cn.get(r.chrom) ?? []
    arr.push({ chrom: r.chrom, pos: r.pos, called: r.genotype !== 'NC', log2R: r.log2R })
    cn.set(r.chrom, arr)
    // Allele dosage is collected ONLY where BOTH parents are known and homozygous for opposite
    // alleles, which is what paternalShare requires and what makes the two directions symmetric.
    // With one parent it is tempting to substitute the opposite homozygote for the missing one;
    // that asserts something false at every marker and is how a one-directional artefact is built.
    // Measured on the one-parent path, the headroom is +0.019 above the sample's centre against
    // -0.090 below, so one direction cannot reach detection and any split is a property of the
    // geometry. Regions in one-parent cells are therefore reported WITHOUT a parental label.
    const fa = ref.gt.get(r.probesetId)
    const mo = mother ? mother.gt.get(r.probesetId) : undefined
    if (fa && mo && r.genotype !== 'NC') {
      const share = paternalShare(fa as 'AA' | 'BB' | 'AB' | 'NC',
        mo as 'AA' | 'BB' | 'AB' | 'NC', r.baf)
      if (share !== null) {
        const d = dosage.get(r.chrom) ?? []
        d.push({ chrom: r.chrom, pos: r.pos, patShare: share })
        dosage.set(r.chrom, d)
      }
    }
  }
  for (const [chrom, ms] of cn) noCall.set(chrom, [ms.length, ms.filter((m) => !m.called).length])
  const all = [...cn.values()].flat().map((m) => m.log2R)
    .filter((x): x is number => x !== null).sort((a, b) => a - b)
  const genomeLrr = all.length ? all[all.length >> 1] : 0
  const segs = [...cn].flatMap(([chrom, ms]) =>
    scanCopyNumber(ms, externalNull(noCall, chrom), genomeLrr))

  const allShare = [...dosage.values()].flat()
  const centre = allShare.length ? recentre(allShare) : NaN
  for (const s of segs) {
    const co = segmentCoords(s)
    const inside = (dosage.get(s.chrom) ?? []).filter((d) => d.pos >= co.start && d.pos <= co.end)
    const call = !mother
      ? { origin: 'unclear',
          why: 'one parent loaded; a direction from a single reference is a property of the '
            + 'measurement geometry rather than of the cell, so none is reported' }
      : s.kind === 'copy-loss'
        ? callLossOrigin(inside, centre, MIN_INFORMATIVE_DEFAULT)
        : callGainOrigin(inside, centre, MIN_INFORMATIVE_DEFAULT)
    out.push({
      cell: c.name, parent: c.parent, chrom: s.chrom, startBp: co.start, endBp: co.end,
      kind: s.kind, informative: inside.length, origin: call.origin, why: call.why,
    })
  }
}
const called = out.filter((r) => (r as { origin: string }).origin !== 'unclear')
console.log(`\n${out.length} regions in linked cells, ${called.length} carry a parental direction`)
writeFileSync(OUT, JSON.stringify({ refs: refs.map((r) => r.name), linked, regions: out }, null, 2))
console.log(`wrote ${OUT}`)

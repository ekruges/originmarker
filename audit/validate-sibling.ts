/**
 * Does the sibling-referenced call actually work? Measured against an INDEPENDENT channel.
 *
 * The two channels share no arithmetic. The intensity channel finds a region by log2R against a
 * per-chromosome external null and never looks at a genotype. The sibling channel reads genotypes
 * at markers the embryo's other cells established as heterozygous and never looks at intensity.
 * So their agreement is evidence rather than a restatement.
 *
 * POSITIVE CASES are regions the intensity channel called a copy-loss in a cell that has sibling
 * cells in the same embryo. The sibling call should report a missing copy.
 * NEGATIVE CASES are regions of identical marker count drawn on chromosomes of the SAME CELL where
 * the intensity channel found nothing. The sibling call should report no deletion. Matching on
 * marker count matters: a longer region is easier, so unmatched negatives would flatter the result.
 *
 * usage: node --experimental-strip-types audit/validate-sibling.ts <data-dir> [out.json]
 */
import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs'

const W = new URL('../web/src/', import.meta.url).pathname
const { headerMap, parseRow } = await import(`${W}ingest.ts`)
const { isAutosome } = await import(`${W}parentage.ts`)
const { scanCopyNumber, externalNull, segmentCoords } = await import(`${W}segments.ts`)
const { callSiblingOrigin, hetRule } = await import(`${W}siblingOrigin.ts`)

const DIR = process.argv[2]
const OUT = process.argv[3] ?? 'audit/validate-sibling.json'
if (!DIR) throw new Error('usage: validate-sibling.ts <data-dir> [out.json]')

const walk = (d: string): string[] => readdirSync(d, { withFileTypes: true })
  .flatMap((e) => (e.isDirectory() ? walk(`${d}/${e.name}`) : [`${d}/${e.name}`]))
const files = walk(DIR).filter((f) => f.endsWith('.probes') && statSync(f).size > 1_000_000)

type Row = { probesetId: string, chrom: string, pos: number, genotype: string, log2R: number | null }
const read = (p: string): Row[] => {
  const lines = readFileSync(p, 'utf8').split('\n')
  const map = headerMap(lines[0])
  const out: Row[] = []
  for (let i = 1; i < lines.length; i += 1) {
    if (!lines[i]) continue
    const r = parseRow(lines[i], map)
    if (r && isAutosome(r.chrom)) out.push(r as Row)
  }
  return out
}

/**
 * Embryo identity from the filename, and it must be the EMBRYO rather than the experiment.
 *
 * The construction needs cells of ONE embryo, which share a single genome. Cells of different
 * embryos from the same couple are full siblings, which is not the same thing at all: they are
 * independent draws from the parental haplotypes, so their heterozygous sites differ and using one
 * to establish het for another would be wrong rather than merely noisy.
 *
 * Only one naming scheme here encodes the embryo. Files like 01_chr16pcq_Z10__49 carry it as a
 * letter-plus-digits token, so several biopsies of embryo Z10 group correctly. Files like
 * 13605-23_51 carry an EXPERIMENT prefix and a biopsy index, and different indices are different
 * embryos, so no grouping is possible from the name and those arrays are excluded rather than
 * grouped wrongly. Returning null is how that exclusion is expressed.
 */
const embryoOf = (p: string): string | null => {
  const n = p.split('/').pop()!.replace(/\.CEL\.probes$/i, '')
  const m = n.match(/_([A-Za-z]\d+)_/)
  return m ? m[1].toUpperCase() : null
}

const groups = new Map<string, string[]>()
let unnamed = 0
for (const f of files) {
  const e = embryoOf(f)
  if (!e) { unnamed += 1; continue }
  groups.set(e, [...(groups.get(e) ?? []), f])
}
if (unnamed) {
  console.log(`${unnamed} arrays excluded: the filename does not identify an embryo, and cells of `
    + 'different embryos are not siblings for this purpose')
}
const usable = [...groups].filter(([, fs]) => fs.length >= 3).slice(0, 12)
console.log(`${files.length} arrays, ${groups.size} embryo groups, ${usable.length} with >=3 cells`)

let tp = 0; let fn = 0; let tn = 0; let fp = 0; let refusedPos = 0; let refusedNeg = 0
const rows: unknown[] = []

for (const [embryo, fs] of usable) {
  const loaded = fs.slice(0, 8).map((f) => ({ f, rows: read(f) })).filter((x) => x.rows.length > 100_000)
  if (loaded.length < 3) continue
  for (const me of loaded) {
    const sibs = loaded.filter((x) => x.f !== me.f)
    const need = hetRule(sibs.length)
    // het established from SIBLINGS only
    const sibHet = new Set<string>()
    const counts = new Map<string, number>()
    for (const s of sibs) {
      for (const r of s.rows) {
        if (r.genotype === 'AB') {
          const c = (counts.get(r.probesetId) ?? 0) + 1
          counts.set(r.probesetId, c)
          if (c >= need) sibHet.add(r.probesetId)
        }
      }
    }
    // intensity channel, independent of every genotype above
    const cn = new Map<string, { chrom: string, pos: number, called: boolean, log2R: number | null }[]>()
    for (const r of me.rows) {
      const a = cn.get(r.chrom) ?? []
      a.push({ chrom: r.chrom, pos: r.pos, called: r.genotype !== 'NC', log2R: r.log2R })
      cn.set(r.chrom, a)
    }
    const noCall = new Map<string, [number, number]>()
    for (const [c, ms] of cn) noCall.set(c, [ms.length, ms.filter((m) => !m.called).length])
    const lrr = me.rows.map((r) => r.log2R).filter((x): x is number => x !== null).sort((a, b) => a - b)
    const genomeLrr = lrr.length ? lrr[lrr.length >> 1] : 0
    const segs = [...cn].flatMap(([c, ms]) => scanCopyNumber(ms, externalNull(noCall, c), genomeLrr))
      .filter((sg) => sg.kind === 'copy-loss')
    const hitChroms = new Set(segs.map((s) => s.chrom))

    const gt = new Map<string, string>()
    const pos = new Map<string, { chrom: string, pos: number }>()
    for (const r of me.rows) { gt.set(r.probesetId, r.genotype); pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos }) }

    const score = (chrom: string, a: number, b: number): ReturnType<typeof callSiblingOrigin> => {
      const obs: string[] = []
      for (const p of sibHet) {
        const q = pos.get(p)
        if (!q || q.chrom !== chrom || q.pos < a || q.pos > b) continue
        const g = gt.get(p)
        if (g) obs.push(g)
      }
      return callSiblingOrigin(obs as never, sibs.length, noCall.get(chrom)![1] / noCall.get(chrom)![0])
    }

    for (const sg of segs) {
      const co = segmentCoords(sg)
      const c = score(sg.chrom, co.start, co.end)
      const missing = c.why.includes('a copy is missing')
      if (c.hypothesis === 'refused' && missing) tp += 1
      else if (c.hypothesis === 'no-deletion') fn += 1
      else refusedPos += 1
      rows.push({ embryo, cell: me.f.split('/').pop(), kind: 'positive', chrom: sg.chrom,
        startBp: co.start, endBp: co.end, markers: c.markers, verdict: missing ? 'copy missing' : c.hypothesis })

      // matched negative: same marker count, a chromosome with no called segment
      const free = [...cn.keys()].filter((k) => !hitChroms.has(k))
      if (!free.length) continue
      const ch = free[(sg.chrom.length + free.length) % free.length]
      const ps = [...sibHet].map((p) => pos.get(p)).filter((q) => q && q.chrom === ch)
        .map((q) => q!.pos).sort((x, y) => x - y)
      if (ps.length < c.markers + 1 || c.markers < 20) continue
      const i = Math.floor((ps.length - c.markers) / 2)
      const n = score(ch, ps[i], ps[i + c.markers - 1])
      const nMissing = n.why.includes('a copy is missing')
      if (n.hypothesis === 'no-deletion') tn += 1
      else if (nMissing) fp += 1
      else refusedNeg += 1
      rows.push({ embryo, cell: me.f.split('/').pop(), kind: 'negative', chrom: ch,
        startBp: ps[i], endBp: ps[i + c.markers - 1], markers: n.markers,
        verdict: nMissing ? 'copy missing' : n.hypothesis })
    }
  }
}
const sens = tp / Math.max(tp + fn, 1)
const spec = tn / Math.max(tn + fp, 1)
console.log(`\nPOSITIVES (intensity called a copy-loss)`)
console.log(`  sibling agrees a copy is missing : ${tp}`)
console.log(`  sibling says no deletion         : ${fn}`)
console.log(`  refused                          : ${refusedPos}`)
console.log(`NEGATIVES (matched on marker count, chromosomes with no called segment)`)
console.log(`  sibling says no deletion         : ${tn}`)
console.log(`  sibling says a copy is missing   : ${fp}`)
console.log(`  refused                          : ${refusedNeg}`)
console.log(`\nsensitivity ${sens.toFixed(4)}   specificity ${spec.toFixed(4)}`)
writeFileSync(OUT, JSON.stringify({ tp, fn, tn, fp, refusedPos, refusedNeg, sens, spec, rows }, null, 2))
console.log(`wrote ${OUT}`)

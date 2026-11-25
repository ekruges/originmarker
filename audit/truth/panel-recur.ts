/**
 * Do bulk ES lines of the same two parents read homozygous at the same opposite-homozygote sites?
 * If wrong calls in the parents' arrays make a bulk child read homozygous there, the sites recur
 * between sibling lines; independent failure would not make them recur. See
 * joy-tables/dev/MARKER-PANEL-plan-item.md, section 8, for the prediction and its bar.
 *
 * Counted on the same terms as audit/truth/panel-reads.ts: parents chosen the same way, lines the
 * tool reads as carrying both parents, chromosomes with no change on any line's result card.
 * Parent pairs and lines are named by index only.
 *
 *   OM_TRIOS=<dir> node --experimental-strip-types audit/truth/panel-recur.ts <records.jsonl>
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const stage = await import(`${W}stage.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('panel-recur: OM_TRIOS is not set to a readable directory. NOT RUN.')
  process.exit(0)
}

type AB = 'AA' | 'AB' | 'BB' | 'NC'
function load(gsm: string): { gt: Map<string, AB>; chrom: Map<string, string>; het: number } {
  const lines = gunzipSync(readFileSync(`${DIR}/${gsm}.probes.gz`)).toString('utf8').split('\n')
  let h = 0
  while (lines[h].startsWith('#')) h += 1
  const map = ingest.headerMap(lines[h])
  const gt = new Map<string, AB>(); const chrom = new Map<string, string>()
  let het = 0; let called = 0
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (!r || !/^\d+$/.test(r.chrom)) continue
    gt.set(r.probesetId, r.genotype as AB); chrom.set(r.probesetId, r.chrom)
    if (r.genotype !== 'NC') { called += 1; if (r.genotype === 'AB') het += 1 }
  }
  return { gt, chrom, het: called ? het / called : NaN }
}

const rows = readFileSync(`${DIR}/trio_manifest_full.csv`, 'utf8').trim().split('\n')
const head = rows[0].split(',')
const recs = rows.slice(1).map((l) => {
  const f = l.split(',')
  return Object.fromEntries(head.map((k, i) => [k, f[i] ?? ''])) as Record<string, string>
})
const present = (g: string) => existsSync(`${DIR}/${g}.probes.gz`)
const parents = recs.filter((r) => r.role === 'parent' && present(r.gsm)).map((r) => r.gsm)
const hetOf = new Map(parents.map((g) => [g, load(g).het]))
const panelHet = Math.min(...[...hetOf.values()].filter((x) => Number.isFinite(x)))
const over = new Set(parents.filter((g) => stage.parentOverCeiling(hetOf.get(g) as number, panelHet)))
const first = (csv: string) =>
  csv.split(/[;| ]/).map((x) => x.trim()).filter((x) => x && !over.has(x) && present(x))[0]

const verdicts = new Map<string, { both: boolean; changed: string[] }>()
for (const l of readFileSync(process.argv[2], 'utf8').split('\n')) {
  if (!l.startsWith('{')) continue
  const d = JSON.parse(l)
  if (d.tag !== 'two') continue
  verdicts.set(d.id, {
    both: d.parentGenome?.verdict === 'carries' && d.parentGenomeOther?.verdict === 'carries',
    changed: (d.card?.counted ?? []).map((c: { chrom: string }) => String(c.chrom)),
  })
}

const groups = new Map<string, string[]>()
for (const r of recs) {
  if (r.material !== 'esc-line' || r.complete_trio.toLowerCase() !== 'true' || !present(r.gsm)) continue
  const fa = first(r.father_gsms); const mo = first(r.mother_gsms)
  if (!fa || !mo || !verdicts.get(r.gsm)?.both) continue
  const k = `${fa}\t${mo}`
  groups.set(k, [...(groups.get(k) ?? []), r.gsm])
}

/** P(at least 2 of the lines fail) when line i fails on its own with chance p[i]. */
function atLeastTwo(p: number[]): number {
  let none = 1; let one = 0
  for (const x of p) { one = one * (1 - x) + none * x; none *= 1 - x }
  return 1 - none - one
}

console.log('pair\tlines\tsites counted\thomozygous per line\thomozygous in 2 or more\texpected if independent'
  + '\tmean share of a line\'s homozygous sites homozygous in a sibling\tmean sibling homozygous rate\tbar\tlowest pairwise share\tholds')
let pi = 0
for (const [k, lines] of groups) {
  if (lines.length < 2) continue
  pi += 1
  const [faId, moId] = k.split('\t')
  const fa = load(faId); const mo = load(moId)
  const kids = lines.map(load)
  const changed = new Set(lines.flatMap((g) => verdicts.get(g)!.changed))
  const homs = kids.map(() => new Set<string>())
  let n = 0
  for (const [id, f] of fa.gt) {
    const m = mo.gt.get(id)
    if (!((f === 'AA' && m === 'BB') || (f === 'BB' && m === 'AA'))) continue
    if (changed.has(kids[0].chrom.get(id) ?? '')) continue
    const calls = kids.map((c) => c.gt.get(id))
    if (calls.some((c) => c === undefined || c === 'NC')) continue
    n += 1
    calls.forEach((c, i) => { if (c !== 'AB') homs[i].add(id) })
  }
  const all = new Map<string, number>()
  for (const s of homs) for (const id of s) all.set(id, (all.get(id) ?? 0) + 1)
  const multi = [...all.values()].filter((v) => v >= 2).length
  const rate = homs.map((s) => s.size / n)
  const expected = n * atLeastTwo(rate)
  const shares: number[] = []; const rates: number[] = []
  for (let a = 0; a < homs.length; a += 1) {
    for (let b = 0; b < homs.length; b += 1) {
      if (a === b || !homs[a].size) continue
      let both = 0
      for (const id of homs[a]) if (homs[b].has(id)) both += 1
      shares.push(both / homs[a].size); rates.push(rate[b])
    }
  }
  const mean = (x: number[]) => x.reduce((s, v) => s + v, 0) / x.length
  const share = mean(shares); const sib = mean(rates); const bar = (sib + 0.5) / 2
  console.log([pi, lines.length, n, homs.map((s) => s.size).join(' '), multi, expected.toFixed(1),
    share.toFixed(4), sib.toFixed(4), bar.toFixed(4), Math.min(...shares).toFixed(4), share > bar ? 'yes' : 'no'].join('\t'))
}
console.log(`parent pairs with two or more bulk lines: ${pi}`)

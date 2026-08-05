// The TypeScript half of the cross-implementation check. Emits every quantity the browser
// builder reports, on real public arrays, for `tools/parity_check.py` to compare against the
// Python. Any divergence between the two is a divergence a user would see, in whichever half
// they happen to be using.
//
// Run: node --experimental-strip-types tools/parity_dump.ts <geo-dir> [out.json]
//      python tools/parity_check.py <geo-dir> [out.json]
//
// <geo-dir> holds the GSE148488 .CEL.txt.gz files, downloadable from the accession.
import { createReadStream, readdirSync, writeFileSync } from 'node:fs'
import { createGunzip } from 'node:zlib'
import { createInterface } from 'node:readline'
import { headerMap, parseRow } from '../web/src/ingest.ts'
import {
  ProductSet, groupByParent, MIN_ASCERTAINMENT, MIN_PRODUCTS,
  SAME_PARENT_MAX, DIFFERENT_PARENT_MIN,
} from '../web/src/inferredReference.ts'

const G = process.argv[2]
if (!G) throw new Error('usage: node tools/parity_dump.ts <geo-dir> [out.json]')
const all = readdirSync(G)
const pick = (g: string) => all.find((f) => f.startsWith(g + '_') && f.endsWith('.CEL.txt.gz'))!
const PRODUCTS = ['GSM4774680','GSM4774681','GSM4774682','GSM4774683','GSM4774685']

const ps = new ProductSet()
for (const g of PRODUCTS) {
  const slot = ps.begin(g)
  const baf = { inBand: 0, total: 0 }
  let map: ReturnType<typeof headerMap> = null
  const gz = createReadStream(`${G}/${pick(g)}`).pipe(createGunzip())
  for await (const line of createInterface({ input: gz, crlfDelay: Infinity })) {
    if (!map) { map = headerMap(line); continue }
    const r = parseRow(line, map)
    if (r) ps.add(slot, r, baf)
  }
  ps.end(slot, baf)
}
const out: Record<string, unknown> = { products: ps.ids, band: ps.band, called: ps.called }
const { mMin, ratios } = ps.chooseM()
out.chosenM = mMin
out.ascertainment = Object.fromEntries(ratios)
out.byM = {}
for (let m = 2; m <= PRODUCTS.length; m += 1) {
  const r = ps.build(m)
  ;(out.byM as Record<number, unknown>)[m] = {
    markers: r.markers, meanM: r.meanM, hRetained: r.hRetained, contamination: r.contamination,
  }
}
out.loo = {}
for (const id of PRODUCTS) {
  const r = ps.build(Math.min(mMin, PRODUCTS.length - 1), [id])
  ;(out.loo as Record<string, unknown>)[id] = { n: r.nProducts, markers: r.markers,
    contamination: r.contamination }
}
const pairs: Record<string, number> = {}
for (let a = 0; a < ps.size; a += 1) for (let b = a + 1; b < ps.size; b += 1) {
  pairs[`${ps.ids[a]}|${ps.ids[b]}`] = ps.opposite(a, b).rate
}
out.pairs = pairs
out.constants = {
  minAscertainment: MIN_ASCERTAINMENT, minProducts: MIN_PRODUCTS,
  sameParentMax: SAME_PARENT_MAX, differentParentMin: DIFFERENT_PARENT_MIN,
}
out.groups = groupByParent(ps.size, (a, b) => ps.opposite(a, b).rate).map((g) => g.map((i) => ps.ids[i]))
writeFileSync(
  process.argv[3] ?? 'ts_parity.json', JSON.stringify(out, null, 1))
console.log('TS done: chosenM', mMin, 'groups', (out.groups as string[][]).map(g => g.length))

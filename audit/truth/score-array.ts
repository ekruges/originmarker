/**
 * The scoring half of the audit runner (run.ts): an array read, its parents read, and the shipped
 * scoreSample run over it, as the page scores it. Nothing here builds a record.
 *
 * Every statement of this file is the root of the scoring closure (audit/gate/closure.ts), so its
 * imports are named: a name used from a module brings in that declaration, not the module.
 */
import { readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import {
  headerMap, parseRow, emptyBuildSums, accumulateBuild, buildVerdict, accumulate, accumulateBaf, emptyBafSums,
  finishProfile,
} from '../../web/src/ingest.ts'
import { noParent, emptyParent, collectParentRow, finishParent, emptyCollected, collectRow, scoreSample } from '../../web/src/scoreSample.ts'
import { parsePanel, FOUNDER_PANEL_PATH, emptySameHom, addSameHom } from '../../web/src/answerKeys.ts'
import { parseInferredBanner, MARK_BYTES } from '../../web/src/inferredArray.ts'
import { emptyHomTally, addHomozygous } from '../../web/src/informationRatio.ts'
import { unreadableLines, missingRequired, recogniseRows } from '../../web/src/platform.ts'

export const text = (p: string) => { const raw = readFileSync(p); return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8') }
export function rows(path: string): any[] {
  const lines = text(path).split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = headerMap(lines[h])
  if (!map) throw new ChipRefused(path, { text: unreadableLines(lines, panel.freqB),
    causes: missingRequired(lines[h]) })
  const o: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = parseRow(lines[i], map); if (r) o.push(r) }
  checkChip(path, o, map)
  return o
}
export const panel = parsePanel(text(new URL(`../../web/public/${FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
export class ChipRefused extends Error {
  file: string; call: any
  constructor(file: string, call: any) { super(call.text); this.file = file; this.call = call }
}
/** The same recognition the page runs, from the rows already read. Throws on another chip. */
function checkChip(path: string, rs: any[], map: any) {
  const b = emptyBuildSums()
  for (const r of rs) accumulateBuild(r, b)
  const call = recogniseRows(rs, buildVerdict(b), map, panel.freqB)
  if (!call.proven) throw new ChipRefused(path, call)
}
const parents = new Map<string, any>()
export function parentIndex(path: string) {
  if (path === '-') return noParent()
  if (parents.has(path)) return parents.get(path)
  const banner = parseInferredBanner(text(path).slice(0, MARK_BYTES)) ?? undefined
  const acc = emptyParent()
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const r of rows(path)) {
    if (!first) first = r.probesetId
    accumulate(r, byChrom); accumulateBaf(r, bafSums); collectParentRow(r, acc)
  }
  const idx = finishParent(acc, finishProfile(path, byChrom, bafSums, first).build.build, banner)
  parents.set(path, idx)
  return idx
}

// The called autosomal genotypes of an array of the same run, as the page keeps them (acc.myGt): read
// by probe id, held as one byte per probe so a whole run's arrays fit in memory together. Iterable
// with a size, so parentage.agreement walks it as it walks the page's map.
const probeIndex = new Map<string, number>()
const probeChrom: string[] = []
const CODES = ['', 'AA', 'AB', 'BB']
interface Calls { size: number; get(id: string): string | undefined; [Symbol.iterator](): Iterator<[string, string]> }
const callsOf = new Map<string, Calls>()
export function calls(path: string): Calls {
  let m = callsOf.get(path)
  if (!m) {
    const rs = rows(path).filter((r) => /^\d+$/.test(r.chrom))
    for (const r of rs) {
      if (!probeIndex.has(r.probesetId)) { probeIndex.set(r.probesetId, probeIndex.size); probeChrom.push(r.chrom) }
    }
    const code = new Uint8Array(probeIndex.size)
    for (const r of rs) if (r.genotype !== 'NC') code[probeIndex.get(r.probesetId)!] = CODES.indexOf(r.genotype)
    const at = (i: number | undefined) => (i === undefined || i >= code.length ? undefined : CODES[code[i]] || undefined)
    let size = 0
    for (let i = 0; i < code.length; i++) if (at(i)) size++
    m = {
      size,
      get: (id: string) => at(probeIndex.get(id)),
      *[Symbol.iterator]() { for (const [id, i] of probeIndex) { const g = at(i); if (g) yield [id, g] as [string, string] } },
    }
    callsOf.set(path, m)
  }
  return m
}
// Where each probe sits, as the page's acc.markerPos answers it for the same-genome comparison.
export const chromOf = { get: (id: string) => { const i = probeIndex.get(id); return i === undefined ? undefined : { chrom: probeChrom[i] } } }

/** One array read and collected as the page collects it, short of scoring. */
export function collectArray(id: string, path: string, parentPath: string, otherPath: string | undefined) {
  const pat = parentIndex(parentPath)
  const mat = otherPath ? parentIndex(otherPath) : null
  const acc = emptyCollected(pat, mat, panel as never)
  const byChrom = new Map(); const bafSums = emptyBafSums(); const builds = emptyBuildSums(); let first = ''
  const sameHom = emptySameHom(); const homTally = emptyHomTally()
  for (const r of rows(path)) {
    if (!first) first = r.probesetId
    accumulate(r, byChrom); accumulateBaf(r, bafSums); collectRow(r, pat, mat, acc as never)
    accumulateBuild(r as never, builds)
    if (!/^\d+$/.test(r.chrom)) continue
    const p = panel.get(r.probesetId)
    if (p) { addSameHom(sameHom, r.genotype, p, p); addHomozygous(homTally, r.genotype, p) }
  }
  const profile = finishProfile(id, byChrom, bafSums, first, builds)
  return { pat, mat, acc, profile, sameHom, homTally }
}

/** One array scored as the page scores it. */
export async function scoreArray(id: string, path: string, parentPath: string, otherPath: string | undefined,
  role: string | undefined, sibPaths: string[]) {
  const { pat, mat, acc, profile, sameHom, homTally } = collectArray(id, path, parentPath, otherPath)
  const res = await scoreSample({ acc, profile, pat, mat, soloRole: role === 'maternal' ? 'maternal' : 'paternal',
    sibs: sibPaths.map(calls),
    sampleName: id, log: () => {} } as never) as any
  return { pat, mat, acc, profile, sameHom, homTally, res }
}

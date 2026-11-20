/**
 * The sister-cell gain reading (web/src/sisterGain.ts) on every open embryo group, and its null.
 *
 * For each cell, every interval another cell of the same embryo shows lost is read as sisterGain reads
 * it: those are the `sister` lines. The same intervals are then read in cells of OTHER embryos, where
 * nothing reciprocal can be: the `null` lines. Every reading is written, called or not, so a threshold
 * is scored offline.
 *
 * Inputs: run.ts records (what each cell shows lost, and its reference chromosomes), features.ts
 * records (the 5 Mb window medians), and a tab-separated file of array id and embryo.
 *
 * Run: node --experimental-strip-types audit/truth/sisters.ts <runs.jsonl> <tag> <features.jsonl.gz>
 *        <embryos.tsv> <out.jsonl>
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const S = await import(new URL('../../web/src/sisterGain.ts', import.meta.url).pathname)
const [runsFile, tag, featFile, embryoFile, out] = process.argv.slice(2)
const lines = (p: string) => (p.endsWith('.gz') ? gunzipSync(readFileSync(p)) : readFileSync(p)).toString('utf8').split('\n')

const AUTOSOMES = Array.from({ length: 22 }, (_, i) => String(i + 1))
interface Unit { id: string; embryo: string; windows: any[]; reference: string[]; losses: any[]; changed: Set<string> }

const windows = new Map<string, any[]>()
for (const l of lines(featFile)) {
  if (!l.startsWith('{')) continue
  const r = JSON.parse(l)
  if (r.error) continue
  windows.set(r.id, r.windows
    .filter((w: any[]) => w[1] >= S.MIN_WINDOW_MARKERS && w[4] !== null)
    .map((w: any[]) => { const [c, i] = String(w[0]).split(':'); return [c, Number(i), w[4]] }))
}
const embryo = new Map<string, string>()
for (const l of lines(embryoFile)) { const [id, e] = l.split('\t'); if (id && e) embryo.set(id, e) }

const units: Unit[] = []
for (const l of lines(runsFile)) {
  if (!l.startsWith('{')) continue
  const r = JSON.parse(l)
  if (r.error || r.tag !== tag || !windows.has(r.id) || !embryo.has(r.id)) continue
  const whole = r.chroms.filter((c: any) => c.aneuploidy && /^\d+$/.test(c.chrom))
  const segs = r.segments.filter((s: any) => /^\d+$/.test(s.chrom) && (s.kind === 'copy-loss' || s.kind === 'copy-gain'))
  const changed = new Set<string>([...whole.map((c: any) => c.chrom), ...segs.map((s: any) => s.chrom)])
  units.push({
    id: r.id, embryo: embryo.get(r.id) as string, windows: windows.get(r.id) as any[],
    reference: (r.copyReference ?? AUTOSOMES).filter((c: string) => !changed.has(c)),
    losses: [
      ...whole.filter((c: any) => c.aneuploidy === 'loss').map((c: any) => S.WHOLE(c.chrom)),
      ...segs.filter((s: any) => s.kind === 'copy-loss')
        .map((s: any) => ({ chrom: s.chrom, startBp: s.start, endBp: s.end, whole: false })),
    ],
    changed,
  })
}

const rows: string[] = []
for (const a of units) {
  for (const loss of a.losses) {
    for (const b of units) {
      if (b.id === a.id) continue
      const r = S.readInterval(b.windows, b.reference, loss.chrom, loss.startBp, loss.endBp)
      if (!r) continue
      rows.push(JSON.stringify({
        kind: b.embryo === a.embryo ? 'sister' : 'null', cell: b.id, from: a.id, chrom: loss.chrom,
        whole: loss.whole, startBp: loss.startBp, endBp: loss.whole ? null : loss.endBp,
        changed: b.changed.has(loss.chrom), shift: +r.shift.toFixed(4), z: +r.z.toFixed(3), windows: r.windows,
      }))
    }
  }
}
writeFileSync(out, rows.join('\n') + '\n')
console.log(`${units.length} cells in ${new Set(units.map((u) => u.embryo)).size} embryos; ${rows.length} readings`)

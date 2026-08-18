// Fit the isotonic recalibration maps from THIS implementation's own injections.
//
// The methods review fitted maps and did not deliver them, and its maps could not be reproduced
// here because the corpus carries no material labels. But its maps are not the only honest ones:
// a map fitted on OUR posterior, over OUR noise, keyed on the material THIS TOOL ASSIGNS AT
// RUNTIME, is self-consistent by construction. It calibrates the assignment the tool actually
// makes rather than one it would make with labels it does not have.
//
// LEAVE-ONE-ARRAY-OUT, which is the whole difference between demonstrating calibration and
// asserting it. Every accuracy quoted for a map comes from rows the map never saw.
//
// Run: node audit/calibration/fit_maps.ts > web/src/calibrationMaps.json
import { readFileSync } from 'node:fs'
const W = new URL('../../web/src/', import.meta.url).pathname
const { originPosterior, shiftMean } = await import(`${W}originPosterior.ts`)

const rows = readFileSync(new URL('./null_shifts.csv', import.meta.url), 'utf8')
  .trim().split('\n').slice(1)
const byArray = new Map<string, number[]>()
for (const l of rows) {
  const p = l.split(','); const v = +p[3]
  if (Number.isFinite(v)) { if (!byArray.has(p[0])) byArray.set(p[0], []); byArray.get(p[0])!.push(v) }
}
const sd = (xs: number[]) => {
  const m = xs.reduce((a, b) => a + b, 0) / xs.length
  return Math.sqrt(xs.reduce((a, b) => a + (b - m) * (b - m), 0) / (xs.length - 1))
}
const arrays = [...byArray].filter(([, v]) => v.length >= 15).map(([n, v]) => ({ n, s: sd(v), nulls: v }))
arrays.sort((a, b) => a.s - b.s)

// The tool assigns four materials. Arrays are split into four noise tiers as a STAND-IN, stated as
// such: the corpus has no material labels, so this maps the noise each material is characterised by
// rather than the material itself.
const per = Math.ceil(arrays.length / 4)
const tiers: Record<string, typeof arrays> = {
  bulk: arrays.slice(0, per),
  trophectoderm: arrays.slice(per, 2 * per),
  'esc-single': arrays.slice(2 * per, 3 * per),
  blastomere: arrays.slice(3 * per),
}

let seed = 20260818
const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff
const CLASSES = ['loss', 'gain', 'cnn-loh'] as const

/** Pool-adjacent-violators: the standard isotonic fit, monotone by construction. */
function pava(xs: { x: number; y: number }[]): { raw: number[]; calibrated: number[] } {
  const pts = xs.slice().sort((a, b) => a.x - b.x)
  const v: { sum: number; n: number; x: number }[] = []
  for (const p of pts) {
    v.push({ sum: p.y, n: 1, x: p.x })
    while (v.length > 1 && v[v.length - 2].sum / v[v.length - 2].n > v[v.length - 1].sum / v[v.length - 1].n) {
      const b = v.pop()!; const a = v.pop()!
      v.push({ sum: a.sum + b.sum, n: a.n + b.n, x: b.x })
    }
  }
  const raw: number[] = []; const calibrated: number[] = []
  let i = 0
  for (const blk of v) {
    const val = blk.sum / blk.n
    raw.push(pts[Math.min(pts.length - 1, i)].x); calibrated.push(val)
    i += blk.n
  }
  // Thin to at most 24 knots so the shipped file stays small and the curve stays readable.
  const step = Math.max(1, Math.ceil(raw.length / 24))
  const rk: number[] = []; const ck: number[] = []
  for (let k = 0; k < raw.length; k += step) { rk.push(raw[k]); ck.push(calibrated[k]) }
  if (rk[rk.length - 1] !== raw[raw.length - 1]) {
    rk.push(raw[raw.length - 1]); ck.push(calibrated[calibrated.length - 1])
  }
  return { raw: rk, calibrated: ck }
}

const out: Record<string, { marginal: { raw: number[]; calibrated: number[] } }> = {}
const report: string[] = []
for (const [material, group] of Object.entries(tiers)) {
  if (!group.length) continue
  // Rows tagged with the array that produced them, so a map can be fitted leaving each one out.
  const all: { conf: number; ok: boolean; arr: string }[] = []
  for (const a of group) {
    for (let i = 0; i < 6000; i += 1) {
      const cls = CLASSES[Math.floor(rnd() * 3)]
      const f = 0.02 + rnd() * 0.68
      const aff = rnd() < 0.5 ? 'loaded' : 'other'
      const obs = shiftMean(cls, f, aff) + a.nulls[Math.floor(rnd() * a.nulls.length)]
      const p = originPosterior({ shift: obs, shiftSd: a.s, material: material as never, markers: 800 })
      if (Number.isFinite(p.confidence)) all.push({ conf: p.confidence, ok: p.parent === aff, arr: a.n })
    }
  }
  // Binned reliability, then isotonic through the bins.
  const bin = (rs: typeof all) => {
    const B = 40
    const acc: { x: number; y: number }[] = []
    for (let b = 0; b < B; b += 1) {
      const lo = 0.5 + (0.5 * b) / B; const hi = 0.5 + (0.5 * (b + 1)) / B
      const inb = rs.filter((r) => r.conf >= lo && r.conf < hi)
      if (inb.length >= 40) {
        acc.push({ x: inb.reduce((s, r) => s + r.conf, 0) / inb.length,
          y: inb.filter((r) => r.ok).length / inb.length })
      }
    }
    return acc
  }
  out[material] = { marginal: pava(bin(all)) }

  // LEAVE-ONE-ARRAY-OUT: fit without an array, score that array's rows through it.
  let ece = 0; let n = 0
  for (const a of group) {
    const fit = pava(bin(all.filter((r) => r.arr !== a.n)))
    const held = all.filter((r) => r.arr === a.n)
    const apply = (x: number) => {
      const { raw, calibrated } = fit
      if (!raw.length) return x
      if (x <= raw[0]) return calibrated[0]
      if (x >= raw[raw.length - 1]) return calibrated[calibrated.length - 1]
      let i = 1; while (i < raw.length && raw[i] < x) i += 1
      const t = (x - raw[i - 1]) / (raw[i] - raw[i - 1])
      return calibrated[i - 1] + t * (calibrated[i] - calibrated[i - 1])
    }
    const B = 10
    for (let b = 0; b < B; b += 1) {
      const lo = 0.5 + (0.5 * b) / B; const hi = 0.5 + (0.5 * (b + 1)) / B
      const inb = held.filter((r) => apply(r.conf) >= lo && apply(r.conf) < hi)
      if (inb.length < 30) continue
      const stated = inb.reduce((s, r) => s + apply(r.conf), 0) / inb.length
      const actual = inb.filter((r) => r.ok).length / inb.length
      ece += inb.length * Math.abs(stated - actual); n += inb.length
    }
  }
  report.push(`  ${material.padEnd(15)} ${group.length} arrays, ${all.length} rows, `
    + `${out[material].marginal.raw.length} knots, leave-one-out ECE ${n ? (ece / n).toFixed(4) : 'n/a'}`)
}
process.stderr.write(`fitted maps\n${report.join('\n')}\n`)
process.stdout.write(`${JSON.stringify(out, null, 1)}\n`)

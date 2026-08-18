// Does the shipped posterior treat the two ARMS of the inference equally?
//
// "Call M by absence of P" is the workhorse for a one-parent cohort: the father is loaded, and a
// maternal-origin event is named because HIS contribution is intact while a change is present. The
// methods review measured that arm as materially weaker on amplified material, a maternal gain on
// trophectoderm called correctly 0.402 of the time against 0.727 for the paternal arm, and
// concluded the two must be calibrated separately or a pooled number is over-confident on exactly
// the maternal side.
//
// That was measured on THEIR pipeline. This measures the same split on ours, with real array noise,
// because a correction fitted to someone else's asymmetry is not a correction.
import { readFileSync } from 'node:fs'
const W = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/src/'
const { originPosterior, shiftMean, logRMean } = await import(`${W}originPosterior.ts`)
const HERE = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/audit/calibration/'

const rows = readFileSync(HERE + 'null_shifts.csv', 'utf8').trim().split('\n').slice(1)
const byArray = new Map<string, number[]>()
for (const l of rows) {
  const p = l.split(','); const v = +p[3]
  if (Number.isFinite(v)) { if (!byArray.has(p[0])) byArray.set(p[0], []); byArray.get(p[0])!.push(v) }
}
const arrays = [...byArray].filter(([, v]) => v.length >= 15)
const sd = (xs: number[]) => { const m = xs.reduce((a, b) => a + b, 0) / xs.length
  return Math.sqrt(xs.reduce((a, b) => a + (b - m) * (b - m), 0) / (xs.length - 1)) }
let seed = 4242
const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff
const CLASSES = ['loss', 'gain', 'cnn-loh'] as const

type R = { arm: 'loaded' | 'other'; cls: string; correct: boolean; conf: number; f: number }
const res: R[] = []
for (const [, nulls] of arrays) {
  const s = sd(nulls)
  for (let i = 0; i < 6000; i++) {
    const cls = CLASSES[Math.floor(rnd() * 3)]
    const f = 0.02 + rnd() * 0.68
    // The ARM: whose copy was actually affected. 'other' is the un-genotyped parent, which is the
    // maternal-by-absence case when the father is the one loaded.
    const arm = rnd() < 0.5 ? 'loaded' : 'other'
    const observed = shiftMean(cls, f, arm) + nulls[Math.floor(rnd() * nulls.length)]
    const withI = i % 2 === 0
    const logRSd = 0.02
    const logR = withI ? logRMean(cls, f) + (rnd() - 0.5) * logRSd * 2 : undefined
    const p = originPosterior({ shift: observed, shiftSd: s, logR, logRSd: withI ? logRSd : undefined,
      material: 'trophectoderm', markers: 800 })
    if (!Number.isFinite(p.confidence)) continue
    res.push({ arm, cls, correct: p.parent === arm, conf: p.confidence, f })
  }
}
const acc = (rs: R[]) => (rs.length ? rs.filter((r) => r.correct).length / rs.length : NaN)
const show = (label: string, f: (r: R) => boolean) => {
  const L = res.filter((r) => f(r) && r.arm === 'loaded')
  const O = res.filter((r) => f(r) && r.arm === 'other')
  if (L.length < 100 || O.length < 100) return
  const gap = acc(L) - acc(O)
  console.log(`  ${label.padEnd(30)} loaded ${acc(L).toFixed(4)}  other ${acc(O).toFixed(4)}  `
    + `gap ${gap >= 0 ? '+' : ''}${gap.toFixed(4)}`)
}
console.log(`${res.length} injections, split by which parent's copy was actually affected\n`)
console.log('  "other" is the un-genotyped parent: the call-M-by-absence-of-P arm\n')
show('all', () => true)
for (const c of CLASSES) show(`class ${c}`, (r) => r.cls === c)
console.log('')
for (const c of CLASSES) show(`class ${c}, low fraction`, (r) => r.cls === c && r.f <= 0.15)
console.log('')
show('top band only (conf >= 0.985)', (r) => r.conf >= 0.985)
show('band B (0.90-0.985)', (r) => r.conf >= 0.90 && r.conf < 0.985)

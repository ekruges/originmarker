/**
 * The cell-level split of the pair likelihood, scored exactly as
 * audit/gain-origin-split-preregistration.txt fixes it before the held-out half is seen.
 *
 * Reads the constructed builds written by gain-origin-truth.ts with OM_OUT. Splits the cells by
 * accession, even to tuning and odd to held-out, picks both thresholds on the tuning half by the
 * pre-registered rules, and applies them unchanged to the held-out half.
 *
 * Run: OM_WIDE=<dir of harness JSON> node --experimental-strip-types audit/gain-origin-split.ts
 */
import { readFileSync, readdirSync } from 'node:fs'
import { join } from 'node:path'

interface Row {
  gsm: string; chrom: string; extra: string; kind?: string; f?: number; verdict: string
  parent?: string; scale: number; shift?: number; lrr?: number
  pairState?: string; pairP?: number; pairTilt?: string; pairTiltP?: number
}

const DIR = process.env.OM_WIDE ?? '.'
const rows: Row[] = []
for (const f of readdirSync(DIR).filter((x) => x.endsWith('.json')).sort()) {
  rows.push(...(JSON.parse(readFileSync(join(DIR, f), 'utf8')).built ?? []))
}
// A build is one cell, one chromosome, one parent, one kind at one fraction, and a shard covers a
// disjoint set of cells. The same build twice means the directory holds more than one run, and
// averaging two code versions over the same arms would read as a larger sample than exists.
{
  const seen = new Set<string>()
  const dup = rows.filter((r) => {
    const k = `${r.gsm}|${r.chrom}|${r.extra}|${r.kind}|${r.f}`
    if (seen.has(k)) return true
    seen.add(k); return false
  })
  if (dup.length) {
    console.error(`${DIR} holds ${dup.length} duplicate builds of ${rows.length}: more than one run. `
      + `Point OM_WIDE at one run's output only. First: ${dup[0].gsm} chr${dup[0].chrom} `
      + `${dup[0].kind} ${dup[0].f} ${dup[0].extra}`)
    process.exit(1)
  }
}

/** Which parent the build leaves over-represented: a gain of one is a loss of the other. */
const truthTilt = (r: Row): string =>
  r.kind === 'gain' ? r.extra : (r.extra === 'paternal' ? 'maternal' : 'paternal')
const named = (r: Row): boolean => !!r.pairState && r.pairState.includes('-')
const dirOf = (s: string): string => s.split('-')[0]

// Harness output written before callPair returned the tilt carries only the state. The tilt is a
// pure function of the state, and the state's own posterior is a lower bound on the summed one, so
// an older file scores conservatively rather than being unreadable.
for (const r of rows) {
  if (!named(r)) continue
  if (!r.pairTilt) {
    const [k, p] = r.pairState!.split('-')
    r.pairTilt = k === 'gain' ? p : (p === 'paternal' ? 'maternal' : 'paternal')
  }
  if (r.pairTiltP == null) r.pairTiltP = r.pairP
}

const cells = [...new Set(rows.map((r) => r.gsm))].sort()
const tune = new Set(cells.filter((_, i) => i % 2 === 0))
const halfOf = (r: Row): 'tune' | 'held' => (tune.has(r.gsm) ? 'tune' : 'held')

const pct = (n: number, d: number): string => (d ? `${((100 * n) / d).toFixed(1)}%` : '-')
const fmt = (n: number, d: number): string => `${n} of ${d} (${pct(n, d)})`

// --- the two curves, on whichever half is asked for -------------------------------------------
function tiltCurve(set: Row[], t: number) {
  const hit = set.filter((r) => named(r) && (r.pairTiltP ?? 0) >= t)
  const wrong = hit.filter((r) => r.pairTilt !== truthTilt(r))
  return { named: hit.length, wrong: wrong.length, right: hit.length - wrong.length, of: set.length }
}
function dirCurve(set: Row[], t: number) {
  const hit = set.filter((r) => named(r) && (r.pairP ?? 0) >= t)
  const wrong = hit.filter((r) => dirOf(r.pairState!) !== r.kind)
  return { named: hit.length, wrong: wrong.length, right: hit.length - wrong.length, of: set.length }
}

const GRID = [0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9999]
const T = rows.filter((r) => halfOf(r) === 'tune')
const H = rows.filter((r) => halfOf(r) === 'held')

console.log(`${rows.length} constructed builds over ${cells.length} cells`)
console.log(`  tuning   ${tune.size} cells, ${T.length} builds`)
console.log(`  held out ${cells.length - tune.size} cells, ${H.length} builds`)

// Is the split comparing like with like? A secondary endpoint, read before the thresholds.
const med = (xs: number[]): number => xs.length ? xs.slice().sort((a, b) => a - b)[xs.length >> 1] : NaN
for (const [n, s] of [['tuning', T], ['held out', H]] as [string, Row[]][]) {
  console.log(`  ${n.padEnd(8)} median cell share spread ${med(s.map((r) => r.scale)).toFixed(4)}`)
}

console.log('\nTUNING HALF, TILT')
console.log(`  ${'t'.padStart(7)} ${'named'.padStart(16)} ${'wrong tilt'.padStart(14)}`)
for (const t of GRID) {
  const c = tiltCurve(T, t)
  console.log(`  ${t.toFixed(4).padStart(7)} ${fmt(c.named, c.of).padStart(16)} ${fmt(c.wrong, c.named).padStart(14)}`)
}

console.log('\nTUNING HALF, DIRECTION')
console.log(`  ${'t'.padStart(7)} ${'named'.padStart(16)} ${'wrong direction'.padStart(16)}`)
for (const t of GRID) {
  const c = dirCurve(T, t)
  console.log(`  ${t.toFixed(4).padStart(7)} ${fmt(c.named, c.of).padStart(16)} ${fmt(c.wrong, c.named).padStart(16)}`)
}

// --- the pre-registered rules, applied to the tuning half only ---------------------------------
// T_tilt: the lowest t with no wrong tilt; failing that, the lowest with at most 1 percent wrong.
// T_dir:  the lowest t with at most 5 percent of named rows in the wrong direction.
const pick = (
  curve: (s: Row[], t: number) => { named: number; wrong: number }, cap: number,
): { t: number; relaxed: boolean } => {
  // A threshold that names nothing has no wrong calls either, so it has to be excluded or it
  // would satisfy the rule by refusing everything.
  for (const t of GRID) { const c = curve(T, t); if (c.named && c.wrong === 0) return { t, relaxed: false } }
  for (const t of GRID) {
    const c = curve(T, t)
    if (c.named && c.wrong / c.named <= cap) return { t, relaxed: true }
  }
  return { t: Infinity, relaxed: true }
}
const tTilt = pick(tiltCurve, 0.01)
const tDir = pick(dirCurve, 0.05)
console.log(`\nT_tilt = ${tTilt.t}${tTilt.relaxed ? ' (no threshold reached zero; the 1 percent rule applied)' : ' (no wrong tilt on the tuning half)'}`)
console.log(`T_dir  = ${tDir.t}${tDir.relaxed ? ' (the 5 percent rule applied, as expected)' : ' (no wrong direction on the tuning half)'}`)

// --- and unchanged to the held-out half --------------------------------------------------------
console.log('\nHELD-OUT HALF AT THE THRESHOLDS PICKED ON THE TUNING HALF')
for (const [label, c] of [['tilt', tiltCurve(H, tTilt.t)], ['direction', dirCurve(H, tDir.t)]] as [string, ReturnType<typeof tiltCurve>][]) {
  console.log(`  ${label.padEnd(10)} named ${fmt(c.named, c.of).padEnd(22)} right ${fmt(c.right, c.named).padEnd(22)} wrong ${fmt(c.wrong, c.named)}`)
}

console.log('\nPRIMARY ENDPOINT: held-out tilt by built fraction')
console.log(`  ${'f'.padStart(5)} ${'builds'.padStart(7)} ${'tilt named'.padStart(18)} ${'of those right'.padStart(18)} ${'shipped detects'.padStart(17)}`)
for (const f of [...new Set(rows.map((r) => r.f))].filter((x): x is number => x != null).sort()) {
  for (const kind of ['gain', 'loss']) {
    const s = H.filter((r) => r.f === f && r.kind === kind)
    if (!s.length) continue
    const c = tiltCurve(s, tTilt.t)
    const shipped = s.filter((r) => r.verdict === 'gain' && r.parent === r.extra && r.kind === 'gain').length
      + s.filter((r) => r.verdict === 'gain' && r.kind === 'loss').length // the latter would be false
    console.log(`  ${String(f).padStart(5)} ${(kind + ' ' + s.length).padStart(7)} ${fmt(c.named, c.of).padStart(18)} `
      + `${fmt(c.right, c.named).padStart(18)} ${fmt(shipped, s.length).padStart(17)}`)
  }
}

console.log('\nHELD-OUT DIRECTION AT T_dir, BY FRACTION')
console.log(`  ${'f'.padStart(5)} ${'kind'.padStart(5)} ${'named'.padStart(16)} ${'wrong direction'.padStart(16)}`)
for (const f of [...new Set(rows.map((r) => r.f))].filter((x): x is number => x != null).sort()) {
  for (const kind of ['gain', 'loss']) {
    const s = H.filter((r) => r.f === f && r.kind === kind)
    if (!s.length) continue
    const c = dirCurve(s, tDir.t)
    console.log(`  ${String(f).padStart(5)} ${kind.padStart(5)} ${fmt(c.named, c.of).padStart(16)} ${fmt(c.wrong, c.named).padStart(16)}`)
  }
}

console.log('\nHELD-OUT CURVES, so the operating point can be read against what was available')
console.log(`  ${'t'.padStart(7)} ${'tilt named'.padStart(16)} ${'tilt wrong'.padStart(14)} ${'dir named'.padStart(16)} ${'dir wrong'.padStart(14)}`)
for (const t of GRID) {
  const a = tiltCurve(H, t); const b = dirCurve(H, t)
  console.log(`  ${t.toFixed(4).padStart(7)} ${fmt(a.named, a.of).padStart(16)} ${fmt(a.wrong, a.named).padStart(14)} `
    + `${fmt(b.named, b.of).padStart(16)} ${fmt(b.wrong, b.named).padStart(14)}`)
}

// --- the verdict the pre-registration fixed ----------------------------------------------------
const h = tiltCurve(H, tTilt.t); const d = dirCurve(H, tDir.t)
const shippedDetects = H.filter((r) => r.verdict === 'gain' && r.kind === 'gain' && r.parent === r.extra).length
console.log('\nAGAINST THE PRE-REGISTERED BAR')
console.log(`  tilt may decide a row:      ${h.named && h.wrong / h.named <= 0.01 && h.named > shippedDetects ? 'MET' : 'NOT MET'}`
  + ` (wrong ${pct(h.wrong, h.named)} against a 1% bar, named ${h.named} against ${shippedDetects} the thresholds detect)`)
console.log(`  direction may be a call:    ${d.named && d.wrong / d.named <= 0.05 ? 'MET' : 'NOT MET'}`
  + ` (wrong ${pct(d.wrong, d.named)} against a 5% bar)`)

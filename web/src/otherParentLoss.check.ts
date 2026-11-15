// Self-check for the one-parent scan for the other parent's copy. Run: node src/otherParentLoss.check.ts
import assert from 'node:assert/strict'
import {
  otherRates, scanOtherParent, coveredShare, OTHER_LRT, OTHER_MIN_SPAN_BP, type OtherMarker,
} from './otherParentLoss.ts'

/** A deterministic stream, so a failure is the code changing and never the draw changing. */
function rng(seed: number): () => number {
  let x = seed >>> 0
  return () => {
    x ^= x << 13; x >>>= 0
    x ^= x >> 17
    x ^= x << 5; x >>>= 0
    return x / 0x100000000
  }
}

const B = 0.03; const A = 0.7; const H = 0.3
const SPACING = 3_600

interface Stretch { from: number; to: number; kind: 'gone' | 'parent-gone' | 'opp' }

/**
 * One chromosome of `n` markers. `clusters` makes calls fail together: in a cluster no second allele
 * is seen, whatever is there.
 */
function chromosome(
  r: () => number, n: number, stretches: Stretch[] = [], clusters = 0, clusterMarkers = 300,
): OtherMarker[] {
  const out: OtherMarker[] = []
  let inCluster = 0
  for (let i = 0; i < n; i += 1) {
    const pos = i * SPACING
    if (clusters && inCluster === 0 && r() < clusters) inCluster = clusterMarkers
    const failed = inCluster > 0
    if (inCluster) inCluster -= 1
    const st = stretches.find((s) => pos >= s.from && pos < s.to)
    const het = r() < 0.16
    const q = het ? null : Math.max(0.001, r() ** 2 * 0.9)
    let seen: boolean; let opp = false
    if (st?.kind === 'parent-gone') {
      // P's own copy is the missing one: the other parent's allele alone, as the opposite homozygote.
      seen = q !== null && r() < A * q
      opp = seen
    } else if (st?.kind === 'opp') {
      seen = q !== null && r() < 0.05
      opp = seen
    } else {
      const p = st ? B : q === null ? H : B + A * q
      seen = !failed && r() < p
      opp = seen && q !== null && r() < 0.02
    }
    out.push({ pos, q, seen, opp })
  }
  return out
}

const r = rng(11)
const reference = Array.from({ length: 6 }, () => chromosome(r, 20_000))
const rates = otherRates(reference)!

// --- 1. the array's own rates come back, and independent calls are not clustered --------------
{
  assert.ok(Math.abs(rates.b - B) < 0.01, `b ${rates.b}`)
  assert.ok(Math.abs(rates.a - A) < 0.05, `a ${rates.a}`)
  assert.ok(Math.abs(rates.h - H) < 0.02, `h ${rates.h}`)
  assert.ok(rates.inflation < 1.5, `independent calls read clustered: ${rates.inflation}`)
  assert.equal(otherRates(reference.slice(0, 1)), null, 'one reference chromosome is not a reference')
}

// --- 2. a chromosome carrying both copies throughout reports nothing -------------------------
for (let k = 0; k < 5; k += 1) {
  assert.deepEqual(scanOtherParent('1', chromosome(r, 20_000), rates), [])
}

// --- 3. a stretch missing the other parent's copy is found, where it is -----------------------
{
  const segs = scanOtherParent('5', chromosome(r, 20_000, [{ from: 30e6, to: 36e6, kind: 'gone' }]), rates)
  assert.equal(segs.length, 1)
  const [sg] = segs
  assert.equal(sg.kind, 'parental-absence')
  assert.ok(sg.score >= OTHER_LRT)
  assert.ok(Math.abs(sg.startBp - 30e6) < 0.3e6 && Math.abs(sg.endBp - 36e6) < 0.3e6,
    `edges ${sg.startBp}-${sg.endBp}`)
  assert.ok(sg.rate > sg.nullRate, 'the rate reported is the one that rose')
}

// --- 4. under the span floor nothing is reported, however clear ------------------------------
{
  const short = OTHER_MIN_SPAN_BP - 1e6
  const c = chromosome(r, 20_000, [{ from: 30e6, to: 30e6 + short, kind: 'gone' }])
  assert.deepEqual(scanOtherParent('5', c, rates), [])
  // Asked for from 1 Mb it comes back, for the caller to withhold and name.
  const held = scanOtherParent('5', c, rates, 1e6)
  assert.equal(held.length, 1)
  assert.ok(held[0].spanBp < OTHER_MIN_SPAN_BP && held[0].spanBp > short - 0.6e6)
}

// --- 5. P's own copy missing is not this event ------------------------------------------------
{
  assert.deepEqual(scanOtherParent('5', chromosome(r, 20_000, [{ from: 30e6, to: 40e6, kind: 'parent-gone' }]), rates), [])
  // No second allele seen, yet one homozygous call in twenty is the homozygote P cannot give: the
  // copy that remains is not shown to be P's.
  assert.deepEqual(scanOtherParent('5', chromosome(r, 20_000, [{ from: 30e6, to: 40e6, kind: 'opp' }]), rates), [])
}

// --- 6. calls that fail in clusters raise the bar and report nothing --------------------------
{
  const rc = rng(5)
  const noisy = (): OtherMarker[] => chromosome(rc, 20_000, [], 0.001, 500)
  const clustered = otherRates(Array.from({ length: 6 }, noisy))!
  assert.ok(clustered.inflation > 2, `clustering not measured: ${clustered.inflation}`)
  // Read as independent calls the same chromosomes report runs, which is what the division is for.
  let divided = 0; let undivided = 0
  for (let k = 0; k < 10; k += 1) {
    const c = noisy()
    divided += scanOtherParent('1', c, clustered).length
    undivided += scanOtherParent('1', c, { ...clustered, inflation: 1 }).length
  }
  assert.equal(divided, 0)
  assert.ok(undivided > 0, 'the clustered construction must be one the undivided score falls for')
}

// --- 7. a real loss is still found on a clustered array ---------------------------------------
{
  const rc = rng(9)
  const clustered = otherRates(Array.from({ length: 6 }, () => chromosome(rc, 20_000, [], 0.0005)))!
  const segs = scanOtherParent('5', chromosome(rc, 20_000, [{ from: 20e6, to: 50e6, kind: 'gone' }], 0.0005), clustered)
  assert.equal(segs.length, 1)
  assert.ok(segs[0].startBp < 21e6 && segs[0].endBp > 49e6)
}

// --- 8. how much of a run other intervals cover: overlaps counted once ------------------------
{
  const iv = (startBp: number, endBp: number) => ({ startBp, endBp })
  assert.equal(coveredShare(0, 100, []), 0)
  assert.equal(coveredShare(0, 100, [iv(-50, 50)]), 0.5)
  assert.equal(coveredShare(0, 100, [iv(10, 40), iv(30, 60), iv(90, 200)]), 0.6)
  assert.equal(coveredShare(0, 100, [iv(200, 300)]), 0)
}

console.log('otherParentLoss ok')

// Self-check for the relationship statistics. Run: node src/relatedness.check.ts
//
// These moved here from the command-line file, where they were the only copy, so the browser had
// no way to notice that the array in the parental slot was not a parent. The arithmetic is
// unchanged and is exercised end to end by cli/om.check.ts; what is pinned here is the part that
// only became reachable once the scorer started calling it with a position map built from a
// DIFFERENT file than the genotypes.
import assert from 'node:assert/strict'
import {
  oppositeHom, windowedIbs0, relate, REL, OPPOSITE_HOM_MAX, IBD0_SPREAD_MIN, IBD0_FLOOR_MAX,
} from './relatedness.ts'
import type { AB } from './informativity.ts'

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

/** A parent and a transmitted child, at a given number of markers on two chromosomes. */
function trio(seed: number, n: number) {
  const r = rng(seed)
  const pg = new Map<string, AB>()
  const cg = new Map<string, AB>()
  const pos = new Map<string, { chrom: string; pos: number }>()
  for (let i = 0; i < n; i += 1) {
    const id = `p${i}`
    // A parent that is heterozygous everywhere gives no informative marker to this statistic, so
    // the draw matches the real panel: about 17 percent heterozygous.
    const p: AB = r() < 0.17 ? 'AB' : (r() < 0.5 ? 'AA' : 'BB')
    pg.set(id, p)
    // The child gets one allele from this parent and one from an unrelated other parent. It can
    // never be the OPPOSITE homozygote, which is the whole basis of the test.
    const fromP = p === 'AB' ? (r() < 0.5 ? 'A' : 'B') : p[0]
    const fromOther = r() < 0.5 ? 'A' : 'B'
    const c = fromP === fromOther ? (fromP + fromP) as AB : 'AB'
    cg.set(id, c)
    pos.set(id, { chrom: i < n / 2 ? '1' : '2', pos: i * 1000 })
  }
  return { pg, cg, pos }
}

// --- 1. A PROBE WITH NO POSITION IS SKIPPED, NOT SORTED -----------------------------------------
//
// THE CRASH THIS PINS WAS REAL AND REACHED AN OPERATOR AS A STACK TRACE. The scorer calls this with
// the parent's genotypes and the SAMPLE's position map, so any probe the parent carries and the
// sample does not has no position. Sorting on it threw "Cannot read properties of undefined" from
// inside Array.sort, on a pair of perfectly good files. Two arrays need not carry the same probe
// set: panel versions differ, and a lab running hundreds of samples will mix them.
{
  // 40,000 markers: at 17 percent parental heterozygosity, and a child homozygous about half the
  // time, that leaves roughly 16,000 informative markers and 80 windows of 200.
  const { pg, cg, pos } = trio(11, 40_000)
  // The parent carries 400 probes the position map has never heard of.
  for (let i = 0; i < 400; i += 1) pg.set(`orphan${i}`, i % 2 ? 'AA' : 'BB')
  const out = windowedIbs0({ gt: pg, pos }, cg, 200)
  assert.ok(Number.isFinite(out.spread),
    'a probe with no position must be skipped rather than crash the sort')
  assert.ok(out.windows > 0, 'and the windows that DO have positions must still be measured')
  // And it must not silently measure nothing: skipping everything would look like a pass here.
  assert.ok(out.windows >= 10,
    `the orphans must not take the real markers with them, got ${out.windows} windows`)
}

// --- 2. THE SYMMETRY THAT MAKES DIRECTION IMPOSSIBLE --------------------------------------------
//
// Opposite homozygotes are symmetric BY CONSTRUCTION, which is why no verdict here names a
// direction. Pinned so nobody reintroduces "child" as an answer.
{
  const { pg, cg } = trio(23, 20_000)
  const ab = oppositeHom(pg, cg)
  const ba = oppositeHom(cg, pg)
  assert.equal(ab.rate, ba.rate, 'swapping the two files must change nothing at all')
  assert.equal(ab.n, ba.n)
  // The verdict that DOES assert a parent-child relationship has to say the direction is open. The
  // others may mention "child" freely: "not parent and child" is a denial, not a direction.
  assert.ok(/direction not resolved/.test(REL.parentChild),
    `the parent-child verdict must not claim a direction: "${REL.parentChild}"`)
  assert.ok(!Object.values(REL).some((v) => v === 'child' || v === 'parent'),
    'no verdict may be a bare "child" or "parent", which is what the symmetric statistic cannot '
    + `support: ${JSON.stringify(REL)}`)
}

// --- 3. A TRANSMITTED CHILD READS AS A FIRST-DEGREE RELATIVE ------------------------------------
{
  const { pg, cg } = trio(31, 40_000)
  const o = oppositeHom(pg, cg)
  assert.ok(o.rate < OPPOSITE_HOM_MAX,
    `a parent and its child cannot be opposite homozygotes, so the rate must sit under the `
    + `${OPPOSITE_HOM_MAX} gate. Got ${o.rate}`)
  assert.ok(o.n > 10_000, 'and there must be enough informative markers to say so')
}

// --- 4. THE THRESHOLDS ARE THE MEASURED ONES ----------------------------------------------------
//
// 181 verified parent-child pairs reach a window spread of at most 0.0238 and 166 constructed
// full-sibling pairs at least 0.0312. If either constant drifts outside that gap this fails.
{
  assert.ok(IBD0_SPREAD_MIN > 0.0238 && IBD0_SPREAD_MIN < 0.0312,
    'the separator must sit between the two measured populations and not inside either: '
    + `got ${IBD0_SPREAD_MIN}`)
  assert.ok(IBD0_FLOOR_MAX > 0 && IBD0_FLOOR_MAX < IBD0_SPREAD_MIN,
    'over the error floor the segmental test has no power, so the floor sits below the separator')
  assert.ok(OPPOSITE_HOM_MAX > 0 && OPPOSITE_HOM_MAX < 0.05,
    'the first-degree gate is a rate, not a spread')
}

// --- 5. TWO UNRELATED GENOMES ARE CALLED UNRELATED ----------------------------------------------
//
// The scorer withholds every named parent on this verdict, so it has to be reached on real
// unrelatedness and not on a thin panel.
{
  const a = trio(41, 30_000)
  const b = trio(97, 30_000)
  const v = relate({ gt: a.pg, pos: a.pos }, { gt: b.pg, pos: a.pos }, OPPOSITE_HOM_MAX)
  assert.equal(v.relationship, REL.unrelated,
    `two independently drawn genomes must read as unrelated, got "${v.relationship}" at rate `
    + `${v.opp.rate}`)
  // And a real parent-child pair must NOT, or the scorer refuses every correct run.
  const t = trio(53, 30_000)
  const ok = relate({ gt: t.pg, pos: t.pos }, { gt: t.cg, pos: t.pos }, OPPOSITE_HOM_MAX)
  assert.notEqual(ok.relationship, REL.unrelated,
    `a transmitted child must not read as unrelated, got "${ok.relationship}"`)
}

// --- 6. TOO LITTLE DATA REFUSES RATHER THAN GUESSING --------------------------------------------
{
  const t = trio(67, 500)
  const v = relate({ gt: t.pg, pos: t.pos }, { gt: t.cg, pos: t.pos }, OPPOSITE_HOM_MAX)
  assert.equal(v.relationship, REL.refused,
    'under ten thousand shared homozygous markers there is nothing to say')
}

console.log('relatedness.check.ts: a probe with no position is skipped rather than crashing the '
  + 'sort, the statistic stays symmetric, and the thresholds are the measured ones')

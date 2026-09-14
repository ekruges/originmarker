/**
 * EVERY HEADLINE RATE WITH ITS INTERVAL, BECAUSE A POINT ESTIMATE AT THIS n IS NOT A NUMBER.
 *
 * WHY THIS FILE EXISTS. This project reports figures like 8 of 8 and 0 of 60. Both are honest and
 * neither is a rate: 8 of 8 is consistent with a true accuracy anywhere from about 63 percent
 * upward, and 0 of 60 does not mean the tool never detects anything. Publishing either bare invites
 * the first question a reviewer asks and answers it badly.
 *
 * WHICH INTERVAL AND WHY. Wilson score intervals, not the normal approximation. At a proportion of
 * 0 or 1 the normal approximation produces an interval of zero width, which is exactly the case
 * this project keeps landing on, and a zero-width interval on 60 observations is a false claim of
 * certainty. Wilson stays sensible at the boundaries.
 *
 * WHERE CLUSTERING MATTERS AND IS SAID SO. Several of these counts are not independent draws.
 * Chromosomes within one array share that array's amplification, and blastomeres within one embryo
 * share a genome. A binomial interval on clustered data is too narrow, measured at 3.5 to 5.8 times
 * too narrow elsewhere in this corpus. Rows below are marked for it rather than silently corrected,
 * because the correction needs the cluster sizes and those differ per arm.
 *
 * Run: node --experimental-strip-types audit/intervals.ts
 */

/** Wilson score interval for a binomial proportion. */
function wilson(k: number, n: number, z = 1.96): [number, number] {
  if (n === 0) return [NaN, NaN]
  const p = k / n
  const d = 1 + (z * z) / n
  const centre = p + (z * z) / (2 * n)
  const half = z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))
  return [Math.max(0, (centre - half) / d), Math.min(1, (centre + half) / d)]
}

/** The rule of three: with zero events in n trials the upper bound is about 3/n. */
const ruleOfThree = (n: number) => 3 / n

interface Row {
  claim: string
  k: number
  n: number
  clustered: boolean
  note?: string
}

/**
 * Every rate this project reports as evidence. The counts are the measured ones; nothing here
 * recomputes them, it only puts an interval around what the harnesses returned.
 */
const rows: Row[] = [
  { claim: 'parent named correctly, trio-resolved breakages', k: 8, n: 8, clustered: true,
    note: '8 events from 4 arrays of one experiment, one sperm donor' },
  { claim: 'parent named correctly, all calls to date', k: 263, n: 263, clustered: true },
  { claim: 'genome-level origin, pronuclei', k: 11, n: 14, clustered: false,
    note: '3 refused rather than wrong' },
  { claim: 'reconstruction verdicts correct', k: 35, n: 41, clustered: false,
    note: '6 refused, 0 incorrect' },
  { claim: 'mirrored losses, paternal direction', k: 12, n: 12, clustered: true },
  { claim: 'mirrored losses, maternal direction', k: 10, n: 12, clustered: true },
  { claim: 'whole-chromosome false positives, adult donors', k: 0, n: 242, clustered: true,
    note: '11 arrays x 22 autosomes, so 11 clusters not 242 draws' },
  { claim: 'false calls, euploid lines and donors, both gates', k: 1, n: 1584, clustered: true,
    note: '72 arrays x 22 autosomes' },
  { claim: 'sensitivity, trisomy 21 on single cells, shipped point', k: 0, n: 60, clustered: false,
    note: 'one chromosome per array, so these are independent' },
  { claim: 'replicate agreement, zygosity', k: 14, n: 14, clustered: false },
  { claim: 'replicate agreement, whole-chromosome set', k: 14, n: 14, clustered: false },
  { claim: 'sex call, hydatidiform moles', k: 13, n: 13, clustered: false },
  { claim: 'sex call, maternal pronuclei carrying no Y', k: 7, n: 7, clustered: false },
  { claim: 'across-family pairs called relatives', k: 0, n: 12, clustered: false,
    note: 'a false positive here would be unambiguous' },
]

const pct = (x: number) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : '  -  ')

console.log('Wilson 95% intervals on every rate this project reports as evidence.')
console.log('A rate at 0 or 1 has a one-sided story: the interval is what bounds it.')
console.log('')
console.log('claim                                                    k/n        point     95% interval')
for (const r of rows) {
  const [lo, hi] = wilson(r.k, r.n)
  const p = r.n ? r.k / r.n : NaN
  console.log(`  ${r.claim.padEnd(54)}${`${r.k}/${r.n}`.padStart(9)}  ${pct(p).padStart(7)}  `
    + `[${pct(lo)}, ${pct(hi)}]${r.clustered ? '  clustered' : ''}`)
}

console.log('')
console.log('=== THE TWO THAT CARRY THE MOST WEIGHT')
{
  const [lo] = wilson(8, 8)
  console.log(`  8 of 8 correct attributions is consistent with a true accuracy as low as ${pct(lo)}.`)
  console.log('  That is why the twelve-family arm matters: it replaces 8 events with about 146,000')
  console.log('  informative markers in each direction.')
}
{
  const [, hi] = wilson(0, 60)
  console.log('')
  console.log(`  0 of 60 detected does NOT mean the tool detects nothing: the upper bound is ${pct(hi)},`)
  console.log(`  and the rule of three puts it near ${pct(ruleOfThree(60))}. What can be said is that`)
  console.log('  sensitivity on this material is below roughly 6 percent, not that it is zero.')
}
{
  const [, hi] = wilson(0, 242)
  console.log('')
  console.log(`  0 false positives in 242 bounds the per-chromosome rate at ${pct(hi)} by Wilson, or`)
  console.log(`  about ${pct(ruleOfThree(242))} by the rule of three. Both assume independent draws, and these`)
  console.log('  are 11 arrays rather than 242, so the honest bound is wider than either.')
}
console.log('')
console.log('CLUSTERING, stated rather than corrected. Chromosomes within an array share that')
console.log('array\'s amplification and are not independent draws. A cluster-robust interval on this')
console.log('corpus measured 3.5 to 5.8 times wider than the naive one. Rows marked clustered above')
console.log('should be read with that factor in mind; the interval printed is the optimistic one.')

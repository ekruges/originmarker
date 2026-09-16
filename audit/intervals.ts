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
 * CLUSTERING, BOUNDED RATHER THAN ASSERTED. Several counts are not independent draws: 242 autosome
 * observations are 11 arrays of 22, and chromosomes within one array share that array's
 * amplification. The dependence is not measurable from a count alone, so each clustered row is
 * printed with BOTH ends of what it could be:
 *
 *   independent      Wilson on all n observations. The optimistic end, true only at zero
 *                    within-array correlation.
 *   fully clustered  Wilson on the m arrays, each contributing one observation. The conservative
 *                    end, true when everything within an array moves together.
 *
 * The real interval lies between them, and a cluster bootstrap would place it. That needs each
 * array's own count, which the harnesses do not yet emit for every arm; rows where the cluster
 * count is not recorded say so rather than borrowing one.
 *
 * A ROW WHOSE OBSERVATIONS ARE ONE PER ARRAY IS NOT CLUSTERED AT ALL, and two rows here were
 * marked as though they were. The mirrored-loss arms are twelve arrays carrying one constructed
 * event each, so their Wilson intervals stand as printed.
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
  /** Independent units behind the n observations. Equal to n where each unit contributes one. */
  clusters?: number
  /** Units carrying the counted outcome, for the conservative end. */
  clusterK?: number
  note?: string
}

/**
 * Every rate this project reports as evidence. The counts are the measured ones; nothing here
 * recomputes them, it only puts an interval around what the harnesses returned.
 */
const rows: Row[] = [
  { claim: 'parent named correctly, trio-resolved breakages', k: 8, n: 8, clusters: 4, clusterK: 4,
    note: '8 events on 4 arrays of one experiment, one sperm donor' },
  { claim: 'parent named correctly, all calls to date', k: 263, n: 263,
    note: 'clustered by array; the harnesses do not record how many arrays' },
  { claim: 'genome-level origin, pronuclei', k: 11, n: 14, clusters: 14, clusterK: 11,
    note: '3 refused rather than wrong' },
  { claim: 'reconstruction verdicts correct', k: 35, n: 41, clusters: 41, clusterK: 35,
    note: '6 refused, 0 incorrect; all 41 share one reconstruction as their reference' },
  { claim: 'mirrored losses, paternal direction', k: 12, n: 12, clusters: 12, clusterK: 12 },
  { claim: 'mirrored losses, maternal direction', k: 10, n: 12, clusters: 12, clusterK: 10 },
  { claim: 'whole-chromosome false positives, adult donors', k: 0, n: 242, clusters: 11, clusterK: 0,
    note: '11 arrays x 22 autosomes' },
  { claim: 'false calls, euploid lines and donors, both gates', k: 1, n: 1584, clusters: 72, clusterK: 1,
    note: '72 arrays x 22 autosomes, the one false call on a single array' },
  { claim: 'heterozygosity-loss findings, one-complement genomes', k: 0, n: 144, clusters: 144, clusterK: 0,
    note: 'moles, sperm and pronuclei, driven with the ploidy call withheld' },
  { claim: 'isodisomy flagged on the stated chromosome, accepted arrays', k: 4, n: 4, clusters: 4, clusterK: 4,
    note: 'two cell lines, two arrays each' },
  { claim: 'moles read as one parental contribution', k: 100, n: 100, clusters: 100, clusterK: 100 },
  { claim: 'sensitivity, trisomy 21 on single cells, shipped point', k: 0, n: 60, clusters: 60, clusterK: 0,
    note: 'one chromosome per array' },
  { claim: 'two-person mixtures tripping any signal', k: 3, n: 4, clusters: 4, clusterK: 3 },
  { claim: 'replicate agreement, zygosity', k: 14, n: 14, clusters: 14, clusterK: 14 },
  { claim: 'replicate agreement, whole-chromosome set', k: 14, n: 14, clusters: 14, clusterK: 14 },
  { claim: 'sex call, hydatidiform moles', k: 13, n: 13, clusters: 13, clusterK: 13 },
  { claim: 'sex call, maternal pronuclei carrying no Y', k: 7, n: 7, clusters: 7, clusterK: 7 },
  { claim: 'across-family pairs called relatives', k: 0, n: 12, clusters: 12, clusterK: 0,
    note: 'a false positive here would be unambiguous' },
]

const pct = (x: number) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : '  -  ')
const band = (k: number, n: number) => {
  const [lo, hi] = wilson(k, n)
  return `[${pct(lo)}, ${pct(hi)}]`
}

console.log('Wilson 95% intervals on every rate this project reports as evidence.')
console.log('A rate at 0 or 1 has a one-sided story: the interval is what bounds it.')
console.log('')
console.log('Two ends are printed where observations share an array. INDEPENDENT assumes none of')
console.log('that sharing; CLUSTERED assumes an array moves as one. The truth is between them.')
console.log('')
console.log('claim                                                    k/n      independent          clustered')
for (const r of rows) {
  const same = r.clusters === r.n
  const right = r.clusters === undefined
    ? 'cluster count not recorded'
    : same ? 'one per array, so the same' : `${r.clusterK}/${r.clusters}  ${band(r.clusterK as number, r.clusters)}`
  console.log(`  ${r.claim.padEnd(58)}${`${r.k}/${r.n}`.padStart(9)}  ${band(r.k, r.n).padEnd(18)}  ${right}`)
}

console.log('')
console.log('=== THE ONES THAT CARRY THE MOST WEIGHT')
{
  const [lo] = wilson(8, 8)
  const [cLo] = wilson(4, 4)
  console.log(`  8 of 8 correct attributions is consistent with a true accuracy as low as ${pct(lo)}`)
  console.log(`  treating the events as independent, and ${pct(cLo)} treating the four arrays as the unit.`)
  console.log('  That is why the twelve-family arm matters: it replaces 8 events with about 146,000')
  console.log('  informative markers in each direction.')
}
{
  const [, hi] = wilson(0, 60)
  console.log('')
  console.log(`  0 of 60 detected does NOT mean the tool detects nothing: the upper bound is ${pct(hi)},`)
  console.log(`  and the rule of three puts it near ${pct(ruleOfThree(60))}. What can be said is that`)
  console.log('  sensitivity on this material is below roughly 6 percent, not that it is zero.')
  console.log('  Each of the 60 is one chromosome on its own array, so no clustering applies.')
}
{
  const [, hi] = wilson(0, 242)
  const [, cHi] = wilson(0, 11)
  console.log('')
  console.log(`  0 false positives in 242 autosome observations bounds the per-chromosome rate at`)
  console.log(`  ${pct(hi)} if those observations are independent. They are 11 arrays of 22, and if an`)
  console.log(`  array moves as one the bound is ${pct(cHi)} per array. Both are true statements about`)
  console.log('  different quantities, and the per-array one is what a reader of a single report wants.')
}
{
  const [, hi] = wilson(0, 144)
  console.log('')
  console.log(`  0 heterozygosity-loss findings on 144 one-complement genomes, each its own array, bounds`)
  console.log(`  that rate at ${pct(hi)} per genome. Before 5.30.0 the same 144 arrays returned 2,747.`)
}
console.log('')
console.log('WHAT WOULD NARROW THE CLUSTERED ROWS. A cluster bootstrap over arrays, which needs each')
console.log('array\'s own numerator and denominator rather than the pooled total. The harnesses emit')
console.log('the pooled figure today, so the two ends above are what can be stated without inventing')
console.log('a within-array correlation.')

/**
 * What two arrays are to each other, and whether one can be the other's parent.
 *
 * MOVED HERE FROM THE COMMAND LINE, WHERE IT WAS THE ONLY COPY. `om link` could tell a parent from
 * a sibling and from an unrelated adult, and the browser could not, because these functions lived
 * in the CLI file. The surface most lab staff use had no way to notice that the array in the
 * parental slot was not a parent at all. Every measurement below is the one the command line has
 * always run; nothing here is new arithmetic.
 *
 * DIRECTION IS NOT IN HERE, and no amount of arithmetic on two genotype files will put it there.
 * Under Hardy-Weinberg the likelihood of a pair factors as P(a)P(b|a) = P(b)P(a|b), so parent-child
 * and child-parent are the same hypothesis, and a tool that says "child" is reading the order of
 * its own arguments. Direction needs something outside the genotypes: the age of the material, a
 * third relative, or the operator's own knowledge of which sample is the embryo.
 */
import type { AB } from './informativity.ts'
import { emptyHet, addOneParent, hetCall } from './obligateHet.ts'

/** A parent array with marker positions, which the windowed test needs to walk in genomic order. */
export interface Positioned {
  gt: Map<string, AB>
  pos: Map<string, { chrom: string; pos: number }>
}

/** Opposite homozygotes. A parent and child cannot be AA and BB at the same marker. */
export function oppositeHom(a: Map<string, AB>, b: Map<string, AB>): { rate: number, n: number } {
  let n = 0
  let opp = 0
  for (const [probe, ga] of a) {
    if (ga !== 'AA' && ga !== 'BB') continue
    const gb = b.get(probe)
    if (gb !== 'AA' && gb !== 'BB') continue
    n += 1
    if (ga !== gb) opp += 1
  }
  return { rate: n ? opp / n : NaN, n }
}

/**
 * Opposite homozygotes in genomic windows, and the spread between the top decile and the median.
 *
 * WHY A SPREAD AND NOT A RATE. Opposite homozygotes are symmetric BY CONSTRUCTION: swap the two
 * files and every count is the same, so nothing built from them can say which array is the parent.
 * What they can say is whether the two genomes share an allele EVERYWHERE. A parent and a child
 * share one allele at every marker, so the only opposite homozygotes between them are genotyping
 * error, and error is spread evenly along the genome. Full siblings share no allele over about a
 * quarter of the genome, and those stretches carry the rate of two unrelated people. So the
 * top-decile window minus the median window is small for a parent-child pair whatever the error
 * rate, and large for siblings. The overall rate cannot make that split: measured over the 106
 * published trios, sibling pairs sit at 0.0114 and parent-child pairs at 0.0037, both under the
 * 0.020 gate, which is why the gate alone called siblings children.
 *
 * The median is the error floor because a pair only reaches this test after passing that gate, and
 * the gate already excludes everything sharing no allele over more than a third of the genome.
 */
export function windowedIbs0(ref: Positioned, s: Map<string, AB>, win = 200) {
  // A PROBE WITH NO POSITION IS SKIPPED RATHER THAN SORTED. The two arrays need not carry the same
  // probe set: panel versions differ, and the position map may come from the sample while the
  // genotypes come from the parent. Sorting on a position that is not there threw
  // `Cannot read properties of undefined` from inside Array.sort, which reaches an operator as a
  // stack trace on a pair of perfectly good files.
  const order = [...ref.gt].filter(([probe, g]) => g !== 'AB' && ref.pos.has(probe))
  // Genomic order, not file order. The statistic is about NEIGHBOURING markers, so a file whose
  // rows are shuffled would spread the sibling signal evenly across windows and read as a parent.
  order.sort(([a], [b]) => {
    const pa = ref.pos.get(a)!
    const pb = ref.pos.get(b)!
    return pa.chrom === pb.chrom ? pa.pos - pb.pos : Number(pa.chrom) - Number(pb.chrom)
  })
  const rates: number[] = []
  let chrom = ''
  let n = 0
  let opp = 0
  for (const [probe, ga] of order) {
    const at = ref.pos.get(probe)!  // guaranteed by the filter above
    // A window never spans two chromosomes, and its unfinished tail is dropped rather than
    // counted at a smaller denominator.
    if (at.chrom !== chrom) { chrom = at.chrom; n = 0; opp = 0 }
    const gb = s.get(probe)
    if (gb !== 'AA' && gb !== 'BB') continue
    n += 1
    if (ga !== gb) opp += 1
    if (n === win) { rates.push(opp / n); n = 0; opp = 0 }
  }
  rates.sort((a, b) => a - b)
  const q = (p: number) => rates[Math.min(rates.length - 1, Math.floor(p * rates.length))]
  const floor = rates.length ? q(0.5) : NaN
  const top = rates.length ? q(0.9) : NaN
  return { windows: rates.length, floor, top, spread: top - floor }
}

/** The one-parent heterozygosity call: is there a second parental contribution, or not. */
export function secondParent(ref: Map<string, AB>, s: Map<string, AB>) {
  const t = emptyHet()
  for (const [probe, pg] of ref) {
    const cg = s.get(probe)
    if (cg) addOneParent(pg as never, cg as never, t as never)
  }
  return hetCall(t as never, 1) as {
    ploidy: string, fraction: number, informative: number, why: string
  }
}

/**
 * The segmental test's constants, measured over the published trios.
 *
 * Of the pairs that pass the opposite-homozygote gate and whose window floor is under
 * IBD0_FLOOR_MAX, the 181 verified parent-child pairs reach a spread of at most 0.0238 and the 166
 * constructed full-sibling pairs a spread of at least 0.0312. The threshold sits between them.
 * Over the floor the test has no power at all: error alone then moves a window further than a
 * missing quarter-genome does.
 *
 * WHAT ELSE CAN PRODUCE THE SIGNAL, and it is not parentage. Two cells of ONE embryo that lost
 * opposite parental copies of a region read as opposite homozygotes across that whole region, and
 * four of 122 same-embryo blastomere pairs do exactly that, at 0.053-0.056. So this test separates
 * a parent-child pair from a pair that is not one; it does not certify that the other pair is
 * siblings. A parent and a child cannot produce it in either direction: whatever copy a child
 * keeps in a region it lost, that copy carries a parental allele.
 */
/**
 * Opposite-homozygote rate above which two arrays are not first-degree relatives.
 *
 * The command line has exposed this as `--opposite-hom-max` with this default since the linkage
 * command existed. It is here so the scorer and the linkage command share one number.
 */
export const OPPOSITE_HOM_MAX = 0.020

export const IBD0_SPREAD_MIN = 0.028
export const IBD0_FLOOR_MAX = 0.012
export const IBD0_MIN_WINDOWS = 100

/** Verdicts `relate` can return. Kept as constants because two commands branch on them. */
export const REL = {
  unrelated: 'unrelated',
  parentChild: 'parent and child, direction not resolved',
  sibling: 'not parent and child: stretches sharing no allele (full siblings)',
  firstDegree: 'first-degree, parent-child and sibling not separated',
  oneParent: 'this parent only, a duplicate or a haploid product',
  ambiguous: 'ambiguous',
  refused: 'refused, too few homozygous markers in common',
} as const

/**
 * What two arrays are to each other, as far as two arrays can say.
 *
 * DIRECTION IS NOT IN HERE, and no amount of arithmetic on two genotype files will put it there.
 * Under Hardy-Weinberg the likelihood of a pair factors as P(a)P(b|a) = P(b)P(a|b), so parent-child
 * and child-parent are the same hypothesis; the tool that says "child" is reading the order of its
 * own arguments. Direction needs something outside the genotypes: the age of the material, a third
 * relative, or the caller's own knowledge of which sample is the embryo.
 */
export function relate(ref: Positioned, s: Positioned, oppMax: number) {
  const opp = oppositeHom(ref.gt, s.gt)
  const link = secondParent(ref.gt, s.gt)
  const win = windowedIbs0(ref, s.gt)
  const verdict = (): string => {
    if (opp.n < 10_000) return REL.refused
    if (!(opp.rate <= oppMax)) return REL.unrelated
    if (link.ploidy === 'uniparental') return REL.oneParent
    if (link.ploidy !== 'biparental') return REL.ambiguous
    if (win.windows < IBD0_MIN_WINDOWS || !(win.floor <= IBD0_FLOOR_MAX)) return REL.firstDegree
    return win.spread >= IBD0_SPREAD_MIN ? REL.sibling : REL.parentChild
  }
  return { opp, link, win, relationship: verdict() }
}

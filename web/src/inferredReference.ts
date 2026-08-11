/**
 * Reconstruct a parental genotype from haploid meiotic products, in the browser.
 *
 * A haploid pronucleus is one meiotic product. Where the parent is homozygous every product
 * carries that allele; where he is heterozygous each product carries one at random, so a
 * heterozygous site is recognised only by observing both alleles across products. Products that
 * happen to agree at a heterozygous site promote it into the reference as homozygous, and a true
 * offspring then reads as missing the parental allele at half of those.
 *
 * Two independent quantities control that error and must not be conflated:
 *
 *   m   how many products were called at a marker, which sets how often a heterozygous site
 *       survives unrecognised
 *   n   how many products exist, which sets how many markers reach any given m
 *
 * This is the browser half of `tools/inferred_reference.py` and must agree with it numerically:
 * both are pinned by the cross-implementation fixture, which is what has caught every divergence
 * in this codebase so far.
 *
 * The storage differs from the Python deliberately. Python holds a dict of per-marker
 * observation lists, which is fine on a workstation and ruinous in a tab: eight arrays of
 * 825,657 probes would be tens of millions of boxed tuples. Here the probe identifiers are
 * mapped to integers once and each product becomes one Uint8Array of allele codes, so eight
 * products cost about 6.6 MB in total and every operation is a linear scan over typed arrays.
 */
import { isAutosome } from './parentage.ts'
import type { AB } from './informativity.ts'
import type { ProbeRow } from './ingest.ts'

/** Allele codes. 0 means the product said nothing usable here. */
const NONE = 0
const A = 1
const B = 2

/**
 * Measured excess of P(all m products agree | parent heterozygous) over the independent-draw
 * value 2^(1-m). Products are not independent at a heterozygous site: a probe with
 * allele-specific dropout calls the same homozygote in every product, and those are exactly the
 * low-m probes.
 *
 * Measured against a real parental array over 133,631 heterozygous autosomal markers, using six
 * haploid products of that parent:
 *
 *     m      markers   P(all agree | het)   2^(1-m)   ratio
 *     2        4,332              0.61958   0.50000   1.24
 *     3       10,588              0.35304   0.25000   1.41
 *     4       23,830              0.17709   0.12500   1.42
 *     5       42,153              0.08023   0.06250   1.28
 *     6       49,196              0.03132   0.03125   1.00
 *
 * The excess vanishes by m=6, so a reference deep enough to be usable needs no correction.
 */
const AGREEMENT_EXCESS: Record<number, number> = { 2: 1.24, 3: 1.41, 4: 1.42, 5: 1.28 }

/** P(all m products show the same allele | the parent is heterozygous). */
export const pAgree = (m: number): number =>
  Math.min(1, 2 ** (1 - m) * (AGREEMENT_EXCESS[m] ?? 1))

/**
 * Largest tolerated drift of the retained marker set away from the genome, as a fraction of the
 * parent's heterozygosity.
 *
 * Raising m lowers contamination and narrows the marker set toward probes nearly every product
 * called, which are enriched for sites where the parent is homozygous because the minor allele
 * is rare. Unrelated people carry the common allele there too, so past a point the reference
 * loses its grip on them. Measured on a father whose real array gives 16.66% heterozygosity:
 *
 *     m>=   markers   h(retained)   of true   contamination   products back   closest negative
 *      2    663,542       16.93%      102%         4.288%          0 of 5           5.96x
 *      3    596,687       16.48%       99%         3.402%          2 of 5           5.65x
 *      4    457,506       15.47%       93%         2.286%          3 of 5           5.25x
 *      5    231,064       13.74%       82%         1.258%          3 of 5           4.67x
 *
 * Recovery plateaus at m=4 while everything else keeps degrading, so 0.90 sits at the knee.
 */
export const MIN_ASCERTAINMENT = 0.90

/** Below this many products the method inverts: true offspring read as decisively absent. */
export const MIN_PRODUCTS = 5

export interface ReferenceStats {
  nProducts: number
  mMin: number
  markers: number
  /** The PARENT's heterozygosity over markers passing the m filter, from product disagreement.
   *  Measures ascertainment. NOT the reference's own heterozygosity, which is `contamination`. */
  hRetained: number
  /** Fraction of the reference that is a heterozygous site mistaken for homozygous. */
  contamination: number
  /** Absence a HAPLOID true offspring reads at those markers. A diploid offspring must also
   *  receive the minor allele from its other parent, so it lands well short of half. */
  spuriousAbsence: number
  meanM: number
  excluded: string[]
}

/** A reconstructed genotype, plus everything a reader needs to judge it. */
export interface Reference extends ReferenceStats {
  /** probe id -> the homozygous call the reference asserts there. */
  genotype: Map<string, AB>
}

export class ProductSet {
  readonly ids: string[] = []
  /** probe id -> dense index, assigned once and shared by every product. */
  private index = new Map<string, number>()
  private chrom: string[] = []
  /** Kept alongside the chromosome so the reference can be written out as an array file. The
   *  reconstruction itself never needs a coordinate; every consumer of the export does. */
  private pos: number[] = []
  /** Probe id by dense index, the inverse of `index`. Kept so a build can walk probes by
   *  integer instead of materialising `[...index.entries()]`, which allocated 825,657 pairs
   *  per build and there are a dozen builds in a run. */
  private probeId: string[] = []
  private alleles: Uint8Array[] = []
  private capacity = 0
  /** Per product diagnostics, in the order of `ids`. */
  readonly called: number[] = []
  readonly hetCalls: number[] = []
  readonly band: number[] = []

  private grow(need: number): void {
    if (need <= this.capacity) return
    const next = Math.max(need, this.capacity ? this.capacity * 2 : 1 << 20)
    for (let i = 0; i < this.alleles.length; i += 1) {
      const wider = new Uint8Array(next)
      wider.set(this.alleles[i])
      this.alleles[i] = wider
    }
    this.capacity = next
  }

  /** Start a product. Rows are then fed one at a time, as the file streams. */
  begin(id: string): number {
    this.ids.push(id)
    this.alleles.push(new Uint8Array(this.capacity))
    this.called.push(0)
    this.hetCalls.push(0)
    this.band.push(NaN)
    return this.ids.length - 1
  }

  /**
   * One row of the product currently being read. Only homozygous autosomal calls become
   * evidence: a heterozygous call in a haploid genome is an error by construction, and counting
   * it as evidence of parental heterozygosity inflates the recovered figure several-fold.
   */
  add(slot: number, row: ProbeRow, baf: { inBand: number; total: number }): void {
    if (!isAutosome(row.chrom)) return
    if (row.baf !== null) {
      baf.total += 1
      if (row.baf >= 0.35 && row.baf <= 0.65) baf.inBand += 1
    }
    if (row.genotype === 'NC') return
    this.called[slot] += 1
    if (row.genotype === 'AB') { this.hetCalls[slot] += 1; return }
    if (row.genotype !== 'AA' && row.genotype !== 'BB') return
    let i = this.index.get(row.probesetId)
    if (i === undefined) {
      i = this.index.size
      this.index.set(row.probesetId, i)
      this.chrom.push(row.chrom)
      this.pos.push(row.pos)
      this.probeId.push(row.probesetId)
      this.grow(i + 1)
    }
    this.alleles[slot][i] = row.genotype === 'AA' ? A : B
  }

  /** Close a product once its file has streamed. */
  end(slot: number, baf: { inBand: number; total: number }): void {
    this.band[slot] = baf.total ? baf.inBand / baf.total : NaN
  }

  get size(): number { return this.ids.length }

  /** Distinct markers any product called homozygous, which is what every pass walks. */
  get markerCount(): number { return this.index.size }

  /** Where a probe sits, for writing the reference out. Empty for a probe never seen homozygous
   *  in any product, which is also a probe the reference has no call at. */
  locus(id: string): { chrom: string; pos: number } | null {
    const i = this.index.get(id)
    return i === undefined ? null : { chrom: this.chrom[i], pos: this.pos[i] }
  }

  private keep(exclude: readonly string[]): number[] {
    const drop = new Set(exclude)
    return this.ids.map((_, i) => i).filter((i) => !drop.has(this.ids[i]))
  }

  /** A reference from every marker where at least mMin products spoke and all agreed. */
  /**
   * The parent's heterozygosity over the markers a given depth retains, and nothing else.
   *
   * `chooseM` needs this scalar at every candidate depth and used to get it by running a full
   * build, which meant constructing and discarding a 600,000-entry map per candidate. Choosing
   * the threshold cost several times the build it was choosing for.
   */
  hRetainedAt(mMin: number, exclude: readonly string[] = []): number {
    const use = this.keep(exclude)
    const n = this.index.size
    let disagreed = 0
    let surviving = 0
    let detect = 0
    for (let i = 0; i < n; i += 1) {
      let m = 0
      let first = NONE
      let mixed = false
      for (const p of use) {
        const v = this.alleles[p][i]
        if (v === NONE) continue
        m += 1
        if (first === NONE) first = v
        else if (v !== first) mixed = true
      }
      if (m < mMin) continue
      surviving += 1
      detect += 1 - pAgree(m)
      if (mixed) disagreed += 1
    }
    return surviving && detect > 0 ? (disagreed / surviving) / (detect / surviving) : NaN
  }

  build(mMin: number, exclude: readonly string[] = [], hPrior?: number): Reference {
    const use = this.keep(exclude)
    const n = this.index.size
    const genotype = new Map<string, AB>()
    let mSum = 0
    let disagreed = 0
    let surviving = 0
    let detect = 0
    // How many products spoke at each marker the reference kept, in marker order. Contamination
    // used to be a SECOND full pass over every probe with a map lookup per probe to find these
    // again; they are free to record here and the pass is over the kept markers alone.
    const keptM: number[] = []

    for (let i = 0; i < n; i += 1) {
      let m = 0
      let first = NONE
      let mixed = false
      for (const p of use) {
        const v = this.alleles[p][i]
        if (v === NONE) continue
        m += 1
        if (first === NONE) first = v
        else if (v !== first) mixed = true
      }
      if (m < mMin) continue
      surviving += 1
      detect += 1 - pAgree(m)
      if (mixed) { disagreed += 1; continue }
      genotype.set(this.probeId[i], first === A ? 'AA' : 'BB')
      keptM.push(m)
      mSum += m
    }

    const hRetained = surviving && detect > 0
      ? (disagreed / surviving) / (detect / surviving) : NaN
    const h = hPrior !== undefined && Number.isFinite(hPrior) ? hPrior : hRetained

    // Contamination is computed per marker from its own m, then averaged, because a marker seen
    // by ten products is far safer than one seen by the minimum, and averaging the exponent
    // instead of the probability would understate the tail.
    let contamination = NaN
    if (genotype.size && Number.isFinite(h)) {
      let acc = 0
      // Same values in the same order the two-pass version summed them in, so the float result
      // is identical and the cross-implementation fixture still pins.
      for (const m of keptM) {
        const odds = h * pAgree(m)
        acc += odds / (1 - h + odds)
      }
      contamination = acc / genotype.size
    }

    return {
      genotype,
      nProducts: use.length,
      mMin,
      markers: genotype.size,
      hRetained,
      contamination,
      spuriousAbsence: contamination / 2,
      meanM: genotype.size ? mSum / genotype.size : NaN,
      excluded: [...exclude],
    }
  }

  /**
   * The largest m whose marker set still resembles the genome.
   *
   * Needs no real parental array. h at m=2 is measured unbiased against a known father (16.93%
   * against 16.66%), so it serves as the baseline the stricter settings are judged against, and
   * the whole rule is computable from the products alone.
   *
   * An earlier rule fixed m at n-1. That happens to be right at n=5 and stops being right as n
   * grows, because what governs the drift is the absolute number of products required to agree,
   * not the slack below n: at n=8 it demanded seven and left the retained set at roughly 62%
   * genome's heterozygosity, which cost the reference its grip on unrelated diploids.
   */
  chooseM(exclude: readonly string[] = []): { mMin: number; ratios: Map<number, number> } {
    const n = this.keep(exclude).length
    const base = this.hRetainedAt(2, exclude)
    const ratios = new Map<number, number>()
    let best = 2
    for (let m = 2; m <= n; m += 1) {
      const r = this.hRetainedAt(m, exclude) / base
      ratios.set(m, r)
      if (Number.isFinite(r) && r >= MIN_ASCERTAINMENT) best = m
    }
    return { mMin: best, ratios }
  }

  /**
   * Opposite-homozygote rate between two products, which settles whether they share a parent
   * with no reference at all.
   *
   * Two haploid products of one man differ only where he is heterozygous and they drew
   * differently, so this needs no reference and is what decides membership. It does NOT decide
   * it one pair at a time: the ranges overlap by 0.18 of a percentage point and a real
   * cross-father pair sits under the same-parent cut. See `SAME_PARENT_MAX` below for the
   * measured ranges and `groupByParent` for the all-pairs rule that follows from them.
   */
  opposite(a: number, b: number): { shared: number; opposite: number; rate: number } {
    const x = this.alleles[a]
    const y = this.alleles[b]
    let shared = 0
    let opp = 0
    const n = this.index.size
    for (let i = 0; i < n; i += 1) {
      const u = x[i]
      const v = y[i]
      if (u === NONE || v === NONE) continue
      shared += 1
      if (u !== v) opp += 1
    }
    return { shared, opposite: opp, rate: shared ? opp / shared : NaN }
  }
}

/**
 * Opposite-homozygote thresholds, and a warning about how much weight one pair can carry.
 *
 * Measured across two experiments: 4.68% to 9.70% between products of one father over 46 pairs,
 * and 9.88% to 16.10% between products of different fathers over 45. The separation is real but
 * the margin is 0.18 of a percentage point, not the empty void a first look suggested, and one
 * genuine cross-father pair sits at 9.88%, under the cut below.
 *
 * So a per-pair label is not evidence of group membership. Grouping uses all-pairs consistency
 * instead, which turns that one misread pair into no group-level error at all.
 */
export const SAME_PARENT_MAX = 0.105
export const DIFFERENT_PARENT_MIN = 0.125

export type Kinship = 'same parent' | 'different parents' | 'ambiguous'

export const kinship = (rate: number): Kinship =>
  (rate < SAME_PARENT_MAX ? 'same parent'
    : rate >= DIFFERENT_PARENT_MIN ? 'different parents' : 'ambiguous')

/**
 * Split products into groups sharing one parent. EVERY pair inside a group must read same
 * parent; a chain of pairwise links is not enough.
 *
 * Single linkage fails on real data. Across the two experiments one genuine cross-father pair
 * reads 9.88%, under the same-parent cut, so connected components chain two men into one group
 * through that single edge. Requiring all pairs rejects it, because the same product reads 12.1%
 * to 12.4% against the rest of that group: one misread pair in 45 becomes zero group errors.
 *
 * Products must already be quality-gated and haploid before this runs. A diploid compared with a
 * haploid has a different expected rate entirely and bridges unrelated groups, and an array below
 * the call-rate floor reads high against everyone, which splits one father into several.
 */
export function groupByParent(
  n: number, rate: (a: number, b: number) => number,
): number[][] {
  const same = (a: number, b: number): boolean => kinship(rate(a, b)) === 'same parent'
  const adj: Set<number>[] = Array.from({ length: n }, (_, a) => {
    const s = new Set<number>()
    for (let b = 0; b < n; b += 1) if (b !== a && same(a, b)) s.add(b)
    return s
  })

  /**
   * The largest set in which every pair reads same parent, found exactly.
   *
   * Greedy clique-building is not sufficient and the failure is not hypothetical. With two
   * parents of three products each and one genuine cross-parent pair reading under the cut, a
   * greedy pass seeded on either member of that pair returns the pair itself as a group and
   * splits both parents in half. Exact search returns the two real groups. Bron-Kerbosch with
   * pivoting, over at most a few dozen products, which is what an experiment holds.
   */
  const best: { members: number[] } = { members: [] }
  const expand = (R: number[], P: Set<number>, X: Set<number>): void => {
    if (!P.size && !X.size) {
      if (R.length > best.members.length) best.members = [...R]
      return
    }
    // Pivot on the candidate with the most neighbours, which prunes the branches that cannot
    // beat it.
    let pivot = -1
    let deg = -1
    for (const u of [...P, ...X]) {
      const d = [...P].filter((v) => adj[u].has(v)).length
      if (d > deg) { deg = d; pivot = u }
    }
    for (const v of [...P]) {
      if (pivot >= 0 && adj[pivot].has(v)) continue
      expand(
        [...R, v],
        new Set([...P].filter((u) => adj[v].has(u))),
        new Set([...X].filter((u) => adj[v].has(u))),
      )
      P.delete(v)
      X.add(v)
    }
  }

  const remaining = new Set(Array.from({ length: n }, (_, i) => i))
  const groups: number[][] = []
  while (remaining.size) {
    best.members = []
    expand([], new Set(remaining), new Set())
    // A product sharing no parent with anything else is its own group of one.
    const clique = best.members.length ? best.members : [[...remaining][0]]
    clique.forEach((x) => remaining.delete(x))
    groups.push([...clique].sort((a, b) => a - b))
  }
  return groups.sort((a, b) => b.length - a.length || a[0] - b[0])
}


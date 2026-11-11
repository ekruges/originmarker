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
import { emptyShareTally, addShare, shareCall, type ShareCall } from './oneParentState.ts'

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
  /** Markers dropped because their window is one where every product carries the same haplotype.
   *  Zero when no founder frequencies were given, and then no window is masked. */
  masked: number
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
  /** 1 where a product called the marker heterozygous. Never evidence for the reconstruction, which
   *  reads `alleles` alone; read only by `productAgainst`, where a diploid shows itself by them. */
  private hets: Uint8Array[] = []

  private grow(need: number): void {
    if (need <= this.capacity) return
    const next = Math.max(need, this.capacity ? this.capacity * 2 : 1 << 20)
    for (let i = 0; i < this.alleles.length; i += 1) {
      const wider = new Uint8Array(next)
      wider.set(this.alleles[i])
      this.alleles[i] = wider
      const h = new Uint8Array(next)
      h.set(this.hets[i])
      this.hets[i] = h
    }
    this.capacity = next
  }

  /** Start a product. Rows are then fed one at a time, as the file streams. */
  begin(id: string): number {
    this.ids.push(id)
    this.alleles.push(new Uint8Array(this.capacity))
    this.hets.push(new Uint8Array(this.capacity))
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
    if (row.genotype === 'AB') this.hetCalls[slot] += 1
    else if (row.genotype !== 'AA' && row.genotype !== 'BB') return
    let i = this.index.get(row.probesetId)
    if (i === undefined) {
      i = this.index.size
      this.index.set(row.probesetId, i)
      this.chrom.push(row.chrom)
      this.pos.push(row.pos)
      this.probeId.push(row.probesetId)
      this.grow(i + 1)
    }
    if (row.genotype === 'AB') this.hets[slot][i] = 1
    else this.alleles[slot][i] = row.genotype === 'AA' ? A : B
  }

  /** Close a product once its file has streamed. */
  end(slot: number, baf: { inBand: number; total: number }): void {
    this.band[slot] = baf.total ? baf.inBand / baf.total : NaN
  }

  get size(): number { return this.ids.length }

  /** Distinct markers any product called, which is what every pass walks. A marker only ever called
   *  heterozygous has no allele code anywhere, so every reconstruction pass skips it. */
  get markerCount(): number { return this.index.size }

  /** Where a probe sits, for writing the reference out. Empty for a probe no product called. */
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

  build(
    mMin: number, exclude: readonly string[] = [], hPrior?: number,
    freqB?: ReadonlyMap<string, number>,
  ): Reference {
    const use = this.keep(exclude)
    const masked = freqB ? this.uniformWindows(use, freqB) : null
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
    let maskedN = 0

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
      if (masked && masked.has(this.windowOf(i))) { maskedN += 1; continue }
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
      masked: maskedN,
    }
  }

  /** The 5 Mb window a marker sits in, as `chrom:index`. */
  private windowOf(i: number): string {
    return `${this.chrom[i]}:${Math.floor(this.pos[i] / UNIFORM_WINDOW_BP)}`
  }

  /**
   * Windows where every product carries the same haplotype, which the reconstruction cannot see.
   *
   * WHY. With n products, 2 x 0.5^n of any stretch where the parent is heterozygous is expected to be
   * carried identically by all of them: 6.25% at five products. There the products agree, so the
   * reconstruction asserts a homozygote at the parent's heterozygous sites, and every true offspring
   * reads the parent's allele absent there. The parity consult measured 5 Mb windows of a five-product
   * rebuild whose contamination against the real father reached 0.15 against 0.019 genome-wide, and a
   * published maternal monosomy read wrongly because of one.
   *
   * THE TEST NEEDS NO PARENT. Per window, how often the products disagree, over how heterozygous the
   * founders are there, relative to the genome's median of that ratio. A window where the parent is
   * heterozygous and the products happened to split reads near the median; one where they all drew
   * the same homozygous haplotype reads near zero. Under UNIFORM_RELATIVE_MAX of the median it is
   * masked: the reference makes no claim there. Measured: at 0.6 it flags 42 of 533 windows and catches
   * all 21 whose disagreement, against the real father, is under 0.35 of what it should be.
   */
  uniformWindows(use: readonly number[], freqB: ReadonlyMap<string, number>): Set<string> {
    const n = this.index.size
    const seen = new Map<string, { markers: number; disagreed: number; het: number }>()
    for (let i = 0; i < n; i += 1) {
      const q = freqB.get(this.probeId[i])
      if (q === undefined || q <= 0 || q >= 1) continue
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
      if (m < 2) continue
      const w = this.windowOf(i)
      let t = seen.get(w)
      if (!t) seen.set(w, t = { markers: 0, disagreed: 0, het: 0 })
      t.markers += 1
      if (mixed) t.disagreed += 1
      t.het += 2 * q * (1 - q)
    }
    const ratio = new Map<string, number>()
    for (const [w, t] of seen) {
      if (t.markers < UNIFORM_MIN_MARKERS || t.het <= 0) continue
      ratio.set(w, (t.disagreed / t.markers) / (t.het / t.markers))
    }
    const sorted = [...ratio.values()].sort((a, b) => a - b)
    const median = sorted.length ? sorted[sorted.length >> 1] : NaN
    const out = new Set<string>()
    if (!(median > 0)) return out
    for (const [w, r] of ratio) if (r < UNIFORM_RELATIVE_MAX * median) out.add(w)
    return out
  }

  /**
   * How one product reads against a reference built from the OTHERS: the one-parent copy state taken
   * over its whole genome, from its own homozygous calls.
   *
   * WHY THIS AND NOT HETEROZYGOSITY. A haploid product of the parent carries that parent's alleles
   * alone, so against a reference of the parent the allele the parent lacks appears only by error,
   * whatever its drop-in. Heterozygosity cannot say that: drop-in raises it, and three of the eight
   * paternal pronuclei of GSE148488, haploid by micromanipulation, read 0.25 to 0.43 heterozygous and
   * were refused as not a genome. Read this way the consult found all eight at 0.061 to 0.125, every
   * diploid offspring at 0.355 or more, products of the other parent at 0.962 or more, and arrays of
   * noise reading no allele contrast at all.
   */
  productAgainst(slot: number, ref: Reference, freqB: ReadonlyMap<string, number>,
    monomorphic: (id: string, g: 'AA' | 'BB') => boolean): ShareCall {
    const t = emptyShareTally()
    const alleles = this.alleles[slot]
    const hets = this.hets[slot]
    for (let i = 0; i < this.index.size; i += 1) {
      const v = alleles[i]
      if (v === NONE && !hets[i]) continue
      const id = this.probeId[i]
      const g = ref.genotype.get(id)
      if (g !== 'AA' && g !== 'BB') continue
      // A call is read as the B allele frequency it stands for: 0 or 1 homozygous, 0.5 heterozygous.
      // The heterozygous calls carry the test: a diploid offspring is heterozygous wherever its
      // other parent gave the other allele, and that is r near 0.5 where a product reads near 0.
      addShare(t, g, v === NONE ? 0.5 : v === B ? 1 : 0, freqB.get(id), monomorphic(id, g))
    }
    return shareCall(t)
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
/** The window, in base pairs, over which the reconstruction's uniform stretches are found. */
export const UNIFORM_WINDOW_BP = 5_000_000
/** Markers a window needs before its disagreement is a measurement. */
export const UNIFORM_MIN_MARKERS = 50
/** A window disagreeing at under this share of the genome's median ratio is masked. */
export const UNIFORM_RELATIVE_MAX = 0.6

/** Products a haploid-staged group needs before it can judge the arrays outside it. One under the
 *  build floor: the depth leave-one-out already judges every member of a group at the floor. */
export const MIN_SEED = MIN_PRODUCTS - 1

export type ProductVerdict = 'product' | 'not this parent' | 'refused'

export interface JudgedProduct {
  id: string
  verdict: ProductVerdict
  call: ShareCall
  /** Judged against a rebuild from the group without it, or, for an array outside the group,
   *  against the group's whole rebuild. */
  against: 'the others' | 'the group'
}

const verdictOf = (c: ShareCall): ProductVerdict =>
  c.state === 'other-lost' ? 'product' : c.state === 'parent-lost' ? 'not this parent' : 'refused'

/**
 * Which arrays are products of the parent a group was built from. Each member is judged against a
 * reference built from the OTHER members; each candidate, an array outside the group, against the
 * whole group's. A product is one when it reads the other parent's copy absent across its genome
 * (the copy state's `other-lost`, r under 0.20); not this parent's when it reads the parent's own
 * absent (`parent-lost`, over 0.75); refused otherwise, which is a diploid, a mixture, or calls with
 * no allele contrast (`no-template`).
 *
 * WHY CANDIDATES. The ploidy triage refuses on heterozygosity, which drop-in raises on a true haploid
 * product. The test here does not move with drop-in, so a product the triage set aside comes back in
 * on its alleles. Grouping still runs on the triaged products alone: a diploid admitted to the pair
 * table splits one parent's group in two.
 */
export function judgeProducts(
  ps: ProductSet, members: readonly string[], candidates: readonly string[], mMin: number,
  freqB: ReadonlyMap<string, number>, monomorphic: (id: string, g: 'AA' | 'BB') => boolean,
): JudgedProduct[] {
  const inGroup = new Set(members)
  const outside = ps.ids.filter((id) => !inGroup.has(id))
  const slot = (id: string): number => ps.ids.indexOf(id)
  const out: JudgedProduct[] = []
  for (const id of members) {
    const ref = ps.build(Math.min(mMin, members.length - 1), [...outside, id], undefined, freqB)
    const call = ps.productAgainst(slot(id), ref, freqB, monomorphic)
    out.push({ id, verdict: verdictOf(call), call, against: 'the others' })
  }
  if (candidates.length) {
    const ref = ps.build(mMin, outside, undefined, freqB)
    for (const id of candidates) {
      const call = ps.productAgainst(slot(id), ref, freqB, monomorphic)
      out.push({ id, verdict: verdictOf(call), call, against: 'the group' })
    }
  }
  return out
}

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


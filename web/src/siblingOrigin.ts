/**
 * Which parental copy is missing, using the embryo's OWN unaffected cells as the reference.
 *
 * This is the construction that removes the parent requirement. Multiple blastomeres of one embryo
 * share both parents exactly, so where an event is present in some cells and absent in others, the
 * unaffected siblings define that embryo's heterozygous sites and the affected cell is read against
 * them. No parental array is involved at any point.
 *
 * THREE THINGS THAT LOOK LIKE DETAIL AND ARE NOT.
 *
 * 1. HETEROZYGOSITY MUST COME FROM THE SIBLINGS, NEVER FROM THE CELL BEING SCORED. Selecting on the
 *    affected cell's own heterozygous calls selects exactly the markers where both alleles survived,
 *    which forces the share to 0.5 at every retained marker and leaves nothing to read: simulated
 *    mean 0.5000 with zero markers at either tail. Sibling-established heterozygosity keeps both
 *    tails, 15.3% at each. It also avoids re-deriving the region from heterozygosity, which is the
 *    failure behind every artefact in this repo's audit.
 *
 * 2. THERE IS NO CENTRE HERE, DELIBERATELY. The one-directional results this project produced came
 *    from centring a share statistic on the median of a distribution whose median IS its boundary.
 *    This module never forms a share and never centres: each marker contributes a likelihood under
 *    three hypotheses and the region's evidence is their product. A quantity that is never centred
 *    cannot be mis-centred.
 *
 * 3. THIS MODULE DOES NOT PHASE, AND CANNOT NAME A PARENT ON ITS OWN. `observations` must already
 *    be oriented to a CONSISTENT HAPLOTYPE, meaning 'AA' denotes the same parental side at every
 *    marker in the region. That orientation is not free: the reference allele is defined per
 *    marker, so a region that lost one parent's copy retains the A allele at some markers and the
 *    B allele at others, and counting AA against BB over unphased markers cancels to noise. Given
 *    phased input the two hypotheses below are the two parental sides; given unphased input they
 *    are not, and the result is meaningless rather than merely weak.
 *
 *    So the honest decomposition is: this module answers DELETION versus DROPOUT, which needs no
 *    orientation, and answers WHICH SIDE only when the caller has supplied phase. Naming that side
 *    maternal or paternal needs one anchor per donor group on top of that.
 *
 * 4. FALSE HETEROZYGOSITY IS THE ONLY ONE-DIRECTIONAL TERM. A marker that is truly homozygous but
 *    called heterozygous by the siblings always retains the reference allele, so it always reads as
 *    loss of the OTHER copy. Its rate, phi, is therefore reported with every call and bounds it.
 *    Above PHI_WALL the call is refused outright rather than downgraded, because past that point the
 *    false-het floor exceeds the separation between the hypotheses and no marker count rescues it.
 */

/** Genotype as this codebase writes it elsewhere. */
export type AB = 'AA' | 'AB' | 'BB' | 'NC'

/**
 * Drop-in on this platform: a heterozygous call at a truly homozygous marker.
 *
 * MEASURED, on 113 same-genome array pairs across four experiments, being replicate and clonal
 * material where a marker the reference calls homozygous is homozygous in truth:
 *
 *     ROBLES  n=40  median 0.0525      DIETER  n=25  median 0.0435
 *     JENNA   n=16  median 0.0410      TREFF   n=32  median 0.0390
 *
 * The pooled median is used. An external assessment assumed 0.0854 from 7 arrays; every experiment
 * here sits below that, and the worst single pair observed was 0.0646.
 */
export const DROP_IN = 0.0435

/** Genotyping error other than drop-in: a homozygote called as the wrong homozygote. */
export const GENOTYPE_ERROR = 0.02

/**
 * False-het fraction above which no call is made.
 *
 * This is a wall rather than a slope. Past roughly 0.15 the reference-copy-lost hypothesis becomes
 * unresolvable at any marker count, because the floor contributed by false heterozygotes exceeds
 * the separation between hypotheses. Refuse rather than report weakly.
 */
export const PHI_WALL = 0.15

/** Posterior a hypothesis must reach before the region is called. */
export const CALL_POSTERIOR = 0.95

/** Fewest informative markers before a region is scored at all. */
export const MIN_MARKERS = 20

/** Fewest event-free sibling cells that can establish heterozygosity. */
export const MIN_SIBLINGS = 2

/**
 * How many siblings must agree before a marker counts as heterozygous.
 *
 * SCALES WITH THE PANEL, and must. At a fixed rule of two, phi RISES as siblings are added, because
 * more cells give a spurious heterozygous call more chances to appear somewhere: at the measured
 * drop-in it runs 0.051 at two clean siblings and 0.339 at eight, crossing the wall at five. So
 * adding arrays would make the bias worse while appearing to make the call more reliable, which is
 * exactly the shape of error that produced this project's one-directional results.
 *
 * A majority rule removes it: the same panel reads 0.023 at five siblings and 0.006 at seven.
 */
export const hetRule = (siblings: number): number => Math.max(MIN_SIBLINGS, Math.ceil(siblings / 2))

/**
 * Expected false-het fraction for a panel of `siblings` cells under a majority rule.
 *
 * A truly homozygous marker is called heterozygous only if at least `hetRule` of the siblings
 * independently drop in. Binomial upper tail.
 */
export function expectedPhi(siblings: number, dropIn = DROP_IN): number {
  const need = hetRule(siblings)
  if (siblings < need) return 0
  let tail = 0
  for (let k = need; k <= siblings; k += 1) {
    let c = 1
    for (let i = 0; i < k; i += 1) c = (c * (siblings - i)) / (i + 1)
    tail += c * dropIn ** k * (1 - dropIn) ** (siblings - k)
  }
  return tail
}

/**
 * The two loss hypotheses are HAPLOTYPE SIDES, not parents, and are only distinguishable when the
 * caller has phased the region. `no-deletion` is meaningful either way.
 */
export type Hypothesis = 'reference-copy-lost' | 'other-copy-lost' | 'no-deletion'

export interface SiblingCall {
  /** True when the caller supplied phased observations, so a side call means a parental side. */
  phased: boolean
  /** The parental copy missing, expressed relative to the reference allele, or a refusal. */
  hypothesis: Hypothesis | 'refused'
  posterior: number
  /** Informative markers used: embryo heterozygous by sibling consensus, sample called. */
  markers: number
  /** Sibling cells that established heterozygosity, and how many had to agree. */
  siblings: number
  agreement: number
  /** Expected fraction of the informative set that is falsely heterozygous. */
  phi: number
  why: string
}

/**
 * One marker's likelihood under each hypothesis.
 *
 * `refAllele` is whichever allele the region is being scored relative to; the caller fixes it once
 * per region so the two loss hypotheses stay distinguishable. A falsely heterozygous marker is
 * truly homozygous for it, and so reads as the other copy having been lost whatever the truth is.
 */
function likelihoods(observed: AB, phi: number, ado: number, eps: number): [number, number, number] {
  // Under a real deletion the cell carries ONE allele, so a heterozygous observation is error.
  const lostRef = observed === 'BB' ? 1 - eps : observed === 'AA' ? eps / 2 : eps / 2
  const lostOther = observed === 'AA' ? 1 - eps : observed === 'BB' ? eps / 2 : eps / 2
  // With both copies present the cell is heterozygous unless an allele dropped.
  const noDel = observed === 'AB' ? 1 - ado : ado / 2
  // A falsely heterozygous marker is homozygous for the reference allele, so it looks like the
  // other copy was lost regardless of which hypothesis is true. This is the one-directional term.
  const asRef = observed === 'AA' ? 1 - eps : eps / 2
  return [
    (1 - phi) * lostRef + phi * asRef,
    (1 - phi) * lostOther + phi * asRef,
    (1 - phi) * noDel + phi * asRef,
  ]
}

/**
 * Call which parental copy is missing across a region.
 *
 * `observations` are the affected cell's genotypes at markers the SIBLINGS established as
 * heterozygous, oriented so that 'AA' is homozygous for the reference allele. `ado` is the
 * amplification dropout rate for the affected cell's stage.
 */
export function callSiblingOrigin(
  observations: readonly AB[],
  siblings: number,
  ado: number,
  dropIn = DROP_IN,
  eps = GENOTYPE_ERROR,
  /** Whether `observations` are oriented to a consistent haplotype. Unphased input can still
   *  separate a deletion from a dropout cluster; it cannot say which side was lost. */
  phased = false,
): SiblingCall {
  const used = observations.filter((g) => g !== 'NC')
  const phi = expectedPhi(siblings, dropIn)
  const agreement = hetRule(siblings)
  const base = { markers: used.length, siblings, agreement, phi, phased }

  if (siblings < MIN_SIBLINGS) {
    return {
      ...base, hypothesis: 'refused', posterior: NaN,
      why: `${siblings} event-free sibling cells is under the ${MIN_SIBLINGS} needed to establish `
        + 'heterozygosity without using the cell being scored',
    }
  }
  if (phi > PHI_WALL) {
    return {
      ...base, hypothesis: 'refused', posterior: NaN,
      why: `expected false-heterozygous fraction ${phi.toFixed(3)} is over the ${PHI_WALL} wall, `
        + 'past which the reference-copy-lost hypothesis is unresolvable at any marker count',
    }
  }
  if (used.length < MIN_MARKERS) {
    return {
      ...base, hypothesis: 'refused', posterior: NaN,
      why: `${used.length} informative markers is under the ${MIN_MARKERS} this needs`,
    }
  }

  // Log space: a region can carry thousands of markers and the product underflows immediately.
  const logs: [number, number, number] = [0, 0, 0]
  for (const g of used) {
    const l = likelihoods(g, phi, ado, eps)
    logs[0] += Math.log(l[0])
    logs[1] += Math.log(l[1])
    logs[2] += Math.log(l[2])
  }
  const top = Math.max(...logs)
  const w = logs.map((x) => Math.exp(x - top))
  const total = w[0] + w[1] + w[2]
  const post = w.map((x) => x / total)
  const names: Hypothesis[] = ['reference-copy-lost', 'other-copy-lost', 'no-deletion']
  // Without phase the two side hypotheses are not distinguishable in principle, so they are pooled
  // and only the deletion/no-deletion contrast is reported.
  if (!phased) {
    const deletion = post[0] + post[1]
    if (Math.max(deletion, post[2]) < CALL_POSTERIOR) {
      return {
        ...base, hypothesis: 'refused', posterior: Math.max(deletion, post[2]),
        why: 'evidence separates neither deletion nor an intact region',
      }
    }
    if (deletion > post[2]) {
      return {
        ...base, hypothesis: 'refused', posterior: deletion,
        why: `a copy is missing here at posterior ${deletion.toFixed(4)}, but the observations are `
          + 'not phased, so WHICH side was lost is not determinable from this region alone. '
          + 'Naming it needs phase, and naming the parent needs one anchor per donor group',
      }
    }
    return {
      ...base, hypothesis: 'no-deletion', posterior: post[2],
      why: `both copies present at posterior ${post[2].toFixed(4)}; the region is a dropout `
        + 'cluster rather than a lost copy',
    }
  }
  let best = 0
  for (let i = 1; i < 3; i += 1) if (post[i] > post[best]) best = i

  if (post[best] < CALL_POSTERIOR) {
    return {
      ...base, hypothesis: 'refused', posterior: post[best],
      why: `best hypothesis reaches ${post[best].toFixed(3)}, under the ${CALL_POSTERIOR} needed; `
        + 'the evidence does not separate them',
    }
  }
  return {
    ...base,
    hypothesis: names[best],
    posterior: post[best],
    why: `${names[best]} at posterior ${post[best].toFixed(4)} over ${used.length} markers, `
      + `${siblings} siblings agreeing ${agreement}, false-het fraction ${phi.toFixed(3)}`,
  }
}

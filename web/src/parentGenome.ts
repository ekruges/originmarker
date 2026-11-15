/**
 * Whether an array carries the genome of the adult loaded as its parent, read before any allele is
 * read against that adult.
 *
 * Read at markers where the adult is homozygous for an allele the founder panel holds at RARE_MAX or
 * under, from the sample's calls there:
 *
 *   D = (same homozygote - other homozygote) / called. A child reads above zero and anyone else below,
 *       at any dropout that takes either allele alike; the adult's own haploid product reads near 1.
 *   E = heterozygous share there / heterozygous share where the adult is the common homozygote. A
 *       heterozygote needs both alleles, so E holds where dropout favours one parent's allele. Drop-in
 *       pulls it toward 1, so a low E says nothing.
 *
 * Not read against a reconstructed parent, whose rare homozygotes include the adult's heterozygotes.
 * Does not tell a parent from a child or a full sibling, and does not count the genomes of an array
 * that carries none of the adult's. Measurements: audit/truth/v6-2026-11-15.txt.
 */

export interface CarryTally {
  /** Adult homozygous for a rare allele: the sample the same homozygote, heterozygous, the other. */
  same: number; het: number; opp: number
  /** Adult homozygous for the common allele: the sample the same homozygote, heterozygous, the other. */
  ctlSame: number; ctlHet: number; ctlOpp: number
}
export const emptyCarryTally = (): CarryTally => ({ same: 0, het: 0, opp: 0, ctlSame: 0, ctlHet: 0, ctlOpp: 0 })

/** Founder frequency at or under which the allele the adult is homozygous for is rare. */
export const RARE_MAX = 0.2

/** One marker into the tally. `freqB` is the founders' B allele frequency. */
export function addCarry(
  t: CarryTally, parent: string | null | undefined, sample: string | null | undefined,
  freqB: number | undefined,
): void {
  if ((parent !== 'AA' && parent !== 'BB') || freqB === undefined) return
  if (sample !== 'AA' && sample !== 'BB' && sample !== 'AB') return
  const own = parent === 'BB' ? freqB : 1 - freqB
  if (own <= RARE_MAX) {
    if (sample === parent) t.same += 1
    else if (sample === 'AB') t.het += 1
    else t.opp += 1
  } else if (own >= 1 - RARE_MAX) {
    if (sample === parent) t.ctlSame += 1
    else if (sample === 'AB') t.ctlHet += 1
    else t.ctlOpp += 1
  }
}

/** Calls needed in each class before a chromosome is read. Chromosome 21 offers a few hundred. */
export const MIN_CARRY_MARKERS = 50
/** Share of rare homozygotes, where the adult is the common one, over which the calls are noise. */
export const CONTROL_MAX = 0.15
/** Standard errors from zero before a chromosome is said to carry the adult's genome, or not to. */
export const CARRY_Z = 2.5

export type CarryState = 'present' | 'absent' | 'unclear' | 'noise' | 'thin'

export interface CarryCall { state: CarryState; d: number; z: number; n: number }

export function carryCall(t: CarryTally): CarryCall {
  const n = t.same + t.het + t.opp
  const ctl = t.ctlSame + t.ctlOpp   // homozygous calls only: the control asks which homozygote
  if (n < MIN_CARRY_MARKERS || ctl < MIN_CARRY_MARKERS) return { state: 'thin', d: NaN, z: NaN, n }
  const d = (t.same - t.opp) / n
  if (t.ctlOpp / ctl > CONTROL_MAX) return { state: 'noise', d, z: NaN, n }
  const se = Math.sqrt(Math.max(1e-12, ((t.same + t.opp) / n - d * d) / n))
  const z = d / se
  return { state: z >= CARRY_Z ? 'present' : z <= -CARRY_Z ? 'absent' : 'unclear', d, z, n }
}

/** Chromosomes that must carry the adult's genome before the array is said to. */
export const PRESENT_MIN = 3
/** Enrichment at or over which the array carries the adult's genome whatever D reads. */
export const ENRICHMENT_MIN = 5
/** Heterozygous calls at the adult's rare homozygotes before the enrichment is read. */
export const MIN_ENRICHMENT_HETS = 200
/** Chromosomes that must read without it, and none with it, before the array is said to hold none.
 *  Also the chromosomes that must be read before the genome is judged as a whole. */
export const ABSENT_MIN = 10

/**
 * `carries`: the adult's genome is in the array. `none`: it is on no readable chromosome, which is
 * the other parent's pronucleus or polar body, a cell that lost every chromosome of this parent, or
 * an adult who is not a parent; the alleles cannot tell those apart. `unshown`: too few chromosomes
 * read to say either. `unread`: no founder frequencies, or a reconstructed parent.
 */
export type ParentGenomeVerdict = 'carries' | 'none' | 'unshown' | 'unread'

export interface ParentGenome {
  verdict: ParentGenomeVerdict
  present: number
  absent: number
  /** Chromosomes read: neither thin nor noise. */
  read: number
  noise: number
  /** Median D over the chromosomes read. */
  d: number
  /** Heterozygote enrichment over the chromosomes read; NaN where too few calls carry it. */
  enrichment: number
}

export function parentGenome(tallies: readonly CarryTally[]): ParentGenome {
  const all = tallies.map((t) => ({ t, c: carryCall(t) }))
  const read = all.filter(({ c }) => c.state === 'present' || c.state === 'absent' || c.state === 'unclear')
  const present = read.filter(({ c }) => c.state === 'present').length
  const absent = read.filter(({ c }) => c.state === 'absent').length
  const ds = read.map(({ c }) => c.d).sort((a, b) => a - b)
  const sum = (f: (t: CarryTally) => number): number => read.reduce((a, { t }) => a + f(t), 0)
  const het = sum((t) => t.het)
  const ctlHet = sum((t) => t.ctlHet)
  const enrichment = het >= MIN_ENRICHMENT_HETS && ctlHet > 0
    ? (het / sum((t) => t.same + t.het + t.opp)) / (ctlHet / sum((t) => t.ctlSame + t.ctlHet + t.ctlOpp))
    : NaN
  // The genome as a whole: D pooled over the chromosomes read, admitted only where none is absent.
  const n = sum((t) => t.same + t.het + t.opp)
  const homs = sum((t) => t.same + t.opp)
  const pooled = n ? (sum((t) => t.same) - sum((t) => t.opp)) / n : NaN
  const pooledZ = n ? pooled / Math.sqrt(Math.max(1e-12, (homs / n - pooled * pooled) / n)) : NaN
  const asAWhole = absent === 0 && read.length >= ABSENT_MIN && pooledZ >= CARRY_Z
  const verdict: ParentGenomeVerdict = present >= PRESENT_MIN || enrichment >= ENRICHMENT_MIN || asAWhole
    ? 'carries'
    : present === 0 && absent >= Math.max(ABSENT_MIN, read.length / 2) ? 'none' : 'unshown'
  return {
    verdict, present, absent, read: read.length,
    noise: all.filter(({ c }) => c.state === 'noise').length,
    d: ds.length ? ds[Math.floor(ds.length / 2)] : NaN,
    enrichment,
  }
}

export const UNREAD: ParentGenome = {
  verdict: 'unread', present: 0, absent: 0, read: 0, noise: 0, d: NaN, enrichment: NaN,
}

/** What to tell the reader where the adult's genome is not shown. `who` is 'father' or 'mother'. */
export function parentGenomeNotice(g: ParentGenome, who: string): { headline: string; body: string } | null {
  if (g.verdict === 'none') {
    return {
      headline: `None of the loaded ${who}'s genome is in this array`,
      body: `${g.absent} of ${g.read} chromosomes read are without it and none carries it. That is what the `
        + `other parent's pronucleus or polar body reads, what a cell that lost every chromosome of this `
        + `parent reads, and what any cell reads against an adult who is not its parent. The alleles cannot `
        + `tell these apart, so how many genomes the array holds is read without this parent.`,
    }
  }
  if (g.verdict === 'unshown') {
    return {
      headline: `The loaded ${who} is not shown to be a parent of this array`,
      body: `${g.present} of ${g.read} chromosomes read carry that genome`
        + `${g.noise ? `, and ${g.noise} could not be read because their calls are noise` : ''}. `
        + `How many genomes the array holds is read without this parent.`,
    }
  }
  return null
}

/** The same in one sentence, for the run log and the limits. */
export function parentGenomeSentence(g: ParentGenome, who: string): string {
  const n = parentGenomeNotice(g, who)
  return n ? `${n.headline}: ${n.body}` : ''
}

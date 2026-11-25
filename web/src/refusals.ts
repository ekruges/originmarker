/**
 * Why a claim was withheld, as a record rather than a sentence.
 *
 * THE RULE THIS SERVES. A claim is withheld only when the cell's own data shows the claim cannot be
 * measured: the DNA is not there, the calls are at chance, or a second genome is present where the
 * claim needs one. Each of those is a PHYSICAL cause, named in `PHYSICAL_CAUSES`, and a refusal for
 * one of them is final. Anything else that currently withholds a claim is a PROXY standing in for a
 * physical cause it does not measure, and is recorded as `proxy` with the gate that fired, so the
 * proxies can be counted, scored against published truth, and removed one at a time.
 *
 * WHY A RECORD AND NOT TEXT. The reasons this tool printed were accurate and still left a reader
 * unable to say which refusals had thrown away calls the published data confirms. A record keeps the
 * value that was withheld, so it can be scored, and the numbers the gate fired on, so the reason a
 * reader sees is always the measurement and never a paraphrase of it.
 *
 * `withheld` is kept for scoring and is never rendered. A display shows `cause` and `measured`.
 */

/** Physical failures of the data itself. The closed list; nothing outside it is a final refusal. */
export const PHYSICAL_CAUSES = [
  /** The DNA is not there, on the whole array or on this chromosome. */
  'no-template',
  /** Most of the calls are random, so they say more about the noise than about any genome. */
  'no-genotype-signal',
  /** Another person's DNA is in the reaction, where the claim needs a single genome. */
  'second-genome',
  /** Nothing in this genome can anchor the intensity, so copy number has no reference. */
  'no-intensity-reference',
  /** The family offers too few informative markers here, before this sample is read at all. */
  'parents-uninformative',
  /** The claim needs a parent array that was not supplied. */
  'parent-not-loaded',
  /** Drop-in moved the B allele frequencies, so where they band no longer says how many copies
   *  there are: two haploid pronuclei of GSE148488 read their band as diploid at drop-in of 0.30. */
  'baf-displaced',
] as const

export type PhysicalCause = (typeof PHYSICAL_CAUSES)[number]

/**
 * What a record's cause is when it is not physical. `proxy`: a gate standing in for a physical cause
 * it does not measure. `none`: a label that withholds nothing, so there is no cause to name (a
 * sample whose changes were measured against an undisturbed remainder).
 */
export type Cause = PhysicalCause | 'proxy' | 'none'

/** Each cause in the words every surface prints. `none` withholds nothing and is never printed. */
export const CAUSE_SAYS: Record<Exclude<Cause, 'none'>, string> = {
  'no-template': 'the DNA is not there',
  'no-genotype-signal': 'most of the genotype calls are random',
  'second-genome': 'a second person\'s DNA is in the reaction',
  'no-intensity-reference': 'no normal level is left to read intensity against',
  'parents-uninformative': 'the parents offer too few informative markers here',
  'parent-not-loaded': 'it needs a parent array that was not loaded',
  'baf-displaced': 'drop-in moved the allele fractions',
  proxy: 'a gate that stands in for a physical cause',
}

/** Every claim kind a refusal can withhold. Derived rows (gains, losses, dosage, Mendelian) follow
 *  from these, so withholding at the source withholds them too. */
export type ClaimKind =
  | 'aneuploidy' | 'lostParent' | 'aneuploidyParent' | 'mixed' | 'alleleGain'
  | 'segment' | 'findings' | 'triploidy' | 'complex' | 'relatedness' | 'sample'

export interface Refusal {
  claim: ClaimKind
  /** Absent for a refusal of the whole sample or a whole channel. */
  chrom?: string
  /** The code path that withheld it, so a proxy can be traced to the line that fired. */
  gate: string
  /** A physical cause, `proxy` for a gate that stands in for one without measuring it, or `none`
   *  for a label that withholds nothing. */
  cause: Cause
  /** The numbers the gate fired on. Never empty: a refusal without its measurement is a sentence. */
  measured: Record<string, number>
  /** What would have been reported. Kept for scoring, never shown. */
  withheld?: unknown
  /** True where the claim is still on the result and only labelled unreliable, rather than removed. */
  soft?: boolean
}

export const isPhysical = (c: string): c is PhysicalCause =>
  (PHYSICAL_CAUSES as readonly string[]).includes(c)

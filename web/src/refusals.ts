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
  /** The calls carry no information about any genome: hybridisation or clustering failed. */
  'no-genotype-signal',
  /** Another person's DNA is in the reaction, where the claim needs a single genome. */
  'second-genome',
  /** Nothing in this genome can anchor the intensity, so copy number has no reference. */
  'no-intensity-reference',
  /** The family offers too few informative markers here, before this sample is read at all. */
  'parents-uninformative',
  /** The claim needs a parent array that was not supplied. */
  'parent-not-loaded',
] as const

export type PhysicalCause = (typeof PHYSICAL_CAUSES)[number]

/** Every claim kind a refusal can withhold. Derived rows (gains, losses, dosage, Mendelian) follow
 *  from these, so withholding at the source withholds them too. */
export type ClaimKind =
  | 'aneuploidy' | 'lostParent' | 'aneuploidyParent' | 'mixed' | 'alleleGain'
  | 'segment' | 'findings' | 'relatedness' | 'sample'

export interface Refusal {
  claim: ClaimKind
  /** Absent for a refusal of the whole sample or a whole channel. */
  chrom?: string
  /** The code path that withheld it, so a proxy can be traced to the line that fired. */
  gate: string
  /** A physical cause, or `proxy` for a gate that stands in for one without measuring it. */
  cause: PhysicalCause | 'proxy'
  /** The numbers the gate fired on. Never empty: a refusal without its measurement is a sentence. */
  measured: Record<string, number>
  /** What would have been reported. Kept for scoring, never shown. */
  withheld?: unknown
  /** True where the claim is still on the result and only labelled unreliable, rather than removed. */
  soft?: boolean
}

export const isPhysical = (c: string): c is PhysicalCause =>
  (PHYSICAL_CAUSES as readonly string[]).includes(c)

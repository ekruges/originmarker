/**
 * What a sample withheld, and why, in the words every surface prints.
 *
 * One row for each refusal on the result (refusals.ts) and one for each row of a sample not
 * reportable (pageSays.changeList's notReported). The count of rows is the count the record keeps:
 * refusals plus rows not reported. The page and the report both print this list.
 */
import type { Refusal } from './refusals.ts'
import { NOT_REPORTED_SAYS, whereText, whoseText, type ChangeList } from './pageSays.ts'
import { OTHER_MIN_SPAN_BP } from './otherParentLoss.ts'

export const CLAIM_SAYS: Record<string, string> = {
  aneuploidy: 'whole-chromosome copy number',
  lostParent: 'whose copy a lost chromosome lacks',
  aneuploidyParent: 'the parent of a whole-chromosome change',
  mixed: 'a mixed signal on a chromosome',
  alleleGain: 'an extra copy read from the alleles',
  segment: 'a segment',
  findings: 'a finding',
  triploidy: 'triploidy',
  complex: 'a disturbed genome',
  relatedness: 'whether the loaded parent is this sample\'s parent',
  sample: 'the sample as a whole',
}

export const CAUSE_SAYS: Record<string, string> = {
  'no-template': 'the DNA is not there',
  'no-genotype-signal': 'most of the genotype calls are random',
  'second-genome': 'a second person\'s DNA is in the reaction',
  'no-intensity-reference': 'no normal level is left to read intensity against',
  'parents-uninformative': 'the parents offer too few informative markers here',
  'parent-not-loaded': 'it needs a parent array that was not loaded',
  'baf-displaced': 'drop-in moved the allele fractions',
  proxy: 'a gate that stands in for a physical cause',
}

const GATE_SAYS: Record<string, string> = {
  'normal-majority': 'most autosomes are not in one copy state',
  'drop-in.no-alleles': 'drop-in on this array, and no alleles to read the segment by',
  'stage.drop-in': 'drop-in on this array',
  'stage.no-genotype-signal': 'the stage found no genotype signal',
  'stage.second-genome': 'the stage found a second genome',
  'stage.heterozygosity-ceiling': 'heterozygosity over what one genome can carry',
  'stage.declared': 'the array did not resolve to any material',
  'copy-scan.chaotic-genome': 'the copy-number scan found the genome too disturbed',
  'other-parent.under-span': `under ${OTHER_MIN_SPAN_BP / 1e6} Mb, no second line`,
}

export interface WithheldRow { what: string; why: string; measured: string }

const num = (v: number): string => (Number.isFinite(v) ? String(+v.toPrecision(3)) : 'not measured')

/** Every refusal and every row not reported of one sample, in that order. */
export function withheldRows(
  r: { refusals?: readonly Refusal[]; integrity?: { headline: string } },
  list: Pick<ChangeList, 'notReported'>,
): WithheldRow[] {
  const refused = (r.refusals ?? []).map((x): WithheldRow => ({
    what: `${CLAIM_SAYS[x.claim] ?? x.claim}${x.chrom ? `, chr${x.chrom}` : ''}${x.soft ? ' (labelled, not removed)' : ''}`,
    why: x.claim === 'sample' && r.integrity ? r.integrity.headline
      : `${GATE_SAYS[x.gate] ?? x.gate}: ${CAUSE_SAYS[x.cause] ?? x.cause}`,
    measured: Object.entries(x.measured ?? {}).map(([k, v]) => `${k} ${num(v)}`).join(', '),
  }))
  const seen = (list.notReported?.rows ?? []).map((g): WithheldRow => ({
    what: `${whereText(g)}: ${whoseText(g)}`, why: NOT_REPORTED_SAYS, measured: '',
  }))
  return [...refused, ...seen]
}

/** The heading over the list, with its count. */
export const withheldHeading = (n: number): string => `Withheld, and why: ${n}`

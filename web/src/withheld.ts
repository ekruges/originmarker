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
  'stage.unknown': 'the stage was not read',
  'copy-scan.chaotic-genome': 'the copy-number scan found the genome too disturbed',
  'other-parent.under-span': `under ${OTHER_MIN_SPAN_BP / 1e6} Mb, no second line`,
}

const pct = (v: number): string => `${(100 * v).toFixed(1)}%`
const num = (v: number): string => String(+v.toPrecision(3))
const int = (v: number): string => Math.round(v).toLocaleString('en-US')

/** Every number a refusal records, in plain words with its unit. A key with no entry here throws. */
export const MEASURE_SAYS: Record<string, { label: string; fmt: (v: number) => string }> = {
  callRate: { label: 'call rate', fmt: pct },
  heterozygosity: { label: 'heterozygosity', fmt: pct },
  dropIn: { label: 'drop-in', fmt: pct },
  heterozygositySpread: { label: 'spread of heterozygosity between chromosomes', fmt: (v) => `${num(v)} (coefficient of variation)` },
  normalFraction: { label: 'autosomes in the commonest copy state', fmt: pct },
  deviantFraction: { label: 'share of autosomes read as changed', fmt: pct },
  affectedAutosomes: { label: 'autosomes read as changed', fmt: int },
  totalAutosomes: { label: 'autosomes read', fmt: int },
  spanBp: { label: 'span', fmt: (v) => `${(v / 1e6).toFixed(2)} Mb` },
  score: { label: 'scan score', fmt: (v) => `${num(v)} (log-likelihood ratio)` },
  markers: { label: 'markers', fmt: int },
  lrrZ: { label: 'intensity', fmt: (v) => `${num(v)} spreads from the cell's own level` },
}

const known = <T>(map: Record<string, T>, key: string, what: string): T => {
  const v = map[key]
  if (v === undefined) throw new Error(`withheld: no words for the ${what} ${key}`)
  return v
}

/** The numbers a refusal fired on, each with its label and unit. */
export const measuredText = (m: Record<string, number> | undefined): string =>
  Object.entries(m ?? {}).map(([k, v]) => {
    const say = known(MEASURE_SAYS, k, 'measurement')
    return `${say.label} ${Number.isFinite(v) ? say.fmt(v) : 'not measured'}`
  }).join(', ')

export interface WithheldRow { what: string; why: string; measured: string }

/**
 * Every refusal and every row not reported of one sample, in that order. A refusal of the whole
 * sample is said in the integrity call's own reason, which counts no change, so it cannot disagree
 * with the card's count.
 */
export function withheldRows(
  r: { refusals?: readonly Refusal[]; integrity?: { why: string } },
  list: Pick<ChangeList, 'notReported'>,
): WithheldRow[] {
  const refused = (r.refusals ?? []).map((x): WithheldRow => ({
    what: `${known(CLAIM_SAYS, x.claim, 'claim')}${x.chrom ? `, chr${x.chrom}` : ''}${x.soft ? ' (labelled, not removed)' : ''}`,
    why: x.claim === 'sample' && r.integrity ? r.integrity.why
      : `${known(GATE_SAYS, x.gate, 'gate')}: ${known(CAUSE_SAYS, x.cause, 'cause')}`,
    measured: measuredText(x.measured),
  }))
  const seen = (list.notReported?.rows ?? []).map((g): WithheldRow => ({
    what: `${whereText(g)}: ${whoseText(g)}`, why: NOT_REPORTED_SAYS, measured: '',
  }))
  return [...refused, ...seen]
}

/** The heading over the list, with its count. */
export const withheldHeading = (n: number): string => `Withheld, and why: ${n}`

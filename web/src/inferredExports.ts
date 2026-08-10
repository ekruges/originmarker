/**
 * Exports for an inferred-reference run, in forms another program can actually read.
 *
 * Three shapes, because they have three different consumers:
 *
 *   long CSV     one row per pair, for plotting and clustering. The shape R and pandas want.
 *   matrix CSV   the square matrix, for a heatmap or a supplementary table as printed.
 *   JSON         the whole run including provenance, for a pipeline or a LIMS.
 *
 * Provenance travels with every one of them. In the CSVs it is a `#` comment block, which
 * `pandas.read_csv(comment='#')` and `read.csv(comment.char='#')` both skip by default, so the
 * file stays loadable in one call while still carrying what it was made from. The JSON carries
 * the same fields as real keys.
 *
 * The reconstructed genotype itself is written by `inferredArray.ts`, which carries its own note
 * on why a file of homozygous calls belonging to an identifiable person is marked the way it is.
 */

export interface RunProvenance {
  tool: string
  generatedAt: string
  reportId: string
  experiment: string
  products: string[]
  mMin: number
  markers: number
  meanM: number
  ascertainment: number
  contamination: number
  spuriousAbsence: number
  sameParentMax: number
  differentParentMin: number
}

export interface PairRate {
  a: string
  b: string
  rate: number
  shared?: number
  opposite?: number
  sameGroup: boolean
}

export interface SampleResult {
  sample: string
  role: string
  group: string
  /** The parental-origin call, which is the answer the run exists to produce. */
  origin?: string
  callRate?: number
  hetBand?: number
  absence?: number
  ceiling?: number
  ratio?: number
  verdict?: string
  excludedBecause?: string
}

const HEADER = [
  'This reference was RECONSTRUCTED from haploid meiotic products.',
  'No array of the parent was used. Every rate here was measured against that inference',
  'and must not be compared against a run that used a real parental array: the denominators',
  'differ, because a measured array contributes its heterozygous markers and this does not.',
]

/** `#` comments so the file still loads in one call: pandas comment='#', R comment.char='#'. */
function provenanceComments(p: RunProvenance): string[] {
  return [
    ...HEADER.map((l) => `# ${l}`),
    '#',
    `# tool: ${p.tool}`,
    `# generated: ${p.generatedAt}`,
    `# report_id: ${p.reportId}`,
    `# experiment: ${p.experiment}`,
    `# products: ${p.products.join(' ')}`,
    `# agreement_rule: at least ${p.mMin} products called and agreeing at a marker`,
    `# markers: ${p.markers}`,
    `# mean_products_per_marker: ${p.meanM.toFixed(3)}`,
    `# ascertainment: ${p.ascertainment.toFixed(4)} of the parent's heterozygosity`,
    `# contamination: ${p.contamination.toFixed(6)}`,
    `# spurious_absence_haploid: ${p.spuriousAbsence.toFixed(6)}`,
    `# same_parent_below: ${p.sameParentMax}`,
    `# different_parent_at_or_above: ${p.differentParentMin}`,
    '#',
  ]
}

/** Quote only when a field would otherwise break the row, so the output stays diffable. */
const cell = (v: string | number | undefined): string => {
  if (v === undefined || v === null) return ''
  const s = String(v)
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s
}

const rowOf = (xs: (string | number | undefined)[]): string => xs.map(cell).join(',')

/** One row per pair. The shape a plotting or clustering script wants. */
export function concordanceLongCsv(pairs: PairRate[], p: RunProvenance): string {
  const cols = ['sample_a', 'sample_b', 'opposite_homozygote_rate', 'shared_markers',
    'opposite_markers', 'same_group']
  return [
    ...provenanceComments(p),
    '# One row per pair of products. opposite_homozygote_rate is the fraction of markers where',
    '# both are called homozygous and they disagree. Two haploid products of one parent differ',
    '# only where that parent is heterozygous and they drew differently.',
    '#',
    cols.join(','),
    ...pairs.map((r) => rowOf([r.a, r.b, r.rate.toFixed(6), r.shared, r.opposite,
      r.sameGroup ? 'yes' : 'no'])),
  ].join('\n') + '\n'
}

/** The square matrix as printed, for a heatmap or a supplementary table. */
export function concordanceMatrixCsv(
  names: string[], rate: (a: number, b: number) => number, p: RunProvenance,
): string {
  return [
    ...provenanceComments(p),
    '# Square matrix of the opposite-homozygote rate. The diagonal is empty. Row and column',
    '# order is by parent group, so a group appears as a contiguous block.',
    '#',
    rowOf(['', ...names]),
    ...names.map((n, i) => rowOf([n, ...names.map((_, j) =>
      (i === j ? '' : rate(i, j).toFixed(6)))])),
  ].join('\n') + '\n'
}

/** One row per array, including the ones excluded and why. */
export function sampleResultsCsv(rows: SampleResult[], p: RunProvenance): string {
  const cols = ['sample', 'origin', 'role', 'parent_group', 'call_rate', 'het_baf_band',
    'absence', 'noise_ceiling', 'ratio_to_ceiling', 'verdict', 'excluded_because']
  return [
    ...provenanceComments(p),
    '# One row per array. origin is the parental-origin call: an array carrying the reconstructed',
    '# parent\'s genome came from that parent, one decisively lacking it came from the other.',
    '# Excluded arrays are kept with the reason, so the table accounts for',
    '# every file that went in. ratio_to_ceiling is absence divided by that sample\'s own noise',
    '# ceiling; at or below 1 reads present, at or above 3 reads absent, between is uncalled.',
    '#',
    cols.join(','),
    ...rows.map((r) => rowOf([r.sample, r.origin, r.role, r.group,
      r.callRate?.toFixed(4), r.hetBand?.toFixed(4), r.absence?.toFixed(6),
      r.ceiling?.toFixed(6), r.ratio?.toFixed(4), r.verdict, r.excludedBecause])),
  ].join('\n') + '\n'
}

/** The whole run, for a pipeline that would rather not parse comments. */
export function runManifestJson(
  p: RunProvenance, groups: string[][], samples: SampleResult[], pairs: PairRate[],
  ascertainmentLadder: Record<string, number>, limits: { what: string; why: string }[],
): string {
  return `${JSON.stringify({
    reference_kind: 'inferred',
    warning: HEADER.join(' '),
    tool: p.tool,
    generated_at: p.generatedAt,
    report_id: p.reportId,
    experiment: p.experiment,
    reference: {
      products: p.products,
      m_min: p.mMin,
      markers: p.markers,
      mean_products_per_marker: p.meanM,
      ascertainment: p.ascertainment,
      ascertainment_by_m: ascertainmentLadder,
      contamination: p.contamination,
      spurious_absence_haploid: p.spuriousAbsence,
    },
    membership: {
      method: 'all pairs within a group must read same parent',
      same_parent_below: p.sameParentMax,
      different_parent_at_or_above: p.differentParentMin,
      groups: groups.map((g, i) => ({ group: i + 1, members: g })),
    },
    samples,
    pairs,
    not_reported: limits,
  }, null, 2)}\n`
}

/** Hand a blob to the browser. The only way anything leaves the page. */
export function download(text: string, filename: string, mime: string): void {
  const url = URL.createObjectURL(new Blob([text], { type: mime }))
  const a = document.createElement('a')
  a.href = url
  a.download = filename
  document.body.appendChild(a)
  a.click()
  a.remove()
  URL.revokeObjectURL(url)
}

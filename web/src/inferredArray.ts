/**
 * Write a reconstructed parental genotype out as an array file.
 *
 * The tool withheld this for a while, on the reasoning that a file of homozygous calls belonging
 * to an identifiable person could be re-imported as though it had been measured. The reasoning
 * was right about the hazard and wrong about the gain: the whole point of reconstructing a parent
 * is to have something to call parental origin AGAINST, and an experiment with no array of that
 * parent has nothing else to use.
 *
 * So it is written, and the hazard is handled where it lives. Every file opens with a banner that
 * any reader sees first, and carries `INFERRED_MARK` in a machine-readable line. `inferredMark()`
 * finds it, and the tools that ingest arrays refuse to treat a file carrying it as a measured
 * one. The header lines are comments the ingest header sniffer skips, so the file still loads
 * anywhere a real export does.
 *
 * WHAT IT IS NOT: a genotype of the parent. It is the subset of that parent's genome the products
 * could establish, every call homozygous, heterozygous sites either detected and dropped or
 * missed and asserted as homozygous. `contamination` is the fraction that is the latter.
 */
import type { AB } from './informativity.ts'

/** Present in every inferred file, on its own line, in a form no measured export produces. */
export const INFERRED_MARK = 'OriginMarker-inferred-reference'

export interface InferredArrayInput {
  /** probe id -> the homozygous call the reference asserts there. */
  genotype: Map<string, AB>
  /** Where a probe sits. Probes it cannot place are dropped rather than written at position 0. */
  locus: (id: string) => { chrom: string; pos: number } | null
  products: readonly string[]
  mMin: number
  contamination: number
  spuriousAbsence: number
  hRetained: number
  /** 'paternal', 'maternal', or '' when the products did not say which parent this is. */
  side: string
  reportId: string
  generatedAt: string
  build: string | null
}

/** True if this text is an inferred reference rather than a measured array. Cheap by design:
 *  callers pass the first few kilobytes of a file, not the whole thing. */
export const inferredMark = (head: string): boolean => head.includes(INFERRED_MARK)

/** How much of a file has to be read to answer. The banner is the first thing written, so this
 *  is generous by an order of magnitude and still costs one small read. */
export const MARK_BYTES = 4096

/** Whether a dropped file is an inferred reference. Reads only the head. */
export async function isInferredFile(f: File): Promise<boolean> {
  try {
    return inferredMark(await f.slice(0, MARK_BYTES).text())
  } catch {
    // A file that cannot be read at all fails later and more informatively than here.
    return false
  }
}

export function inferredArrayText(i: InferredArrayInput): string {
  const nl = String.fromCharCode(10)
  const head = [
    `# ${INFERRED_MARK}`,
    '#',
    '# THIS IS NOT A MEASURED ARRAY. It is a parental genotype reconstructed from the haploid',
    '# cells that parent produced, and every call in it is homozygous by construction. Treating',
    `# it as an array of a real person overstates what was measured by ${
      (i.contamination * 100).toFixed(2)}% of markers.`,
    '#',
    `# side                  ${i.side || 'not established'}`,
    `# products              ${i.products.length} (${i.products.join(', ')})`,
    `# markers               ${i.genotype.size}`,
    `# min products / marker ${i.mMin}`,
    `# contamination         ${(i.contamination * 100).toFixed(3)}%  heterozygous sites asserted `
      + 'as homozygous',
    `# spurious absence      ${(i.spuriousAbsence * 100).toFixed(3)}%  a true haploid offspring `
      + 'reads this much absence against it for that reason alone',
    `# parent heterozygosity ${(i.hRetained * 100).toFixed(3)}%  over the markers retained`,
    `# assembly              ${i.build ?? 'undetermined'}`,
    `# report                ${i.reportId}`,
    `# generated             ${i.generatedAt}`,
    '#',
    'probeset_id\tchr\tposition\tgenotype',
  ].join(nl)

  // Sorted, so the file is byte-identical between runs of the same products and a diff between
  // two references is readable.
  const rows: string[] = []
  for (const [id, gt] of i.genotype) {
    const at = i.locus(id)
    if (!at) continue
    rows.push(`${id}\t${at.chrom}\t${at.pos}\t${gt}`)
  }
  rows.sort()
  return `${head}${nl}${rows.join(nl)}${nl}`
}

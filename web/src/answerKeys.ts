/**
 * Answer keys: markers whose right answer is known before the sample is read, measured per cell.
 *
 * WHY. Every quality gate this tool had compared a number to a line chosen once, on other material:
 * a heterozygosity ceiling derived from bulk drop-in, a heterozygous fraction under which a chromosome
 * counts as one copy. On a cell whose error rates differ from the material the line was drawn on, the
 * line refuses real calls or admits false ones, and nothing in the output says which. An answer key
 * measures the error on THIS cell, at markers where the parents have already said what the right
 * answer is, so a claim can be judged against what this array actually does.
 *
 * SAME-HOMOZYGOTE MARKERS. Where both measured parents carry the same homozygote, the sample can only
 * be that homozygote: no parent, no relative and no recombination can put another allele there. So a
 * heterozygous call is DROP-IN, and the opposite homozygote is an outright ERROR, the rate at which an
 * allele absent from the template appears as a homozygous call. That second rate is what a lost
 * parent's allele should read on a chromosome that has lost that parent's copy.
 *
 * Measured by the over-refusal consult on GSE148488: drop-in 0.220 to 0.307 on the seven cells refused
 * by the heterozygosity ceiling against a 0.020 median on admitted cells, wrong homozygote 0.007 to
 * 0.039 on the same seven. Both need two MEASURED bulk parents; a reconstructed parent is homozygous at
 * markers where the real one is not, and those markers would read as errors that are not.
 */
import type { AB } from './informativity.ts'

export interface SameHomTally {
  /** Called sample genotypes at markers where both parents are the same homozygote. */
  n: number
  /** Of those, heterozygous: drop-in. */
  het: number
  /** Of those, the opposite homozygote: an allele absent from the template, called homozygous. */
  wrong: number
}

export const emptySameHom = (): SameHomTally => ({ n: 0, het: 0, wrong: 0 })

export const addSameHom = (t: SameHomTally, sample: AB | null, fa: AB | null, mo: AB | null): void => {
  if (fa !== mo || (fa !== 'AA' && fa !== 'BB')) return
  if (sample !== 'AA' && sample !== 'AB' && sample !== 'BB') return
  t.n += 1
  if (sample === 'AB') t.het += 1
  else if (sample !== fa) t.wrong += 1
}

/** Heterozygous calls where there can be none. */
export const dropIn = (t: SameHomTally): number => (t.n ? t.het / t.n : NaN)

/** Same-homozygote markers needed before drop-in is a measurement. A cell carries hundreds of
 *  thousands; under this the parents barely overlap and the rate is not used. */
export const MIN_SAME_HOM = 10_000

/**
 * Put a panel drop-in on the scale of the drop-in read between measured parents.
 *
 * The panel reads slightly HIGH: a founder panel's shared homozygote is occasionally a variant every
 * founder lacks and the sample carries, so a real heterozygote is counted as drop-in. Measured over
 * the 31 GSE148488 cells with a published karyotype, each against a panel of four or five unrelated
 * founders built without its own parents: panel = 0.932 x parents + 0.0424, r = 0.9983. Inverted
 * here, so a reading made with no parent at all is the reading two parents would have given.
 */
export const PANEL_DROPIN_SLOPE = 0.932
export const PANEL_DROPIN_INTERCEPT = 0.0424
export const onParentScale = (panelDropIn: number): number =>
  Number.isFinite(panelDropIn)
    ? Math.min(1, Math.max(0, (panelDropIn - PANEL_DROPIN_INTERCEPT) / PANEL_DROPIN_SLOPE))
    : NaN

/**
 * A founder panel as shipped: a header, then one line per marker, its probeset and the allele every
 * founder carries on both chromosomes. Any columns between are ignored, so the builder may carry the
 * chromosome for its own checks. One parser for the browser and the command line, so the two cannot
 * read the same file differently.
 */
export function parsePanel(text: string): Map<string, 'AA' | 'BB'> {
  const out = new Map<string, 'AA' | 'BB'>()
  const lines = text.split('\n')
  for (let i = 1; i < lines.length; i += 1) {
    const l = lines[i]
    if (!l) continue
    const tab = l.indexOf('\t')
    if (tab < 0) continue
    const allele = l.slice(l.lastIndexOf('\t') + 1).trim()
    if (allele === 'A') out.set(l.slice(0, tab), 'AA')
    else if (allele === 'B') out.set(l.slice(0, tab), 'BB')
  }
  return out
}

/** Where the shipped panel lives, relative to the site root and to the repository's web/public. */
export const FOUNDER_PANEL_PATH = 'panel/founders.tsv.gz'

/**
 * The shipped panel, fetched once. Inflated from the magic number rather than the filename, as the
 * examples are, because a host may already have removed the gzip encoding in transit. Null on any
 * failure: without a panel the stage behaves exactly as it did before panels existed.
 */
let panelPromise: Promise<Map<string, 'AA' | 'BB'> | null> | null = null
export function loadFounderPanel(base = ''): Promise<Map<string, 'AA' | 'BB'> | null> {
  panelPromise ??= (async () => {
    try {
      const res = await fetch(base + FOUNDER_PANEL_PATH)
      if (!res.ok) return null
      let blob = await res.blob()
      const head = new Uint8Array(await blob.slice(0, 2).arrayBuffer())
      if (head[0] === 0x1f && head[1] === 0x8b) {
        blob = await new Response(blob.stream().pipeThrough(new DecompressionStream('gzip'))).blob()
      }
      return parsePanel(await blob.text())
    } catch {
      return null
    }
  })()
  return panelPromise
}

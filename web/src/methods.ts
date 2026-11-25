/**
 * The methods paragraph for the run that is loaded, in the past tense.
 *
 * Every number in it is a bar applied to this run (bars.ts, read by name through barIn, which
 * throws for a bar the run did not apply) or a count of the run itself. A sentence is written only
 * for a step this run did. A citation is taken only from CHECKED, and only where the code reads the
 * same measurement that paper produced; a step with no such paper carries none.
 */
import { barIn, barValue, type Bar } from './bars.ts'

export interface Cite { key: string; ref: string; doi: string; inText: string; produced: string }

/** Copied from the literature table checked on PubMed and the publishers' pages, with DOIs. */
export const CHECKED: Record<'peiffer2006' | 'verdyck2022', Cite> = {
  peiffer2006: {
    key: 'peiffer2006', inText: 'Peiffer et al., 2006', doi: '10.1101/gr.5402306',
    ref: 'Peiffer DA et al. High-resolution genomic profiling of chromosomal aberrations using Infinium '
      + 'whole-genome genotyping. Genome Res 2006;16:1136-1148. doi:10.1101/gr.5402306.',
    produced: 'SNP arrays read both intensity and allelic ratio, so they detect copy changes and copy-neutral LOH',
  },
  verdyck2022: {
    key: 'verdyck2022', inText: 'Verdyck et al., 2022', doi: '10.1016/j.rbmo.2021.10.023',
    ref: 'Verdyck P et al. Analysis of parental contribution for aneuploidy detection (APCAD): a novel '
      + 'method to detect aneuploidy and mosaicism in preimplantation embryos. Reprod Biomed Online '
      + '2022;44:459-468. doi:10.1016/j.rbmo.2021.10.023.',
    produced: 'the share of DNA of one parent\'s origin per chromosome, from SNP genotypes of both parents, '
      + 'off one half flagging an anomaly and its parent',
  },
}

export interface MethodsRun {
  /** The tool and its version, as the report's footer prints it. */
  tool: string
  /** The chip every file was read as. */
  chip: string
  /** Genome builds of the files, 'undetermined' where none was found. */
  builds: readonly string[]
  samples: number
  father: 'none' | 'measured' | 'rebuilt'
  mother: boolean
  /** The largest group of arrays read as units of one embryo, 0 where none was grouped. */
  grouped: number
  /** Samples whose change list holds a sex-chromosome row (pageSays, step 'sex'). */
  sexRows: number
  /** Per sample with a profile: markers in the file, and markers called on chromosomes 1 to 22. */
  markers: readonly { file: number; autosomal: number }[]
  applied: readonly Bar[]
}

const int = (n: number): string => n.toLocaleString('en-US')
const span = (xs: readonly number[]): string => {
  const lo = Math.min(...xs); const hi = Math.max(...xs)
  return lo === hi ? int(lo) : `${int(lo)} to ${int(hi)}`
}

export interface Methods { text: string; cites: Cite[] }

export function methodsParagraph(run: MethodsRun): Methods {
  const v = (name: string, file?: string) => barValue(barIn(run.applied, name, file))
  const cites: Cite[] = []
  const cite = (c: Cite) => { if (!cites.includes(c)) cites.push(c); return c.inText }
  const parents = (run.father !== 'none' ? 1 : 0) + (run.mother ? 1 : 0)
  const build = run.builds.length === 1 && run.builds[0] !== 'undetermined'
    ? `on genome build ${run.builds[0]}, as read from the marker positions`
    : 'with the genome build not determined from the marker positions'
  const s: string[] = []
  s.push(`${run.samples === 1 ? 'One sample array was' : `${run.samples} sample arrays were`} read with `
    + `${run.tool}, ${run.samples === 1 ? '' : 'each '}as a ${run.chip} file ${build}.`)
  if (run.markers.length) {
    const one = run.markers.length === 1
    s.push(`${one ? 'The file held' : 'The files held'} ${span(run.markers.map((m) => m.file))} markers, `
      + `${span(run.markers.map((m) => m.autosomal))} of them called on chromosomes 1 to 22${one ? '' : ', by file'}.`)
  }
  s.push(parents === 0 ? 'No parental array was loaded.'
    : run.father === 'rebuilt' && !run.mother
      ? 'The father\'s genotype, rebuilt from haploid cells he produced, was loaded in place of his array; no array of the mother was loaded.'
      : run.father === 'measured' && !run.mother ? 'The father\'s array was loaded; no array of the mother was loaded.'
        : run.father === 'none' ? 'The mother\'s array was loaded; no array of the father was loaded.'
          : `Both parents' arrays were loaded${run.father === 'rebuilt' ? ', the father\'s rebuilt from haploid cells he produced' : ''}.`)
  s.push('Copy number was read on each array against its own other chromosomes, from its call rate, its '
    + `intensity and its allele fractions (log R ratio and B allele frequency; ${cite(CHECKED.peiffer2006)}). `
    + `A whole chromosome was read from its call rate against the genome's median (under ${v('CALL_COLLAPSE')} `
    + `of it, not there to genotype) and from its intensity against the array's own null (a calibrated z of `
    + `${v('Z_CHROMOSOME')}, and a shift of at least ${v('COPY_SHIFT_FLOOR')} for a loss or `
    + `${v('COPY_SHIFT_GAIN_FLOOR')} for a gain). Segments were scanned in windows of at least `
    + `${v('MIN_SEGMENT_MARKERS')}, kept where the log-likelihood ratio cleared ${v('SEGMENT_LRT')} and the `
    + `intensity moved by ${v('SEGMENT_LRR_SHIFT')} for a loss or ${v('SEGMENT_GAIN_SHIFT')} for a gain, and no `
    + `segment shorter than ${v('MIN_SEGMENT_SPAN_BP')} (a gain piece ${v('MIN_GAIN_SPAN_BP')}) was printed.`)
  s.push(`A sample was not reported where its call rate was under ${v('COMPLEX_CALL_RATE')}, where a share of `
    + `${v('COMPLEX_DEVIANT_FRACTION')} or more of its autosomes was changed, or where its stage did not `
    + 'resolve; what was seen on it was listed apart, and every claim withheld was listed with its cause.')
  s.push('The stage was read from the heterozygosity and the spread of the allele fractions, a cell under '
    + `${v('HAPLOID_MAX_HET')} heterozygous carrying one genome.`)
  s.push('The sex chromosomes were read from chrX and chrY against the array\'s own autosomes, from the chrY '
    + `call rate, intensity and allele fractions, with no Y read as absent on fewer than ${v('Y_MIN_PROBES')}.`)
  if (run.sexRows) {
    s.push(`In ${run.sexRows === 1 ? 'one sample' : `${run.sexRows} samples`} that reading entered the change list `
      + 'as a row: a loss or a gain where it states a sex chromosome missing or extra, a finding where it names '
      + 'two possibilities; such a row is not drawn on the figures, which draw chromosomes 1 to 22.')
  }
  if (parents) {
    s.push('Parent of origin was read from the loaded parent\'s genotypes: at markers where that parent is '
      + 'homozygous, an allele the parent lacks shows that the parent\'s copy is absent. A sample was called '
      + `as lacking that parent's contribution only above ${v('ABSENCE_MARGIN')} its own noise bound, its `
      + `no-call rate times its heterozygous fraction plus ${v('ABSENCE_ERROR_FLOOR')}.`)
    if (parents === 2) {
      s.push('With both parents loaded, each parent\'s share of the DNA was read at the markers where the '
        + `parents are opposite homozygotes (${cite(CHECKED.verdyck2022)}), and a chromosome reading one parent `
        + `alone at a share of ${v('SURVIVOR_MIN')} or more was called as missing the other parent's copy.`)
    } else {
      s.push('With one parent loaded, a stretch missing the other parent\'s copy was read from the same '
        + `markers and printed only over ${v('OTHER_MIN_SPAN_BP')}.`)
    }
    s.push(`Each parent named carries a grade from a calibrated probability: A at ${v('BAND_A_MIN')} or more, `
      + `B at ${v('BAND_B_MIN')}, C at ${v('BAND_C_MIN')}, and F under ${v('BAND_F_MIN')}.`)
  }
  if (run.grouped > 1) {
    s.push(`The ${run.grouped} arrays read as units of one embryo, by genotype agreement over `
      + `${v('SAME_EMBRYO_AGREEMENT')}, were compared on each stretch between the breakpoints of all of them, `
      + `breakpoints within ${v('SAME_BREAK_BP')} taken as one: a change every unit carries was read as from the `
      + 'gamete with three units or more and as consistent with it with two, and a change some units carry '
      + 'while another reads that stretch normal as arising after fertilisation.')
  }
  s.push(`The ${run.applied.length} bars applied to this run are listed with their values in the report.`)
  return { text: s.join(' '), cites }
}

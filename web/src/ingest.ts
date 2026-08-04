/**
 * ingest - what is this file, whose sample is it, and can anything be claimed from it.
 *
 * Phase 0 of OriginMarker 2.0: the drop box. Files arrive unlabelled, with no declared
 * platform, build, role or amplification regime, and this module works out what it can and
 * refuses cleanly when it cannot.
 *
 * The one structural result that shapes everything here: two parents of the same embryo are
 * BOTH exactly first-degree to it, with identical kinship 0.25 and IBD1 = 1. No autosomal
 * statistic can order them, so father-versus-mother is a SEX question, not a relatedness
 * question. And the pairwise sex test inverts under mixed chemistry: with bulk parents and an
 * amplified embryo, IBS0(father, embryo) exceeds IBS0(mother, embryo), so the mother reads as
 * LESS related and the sign follows the chemistry rather than the biology.
 *
 * Therefore: roles are DECLARED by the user, and inference is used to CONTRADICT a declaration,
 * never to make one. A contradicted declaration is a robust finding; an assigned role is not.
 *
 * Sex is called per sample from that sample's own chrX/autosome heterozygosity ratio, which is
 * dropout-robust because dropout scales numerator and denominator together and cancels.
 *
 * Self-check:  node src/ingest.check.ts
 */

import type { AB } from './informativity.ts'
import { CHROM_LEN, GAPS, type Build } from './buildref.ts'

export type { Build }

// --- the file format ---------------------------------------------------------------------

/**
 * A `.CEL.probes` row. Axiom post-ps-classification output, one file per sample, sample ID
 * carried only by the filename.
 *
 * Columns are located BY HEADER NAME rather than by position: a different export ordering is
 * a plausible variation and reading column 6 blindly would silently swap copy_number for
 * genotype.
 */
export interface ProbeRow {
  probesetId: string
  chrom: string
  pos: number
  /** null when the file left the field empty, which is a different fact from zero. */
  log2R: number | null
  baf: number | null
  copyNumber: number | null
  genotype: AB
  bestProbeset: boolean
  /** The file gave a genotype this module cannot read. Counted rather than swallowed: read as
   *  a no-call it looks like a failed measurement, which is a different fact from a dialect
   *  nothing here speaks. */
  unreadableGenotype: boolean
  /** The genotype exactly as the file wrote it, so the dialect can be named. */
  rawGenotype: string
}

const REQUIRED = ['probeset_id', 'chr', 'position', 'genotype'] as const

/**
 * Column aliases, because the same data family ships under different names.
 *
 * Both of these are real, on the same UK Biobank Axiom array, from two papers by the same lab:
 * Turocy 2026 (GSE186407) writes `log2R` in a TAB-separated `.CEL.probes.txt`, and Zuccaro 2020
 * (GSE148488) writes `normalized_intensity` in a COMMA-separated `.CEL.txt` with the columns in a
 * different order. Matching one spelling would have silently dropped the entire intensity channel
 * on the other, and phase 3 would then report copy-neutral LOH as untestable for a reason that
 * had nothing to do with the sample.
 */
const ALIASES: Record<string, readonly string[]> = {
  log2r: ['log2r', 'normalized_intensity', 'lrr', 'log_r_ratio', 'logrratio'],
  baf: ['baf', 'b_allele_freq', 'b allele freq', 'ballelefreq'],
  copy_number: ['copy_number', 'copynumber', 'cn'],
  bestprobeset: ['bestprobeset', 'best_probeset'],
  probe_classification: ['probe_classification', 'snp_classification'],
}

/** Delimiters seen in the wild for this family. Sniffed, never assumed. */
const DELIMITERS = ['\t', ','] as const

/**
 * Axiom's numeric genotype coding. 0/1/2/-1 is the "Numeric Call Codes" export option, and
 * which numeral means which homozygote is a convention, not a fact about the file.
 *
 * It is however a DETECTABLE convention, unlike the nucleotide one. See `verifyCoding`.
 */
const CODE: Record<string, AB> = {
  // Axiom's "Numeric Call Codes" export.
  '0': 'AA', '1': 'AB', '2': 'BB', '-1': 'NC',
  // AB space, which is what most exports write and what this module works in. Missing until
  // 2.2.3, so an ordinary AB-space file read as 100% no-call and was excluded on its call rate
  // with nothing saying why.
  'AA': 'AA', 'AB': 'AB', 'BA': 'AB', 'BB': 'BB',
  // No-call, spelled the several ways vendors spell it.
  'NC': 'NC', '--': 'NC', '---': 'NC', 'NoCall': 'NC', 'nocall': 'NC', '?': 'NC', '.': 'NC',
}

/** A pair of nucleotides rather than AB space: A/C/G/T, in any order, excluding AA which is
 *  also legal AB space. Detected so it can be REFUSED rather than silently read as no-call. */
const NUCLEOTIDE = /^[ACGT][ACGT]$/

export interface ColumnMap {
  [name: string]: number | string | undefined
  /** The delimiter this header was parsed with, carried so rows are split the same way. */
  __delim?: string
}

/**
 * Locate columns by name, sniffing the delimiter and resolving aliases.
 *
 * Returns null when no delimiter yields all four required columns, which is the honest answer
 * for a file this module does not understand - better than parsing one giant column and
 * reporting zero usable markers.
 */
export function headerMap(headerLine: string): ColumnMap | null {
  for (const delim of DELIMITERS) {
    const cols = headerLine.trim().split(delim).map((c) => c.trim().toLowerCase())
    if (cols.length < REQUIRED.length) continue
    const map: ColumnMap = {}
    cols.forEach((c, i) => { map[c] = i })
    // Fold each alias onto its canonical name so the rest of the module sees one spelling.
    for (const [canonical, spellings] of Object.entries(ALIASES)) {
      if (canonical in map) continue
      for (const alt of spellings) {
        if (alt in map) { map[canonical] = map[alt]; break }
      }
    }
    if (REQUIRED.every((r) => typeof map[r] === 'number')) {
      map.__delim = delim
      return map
    }
  }
  return null
}

const at = (f: string[], map: ColumnMap, key: string): string | undefined => {
  const i = map[key]
  return typeof i === 'number' ? f[i] : undefined
}

const num = (s: string | undefined): number | null => {
  if (s === undefined) return null
  const t = s.trim()
  if (t === '') return null
  const v = Number(t)
  return Number.isFinite(v) ? v : null
}

export function parseRow(line: string, map: ColumnMap): ProbeRow | null {
  const f = line.split(typeof map.__delim === 'string' ? map.__delim : '\t')
  const pos = num(at(f, map, 'position'))
  const id = (at(f, map, 'probeset_id') ?? '').trim()
  if (!id || pos === null) return null
  const best = at(f, map, 'bestprobeset')
  const raw = (at(f, map, 'genotype') ?? '').trim()
  return {
    probesetId: id,
    chrom: (at(f, map, 'chr') ?? '').trim(),
    pos,
    log2R: num(at(f, map, 'log2r')),
    baf: num(at(f, map, 'baf')),
    copyNumber: num(at(f, map, 'copy_number')),
    genotype: CODE[raw] ?? 'NC',
    unreadableGenotype: raw !== '' && CODE[raw] === undefined,
    rawGenotype: raw,
    // Absent means unfiltered rather than failed: a file with no such column has not told us
    // anything about probe quality, and treating that as "excluded" would drop every marker.
    bestProbeset: best === undefined ? true : best.trim() === '1',
  }
}

// --- per-sample profile ------------------------------------------------------------------

export interface ChromStats {
  markers: number; called: number; het: number; nocall: number
  /** Genotypes this module could not read. A nucleotide-space file lands entirely here. */
  unreadable: number
  /** Of those, ones that look like a nucleotide pair, which names the likely dialect. */
  nucleotide: number
}

export type Sex = 'male' | 'female' | 'ambiguous'

export interface SampleProfile {
  id: string
  markers: number
  called: number
  callRate: number
  /** Heterozygous fraction of CALLED markers. Dropout suppresses this and nothing raises it. */
  hetRate: number
  nocallRate: number
  byChrom: Map<string, ChromStats>
  /** chrX het rate divided by mean autosomal het rate. Dropout cancels in the ratio. */
  chrXHetRatio: number | null
  sex: Sex
  product: string
  idPrefix: string
  /** Whether BAF corroborates the numeric genotype coding, or contradicts it. */
  coding: CodingCheck
  /** Which assembly the positions belong to, from illegal-placement counting. */
  build: BuildCall
}

const AUTOSOMES = new Set(Array.from({ length: 22 }, (_, i) => String(i + 1)))

export const emptyChrom = (): ChromStats =>
  ({ markers: 0, called: 0, het: 0, nocall: 0, unreadable: 0, nucleotide: 0 })

/** Accumulate one row into a running profile. Streaming, so a 42 MB file costs one pass. */
export function accumulate(row: ProbeRow, byChrom: Map<string, ChromStats>): void {
  const c = row.chrom.replace(/^chr/i, '')
  let s = byChrom.get(c)
  if (!s) { s = emptyChrom(); byChrom.set(c, s) }
  s.markers++
  if (row.unreadableGenotype) {
    s.unreadable++
    if (NUCLEOTIDE.test(row.rawGenotype)) s.nucleotide++
  }
  if (row.genotype === 'NC') s.nocall++
  else { s.called++; if (row.genotype === 'AB') s.het++ }
}

/**
 * Sex from the sample's OWN chrX/autosome heterozygosity ratio.
 *
 * Dropout-robust by construction: it suppresses heterozygotes on chrX and the autosomes in the
 * same proportion, so the ratio survives where an absolute het rate does not. A male has one X
 * and therefore no chrX heterozygotes outside the pseudoautosomal region, which is about 7% of
 * the chromosome, so a male ratio sits near zero at any dropout rate.
 *
 * MEASURED, not assumed. On this array (UK Biobank Axiom) real known-sex bulk samples give:
 * male 0.113 and 0.115 (Zuccaro sperm donor, two replicates), female 1.060 and 1.144 (egg donors
 * A and C, GSE148488). So the female expectation is ABOUT 1.0 here, not the 0.7-0.8 predicted
 * from the X's smaller effective population size - that prediction does not survive contact with
 * this array's ascertainment. The bands below are wide enough to hold either way, and the raw
 * ratio is always reported so a borderline call is inspected rather than trusted.
 */
export function callSex(byChrom: Map<string, ChromStats>): { sex: Sex; ratio: number | null } {
  const x = byChrom.get('X')
  let autoCalled = 0, autoHet = 0
  for (const [c, s] of byChrom) if (AUTOSOMES.has(c)) { autoCalled += s.called; autoHet += s.het }
  if (!x || x.called === 0 || autoCalled === 0 || autoHet === 0) return { sex: 'ambiguous', ratio: null }
  const ratio = (x.het / x.called) / (autoHet / autoCalled)
  // A male carries heterozygotes only in the pseudoautosomal region. That is 6.8% of the chrX
  // markers on this array (178 of 2,617, GRCh37 PAR1 60001-2699520 and PAR2 154931044-155260560),
  // not the 1.99% of physical length this comment used to claim: the array is PAR-enriched. The
  // bands below rest on the measured male ratio of ~0.02-0.11, not on that figure. So his
  // ratio sits near 0.02; a female runs 0.7-0.8. The wide gap is deliberate: anything between is
  // a candidate chrX abnormality, a pooled file or contamination, each a finding rather than a
  // female with noise.
  return { sex: ratio < 0.15 ? 'male' : ratio > 0.55 ? 'female' : 'ambiguous', ratio }
}

/**
 * Product from marker count and ID prefix.
 *
 * 825,927 markers is the UK Biobank Axiom Array's delivered count (Bycroft 2018). Nothing else
 * in the Axiom catalogue is close, which is why the band works here and does NOT generalise:
 * Thermo publishes only inequalities (">900,000") for the rest of the range.
 */
export function detectProduct(n: number, idPrefix: string): string {
  if (idPrefix === 'AX-' && n >= 820_000 && n <= 826_500) return 'UK Biobank Axiom Array'
  if (idPrefix === 'AX-') return `Axiom array, product unidentified (${n} markers)`
  if (idPrefix === 'rs') return `Illumina or generic rsID panel (${n} markers)`
  if (idPrefix === 'cg' || idPrefix === 'ch') return 'REJECT: methylation array, not genotyping'
  return `unidentified platform (${n} markers, IDs like "${idPrefix}")`
}

export const idPrefixOf = (id: string): string => {
  const m = /^(AX-|Affx-|rs|exm|kgp|cg|ch|JHU_)/.exec(id)
  return m ? m[1] : id.slice(0, 2)
}

// --- the numeric coding IS detectable, unlike the nucleotide convention -------------------

export interface CodingCheck {
  /** 'ok' when BAF corroborates 0=AA/2=BB; 'inverted' when it contradicts it. */
  verdict: 'ok' | 'inverted' | 'untestable'
  meanBafHom0: number | null
  meanBafHom2: number | null
  meanBafHet: number | null
  note: string
}

/**
 * Check the numeric genotype coding against BAF.
 *
 * The NUCLEOTIDE convention (which base is called A) is unrecoverable from a file carrying no
 * nucleotides, and no detector should be attempted for it. But the numeric CODING - whether 0
 * or 2 is the B-homozygote - is a different question and BAF answers it: whichever code pairs
 * with BAF near 0 carries no B alleles.
 *
 * This matters because it is the one convention error that would silently invert results while
 * every internal consistency check still passed. Note that the informativity logic itself is
 * invariant under a consistent A/B relabel, so a uniformly inverted file is harmless; what this
 * catches is a file inverted relative to its own BAF column, which means the two channels
 * disagree and one of them is wrong.
 */
export function verifyCoding(
  sums: { hom0: number; n0: number; hom2: number; n2: number; het: number; nHet: number },
): CodingCheck {
  const m0 = sums.n0 ? sums.hom0 / sums.n0 : null
  const m2 = sums.n2 ? sums.hom2 / sums.n2 : null
  const mh = sums.nHet ? sums.het / sums.nHet : null
  if (m0 === null || m2 === null) {
    return { verdict: 'untestable', meanBafHom0: m0, meanBafHom2: m2, meanBafHet: mh,
             note: 'no BAF column, or one homozygous class absent' }
  }
  if (m0 < 0.25 && m2 > 0.75) {
    return { verdict: 'ok', meanBafHom0: m0, meanBafHom2: m2, meanBafHet: mh,
             note: 'BAF corroborates 0=AA, 2=BB' }
  }
  if (m2 < 0.25 && m0 > 0.75) {
    return { verdict: 'inverted', meanBafHom0: m0, meanBafHom2: m2, meanBafHet: mh,
             note: 'BAF says code 2 carries no B allele: the numeric coding is inverted' }
  }
  return { verdict: 'untestable', meanBafHom0: m0, meanBafHom2: m2, meanBafHet: mh,
           note: 'homozygous BAF means are not separated; BAF may be unreliable on this sample' }
}

// --- dropout, fitted per sample ----------------------------------------------------------

export interface DropoutEstimate { d: number; se: number; n: number; hetObserved: number }

/**
 * Per-sample dropout from the heterozygote deficit at father-het x mother-homozygous markers.
 *
 * At those markers Mendel gives exactly half AA and half AB, INDEPENDENT of allele frequency,
 * which is what makes this admissible client-side: it needs no reference panel and no assumed
 * population het rate. Dropout converts AB to a homozygous call without changing the
 * denominator, so observed het h = (1/2)(1 - d) and therefore
 *
 *     d_hat = 1 - 2h        se(d_hat) = 2 * sqrt(h(1-h)/n)
 *
 * The precision is load-bearing rather than decorative: a dropout mismatch of only 2
 * percentage points between two samples inverts genome-wide IBS0, so the estimator must resolve
 * dropout to better than 2 pp before any count-based relatedness statistic can be trusted. At
 * n = 10,000 informative markers se is about 0.008, comfortably inside that.
 */
export function estimateDropout(
  father: Map<string, AB>,
  mother: Map<string, AB>,
  subject: Map<string, AB>,
): DropoutEstimate | null {
  let n = 0, het = 0
  for (const [id, fa] of father) {
    if (fa !== 'AB') continue
    const mo = mother.get(id)
    if (mo !== 'AA' && mo !== 'BB') continue
    const s = subject.get(id)
    if (s === undefined || s === 'NC') continue      // a no-call is missing, not a homozygote
    n++
    if (s === 'AB') het++
  }
  if (n === 0) return null
  const h = het / n
  return {
    d: Math.min(1, Math.max(0, 1 - 2 * h)),
    se: 2 * Math.sqrt((h * (1 - h)) / n),
    n,
    hetObserved: h,
  }
}

// --- QC gates, re-derived for amplified material -------------------------------------------

export type GateVerdict = 'usable' | 'marginal' | 'exclude' | 'report_only'

export interface Gate {
  name: string
  verdict: GateVerdict
  value: number | null
  detail: string
}

/**
 * Sample gates.
 *
 * The vendor Axiom sample gate is a call rate of 97% or better, established on unamplified
 * bulk DNA. Amplified single-cell and few-cell material sits well below it, so a
 * vendor-compliant pipeline refuses the entire target cohort. The only published exclusion
 * threshold measured ON amplified material is Natesan 2014's: 75-95% call rate is the usable
 * band and below 60% the sample is excluded.
 *
 * Nothing here invents a threshold. Where no amplified-material threshold exists in the
 * literature the gate is `report_only`, which is the honest verdict rather than a guess.
 */
export function gates(p: SampleProfile): Gate[] {
  const out: Gate[] = []
  const cr = p.callRate

  out.push({
    name: 'call rate',
    value: cr,
    verdict: cr < 0.60 ? 'exclude' : cr < 0.75 ? 'marginal' : 'usable',
    detail: cr < 0.60
      ? 'below 60%: excluded (Natesan 2014, the only published threshold measured on amplified material)'
      : cr < 0.75
        ? 'between 60% and 75%: below the published usable band, interpret with care'
        : cr > 0.97
          ? 'inside the vendor bulk-DNA band'
          : 'inside the 75-95% amplified usable band; the Axiom vendor gate of 97% would reject this '
            + 'sample, and that gate was established on unamplified bulk DNA',
  })

  // The het-to-hom asymmetry key SNPs rest on is an approximation on amplified material, not a
  // theorem: erroneous heterozygous calls are common below 60% call rate.
  out.push({
    name: 'het-to-hom asymmetry valid',
    value: cr,
    verdict: cr < 0.60 ? 'exclude' : cr < 0.75 ? 'marginal' : 'usable',
    detail: cr < 0.60
      ? 'SUSPENDED: below 60% call rate erroneous heterozygous calls are common, so a het call is '
        + 'no longer robust and the key/non-key partition loses its guarantee'
      : 'holds within the usable call-rate band',
  })

  // Turocy 2026 excluded samples showing LOH along all chromosomes as abnormal fertilisation.
  // That criterion is qualitative and embryo-clustered, which an unlabelled drop box cannot
  // reproduce, so this flags a candidate rather than reproducing the filter.
  out.push({
    name: 'genome-wide LOH',
    value: p.hetRate,
    verdict: p.hetRate < 0.02 ? 'exclude' : 'report_only',
    detail: p.hetRate < 0.02
      ? 'heterozygosity near zero genome-wide: candidate abnormal fertilisation or failed genome '
        + 'unification (Turocy 2026 excluded such samples). The published criterion is qualitative '
        + 'and embryo-clustered, which a single unlabelled sample cannot reproduce.'
      : 'no genome-wide LOH signature',
  })

  // A diploid human cannot be 56% heterozygous. A real WGA'd blastomere measured exactly that at
  // 67% call rate, passed every other gate as merely marginal, and gave a sex call contradicting
  // the other blastomere of the same embryo.
  // Per CHROMOSOME, not genome-wide. A real blastomere came through with ten of its twenty-two
  // autosomes carrying no homozygous BAF population at all, 68.9% heterozygous, while the twelve
  // intact ones diluted the genome-wide figure to 32%: under both thresholds below, so nothing
  // fired. The failure is per-chromosome and the test has to be too.
  const perChrom = [...p.byChrom]
    .filter(([c, s2]) => /^(?:[1-9]|1[0-9]|2[0-2])$/.test(c) && s2.called >= 200)
    .map(([c, s2]) => [c, s2.het / s2.called] as const)
  const worst = perChrom.reduce<readonly [string, number]>(
    (a, b) => (b[1] > a[1] ? b : a), ['', 0])
  if (perChrom.length && worst[1] > 0.50) {
    out.push({
      name: 'heterozygosity plausible, per chromosome',
      value: worst[1],
      verdict: 'exclude',
      detail: `chr${worst[0]} is ${(100 * worst[1]).toFixed(1)}% heterozygous, which one diploid `
        + 'genome cannot be. The genome-wide figure can stay inside its bounds while individual '
        + 'chromosomes sit far outside them, so this is tested per chromosome: a chromosome whose '
        + 'B-allele frequencies carry no homozygous population is a clustering failure on that '
        + 'chromosome, not biology, and both the noise ceiling and the zygosity call read it as '
        + 'signal.',
    })
  }

  out.push({
    name: 'heterozygosity plausible',
    value: p.hetRate,
    verdict: p.hetRate > 0.50 ? 'exclude' : p.hetRate > 0.40 ? 'marginal' : 'usable',
    detail: p.hetRate > 0.50
      ? `heterozygosity ${(100 * p.hetRate).toFixed(1)}% is not attainable by a diploid genome: `
        + 'spurious heterozygous calls, contamination, or a chimeric sample. Excluded.'
      : p.hetRate > 0.40
        ? `heterozygosity ${(100 * p.hetRate).toFixed(1)}% is above the plausible range for a `
          + 'single diploid genome on this array; treat allelic calls with suspicion'
        : 'within the plausible range for one diploid genome',
  })

  // A file this module cannot read at all is a different finding from a file that failed to
  // call, and reporting it as a 0% call rate sends the reader to look at their chemistry.
  const unreadable = [...p.byChrom.values()].reduce((a, c) => a + c.unreadable, 0)
  const nucleotide = [...p.byChrom.values()].reduce((a, c) => a + c.nucleotide, 0)
  if (unreadable > 0) {
    const share = p.markers ? unreadable / p.markers : 0
    out.push({
      name: 'genotype format',
      value: share,
      verdict: share > 0.5 ? 'exclude' : 'marginal',
      detail: nucleotide > unreadable / 2
        ? `${(100 * share).toFixed(1)}% of genotypes are nucleotide pairs (A/C/G/T) rather than `
          + 'AB space. Which nucleotide is allele A is a per-marker convention that has to be '
          + 'resolved by pooling across every sample in the run, which this page does not do, so '
          + 'they are refused rather than assigned. Use an AB-space or numeric export, or the '
          + 'command line, which pools.'
        : `${(100 * share).toFixed(1)}% of genotypes are in a spelling this module does not `
          + 'read. AA/AB/BB, 0/1/2/-1 and the usual no-call tokens are understood.',
    })
  }

  out.push({
    name: 'numeric genotype coding',
    value: null,
    verdict: p.coding.verdict === 'inverted' ? 'exclude'
      : p.coding.verdict === 'untestable' ? 'report_only' : 'usable',
    detail: p.coding.note,
  })

  out.push({
    name: 'sex call',
    value: p.chrXHetRatio,
    verdict: p.sex === 'ambiguous' ? 'report_only' : 'usable',
    detail: p.sex === 'ambiguous'
      ? `chrX/autosome heterozygosity ratio ${p.chrXHetRatio?.toFixed(3) ?? 'n/a'} falls between `
        + 'the male band (near 0.02, pseudoautosomal only) and the female band (0.7-0.8). Three '
        + 'things land here and they are different findings: a candidate chrX copy-number '
        + 'abnormality, a sex-mixed or pooled file, or contamination. Not resolved, because '
        + 'guessing here is the one place this module could invert a declared role.'
      : `chrX/autosome heterozygosity ratio ${p.chrXHetRatio?.toFixed(3) ?? 'n/a'} is inside the `
        + `${p.sex} band. Dropout cancels in the ratio, so this survives amplified material.`,
  })

  return out
}

// --- local no-call clustering --------------------------------------------------------------

/** log of the binomial pmf, by recurrence, so nothing underflows at any n. */
function logPmfSeries(n: number, p: number): Float64Array {
  const out = new Float64Array(n + 1)
  const lp = Math.log(p), lq = Math.log1p(-p)
  out[0] = n * lq
  for (let k = 0; k < n; k++) out[k + 1] = out[k] + Math.log(n - k) - Math.log(k + 1) + lp - lq
  return out
}

/** P(X >= k) for X ~ Binomial(n, p), summed from the tail so small values keep precision. */
export function binomUpperTail(k: number, n: number, p: number): number {
  if (k <= 0) return 1
  if (k > n) return 0
  const lg = logPmfSeries(n, p)
  let s = 0
  for (let i = n; i >= k; i--) s += Math.exp(lg[i])
  return Math.min(1, s)
}

/** Smallest no-call count in a window of n whose upper tail falls below alpha. */
export function criticalNocalls(n: number, p0: number, alpha: number): number {
  const lg = logPmfSeries(n, p0)
  let s = 0
  for (let k = n; k >= 0; k--) {
    s += Math.exp(lg[k])
    if (s > alpha) return k + 1
  }
  return 0
}

export interface NocallCluster {
  startIndex: number
  endIndex: number
  nocalls: number
  windowMarkers: number
  rate: number
  p: number
}

/**
 * Scan for a LOCAL no-call excess against the sample's OWN genome-wide rate.
 *
 * This replaces an absolute ceiling, which is vacuous on amplified material: a fixed 5% limit
 * rejects every window of a sample running 8-13% genome-wide while detecting nothing. The real
 * signal was always a local deficit against that sample's baseline, because a non-random
 * no-call cluster is itself evidence of a structural event in that genome.
 *
 * Exact binomial, null p0 = the sample's own genome-wide no-call rate, Bonferroni-corrected
 * across markers. Power against a true local rate of 0.30 is about 0.54 at a 100-marker window
 * and 0.98 at 200, so 200 is the smallest window with usable power against a partial event; a
 * near-complete dropout event is certain from 20 markers.
 *
 * No published method exists for this, so it is assembled from scan-statistic primitives and
 * needs validation. It deliberately avoids Monte-Carlo inference, needing only the exact
 * binomial tail above.
 */
export function scanNocallClusters(
  nocallFlags: ArrayLike<0 | 1> | boolean[],
  p0: number,
  windowMarkers = 200,
  alpha = 0.05,
): NocallCluster[] {
  const n = nocallFlags.length
  if (n < windowMarkers || p0 <= 0 || p0 >= 1) return []
  const crit = criticalNocalls(windowMarkers, p0, alpha / Math.max(1, n))
  const flag = (i: number): number => (nocallFlags[i] ? 1 : 0)

  let count = 0
  for (let i = 0; i < windowMarkers; i++) count += flag(i)

  const hits: NocallCluster[] = []
  const push = (start: number, c: number) => {
    hits.push({
      startIndex: start,
      endIndex: start + windowMarkers - 1,
      nocalls: c,
      windowMarkers,
      rate: c / windowMarkers,
      p: binomUpperTail(c, windowMarkers, p0),
    })
  }
  if (count >= crit) push(0, count)
  for (let i = windowMarkers; i < n; i++) {
    count += flag(i) - flag(i - windowMarkers)
    if (count >= crit) {
      const start = i - windowMarkers + 1
      const last = hits[hits.length - 1]
      // Merge overlapping hits: one event should be one finding, not one per offset.
      if (last && start <= last.endIndex) {
        last.endIndex = i
        last.nocalls = Math.max(last.nocalls, count)
        last.rate = last.nocalls / windowMarkers
        last.p = binomUpperTail(last.nocalls, windowMarkers, p0)
      } else push(start, count)
    }
  }
  return hits
}

// --- genome build, from positions alone ----------------------------------------------------

export interface BuildCall {
  build: Build | null
  /** Illegal placements under each candidate: position past the chromosome end, or inside an
   *  assembly N-gap. The correct build gives exactly zero, which is what makes this decisive. */
  illegal: Record<Build, number>
  tested: number
  note: string
}

const BUILDS: Build[] = ['GRCh37', 'GRCh38']

/** Is a 1-based position inside a half-open 0-based gap interval. Binary search, so a whole
 *  array costs one pass rather than one scan per marker. */
function inGap(intervals: ReadonlyArray<readonly [number, number]>, pos1: number): boolean {
  const p = pos1 - 1                                   // 1-based position to 0-based coordinate
  let lo = 0, hi = intervals.length - 1
  while (lo <= hi) {
    const mid = (lo + hi) >> 1
    const [s, e] = intervals[mid]
    if (p < s) hi = mid - 1
    else if (p >= e) lo = mid + 1
    else return true
  }
  return false
}

/**
 * Determine the genome build from marker positions, with no rsIDs, no manifest, no chain file.
 *
 * A marker cannot sit inside an assembly N-gap or past the end of its chromosome, so under the
 * correct build the illegal-placement count is exactly zero. Under the wrong build the same
 * positions land in gaps at around 1-1.5%, so at array marker counts the two are separated by a
 * margin no noise can cross: at 825,656 markers the probability of zero illegal placements under
 * the wrong build is astronomically small.
 *
 * The asymmetry of the test is worth knowing: a NONZERO count EXCLUDES a build at any marker
 * count, because the correct-build rate is exactly zero rather than merely small. A zero count
 * only supports a build once enough markers have been examined to have expected a nonzero one.
 * Below roughly a thousand markers, therefore, only exclusions are trustworthy.
 */
export function detectBuild(rows: ReadonlyArray<{ chrom: string; pos: number }>): BuildCall {
  const sums = emptyBuildSums()
  for (const r of rows) accumulateBuild(r, sums)
  return buildVerdict(sums)
}

/** The verdict, split out so the streaming and whole-array routes cannot disagree. */
export function buildVerdict(sums: BuildSums): BuildCall {
  const illegal: Record<Build, number> = { GRCh37: sums.GRCh37, GRCh38: sums.GRCh38 }
  const tested = sums.tested
  const clean = BUILDS.filter((b) => illegal[b] === 0)

  if (tested === 0) {
    return { build: null, illegal, tested, note: 'no markers on a recognised primary chromosome' }
  }
  if (clean.length === 1) {
    const other = BUILDS.find((b) => b !== clean[0])!
    // A single clean build is only a positive call once the alternative was given a real chance
    // to look clean too. Below that, say undetermined rather than claim a coin flip.
    if (tested < 1_000 && illegal[other] === 0) {
      return { build: null, illegal, tested, note: `only ${tested} markers tested: too few to call` }
    }
    return {
      build: clean[0], illegal, tested,
      note: `${clean[0]}: 0 illegal placements of ${tested} markers, against ${illegal[other]} `
        + `(${((100 * illegal[other]) / tested).toFixed(3)}%) under ${other}`,
    }
  }
  if (clean.length === 0) {
    return {
      build: null, illegal, tested,
      note: 'both builds place markers illegally, so neither is right: a third assembly, a '
        + 'coordinate convention this does not model, or corrupt positions',
    }
  }
  return {
    build: null, illegal, tested,
    note: tested < 1_000
      ? `only ${tested} markers tested: both builds are clean, which at this count is expected `
        + 'under either. Undetermined.'
      : 'both builds are clean at a marker count where one should not be. Positions may not be '
        + 'genomic coordinates at all.',
  }
}

// --- relatedness: verification of a DECLARED role, never assignment ------------------------

export type Relationship =
  | 'same_individual'
  /** Flat per-window dispersion. Consistent with parent-offspring AND with unrelated, which
   *  this statistic cannot separate once there is dropout - see `assessRelatedness`. */
  | 'flat_dispersion_parent_offspring_or_unrelated'
  /** Variable per-window dispersion: mixed IBD states, so NOT a parent-offspring pair. */
  | 'sibling_or_more_distant'
  | 'indeterminate'

export interface RelatednessCall {
  relationship: Relationship
  /** Fraction of co-called markers that are opposite homozygotes. Dropout INFLATES this, in both
   *  classes and by different amounts when the two samples differ in chemistry, which is why it
   *  is reported but never used to separate parent-offspring from sibling. */
  ibs0Rate: number
  /** Standard deviation of the per-window IBS0 rate. This is the discriminator: a parent and
   *  child share one allele by descent in EVERY window, so the rate is uniform; siblings mix
   *  IBD0, IBD1 and IBD2 windows, so it varies. Dropout is independent across markers, so it
   *  shifts every window's mean together and leaves the dispersion intact. */
  windowedP0Sd: number | null
  concordance: number
  informativeMarkers: number
  windows: number
  note: string
}

/** READv2 refuses below this many expected mismatches, and so does this. */
const DEGREE_CALL_FLOOR = 10_000

/**
 * Maximum fitted-dropout difference between two samples before the degree call is refused.
 *
 * Deliberately low. The ingestion research measured a 2 percentage point mismatch inverting
 * genome-wide IBS0; real data here shows a 33 point mismatch inverting the windowed dispersion
 * as well, which was supposed to be the dropout-robust statistic. 0.05 keeps the comparison to
 * pairs whose chemistry is genuinely comparable.
 */
const DROPOUT_MISMATCH_CEILING = 0.05

/**
 * Per-window IBS0 dispersion above which a pair has mixed IBD states and is therefore not
 * parent-offspring. [N] - derived from simulation with linkage and random phase, not from
 * published calibration, and it needs validation on real family data. Measured separation, sd:
 *
 *     dropout      parent-offspring   full sibling   unrelated
 *     0%           0.0000             0.0490         0.0258
 *     30%          0.0198             0.0653         0.0317
 *     63%          0.0279             0.0826         0.0337
 *
 * The threshold sits below every sibling value and above every non-sibling value across that
 * range. The margin is roughly 30% at the tightest point, which is thin enough that the raw
 * dispersion is always reported so a caller can compare pairs within their own cohort rather
 * than rely on this constant.
 */
const SIBLING_DISPERSION = 0.045

/**
 * Compare two samples. Verification only.
 *
 * This can never say which of two parents is the father. Both parents of one embryo are exactly
 * first-degree to it, with identical expected kinship 0.25 and IBD1 = 1, so no autosomal
 * statistic orders them - a structural identity, not a power limitation. What it can do is
 * contradict a declaration: a pair declared parent-and-child that reads unrelated is a sample
 * swap, and that finding is robust in a way an assignment would not be.
 *
 * Markers must be supplied in genomic order for the windowed statistic to mean anything, since
 * the windows are what carry the linkage the dispersion depends on.
 */
export function assessRelatedness(
  a: Map<string, AB>,
  b: Map<string, AB>,
  orderedMarkerIds: readonly string[],
  windowMarkers = 500,
  /** Fitted dropout for each sample, when known. Supplying them is what lets this function
   *  refuse a degree call it cannot make; omitting them means the refusal cannot fire. */
  dropout?: { a: number; b: number },
): RelatednessCall {
  let coCalled = 0, agree = 0, ibs0 = 0
  const perWindow: number[] = []
  let winN = 0, winIbs0 = 0

  for (const id of orderedMarkerIds) {
    const x = a.get(id), y = b.get(id)
    if (x === undefined || y === undefined || x === 'NC' || y === 'NC') continue
    coCalled++
    if (x === y) agree++
    const opposite = (x === 'AA' && y === 'BB') || (x === 'BB' && y === 'AA')
    if (opposite) ibs0++
    winN++
    if (opposite) winIbs0++
    if (winN === windowMarkers) { perWindow.push(winIbs0 / winN); winN = 0; winIbs0 = 0 }
  }

  const concordance = coCalled ? agree / coCalled : 0
  const rate = coCalled ? ibs0 / coCalled : 0
  const mean = perWindow.length ? perWindow.reduce((s, v) => s + v, 0) / perWindow.length : 0
  const sd = perWindow.length > 1
    ? Math.sqrt(perWindow.reduce((s, v) => s + (v - mean) ** 2, 0) / (perWindow.length - 1))
    : null

  const base = {
    ibs0Rate: rate, windowedP0Sd: sd, concordance, informativeMarkers: coCalled,
    windows: perWindow.length,
  }

  if (coCalled < 1_000) {
    return { ...base, relationship: 'indeterminate', note: `only ${coCalled} co-called markers` }
  }
  // Excluded before any role logic, or a duplicate reads as a first-degree relative. 0.90, not
  // 0.99: real replicates of the same bulk DNA concord at 95.8%, while a parent-offspring pair on
  // the same data concords at 54.9%.
  if (concordance > 0.90 && rate < 0.02) {
    return { ...base, relationship: 'same_individual',
             note: `concordance ${(100 * concordance).toFixed(2)}% with essentially no opposite `
               + 'homozygotes: the same individual, a replicate, or monozygotic twins' }
  }
  // Measured on real data, and it defeats the statistic outright: a bulk parent against a WGA'd
  // embryo gives parent-offspring dispersion 0.073-0.098 while real full siblings give 0.041. The
  // ordering inverts, so the band does not merely miscalibrate, it points the wrong way.
  if (dropout && Math.abs(dropout.a - dropout.b) > DROPOUT_MISMATCH_CEILING) {
    return { ...base, relationship: 'indeterminate',
             note: `fitted dropout differs by `
               + `${Math.abs(dropout.a - dropout.b).toFixed(3)} between the two samples, above the `
               + `${DROPOUT_MISMATCH_CEILING} ceiling. Degree REFUSED: on real mixed-chemistry `
               + 'pairs the dispersion ordering inverts, so parent-offspring reads as sibling. '
               + `Statistics are reported (sd ${sd?.toFixed(4) ?? 'n/a'}, IBS0 `
               + `${(100 * rate).toFixed(2)}%) but no degree can be inferred from them.` }
  }
  if (coCalled < DEGREE_CALL_FLOOR || sd === null) {
    return { ...base, relationship: 'indeterminate',
             note: `${coCalled} co-called markers is below the ${DEGREE_CALL_FLOOR} floor for a `
               + 'degree call (READv2 refuses here too)' }
  }
  // The dispersion, NOT the genome-wide rate, which is reported for context and unsafe to
  // threshold on. Across 0%, 30% and 63% dropout, sibling dispersion stayed 2.5-3x
  // parent-offspring while the genome-wide rate moved from 3% to 22% and separated nothing:
  // dropout is independent across markers, so it lifts every window's mean and leaves the spread.
  if (sd > SIBLING_DISPERSION) {
    return { ...base, relationship: 'sibling_or_more_distant',
             note: `per-window IBS0 dispersion ${sd.toFixed(4)} exceeds ${SIBLING_DISPERSION}: `
               + 'windows are in mixed IBD states, so this is NOT a parent-offspring pair. A pair '
               + 'declared parent-and-child that lands here is a contradiction worth chasing.' }
  }
  // Flat dispersion is where parent-offspring lives, and unrelated too: neither has IBD structure
  // to vary. The mean rate separates them at zero dropout and not under it, so naming either here
  // could confirm a swapped sample instead of catching it.
  return {
    ...base, relationship: 'flat_dispersion_parent_offspring_or_unrelated',
    note: `per-window IBS0 dispersion ${sd.toFixed(4)} is flat, which is consistent with `
      + 'parent-offspring AND with unrelated: neither has IBD structure along the genome. '
      + `Genome-wide IBS0 is ${(100 * rate).toFixed(2)}%, and near 0% would indicate `
      + 'parent-offspring - but that figure is inflated by dropout in both samples, by different '
      + 'amounts when their chemistry differs, so it is reported and not thresholded. What IS '
      + 'ruled out here is a sibling or more distant relationship.',
  }
}

// --- assembling a profile ------------------------------------------------------------------

/** Illegal-placement tallies, accumulated in the same pass as everything else. */
export interface BuildSums { GRCh37: number; GRCh38: number; tested: number }

export const emptyBuildSums = (): BuildSums => ({ GRCh37: 0, GRCh38: 0, tested: 0 })

export function accumulateBuild(row: { chrom: string; pos: number }, s: BuildSums): void {
  const c = row.chrom.replace(/^chr/i, '')
  if (!(c in CHROM_LEN.GRCh37) || !(c in CHROM_LEN.GRCh38)) return
  s.tested++
  for (const b of BUILDS) {
    if (row.pos < 1 || row.pos > CHROM_LEN[b][c] || inGap(GAPS[b][c] ?? [], row.pos)) s[b]++
  }
}

export interface BafSums { hom0: number; n0: number; hom2: number; n2: number; het: number; nHet: number }

export const emptyBafSums = (): BafSums => ({ hom0: 0, n0: 0, hom2: 0, n2: 0, het: 0, nHet: 0 })

export function accumulateBaf(row: ProbeRow, s: BafSums): void {
  if (row.baf === null) return
  if (row.genotype === 'AA') { s.hom0 += row.baf; s.n0++ }
  else if (row.genotype === 'BB') { s.hom2 += row.baf; s.n2++ }
  else if (row.genotype === 'AB') { s.het += row.baf; s.nHet++ }
}

export function finishProfile(
  id: string,
  byChrom: Map<string, ChromStats>,
  bafSums: BafSums,
  firstId: string,
  buildSums: BuildSums = emptyBuildSums(),
): SampleProfile {
  let markers = 0, called = 0, het = 0, nocall = 0
  for (const s of byChrom.values()) {
    markers += s.markers; called += s.called; het += s.het; nocall += s.nocall
  }
  const { sex, ratio } = callSex(byChrom)
  const prefix = idPrefixOf(firstId)
  return {
    id,
    markers,
    called,
    callRate: markers ? called / markers : 0,
    hetRate: called ? het / called : 0,
    nocallRate: markers ? nocall / markers : 0,
    byChrom,
    chrXHetRatio: ratio,
    sex,
    product: detectProduct(markers, prefix),
    idPrefix: prefix,
    coding: verifyCoding(bafSums),
    build: buildVerdict(buildSums),
  }
}

/** Profile a whole file held as text. The browser path streams instead; the accumulators are
 *  shared so both routes produce the same profile. */
export function profileText(id: string, text: string): SampleProfile | null {
  const lines = text.split('\n')
  let map: ColumnMap | null = null
  let i = 0
  for (; i < lines.length; i++) {
    map = headerMap(lines[i])
    if (map) { i++; break }
  }
  if (!map) return null
  const byChrom = new Map<string, ChromStats>()
  const baf = emptyBafSums()
  const build = emptyBuildSums()
  let firstId = ''
  for (; i < lines.length; i++) {
    const line = lines[i]
    if (!line.trim()) continue
    const row = parseRow(line, map)
    if (!row) continue
    if (!firstId) firstId = row.probesetId
    accumulate(row, byChrom)
    accumulateBaf(row, baf)
    accumulateBuild(row, build)
  }
  return byChrom.size ? finishProfile(id, byChrom, baf, firstId, build) : null
}

/**
 * Public arrays that load with one click, so the tool can be tried without a family's data.
 *
 * Every one is a GEO release file, from two series on the same Affymetrix Axiom platform:
 * GSE148488, an Egli lab study of human zygotes, and GSE290961, a later deposit from the same
 * laboratory that carries a second family. They are subsets at every eighth marker, except the
 * three arrays of the trio that demonstrates the allele count, which are at every second: a
 * per-chromosome count needs 200 informative markers, and at an eighth no chromosome reaches it.
 * No value is rounded or altered.
 *
 * WHAT THE SUBSET COSTS, since it is not nothing. Linkage, stage inference and whole-chromosome
 * events all survive it, because each uses either the whole genome or a whole chromosome. A
 * SEGMENTAL change does not: the caller's resolution floor is 2,400 informative markers, and a
 * 12 Mb region holds about 390 at this density. So no example here shows a segmental call, and
 * that is a property of the download size rather than of the samples.
 */
export interface Example {
  /** GEO sample accession, which is also the filename stem. */
  gsm: string
  file: string
  role: 'donor' | 'oocyte' | 'sample'
  what: string
  expect: string
}

export const EXAMPLE_SERIES = 'GSE148488'
/** The second series, which supplies the donor-and-two-children family. */
export const EXAMPLE_SERIES_2 = 'GSE290961'
export const EXAMPLE_CITATION =
  'Zuccaro MV et al. Allele-specific chromosome removal after Cas9 cleavage in human embryos. '
  + 'Cell 2020;183(6):1650-1664. GEO GSE148488, Affymetrix Axiom, GRCh37. '
  + 'Second family from GEO GSE290961, same platform and laboratory, cited by accession because '
  + 'no associated publication was confirmed for it.'
export const EXAMPLE_STRIDE = 8
export const EXAMPLE_MARKERS = 103_207
/** The trio for the allele count: every second marker. */
export const EXAMPLE_TRIO_STRIDE = 2
export const EXAMPLE_TRIO_MARKERS = 412_828
export const EXAMPLE_DENSITY = `Every eighth marker, ${EXAMPLE_MARKERS.toLocaleString('en-US')} per file, `
  + `except the sperm donor, oocyte donor D and embryo C at every second, `
  + `${EXAMPLE_TRIO_MARKERS.toLocaleString('en-US')} per file. Values unaltered.`

export const EXAMPLES: Example[] = [
  {
    gsm: 'GSM4472397',
    file: 'GSM4472397_sperm_DNA_71.half.csv.gz',
    role: 'donor',
    what: 'sperm donor, bulk sperm DNA',
    expect: 'the paternal genome every sample below is measured against',
  },
  {
    gsm: 'GSM4472409',
    file: 'GSM4472409_A8_45.subset.csv.gz',
    role: 'sample',
    what: 'embryo A8, a child of both donors',
    expect: 'refused, and the refusal is the point. Both parents are there, 9.88% paternal and '
      + '9.27% maternal absence under a 10.91% ceiling, but at a 53% call rate the array does not '
      + 'resolve to any material, so nothing on it is reportable and no chromosome is read',
  },
  {
    gsm: 'GSM4472408',
    file: 'GSM4472408_donor_D.half.csv.gz',
    role: 'oocyte',
    what: 'oocyte donor D, the mother of embryo C below, loaded as the oocyte donor',
    expect: 'her array is over the diploid heterozygosity ceiling against the sperm donor, 34.5% '
      + 'heterozygous where one person reads 14 to 17%, and the run says so as a limit. Nothing is '
      + 'withheld on that account: her homozygous calls still make markers informative, and the '
      + 'allele count reads only those. Every sample here is read against both parents',
  },
  {
    gsm: 'GSM4552385',
    file: 'GSM4552385_embryoC_blastomere3.half.csv.gz',
    role: 'sample',
    what: 'embryo C, blastomere 3, a child of both donors, Cas9 injected at fertilisation',
    expect: 'the case both parents exist for. The source study\'s hand karyotype records twelve '
      + 'monosomies on this cell. At the markers where the parents are opposite homozygotes, '
      + 'eleven of them read one parent\'s allele alone: the maternal copy absent on 1, 2, 4, 7, '
      + '8, 17 and 20, the paternal on 5, 6, 10 and 19, each confirmed by the Mendelian channel. '
      + 'Chromosome 22 has too few markers at this density. The intensity reaches none of the '
      + 'eleven and flags 12 and 13 instead, where the file\'s copy number collapses to nothing '
      + 'and the alleles still read both parents',
  },
  {
    gsm: 'GSM4472407',
    file: 'GSM4472407_donor_A_47.subset.csv.gz',
    role: 'sample',
    what: 'oocyte donor A, unrelated to both donors above',
    expect: 'neither parent: 5.53% paternal absence against a 0.72% ceiling, and donor D is not '
      + 'her either',
  },
  {
    gsm: 'GSM4472415',
    file: 'GSM4472415_donor_C_70.subset.csv.gz',
    role: 'sample',
    what: 'oocyte donor C, unrelated to the sperm donor and the mother of embryo 4',
    expect: 'neither parent: 5.22% paternal absence against a 1.07% ceiling, and donor D is not '
      + 'her either. As a sample she is the negative control; as the oocyte she is embryo 4\'s mother',
  },
  {
    gsm: 'GSM4472424',
    file: 'GSM4472424_embryo4_TE.subset.csv.gz',
    role: 'sample',
    what: 'embryo 4, a trophectoderm biopsy, the sperm donor\'s child by oocyte donor C',
    expect: 'his child at a 0.66% opposite-homozygote rate with a second parental contribution at '
      + '14.7% of 74,399 informative markers, and no chromosomal change anywhere. Its mother is '
      + 'donor C, loaded here as a sample, so against donor D as the oocyte it reads unclear; make '
      + 'donor C the oocyte and it is a clean pass against both parents',
  },
  {
    gsm: 'GSM8826446',
    file: 'GSM8826446_donor_bulk.subset.csv.gz',
    role: 'sample',
    what: 'a second family: bulk DNA of a different donor, from GSE290961',
    expect: 'neither parent, as it should be: unrelated to both donors above. Load it as the sperm '
      + 'donor instead and the two samples below become its children. It is here so the set '
      + 'contains a second family rather than one donor and a row of negatives',
  },
  {
    gsm: 'GSM8826445',
    file: 'GSM8826445_child_bulk.subset.csv.gz',
    role: 'sample',
    what: 'that donor\'s child, as bulk DNA',
    expect: 'against GSE290961\'s donor: a child at a 1.11% opposite-homozygote rate with a second '
      + 'parental contribution at 9.0% of 81,217 informative markers. Read as bulk material, '
      + '98.5% called and 17.5% heterozygous, so the dropout inferred from it is 0.008',
  },
  {
    gsm: 'GSM8826436',
    file: 'GSM8826436_child_singlecell.subset.csv.gz',
    role: 'sample',
    what: 'the same donor\'s child, as a single amplified cell',
    expect: 'the same relationship read from far worse material, which is the comparison worth '
      + 'making. A child at 1.82% opposite-homozygote with a 8.7% second contribution, but 91.4% '
      + 'called and 13.6% heterozygous, so the inferred dropout is 0.193 rather than 0.008. '
      + 'Nothing was declared: the tool reads the amplification off the array itself, and every '
      + 'threshold downstream moves with it',
  },
  {
    gsm: 'GSM4472398',
    file: 'GSM4472398_sperm_DNA_79.subset.csv.gz',
    role: 'sample',
    what: 'a second array of the same sperm donor',
    expect: 'every allele traces back: 0.15% absence against a 1.45% ceiling. Reported as '
      + 'androgenetic because that is what a genome carrying only this donor\'s alleles is; '
      + 'on a bulk sperm sample read it as an identity match, not an embryo',
  },
]

/**
 * Fetch one example and hand back a File the drop zone cannot tell from a dropped one.
 *
 * Stored gzipped, and inflated from the magic number rather than from the filename: a host that
 * serves these with Content-Encoding: gzip has already inflated them by the time the bytes
 * arrive, and one that does not has not. Both are correct, so the bytes decide.
 */
export async function loadExample(ex: Example, base = 'examples/'): Promise<File> {
  const res = await fetch(base + ex.file)
  if (!res.ok) throw new Error(`${ex.file}: HTTP ${res.status}`)
  let blob = await res.blob()
  const head = new Uint8Array(await blob.slice(0, 2).arrayBuffer())
  if (head[0] === 0x1f && head[1] === 0x8b) {
    blob = await new Response(blob.stream().pipeThrough(new DecompressionStream('gzip'))).blob()
  }
  return new File([blob], ex.file.replace(/\.gz$/, ''), { type: 'text/csv' })
}

/**
 * Public arrays that load with one click, so the tool can be tried without a family's data.
 *
 * Every one is a GEO release file from Zuccaro et al. 2020 (GSE148488), an Egli lab study of
 * human zygotes. They are subsets: every eighth marker, and the two columns this tool never
 * reads are dropped. No value is rounded or altered, and the calls below are the ones the full
 * 825,657-marker files produce, so the demonstration is the real one at a tenth of the download.
 */
export interface Example {
  /** GEO sample accession, which is also the filename stem. */
  gsm: string
  file: string
  role: 'donor' | 'sample'
  what: string
  expect: string
}

export const EXAMPLE_SERIES = 'GSE148488'
export const EXAMPLE_CITATION =
  'Zuccaro MV et al. Allele-specific chromosome removal after Cas9 cleavage in human embryos. '
  + 'Cell 2020;183(6):1650-1664. GEO GSE148488, Affymetrix Axiom, GRCh37.'
export const EXAMPLE_STRIDE = 8
export const EXAMPLE_MARKERS = 103_207

export const EXAMPLES: Example[] = [
  {
    gsm: 'GSM4472397',
    file: 'GSM4472397_sperm_DNA_71.subset.csv.gz',
    role: 'donor',
    what: 'sperm donor, bulk sperm DNA',
    expect: 'the paternal genome every sample below is measured against',
  },
  {
    gsm: 'GSM4472409',
    file: 'GSM4472409_A8_45.subset.csv.gz',
    role: 'sample',
    what: 'embryo A8, this donor\'s own zygote',
    expect: 'unclear, and the refusal is the point. Its paternal absence of 9.88% sits under a '
      + '10.91% ceiling, so the sperm donor\'s genome is there, but at a 53% call rate zygosity '
      + 'is withheld and the androgenetic/biparental split cannot be made. Earlier versions '
      + 'called it biparental on a heterozygous band that a call rate this low does not support',
  },
  {
    gsm: 'GSM4472407',
    file: 'GSM4472407_donor_A_47.subset.csv.gz',
    role: 'sample',
    what: 'oocyte donor A, unrelated to the sperm donor',
    expect: 'no paternal contribution. 5.53% absence against a 0.72% ceiling, 7.7x over',
  },
  {
    gsm: 'GSM4472415',
    file: 'GSM4472415_donor_C_70.subset.csv.gz',
    role: 'sample',
    what: 'oocyte donor C, unrelated to the sperm donor',
    expect: 'no paternal contribution. 5.22% absence against a 1.07% ceiling, 4.9x over',
  },
  {
    gsm: 'GSM4472424',
    file: 'GSM4472424_embryo4_TE.subset.csv.gz',
    role: 'sample',
    what: 'embryo 4, a trophectoderm biopsy, this donor\'s own child',
    expect: 'a clean pass, and the only example here that is one. Confirmed as his child at a '
      + '0.66% opposite-homozygote rate with a second parental contribution at 14.7% of 74,399 '
      + 'informative markers, and no chromosomal change anywhere. Every other sample below is a '
      + 'refusal or a negative, which made the set read as though the tool only ever says no',
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

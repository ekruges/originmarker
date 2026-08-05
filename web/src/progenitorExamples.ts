/**
 * Five haploid products of one man, public, that load with one click.
 *
 * They are the paternal pronuclei of Zuccaro et al. 2020, the same series the rest of the tool
 * is validated against, and they are a POSITIVE case on purpose: five products of one parent
 * that pass every gate, group as one, and reconstruct. A user seeing the feature for the first
 * time should see it work before they see it refuse.
 *
 * Subsets at every sixteenth marker. Every rate the tool reports is a rate, so the subset gives
 * the same answer as the full array: reconstructing from these gives 15.72% ascertainment and
 * 2.33% contamination against 15.47% and 2.29% on the full 825,657-marker files.
 *
 * The series has only one father, so it cannot demonstrate the split that a second father
 * produces. That case is real and is documented rather than staged.
 */
export interface ProductExample {
  gsm: string
  file: string
  label: string
  what: string
}

export const PROGENITOR_CITATION =
  'Zuccaro MV et al. Allele-specific chromosome removal after Cas9 cleavage in human embryos. '
  + 'Cell 2020;183(6):1650-1664. GEO GSE148488, Affymetrix Axiom, GRCh37.'

export const PROGENITOR_STRIDE = 16
export const PROGENITOR_MARKERS = 51_604

export const PROGENITOR_EXAMPLES: ProductExample[] = [
  {
    gsm: 'GSM4774680',
    file: 'GSM4774680_pMII-1.subset.csv.gz',
    label: 'pMII-1',
    what: 'paternal pronucleus, zygote 1',
  },
  {
    gsm: 'GSM4774681',
    file: 'GSM4774681_pMII-2.subset.csv.gz',
    label: 'pMII-2',
    what: 'paternal pronucleus, zygote 2. The noisiest of the five at a 7.9% band, and still a '
      + 'usable product',
  },
  {
    gsm: 'GSM4774682',
    file: 'GSM4774682_pMII-3.subset.csv.gz',
    label: 'pMII-3',
    what: 'paternal pronucleus, zygote 3',
  },
  {
    gsm: 'GSM4774683',
    file: 'GSM4774683_pMII-4.subset.csv.gz',
    label: 'pMII-4',
    what: 'paternal pronucleus, zygote 4',
  },
  {
    gsm: 'GSM4774685',
    file: 'GSM4774685_pMII-6.subset.csv.gz',
    label: 'pMII-6',
    what: 'paternal pronucleus, zygote 6',
  },
]

/**
 * Fetch one example as a File the drop zone cannot tell from a dropped one.
 *
 * Inflated from the magic number rather than the filename: a host serving these with
 * Content-Encoding: gzip has already inflated them by the time the bytes arrive, and one that
 * does not has not. Both are correct, so the bytes decide.
 */
export async function loadProductExample(
  ex: ProductExample, base = 'examples/',
): Promise<File> {
  const res = await fetch(base + ex.file)
  if (!res.ok) throw new Error(`${ex.file}: HTTP ${res.status}`)
  let blob = await res.blob()
  const head = new Uint8Array(await blob.slice(0, 2).arrayBuffer())
  if (head[0] === 0x1f && head[1] === 0x8b) {
    blob = await new Response(blob.stream().pipeThrough(new DecompressionStream('gzip'))).blob()
  }
  return new File([blob], ex.file.replace(/\.gz$/, ''), { type: 'text/csv' })
}

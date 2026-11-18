/**
 * Public arrays that load with one click, so the tool can be tried without a family's data.
 *
 * Three runs, each one button: both parents, the father alone, and no parent. The cells are of
 * GSE148488 (Zuccaro et al. 2020) and were picked by scoring every open cell of that study's
 * karyotype table (Table S3) on this code, in all three setups, and keeping the cells whose card
 * states what the table states with nothing else added. Every card names a result; the sets hold
 * no negative control and no refusal.
 *
 * The second file labelled sperm is an unrelated man of GSE186441, the public series of adults the
 * founder panel is built from. The run reads each candidate against every cell and names the
 * father, which is the laboratory's own question when several men are on file.
 *
 * Every second marker, the same rows in every file, values unaltered (audit/make-examples.py).
 * That is the thinnest density at which each card prints what the full arrays print. At every
 * third marker the ten-change cell has lost too many chromosomes for the oocyte donor's check to
 * read her genome, and the card falls to unclear.
 */
export interface Example {
  /** GEO sample accession, which is also the filename stem. */
  gsm: string
  file: string
  role: 'donor' | 'oocyte' | 'sample'
  what: string
  /** One plain line, true of what the run prints for this file. */
  expect: string
  /** The published statement the example rests on. */
  source: string
}

export interface ExampleSet {
  id: string
  /** The button. */
  label: string
  /** What this run shows, one line. */
  line: string
  examples: Example[]
}

export const EXAMPLE_SERIES = 'GSE148488'
export const EXAMPLE_CITATION =
  'Zuccaro MV et al. Allele-specific chromosome removal after Cas9 cleavage in human embryos. '
  + 'Cell 2020;183(6):1650-1664. GEO GSE148488, Affymetrix Axiom, GRCh37. The second file labelled '
  + 'sperm is an unrelated adult of GEO GSE186441, the series the founder panel is built from.'
export const EXAMPLE_STRIDE = 2
export const EXAMPLE_MARKERS = 412_828
export const EXAMPLE_DENSITY = `Every second marker, ${EXAMPLE_MARKERS.toLocaleString('en-US')} per file. `
  + 'Values unaltered.'

const TABLE_S3 = 'Zuccaro et al. 2020, Table S3'

const FATHER: Example = {
  gsm: 'GSM4472397',
  file: 'GSM4472397_sperm_donor.half.csv.gz',
  role: 'donor',
  what: 'the sperm donor, bulk DNA',
  expect: 'the father of every cell in this run.',
  source: 'GEO GSE148488: the study\'s one sperm donor',
}

const OTHER_MAN: Example = {
  gsm: 'GSE186441',
  file: 'GSE186441_Sample_003_man.half.csv.gz',
  role: 'donor',
  what: 'an unrelated man, Sample_003 of GSE186441',
  expect: 'a second candidate father. The run reads both candidates against every cell and names '
    + 'the sperm donor; this man\'s genome is in none of them.',
  source: 'GEO GSE186441: constitutional blood DNA of unrelated adults, same chip',
}

const MOTHER: Example = {
  gsm: 'GSM4472408',
  file: 'GSM4472408_egg_donor_D.half.csv.gz',
  role: 'oocyte',
  what: 'egg donor D, bulk DNA',
  expect: 'the mother of every cell in this run.',
  source: 'GEO GSE148488: the three cells below are egg donor D\'s',
}

const NORMAL: Example = {
  gsm: 'GSM4472410',
  file: 'GSM4472410_embryoD_TE.half.csv.gz',
  role: 'sample',
  what: 'embryo D, a trophectoderm biopsy',
  expect: 'BOTH PARENTS, no change: both genomes are in the cell and every chromosome is whole.',
  source: `${TABLE_S3}: normal on all 22 autosomes`,
}

const FATHER_LOSS: Example = {
  gsm: 'GSM4472412',
  file: 'GSM4472412_embryoE_TE.half.csv.gz',
  role: 'sample',
  what: 'embryo E, a trophectoderm biopsy',
  expect: 'BOTH PARENTS, 1 change: the father\'s copy of the long arm of chromosome 6 is missing. '
    + 'Two readings agree.',
  source: `${TABLE_S3}: del 6q, paternal`,
}

const TEN_CHANGES: Example = {
  gsm: 'GSM4552364',
  file: 'GSM4552364_embryoA_biopsy1.half.csv.gz',
  role: 'sample',
  what: 'embryo A, a cleavage-stage blastomere',
  expect: 'BOTH PARENTS, 10 changes, each naming the parent whose copy is missing: five of the '
    + 'mother\'s and five of the father\'s.',
  source: `${TABLE_S3}: ten changes on this cell, each with its parent`,
}

/** The same cells against the sperm donor alone, which is the laboratory's usual case. */
const ONE_PARENT_NORMAL = 'FATHER + SECOND PARENT, no change: his genome is in the cell, a second '
  + 'parent\'s copy is there with it, and every chromosome is whole.'
const ONE_PARENT_FATHER_LOSS = 'FATHER + SECOND PARENT, 1 change: the father\'s copy of the long arm '
  + 'of chromosome 6 is missing, named from his array alone.'
const ONE_PARENT_TEN = 'FATHER + SECOND PARENT, 8 changes: five on the father\'s copy and three on '
  + 'the other parent\'s, read against his array alone.'

/** The same two cells read with nothing loaded, where the headline is the material. */
const NO_PARENT_CHANGE: Example = {
  ...FATHER_LOSS,
  expect: 'BULK, 1 change: chromosome 6 loses copy over its long arm. With no parent the position '
    + 'is read and no parent is named.',
}
const NO_PARENT_NORMAL: Example = {
  ...NORMAL,
  expect: 'TROPHECTODERM, no change: the material, XX, and every chromosome at two copies.',
}

export const EXAMPLE_SETS: ExampleSet[] = [
  {
    id: 'two',
    label: 'Both parents',
    line: 'Both parents loaded: a cell with nothing changed, a cell missing the father\'s copy of '
      + 'a chromosome arm, and a cell with ten changes each named to a parent.',
    examples: [FATHER, OTHER_MAN, MOTHER, NORMAL, FATHER_LOSS, TEN_CHANGES],
  },
  {
    id: 'father',
    label: 'One parent',
    line: 'The laboratory\'s usual case: only the sperm donor on file. The run names him among the '
      + 'candidates and reads each cell against him alone.',
    examples: [
      FATHER, OTHER_MAN,
      { ...NORMAL, expect: ONE_PARENT_NORMAL },
      { ...FATHER_LOSS, expect: ONE_PARENT_FATHER_LOSS },
      { ...TEN_CHANGES, expect: ONE_PARENT_TEN },
    ],
  },
  {
    id: 'none',
    label: 'No parent',
    line: 'No parent on file: the material, the sex chromosomes and copy number, with no parent named.',
    examples: [NO_PARENT_CHANGE, NO_PARENT_NORMAL],
  },
]

/** Every file any set loads, once. */
export const EXAMPLE_FILES: Example[] = [...new Map(
  EXAMPLE_SETS.flatMap((s) => s.examples).map((e) => [e.file, e]),
).values()]

/** The example a loaded file came from, by its name once inflated. The set that was pressed
 *  decides the line, because the no-parent set reads two of the same files. */
export const exampleFor = (set: ExampleSet | null, name: string): Example | undefined =>
  set?.examples.find((x) => x.file.replace(/\.gz$/, '') === name)

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

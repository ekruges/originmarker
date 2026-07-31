// Every reference below was resolved against the Crossref API before being written
// here - none is from recall. That is the same rule the app applies to genomic
// coordinates (see R1 in the README): a citation asserted from memory is exactly the
// kind of fluent, plausible, unverifiable claim this tool exists to avoid.
// Regenerate/re-verify with the script in deploy/verify-citations.py.

export interface Citation {
  id: string
  authors: string
  year: number
  title: string
  journal: string
  volume: string
  page: string
  doi: string | null
  url: string | null
  note?: string
}

export const CITATIONS: Record<string, Citation> = {
  karyomapping: {
    id: 'karyomapping',
    authors: 'Handyside et al.',
    year: 2009,
    title: 'Karyomapping: a universal method for genome wide analysis of genetic disease based on mapping crossovers between parental haplotypes',
    journal: 'Journal of Medical Genetics',
    volume: '47',
    page: '651-658',
    doi: '10.1136/jmg.2009.069971',
    url: 'https://doi.org/10.1136/jmg.2009.069971',
  },
  eshre_pgt_m: {
    id: 'eshre_pgt_m',
    authors: 'Carvalho et al.',
    year: 2020,
    title: 'ESHRE PGT Consortium good practice recommendations for the detection of monogenic disorders†',
    journal: 'Human Reproduction Open',
    volume: '2020',
    page: '',
    doi: '10.1093/hropen/hoaa018',
    url: 'https://doi.org/10.1093/hropen/hoaa018',
  },
  thousand_g: {
    id: 'thousand_g',
    authors: 'Auton et al.',
    year: 2015,
    title: 'A global reference for human genetic variation',
    journal: 'Nature',
    volume: '526',
    page: '68-74',
    doi: '10.1038/nature15393',
    url: 'https://doi.org/10.1038/nature15393',
  },
  gnomad_v4: {
    id: 'gnomad_v4',
    authors: 'Chen et al.',
    year: 2023,
    title: 'A genomic mutational constraint map using variation in 76,156 human genomes',
    journal: 'Nature',
    volume: '625',
    page: '92-100',
    doi: '10.1038/s41586-023-06045-0',
    url: 'https://doi.org/10.1038/s41586-023-06045-0',
  },
  clinvar: {
    id: 'clinvar',
    authors: 'Landrum et al.',
    year: 2019,
    title: 'ClinVar: improvements to accessing data',
    journal: 'Nucleic Acids Research',
    volume: '48',
    page: 'D835-D844',
    doi: '10.1093/nar/gkz972',
    url: 'https://doi.org/10.1093/nar/gkz972',
  },
  ensembl: {
    id: 'ensembl',
    authors: 'Harrison et al.',
    year: 2023,
    title: 'Ensembl 2024',
    journal: 'Nucleic Acids Research',
    volume: '52',
    page: 'D891-D899',
    doi: '10.1093/nar/gkad1049',
    url: 'https://doi.org/10.1093/nar/gkad1049',
  },
  ensembl_rest: {
    id: 'ensembl_rest',
    authors: 'Yates et al.',
    year: 2014,
    title: 'The Ensembl REST API: Ensembl Data for Any Language',
    journal: 'Bioinformatics',
    volume: '31',
    page: '143-145',
    doi: '10.1093/bioinformatics/btu613',
    url: 'https://doi.org/10.1093/bioinformatics/btu613',
  },
  dbsnp: {
    id: 'dbsnp',
    authors: 'Sherry',
    year: 2001,
    title: 'dbSNP: the NCBI database of genetic variation',
    journal: 'Nucleic Acids Research',
    volume: '29',
    page: '308-311',
    doi: '10.1093/nar/29.1.308',
    url: 'https://doi.org/10.1093/nar/29.1.308',
  },
  hgvs: {
    id: 'hgvs',
    authors: 'den Dunnen et al.',
    year: 2016,
    title: 'HGVS Recommendations for the Description of Sequence Variants: 2016 Update',
    journal: 'Human Mutation',
    volume: '37',
    page: '564-569',
    doi: '10.1002/humu.22981',
    url: 'https://doi.org/10.1002/humu.22981',
  },
  decode_map: {
    id: 'decode_map',
    authors: 'Halldorsson et al.',
    year: 2019,
    title: 'Characterizing mutagenic effects of recombination through a sequence-level genetic map',
    journal: 'Science',
    volume: '363',
    page: '',
    doi: '10.1126/science.aau1043',
    url: 'https://doi.org/10.1126/science.aau1043',
  },
  haldane: {
    id: 'haldane',
    authors: 'Haldane',
    year: 1919,
    title: 'The combination of linkage values, and the calculation of distances between the loci of linked factors',
    journal: 'Journal of Genetics',
    volume: '8',
    page: '299-309',
    doi: null,
    url: null,
    note: 'Predates DOI assignment; no persistent identifier exists for the original.',
  },
  kosambi: {
    id: 'kosambi',
    authors: 'Kosambi',
    year: 1943,
    title: 'The estimation of map distances from recombination values',
    journal: 'Annals of Eugenics',
    volume: '12',
    page: '172-175',
    doi: '10.1111/j.1469-1809.1943.tb02321.x',
    url: 'https://doi.org/10.1111/j.1469-1809.1943.tb02321.x',
  },
  ldlink: {
    id: 'ldlink',
    authors: 'Machiela & Chanock',
    year: 2015,
    title: 'LDlink: a web-based application for exploring population-specific haplotype structure and linking correlated alleles of possible functional variants',
    journal: 'Bioinformatics',
    volume: '31',
    page: '3555-3557',
    doi: '10.1093/bioinformatics/btv402',
    url: 'https://doi.org/10.1093/bioinformatics/btv402',
  },
  ado: {
    id: 'ado',
    authors: 'Thornhill et al.',
    year: 2001,
    title: 'A comparison of different lysis buffers to assess allele dropout from single cells for preimplantation genetic diagnosis',
    journal: 'Prenatal Diagnosis',
    volume: '21',
    page: '490-497',
    doi: '10.1002/pd.109',
    url: 'https://doi.org/10.1002/pd.109',
  },
  abcc8_chi: {
    id: 'abcc8_chi',
    authors: 'Arnoux et al.',
    year: 2011,
    title: 'Congenital hyperinsulinism: current trends in diagnosis and therapy',
    journal: 'Orphanet Journal of Rare Diseases',
    volume: '6',
    page: '',
    doi: '10.1186/1750-1172-6-63',
    url: 'https://doi.org/10.1186/1750-1172-6-63',
  },
  ma_2017: {
    id: 'ma_2017',
    authors: 'Ma et al.',
    year: 2017,
    title: 'Correction of a pathogenic gene mutation in human embryos',
    journal: 'Nature',
    volume: '548',
    page: '413-419',
    doi: '10.1038/nature23305',
    url: 'https://doi.org/10.1038/nature23305',
    note: 'Reported inter-homologue repair of a paternal mutation in human embryos. The claim rests on which parental allele the embryo carries, which is the question Syngamy measures.',
  },
  egli_2018: {
    id: 'egli_2018',
    authors: 'Egli et al.',
    year: 2018,
    title: 'Inter-homologue repair in fertilized human eggs?',
    journal: 'Nature',
    volume: '560',
    page: 'E5-E7',
    doi: '10.1038/s41586-018-0379-5',
    url: 'https://doi.org/10.1038/s41586-018-0379-5',
    note: 'Argues the apparent correction is equally explained by loss of the paternal allele, so an assay that cannot separate "corrected" from "absent" cannot settle it. This is the failure mode Syngamy is built to detect.',
  },
  zuccaro_2020: {
    id: 'zuccaro_2020',
    authors: 'Zuccaro et al.',
    year: 2020,
    title: 'Allele-Specific Chromosome Removal after Cas9 Cleavage in Human Embryos',
    journal: 'Cell',
    volume: '183',
    page: '1650-1664.e15',
    doi: '10.1016/j.cell.2020.10.025',
    url: 'https://doi.org/10.1016/j.cell.2020.10.025',
    note: 'Cas9 cleavage in human embryos frequently removes the cut chromosome rather than repairing it, in segments up to a whole chromosome arm. GEO GSE148488 supplies the bundled example arrays.',
  },
  kothiyal_2019: {
    id: 'kothiyal_2019',
    authors: 'Kothiyal et al.',
    year: 2019,
    title: 'Mendelian Inconsistent Signatures from 1314 Ancestrally Diverse Family Trios',
    journal: 'Journal of Computational Biology',
    volume: '26',
    page: '405-419',
    doi: '10.1089/cmb.2018.0253',
    url: 'https://doi.org/10.1089/cmb.2018.0253',
    note: 'Mendelian inconsistencies at 0.63% of variant sites across 1,314 trios. Used as a lower bound on the spurious-violation rate, never as the rate itself.',
  },
  natesan_2014: {
    id: 'natesan_2014',
    authors: 'Natesan et al.',
    year: 2014,
    title: 'Genome-wide karyomapping accurately identifies the inheritance of single-gene defects in human preimplantation embryos',
    journal: 'Genetics in Medicine',
    volume: '16',
    page: '838-845',
    doi: '10.1038/gim.2014.45',
    url: 'https://doi.org/10.1038/gim.2014.45',
    note: 'The call-rate band this tool gates on, and the only published threshold measured on amplified embryo biopsy material rather than bulk DNA.',
  },
}

/** Vancouver-ish one-liner: 'Handyside et al. Karyomapping... J Med Genet. 2009;47:651-658.' */
export function formatCitation(c: Citation): string {
  const vol = c.volume ? `;${c.volume}` : ''
  const pg = c.page ? `:${c.page}` : ''
  // 'Handyside et al.' already ends in a period; blindly appending another gives
  // 'Handyside et al..', which is precisely the kind of tell that makes a reference
  // list look generated rather than typeset.
  const authors = c.authors.endsWith('.') ? c.authors : `${c.authors}.`
  const title = c.title.endsWith('.') ? c.title : `${c.title}.`
  return `${authors} ${title} ${c.journal}. ${c.year}${vol}${pg}.`
}

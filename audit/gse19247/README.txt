GSE19247 - EXTERNAL VALIDATION ON CONFIRMED TRISOMIES
=====================================================

The audit's largest gap was that the gain path had only ever run on constructed events, because
every dosage event in the laboratory series is a loss. This closes that for DETECTION, on all 330
karyotype-confirmed trisomy 21 single cells in a public series: AUC 0.9779 against 114 confirmed
euploid controls on GPL8855, 0.9991 against 27 on GPL6985, with chr21 ranking first of 22
autosomes on both inside the trisomy cells themselves. Results and caveats are in
../BREAKPOINTS-AUDIT.txt section E. Do not quote numbers from here without reading E.

RUN IT

  python3 -m venv .venv && ./.venv/bin/pip install numpy      # numpy is the only dependency
  ./.venv/bin/python cluster.py                               # ~40s, writes clusters_gpl6985.npz
  ./.venv/bin/python trisomy.py                               # GPL6985, 83 trisomy cells
  ./.venv/bin/python trisomy8855.py                           # GPL8855, 247 trisomy, 114 euploid
  ./.venv/bin/python obligate.py                              # the negative result, section E2

  idat.py     Minimal Illumina IDAT reader. Self-checks in __main__: run it on any .idat.gz.
  cluster.py  Per-marker genotype cluster positions, estimated from the 255 arrays in the series.
  geno.py     Genotypes and BAF from those clusters, plus three concordance checks.
  trisomy.py  The dosage channel on confirmed trisomies, with every autosome scored as control.
  obligate.py The one-parent obligate-het boundaries, measured. They do not hold here; see E2.
  trisomy8855.py  The same dosage test on the second platform, which is where the euploid
            controls live. Intensity only, never genotypes: that platform's GEO table lists one
            Address per marker and does not mark Infinium I, which would break a genotype caller
            but leaves a total intensity valid whichever chemistry the marker uses.

GETTING THE DATA. All public, and no Illumina file is involved.

  manifest  GEO's own platform record carries the address-to-marker map, which is the part that
            would otherwise need Illumina's proprietary .bpm:
            https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6985&targ=self&form=text&view=data
            Save as gpl6985.txt beside these scripts. 373,397 markers; 342,946 are Infinium II
            and only those are used, since Infinium I needs a per-marker channel the table lacks.
  arrays    ftp.ncbi.nlm.nih.gov/geo/samples/GSM477nnn/<GSM>/suppl/ , Grn and Red .idat.gz.
            255 arrays are on GPL6985, about 1.1 GB. Put them in ./idat/ .
            367 arrays are on GPL8855, about 1.4 GB. Put them in ./idat8/ , with that platform's
            table saved as gpl8855.txt.
  index     gpl6985_samples.json, one record per GPL6985 sample: gsm, title, src, characteristics
            and its two IDAT filenames. Built from the series metadata (targ=gsm&form=text).

  LABELS ARE IN THE GEO source_name FIELD. "lymphoblast from Coriell family 1990" is a confirmed
  trisomy 21 cell. That is the ground truth these scripts use, and it matches the paper's counts.

TWO NON-OBVIOUS PROPERTIES OF THIS DATA

  1. The channels must be normalised before any angle is taken. Green and red differ in brightness
     and background by a different amount on every array; on one sample here the 98th percentiles
     are 32,170 and 20,813. Left raw, the mode of theta lands in the MIDDLE of the range, which
     reads as a 50% heterozygous sample and is impossible in a real person. cluster.py does this
     per channel per sample.
  2. There is NO BULK DNA in this series. All 622 arrays are single cells: the four buccal arrays
     are single buccal cells, not parental gDNA. Anything needing a reference parent is therefore
     limited here, and section E2 measures how much.

WHAT THIS SET CANNOT DO. Name the parent of a gain. That needs a genotyped parent beside a
confirmed gain, and this series has neither bulk parents nor a trio. Detection is validated; origin
is not, and a named parent stays provisional until a series with both is found.

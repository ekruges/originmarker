<!-- Laboratory cells are pseudonymised per audit/local/barcode-map.tsv (S12, S13, S20, S23).
     GSM accessions are public GEO identifiers and are given as they are. -->

# Consult: one reconstructed parent must give what two measured parents give

Written 11 November 2026 against version 5.51.0.

## The requirement

The laboratory never has both parents arrayed. It has one, and that one is RECONSTRUCTED from the
parent's haploid products (pronuclei, polar bodies) by `inferredReference.ts`, homozygous at every
marker by construction. The developer's requirement, in his words: the results "need to be 1:1 with
and without" both parents. This consult is about where they are not, and whether they can be.

## What is already 1:1 (5.51.0)

Measured over the 31 GSE148488 cells with a published Zuccaro 2020 Table S3 karyotype, scored once
against both measured parents and once against the father reconstructed from five of his pronuclei
alone (`audit/refusal-census.ts`, `OM_FATHER`):

- the stage, its refusal and the drop-in it states: identical on 31 of 31, because they are read
  against a founder panel (`web/public/panel/founders.tsv.gz`) and never against a parent
- whole-chromosome copy number: identical, 19 calls each
- the founder-panel answer keys: identical, they never see a parent

## What is not

| claim | two measured parents | one reconstructed father |
|---|---|---|
| which parent lost a whole chromosome (`lostParent`) | 87 (82 confirmed, 1 wrong) | 0: the channel requires two bulk parents |
| which parent supplied an extra copy (`alleleGain`) | 484 readings | 0: same requirement |
| segments | 54 | 115, of which 40 shared; most of the difference is absence segments against different marker sets |

Plus one dependence pinned in `scoreSample.check.ts`: the copy-number scan skips any chromosome the
alleles name as missing a parent, and those names exist only with two parents.

## What a reconstruction can and cannot carry

The reconstructed father here: 457,506 markers, heterozygosity 15.47% established from the products,
2.286% of its asserted homozygotes estimated to be heterozygous sites mistaken for homozygous, and a
spurious absence rate of 1.143% for a haploid true offspring. His heterozygous sites are not usable as
"this parent is heterozygous here" markers, so any test that needs them cannot run against him.

Separately, the reconstruction itself refuses half this father's genuine pronuclei. Of 8, haploid by
micromanipulation: 4 stage haploid (heterozygosity 0.076 to 0.095), 1 reads blastomere (0.121), 3 were
refused as not a genome (0.250, 0.330, 0.430). The rebuild needs 5. On every one of the 8 the panel
drop-in equals the heterozygosity to within 0.011, because a haploid has no real heterozygote: all of
its heterozygosity is drop-in.

## Questions

**Q1. Both directions of a whole-chromosome loss from one reconstructed parent P.** At markers where P is
homozygous, P's allele present or absent says whether P's copy is there. Intensity says one copy or
two. So "one copy, P present" is the other parent's copy lost, and "one copy, P absent" is P's. Is that
complete and sound on amplified single cells, with intensity unreliable on chaotic genomes (AUC 0.67 in
the over-refusal consult) and drop-in of up to 0.3? What should replace intensity where it fails?

**Q2. A normal-majority test from the alleles with one reconstructed parent.** 5.51.0 withholds copy
number on drop-in arrays because no normal majority can be shown, and that is what still hides 11 of 26
hand-called gains on the laboratory's cells. The two-parent version (the share of autosomes the
alleles read as disomic) separated chaotic from non-chaotic cells at 0.000 to 0.727 against 0.909 to
1.000. What is the one-parent analogue, and does it need P's heterozygous sites?

**Q3. Accepting a haploid product under heavy drop-in.** Heterozygosity minus panel drop-in is near
zero on haploid products (-0.011 to +0.010 here) and +0.024 to +0.045 on the diploid single cells with
drop-in. A two-equation reading, genome-wide heterozygosity d + (1-d)t and panel heterozygosity
d + (1-d)t*rho with rho the panel's residual, recovers d exactly on haploids and to within about 0.005
of the two-parent reading on diploids in a rough check. Is this a sound ploidy test for a product the
operator declares haploid, and what bar would you fix before scoring it?

**Q4. The copy-number scan's dependence on allele names.** Should the scan skip chromosomes the alleles
name as missing a parent? Not skipping makes it parent-independent and reports copy number and
parental composition as two channels on one chromosome; skipping avoids describing one event twice.

**Q5. The panel.** Four public founders give a residual offset of about 0.065 on the laboratory's
population and 0.042 on GSE148488 with four or five founders. What would a better public panel be
built from, given the platform is UK Biobank Axiom and public bulk arrays on it are few?

## Material

- 31 GSE148488 cells with published Table S3 karyotypes and both parents; 15 pronuclei; the
  reconstructed father (`/root/recon-father.probes` on tesla, five pronuclei)
- founder panels: the shipped public one, and four built without each family's parents
- the laboratory's 93 arrays, no parents, of which 26 are refused today and 4 now staged with drop-in
- two unrelated laboratory bulk arrays (S20, S23)

## Constraints

- One reconstructed parent is the target configuration; two measured parents are the reference.
- The 72 laboratory arrays with no hand call stay the specificity arm; the 477 published-normal
  autosomes stay the published arm; zero tolerance for a wrong parent.
- Every change pre-registered with its bars fixed before it is scored.

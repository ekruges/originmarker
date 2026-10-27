# Consult: rebuilding the missing parent from diploid siblings, and the circularity in it

Fourth consult in this thread. Public data for every figure. No laboratory sample is named.

This is a method question with a circularity at its centre, and I want the circularity addressed
before any code. The engineering is small; I do not think the science is.

## 1. Why this is the binding problem for the laboratory

Every validated figure this tool has is on public data, because no experiment in the laboratory's
own set has both parents arrayed.

| experiment | samples | father | mother |
|---|---|---|---|
| Experiment 2 | 18 | sperm array measured | **none** |
| Experiment 1 | 9 | **never arrayed** | egg donor measured |
| Experiment 5 | 18 | 3 products, 5 needed | rebuilt from 6 of her own products |
| MZ control | 15 | measured | measured, and it is public Zuccaro data |

With one parent the tool can exclude that parent and little else: measured at 0.854 correct about
the loaded parent's own copy and **0.289** about the other parent's. Not matching the mother is not
the same as matching the father.

So the single highest-value thing available is to recover the missing parent from the material that
exists. Experiment 2 is the candidate: a measured father and 18 diploid embryos.

## 2. What exists, and why it does not apply

`inferredReference.ts` reconstructs a parent and is validated: a father rebuilt from his own haploid
products read 0.03 percent against his measured array, and a held-out product of his read 0.71
percent against the rebuild against 0.63 measured. `MIN_PRODUCTS = 5`, and below that the method
inverts so that true offspring read as decisively absent.

It takes **haploid meiotic products only**. A haploid pronucleus gives one allele per marker
cleanly, and a heterozygous site in the parent is recognised by seeing both alleles across products.
Experiment 2 has no pronuclei. It has diploid embryos.

## 3. The proposed extension, and the arithmetic that makes it look easy

At a marker where the **father is homozygous**, say AA, a diploid embryo's genotype names the
maternal allele directly: AB means she transmitted B, AA means she transmitted A. So each embryo
contributes one maternal allele per father-homozygous marker, which is formally what one haploid
product contributes. Eighteen embryos would be eighteen products, well past the floor of five.

Father-homozygous markers are the majority of the array, so the marker budget is not the issue.

## 4. The circularity, which the haploid case does not have

An embryo that has **lost the maternal copy** of a chromosome carries only the father's allele
there. At a father-AA marker it reads AA, and the method records that the mother transmitted A. She
may be BB. So an aneuploid embryo does not merely fail to contribute; it contributes the father's
allele **as though it were hers**, and it does so across the whole affected chromosome.

Which means the rebuild is corrupted exactly on the chromosomes carrying the events the study is
about, and in the direction that makes those events harder to see afterwards: the rebuilt mother
drifts toward the father on precisely the chromosomes where embryos lost her copy.

And the circularity closes: deciding which embryos to exclude requires knowing which lost a
maternal copy, which is the question the rebuilt mother is being built to answer.

A haploid product has none of this. It is one meiotic product by dissection, and a chromosome it
lacks is absent rather than substituted.

## 5. A second asymmetry, smaller but real

A haploid product reading AA at a father-homozygous marker means that parent transmitted A. A
diploid embryo reading AA at a father-AA marker means the mother transmitted A **or** her B dropped
out, because the father's A is present either way and masks it. So per marker the diploid case
carries strictly less information than the haploid case, and `MIN_PRODUCTS = 5` was calibrated on
the haploid case. Whatever the diploid equivalent is, it is larger, and I do not know by how much.

Measured dropout on this material, for scale: amplified single cells read heterozygosity three to
six times below their own unamplified control, with dropout of 0.30 to 0.85 per allele depending on
chemistry.

## 6. What I am asking

### 6.1 Is the circularity breakable, and how?

The routes I can see, none of which I trust:

1. **Iterate.** Rebuild from all embryos, call the losses, exclude the affected chromosome of each
   affected embryo, rebuild again. I do not know whether that converges or whether it converges on
   the wrong answer, and the first pass is corrupted in a direction that would suppress the very
   events the second pass needs to find.
2. **Per chromosome rather than per embryo.** Rebuild each chromosome from only the embryos that
   read biparental there, using a statistic that needs no mother. The obligate-heterozygote test
   needs both parents, so the available one is the sample's own heterozygosity on the
   father-homozygous set, which is confounded by dropout at the rates above.
3. **Use the majority.** At a father-AA marker, if most embryos read AB, she carries B. An embryo
   that lost her copy is outvoted. This is the one I find most plausible and I can see it failing
   where many embryos share an event, which on sibling embryos of one cohort is exactly what
   happens: in the public corpus all four `del12p (pat)` rows are one embryo's blastomeres, and 14
   of 31 cells are siblings of 4 embryos.

Is there a formulation that does not have the problem rather than managing it?

### 6.2 What is the diploid equivalent of `MIN_PRODUCTS = 5`?

Derived rather than searched, please, given the information asymmetry in section 5 and the measured
dropout. The existing floor has a stated mechanism behind it and I would want the same.

### 6.3 What is the right validation, and the right bar?

This can be validated on public data before it touches the laboratory's material, which is the part
I am confident about. GSE148488 carries embryos with **both** parents arrayed, so the mother can be
rebuilt from the embryos using only the father and then scored against her actual array.

What I do not know is what bar means success. The haploid precedent is 0.03 percent against the
measured parent, but that was from clean pronuclei. And the leave-one-out discipline is mandatory
here for the reason this project already measured: a sample that helped build a reference reads 0.00
percent against it where a held-out sample reads about 1.06.

### 6.4 Is rebuilding the parent even the right move?

The alternative is to stop trying to manufacture the second parent and instead get more out of the
one that exists. The measured figures are 0.854 about the loaded parent's own copy and 0.289 about
the other's. If the second number can be raised directly, that serves the laboratory without any
reconstruction and without this circularity. I do not know whether 0.289 is a calibration problem or
an identifiability limit, and that distinction decides which of these two is worth building.

## 7. Context you may want, from the rounds before this

Settled with pre-registered measurements, so not worth your effort:

- A true single-parent loss reads absence at 0.1153 to 0.1244 against the parent whose copy went,
  computed from parental genotypes with no spectrum assumed.
- An intact chromosome on this material reads 0.0210 at the median, ten times the null floor.
- The segment scan's null is a median over the sample's other autosomes and no replacement helps:
  taking it over only the hand-normal chromosomes still leaves 8 of 20 cells unable to detect a real
  loss.
- Amplification lowers heterozygosity rather than raising it, median spurious rate 0.033 against a
  bulk 0.0525, measured on 20 single cells of a known diploid genome against unamplified controls.
- Attribution, where an event is detected, is 15 of 15 against the published karyotype. Detection is
  15 of 22 deletions and 0 of 6 duplications, measured positionally.

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

---

# What the measurements said, added after the response

Run on both corpora. The response's order was followed: step 1 first, then the floor, then whether
Experiment 2 can be rebuilt at all.

## Step 3 came back decisive, and it ends the plan for Experiment 2

The shipped `om reconstruct` refuses, and the refusal is the finding. Run on the seven haploid
products of the series that carry a genome which is not the measured parent's:

```
group 1 (4 products)   positions 03, 06, 10, 14
group 2 (2 products)   positions 07, 18
group 3 (1 product)    position 01
```

`groupByParent` finds three parents where the pairs file records one, which is the review's "at
least three families are present", confirmed independently by the shipped grouping rather than by
all-pairs rates. The largest group holds 4 products against a floor of 5.

So Experiment 2 cannot have its egg donor rebuilt. Not because of the circularity, which the
response solves, and not because of the floor, which the response derives. Because the material is
three donors and no one of them has five products.

**And Experiment 1 is in the same position.** Its four paternal products group together and the
fifth sample sits alone, which also confirms this project's earlier finding that that sample matches
nothing in the set:

```
group 1 (4 products)   and one product in a group of its own
```

Both experiments are **exactly one haploid product short**. That is the laboratory request, and it
is far more specific than an array: one more product from the donor of a named group.

## Step 1 did not reproduce, except where it did

The marker-class claim fails on the public corpus, where the truth is the authors' published
karyotype over 1,063 chromosome-observations rather than 175 with self-assigned labels. Both classes
overlap. The excluded class does separate better, a ratio of medians of 2.25 against 1.45, which is
the predicted direction, but neither gives a boundary.

The figures on the laboratory series reproduce exactly, including 111 of 154. The reading differs:
those seven samples read heterozygosity of 0.040 to 0.064, so they are single pronuclei. A haploid
product has one genome, so calling its chromosomes uniparental is correct and 72 percent is the
right answer. On the excluded class the same chromosomes give 90 of 154, which is worse.

## Step 2 reproduced to four decimals and is fixed

312 parent-heterozygous markers, 0.5514 and 0.7147, clearing the absolute floor. Exactly as
reported.

The denominator took two attempts. A fraction of the sample's median autosome refuses chr22 of a
public cell for being chr22, so a count compared across chromosomes of different sizes cannot work.
The chromosome's own call rate separates it with nothing between 0.175 and 0.813, and
`CALL_COLLAPSE = 0.60` already existed for the same question, so the fix adds no constant.

## Step 4 is not built, and the premise for it is gone

The response said to build the het-only rebuild only if step 1 leaves a gap. Step 1 leaves a gap:
the marker class is not the answer to 0.289. But step 3 says there is nothing in the laboratory's
material to rebuild from, so the rebuild would have no input even if it were built. The
one-sided formulation is the right design and it is recorded rather than implemented.

## What of the response survives in full

- The one-sided reformulation, asserting heterozygosity only. It is the correct answer to the
  circularity and would be what gets built the moment a group reaches five products.
- The demonstration that majority vote fails on Mendel alone, 35 to 48 percent of embryos, with no
  events involved. That retires the route this project was going to take.
- The informative floor, which was real and is fixed.
- The warning that this class answers "is there a second contribution" and not "whose". Nothing
  here contradicts it.

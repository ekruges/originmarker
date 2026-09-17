# Stage, abnormality, parent and mechanism across three public embryo series

What the tool reports over every publicly deposited human embryo array it can read, grouped the way
a batch is usually asked about: how many samples at each developmental stage, how many chromosome
abnormalities at each, which parent each one came from, and whether the change was carried by the
gamete or arose after fertilisation.

Run with `audit/embryo/summary.ts`. Stage and embryo identity come from each series' own GEO
metadata, read by `audit/embryo/geo_sheet.py`; they are the submitters' statements, not anything
this tool inferred. No series states either in a structured field, so each is parsed separately.

| series | arrays | material |
|---|---|---|
| GSE148488 | 135 | pronuclei, cleavage blastomeres, trophectoderm biopsies, both parents arrayed |
| GSE186407 | 166 | cleavage blastomeres, polar bodies, cell fragments, no parents |
| GSE290961 | 158 | cleavage blastomeres, whole embryos, no parents |

445 sample arrays, 355 of which resolve to a stage and are scored. The rest are refused: the
material did not produce an array this tool will read, which is itself the result for that array.

## Samples and abnormalities by stage

| stage | arrays | scored | with any change | rate | whole chromosomes | segments |
|---|---|---|---|---|---|---|
| pronucleus | 15 | 13 | 9 | 0.69 | 3 | 24 |
| whole zygote | 6 | 5 | 2 | 0.40 | 1 | 1 |
| cleavage blastomere | 320 | 265 | 155 | 0.58 | 226 | 377 |
| cell fragment | 7 | 2 | 0 | 0.00 | 0 | 0 |
| polar body | 14 | 14 | 7 | 0.50 | 5 | 2 |
| trophectoderm | 52 | 25 | 11 | 0.44 | 14 | 15 |
| whole embryo | 6 | 6 | 2 | 0.33 | 1 | 2 |
| stem cell line | 25 | 25 | 11 | 0.44 | 0 | 11 |

Whole chromosomes and segments are counted apart because their error rates are not the same thing.
The whole-chromosome call is measured against stated karyotypes on bulk DNA (96.5% of stated gains,
no false call on 58 stated-normal arrays; see `docs/validation.md`). A segment call on an amplified
single cell has no such corpus behind it, and the segment counts here should be read as a rate of
reported events rather than a rate of real ones.

## Meiotic or post-zygotic

An error carried by the egg or the sperm is in every cell of the embryo. One that happens at a
division after fertilisation is not. Neither is visible on a single array, so this reads two or
more independently sampled units of the SAME embryo (`web/src/timing.ts`).

| | whole chromosomes | segments |
|---|---|---|
| reciprocal: one unit gains what another loses | 55 | 12 |
| present in every unit (meiotic) | 5 | 17 |
| present in some units and not others | 187 | 393 |
| one unit only, nothing to compare | 3 | 10 |

**The reciprocal row is the decisive one.** A gain in one blastomere against a loss in its sister
over the same interval is one chromosome split between two daughter cells, which can only have
happened at a division after fertilisation. A change carried by the gamete is in every cell and
cannot read that way. 55 whole-chromosome events of this corpus are reciprocal.

The "some units and not others" row is not the same kind of evidence. It is consistent with a
post-zygotic error, but it is equally consistent with the event being missed in the other units,
and on amplified single cells that happens. Read it as an upper bound, and the reciprocal row as
the count that stands on its own.

## Which parent

Naming a parent needs both parents genotyped. Of these three series only GSE148488 arrayed them, so
the parental columns are filled for it and left blank for the other two rather than guessed.

| | whole chromosomes |
|---|---|
| paternal | 7 |
| maternal | 2 |
| not called | 25 |

A loss is where this question is usually answerable: the surviving alleles say whose copy went. On a
gain it needs both parents and enough markers where they are homozygous for different alleles, which
is why most events here are not called rather than called wrongly.

## Chromosome 21

11 events on chromosome 21 across the corpus, all losses: 5 whole chromosome and 6 segmental. No
chromosome 21 gain is reported on any array of any of the three series.

## What this does not cover

- **The lab's own 884 arrays are not in it.** They carry no sample sheet stating stage, embryo or
  parents, and the stage inference alone cannot supply embryo identity: two cells of one embryo look
  exactly like two cells of two embryos unless something outside the array says otherwise.
- **Parent of origin is answerable on 69 arrays**, the ones with both parents on the same panel.
- **Mechanism is answerable only where an embryo contributed two or more units**: of the 93 embryos
  represented among the scored arrays, 74 contributed two or more, covering 328 arrays.

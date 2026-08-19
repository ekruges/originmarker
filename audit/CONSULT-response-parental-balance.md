# Consult response: parental balance, the zygosity channel, and lab fit

Everything below is measured against your repo at
/Users/ezrakruger/claudecodegeneralworkspace/originmarker and the audit tables in it. Where I could
not measure something I say so. Two defects are hard blockers and I lead with them.

---

## Summary of rulings

| # | Question | Ruling |
|---|---|---|
| Q1a | Is the uniparental inference sound? | Sound in principle, BROKEN in implementation. The gate is one-sided: it fired on 7/7 gynogenetic and 0/5 androgenetic samples in your own motivating run. |
| Q1b | Is an uncalibrated bounded number right? | Emit no number. The function is a reparameterisation of fold, and it cannot leave band D below ~13,000x fold. It is a verdict wearing a number's clothes. |
| Q1c | Is row text enough to stop 1,272 reading as 1,272? | No. Emit ONE genome-level row, not per-event rows. |
| Q2a | What is the right unit? | Not regions. Use merged maximal intervals for the rate and per-genome binary as the headline. But see the LOH defect: the 183 findings are mostly detector artefact, not slicing. |
| Q2b | What denominator? | Not informative markers. It varies 1.8x while artefact propensity varies 20x, and normalising by it AMPLIFIES the confounder 1.2x. |
| Q2c | Is label permutation right? | The permutation is right; the STATISTIC is wrong. Pooled rate is not exchangeable-safe. Use rank of per-sample rate. |
| Q2d | MAX_POWER_SKEW and MIN_PER_GROUP? | 1.5 and 3 are both wrong. Measured: skew 1.1 (pooled) or 1.4 (rank); MIN_PER_GROUP 5 as a floor, 20 for a usable test. |
| Q2e | Is gyno vs andro the right comparison? | It answers a different question than the lab asked. See Q2e below. |

---

## Q1a. The inference is sound. The implementation is one-sided.

The logic is correct and I have no objection to it: a gynogenetic genome has no paternal complement,
so a change in it is maternal by construction, and that is Mendelian rather than metric. It is the
same reasoning as the obligate-het channel and it correctly needs no detection floor.

The implementation does not do this. `uniparentalOrigin()` gates on

    fold = genomeRate / explainable;  if (!(fold > 1)) return null

where `genomeRate` is absence of the LOADED parent's contribution. That quantity is large when the
loaded parent is ABSENT and small when the loaded parent is PRESENT. In a run that loads the sperm
donor, a gynogenetic sample has high absence and an androgenetic sample has near-zero absence. So
the gate can only ever fire on the class that lacks the loaded parent.

Measured on your own motivating run (audit/results.json, GSE148488, sperm donor loaded):

| class | n | fold range | gate fires |
|---|---|---|---|
| gynogenetic (maternal pronuclei) | 7 | 8.64 to 17.44 | 7 of 7 |
| androgenetic (paternal pronuclei) | 5 | 0.47 to 0.62 | 0 of 5 |
| unclear | 3 | 0.44 to 0.49 | 0 of 3, correctly |

Every androgenetic sample has fold < 1 because the loaded parent IS the parent it carries. The
channel returns null for all of them. This is not a marginal case: it is half your samples, and it
is the half the channel was written to rescue.

The consequence for `parentalBalance` is severe and it runs the same direction: if origin rows are
what feed the aggregation, maternal genomes get 7 samples' worth of called events and paternal
genomes get zero, and the "differential" verdict is guaranteed regardless of biology.

Fix: for an androgenetic sample with the paternal array loaded, the evidence is not absence, it is
PRESENCE at full dosage plus uniparental zygosity. The gate needs to branch on whether the
genome-level class matches the loaded role:

    genomeIs === role   -> margin comes from the PRESENCE call (parent_genome_present) and the
                           zygosity call, not from genomeRate/explainable
    genomeIs !== role   -> margin is genomeRate/explainable, as now

Both branches need their own margin definition. Until that exists, the channel should refuse the
matching-role case explicitly rather than returning null, so the silence is visible as a refusal
instead of reading as "not evaluable".

Second point on Q1a: your three named mechanisms are the right ones and the literature section
below covers them. But the mechanism that most threatens THIS run is not biological. The channel
inherits from `originClass`, which is assigned in parentage.ts from `zygosity`, which is read from
`hetBand > HET_BAND_DIPLOID` (0.08). All 15 pronuclei have hetBand 0.016 to 0.268, and three of
them sit ABOVE 0.08 yet still land in the uniparental branch through the absence route. A
misassigned genome-level class propagates to every event in the sample with no per-event evidence
able to contradict it. That is the real single point of failure: one wrong genome call is 1,272
wrong rows, and nothing downstream can catch it.

## Q1b. Emit no number.

Emit the verdict and the fold. Drop the confidence.

`reportedConfidence(log(fold), bound=0.206, halfEvidence=1.0)` is monotone in fold and nothing
else. It contains no information that `fold` does not already contain; it is a strictly increasing
reparameterisation of a quantity you already report. A reader who has the fold learns nothing from
the confidence, and a reader who has only the confidence has lost the interpretable quantity.

Worse, its range makes it misleading. Measured:

| fold | confidence | band |
|---|---|---|
| 1.5 | 0.466 | D |
| 5 | 0.617 | D |
| 12.9 | 0.664 | D |
| 17.4 (best in your run) | 0.675 | D |
| 1,000 | 0.736 | D |
| infinity | 0.794 | C |

Every value the channel can ever emit is band D until fold exceeds about 13,000x. Your entire
observed range, 8.6x to 17.4x, compresses into 0.648 to 0.675. So the number varies by 0.03 across
the whole plausible input range while claiming three decimal places, and it labels a Mendelian
certainty with the same band the dosage channel uses for its weakest measured guesses, where band D
accuracy is 0.60 to 0.64.

That last point is the serious one. The bands are shared across channels precisely so a reader can
compare rows. A gynogenetic inference is not 60% accurate; conditional on the genome-level call it
is essentially certain. Putting it in band D says the opposite of what is true, and the band is
what a reader will actually look at.

I checked the bound: 1 - 0.05^(1/13) = 0.2058, so SYSTEMATIC_ERROR_BOUND = 0.206 is a correct rule
of three on 13 units with zero events. That is a sound bound on the GENOME-LEVEL call. It is not a
bound on the inheritance step, which has no error of its own: given the class, the parent follows
deductively. You are applying a bound to the wrong link in the chain.

What to emit instead:

    parent: maternal
    basis:  genome-level call, gynogenetic at 11.3x the explainable ceiling
    scope:  applies to every change in this sample; not measured per event
    caveat: rests entirely on the genome-level call; if that is wrong, all of these are wrong

That is more informative than 0.659 and it cannot be quoted as a per-event accuracy.

On the lab's requirement that every abnormality carry a confidence-scored parental call: the
requirement is satisfiable only where a per-event measurement exists. Where it does not, inventing
a number to satisfy a reporting schema is how uncalibrated numbers enter clinical conversation.
Push back on the requirement rather than on the arithmetic.

## Q1c. Row text is not enough. Emit one row.

Repeating a caveat 1,272 times does not stop the count reading as 1,272. Readers count rows;
they do not read repeated boilerplate. A per-row caveat that appears on every row carries zero
information by construction, and both PDF and viewer let a reader sort or filter the table in ways
that strip the text from the number.

The structural fix is to stop emitting per-event parental rows for this channel entirely:

- one genome-level statement: "this sample is gynogenetic at 11.3x; all changes in it are maternal"
- per-event rows keep their location and class, and their ORIGIN column reads "inherited from
  genome-level call" as a categorical value, not a confidence
- the aggregation counts the SAMPLE once, never the rows

If you keep any per-event confidence at all, the count of independent parental determinations in
this run is 12, not 1,272 - one per uniparental sample.

---

## Q2a. The unit. Ruling: the framing is not wrong, but your named threat is not the real one.

You asked whether regions, breakpoints, megabases, per-genome binary, or merged intervals is the
right unit, and whether the rate framing is wrong. Measured on the 214 cells in
audit/features/regions-*.tsv, 1,513 regions total:

| unit | total | per cell median | per cell max | spearman vs regions |
|---|---|---|---|---|
| regions (shipped) | 1,513 | 4 | 59 | 1.000 |
| merged maximal intervals (1 Mb join) | 1,242 | 3 | 52 | 0.966 |
| breakpoints | 2,484 | 6 | 104 | 0.966 |
| affected Mb | 56,885 | 132 | 991 | 0.721 |
| chromosomes touched | - | 3 | 21 | 0.933 |
| any change (binary) | 214 | - | - | undefined |

Merging removes only 17.9% of rows and the ranking of samples is essentially unchanged (Spearman
0.966). Your worst cell goes from 59 regions to 52 merged intervals, not to 1. So on THIS corpus,
window-slicing is a real but modest inflation, not the dominant effect. The one unit that reorders
samples materially is affected megabases (Spearman 0.72), because it weights one 245 Mb whole-
chromosome loss above twelve small segments.

That said, the 183-finding sample you describe is a copy-neutral LOH case, and those are not in
these tables (`kind` is loss/gain only). For cnn-LOH the slicing problem is worse than merging can
fix, and it is a detector defect rather than a unit choice. `detectLoh` computes depletion against
the array's OWN mean heterozygosity:

    background = sum(het) / sum(called)   over the whole array
    depletion  = 1 - (het_window/called_window) / background
    call if depletion >= 0.65

On a uniparental homozygous genome the background is near zero, so the ratio is unstable and pure
counting noise trips the threshold. Simulated on an event-free genome, 800 windows of 2,000 called
markers, no event anywhere:

| true heterozygosity | false cnn-LOH windows |
|---|---|
| 0.030 | 0 |
| 0.010 | 8 |
| 0.005 | 42 |
| 0.002 | 178 |
| 0.001 | 307 |

Your uniparental samples sit at hetBand 0.016 to 0.036 genome-wide, and residual het inside a
homozygous genome is lower still. 178 false windows at het 0.002 is within noise of the 183 you
observed. I cannot prove your 183 are artefact without the array, but the detector reproduces that
count from nothing, and that is enough to stop the number being used as an effect size.

So the ruling on Q2a has two parts:

1. cnn-LOH must not enter the aggregation on uniparental genomes at all. A genome that is
   homozygous by construction has no heterozygosity to lose; "copy-neutral LOH" is not a
   well-defined event there. Gate `detectLoh` on `zygosity !== 'uniparental_homozygous'`, or
   require an absolute het floor (background > ~0.05) rather than a relative depletion. This
   single change probably removes most of your 183 and much of the 1,272.

2. For loss and gain, use merged maximal intervals as the countable unit and report per-genome
   binary alongside it. Not because merging changes the answer much (it does not, 0.966), but
   because the merged count is defensible as "distinct events" while a window count is defensible
   only as "windows".

The headline number should be the per-genome binary rate: what fraction of maternal genomes carry
any change, against what fraction of paternal genomes. It is the only unit that is immune to both
slicing and detector sensitivity, it is what the lab's question actually asks, and it is
interpretable without a denominator. Its cost is power, quantified in Q2d.

## Q2b. The denominator. Informative markers is the wrong one and it makes things worse.

Measured on the 26 arrays in audit/results.json that carry both quality fields:

- informative markers range 348,211 to 632,618, a spread of 1.82x
- the explainable-noise ceiling (the array's own artefact propensity) ranges 0.0065 to 0.1288, a
  spread of 19.7x
- Spearman(informative, ceiling) = -0.533, p = 0.005. They are NEGATIVELY correlated: arrays with
  more informative markers are cleaner.

That negative correlation is what breaks the normalisation. Comparing the 8 cleanest against the 8
noisiest arrays, and assuming false regions scale with artefact propensity:

| | bias between clean and noisy groups |
|---|---|
| raw counts | 5.89x |
| divided by informative markers (shipped) | 7.06x |

Dividing by informative markers amplifies the confounder by 1.20x rather than removing it, because
the noisy arrays have both more artefact AND fewer markers, so the division pushes the same
direction as the bias.

What would actually flatten it, tested the same way:

| denominator | residual bias |
|---|---|
| none (raw counts) | 5.89x |
| informative markers | 8.65x |
| called markers | 13.09x |
| explainable-noise ceiling | 1.00x |

The right denominator is the array's own artefact propensity, which you already compute per sample
as `explainable`. That is the quantity that predicts how many false regions an array will produce.
Marker count does not.

Recommendation: do not normalise by markers. Either
(a) use the per-genome binary unit, which needs no denominator, or
(b) keep counts raw and put the explainable-noise ceiling in the model as a per-sample covariate or
    offset, so a noisy array is expected to carry more findings and is not credited for them.

Your parenthetical about call-rate collapse being both a detection target and a hole in coverage is
exactly right and it is the reason (b) is preferable to any fixed denominator: the same physical
event both creates a finding and destroys the coverage you would measure it against. A marker-based
denominator cannot represent that; a per-sample noise term can.

## Q2c. The permutation is right. The statistic is not.

Label permutation is the correct frame here and I would not replace it with a parametric model.
Your reasoning in the code comment is sound: counts are neither normal nor Poisson, and permuting
labels keeps each sample's event and marker counts together.

Measured overdispersion on the 214-cell corpus: mean 7.07 regions per cell, variance 69.0, so
variance/mean = 9.8x Poisson, negative-binomial k = 0.81. A Poisson or quasi-Poisson model would be
badly wrong; a negative binomial with sample as a random effect is defensible but at 12 samples the
variance components are not estimable, and you would be trading a distribution-free test for a
model whose assumptions you cannot check. Stay with the permutation.

The defect is the test STATISTIC. `rateOf` pools: `sum(events)/sum(informative)` over the group.
Under label permutation the samples are exchangeable but their DENOMINATORS travel with them, so a
pooled ratio is not a symmetric function of the labels and the permutation null does not have the
right size. Measured, 7 vs 5, true null (both groups drawn from the same distribution), alpha 0.05:

| marker skew between groups | pooled (shipped) | mean of per-sample rates | rank of per-sample rate |
|---|---|---|---|
| 1.0 | 0.044 | 0.046 | 0.049 |
| 1.1 | 0.070 | 0.062 | 0.049 |
| 1.2 | 0.099 | 0.083 | 0.058 |
| 1.3 | 0.128 | 0.094 | 0.059 |
| 1.5 (the shipped guard limit) | 0.164 | 0.119 | 0.077 |
| 2.0 | 0.306 | 0.172 | 0.126 |

At skew 1.0 all three are correct. The shipped statistic degrades fastest: at the guard's own limit
of 1.5x it runs at 0.164, more than three times nominal. The rank statistic holds near nominal out
to about 1.4x.

Recommendation: replace the pooled-rate difference with the rank-sum of per-sample rates
(Wilcoxon statistic, permutation null). It is the smallest change that fixes the size, it needs no
distributional assumption, it keeps your seeded reproducibility, and it is robust to the one
disturbed genome that dominates a pooled count.

Also correct `minAchievableP`. The shipped formula returns 0.1429 at 3 vs 3 where the exact
two-sided minimum is 2/C(6,3) = 0.100, and 0.0423 at 4 vs 4 where the exact is 0.0286. It is
conservative rather than wrong, so it over-declares `underpowered`, but the exact value is one line:
2/C(nA+nB, min(nA,nB)).

## Q2d. The thresholds. Both are wrong, and here are measured replacements.

MAX_POWER_SKEW = 1.5. Derived from the type-I table above: at 1.5x the shipped statistic runs at
0.164 against a nominal 0.05. If you keep the pooled statistic, the largest skew holding
false positives at or below 1.5x nominal (0.075) is 1.1. If you adopt the rank statistic it is 1.4.

  keep pooled rate  -> MAX_POWER_SKEW = 1.1
  adopt rank        -> MAX_POWER_SKEW = 1.4

I would take the second: the guard is doing work the statistic should be doing, and a statistic that
is robust to the confounder is better than a threshold that refuses whenever the confounder appears.

Note the skew guard did not bind in your run either way: median informative was 551,558 maternal vs
518,853 paternal, skew 1.06. But that is an accident of which samples got a class. The three arrays
that fell out as "unclear" were the three worst (call rates 0.538, 0.556, 0.591, all paternal), so
the exclusion silently repaired the skew from 1.06 to 1.01 by dropping the paternal group's worst
members. The guard measures the surviving samples, not the ones the pipeline dropped, so a
quality-correlated exclusion is invisible to it. Add a second check: report the skew BEFORE the
class gate, and refuse if the EXCLUSION rate differs between groups (here 3/8 paternal vs 0/7
maternal).

MIN_PER_GROUP = 3. At 3 vs 3 the exact minimum achievable p is 0.100, which cannot clear alpha 0.05
at all, so no result of any kind can be reported. The current value permits a test that is
incapable of producing a finding. Measured floors and power (negative binomial k = 0.81 from your
corpus, alpha 0.05):

| n per group | min achievable p | clears 0.05 | clears 0.05/3 | power at 3x | power at 5x |
|---|---|---|---|---|---|
| 3 | 0.1000 | no | no | 0.00 | 0.00 |
| 4 | 0.0286 | yes | no | 0.14 | 0.22 |
| 5 | 0.0079 | yes | yes | 0.15 | 0.36 |
| 7 | 0.0006 | yes | yes | 0.28 | 0.55 |
| 8 | 0.0002 | yes | yes | 0.38 | 0.65 |
| 10 | <0.0001 | yes | yes | 0.45 | 0.77 |
| 12 | <0.0001 | yes | yes | 0.58 | 0.82 |
| 20 | <0.0001 | yes | yes | 0.82 | 0.99 |

  MIN_PER_GROUP = 5   (hard floor: below this, per-class testing cannot clear a Bonferroni
                       threshold for 3 classes even in principle)
  MIN_PER_GROUP = 20  (reporting floor for a headline claim: the smallest group size with 80%
                       power at a 3x effect)

Report both. Between 5 and 20 the test runs but should be labelled exploratory, and the report
should print the effect size it COULD have detected at 80% power rather than only the p.

For your actual run, 7 maternal vs 5 paternal: power is 0.17 at 3x, 0.30 at 5x, 0.48 at 10x. A
planted 7.5x effect is missed about two times in three. Your own note says 4 vs 4 admits no p below
0.042 and 8 vs 7 resolves it; that is right about the FLOOR but it conflates floor with power. 8 vs
7 can report a small p, but it only finds a 3x effect about 38% of the time.

On the per-genome binary unit, power is lower still: at a 1.5x relative risk you need roughly 40
genomes per group. That is the honest cost of the unit that is immune to slicing, and it is worth
stating in the report so the choice is visible.

## Q2e. Gyno vs andro answers a different question than the lab asked.

Your instinct is right and I would go further than the caveat you propose.

The lab's question is whether breakage falls differentially on the two parental genomes in embryos.
Gynogenetic and androgenetic conceptuses are not a sample of embryos; they are a sample of
fertilisation failures, and each class arises by a different mechanism. Whatever difference you
measure between them is confounded with the mechanism that produced the class, not just with the
parental genome. Two specific confounds, both of which run the same direction:

- a gynogenetic genome has never been through a sperm-derived protamine-to-histone exchange or the
  paternal pronucleus's first S phase, and an androgenetic one has never been through the oocyte's
  meiosis II completion. The two classes differ in replication history before any parental genome
  effect is considered.
- the classes have different quality distributions in your own data: median call rate 0.858
  maternal vs 0.807 paternal across all 15 arrays, and all three excluded arrays were paternal.

Additionally, and this is specific to the run in question: audit/results.json describes these
samples as "maternal nucleus isolated from 2PN zygote" and "paternal pronucleus". If they are
dissected pronuclei rather than conceptuses, they are haploid single genomes, not embryos, and
statements about parental balance in embryos do not transfer at all. The literature section below
establishes what the material actually is.

On using the biparental samples: yes, there is a defensible route, and it is not per-event calling.
Two options, in order of preference:

1. Restrict to WHOLE-CHROMOSOME events on biparental samples, where the dosage channel does have a
   detection floor on your material (trophectoderm 0.186 loss, 0.135 cnn-LOH) and can reach band A
   or B. You lose all segmental events, which is most of them, but the events you keep carry a
   measured per-event origin rather than an inherited one. Then the comparison is within-embryo,
   maternal against paternal events in the SAME genome, which removes every between-sample quality
   confound at a stroke. A paired test on within-embryo counts is far more powerful per sample than
   a between-group test, and it is the design the question deserves.

2. Report the biparental events with origin unresolved, and use them only for the location question
   (where breakpoints fall against fragile sites, long genes, replication timing), which your
   features track already does without needing a parent.

Option 1 is the one I would build. It answers the lab's actual question, it uses the material they
actually have, and it does not require the uniparental samples at all.

---

## Correction that changes the framing of the whole run

Before Q3, one factual correction, because it affects Q1 and Q2 as well.

GSE148488 is your own lab's data. It is the Zuccaro 2020 Cell paper, "Allele-Specific Chromosome
Removal after Cas9 Cleavage in Human Embryos", DOI 10.1016/j.cell.2020.10.025, PMID 33125898,
platform GPL28377 Affymetrix Axiom UK Biobank. Verified from the GEO record via eutils.

GSM4774673-4774679 are titled "maternal nucleus isolated from 2PN zygote N" and GSM4774680-4774687
"paternal nucleus isolated from 2PN zygote N". The source_name field reads "single nucleus isolated
from 2PN zygote at 20h post Cas9 RNP injection and fertilization" and characteristics reads
"tissue: 2PN zygote".

So these are dissected single pronuclei from normally fertilised BIPARENTAL zygotes, not
gynogenetic and androgenetic embryos. Three consequences:

1. They are a better validation set than embryos for the absence-detection step, because
   uniparental status is established by dissection rather than inferred. Nothing wrong with using
   them that way.
2. They cannot validate anything post-zygotic. A dissected pronucleus has not been through syngamy
   or a single mitosis, so every mechanism that would break the uniparental inference in real
   material - mosaicism, post-zygotic genome loss, lineage segregation, chimaerism - is
   structurally absent. Sensitivity measured here is an upper bound on embryo performance, not an
   estimate of it.
3. Every zygote in the series was Cas9 RNP treated, and the paper's finding is that unrepaired
   breaks produce segmental and whole-chromosome loss. The losses in this dataset are largely
   nuclease-induced. Using them to estimate the parental balance of spontaneous segmental loss is
   circular. Using them to test whether the caller detects known induced losses is legitimate.

This also resolves Q2e more sharply than my own answer above: the run that motivated
`parentalBalance` cannot speak to parental balance in embryos at all, whatever unit or test you
choose. It is a positive control, and it should be labelled as one in the report.

Minor: your accession range for the donor arrays is slightly off. The sperm-donor replicates start
at GSM4472397, not GSM4472398. Also, GEO does not state a marker count for GPL28377; the 825,656
figure is your own measurement off the files, and Rana 2023 quotes 820,967 for nominally the same
array, so the two should not be used interchangeably.

The one thing that could not be established: the WGA kit per sample. GEO says only "amplification
kit as indicated for each sample" and the indication is absent from the SOFT metadata. Since MDA
and PCR-based WGA differ in exactly the regional-bias property that drives spurious LOH, this is
load-bearing for the false-positive question in Q2a, and you can settle it from your own STAR
Methods.

---

## Q3a. Who runs this, and when

Your assumption is right, and the literature supports keeping it. This is a research tool for a
group re-analysing arrays after the fact, and the person at the keyboard is a bioinformatician or a
research-track embryologist, not a genetic counsellor and not a clinical lab director.

What makes that more than an assumption: ESHRE 2020 organisation good practice requires the
examination process and reporting to comply with local guidelines or ISO 15189, with documented
validation including "a summary of validation results" before any clinical PGT cycle. A browser
tool with no LIS integration, no sample-accessioning trail, and no validation dossier cannot enter
that path, and should not try to. The ceiling on what the output can be used for is set by that,
not by the quality of the statistics.

What it changes about the output: a research user wants the intermediate quantities and the
refusals, which is what you already produce, and does not want a clinical-report format. The
12-page PDF is arguably the wrong artefact for this user. A machine-readable table plus the viewer
would serve better, because the actual next step for a bioinformatician is to join your calls
against their own metadata. Keep the PDF for the record; make the TSV the primary output.

## Q3b. What decision does a parental-origin call feed

Blunt answer, since you asked for one: in a clinical path, nothing. There is no decision in PGT-A
or PGT-M that currently takes parent-of-origin as an input. No professional-body document addresses
reporting it - not the four ESHRE 2020 good-practice documents, not ESHRE 2022 mosaicism, not the
PGDIS statements, not the ACMG standards. That is a genuine void rather than an oversight I
happened not to find.

In a research path it feeds three real questions, and these are worth building for:

- mechanism attribution: whether an event is meiotic or mitotic, which is what haplarithmisis and
  karyomapping deliver and what determines whether an error is recurrent for a couple
- your lab's own question, on which see below
- counselling-adjacent research on recurrence risk, where a maternal meiotic error carries a
  different recurrence expectation than a post-zygotic one

But there is a sharper framing available for your driving question, and I would change the pitch.
The literature already answers the question as posed, in both directions:

| | answer | source |
|---|---|---|
| whole-chromosome aneuploidy | overwhelmingly maternal, strongly maternal-age dependent | McCoy 2015: maternal BPH trisomy 4.57% of blastomeres vs paternal BPH 0.16% |
| segmental imbalance | predominantly PATERNAL, carried by deletions | Tsuiko 2021: 61.5% paternal, n=162, P=0.004 |
| segmental mechanism | ~70% post-zygotic mitotic | Girardi 2020 |

So "do breakpoints and losses fall differentially on the parental genomes" is answered, and the
answer differs by event class in opposite directions. Re-deriving it validates the tool but is not
new biology.

What is genuinely open: the parental split of MITOTIC segmental events specifically. The reason it
is open is structural rather than incidental. Haplotype methods DEFINE a mitotic event as one
showing a normal biparental haplotype pattern, so the method that yields the meiotic parental split
cannot in principle yield the mitotic one. Tsuiko 2021 reports that post-zygotic aneuploidy does
not discriminate between parental homologs; Handyside 2025 counted mitotic segmental events without
assigning them a parent.

A dosage-based method that assigns a parent to a mitotic segmental event without relying on the
haplotype definition would be new. That is the defensible contribution, it is narrower than the
current pitch, and it happens to be exactly what a self-referenced allele-dosage channel is for.
I would rewrite the tool's stated purpose around it.

## Q3c. Bands C and D. Do not show the number.

Show the verdict, dim it, and drop the digits. The field has already run this experiment.

The PGT-A mosaicism episode is the precedent and it is directly on point. ESHRE's 2022 survey of
239 centres found 80% of centres that biopsy three or more cells report mosaicism, while only 66.9%
of all centres had validated their technology and only 61.8% of those had validated specifically
for calling mosaicism. Criteria for designating and reporting it "vary significantly across the
centres". Campos 2023 argues from the concordance data that intermediate copy-number profiles do
not represent strong evidence of mosaicism but only an inaccurate and misleading assumption, and
that laboratories should limit predictions to euploid and aneuploid and stop reporting mosaicism
entirely.

The lesson is not that hedging language failed. It is that a low-specificity intermediate category
gets reported and quoted regardless of how it is labelled, and the field's corrective was to
propose withdrawing the category, not to improve the caveat. A parental call at 0.60 accuracy is a
low-specificity intermediate category, and "not for reporting" is exactly the kind of label that
episode showed does not survive contact with a downstream reader.

Concretely, for band D at accuracy 0.60: a coin flip is 0.50. You are offering four percentage
points of information over guessing, with three decimal places of apparent precision attached.
Render it as "not evaluable" with the measurement available on hover or in an appendix table, and
keep the number out of the row. Band C at 0.75 I would show as a verdict with no number. Bands A
and B carry numbers.

This is the same argument as Q1b arriving from a different direction, and the two should be settled
together: the tool currently has one confidence scale spanning a Mendelian certainty and a 0.60
guess, and neither end is served by it.

## Q3d. The run log

What an operator needs while a run proceeds is: which file is being processed, how many remain, and
any refusal that will affect the final answer. That is three lines, updated in place.

Everything else belongs in the report. The specific failure you found - 319 identical warnings - is
worth a general rule rather than a fix: a warning that fires per marker or per window is a property
of the ARRAY, not of the run, so it should be counted and emitted once with its count at the point
where it becomes a refusal. "319 windows below the marker floor" is information; 319 copies of the
same string is noise that hides the one line that mattered.

The one thing that should interrupt an operator mid-run is a refusal that makes the rest of the run
pointless, for example call rate below the floor on the parent array, since every downstream
channel depends on it. That should stop and say so rather than proceed to produce a 12-page PDF of
refusals.

## Q3e. How it sits beside PGT-A and PGT-M

Not redundant, but the niche is narrower than it looks, and for a structural reason.

Every published parent-of-origin method requires BOTH parents genotyped, and most also require a
phasing reference:

| method | input | both parents | extra reference | parent for an aneuploidy | meiotic vs mitotic |
|---|---|---|---|---|---|
| karyomapping (Handyside 2009) | SNP array, single cell | yes | yes, sibling or relative | yes, parental and grandparental | partially |
| haplarithmisis / siCHILD (Zamani Esteki 2015) | SNP array BAF + genotypes | yes, phased | yes, grandparents or sibling | yes | yes, except monosomies |
| Parental Support (Johnson 2010, McCoy 2015) | SNP array, blastomere | yes | no | yes | yes, via BPH/SPH |
| MeioMap (Ottolini 2015) | all products of one meiosis | mother reconstructed | the other products are the reference | maternal by construction | yes |
| low-pass WGS PGT-A | 0.01-0.05x WGS | no | no | NO in standard practice | no |

Two corrections to your framing. Karyomapping is a three-genome minimum, not two: parents plus a
sibling or other relative to establish phase. And haplarithmisis excludes monosomies from its
meiotic/mitotic tracing, which means the field's most complete method declines to call mechanism on
exactly the event class your tool primarily targets.

So one parent plus the sample is not a weaker version of a standard, because there is no standard
operating in that configuration. But the reason is structural: an informative site is conventionally
defined as one parent homozygous reference and the other homozygous alternate, so with one parent a
site where the sample is homozygous for the known parent's allele is ambiguous between "inherited
that allele from the absent parent too" and "the absent parent's allele is missing". Recovering the
absent contribution needs either a population LD panel, which is the LD-PGTA route (Ariad 2021), or
an assumption about the absent genome. Your channels do statistical inference the two-parent
methods do not have to do, and that is where the detection floors in your own table come from.

The closest published number to your problem: haplarithmisis detected segmental events above
4.4 Mb using about 88,000 informative SNPs genome-wide, with breakpoint accuracy averaging 0.51 Mb
(SD 0.48). That is a two-parents-plus-reference number, and it is the bar a one-parent method is
measured against.

The genuinely non-redundant part of your tool is the uniparental case, where absence itself carries
the information and no second parent is needed. That is the defensible core, subject to Q1a's
mechanisms and the implementation defect.

## Q3f. What would make a lab distrust it on sight

In the order they would hit it:

1. A parental call on every one of 1,272 rows from a single genome-level fact. Anyone who has
   worked with these data will ask what the independent unit is within a minute, and if the answer
   is "one", the row count reads as inflation. This is the biggest one and it is Q1c.
2. Three decimal places on an uncalibrated number. 0.659 asserts a precision that the underlying
   quantity does not have; a reader who checks will find the whole observed range spans 0.03 and
   will discount everything else in the report.
3. A confidence band shared between a Mendelian certainty and a 0.60 coin flip. When a reader
   notices that the gynogenetic inference and the weakest dosage guess are both band D, they will
   conclude the bands are decorative.
4. Copy-neutral LOH called across a homozygous genome. A reviewer who knows the material will ask
   how you call loss of heterozygosity in a genome that has none, and the answer needs to be that
   you do not.
5. "Research use only" on an output formatted as a clinical report. The 12-page PDF signals a
   deliverable; the disclaimer says it is not one. Pick one.
6. A validation set that is the authors' own Cas9 experiment, described as embryos when the GEO
   record says dissected pronuclei. This one is recoverable simply by describing it accurately, and
   it costs nothing to fix, but if a reader finds it themselves it costs the rest of the report.

None of these are about the statistics being wrong. They are about the output claiming more
resolution than the measurement has, which is the same failure the mosaicism episode turned on.

---

## What I would do next, in order

1. Fix the one-sided gate in `uniparentalOrigin` or make it refuse explicitly. Nothing else in the
   parental-balance chain means anything until this is done, because it guarantees a maternal
   excess. Measured: 7/7 vs 0/5.
2. Gate `detectLoh` on zygosity, or on an absolute heterozygosity floor. This is probably where
   most of the 1,272 comes from, and it is a three-line change.
3. Drop the confidence number from the zygosity channel; emit verdict plus fold.
4. Collapse per-event parental rows to one genome-level statement per uniparental sample.
5. Replace the pooled-rate statistic with the rank of per-sample rates; set MAX_POWER_SKEW to 1.4,
   MIN_PER_GROUP to 5 as a hard floor and 20 as a reporting floor; use the exact
   2/C(n, min(nA,nB)) for the achievable floor.
6. Relabel the GSE148488 run as a positive control, not a parental-balance result.
7. Re-pitch the tool at the parental split of mitotic segmental events, which is the open question
   the existing methods structurally cannot answer.

## The measurement I would make before any of this

There is one empirical number nobody has published and you are unusually well placed to produce:
the tract-scale false-positive rate for WGA-induced apparent LOH or apparent parental absence on a
SNP array, as a function of tract length. Per-site allele dropout is measured (10.78% of loci on
MDA trophectoderm biopsies, Zhu 2026), but the independent-versus-correlated dropout distinction is
what actually sets your error budget, and your own segments.ts comment already documents that runs
of two occur 6 to 10 times more often than independence predicts.

GSE148488 contains the right controls: the sperm-donor and egg-donor replicate arrays are diploid
and biparental by construction, so any apparent uniparental tract in them is an artefact. That
yields the false-tract-rate-versus-length curve on your exact platform, which is the number that
would let you state a real detection floor for segments instead of "none". It would also settle
whether the 183 findings are what I think they are.


---

## References

GSE148488 series record - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148488
Zuccaro 2020, allele-specific chromosome removal after Cas9 cleavage in human embryos - https://doi.org/10.1016/j.cell.2020.10.025
Handyside 2009, karyomapping - https://doi.org/10.1136/jmg.2009.069971
Zamani Esteki 2015, haplarithmisis / siCHILD - https://doi.org/10.1016/j.ajhg.2015.04.011
McCoy 2015, common variants spanning PLK4 and mitotic aneuploidy - https://doi.org/10.1126/science.aaa3337
Tsuiko 2021, karyotype of the blastocoel fluid and parental origin of segmental imbalance - https://doi.org/10.1038/s41525-021-00246-0
Girardi 2020, segmental aneuploidies are predominantly mitotic - https://doi.org/10.1016/j.ajhg.2020.03.005
Ottolini 2015, MeioMap, genome-wide maps of recombination in human oocytes - https://doi.org/10.1038/ng.3306
Gruhn 2019, chromosome errors across the female reproductive lifespan - https://doi.org/10.1126/science.aav7321
Destouni 2016, zygotes segregate entire parental genomes into distinct blastomere lineages - https://doi.org/10.1101/gr.200527.115
Sunde 2011, biparental-marker diploid moles and androgenetic/biparental mosaicism - https://doi.org/10.1038/ejhg.2011.93
Ariad 2021, LD-PGTA, parent-free recovery of BPH/SPH from low-pass sequencing - https://doi.org/10.1126/sciadv.aaz7602
ESHRE PGT Consortium 2020, good practice recommendations for the organisation of PGT - https://doi.org/10.1093/hropen/hoaa021
ESHRE PGT Consortium 2020, good practice recommendations for PGT-SR and PGT-A - https://doi.org/10.1093/hropen/hoaa017
ESHRE PGT Consortium 2020, good practice recommendations for PGT-M - https://doi.org/10.1093/hropen/hoaa020
Campos 2023, argument against reporting PGT-A mosaicism - https://doi.org/10.1007/s10815-023-02936-3
Gonzales 2021, ACMG revision on regions of homozygosity and suspected UPD - https://doi.org/10.1038/s41436-021-01203-z
Riggs 2019, ACMG/ClinGen CNV interpretation standards; erratum Genet Med 2021;23(11):2230 - https://doi.org/10.1038/s41436-019-0686-8

Full literature record, including the 55-row retraction and correction audit and 10 stated gaps:
LIT-parental-balance-consult.md

Retraction and correction status: 48 DOIs checked against OpenAlex is_retracted and PubMed
PublicationType plus CommentsCorrections. Zero retractions, zero expressions of concern. Three
errata exist and were not read: Masset 2022, Deveault 2008, and Riggs 2019 (correction at Genet Med
2021;23(11):2230, which should be pulled before implementing any CNV scoring rule).

---

## Measurement provenance

Every number in this document is computed from your repo or from a named primary source. The
computations are in fig_balance_diagnostics.png and the table below.

source | what was measured
--- | ---
audit/results.json | 26 arrays with absence, ceiling, margin, informative, callRate, hetBand. Source of the 7/7 vs 0/5 gate result, the skew and exclusion analysis, and the denominator correlation.
audit/features/regions-*.tsv | 214 cells, 1,513 regions across 5 cohorts. Source of the unit table, the dispersion estimate k=0.81, and the merge comparison.
web/src/parentalBalance.ts | reimplemented faithfully in Python, including the LCG, to measure type-I error and power. Confirmed to have no callers anywhere in the repo.
web/src/uniparentalOrigin.ts, oneParentOrigin.ts | reportedConfidence reimplemented; band boundaries and the 0.206 bound verified (rule of three on 13 units = 0.2058).
web/src/abnormalities.ts | LOH_DEPLETION = 0.65 and the self-referenced background; simulated to produce the false-window table.

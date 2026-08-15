# Handover: answer the Egli question with OriginMarker, and map the events onto genomic factors

Everything below runs on this machine. The tool, the arrays and the harnesses are already here.
Nothing needs downloading except the factor tracks in part 3.

## The question

Dr. Dieter Egli, Columbia: **"The insight we are looking for is if breakpoints/losses affect
paternal or maternal genomes differentially, or equally. Either answer is very valuable and can be
published. It reflects a different biology."**

So the unit is the EVENT, deduplicated within embryo, and the deliverable is a per-event parental
label with a defensible error rate plus an ascertainment argument. A genome-wide summary statistic
is not the deliverable.

## Where everything is

**Run the corpus work on the NAS, not on this Mac.** The 884 arrays used to live on a USB stick,
which is why they were unreachable: absent whenever the stick is unplugged and invisible to any
process that cannot see that volume. They now live on the NAS, which is always up.

    wrapper         scripts/om-nas.sh <om subcommand and flags...>
    repo (Mac)      /Users/ezrakruger/claudecodegeneralworkspace/originmarker
    repo (NAS)      /volume1/docker/originmarker         updated by every deploy, same code
    lab arrays      /volume1/docker/originmarker-data/"SNP array data"
                    884 .probes, 4 groups: DIETER 93, JENNA 178, ROBLES 264, TREFF 349
    public series   GSE148488 (135 arrays, 1 sperm donor + 4 egg donors + embryos)
                    GSE290961 (158 arrays, 9 bulk-quality candidates, 6 with confirmed children)
    called regions  audit/features/regions-*.tsv   1,189 regions, 318 gains, 871 losses,
                    on 175 of 884 arrays. Columns: cell, chr, start_bp, end_bp, log2R_dev, kind.
                    NO parental label on any of them. That is the gap.

If the GEO series are not still in the session scratchpad, re-fetch: the accession list is in
`audit/asymmetry/`, and the fetch pattern is in the git history of `scripts/`.

## 1. How to run the tool

Node 22+. No install step for the CLI; it imports the same modules the web app runs, so it cannot
drift from what ships.

    # On the NAS, against the corpus. This is the one to use for anything touching the arrays.
    ./scripts/om-nas.sh help
    ./scripts/om-nas.sh census /data
    ./scripts/om-nas.sh stage "/data/SNP array data/DIETER/<file>.probes"

    # Locally, for the public series or for development.
    cd /Users/ezrakruger/claudecodegeneralworkspace/originmarker
    node --experimental-strip-types cli/om.ts help
    node --experimental-strip-types cli/om.ts constants     # every tunable + its provenance

Inside the container `/w` is the repo and `/data` is the corpus, so corpus paths start `/data/`.
The NAS host runs node v18 and this needs 22, so the wrapper uses a `node:22-slim` container: the
right runtime, no change to a host that also serves the live site.

**Memory is the one real constraint.** The NAS has roughly 1 GB free of 7 GB. Use `--stride 8` for
anything sweeping the whole corpus, which is what the census and screening paths already do and
costs nothing statistically for a proportion measured over 100,000 markers. `OM_NAS_HEAP=6144
scripts/om-nas.sh ...` raises the heap if a single full-density array needs it, but do not run two
at once. If it is not enough, the alternative host is `tesla` (12 cores, 31 GB, node 22 already in
/opt/node22), which was down when this was written and needs a physical poke.

    om stage       <array>                     material, dropout, marker floor
    om link        <parent> <sample>...        child / duplicate / unrelated / refused
    om origin      <parent> <sample>           which parent's copy is missing, per region
    om cohort      <dir> --ref <array>         the same for every confirmed child in a directory
    om census      <dir>                       haploid products per group
    om reconstruct <product>...                a parent's genotypes from her own haploid cells
    om enrich      <regions.tsv> --track T --markers A

Flags that matter here: `--json`, `--stride N` (screening only), `--channel
auto|genotype|dosage|both`, `--state loss|gain|cnn-loh`, `--parents 1|2`, `--region chr6:1-9000000`,
`--out FILE`. Every constant in `om constants` is a flag, lower-cased with dashes. **Omitting all
flags is exactly the configuration the web tool runs and the audits measured**; any override is
printed in the output because a number produced under a changed constant is not comparable with the
validation figures.

Before and after any change: `./scripts/checks.sh` (33 self-checks), `./scripts/release-check.sh`.

Harnesses already written, in `audit/asymmetry/`:

    one_parent_cohort.ts   quality gate -> bulk reference -> linkage established -> regions from
                           intensity -> origin. Three bounded passes; at most two arrays resident.
    haploid_census.ts      haploid products per donor group
    find_examples.ts       bulk-quality candidates and their confirmed children
    residual_corr.ts       the dosage/intensity channel correlation

## 2. What is established, so you do not re-derive it

Read `audit/ORIGIN-CONSULT.md` with `boundary_of_certainty.csv` and `refusal_taxonomy.csv` first.
It is a prior methods review on 135 arrays of this platform and it is the authority for the floors.

Load-bearing facts:

- **Origin is defined by reference to a parent.** With none, the false origin-call rate is 1.0000
  against 0.0790 with a real parental array. There is no parent-free path.
- **Copy-neutral LOH is the EASIEST class**, not the hardest: deviation f/2, the largest of the
  three, because copy number stays at 2 so the caller never degrades. On trophectoderm with ONE
  parent its origin floor is 0.186, callable, where a loss on the same array is 0.625.
- **A single blastomere loss with one parent has no floor at any fraction to 0.70**, and still none
  at eight cells. Material limit; no file the user can load moves it.
- **A second parent is worth 3.36x**, more than four extra cells of the same embryo (1.38x).
- **A parent reconstructed from >=5 of her own HAPLOID products is indistinguishable from her bulk
  array.** Reconstruction from diploid SIBLINGS is disqualifying and asymmetric: it flips the sign
  on 29-38% of chromosomes and does not converge.
- **Ascertainment is the whole risk.** If maternal and paternal events are detected with unequal
  sensitivity, true parity reads as asymmetry: on one-cell-per-embryo material a 7.7% maternal
  shortfall alone manufactures an apparent paternal fraction of 0.52. Claimable floor: nothing
  below a paternal fraction of 0.60, nothing on one-cell-per-embryo material, roughly 20-30 oocyte
  donors before a null means anything.
- Two channels exist. Obligate-het (parent homozygous) and untransmitted-haplotype (parent
  heterozygous, sample homozygous, disjoint set, 1.40-1.98x SNR, the only channel giving a
  blastomere any floor at all).

## 3. What to do

### Step 1, the one that turns "unknown" into a number

Sort the haploid products by parent. `om census` says six donor groups hold five or more haploid
products, but a group holds products of BOTH parents and a reconstruction needs five from ONE. In
the single group already sorted, six became four. So the real count of anchorable groups is unknown
and everything downstream depends on it.

    om census "/Volumes/SANDISK USB/SNP array data"
    om reconstruct <the haploid products of one group>... --group-only --pairs
    om reconstruct <the products of the largest single-parent group>... --out parent.probes

Products of one parent share half their haplotypes and are unrelated to the other parent's, so they
separate on mutual relatedness with no reference. Report: per group, how many products, how many
parents they resolve into, how many products per parent, and which groups clear five.

### Step 2, the cohort run

For every group with an anchor, run the cohort and attribute what is attributable.

    om cohort "<group dir>" --ref parent.probes --role maternal --json --out group.json

Then filter the 1,189 regions to those on arrays that (a) have an anchor, (b) sit on a callable
material/state/parent combination, (c) pass the array gate. Report the count that survives each
filter separately: the attrition is itself the result, and it is what tells Dr. Egli what a bigger
study would need.

**Look specifically for copy-neutral events on trophectoderm.** That is the one class callable with
a single parent on this material, at 0.186, and nobody has searched this corpus for it. Detection
there must NOT use heterozygosity depletion (floor 0.65-0.90, and 1% of clean windows already read
40% depleted from clustered dropout); use the allelic deviation as the trigger and report depletion
as a consequence.

### Step 3, the factor mapping, which is the part I most want widened

`audit/features/README.txt` has the existing analysis: 1,189 regions against four tracks with a
marker-matched permutation null. Result was negative everywhere, every ratio between 0.77 and 1.04,
including no fragile-site enrichment at n=1,189.

    om enrich audit/features/regions-ROBLES.tsv --track <features.json> --markers <an array>

**The null is the entire difficulty and it is already solved: do not change it.** A region can only
be called where the array carries markers AND where amplification produced calls, so a uniform null
reports enrichment for anything correlated with marker density, which gene density is and fragile
sites therefore are. Every null interval is drawn on the SAME chromosome with the SAME number of
informative markers.

Current tracks: common fragile sites (24, aphidicolin-induced), genes over 500 kb, centromere and
telomere gaps, late-replication valleys (ES and constitutive). All hg19, all public, sources named
in the README.

**Add every factor you can justify.** The ask is explicitly "all of them". Candidates, and please
add any I have missed:

- CpG islands, shores and shelves
- Imprinted loci and imprinting control regions. This one is the most interesting for a
  parent-of-origin question: if breakage tracked imprinted domains, the parental asymmetry would
  have a mechanism rather than just a number.
- Methylation: germline differentially methylated regions, oocyte-specific and sperm-specific
  domains, partially methylated domains. Note the material is preimplantation, so use
  gamete/early-embryo maps rather than somatic ones where they exist.
- Replication timing beyond the valleys already used: full Repli-seq domains, early/late
  transitions, replication origins
- Recombination hotspots and the deCODE map (the repo already ships deCODE maps in `data/maps`)
- Topologically associating domain boundaries, CTCF sites, cohesin binding
- Segmental duplications, repeats, palindromes, satellite and pericentromeric regions
- G-quadruplex motifs and other secondary-structure predictors
- Chromatin state in early embryo or ESC: open chromatin, histone marks, lamina-associated domains
- Nuclear position and chromosome territory proxies
- Gene expression level in early embryo, and zygotic genome activation timing
- Structural variant hotspots from population catalogues

For each: state the source, the build, and whether it is a somatic map being applied to embryonic
material, which is a real transfer risk worth flagging in the output rather than in a footnote.

**Two disciplines from the existing work.** Read effect sizes, not p-values: with 1,189 regions the
permutation null is narrow enough that a four percent difference reaches p 0.0005 while meaning
nothing. And split by event class, since gains and losses may differ and a pooled null can be an
artefact of pooling.

**And the split that has never been done:** once step 2 produces parental labels, run the factor
comparison SEPARATELY for paternal and maternal events. That is the direct form of Dr. Egli's
question. Be aware of the caveat already on record: the late-replicating fragile compartment is
reported as established symmetrically on both parental genomes from the first cell cycle, so a
positional result may well not license a parental one. Say so if it does not.

## 4. What to report

- How many donor groups can be anchored, and by what route
- How many of the 1,189 regions survive each filter, with the attrition shown
- The parental split of what survives, with its denominator and its ascertainment caveat stated in
  the same sentence as the number
- The factor comparison, all tracks, effect sizes foremost, split by event class and by parent
- What a study designed to answer this properly would need, in samples and in genotyped parents

Refuse rather than reach. A wrong parental call is indistinguishable from a right one to a reader,
and this project has already produced six one-directional results that were all wrong. If the
honest answer is that the corpus cannot support the general claim, that is a publishable finding
about study design and should be written as one.

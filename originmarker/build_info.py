"""
Build identity: version + codename.

Releases are named after the biology of crossing over. RELEASES is ordered by how good the
name is, not by the meiotic sequence.

Bump VERSION and CODENAME together. CODENAME must be the next unused entry in RELEASES.
"""
from __future__ import annotations

from typing import NamedTuple


class Release(NamedTuple):
    name: str
    gloss: str


# Each gloss ships in the docs and in as_dict(), so a name without a definition is not
# enough. Add new names at the position their quality earns.
RELEASES: tuple[Release, ...] = (
    Release("Chiasma",
            "The X-shaped point where two chromatids have crossed over - the visible "
            "evidence of recombination, and the event every flanking marker is chosen to "
            "bracket."),
    Release("Synapsis",
            "The lengthwise pairing of homologous chromosomes that has to happen before "
            "they can exchange anything."),
    Release("Pachytene",
            "The substage of prophase I in which crossing over actually occurs."),
    Release("Holliday",
            "The four-armed branched junction intermediate that resolves into a crossover."),
    Release("Bivalent",
            "A pair of synapsed homologues, held together by the chiasmata between them."),
    Release("Diplotene",
            "The substage where the homologues pull apart and the chiasmata become visible."),
    Release("Synaptonemal",
            "The protein scaffold that zips paired homologues together along their length."),
    Release("Leptotene",
            "The beginning of prophase I, where chromosomes first condense into threads."),
    Release("Zygotene",
            "The substage where homologues find each other and begin to pair."),
    Release("Diakinesis",
            "The final condensation of prophase I, chiasmata still holding the homologues."),
    Release("Tetrad",
            "The four chromatids of a paired chromosome pair, the unit crossover acts on."),
    Release("Cohesin",
            "The ring complex holding sister chromatids together, without which crossovers "
            "could not be resolved correctly."),
    Release("Kinetochore",
            "The structure that couples a chromosome to the spindle and pulls it to a pole."),
    Release("Recombinase",
            "The enzyme class that catalyses strand invasion and exchange."),
    Release("Crossover",
            "The reciprocal exchange itself - the reason a marker on one side of a locus "
            "can stop predicting the other."),
    Release("Anaphase",
            "The moment the two copies are pulled to opposite poles - the physical event "
            "behind the question of which parent a copy came from."),
    Release("Disjunction",
            "The clean separation of homologues, whose failure is most of what this tool "
            "is looking for."),
    Release("Interference",
            "The rule that one crossover suppresses another nearby, which is why markers "
            "cannot be treated as independent."),
    Release("Centromere",
            "The constriction the spindle grips, and the region whose retained "
            "heterozygosity was long read as evidence of a meiosis I error."),
    Release("Chromatid",
            "One of the two copies a replicated chromosome is made of."),
    Release("Telomere",
            "The chromosome end, and the boundary a terminal event runs to."),
    Release("Bouquet",
            "The prophase I arrangement where telomeres gather at one pole of the nucleus - "
            "position in the nucleus, rather than sequence, deciding what meets what."),
    Release("Synizesis",
            "The tight clustering of chromatin that accompanies the bouquet."),
    Release("Dictyate",
            "The arrested prophase an oocyte holds for decades before it finishes dividing."),
    Release("Pronucleus",
            "Either parental genome after fertilisation and before they combine, which is the "
            "one moment the two are separable by looking."),
    Release("Karyogamy",
            "The fusion of the two pronuclei, after which parent of origin can only be inferred."),
    Release("Metaphase",
            "The alignment before separation, where a chromosome is most visible."),
    Release("Spindle",
            "The apparatus that pulls the homologues apart, and misattaches when it fails."),
    Release("Nondisjunction",
            "The failure of a pair to separate, which is the event most of what this tool "
            "detects actually is."),
    Release("Karyotype",
            "The complete set of chromosomes taken together, which is what the genome map "
            "draws rather than one chromosome at a time."),
    Release("Univalent",
            "A chromosome with no partner at metaphase I, which segregates at random."),
    Release("Syngamy",
            "The fusion the tool is named for, and the moment after which the two parental "
            "genomes are one and every question about them is asked indirectly."),
    Release("Cleavage",
            "The divisions from zygote to morula, where a blastomere still carries both parental "
            "genomes and an absent allele is therefore Mendelian rather than statistical."),
    Release("Interphase",
            "The state between divisions, where a genome is being copied rather than segregated "
            "and much of what looks like structure is flux."),
    Release("Amphimixis",
            "The union of the paternal and maternal genomes into one nucleus, after which "
            "there is a single genome rather than two that have to be kept in step."),
    Release("Anaphase II",
            "The second separation, where sister chromatids part and a premature one is the "
            "commonest source of an aneuploid egg."),
    Release("Blastomere",
            "One cell of a cleaving embryo, carrying both parental genomes in a single copy of "
            "the material, which is where parent of origin becomes answerable and where the "
            "assay is at its noisiest."),
    Release("Reciprocal",
            "The pair of products of one exchange, which carry the same event read from opposite "
            "sides and are only fully described together."),
    Release("Polar Body",
            "The discarded sister product of a meiotic division, biopsied in clinical work to "
            "infer indirectly what the oocyte kept. A second reading of the same event, from the "
            "side that was thrown away."),
    Release("Zona Pellucida",
            "The coat around the oocyte that decides what is allowed to fuse with it, and blocks "
            "everything after the first. A check on the identity of what has been presented, "
            "before anything downstream depends on it being what it claims."),
    Release("Parthenote",
            "An embryo developing from an egg alone, carrying a maternal genome and no paternal "
            "one. The case the genome-level class named backwards."),
    Release("Heteroduplex",
            "A stretch where two strands from different sources are paired, and mismatches show "
            "only once both are read. What a second study is for."),
    Release("Resolvase",
            "The enzyme that cuts a junction one way or the other. Leaving it uncut is also an "
            "outcome, and the safer one when neither cut can be justified."),
    Release("Gene Conversion",
            "Where one sequence is overwritten by the copy it was paired with, and the record of "
            "what was there is gone. What a reference built from the wrong majority does."),
    Release("Plasmogamy",
            "The step where two cells join but their nuclei have not yet met. A handoff that has "
            "to complete before anything downstream can happen at all."),
    Release("Capacitation",
            "The conditioning a sperm undergoes in the female tract before it can fertilise "
            "anything: the same cell, competent only once it has adjusted to where it is."),
    Release("Isodisomy",
            "Both copies of a chromosome from one homologue of one parent, visible only as "
            "heterozygosity gone. One copy reads the same way, and only intensity or the array's "
            "own background can say which."),
    Release("Monosomy",
            "One copy of a chromosome where there should be two. Whose copy survived is what the "
            "alleles say; whether one copy remains or two is what only the intensity can."),
    Release("Gamete",
            "The haploid cell a parent's genome is read back from, one meiotic product at a "
            "time."),
    Release("Homologue",
            "The partner chromosome a copy is measured against. Whose homologue remains is what "
            "the alleles say; how many copies of it is what the intensity says."),
    Release("Trisomy",
            "Three copies of a chromosome where there should be two. Which parent supplied the "
            "third is what the alleles say, and whether there is a third at all is what only the "
            "intensity can confirm."),
    Release("Mosaic",
            "A change carried by some cells of an embryo and not others, which moves every "
            "measurement part of the way and is the hardest thing on an array to be sure of."),
    Release("Equational",
            "The second meiotic division, which parts sister chromatids and leaves the "
            "chromosome count unchanged. Two divisions reaching the same count, told apart by "
            "what the alleles say and not by how much material is there."),
    Release("Reductional",
            "The first meiotic division, which parts the two homologues and halves the "
            "chromosome count. Each parent contributed one of them, so which homologue has gone "
            "is a question that has to be put to both and not only to the first."),
    Release("Hemizygous",
            "One copy of a region where there should be two, so only one parent's alleles are "
            "there to read. Which parent is a question each parent's own markers answer "
            "separately, and a region that reads absent against both is answering neither."),
    Release("Achiasmate",
            "A chromosome pair that formed no crossover, so the two homologues separate exactly as "
            "they arrived. There is nothing in the result to say a crossover was looked for and "
            "not found, which is the same silence a chromosome returns when it was barely "
            "genotyped at all."),
    Release("Asynapsis",
            "Homologues that fail to pair, so neither has the partner it would be read against. "
            "A chromosome judged by the median of its neighbours has the same difficulty once most "
            "of them carry events: the reference it needs is not there, and nothing in the result "
            "says so."),
    Release("Dyad",
            "One chromosome of a pair with its two chromatids still joined, which is what a cell "
            "carries between the two meiotic divisions. One copy where there should be two has no "
            "heterozygous sites left to call, and a region that failed to amplify has none "
            "either, which is the whole of the difficulty in telling them apart."),
    Release("Chimaera",
            "Two genomes in one body, which raises the heterozygous rate at every marker because "
            "two people differ everywhere. A genome whose own chromosomes are at the wrong number "
            "raises it on those and nowhere else, and that difference is the whole of how the two "
            "are told apart."),
    Release("Nullisomy",
            "Both parental copies of a region gone, rather than one. Two absence scans, one per "
            "parent, both read it as missing, which is the one case where neither of them is "
            "naming a parent and both are right."),
    Release("Heterodisomy",
            "Both homologues of ONE parent, rather than one from each. It is heterozygous wherever "
            "a crossover separated them, so it reads like a biparental genome on the marker class "
            "that counts heterozygosity, and like a single contribution on the class that asks "
            "whose allele is present. The case a one-parent call has to decline."),
    Release("Euploidy",
            "The right number of chromosomes. Naming it is a separate question from naming the "
            "material an array came from, and until now only the material had a word: a diploid "
            "cell was described as a trophectoderm and nothing said it carried two genomes."),
    Release("Obligate Crossover",
            "Every bivalent must have at least one, or the homologues do not segregate. The "
            "crossover is not optional and it is not a property of the cell trying harder: either "
            "the pairing made one available or it did not. The same split applies to a marker: a "
            "chromosome where the parents are rarely opposite homozygotes offers nothing to read, "
            "however good the array, and that is a different refusal from an array that was handed "
            "thousands of such markers and returned ninety."),
    Release("Free Recombination",
            "Two loci far enough apart that crossovers between them randomise completely. The "
            "recombination fraction reaches 0.5 and the pair carries no linkage information at "
            "all, which no amount of further data recovers: the markers are not broken, there is "
            "nothing there to measure. An array whose calls agree with a parent's at exactly the "
            "rate chance predicts is the same statement about genotypes rather than about "
            "distance, and a two-person mixture reads below that rate rather than at it."),
    Release("Unequal Exchange",
            "A crossover between misaligned repeats, which hands one chromatid a duplication and "
            "the other the reciprocal deletion. The two products are one event and they do not "
            "announce themselves alike: the deletion takes the DNA away and the array stops "
            "calling there, while the duplication leaves every marker present and callable. A "
            "scan that looks for where genotyping stopped therefore finds one of them and is "
            "blind to the other, however it is thresholded."),
    Release("Mismatch Repair",
            "Where a crossover leaves two strands disagreeing, the cell does not guess which base is "
            "right: it reads the strand it trusts and corrects the other against it. The same move "
            "applies to an array. Where both parents carry the same homozygote the right answer is "
            "known before the sample is read, so every disagreement there is this array's own error, "
            "measured on it, and the calls that matter are read with that error removed rather than "
            "against a line drawn once on other material."),
    Release("Linkage Equilibrium",
            "In unrelated people, alleles at one marker carry no information about who anyone's "
            "parents are, because no shared pedigree ties them together. That is what lets a panel "
            "of strangers stand where a family's parents are missing: where every one of them "
            "carries the same allele on both chromosomes, a sample almost surely does too, so a "
            "heterozygous call there is the array's error, measured without a single relative."),
    Release("Haplotype Block",
            "Recombination is rare enough that a chromosome passes to a child in long stretches, "
            "each copied whole from one of the parent's two homologues. Two of a parent's gametes "
            "carry the same stretch half the time, and all of a handful of them some of the time, "
            "and there a parent rebuilt from its gametes cannot see that it had two."),
    Release("Coupling Phase",
            "Alleles that sit on the same parental chromosome travel together, in coupling, and an "
            "array reads a cell's genome through them as well as through its heterozygosity. Where "
            "dropout empties the heterozygous calls, the other parent's alleles are still there, "
            "and they say whether two genomes are present."),
    Release("Repulsion Phase",
            "Alleles on opposite homologues are in repulsion: one parent's chromosome carries one, the "
            "other parent's the other. Count the copies of each and an extra chromosome shows whose "
            "it is, with only one of the two parents there to read them against."),
    Release("Identity by Descent",
            "A child holds one of each parent's two chromosomes, whole between the crossovers: "
            "identical by descent to that parent along its length. An array that carries none of an "
            "adult's genome is not that adult's child, whichever slot the file was loaded into."),
    Release("Identity by State",
            "Two genomes can show the same allele without either having it from the other: alike in "
            "state, not by descent. A miscall that amplified arrays repeat is alike in state, and a "
            "parent rebuilt from amplified cells is now read without them."),
    Release("Loss of Heterozygosity",
            "Where one parent's copy of a stretch is gone, the cell stops showing that parent's "
            "alleles there. With only the other parent on an array the missing alleles are still "
            "countable, by how often unrelated people carry them, down to where two parents who "
            "share a stretch read the same."),
    Release("Double-Strand Break",
            "A crossover begins with a break through both strands of one chromosome. A break left "
            "unrepaired loses what lies beyond it, and the piece lost can be short. Read marker by "
            "marker against the parent whose copy it was, it shows where a window of two thousand "
            "markers does not."),
    Release("Sister Chromatid",
            "After replication a chromosome is two sister chromatids, each a copy of the other. The "
            "cells of one embryo are copies in the same way: what the parents gave is in every one "
            "of them. So a stretch missing from one cell alone is a loss, and a stretch missing from "
            "all of them is what the parents share."),
    Release("Proofreading",
            "A polymerase checks each base it has just laid down and takes back one that does not "
            "pair, before it goes on. Nothing the tool reads or reports changes here: a check that "
            "failed without cause now measures what it claims to, and the libraries under the "
            "tool are brought up to date."),
    Release("Ideogram",
            "An ideogram is a chromosome drawn to scale, with its landmarks marked along its "
            "length. Every call now has one: the chromosome's own markers from end to end, allele "
            "fraction above and intensity below, with a box over the stretch the call covers."),
    Release("Pedigree",
            "A pedigree is the record of who descends from whom. Given several men's arrays, the "
            "tool now reads which one the cells descend from: the one whose genome is in them. "
            "Where none is, or where two are, it says so and names nobody."),
    Release("Premature Separation",
            "A chromatid that parts from its sister too early can travel to the wrong pole and "
            "leave a cell with an extra copy. Which parent that copy came from is in the alleles. "
            "With one parent on an array, an extra copy the intensity sees is now put to that "
            "question even where the alleles alone read it as a missing parent."),
    Release("Cortical Reaction",
            "When one sperm has fused, the egg's cortical reaction hardens its coat and nothing "
            "else gets in. A file from a chip the tool is not proven on is now stopped at the door "
            "the same way, with the reason, and is never scored with another chip's numbers."),
    Release("Segregation",
            "A parent's two copies of a chromosome segregate among her children, each child "
            "drawing one. Enough draws show what she carried. A run that holds other cells of the "
            "same mother now reads her from them, and with one parent on an array the parent of "
            "an extra chromosome is named as often as with two."),
    Release("Pseudoautosomal",
            "The X and the Y pair and cross over only at their tips, the pseudoautosomal regions, "
            "which behave like any autosome. Everywhere else a cell holds X by copies of its own. "
            "With no parent on an array those copies are now read for what they show, and a line "
            "that depends on how many genomes the cell holds says so."),
    Release("Barr Body",
            "A cell with two X chromosomes shuts one down, and it shows as a dense spot at the "
            "edge of the nucleus, the Barr body. That spot is how two X were first told from "
            "one. With a parent on an array, a second X is now printed only where the cell "
            "shows it."),
    Release("Karyogram",
            "A karyogram lays a cell's chromosomes out in order, so a missing or extra one is "
            "seen at a glance. The result card now draws each parent's copy across the "
            "chromosomes that way, and its headline says what was found."),
    Release("Deletion",
            "A deletion takes part of a chromosome and leaves the rest. A chromosome whose "
            "median the whole-chromosome test pushed over the line was called lost entire, "
            "and the cut was never looked for. It is now read again in pieces."),
    Release("Trophectoderm",
            "The outer cells of a blastocyst, the few a clinic takes for testing and puts on an "
            "array. The examples are now real trophectoderm biopsies and the donors beside them, "
            "read three ways: both parents, one parent, and none."),
    Release("Null Allele",
            "An allele that makes no product at all, and shows only by what is missing. A reading "
            "taken without the reference panel is the same thing: an answer with nothing behind "
            "it. The panel is now read, or the file is refused."),
    Release("Concordance",
            "Two readings of the same cell that agree. A card is one reading made of many parts, "
            "and where its headline, its badge, its notice and the figure under it disagreed, one "
            "of them was wrong. They are made to agree, or the card says it cannot tell."),
    Release("Dosage",
            "How many copies of a stretch of chromosome are present. It is the quantity the "
            "intensity channel reads, and splitting it by parent is what the figure draws."),
    Release("Duplication",
            "An extra copy of a stretch of chromosome. The counterpart of a deletion, and the "
            "harder of the two to see, because an extra copy raises a reading less than a "
            "missing one lowers it."),
)

VERSION = "6.17.0"
CODENAME = "Dosage"

BUILD = f'Build {VERSION} "{CODENAME}"'


def current() -> Release:
    """The Release record for CODENAME. Raises if the codename is not on the ladder."""
    for r in RELEASES:
        if r.name == CODENAME:
            return r
    raise ValueError(f"CODENAME {CODENAME!r} is not in RELEASES - add it or fix the typo")


def as_dict() -> dict:
    r = current()
    return {"version": VERSION, "codename": r.name, "gloss": r.gloss, "build": BUILD}


if __name__ == "__main__":
    r = current()
    assert r.name == CODENAME
    names = [x.name for x in RELEASES]
    assert len(names) == len(set(names)), "duplicate codename on the ladder"
    assert all(x.gloss.strip().endswith(".") for x in RELEASES), "every gloss is a sentence"
    assert BUILD == f'Build {VERSION} "{CODENAME}"'
    d = as_dict()
    assert d["build"] == BUILD and d["gloss"] == r.gloss
    print(f"{BUILD}\n  {r.gloss}\n  ladder: {len(RELEASES)} names, next up: "
          f"{names[names.index(CODENAME) + 1] if names.index(CODENAME) + 1 < len(names) else '(end)'}")

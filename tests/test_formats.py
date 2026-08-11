"""Normalisation: the format families a genotyping lab actually produces, into one shape.

Every case here is modelled on a real public file that was parsed during development, and the
header lines are copied from those files rather than invented.
"""

from __future__ import annotations

import collections
from pathlib import Path

import pytest

from originmarker import origin
def w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_a_genomestudio_preamble_does_not_hide_the_columns(tmp_path):
    """GenomeStudio writes a [Header] block before [Data], so the columns are on line ten. The
    reader takes the first line that RESOLVES to the fields it needs, which handles that marker,
    a comment preamble, and a plain header, without knowing which it is looking at."""
    p = w(tmp_path, "fr.txt", "\n".join([
        "[Header]", "BSGT Version\t3.1.12", "Num SNPs\t561466", "Num Samples\t120", "[Data]",
        "Sample ID\tSNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\t"
        "B Allele Freq\tLog R Ratio",
        "NA06985\trs1000000\t12\t125456933\tT\tC\t0.5416\t0.3096",
        "NA06985\trs10000010\t4\t21227772\tC\tC\t1.0000\t0.1894",
        "NA06989\trs1000000\t12\t125456933\tT\tT\t0.0100\t0.0200",
        "NA06989\trs10000010\t4\t21227772\tC\tT\t0.4900\t0.0100",
    ]))
    got = origin.read_samples(p)
    assert sorted(got) == ["NA06985", "NA06989"], "long format: one row per sample and SNP"
    assert got["NA06985"].probes["rs1000000"].gt == "TC", "nucleotides are kept until harmonised"
    assert got["NA06985"].probes["rs1000000"].baf == pytest.approx(0.5416)
    assert got["NA06985"].probes["rs10000010"].lrr == pytest.approx(0.1894)


def test_alleles_split_across_two_columns_are_joined(tmp_path):
    p = w(tmp_path, "a.txt", "\n".join([
        "SNP Name\tChr\tPosition\tAllele1 - Top\tAllele2 - Top",
        "rs1\t1\t100\tA\tG",
        "rs2\t1\t200\t-\t-",
    ]))
    s = origin.read_sample(p)
    assert s.probes["rs1"].gt == "AG"
    assert s.probes["rs2"].gt == "NC", "a dashed no-call is a no-call, not an empty genotype"


def test_columns_named_only_for_the_samples_are_recognised(tmp_path):
    """HapMap and several public releases write one column per sample, headed by nothing but the
    sample id. Measured: 174 samples out of one chr22 file."""
    p = w(tmp_path, "hm.txt", "\n".join([
        "rs# alleles chrom pos strand NA12891 NA12892 NA12878",
        "rs1 A/G chr22 100 + AA AG AG",
        "rs2 C/T chr22 200 + CC CT CC",
    ]))
    got = origin.read_samples(p)
    assert sorted(got) == ["NA12878", "NA12891", "NA12892"]
    assert got["NA12892"].probes["rs1"].gt == "AG"
    assert got["NA12891"].probes["rs2"].chrom == "22", "the chr prefix is stripped"


def test_nucleotides_are_harmonised_across_samples_never_within_one(tmp_path):
    """A marker whose alleles are A and G shows only "AA" in a sample homozygous for the first,
    so that sample alone cannot know G exists. Two files mapping independently would disagree
    about which allele is A, so the mapping is agreed across every sample being compared."""
    p = w(tmp_path, "n.txt", "\n".join([
        "rs# alleles chrom pos DAD KID",
        "rs1 A/G 1 100 AA GG",
        "rs2 C/T 1 200 CC CT",
        "rs3 A/G 1 300 GG GG",
    ]))
    got = origin.read_samples(p)
    dad, kid = origin.harmonise([got["DAD"], got["KID"]])
    # rs1: alleles {A,G} pooled across BOTH samples, so the father is one homozygote and the
    # child the other. Harmonising the father alone would have made them both "AA".
    assert (dad.probes["rs1"].gt, kid.probes["rs1"].gt) == ("AA", "BB")
    assert (dad.probes["rs2"].gt, kid.probes["rs2"].gt) == ("AA", "AB")
    assert (dad.probes["rs3"].gt, kid.probes["rs3"].gt) == ("AA", "AA"), "one allele: both hom"
    solo, = origin.harmonise([got["DAD"]])
    assert solo.probes["rs1"].gt == "AA", "alone, the father cannot see the other allele"


def test_harmonising_ab_space_data_changes_nothing(tmp_path):
    """It has to be safe to run on everything, so the caller never has to know which space a file
    is in. In AB space the alleles are A and B and sorted order maps each to itself."""
    p = w(tmp_path, "ab.txt", "\n".join([
        "Name,Chr,Position,s1.GType,s2.GType",
        "rs1,1,100,AA,AB",
        "rs2,1,200,BB,AB",
    ]))
    got = origin.read_samples(p)
    before = {k: {r: x.gt for r, x in v.probes.items()} for k, v in got.items()}
    after = {s.sample_id: {r: x.gt for r, x in s.probes.items()}
             for s in origin.harmonise(list(got.values()))}
    assert after == before


def test_a_file_with_markers_but_no_genotypes_says_which_layouts_it_looked_for(tmp_path):
    p = w(tmp_path, "x.txt", "Name,Chr,Position,B Allele Freq\nrs1,1,100,0.5\n")
    with pytest.raises(ValueError, match="no genotypes"):
        origin.read_samples(p)


def test_delimiters_are_sniffed_including_whitespace(tmp_path):
    for name, sep in (("t.txt", "\t"), ("c.txt", ","), ("s.txt", " ")):
        p = w(tmp_path, name, sep.join(["Name", "Chr", "Position", "GType"]) + "\n"
              + sep.join(["rs1", "1", "100", "AB"]) + "\n")
        assert origin.read_sample(p).probes["rs1"].gt == "AB", name


def test_one_sample_can_be_pulled_out_of_a_multi_sample_file(tmp_path):
    p = w(tmp_path, "m.txt", "\n".join([
        "rs# chrom pos A B C", "rs1 1 100 AA AG GG"]))
    assert origin.read_sample(p, "B").probes["rs1"].gt == "AG"
    with pytest.raises(ValueError, match="no sample"):
        origin.read_sample(p, "Z")

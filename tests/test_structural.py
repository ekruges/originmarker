"""PLINK input, both-parent classification, and the read-level evidence an array cannot carry."""

from __future__ import annotations

from pathlib import Path

import pytest

from originmarker import origin
from originmarker import structural


# --- PLINK ------------------------------------------------------------------------------------


def write_bed(tmp_path, ids, snps, codes, mode=1):
    (tmp_path / "t.bim").write_text(
        "\n".join(f"{c}\t{r}\t0\t{p}\t{a1}\t{a2}" for r, c, p, a1, a2 in snps) + "\n")
    (tmp_path / "t.fam").write_text("\n".join(f"F {i} 0 0 1 -9" for i in ids) + "\n")
    buf = bytearray(b"\x6c\x1b" + bytes([mode]))
    for row in codes:
        packed = bytearray((len(ids) + 3) // 4)
        for j, c in enumerate(row):
            packed[j >> 2] |= (c & 3) << ((j & 3) * 2)
        buf += packed
    (tmp_path / "t.bed").write_bytes(bytes(buf))
    return tmp_path / "t"


SNPS = [("rs1", "1", 100, "A", "G"), ("rs2", "1", 200, "C", "T"), ("rs3", "2", 300, "A", "T")]


def test_plink_binary_decodes_the_documented_two_bit_codes(tmp_path):
    """00 homozygous A1, 01 MISSING, 10 heterozygous, 11 homozygous A2. The missing code sits in
    the middle of the run, so reading it as a genotype invents a call for every failed sample."""
    pre = write_bed(tmp_path, ["DAD", "MUM", "KID"], SNPS,
                    [[0, 2, 3], [0, 2, 1], [2, 2, 2]])
    got = origin.read_plink(pre)
    assert sorted(got) == ["DAD", "KID", "MUM"]
    assert [got[s].probes["rs1"].gt for s in ("DAD", "MUM", "KID")] == ["AA", "AB", "BB"]
    assert got["KID"].probes["rs2"].gt == "NC", "code 01 is missing, not a genotype"
    assert (got["KID"].probes["rs3"].chrom, got["KID"].probes["rs3"].pos) == ("2", 300)


def test_the_ab_mapping_comes_from_the_bim_rather_than_being_pooled(tmp_path):
    """`.bim` names A1 and A2 per marker, so the mapping is exact here. Everywhere else it has to
    be agreed across samples, because a homozygous sample shows only one of the two alleles."""
    pre = write_bed(tmp_path, ["ONE"], SNPS, [[0], [3], [2]])
    s = origin.read_plink(pre)["ONE"]
    assert [s.probes[r].gt for r in ("rs1", "rs2", "rs3")] == ["AA", "BB", "AB"]


def test_individual_major_plink_is_refused_rather_than_misread(tmp_path):
    pre = write_bed(tmp_path, ["A", "B"], SNPS, [[0, 0], [0, 0], [0, 0]], mode=0)
    with pytest.raises(ValueError, match="individual-major"):
        origin.read_plink(pre)


def test_a_bad_magic_number_is_refused(tmp_path):
    write_bed(tmp_path, ["A"], SNPS, [[0], [0], [0]])
    (tmp_path / "t.bed").write_bytes(b"\x00\x00\x01" + b"\x00" * 3)
    with pytest.raises(ValueError, match="magic bytes"):
        origin.read_plink(tmp_path / "t")


def test_plink_text_is_read_and_harmonised(tmp_path):
    (tmp_path / "t.map").write_text("1\trs1\t0\t100\n1\trs2\t0\t200\n")
    (tmp_path / "t.ped").write_text("F DAD 0 0 1 -9 A A C C\nF KID 0 0 2 -9 G G C T\n")
    got = origin.read_plink(tmp_path / "t")
    # .ped carries nucleotides with no A1/A2 column, so the mapping is pooled across samples.
    assert (got["DAD"].probes["rs1"].gt, got["KID"].probes["rs1"].gt) == ("AA", "BB")
    assert (got["DAD"].probes["rs2"].gt, got["KID"].probes["rs2"].gt) == ("AA", "AB")


def test_read_sample_dispatches_on_a_plink_extension(tmp_path):
    write_bed(tmp_path, ["ONE"], SNPS, [[0], [3], [2]])
    assert origin.read_sample(tmp_path / "t.bed", "ONE").probes["rs2"].gt == "BB"


# --- both parents -----------------------------------------------------------------------------


def trio(n=4_000, *, pat_absent=0.002, mat_absent=0.002, parent_het=0.17):
    """A trio where both parents carry the heterozygosity a real person does.

    `parent_het` is load-bearing rather than decorative: the second-parent axis is derived from
    it, so a parent homozygous at every marker drives the expectation to zero and every sample
    reads biparental whatever it contains. Each parent's heterozygous markers are assigned
    independently of `scatter`, because carving them out of the region where the child carries
    both alleles would shrink the informative denominator without touching the numerator and
    inflate the absence rate by 1/(1 - parent_het).
    """
    fa, mo, kid = {}, {}, {}
    for i in range(n):
        k, pos = f"1:{i}", 1_000 + i * 1_000
        scatter = (i * 7919) % 1000
        f_het = ((i * 2654435761) >> 16) % 1000 < parent_het * 1000
        m_het = ((i * 2246822519) >> 13) % 1000 < parent_het * 1000
        fa[k] = origin.Probe("1", pos, "AB" if f_het else "AA", baf=0.5 if f_het else 0.0)
        mo[k] = origin.Probe("1", pos, "AB" if m_het else "BB", baf=0.5 if m_het else 1.0)
        # The child should be AB. Losing one parent's allele shows as that parent's homozygote.
        gt = "AB"
        if scatter < pat_absent * 1000:
            gt = "BB"
        elif scatter >= 1000 - mat_absent * 1000:
            gt = "AA"
        # A real diploid genome reads 15%-16% in the heterozygous BAF band. Setting every
        # marker to 0.5 gives 100%, which no genome produces and which the plausibility gate in
        # parental_origin now refuses, correctly. The band is a separate channel from the
        # genotype here, as it is for the parents above.
        in_band = ((i * 2246822519) >> 11) % 1000 < 160
        kid[k] = origin.Probe("1", pos, gt, baf=0.5 if in_band else 1.0)
    return (origin.Sample("dad", fa, []), origin.Sample("mum", mo, []),
            origin.Sample("kid", kid, []))


def test_both_parents_classify_from_two_measurements_not_one_and_an_inference():
    """With only the father, whether a second parent contributed is read off the rate of alleles
    he lacks, which separates by about 1.6x. With the mother present each parent is tested on the
    axis that separates by thirty-fold, and the class falls out of the pair."""
    dad, mum, kid = trio()
    r = origin.both_parents(dad, mum, kid)
    assert r.origin_class == "biparental"
    assert r.paternal.verdict == r.maternal.verdict == "parent_genome_present"


def test_a_missing_maternal_complement_is_confirmed_rather_than_inferred():
    dad, mum, kid = trio(mat_absent=0.10)
    r = origin.both_parents(dad, mum, kid)
    assert r.origin_class == "androgenetic"
    assert any("confirms it directly" in n for n in r.notes)


def test_neither_parent_accounting_for_a_genome_points_at_the_pairing_first():
    dad, mum, kid = trio(pat_absent=0.10, mat_absent=0.10)
    r = origin.both_parents(dad, mum, kid)
    assert r.origin_class == "neither_parent"
    assert any("check the pairing" in n for n in r.notes)


def test_two_declared_parents_that_are_one_person_are_flagged():
    """Both tests pass against the same genome, so it reads as a confident biparental result."""
    dad, _mum, kid = trio()
    r = origin.both_parents(dad, dad, kid)
    assert any("duplicate or a relabelled file" in n for n in r.notes)


# --- structural evidence ------------------------------------------------------------------------


SV_VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "6\t65399000\tdel1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=65402000;SVLEN=-3000\tGT:SR\t0/1:12,7",
    "6\t65400500\tins1\tN\t<INS>\t50\tPASS\tSVTYPE=INS;SVLEN=340\tGT:SR\t0/1:9,4",
    "6\t70000000\tdel2\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=70010000\tGT\t0/0",
    "6\t65401000\tlow\tN\t<DEL>\t50\tLowQual\tSVTYPE=DEL;END=65402000\tGT\t1/1",
])


def test_structural_calls_carry_the_read_support_behind_them(tmp_path):
    p = tmp_path / "sv.vcf"
    p.write_text(SV_VCF + "\n")
    svs = structural.read_sv_vcf(p)
    d = next(v for v in svs if v.variant_id == "del1")
    assert (d.svtype, d.span, d.genotype, d.split_reads) == ("DEL", 3000, "AB", 7)


def test_a_call_in_the_file_is_not_a_call_in_the_sample(tmp_path):
    p = tmp_path / "sv.vcf"
    p.write_text(SV_VCF + "\n")
    ev = structural.at_locus(structural.read_sv_vcf(p), "6", 65_400_000)
    assert [v.variant_id for v in ev.deletions] == ["del1"], "0/0 is not carried, LowQual fails"
    assert [v.variant_id for v in ev.insertions] == ["ins1"]


def test_absence_of_a_call_is_not_reported_as_absence_of_an_event(tmp_path):
    """A caller run with a size threshold above the event, or without the type enabled, gives the
    same empty result as a genome with nothing there."""
    p = tmp_path / "sv.vcf"
    p.write_text(SV_VCF + "\n")
    ev = structural.at_locus(structural.read_sv_vcf(p), "6", 90_000_000)
    assert not ev.spans_locus
    assert any("not absence of an event" in n for n in ev.notes)
    assert "excludes an event the caller could detect" in ev.summary()


def test_no_structural_input_at_all_says_the_mechanisms_are_untested():
    ev = structural.at_locus([], "6", 65_400_000, supplied=False)
    assert not ev.supplied and "NOT TESTED" in ev.summary()


def test_reads_make_the_two_mechanisms_an_array_cannot_see_testable(tmp_path):
    """An array cannot see an insertion at all, and cannot see a deletion shorter than r_min
    markers however clean the data. Reads have neither limit."""
    from originmarker import normalize as nz
    p = tmp_path / "sv.vcf"
    p.write_text(SV_VCF + "\n")
    svs = structural.read_sv_vcf(p)
    cal = nz.Calibration(center=0.0, lrr_compression=0.594, sigma_lrr=0.11, route="declared")
    ms = [origin.Marker(f"m{i}", "6", 65_390_000 + i * 2_000, "AA", "BB", "AB", baf=0.5, lrr=0.0)
          for i in range(120)]

    def ledger(ev):
        v = origin.analyse_embryo("E", ms, "6", 65_400_000, calibration=cal, dropout=0.01,
                                  structural_evidence=ev)
        return {x.mechanism: x.verdict for x in v.windows[-1].ledger}

    bare = ledger(None)
    assert bare["H3f insertion or balanced rearrangement"] == "not_tested"
    with_reads = ledger(structural.at_locus(svs, "6", 65_400_000))
    assert with_reads["H3f insertion or balanced rearrangement"] == "not_excluded"
    assert with_reads["H3b large on-target deletion"] == "excluded"

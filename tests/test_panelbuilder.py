"""
Golden integration test for panelbuilder (ABCC8 c.3989-9G>A).

Offline and deterministic: PANELBUILDER_CACHE points at fixtures/, the recorded API
responses from a real run. A single live smoke test is gated behind RUN_LIVE=1.

    pytest test_panelbuilder.py
    RUN_LIVE=1 pytest test_panelbuilder.py::test_live_smoke   # optional, hits APIs
"""
import io
import json
import os
import random
import urllib.error
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
os.environ.setdefault("PANELBUILDER_CACHE", str(FIXTURES))
# 0 = never expire. The fixtures are the recorded source of truth, not a cache; letting
# them age out would quietly turn this offline suite into a live one.
os.environ.setdefault("PANELBUILDER_CACHE_TTL", "0")

import panelbuilder as pb  # noqa: E402  (env set before import)

GOLDEN_QUERY = "NM_000352.6(ABCC8):c.3989-9G>A"


@pytest.fixture(scope="module")
def result():
    if not FIXTURES.exists():
        pytest.skip("fixtures/ not present - run `python panelbuilder.py` once, then "
                    "copy .panelbuilder_cache -> fixtures/ to record them")
    return pb.build(GOLDEN_QUERY)


def test_variant_identity(result):
    v = result.variant
    assert v.rsid == "rs151344623"
    assert v.gene == "ABCC8"
    assert v.chrom == "11"
    assert v.pos_grch38 == 17_397_055
    assert (v.vcf_ref, v.vcf_alt) == ("C", "T")
    assert v.strand == -1
    assert v.clinical_significance == "Pathogenic"
    assert v.clinvar_accession == "VCV000009088"


@pytest.mark.parametrize("prose", [
    "ABCC8 splice mutation",
    "the cystic fibrosis one",
    "a pathogenic variant",
    "BRCA1",
    "pathogenic",
    "",
    "   ",
])
def test_prose_never_resolves_to_a_variant(prose):
    """Free text must not become a locus.

    ClinVar's esearch is a full-text search, not a lookup: it does not fail on prose, it
    returns its best match, and a coordinate for a variant nobody asked about is worse
    than an error because an error is visible. The gate rejects before any request.
    """
    with pytest.raises(pb.ApiError, match="not a variant identifier"):
        pb.resolve_variant(prose)


def test_identifier_shapes_are_accepted(result):
    """The gate must not be so tight that it refuses real input."""
    for ok in ["rs151344623", "NM_000352.6(ABCC8):c.3989-9G>A", "NM_000352.6:c.3989-9G>A",
               "VCV000009088", "NC_000011.10:g.17397055C>T"]:
        assert pb._looks_like_identifier(ok), ok


def test_record_must_match_the_query():
    """A hit that is not the variant asked for is refused, not reported."""
    vset = {"variation_xrefs": [{"db_source": "dbSNP", "db_id": "999"}]}
    summ = {"title": "NM_000000.1(WRONG):c.1A>T", "accession": "VCV000000001"}
    # rsID that the record does not carry
    with pytest.raises(pb.ApiError, match="was not asked for"):
        pb._assert_record_matches("rs151344623", summ, vset)
    # HGVS naming a different variant
    with pytest.raises(pb.ApiError, match="not the same variant"):
        pb._assert_record_matches("NM_000352.6(ABCC8):c.3989-9G>A", summ, vset)
    # the real pairing passes, with and without the gene parenthetical
    good_v = {"variation_xrefs": [{"db_source": "dbSNP", "db_id": "151344623"}]}
    good_s = {"title": "NM_000352.6(ABCC8):c.3989-9G>A", "accession": "VCV000009088"}
    pb._assert_record_matches("rs151344623", good_s, good_v)
    pb._assert_record_matches("NM_000352.6(ABCC8):c.3989-9G>A", good_s, good_v)
    pb._assert_record_matches("NM_000352.6:c.3989-9G>A", good_s, good_v)
    pb._assert_record_matches("VCV000009088", good_s, good_v)


def test_strand_aware_transcript_change(result):
    # minus strand: genomic C>T == transcript-sense G>A
    assert result.variant.transcript_sense_change().startswith("G>A")


def test_variant_is_not_its_own_marker():
    """A variant does not flank itself.

    The golden case is blind to this: its AF never clears the MAF floor, so it can never
    be recruited into its own panel. A common pathogenic allele can, and lands at dist 0,
    on a side it is not on, counted in neither lower_count nor higher_count.
    """
    v = pb.VariantRecord(query="q", rsid="rsSELF", gene="G", strand=1, chrom="1",
                         pos_grch38=1000, vcf_ref="C", vcf_alt="T",
                         clinical_significance="Pathogenic", review_status=None,
                         clinvar_accession=None)
    common = {"genome": {"af": 0.34, "an": 100_000, "populations": []}}
    raw = {
        # the variant itself, common enough to pass the MAF floor
        "self": {"variant_id": "1-1000-C-T", "rsid": "rsSELF", "pos": 1000,
                 "ref": "C", "alt": "T", **common},
        # a genuine flanking marker
        "flank": {"variant_id": "1-1500-A-G", "rsid": "rsFLANK", "pos": 1500,
                  "ref": "A", "alt": "G", **common},
    }
    markers = pb.annotate(raw, v)
    rsids = {m.rsid for m in markers}
    assert "rsSELF" not in rsids, "the variant was recruited as a marker of itself"
    assert "rsFLANK" in rsids, "a real flanking marker was dropped"
    assert all(m.dist != 0 for m in markers)


def test_gene_pick_ignores_alphabetical_first_gene():
    """ClinVar orders genes[] alphabetically, so genes[0] is not the curated gene.

    The record below is the real shape of VCV000002556: the LOC locus sorts first and sits
    on the opposite strand, so picking it flips the strand as well as the gene.
    """
    summ = {"title": "NM_000243.3(MEFV):c.2078TGA[1] (p.Met694del)",
            "gene_sort": "LOC126862264",
            "genes": [{"symbol": "LOC126862264", "strand": "+"},
                      {"symbol": "MEFV", "strand": "-"}]}
    assert pb._gene_and_strand(summ, {}) == ("MEFV", -1)

    # single-gene record: strand still comes off the ClinVar entry, no Ensembl call
    hexa = {"title": "NM_000520.6(HEXA):c.1274_1277dup (p.Tyr427fs)",
            "genes": [{"symbol": "HEXA", "strand": "-"}]}
    assert pb._gene_and_strand(hexa, {}) == ("HEXA", -1)


def test_revcomp_handles_indels():
    # A multi-base allele must reverse complement, not complement base by base.
    assert pb._revcomp("GGATA") == "TATCC"
    assert pb._revcomp("C") == "G"          # reversal is a no-op for an SNV
    assert pb._revcomp("") is None
    assert pb._revcomp("<DEL>") is None     # symbolic allele has no transcript sense

    v = pb.VariantRecord(query="VCV000003889", rsid=None, gene="HEXA", strand=-1,
                         chrom="15", pos_grch38=72346579, vcf_ref="G", vcf_alt="GGATA",
                         clinical_significance="Pathogenic", review_status=None,
                         clinvar_accession="VCV000003889")
    assert v.transcript_sense_change() == "C>TATCC (transcript sense; minus strand)"


def test_qc_failed_sites_are_never_offered():
    """gnomAD's own QC verdict is honoured: a site flagged AC0 (no high-quality genotype
    survived) or AS_VQSR (looks like an artifact) is not offered as a marker.

    Not fastidiousness: a mismapped site collects spurious heterozygotes, that inflates
    2pq, and 2pq is the ranking key, so QC-rejected sites sort to the top of the panel.
    """
    v = pb.VariantRecord(query="q", rsid="rsX", gene="G", strand=1, chrom="1",
                         pos_grch38=1000, vcf_ref="C", vcf_alt="T",
                         clinical_significance=None, review_status=None,
                         clinvar_accession=None)
    common = {"af": 0.5, "an": 100_000, "populations": [{"id": "nfe", "ac": 5000, "an": 10_000}]}
    raw = {
        "pass": {"variant_id": "1-2000-A-G", "rsid": "rsPASS", "pos": 2000, "ref": "A",
                 "alt": "G", "genome": dict(common, filters=[])},
        "ac0": {"variant_id": "1-3000-A-G", "rsid": "rsAC0", "pos": 3000, "ref": "A",
                "alt": "G", "genome": dict(common, filters=["AC0"])},
        "vqsr": {"variant_id": "1-4000-A-G", "rsid": "rsVQSR", "pos": 4000, "ref": "A",
                 "alt": "G", "genome": dict(common, filters=["AS_VQSR"])},
    }
    got = {m.rsid for m in pb.annotate(raw, v)}
    assert got == {"rsPASS"}, got


def test_thin_population_frequencies_are_not_reported():
    """A frequency from four people is not a frequency.

    Per-population MAF needs an AN floor, or an AN of 8 yields a MAF that looks like any
    other and can be selected on by the ancestry ranking. AN is carried alongside so a
    reader can see what each figure rests on.
    """
    maf, ans = pb._per_pop_maf({"genome": {"populations": [
        {"id": "nfe", "ac": 5_000, "an": 10_000},   # solid
        {"id": "mid", "ac": 4, "an": 8},            # four people
    ]}})
    assert "NFE" in maf and maf["NFE"] == 0.5
    assert "MID" not in maf, "a MAF from 4 people was reported"
    assert ans["NFE"] == 10_000


def test_failed_1000g_join_is_unknown_not_zero():
    """A failed allele match is not a count of zero.

    Ensembl reports indel alleles in its own notation ('-', 'A') where ClinVar and gnomAD
    say 'C'>'CA'. No string comparison bridges those, so the join misses routinely and the
    miss must read as unavailable, never as a measured count.
    """
    import unittest.mock as mock
    v = pb.VariantRecord(query="x", rsid="rs999", gene="G", strand=1, chrom="1",
                         pos_grch38=1, vcf_ref="C", vcf_alt="CA",
                         clinical_significance=None, review_status=None,
                         clinvar_accession=None)
    ens = json.dumps({"name": "rs999", "mappings": [], "populations": [
        {"population": "1000GENOMES:phase_3:ALL", "allele": "-", "allele_count": 1200},
        {"population": "1000GENOMES:phase_3:ALL", "allele": "A", "allele_count": 3800}]})
    gnom = {"data": {"variant": {"genome": {"ac": 50, "an": 100_000, "af": 5e-4}}}}
    with mock.patch.object(pb, "_graphql", return_value=gnom), \
         mock.patch.object(pb, "_get", return_value=ens):
        r = pb.assess_rarity(v)
    assert r.thousand_genomes_ac is None, "a failed join was reported as a count"
    # The verdict is still soundly "too rare" on gnomAD's own AF. What must not happen is
    # the missing 1000G figure rendering as a measured zero.
    assert "count 0" not in r.reason, "a failed join printed as a count of zero"
    assert "unavailable" in r.reason, "the missing figure is not disclosed"
    assert r.population_LD_usable is False


def test_rarity_is_never_asserted_from_a_gap():
    """A source that did not answer is not a source that answered "rare".

    Mocked rather than golden, and this is why the mock table has to exist: both of
    ABCC8's sources answer and both say rare, so the fixture only ever takes the one
    branch that is always correct. The defect class lives in the branches it never
    reaches.
    """
    import unittest.mock as mock
    v = pb.VariantRecord(query="x", rsid="rs1801133", gene="MTHFR", strand=1, chrom="1",
                         pos_grch38=11_796_321, vcf_ref="C", vcf_alt="CA",
                         clinical_significance=None, review_status=None,
                         clinvar_accession=None)

    def kg(*pairs):
        return json.dumps({"name": "rs1801133", "populations": [
            {"population": "1000GENOMES:phase_3:ALL", "allele": a, "allele_count": c}
            for a, c in pairs]})

    def rarity(ens, gnom):
        """ens/gnom None means that source raised: an outage rather than an answer."""
        boom = mock.Mock(side_effect=Exception("outage"))
        with mock.patch.object(pb, "_get", mock.Mock(return_value=ens) if ens else boom), \
             mock.patch.object(pb, "_graphql", mock.Mock(return_value=gnom) if gnom else boom):
            return pb.assess_rarity(v)

    common = {"data": {"variant": {"genome": {"ac": 1377, "an": 5008, "af": 0.275}}}}
    indel_miss = kg(("-", 1200), ("A", 3800))    # notation mismatch: joins to nothing
    seen = []

    # The two ways the 1000G count goes missing on a common SNP. The second needs no
    # outage: it fires on healthy infrastructure, every time, for this indel.
    for label, ens in (("ensembl outage", None), ("indel join failure", indel_miss)):
        r = rarity(ens, common)
        seen.append(r)
        assert r.thousand_genomes_ac is None, label
        assert r.ld_status == "unknown", (label, r.ld_status)
        assert "too rare" not in r.reason.lower(), (label, r.reason)
        # One source quiet, the other answered 0.275: the verdict must disclose the figure
        # it does have, rather than claim it could not retrieve any.
        assert "2.75e-1" in r.reason, (label, r.reason)
        assert "could not be retrieved" not in r.reason, (label, r.reason)
        assert "unavailable" in r.reason, (label, r.reason)

    # Both sources silent: this one really could not look, and is allowed to say so.
    r = rarity(None, None)
    seen.append(r)
    assert r.ld_status == "unknown", r.ld_status
    assert "too rare" not in r.reason.lower(), r.reason
    assert "could not be retrieved" in r.reason, r.reason

    # The mirror image: a source that actually returned a low figure still reads as too
    # rare from either side alone. A true rarity finding is not dropped because the other
    # source was quiet.
    for label, ens, gnom in (
            ("gnomAD alone says rare", indel_miss, {"data": {"variant": {
                "genome": {"ac": 50, "an": 100_000, "af": 5e-4}}}}),
            ("1000G alone says rare", kg(("CA", 1)), None)):
        r = rarity(ens, gnom)
        seen.append(r)
        assert r.ld_status == "undefined", (label, r.ld_status)
        assert r.population_LD_usable is False, label
        assert "too rare" in r.reason.lower(), (label, r.reason)

    # Confirmed common: 1000G answered with a real count, above the threshold.
    r = rarity(kg(("CA", 1377)), common)
    seen.append(r)
    assert r.ld_status == "defined" and r.population_LD_usable is True, r.ld_status

    for r in seen:
        # An Optional interpolated into prose reads to a scientist as a measurement, so
        # every figure goes through fmt_af or an explicit "unavailable", on every branch.
        assert "None" not in r.reason, r.reason
        # The verdict and its evidence, and nothing else. It renders on a one-line card, so
        # a branch that grows a standing caveat back into itself is a bug: R3 rides on
        # DISCLAIMER, which every surface serving this string also serves.
        assert len(r.reason) <= 160, f"{len(r.reason)} chars is too long for the card: {r.reason}"
        assert r.ld_status in ("defined", "undefined", "unknown"), r.ld_status
        # LD is never asserted as usable without a confirmed count behind it.
        assert r.population_LD_usable is (r.ld_status == "defined")


def test_classification_is_the_aggregate_not_one_rcv():
    """One condition's opinion is not the variant's classification.

    A VCV carries a GermlineClassification per RCV (variant+condition pairing) and those
    disagree by design, so only the variation-level aggregate speaks for the variant: not
    the first <Description> in document order. The review status must come from the same
    block as the description it is paired with, or the two can describe different records.
    """
    xml = """<ClinVarResult-Set><VariationArchive><ClassifiedRecord>
      <Classifications>
        <GermlineClassification>
          <ReviewStatus>reviewed by expert panel</ReviewStatus>
          <Description>drug response</Description>
        </GermlineClassification>
      </Classifications>
      <RCVList>
        <RCVAccession><RCVClassifications><GermlineClassification>
          <ReviewStatus>no assertion criteria provided</ReviewStatus>
          <Description>Conflicting classifications of pathogenicity</Description>
        </GermlineClassification></RCVClassifications></RCVAccession>
        <RCVAccession><RCVClassifications><GermlineClassification>
          <ReviewStatus>criteria provided, single submitter</ReviewStatus>
          <Description>Benign</Description>
        </GermlineClassification></RCVClassifications></RCVAccession>
      </RCVList>
    </ClassifiedRecord></VariationArchive></ClinVarResult-Set>"""
    import xml.etree.ElementTree as ET
    sig, rev = pb._aggregate_classification(ET.fromstring(xml))
    assert sig == "drug response", sig
    assert rev == "reviewed by expert panel", rev

    # No aggregate block: say nothing rather than promote an RCV's opinion.
    bare = "<ClinVarResult-Set><VariationArchive><ClassifiedRecord><RCVList>" \
           "<RCVAccession><RCVClassifications><GermlineClassification>" \
           "<Description>Pathogenic</Description></GermlineClassification>" \
           "</RCVClassifications></RCVAccession></RCVList></ClassifiedRecord>" \
           "</VariationArchive></ClinVarResult-Set>"
    assert pb._aggregate_classification(ET.fromstring(bare)) == (None, None)


def test_rarity_forces_phasing(result):
    r = result.rarity
    assert r.population_LD_usable is False      # single copy in 1000G
    assert r.thousand_genomes_ac in (0, 1)
    assert r.gnomad_af_genome is not None and r.gnomad_af_genome < 1e-3


def test_candidate_pool(result):
    assert len(result.candidates) > 1000
    for m in result.candidates:
        assert m.maf >= pb.COMMON_MAF
        assert len(m.ref) == 1 and len(m.alt) == 1


def test_nearest_marker_is_rs757110(result):
    nearest = min(result.candidates, key=lambda m: abs(m.dist))
    assert nearest.rsid == "rs757110"
    assert nearest.pos == 17_396_930
    assert nearest.dist == -125


def test_panel_covers_both_sides(result):
    tel = [m for m in result.recommended if m.dist < 0]
    cen = [m for m in result.recommended if m.dist > 0]
    assert len(tel) >= 2 and len(cen) >= 2
    # Both sides carry close markers, so no thin-coverage flag. The hotspots this panel
    # does cross are a separate flag, asserted below.
    assert not [f for f in result.coverage["flags"] if "30 kb" in f or "No markers" in f]


def _marker(pos, het, *, rsid=None, variant_id=None, hot=False, approx=False,
            het_max_pop=None):
    """A marker at `pos` for a variant at 1000. Only the selection inputs are meaningful.

    het_max_pop defaults to a value deliberately different from het, because a max over
    populations is >= the global 2pq and strictly greater in real data. Defaulting it to
    het would make the two quantities equal by construction, and no test built on this
    factory could tell them apart. Pass it explicitly when a test cares about the value.
    """
    dist = pos - 1000
    hmp = het_max_pop if het_max_pop is not None else min(0.5, round(het * 1.2 + 0.04, 4))
    return pb.Marker(
        rsid=rsid or f"rs{pos}", variant_id=variant_id or f"1-{pos}-A-G", chrom="1",
        pos=pos, ref="A", alt="G", af=0.5, maf=0.5, het=het, het_max_pop=hmp, dist=dist,
        side="higher coord" if dist > 0 else "lower coord", tier=pb._tier(dist),
        cm=0.0, recomb_fraction=0.0, hotspot_between=hot, map_approx=approx)


def test_hotspot_markers_are_flagged_not_dropped(result):
    """A flanking marker is only useful while it stays in phase with the locus, and a
    recombination hotspot between the two is the one thing that breaks that. Such markers
    are higher risk, not worthless: they stay in the panel with their per-marker theta,
    and the coverage judgment has to say how many there are.
    """
    hot = [m for m in result.recommended if m.hotspot_between]
    assert hot, "golden panel crosses hotspots on the higher-coordinate side"
    assert all(m.recomb_fraction is not None for m in hot)
    flag = [f for f in result.coverage["flags"] if "hotspot" in f]
    assert flag, result.coverage["flags"]
    assert f"{len(hot)} of the {len(result.recommended)}" in flag[0]


def test_close_markers_behind_a_hotspot_do_not_count_as_cover():
    """Six markers within 30 kb, every one of them across a hotspot, is not coverage."""
    markers = ([_marker(1000 - d, 0.5, hot=True) for d in (500, 1_500, 5_000, 20_000)]
               + [_marker(1000 + d, 0.5) for d in (500, 1_500, 5_000, 20_000)])
    _, recommended, cov = pb.select_panel(markers)
    assert cov["lower_core_near"] == 4          # they are still in the panel
    assert len([m for m in recommended if m.dist < 0]) == 4
    flags = " ".join(cov["flags"])
    assert "lower-coordinate side are clear of an intervening recombination hotspot" in flags
    assert "(0 of 4)" in flags
    assert "higher-coordinate" not in flags   # the clear side is not accused


def test_unmapped_chromosome_is_not_reported_as_hotspot_free():
    """No map is not the same claim as no hotspot.

    The 1 cM/Mb fallback can never reach the hotspot threshold, so hotspot_between is
    False for every marker on a chromosome with no bundled map. Counting that as clear
    turns "not assessed" into a clean bill of health.
    """
    markers = [_marker(1000 + d, 0.5, approx=True) for d in
               (-20_000, -1_000, 1_000, 20_000)]
    _, _, cov = pb.select_panel(markers)
    flags = " ".join(cov["flags"])
    assert "not assessed for 4 of the 4" in flags, cov["flags"]
    assert "clear of an intervening recombination hotspot" not in flags


def test_empty_side_is_not_the_same_as_a_thin_side():
    """Zero markers on a side cannot detect recombination there at all; one can."""
    one_side = [_marker(1000 - d, 0.5) for d in (500, 1_500)]
    _, _, cov = pb.select_panel(one_side)
    assert any("No markers at all on the higher-coordinate side" in f
               for f in cov["flags"]), cov["flags"]

    thin = one_side + [_marker(1000 + 500, 0.5)]
    _, _, cov = pb.select_panel(thin)
    assert any(f == "Fewer than 2 markers within 30 kb on the higher-coordinate side."
               for f in cov["flags"]), cov["flags"]
    assert not any("No markers at all" in f for f in cov["flags"])


def test_multiallelic_site_keeps_the_better_allele():
    """One assay reads a site once.

    Multi-allelic sites share an rsID, so selection keyed on rsid lets one alt evict the
    other and gnomAD's response order decides the survivor.
    """
    better = _marker(1500, 0.49, rsid="rsMULTI", variant_id="1-1500-A-G")
    worse = _marker(1500, 0.10, rsid="rsMULTI", variant_id="1-1500-A-C")
    filler = [_marker(1000 - d, 0.5) for d in (500, 1_500)]
    for order in ([better, worse], [worse, better]):
        _, recommended, _ = pb.select_panel(order + filler)
        picked = [m for m in recommended if m.rsid == "rsMULTI"]
        assert len(picked) == 1, "one assay reads the site once"
        assert picked[0].variant_id == "1-1500-A-G", "the worse allele won the slot"


def test_selection_does_not_depend_on_api_response_order(result):
    """Two labs running the same query must get the same panel back."""
    markers = list(result.candidates)
    expected = [m.variant_id for m in pb.select_panel(markers)[1]]
    for seed in range(10):
        random.Random(seed).shuffle(markers)
        assert [m.variant_id for m in pb.select_panel(markers)[1]] == expected


def test_monomorphic_sites_are_not_markers():
    """common_maf=0.0 is a legal floor, so 2pq=0 sites reach the pool. A site heterozygous
    in nobody cannot show which parental haplotype an embryo inherited, and must not be
    offered as a marker or counted as cover for a side.
    """
    real = [_marker(1000 - d, 0.5) for d in (500, 1_500)]
    # het_max_pop=0.0 explicitly: "heterozygous in nobody" means zero in EVERY
    # population. The factory's default is deliberately nonzero, so omitting it would
    # build a site that is monomorphic globally and polymorphic somewhere: the test below.
    dead = [_marker(1000 + d, 0.0, het_max_pop=0.0) for d in (500, 1_500, 5_000)]
    candidates, recommended, cov = pb.select_panel(real + dead)
    assert all(m.het > 0 for m in candidates + recommended)
    assert cov["higher_count"] == 0
    assert any("No markers at all on the higher-coordinate side" in f
               for f in cov["flags"]), cov["flags"]


def test_globally_monomorphic_but_polymorphic_somewhere_is_kept():
    """The complement of the test above, and the reason the guard reads max() rather than
    het: a site with a global 2pq of 0 can still be common in one population, where it is
    informative for a family. "Informative for nobody" is the bar for exclusion, not
    "informative for the average of everybody".
    """
    live = [_marker(1000 - 500, 0.0, het_max_pop=0.42),
            _marker(1000 + 500, 0.0, het_max_pop=0.38)]
    candidates, recommended, _ = pb.select_panel(live)
    assert len(candidates) == 2, "a site common in one population is not a dead site"
    assert recommended, "and it is eligible for the shortlist"


def test_ensembl_cross_check_passes(result):
    checked = [m for m in result.recommended if m.ensembl_pos_check is not None]
    assert checked, "expected some markers cross-checked"
    assert all(m.ensembl_pos_check == "ok" for m in checked)


def test_ranking_never_uses_ld(result):
    # Sanity: recommended markers are high-het and near, not LD-selected.
    near = [m for m in result.recommended if abs(m.dist) < 30_000]
    assert max(m.het_max_pop for m in near) > 0.45


def test_provenance_and_disclaimer(result):
    p = result.provenance
    assert p["build"] == "GRCh38"
    assert "Research use only" in p["disclaimer"]
    assert "phasing" in p["disclaimer"].lower()


@pytest.mark.skipif(os.environ.get("RUN_LIVE") != "1", reason="live APIs; set RUN_LIVE=1")
def test_live_smoke(tmp_path):
    os.environ["PANELBUILDER_CACHE"] = str(tmp_path)   # fresh cache -> real calls
    r = pb.build(GOLDEN_QUERY)
    assert r.variant.rsid == "rs151344623"
    assert r.rarity.population_LD_usable is False


# --------------------------------------------------------------------------- #
# resolve_variant: esearch is a relevance ranking, not a lookup
#
# The golden fixture's esearch returns count=1, so it exercises the reconcile loop with a
# single element and is blind to everything below. These drive the loop with a synthetic
# esearch/esummary pair instead: offline, no fixture, no network.
# --------------------------------------------------------------------------- #

def _fake_clinvar(monkeypatch, ids, records, seen=None):
    """Stub _get for the esearch -> esummary pair; efetch (clinical detail) is optional."""
    def fake_get(base, path, params=None, **kw):
        if params and seen is not None:
            seen.append(dict(params))
        if "esearch" in path:
            return json.dumps({"esearchresult": {"count": str(len(ids)),
                                                 "idlist": ids[:int(params.get("retmax", 20))]}})
        if "esummary" in path:
            return json.dumps({"result": {i: records[i] for i in params["id"].split(",")
                                          if i in records}})
        raise pb.ApiError("efetch not stubbed")   # best-effort path, caught by resolve
    monkeypatch.setattr(pb, "_get", fake_get)


def _cv_record(title, accession, rsid_digits):
    return {"title": title, "accession": accession, "obj_type": "single nucleotide variant",
            "genes": [{"symbol": "HBB", "strand": "-"}],
            "variation_set": [{"variation_name": title,
                               "variation_xrefs": [{"db_source": "dbSNP",
                                                    "db_id": rsid_digits}],
                               "variation_loc": [{"assembly_name": "GRCh38", "chr": "11",
                                                  "start": "5227002", "ref": "T",
                                                  "alt": "A"}]}]}


SICKLE = "NM_000518.5(HBB):c.20A>T"


def test_reconciles_against_every_id_esearch_returned():
    """The matching record may sit anywhere in the relevance order, including last.

    A prefix of a relevance ranking is not a shortlist: relevance does not rank "is this
    the variant you named". Truncating the candidate ids refuses common indications on the
    exact expression ClinVar itself displays.
    """
    seen = []
    ids = [str(9000 + i) for i in range(12)] + ["15333"]
    records = {i: _cv_record(f"NM_000518.5(HBB):c.{380 + n}_396del (p.Val127fs)",
                             f"VCV00000{i}", "999")
               for n, i in enumerate(ids[:12])}
    records["15333"] = _cv_record(f"{SICKLE} (p.Glu7Val)", "VCV000015333", "334")

    with pytest.MonkeyPatch.context() as mp:
        _fake_clinvar(mp, ids, records, seen)
        v = pb.resolve_variant(SICKLE)

    # Only the last record carries rs334; every decoy carries rs999.
    assert v.rsid == "rs334", v
    assert (v.chrom, v.pos_grch38) == ("11", 5_227_002)
    assert v.gene == "HBB" and v.strand == -1
    # retmax must be explicit: the esearch default is 20, so relying on it silently caps
    # the candidate set at a number nobody chose.
    assert int(seen[0]["retmax"]) >= len(ids), seen[0]


def test_unmatched_search_does_not_blame_the_identifier():
    """No reconciling record is a statement about the SEARCH, not about the user.

    Re-raising the first candidate's mismatch names whatever esearch ranked first as a
    best match and tells the user their canonical expression is suspect. Both are false.
    """
    ids = [str(9000 + i) for i in range(3)]
    records = {i: _cv_record("NM_000518.5(HBB):c.380_396del (p.Val127fs)",
                             f"VCV00000{i}", "999") for i in ids}

    with pytest.MonkeyPatch.context() as mp:
        _fake_clinvar(mp, ids, records)
        with pytest.raises(pb.ApiError) as ei:
            pb.resolve_variant(SICKLE)

    msg = str(ei.value)
    assert "Check the identifier" not in msg, msg      # it was not the user's mistake
    assert "best match" not in msg, msg                # do not dignify an unrelated record
    assert "VCV000009000" not in msg, msg              # ... or name it at all
    assert "3 record(s)" in msg, msg
    # and it must still refuse, rather than build a panel on the near miss
    assert "is this variant" in msg, msg


def test_ensembl_outage_is_not_reported_as_absence():
    """A 503 means the question was never asked, so it cannot answer it.

    Ensembl signals genuine absence with 400 and {"error": "<id> not found ..."}, NOT with
    404. Both directions are pinned here: keying the absence claim on 404 alone reports
    every real typo as an outage, which is the same bug reversed.
    """
    def raiser(msg):
        def fake_get(*a, **kw):
            raise pb.ApiError(msg)
        return fake_get

    outage = "GET https://rest.ensembl.org/variation/homo_sapiens/rs334 failed after 2 tries: HTTP Error 503: Service Unavailable"
    absent = "GET https://rest.ensembl.org/variation/homo_sapiens/rs9 failed after 2 tries: HTTP Error 400: Bad Request"

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(pb, "_get", raiser(outage))
        with pytest.raises(pb.ApiError) as ei:
            pb._resolve_via_ensembl("rs334")
    msg = str(ei.value)
    assert "not found" not in msg, msg                 # no claim about existence
    assert "Check the identifier" not in msg, msg      # no blame
    assert "503" in msg and "Retry" in msg, msg
    assert "NOT a statement about whether" in msg, msg

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(pb, "_get", raiser(absent))
        with pytest.raises(pb.ApiError, match="was not found"):
            pb._resolve_via_ensembl("rs9")             # 400 = Ensembl's own "no"

    # a timeout carries no status at all -> still an outage, never an absence
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(pb, "_get", raiser("GET https://rest.ensembl.org/... failed after 2 tries: timed out"))
        with pytest.raises(pb.ApiError, match="could not be reached"):
            pb._resolve_via_ensembl("rs334")


def test_ensembl_status_survives_the_trip_through_http(tmp_path):
    """The status discriminator must survive _http's real message, not a hand-written one.

    _resolve_via_ensembl recovers the HTTP status by regexing 'HTTP Error (\\d{3})' back
    out of the prose _http formats with {last}. Nothing makes _http promise to keep
    saying that: str(HTTPError) is 'HTTP Error 400: Bad Request' where repr is
    '<HTTPError 400: ...>', so {last} -> {last!r} would drop the status to None and turn
    every genuine typo into "could not be reached, retry". Patching at urlopen, the real
    network boundary, lets the real _http build the message: that is what pins the two
    functions together.
    """
    def urlopen_raising(code, reason):
        def fake(req, timeout=None):
            raise urllib.error.HTTPError(req.full_url, code, reason, {}, io.BytesIO(b""))
        return fake

    for code, reason in [(400, "Bad Request"), (404, "Not Found")]:
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(pb, "CACHE_DIR", tmp_path)          # never read the fixtures
            mp.setattr(pb.time, "sleep", lambda *_: None)  # don't pay the retry backoff
            mp.setattr(pb.urllib.request, "urlopen", urlopen_raising(code, reason))
            with pytest.raises(pb.ApiError) as ei:
                pb._resolve_via_ensembl("rs99999999999")
        assert "was not found" in str(ei.value), f"{code} lost its status: {ei.value}"

    for code, reason in [(500, "Internal Server Error"), (503, "Service Unavailable")]:
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(pb, "CACHE_DIR", tmp_path)
            mp.setattr(pb.time, "sleep", lambda *_: None)
            mp.setattr(pb.urllib.request, "urlopen", urlopen_raising(code, reason))
            with pytest.raises(pb.ApiError) as ei:
                pb._resolve_via_ensembl("rs334")
        msg = str(ei.value)
        assert "could not be reached" in msg, msg
        assert "not found" not in msg, msg
        assert str(code) in msg, msg      # the real status reached the sentence


# --------------------------------------------------------------------------- #
# The fetch ledger is per-build, not per-process
# --------------------------------------------------------------------------- #

def test_concurrent_builds_do_not_share_a_data_date():
    """The fetch ledger is per-build, not per-process.

    app/jobs.py runs MAX_CONCURRENT builds in one process, and queried_utc is the age of
    the data on a filed artifact. A module global would be thread-safe and build-unsafe:
    a lock cannot help, because the builds do not race on the dict, they race for it.
    """
    import concurrent.futures as cf
    import threading

    barrier = threading.Barrier(2)

    def one(epochs):
        pb._reset_fetch_log()
        barrier.wait()                      # force the interleave the race needs
        for e in epochs:
            pb._note_fetch(e, cached=True)
        barrier.wait()
        log = pb._fetch_log()
        return pb._utc(log["oldest"]), log["from_cache"]

    old = [1_600_000_000 + i for i in range(5)]   # 2020 sources
    new = [1_780_000_000 + i for i in range(5)]   # 2026 sources
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        a, b = list(ex.map(one, [old, new]))

    assert a[0].startswith("2020"), f"build A took the other build's data date: {a}"
    assert b[0].startswith("2026"), f"build B took the other build's data date: {b}"
    assert a[1] == b[1] == 5, f"counts leaked across builds: {a} {b}"


def test_pool_workers_tally_into_the_builds_ledger():
    """The counterpart to the test above, and the reason the ledger is a ContextVar rather
    than a thread-local: _note_fetch runs inside the enumerate_candidates worker pool,
    so a thread-local would hand every worker its own empty ledger and build() would read
    back zeros. Isolation is worthless if it isolates the build from its own data.
    """
    import contextvars
    import concurrent.futures as cf

    pb._reset_fetch_log()
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        # copy_context() must be evaluated on the submitting thread, as
        # enumerate_candidates does it: inside the worker it would copy the worker's empty
        # context and silently tally nowhere.
        futs = [ex.submit(contextvars.copy_context().run, pb._note_fetch, 1_700_000_000, False)
                for _ in range(8)]
        [f.result() for f in futs]
    assert pb._fetch_log()["from_network"] == 8, pb._fetch_log()


def test_efetch_failure_is_not_diagnosed_as_a_structural_variant(monkeypatch):
    """A failed lookup must never be reported as a property of the variant.

    efetch supplies the VCF alleles. If a failure falls through to the structural-variant
    guard, that guard reads the empty alleles and refuses a routine SNV as having no
    single reference and alternate allele.
    """
    real_get = pb._get

    def broken_efetch(base, path, params=None, **kw):
        if "efetch" in path:
            raise pb.ApiError("HTTP Error 502: Bad Gateway")
        return real_get(base, path, params, **kw)

    monkeypatch.setattr(pb, "_get", broken_efetch)
    with pytest.raises(pb.ApiError) as ei:
        pb.resolve_variant("NM_000352.6(ABCC8):c.3989-9G>A", build="GRCh38")
    msg = str(ei.value)
    assert "no single reference and alternate allele" not in msg, \
        f"an efetch outage is still being diagnosed as a CNV: {msg}"
    assert "not a statement about the variant" in msg, msg
    assert "502" in msg, "the actual cause should reach the person who has to retry"


def test_rsid_resolves_to_the_allele_not_a_haplotype_containing_it():
    """The rsID and the HGVS for one allele must agree. They are the same variant.

    ClinVar's relevance order returns haplotype records ahead of the allele itself, and a
    haplotype record still resolves: it hands back one constituent's position with the
    haplotype's own classification attached, so the right coordinate arrives carrying the
    wrong significance.
    """
    by_rsid = pb.resolve_variant("rs334", build="GRCh38")
    by_hgvs = pb.resolve_variant("NM_000518.5(HBB):c.20A>T", build="GRCh38")

    assert by_rsid.clinvar_accession == "VCV000015333", \
        f"rs334 resolved to {by_rsid.clinvar_accession} ({by_rsid.clinical_significance!r})"
    assert by_rsid.clinical_significance == "Pathogenic"
    assert (by_rsid.chrom, by_rsid.pos_grch38) == (by_hgvs.chrom, by_hgvs.pos_grch38)
    assert by_rsid.clinvar_accession == by_hgvs.clinvar_accession, \
        "the rsID and the HGVS for one allele disagree about which record it is"


def test_combination_records_are_recognised():
    """The guard is on obj_type, not on a hardcoded accession list."""
    assert pb._is_combination({"obj_type": "Haplotype"})
    assert pb._is_combination({"obj_type": "CompoundHeterozygote"})
    assert pb._is_combination({"obj_type": "Diplotype"})
    assert not pb._is_combination({"obj_type": "single nucleotide variant"})
    assert not pb._is_combination({"obj_type": "Deletion"})
    assert not pb._is_combination({})


def test_position_disagreement_reaches_the_coverage_card(monkeypatch):
    """Two sources disagreeing about where a marker IS is not a cell value.

    The cross-check costs a request per marker, so it is spent on the shortlist and runs
    after selection: a disputed marker is already shortlisted when the disagreement
    surfaces, and the panel has to say so rather than leave it in its own row. Injected
    because the golden panel cross-checks clean, so building it would pass either way.
    """
    real = pb.cross_check_ensembl

    def one_disputed(markers, top_n=8):
        real(markers, top_n=top_n)
        markers[0].ensembl_pos_check = "MISMATCH:99999999"

    monkeypatch.setattr(pb, "cross_check_ensembl", one_disputed)
    r = pb.build(pb.StructuredQuery(variant=GOLDEN_QUERY, build="GRCh38"))

    flags = r.coverage["flags"]
    hit = [f for f in flags if "placed differently" in f]
    assert hit, f"a disputed marker produced no panel-level warning: {flags}"
    assert "99,999,999" in hit[0] or "99999999" in hit[0], hit[0]
    assert r.recommended[0].rsid in hit[0], "the flag should name the marker to distrust"


def test_clean_cross_check_produces_no_disagreement_flag():
    """The other half: the flag must not cry wolf on a panel that checks out."""
    r = pb.build(pb.StructuredQuery(variant=GOLDEN_QUERY, build="GRCh38"))
    assert not [f for f in r.coverage["flags"] if "placed differently" in f]


def test_a_deletion_overhanging_the_template_is_still_masked():
    """A deletion is anchored at the base BEFORE the bases it removes.

    So one straddling the template's left edge has pos < region_start while its reference
    footprint still lies under the forward primer. Filtering mask sites on the anchor
    dropped it, primer3 was never told, and the note told the reader the region was clear.
    primers.design already clamps an overhanging site; the caller could simply never hand
    it one, so the clamp's own check exercised a path production could not reach.
    """
    import primers

    flank = 400
    marker_pos = 17_397_055
    lo, hi = marker_pos - flank, marker_pos + flank

    # anchored 4 bp before the template, deleting 12 bp: covers the first 8 template bases
    overhang = primers.MaskSite(pos=lo - 4, ref="A" * 13, alt="A", maf=0.30)
    inside = primers.MaskSite(pos=marker_pos + 120, ref="C", alt="T", maf=0.30)
    hazards = [overhang, inside]

    near = [s for s in hazards if s.pos + s.span - 1 >= lo and s.pos <= hi]
    assert overhang in near, "a deletion overhanging the template edge was dropped"
    assert inside in near

    # and the anchor test, which is what shipped, drops it
    anchored = [s for s in hazards if lo <= s.pos <= hi]
    assert overhang not in anchored, "this test no longer reproduces the bug it pins"

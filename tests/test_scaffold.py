"""Which paternal homologue each embryo received, from co-inheritance between siblings."""

from __future__ import annotations

import math

import pytest

from originmarker import scaffold
def trio_markers(n=40, *, embryos=("A", "B", "C"), carries=None, start=17_200_000, step=5_000,
                 mother="AA", flip_after=None):
    """Markers with a heterozygous father and a homozygous mother, so every allele is deducible.

    `carries` maps embryo -> which of the father's homologues it takes (0 or 1). `flip_after`
    maps embryo -> marker index beyond which it switches, standing in for a crossover.
    """
    carries = carries or {"A": 0, "B": 1, "C": 0}
    hap = [("A", "B") if i % 3 else ("B", "A") for i in range(n)]
    out = []
    for i in range(n):
        gts = {}
        for e in embryos:
            side = carries[e]
            if flip_after and e in flip_after and i > flip_after[e]:
                side ^= 1
            gts[e] = "".join(sorted(hap[i][side] + mother[0]))
        out.append(scaffold.Marker(f"m{i}", "11", start + i * step, "AB", mother, gts))
    return out


def test_the_paternal_allele_is_deducible_only_from_a_heterozygous_father():
    """The exact complement of the presence test in `origin`, which needs him homozygous. A
    homozygous father transmits a known allele but his two chromosomes are indistinguishable
    through it, so it says nothing about WHICH one came through."""
    assert scaffold.paternal_allele("AB", "AA", "AB") == "B"
    assert scaffold.paternal_allele("AB", "AA", "AA") == "A"
    assert scaffold.paternal_allele("AA", "AB", "AB") is None
    assert scaffold.paternal_allele("AB", "AB", "AB") is None, "both het: cannot say"
    assert scaffold.paternal_allele("AB", "AB", "AA") == "A", "a homozygous embryo resolves"
    assert scaffold.paternal_allele("AB", "AA", "NC") is None
    # The mother helps but is not required, which is why this half degrades more gracefully
    # than the presence half does.
    assert scaffold.paternal_allele("AB", None, "AA") == "A"
    assert scaffold.paternal_allele("AB", None, "AB") is None


def test_embryos_split_into_two_groups_without_knowing_which_is_which():
    s = scaffold.build_scaffold(trio_markers(), "11", 17_300_000, dropout=0.01)
    assert {c.embryo_id: c.relationship for c in s.calls} == {
        "A": "same", "B": "opposite", "C": "same"}
    assert any("No anchor" in r for r in s.refusals)
    assert all(c.inherited is None for c in s.calls), "unlabelled groups stay unlabelled"


def test_one_external_anchor_names_both_groups_at_once():
    s = scaffold.build_scaffold(trio_markers(), "11", 17_300_000,
                                anchor="A", anchor_carries="mutant", dropout=0.01)
    assert {c.embryo_id: c.inherited for c in s.calls} == {
        "A": "H2_mutant_inherited",
        "B": "H1_wildtype_inherited",
        "C": "H2_mutant_inherited"}
    # And the opposite anchor inverts every call, which is what makes the anchor load-bearing.
    t = scaffold.build_scaffold(trio_markers(), "11", 17_300_000,
                                anchor="A", anchor_carries="wildtype", dropout=0.01)
    assert {c.embryo_id: c.inherited for c in t.calls} == {
        "A": "H1_wildtype_inherited",
        "B": "H2_mutant_inherited",
        "C": "H1_wildtype_inherited"}


def test_confidence_saturates_with_the_map_rather_than_the_marker_count():
    """The regression guard for the flaw this module was rewritten to fix.

    Markers on the same side of a crossover all report one fact, not many. Summing their
    log-odds multiplies a single observation by the marker count and manufactures certainty:
    559 real markers inside 1 cM produced odds near 10^50 for a claim whose true error rate is
    set by whether a crossover fell in between, which no marker density can reduce.
    """
    sparse = scaffold.build_scaffold(trio_markers(10, step=20_000), "11", 17_300_000,
                                     anchor="A", anchor_carries="mutant", dropout=0.01)
    dense = scaffold.build_scaffold(trio_markers(200, step=1_000), "11", 17_300_000,
                                    anchor="A", anchor_carries="mutant", dropout=0.01)
    lo = {c.embryo_id: abs(c.log10_odds) for c in sparse.calls if c.embryo_id != "A"}
    hi = {c.embryo_id: abs(c.log10_odds) for c in dense.calls if c.embryo_id != "A"}
    assert dense.informative_markers == 20 * sparse.informative_markers
    for e in lo:
        assert hi[e] < 3.0 * lo[e], (
            f"{e}: twentyfold the markers moved the odds from 10^{lo[e]:.1f} to 10^{hi[e]:.1f}, "
            "which is marker count masquerading as evidence")


def test_a_crossover_moves_the_answer_and_the_chain_follows_it():
    """An embryo that switches homologue partway must be called by the state AT THE VARIANT,
    not by whichever side of the crossover happens to carry more markers."""
    # C matches A over the first half and B over the second. The variant sits in the first half.
    ms = trio_markers(40, flip_after={"C": 20})
    early = scaffold.build_scaffold(ms, "11", 17_220_000, anchor="A", anchor_carries="mutant",
                                    dropout=0.01)
    late = scaffold.build_scaffold(ms, "11", 17_380_000, anchor="A", anchor_carries="mutant",
                                   dropout=0.01)
    got_early = {c.embryo_id: c.relationship for c in early.calls}["C"]
    got_late = {c.embryo_id: c.relationship for c in late.calls}["C"]
    assert got_early == "same", got_early
    assert got_late == "opposite", got_late


def test_a_homozygous_father_yields_nothing_however_many_embryos():
    ms = [scaffold.Marker(m.rsid, m.chrom, m.pos, "AA", m.mother, m.embryos)
          for m in trio_markers()]
    s = scaffold.build_scaffold(ms, "11", 17_300_000, anchor="A", anchor_carries="mutant")
    assert s.informative_markers == 0
    assert any("heterozygous" in r for r in s.refusals)
    assert not s.calls


def test_a_single_embryo_cannot_be_scaffolded_at_all():
    ms = trio_markers(embryos=("A",), carries={"A": 0})
    s = scaffold.build_scaffold(ms, "11", 17_300_000, anchor="A", anchor_carries="mutant")
    assert any("Fewer than two embryos" in r for r in s.refusals)
    assert not s.calls


def test_an_anchor_that_is_not_among_the_embryos_is_refused():
    s = scaffold.build_scaffold(trio_markers(), "11", 17_300_000,
                                anchor="Z", anchor_carries="mutant", dropout=0.01)
    assert any("not among the embryos" in r for r in s.refusals)
    assert all(c.inherited is None for c in s.calls)


def test_dropout_raises_the_deduction_error_rather_than_being_ignored():
    """A dropped allele does not merely lose a marker here. A heterozygous embryo reading as
    homozygous gives a confident deduction of the WRONG paternal allele, at about half the
    dropout rate, so the error term has to carry it."""
    clean = scaffold.build_scaffold(trio_markers(), "11", 17_300_000,
                                    anchor="A", anchor_carries="mutant", dropout=0.01)
    noisy = scaffold.build_scaffold(trio_markers(), "11", 17_300_000,
                                    anchor="A", anchor_carries="mutant", dropout=0.40)
    assert noisy.deduction_error == pytest.approx(0.20)
    assert clean.deduction_error == scaffold.DEDUCTION_FLOOR
    sharp = {c.embryo_id: abs(c.log10_odds) for c in clean.calls if c.embryo_id != "A"}
    blunt = {c.embryo_id: abs(c.log10_odds) for c in noisy.calls if c.embryo_id != "A"}
    for e in sharp:
        assert blunt[e] < sharp[e], e


def test_the_two_state_chain_behaves_at_its_limits():
    # No observations anywhere: the posterior must stay at the prior.
    assert scaffold.relationship_posterior([None], [], 0.01) == [0.5]
    assert scaffold.relationship_posterior([], [], 0.01) == []
    # A free crossover between two nodes carries nothing across.
    post = scaffold.relationship_posterior([True, None], [scaffold.flip_probability(0.5)], 0.01)
    assert post[1] == pytest.approx(0.5)
    # No recombination: the observation transfers intact, and ONE marker is worth exactly what
    # one marker is worth. Two independent meioses each mis-deduced at 0.01 agree by accident
    # when both are wrong, so P(agree | same) = 0.99^2 + 0.01^2 and the posterior is that.
    post = scaffold.relationship_posterior([True, None], [0.0], 0.01)
    assert post[1] == pytest.approx(0.99 ** 2 + 0.01 ** 2)
    assert scaffold.flip_probability(0.0) == 0.0
    assert scaffold.flip_probability(0.5) == 0.5


def test_the_map_approximation_is_reported_not_hidden():
    ms = [scaffold.Marker(m.rsid, "Y", m.pos, m.father, m.mother, m.embryos)
          for m in trio_markers()]
    s = scaffold.build_scaffold(ms, "Y", 17_300_000, anchor="A", anchor_carries="mutant")
    assert s.map_approx, "chrY has no bundled map; the reader must be told"
    assert "APPROXIMATE" in s.summary()


def test_a_build_mismatch_stops_the_map_being_read_at_foreign_coordinates():
    """Confirmed necessary on real files: the lab's Axiom exports are GRCh37 and the bundled
    deCODE maps are GRCh38. Reading a real map at another assembly's positions returns a
    plausible number from the wrong locus, which is worse than an approximation that says so.
    There is no liftOver by design: UCSC chain files carry a non-commercial field-of-use
    restriction Apache 2.0 cannot sublicense."""
    ms = trio_markers()
    matched = scaffold.build_scaffold(ms, "11", 17_300_000, anchor="A",
                                      anchor_carries="mutant", dropout=0.01, build="GRCh38")
    foreign = scaffold.build_scaffold(ms, "11", 17_300_000, anchor="A",
                                      anchor_carries="mutant", dropout=0.01, build="GRCh37")
    assert not any("bundled genetic map" in r for r in matched.refusals)
    assert any("was NOT consulted" in r for r in foreign.refusals)
    assert foreign.map_approx and not matched.map_approx
    # The calls survive; it is the stated confidence that must not.
    assert ({c.embryo_id: c.inherited for c in foreign.calls}
            == {c.embryo_id: c.inherited for c in matched.calls})
    assert abs(foreign.calls[1].log10_odds) != abs(matched.calls[1].log10_odds)

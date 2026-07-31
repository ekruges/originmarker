"""The integrator: files to one verdict, with the refusals attached.

The run-length arithmetic here is a second implementation of what `runlength.ts` does, so it is
checked against the same golden fixture. That is the whole point of the duplication being
acceptable: a divergence between the two languages fails a test rather than reaching a result.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path
from typing import Optional, Sequence

import pytest

import normalize as nz
import origin

FIXTURES = Path(__file__).resolve().parent / "fixtures"

#: The two header lines below are copied verbatim from published GEO supplementary files, so they
#: are the format contract rather than a guess at it: GSE148488 (Zuccaro, comma-delimited, the
#: intensity column named `normalized_intensity`) and GSE186407 (Turocy, tab-delimited, `log2R`).
ZUCCARO_HEADER = "probeset_id,copy_number,chr,position,probe_classification,baf,genotype,normalized_intensity"
TUROCY_HEADER = "probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset"


# --- the run statistic, against the same fixture the TypeScript side uses -------------------


def test_run_length_matches_the_golden_fixture():
    """All twelve acceptance cases, so the Python and TypeScript implementations cannot drift."""
    vectors = list(csv.DictReader((FIXTURES / "golden_test_vectors.csv").open()))
    expected = list(csv.DictReader((FIXTURES / "golden_test_expected_outputs.csv").open()))
    assert len(expected) == 12

    for exp in expected:
        rows = [v for v in vectors if v["case"] == exp["case"]]
        assert len(rows) == 20
        l3 = [r for r in rows if r["father_gt"] in ("AA", "BB")]
        flags = [origin.score_paternal(r["father_gt"], r["mother_gt"], r["embryo_gt"]) == 0
                 for r in l3]
        best = cur = 0
        for f in flags:
            cur = cur + 1 if f else 0
            best = max(best, cur)
        n = len(flags)
        assert n == int(exp["n_L3"]), f"n_L3: {exp['case']}"
        assert sum(flags) == int(exp["z_sum"]), f"z_sum: {exp['case']}"
        assert best == int(exp["longest_run"]), f"longest_run: {exp['case']}"
        assert origin.r_min(n, origin.KOTHIYAL_FLOOR) == int(exp["r_min"]), f"r_min: {exp['case']}"
        want_p = float(exp["run_p"])
        if want_p > 0:
            got = origin.run_length_p(best, n, origin.KOTHIYAL_FLOOR)
            assert abs(got - want_p) <= 0.02 * want_p, f"run_p: {exp['case']}"


def test_the_erdos_renyi_form_is_computed_stably():
    """The direct expression cancels to exactly zero at r = 10; this one does not."""
    q = origin.KOTHIYAL_FLOOR
    assert 1 - (1 - q ** 10) ** 1 == 0.0, "the naive form really does return zero"
    assert origin.run_length_p(10, 10, q) == pytest.approx(9.849e-23, abs=1e-25)
    # And where the naive form is fine, the two agree.
    for r, n in ((2, 100), (4, 10)):
        naive = 1 - (1 - q ** r) ** (n - r + 1)
        assert origin.run_length_p(r, n, q) == pytest.approx(naive, abs=1e-12)


def test_q_composes_as_max_not_sum():
    assert origin.q_hat(None) == (origin.KOTHIYAL_FLOOR, "kothiyal_floor")
    assert origin.q_hat(0.001) == (origin.KOTHIYAL_FLOOR, "kothiyal_floor")
    assert origin.q_hat(0.33) == (0.33, "empirical")
    assert origin.q_hat(float("nan"))[0] == origin.KOTHIYAL_FLOOR


def test_paternal_scoring_is_defined_only_where_the_father_is_homozygous():
    assert origin.score_paternal("AA", "BB", "AB") == 1
    assert origin.score_paternal("AA", "BB", "BB") == 0, "paternal allele absent"
    assert origin.score_paternal("BB", "AA", "AA") == 0
    assert origin.score_paternal("AB", "AA", "AA") is None, "father het: cannot say"
    assert origin.score_paternal("AA", "BB", "NC") is None
    assert origin.score_paternal("NC", "BB", "AB") is None
    # And it survives the mother being unavailable, which is why the degraded mode is worth running.
    assert origin.score_paternal("AA", None, "BB") == 0


def test_presence_needs_a_mother_who_could_not_have_supplied_the_allele():
    """The defect the PennCNV trio exposed. Absence is Mendelian and needs no mother; presence is
    an identity claim and does. Scoring presence regardless breaks a real run at every marker
    where a heterozygous mother transmitted the same allele the father carries, which turned a
    whole chromosome of paternal absence into a longest run of 21."""
    assert origin.score_paternal("AA", "BB", "AB") == 1, "mother cannot supply A: it is paternal"
    assert origin.score_paternal("AA", "AB", "AB") is None, "mother could have supplied that A"
    assert origin.score_paternal("AA", "AA", "AA") is None
    assert origin.score_paternal("AA", None, "AB") is None, "no mother: presence unprovable"
    # Absence is unaffected by any of it.
    for mother in ("BB", "AB", "AA", None, "NC"):
        assert origin.score_paternal("AA", mother, "BB") == 0, mother


def test_a_marker_that_cannot_say_is_excluded_not_counted_as_presence():
    """Folding None into the flag list both breaks runs and inflates n. A hemizygous embryo whose
    mother is heterozygous produces exactly this marker, so it is the common case, not an edge."""
    ms = []
    for i in range(60):
        # A clean paternal deletion, but the mother is heterozygous at every third marker, so the
        # embryo's hemizygous call there sometimes still contains the father's allele.
        mother = "AB" if i % 3 == 0 else "BB"
        embryo = "AA" if (i % 3 == 0 and i % 2 == 0) else "BB"
        ms.append(origin.Marker(f"m{i}", "11", 1_000_000 + i * 2_000, "AA", mother, embryo,
                                baf=1.0, lrr=-0.594))
    v = origin.analyse_embryo("E", ms, "11", 1_500_000, calibration=GOOD, dropout=0.01)
    w = v.windows[-1]
    assert w.n_l3 == w.longest_run, "every informative marker shows absence, so the run is all of them"
    assert w.n_l3 < len(ms), "the markers that cannot say are excluded from n, not scored"


# --- dropout, and the mistake worth pinning -------------------------------------------------


def test_dropout_is_recovered_exactly_and_needs_no_population_prior():
    n = 10_000
    fa = {f"r{i}": "AB" for i in range(n)}
    mo = {f"r{i}": "AA" for i in range(n)}
    # 18.5% heterozygous where Mendel predicts 50%: d = 1 - 2(0.185) = 0.63.
    emb = {f"r{i}": ("AB" if i < 1_850 else "AA") for i in range(n)}
    d, se, k = origin.estimate_dropout(fa, mo, emb)
    assert d == pytest.approx(0.63, abs=1e-9)
    assert k == n
    assert se < 0.02, "must resolve better than the 2 pp that inverts count statistics"


def test_handing_the_estimator_a_parent_returns_a_clamp_not_a_measurement():
    """A mistake I made and reported as a negative control. The father IS heterozygous at
    father-het markers, so h = 1.0 and 1 - 2h is negative; the clamp returns 0.0, which is not
    evidence of anything. Pinned so it cannot be mistaken for a validation again."""
    fa = {f"r{i}": "AB" for i in range(5_000)}
    mo = {f"r{i}": "AA" for i in range(5_000)}
    d, _, _ = origin.estimate_dropout(fa, mo, dict(fa))
    assert d == 0.0
    assert origin.estimate_dropout({}, {}, {}) is None


# --- the verdict, and its two structural refusals -------------------------------------------


def markers(n: int, *, absent_from: int = -1, absent_to: int = -1) -> list[origin.Marker]:
    out = []
    for i in range(n):
        absent = absent_from <= i <= absent_to
        out.append(origin.Marker(
            rsid=f"m{i}", chrom="11", pos=1_000_000 + i * 2_000,
            father="AA", mother="BB",
            embryo="BB" if absent else "AB",
            baf=1.0 if absent else 0.5,
            lrr=-0.6 if absent else 0.0,
        ))
    return out


GOOD = nz.Calibration(center=0.0, lrr_compression=0.598, sigma_lrr=0.417, route="chrx_male")
UNSCALED = nz.Calibration(center=0.0, lrr_compression=math.nan, sigma_lrr=1.363, route="none")


def test_h1_versus_h2_is_refused_without_a_phase_source():
    """The honest boundary of the whole tool. Those two differ ONLY in which paternal homologue
    was transmitted, and no genotype file says which of the father's chromosomes carried the
    mutation."""
    v = origin.analyse_embryo("E1", markers(60), "11", 1_060_000, calibration=GOOD, dropout=0.01)
    assert any("REFUSED" in r or "refused" in r.lower() for r in v.refusals)
    assert any("differ ONLY in which paternal homologue" in r for r in v.refusals)
    # Declaring one removes that refusal, and only that one.
    v2 = origin.analyse_embryo("E1", markers(60), "11", 1_060_000, calibration=GOOD, dropout=0.01,
                               phase=origin.PhaseSource(mutant_allele="A", route="long_read_father"))
    assert not any("differ ONLY" in r for r in v2.refusals)


def test_an_unscaled_channel_refuses_every_copy_number_mechanism():
    """Real WGA'd embryo behaviour: sigma 1.363 puts the intensity channel past the usable
    ceiling, so copy-number states must be NOT TESTED however clean the genotypes look."""
    v = origin.analyse_embryo("E1", markers(60, absent_from=20, absent_to=40),
                              "11", 1_060_000, calibration=UNSCALED, dropout=0.01)
    assert any("UNSCALED" in r for r in v.refusals)
    ledger = {e.mechanism: e.verdict for e in v.windows[-1].ledger}
    for m in ("H3c segmental paternal loss", "H3d whole-chromosome loss",
              "H3e copy-neutral LOH", "H3f insertion or balanced rearrangement"):
        assert ledger[m] == "not_tested", m


def test_a_real_run_is_found_and_localised_with_a_calibrated_channel():
    v = origin.analyse_embryo("E1", markers(60, absent_from=20, absent_to=40),
                              "11", 1_060_000, calibration=GOOD, dropout=0.01,
                              phase=origin.PhaseSource(mutant_allele="A", route="single_sperm"))
    whole = v.windows[-1]
    assert whole.n_l3 == 59, "the variant position itself is excluded"
    # 21 markers are absent, but one of them IS the variant position, and dropping it from the
    # subsequence joins the two halves rather than splitting the run.
    assert whole.longest_run == 20
    assert whole.r_min is not None and whole.longest_run >= whole.r_min
    assert whole.run_p < 1e-20
    assert whole.segment_state == "PAT0_MAT1", whole.segment_state
    # Separation is quoted over the SEGMENT, which is what makes it callable at all.
    assert whole.separation_sigma > 5.0
    ledger = {e.mechanism: e.verdict for e in whole.ledger}
    assert ledger["H3a stochastic dropout"] == "excluded", "contiguity excludes dropout"
    assert ledger["H3e copy-neutral LOH"] == "excluded"
    # An insertion stays untestable on array input no matter how good the rest is.
    assert ledger["H3f insertion or balanced rearrangement"] == "not_tested"


def test_the_segment_reported_is_the_one_at_the_variant_not_the_longest_anywhere():
    """Two events on one chromosome. Quoting the longer, unrelated one attaches a real sigma to
    the wrong locus, and counting every marker carrying that label anywhere in the window inflates
    it further: a normal chromosome 3 came back as 7,499 markers and 359 sigma."""
    ms = markers(400, absent_from=20, absent_to=30)              # the event at the variant
    for i in range(200, 300):                                    # a longer, unrelated one
        ms[i] = origin.Marker(f"m{i}", "11", 1_000_000 + i * 2_000, "AA", "BB", "BB",
                              baf=1.0, lrr=-0.6)
    v = origin.analyse_embryo("E", ms, "11", 1_050_000, calibration=GOOD, dropout=0.01)
    w = v.windows[-1]
    assert w.segment_markers <= 20, f"reported the far segment: {w.segment_markers} markers"
    assert w.separation_sigma < 30


def test_chrx_copy_states_are_refused_without_a_measured_probe_offset():
    """chrX probes carry a per-array intensity offset unrelated to copy number. Without it a
    male's single X reads as two copies, so the copy layer is withheld while the genotype
    statistics, which never touch intensity, are reported as normal."""
    ms = [origin.Marker(f"x{i}", "X", 1_000_000 + i * 2_000, "AA", "BB", "BB", baf=1.0, lrr=-0.3)
          for i in range(60)]
    no_offset = nz.Calibration(center=0.0, lrr_compression=0.594, sigma_lrr=0.11,
                               route="declared", chrx_offset=None)
    v = origin.analyse_embryo("E", ms, "X", 1_060_000, calibration=no_offset, dropout=0.01)
    w = v.windows[-1]
    assert any("chrX copy-number states NOT TESTED" in r for r in v.refusals)
    assert w.longest_run == w.n_l3 == 59, "the genotype claim still stands"
    assert w.segment_state is None
    ledger = {e.mechanism: e.verdict for e in w.ledger}
    assert ledger["H3e copy-neutral LOH"] == "not_tested"
    # With the offset on file the copy layer runs.
    with_offset = nz.Calibration(center=0.0, lrr_compression=0.594, sigma_lrr=0.11,
                                 route="chrx_male", chrx_offset=0.19)
    v2 = origin.analyse_embryo("E", ms, "X", 1_060_000, calibration=with_offset, dropout=0.01)
    assert not any("chrX copy-number states" in r for r in v2.refusals)
    assert v2.windows[-1].segment_state is not None


def test_a_clean_embryo_claims_nothing():
    v = origin.analyse_embryo("E1", markers(60), "11", 1_060_000, calibration=GOOD, dropout=0.01,
                              phase=origin.PhaseSource(mutant_allele="A", route="single_sperm"))
    whole = v.windows[-1]
    assert whole.longest_run == 0
    assert {e.verdict for e in whole.ledger} <= {"not_excluded", "not_tested"}
    assert whole.segment_state is None, "no candidate event means no state call"


def test_high_dropout_raises_the_bar_rather_than_lowering_confidence():
    """At q = 0.33 a run that would be decisive on clean data is no longer significant. The tool
    must decline rather than report a weaker version of the same claim. Five is the interesting
    length: r_min is 2 on clean data and 8 at 33% dropout, so the same evidence lands either side."""
    clean = origin.analyse_embryo("E", markers(60, absent_from=20, absent_to=24),
                                  "11", 1_060_000, calibration=GOOD, dropout=0.01)
    noisy = origin.analyse_embryo("E", markers(60, absent_from=20, absent_to=24),
                                  "11", 1_060_000, calibration=GOOD, dropout=0.33)
    assert clean.windows[-1].longest_run == noisy.windows[-1].longest_run == 5
    assert clean.windows[-1].r_min < noisy.windows[-1].r_min
    assert clean.windows[-1].longest_run >= clean.windows[-1].r_min
    assert noisy.windows[-1].longest_run < noisy.windows[-1].r_min
    assert {e.mechanism: e.verdict for e in noisy.windows[-1].ledger}[
        "H3a stochastic dropout"] == "not_excluded"


# --- the two gates: is this a real loss, or amplification wearing its shape -------------------


def scored(pattern: Sequence[tuple[int, Optional[int]]]) -> list[tuple[origin.Marker, int]]:
    return [(origin.Marker(f"m{i}", "11", pos, "AA", "BB", "BB"), s)
            for i, (pos, s) in enumerate(pattern)]


def test_a_run_splits_where_the_chromosome_was_never_interrogated():
    """Three violations spread over megabases are three observations, not a deletion. On the real
    trio one such 'run of 3' spanned 5,580 kb, which as a hemizygous deletion in a healthy adult
    would be a major clinical finding rather than a quiet data point."""
    tight = scored([(1_000 + i * 1_000, 0) for i in range(5)])
    assert [r.length for r in origin.absence_runs(tight, 10_000)] == [5]
    spread = scored([(1_000, 0), (2_000, 0), (5_000_000, 0)])
    runs = origin.absence_runs(spread, 10_000)
    assert [r.length for r in runs] == [2, 1], "the far marker is not part of the same event"
    assert (runs[0].lo, runs[0].hi) == (1_000, 2_000)
    # With no cap the old behaviour returns, which is what the constraint exists to prevent.
    assert [r.length for r in origin.absence_runs(spread, None)] == [3]


def test_runs_end_at_the_markers_that_actually_carry_them():
    seq = scored([(1_000, 0), (2_000, 0), (3_000, 1), (4_000, 0), (5_000, None)])
    runs = origin.absence_runs(seq, 10_000)
    assert [(r.length, r.lo, r.hi) for r in runs] == [(2, 1_000, 2_000), (1, 4_000, 4_000)]
    assert origin.absence_runs(scored([(1_000, 1)]), 10_000) == []


def test_the_maternal_side_separates_dropout_from_deletion():
    """A deletion removes one parent's copy. Dropout removes alleles without caring whose."""
    # Violations at roughly the dropout rate: dropout explains it, a deletion cannot.
    assert origin.bidirectional_test(60, 20, 0.33)[0] == "artefact"
    # None at all, with enough markers to have seen some: only a real loss looks like this.
    v, lr = origin.bidirectional_test(30, 0, 0.33)
    assert v == "one_directional" and lr < -1
    # Too few markers to distinguish, and no markers at all.
    assert origin.bidirectional_test(2, 0, 0.33)[0] == "insufficient"
    assert origin.bidirectional_test(0, 0, 0.33) == ("insufficient", 0.0)


def test_the_test_has_no_power_on_clean_dna_and_that_is_correct():
    """When dropout falls to the error floor the two hypotheses coincide, so the ratio goes to
    one and the answer is 'cannot say'. A sample with no excess dropout cannot be faked by
    dropout, so there is nothing to rule out. Measured on the HapMap trio at q = 0.0076."""
    v, lr = origin.bidirectional_test(200, 0, origin.KOTHIYAL_FLOOR)
    assert v == "insufficient" and abs(lr) < 1e-9
    assert origin.bidirectional_test(200, 0, 0.0076)[0] == "insufficient"


def interleaved(*, maternal_also_absent: bool) -> list[origin.Marker]:
    """Half the markers build a paternal-absence run; the other half can only speak for the
    mother, because a heterozygous father says nothing about paternal transmission."""
    ms = []
    for i in range(60):
        pos = 1_000_000 + i * 2_000
        if i % 2 == 0:
            ms.append(origin.Marker(f"p{i}", "11", pos, "AA", "BB", "BB", baf=1.0, lrr=-0.6))
        else:
            embryo = "BB" if maternal_also_absent else "AA"
            ms.append(origin.Marker(f"q{i}", "11", pos, "AB", "AA", embryo, baf=1.0, lrr=-0.6))
    return ms


def test_a_dropout_shaped_run_is_refuted_and_stops_driving_the_verdict():
    v = origin.analyse_embryo("E", interleaved(maternal_also_absent=True), "11", 1_059_000,
                              calibration=GOOD, dropout=0.33)
    w = v.windows[-1]
    assert w.longest_run == 30, w.longest_run
    assert w.bidirectional == "artefact", (w.bidirectional, w.bidirectional_log10_lr)
    assert w.bidirectional_log10_lr > 10
    assert any("REFUTED as amplification dropout" in r for r in w.refusals)
    assert w.segment_state is None, "a refuted run must not go on to drive a state call"
    assert {e.mechanism: e.verdict for e in w.ledger}["H3a stochastic dropout"] == "not_excluded"


def test_a_one_directional_run_of_the_same_length_survives():
    # Variant placed at a paternal-informative marker: the father-heterozygous ones carry the
    # maternal evidence but look ordinary to the state layer, which is the point of them.
    v = origin.analyse_embryo("E", interleaved(maternal_also_absent=False), "11", 1_061_000,
                              calibration=GOOD, dropout=0.33)
    w = v.windows[-1]
    assert w.longest_run == 30
    assert w.bidirectional == "one_directional", w.bidirectional
    assert not any("REFUTED" in r for r in w.refusals)
    assert w.segment_state == "PAT0_MAT1"


def test_without_a_mother_no_significance_is_claimed_at_all():
    """Father plus embryo cannot support a run statistic. Presence is unobservable without
    someone who could not have supplied the allele, so the informative set holds absences only,
    nothing can break a run, and the length is a violation count wearing a p-value. On a normal
    chromosome 20 that produced a run of 3 across 35 Mb at p = 2.5e-07."""
    ms = [origin.Marker(f"m{i}", "11", 1_000_000 + i * 2_000, "AA", None,
                        "BB" if 20 <= i <= 40 else "AB") for i in range(60)]
    v = origin.analyse_embryo("E", ms, "11", 1_061_000, calibration=GOOD, dropout=None)
    w = v.windows[-1]
    assert any("No maternal genotypes" in r for r in v.refusals)
    assert math.isnan(w.run_p), "no p-value may be quoted in this mode"
    assert w.q_source == "uncalibrated_no_mother"
    assert w.segment_state is None, "and nothing downstream may run off it"
    assert not w.hypotheses
    # The absences are still counted and shown; only the significance claim is withheld.
    assert w.z_sum == w.n_l3 > 0


def test_whole_chromosome_absence_survives_marker_deserts():
    """The compactness rule splits at centromeres and heterochromatin, so a male's X fragments
    even though it is genuinely absent end to end. Fragmentation describes where markers are,
    not what is missing, so the whole-chromosome claim rests on every informative marker."""
    ms = [origin.Marker(f"m{i}", "11", p, "AA", "BB", "BB", baf=1.0, lrr=-0.6)
          for i, p in enumerate(list(range(1_000_000, 1_060_000, 2_000))
                                + list(range(9_000_000, 9_060_000, 2_000)))]
    v = origin.analyse_embryo("E", ms, "11", 5_000_000, calibration=GOOD, dropout=0.01)
    w = v.windows[-1]
    assert w.longest_run < w.n_l3, "the desert splits the run"
    assert w.z_sum == w.n_l3, "but every informative marker still shows absence"
    assert "H3_paternal_chromosome_absent" in w.hypotheses


# --- parental origin from sperm and sample alone ----------------------------------------------


def sample_pair(rate_by_chrom, *, n_per_chrom=4_000, y_called=False, chroms=None,
                nonpaternal=0.0, het_band=0.0):
    """A father homozygous everywhere, and a sample built to controlled rates on three axes.

    `rate_by_chrom` sets how often the father's obligate allele is missing. `nonpaternal` sets
    how often the sample carries an allele he cannot supply, which is what separates a
    paternal-only genome from a biparental one. `het_band` sets the fraction of BAF in the
    heterozygous band, which says whether the genome is diploid at all.
    """
    chroms = chroms or [str(i) for i in range(1, 23)]
    fp, sp = {}, {}
    for c in chroms:
        r = rate_by_chrom.get(c, rate_by_chrom.get("*", 0.002))
        for i in range(n_per_chrom):
            k = f"{c}:{i}"
            pos = 1_000 + i * 1_000
            # The father is AA, so "BB" loses his allele and "AB" carries one he lacks.
            if (i % 1000) < r * 1000:
                gt = "BB"
            elif (i % 1000) >= 1000 - nonpaternal * 1000:
                gt = "AB"
            else:
                gt = "AA"
            fp[k] = origin.Probe(c, pos, "AA", baf=0.0)
            sp[k] = origin.Probe(c, pos, gt,
                                 baf=0.5 if (i % 1000) < het_band * 1000 else 1.0)
    for i in range(812):
        fp[f"Y:{i}"] = origin.Probe("Y", 1_000 + i * 1_000, "AA")
        sp[f"Y:{i}"] = origin.Probe("Y", 1_000 + i * 1_000, "AA" if y_called else "NC")
    return (origin.Sample("sperm", fp, []), origin.Sample("sample", sp, []))


AUTO_AND_X = [str(i) for i in range(1, 23)] + ["X"]


AXIOM = "UK Biobank Axiom Array"


def test_a_paternal_genome_is_recognised_without_any_mother():
    """The question a lab asks first, and it needs no mother: attributing one allele needs
    someone who could not have supplied it, but measuring how often the father's obligate allele
    is missing does not."""
    fa, s = sample_pair({"*": 0.002})
    rep = origin.parental_origin(fa, s, product=AXIOM)
    assert rep.verdict == "paternal_genome_present"
    assert rep.genome_rate < 0.01
    assert all(c.verdict == "paternal_present" for c in rep.chroms)


def test_an_unrelated_genome_is_recognised_and_its_chromosomes_are_not_dissected():
    """Measured on the real lab data at 6.81%, against 0.82% that its own noise could explain."""
    fa, s = sample_pair({"*": 0.068})
    rep = origin.parental_origin(fa, s, product=AXIOM)
    assert rep.verdict == "no_paternal_contribution"
    assert math.isnan(rep.baseline), "no baseline exists, so none may be quoted"
    assert not rep.chroms, "per-chromosome loss is meaningless without a genome to lose it from"
    assert rep.genome_rate > rep.explainable * origin.ABSENCE_MARGIN


def test_one_lost_chromosome_is_found_against_a_present_genome():
    fa, s = sample_pair({"*": 0.002, "7": 0.055})
    rep = origin.parental_origin(fa, s, product=AXIOM)
    assert rep.verdict == "paternal_genome_present"
    lost = [c for c in rep.chroms if c.verdict == "paternal_absent"]
    assert [c.chrom for c in lost] == ["7"], [(c.chrom, c.rate) for c in rep.chroms]
    assert lost[0].log10_lr < -3


def test_a_male_offspring_loses_the_paternal_X_by_biology_not_by_defect():
    """His father sent a Y instead. Reporting ordinary sex determination as a chromosome loss
    would be a false alarm on roughly half of all samples."""
    fa, s = sample_pair({"*": 0.002, "X": 0.315},
                        chroms=[str(i) for i in range(1, 23)] + ["X"], y_called=True)
    rep = origin.parental_origin(fa, s, product=AXIOM)
    assert rep.y_bearing_sperm
    x = next(c for c in rep.chroms if c.chrom == "X")
    assert x.verdict == "expected_absent", x
    assert "carries the father's Y" in x.note
    # Without a Y the same X rate is a real loss, not sex determination.
    fa2, s2 = sample_pair({"*": 0.002, "X": 0.315},
                          chroms=[str(i) for i in range(1, 23)] + ["X"], y_called=False)
    x2 = next(c for c in origin.parental_origin(fa2, s2, product=AXIOM).chroms if c.chrom == "X")
    assert x2.verdict == "paternal_absent"


def test_the_case_the_chrY_test_cannot_reach_is_called_from_snps_alone():
    """An X-bearing sperm leaves no Y, so chrY cannot separate a full paternal genome from none.
    On the real lab samples that is 04 and 05 (present) against 03 (absent): all three are
    Y-negative and the SNP rate separates them by thirty to forty fold."""
    fa, present = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, y_called=False)
    fa2, absent = sample_pair({"*": 0.068}, chroms=AUTO_AND_X, y_called=False)
    a = origin.parental_origin(fa, present, product=AXIOM)
    b = origin.parental_origin(fa2, absent, product=AXIOM)
    assert not a.y_bearing_sperm and not b.y_bearing_sperm, "neither has a Y to go on"
    assert a.verdict == "paternal_genome_present"
    assert b.verdict == "no_paternal_contribution"
    # The sperm type comes from the chrX SNP rate, so chrY is a cross-check and not the method.
    assert a.sperm_type == "X_bearing"
    assert any("X-bearing sperm" in n for n in a.notes)
    assert any("chrY" in n for n in b.notes)


def test_the_three_axes_separate_androgenetic_gynogenetic_and_biparental():
    """Paternal absence alone cannot distinguish a paternal-only genome from a biparental one:
    the father's alleles are present in both. The second axis, alleles he cannot supply, is what
    separates them, and the BAF band says whether the genome is diploid at all."""
    fa, andro = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, nonpaternal=0.04, het_band=0.02)
    fa2, gyno = sample_pair({"*": 0.068}, chroms=AUTO_AND_X, nonpaternal=0.12, het_band=0.03)
    fa3, bip = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, nonpaternal=0.30, het_band=0.16)

    a = origin.parental_origin(fa, andro, product=AXIOM)
    g = origin.parental_origin(fa2, gyno, product=AXIOM)
    b = origin.parental_origin(fa3, bip, product=AXIOM)

    assert a.origin_class == "androgenetic", (a.nonpaternal_rate, a.het_band)
    assert g.origin_class == "gynogenetic", (g.nonpaternal_rate, g.het_band)
    assert b.origin_class == "biparental", (b.nonpaternal_rate, b.het_band)
    assert a.zygosity == "uniparental_homozygous"
    assert b.zygosity == "diploid"
    # An androgenote and a biparental embryo look identical on the paternal-absence axis, which
    # is exactly why the second axis has to exist.
    assert abs(a.genome_rate - b.genome_rate) < 0.01
    assert b.nonpaternal_rate > 4 * a.nonpaternal_rate


def test_dropout_is_not_mistaken_for_a_uniparental_genome():
    """Dropout smears the heterozygous cluster; it does not empty it. A biparental embryo at 33%
    dropout still reads 22.2% in the BAF band against 15-16% for a clean diploid, while a truly
    uniparental genome reads 1.3-3.4%."""
    fa, noisy = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, nonpaternal=0.30, het_band=0.22)
    r = origin.parental_origin(fa, noisy, product=AXIOM)
    assert r.zygosity == "diploid"
    assert r.origin_class == "biparental"


# --- reading the two real formats -----------------------------------------------------------


def write(tmp_path, name, header, rows):
    p = tmp_path / name
    body = header + "\n" + "\n".join(rows) + "\n"
    if name.endswith(".gz"):
        p.write_bytes(gzip.compress(body.encode()))
    else:
        p.write_text(body)
    return p


def test_both_published_formats_parse_to_the_same_thing(tmp_path):
    """Different delimiter, different column order, different intensity column name, one result."""
    z = write(tmp_path, "z.CEL.txt", ZUCCARO_HEADER, [
        "AX-13216142,2.0,1,86028,Other,0.009244147,0,0.26140523",
        "AX-37361813,2.0,X,727841,NoMinorHom,0.9589331,2,-0.37207379",
    ])
    t = write(tmp_path, "t.CEL.probes.txt.gz", TUROCY_HEADER, [
        "AX-13216142\t1\t86028\t0.26140523\t0.009244147\t2.0\t0\t1",
        "AX-37361813\tX\t727841\t-0.37207379\t0.9589331\t2.0\t2\t1",
    ])
    a, b = origin.read_sample(z), origin.read_sample(t)
    assert a.probes == b.probes
    assert a.probes["AX-13216142"] == origin.Probe(
        chrom="1", pos=86028, gt="AA", baf=pytest.approx(0.009244147),
        lrr=pytest.approx(0.26140523), cn=2.0)
    assert a.sample_id == "z" and b.sample_id == "t"


def test_unparseable_rows_become_nothing_rather_than_something(tmp_path):
    p = write(tmp_path, "x.CEL.txt", ZUCCARO_HEADER, [
        "AX-1,2.0,1,100,Other,0.5,1,0.0",
        "AX-2,2.0,1,NA,Other,0.5,1,0.0",          # no position: dropped entirely
        "AX-3,2.0,1,300,Other,NaN,-1,---",        # unusable BAF/LRR and a no-call genotype
        "AX-4,2.0,1,400,Other,0.5",               # truncated: dropped
    ])
    s = origin.read_sample(p)
    assert set(s.probes) == {"AX-1", "AX-3"}
    assert s.probes["AX-3"] == origin.Probe(chrom="1", pos=300, gt="NC", baf=None, lrr=None, cn=2.0)


def test_a_file_with_no_genotype_column_is_rejected_by_name(tmp_path):
    p = write(tmp_path, "x.txt", "probeset_id,chr,position,baf", ["AX-1,1,100,0.5"])
    with pytest.raises(ValueError, match="no genotypes"):
        origin.read_sample(p)
    with pytest.raises(ValueError, match="empty file"):
        origin.read_sample(write(tmp_path, "e.txt", "", []))


def test_sex_from_the_heterozygosity_RATIO_survives_dropout(tmp_path):
    """The ratio is the point. Dropout deflates chrX and autosomal heterozygosity by the same
    factor, so it cancels, which an absolute chrX-het threshold would not do. Real WGA samples in
    this dataset reach 47% no-call and still call sex correctly."""
    def sample(dropout: float) -> origin.Sample:
        probes = {}
        for i in range(2000):
            lost = (i % 100) < dropout * 100
            probes[f"a{i}"] = origin.Probe("1", i, "AA" if lost else "AB")
            probes[f"x{i}"] = origin.Probe("X", i, "AA" if lost else "AB")
        return origin.Sample("s", probes, [])

    for d in (0.0, 0.3, 0.6):
        assert origin.call_sex(sample(d))[0] == "female", d
    male = origin.Sample("m", {f"a{i}": origin.Probe("1", i, "AB" if i % 2 else "AA")
                               for i in range(500)}
                         | {f"x{i}": origin.Probe("X", i, "AA") for i in range(500)}, [])
    assert origin.call_sex(male)[0] == "male"
    # An intermediate ratio names neither, because a wrong male call silently doubles the scale.
    ambiguous = origin.Sample("q", {f"a{i}": origin.Probe("1", i, "AB" if i % 2 else "AA")
                                    for i in range(500)}
                              | {f"x{i}": origin.Probe("X", i, "AB" if i % 6 == 0 else "AA")
                                 for i in range(500)}, [])
    assert origin.call_sex(ambiguous)[0] == "unknown"
    assert origin.call_sex(origin.Sample("empty", {}, []))[0] == "unknown"


def test_the_join_is_an_inner_join_on_probe_id(tmp_path):
    fa = origin.Sample("f", {"r1": origin.Probe("1", 10, "AA"),
                             "r2": origin.Probe("1", 20, "AA")}, [])
    em_ = origin.Sample("e", {"r2": origin.Probe("1", 20, "AB"),
                              "r3": origin.Probe("1", 30, "AB")}, [])
    got = origin.build_markers(fa, em_)
    assert [m.rsid for m in got] == ["r2"], "r1 and r3 are not shared"
    assert got[0].mother is None, "no mother supplied means None, not a guessed genotype"


def test_markers_need_not_arrive_sorted_and_off_chromosome_ones_do_not_join_a_run():
    ms = markers(60, absent_from=20, absent_to=40)
    rev = list(reversed(ms))
    a = origin.analyse_embryo("E", ms, "11", 1_060_000, calibration=GOOD, dropout=0.01)
    b = origin.analyse_embryo("E", rev, "11", 1_060_000, calibration=GOOD, dropout=0.01)
    assert a.windows[-1].longest_run == b.windows[-1].longest_run
    off = [origin.Marker(m.rsid, "7" if i % 2 else "11", m.pos, m.father, m.mother, m.embryo,
                         m.baf, m.lrr) for i, m in enumerate(ms)]
    c = origin.analyse_embryo("E", off, "11", 1_060_000, calibration=GOOD, dropout=0.01)
    assert c.windows[-1].longest_run <= a.windows[-1].longest_run


# --- the whole plate: parentage gates the homologue call --------------------------------------


def plate_markers(sample_id, *, hap, chrom="11", n=2_000, start=17_200_000, step=5_000):
    """One androgenote: homozygous for a single sperm's haplotype, which is what makes the
    scaffold work so well on them. Every fifth marker has a heterozygous father, so the plate
    carries both marker classes: homozygous ones drive parentage, heterozygous ones the scaffold.
    """
    fp, sp = {}, {}
    for i in range(n):
        k, pos = f"{chrom}:{i}", start + i * step
        if i % 5 == 0:
            allele = ("A", "B")[hap] if (i // 5) % 3 else ("B", "A")[hap]
            fp[k] = origin.Probe(chrom, pos, "AB", baf=0.5)
        else:
            allele = "A"
            fp[k] = origin.Probe(chrom, pos, "AA", baf=0.0)
        sp[k] = origin.Probe(chrom, pos, allele * 2, baf=0.0 if allele == "A" else 1.0)
    return origin.Sample("sperm", fp, []), origin.Sample(sample_id, sp, [])


def test_the_plate_runs_parentage_first_and_uses_it_as_a_gate():
    fa, a = plate_markers("andro_A", hap=0)
    _, b = plate_markers("andro_B", hap=1)
    _, gyno = sample_pair({"*": 0.068}, chroms=["11"], nonpaternal=0.12, het_band=0.03)
    gyno.sample_id = "gyno"
    rep = origin.run_experiment(fa, [a, b, gyno], variant_chrom="11", variant_pos=17_300_000,
                                anchor="andro_A", anchor_carries="mutant")
    by = {s.sample_id: s for s in rep.samples}
    assert by["gyno"].parentage.origin_class == "gynogenetic"
    assert by["gyno"].relationship is None, "a sample with no paternal genome is set aside"
    assert any("no paternal genome" in n for n in by["gyno"].notes)
    # The two androgenotes carry opposite paternal homologues, and the anchor names both.
    assert by["andro_A"].inherited == "H2_mutant_inherited"
    assert by["andro_B"].inherited == "H1_wildtype_inherited"


def test_the_runner_reads_the_gated_significance_and_never_recomputes_it():
    """The regression guard for a bug introduced in the runner itself. Re-deriving 'significant'
    from longest_run and r_min bypasses every gate above it, including the no-mother refusal, and
    reported a paternal loss on three real samples whose paternal genome is present at 0.2%."""
    fa, a = plate_markers("andro_A", hap=0)
    _, b = plate_markers("andro_B", hap=1)
    rep = origin.run_experiment(fa, [a, b], variant_chrom="11", variant_pos=17_300_000)
    for sv in rep.samples:
        if sv.variant is None:
            continue
        for w in sv.variant.windows:
            assert not w.significant, "no mother means no significance may be claimed"
        assert sv.relationship != "void", "and nothing downstream may act as if there were"
    assert "not tested" in rep.summary()


def test_an_absent_paternal_segment_voids_the_homologue_call():
    """You cannot inherit a chromosome that is not there. Reporting a homologue anyway is how a
    deletion gets read as a correction, which is the ordering this exists to enforce."""
    ms = [origin.Marker(f"m{i}", "11", 1_000_000 + i * 2_000, "AA", "BB",
                        "BB" if 20 <= i <= 60 else "AB", baf=1.0, lrr=-0.6) for i in range(120)]
    v = origin.analyse_embryo("E", ms, "11", 1_060_000, calibration=GOOD, dropout=0.01)
    assert any(w.significant for w in v.windows), "the deletion is found"
    sv = origin.SampleVerdict("E", parentage=None, variant=v, relationship="same",
                              inherited="H2_mutant_inherited")
    absent = any(w.significant for w in sv.variant.windows)
    assert absent, "and that is the condition the runner voids on"


def test_the_decision_boundary_is_this_samples_own_noise_not_a_population_constant():
    """An "unrelated rate" cannot be a decision boundary. Measured across 50+ pairs on one array
    it ran from 6.8% to 50%, because it depends on the ancestry the two people share (same-study
    pairs near 7%, cross-study 8-15%) and on how noisy the sample used as the father is."""
    # Same absence rate, different sample quality: the call must differ, because what counts as
    # surprising depends on what this sample's own noise can manufacture.
    fa, clean = sample_pair({"*": 0.03}, chroms=AUTO_AND_X, het_band=0.30)
    r_clean = origin.parental_origin(fa, clean, product=AXIOM)
    assert r_clean.verdict == "no_paternal_contribution", r_clean.explainable

    fa2, noisy = sample_pair({"*": 0.03}, chroms=AUTO_AND_X, het_band=0.30)
    for k, pr in list(noisy.probes.items())[::2]:      # half the calls lost
        noisy.probes[k] = origin.Probe(pr.chrom, pr.pos, "NC", baf=pr.baf)
    r_noisy = origin.parental_origin(fa2, noisy, product=AXIOM)
    assert r_noisy.explainable > r_clean.explainable * 5
    assert r_noisy.verdict != "no_paternal_contribution", (
        "at this noise level a true father-offspring pair can reach this rate, so absence "
        "must not be claimed")


def test_the_noise_bound_needs_both_dropout_and_heterozygosity():
    """Dropout only manufactures paternal absence by turning a HETEROZYGOUS call homozygous and
    discarding the paternal allele. A homozygous genome is immune however much drops out, which
    is why androgenotes hold at 0.16% absence despite a 13.6% no-call rate while a diploid embryo
    at 46.9% reaches 9.69% against its own father."""
    assert origin.absence_explainable(0.50, 0.0) == 0.0, "homozygous: immune to dropout"
    assert origin.absence_explainable(0.0, 0.50) == 0.0, "no dropout: nothing to inflate"
    assert origin.absence_explainable(0.136, 0.013) == pytest.approx(0.0018, abs=1e-4)
    assert origin.absence_explainable(0.469, 0.222) == pytest.approx(0.104, abs=1e-3)
    assert math.isnan(origin.absence_explainable(float("nan"), 0.2))


def test_a_true_father_offspring_pair_is_never_called_unrelated_however_noisy():
    """The adversarial case, from real data: Zuccaro A8 against its own father shows 9.69%
    absence, ABOVE genuinely unrelated pairs at 5.1-5.6%, because it has a 46.9% no-call rate.
    Any fixed threshold between those two numbers calls a real father-offspring pair unrelated."""
    fa, s = sample_pair({"*": 0.097}, chroms=AUTO_AND_X, het_band=0.222)
    for k, pr in list(s.probes.items()):
        if hash(k) % 100 < 47:
            s.probes[k] = origin.Probe(pr.chrom, pr.pos, "NC", baf=pr.baf)
    r = origin.parental_origin(fa, s, product=AXIOM)
    assert r.explainable > r.genome_rate, (r.explainable, r.genome_rate)
    assert r.verdict == "paternal_genome_present"


def test_the_second_parent_signal_is_derived_from_the_father_not_fitted():
    """A second parent supplies an allele the father lacks only where he is homozygous, at rate
    sum(p^2 q + q^2 p) = sum(pq). His own heterozygosity is sum(2pq) over the same markers, so the
    quantity is exactly half of it, measurable from his array with no reference panel.

    Checked against six pairs: 8.48% predicted vs 3.4-5.2% observed for three paternal-only
    genomes, and 8.33/16.10/9.72% predicted vs 34.66/16.35/10.65% for three biparental ones."""
    assert origin.second_parent_signal(0.170) == pytest.approx(0.085)
    assert origin.second_parent_signal(0.322) == pytest.approx(0.161)
    assert math.isnan(origin.second_parent_signal(float("nan")))


def test_a_homozygous_genome_does_not_consult_the_second_axis_at_all():
    """One allele per locus. If those alleles are the father's there is no room for a maternal
    complement, so the residual on axis 2 is error by construction. This matters because axis 2
    separates by about 1.6x where axis 1 separates by 30x, and a noisy paternal-only genome
    drifts up toward the boundary."""
    fa, s = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, nonpaternal=0.30, het_band=0.02)
    r = origin.parental_origin(fa, s, product=AXIOM)
    assert r.zygosity == "uniparental_homozygous"
    assert r.nonpaternal_rate > r.second_parent_expected, "axis 2 would say biparental"
    assert r.origin_class == "androgenetic", "but zygosity settles it without consulting axis 2"


def test_a_diploid_paternal_only_genome_is_called_and_the_thin_margin_declared():
    """Two sperm rather than a duplicated one: paternal-only but heterozygous. This is the only
    case where axis 2 has to carry the decision, so the margin is reported with it."""
    fa, s = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, nonpaternal=0.02, het_band=0.20)
    # The boundary is HALF THE FATHER'S HETEROZYGOSITY, so he has to have some.
    for i, (k, pr) in enumerate(list(fa.probes.items())):
        if i % 5 == 0 and pr.chrom != "Y":
            fa.probes[k] = origin.Probe(pr.chrom, pr.pos, "AB", baf=0.5)
    r = origin.parental_origin(fa, s, product=AXIOM)
    assert r.zygosity == "diploid"
    assert r.second_parent_expected > 0.05, r.second_parent_expected
    assert r.origin_class == "androgenetic"
    assert any("two sperm" in n for n in r.notes)


def test_zygosity_falls_back_to_genotypes_when_no_intensities_are_present():
    """HapMap-style genotype-only exports have no BAF at all. The same derived reference works:
    a biparental sample resembles the father, a uniparental one tends to zero."""
    fa, s = sample_pair({"*": 0.002}, chroms=AUTO_AND_X, nonpaternal=0.30, het_band=0.0)
    for k, pr in list(s.probes.items()):
        s.probes[k] = origin.Probe(pr.chrom, pr.pos, pr.gt, baf=None)
    for k, pr in list(fa.probes.items()):
        fa.probes[k] = origin.Probe(pr.chrom, pr.pos, "AB" if hash(k) % 5 else "AA", baf=None)
    r = origin.parental_origin(fa, s, product=AXIOM)
    assert math.isnan(r.het_band)
    assert r.zygosity in ("diploid", "uniparental_homozygous"), "a call is still made"
    assert any("No B-allele frequencies" in n for n in r.limits), "and the weaker basis is stated"


# --- segmental parent of origin ---------------------------------------------------------------


def segmental_pair(absent_span=None, *, n=6_000, chrom="4", noise=0.002, absent_rate=0.07):
    """A genome carrying the paternal contribution, optionally losing it over one span."""
    fp, sp = {}, {}
    lo, hi = absent_span or (-1, -1)
    for i in range(n):
        k, pos = f"{chrom}:{i}", 1_000 + i * 5_000
        fp[k] = origin.Probe(chrom, pos, "AA", baf=0.0)
        rate = absent_rate if lo <= pos <= hi else noise
        # Scattered, not blocked: `i % 1000` would put every absence in one contiguous run and
        # make the window rates bimodal by construction rather than by biology.
        scatter = (i * 7919) % 1000
        sp[k] = origin.Probe(chrom, pos, "BB" if scatter < rate * 1000 else "AA", baf=0.0)
    return origin.Sample("sperm", fp, []), origin.Sample("sample", sp, [])


def test_a_uniform_genome_is_not_carved_into_segments():
    """The failure this was built with. Estimating the two rates as percentiles of the observed
    window rates finds two populations even in a uniform genome, because those are just the tails
    of one noise distribution: a sample whose paternal genome is present everywhere at 0.16% came
    back as 28 'absent' segments running at 0.3-1.5%. The rates must be anchored to the
    calibrated noise bound instead."""
    fa, s = segmental_pair()
    segs, info = origin.segmental_origin(fa, s, explainable=0.007)
    assert not segs, [(x.chrom, x.rate) for x in segs]
    assert "no window exceeds" in info["reason"]


def test_a_genome_missing_the_paternal_contribution_throughout_is_not_segmented():
    fa, s = segmental_pair(noise=0.07)
    segs, info = origin.segmental_origin(fa, s, explainable=0.007)
    assert not segs
    assert "absent throughout" in info["reason"]


def test_one_lost_stretch_is_found_with_its_boundaries():
    lo, hi = 6_000_000, 14_000_000
    fa, s = segmental_pair((lo, hi))
    segs, _ = origin.segmental_origin(fa, s, explainable=0.007)
    lost = [x for x in segs if x.state == "paternal_absent"]
    assert len(lost) == 1, [(x.start, x.end, x.rate) for x in lost]
    got = lost[0]
    assert abs(got.start - lo) < 500_000 and abs(got.end - hi) < 500_000
    assert got.confident and got.margin > 3


def test_segment_margin_separates_a_real_loss_from_a_hard_region():
    """On real data the one true finding, a whole non-PAR X, sat at 18x the noise floor while
    every spurious subtelomeric stretch sat at 1-2.5x. Reporting them alike would bury it."""
    real = origin.OriginSegment("X", 2_700_151, 155_233_115, 4_880, 1_796, 0.368,
                                "paternal_absent", margin=18.0)
    hard = origin.OriginSegment("16", 85_667, 2_290_138, 1_123, 46, 0.041,
                                "paternal_absent", margin=2.0)
    assert real.confident and not hard.confident
    assert real.span_mb > 150


# --- sequencing input -------------------------------------------------------------------------


VCF_HEAD = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDAD\tKID"""


def write_vcf(tmp_path, rows, name="t.vcf"):
    p = tmp_path / name
    p.write_text(VCF_HEAD + "\n" + "\n".join(rows) + "\n")
    return p


def test_a_vcf_supplies_read_level_evidence_an_array_cannot():
    """AD makes the alternate fraction a MEASUREMENT rather than a three-way call, and DP is
    linear in copy number so log2(DP/median) is an LRR with compression exactly one. None of the
    per-array calibration that intensity needs arises."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = write_vcf(Path(d), [
            "chr1\t100\trs1\tA\tG\t50\tPASS\t.\tGT:DP:AD\t0/0:40:40,0\t0/1:40:20,20",
            "chr1\t200\trs2\tC\tT\t50\tPASS\t.\tGT:DP:AD\t1/1:40:0,40\t0/1:20:10,10",
            "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT:DP:AD\t0/1:40:20,20\t./.:0:0,0",
            "chr1\t400\trs4\tT\tC,G\t50\tPASS\t.\tGT:DP:AD\t0/1:40:20,20\t0/1:40:20,20",
        ] + [  # padding so the median depth is stable, as it is on real WGS
            f"chr1\t{1000+i}\tp{i}\tA\tG\t50\tPASS\t.\tGT:DP:AD\t0/0:40:40,0\t0/0:40:40,0"
            for i in range(20)
        ])
        dad = origin.read_sample(p, "DAD")
        kid = origin.read_sample(p, "KID")
    assert len(dad.probes) == 23, "the multi-allelic site is skipped, not guessed at"
    assert [dad.probes[k].gt for k in ("rs1", "rs2", "chr1:300")] == ["AA", "BB", "AB"]
    assert kid.probes["chr1:300"].gt == "NC"
    assert dad.probes["rs1"].baf == 0.0 and dad.probes["rs2"].baf == 1.0
    assert kid.probes["rs1"].baf == pytest.approx(0.5), "alt fraction, straight from the reads"
    # Half the depth is exactly one copy, with no compression constant anywhere.
    assert kid.probes["rs2"].lrr == pytest.approx(-1.0)
    assert kid.probes["rs1"].lrr == pytest.approx(0.0)


def test_a_variant_only_vcf_is_flagged_because_its_marker_set_is_conditioned():
    """The trap in real data: a single-sample benchmark VCF lists only sites where the sample
    differs from the reference, so its homozygous-REFERENCE positions are absent from the file.
    Those are exactly the markers a paternal-absence test needs, and the analysis would still run,
    silently, on a marker set conditioned on which sample happened to vary. Measured on GIAB
    HG001: 0 of 89,974 records homozygous reference."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = write_vcf(Path(d), [
            f"chr1\t{100+i}\trs{i}\tA\tG\t50\tPASS\t.\tGT:DP:AD\t0/1:40:20,20\t1/1:40:0,40"
            for i in range(200)
        ])
        s = origin.read_sample(p, "DAD")
    assert any("variant-only VCF" in n for n in s.notes)
    assert any("homozygous REFERENCE" in n for n in s.notes)


def test_vcf_and_array_input_reach_the_same_verdict():
    """Round-trip: the same genotypes through both readers must classify identically, so the
    sequencing path is not a second implementation with its own behaviour."""
    import tempfile
    fa, s = sample_pair({"*": 0.002}, chroms=["1"], n_per_chrom=3_000, het_band=0.02)
    rows = []
    for i, (k, fp) in enumerate(sorted(fa.probes.items())):
        if fp.chrom != "1":
            continue
        sp = s.probes[k]
        code = {"AA": "0/0", "AB": "0/1", "BB": "1/1", "NC": "./."}
        rows.append(f"chr1\t{fp.pos}\t{k}\tA\tG\t50\tPASS\t.\tGT:DP\t"
                    f"{code[fp.gt]}:40\t{code[sp.gt]}:40")
    with tempfile.TemporaryDirectory() as d:
        p = write_vcf(Path(d), rows)
        v_fa, v_s = origin.read_sample(p, "DAD"), origin.read_sample(p, "KID")
    a = origin.parental_origin(fa, s, product=AXIOM)
    b = origin.parental_origin(v_fa, v_s, product=AXIOM)
    assert b.verdict == a.verdict == "paternal_genome_present"
    assert b.genome_rate == pytest.approx(a.genome_rate, abs=1e-6)

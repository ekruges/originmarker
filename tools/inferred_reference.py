"""Reconstruct a parental genotype from haploid meiotic products, and score against it.

A haploid pronucleus is one meiotic product. Where the parent is homozygous every product
carries that allele; where the parent is heterozygous each product carries one at random, so a
heterozygous site is recognised only by observing both alleles across products. Products that
all happen to agree at a heterozygous site promote it into the reference as homozygous, and a
true offspring then reads as missing the parental allele at half of those.

Two independent quantities control that error, and an earlier design conflated them by
requiring EVERY product to be called at a marker:

  m   how many products were called there, which sets how often a heterozygous site survives
      unrecognised, at 2^(1-m)
  n   how many products exist, which sets how many markers reach any given m

Requiring m = n ties marker retention to the product with the worst call rate. At five products
of 63% to 88% call rate only a quarter of markers survive, and those are the cleanest probes,
which are enriched for sites where the parent is homozygous because the minor allele is rare.
Unrelated people carry the common allele there too, so the reference loses its ability to
exclude them: measured, an unrelated adult fell from 4.15% absence against a real array to
2.24% against such a reference, while true offspring rose. Both sides move toward each other.

Separating m from n lets contamination fall with m while retention stays governed by n.

Nothing here decides parentage on its own. It builds a reference; origin.parental_origin makes
the call, unchanged, so a reconstructed reference is judged by the same rule as a real array.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from originmarker import origin  # noqa: E402

#: Products contribute homozygous calls only. A heterozygous call in a haploid genome is an
#: error by construction, and counting it as evidence of parental heterozygosity inflates the
#: recovered figure several-fold, so it is discarded and only tallied for diagnostics.
HOM = ("AA", "BB")


@dataclass
class Products:
    """Per-marker allele observations across the haploid products, built in one pass.

    Held as observations rather than as a finished reference so that leaving a product out
    costs a filter rather than a re-read: a product scored against a reference containing
    itself shows zero absence, which is a bias equal to the entire signal.
    """

    ids: list[str] = field(default_factory=list)
    #: probe -> [(product index, allele)]
    obs: dict[str, list[tuple[int, str]]] = field(default_factory=dict)
    #: probe -> (chromosome, position)
    meta: dict[str, tuple[str, int]] = field(default_factory=dict)
    het_calls: dict[str, int] = field(default_factory=dict)
    called: dict[str, int] = field(default_factory=dict)
    #: Heterozygous BAF band per product. The m-of-n model assumes each product carries ONE
    #: allele drawn at random, so a diploid genome breaks it even though its homozygous calls
    #: are still real parental alleles: a homozygous call in a diploid means both copies agreed,
    #: which is a different sampling process and not a single meiotic draw.
    #:
    #: Gating on HET_BAND_DIPLOID separates the clear cases and nothing more. There is no empty
    #: gap to sit a threshold in: confirmed clean haploid products measure 1.0% to 7.8%, and the
    #: samples excluded here measure 10.0% and 10.7%, so the margin is a couple of points on a
    #: statistic that a low call rate inflates. Where admitting or excluding one product decides
    #: whether a group clears the depth floor, that decision needs a better ploidy call than
    #: this, not a tighter cutoff.
    band: dict[str, float] = field(default_factory=dict)

    def add(self, sample_id: str, sample: origin.Sample) -> None:
        idx = len(self.ids)
        self.ids.append(sample_id)
        het = tot = 0
        baf_in = baf_tot = 0
        for probe, p in sample.probes.items():
            if not origin._is_autosome(p.chrom):
                continue
            if p.baf is not None:
                baf_tot += 1
                baf_in += 0.35 <= p.baf <= 0.65
            if p.gt == "NC":
                continue
            tot += 1
            if p.gt == "AB":
                het += 1
                continue
            if p.gt not in HOM:
                continue
            self.obs.setdefault(probe, []).append((idx, p.gt[0]))
            self.meta[probe] = (p.chrom, p.pos)
        self.het_calls[sample_id] = het
        self.called[sample_id] = tot
        self.band[sample_id] = baf_in / baf_tot if baf_tot else float("nan")


@dataclass
class Reference:
    sample: origin.Sample
    n_products: int
    m_min: int
    markers: int
    #: The PARENT's heterozygous fraction over the markers that passed the m_min filter,
    #: recovered from disagreement between products. It measures ascertainment: how far the
    #: marker set the reference draws from has drifted from the genome, which is what decides
    #: whether the reference can still exclude an unrelated genome.
    #:
    #: This is NOT the reference's own heterozygosity. The reference is homozygous everywhere by
    #: construction, and the fraction of it that is secretly heterozygous is `contamination`.
    #: The two differ by an order of magnitude (16.4% against 1.6% at n=6, m>=4) and conflating
    #: them misreads a headline figure, so both are reported.
    h_retained: float
    #: Expected fraction of the reference that is a heterozygous site mistaken for homozygous.
    contamination: float
    #: Absence a HAPLOID true offspring reads at contaminated markers: it carries the other
    #: allele at exactly half of them. A diploid offspring must also receive the minor allele
    #: from its other parent, so it lands well short of half. Measured on a real diploid
    #: offspring against a contaminated reference: 25.7% of contaminated markers, not 50%, so
    #: this figure overstates the damage to a diploid by roughly 1.9x and understates nothing.
    #: It is the haploid case, which is the case this builder is for.
    spurious_absence: float
    mean_m: float
    excluded: Optional[str]



#: Measured excess of P(all m products agree | parent heterozygous) over the independent-draw
#: value 2^(1-m). Products are not independent at a heterozygous site: a probe with
#: allele-specific dropout calls the same homozygote in every product, and those are exactly
#: the low-m probes. Without this, reported contamination runs 20%-29% below truth at m<=5.
#:
#: Measured against a real parental array over 133,631 heterozygous autosomal markers, using
#: six haploid products of that parent:
#:
#:     m      markers   P(all agree | het)   2^(1-m)   ratio
#:     2        4,332              0.61958   0.50000   1.24
#:     3       10,588              0.35304   0.25000   1.41
#:     4       23,830              0.17709   0.12500   1.42
#:     5       42,153              0.08023   0.06250   1.28
#:     6       49,196              0.03132   0.03125   1.00
#:
#: The excess vanishes by m=6, so a reference deep enough to be usable needs no correction.
AGREEMENT_EXCESS = {2: 1.24, 3: 1.41, 4: 1.42, 5: 1.28}


def _p_agree(m: int) -> float:
    """P(all m products show the same allele | the parent is heterozygous)."""
    return min(1.0, 2.0 ** (1 - m) * AGREEMENT_EXCESS.get(m, 1.0))


def _h_from_disagreement(obs: Iterable[list[tuple[int, str]]], m_min: int) -> tuple[float, float]:
    """Recover the parent's heterozygous fraction from how often the products disagree.

    With m products called at a marker, a heterozygous parent is detected unless all m agree,
    which happens at 2^(1-m). So the observed disagreement rate estimates h times the mean
    detection probability, and dividing by that mean inverts it. m varies per marker, so the
    mean is taken over markers rather than assumed constant.
    """
    disagreed = 0
    total = 0
    detect = 0.0
    for v in obs:
        m = len(v)
        if m < m_min:
            continue
        total += 1
        detect += 1.0 - _p_agree(m)
        if len({a for _, a in v}) > 1:
            disagreed += 1
    if not total or detect <= 0:
        return math.nan, math.nan
    mean_detect = detect / total
    return (disagreed / total) / mean_detect, mean_detect


#: Largest tolerated drift of the retained marker set away from the genome, as a fraction of the
#: parent's heterozygosity. Raising m lowers contamination and narrows the marker set toward
#: probes nearly every product called, which are enriched for sites where the parent is
#: homozygous because the minor allele is rare. Unrelated people carry the common allele there
#: too, so past a point the reference loses its grip on them.
#:
#: Measured on a father whose real array gives 16.66% heterozygosity, five haploid products:
#:
#:     m>=   markers   h(retained)   of true   contamination   products back   closest negative
#:      2    663,542       16.93%      102%         4.288%          0 of 5           5.96x
#:      3    596,687       16.48%       99%         3.402%          2 of 5           5.65x
#:      4    457,506       15.47%       93%         2.286%          3 of 5           5.25x
#:      5    231,064       13.74%       82%         1.258%          3 of 5           4.67x
#:
#: Recovery plateaus at m=4 while everything else keeps degrading, so 0.90 sits at the knee.
MIN_ASCERTAINMENT = 0.90

#: Below this many products the method inverts: at three and fewer, every true offspring
#: tested came back a decisive wrong answer rather than a refusal, 24 of 24 across two
#: experiments. Mirrors `MIN_PRODUCTS` in web/src/inferredReference.ts, which is where the
#: refusal is enforced; here it is the number the validation scripts report against.
MIN_PRODUCTS = 5


#: Opposite-homozygote thresholds. Measured across two experiments: 4.68% to 9.70% between
#: products of one parent over 46 pairs, and 9.88% to 16.10% between products of different
#: parents over 45. The separation is real but the margin is 0.18 of a percentage point, and one
#: genuine cross-parent pair sits at 9.88%, under the cut, so a per-pair label is not evidence of
#: group membership. Grouping requires every pair, which turns that one misread into no error.
SAME_PARENT_MAX = 0.105
DIFFERENT_PARENT_MIN = 0.125


def kinship(rate: float) -> str:
    if rate < SAME_PARENT_MAX:
        return "same parent"
    return "different parents" if rate >= DIFFERENT_PARENT_MIN else "ambiguous"


def group_by_parent(n: int, rate) -> list[list[int]]:
    """Split products into groups sharing one parent. EVERY pair inside a group must agree.

    The largest such set is found exactly, because greedy clique-building is not sufficient and
    the failure is not hypothetical: with two parents of three products each and one genuine
    cross-parent pair reading under the cut, a greedy pass seeded on either member of that pair
    returns the pair itself and splits both parents in half. Bron-Kerbosch with pivoting, over at
    most a few dozen products, which is what an experiment holds.

    Products must already be quality-gated and haploid. A diploid compared with a haploid has a
    different expected rate entirely and bridges unrelated groups, and an array below the
    call-rate floor reads high against everyone, which splits one parent into several.

    Mirrored in web/src/inferredReference.ts and pinned against it by tools/parity_check.py.
    """
    adj = [{b for b in range(n) if b != a and kinship(rate(a, b)) == "same parent"}
           for a in range(n)]
    best: list[int] = []

    def expand(r: list[int], pset: set[int], x: set[int]) -> None:
        nonlocal best
        if not pset and not x:
            if len(r) > len(best):
                best = list(r)
            return
        # Pivot on the candidate with the most neighbours, which prunes branches that cannot win.
        pivot, deg = -1, -1
        for u in pset | x:
            d = len(pset & adj[u])
            if d > deg:
                pivot, deg = u, d
        for v in list(pset):
            if pivot >= 0 and v in adj[pivot]:
                continue
            expand(r + [v], pset & adj[v], x & adj[v])
            pset.discard(v)
            x.add(v)

    remaining = set(range(n))
    groups: list[list[int]] = []
    while remaining:
        best = []
        expand([], set(remaining), set())
        # A product sharing no parent with anything else is its own group of one.
        clique = best if best else [min(remaining)]
        remaining -= set(clique)
        groups.append(sorted(clique))
    return sorted(groups, key=lambda g: (-len(g), g[0]))


def choose_m(p: Products,
             exclude: Optional[str | Iterable[str]] = None) -> tuple[int, dict[int, float]]:
    """The largest m whose marker set still resembles the genome.

    m = n - 1 was the earlier rule and it happens to be right at n=5. It stops being right as n
    grows, because what governs the drift is the absolute number of products required to agree,
    not the slack below n: at n=8 the same rule demands seven and the retained set falls to
    roughly 62% of the genome's heterozygosity.

    Needs no real parental array. h_retained at m=2 is measured unbiased against a known father
    (16.93% against 16.66%), so it serves as the baseline the stricter settings are judged
    against, and the whole rule is computable from the products alone.
    """
    names = ([exclude] if isinstance(exclude, str) else list(exclude or []))
    n = len(p.ids) - len(names)
    base = build(p, 2, exclude=exclude).h_retained
    ratios: dict[int, float] = {}
    best = 2
    for m in range(2, n + 1):
        h = build(p, m, exclude=exclude).h_retained
        ratios[m] = h / base if base else math.nan
        if math.isfinite(ratios[m]) and ratios[m] >= MIN_ASCERTAINMENT:
            best = m
    return best, ratios


def build(p: Products, m_min: int, exclude: Optional[str | Iterable[str]] = None,
          h_prior: Optional[float] = None) -> Reference:
    """A reference from every marker where at least m_min products spoke and all agreed.

    `exclude` takes one product id or several, so subsets of a product set that was read once
    cost a filter rather than a re-read.
    """
    names = ([exclude] if isinstance(exclude, str) else list(exclude or []))
    dropped = {p.ids.index(x) for x in names}
    kept: dict[str, origin.Probe] = {}
    surviving: list[list[tuple[int, str]]] = []
    m_sum = 0
    for probe, v in p.obs.items():
        w = [x for x in v if x[0] not in dropped] if dropped else v
        if len(w) < m_min:
            continue
        surviving.append(w)
        if len({a for _, a in w}) > 1:
            continue
        chrom, pos = p.meta[probe]
        kept[probe] = origin.Probe(chrom, pos, w[0][1] * 2)
        m_sum += len(w)

    h_ret, _ = _h_from_disagreement(surviving, m_min)
    h = h_prior if h_prior is not None and math.isfinite(h_prior) else h_ret

    # Contamination is computed per marker from its own m, then averaged, because a marker
    # seen by ten products is far safer than one seen by the minimum and averaging the
    # exponent instead of the probability would understate the tail.
    # Not zero: an empty or unmeasurable reference has UNKNOWN contamination, and reporting
    # 0.0 there makes a failed build read as a perfect one in any aggregation of the column.
    contam = math.nan
    if kept and math.isfinite(h):
        acc = 0.0
        for probe in kept:
            v = p.obs[probe]
            m = len([x for x in v if x[0] not in dropped]) if dropped else len(v)
            odds = h * _p_agree(m)
            acc += odds / (1.0 - h + odds)
        contam = acc / len(kept)

    n = len(p.ids) - len(dropped)
    return Reference(
        sample=origin.Sample(f"inferred-n{n}-m{m_min}", kept, []),
        n_products=n, m_min=m_min, markers=len(kept), h_retained=h_ret,
        contamination=contam, spurious_absence=contam / 2.0,
        mean_m=m_sum / len(kept) if kept else math.nan,
        excluded=";".join(names) if names else None,
    )


if __name__ == "__main__":
    # The browser half is checked by web/src/inferredReference.check.ts and the two are pinned to
    # each other by tools/parity_check.py, which needs real arrays. This is the part that runs
    # with no data at all: the arithmetic, and the two rules that decide what gets built.
    import math

    # The agreement probability is the MEASURED one below m=6, not the independent draw.
    assert abs(_p_agree(2) - 0.62) < 1e-12
    assert abs(_p_agree(4) - 0.1775) < 1e-12
    assert _p_agree(6) == 2 ** -5, "the excess vanishes by six, so no correction applies"
    assert _p_agree(10) == 2 ** -9

    # Thresholds are the measured ones. A same-parent pair at 9.88% is real and sits under the
    # cut, which is why membership is all-pairs rather than per-pair.
    assert (SAME_PARENT_MAX, DIFFERENT_PARENT_MIN) == (0.105, 0.125)
    assert kinship(0.0468) == "same parent"
    assert kinship(0.0988) == "same parent"
    assert kinship(0.161) == "different parents"
    assert kinship(0.115) == "ambiguous", "the gap is 0.18 points wide, so this exists"

    # Grouping must not chain two parents through one bridging pair under the cut. Measured:
    # one genuine cross-father pair reads 9.88% while reading 12.1% to 12.4% against the rest.
    A, B = {0, 1, 2}, {3, 4, 5}

    def _rate(a: int, b: int) -> float:
        if a == b:
            return 0.0
        if (a in A and b in A) or (a in B and b in B):
            return 0.07
        return 0.0988 if {a, b} == {2, 3} else 0.13

    assert group_by_parent(6, _rate) == [[0, 1, 2], [3, 4, 5]]
    assert group_by_parent(5, lambda a, b: 0.0 if a == b else 0.07) == [[0, 1, 2, 3, 4]]

    # Contamination is a POSTERIOR, so it carries the prior that most markers are homozygous:
    # at 17% heterozygosity with four products agreeing, only 3.5% of the reference is wrong,
    # not the 17.75% that P(agree | het) alone would suggest.
    def _contam(h: float, m: int) -> float:
        odds = h * _p_agree(m)
        return odds / (1 - h + odds)

    assert abs(_contam(0.17, 4) - 0.0351) < 5e-4
    assert _contam(0.0, 4) == 0.0, "a homozygous parent contaminates nothing"
    assert math.isclose(_contam(1.0, 4), 1.0), "an all-het parent is all contamination"
    for m in range(2, 9):
        assert _contam(0.17, m + 1) < _contam(0.17, m), "deeper must be cleaner"
    assert MIN_ASCERTAINMENT == 0.90 and MIN_PRODUCTS == 5

    print("tools/inferred_reference.py self-check OK")

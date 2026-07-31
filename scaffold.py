"""scaffold - which of the father's two homologues each embryo received.

This is the half of the question that genotypes alone cannot answer, and the reason is worth
stating plainly. "Was the disease allele inherited and corrected, or was the wild-type
inherited" is a question about WHICH of the father's chromosomes came through. Both of his
homologues are simply "his alleles" in any genotype file; nothing in the data says which one
carries the mutation.

What breaks the deadlock is having several embryos from the same father. Each embryo received
exactly one paternal homologue, so embryos that received the same one agree at every marker
where the father is heterozygous, and embryos that received different ones disagree at every
such marker. That splits the embryos into two groups without needing to know which group is
which. One embryo whose variant-site genotype was established outside the array then labels
both groups at once, and the label propagates to every sibling.

So the pieces are: co-inheritance to build the partition, and a single external anchor to name
it. No long reads, no reference panel, and no new sequencing.

Two things limit it, and both are reported rather than assumed away:

  crossover   a recombination between the markers and the variant moves the answer to the other
              homologue. Distance is measured in centimorgans from the bundled deCODE map, not
              in base pairs, because base pairs do not recombine at a constant rate.
  deduction   working out which paternal allele an embryo received can be wrong, not merely
              absent, when an allele drops out. A heterozygous embryo reading as homozygous
              yields a confident deduction of the wrong allele, at roughly half the dropout
              rate, which is why bulk DNA matters so much more here than sample count.

Self-check:  python -m scaffold
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal, Optional, Sequence

import genetic_map

#: Beyond this the recombination fraction is close enough to one half that a marker carries
#: almost nothing about the variant, and including it only adds error.
DEFAULT_WINDOW_CM = 10.0

#: Residual per-marker deduction error when nothing else is known: Kothiyal's Mendelian floor.
DEDUCTION_FLOOR = 0.0063


def paternal_allele(father: str, mother: Optional[str], embryo: str) -> Optional[str]:
    """Which allele the embryo received FROM THE FATHER, or None if more than one fits.

    Needs a heterozygous father: a homozygous one transmits a known allele, but his homologues
    are indistinguishable through it, so it says nothing about which chromosome came through.
    That is the exact complement of the presence test in `origin`, which needs him homozygous.
    The two halves of this tool read opposite marker classes, which is why a panel informative
    for one can be useless for the other.
    """
    if father != "AB" or embryo in ("NC", "", None):
        return None
    maternal = ["A", "B"] if mother in (None, "NC", "") else sorted(set(mother))
    fits = {p for p in ("A", "B") for m in maternal
            if sorted(p + m) == sorted(embryo)}
    return fits.pop() if len(fits) == 1 else None


@dataclass
class Marker:
    rsid: str
    chrom: str
    pos: int
    father: str
    mother: Optional[str]
    #: embryo id -> genotype
    embryos: dict[str, str] = field(default_factory=dict)


Relationship = Literal["same", "opposite", "undetermined"]


def flip_probability(theta: float) -> float:
    """Chance the RELATIONSHIP between two embryos flips across an interval.

    It flips when exactly one of the two meioses recombined there, so 2*theta*(1-theta). At
    theta = 0.5 this is 0.5: the two embryos become independent and nothing can be carried
    across.
    """
    return 2.0 * theta * (1.0 - theta)


def relationship_posterior(
    agreements: Sequence[Optional[bool]], flips: Sequence[float], err: float
) -> list[float]:
    """P(same homologue) at each position, by forward-backward over a two-state chain.

    Summing per-marker log-odds instead would be wrong, and wrong in the dangerous direction.
    Markers on the same side of a crossover all report the same relationship, so they are one
    observation repeated, not many independent ones. Adding them multiplies a single fact by the
    marker count and manufactures certainty: 559 markers within 1 cM produced odds around 10^50
    for a claim whose real error rate is set by whether a crossover fell in between, which no
    amount of marker density reduces.

    Modelling the relationship as a chain with map-derived transitions fixes both halves. The
    state can switch where a real crossover happened, and confidence at the variant stays bounded
    by the recombination distance to the markers that support it.

    `agreements[i]` is None at a position carrying no observation, which is how the variant
    itself is included and read off.
    """
    n = len(agreements)
    if n == 0:
        return []
    # Two independent meioses, each mis-deduced with probability err: they agree when both are
    # right or both are wrong.
    p_agree_same = (1.0 - err) ** 2 + err ** 2
    p_agree = (p_agree_same, 1.0 - p_agree_same)

    def emit(s: int, i: int) -> float:
        a = agreements[i]
        if a is None:
            return 1.0
        return p_agree[s] if a else 1.0 - p_agree[s]

    fwd = [[0.0, 0.0] for _ in range(n)]
    fwd[0] = [0.5 * emit(0, 0), 0.5 * emit(1, 0)]
    for i in range(1, n):
        f = flips[i - 1]
        for s in (0, 1):
            fwd[i][s] = (fwd[i - 1][s] * (1.0 - f) + fwd[i - 1][1 - s] * f) * emit(s, i)
        tot = fwd[i][0] + fwd[i][1]
        if tot > 0.0:
            fwd[i] = [fwd[i][0] / tot, fwd[i][1] / tot]

    bwd = [[1.0, 1.0] for _ in range(n)]
    for i in range(n - 2, -1, -1):
        f = flips[i]
        for s in (0, 1):
            bwd[i][s] = ((1.0 - f) * emit(s, i + 1) * bwd[i + 1][s]
                         + f * emit(1 - s, i + 1) * bwd[i + 1][1 - s])
        tot = bwd[i][0] + bwd[i][1]
        if tot > 0.0:
            bwd[i] = [bwd[i][0] / tot, bwd[i][1] / tot]

    out = []
    for i in range(n):
        a, b = fwd[i][0] * bwd[i][0], fwd[i][1] * bwd[i][1]
        out.append(a / (a + b) if (a + b) > 0.0 else 0.5)
    return out


@dataclass
class PairVerdict:
    other: str
    relationship: Relationship
    log10_odds: float
    n_markers: int
    agreements: int
    nearest_cm: float


@dataclass
class EmbryoCall:
    embryo_id: str
    #: Relative to the anchor embryo.
    relationship: Relationship
    log10_odds: float
    n_markers: int
    nearest_cm: float
    #: Only once an anchor names the groups.
    homologue: Optional[Literal["mutant", "wildtype"]] = None
    inherited: Optional[Literal["H2_mutant_inherited", "H1_wildtype_inherited"]] = None
    notes: list[str] = field(default_factory=list)


@dataclass
class Scaffold:
    chrom: str
    variant_pos: int
    embryo_ids: list[str]
    informative_markers: int
    window_cm: float
    deduction_error: float
    anchor: Optional[str] = None
    anchor_carries: Optional[str] = None
    calls: list[EmbryoCall] = field(default_factory=list)
    pairs: list[PairVerdict] = field(default_factory=list)
    refusals: list[str] = field(default_factory=list)
    map_approx: bool = False

    def summary(self) -> str:
        lines = [f"paternal homologue scaffold, chr{self.chrom}:{self.variant_pos:,}",
                 f"  {len(self.embryo_ids)} embryos, {self.informative_markers} markers with a "
                 f"heterozygous father inside {self.window_cm:g} cM",
                 f"  deduction error rate {self.deduction_error:.4f}"]
        if self.map_approx:
            lines.append("  genetic distances are APPROXIMATE for this region (see map source)")
        if self.anchor:
            lines.append(f"  anchored on {self.anchor}, externally typed as "
                         f"{self.anchor_carries}")
        for c in self.calls:
            got = c.inherited or f"group {c.relationship} (unlabelled)"
            lines.append(
                f"  {c.embryo_id:<20s} {got:<28s} "
                f"10^{c.log10_odds:+.1f} on {c.n_markers} markers, "
                f"nearest informative {c.nearest_cm:.3f} cM"
            )
            for n in c.notes:
                lines.append(f"      {n}")
        for r in self.refusals:
            lines.append(f"  REFUSED: {r}")
        return "\n".join(lines)


def build_scaffold(
    markers: Sequence[Marker],
    chrom: str,
    variant_pos: int,
    *,
    anchor: Optional[str] = None,
    anchor_carries: Optional[Literal["mutant", "wildtype"]] = None,
    dropout: Optional[float] = None,
    window_cm: float = DEFAULT_WINDOW_CM,
    decisive_log10: float = 2.0,
    build: Optional[str] = None,
) -> Scaffold:
    """Group embryos by which paternal homologue they carry, and label the groups if anchored."""
    ids: list[str] = []
    for m in markers:
        for e in m.embryos:
            if e not in ids:
                ids.append(e)

    # An allele that drops out does not merely go missing here. A heterozygous embryo reading as
    # homozygous produces a CONFIDENT deduction of the wrong paternal allele, which happens at
    # about half the dropout rate. Silence would be harmless; a wrong answer is not.
    err = max(DEDUCTION_FLOOR, (dropout or 0.0) / 2.0)

    want = chrom.replace("chr", "").lower()
    gm = genetic_map.load(want)
    approx = gm.approx or not gm.covers(variant_pos)

    # Coordinates from another assembly index the map at the wrong locations and return a
    # plausible number. Where the build is known and does not match, the map is not consulted: a
    # uniform rate is a stated approximation, a misread map is a silent error.
    mismatch = build is not None and build != genetic_map.MAP_BUILD

    def cm_between(a: int, b: int) -> float:
        if mismatch:
            return abs(b - a) / 1e6 * genetic_map.APPROX_CM_PER_MB
        return gm.cm_between(a, b)

    deduced: list[tuple[int, float, dict[str, str]]] = []
    for m in markers:
        if m.chrom.replace("chr", "").lower() != want or m.father != "AB":
            continue
        cm = cm_between(variant_pos, m.pos)
        if cm > window_cm:
            continue
        row = {e: a for e in ids
               if (a := paternal_allele(m.father, m.mother, m.embryos.get(e, "NC"))) is not None}
        if len(row) >= 2:
            deduced.append((m.pos, cm, row))
    # Ordered along the chromosome, because the chain's transitions are physical.
    deduced.sort(key=lambda t: t[0])

    sc = Scaffold(chrom=want, variant_pos=variant_pos, embryo_ids=ids,
                  informative_markers=len(deduced), window_cm=window_cm,
                  deduction_error=err, anchor=anchor, anchor_carries=anchor_carries,
                  map_approx=approx or mismatch)
    if mismatch:
        sc.refusals.append(
            f"Coordinates are {build} and the bundled genetic map is {genetic_map.MAP_BUILD}. "
            "The map was NOT consulted: reading it at another assembly's positions returns a "
            "plausible number from the wrong locus, which is worse than an approximation that "
            "says so. Distances here are a uniform 1 cM/Mb, which ignores hotspots, so the odds "
            "below are indicative rather than calibrated. Supplying build-matched positions "
            "restores the real map. There is no liftOver here by design: UCSC chain files carry "
            "a non-commercial field-of-use restriction Apache 2.0 cannot sublicense."
        )

    if len(ids) < 2:
        sc.refusals.append(
            "Fewer than two embryos. The partition is built from co-inheritance BETWEEN embryos, "
            "so a single embryo has nothing to be compared with and its homologue stays unnamed."
        )
        return sc
    if not deduced:
        sc.refusals.append(
            f"No markers within {window_cm:g} cM where the father is heterozygous and at least "
            "two embryos yield a deducible paternal allele. A homozygous father transmits a known "
            "allele but tells you nothing about WHICH of his chromosomes carried it."
        )
        return sc

    reference = anchor if anchor in ids else ids[0]
    for other in ids:
        if other == reference:
            continue
        # Positions along the chromosome, with the variant inserted as an unobserved node so its
        # posterior can be read directly rather than borrowed from the nearest marker.
        nodes: list[tuple[int, Optional[bool]]] = [
            (pos, row[reference] == row[other])
            for pos, _cm, row in deduced if reference in row and other in row
        ]
        n = len(nodes)
        agree_n = sum(1 for _, a in nodes if a)
        nearest = min((abs(p - variant_pos) for p, _ in nodes), default=None)
        nearest_cm = cm_between(variant_pos, variant_pos + nearest) if nearest else 0.0
        nodes.append((variant_pos, None))
        nodes.sort(key=lambda t: t[0])
        at = next(i for i, (p, a) in enumerate(nodes) if p == variant_pos and a is None)

        flips = [flip_probability(genetic_map.haldane_theta(
            cm_between(nodes[i][0], nodes[i + 1][0]))) for i in range(len(nodes) - 1)]
        post = relationship_posterior([a for _, a in nodes], flips, err)
        p_same = min(max(post[at], 1e-300), 1.0 - 1e-16)
        log10 = math.log10(p_same / (1.0 - p_same))
        rel: Relationship = ("same" if log10 >= decisive_log10
                             else "opposite" if log10 <= -decisive_log10
                             else "undetermined")
        sc.pairs.append(PairVerdict(other, rel, log10, n, agree_n, nearest_cm))

    by_id = {p.other: p for p in sc.pairs}
    for e in ids:
        if e == reference:
            sc.calls.append(EmbryoCall(e, "same", math.inf, len(deduced), 0.0,
                                       notes=["reference embryo for the partition"]))
            continue
        p = by_id[e]
        call = EmbryoCall(e, p.relationship, p.log10_odds, p.n_markers, p.nearest_cm)
        if p.relationship == "undetermined":
            call.notes.append(
                f"co-inheritance with {reference} is not resolved: {p.agreements}/{p.n_markers} "
                "markers agree, which neither hypothesis explains well. Most often a crossover "
                "inside the window, or deduction errors from dropout."
            )
        sc.calls.append(call)

    if anchor is None or anchor_carries is None:
        sc.refusals.append(
            "No anchor. Co-inheritance splits the embryos into two groups but cannot name them: "
            "nothing in the array says which of the father's homologues carries the mutation. "
            "Supply one embryo whose variant-site genotype was established outside the array "
            "and both groups are named at once. Until then H1 and H2 stay unseparated."
        )
        return sc
    if anchor not in ids:
        sc.refusals.append(f"Anchor {anchor!r} is not among the embryos: {', '.join(ids)}.")
        return sc

    other_label = "wildtype" if anchor_carries == "mutant" else "mutant"
    for c in sc.calls:
        if c.relationship == "undetermined":
            continue
        c.homologue = anchor_carries if c.relationship == "same" else other_label
        c.inherited = ("H2_mutant_inherited" if c.homologue == "mutant"
                       else "H1_wildtype_inherited")
    return sc


if __name__ == "__main__":
    # Two paternal homologues over 40 markers around a variant, three embryos: A and C receive
    # homologue 1, B receives homologue 2. The mother is homozygous throughout so every paternal
    # allele is deducible.
    hap1 = ["A" if i % 3 else "B" for i in range(40)]
    hap2 = ["B" if a == "A" else "A" for a in hap1]
    ms = []
    for i in range(40):
        mat = "AA"
        ms.append(Marker(f"m{i}", "11", 17_200_000 + i * 5_000, "AB", mat, {
            "A": "".join(sorted(hap1[i] + "A")),
            "B": "".join(sorted(hap2[i] + "A")),
            "C": "".join(sorted(hap1[i] + "A")),
        }))

    s = build_scaffold(ms, "11", 17_300_000, dropout=0.01)
    assert s.informative_markers == 40, s.informative_markers
    rel = {c.embryo_id: c.relationship for c in s.calls}
    assert rel == {"A": "same", "B": "opposite", "C": "same"}, rel
    assert any("No anchor" in r for r in s.refusals), s.refusals
    assert all(c.inherited is None for c in s.calls), "unanchored groups must stay unnamed"

    anchored = build_scaffold(ms, "11", 17_300_000, anchor="A", anchor_carries="mutant",
                              dropout=0.01)
    got = {c.embryo_id: c.inherited for c in anchored.calls}
    assert got == {"A": "H2_mutant_inherited",
                   "B": "H1_wildtype_inherited",
                   "C": "H2_mutant_inherited"}, got
    assert not any("No anchor" in r for r in anchored.refusals)

    # A homozygous father is silent here, however many embryos there are.
    flat = [Marker(m.rsid, m.chrom, m.pos, "AA", m.mother, m.embryos) for m in ms]
    assert build_scaffold(flat, "11", 17_300_000).informative_markers == 0

    assert paternal_allele("AB", "AA", "AB") == "B"
    assert paternal_allele("AB", "AA", "AA") == "A"
    assert paternal_allele("AB", "AB", "AB") is None, "both parents het: cannot say"
    assert paternal_allele("AB", "AB", "AA") == "A", "a homozygous embryo still resolves"
    assert paternal_allele("AA", "AB", "AB") is None, "homozygous father says nothing"
    assert paternal_allele("AB", None, "AA") == "A", "works without the mother, less often"
    assert paternal_allele("AB", None, "AB") is None

    print("scaffold self-check OK")

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
)

VERSION = "1.3.2"
CODENAME = "Diakinesis"

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

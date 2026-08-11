"""
OriginMarker's library modules: everything the API, the tools and the tests import.

Flat on purpose. These twelve are one dependency chain rather than two subsystems that happen to
share a repo - `origin` reaches `hmm` reaches `emissions`, `panelbuilder` reaches `genetic_map`
reaches the bundled maps - so splitting them into sub-packages would put a directory boundary
through the middle of a call path and buy nothing.

Each module carries a self-check named in its own docstring. They run as
`python -m originmarker.<module>`, not as a path: the package has to be importable for a module
to find its neighbours.

`app/` is the HTTP layer over this. `tools/` and `audit/` are scripts that import it. `data/`
stays at the repo root rather than in here: it is 23MB of bundled recombination maps that the
image copies as its own layer, and `genetic_map` resolves it relative to the repo, not to itself.
"""

<img src="docs/originmarker.svg" alt="OriginMarker" width="340">

Builds a ranked, downloadable menu of candidate flanking SNP markers around a pathogenic variant,
for PGT-M linkage and karyomapping. Syngamy then reads SNP arrays after an experiment and reports
which parental genome is present, on which chromosomes, and across which regions within them.
Progenitor builds a SNP array for a parent nobody genotyped, out of the haploid cells that parent
produced. It works out which of the two parents it built from the Y those cells carry, writes the
array out as a file usable anywhere an array goes, and calls the parental origin of every input
against it.

Developed by & for the [Egli Lab](https://eglilab.com) at Columbia University Irving Medical
Center.

> **Research use only. Candidate markers require validation and per-family phasing in a
> qualified genetics laboratory. Not a clinical diagnostic.**

Live at **[originmarker.app](https://originmarker.app/)**, and in parallel at
**[ezrakruger.cc/originmarker](https://ezrakruger.cc/originmarker/)**.

Genotype files are read in the browser. There is no endpoint that receives one.

## Running it

```sh
docker compose up -d --build
```

Then <http://localhost:8091>. No API keys are required; optional ones are listed in
`deploy/README-deploy.md` and the app degrades gracefully without them.

Without Docker:

```sh
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --reload
cd web && npm install && npm run dev
```

The Syngamy analysis also runs standalone, with more than the page exposes:

```sh
python -m originmarker.origin --father sperm.txt --samples embryo1.txt embryo2.txt
```

## Tests

```sh
export PANELBUILDER_CACHE=tests/fixtures PANELBUILDER_CACHE_TTL=0
.venv/bin/python -m pytest tests/
cd web && npm run build && for f in src/*.check.ts; do node "$f"; done
```

Offline, against recorded API responses in `tests/fixtures/`. `npm run build` is `tsc -b && vite
build`, and the `-b` is the typecheck: bare `tsc` has no inputs here and exits 0 regardless.

## Layout

| | |
|---|---|
| `originmarker/` | the library: panel building, array reading, the HMM, the report writer. Every module carries a self-check, run as `python -m originmarker.<module>` |
| `app/` | FastAPI over that library. `app/main.py` is the only entry point |
| `web/` | the SPA. Syngamy and Progenitor run entirely in the browser; nothing is uploaded |
| `data/` | 23MB of bundled deCODE recombination maps, read from disk and never re-downloaded |
| `tests/` | offline, against recorded API responses in `tests/fixtures/` |
| `tools/` | one-off analysis scripts, kept because their numbers are cited in the docs |
| `audit/` | the accuracy audit and its record |
| `docs/` | validation write-ups per major version |
| `deploy/` | the deployment runbook |

## Documentation

Method, scope, data sources and limitations: [Documentation](https://originmarker.app/#/docs),
[Syngamy](https://originmarker.app/#/syngamy-docs) and
[Progenitor](https://originmarker.app/#/progenitor-docs).

Version history and the bugs each release fixed: [CHANGELOG.md](CHANGELOG.md).

Accuracy audit, 27 public arrays with bench-established answers: [audit/](audit/).

## Command line

`om` exposes the same modules the browser runs, with every constant as a flag. The web tool infers
them all and offers no knobs on purpose; this does the opposite.

```bash
node --experimental-strip-types cli/om.ts constants
```

| | |
|---|---|
| `om stage <array>` | material, dropout, marker floor |
| `om link <parent> <sample>...` | child, duplicate, haploid product, or unrelated |
| `om origin <parent> <sample>` | which parent's copy is missing, per region |
| `om cohort <dir> --ref <array>` | the same for every confirmed child under a directory |
| `om census <dir>` | haploid products per donor group |
| `om reconstruct <product>...` | a parent's genotypes from that parent's haploid cells |
| `om enrich <regions.tsv>` | positional enrichment against a marker-matched null |
| `om constants` | every tunable, its value, and why it is that value |

Omitting every flag is exactly the configuration the web tool runs and the audits measured. Moving
one is printed in the output, since a number produced under a changed constant is not comparable
with the validation figures.

## License

[Apache 2.0](LICENSE).

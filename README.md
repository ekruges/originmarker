<img src="docs/originmarker.svg" alt="OriginMarker" width="340">

Builds a ranked, downloadable menu of candidate flanking SNP markers around a pathogenic variant,
for PGT-M linkage and karyomapping. Syngamy then reads SNP arrays after an experiment and reports
which parental genome is present, on which chromosomes, and across which regions within them.
Progenitor reconstructs a parent's genotype from the haploid cells that parent produced, for when
no array of the parent exists.

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
python origin.py --father sperm.txt --samples embryo1.txt embryo2.txt
```

## Tests

```sh
export PANELBUILDER_CACHE=tests/fixtures PANELBUILDER_CACHE_TTL=0
.venv/bin/python -m pytest tests/
cd web && npm run build && for f in src/*.check.ts; do node "$f"; done
```

Offline, against recorded API responses in `tests/fixtures/`. `npm run build` is `tsc -b && vite
build`, and the `-b` is the typecheck: bare `tsc` has no inputs here and exits 0 regardless.

## Documentation

Method, scope, data sources and limitations: [Documentation](https://originmarker.app/#/docs),
[Syngamy](https://originmarker.app/#/syngamy-docs) and
[Progenitor](https://originmarker.app/#/progenitor-docs).

Version history and the bugs each release fixed: [CHANGELOG.md](CHANGELOG.md).

Accuracy audit, 27 public arrays with bench-established answers: [audit/](audit/).

## License

[Apache 2.0](LICENSE).

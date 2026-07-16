<img src="docs/originmarker.svg" alt="OriginMarker" width="340">

Builds a ranked, downloadable menu of candidate flanking SNP markers around a pathogenic
variant, for determining which parental allele an embryo inherited (PGT-M linkage /
karyomapping).

> **Research use only. Candidate markers require validation and per-family phasing in a
> qualified genetics laboratory. Not a clinical diagnostic.**

Live at **[ezrakruger.cc/originmarker](https://ezrakruger.cc/originmarker/)**.

## What it does

Give it a gene and a pathogenic variant (HGVS or rsID). It resolves the variant against
ClinVar and Ensembl, pulls common SNPs from gnomAD in a window either side, ranks them by
expected heterozygosity and proximity, annotates genetic distance from the bundled deCODE
map, and selects a panel covering both flanks. Exports to CSV, JSON, XLSX and PDF.

It proposes candidates. It cannot phase them: that needs the family.

## Running it

```sh
docker compose up -d --build
```

Then <http://localhost:8091>. No API keys are required. Optional keys are listed in
`deploy/README-deploy.md`; the app degrades gracefully without them.

Without Docker:

```sh
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --reload
cd web && npm install && npm run dev
```

## Tests

```sh
export PANELBUILDER_CACHE=tests/fixtures PANELBUILDER_CACHE_TTL=0
.venv/bin/python -m pytest tests/
cd web && npx tsc --noEmit && npm run build
```

The suite runs offline against recorded API responses in `tests/fixtures/` (gzipped).
`tests/genome_sweep.py` runs against live APIs and is not part of the suite.

## Documentation

Method, scope, data sources and known limitations are documented in the app:
[Documentation](https://ezrakruger.cc/originmarker/#/docs).

Version history, and the bugs each release fixed: [CHANGELOG.md](CHANGELOG.md).

## License

[Apache 2.0](LICENSE).

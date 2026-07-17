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

## Using it

1. **Type the variant.** An rsID (`rs334`), HGVS (`NM_000352.6(ABCC8):c.3989-9G>A`), a
   ClinVar accession (`VCV000009088`), or plain words (`the sickle cell mutation, in
   Africans`). Options go in the same line: `rs6025 in Europeans`, `rs334 with a 500kb
   window`. The `ex` button fills in a working example.
2. **Check what it resolved to.** Gene, coordinate, ClinVar's verdict. This is the step
   nobody can do for you: if it is the wrong variant, the panel will be a correct answer to
   a question you did not ask.
3. **Build the panel.** 20 to 60 seconds. Open the build log to watch it work.
4. **Read the coverage warnings, then download** CSV, XLSX, JSON or PDF.

The result is a list of candidates to genotype. You still have to genotype the carrier,
drop the markers where they are not heterozygous, and phase the rest against a relative.

Full walkthrough: [Using the site](https://ezrakruger.cc/originmarker/#/docs/using).

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
cd web && npm run build && for f in src/*.check.ts; do node "$f"; done
```

The suite runs offline against recorded API responses in `tests/fixtures/` (gzipped).
`tests/genome_sweep.py` runs against live APIs and is not part of the suite.

`npm run build` is `tsc -b && vite build`, and the `tsc -b` is the typecheck. This line used
to read `npx tsc --noEmit`, which checks **nothing**: `tsconfig.json` is a solution file with
`"files": []` and two references, so bare tsc has no inputs and exits 0 over a codebase that
does not compile. Only `-b` follows the references. Nothing shipped broken, because the build
ran the real check either way, but the reassurance was empty.

The `*.check.ts` files are assert-based self-checks over the pieces that are worth pinning:
the docs numbering, the log tags, the primer UI's honesty rules. They are plain node scripts.

## Documentation

Method, scope, data sources and known limitations are documented in the app:
[Documentation](https://ezrakruger.cc/originmarker/#/docs).

Version history, and the bugs each release fixed: [CHANGELOG.md](CHANGELOG.md).

## License

[Apache 2.0](LICENSE).

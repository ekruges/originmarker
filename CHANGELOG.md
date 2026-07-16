# Changelog

Releases are named after the biology of crossing over. All dates 2026.

The entries below are unusually detailed about bugs. That is deliberate: this tool's output
informs which parental allele an embryo inherited, so a plausible wrong answer is worse than
an error. An error gets investigated; a plausible wrong answer gets used. Anyone deciding
whether to trust a panel from an older build deserves to know exactly what it got wrong.

---

## 1.0.0 "Synaptonemal"

The protein scaffold that zips paired homologues together along their length.

**First public release.** Open source under Apache 2.0. Versioning restarts at semver here:
the 1.x/2.x numbers below were pre-release build counters, not compatibility promises.

### Added

- **Free-text input.** Two paths, and the difference is not cosmetic. Name an identifier
  (rsID, HGVS, ClinVar accession) and a regex reads it: no model runs, and it costs
  nothing. Describe the variant in words and a model is asked what you meant, answering
  from its own knowledge.
- **The model cannot touch a coordinate.** The typed query it fills has no field for a
  chromosome, position, strand or allele, so the coordinate on every panel came from a live
  lookup regardless of how the variant was named. This is a property of the code rather
  than a promise, and it was attacked specifically: extra coordinate keys are ignored,
  coordinate-shaped identifiers are refused by an allow-list.
- **A gene cross-check.** Any gene symbol in your own text is compared against the gene of
  the record that actually resolved. Ask for "the ABCC8 splice mutation" and get a variant
  in HBB, and it refuses rather than warns. It reads the resolved record, never the model's
  claim about itself: a model wrong about the variant can be wrong about its gene in the
  same breath.
- **Model provenance in every export.** A panel whose variant a model chose says so, in the
  CSV, JSON, XLSX and on the PDF variant card, and it says the two things separately: the
  model chose which variant, and it did not supply the coordinate. Previously a panel built
  from prose was byte-identical to one typed by hand.
- **An intent cache**, so identical prose is not billed twice.
- Documentation of all of the above, including what it cannot do.
- The monogram: a geometric sans ring around the wordmark's own serif M, outlined from the
  Merriweather glyph rather than redrawn. It is the site icon and the corner mark.
- Rotating input examples cycling the accepted forms.

### Fixed

- **The search box accepted one character at a time.** A component declared inside the
  render body took a fresh identity every render, so React remounted the whole subtree on
  each keystroke and destroyed the input's focus.
- **The gene extractor could not read `C9orf72`.** HGNC's Cxorfy convention carries
  lowercase, which a capitals-only pattern reads straight past. The commonest genetic cause
  of ALS was invisible to the safeguard.
- `MB` was excluded as a unit, but it is also myoglobin. Units are now stripped before gene
  extraction, so `500MB` is a window and `MB` is a gene.
- Ordinary report jargon (NIPT, CVS, WES, MLPA and others) read as gene symbols and refused
  correct answers, after billing for them.
- Legitimate aliases (`SUR1` for ABCC8, `ND1` for MT-ND1) were refused with no way past it.

### Known limitations

A gene symbol you lowercase reads as prose, and the cross-check then stays silent exactly
as if you had named no gene. Capitalisation is what separates a symbol from an ordinary
word here, and relaxing that makes every English word a symbol. Documented rather than
papered over.

---

## 2.4 "Bivalent"

A pair of synapsed homologues, held together by the chiasmata between them.

### Fixed

- **Default ranking used a best-case-over-ancestries statistic.** With no ancestry selected
  (the default), the primary sort key was the maximum expected heterozygosity across the
  eight gnomAD populations. That is not a prior for any carrier: it is the order statistic
  E[max], upward-biased by construction (+0.09 absolute, +26% relative on a representative
  locus), and each marker's figure came from a different population, so the ranking assumed
  a different ancestry for every row. One marker held a core slot advertising 2pq = 0.4884,
  the African value; for a Northern European family the same marker is 0.1418, and for an
  East Asian family 0.0241. Ranking is now on the global 2pq. On the reference panel this
  changed 15 of 20 shortlisted markers, took the count of markers with 2pq < 0.10 from 2 to
  0, and improved the worst-case marker in all eight populations.
- **Every export declared a ranking basis the engine did not use.** Exports stated "global
  2pq prior" while the sort keyed on the max across populations. The engine now names the
  quantity that produced the order, and exports render that name verbatim rather than
  restating it, so the two cannot drift apart again.
- **The PDF still labelled flanks `tel`/`cen`.** The engine moved to `lower coord` /
  `higher coord` in 2.3, but the PDF recomputed the label from the sign of the distance.
  This is the artifact that gets printed and filed, and it was anatomically inverted for
  every gene on a q arm.
- **`rs334` (sickle cell) resolved to a haplotype record.** ClinVar's relevance ranking
  returns eight haplotypes ahead of the HbS allele, each pairing it with a second HBB
  variant. A haplotype record resolves cleanly by returning one constituent's position with
  the haplotype's own classification, so the canonical sickle cell rsID produced a panel
  labelled "other" instead of "Pathogenic". Records describing combinations of variants
  (haplotype, compound heterozygote, diplotype) are now recognised by object type and
  refused with an explanation.
- **The client-side ranking mirror diverged from the engine.** `rank.ts` still fell back to
  the max-across-populations statistic, so the table headed "Engine rank" was ordered on it
  for all candidates while the PDF reported a different basis. Screen and export now agree
  for every row.
- **Concurrent builds shared one fetch ledger.** The ledger records the age of the data and
  supplies `queried_utc`, which is the data date printed on a filed export. It was a module
  global: thread-safe, but not build-safe. Two builds in one process shared it, so an export
  could carry the other panel's data date. Now scoped per build.
- **A failed ClinVar efetch was diagnosed as a structural variant.** The failure was
  swallowed, and the empty alleles fell through to the copy-number guard, refusing a routine
  SNV with a message that contradicted itself: "is single nucleotide variant (...), which
  has no single reference and alternate allele". An efetch failure now reports as a failed
  lookup.
- **A position disagreement between sources produced no panel-level warning.** The
  cross-check runs after selection, so a marker gnomAD and Ensembl place differently was
  already shortlisted when the disagreement surfaced, and it appeared only in that marker's
  own cell. It now raises a coverage flag, which renders in the UI and in all four exports.
- **Panels reported the live server's Ensembl release, not their own.** One frozen panel
  could render two different release numbers as the server moved on beneath it. The release
  is captured at build time; when it is unknown the panel says so rather than borrowing the
  current value.
- **PDF printed an allele frequency at full float precision** (`0.012718820176763365`) four
  lines above its own banner rendering the same quantity as `1.27e-2`.

### Changed

- Module self-checks now run under `pytest`. They existed and were never executed: with the
  ranking defect reintroduced, the suite still reported 54 passed. It now fails.
- Both test factories (Python and TypeScript) defaulted the global and max-across-population
  heterozygosities to the same number, making the two quantities indistinguishable by
  construction. They now default apart.
- Tests: 49 to 69.

---

## 2.3 "Holliday"

The four-armed branched junction intermediate that resolves into a crossover.

### Fixed

- **Genetic distances beyond the end of a chromosome map were clamped and reported as
  measurements.** The bundled maps stop short of the telomeres (chr11 by ~11 kb), and
  positions past the last map point collapsed onto a single value, reporting cM = 0 and
  recombination fraction = 0 sourced as deCODE. A fabricated theta of zero does not merely
  mislead: it ranks as a perfect marker. Positions off the map are now extrapolated at the
  documented 1 cM/Mb fallback and flagged as approximate.
- **`chrX` was labelled sex-averaged and is not.** The bundled chrX reads 1.297 cM/Mb
  against 1.175 (chr7) and 1.160 (chr8). A genuinely sex-averaged X must sit well below the
  autosomes, since paternal meioses are a third of X transmissions and contribute no
  crossovers outside the pseudoautosomal regions. It is the female map: correct for a female
  carrier, and an overstatement for a male one. It now says which.
- **Hotspot detection was off by one on both bounds**, so a hotspot outside the queried span
  could be reported for it.
- **Coverage was blind to recombination.** Eight of twenty shortlisted markers on the
  reference panel have a recombination hotspot between them and the variant, and coverage
  reported no flags. Sides are now judged on markers clear of an intervening hotspot, with
  "not assessed" kept distinct from "clear": an unmapped chromosome can never reach the
  hotspot threshold, so a negative there means unmeasured, not safe.
- **`tel`/`cen` flank labels were anatomically inverted for every gene on a q arm** (roughly
  60% of the genome, including CFTR, BRCA2, SMN1, MECP2, FMR1, and every acrocentric, whose
  p arm carries no genes). No centromere table exists anywhere in the tool, so it was
  asserting an anatomical fact it never looked up. Sides are now named by coordinate.
- **gnomAD quality filters were never requested**, so 33 QC-failed sites in the reference
  window were offered as candidates. Artifacts have inflated heterozygosity, so the ranking
  actively preferred them.
- **A failed 1000 Genomes join was reported as a count of zero**, i.e. a lie about the data
  rather than an admission that the lookup failed. Now reported as unavailable.
- **ClinVar significance was read from the first RCV rather than the variation-level
  aggregate**, and the significance and review status were scraped by independent passes, so
  they could describe different submissions.
- **`queried_utc` was wall-clock "now" even for a panel assembled entirely from cache.** The
  panel now reports both when it was built and how old the data in it is.
- **The natural-language rate limiter billed for refusals.** The global cap was tested after
  the model call, so a rejected request still cost money.
- **The free-text parser took the first identifier and discarded the rest**, so "not
  rs1801133, I mean rs151344623" built the wrong panel. It now refuses and names both.
- **The locus track's y-axis named an ancestry while plotting the global statistic.**
- **Rate limiting could be bypassed** by supplying an `X-Forwarded-For` header, since
  proxies append to it and the first entry is client-controlled.
- Error paths no longer leak module paths or raw provider errors to the browser.

### Added

- Mobile support.
- Recombination hotspot annotation and per-marker genetic distance from the bundled deCODE
  2019 map.

---

## 2.2 "Pachytene"

The substage of prophase I in which crossing over actually occurs.

### Fixed

- **Free text resolved to an arbitrary real variant.** ClinVar's esearch is a full-text
  search, and the code took the first hit on trust, so the phrase "a pathogenic variant"
  returned a real gene at real coordinates with no indication anything was wrong. Only
  identifiers are accepted now, and every candidate record must reconcile against the query
  before it is used.
- **The first ClinVar hit was not always the right one.** `rs1801133` (MTHFR, chr1) returned
  four records, the first of which is a CPS1 variant on chr2. All candidates are now fetched
  in one call and the one that reconciles is used.
- A copy-number variant resolved to empty alleles and a nonsense gnomAD identifier rather
  than being refused.
- The variant could appear as a marker for itself.
- Strand was reported as plus when it was unknown.
- A failed gnomAD chunk returned an empty list instead of raising.

### Added

- Genome-wide verification sweep against live data across all 23 chromosomes, both strands,
  and SNV/insertion/deletion/duplication classes, cross-checking every coordinate against a
  second source.

---

## 2.1 "Synapsis"

The lengthwise pairing of homologous chromosomes that has to happen before they can exchange
anything.

### Changed

- Two-tone wordmark.
- Removed the decorative helix.

---

## 1.0 "Chiasma"

The X-shaped point where two chromatids have crossed over: the visible evidence of
recombination, and the event every flanking marker is chosen to bracket.

First complete release.

### Added

- Panel builder: resolve a variant, enumerate common flanking SNPs from gnomAD, rank by
  expected heterozygosity and proximity, select a balanced panel covering both sides.
- Coordinates from live APIs only (ClinVar, Ensembl, gnomAD). No hardcoded positions.
- CSV, JSON, XLSX and PDF exports, each self-describing: build, both variant forms, source
  versions, timestamps, the wet-lab protocol, and the disclaimer.
- Progress streaming, docs with verified citations, optional LDlink annotation.
- Landing page, documentation.

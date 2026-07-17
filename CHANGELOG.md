# Changelog

Releases are named after the biology of crossing over. All dates 2026.

The entries below are unusually detailed about bugs. That is deliberate: this tool's output
informs which parental allele an embryo inherited, so a plausible wrong answer is worse than
an error. An error gets investigated; a plausible wrong answer gets used. Anyone deciding
whether to trust a panel from an older build deserves to know exactly what it got wrong.

---

## 1.4.0 "Tetrad"

The four chromatids of a paired chromosome, the unit crossover acts on.

Responds to an independent scientific audit of the methodology. The audit found no errors in
the statistics or genetics that is implemented; every change here closes a gap between what
the method assumes and the inputs it silently accepted, or improves a displayed figure.

### Added

- **mtDNA variants are refused at resolve time.** Mitochondrial DNA is maternally inherited,
  does not recombine, and is often heteroplasmic, so "which parental chromosome an embryo
  inherited" is undefined and flanking-SNP linkage does not apply. A chrM variant used to
  resolve like any point variant and hand back a thin panel that looked valid; it now returns
  a clear statement that the method does not apply.
- **A "cases this tool does not handle, or cannot detect" section** in Scope and limits:
  de novo variants, uniparental disomy, mosaicism, consanguinity, mtDNA and repeat expansions,
  each with the consequence for a panel. The tool never meets the family, so most of these it
  cannot detect and names as the reader's to rule out.

### Fixed

- **gnomAD exome frequencies were fetched and then dropped.** The single-variant query asked
  for both the genome and exome callsets but read only the 76k genomes. A coding pathogenic
  variant, which is most of this tool's input, carries its frequency in the 730k exomes and
  can be sparse or absent from the genomes, so the rarity card could read blank while a good
  exome answer sat unused. The verdict and the headline frequency now take whichever callset
  observed more chromosomes at the site, and both callsets are shown on the card and in every
  export. The LD-usable verdict is unchanged: it remains gated on the 1000 Genomes count.

### Changed

- The primer field reference notes that the 69 C Tm default is deliberately stringent, higher
  than routine genotyping, and can be lowered toward 60 C for a standard single-anneal assay.
- The star's 1 Mb clause is explained as inert at the default window (every candidate is
  already well inside 1 Mb) and binding only if the window is widened; ESHRE's 2 Mb
  "acceptable but not advisable" allowance is noted for that case.

---

## 1.3.4 "Diakinesis"

### Added

- A dismissable nudge at the top of the candidate list recommending a genome check on the
  primer pairs. Its button reopens the build log and streams the verification into it.
- A download button on the build log that exports it as a .txt, with a footer line carrying
  the release, job, instance URL, timestamp and data versions.

---

## 1.3.3 "Diakinesis"

### Fixed

- **Free-text descriptions never resolved.** The intent prompt forbade the model to recall
  an rsID, so "the sickle cell mutation, in Africans" returned no variant at all rather than
  rs334. It now names the standard identifier for a described variant, while still never
  emitting a coordinate; the identifier is confirmed by live lookup as before.

---

## 1.3.2 "Diakinesis"

Terms of use, brought up to what the app actually does. No behaviour changed.

### Changed

- **The terms disclose the two things that leave the server carrying your input**, which they
  did not. Section 5 was called "Third-party data" and described only what is retrieved:
  - **Free text with no identifier in it is sent verbatim to Anthropic's API**, to be read by
    a small model. Text containing an rsID, HGVS or a ClinVar accession is read here by a
    regular expression and reaches no model at all. The distinction is the whole privacy
    story of that box, the terms already asked readers to keep identifying information out of
    it, and they never said why.
  - **Primer sequences are sent to UCSC** when the check is asked for. They come from the
    reference genome rather than any sample, and nothing else from the query goes with them.
- **Primers are in the terms at all now.** Section 1 described a tool that proposes markers,
  which stopped being the whole of it in 1.3.0. A primer pair is a candidate in the same
  sense a marker is: nothing here has run a PCR, a reference is not a patient's genome, and a
  private variant under a primer site causes exactly the dropout the design cannot see.
- **Section 6 no longer lets "no tracking" imply "nothing leaves".** They are different
  claims and only the first was true. It says which is which, and that browser history stays
  in the browser.
- **Section 7 states the one dependency that is not Apache 2.0.** primer3-py is GPLv2, which
  is why it is optional and absent from the default image. Redistributing an image built with
  it switched on makes the combined work GPLv2. Previously the terms said "redistribute,
  including commercially" with no qualification.
- **Section 4 says whose UCSC quota the check spends.** It is the operator's key answering for
  UCSC's published limit, which is why it is rate limited per client.

### Added

- A check that fails when the app reaches a host the terms do not name, and when the terms
  stop saying a load-bearing thing ("not a clinical diagnostic", "candidate", "wet-lab",
  "GPLv2"). Prose drifts from code silently and in one direction: someone adds an outbound
  call, every test passes because the call works, and the terms keep saying it does not
  happen. It reads the rendered prose with the hrefs stripped, because the first version of
  it passed a gutted disclosure on the strength of a leftover URL.

---

## 1.3.1 "Diakinesis"

A patch, so it keeps 1.3.0's name.

### Added

- **The primer design is now reachable from Manual input**, which is where it was always
  meant to be. Pick which markers get a pair, and set every constraint they are designed
  under: the melting temperatures, the lengths, the composition, the product window, the
  reaction conditions the Tm is only meaningful beside, and the mask. The form seeds from the
  server's own numbers, so what is on screen is what the build is asked for.
- **A checkbox to check every pair against the genome as part of the build**, for anyone who
  would rather wait once than build and then press a second button. It states the cost beside
  itself: UCSC allows one query every 15 seconds, so it adds about 15 s per designed pair
  while the panel alone takes 20 to 60 s. The build log names each verdict as it lands, so a
  long run is not a blank wait, and it is off unless ticked, every time.

### Fixed

- **The Manual input primer form has never rendered for anyone.** It draws only against a
  server that states its defaults, which is right: the numbers must come from the engine that
  will use them rather than a copy in the browser. `/api/health` never sent them. So the
  section was gated on a field that did not exist, the flag beside it said primers were
  enabled, and the form was simply absent with nothing anywhere reporting a problem.
- **The build log claimed no pair had been checked, in builds that then checked them.** It
  said "none has been checked against the genome by this build", which was true when the
  design emitted it and false three lines later once the bundled check ran in the same job,
  about the same pairs, into the same console. It says "not yet" now, which is true on both
  paths.
- **A bundled check held a build slot for the whole of its wait**, found by watching the live
  site refuse an ordinary build while two ticked boxes sat waiting on UCSC. The slot bounds
  BUILDS, as its `MAX_CONCURRENT_BUILDS` name says, and the default is 2: two people ticking
  the box therefore blocked panel builds for everyone, for minutes, over work that had
  finished. The slot goes back when the build ends now. Verification stays bounded by the
  per-client budget and by the process-wide gate, which is what was keeping it polite anyway.
- **`npx tsc --noEmit` type-checks nothing, and the README told you to run it.** `tsconfig.json`
  is a solution file with `"files": []` and two references, so bare tsc has no inputs and
  exits 0 over a codebase that does not compile: only `tsc -b` follows the references. Proven
  rather than reasoned: a deliberate type error exits 0 under `--noEmit` and 2 under `-b`.
  Nothing shipped broken, because `npm run build` runs `tsc -b` and was in the same gate, but
  the line above it was pure reassurance. The README now documents the command that checks.

### Security

- **The bundled check spends the verification budget, not the build budget.** Both routes to
  UCSC reach one published daily quota, and the budgets differ: 20 builds a client per window
  against 4 verification runs. Charged only to the build budget, the checkbox would have been
  a five-fold rate-limit bypass on someone else's server. It is charged to the same key the
  button spends, and a test fails if that stops being true.

---

## 1.3.0 "Diakinesis"

The final condensation of prophase I, chiasmata still holding the homologues.

### Added

- **Primer design for the markers that meet the flanking criteria.** Each one gets a
  candidate FWD/REV pair for genotyping it by PCR, designed by primer3 against a reference
  template fetched around the marker, and folded into the PDF, CSV, XLSX and JSON. Defaults
  are 20 to 35 bases, GC 40 to 60%, Tm 69 C, product under 600 bp, GRCh38. Every field is
  settable, from the panel's own primer box or through the API, and the settings a panel was
  built under travel with it as provenance.
- **Common variants are masked out of both primer sites.** This is the part that matters. A
  primer sitting on a common SNP fails to bind in exactly the carriers who have it, their
  allele goes unamplified, and a heterozygote is read as a homozygote: allele dropout, which
  is silent and yields a genotype that looks clean and is wrong. The pool of markers is also
  the pool of hazards, so every gnomAD variant at or above the mask floor is excluded from
  under both primers, and the marker itself sits in the product under neither of them.
- **Optional verification against the whole genome**, through UCSC In-Silico PCR, behind a
  button and never part of a build: UCSC publishes one request every 15 seconds, so a public
  URL that verified on its own would spend the owner's quota on visitors who never looked at
  the result. It needs a UCSC API key (`deploy/README-deploy.md` says where to get one) and
  reports itself unavailable without one rather than offering a button that cannot run.
- **A primers chapter in the documentation**: every field with its bounds and what it
  constrains, what each warning means, what a clean in-silico result is and is not worth, and
  the steps to obtain a key. The field table is generated from the form's own field list, so
  a knob that reaches the form and not the docs fails the typecheck.

### Fixed

- **A verified pair kept the warning saying it had never been verified.** The verdict was
  written onto the pair as a bare state code, so a pair UCSC had called dangerous still
  carried "NOT CHECKED AGAINST THE GENOME" on the same row, and the reader was being asked to
  trust exactly one of two statements the same document made. The finding itself, naming the
  loci and their positions, reached only the build log: never the table, never the PDF. The
  verdict's own words are welded onto the pair now, and the caveat with them.
- **The check asked UCSC a narrower question than UCSC asks itself, and reported the answer
  as genome-wide.** Max product size went out at 1000 bp against hgPcr's own default of 4000.
  Measured against the live endpoint with a pair known to give one 549 bp product: at 400 the
  product is not reported at all. So the field bounds the search rather than filtering its
  result, and a pair whose second locus amplified between 1001 and 4000 bp came back holding
  one product, classified clean, and printed VERIFIED CLEAN (in silico) on a filed PDF. That
  is the multi-locus pass this lane exists to prevent, reached without a single component
  behaving incorrectly. It cut the other way too, since a design may ask for a 3000 bp
  product: that product could not be reported, and its absence classified as "found no
  product, do not order", which was our own request accusing a good pair. The question is
  UCSC's default now, the design has a server-side ceiling below it, and a check fails if
  either moves.
- **A page listing two loci could come back holding one.** The parser decided a line was a
  FASTA header before stripping the HTML around it, so a header whose ">" UCSC had wrapped or
  escaped was skipped as quietly as a sequence line, without being marked unreadable. The
  remaining product then classified as a clean single band. Tags come off and entities decode
  before that test now, so the test and the regex read the same text.
- **A UCSC timeout was rendered as a UCSC finding.** The page had no name for the state that
  means "asked, and the answer could not be read", so a timeout or a spent quota drew a red
  DANGER badge and a banner reading "In-silico PCR contradicts N primer pairs in this panel",
  over pairs UCSC had never answered about. The PDF from the same job id said NOT VERIFIED
  for all of them: one job, two documents, and two different instructions. Not verified is
  now neither dangerous nor clean, which are separate questions and were fixed together, as
  answering only the first renders a quota stop green.
- **The daily verdict on alt scaffolds overclaimed.** A hit on chr6 and one on a chr6 alt
  haplotype are usually one locus reported twice, and the note asserted flatly that the pair
  amplified more than one locus and must be redesigned. It says products rather than loci
  now, and names the ambiguity where alt or fix scaffolds are among the hits. The state stays
  DANGER: hgPcr cannot separate a redundant alt copy from a real second locus on that
  haplotype, and guessing toward clean is this tool's worst direction to guess in.
- **A UCSC key in `.env` never reached the container.** Compose reads that file for variable
  substitution and passes nothing it was not asked to, and the key was not named in the
  environment block. The key could be correctly generated, correctly stored, and the feature
  would still report itself unavailable, with no error anywhere: every layer behaving exactly
  as written. A test now fails if any key-shaped variable the app reads is not forwarded.
- **The JSON export shipped the shortlist twice in two different shapes.** The file writes it
  once inside `candidates` and again whole under `recommended`, and only the first was given
  its primers and its ancestry column, so a reader who took the shortlist rather than
  filtering the candidates got the lesser copy of the same marker.

### Changed

- **Every repeated block on the panel is one sentence and a link now.** Each warning carries
  a short form and a long one, written beside each other in the engine so they cannot drift.
  The table takes the short one; the exports keep the full wording, because a filed page
  cannot follow a link. The not-checked warning goes from 242 characters to 67, the mask note
  from 244 to 71, and the primer form's six paragraphs of constraint move to the docs.
- **The primer box is collapsed behind a thin line.** Four lines of detail under every marker
  is a table nobody reads. What collapsing never does is hide a finding: a dangerous pair and
  a failed design open themselves, the verdict sits on the summary line either way, and the
  panel-level alert lists every dangerous pair regardless of what is open.
- **The star's key is a key again.** It said what the criteria were, next to every star,
  beside a hover that already said it; it names the claim now and links to the chapter.
- Documentation section numbers derive from the section list, the way citation numbers always
  have. Inserting a chapter used to renumber eight headings and strand every cross-reference
  that named one.
- primer3-py stays an optional dependency and is not in `requirements.txt`: it is GPLv2 and
  this repo is Apache 2.0. Absent, panels build exactly as before and simply carry no
  primers. Nothing from kent/isPcr is vendored; verification calls UCSC's hosted service.

### Known limitations

- Verification is bounded at 4000 bp, which is UCSC's own bound and not an exhaustive search:
  a second locus amplifying wider than that is still invisible to it.
- A clean in-silico result is not a wet-lab validation. It does not model cycling conditions,
  and it cannot see a carrier's private variants under a primer site, which cause dropout
  exactly where the reference cannot show it.

---

## 1.2.0 "Zygotene"

The substage where homologues find each other and begin to pair.

### Added

- **Local query history.** Click the search box and your previous queries drop down, each
  with its candidate count and an x to forget it, plus Clear all. It lives in your browser
  and goes nowhere: no account, no server, nothing about anyone leaves the page. Storage
  that is full, disabled, or holding another version's garbage degrades to an empty list
  rather than taking the search box down.
- **A star on markers that meet ESHRE's flanking criteria**, in the table and in all four
  exports. What it means, exactly: within 1 Mb of the variant, no recombination hotspot in
  between on an assessed map, and no position disagreement between gnomAD and Ensembl. It
  is a structural check, not a ranking, and not a claim about any carrier's genotype.

### Why the star is a predicate and not a top-three

**There is no convention for a "top 3".** ESHRE's PGT-M recommendations (doi:
10.1093/hropen/hoaa018) do give real numbers, but a different shape: *at least three SNPs
proximal and three SNPs distal*, within 1 Mb of the variant, avoiding known hotspots. That
is six, per side, as a minimum count, not a ranking. ESHRE's own informativity rank is a
function of the couple's actual genotypes, which this tool structurally cannot have: it
proposes candidates and has no genotypes (R3), and 2pq is a population prior, not a
genotype (R4). So "strongest" in the only sense the convention defines is not computable
here, and claiming it would be an invented recommendation printed on a filed PDF.

**Starring the top 3 by the existing rank would have been actively harmful.** On the
reference ABCC8 panel, the three top-ranked markers are all on the same side and all three
have a recombination hotspot between them and the variant. Starring them would endorse
exactly the markers the tool's own coverage flag warns about, and violate the both-sides
rule on the same page that states it.

**So the star is a predicate, and it is not capped.** Every marker meeting the criteria
gets one; the count per side is reported against ESHRE's minimum of three. A cap would
force an ordering, and no source gives an exchange rate between heterozygosity and
distance: on the reference panel, capping at three would star a marker 24 kb away over one
125 bp away on a 0.01 difference in a population prior. The five unstarred qualifiers are
not worse than the three starred ones, and a star saying so would be a fabrication.

The label says "meets ESHRE's structural criteria". It never says strongest, best, or
recommended.

### Fixed

- **The history dropdown could submit a variant you never named.** Highlight a row, type a
  character that filters it out, press Enter, and the panel built for whichever row had
  slid into that index: a plausible wrong answer, produced silently. The selection index is
  a reference into a list that re-filters as you type, and resetting it in an effect runs
  after paint, which is too late for the keypress that follows. It resets in the same
  handler that changes the list now.
- **The arrow keys did nothing.** The input's own key handler and the combobox's collided,
  and the result was neither: navigation silently never ran. All keys are handled in one
  place now.
- **The export contradicted itself on 349 rows.** The criteria are evaluated over the
  shortlist, and the column is written for every candidate, so ~1,200 markers nobody judged
  shipped a False verdict. Several were nearer the variant than the starred ones, and the
  printed note routes a reader looking for a fourth marker straight into them. The field is
  None where it was not assessed, and the column prints empty: False is an assertion,
  absence is not.
- Four documentation cross-references named the wrong section number after the previous
  release inserted a section. A check now derives them from the section list.

---

## 1.1.0 "Leptotene"

The beginning of prophase I, where chromosomes first condense into threads.

### Added

- **A build log you can watch.** A dropdown beside the progress bar streams what the engine
  is actually doing, one tagged line per event: `[FETCH]` a request going out, `[CACHE]` one
  answered from disk and how old that answer is, `[INFO]` a count, `[SKIP]` sites dropped
  and why, `[WARN]` something worth reading, `[DONE]` the summary. It survives the build
  finishing, which is when it is most useful, and a late subscriber replays the whole log
  rather than joining midway.
- **The monogram on the PDF and XLSX masthead**, rendered from `web/public/favicon.svg` at
  export time. Page 1 only.

### Fixed

- **A provider URL, with its query string, reached the browser.** A failed call raised an
  error built from the raw URL, and the job's error text is both shown to the user and
  appended to the build log. The NCBI api_key travels in that query string. It did not in
  fact leak, but only because the URL was truncated at 80 characters and the key happens to
  be appended last: two couplings nothing asserted, one reordered dict away from publishing
  the key. The error is now built from the same scrubbed label the log uses, and a check
  drives a real failing call and asserts the key is absent.
- **`NM_000518.5(HBB):c.20A>T, MAF at least 0.1` did not resolve.** The search box decided
  client-side whether text was an identifier, and its HGVS test was a substring search, so
  any text containing `:c.` was posted whole, modifiers included, as the variant. The
  predicate is deleted: everything goes to the one parser, which is the server's. This also
  fixes free text being switched off refusing `rs334 in Europeans`, which needs no model.

### Changed

- The monogram's geometry is checked across its three copies. `favicon.svg` feeds the tab
  icon, the PDF and the spreadsheet; `Mark.tsx` holds a copy because its colours are CSS
  variables and the favicon's must be literal. Changing one used to leave the other alone,
  silently, with every gate green. The colours are still allowed to differ; the geometry is
  not, and a check says so.
- The log tag set is likewise pinned across Python and TypeScript.

### Performance

**No speedup shipped, and the measurements say why.** A cold build is 6.8 s, about 70% of it
waiting on the network. The candidates worth trying were tried and rejected on evidence:

- **The fan-out is CPU-bound, not network-bound.** Against a warm cache, one worker takes
  1,027 ms and eight take 887 ms: eight times the workers buys 1.16x. `json.loads` is 40 ms
  per 3.7 MB chunk and does not release the GIL, so ~890 ms is a floor no worker count
  removes. Sixteen workers is slower than eight.
- **Bigger requests are not fewer costs.** 40 kb chunks (7 requests instead of 13) measured
  *slower* than 20 kb (1.35 s vs 1.07 s): the cost is bytes and parsing, not round trips.
  10 kb was far worse at 5.97 s.
- **Nothing unread is worth removing.** The region query reads every field it asks for.
  Populations are 68% of the payload and gnomAD cannot filter them server-side.
- `annotate` is linear, and gnomAD's own run-to-run variance on an identical request is 8x
  (1.07 s vs 8.67 s), which is the noise floor any future claim here has to clear.

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

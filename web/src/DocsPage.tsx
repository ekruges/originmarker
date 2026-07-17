import { useEffect, type ReactNode } from 'react'
import { Alert, Anchor, Button, Code, Group, List, Paper, Table, Text, Title } from '@mantine/core'
import { CITATIONS, formatCitation } from './citations'
import { PRIMER_SIZE_CAP, type Health } from './api'
import {
  FIELD_KEYS, PRIMER_FIELDS, PRIMER_GROUPS, type FieldKey, type GroupKey,
} from './PrimerOptions'

export const REPO_URL = 'https://github.com/ekruges/originmarker'
export const HOME_URL = 'https://ezrakruger.cc/'
const CONTACT = 'kruger.ezra.s@gmail.com'

/** GitHub's mark, inlined. currentColor so one CSS rule drives its resting and hover grey. */
const GithubMark = () => (
  <svg viewBox="0 0 16 16" width="17" height="17" fill="currentColor" aria-hidden focusable="false">
    <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38
      0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01
      1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95
      0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27
      2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82
      2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0
      .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z" />
  </svg>
)

/** Person in a circle: the link home. Same 16-unit box and currentColor as the mark beside
 *  it, so one CSS rule greys both. The shoulders are clipped by the ring rather than drawn
 *  to meet it, which is what keeps the silhouette reading at 17px. */
const AvatarMark = () => (
  <svg viewBox="0 0 16 16" width="17" height="17" aria-hidden focusable="false">
    <defs>
      <clipPath id="om-av-clip">
        <circle cx="8" cy="8" r="7.25" />
      </clipPath>
    </defs>
    <circle cx="8" cy="8" r="7.25" fill="none" stroke="currentColor" strokeWidth="1.5" />
    <g clipPath="url(#om-av-clip)" fill="currentColor">
      <circle cx="8" cy="6.1" r="2.6" />
      <ellipse cx="8" cy="14.6" rx="4.7" ry="3.6" />
    </g>
  </svg>
)

// Numbering comes from CITATIONS key order, so inline markers and the reference list
// cannot drift apart: there is one ordering, not two.
const ORDER = Object.keys(CITATIONS)
const numberOf = (id: string) => ORDER.indexOf(id) + 1

/** In-page anchor href. The app routes on the hash, so a bare `href="#sources"` would
 *  replace the route and throw the reader back to the landing page: anchors stay
 *  namespaced under '#/docs/...'. */
const docHref = (id: string) => `#/docs/${id}`

/** The section id from a '#/docs/<id>' hash, or '' on the bare docs route. */
export function docSectionFromHash(hash: string): string {
  const m = /^#\/docs\/([\w:.-]+)$/.exec(hash)
  return m ? m[1] : ''
}

/** Inline citation marker: [n], linking to its entry in the reference list. */
function Ref({ id }: { id: string }) {
  const n = numberOf(id)
  if (n < 1) return null // an unknown id renders as nothing rather than a fake number
  return (
    <Anchor href={docHref(`ref-${id}`)} className="om-mono" style={{ fontSize: 11 }} aria-label={`reference ${n}`}>
      [{n}]
    </Anchor>
  )
}

const SECTIONS = [
  { id: 'what', label: 'What this is' },
  { id: 'using', label: 'Using the site' },
  { id: 'not', label: 'Scope and limits' },
  { id: 'pipeline', label: 'How a panel is built' },
  { id: 'star', label: 'The star: ESHRE flanking criteria' },
  { id: 'ld', label: 'Linkage disequilibrium and rare variants' },
  { id: 'prior', label: 'Expected heterozygosity is a prior' },
  { id: 'recomb', label: 'Recombination and the genetic map' },
  { id: 'layerb', label: 'Using the panel in the lab' },
  { id: 'primers', label: 'Primers: design, settings and checking' },
  { id: 'sources', label: 'Data sources and versions' },
  { id: 'conventions', label: 'Conventions' },
  { id: 'freetext', label: 'Free text and the model' },
  { id: 'limits', label: 'Known limitations' },
  { id: 'example', label: 'Worked example' },
  { id: 'api', label: 'API' },
  { id: 'references', label: 'References' },
]

/** A section's number, from its position in SECTIONS. The one place the ordering is stated,
 *  so headings, the nav and every cross-reference all read it rather than restate it. */
const sectionNo = (id: string) => SECTIONS.findIndex((s) => s.id === id) + 1

/** The constraint each primer group carries that its numbers cannot show. Lifted out of the
 *  form itself: they are worth reading once, not above every row of boxes forever. */
const GROUP_NOTE: Record<GroupKey, string> = {
  tm: 'A high target is a deliberate choice for a single band, and nothing relaxes it: a '
    + 'marker where no primer reaches it fails and names the knob, because a pair handed back '
    + 'below the target you set anneals at a temperature you will not run. The pair is held '
    + 'together because both primers anneal in the same tube.',
  size: 'Length floats to reach the Tm target: within a given GC range a fixed short oligo '
    + `cannot reach a high Tm, so a target and a fixed length cannot both hold. Capped at `
    + `${PRIMER_SIZE_CAP}, past which primer3's Tm model is no longer defined, and a Tm past `
    + 'it would be computed, plausible and out of model.',
  gc: 'GC content drives Tm, and a window whose own GC sits outside this range is where the '
    + 'design fails rather than compromises.',
  product: 'The template fetched around each marker widens with the maximum product. The '
    + 'server derives it from these, so it is not a separate knob: a template that did not '
    + 'widen with the product would starve the design instead of failing.',
  salt: 'Stated, never inherited. Tm depends on all four, and primer3 ships two default sets '
    + 'whose divalent and dNTP values differ and which compute a different Tm for the same '
    + 'oligo, several degrees apart. A Tm is only reproducible beside the conditions it was '
    + 'computed under, so these travel with every panel.',
  mask: 'Deliberately lower than the marker MAF floor: a variant too rare to be a useful '
    + 'marker still stops a primer binding in the carriers who have it, and the allele it '
    + 'should have amplified is then read as absent.',
}

/**
 * Every knob in the primer form, in one line each.
 *
 * Record<FieldKey, string>, so this is exhaustive by typecheck: a knob added to PRIMER_FIELDS
 * and not described here fails the build. That is the point. A field reference that silently
 * omits the field you are looking at is worse than none.
 */
const FIELD_DOC: Record<FieldKey, string> = {
  min_tm: 'No primer below this is accepted. The floor, not a preference.',
  opt_tm: 'What primer3 aims for and scores against.',
  max_tm: 'No primer above this is accepted.',
  max_pair_diff_tm: 'How far apart the two primers\' Tm may be. They anneal in the same tube, '
    + 'so a mismatched pair has no single annealing temperature that suits both.',
  min_size: 'Shortest oligo allowed.',
  opt_size: 'Length primer3 aims for, and drifts from to reach the Tm target.',
  max_size: `Longest oligo allowed. ${PRIMER_SIZE_CAP} is the ceiling: primer3's Tm model is `
    + 'not defined past it.',
  min_gc: 'Lowest GC fraction accepted in an oligo.',
  max_gc: 'Highest GC fraction accepted.',
  gc_clamp: 'How many of the 3\' end bases must be G or C. The 3\' end is where extension '
    + 'starts, so a clamp there stabilises the end that matters.',
  max_poly_x: 'Longest run of one base allowed anywhere in the oligo, e.g. 4 forbids AAAAA. '
    + 'Long runs misprime against other runs.',
  min_product: 'Shortest amplicon accepted.',
  max_product: 'Longest amplicon accepted. Also sets how much template is fetched.',
  salt_monovalent: 'Na+ / K+ concentration, for the salt correction to the Tm.',
  salt_divalent: 'Mg2+ concentration. Part of the same correction, and one of the two values '
    + 'primer3\'s own defaults disagree about.',
  dntp_conc: 'dNTP concentration. It binds Mg2+, so it enters the divalent correction.',
  dna_conc: 'Annealing oligo concentration, which enters the Tm directly.',
  mask_maf: 'Every gnomAD variant at or above this frequency in the window is kept out from '
    + 'under both primers.',
  target_pad: 'How far the marker itself is kept from either primer, so it sits in the '
    + 'product by a margin rather than by luck.',
}

// R8: a byte-for-byte copy of pb.DISCLAIMER, never a paraphrase. Used only when
// /api/health is unreachable.
const FALLBACK_DISCLAIMER =
  'Research use only. Candidate markers require validation and per-family phasing in a qualified genetics laboratory. Not a clinical diagnostic.'

const FALLBACK_STEPS = [
  'Genotype the carrier parent at these candidate markers.',
  'Keep only the markers where that carrier is actually heterozygous. Expected heterozygosity is a population average, not this individual’s genotype.',
  'Phase the retained markers against an informative relative (an affected or unaffected child, a proband, or a grandparent), or by read-based sequencing. OriginMarker cannot determine phase.',
  'Genotype the embryo biopsy at the phased markers.',
  'Require at least two concordant markers per side before calling parental origin.',
  'Use markers on both sides: one side alone cannot reveal a recombination between the marker and the locus, and per-side redundancy guards against allele dropout.',
]

const ENDPOINTS: [string, string][] = [
  ['GET /api/health', 'Versions, build, gnomAD dataset, genetic-map source, the disclaimer and the lab protocol steps. Feature flags for the optional LD annotation and free-text search.'],
  ['POST /api/resolve', 'Resolve an HGVS expression or rsID to a canonical variant record: rsID, gene, strand, GRCh38 coordinate, VCF ref/alt, clinical significance. Also returns the rarity verdict. Returns 400 on unresolvable input.'],
  ['POST /api/panel', 'Start a panel build. Returns 202 and a job_id; the build itself runs off the request path.'],
  ['GET /api/panel/{job_id}/stream', 'Server-sent events: progress frames per pipeline stage, then done or error. Heartbeats keep the stream alive behind a proxy.'],
  ['GET /api/panel/{job_id}', 'Poll a job. Returns status, current stage and fraction, and the full result once done. Works as a fallback where SSE is buffered.'],
  ['GET /api/export/{job_id}.{csv|xlsx|json|pdf}', 'Download a finished panel. Every export carries the build, the provenance stamp, the disclaimer and the lab protocol.'],
  ['POST /api/nl', 'Parse free text to a typed query: which variant, how wide a window, which ancestry. Intent only; it has no coordinate fields.'],
  ['GET /api/genes', 'Exact gene-symbol lookup for the search box.'],
  ['GET /api/ld', 'Optional LD annotation between two common SNPs. Refuses rare-variant queries.'],
]

export function DocsPage({ health }: { health: Health | null }) {
  const ensembl = health?.ensembl_release ?? '-'
  const gnomad = health?.gnomad_dataset ?? 'gnomad_r4'
  const build = health?.build ?? 'GRCh38'
  const mapSource = health?.map_source ?? 'deCODE 2019 sex-averaged (Beagle GRCh38 liftover, plink format)'
  const steps = health?.layer_b_steps?.length ? health.layer_b_steps : FALLBACK_STEPS

  // The hash is a route, not an element id, so the browser will not scroll for us.
  useEffect(() => {
    const jump = () => {
      const id = docSectionFromHash(window.location.hash)
      if (!id) return
      // Wait a frame: on a cold deep-link the section may not be mounted yet.
      requestAnimationFrame(() => {
        document.getElementById(id)?.scrollIntoView({ block: 'start', behavior: 'smooth' })
      })
    }
    jump()
    window.addEventListener('hashchange', jump)
    return () => window.removeEventListener('hashchange', jump)
  }, [])

  return (
    <div className="om-docs-wrap" style={{ display: 'flex', gap: 24, maxWidth: 1100, margin: '0 auto', padding: 12, alignItems: 'flex-start' }}>
      <nav
        className="om-docs-nav"
        aria-label="Documentation sections"
        style={{ position: 'sticky', top: 12, flex: '0 0 200px', alignSelf: 'flex-start' }}
      >
        <ol style={{ listStyle: 'none', margin: 0, padding: 0 }}>
          {SECTIONS.map((s, i) => (
            <li key={s.id}>
              <a href={docHref(s.id)}>
                <span className="om-mono" style={{ marginRight: 6 }}>
                  {i + 1}
                </span>
                {s.label}
              </a>
            </li>
          ))}
        </ol>
        <Text size="xs" c="dimmed" mt={10} pl={8} className="om-mono">
          {build} · {gnomad} · Ensembl {ensembl}
        </Text>

        <div className="om-docs-links">
          <a href={REPO_URL} target="_blank" rel="noreferrer" aria-label="Source on GitHub" title="Source on GitHub">
            <GithubMark />
          </a>
          <a href={HOME_URL} aria-label="ezrakruger.cc" title="ezrakruger.cc">
            <AvatarMark />
          </a>
        </div>
      </nav>

      <article className="om-docs-body" style={{ flex: 1, minWidth: 0 }}>
        <Title order={1} mb={4}>
          OriginMarker documentation
        </Title>
        <Text size="xs" c="dimmed" mb="md">
          Candidate flanking-SNP panels for PGT-M linkage, built from population data. Research decision
          support, not a diagnostic.
        </Text>

        <Section id="what" title="What this is">
          <Text mb={8}>
            A carrier parent transmits either the wild-type or the mutant allele. After an embryo is
            edited, a wild-type read at the variant site is ambiguous: the egg may have been
            wild-type-fertilized, in which case the editor did nothing, or mutant-fertilized and
            corrected. Reading the variant site alone cannot separate those two histories.
          </Text>
          <Text mb={8}>
            SNPs flanking the variant do separate them. They were not the editing target, so genotyping
            them reports which parental chromosome came in regardless of what happened at the variant
            itself. That is established PGT-M linkage and karyomapping methodology{' '}
            <Ref id="karyomapping" /> <Ref id="eshre_pgt_m" />, repurposed here for editing attribution.
          </Text>
          <Text>
            OriginMarker does the part population data can answer: pulling candidate SNPs from gnomAD,
            annotating their frequencies, ranking them and laying them out against the locus, which is
            work a scientist otherwise does by hand across gnomAD, Ensembl and ClinVar. Genotyping the
            carrier and phasing the markers in the family is bench work, and no app can do it. The output
            is a menu of candidates.
          </Text>
        </Section>

        <Section id="using" title="Using the site">
          <Text mb={10}>
            Type the variant, check what it resolved to, build the panel, read the warnings,
            download it. The one step nobody can do for you is the second: confirming that
            the record on screen is the variant you actually meant.
          </Text>

          <Title order={3} mt={14} mb={4}>1. Name the variant</Title>
          <Text mb={6}>
            One box, five accepted forms. The first four are read by a regex, cost nothing,
            and are exact:
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Tbody>
                <ExRow k="rsID" v="rs334" mono />
                <ExRow k="HGVS, with a transcript" v="NM_000352.6(ABCC8):c.3989-9G>A" mono />
                <ExRow k="HGVS, gene omitted" v="NM_000352.6:c.3989-9G>A" mono />
                <ExRow k="ClinVar accession" v="VCV000009088" mono />
                <ExRow k="Plain words" v="the sickle cell mutation, in Africans" />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={6}>
            The last one is different in kind: no identifier is in the text, so a language
            model is asked what you meant and answers from its own knowledge. That is a
            claim, not a reading, and it is fenced accordingly. See{' '}
            <Anchor href={docHref('freetext')}>Free text and the model</Anchor> before you
            rely on it. If you know the rsID, type it: exact, instant, free.
          </Text>
          <Text mb={6}>
            Options go in the same line, and are always read locally whichever form you use:
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Tbody>
                <ExRow k="Ancestry" v="rs6025 in Europeans" mono />
                <ExRow k="Window" v="rs334 with a 500kb window" mono />
                <ExRow k="MAF floor" v="NM_000518.5(HBB):c.20A>T, MAF at least 0.1" mono />
                <ExRow k="Together" v="rs113993960 in East Asians, 100kb, MAF 0.1" mono />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            The <Code>ex</Code> button fills in a random working example. Prefer knobs to
            prose? <b>Manual input</b> exposes every parameter as its own field: window, MAF
            floor, ancestry, build, whether to cross-check positions against Ensembl, and a{' '}
            <b>Primers</b> section with the whole design in reach: which markers get a pair,
            every constraint it is designed under, and a checkbox to check each pair against
            the genome as part of the build rather than afterwards. That check is much slower
            than the panel, so it is off unless you tick it (<SecRef id="primers" />).
          </Text>

          <Title order={3} mt={14} mb={4}>2. Check the resolved variant. This step is yours</Title>
          <Text mb={6}>
            Resolve looks the identifier up live and shows you the record: gene and strand,
            the GRCh38 coordinate, the variant in both genomic and transcript-sense form,
            ClinVar's classification and review status, the gnomAD frequency, and whether
            population LD is defined for it. Nothing is built until you press the button, so
            this is the checkpoint.
          </Text>
          <Text mb={6}>
            <b>Read the gene and the coordinate and confirm they are the variant you mean.</b>{' '}
            Everything downstream is arithmetic on that record: if it is the wrong variant,
            the panel will be a technically perfect answer to a question you did not ask.
            That is the failure this tool cannot catch for you, and it is likelier than it
            sounds when a gene has many pathogenic variants.
          </Text>
          <Text mb={8}>
            Refusals here are deliberate. A haplotype record, a structural variant, prose
            that names no identifier, or a genomic position typed as a coordinate are all
            turned away rather than approximated: the message says what to type instead.
          </Text>

          <Title order={3} mt={14} mb={4}>3. Build the panel</Title>
          <Text mb={6}>
            Typically 20 to 60 seconds, because it pulls every common SNP in the window from
            gnomAD, which for a 250 kb window is a few hundred thousand rows. Open{' '}
            <b>Build log</b> to watch it work. Each line carries a tag:
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Tbody>
                <ExRow k="[FETCH]" v="a request going out to ClinVar, Ensembl or gnomAD" />
                <ExRow k="[CACHE]" v="answered from disk instead, and how old that answer is" />
                <ExRow k="[INFO]" v="a count: candidates, markers shortlisted, positions confirmed" />
                <ExRow k="[SKIP]" v="sites dropped, and the reason: MAF floor, QC filters, monomorphic" />
                <ExRow k="[WARN]" v="something worth reading before you use the panel" />
                <ExRow k="[DONE]" v="the summary: markers, candidates, elapsed, how much came off the network" />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            The log is worth opening when a build is slow: it tells you whether you are
            waiting on gnomAD or on your own window size. It stays available after the build
            finishes.
          </Text>

          <Title order={3} mt={14} mb={4}>4. Read the panel, warnings first</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              <b>Flanking coverage</b> is the first thing to read. It says how many markers
              sit each side, how many are close, and it raises a flag when a side is thin,
              empty, or separated from the variant by a recombination hotspot. A red flag
              here is not cosmetic: it is the tool telling you this panel may not do the job.
            </List.Item>
            <List.Item>
              <b>The locus track</b> plots each candidate against the variant: height is the
              expected heterozygosity, colour is the side.
            </List.Item>
            <List.Item>
              <b>The table</b> lists every candidate, shortlisted ones marked. Sort it, filter
              it, or apply a preset. <i>Engine rank</i> is the order the panel was actually
              selected in.
            </List.Item>
            <List.Item>
              <b>Ancestry</b> re-orders the page by that population's own frequencies. It does
              not re-select the shortlist: to do that, rebuild.
            </List.Item>
            <List.Item>
              <b>Provenance</b> records when the panel was built, how old its data is, which
              sources answered, and what the ranking keyed on. If <Code>queried_utc</Code> is
              much older than <Code>built_utc</Code>, the panel came from cache and rebuilding
              it is not a re-check.
            </List.Item>
          </List>

          <Title order={3} mt={14} mb={4}>5. Download it</Title>
          <Text mb={8}>
            <b>CSV</b> for a pipeline, <b>XLSX</b> for a spreadsheet, <b>JSON</b> for code,{' '}
            <b>PDF</b> to read or file. All four carry the same facts: the variant in both
            forms, the sources and their versions, the timestamps, the wet-lab protocol and
            the disclaimer. One caveat worth knowing:{' '}
            <Code>pandas.read_csv(comment='#')</Code> strips the CSV's entire header block,
            caveats included. The JSON and XLSX carry those as data.
          </Text>

          <Alert color="gray" variant="light" title="What you have at the end">
            <Text size="xs">
              A list of candidates to genotype, not an answer. The panel says which SNPs are
              worth typing in this family; it cannot say which parental allele an embryo
              inherited. That needs the carrier genotyped, the uninformative markers dropped,
              and the rest phased against a relative.{' '}
              <Anchor href={docHref('layerb')}>Using the panel in the lab</Anchor> is the
              protocol.
            </Text>
          </Alert>
        </Section>

        <Section id="not" title="Scope and limits">
          <List spacing={4} mb={8}>
            <List.Item>No editing-reagent, guide-RNA or gamete design.</List.Item>
            <List.Item>
              No phase determination. Phase requires family samples or read-based sequencing; the app
              cannot infer it.
            </List.Item>
            <List.Item>No clinical interpretation and no diagnosis.</List.Item>
            <List.Item>No patient data storage. Nothing about a family is submitted or retained.</List.Item>
            <List.Item>
              No claim that a marker is informative in a specific family. Population frequency cannot
              establish that; only genotyping the actual carrier can.
            </List.Item>
          </List>
          <Alert color="gray" variant="light" title="Disclaimer">
            <Text size="xs">{health?.disclaimer ?? FALLBACK_DISCLAIMER}</Text>
          </Alert>
        </Section>

        <Section id="pipeline" title="How a panel is built">
          <List spacing={5} type="ordered" mb={10}>
            <List.Item>
              <b>Resolve.</b> The variant is looked up live against ClinVar, falling back to Ensembl for a
              bare rsID, producing a canonical record: rsID, gene, strand, {build} coordinate, VCF ref/alt,
              clinical significance, review status, accession. Input that cannot be resolved returns an
              error; no panel is built on an approximation.
            </List.Item>
            <List.Item>
              <b>Rarity.</b> The gnomAD single-variant record and the 1000 Genomes allele count decide
              whether population LD with this allele is even defined. For a pathogenic variant the answer
              is almost always no (see <SecRef id="ld" />).
            </List.Item>
            <List.Item>
              <b>Enumerate.</b> A ±250 kb window by default, pulled from gnomAD as 20 kb region chunks
              fetched concurrently. Only SNPs are kept, and only where gnomAD genome AN ≥ 10,000, which
              discards sites whose frequency rests on a poor call rate, and where MAF ≥ 0.05, which is the
              floor for a marker worth genotyping.
            </List.Item>
            <List.Item>
              <b>Annotate.</b> Global and per-ancestry MAF; expected heterozygosity 2pq; signed distance
              from the variant; genomic side; tier; and cM, recombination fraction and a hotspot flag
              interpolated from the bundled genetic map.
            </List.Item>
            <List.Item>
              <b>Rank and select.</b> Ranked by expected heterozygosity, then by proximity. A balanced set
              is chosen from distance bands on each side independently, so a dense cluster on one side
              cannot crowd the other out, plus distant sentinels near the window edge.
            </List.Item>
            <List.Item>
              <b>Cross-check.</b> The nearest markers' positions are re-verified against Ensembl, an
              independent source from the one that supplied them. A disagreement is flagged.
            </List.Item>
          </List>
          <Text mb={8}>
            <b>Tiers</b> group candidates by absolute distance: <b>A</b> core, under 2 kb; <b>B</b> near,
            2–30 kb; <b>C</b> flank, 30 kb and beyond. Closer is better, since there is less chance of a
            recombination between marker and locus, but a panel of only the closest markers is fragile:
            one nearby crossover can mislead all of them at once. Distant sentinels catch that. A marker
            that disagrees with its near neighbours on the same side localizes the crossover.
          </Text>
          <Text>
            Selection is balanced across both genomic sides. A one-sided panel cannot detect a
            recombination between the variant and the markers, so a side with fewer than two markers within
            30 kb raises a coverage flag.
          </Text>
        </Section>

        <Section id="star" title="The star: ESHRE flanking criteria">
          <Text mb={8}>
            A star beside a shortlisted marker means it meets three structural criteria drawn from
            the ESHRE PGT Consortium's good practice recommendations{' '}
            <Ref id="eshre_pgt_m" />. It is this tool's own check of those criteria. ESHRE has not
            reviewed this tool.
          </Text>

          <Title order={3} mt={12} mb={4}>What a star means, exactly</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              <b>Within 1 Mb of the variant.</b> ESHRE recommends staying inside 1 Mb (about 1 cM)
              of the pathogenic variant, because loci 1 cM apart are expected to recombine about 1%
              of the time.
            </List.Item>
            <List.Item>
              <b>No recombination hotspot between the marker and the variant</b>, on the bundled
              deCODE map. A marker separated from the locus by a hotspot is the one most likely to
              have lost phase with it, which is the failure this criterion exists to avoid.
            </List.Item>
            <List.Item>
              <b>Its GRCh38 position is not disputed</b> between gnomAD and Ensembl. Note the
              difference between disputed and unchecked: only the nearest few markers are ever
              cross-checked, so most carry no verdict, and no verdict is not a mark against them.
            </List.Item>
          </List>

          <Title order={3} mt={12} mb={4}>What a star does not mean</Title>
          <Text mb={6}>
            <b>It is not an informativity rank, and it is not a ranking of any kind.</b> ESHRE's
            marker informativity (Tables I and II of that paper) is computed from the couple's and
            their relatives' actual genotypes. This tool has no genotypes: it proposes candidates
            and cannot phase them, and 2pq is a population prior rather than a claim about your
            carrier. So the only sense in which ESHRE defines a "best" marker is one this tool
            structurally cannot compute, and it does not pretend to.
          </Text>
          <Text mb={6}>
            <b>Starred markers are not ordered</b>, and no starred marker is preferred over another.
            There is no cap: every marker meeting the criteria is starred. A cap would force an
            ordering, and no published source gives an exchange rate between heterozygosity and
            distance, so any such order would be this tool's invention presented as a
            recommendation.
          </Text>
          <Text mb={6}>
            <b>An unstarred marker is not unusable.</b> ESHRE recommends at least three SNPs on each
            side of the locus, and unstarred markers may well be needed to reach that. A star is a
            property of one marker; sufficiency is a property of the panel, and that is what{' '}
            <Anchor href={docHref('layerb')}>the lab protocol</Anchor> and the coverage warnings
            are for.
          </Text>

          <Alert color="gray" variant="light" title="In short">
            <Text size="xs">
              A star says the marker is structurally well placed: close enough, no hotspot in the
              way, and its position agreed on. It says nothing about whether your carrier is
              heterozygous there, which only genotyping can answer.
            </Text>
          </Alert>
        </Section>

        <Section id="ld" title="Linkage disequilibrium and rare variants">
          <Text mb={8}>
            Linkage disequilibrium (r² or D′) between a marker and the pathogenic variant is{' '}
            <b>undefined</b> when the pathogenic variant is rare. LD is a property of a haplotype
            frequency distribution estimated from a reference panel <Ref id="thousand_g" />. If the
            pathogenic allele appears once in that panel, as in the worked example below, there is no
            distribution to estimate from. The statistic is not weak or noisy; it does not exist.
          </Text>
          <Text mb={8}>
            So LD cannot rank candidate markers here and cannot call parental origin. The ranking key
            never sees an LD value. Origin comes from per-family phasing.
          </Text>
          <Text>
            The optional LDlink annotation <Ref id="ldlink" /> is a different measurement: a labelled
            prior <b>between two common SNPs</b>, which says only whether two markers carry redundant
            information. It is never computed against the pathogenic variant; the query is refused when
            either allele is rare.
          </Text>
        </Section>

        <Section id="prior" title="Expected heterozygosity is a prior">
          <Text mb={8}>
            Expected heterozygosity, 2pq under Hardy-Weinberg, is the probability that a randomly drawn
            individual from a reference population is heterozygous at the site. That is the entire content
            of the number.
          </Text>
          <Text mb={8}>
            It is <b>not</b> a statement about your carrier. A marker with 2pq = 0.5 is not half
            heterozygous in one person; that person is heterozygous or homozygous, and population data
            cannot say which. A marker is informative for a family only if the carrier is in fact
            heterozygous at it. Otherwise the marker reports nothing about which chromosome was
            transmitted, whatever its frequency in gnomAD.
          </Text>
          <Text>
            Ranking by 2pq maximizes the chance a candidate survives carrier genotyping. It does not
            substitute for it.
          </Text>
        </Section>

        <Section id="recomb" title="Recombination and the genetic map">
          <Text mb={8}>
            A crossover between a marker and the pathogenic allele makes the marker report the wrong
            parental chromosome. Each candidate carries its cM distance from the variant, interpolated
            from the bundled map <Ref id="decode_map" />. Physical distance in bp is a poor proxy for that
            risk: recombination rate varies by orders of magnitude along the genome.
          </Text>
          <Text mb={8}>
            The recombination fraction θ is derived from cM by Haldane's map function <Ref id="haldane" />,
            θ = 0.5(1 − e<sup>−2d</sup>) with d in Morgans. Haldane assumes no interference, so at the same
            cM it returns a slightly higher θ than Kosambi <Ref id="kosambi" />. That is the conservative
            direction: it overstates the risk of a marker rather than understating it.
          </Text>
          <Text mb={8}>
            A hotspot between the marker and the variant is flagged explicitly rather than left implicit
            in the cM value.
          </Text>
          <Text>
            The map is <b>population-averaged</b>. It describes expected recombination across many meioses
            in a reference cohort, not what happened in the meiosis that produced the embryo being tested.
            Like 2pq, it is a prior. A low recombination fraction means a marker is likely to have stayed
            in phase; it is not a guarantee that it did. That is why the protocol below requires two
            concordant markers per side, and both sides.
          </Text>
        </Section>

        <Section id="layerb" title="Using the panel in the lab">
          <Text mb={8}>
            The panel is a starting point. To actually call parental origin:
          </Text>
          <List spacing={5} type="ordered" mb={10}>
            {steps.map((s, i) => (
              <List.Item key={i}>{s}</List.Item>
            ))}
          </List>
          <Text mb={8}>
            Two markers per side means a genotyping error, a mis-assigned phase or a crossover shows up as
            disagreement instead of a confident wrong answer. If the two sides disagree, a crossover
            occurred between them.
          </Text>
          <Text>
            Per-side redundancy also covers allele dropout, in which one allele fails to amplify from a
            single-cell or few-cell biopsy and a heterozygote is read as a homozygote <Ref id="ado" />. A
            dropout is silent: the genotype it yields looks clean and is wrong. Only a second marker on the
            same side, disagreeing with it, catches that.
          </Text>
        </Section>

        <Section id="primers" title="Primers: design, settings and checking">
          <Text mb={8}>
            Where a panel carries primers, each one is a candidate FWD/REV pair for genotyping
            that marker by PCR. Primer3 designs them against a reference template fetched
            around the marker. They are candidates in the same sense the markers are: nothing
            here has run a PCR, and every pair needs validating at the bench before use.
          </Text>
          <Text mb={8}>
            <b>Common variants are masked out of both primer sites.</b> This is the part worth
            reading. A primer sitting on a common SNP fails to bind in exactly the carriers who
            have that SNP, so their allele goes unamplified and a heterozygote is read as a
            homozygote. That is allele dropout <Ref id="ado" />, it is silent, and it is the
            worst error this tool can contribute to: the genotype it yields looks clean and is
            wrong. The pool of markers is also the pool of hazards, so every gnomAD variant at
            or above the mask MAF floor in the window is excluded from under both primers, and
            the marker itself sits in the product under neither of them.
          </Text>

          <Title order={3} mt={16} mb={4}>What the warnings mean</Title>
          <Text mb={8}>
            Every pair carries at least one note, always. A pair with nothing said about it
            would read as a pair with nothing wrong with it, so the design never hands one
            back. On the panel each note is one line; the full wording is in the PDF, CSV,
            XLSX and JSON exports, which are read away from this page.
          </Text>
          <Wide>
          <Table className="om-table" withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>Note</Table.Th>
                <Table.Th>What it means</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              <Table.Tr>
                <Table.Td><b>Not checked against the genome</b></Table.Td>
                <Table.Td>
                  The default, and the state of every pair until you ask for the check below.
                  Primer3 saw one 800 bp template and cannot see the other 3.1 Gb, so the pair
                  may prime a second locus somewhere this build has not looked. Not a verdict,
                  and not a pass.
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td><b>DANGER</b></Table.Td>
                <Table.Td>
                  UCSC was asked and the answer was bad: several products, none at all, one on
                  the wrong chromosome, or one at the wrong size. A pair that amplifies more
                  than one locus cannot genotype the marker, because the trace is a mixture and
                  a heterozygote is indistinguishable from two loci differing at that base. The
                  pair is still shown, with the finding: a hidden primer is a decision made for
                  you, a warned one is information.
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td><b>One product</b></Table.Td>
                <Table.Td>
                  UCSC found exactly one product, on the right chromosome, at the designed
                  size. This is the best result available here and it is still not a validation:
                  see the caveat below.
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td><b>Still unverified</b></Table.Td>
                <Table.Td>
                  UCSC was asked and the answer could not be read: an unreadable page, a
                  timeout, a spent quota. An unreadable answer is not a clean answer, so it
                  stays unverified rather than becoming a pass.
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td><b>No pair for this marker</b></Table.Td>
                <Table.Td>
                  The design ran and found nothing meeting your settings, and it says which
                  constraint it could not satisfy. Nothing is relaxed to force a pair: a pair
                  handed back below the Tm you set anneals at a temperature you will not run.
                  A window at 63% GC will not yield a 40-60% GC primer, and that is the honest
                  answer, not a failure to try.
                </Table.Td>
              </Table.Tr>
            </Table.Tbody>
          </Table>
          </Wide>

          <Title order={3} mt={16} mb={4}>Choosing what gets a pair, and tuning the design</Title>
          <Text mb={8}>
            By default only the starred markers get primers, which keeps a build fast and
            covers the markers most likely to be ordered. Switch the query to{' '}
            <b>Manual input</b> and a <b>Primers</b> section appears with the whole design
            exposed: the scope, and every field in the reference table below.
          </Text>
          <List spacing={4} mb={10}>
            <List.Item>
              <b>Design primers for</b> chooses the scope. <i>Markers meeting the flanking
              criteria</i> is the default. <i>Every shortlisted marker</i> designs for the
              whole shortlist, which is more pairs and a longer build. <i>No primers</i>{' '}
              builds the panel alone.
            </List.Item>
            <List.Item>
              <b>Every field is yours.</b> The numbers the form opens with are the server's
              own, fetched from the engine that will use them rather than kept in the browser,
              so what you see is what a build is asked for. Whatever you set travels with the
              panel as provenance, and the primer box on the result restates it.
            </List.Item>
            <List.Item>
              <b>Nothing is relaxed to force a pair.</b> Ask for something a window cannot
              give, and that marker fails and names the constraint it could not meet. A pair
              handed back below the Tm you set would anneal at a temperature you will not run.
            </List.Item>
          </List>
          <Text mb={8}>
            The same form is on every finished panel, under the <b>Primers</b> chip: read-only
            there if the page cannot rebuild, and a <b>Rebuild</b> otherwise.
          </Text>

          <Title order={3} mt={16} mb={4}>In-silico PCR: what a pass is worth</Title>
          <Text mb={8}>
            Checking is optional and off by default. Run it two ways: the <b>Check pairs</b>{' '}
            button in the primer box on a finished panel, or the checkbox in Manual input to
            fold it into the build itself. Either way it sends each pair to UCSC In-Silico PCR,
            which aligns the two primers against the whole GRCh38 reference and reports where
            they would amplify it.
          </Text>
          <Text mb={8}>
            <b>It is much slower than the build.</b> UCSC publishes one request every 15
            seconds for programmatic use, so each pair adds about 15 seconds: eleven pairs is
            roughly three minutes, against 20 to 60 seconds for the panel itself. Bundling it
            means one wait instead of two; leaving it off means the panel is on screen while
            you decide. The build log names each verdict as it lands, so a bundled run is not
            a blank wait. Nothing is ever checked unless you ask: a public URL that verified
            on its own would spend the owner's daily quota on visitors who never read the
            result, and both routes draw on the same per-client budget.
          </Text>
          <Alert color="yellow" p={8} mb={8}>
            <Text size="sm">
              <b>A clean result is not a wet-lab validation.</b> It is not a PCR. It does not
              see this carrier's genome, so a private variant under a primer site will not
              appear and can still cause dropout. It does not model your cycling conditions, so
              one reference product is not a promise of one band on a gel. It means the pair is
              not obviously multi-locus against one reference sequence, and nothing more.
            </Text>
          </Alert>

          <Title order={3} mt={16} mb={4}>Turning the check on: a UCSC API key</Title>
          <Text mb={8}>
            {health?.insilico_pcr_enabled
              ? 'This instance has a key configured, so the Check pairs button is live.'
              : 'This instance has no key configured, so pairs stay marked as not checked. '
                + 'If you are running your own copy:'}
          </Text>
          <List spacing={4} type="ordered" mb={10}>
            <List.Item>
              Make a UCSC Genome Browser account, or log in to yours, at{' '}
              <Anchor href="https://genome.ucsc.edu/cgi-bin/hgLogin" target="_blank" rel="noreferrer">
                genome.ucsc.edu/cgi-bin/hgLogin
              </Anchor>.
            </List.Item>
            <List.Item>
              Open the{' '}
              <Anchor href="https://genome.ucsc.edu/cgi-bin/hgHubConnect#dev" target="_blank" rel="noreferrer">
                Hub Development page
              </Anchor>{' '}
              and generate a key in the <b>API key</b> section at the bottom.
            </List.Item>
            <List.Item>
              Set it as <Code>UCSC_API_KEY</Code> in the server's environment and restart. It
              is read server-side only and never reaches the browser.
            </List.Item>
          </List>
          <Text mb={8} size="sm" c="dimmed">
            Keys work only on the primary site, not the genome-euro or genome-asia mirrors.
            Without a key, hgPcr answers with a CAPTCHA, which this app classifies as
            unverified and leaves that way; it does not attempt to solve it. UCSC's published
            limit for programmatic use is one request every 15 seconds and 5,000 a day, and
            those are UCSC's numbers rather than tuning knobs, so they are not settable here.
          </Text>

          <Title order={3} mt={16} mb={4}>Field reference</Title>
          <Text mb={8}>
            Every knob the primer box draws, with the bounds it accepts. Defaults come from the
            server that will use them, so the numbers in the form are the engine's, not a copy.
          </Text>
          {PRIMER_GROUPS.map((g) => (
            <div key={g.key}>
              <Title order={4} mt={12} mb={3}>{g.title}</Title>
              <Text mb={6} size="sm">{GROUP_NOTE[g.key]}</Text>
              <Wide>
              <Table className="om-table" withTableBorder>
                <Table.Thead>
                  <Table.Tr>
                    <Table.Th>Field</Table.Th>
                    <Table.Th>Range</Table.Th>
                    <Table.Th>What it does</Table.Th>
                  </Table.Tr>
                </Table.Thead>
                <Table.Tbody>
                  {FIELD_KEYS.filter((k) => PRIMER_FIELDS[k].group === g.key).map((k) => (
                    <Table.Tr key={k}>
                      <Table.Td className="om-mono">{PRIMER_FIELDS[k].label}</Table.Td>
                      <Table.Td className="om-mono" style={{ whiteSpace: 'nowrap' }}>
                        {PRIMER_FIELDS[k].min} to {PRIMER_FIELDS[k].max}
                      </Table.Td>
                      <Table.Td>{FIELD_DOC[k]}</Table.Td>
                    </Table.Tr>
                  ))}
                </Table.Tbody>
              </Table>
              </Wide>
            </div>
          ))}
        </Section>

        <Section id="sources" title="Data sources and versions">
          <Text mb={8}>
            ClinVar, Ensembl and gnomAD are queried when you build a panel, so a panel is a snapshot of
            those databases at that moment, and every result carries its pull timestamp. The genetic map
            is bundled and versioned below.
          </Text>
          <Wide>
          <Table className="om-table" withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>Source</Table.Th>
                <Table.Th>What it provides</Table.Th>
                <Table.Th>Version</Table.Th>
                <Table.Th>Ref</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              <Table.Tr>
                <Table.Td>ClinVar</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  Variant identity, clinical significance, review status, accession, and the GRCh37/GRCh38
                  mapping used when a build conversion is needed.
                </Table.Td>
                <Table.Td className="om-mono">live (E-utilities)</Table.Td>
                <Table.Td>
                  <Ref id="clinvar" />
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td>Ensembl</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  Gene model and strand, 1000 Genomes frequencies, gene lookup, and the independent
                  position cross-check on the top markers.
                </Table.Td>
                <Table.Td className="om-mono">release {ensembl}</Table.Td>
                <Table.Td>
                  <Ref id="ensembl" /> <Ref id="ensembl_rest" />
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td>gnomAD</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  Allele frequencies, global and per-ancestry, and the region enumeration that produces the
                  candidate pool.
                </Table.Td>
                <Table.Td className="om-mono">{gnomad}</Table.Td>
                <Table.Td>
                  <Ref id="gnomad_v4" />
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td>dbSNP</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  The rsID namespace every marker is reported in, and the identifier the outbound record
                  links resolve against.
                </Table.Td>
                <Table.Td className="om-mono">via rsID</Table.Td>
                <Table.Td>
                  <Ref id="dbsnp" />
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td>1000 Genomes</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  The reference panel whose allele count decides whether population LD with the pathogenic
                  variant is defined at all.
                </Table.Td>
                <Table.Td className="om-mono">via Ensembl</Table.Td>
                <Table.Td>
                  <Ref id="thousand_g" />
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td>Genetic map</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  cM distance, recombination fraction and hotspot flags. Bundled and offline, so a build
                  needs no network for it and stays reproducible.
                </Table.Td>
                <Table.Td className="om-mono" style={{ whiteSpace: 'normal' }}>
                  {mapSource}
                </Table.Td>
                <Table.Td>
                  <Ref id="decode_map" />
                </Table.Td>
              </Table.Tr>
              <Table.Tr>
                <Table.Td>LDlink</Table.Td>
                <Table.Td style={{ whiteSpace: 'normal' }}>
                  Optional labelled LD prior between two common SNPs. Disabled without a token, which
                  removes the annotation and nothing else.
                </Table.Td>
                <Table.Td className="om-mono">{health ? (health.ldlink_enabled ? 'enabled' : 'not configured') : '-'}</Table.Td>
                <Table.Td>
                  <Ref id="ldlink" />
                </Table.Td>
              </Table.Tr>
            </Table.Tbody>
          </Table>
          </Wide>
          <Text size="xs" c="dimmed" mt={8}>
            API {health ? `v${health.version}` : '-'} · genome build {build} · every coordinate is
            stamped with its build, in the interface and in every export.
          </Text>

          {health?.release && (
            <Paper mt={12} p={8} withBorder>
              <Text size="xs" fw={600}>
                {health.release}
              </Text>
              <Text size="xs" c="dimmed" mt={2}>
                {health.release_gloss}
              </Text>
            </Paper>
          )}
          <Text size="xs" c="dimmed" mt={4}>
            Where a chromosome has no map data, cM is approximated at a uniform 1 cM/Mb and marked as an
            approximation everywhere it appears.
          </Text>
        </Section>

        <Section id="conventions" title="Conventions">
          <Text mb={8}>
            <b>Nomenclature.</b> Variants are written per the HGVS recommendations{' '}
            <Ref id="hgvs" />, with the transcript accession included, as in{' '}
            <span className="om-mono">NM_000352.6(ABCC8):c.3989-9G&gt;A</span>. A c. description without a
            transcript is ambiguous and is not accepted as sufficient.
          </Text>
          <Text mb={8}>
            <b>Strand.</b> HGVS c. descriptions are in transcript sense; VCF is in genomic sense, on the
            plus strand. For a minus-strand gene these disagree by complementation, and the app shows both
            with an explicit label rather than assuming they match. ABCC8 is on the minus strand, so the
            transcript-sense change <span className="om-mono">G&gt;A</span> is the genomic change{' '}
            <span className="om-mono">C&gt;T</span>. The complement is computed from the resolved strand,
            never assumed.
          </Text>
          <Text mb={8}>
            <b>Build.</b> Everything is computed on {build}. Builds are never silently mixed. GRCh37 input
            is converted explicitly through ClinVar's own assembly mapping and the conversion is labelled
            on the record; a GRCh37 position, when shown, is display-only and marked as such.
          </Text>
          <Text>
            <b>Sides.</b> The two flanks are named by <b>coordinate</b>, not by strand and not by
            chromosome arm: <b>lower coordinate</b> and <b>higher coordinate</b>, counted in the API
            payload as <span className="om-mono">lower_count</span> and{' '}
            <span className="om-mono">higher_count</span>. Upstream and downstream would flip meaning
            with gene orientation; a signed bp distance and a coordinate-based side name mean the same
            thing on every gene. No side name here implies a chromosome arm: which flank runs toward the
            telomere depends on where the centromere sits relative to the locus, and nothing in this tool
            looks that up.
          </Text>
        </Section>

        <Section id="freetext" title="Free text and the model">
          <Text mb={8}>
            The search box takes free text. What happens to it turns on one thing: whether your text
            contains a variant identifier. The two paths are not the same feature with different
            wording, and the difference is worth knowing before you trust an answer from either.
          </Text>

          <Title order={3} mt={12} mb={4}>You named it: a regex reads it, and no model runs</Title>
          <Text mb={8}>
            If your text contains an rsID, an HGVS expression or a ClinVar accession, a regular
            expression lifts it out verbatim and nothing is sent to a model. The same local reading
            picks up your modifiers: window size, ancestry, MAF floor. The identifier that reaches the
            pipeline is the one you typed, character for character, and the request costs nothing.
            Text naming two identifiers is refused rather than settled by guesswork, since{' '}
            <i>not rs1, use rs2</i> and <i>rs1, or maybe rs2</i> read alike to a regex. This is the
            path nearly every real request takes.
          </Text>

          <Title order={3} mt={12} mb={4}>You described it: a model is asked, and it answers from memory</Title>
          <Text mb={8}>
            If no identifier is present, the text goes to a small model (Haiku, at temperature 0) with
            one question: which variant did this person mean. The model is not extracting the answer
            from your text, because the answer is not in your text. It answers from its own training.
            That is a knowledge claim rather than a reading, and it is the only place in this app where
            something other than a database decides what you are looking at.
          </Text>

          <Title order={3} mt={12} mb={4}>What the model can and cannot touch</Title>
          <Text mb={8}>
            It cannot hand you a coordinate. Not a chromosome, not a position, not an allele, not a
            strand. The typed query it fills, <Code>pb.StructuredQuery</Code>, has no field for any of
            them (see <SecRef id="api" />). That is a property of the
            code, not a rule the model is asked to observe: there is nowhere to put a coordinate, so a
            recalled one has no route into a panel. Whichever path you arrive by, the coordinate on
            your panel came from the same live lookup described in{' '}
            <SecRef id="pipeline" />.
          </Text>
          <Text mb={8}>
            What it does decide is <b>which variant you meant</b>, and it can be wrong. A gene may
            carry hundreds of pathogenic variants, and <i>the ABCC8 splice mutation</i> does not name
            one of them uniquely. So read the result in two parts: the panel around the variant is as
            trustworthy as any other panel this tool builds, and <b>which variant it is a panel about</b>
            {' '}is the part that needs your eyes.
          </Text>

          <Title order={3} mt={12} mb={4}>Safeguards</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              <b>The query type has no coordinate fields.</b> The strongest one, because it is
              structural rather than a promise: there is no field to carry a fabricated position, so
              nothing downstream has to be careful.
            </List.Item>
            <List.Item>
              <b>Coordinate-shaped output is refused.</b> What the model returns must match an
              allow-list of identifier shapes: rsID, HGVS, ClinVar accession. An allow-list, not a
              list of forbidden coordinate formats, because a list of the ways a position can be
              written is never complete, whereas anything not recognisable as an identifier is refused
              by construction. <span className="om-mono">chr11:17397055</span> and{' '}
              <span className="om-mono">11-17397055-C-T</span> do not pass; they raise.
            </List.Item>
            <List.Item>
              <b>The identifier is looked up live</b>, by the same code path a typed one takes, and the
              record that comes back is reconciled against that identifier. A ClinVar hit whose dbSNP
              ids, accession or variant name are not the ones asked for is refused rather than
              reported.
            </List.Item>
            <List.Item>
              <b>The gene is cross-checked against a gene you named.</b> A gene symbol written the way
              HGNC writes it, capitalised, as <Code>ABCC8</Code>, is read out of your text by the same
              local regex, never by the model, and compared against the gene the lookup returns. If you asked about ABCC8 and the identifier turns out to be
              a variant in another gene, that is not a note to read past: Build is withheld until you
              tick a box saying you meant it. The check is worth exactly what you put into it, since
              the gene it compares against is one you typed rather than one the model chose.
            </List.Item>
            <List.Item>
              <b>Nothing is built without you.</b> The resolved variant appears as a card (rsID, gene,
              strand, {build} coordinate, clinical significance) and no panel exists until you click
              Build. When the model supplied the identifier, that card carries a caveat saying so,
              naming the identifier and stating plainly that you did not type it.
            </List.Item>
            <List.Item>
              <b>Every export says so.</b> A model-chosen panel carries the model's id and the text it
              was given into the CSV, XLSX, JSON and PDF, so the caveat survives the download and
              reaches whoever opens the file next. In the CSV it is a column on every row rather than a
              header comment, which is what makes it survive{' '}
              <span className="om-mono">read_csv(comment='#')</span>.
            </List.Item>
            <List.Item>
              <b>The model path is metered before it runs.</b> Per-client and global caps are checked
              ahead of the call, and only text that would actually reach a model is counted against
              them.
            </List.Item>
          </List>

          <Title order={3} mt={12} mb={4}>What this cannot do</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              The model is small, and it is being asked to recall one specific variant out of
              everything it has ever read. It can recall the wrong one, and the risk is worst for
              exactly the genes you are most likely to ask about: the well-studied ones with many
              known pathogenic variants.
            </List.Item>
            <List.Item>
              <b>It cannot tell you when it is unsure.</b> A confidently wrong identifier looks exactly
              like a right one, and there is no score here that separates them. The live lookup is no
              help: it faithfully returns the correct record for the wrong variant. The reconciliation
              above confirms the record matches the identifier, but it cannot know whether the
              identifier was what you meant, and the gene check above only fires when the gene differs
              too.
            </List.Item>
            <List.Item>
              <b>Name no gene and the cross-check has nothing to check against.</b> It compares the
              answer to a gene you typed, so if you typed none it stays silent, and silence here is not
              approval. That is reason enough to name the gene even when the variant is obvious to you.
              Name it as <Code>ABCC8</Code> and not as <Code>abcc8</Code>: capitalisation is what
              distinguishes a symbol from an ordinary word here, so a lowercased symbol reads as prose
              and the check stays silent exactly as if you had named no gene. Symbols that carry
              lowercase by convention, as <Code>C9orf72</Code> and <Code>MT-ND1</Code>, are read
              correctly.
            </List.Item>
            <List.Item>
              The cross-check only catches the wrong <i>gene</i>. Ask for the wrong variant{' '}
              <i>within</i> the gene you named and every check in this section passes: the gene agrees,
              the record is real, the coordinate is live and correct. That is the failure this section
              exists to warn you about, and you are the only one who can catch it.
            </List.Item>
          </List>
          <Text mb={8}>
            The remedy is one line long: <b>if you know the rsID or the HGVS, type it.</b> It is free,
            it is exact, and it skips this entire section.
          </Text>

          <Title order={3} mt={12} mb={4}>Cost, and switching it off</Title>
          <Text>
            An identifier never reaches a model, so ordinary use of this box costs nothing to run. The
            feature is optional and the app is fully functional without it: with no API key configured
            the box asks for an identifier instead of prose, and nothing else changes. That is why an
            absent key degrades the search box rather than breaking the tool. Free text is currently{' '}
            <b>{health ? (health.nl_enabled ? 'enabled' : 'not configured') : 'unknown'}</b> on this
            deployment.
          </Text>
        </Section>

        <Section id="limits" title="Known limitations">
          <Text mb={8}>
            Everything below is a real property of this tool or its data, not a hypothetical. It is
            listed because a plausible wrong answer is worse than an error: an error gets
            investigated, a plausible wrong answer gets used. This list is not exhaustive.
          </Text>

          <Title order={3} mt={12} mb={4}>The genetic map describes a population, not your carrier</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              The autosomal map is <b>sex-averaged</b>, so it describes neither parent. Linkage runs
              through one carrier, of known sex. Female recombination exceeds male genome-wide and the
              ratio varies by region, so cM and θ here can err in <b>either direction</b> for the
              carrier being tested.
            </List.Item>
            <List.Item>
              <b>chrX is the female map</b>, not a sex-averaged one: deCODE reports X distances in
              female meioses. It is correct for a female carrier and overstates recombination for a
              male one, whose X passes to each daughter as a single haplotype outside the
              pseudoautosomal regions.
            </List.Item>
            <List.Item>
              The maps stop short of the telomeres. Beyond either end, cM is <b>extrapolated</b> at a
              uniform 1 cM/Mb and flagged <span className="om-mono">map_approx</span>. chrY and chrM
              have no map at all and use the same fallback throughout.
            </List.Item>
            <List.Item>
              A hotspot verdict on an approximated stretch means <b>not assessed</b>, never
              <i> no hotspot</i>. The uniform fallback can never reach the hotspot threshold.
            </List.Item>
          </List>

          <Title order={3} mt={12} mb={4}>Frequencies are priors, and some are thin</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              2pq is the chance a <i>random</i> member of a population is heterozygous. Your carrier
              either is or is not. A 2pq of 0.5 marker is a coin flip, not a guarantee, and the only
              way to know is to genotype them.
            </List.Item>
            <List.Item>
              gnomAD population sizes differ by orders of magnitude. A per-population floor of AN 200
              (100 people) applies, but at that floor a four-decimal frequency rests on very few
              individuals. Every export carries the <span className="om-mono">an_*</span> columns so
              you can see which figures to lean on.
            </List.Item>
            <List.Item>
              With no ancestry selected, ranking uses the <b>global</b> 2pq. That is a real quantity
              for the whole cohort and a poor one for any particular family: if you know the ancestry,
              select it.
            </List.Item>
          </List>

          <Title order={3} mt={12} mb={4}>Resolution and data</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              ClinVar's search is a <b>relevance ranking</b>, not a lookup. Every candidate record is
              reconciled against the query before use, and unreconcilable input is refused rather than
              approximated, but a variant whose record ClinVar does not surface will not resolve. The
              rsID or the VCV accession always resolves directly.
            </List.Item>
            <List.Item>
              An rsID can name <b>more than one allele</b> at a position (rs334 covers both HbS and
              HbC). Records describing combinations of variants (haplotypes, compound heterozygotes)
              are refused: this tool anchors on one point variant.
            </List.Item>
            <List.Item>
              Structural and copy-number variants are <b>not supported</b>. They have no single
              reference and alternate allele to anchor on.
            </List.Item>
            <List.Item>
              The 1000 Genomes allele count is joined by notation, and Ensembl writes indels
              differently from ClinVar and gnomAD. The join can fail; when it does the count reads
              <i> unavailable</i>, never zero.
            </List.Item>
            <List.Item>
              Responses are cached. Every panel reports both when it was built and how old its data
              is: check <span className="om-mono">queried_utc</span> against
              <span className="om-mono"> built_utc</span> before treating a rebuild as a re-check.
              Classifications are revised and frequencies are re-called upstream.
            </List.Item>
            <List.Item>
              The independent position cross-check runs only on the nearest few shortlisted markers,
              because it costs a request each. A blank cross-check means <b>not checked</b>, not
              agreement.
            </List.Item>
          </List>

          <Title order={3} mt={12} mb={4}>Method</Title>
          <List spacing={4} mb={8}>
            <List.Item>
              Markers are ranked by heterozygosity and proximity only. Genetic distance and hotspot
              status are <b>reported but not ranked on</b>, so a well-placed marker can still sit
              behind a hotspot. The coverage card says how many do.
            </List.Item>
            <List.Item>
              Linkage disequilibrium is never a ranking key, by design. LD is a population statistic
              and cannot substitute for phasing the family.
            </List.Item>
            <List.Item>
              No side name implies a chromosome arm. Nothing here knows where the centromere is.
            </List.Item>
            <List.Item>
              Loading a CSV with <span className="om-mono">pandas.read_csv(comment='#')</span> strips
              every caveat in the header, including the disclaimer. The JSON and XLSX exports carry
              the same fields as data.
            </List.Item>
          </List>
        </Section>

        <Section id="example" title="Worked example">
          <Text mb={8}>The reference case, end to end.</Text>
          <Wide>
          <Table className="om-table" withTableBorder>
            <Table.Tbody>
              <ExRow k="Query" v="NM_000352.6(ABCC8):c.3989-9G>A" mono />
              <ExRow k="Resolves to" v="rs151344623" mono />
              <ExRow k="Gene" v="ABCC8, chr11, minus strand, 11p15.1" />
              <ExRow k="Coordinate" v="GRCh38 chr11:17,397,055" mono />
              <ExRow k="VCF" v="11-17397055-C-T" mono />
              <ExRow k="Transcript sense" v="G>A (genomic C>T; minus-strand gene)" mono />
              <ExRow k="ClinVar" v="VCV000009088. Pathogenic, multiple submitters, no conflicts" />
              <ExRow k="gnomAD genome AF" v="≈1.8e-4" mono />
              <ExRow k="1000 Genomes AC" v="1 (LD with this allele is undefined)" mono />
              <ExRow k="Candidate pool" v="~~1,200 SNPs with MAF ≥ 0.05 in ±250 kb" />
              <ExRow k="Nearest strong candidate" v="rs757110 at chr11:17,396,930 (−125 bp)" mono />
              <ExRow k="Distant sentinels" v="rs615358 (~−208 kb), rs1476699 (~+209 kb)" mono />
              <ExRow k="Coverage" v="both sides covered; no flags" />
            </Table.Tbody>
          </Table>
          </Wide>
          <Text mb={8}>
            An allele count of 1 puts this variant in the undefined-LD case of{' '}
            <SecRef id="ld" />, so ranking is on 2pq and proximity.
            rs757110 sits 125 bp from the variant, so recombination between them is negligible. It still
            has to be heterozygous in the carrier, and it still needs a partner on its own side and
            coverage on the other.
          </Text>
          <Text>
            For context only: ABCC8 variants cause congenital hyperinsulinism <Ref id="abcc8_chi" />. That
            is background, not an interpretation of this variant in any individual.
          </Text>
        </Section>

        <Section id="api" title="API">
          <Wide>
          <Table className="om-table" withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>Endpoint</Table.Th>
                <Table.Th>Purpose</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              {ENDPOINTS.map(([ep, what]) => (
                <Table.Tr key={ep}>
                  <Table.Td className="om-mono" style={{ verticalAlign: 'top' }}>
                    {ep}
                  </Table.Td>
                  <Table.Td style={{ whiteSpace: 'normal' }}>{what}</Table.Td>
                </Table.Tr>
              ))}
            </Table.Tbody>
          </Table>
          </Wide>
          <Text mb={8}>
            The engine, <Code>panelbuilder</Code>, is importable on its own, with no web framework:
          </Text>
          <Paper p={8} mb={10} style={{ background: 'var(--om-zebra)' }}>
            <pre className="om-mono" style={{ margin: 0, whiteSpace: 'pre-wrap' }}>
              {`import panelbuilder as pb
r = pb.build("NM_000352.6(ABCC8):c.3989-9G>A", window=250_000, ancestry="NFE")
r.rarity.population_LD_usable   # False -> per-family phasing required
r.recommended                   # balanced, both-sided candidate panel`}
            </pre>
          </Paper>
          <Text>
            The free-text endpoint is typed as <Code>pb.StructuredQuery</Code>, which has no chromosome,
            position, ref, alt or strand field: it carries the variant as an identifier that is looked up
            live, exactly as the manual path does. Other consumers (a CLI, a bot, an MCP server) can be
            built against the same contract without touching the genetics code.
          </Text>
        </Section>

        <Section id="references" title="References">
          <ol style={{ margin: 0, paddingLeft: 22 }}>
            {ORDER.map((id) => {
              const c = CITATIONS[id]
              return (
                <li key={id} id={`ref-${id}`} style={{ marginBottom: 5, fontSize: 12, lineHeight: 1.45 }}>
                  {formatCitation(c)}{' '}
                  {c.doi ? (
                    <Anchor href={`https://doi.org/${c.doi}`} target="_blank" rel="noopener noreferrer" size="xs" className="om-mono">
                      doi:{c.doi}
                    </Anchor>
                  ) : (
                    <Text span size="xs" c="dimmed">
                      {c.note}
                    </Text>
                  )}
                </li>
              )
            })}
          </ol>
        </Section>

        <Group justify="space-between" align="center" mt="xl" pt={12}
               style={{ borderTop: '1px solid var(--om-border)' }} wrap="wrap" gap={8}>
          <Text size="xs" c="dimmed">
            Questions, corrections, or a bug in a panel?
          </Text>
          <Button component="a" href={`mailto:${CONTACT}?subject=OriginMarker`}
                  variant="default" size="xs">
            Contact
          </Button>
        </Group>

        <Text size="xs" c="dimmed" mt="lg" pt={8} style={{ borderTop: '1px solid var(--om-border)' }}>
          {health?.disclaimer ?? FALLBACK_DISCLAIMER}
        </Text>
      </article>
    </div>
  )
}

/**
 * One numbered section. The number comes from SECTIONS order, exactly as a citation's comes
 * from CITATIONS order: there is one ordering, not two, so the nav and the heading cannot
 * disagree and inserting a section renumbers everything below it for free.
 */
function Section({ id, title, children }: { id: string; title: string; children: ReactNode }) {
  return (
    <section id={id} style={{ scrollMarginTop: 12, marginBottom: 22 }}>
      <Title order={2} mb={6} pb={3} style={{ borderBottom: '1px solid var(--om-border)' }}>
        {sectionNo(id)} · {title}
      </Title>
      {children}
    </section>
  )
}

/** A cross-reference in prose. Written as an id, rendered as whatever number that section
 *  currently holds, so a renumber carries it rather than stranding it. */
const SecRef = ({ id }: { id: string }) => (
  <Anchor href={docHref(id)}>section {sectionNo(id)}</Anchor>
)

/** Scroll box for tables, which do not fit the prose column's measure. */
const Wide = ({ children }: { children: ReactNode }) => (
  <div style={{ overflowX: 'auto', marginBottom: 10 }}>{children}</div>
)

const ExRow = ({ k, v, mono }: { k: string; v: string; mono?: boolean }) => (
  <Table.Tr>
    <Table.Td style={{ width: 190, color: 'var(--om-text-dim)' }}>{k}</Table.Td>
    <Table.Td className={mono ? 'om-mono' : undefined} style={{ whiteSpace: 'normal' }}>
      {v}
    </Table.Td>
  </Table.Tr>
)

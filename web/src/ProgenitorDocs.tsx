import { Alert, Anchor, Code, List, Paper, Table, Text } from '@mantine/core'
import { CITATIONS, formatCitation } from './citations'
import type { Health } from './api'
import { DocsShell, docsHelpers, type DocSection } from './DocsShell'
import { LIMITS } from './InferredLimits'
import {
  MIN_ASCERTAINMENT, MIN_PRODUCTS, SAME_PARENT_MAX, DIFFERENT_PARENT_MIN, pAgree,
} from './inferredReference'
import { HET_RATE_DIPLOID } from './ProductPanel'
import { ABSENCE_MARGIN, CALL_RATE_FLOOR, HET_BAND_DIPLOID, pct } from './parentage'
import {
  PROGENITOR_CITATION, PROGENITOR_EXAMPLES, PROGENITOR_MARKERS, PROGENITOR_STRIDE,
} from './progenitorExamples'

/**
 * Progenitor's documentation, on the same furniture as the Syngamy and panel docs.
 *
 * A reader arrives here holding files and one question: can a reference built without the
 * parent's array be trusted. So the page is organised around the error rather than the feature.
 * Every threshold below is quoted from the constant the code actually uses, not retyped, so the
 * page cannot drift from the tool the way a hand-maintained number would.
 */

const ORDER = Object.keys(CITATIONS)
const PREFIX = 'progenitor-docs'

const SECTIONS: DocSection[] = [
  { id: 'what', label: 'What this answers' },
  { id: 'procedure', label: 'Methods text, for a paper' },
  { id: 'method', label: 'The method, and its one error' },
  { id: 'products', label: 'What counts as a product' },
  { id: 'using', label: 'Using the page' },
  { id: 'formats', label: 'File formats it accepts' },
  { id: 'quality', label: 'Stage 1: quality and ploidy' },
  { id: 'membership', label: 'Stage 2: which products share a parent' },
  { id: 'naming', label: 'Which group is the father\u2019s' },
  { id: 'origin', label: 'Calling parental origin' },
  { id: 'depth', label: 'Stage 3: how deep to build' },
  { id: 'contamination', label: 'Contamination, and the absence it adds' },
  { id: 'verify', label: 'Leave-one-out verification' },
  { id: 'floor', label: 'The five-product floor' },
  { id: 'refusals', label: 'What it refuses to report' },
  { id: 'exports', label: 'Exports and the report' },
  { id: 'examples', label: 'The bundled example data' },
  { id: 'validation', label: 'Validation record' },
  { id: 'parity', label: 'Two implementations, one answer' },
  { id: 'limits', label: 'Scope and limits' },
  { id: 'privacy', label: 'Privacy and data handling' },
  { id: 'references', label: 'References' },
]

const { Section, SecRef } = docsHelpers(PREFIX, SECTIONS)

export function progenitorSectionFromHash(hash: string): string {
  const m = new RegExp(`^#/${PREFIX}/([\\w:.-]+)$`).exec(hash)
  return m ? m[1] : ''
}

function Ref({ id }: { id: string }) {
  const n = ORDER.indexOf(id) + 1
  if (n < 1) return null
  return (
    <Anchor
      href={`#/${PREFIX}/ref-${id}`} className="om-mono" style={{ fontSize: 11 }}
      aria-label={`reference ${n}`}
    >
      [{n}]
    </Anchor>
  )
}

const Wide = ({ children }: { children: React.ReactNode }) => (
  <div style={{ overflowX: 'auto', marginBottom: 10 }}>{children}</div>
)

/** A row of figures, right-aligned and monospaced so a column reads down. */
const NumRow = ({ c }: { c: string[] }) => (
  <Table.Tr>
    {c.map((v, i) => (
      <Table.Td
        key={v + String(i)}
        className={i ? 'om-mono' : undefined}
        style={i ? { textAlign: 'right', fontVariantNumeric: 'tabular-nums' } : undefined}
      >
        {v}
      </Table.Td>
    ))}
  </Table.Tr>
)

export function ProgenitorDocsPage({ health }: { health: Health | null }) {
  return (
    <DocsShell
      prefix={PREFIX}
      sections={SECTIONS}
      health={health}
      title="Progenitor documentation"
      subtitle={
        <>
          Reconstructing a parent&rsquo;s genotype from the haploid cells that parent produced,
          with no array of the parent. Research use only, not a clinical diagnostic.
        </>
      }
      siblings={[
        { label: 'Syngamy documentation', href: '#/syngamy-docs' },
        { label: 'Panel documentation', href: '#/docs' },
      ]}
    >
      {/* --- 1 --------------------------------------------------------------------------- */}
      <Section id="what" title="What this answers">
        <Text size="sm" mb={10}>
          You have several haploid cells from one parent and no array of the parent. Progenitor
          reconstructs that parent&rsquo;s genotype from the cells alone, then reports how much
          of the parent it recovered, how much of the result is wrong, and which of the submitted
          arrays actually belong to that parent.
        </Text>
        <Text size="sm" mb={10}>
          The case it exists for is the ordinary one: the man is not available for an array, or
          consent covers the embryos and not him, or the sample predates the question. Syngamy
          needs his array. This does not.
        </Text>
        <Text size="sm" mb={10}>
          It answers three things and refuses several others. It answers whether the submitted
          arrays come from one parent or more than one (<SecRef id="membership" />), what the
          reconstructed reference contains and how contaminated it is
          (<SecRef id="depth" />, <SecRef id="contamination" />), and whether a given sample
          carries that parent&rsquo;s genome (<SecRef id="verify" />). What it will not report,
          and why each refusal was measured rather than assumed, is <SecRef id="refusals" />.
        </Text>
        <Alert color="orange" p="xs">
          <Text size="xs">
            A reconstructed reference is not a measured array and the two are not interchangeable.
            Rates from a run against a real parental array and rates from a run against a
            reconstruction have different denominators, so they must not be compared or pooled.
            Every export repeats this, on every page, for that reason.
          </Text>
        </Alert>
      </Section>

      {/* --- 2 --------------------------------------------------------------------------- */}
      <Section id="procedure" title="Methods text, for a paper">
        <Text size="sm" mb={10}>
          The whole procedure in prose, written to be quoted or adapted in a methods section
          rather than read as documentation. It states what the software does and what it
          refuses to do; it does not state what your samples are, so the bracketed parts are
          yours to fill in. Every threshold named here is the shipped value and is derived in
          the sections below.
        </Text>
        <Paper withBorder p="sm" style={{ background: 'var(--om-zebra)' }}>
          <Text size="sm" style={{ lineHeight: 1.7 }}>
            Parental origin was determined with OriginMarker (Progenitor,
            v{health?.version ?? 'x.y.z'}), which reconstructs a parental genotype from
            haploid meiotic products and requires no array of either parent. SNP array genotypes
            for [N] [pronuclei / polar bodies / single sperm] were supplied as [platform] exports
            and processed entirely client-side.
          </Text>
          <Text size="sm" mt={8} style={{ lineHeight: 1.7 }}>
            Arrays were first gated on call rate and ploidy: any array below a{' '}
            {pct(CALL_RATE_FLOOR, 0)} call rate, or whose B-allele frequency band and genotype
            heterozygosity were inconsistent with a single haploid genome, was excluded from
            reconstruction. Remaining products were compared pairwise at markers where both were
            called homozygous, and the rate of opposite homozygous calls was used to partition
            them into groups sharing one parent. Grouping required every pair within a group to
            fall below {pct(SAME_PARENT_MAX, 1)}; a chain of pairwise links was not accepted,
            and groups were found as exact maximum cliques rather than by greedy assignment.
          </Text>
          <Text size="sm" mt={8} style={{ lineHeight: 1.7 }}>
            Which group derived from the father was established from chromosome Y, on the basis
            that a maternal cell cannot carry one. A product was called Y-bearing only where its
            Y call rate and its median Y intensity, each relative to its own autosomes, both
            indicated a whole chromosome; either signal alone is insufficient. Where no group
            carried a Y, or more than one did, the assignment of groups to parents was withheld
            and the two sides are reported without parental labels.
          </Text>
          <Text size="sm" mt={8} style={{ lineHeight: 1.7 }}>
            The parental genotype was reconstructed from the largest group meeting a minimum of{' '}
            {MIN_PRODUCTS} products. At each marker, the parental allele was called where at
            least m products were called and all agreed; m was chosen as the deepest threshold
            retaining at least {pct(MIN_ASCERTAINMENT, 0)} of the parent&rsquo;s genome-wide
            heterozygosity. Parental heterozygosity was estimated from the rate at which products
            disagreed. Sites at which the parent is heterozygous but all products drew the same
            allele are asserted homozygous by construction; this residual is reported as
            contamination and propagates to a known floor of spurious absence.
          </Text>
          <Text size="sm" mt={8} style={{ lineHeight: 1.7 }}>
            Every array was then scored against the reconstructed genotype as the rate at which
            the parent&rsquo;s obligate allele was absent where that parent is homozygous, and
            compared with a per-array noise ceiling computed from that array&rsquo;s own no-call
            rate and heterozygous fraction. An array at or below its ceiling was called as
            carrying that parent&rsquo;s genome, at or above {ABSENCE_MARGIN.toFixed(0)}x its
            ceiling as lacking it, and between the two as uncalled. Products contributing to the
            reconstruction were scored against a reference rebuilt without them. No confidence
            percentage is reported.
          </Text>
          <Text size="sm" mt={8} style={{ lineHeight: 1.7 }}>
            Arrays excluded by the ploidy and call-rate gates were scored and reported, since
            exclusion from the reconstruction is not a reason to withhold a call about the array
            itself. The tool does not report [list what you are not claiming, for example
            androgenetic versus biparental status, results on chromosome X, or discrimination of
            the reconstructed parent from a close relative], all of which it withholds against an
            inferred reference.
          </Text>
        </Paper>
        <Text size="xs" c="dimmed" mt={10} style={{ lineHeight: 1.55 }}>
          Cite the software and the version, which is on every export and on the report. The
          reconstruction is deterministic: the same product files at the same version reproduce
          the same genotype byte for byte, so the input files and their checksums, which the run
          manifest records, are what makes it reproducible.
        </Text>
      </Section>

      <Section id="method" title="The method, and its one error">
        <Text size="sm" mb={10}>
          A haploid cell is one meiotic product: it carries one of the parent&rsquo;s two
          chromosome sets, recombined. So at a marker where the parent is homozygous, every
          product carries that allele, and the parent&rsquo;s genotype there is recovered exactly
          by observing any product at all. That is most of the genome and it is free.
        </Text>
        <Text size="sm" mb={10}>
          At a marker where the parent is heterozygous, each product carries one allele or the
          other. The site is recognised as heterozygous only by observing <b>both</b> alleles
          across the products. Observing the same allele in all of them is indistinguishable, from
          the inside, from a parent who is genuinely homozygous.
        </Text>
        <Text size="sm" mb={10}>
          That is the entire error. A heterozygous site where every product happened to agree
          enters the reference as homozygous, and a true offspring then reads as missing the
          parental allele at half of those. The reconstruction itself is deterministic, so re-running
          it on the same files reproduces it byte for byte; what is uncertain is which sites that
          error took, and how far the retained markers have drifted from the genome.
        </Text>
        <Text size="sm" mb={10}>
          Two quantities control that error and conflating them is the mistake this tool was
          rebuilt to avoid:
        </Text>
        <List size="sm" spacing={4} mb={10}>
          <List.Item>
            <Code>m</Code>, how many products were <i>called</i> at a marker. This sets how often a
            heterozygous site survives unrecognised, at roughly 2<sup>1&minus;m</sup>.
          </List.Item>
          <List.Item>
            <Code>n</Code>, how many products <i>exist</i>. This sets how many markers reach any
            given m, and so how much of the genome the reference covers.
          </List.Item>
        </List>
        <Text size="sm" mb={10}>
          An earlier version fixed the agreement threshold at <Code>m = n &minus; 1</Code>. That is
          right at five products and wrong at eight, where it demanded seven agreeing calls and
          left the retained marker set holding roughly 62% of the genome&rsquo;s heterozygosity.
          The
          threshold is now chosen by measuring the marker set instead. See <SecRef id="depth" />.
        </Text>
        <Text size="sm" mb={10}>
          The agreement probability is not the textbook one. Products are not independent at a
          heterozygous site, because a probe with allele-specific dropout calls the same
          homozygote in every product, and those are exactly the probes with low m. Measured
          against a real parental array over 133,631 heterozygous autosomal markers using six
          products of that parent:
        </Text>
        <Wide>
          <Table striped withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>markers called</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>P(all agree | het)</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>2^(1&minus;m)</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>excess</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>markers</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              <NumRow c={['m = 2', '0.61958', '0.50000', '1.24x', '4,332']} />
              <NumRow c={['m = 3', '0.35304', '0.25000', '1.41x', '10,588']} />
              <NumRow c={['m = 4', '0.17709', '0.12500', '1.42x', '23,830']} />
              <NumRow c={['m = 5', '0.08023', '0.06250', '1.28x', '42,153']} />
              <NumRow c={['m = 6', '0.03132', '0.03125', '1.00x', '49,196']} />
            </Table.Tbody>
          </Table>
        </Wide>
        <Text size="sm" mb={10}>
          The excess vanishes by six, so a reference deep enough to be usable needs no correction
          at all. Below six the measured factor is applied: the tool uses{' '}
          <Code>P(agree | het)</Code> of {pAgree(4).toFixed(4)} at m = 4, not 0.125.
        </Text>
      </Section>

      {/* --- 3 --------------------------------------------------------------------------- */}
      <Section id="products" title="What counts as a product">
        <Text size="sm" mb={10}>
          Anything that is one haploid meiotic product of the parent, arrayed on its own:
        </Text>
        <List size="sm" spacing={4} mb={10}>
          <List.Item>
            a paternal pronucleus removed from a zygote before syngamy, which is what the bundled
            examples are <Ref id="zuccaro_2020" />
          </List.Item>
          <List.Item>a single sperm cell</List.Item>
          <List.Item>
            a polar body, for the maternal side. The method is symmetric and nothing in it is
            specific to the father, though the validation record is paternal
            (<SecRef id="validation" />)
          </List.Item>
        </List>
        <Text size="sm" mb={10}>
          What does not count: an embryo, a blastocyst biopsy, a cell line, or any diploid sample.
          A diploid array is not one meiotic product and admitting it as one is actively harmful
          rather than merely useless, because a diploid compared against a haploid has a different
          expected concordance entirely and bridges unrelated parents into one group. The tool
          refuses them at stage 1 rather than trusting the label on the file
          (<SecRef id="quality" />).
        </Text>
        <Text size="sm">
          Sibling products of one parent are correlated only through that parent. Two products of
          the same zygote are not two products: they carry the same meiotic draw, and the agreement
          arithmetic above would count them twice.
        </Text>
      </Section>

      {/* --- 4 --------------------------------------------------------------------------- */}
      <Section id="using" title="Using the page">
        <List size="sm" spacing={6} type="ordered" mb={10}>
          <List.Item>
            <b>Add files</b>, or drop them on the page. One file per haploid cell. Nothing is
            uploaded: every byte is read in the tab. <b>Examples</b> loads five public products
            instead, which is the fastest way to see the whole flow (<SecRef id="examples" />).
          </List.Item>
          <List.Item>
            Each file is streamed and triaged as it arrives, and appears on the page as it
            finishes rather than when the whole batch does. The <b>run log</b> records what was
            read and what was set aside, in order, and downloads as a text file.
          </List.Item>
          <List.Item>
            <b>Run</b> does the rest in one pass: every pair of usable products is compared and
            the products split into groups sharing a parent, the father&rsquo;s group is named
            from the Y its products carry (<SecRef id="naming" />), the largest group that clears
            the floor is reconstructed, and every array that went in is called against it
            (<SecRef id="origin" />). Members of the reconstructed group are verified against a
            reference built without them.
          </List.Item>
          <List.Item>
            Export. Six artefacts, all written in the browser (<SecRef id="exports" />),
            including the reconstructed genotype as an array file. Under that one is{' '}
            <b>Download and open in Syngamy</b>, which saves the array and opens Syngamy in a new
            tab with it already loaded as the donor, so the two halves of the run do not need a
            file picker between them.
          </List.Item>
        </List>
        <Text size="sm" mb={10}>
          Each stage carries its own headline outcome and a collapsed <b>Show the evidence</b>{' '}
          panel holding the numbers behind it. The headline is the answer; the panel is the
          working, and the first stage is the parental-origin call itself. Adding or removing a
          file invalidates the result rather than leaving a stale one on screen, and Run has to
          be pressed again.
        </Text>
        <Text size="sm">
          The order of the stages is a safety property rather than presentation. Running
          concordance before the call-rate gate splits one parent into several, because a
          degraded array reads high against everyone: a true product at a 53.9% call rate read
          12.7% to 14.4% against its own co-products, which the same-parent threshold calls a
          different man. Running it before the ploidy gate merges unrelated parents. Nothing is
          built until membership is settled, and membership is not attempted until quality is.
        </Text>
      </Section>

      {/* --- 5 --------------------------------------------------------------------------- */}
      <Section id="formats" title="File formats it accepts">
        <Text size="sm" mb={10}>
          The same exports the rest of the tool reads, because the intake is shared code rather
          than a second parser: Axiom <Code>.probes</Code>, or a genotype table as CSV or TSV,
          gzipped or not. A header row is required and the columns are found by name, so column
          order does not matter.
        </Text>
        <Text size="sm" mb={10}>
          What is used: the probe identifier, the chromosome, and the genotype call. B-allele
          frequency is used for the ploidy gate when present. Copy number and log ratio are read
          but contribute nothing here.
        </Text>
        <Text size="sm">
          Markers are matched between products by probe identifier, not by coordinate, so no
          assembly conversion happens and none is needed. Files from different array versions
          intersect on whatever probe identifiers they share, and the marker count in the result
          is that intersection at the chosen depth. Only autosomal markers become evidence:
          chromosome X and Y are read and excluded (<SecRef id="refusals" />).
        </Text>
      </Section>

      {/* --- 6 --------------------------------------------------------------------------- */}
      <Section id="quality" title="Stage 1: quality and ploidy">
        <Text size="sm" mb={10}>
          Two gates, both of which must pass before an array is treated as a product.
        </Text>
        <Text size="sm" mb={10}>
          <b>Call rate</b> must be at least {pct(CALL_RATE_FLOOR, 0)}. Below that an array reads
          high against everyone and fragments a single parent into several groups.
        </Text>
        <Text size="sm" mb={10}>
          <b>Ploidy</b> is judged on two signals, and failing either withholds the array. The
          B-allele frequency band must be at or under {pct(HET_BAND_DIPLOID, 0)}, and the
          heterozygous-call rate under {pct(HET_RATE_DIPLOID, 0)}. Neither alone is sufficient,
          which was measured rather than reasoned: a known unrelated diploid adult measures a
          7.17% BAF band and passes the band gate on its own. Admitted as a product it read 10.9%
          against a genuine product, fragmented one parent&rsquo;s group in two, and produced
          three groups where there are two. Genotype heterozygosity catches it, because a haploid
          genome cannot be heterozygous at all and every such call is error: twelve confirmed
          products ran 5.3% to 13.7% and that adult reads 15.0%.
        </Text>
        <Text size="sm">
          The margin there is 1.3 percentage points, which is thin, so anything between the
          haploid and diploid ranges is set aside as borderline rather than argued about. A
          product set aside is still listed with its reason, and still appears in the per-sample
          export, so the table accounts for every file submitted.
        </Text>
      </Section>

      {/* --- 7 --------------------------------------------------------------------------- */}
      <Section id="membership" title="Stage 2: which products share a parent">
        <Text size="sm" mb={10}>
          Products of one parent are not identical. Each carries a different meiotic draw, so at a
          marker where the parent is heterozygous two products disagree half the time. The rate of
          <i> opposite homozygotes</i> between two products is therefore a direct function of
          relatedness and needs no reference genome to compute. Membership is settled this way,
          before anything is built, because a reference built from two people is worse than no
          reference at all.
        </Text>
        <Text size="sm" mb={10}>
          Measured across two experiments: 4.68% to 9.70% between products of one father over 46
          pairs, and 9.88% to 16.10% between products of different fathers over 45 pairs.
        </Text>
        <Wide>
          <Table striped withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>reads</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>opposite homozygotes</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              <NumRow c={['same parent', `under ${pct(SAME_PARENT_MAX, 1)}`]} />
              <NumRow c={['ambiguous', `${pct(SAME_PARENT_MAX, 1)} to ${pct(DIFFERENT_PARENT_MIN, 1)}`]} />
              <NumRow c={['different parents', `${pct(DIFFERENT_PARENT_MIN, 1)} and above`]} />
            </Table.Tbody>
          </Table>
        </Wide>
        <Alert color="orange" p="xs" mb={10}>
          <Text size="xs">
            The separation is real but the margin is 0.18 of a percentage point, not the empty
            void a first look suggests. One genuine cross-father pair sits at 9.88%, under the
            same-parent cut. A per-pair label is therefore <b>not</b> evidence of group
            membership.
          </Text>
        </Alert>
        <Text size="sm" mb={10}>
          So grouping requires <b>every</b> pair inside a group to read same parent, computed as
          an exact maximum clique rather than by chaining links. Single linkage merges two men
          through that one 9.88% edge; requiring all pairs rejects it, because the same product
          reads 12.1% to 12.4% against the rest of that group. One misread pair in 45 becomes zero
          group-level errors. A greedy all-pairs pass is also not enough and was replaced after it
          returned the wrong partition on a synthetic case that the real products, all of one
          father, could not have caught.
        </Text>
        <Text size="sm">
          When more than one group comes back, the page says so and refuses to build until one is
          chosen. This is not hypothetical: one validation series of fourteen products turned out
          to hold two sperm donors, in groups of nine and five, which the concordance found
          without any reference and without being told to look.
        </Text>
      </Section>

      {/* --- 8 --------------------------------------------------------------------------- */}
      <Section id="naming" title="Which group is the father&rsquo;s">
        <Text size="sm" mb={10}>
          Splitting products into groups that share a parent says nothing about WHICH parent. A
          group of siblings looks the same whether the parent they share is the mother or the
          father, and an experiment with no array of either has nothing to compare against. One
          signal settles it: a maternal cell cannot carry a Y chromosome, so a group with one
          Y-bearing product is the father&rsquo;s.
        </Text>
        <Text size="sm" mb={10}>
          Two measurements are required and neither is sufficient on its own. This is not
          caution; each rule alone inverts a real experiment in the validation set.
        </Text>
        <Text size="sm" mb={10}>
          Call rate alone is wrong on an array that genotypes 86.2% of its Y probes while its Y
          intensity sits a full log2 below its own autosomes, exactly where arrays with no Y sit.
          An absent chromosome still produces calls; those are noise on nothing. Taking that
          array at its word names a maternal group paternal and inverts every call in it.
          Intensity alone is wrong in the other direction, on an array reading &minus;0.10 log2
          while calling not one Y probe: there is nothing there to genotype.
        </Text>
        <Text size="sm" mb={10}>
          Requiring both separates cleanly across 46 arrays. Y-bearing arrays call 93.7% to 97.3%
          of their Y probes at +0.16 to +0.43 log2 against their own autosomes. Every other array
          either calls 0.0% or sits at &minus;0.81 to &minus;1.25. The cuts sit inside both gaps.
        </Text>
        <Text size="sm" mb={10}>
          Naming is withheld rather than guessed in two cases. If no product carries a Y, the
          father is not established: a paternal group of n products is all X-bearing 2
          <sup>&minus;n</sup> of the time, which at five products is 3%. If more than one group
          carries a Y, one sperm donor cannot have produced them, so the input is not what it was
          described as. In both cases the split still holds and the two sides are reported as
          &ldquo;this parent&rdquo; and &ldquo;the other parent&rdquo;. Chromosome Y is read only
          to name the group; no product is selected or dropped by it, which would leave no
          paternal X in the reference by construction. See <SecRef id="limits" />.
        </Text>
      </Section>

      <Section id="origin" title="Calling parental origin">
        <Text size="sm" mb={10}>
          This is what the reconstruction is for. Once one parent&rsquo;s genotype exists, every
          array in the experiment is scored against it, not only the ones it was built from. An
          array carrying that parent&rsquo;s genome came from that parent; one decisively lacking
          it came from the other. The measurement is the same one <SecRef id="verify" /> describes
          and is read the same way: at or below the array&rsquo;s own noise ceiling reads present,
          at or above three times it reads absent, and in between is left uncalled.
        </Text>
        <Text size="sm" mb={10}>
          Which parent gets reconstructed is not chosen. The largest group that clears the
          five-product floor is built, whichever parent it belongs to, because the call does not
          depend on it: an experiment whose paternal products are too few but whose maternal ones
          are not is answered by reconstructing the mother and reading the other side off it.
        </Text>
        <Text size="sm" mb={10}>
          Arrays the gates excluded are called too, and their rows say so. A product has to be a
          clean haploid cell to go INTO a reference, since a diploid compared against a haploid
          has a different expected rate entirely and an array below the call-rate floor reads
          high against everyone. Neither is a reason to withhold an answer about it. A fused
          zygote or a half-failed amplification is the case this gets asked about most.
        </Text>
        <Text size="sm" mb={10}>
          Validation, on an experiment where the answer is known independently: 18 arrays were
          dropped in with the sperm donor&rsquo;s own array held back. Sixteen were callable and
          all sixteen match what that array says, including two the laboratory record has the
          wrong way round. The other two are below the call-rate floor and read unclear against
          both. No inversions. See <SecRef id="validation" />.
        </Text>
      </Section>

      <Section id="depth" title="Stage 3: how deep to build">
        <Text size="sm" mb={10}>
          A marker enters the reference when at least <Code>m</Code> products called it and all of
          them agreed. Raising m lowers contamination and narrows the marker set, and the trade is
          not free in the direction it first appears: the surviving probes are the ones nearly
          every product called, which are enriched for sites where the parent is homozygous
          because the minor allele is rare. Unrelated people carry the common allele there too, so
          past a point the reference loses its grip on genomes that do not belong to the parent.
        </Text>
        <Text size="sm" mb={10}>
          The threshold is chosen by measuring that drift, which needs no ground truth. The
          parent&rsquo;s heterozygosity over the retained markers is recovered from how often the
          products disagree, and the deepest setting still holding{' '}
          {pct(MIN_ASCERTAINMENT, 0)} of the m = 2 baseline is taken. Measured on a father whose
          real array gives 16.66% heterozygosity, with the m = 2 estimate landing at 16.93%
          unbiased against it:
        </Text>
        <Wide>
          <Table striped withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th>threshold</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>markers</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>h retained</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>of true</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>contamination</Table.Th>
                <Table.Th style={{ textAlign: 'right' }}>closest negative</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              <NumRow c={['m ≥ 2', '663,542', '16.93%', '102%', '4.288%', '5.96x']} />
              <NumRow c={['m ≥ 3', '596,687', '16.48%', '99%', '3.402%', '5.65x']} />
              <NumRow c={['m ≥ 4', '457,506', '15.47%', '93%', '2.286%', '5.25x']} />
              <NumRow c={['m ≥ 5', '231,064', '13.74%', '82%', '1.258%', '4.67x']} />
            </Table.Tbody>
          </Table>
        </Wide>
        <Text size="sm">
          Recovery plateaus at m = 4 while everything else keeps degrading, which is where the
          {' '}{pct(MIN_ASCERTAINMENT, 0)} floor sits. The ladder is shown on the page for every
          run, with the chosen row marked and the rows below it labelled too narrow, so the choice
          is visible rather than asserted.
        </Text>
      </Section>

      {/* --- 9 --------------------------------------------------------------------------- */}
      <Section id="contamination" title="Contamination, and the absence it adds">
        <Text size="sm" mb={10}>
          <b>Contamination</b> is the fraction of the reference that is a heterozygous site
          mistaken for a homozygous one. It is computed per marker from that marker&rsquo;s own m
          and then averaged, because a marker seen by ten products is far safer than one seen by
          the minimum, and averaging the exponent rather than the probability would understate the
          tail. It is a posterior, so it carries the prior that most markers are genuinely
          homozygous: at 17% heterozygosity with four products agreeing, only about 3.5% of the
          reference is expected to be wrong even though every agreeing marker is a candidate.
        </Text>
        <Text size="sm" mb={10}>
          <b>Ascertainment</b> is the opposite axis and is often confused with it. Ascertainment
          is how much of the <i>parent&rsquo;s</i> heterozygosity the retained marker set still
          represents; contamination is how much of the <i>reference</i> is wrong. A reference can
          be very clean and very unrepresentative at once, which is exactly the failure at high m.
        </Text>
        <Text size="sm" mb={10}>
          Contamination has a direct, calculable consequence. At a contaminated marker the
          reference says AA while the parent is really AB, so a haploid product carries the other
          allele half the time and reads as an opposite homozygote through no fault of its own.
          That is <b>added absence</b>, equal to half the contamination, and it is reported next to
          the contamination for that reason.
        </Text>
        <Alert color="blue" p="xs" mb={10}>
          <Text size="xs">
            Because it is known rather than guessed, it is admitted into the noise ceiling when a
            haploid product is scored. Leaving it out judges every true product against a ceiling
            roughly its own size, and true products read unclear. A diploid sample carries no such
            term: its added absence needs the minor allele from the other parent as well, and the
            diploid path is calibrated without one.
          </Text>
        </Alert>
        <Text size="sm">
          The estimator has an honest limit, pinned in the self-check rather than papered over.
          Contamination is derived from how often the products <i>disagree</i>. A set that never
          disagrees anywhere is indistinguishable, from the inside, from a parent who is homozygous
          across the whole genome, so the recovered heterozygosity is zero and the reported
          contamination is zero even though such a reference would be nothing but contamination.
          Real product sets disagree at thousands of markers, which is why this is a degenerate
          case and not a live hazard. The tool cannot detect it and does not pretend to.
        </Text>
      </Section>

      {/* --- 10 -------------------------------------------------------------------------- */}
      <Section id="verify" title="Leave-one-out verification">
        <Text size="sm" mb={10}>
          Every product in the chosen group is scored against a reference built <i>without</i> it.
          A product scored against a reference that contains it reads exactly zero absence, which
          is a bias the size of the whole signal, so the exclusion has to be real rather than
          cosmetic.
        </Text>
        <Text size="sm" mb={10}>
          The leave-one-out reference is built from one fewer product and is dirtier for it, so
          each row&rsquo;s ceiling carries that reference&rsquo;s own added absence rather than the
          full build&rsquo;s. When the chosen depth equals the group size the leave-one-out build
          sits one below it, because n &minus; 1 products can never reach a threshold of n.
        </Text>
        <Text size="sm" mb={10}>
          The verdict compares measured absence against that ceiling: at or under it reads present,
          at {ABSENCE_MARGIN} times it or more reads absent, and between the two is left uncalled
          rather than guessed. The band exists because a noisy sample can reach the unrelated range
          against its own parent.
        </Text>
        <Text size="sm">
          This is a check on the reference, not on the products. Every member reading present says
          the reconstruction is self-consistent and usable. A member reading absent says it does
          not belong in the group after all, and the run should be repeated without it.
        </Text>
      </Section>

      {/* --- 11 -------------------------------------------------------------------------- */}
      <Section id="floor" title="The five-product floor">
        <Text size="sm" mb={10}>
          Below {MIN_PRODUCTS} products the tool declines to build at all. This is a hard refusal
          and not a warning, because the failure mode below it is not a weak answer.
        </Text>
        <Text size="sm" mb={10}>
          At three products and fewer, every true offspring tested inverted to a decisive
          <i> wrong</i> answer rather than to a refusal: 24 of 24 across two experiments. Four
          products was not measured separately, so the floor sits at {MIN_PRODUCTS} rather than
          at the last depth known to fail.
        </Text>
        <Text size="sm">
          With few products, agreement at a heterozygous site is common rather than rare, so a large fraction of the reference is
          contaminated, and the added absence of <SecRef id="contamination" /> swamps the real
          signal. The reference stops describing the parent and starts describing whichever alleles
          those few products happened to draw.
        </Text>
      </Section>

      {/* --- 12 -------------------------------------------------------------------------- */}
      <Section id="refusals" title="What it refuses to report">
        <Text size="sm" mb={10}>
          Each of these came out wrong under measurement. None is precautionary, and a refusal
          means the evidence does not reach, which is a different claim from a negative result.
        </Text>
        <List size="sm" spacing={7} mb={10}>
          {LIMITS.map((l) => (
            <List.Item key={l.what}>
              <b>{l.what}</b>, because {l.why}
            </List.Item>
          ))}
        </List>
        <Text size="sm">
          The refusals travel with the results. They are listed on the page, in the report, and in
          the run manifest, so a reader who receives only an export still sees the boundary of the
          evidence.
        </Text>
      </Section>

      {/* --- 13 -------------------------------------------------------------------------- */}
      <Section id="exports" title="Exports and the report">
        <Text size="sm" mb={10}>
          Six artefacts, all written in the browser, all named{' '}
          <Code>progenitor-&lt;artefact&gt;-&lt;report id&gt;</Code> so a directory of them from
          several runs stays sorted and attributable.
        </Text>
        <Wide>
          <Table striped withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th style={{ width: 210 }}>artefact</Table.Th>
                <Table.Th>what it holds</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              {[
                ['inferred (.probes)', 'The reconstructed genotype itself, as an array file in '
                  + 'the same four columns every other export in this family uses. Drop it into '
                  + 'Syngamy as the donor to call arrays this run never saw. It opens with a '
                  + 'banner saying it is not a measured array and carries a machine-readable '
                  + 'mark: Syngamy reads that mark and states it on every artefact of a run made '
                  + 'against it, and Progenitor refuses the file as a product. The button under '
                  + 'this row saves it and hands a copy straight to a new Syngamy tab.'],
                ['report (PDF)', 'Letter, in the same format as the Syngamy report. Every page '
                  + 'states that the reference was inferred, so a page lifted out of context '
                  + 'still says so.'],
                ['concordance (CSV)', 'One row per pair of products with the rate and whether '
                  + 'they grouped together. The shape a plotting or clustering script wants.'],
                ['matrix (CSV)', 'The square concordance matrix as printed, ordered by group. '
                  + 'Loads with index_col=0 for a heatmap or a supplementary table.'],
                ['samples (CSV)', 'One row per array including the excluded ones, each with its '
                  + 'parental-origin call and, where it was excluded, the reason. The table '
                  + 'accounts for every file submitted.'],
                ['run manifest (JSON)', 'The whole run: reference parameters, membership, how the '
                  + 'parent was named and which products carry a Y, every sample with its origin '
                  + 'call, every pair, and every refusal. For a pipeline or a LIMS.'],
              ].map(([k, v]) => (
                <Table.Tr key={k}>
                  <Table.Td className="om-mono" style={{ fontSize: 11 }}>{k}</Table.Td>
                  <Table.Td style={{ whiteSpace: 'normal' }}>{v}</Table.Td>
                </Table.Tr>
              ))}
            </Table.Tbody>
          </Table>
        </Wide>
        <Text size="sm" mb={10}>
          The provenance block repeats in every file. In the CSVs it sits in a leading{' '}
          <Code>#</Code> block, which <Code>pandas.read_csv</Code> and R skip by default, so the
          file still loads in one call with numeric columns typed as numbers.
        </Text>
        <Alert color="orange" p="xs">
          <Text size="xs">
            <b>The reconstructed genotype is exported as an array file</b>, in the same four
            columns every other export in this family uses, so it can be dropped into Syngamy as
            the donor and used to call arrays this run never saw. It was withheld for several
            releases on the reasoning that a file of homozygous calls belonging to an identifiable
            person could be re-imported as though it had been measured. That reasoning was right
            about the hazard and wrong about the gain: reconstructing a parent exists in order to
            have something to call parental origin against. The hazard is handled where it lives.
            Every such file opens with a banner stating it is not a measured array, and carries a
            machine-readable mark saying the same, which the tools that ingest arrays look for.
            The banner lines are comments the header sniffer skips, so the file still loads
            anywhere a real export does.
          </Text>
        </Alert>
      </Section>

      {/* --- 14 -------------------------------------------------------------------------- */}
      <Section id="examples" title="The bundled example data">
        <Text size="sm" mb={10}>
          <b>Examples</b> loads {PROGENITOR_EXAMPLES.length} public paternal pronuclei from the
          series the tool is validated against. They are a positive case on purpose: five products
          of one parent that pass every gate, group as one, reconstruct, and verify. Someone
          seeing the feature for the first time should see it work before they see it refuse.
        </Text>
        <Wide>
          <Table striped withTableBorder>
            <Table.Thead>
              <Table.Tr>
                <Table.Th style={{ width: 120 }}>accession</Table.Th>
                <Table.Th style={{ width: 90 }}>label</Table.Th>
                <Table.Th>what it is</Table.Th>
              </Table.Tr>
            </Table.Thead>
            <Table.Tbody>
              {PROGENITOR_EXAMPLES.map((e) => (
                <Table.Tr key={e.gsm}>
                  <Table.Td className="om-mono" style={{ fontSize: 11 }}>{e.gsm}</Table.Td>
                  <Table.Td className="om-mono" style={{ fontSize: 11 }}>{e.label}</Table.Td>
                  <Table.Td style={{ whiteSpace: 'normal' }}>{e.what}</Table.Td>
                </Table.Tr>
              ))}
            </Table.Tbody>
          </Table>
        </Wide>
        <Text size="sm" mb={10}>
          They are subsets at every {PROGENITOR_STRIDE}th marker, {PROGENITOR_MARKERS.toLocaleString()}{' '}
          markers each, so the download stays under a megabyte. Every quantity the tool reports is
          a rate, so the subset gives the same answer as the full file: reconstructing from these
          gives 15.72% ascertainment and 2.33% contamination against 15.47% and 2.29% on the full
          825,657-marker arrays.
        </Text>
        <Text size="sm" mb={10}>
          The series has only one father, so it cannot demonstrate the split a second father
          produces. That case is real and is documented in <SecRef id="membership" /> rather than
          staged here.
        </Text>
        <Text size="xs" c="dimmed">{PROGENITOR_CITATION} <Ref id="zuccaro_2020" /></Text>
      </Section>

      {/* --- 15 -------------------------------------------------------------------------- */}
      <Section id="validation" title="Validation record">
        <Text size="sm" mb={10}>
          Every threshold on this page came off real arrays. The record below is what was run, not
          a summary of what it is hoped to do.
        </Text>
        <Text size="sm" fw={600} mb={6}>Parental origin, four ways</Text>
        <List size="sm" spacing={6} mb={10}>
          <List.Item>
            <b>Against a held-back parent array, twice.</b> On an experiment of 18 arrays the sperm
            donor&rsquo;s own array was withheld from the run and used only to mark the result: 16
            were callable and all 16 agree, including two the laboratory record has the wrong way
            round. The other two are below the call-rate floor and read unclear against both. On an
            independent public dataset the same test was run from five pronuclei: the
            father&rsquo;s own bulk array reads 0.03x, a SECOND array of the same man 0.09x, and
            his own zygote 0.75x, all present. Two unrelated oocyte donors read 3.06x, absent, and
            2.44x, refused. No wrong answer in either.
          </List.Item>
          <List.Item>
            <b>Specificity, on 26 arrays of other people.</b> A reference built from one
            experiment&rsquo;s father was run against every array of two unrelated experiments.
            Not one read as carrying him. Nineteen are decisive at 1.75x to 10.84x and seven are
            uncalled, five of those being arrays the gates had already excluded.
          </List.Item>
          <List.Item>
            <b>Robustness to depth.</b> Rebuilding from five of the eight products moved no call:
            17 of 18 identical and one degrading to a refusal. No inversion. Contamination rises
            from 0.48% to 2.24% over that drop, which is where the refusal comes from.
          </List.Item>
          <List.Item>
            <b>An independent witness.</b> Y-bearing is read from chromosome Y call rate and
            intensity and shares nothing with the autosomal genotype agreement that produces the
            call. Across the twelve Y-bearing arrays in three experiments the two never disagree:
            seven agree on paternal and five are refused by the autosomal side, all of them
            low-quality arrays.
          </List.Item>
        </List>
        <Text size="sm" mb={10}>
          One thing measured and deliberately NOT acted on: these arrays carry a BestProbeset
          column, and honouring it changes nothing. Filtering to best probes only moves no call of
          18, shifts every ratio by at most 0.07x, raises contamination slightly rather than
          lowering it, and costs 5,795 markers. The flag is parsed and ignored on purpose.
        </Text>
        <List size="sm" spacing={6} mb={10}>
          <List.Item>
            <b>Against a known father.</b> A reference reconstructed from six haploid products,
            with his real array held back and used only to score the result. Ascertainment,
            contamination and the agreement excess in <SecRef id="method" /> and{' '}
            <SecRef id="depth" /> are all measured against that array.
          </List.Item>
          <List.Item>
            <b>A second father, held out entirely.</b> 18 arrays. 16 agree with the answer the
            real array gives, with <b>0 inversions and 0 false relationships</b>. The
            father&rsquo;s own replicate array is recognised at 0.14x by a reference built from
            nothing but his sons&rsquo; pronuclei.
          </List.Item>
          <List.Item>
            <b>Membership found an undeclared second donor.</b> A series of fourteen products
            split into groups of nine and five on concordance alone, with no reference and no
            prompting. That is what forced all-pairs grouping over single linkage
            (<SecRef id="membership" />).
          </List.Item>
          <List.Item>
            <b>The negative controls hold.</b> An unrelated diploid adult reads 3.37x absent
            against the reconstructed reference and 3.35x against the father&rsquo;s real array.
            The reconstruction and the measured array agree on someone who is not the parent.
          </List.Item>
          <List.Item>
            <b>The floor was measured, not chosen.</b> 24 of 24 true products inverted to decisive
            wrong answers at three products and fewer (<SecRef id="floor" />).
          </List.Item>
        </List>
        <Text size="sm">
          Each fix on this page has a check behind it that fails if the defect returns, including
          the ones where the check itself was wrong first: three assertions in the self-check were
          corrected before the code was, because contamination is a posterior, unreachable depths
          give an undefined ascertainment rather than a passing one, and a set that never disagrees
          reports no contamination by construction.
        </Text>
      </Section>

      {/* --- 16 -------------------------------------------------------------------------- */}
      <Section id="parity" title="Two implementations, one answer">
        <Text size="sm" mb={10}>
          The arithmetic exists twice: once in the browser and once in{' '}
          <Code>tools/inferred_reference.py</Code>. That is deliberate, and the storage differs
          because a tab cannot afford what a workstation can. Python holds per-marker observation
          lists; the browser maps probe identifiers to integers once and holds each product as a
          single <Code>Uint8Array</Code> of allele codes, so eight products of 825,657 probes cost
          about 6.6 MB rather than tens of millions of boxed tuples.
        </Text>
        <Text size="sm">
          Separate code paths can drift, and a divergence would be invisible to a user, who only
          ever sees one of them. So a cross-implementation check compares every quantity a user is
          shown, on the public products, at every depth, plus the grouping itself and the synthetic
          case that broke a greedy implementation. It currently reports <b>0 divergences</b>.
        </Text>
      </Section>

      {/* --- 17 -------------------------------------------------------------------------- */}
      <Section id="limits" title="Scope and limits">
        <List size="sm" spacing={6} mb={10}>
          <List.Item>
            Research use only. Candidate markers require validation and per-family phasing in a
            qualified genetics laboratory. This is not a clinical diagnostic.
          </List.Item>
          <List.Item>
            The validation record is paternal. The method is symmetric and polar bodies are
            meiotic products in exactly the same sense, but the numbers above were measured on
            pronuclei and sperm.
          </List.Item>
          <List.Item>
            Autosomes only. See <SecRef id="refusals" /> for why chromosome X cannot be recovered
            here even in principle.
          </List.Item>
          <List.Item>
            Calibrated on common-SNP arrays, where parental heterozygosity runs 15 to 19%. A
            whole-genome variant callset is a different object: most of its sites are rare, they
            cannot show absence in anyone, and they dilute the denominator until related and
            unrelated pairs move together. Restrict such a file to common variants first.
          </List.Item>
          <List.Item>
            Nothing here separates the parent from a close relative. A simulated man sharing half
            the genome landed inside the uncalled band and was never rejected.
          </List.Item>
        </List>
      </Section>

      {/* --- 18 -------------------------------------------------------------------------- */}
      <Section id="privacy" title="Privacy and data handling">
        <Text size="sm" mb={10}>
          Every file is read in the tab. Nothing is uploaded, no request carries genotype data, and
          the reconstruction never leaves the browser. Closing the tab discards it; there is no
          store to clear.
        </Text>
        <Text size="sm">
          The exports are written locally by the same page. Files are identified in them by the
          name you gave them, so a filename carrying a patient or laboratory identifier will carry
          it into the export and the report. Rename before submitting if that matters. See also
          the refusal to export the reconstructed genotype in <SecRef id="exports" />, which exists
          for the same reason.
        </Text>
      </Section>

      {/* --- 19 -------------------------------------------------------------------------- */}
      <Section id="references" title="References">
        <ol style={{ paddingLeft: 20, margin: 0 }}>
          {ORDER.map((id) => {
            const c = CITATIONS[id]
            return (
              <li key={id} id={`ref-${id}`} style={{ marginBottom: 8, scrollMarginTop: 12 }}>
                <Text size="sm">
                  {formatCitation(c)}{' '}
                  {c.url && (
                    <Anchor href={c.url} target="_blank" rel="noreferrer" size="sm">
                      {c.doi ?? 'link'}
                    </Anchor>
                  )}
                </Text>
                {c.note && <Text size="xs" c="dimmed">{c.note}</Text>}
              </li>
            )
          })}
        </ol>
      </Section>
    </DocsShell>
  )
}

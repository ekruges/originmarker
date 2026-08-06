import { type ReactNode } from 'react'
import { Alert, Anchor, Code, List, Table, Text, Title } from '@mantine/core'
import { CITATIONS, formatCitation } from './citations'
import type { Health } from './api'
import { EXAMPLES, EXAMPLE_CITATION, EXAMPLE_MARKERS, EXAMPLE_STRIDE } from './examples'
import { DocsShell, docsHelpers, docsSectionFromHash, type DocSection } from './DocsShell'

/**
 * Syngamy's own documentation, on the same furniture as the panel docs and at its own route.
 *
 * It is a separate page rather than more sections in the other one because it answers a
 * different question. The panel docs are about choosing markers before an experiment; this is
 * about reading arrays after one. A reader arrives here already holding files.
 */

const ORDER = Object.keys(CITATIONS)
const numberOf = (id: string) => ORDER.indexOf(id) + 1

const PREFIX = 'syngamy-docs'
const docHref = (id: string) => `#/${PREFIX}/${id}`

export const synSectionFromHash = (hash: string): string => docsSectionFromHash(PREFIX, hash)

function Ref({ id }: { id: string }) {
  const n = numberOf(id)
  if (n < 1) return null
  return (
    <Anchor
      href={docHref(`ref-${id}`)} className="om-mono" style={{ fontSize: 11 }}
      aria-label={`reference ${n}`}
    >
      [{n}]
    </Anchor>
  )
}

const SECTIONS: DocSection[] = [
  { id: 'what', label: 'What this answers' },
  { id: 'why', label: 'Why the variant site cannot answer it' },
  { id: 'biology', label: 'The biology, from first principles' },
  { id: 'array', label: 'What a SNP array actually measures' },
  { id: 'using', label: 'Using the page' },
  { id: 'formats', label: 'File formats it accepts' },
  { id: 'axes', label: 'The three measurements' },
  { id: 'bound', label: 'The noise bound, derived per sample' },
  { id: 'unrelated', label: 'The constant that had to be abandoned' },
  { id: 'classes', label: 'Reading the verdict' },
  { id: 'band', label: 'The uncalled band' },
  { id: 'chrx', label: 'Sperm type from chrX, not chrY' },
  { id: 'oocyte', label: 'Adding the oocyte donor' },
  { id: 'locus', label: 'The per-locus deletion test' },
  { id: 'confidence', label: 'Confidence, and why no percentage' },
  { id: 'gates', label: 'Quality gates' },
  { id: 'assembly', label: 'Assembly detection, and no liftOver' },
  { id: 'chrom', label: 'Per-chromosome and segmental results' },
  { id: 'segments', label: 'Chromosomal change: segments' },
  { id: 'report', label: 'The report, and how to cite it' },
  { id: 'examples', label: 'The bundled example data' },
  { id: 'backend', label: 'What the backend adds' },
  { id: 'cli', label: 'Command line' },
  { id: 'validation', label: 'Validation record' },
  { id: 'limits', label: 'Scope and limits' },
  { id: 'privacy', label: 'Privacy and data handling' },
  { id: 'references', label: 'References' },
]

const { Section, SecRef } = docsHelpers(PREFIX, SECTIONS)

export function SyngamyDocsPage({ health }: { health: Health | null }) {
  return (
    <DocsShell
      prefix={PREFIX}
      sections={SECTIONS}
      health={health}
      title="Syngamy documentation"
      subtitle={
        <>
          Which parent contributed the genome in this sample, from SNP array genotypes.
          Research use only, not a clinical diagnostic.
        </>
      }
      siblings={[
        { label: 'Progenitor documentation', href: '#/progenitor-docs' },
        { label: 'Panel documentation', href: '#/docs' },
      ]}
    >

        {/* --- 1 ------------------------------------------------------------------------- */}
        <Section id="what" title="What this answers">
          <Text mb={8}>
            You have a sample from an experiment that is already finished: an embryo, a
            blastomere, a cell line. You also have an array from the sperm donor. The question
            is which of the two parental genomes is present in that sample, and on which
            chromosomes. Syngamy answers it from the genotypes alone.
          </Text>
          <Text mb={8}>
            The three outcomes it separates are:
          </Text>
          <List size="sm" spacing={4} mb={10}>
            <List.Item>
              <b>Biparental.</b> Both parents contributed. The ordinary case.
            </List.Item>
            <List.Item>
              <b>Androgenetic.</b> Only the sperm donor&apos;s genome is present. The maternal
              complement is missing, and what is there traces entirely to him.
            </List.Item>
            <List.Item>
              <b>Gynogenetic.</b> No paternal contribution anywhere. Only the maternal genome.
            </List.Item>
          </List>
          <Text mb={8}>
            With an oocyte donor array as well, a fourth outcome becomes reachable:{' '}
            <b>neither parent</b>, meaning the sample belongs to neither of the two people you
            declared. See <SecRef id="oocyte" />.
          </Text>
          <Text>
            It also reports which parts of a genome are missing rather than only whether any of
            it is, so a paternal contribution present on twenty-one chromosomes and absent on one
            reads as that rather than as &quot;present&quot;. See <SecRef id="chrom" />.
          </Text>
        </Section>

        {/* --- 2 ------------------------------------------------------------------------- */}
        <Section id="why" title="Why the variant site cannot answer it">
          <Text mb={8}>
            Suppose an embryo carried a paternal mutation and, after editing, the site reads
            wild-type. Two histories produce that reading and they are not the same result:
          </Text>
          <List size="sm" spacing={4} mb={10} type="ordered">
            <List.Item>
              The mutant allele was <b>repaired</b> to wild-type, and the embryo now carries two
              wild-type copies, one from each parent.
            </List.Item>
            <List.Item>
              The mutant allele was <b>lost</b>. The paternal chromosome carrying it was cut and
              not restored, so what remains is the maternal wild-type copy alone. The site reads
              wild-type because there is only one allele left to read.
            </List.Item>
          </List>
          <Text mb={8}>
            An assay that looks only at the variant site cannot distinguish these, because both
            give the same base at that position. This is not hypothetical: it is the substance of
            a published disagreement about whether human embryos perform inter-homologue repair{' '}
            <Ref id="ma_2017" /> <Ref id="egli_2018" />, and subsequent work found that Cas9
            cleavage in human embryos frequently removes the cut chromosome in segments up to a
            whole arm rather than repairing it <Ref id="zuccaro_2020" />.
          </Text>
          <Text mb={8}>
            Markers away from the cut site separate the two. They were not the editing target, so
            whatever happened at the variant, they still report which parental chromosome is
            physically present. If the paternal chromosome is gone, every paternal-only allele
            along it is gone with it, and that is a signal spread over thousands of markers rather
            than one base.
          </Text>
          <Alert color="blue" p="xs" mb={4}>
            <Text size="xs">
              This is the same logic the panel half of this tool builds panels for, run in the
              opposite direction. There, you choose flanking markers before an experiment. Here,
              you already have genome-wide genotypes and read origin off them directly.
            </Text>
          </Alert>
        </Section>

        {/* --- 3 ------------------------------------------------------------------------- */}
        <Section id="biology" title="The biology, from first principles">
          <Text mb={8}>
            If you already know this, skip to <SecRef id="array" />. Nothing below is unusual, but
            every later section depends on it being precise.
          </Text>
          <Title order={3} mt={14} mb={4}>Two copies, one from each parent</Title>
          <Text mb={8}>
            A human somatic cell carries 22 pairs of autosomes plus two sex chromosomes. At any
            given position you therefore have two alleles. One came from the egg, one from the
            sperm. At a SNP with two variants, conventionally called A and B, your genotype at
            that position is AA, AB or BB. AB is <b>heterozygous</b>: you carry one of each. AA
            and BB are <b>homozygous</b>.
          </Text>
          <Title order={3} mt={14} mb={4}>Syngamy</Title>
          <Text mb={8}>
            Syngamy is the fusion of the two gametic nuclei at fertilisation, the moment the
            single-copy egg genome and the single-copy sperm genome become one two-copy genome.
            The tool is named for it because everything it reports is a statement about whether
            that fusion happened and what survived it.
          </Text>
          <Title order={3} mt={14} mb={4}>Obligate alleles</Title>
          <Text mb={8}>
            The whole method turns on one observation. Take a marker where the sperm donor is
            homozygous, say AA. Every sperm he produces carries an A at that marker, without
            exception, because he has nothing else to give. So if a sample inherited a paternal
            genome, it <b>must</b> carry at least one A there. That A is his <b>obligate
            allele</b>.
          </Text>
          <Text mb={8}>
            Now look at the sample. If it reads BB, it carries no A, and it therefore cannot have
            inherited that stretch of chromosome from him. One such marker is noise. Ten thousand
            of them, on a chromosome where he is homozygous at each, is not.
          </Text>
          <Text mb={8}>
            Markers where he is <b>heterozygous</b> (AB) are useless for this test in the other
            direction: he could have transmitted either allele, so seeing either one in the sample
            proves nothing. They are excluded from the absence measurement entirely, which is why
            the tool reports an &quot;informative&quot; marker count well below the number of
            markers on the chip.
          </Text>
          <Title order={3} mt={14} mb={4}>Androgenote and gynogenote</Title>
          <Text mb={8}>
            An <b>androgenote</b> has only paternal chromosomes. This arises when the egg genome
            fails to participate, and the paternal set is typically then duplicated, giving a
            genome that is two identical paternal copies and therefore homozygous nearly
            everywhere. A <b>gynogenote</b> is the mirror image: only maternal chromosomes.
            Neither is viable, and both are relevant here because both occur in embryo experiments
            and both can be mistaken for something else if you look only at the variant site.
          </Text>
          <Text>
            Homozygosity is why the tool measures zygosity separately from parentage. A genome
            that is homozygous almost everywhere has one allele per position, so it cannot also
            hold a second parent&apos;s complement, whatever the other measurements suggest.
          </Text>
        </Section>

        {/* --- 4 ------------------------------------------------------------------------- */}
        <Section id="array" title="What a SNP array actually measures">
          <Text mb={8}>
            An array does not read bases. It measures fluorescence from two probes per marker, one
            for each allele, and a clustering model turns the pair of intensities into a call. Two
            derived quantities matter here, and knowing what they are makes the rest of this page
            readable.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Tbody>
                <ExRow
                  k="Genotype call"
                  v="AA, AB, BB, or no-call. A no-call means the point fell outside every cluster, so the platform declined to call it. It is not evidence of anything; it is a measurement that failed."
                />
                <ExRow
                  k="B-allele frequency (BAF)"
                  v="The fraction of signal coming from the B probe, between 0 and 1. A true AA sits near 0, BB near 1, and a heterozygote near 0.5. It is continuous, so it survives a genotype call being dropped."
                />
                <ExRow
                  k="Log R ratio"
                  v="Total intensity against a reference, which reports copy number. Read for allele-coding checks and copy-neutral loss of heterozygosity, not used to call parentage."
                />
              </Table.Tbody>
            </Table>
          </Wide>
          <Title order={3} mt={14} mb={4}>Allele dropout, and why it dominates</Title>
          <Text mb={8}>
            An embryo biopsy is a handful of cells, sometimes one. There is not enough DNA to
            hybridise, so it is amplified first, and amplification is uneven: one of the two
            alleles at a position may be amplified far less than the other, or not at all. When
            that happens a genuine AB is measured as AA or BB. This is <b>allele dropout</b>, and
            it is the single largest source of error in this entire workflow{' '}
            <Ref id="ado" /> <Ref id="natesan_2014" />.
          </Text>
          <Text mb={8}>
            Its consequence is precise and worth stating carefully, because the whole calibration
            in <SecRef id="bound" /> follows from it. Dropout can manufacture a false
            paternal-absence call, but only by one route: the sample was genuinely heterozygous,
            the paternal allele dropped, and the call became homozygous for the maternal one. It
            cannot manufacture absence at a marker where the sample is genuinely homozygous,
            because there is no second allele to lose.
          </Text>
          <Alert color="blue" p="xs" mb={4}>
            <Text size="xs">
              BAF survives dropout in a way the genotype call does not. A dropped heterozygote
              still has intensity from both probes, so its BAF still sits mid-band even though its
              genotype was called homozygous or not at all. That is why zygosity is read from BAF,
              and why the BAF band is counted <i>before</i> no-calls are discarded.
            </Text>
          </Alert>
        </Section>

        {/* --- 5 ------------------------------------------------------------------------- */}
        <Section id="using" title="Using the page">
          <Text mb={10}>
            Four steps. Nothing is uploaded at any point, so the only cost of a run is your own
            machine&apos;s time.
          </Text>

          <Title order={3} mt={14} mb={4}>1. Add the files</Title>
          <Text mb={8}>
            Drop array exports on the box, or use <b>Add files</b>. Each becomes a chip. Every
            file is read and profiled the moment it lands, before anything is compared, so the
            chip shows a shape and a quality read immediately: marker count, call rate, sex call,
            and a small per-chromosome heterozygosity track. Nothing in that first pass needs a
            donor, which is why it can run on drop.
          </Text>
          <Text mb={8}>
            <b>Examples</b> loads five public arrays instead. See <SecRef id="examples" />.
          </Text>

          <Title order={3} mt={14} mb={4}>2. Label them</Title>
          <Text mb={8}>
            Each chip carries a three-way control: <Code>sperm</Code>, <Code>oocyte</Code>,{' '}
            <Code>sample</Code>. Exactly one file must be the sperm donor. One may be the oocyte
            donor, which is optional and strongly recommended. Everything else is a sample. The
            first file added is guessed to be the sperm donor and the rest samples; correct it on
            the chip.
          </Text>
          <Alert color="orange" p="xs" mb={10}>
            <Text size="xs">
              Labelling is the one input the tool cannot check for you, and it is the one that
              matters most. A sample labelled as the donor produces a confident, wrong, entirely
              plausible answer. With both parents supplied the tool does check that they are two
              different people (see <SecRef id="oocyte" />), but nothing can tell it that the file
              you called the sperm donor is the right man.
            </Text>
          </Alert>

          <Title order={3} mt={14} mb={4}>3. Run</Title>
          <Text mb={8}>
            The parent files are read once and held as one call per marker. Each sample is then
            streamed against them and discarded, so memory stays flat regardless of how many
            samples you queue. The run log records every stage, and downloads as a text file with
            a build and debug footer.
          </Text>

          <Title order={3} mt={14} mb={4}>4. Read the result</Title>
          <Text mb={8}>
            Each sample gets a card: the class in large type, the measured rates beneath it, and a
            <b> Detail</b> toggle opening the evidence tables, every chromosome, the quality read
            and every gate. <b>Report (PDF)</b> then writes the whole thing out. See{' '}
            <SecRef id="report" />.
          </Text>
        </Section>

        {/* --- 6 ------------------------------------------------------------------------- */}
        <Section id="formats" title="File formats it accepts">
          <Text mb={8}>
            The browser reads delimited text and sniffs the delimiter rather than assuming it.
            Tab and comma are both handled. Four columns are required; the rest are used when
            present and their absence is reported rather than guessed around.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th style={{ width: 170 }}>Column</Table.Th>
                  <Table.Th style={{ width: 90 }}>Needed</Table.Th>
                  <Table.Th>Also accepted as</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <Table.Tr>
                  <Table.Td className="om-mono">probeset_id</Table.Td>
                  <Table.Td>required</Table.Td>
                  <Table.Td>the marker identifier files are joined on</Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td className="om-mono">chr, position</Table.Td>
                  <Table.Td>required</Table.Td>
                  <Table.Td>used for per-chromosome results and assembly detection</Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td className="om-mono">genotype</Table.Td>
                  <Table.Td>required</Table.Td>
                  <Table.Td>
                    AA/AB/BB, Axiom numeric codes 0/1/2/-1, and the usual no-call tokens
                    (NC, --, ---, ?, .). Nucleotide pairs are refused rather than assigned:
                    see below.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td className="om-mono">baf</Table.Td>
                  <Table.Td>strongly wanted</Table.Td>
                  <Table.Td className="om-mono">b_allele_freq, b allele freq, ballelefreq</Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td className="om-mono">log2r</Table.Td>
                  <Table.Td>optional</Table.Td>
                  <Table.Td className="om-mono">
                    normalized_intensity, lrr, log_r_ratio, logrratio
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td className="om-mono">copy_number</Table.Td>
                  <Table.Td>optional</Table.Td>
                  <Table.Td className="om-mono">copynumber, cn</Table.Td>
                </Table.Tr>
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            The alias lists are not decoration. The same UK Biobank Axiom array ships as a
            tab-separated <Code>.CEL.probes.txt</Code> writing <Code>log2R</Code> in one dataset
            and a comma-separated <Code>.CEL.txt</Code> writing{' '}
            <Code>normalized_intensity</Code> in another, from the same lab, with the columns in a
            different order. Matching one spelling would silently drop the entire intensity
            channel on the other.
          </Text>
          <Title order={3} mt={14} mb={4}>Without BAF</Title>
          <Text mb={8}>
            Zygosity falls back to genotype heterozygosity, which is the weaker measure because
            dropout attacks it directly, and the result carries a stated limit saying so. Supply
            BAF where you can.
          </Text>
          <Title order={3} mt={14} mb={4}>Nucleotide genotypes are refused here</Title>
          <Text mb={8}>
            A file writing <Code>AG</Code> rather than <Code>AB</Code> has no inherent A/B
            assignment: which nucleotide is &quot;A&quot; is a per-marker convention. Harmonising
            those must pool across every file in the run, never per file, or two files can be
            assigned opposite conventions at the same marker and every comparison at it becomes
            garbage. This page reads one file at a time and cannot pool, so it refuses such a
            file and says so on the <b>genotype format</b> gate rather than assigning a
            convention it cannot justify. The command line pools and reads them; see{' '}
            <SecRef id="backend" />.
          </Text>
          <Title order={3} mt={14} mb={4}>Beyond the browser</Title>
          <Text>
            The command-line tool additionally reads VCF (GT, AD and DP fields), Illumina
            FinalReport in both wide and long layouts including the bracketed header preamble,
            HapMap raw exports, PLINK binary <Code>.bed/.bim/.fam</Code> and text{' '}
            <Code>.ped/.map</Code>, and structural-variant VCFs. See <SecRef id="cli" />.
          </Text>
        </Section>

        {/* --- 7 ------------------------------------------------------------------------- */}
        <Section id="axes" title="The three measurements">
          <Text mb={10}>
            Every call comes from three numbers, and each is compared against a reference derived
            from the data in front of it rather than a constant fitted somewhere else. That is the
            central design decision of this tool and <SecRef id="unrelated" /> explains what forced
            it.
          </Text>

          <Title order={3} mt={14} mb={4}>1. Paternal absence rate</Title>
          <Text mb={8}>
            Across autosomal markers where the sperm donor is homozygous and the sample is called,
            the fraction at which his obligate allele is missing. This is the axis that carries the
            call: it separates a present paternal genome from an absent one by roughly thirty-fold.
          </Text>
          <Text mb={8}>
            Its reference is the noise bound of <SecRef id="bound" />.
          </Text>

          <Title order={3} mt={14} mb={4}>2. Alleles the donor lacks</Title>
          <Text mb={8}>
            The fraction of called autosomal markers where the sample carries an allele the sperm
            donor does not possess at all. Such an allele came from somewhere else, which means a
            second parent. This is the only axis that can separate a paternal-only genome from a
            biparental one, because his own alleles are present in both cases and absence therefore
            cannot tell them apart.
          </Text>
          <Text mb={8}>
            Its reference is <b>half his heterozygosity</b>, and that is derived rather than
            fitted. A second parent can supply an allele he lacks only at markers where he is
            homozygous, which happens at a rate of{' '}
            <Code>p&sup2;q + q&sup2;p = pq</Code> summed over markers. His own heterozygosity is{' '}
            <Code>2pq</Code> summed over the same markers. So the expected second-parent signal is
            exactly half of a quantity the file already reports.
          </Text>
          <Alert color="orange" p="xs" mb={10}>
            <Text size="xs">
              This axis separates by about 1.6-fold, not thirty. Treat a call that sits near its
              boundary as provisional, and note that the tool says so in the result when it lands
              there.
            </Text>
          </Alert>

          <Title order={3} mt={14} mb={4}>3. Heterozygous BAF band</Title>
          <Text mb={8}>
            The fraction of autosomal markers whose BAF falls between 0.35 and 0.65. A genome with
            two different parental copies is heterozygous at many positions and reads 15 to 16% in
            that band. A uniparental homozygous genome reads 1.3 to 3.4%. The threshold is 8%,
            placed between the two measured populations rather than at a round number.
          </Text>
          <Text>
            Dropout does not empty this band, which is what makes it trustworthy on poor material:
            a biparental embryo at 33% dropout still reads 22.2%.
          </Text>
        </Section>

        {/* --- 8 ------------------------------------------------------------------------- */}
        <Section id="bound" title="The noise bound, derived per sample">
          <Text mb={8}>
            An absence rate means nothing on its own. Against what? Every sample has its own noise,
            and a sample with 47% dropout will show far more spurious absence than one with 0.4%,
            for reasons that have nothing to do with parentage.
          </Text>
          <Text mb={8}>
            The bound this tool uses is:
          </Text>
          <Wide>
            <Table withTableBorder>
              <Table.Tbody>
                <Table.Tr>
                  <Table.Td className="om-mono" style={{ fontSize: 13 }}>
                    explainable absence = no-call rate &times; heterozygous fraction + 0.005
                  </Table.Td>
                </Table.Tr>
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            It follows directly from <SecRef id="array" />. Dropout manufactures paternal absence
            only by turning a genuinely heterozygous call homozygous while discarding the paternal
            allele. So the rate at which it can do so scales with <b>both</b> how much is dropping
            out and how much of the genome is heterozygous in the first place. Neither term alone
            is a bound: a genome that is homozygous everywhere is immune however much drops out,
            and a genome with no dropout has nothing to inflate.
          </Text>
          <Text mb={8}>
            The trailing 0.005 is a residual floor for genotyping error on clean data, measured at
            0.03% and 0.05% on files with essentially no dropout. Without it a near-perfect array
            would have a reference of zero and any single error would exceed it.
          </Text>
          <Title order={3} mt={14} mb={4}>How tightly it holds</Title>
          <Text mb={8}>
            The bound has to cover every true parent-offspring pair without being so loose that it
            covers unrelated ones. Measured against pairs of known relationship:
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th>Sample</Table.Th>
                  <Table.Th ta="right">No-call</Table.Th>
                  <Table.Th ta="right">Het</Table.Th>
                  <Table.Th ta="right">Bound predicts</Table.Th>
                  <Table.Th ta="right">Observed</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <NumRow c={['PennCNV offspring', '0.4%', '31.8%', '0.13%', '0.05%']} />
                <NumRow c={['Androgenote 02', '13.6%', '1.3%', '0.18%', '0.16%']} />
                <NumRow c={['Androgenote 04', '13.7%', '2.0%', '0.27%', '0.22%']} />
                <NumRow c={['Zuccaro A8', '46.9%', '22.2%', '10.4%', '9.69%']} />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            Every one is covered, and none by more than a factor of about two and a half. The last
            row is the important one: a true father-offspring pair scoring 9.69% absence, which is
            higher than some genuinely unrelated pairs score. Only a per-sample bound puts that
            number in the right place.
          </Text>
          <Text>
            A call of &quot;absent&quot; requires exceeding this bound <b>three-fold</b>. See{' '}
            <SecRef id="band" />.
          </Text>
        </Section>

        {/* --- 9 ------------------------------------------------------------------------- */}
        <Section id="unrelated" title="The constant that had to be abandoned">
          <Text mb={8}>
            This section documents a method that was implemented, measured, and discarded. It is
            here because the alternative is the obvious design, and anyone evaluating this tool
            should know it was tried and why it fails.
          </Text>
          <Text mb={8}>
            The obvious approach: two unrelated people mismatch at some characteristic rate, so
            measure that rate once and compare every sample against it. A sample near it is
            unrelated; a sample well below it is the child.
          </Text>
          <Text mb={8}>
            Measured across more than fifty pairs, the &quot;unrelated&quot; rate ranged from{' '}
            <b>6.8% to 50.3%</b>. It is not a constant. It depends on how much ancestry the two
            people happen to share, and on how noisy the file standing in for the father is.
          </Text>
          <Alert color="orange" p="xs" mb={10}>
            <Text size="xs">
              Worse, the two populations overlap. Zuccaro embryo A8 scores 9.69% against its own
              father. Genuinely unrelated pairs in the same dataset score 5.1% to 5.6%. A single
              threshold placed anywhere gets one of those two wrong, and it is not a matter of
              choosing the threshold better: no threshold on that axis separates them.
            </Text>
          </Alert>
          <Text mb={8}>
            The per-sample bound resolves this because it asks a different question. Not &quot;is
            this rate high compared to other people&apos;s rates&quot; but &quot;is this rate
            higher than <i>this file&apos;s own quality</i> can account for&quot;. A8&apos;s 9.69%
            is unremarkable for a sample with 47% dropout and 22% heterozygosity. The same 9.69%
            from a clean array would be forty-fold over its bound and would be called absent.
          </Text>
          <Text>
            With that change the tool called 9 of 9 known relationships correctly, with no
            population constant anywhere in the code. See <SecRef id="validation" />.
          </Text>
        </Section>

        {/* --- 10 ------------------------------------------------------------------------ */}
        <Section id="classes" title="Reading the verdict">
          <Text mb={8}>
            The class is assembled from the three axes in a fixed order, and zygosity is consulted
            first because it constrains what the others can mean.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th style={{ width: 150 }}>Class</Table.Th>
                  <Table.Th>Reached when</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <Table.Tr>
                  <Table.Td><b>androgenetic</b></Table.Td>
                  <Table.Td>
                    The genome is homozygous (BAF band at or below 8%) and the paternal absence
                    rate is at or below the bound. One allele per position, and those alleles are
                    his, so there is no room for a second complement. Zygosity settles this without
                    consulting the narrower axis. Also reached when the genome is diploid, paternal
                    alleles are present, and alleles he lacks fall below the second-parent
                    expectation, which is consistent with two sperm; that case is flagged as
                    provisional in the result.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>gynogenetic</b></Table.Td>
                  <Table.Td>
                    Paternal absence exceeds the bound three-fold. No paternal contribution
                    anywhere. With only the sperm donor supplied this is an inference from his
                    absence, not a measurement of her presence.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>biparental</b></Table.Td>
                  <Table.Td>
                    The genome is diploid, paternal alleles are present, and alleles he lacks
                    exceed half his heterozygosity.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>neither parent</b></Table.Td>
                  <Table.Td>
                    Both parents supplied, and both are absent. Unreachable without an oocyte donor
                    array. See <SecRef id="oocyte" />.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>unclear</b></Table.Td>
                  <Table.Td>
                    The evidence does not separate the options. Always accompanied by a stated
                    limit naming what would settle it. See <SecRef id="band" />.
                  </Table.Td>
                </Table.Tr>
              </Table.Tbody>
            </Table>
          </Wide>
          <Text>
            Alongside the class the card reports zygosity, sperm type, the informative marker
            count, and both the observed rate and the reference it was measured against, so the
            margin is visible rather than implied.
          </Text>
        </Section>

        {/* --- 11 ------------------------------------------------------------------------ */}
        <Section id="band" title="The uncalled band">
          <Text mb={8}>
            Between the bound and three times it, nothing is called. This is a deliberate refusal
            and not an oversight.
          </Text>
          <List size="sm" spacing={4} mb={10}>
            <List.Item>
              At or below the bound: <b>present</b>. The absence observed is no more than this
              sample&apos;s own noise can produce.
            </List.Item>
            <List.Item>
              Above three times the bound: <b>absent</b>.
            </List.Item>
            <List.Item>
              In between: <b>uncalled</b>, with the reason stated.
            </List.Item>
          </List>
          <Text mb={8}>
            The band exists because a noisy diploid sample can reach the unrelated range against
            its own parent, as A8 does. In that region the measurement genuinely does not
            distinguish the hypotheses, and reporting a call would be reporting a coin flip with a
            number attached.
          </Text>
          <Text mb={8}>
            The result says what would resolve it, which is almost always the same thing: the other
            parent&apos;s array. That measures dropout directly instead of bounding it, and
            it turns an inference into a measurement.
          </Text>
          <Alert color="orange" p="xs">
            <Text size="xs">
              Adding data can move a call <i>into</i> the uncalled band. If a sample&apos;s array
              is poor enough that the maternal side cannot be resolved on it, supplying the oocyte
              donor changes a confident single-parent call into an honest &quot;unclear&quot;. That
              is the tool working, not failing, and the per-parent rates in the detail view say
              which side is unresolved and by how much.
            </Text>
          </Alert>
        </Section>

        {/* --- 12 ------------------------------------------------------------------------ */}
        <Section id="chrx" title="Sperm type from chrX, not chrY">
          <Text mb={8}>
            The conventional way to ask whether a sperm genome is present is to look for a Y
            chromosome. It does not work, and understanding why is worth a minute because it is
            the concrete thing this tool does that a chrY test cannot.
          </Text>
          <Text mb={8}>
            A sperm carries either an X or a Y, never both. So a sample containing a complete
            paternal genome delivered by an X-bearing sperm has <b>no Y at all</b>. A chrY test
            reads that identically to a sample with no paternal contribution whatsoever. The two
            are opposite findings and the test cannot separate them.
          </Text>
          <Text mb={8}>
            Syngamy reads chrX instead, using the same obligate-allele logic as everywhere else. If
            the sperm donor&apos;s alleles are present on chrX as well as the autosomes, the sperm
            was X-bearing. If they are present on the autosomes but absent on chrX <i>and</i> the
            sample calls a chrY, it was Y-bearing and the sample simply has no paternal X, which
            is ordinary sex determination rather than a loss. That case is labelled{' '}
            <Code>expected absent</Code> and excluded from the loss count.
          </Text>
          <Text mb={8}>
            The chrY half of that is a measurement, not an inference, and the distinction is
            load-bearing. Reading &quot;he must have sent a Y&quot; from the missing X alone is
            circular, and on a sample carrying no chrY at all it would exempt a genuine paternal X
            deletion as ordinary sex determination while asserting a Y chromosome the file does
            not show. Where the paternal X is absent and no chrY is called, the loss is reported
            as a loss and no sperm type is claimed.
          </Text>
          <Alert color="blue" p="xs" mb={10}>
            <Text size="xs">
              The exemption is strictly paternal. A mother transmits an X to a child of either sex,
              so an absent maternal X is a real finding and is never exempted. The code carries the
              two roles separately for exactly this reason.
            </Text>
          </Alert>
          <Text>
            On one real four-sample set, three of the four had no chrY, so the lab&apos;s existing
            XY test could not separate them. The SNP measurement separated them by 15 to 40-fold:
            one gynogenetic, two androgenetic from X-bearing sperm.
          </Text>
        </Section>

        {/* --- 13 ------------------------------------------------------------------------ */}
        <Section id="oocyte" title="Adding the oocyte donor">
          <Text mb={8}>
            An oocyte donor array is optional. It is also the single highest-value thing you can
            add, and this section is precise about what it buys because &quot;more data is
            better&quot; is not an argument.
          </Text>
          <Text mb={8}>
            With her array, each sample is tallied against both parents in the same pass, by the
            same measurement, and the class is read off the pair. The two rates are comparable only
            because the instrument is identical: the ceiling is a property of the <i>sample</i>, so
            both parents are measured against the same one.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th>Paternal</Table.Th>
                  <Table.Th>Maternal</Table.Th>
                  <Table.Th>Class</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <NumRow plain c={['present', 'present', 'biparental']} />
                <NumRow plain c={['present', 'absent', 'androgenetic']} />
                <NumRow plain c={['absent', 'present', 'gynogenetic']} />
                <NumRow plain c={['absent', 'absent', 'neither parent']} />
                <NumRow plain c={['either unresolved', '', 'unclear']} />
              </Table.Tbody>
            </Table>
          </Wide>
          <Title order={3} mt={14} mb={4}>What only the pair can do</Title>
          <List size="sm" spacing={6} mb={10}>
            <List.Item>
              <b>Separate gynogenetic from neither-parent.</b> With the father alone, a maternal
              origin is the shape left behind by his absence. It is an inference, and a sample
              belonging to some third person produces exactly the same shape. Measuring her
              directly separates them. On the bundled examples an unrelated oocyte donor reads as
              gynogenetic against the sperm donor alone and as <b>neither parent</b> once a second
              parent is measured.
            </List.Item>
            <List.Item>
              <b>Confirm an androgenote rather than infer it.</b> The maternal complement being
              absent becomes a measurement against her own array instead of a deduction from allele
              sharing.
            </List.Item>
            <List.Item>
              <b>Catch a duplicate parent.</b> Two arrays of one person pass both tests and yield a
              confident biparental call. The tool compares the two declared parents against each
              other and says so if they agree at over 99% of shared markers, which is a relabelled
              file rather than a pedigree.
            </List.Item>
            <List.Item>
              <b>Enable the per-locus test in the backend.</b> Without her, paternal presence
              cannot be established at any individual marker at all. See <SecRef id="backend" />.
            </List.Item>
          </List>
          <Alert color="orange" p="xs">
            <Text size="xs">
              A caveat worth knowing. When one parent is confirmed present and the other is
              unresolved, the pair is reported as <b>unclear</b>, which discards the confirmed
              half: the tool knows a paternal genome is there and only the maternal question is
              open. The per-parent rates are all shown, so nothing is hidden, but the class name
              alone does not carry it.
            </Text>
          </Alert>
        </Section>

        {/* --- 14 ------------------------------------------------------------------------ */}
        <Section id="locus" title="The per-locus deletion test">
          <Text mb={8}>
            Everything above is a rate across the genome. This asks a different question about one
            place: is the paternal contribution specifically absent <i>around this variant</i>. It
            appears under the results once a run has finished, as a box taking a position.
          </Text>
          <Text mb={8}>
            Type the variant&apos;s coordinate in the array&apos;s own assembly, as{' '}
            <Code>chr7:117559590</Code>. Each sample is re-read, one chromosome of it, and scored
            marker by marker: 1 where the sperm donor&apos;s allele is provably present, 0 where it
            is provably absent, and nothing at all where the marker cannot say. The statistic is
            the length of the longest <b>run</b> of consecutive absences, against the longest run
            that independent genotyping error would produce at the same marker count.
          </Text>
          <Text mb={8}>
            Contiguity is the whole point. Scattered absences are error; a deletion removes a
            stretch. Three windows are reported, 25 kb, 10 Mb and the whole chromosome, because an
            event smaller than the local marker spacing cannot produce a run at all.
          </Text>
          <Text mb={8}>
            The threshold depends on how often a marker shows paternal absence when the paternal
            genome <i>is</i> present, and that rate is taken from this sample&apos;s own
            genome-wide absence rather than from a published floor. It matters: the floor is a
            lower bound, so using it on amplified material demands a shorter run than the data
            warrants. At a 5% spurious rate the required run is 4 markers where the floor asks 2,
            and at 30% it is 8. The required run is never below 2 in any case, because contiguity
            is the discriminator and a run of one has none to offer.
          </Text>
          <Alert color="orange" p="xs" mb={10}>
            <Text size="xs">
              <b>This test needs the oocyte donor and refuses without her.</b> Absence is Mendelian
              and needs nothing from her: a homozygous father must transmit his allele, so an
              embryo lacking it is missing the paternal contribution whatever she carries. Presence
              is an identity claim, and the embryo carrying his allele only shows it came from{' '}
              <i>him</i> if she could not have supplied it, which needs her homozygous for the
              other allele. Without that, nothing can score as present, so nothing can break a run,
              and the statistic has no null to be significant against. In that mode a normal
              chromosome 20 produced a run of 3 across 35 Mb at p = 2.5e-07.
            </Text>
          </Alert>
          <Text mb={8}>
            Supplying an event size you care about turns an absent run into an explicit{' '}
            <b>below resolution</b> rather than letting it read as no event. That matters: on a
            typical panel the smallest resolvable event is r_min times the local marker spacing,
            which is tens of kilobases, so a 20 kb deletion cannot produce a significant run no
            matter how real it is.
          </Text>
          <Alert color="blue" p="xs">
            <Text size="xs">
              What the run does not tell you: it does not separate copy loss from copy-neutral
              loss of heterozygosity, since both remove paternal alleles contiguously and only the
              intensity channel distinguishes them. The marker count is printed beside every
              verdict because a window holding few informative markers carries little power even
              when it clears the threshold.
            </Text>
          </Alert>
        </Section>

        {/* --- 15 ------------------------------------------------------------------------ */}
        <Section id="confidence" title="Confidence, and why no percentage">
          <Text mb={8}>
            The tool never prints &quot;87% confident&quot;. That is a considered refusal, and
            since it is the first thing a reviewer will ask about, here is the full reasoning.
          </Text>
          <Title order={3} mt={14} mb={4}>The internal number is meaningless as a probability</Title>
          <Text mb={8}>
            The likelihood ratio between &quot;paternal genome present&quot; and &quot;absent&quot;,
            computed across a genome-wide array, runs past 10<sup>10000</sup> on a typical sample.
            Converting that to a percentage gives 100.000...% with thousands of nines. It is not
            wrong, and it is not information: it measures how many markers were read, not whether
            the answer is right.
          </Text>
          <Title order={3} mt={14} mb={4}>The honest number is empirical accuracy</Title>
          <Text mb={8}>
            What a reader actually wants is: when this tool makes a call, how often is it correct?
            That is a property of the method against ground truth, not of one sample. On pairs of
            known relationship it stands at <b>9 of 9 correct</b>. The 95% lower confidence bound
            on 9 of 9 is roughly 72%.
          </Text>
          <Text mb={8}>
            To claim 99% accuracy with the same confidence would take about 300 consecutive correct
            calls; for 95%, about 59. Printing any percentage before then would be asserting
            precision the validation does not support.
          </Text>
          <Title order={3} mt={14} mb={4}>What is reported instead</Title>
          <Text>
            The margin. Every call shows the observed rate, the reference it was measured against,
            and the ratio between them. A call at 30x its ceiling and a call at 3.1x are both
            &quot;absent&quot;, and a reader can see instantly that they are not equally
            comfortable. That is information a single percentage would destroy.
          </Text>
          <Text mt={8}>
            Under the ratio sit the two figures the ceiling is the product of: this sample&apos;s
            no-call rate and its heterozygous fraction. A ratio that fell between two runs can
            mean the signal shrank or the ceiling grew, and those have opposite implications for
            whether to repeat the array. The factors say which one moved.
          </Text>
        </Section>

        {/* --- 15 ------------------------------------------------------------------------ */}
        <Section id="gates" title="Quality gates">
          <Text mb={8}>
            One gate runs per CHROMOSOME rather than per file, and it is the only one that can
            withhold a chromosome whose file passed everything else. If fewer than 40% of a
            chromosome&rsquo;s B-allele frequencies sit at an extreme, outside 0.15 to 0.85, its
            probes clustered badly and its genotype calls are not measuring it. No verdict is
            reported for that chromosome and the reason is printed in its place.
          </Text>
          <Text mb={8}>
            This catches something no allelic statistic can see. When a chromosome mis-clusters,
            every measure built on genotypes agrees that the parental allele is missing, because
            they are all reading the same broken calls, and the result is indistinguishable from a
            real chromosome-scale loss. Measured across the five paternal pronuclei of the bundled
            series, each a meiotic product of the sperm donor and so present on every autosome:
            correctly clustered chromosomes run 75.2% to 97.6% extreme, while one chromosome ran
            13.0% and was reported <b>absent at 18.36%</b> on a sample whose genome-wide verdict is
            present. The gap between those is empty by a factor of 5.8. A genuinely trisomic
            chromosome also fails this gate, since its allelic ratios cluster at a third and two
            thirds; that is the right outcome, because gains are not called on this platform at all.
          </Text>
          <Text mb={8}>
            Every file is also put through six gates on load. They report; they do not silently drop
            anything. Each returns <Code>usable</Code>, <Code>marginal</Code>,{' '}
            <Code>exclude</Code> or <Code>report only</Code>, and all six appear in the detail view
            and the PDF.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th style={{ width: 150 }}>Gate</Table.Th>
                  <Table.Th>What it catches</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <Table.Tr>
                  <Table.Td><b>call rate</b></Table.Td>
                  <Table.Td>
                    Below 60%, excluded, following the only published threshold measured on
                    amplified material <Ref id="natesan_2014" />. 60 to 75% is marginal. Note the
                    Axiom vendor gate of 97% was established on unamplified bulk DNA and would
                    reject nearly every embryo biopsy, which is why it is not used here.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>het-to-hom asymmetry valid</b></Table.Td>
                  <Table.Td>
                    Below 60% call rate, erroneous heterozygous calls become common, so a
                    heterozygous call is no longer robust and any partition resting on it loses its
                    guarantee. The gate suspends that assumption rather than quietly relying on it.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>genome-wide LOH</b></Table.Td>
                  <Table.Td>
                    Heterozygosity near zero across every chromosome: a candidate abnormal
                    fertilisation or failed genome unification. The published criterion is
                    qualitative and clustered by embryo, which a single unlabelled file cannot
                    reproduce, so this flags a candidate rather than reproducing that filter.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>heterozygosity plausible</b></Table.Td>
                  <Table.Td>
                    A diploid human genome cannot be 56% heterozygous. One real amplified
                    blastomere measured exactly that at 67% call rate, passed every other gate as
                    merely marginal, and produced a sex call contradicting the other blastomere of
                    the same embryo. An upper bound catches it directly: above 50% excluded, above
                    40% marginal.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>numeric genotype coding</b></Table.Td>
                  <Table.Td>
                    Axiom&apos;s 0/1/2/-1 export leaves which numeral means which homozygote as a
                    convention rather than a fact. Unlike the nucleotide convention this one is
                    detectable, by checking the numeric calls against mean BAF, and an inverted
                    file is excluded rather than analysed backwards.
                  </Table.Td>
                </Table.Tr>
                <Table.Tr>
                  <Table.Td><b>sex call</b></Table.Td>
                  <Table.Td>
                    From the chrX to autosome heterozygosity ratio, near 0.02 for male
                    (pseudoautosomal only) and 0.7 to 0.8 for female. Dropout cancels in a ratio,
                    so it survives amplified material. A value between the bands is reported and
                    not resolved: a chrX copy-number abnormality, a pooled file and contamination
                    all land there, and guessing is the one place this could invert a declared
                    role.
                  </Table.Td>
                </Table.Tr>
              </Table.Tbody>
            </Table>
          </Wide>
          <Text>
            A gate reading <Code>exclude</Code> does not stop the run. It appears on the sample,
            in the detail view and in the report, and it is your decision what to do about it. The
            tool refuses to make quality judgements silently.
          </Text>
        </Section>

        {/* --- 16 ------------------------------------------------------------------------ */}
        <Section id="assembly" title="Assembly detection, and no liftOver">
          <Text mb={8}>
            A marker position means nothing without knowing which reference assembly it is stated
            against, and array exports frequently do not say. Getting it wrong misplaces every
            marker and silently corrupts anything positional.
          </Text>
          <Text mb={8}>
            The assembly is detected from the positions themselves. Each marker is checked against
            the UCSC <Code>chromInfo</Code> and <Code>gap</Code> tables for GRCh37 and GRCh38: a
            position past the end of its chromosome, or inside an assembly gap, is illegal under
            that build. Counting illegal placements under each and comparing gives the answer,
            with the counts and the number tested reported rather than just a verdict.
          </Text>
          <Alert color="orange" p="xs" mb={10}>
            <Text size="xs">
              <b>No liftOver is performed, and no chain file is shipped.</b> UCSC chain files carry
              a non-commercial field-of-use restriction that an Apache 2.0 project cannot
              sublicense. So a build mismatch is <i>detected and reported</i>, never silently
              corrected. Where the mismatch matters, the affected step refuses rather than
              proceeding on coordinates it cannot trust.
            </Text>
          </Alert>
          <Text>
            This has fired on real data. Lab Axiom exports came back GRCh37 while the bundled
            genetic maps are GRCh38. The scaffold refused to read the map at those coordinates, its
            reported odds dropped from 10<sup>8</sup> to 10<sup>6</sup>, and the calls themselves
            were unchanged. That drop is the honest number: the confidence had been resting on a
            map read at the wrong positions.
          </Text>
        </Section>

        {/* --- 17 ------------------------------------------------------------------------ */}
        <Section id="chrom" title="Per-chromosome and segmental results">
          <Text mb={8}>
            A genome-wide rate hides the finding that matters most. Zuccaro reports Cas9 cleavage
            removing segments up to a whole chromosome arm <Ref id="zuccaro_2020" />, and a sample
            missing the paternal contribution on one chromosome still has a low genome-wide
            absence rate.
          </Text>
          <Text mb={8}>
            So every chromosome is reported separately: informative markers, missing count, rate,
            the multiple of the ceiling, and a verdict. Chromosomes with fewer than 200 informative
            markers are omitted, because the rate there is too noisy to place against the ceiling
            at all.
          </Text>
          <Title order={3} mt={14} mb={4}>The pseudoautosomal region</Title>
          <Text mb={8}>
            Chromosome X is reported in two rows rather than one. PAR1 and PAR2 sit on both the X
            and the Y and recombine between them, so a Y-bearing sperm still delivers the father&apos;s
            PAR alleles while the rest of his X is legitimately gone. Pooled into a single chrX
            figure that distinction is invisible; split, the <code>X:PAR</code> row reads present
            on exactly the samples where <code>X</code> reads expected_absent, which is a positive
            control on the sample rather than an inference about it. The boundaries come from the
            GRC assembly region reports, and the split is skipped when the assembly could not be
            determined: PAR2 is at a different address in GRCh37 and GRCh38, and guessing would
            put ordinary chrX markers into a control. Roughly 976 informative markers fall in the
            PAR on a full Axiom array, so on the one-in-eight example subsets the row falls under
            the 200-marker floor and is omitted.
          </Text>
          <Title order={3} mt={14} mb={4}>Spread between chromosomes</Title>
          <Text mb={8}>
            The report gives the coefficient of variation of the per-chromosome rate alongside the
            cleanest chromosome. Measured on the bundled arrays: 1.11 for two arrays of one man,
            0.76 for a degraded but true parent-offspring pair, 0.11 for an unrelated adult, 0.10
            for a 50:50 blend of two unrelated genomes. A relationship interrupted by a loss is
            patchy and reaches the ceiling somewhere; a genome that differs everywhere by the same
            amount does not.
          </Text>
          <Text mb={8}>
            Where a sample already sits in the uncalled band and its absence is uniform with no
            chromosome near the ceiling, the report says the shape is two genomes blended rather
            than one genome missing pieces, and that no reanalysis separates them. It is worth
            being explicit about what this does not do: uniformity alone does not identify a
            mixture. An unrelated adult measures 0.112 and a blend 0.104, which is no separation
            at all. The statistic distinguishes a genome-wide difference from one confined to part
            of the genome, and it is consulted only where the rate is already ambiguous.
          </Text>
          <Title order={3} mt={14} mb={4}>Below the chromosome</Title>
          <Text mb={8}>
            The backend goes finer, and two constraints keep the finer analysis honest:
          </Text>
          <List size="sm" spacing={6} mb={10}>
            <List.Item>
              <b>Physical compactness.</b> A run of consecutive absent markers is split wherever
              the gap between them exceeds ten times the local informative spacing. Adjacency in
              the list of informative markers is not adjacency on the chromosome, and marker
              deserts would otherwise let a &quot;run&quot; span a region containing no evidence.
            </List.Item>
            <List.Item>
              <b>Whole-chromosome detection keyed to completeness.</b> A whole-chromosome loss is
              recognised by every informative marker being absent, not by the longest run being
              long. Deserts fragment runs, so a genuine whole-chromosome loss can present as
              several runs rather than one.
            </List.Item>
          </List>
          <Text>
            Segment rates are anchored to the calibrated noise bound rather than estimated from the
            distribution of observed window rates. An earlier version estimated them as the 10th
            and 90th percentiles of that distribution, which finds two populations even in a
            uniform genome: it carved a sample present everywhere at 0.16% into 28 &quot;absent&quot;
            segments, and a gynogenetic sample into 720. That is what a percentile-based
            segmenter does when there is nothing to segment.
          </Text>
        </Section>

        {/* --- 18 ------------------------------------------------------------------------ */}
        <Section id="segments" title="Chromosomal change: segments">
          <Text mb={8}>
            A chromosome that is only PARTLY lost is the case the per-chromosome verdict cannot
            report. Its absence rate is the average of a missing region and a present one, so it
            lands between the two references and the chromosome comes back <b>unclear</b>. That
            reads as an absence of information when it is in fact a located event, and a reader
            skimming a table of chromosomes has no reason to look further. So segments are scanned
            for separately and stated at the top of the sample, in the same visual language as a
            warning, before the tables that carry the numbers.
          </Text>
          <Text mb={8}>
            <b>What is reported.</b> Each region carries its chromosome and coordinates, its
            physical span, the number of called informative markers inside it, the local absence
            rate, the rate it was scored against, and its score. It is a <b>loss of the paternal
            contribution</b> over that region. It is not a statement about physical copy number:
            whether the chromosome is deleted or present in two copies from one parent needs
            intensity, which this platform cannot support. And it is not a gain. Gains are not
            called here at all, in either channel.
          </Text>
          <Title order={3} mt={14} mb={4}>The null is external, and that is the whole design</Title>
          <Text mb={8}>
            Each window is scored against the median per-chromosome absence rate of that
            sample&rsquo;s OTHER chromosomes, floored at 0.2%. Scoring against the sample&rsquo;s
            own genome-wide rate fails exactly where it matters: a large event inflates the rate it
            is being tested against, so the biggest losses score lowest. Taking a median rather
            than a mean is what stops several events from doing it collectively, and excluding the
            chromosome under test is what stops one from doing it alone.
          </Text>
          <Text mb={8}>
            A consequence worth stating: a genome where EVERY chromosome is replaced produces no
            segments, because there is nothing clean to form a null from. That is correct rather
            than a gap. A wholly non-paternal genome is what the genome-wide verdict is for, and it
            calls it.
          </Text>
          <Title order={3} mt={14} mb={4}>The threshold is empirical, never a closed-form tail</Title>
          <Text mb={8}>
            A Bonferroni-corrected exact binomial and an Erdos-Renyi run-length scan both fail
            here, for the same reason the per-locus p-value did: absence artefact on amplified
            material is spatially clustered, and a tail computed under independence reports
            artefact as significant. Measured on this platform, runs of two occur 6 to 10 times
            more often than independence predicts on single cells and about twice as often on bulk.
          </Text>
          <Text mb={8}>
            So the threshold comes from genomes known to carry nothing. Five paternal pronuclei of
            the sperm donor, each a meiotic product of his and so present on every autosome, reach
            a maximum score of 139 over 110 chromosome scans. The calling threshold is 250.
            <b> False segments on those genomes: 0 of 110.</b> That threshold is fitted on the same
            null it is evaluated against, so it awaits an out-of-sample clean cohort, and the score
            is printed beside every segment so a marginal call reads as marginal.
          </Text>
          <Title order={3} mt={14} mb={4}>The floor, and why a marker count is not a resolution</Title>
          <Text mb={8}>
            Titrated on real material rather than simulated noise: a block of a maternal
            pronucleus spliced into a clean paternal one of the same series, so the segment is a
            genuine alternative genome carrying real amplification artefact. Twelve constructions
            per size, four backgrounds by three donors.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th>Spliced markers</Table.Th>
                  <Table.Th style={{ textAlign: 'right' }}>Score</Table.Th>
                  <Table.Th style={{ textAlign: 'right' }}>Against null max 139</Table.Th>
                  <Table.Th>Called</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <NumRow c={['1,200', '152 to 197', '1.09x', 'no']} plain />
                <NumRow c={['2,400', '431 to 556', '3.1x', 'yes']} plain />
                <NumRow c={['4,800', '968 to 1319', '7.0x', 'yes']} plain />
                <NumRow c={['9,600', '1807 to 2450', '13x', 'yes']} plain />
                <NumRow c={['19,200', '3489 to 4582', '25x', 'yes']} plain />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            The floor is 2,400 called informative markers, where separation becomes real rather
            than merely nonzero. The previous 200-marker floor was not a weak test, it was no test
            at all.
          </Text>
          <Text mb={8}>
            <b>The span is the resolution, not the event.</b> Marker spacing on this array runs
            from about 1 bp to 21 kb, so a fixed window of 2,400 markers covers anywhere from a
            few hundred kilobases to tens of megabases depending where it falls. A real event
            smaller than the window that found it is reported at the size of that window. Both
            numbers travel with every segment for that reason, and only one of them is a
            resolution.
          </Text>
          <Text mb={8}>
            A chromosome withheld by the allelic-ratio gate of <SecRef id="gates" /> is not
            scanned either. Its genotype calls are not measuring it, and the same broken calls
            would produce a confident segment inside it.
          </Text>

          <Title order={3} mt={16} mb={4}>
            Why there is no gain, and no parent attached to an extra chromosome
          </Title>
          <Text mb={8}>
            The obvious next question is the opposite one: a sample has TWO copies of a chromosome
            where it should have one, so which parent supplied the extra? That was built and
            measured, and it is refused. The measurements are below because a refusal without one
            is just a preference.
          </Text>
          <Text mb={8}>
            Three mechanisms produce two copies in a nominally haploid sample: his one copy
            replicated, both his homologues retained, or a chromosome from someone else kept. The
            first is invisible by construction. Duplicating a chromosome changes no genotype at any
            marker, so allelic evidence cannot see it at all, and the copy count itself needs
            intensity.
          </Text>
          <Text mb={8}>
            The other two were constructed from real arrays: a haploid paternal pronucleus as the
            background, with one chromosome merged either with a second pronucleus of the SAME
            father, which is two independent draws of his genome, or with a maternal pronucleus,
            which is a different person. Heterozygosity at the father&rsquo;s heterozygous sites
            says both his homologues are there; heterozygosity at his HOMOZYGOUS sites says the
            extra material is not his, since he cannot supply a second allele where he has only
            one. Each ratio is self-normalised against the sample&rsquo;s own other chromosomes.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th>Construction</Table.Th>
                  <Table.Th style={{ textAlign: 'right' }}>Het at father-het</Table.Th>
                  <Table.Th style={{ textAlign: 'right' }}>Het at father-hom</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                <NumRow c={['nothing added, 88 chromosome observations',
                  'median 1.01, max 2.03', 'median 1.01, max 1.91']} plain />
                <NumRow c={['his copy duplicated', 'no change at all', 'no change at all']} plain />
                <NumRow c={['both his homologues', '1.46 to 7.25', '1.19 to 2.17']} plain />
                <NumRow c={['a chromosome from someone else', '3.74 to 4.99', '2.18 to 3.17']} plain />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            The classes overlap the noise. A chromosome carrying no extra material at all reaches
            2.03 on the statistic, while the weakest real second-genome construction reaches 1.46,
            so the null runs straight through the signal. And the two mechanisms that ARE both real
            events separate by 2.17 against 2.18, a margin of half a percent. Neither call is
            available.
          </Text>
          <Text mb={8}>
            One further limit sits behind all of that, and it would apply even with clean
            separation. The statistic detects material that is <b>not paternal</b>. It does not
            detect material that is <b>maternal</b>: measured against the true mother and against
            an unrelated adult at matched array quality, the two are indistinguishable. Dispermy, a
            second sperm, contamination from another man in the lab, a swapped sample and a
            genuinely retained maternal chromosome all produce the same signature. So even where an
            extra chromosome is detectable, calling it maternal is an inference from an assumed
            pedigree rather than a measurement, and it would need the oocyte donor&rsquo;s own array
            to become one.
          </Text>
        </Section>

        <Section id="report" title="The report, and how to cite it">
          <Text mb={8}>
            <b>Report (PDF)</b> writes the whole run out as a paginated document intended for a
            supplementary file. It is typeset in the browser, for the same reason the genotypes are
            read there: there is no server to send them to.
          </Text>
          <Text mb={8}>It contains, in order:</Text>
          <List size="sm" spacing={4} mb={10} type="ordered">
            <List.Item>
              Run identity: report id, generation time, run start, software version, both donors
              with marker counts and hash prefixes, donor heterozygosity, assembly, and a statement
              of how the data was handled.
            </List.Item>
            <List.Item>
              A summary table, one row per sample, with each margin against its own ceiling.
            </List.Item>
            <List.Item>
              Per sample: the class in a filled bar, every measured rate, an evidence table per
              parent giving observed value, reference and what the reference is, a full
              per-chromosome table per parent, the complete quality read, all six gates, and any
              findings and limits.
            </List.Item>
            <List.Item>
              Methods, written from the run rather than boilerplate, naming which mode it ran in
              and which limits actually fired.
            </List.Item>
            <List.Item>
              Interpretation, including the confidence reasoning in <SecRef id="confidence" />.
            </List.Item>
            <List.Item>
              Every constant with its provenance, so no number in the document is unattributed.
            </List.Item>
            <List.Item>
              Each input file with its size, marker count and full SHA-256, so a reviewer can
              confirm the inputs without being sent them.
            </List.Item>
            <List.Item>Citations.</List.Item>
          </List>
          <Text mb={8}>
            Every page carries the research-use disclaimer and a footer identifying the tool,
            version, report id, generation time and page number. The report id is derived from the
            file hashes and the run start, so two runs cannot collide, and the filename carries the
            same id as the document header.
          </Text>
          <Title order={3} mt={14} mb={4}>The run log</Title>
          <Text>
            The log is deliberately <i>not</i> in the PDF. It downloads separately as a text file
            with a build and debug footer, which is the right shape for pasting into a bug report
            and the wrong shape for a supplementary file.
          </Text>
        </Section>

        {/* --- 19 ------------------------------------------------------------------------ */}
        <Section id="examples" title="The bundled example data">
          <Text mb={8}>
            <b>Examples</b> loads {EXAMPLES.length} public arrays so the tool can be tried without
            a family&apos;s data. They come from {EXAMPLE_CITATION}
          </Text>
          <Text mb={8}>
            Each is a subset: every {ordinal(EXAMPLE_STRIDE)} marker, {int(EXAMPLE_MARKERS)} per
            file, with the two columns this tool never reads dropped. <b>No value is rounded or
            altered</b>, and the calls below are the ones the full 825,657-marker files produce, so
            the demonstration is the real one at a tenth of the download.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Thead>
                <Table.Tr>
                  <Table.Th style={{ width: 110 }}>Accession</Table.Th>
                  <Table.Th style={{ width: 70 }}>Role</Table.Th>
                  <Table.Th>What it is, and what it should say</Table.Th>
                </Table.Tr>
              </Table.Thead>
              <Table.Tbody>
                {EXAMPLES.map((e) => (
                  <Table.Tr key={e.gsm}>
                    <Table.Td className="om-mono">{e.gsm}</Table.Td>
                    <Table.Td>{e.role === 'donor' ? 'sperm' : e.role}</Table.Td>
                    <Table.Td>
                      <Text size="sm">{e.what}</Text>
                      <Text size="xs" c="dimmed">{e.expect}</Text>
                    </Table.Td>
                  </Table.Tr>
                ))}
              </Table.Tbody>
            </Table>
          </Wide>
          <Title order={3} mt={14} mb={4}>Try the oocyte path on them</Title>
          <Text mb={8}>
            Relabel <Code>GSM4472407_donor_A_47</Code> from sample to oocyte and run again. Two
            results change in instructive ways:
          </Text>
          <List size="sm" spacing={4} mb={10}>
            <List.Item>
              <Code>donor_C_70</Code> moves from <b>gynogenetic</b> to <b>neither parent</b>. She
              is unrelated to both, and only measuring a second parent can say so.
            </List.Item>
            <List.Item>
              <Code>A8_45</Code> stays <b>unclear</b>, and the oocyte donor shows a second reason
              why. Against the sperm donor alone its 53% call rate already withholds zygosity. Add
              her and the maternal side lands at 12.84% against a 10.91% ceiling, 1.2x, inside the
              uncalled band as well. That array has 47% dropout and resolves neither parent. See{' '}
              <SecRef id="band" />.
            </List.Item>
          </List>
          <Alert color="blue" p="xs">
            <Text size="xs">
              One example is labelled androgenetic and is not an embryo. GSM4472398 is a second
              array of the same sperm donor, so every allele traces back to him and the classifier
              says so correctly. On a bulk sperm sample read that as an identity match. It is in
              the set because it demonstrates the absence axis at its cleanest: 0.15% against a
              1.45% ceiling.
            </Text>
          </Alert>
        </Section>

        {/* --- 20 ------------------------------------------------------------------------ */}
        <Section id="backend" title="What the backend adds">
          <Text mb={8}>
            The browser runs the rate-based analysis in full: parent of origin, zygosity, sperm
            type, per-chromosome results, both parents. That is the whole north-star question and
            it needs nothing else.
          </Text>
          <Text mb={8}>
            The Python package goes further in three directions the page does not expose. The
            per-locus deletion test of <SecRef id="locus" /> runs in the browser too, and fits its
            spurious-violation rate the same way.
          </Text>

          <Title order={3} mt={14} mb={4}>Bidirectional cross-check</Title>
          <Text mb={8}>
            With both parents, a further discriminator becomes available. A genuine paternal
            deletion cannot produce maternal-absence calls in the same region, because the maternal
            chromosome is intact there. Dropout produces both. The test is a binomial likelihood
            ratio between those two hypotheses over the same window, with an explicit
            insufficient-evidence outcome when neither is supported.
          </Text>

          <Title order={3} mt={14} mb={4}>Sibling haplotype scaffold</Title>
          <Text mb={8}>
            Across several embryos from one couple, a two-state forward-backward chain over
            co-inheritance assigns which paternal homologue each sample received. Confidence
            saturates on the genetic map rather than on marker count, since two markers a
            centimorgan apart carry less independent information than their count suggests. It
            refuses to read the map at all on an assembly mismatch.
          </Text>

          <Title order={3} mt={14} mb={4}>Read-level and structural evidence</Title>
          <Text mb={8}>
            Structural-variant VCFs are read and matched to a locus, so a called deletion
            corroborates an absence run. Absence of a call is never reported as absence of an
            event, which is a different statement and the file cannot support it.
          </Text>

          <Title order={3} mt={14} mb={4}>Format normalisation</Title>
          <Text>
            The readers listed in <SecRef id="formats" /> live here, along with cross-file
            harmonisation. Nucleotide-to-AB assignment pools across every sample in the run rather
            than being decided per file, because a per-file decision can assign two files opposite
            conventions at the same marker.
          </Text>
        </Section>

        {/* --- 21 ------------------------------------------------------------------------ */}
        <Section id="cli" title="Command line">
          <Text mb={8}>
            The analysis runs standalone with no server and no browser. The parentage scan needs
            only a paternal file and one or more samples.
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Tbody>
                <ExRow k="--father" v="paternal array export. Required." mono />
                <ExRow k="--samples" v="one or more sample exports, scanned for parent of origin" mono />
                <ExRow k="--mother" v="maternal array export. Omit to run the degraded mode, which refuses the per-locus test." mono />
                <ExRow k="--embryo" v="a single embryo or biopsy export, for the per-locus path" mono />
                <ExRow k="--chrom, --pos" v="variant chromosome and position, in the array's build. Omit to scan parent of origin only." mono />
                <ExRow k="--product" v="array product name, which sets the chrX probe offset. Leave unspecified rather than guessing: a wrong product applies the wrong offset." mono />
                <ExRow k="--mutant-allele, --phase-route" v="which allele is the mutant one, and how that was established. The route is recorded verbatim rather than inferred." mono />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            Every module carries a runnable self-check that asserts against real measured values
            rather than against itself, and the browser arithmetic is pinned to the same fixtures
            as the Python so a divergence between them fails a test instead of reaching a result.
          </Text>
          <Text>
            The <Code>--product</Code> flag is worth care. It defaulted to a named array once, and
            that applied one platform&apos;s chrX probe offset to another platform&apos;s data,
            producing a plausible and wrong result. It now defaults to unspecified.
          </Text>
        </Section>

        {/* --- 22 ------------------------------------------------------------------------ */}
        <Section id="validation" title="Validation record">
          <Text mb={8}>
            What follows is what the method has actually been run against. It is deliberately
            small and specific, because a validation section that is vague is worse than none.
          </Text>
          <Title order={3} mt={14} mb={4}>Positive control</Title>
          <Text mb={8}>
            A published trio with a male offspring. His X is entirely maternal, by biology rather
            than by any paper&apos;s claim, so the paternal contribution on chrX must be absent and
            everywhere else present. Measured: chrX 4,324 of 4,324 informative markers absent, and
            autosomal absence at 0.05% against a 0.13% bound.
          </Text>
          <Text mb={8}>
            The same run surfaced two deletions not previously flagged: a de novo 71 kb deletion on
            chromosome 3 and an inherited one on chromosome 11, both confirmed independently by the
            log R ratio channel, which the parentage test never reads.
          </Text>
          <Title order={3} mt={14} mb={4}>The pronucleus audit</Title>
          <Text mb={8}>
            The strongest ground truth available for this question, and the reason it is worth
            more than a larger number of ordinary cases: Zuccaro et al. separated the two
            pronuclei of a fertilised zygote by micromanipulation and arrayed each one alone{' '}
            <Ref id="zuccaro_2020" />. A paternal-nucleus sample contains one parental genome and
            a maternal-nucleus sample the other, physically, before any statistic is applied.
          </Text>
          <Text mb={8}>
            Run over 27 arrays from that series, 26 cases with bench ground truth:
          </Text>
          <Wide>
            <Table striped withTableBorder>
              <Table.Tbody>
                <ExRow k="Correct" v="22" />
                <ExRow k="Refused" v="4, each declining to name a class and saying why" />
                <ExRow k="Incorrect" v="0" />
                <ExRow
                  k="Passing every gate"
                  v="22 of 23 correct, 1 refused, 0 incorrect"
                />
                <ExRow
                  k="Donor identification"
                  v="7 of 7 maternal pronuclei resolved to exactly one of four egg donors, the
                     runner-up 10 to 18 times further away"
                />
              </Table.Tbody>
            </Table>
          </Wide>
          <Text mb={8}>
            The four refusals are the informative part. Three are paternal pronuclei below the
            60% call rate, where the tool excludes the sample on its own gate and withholds
            zygosity; on the axis that matters it had all three right, detecting the paternal
            genome as present at 0.44x, 0.49x and 0.45x of their ceilings. The fourth is the
            noisiest of eight egg-donor arrays, landing at 1.08x inside the uncalled band while
            the other seven, unrelated to the sperm donor in exactly the same way, are called at
            4.8x to 7.7x.
          </Text>
          <Alert color="blue" p="xs" mb={10}>
            <Text size="xs">
              That audit found a real defect and it is worth stating plainly. Its first run
              returned <b>three incorrect calls</b>: haploid paternal pronuclei reported as
              biparental, because spurious heterozygous calls on a poor array pushed the BAF band
              over the diploid threshold. The tool had already excluded all three on its
              call-rate gate and already printed that heterozygous calls were untrustworthy
              there, and the classifier used them anyway. Zygosity is now withheld below that
              floor, in both implementations, and the class falls to unclear with the reason
              stated. The full record, and a report for every case, is in the repository under{' '}
              <Code>audit/</Code>.
            </Text>
          </Alert>
          <Title order={3} mt={14} mb={4}>Earlier known-relationship pairs</Title>
          <Text mb={8}>
            <b>9 of 9 correct</b> across pairs spanning true parent-offspring pairs, unrelated
            pairs, and replicate arrays of one individual, with no population constant anywhere in
            the code. See <SecRef id="confidence" /> for what that does and does not license.
          </Text>
          <Title order={3} mt={14} mb={4}>A correction worth recording</Title>
          <Text mb={8}>
            During validation, 35 absence runs genome-wide in a clean trio were briefly read as a
            failure, on the reasoning that it was far above the per-test analytic rate. That was
            wrong. A normal human genome carries hundreds of deletions, and 32 of the 35 had
            physically plausible spans. The error was comparing a genome-wide count against a
            per-test threshold. It is recorded here because a validation section that only lists
            successes is not evidence of anything.
          </Text>
          <Title order={3} mt={14} mb={4}>Bugs found and fixed by validation</Title>
          <List size="sm" spacing={4}>
            <List.Item>
              Presence was scored without consulting the mother, so at a paternal deletion, wherever
              a heterozygous mother transmitted the same allele, the marker read as present. This
              broke runs at roughly half of maternal heterozygous markers and turned a chrX loss of
              4,324 markers into a reported run of 21.
            </List.Item>
            <List.Item>
              Markers that could not speak were folded in as evidence of presence rather than
              excluded, which both broke runs and inflated the marker count.
            </List.Item>
            <List.Item>
              The segment reported was the longest anywhere on the chromosome rather than the one at
              the variant, so a normal chromosome 3 was reported at 7,499 markers and 359 standard
              deviations.
            </List.Item>
            <List.Item>
              A runner re-derived significance downstream of the gates, bypassing the
              no-mother refusal and reporting paternal loss on three samples whose paternal genome
              was present at 0.2%.
            </List.Item>
            <List.Item>
              The chrX probe offset was applied when deriving one quantity but not when
              interpreting chrX, so a male&apos;s single X read as two maternal copies.
            </List.Item>
            <List.Item>
              A wide multi-sample genotype file was read as its first sample only, dropping the rest
              silently. Found by the tool&apos;s own harmonisation no-op test.
            </List.Item>
          </List>
        </Section>

        {/* --- 23 ------------------------------------------------------------------------ */}
        <Section id="limits" title="Scope and limits">
          <Text mb={8}>
            <b>This is research decision support, not a clinical diagnostic.</b> Nothing here is
            validated for clinical use and no result should inform a transfer decision.
          </Text>
          <List size="sm" spacing={6} mb={10}>
            <List.Item>
              <b>It does not call gains, and does not attach a parent to an extra chromosome.</b>{' '}
              Copy number needs intensity, which whole-genome-amplified material cannot support
              here, and the allelic route was measured and refused: the noise runs through the
              signal, and the two mechanisms that are both real events separate by half a percent.
              See <SecRef id="segments" /> for the measurements. Loss is reported; gain is not.
            </List.Item>
            <List.Item>
              <b>It does not phase.</b> It reports which parent contributed, not which of that
              parent&apos;s two chromosomes. Phase requires an informative relative, reads long
              enough to span the markers, or typed gametes.
            </List.Item>
            <List.Item>
              <b>It cannot detect what the array cannot see.</b> Balanced rearrangements, low-level
              mosaicism below the resolution of a single biopsy, and anything in a region the chip
              does not cover are invisible to it.
            </List.Item>
            <List.Item>
              <b>It assumes your labels.</b> Which file is the sperm donor is your assertion. With
              both parents it checks they are two different people; it can never check they are the
              right two.
            </List.Item>
            <List.Item>
              <b>Consanguinity degrades the second-parent axis.</b> Related parents share alleles,
              so a second parent supplies fewer alleles the first lacks, and that axis already only
              separates by 1.6-fold.
            </List.Item>
            <List.Item>
              <b>A single sample cannot reproduce embryo-clustered criteria.</b> Some published
              quality filters operate across all biopsies of one embryo. An unlabelled drop box
              cannot reconstruct that grouping, so those gates flag candidates rather than
              reproducing the published filter.
            </List.Item>
            <List.Item>
              <b>Contamination is not modelled.</b> A mixed sample can present as biparental. The
              heterozygosity-plausible gate catches gross cases; it is not a contamination assay.
            </List.Item>
            <List.Item>
              <b>No liftOver.</b> A build mismatch is reported, never corrected. See{' '}
              <SecRef id="assembly" />.
            </List.Item>
          </List>
          <Text>
            Confirm any result that matters by an independent method. The tool is built to make
            that confirmation easy to specify: every number it reports carries the reference it was
            measured against, so what would falsify a call is always stated.
          </Text>
        </Section>

        {/* --- 24 ------------------------------------------------------------------------ */}
        <Section id="privacy" title="Privacy and data handling">
          <Text mb={8}>
            Array files from an embryo experiment are among the most identifying artefacts in this
            workflow: they describe a family. The terms promise that nothing about a family is
            submitted or retained, and the architecture is what makes that promise checkable rather
            than merely stated.
          </Text>
          <List size="sm" spacing={6} mb={10}>
            <List.Item>
              <b>There is no endpoint to send them to.</b> Every file is read by the browser using
              the file streaming API. No genotype is transmitted, and no request carrying one
              exists in the code.
            </List.Item>
            <List.Item>
              <b>Nothing is stored.</b> Parent genotypes live in memory for the run and are dropped
              when you clear or reload the page. There is no local database and no cache.
            </List.Item>
            <List.Item>
              <b>The report is built in the tab.</b> The PDF is typeset in the browser rather than
              on a server, specifically so that producing it does not require uploading anything.
            </List.Item>
            <List.Item>
              <b>Hashes travel instead of files.</b> The report carries the SHA-256 of each input,
              which lets a reviewer confirm they hold the same file without a file ever moving.
            </List.Item>
          </List>
          <Text>
            You can verify all of this: open the network panel and run an analysis. The only
            requests are for the page itself and, if you use them, the bundled example files.
          </Text>
        </Section>

        {/* --- 25 ------------------------------------------------------------------------ */}
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

/* --- helpers, matching the panel docs ---------------------------------------------------- */

const Wide = ({ children }: { children: ReactNode }) => (
  <div style={{ overflowX: 'auto', marginBottom: 10 }}>{children}</div>
)

const ExRow = ({ k, v, mono }: { k: string; v: string; mono?: boolean }) => (
  <Table.Tr>
    <Table.Td
      className={mono ? 'om-mono' : undefined}
      style={{ width: 190, color: 'var(--om-text-dim)' }}
    >
      {k}
    </Table.Td>
    <Table.Td style={{ whiteSpace: 'normal' }}>{v}</Table.Td>
  </Table.Tr>
)

/** A row whose leading cell is a label and whose rest are figures, right-aligned and monospaced
 *  so a column of percentages can be read down. `plain` for rows whose cells are words. */
const NumRow = ({ c, plain }: { c: string[]; plain?: boolean }) => (
  <Table.Tr>
    {c.map((v, i) => (
      <Table.Td
        key={v + String(i)}
        className={i && !plain ? 'om-mono' : undefined}
        ta={i && !plain ? 'right' : undefined}
      >
        {v}
      </Table.Td>
    ))}
  </Table.Tr>
)

const int = (n: number): string => n.toLocaleString('en-US')

const ordinal = (n: number): string => {
  const s = ['th', 'st', 'nd', 'rd']
  const v = n % 100
  return `${n}${s[(v - 20) % 10] ?? s[v] ?? s[0]}`
}

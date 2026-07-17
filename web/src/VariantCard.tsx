import { useState } from 'react'
import { Alert, Anchor, Badge, Button, Checkbox, Group, Paper, Table, Text } from '@mantine/core'
import { geneMismatch, links, type NLResponse, type ResolveResponse } from './api'
import { coord, sci, strandLabel } from './fmt'

const Row = ({ k, children }: { k: string; children: React.ReactNode }) => (
  <Table.Tr>
    <Table.Td className="om-kv-k" style={{ color: 'var(--om-text-dim)', verticalAlign: 'top' }}>{k}</Table.Td>
    <Table.Td>{children}</Table.Td>
  </Table.Tr>
)

const sigColor = (s: string | null) => {
  const t = (s ?? '').toLowerCase()
  if (t.includes('pathogenic') && !t.includes('conflicting')) return 'red'
  if (t.includes('benign')) return 'green'
  return 'gray'
}

/**
 * The resolved variant, and the parse that produced it.
 *
 * `nl` is required, not optional: the caveat about a model-chosen variant has to be on
 * screen at the same time as the variant it doubts. Pass null when a human typed the
 * identifier.
 */
export function VariantCard({
  data, nl, onBuild, busy,
}: {
  data: ResolveResponse
  nl: NLResponse | null
  onBuild: () => void
  busy: boolean
}) {
  const { variant: v, rarity, transcript_sense, clinvar_url } = data
  const minus = v.strand === -1

  // The user named a gene and got a variant in another one. Not a warning to read past: the
  // panel would be correct about a locus nobody asked for, so Build is withheld until they
  // say they meant this. Scoped to this card, which is keyed on the variant upstream.
  const named = nl?.named_genes ?? []
  const wrongGene = geneMismatch(v.gene, named)
  const [meantIt, setMeantIt] = useState(false)

  return (
    <Paper mb="sm">
      <Group justify="space-between" className="om-section-title" wrap="nowrap">
        <span>Resolved variant</span>
      </Group>

      <div style={{ padding: 8 }}>
        {nl && (
          <Alert
            mb={8}
            role={nl.used_llm ? 'alert' : undefined}
            color={nl.used_llm ? 'red' : 'gray'}
            variant={nl.used_llm ? 'filled' : 'light'}
            title={nl.used_llm ? 'A language model chose this variant' : 'Parsed from your wording'}
          >
            {/* Two claims, kept apart on purpose. WHICH variant is the model's, recalled
                from its own memory rather than read off anything the user typed. The
                coordinate is not the model's at all. Merging them either accuses the
                position of being invented or lets the choice pass as looked-up. */}
            {nl.used_llm && (
              <Text size="sm" fw={600} mb={6}>
                Nothing you typed named a variant, so a model named one from its own memory. The
                position and alleles below were then looked up live from that identifier, not
                taken from the model: they are right about the variant the model named. Check
                that is the variant you meant.
              </Text>
            )}
            <Text size="xs" className="om-mono">
              variant={nl.query.variant}
              {nl.query.gene ? ` gene=${nl.query.gene}` : ''}
              {nl.query.window_bp ? ` window_bp=${nl.query.window_bp}` : ''}
              {nl.query.ancestry ? ` ancestry=${nl.query.ancestry}` : ''}
              {nl.query.common_maf != null ? ` common_maf=${nl.query.common_maf}` : ''}
            </Text>
            {nl.note && (
              <Text size="xs" mt={2} c={nl.used_llm ? undefined : 'dimmed'}>
                {nl.note}
              </Text>
            )}
          </Alert>
        )}

        {wrongGene && (
          <Alert color="red" variant="filled" mb={8} role="alert" title="That is not the gene you named">
            <Text size="sm" fw={600}>
              You asked about {named.join(' or ')}. That identifier is a variant in {v.gene}.
            </Text>
            <Text size="xs" mt={6}>
              The record below is what {v.query} resolves to, live. Read it before you decide.
              If you did mean a variant in {v.gene}, say so and Build turns on.
            </Text>
            <Checkbox
              mt={8}
              size="xs"
              checked={meantIt}
              onChange={(e) => setMeantIt(e.currentTarget.checked)}
              label={`I meant this variant in ${v.gene}, not ${named.join(' or ')}.`}
            />
          </Alert>
        )}

        {v.build_note && (
          <Alert color="yellow" mb={8} title="Build note" variant="light">
            <Text size="xs">{v.build_note}</Text>
          </Alert>
        )}

        {/* The title says what we can do, not what the variant is: population_LD_usable is
            false both when the allele is too rare and when the frequency lookups came back
            empty. Only `rarity.reason` separates those two, so it is rendered alone: no
            fallback chain, which would let another field mask it. */}
        {!rarity.population_LD_usable && (
          <Alert color="red" mb={8} title="Population LD is not usable for this variant" role="alert">
            <Text size="xs">{rarity.reason}</Text>
          </Alert>
        )}

        <Table className="om-table om-kv" withRowBorders={false}>
          <Table.Tbody>
            <Row k="Gene">
              <Text span fw={600}>{v.gene ?? '-'}</Text>{' '}
              <Text span size="xs" c="dimmed">
                strand {strandLabel(v.strand)}
              </Text>
            </Row>
            <Row k="Query">
              <span className="om-mono">{v.query}</span>
            </Row>
            <Row k="rsID">
              {v.rsid ? (
                <Anchor className="om-mono" href={links.dbsnp(v.rsid)} target="_blank" rel="noreferrer">
                  {v.rsid}
                </Anchor>
              ) : (
                <Text span size="xs" c="dimmed">not in dbSNP</Text>
              )}
            </Row>
            <Row k="Coordinate">
              <span className="om-mono">{coord(v.chrom, v.pos_grch38, 'GRCh38')}</span>
              {v.pos_grch37 != null && (
                <Text span size="xs" c="dimmed" className="om-mono">
                  {'  ·  GRCh37 '}chr{v.chrom.replace(/^chr/, '')}:{v.pos_grch37.toLocaleString('en-US')}
                </Text>
              )}
            </Row>

            <Row k="Genomic (VCF, +strand)">
              <span className="om-mono">
                {v.vcf_ref}&gt;{v.vcf_alt}
              </span>
            </Row>
            <Row k="Transcript sense (HGVS c.)">
              <span className="om-mono">{transcript_sense}</span>{' '}
              {minus && (
                <Badge size="xs" color="orange" variant="light" ml={4}>
                  minus-strand gene: complement of the genomic form
                </Badge>
              )}
            </Row>

            <Row k="Clinical significance">
              <Badge size="sm" variant="light" color={sigColor(v.clinical_significance)}>
                {v.clinical_significance ?? 'not provided'}
              </Badge>
              {v.review_status && (
                <Text span size="xs" c="dimmed" ml={6}>{v.review_status}</Text>
              )}
            </Row>
            <Row k="ClinVar">
              {v.clinvar_accession ? (
                <Anchor
                  className="om-mono"
                  href={clinvar_url || links.clinvar(v.clinvar_accession)}
                  target="_blank"
                  rel="noreferrer"
                >
                  {v.clinvar_accession}
                </Anchor>
              ) : (
                <Text span size="xs" c="dimmed">no accession</Text>
              )}
            </Row>

            <Row k="gnomAD (genomes)">
              <span className="om-mono">AF {sci(rarity.gnomad_af_genome)}</span>
              <Text span size="xs" c="dimmed" className="om-mono">
                {'  ·  AC '}{rarity.gnomad_ac_genome ?? '-'}
                {'  ·  AN '}{rarity.gnomad_an_genome?.toLocaleString('en-US') ?? '-'}
              </Text>
              {rarity.gnomad_callset === 'genome' && (
                <Text span size="xs" c="dimmed"> · rarity read from here</Text>
              )}
            </Row>
            {/* Shown only when the variant is in the exome callset, which for a coding
                pathogenic variant is where its real frequency lives. */}
            {(rarity.gnomad_an_exome ?? 0) > 0 && (
              <Row k="gnomAD (exomes)">
                <span className="om-mono">AF {sci(rarity.gnomad_af_exome)}</span>
                <Text span size="xs" c="dimmed" className="om-mono">
                  {'  ·  AC '}{rarity.gnomad_ac_exome ?? '-'}
                  {'  ·  AN '}{rarity.gnomad_an_exome?.toLocaleString('en-US') ?? '-'}
                </Text>
                {rarity.gnomad_callset === 'exome' && (
                  <Text span size="xs" c="dimmed"> · rarity read from here</Text>
                )}
              </Row>
            )}
            <Row k="1000 Genomes">
              <span className="om-mono">AC {rarity.thousand_genomes_ac ?? '-'}</span>
            </Row>
            <Row k="LD usability">
              <Badge size="sm" color={rarity.population_LD_usable ? 'green' : 'red'} variant="light">
                {rarity.population_LD_usable ? 'LD usable' : 'LD not usable'}
              </Badge>
              {/* The reason is already in the red banner above when LD is not usable. */}
              {rarity.population_LD_usable && (
                <Text size="xs" c="dimmed" mt={2}>{rarity.reason}</Text>
              )}
            </Row>
          </Table.Tbody>
        </Table>

        <Group mt={8} gap={8}>
          {/* Never the default action once the gene disagrees: withheld until acknowledged,
              and still not the primary button afterwards. */}
          <Button
            onClick={onBuild}
            loading={busy}
            disabled={wrongGene && !meantIt}
            variant={wrongGene ? 'default' : 'filled'}
          >
            Build panel
          </Button>
          <Text size="xs" c="dimmed">
            {wrongGene && !meantIt
              ? `Withheld: you named ${named.join(' or ')} and this variant is in ${v.gene}.`
              : 'Pulls every common SNP in the window from gnomAD. Typically 20 to 60 s.'}
          </Text>
        </Group>
      </div>
    </Paper>
  )
}

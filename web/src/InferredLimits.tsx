import { Text } from '@mantine/core'
import { pct } from './parentage'
import { MIN_PRODUCTS } from './inferredReference'

/**
 * What the tool will not report when the reference is reconstructed rather than measured, and
 * the provenance that has to travel with every number.
 *
 * These are not cautions. Each one is a measurement that came out badly, and the list exists so
 * a reader can see the boundary of the evidence without having to find the paper.
 */

export interface Limit {
  what: string
  why: string
}

/** Every entry is measured. Nothing here is precautionary. */
export const LIMITS: Limit[] = [
  {
    what: 'Whether a genome is androgenetic or biparental',
    why: 'that axis is derived from the parent\'s own heterozygosity, and a reconstruction is '
      + 'homozygous everywhere by construction. The expectation collapses to zero and every '
      + 'diploid sample would read biparental on no evidence.',
  },
  {
    what: `Any verdict from fewer than ${MIN_PRODUCTS} products`,
    why: 'at three products and fewer, every true offspring tested inverted to a decisive wrong '
      + 'answer rather than to a refusal: 24 of 24 across two experiments. Four was not measured '
      + `separately, so the floor sits at ${MIN_PRODUCTS} rather than at the last depth known to `
      + 'fail. The tool declines to build below it at all.',
  },
  {
    what: 'Sperm type, and any result on chromosome X',
    why: 'the reconstruction is built from autosomes only. Chromosome Y is read, but only to '
      + 'decide which group is the father\'s; no product is ever selected or dropped by it, '
      + 'which would leave no paternal X in the reference by construction.',
  },
  {
    what: 'Which side is which, when no product carries a Y',
    why: 'a maternal cell cannot carry a Y, so one product that does names its group the '
      + 'father\'s. Nothing names it the other way. A paternal group of n products is all '
      + 'X-bearing 2^-n of the time, which at five products is 3%, so the two sides are reported '
      + 'as "this parent" and "the other parent" rather than guessed at. The split between them '
      + 'is unaffected and stays exact.',
  },
  {
    what: 'A parental-origin call on an array the gates excluded',
    why: 'every array is called, including the ones set aside as too poor or not haploid, '
      + 'because a fused or half-failed zygote is what this gets asked about most. Those rows '
      + 'carry the ratio like any other and say on their face that they were excluded: the '
      + 'reference was not built from them and their own noise ceiling is what it is.',
  },
  {
    what: 'A decisive call on a diploid sample',
    why: 'discrimination against diploid genomes depends on how far the marker set has drifted '
      + 'from the genome, and it degraded measurably at settings that were otherwise fine. '
      + 'Haploid samples are unaffected.',
  },
  {
    what: 'Telling this parent from a close relative',
    why: 'a simulated man sharing half the genome landed inside the uncalled band and was never '
      + 'rejected. Nothing in an inferred reference separates them.',
  },
  {
    what: 'Comparing these rates against a run that used a real array',
    why: 'the denominator differs. A measured array contributes its heterozygous markers to it '
      + 'and a reconstruction contributes none, so the two are not the same quantity.',
  },
]

export function LimitsPanel() {
  return (
    <div style={{ marginTop: 10 }}>
      <Text size="xs" c="dimmed" mb={8} style={{ maxWidth: 780, lineHeight: 1.5 }}>
        Each was measured coming out wrong, not guessed at. A refusal means the evidence does not
        reach, which is a different claim from a negative result.
      </Text>
      <div style={{ display: 'grid', gap: 7 }}>
        {LIMITS.map((l) => (
          <div key={l.what} style={{ display: 'flex', gap: 9, alignItems: 'baseline' }}>
            <span style={{
              color: 'var(--om-higher)', fontFamily: 'var(--om-mono)', fontSize: 11,
              flex: 'none', lineHeight: 1.5,
            }}
            >
              no
            </span>
            <span style={{ fontSize: 11.5, lineHeight: 1.55 }}>
              <b style={{ fontWeight: 600 }}>{l.what}</b>
              <span style={{ color: 'var(--om-text-dim)' }}>, because {l.why}</span>
            </span>
          </div>
        ))}
      </div>
    </div>
  )
}

export interface Provenance {
  products: string[]
  mMin: number
  markers: number
  ascertainment: number
  contamination: number
  spuriousAbsence: number
  tool: string
  generatedAt: string
  reportId: string
}

/**
 * The stamp that makes a reconstructed reference impossible to mistake for a measured one.
 *
 * It has to survive one page being photocopied out of context, so it repeats on every export
 * rather than appearing once at the front.
 */
export function ProvenanceStamp({ p }: { p: Provenance }) {
  const line = (k: string, v: string) => (
    <div style={{ display: 'flex', gap: 6, fontSize: 10.5, lineHeight: 1.6 }}>
      <span style={{ color: 'var(--om-text-dim)', minWidth: 116 }}>{k}</span>
      <span style={{ fontFamily: 'var(--om-mono)', fontVariantNumeric: 'tabular-nums' }}>{v}</span>
    </div>
  )
  return (
    <div style={{
      marginTop: 14, borderLeft: '3px solid var(--om-higher)', background: 'var(--om-zebra)',
      padding: '10px 14px',
    }}
    >
      <Text size="xs" fw={700} mb={4} style={{ color: 'var(--om-higher)' }}>
        Reconstructed reference, not a measured array
      </Text>
      <Text size="xs" c="dimmed" mb={7} style={{ maxWidth: 720, lineHeight: 1.5 }}>
        No array of this parent was used. The genotype below was inferred from {p.products.length}{' '}
        haploid products, and every rate on this page was measured against that inference.
      </Text>
      {line('products', p.products.join(', '))}
      {line('agreement', `at least ${p.mMin} products called and agreeing`)}
      {line('markers', p.markers.toLocaleString())}
      {line('ascertainment', `${pct(p.ascertainment, 1)} of the parent's heterozygosity`)}
      {line('contamination', `${pct(p.contamination, 2)}, adding ${pct(p.spuriousAbsence, 2)}`)}
      {line('tool', p.tool)}
      {line('report', p.reportId)}
      {line('generated', p.generatedAt)}
    </div>
  )
}

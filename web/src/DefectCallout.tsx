/**
 * One chromosomal defect, with its parental origin as the headline and the evidence beneath it.
 *
 * The ordering is the point. A reader wants WHOSE COPY first and the coordinates second, so origin
 * is the headline and everything supporting it is a chip below. Where origin is not determinable
 * the headline says so in the same place and the same size, rather than the box quietly becoming a
 * coordinate list: a defect with no origin is a result, not a missing field.
 *
 * Every chip is evidence someone could check. Nothing here is decorative, and a chip is omitted
 * rather than shown empty, because a chip reading "—" invites the reader to supply a value.
 *
 * The `exclusive` chip is the load-bearing one on a single-parent run. It counts markers where the
 * sample carries an allele the loaded parent does not have, which is Mendelian exclusion rather
 * than a statistic: dropout removes alleles and cannot invent one, so the count is evidence a
 * reader can act on directly. Measured 5,130-5,172 when the loaded parent's copy was removed on a
 * real trio and exactly 0 when the other parent's was.
 */
import { Text } from '@mantine/core'
import type { Defect } from './defects.ts'

const mb = (n: number): string => (n / 1e6).toFixed(2)

const Chip = ({ label, value, mono = false }: {
  label: string, value: string, mono?: boolean
}) => (
  <span style={{
    display: 'inline-flex', alignItems: 'baseline', gap: 5,
    border: '1px solid var(--om-defect)', background: 'var(--om-defect-chip)',
    padding: '2px 8px', fontSize: 12, lineHeight: 1.5, whiteSpace: 'nowrap',
  }}
  >
    <span style={{ color: 'var(--om-defect)', fontWeight: 600 }}>{label}</span>
    <span style={{ fontFamily: mono ? 'var(--om-mono)' : undefined, fontVariantNumeric: 'tabular-nums' }}>
      {value}
    </span>
  </span>
)

/** The headline: whose copy, or that it could not be said. */
const headline = (d: Defect): string => {
  const what = d.kind === 'copy-gain' ? 'extra copy' : 'copy lost'
  if (d.origin === 'unclear') return `chr${d.chrom} ${what} — origin not determined`
  return `${d.origin.toUpperCase()} ${what} · chr${d.chrom}`
}

export function DefectCallout({ defects }: { defects: Defect[] }) {
  if (!defects.length) return null
  const named = defects.filter((d) => d.origin !== 'unclear').length
  return (
    <div style={{
      border: '1px solid var(--om-defect)', borderLeft: '4px solid var(--om-defect)',
      background: 'var(--om-defect-bg)', padding: '11px 14px', margin: '10px 0 4px',
    }}
    >
      <Text style={{ fontSize: 15, fontWeight: 700, color: 'var(--om-defect)', lineHeight: 1.25 }}>
        Chromosomal change: {defects.length} region{defects.length === 1 ? '' : 's'}
        {named ? `, ${named} with a parent named` : ', none with a parent named'}
      </Text>

      <div style={{ marginTop: 9, display: 'grid', gap: 11 }}>
        {defects.map((d) => (
          <div key={`${d.chrom}:${d.startBp}:${d.kind}`}>
            <Text style={{
              fontSize: 14, fontWeight: 700, lineHeight: 1.3,
              color: d.origin === 'unclear' ? 'var(--om-text-dim)' : 'var(--om-defect)',
            }}
            >
              {headline(d)}
            </Text>

            <div style={{ display: 'flex', flexWrap: 'wrap', gap: 5, marginTop: 5 }}>
              <Chip label="at" value={d.locus} mono />
              <Chip label="span" value={`${mb(d.endBp - d.startBp)} Mb`} mono />
              {d.interval ? <Chip label="edge" value={d.interval} mono /> : null}
              <Chip label="event" value={d.kind.replace('copy-', '').replace('-', ' ')} />
              {d.informative !== undefined
                ? <Chip label="markers" value={String(d.informative)} mono /> : null}
              {d.posterior !== undefined && Number.isFinite(d.posterior)
                ? <Chip label="posterior" value={d.posterior.toFixed(3)} mono /> : null}
              {d.exclusive !== undefined
                ? <Chip label="exclusive" value={String(d.exclusive)} mono /> : null}
              {d.basis ? <Chip label="basis" value={d.basis.replace('-', ' ')} /> : null}
              {d.stage ? <Chip label="material" value={d.stage} /> : null}
              {d.dropout !== undefined
                ? <Chip label="dropout" value={d.dropout.toFixed(3)} mono /> : null}
            </div>

            <Text size="sm" mt={5} style={{ maxWidth: 780, lineHeight: 1.5 }}>
              {d.why}
            </Text>
          </div>
        ))}
      </div>

      <Text size="xs" c="dimmed" mt={9} style={{ maxWidth: 780, lineHeight: 1.5 }}>
        Position comes from the intensity channel and needs no parental genotype, so it is reported
        for every region. Naming the parent needs allele dosage, which needs a genotyped parent;
        where the run has none the origin reads not determined and the coordinates stand on their
        own rather than being withheld.
      </Text>
    </div>
  )
}

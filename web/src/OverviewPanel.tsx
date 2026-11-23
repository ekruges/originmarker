/**
 * The run in one flat table (overview.ts), with its download. Plain rows, no chart.
 */
import { Button, Table, Text } from '@mantine/core'
import { overviewTable, overviewTsv, type Arms, type OverviewSample } from './overview.ts'

export function OverviewPanel({ samples, arms }: { samples: readonly OverviewSample[]; arms: Arms | null }) {
  if (!samples.length) return null
  const table = overviewTable(samples, arms ?? new Map())
  const download = () => {
    const blob = new Blob([overviewTsv(table)], { type: 'text/tab-separated-values' })
    const a = document.createElement('a')
    a.href = URL.createObjectURL(blob)
    a.download = `syngamy-overview-${new Date().toISOString().slice(0, 10)}.tsv`
    a.click()
    URL.revokeObjectURL(a.href)
  }
  return (
    <div style={{ border: '1px solid var(--om-line)', padding: '12px 15px', margin: '14px 0 4px' }}>
      <Text style={{ fontSize: 11, fontWeight: 700, letterSpacing: '0.06em', color: 'var(--om-text-dim)' }}>
        THE RUN, COUNTED
      </Text>
      <Text size="xs" c="dimmed" mt={4} style={{ maxWidth: 820, lineHeight: 1.5 }}>
        Changes are the result cards' rows with a gain or a loss; a row with neither is counted as a
        finding. A sample not reportable is counted here and adds no change. An arm runs from one end of an autosome to its centromere (the p/q boundary of the UCSC hg19 band table, the one the whole-chromosome test reads), 5 Mb either way.
        {arms ? '' : ' The chromosome lengths did not load, so every change that is not whole is counted as a piece.'}
      </Text>
      <Table className="om-table" withTableBorder mt={8} style={{ maxWidth: 520 }}>
        <Table.Tbody>
          {table.map(([label, v]) => (
            <Table.Tr key={label}>
              <Table.Td style={{ fontSize: 12 }}>{label}</Table.Td>
              <Table.Td style={{ fontFamily: 'var(--om-mono)', fontSize: 12, textAlign: 'right' }}>{v}</Table.Td>
            </Table.Tr>
          ))}
        </Table.Tbody>
      </Table>
      <Button variant="default" size="xs" mt={8} onClick={download} title="The same table, tab separated">
        Overview (TSV)
      </Button>
    </div>
  )
}

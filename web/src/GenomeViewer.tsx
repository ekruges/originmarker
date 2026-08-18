/**
 * The genome map as something a reader can interrogate, not just look at.
 *
 * WHY A VIEWER AND NOT A PICTURE. The static figure answers "where did the changes land". The next
 * question is always "what is THAT one", and a static figure sends a reader back to a table to find
 * out. Here a region is clickable and answers for itself: its coordinates, every site class it
 * overlaps, and whether any of those classes offers an alternative explanation for it.
 *
 * THE STAR IS NARROW ON PURPOSE. It marks a region that sits in a class whose coincidence with THIS
 * set of regions is more than the matched null produces, not one that merely touches a feature.
 * Late-replication valleys cover most of the genome, so flagging on contact would star nearly
 * everything and mean nothing.
 *
 * The report keeps the static figure. This is the screen's version of the same geometry and the
 * same colours, with the interaction a printed page cannot carry.
 */
import { useMemo, useState } from 'react'
import { Text } from '@mantine/core'
import { regionFlags, scannedTracks, type ComparisonResult } from './comparison.ts'
import { LANE_COLOUR } from './comparisonFigures.ts'
import { FIG } from './figures.ts'

const SHORT: Record<string, string> = {
  'common fragile site': 'fragile site',
  'gene over 500 kb': 'long gene',
  'centromere or telomere': 'centromere/telomere',
  'late-replication valley (ES)': 'late-repl (ES)',
  'late-replication valley (constitutive)': 'late-repl (const.)',
}

export function GenomeViewer({ c }: { c: ComparisonResult }) {
  const [sel, setSel] = useState<number | null>(null)
  const [hidden, setHidden] = useState<Set<string>>(new Set())

  const model = useMemo(() => {
    const tracks = scannedTracks(c.track)
    const lanes = Object.keys(tracks).sort()
    const ends = new Map<string, number>()
    const centro = new Map<string, number>()
    for (const g of c.track.gaps ?? []) {
      ends.set(g.chrom, Math.max(ends.get(g.chrom) ?? 0, g.endBp))
      if ((g as { kind?: string }).kind === 'centromere') {
        centro.set(g.chrom, (g.startBp + g.endBp) / 2)
      }
    }
    const order = [...ends.keys()].filter((k) => /^\d+$/.test(k))
      .sort((a, b) => Number(a) - Number(b))
    return { tracks, lanes, ends, centro, order, max: Math.max(1, ...order.map((k) => ends.get(k) ?? 0)) }
  }, [c])

  const flags = useMemo(() => regionFlags(c), [c])
  const starred = flags.filter((f) => f.related.length).length

  const labelW = 26
  const plotW = 560
  const laneH = 3
  const barH = 7
  const shown = model.lanes.filter((l) => !hidden.has(l))
  const blockH = barH + shown.length * laneH + 9
  const H = 12 + model.order.length * blockH + 6
  const sx = (bp: number) => labelW + (bp / model.max) * plotW

  const chosen = sel === null ? null : flags[sel]

  return (
    <div>
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 10, alignItems: 'center',
        marginBottom: 8 }}
      >
        {model.lanes.map((l) => (
          <button
            key={l}
            type="button"
            onClick={() => setHidden((p) => {
              const n = new Set(p); if (n.has(l)) n.delete(l); else n.add(l); return n
            })}
            style={{
              display: 'flex', alignItems: 'center', gap: 5, background: 'none', cursor: 'pointer',
              border: 'none', padding: 0, fontSize: 11,
              color: hidden.has(l) ? 'var(--om-text-dim)' : 'var(--om-text)',
              opacity: hidden.has(l) ? 0.45 : 1,
            }}
          >
            <span style={{ width: 9, height: 7, background: LANE_COLOUR[l] ?? FIG.line,
              display: 'inline-block' }}
            />
            {SHORT[l] ?? l}
          </button>
        ))}
        <span style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 11 }}>
          <span style={{ width: 9, height: 7, background: FIG.accent, display: 'inline-block' }} />
          detected region
        </span>
        <span style={{ fontSize: 11, color: 'var(--om-text-dim)' }}>
          {`* ${starred} of ${flags.length} have an alternative explanation available`}
        </span>
      </div>

      <svg viewBox={`0 0 ${labelW + plotW + 8} ${H}`} width="100%" role="img"
        aria-label="Detected regions on every autosome, with the site classes scanned beneath each"
      >
        {model.order.map((chrom, i) => {
          const t = 12 + i * blockH
          const end = model.ends.get(chrom) ?? 0
          const cen = model.centro.get(chrom)
          return (
            <g key={chrom}>
              <text x={labelW - 4} y={t + barH - 0.5} fontSize={7} fill={FIG.ink} textAnchor="end">
                {chrom}
              </text>
              <rect x={sx(0)} y={t} width={sx(end) - sx(0)} height={barH} fill="#eceff2" />
              {cen !== undefined && (
                <circle cx={sx(cen)} cy={t + barH / 2} r={2} fill="none" stroke={FIG.axis}
                  strokeWidth={0.8}
                />
              )}
              {shown.map((lane, li) => (
                <g key={lane}>
                  <rect x={sx(0)} y={t + barH + 2 + li * laneH} width={sx(end) - sx(0)}
                    height={laneH - 0.8} fill="#f4f6f8"
                  />
                  {model.tracks[lane].filter((f) => f.chrom === chrom).map((f, k) => (
                    <rect
                      key={k} x={sx(f.startBp)} y={t + barH + 2 + li * laneH}
                      width={Math.max(0.6, sx(f.endBp) - sx(f.startBp))} height={laneH - 0.8}
                      fill={LANE_COLOUR[lane] ?? FIG.line}
                    />
                  ))}
                </g>
              ))}
              {flags.filter((f) => f.chrom === chrom).map((f) => (
                <g key={f.index} style={{ cursor: 'pointer' }}
                  onClick={() => setSel(sel === f.index ? null : f.index)}
                >
                  <rect
                    x={sx(f.startBp)} y={t - 1}
                    width={Math.max(1.6, sx(f.endBp) - sx(f.startBp))} height={barH + 2}
                    fill={FIG.accent}
                    stroke={sel === f.index ? FIG.ink : 'none'} strokeWidth={sel === f.index ? 1 : 0}
                  />
                  {f.related.length > 0 && (
                    <text x={(sx(f.startBp) + sx(f.endBp)) / 2} y={t - 2.5} fontSize={8}
                      fill={FIG.accent} textAnchor="middle" fontWeight={700}
                    >
                      *
                    </text>
                  )}
                  <title>{`${f.name}${f.related.length ? '  *' : ''}`}</title>
                </g>
              ))}
            </g>
          )
        })}
      </svg>

      {chosen ? (
        <div style={{ border: '1px solid var(--om-line)', borderLeft: `3px solid ${FIG.accent}`,
          padding: '8px 12px', marginTop: 8 }}
        >
          <Text style={{ fontSize: 12, fontWeight: 700 }}>
            {chosen.name}
            {chosen.related.length ? '  *' : ''}
          </Text>
          <Text style={{ fontSize: 11, color: 'var(--om-text-dim)', lineHeight: 1.5, marginTop: 3 }}>
            {`${((chosen.endBp - chosen.startBp) / 1e6).toFixed(2)} Mb. `}
            {chosen.overlaps.length
              ? `Overlaps ${chosen.overlaps.map((o) => (SHORT[o] ?? o)).join(', ')}.`
              : 'Overlaps none of the site classes scanned.'}
          </Text>
          <Text style={{ fontSize: 11, lineHeight: 1.5, marginTop: 5 }}>
            {chosen.related.length
              ? `An alternative explanation is available: this region sits in `
                + `${chosen.related.join(' and ')}, whose coincidence with this set of regions is `
                + 'more than the matched null produces. That is a reason to weigh it more '
                + 'cautiously, not a reason to discard it.'
              : 'No site class that cleared the corrected threshold covers this region, so nothing '
                + 'here offers an alternative explanation for it.'}
          </Text>
        </div>
      ) : (
        <Text style={{ fontSize: 11, color: 'var(--om-text-dim)', marginTop: 6 }}>
          Click a region to see what it overlaps. Click a key entry to hide that lane.
        </Text>
      )}
    </div>
  )
}

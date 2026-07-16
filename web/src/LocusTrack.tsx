import { useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react'
import { Button, Checkbox, Group, Paper, Slider, Text, Tooltip } from '@mantine/core'
import { isUpper, TIERS, TIER_LABEL, type Marker, type PanelResult, type Tier } from './api'
import { axisLabel, int, num, shortCoord, sig2, signedBp } from './fmt'
import { ancestryHet, shownHet } from './rank'

const H = 250
// top has to clear the variant's rsID label, drawn ABOVE the plot area at PAD.top - LABEL_DY.
const PAD = { top: 28, right: 12, bottom: 34, left: 46 }
const LABEL_DY = 11
const LOWER = 'var(--om-lower)'
const UPPER = 'var(--om-higher)'
const MAX_HET = 0.5 // 2pq is bounded above by 0.5 at p = 0.5

interface Props {
  result: PanelResult
  ancestry: string | null
}

export function LocusTrack({ result, ancestry }: Props) {
  const { variant, candidates, recommended, provenance } = result
  const center = variant.pos_grch38
  const floor = provenance.common_maf ?? 0.05

  const [tiers, setTiers] = useState<Tier[]>([...TIERS])
  const [minMaf, setMinMaf] = useState(floor)
  const [domain, setDomain] = useState<[number, number] | null>(null)
  const [hover, setHover] = useState<{ m: Marker; x: number; y: number } | null>(null)
  const [brush, setBrush] = useState<{ from: number; to: number } | null>(null)

  const wrapRef = useRef<HTMLDivElement>(null)
  const [w, setW] = useState(900)

  useLayoutEffect(() => {
    const el = wrapRef.current
    if (!el) return
    const ro = new ResizeObserver(([e]) => setW(Math.max(360, e.contentRect.width)))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  // Reset zoom/filters when a different panel loads.
  useEffect(() => {
    setDomain(null)
    setMinMaf(floor)
    setTiers([...TIERS])
  }, [variant.pos_grch38, floor])

  const recSet = useMemo(() => new Set(recommended.map((m) => m.variant_id)), [recommended])

  const full = useMemo<[number, number]>(() => {
    const win = provenance.window_bp ?? 250_000
    let lo = center - win
    let hi = center + win
    for (const m of candidates) {
      lo = Math.min(lo, m.pos)
      hi = Math.max(hi, m.pos)
    }
    return [lo, hi]
  }, [candidates, center, provenance.window_bp])

  const [x0, x1] = domain ?? full
  const span = x1 - x0
  const plotW = Math.max(10, w - PAD.left - PAD.right)
  const plotH = H - PAD.top - PAD.bottom

  const variantLabel = variant.rsid ?? 'variant'
  // ~5px per glyph at 9px in this face. Only has to be close enough to keep the label off
  // the edges; erring wide is harmless.
  const labelHalfW = (variantLabel.length * 5) / 2

  const sx = (pos: number) => PAD.left + ((pos - x0) / span) * plotW
  const sy = (het: number) => PAD.top + plotH - (Math.min(het, MAX_HET) / MAX_HET) * plotH
  const invert = (px: number) => x0 + ((px - PAD.left) / plotW) * span

  // Only ever called for markers in `plotted`, which is filtered on shownHet being non-null.
  const yOf = (m: Marker) => sy(shownHet(m, ancestry) as number)

  const shown = useMemo(
    () => candidates.filter((m) => tiers.includes(m.tier) && m.maf >= minMaf),
    [candidates, tiers, minMaf],
  )
  // A marker with no frequency in the selected population has no height on this axis:
  // its global 2pq belongs to another population, and the baseline would draw "no data"
  // as zero. It leaves the plot instead, and the count above says so.
  const plotted = useMemo(
    () => shown.filter((m) => shownHet(m, ancestry) != null),
    [shown, ancestry],
  )
  const unpriced = shown.length - plotted.length
  const visible = useMemo(
    () => plotted.filter((m) => m.pos >= x0 && m.pos <= x1).sort((a, b) => a.pos - b.pos),
    [plotted, x0, x1],
  )

  const ticks = useMemo(() => {
    const step = Math.pow(10, Math.floor(Math.log10(span / 5)))
    const s = [1, 2, 5, 10].map((k) => k * step).find((k) => span / k <= 6) ?? step
    const out: number[] = []
    for (let t = Math.ceil(x0 / s) * s; t <= x1; t += s) out.push(t)
    return out
  }, [x0, x1, span])

  // Nearest-point hover: binary search over position-sorted visible markers.
  const onMove = (e: React.MouseEvent<SVGSVGElement>) => {
    const rect = e.currentTarget.getBoundingClientRect()
    const px = ((e.clientX - rect.left) / rect.width) * w
    const py = ((e.clientY - rect.top) / rect.height) * H
    if (brush) {
      setBrush({ ...brush, to: Math.min(Math.max(px, PAD.left), PAD.left + plotW) })
      return
    }
    if (!visible.length) return setHover(null)
    const target = invert(px)
    let lo = 0
    let hi = visible.length - 1
    while (lo < hi) {
      const mid = (lo + hi) >> 1
      if (visible[mid].pos < target) lo = mid + 1
      else hi = mid
    }
    const cands = [visible[lo - 1], visible[lo], visible[lo + 1]].filter(Boolean)
    let best: Marker | null = null
    let bestD = Infinity
    for (const m of cands) {
      const d = Math.hypot(sx(m.pos) - px, yOf(m) - py)
      if (d < bestD) { bestD = d; best = m }
    }
    setHover(best && bestD < 22 ? { m: best, x: sx(best.pos), y: yOf(best) } : null)
  }

  const startBrush = (e: React.MouseEvent<SVGSVGElement>) => {
    const rect = e.currentTarget.getBoundingClientRect()
    const px = ((e.clientX - rect.left) / rect.width) * w
    if (px < PAD.left || px > PAD.left + plotW) return
    setHover(null)
    setBrush({ from: px, to: px })
  }

  const endBrush = () => {
    if (!brush) return
    const [a, b] = [brush.from, brush.to].sort((p, q) => p - q)
    // Ignore stray clicks; only a real drag zooms.
    if (b - a > 6) {
      const lo = invert(a)
      const hi = invert(b)
      if (hi - lo > 50) setDomain([lo, hi]) // don't zoom past ~50 bp
    }
    setBrush(null)
  }

  const zoomed = domain != null
  const hetLabel = ancestry ? `${ancestry} 2pq (prior)` : 'global 2pq (prior)'

  return (
    <Paper mb="sm">
      <Group justify="space-between" className="om-section-title" wrap="nowrap">
        <span>
          Locus track: chr{variant.chrom.replace(/^chr/, '')} (GRCh38)
        </span>
        <Text size="xs" c="dimmed">
          {int(visible.length)} of {int(candidates.length)} candidates shown
          {unpriced > 0 ? ` · ${int(unpriced)} not plotted: no gnomAD ${ancestry} frequency` : ''} ·
          drag to zoom
        </Text>
      </Group>

      <Group gap={14} p={8} align="center" wrap="wrap">
        <Group gap={8}>
          <Text size="xs" c="dimmed">Tier</Text>
          {TIERS.map((t) => (
            <Checkbox
              key={t}
              size="xs"
              label={TIER_LABEL[t]}
              checked={tiers.includes(t)}
              onChange={(e) =>
                setTiers(e.currentTarget.checked ? [...tiers, t] : tiers.filter((x) => x !== t))
              }
            />
          ))}
        </Group>

        <Group gap={6} style={{ minWidth: 240, flex: '1 1 240px' }}>
          <Text size="xs" c="dimmed" style={{ whiteSpace: 'nowrap' }}>
            MAF ≥ {minMaf.toFixed(2)}
          </Text>
          <Slider
            style={{ flex: 1 }}
            size="xs"
            min={floor}
            max={0.5}
            step={0.01}
            value={minMaf}
            onChange={setMinMaf}
            label={(v) => v.toFixed(2)}
            aria-label="Minimum minor allele frequency"
          />
        </Group>

        <Button variant="default" size="xs" disabled={!zoomed} onClick={() => setDomain(null)}>
          Reset zoom
        </Button>
      </Group>

      <div ref={wrapRef} style={{ position: 'relative', padding: '0 8px 8px' }}>
        <svg
          className="om-track"
          viewBox={`0 0 ${w} ${H}`}
          height={H}
          role="img"
          aria-label={`Locus track: ${visible.length} candidate SNPs around ${variant.gene ?? 'the variant'} at ${shortCoord(variant.chrom, center)} on GRCh38. Stem height is the ${ancestry ? `expected heterozygosity prior in ${ancestry}` : 'global expected heterozygosity prior'}. The ranked table below carries the same data in accessible form.`}
          onMouseMove={onMove}
          onMouseLeave={() => { setHover(null); setBrush(null) }}
          onMouseDown={startBrush}
          onMouseUp={endBrush}
        >
          {/* y grid + axis */}
          {[0, 0.125, 0.25, 0.375, 0.5].map((h) => (
            <g key={h}>
              <line x1={PAD.left} x2={w - PAD.right} y1={sy(h)} y2={sy(h)} stroke="#eceef0" />
              <text x={PAD.left - 6} y={sy(h) + 3} textAnchor="end" fontSize={9} fill="#8d959e">
                {h.toFixed(2)}
              </text>
            </g>
          ))}
          <text
            transform={`translate(11 ${PAD.top + plotH / 2}) rotate(-90)`}
            textAnchor="middle"
            fontSize={9}
            fill="#6b727b"
          >
            {hetLabel}
          </text>

          {/* x axis */}
          <line
            x1={PAD.left}
            x2={w - PAD.right}
            y1={PAD.top + plotH}
            y2={PAD.top + plotH}
            stroke="#adb4bc"
          />
          {ticks.map((t) => (
            <g key={t}>
              <line x1={sx(t)} x2={sx(t)} y1={PAD.top + plotH} y2={PAD.top + plotH + 4} stroke="#adb4bc" />
              <text x={sx(t)} y={PAD.top + plotH + 14} textAnchor="middle" fontSize={9} fill="#6b727b">
                {axisLabel(t, span)}
              </text>
            </g>
          ))}
          <text
            x={PAD.left + plotW / 2}
            y={H - 4}
            textAnchor="middle"
            fontSize={9}
            fill="#8d959e"
          >
            chr{variant.chrom.replace(/^chr/, '')} position (GRCh38) · lower ◀ variant ▶ higher
          </text>

          {/* lollipops: stem height = the 2pq prior the axis names, colour = side */}
          {visible.map((m) => {
            const rec = recSet.has(m.variant_id)
            const c = isUpper(m) ? UPPER : LOWER
            const X = sx(m.pos)
            const Y = yOf(m)
            return (
              <g key={m.variant_id} opacity={rec ? 1 : 0.45}>
                <line x1={X} x2={X} y1={PAD.top + plotH} y2={Y} stroke={c} strokeWidth={rec ? 1.2 : 0.6} />
                <circle
                  cx={X}
                  cy={Y}
                  r={rec ? 3.4 : 1.7}
                  fill={rec ? c : '#fff'}
                  stroke={c}
                  strokeWidth={rec ? 1 : 0.8}
                />
              </g>
            )
          })}

          {/* the pathogenic variant */}
          {center >= x0 && center <= x1 && (
            <>
              <line
                x1={sx(center)}
                x2={sx(center)}
                y1={PAD.top - 8}
                y2={PAD.top + plotH}
                stroke="#c1272d"
                strokeWidth={1}
                strokeDasharray="3 2"
              />
              {/* Clamped inside the plot: brushed against an edge, a centred label would
                  hang off the side. */}
              <text
                x={Math.min(
                  Math.max(sx(center), PAD.left + labelHalfW),
                  PAD.left + plotW - labelHalfW,
                )}
                y={PAD.top - LABEL_DY}
                textAnchor="middle"
                fontSize={9}
                fill="#c1272d"
              >
                {variantLabel}
              </text>
            </>
          )}

          {hover && (
            <circle cx={hover.x} cy={hover.y} r={5.5} fill="none" stroke="#1d2126" strokeWidth={1} />
          )}

          {brush && (
            <rect
              x={Math.min(brush.from, brush.to)}
              y={PAD.top}
              width={Math.abs(brush.to - brush.from)}
              height={plotH}
              fill="#337ab7"
              opacity={0.12}
              stroke="#337ab7"
              strokeWidth={0.5}
            />
          )}
        </svg>

        {hover && (
          <div
            style={{
              position: 'absolute',
              left: Math.min(Math.max((hover.x / w) * 100, 2), 78) + '%',
              top: hover.y + 6,
              background: '#fff',
              border: '1px solid var(--om-border-strong)',
              padding: '4px 6px',
              fontSize: 11,
              pointerEvents: 'none',
              zIndex: 3,
              whiteSpace: 'nowrap',
            }}
          >
            <div className="om-mono" style={{ fontWeight: 600 }}>{hover.m.rsid}</div>
            <div className="om-mono">{shortCoord(hover.m.chrom, hover.m.pos)} (GRCh38)</div>
            <div className="om-mono">
              {signedBp(hover.m.dist)} bp · {isUpper(hover.m) ? 'higher' : 'lower'} · {hover.m.tier}
            </div>
            <div className="om-mono">MAF {num(hover.m.maf)}</div>
            <div className="om-mono">
              2pq {num(hover.m.het)} <span style={{ color: 'var(--om-text-dim)' }}>(prior)</span>
            </div>
            {ancestry && (
              <div className="om-mono">
                {ancestry} 2pq {num(ancestryHet(hover.m, ancestry))}{' '}
                <span style={{ color: 'var(--om-text-dim)' }}>(prior)</span>
              </div>
            )}
            <div className="om-mono">
              {hover.m.cm == null ? 'cM -' : `${sig2(hover.m.cm)} cM`}
              {hover.m.map_approx ? ' (approx.)' : ''}
            </div>
            {recSet.has(hover.m.variant_id) && (
              <div style={{ color: 'var(--om-blue)' }}>shortlisted</div>
            )}
          </div>
        )}

        {!visible.length && (
          <Text size="xs" c="dimmed" ta="center" py="md">
            {!plotted.length && unpriced > 0
              ? `No candidate here has a gnomAD ${ancestry} frequency, so none can be drawn against this axis. Clear the ancestry to see the global prior.`
              : 'No candidates match these filters. Lower the MAF threshold or re-enable a tier.'}
          </Text>
        )}

        <Group gap={14} mt={4} justify="center">
          <LegendDot color={LOWER} label="lower coordinate" />
          <LegendDot color={UPPER} label="higher coordinate" />
          <Tooltip
            withArrow
            multiline
            w={260}
            label="Filled = shortlisted: the top-scoring markers in each distance band on each side. Hollow = the rest of the candidate pool."
          >
            <Text size="xs" c="dimmed" style={{ cursor: 'help' }}>
              ● shortlisted / ○ other candidate
            </Text>
          </Tooltip>
        </Group>
      </div>
    </Paper>
  )
}

const LegendDot = ({ color, label }: { color: string; label: string }) => (
  <Group gap={4}>
    <svg width={10} height={10} aria-hidden="true">
      <circle cx={5} cy={5} r={3.4} fill={color} />
    </svg>
    <Text size="xs" c="dimmed">{label}</Text>
  </Group>
)

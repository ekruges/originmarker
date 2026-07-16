import { useEffect, useLayoutEffect, useRef, useState } from 'react'
import {
  ActionIcon, Alert, Autocomplete, Button, Checkbox, Group, NumberInput,
  Paper, SegmentedControl, Select, Text, TextInput, Tooltip,
} from '@mantine/core'
import { ANCESTRIES, api, ApiError, type Health, type NLResponse, type StructuredQuery } from './api'

const EXAMPLE = 'NM_000352.6(ABCC8):c.3989-9G>A'

/** Every accepted input form, one per rotation. The last is only offered when free text is
 *  actually enabled: advertising a form the server will refuse is worse than not showing it. */
const PLACEHOLDERS: { text: string; nl?: boolean }[] = [
  { text: EXAMPLE },
  { text: 'rs334' },
  { text: 'NM_000518.5(HBB):c.20A>T' },
  { text: 'VCV000009088' },
  { text: 'NC_000011.10:g.17397055C>T' },
  { text: 'the ABCC8 splice variant, in Europeans', nl: true },
]

const ROTATE_MS = 3600
const FADE_MS = 320

type Mode = 'search' | 'manual'

interface Props {
  health: Health | null
  busy: boolean
  /** `nl` is the parse provenance, passed only when a model produced this query. Handed
   *  up rather than shown here: the caveat has to outlive this box. */
  onResolve: (q: StructuredQuery, nl?: NLResponse) => void
  /** Landing view: one oversized box, chrome stripped, manual entry behind a menu. */
  hero?: boolean
}

/** Detects an identifier we can hand to /api/resolve verbatim; anything else is free text. */
const looksLikeId = (s: string) => /^rs\d+$/i.test(s.trim()) || /:c\.|:g\.|:n\.|>/.test(s)

export function SearchPanel({ health, busy, onResolve, hero = false }: Props) {
  const [mode, setMode] = useState<Mode>('search')
  const [text, setText] = useState('')
  const [nlError, setNlError] = useState<string | null>(null)
  const [parsing, setParsing] = useState(false)

  // Manual mode: every StructuredQuery knob, no parsing in the way.
  const [q, setQ] = useState<Required<Omit<StructuredQuery, 'gene' | 'ancestry'>> & {
    gene: string
    ancestry: string | null
  }>({
    variant: '', gene: '', window_bp: 250_000, build: 'GRCh38',
    ancestry: null, common_maf: 0.05, cross_check: true,
  })

  const [geneOpts, setGeneOpts] = useState<string[]>([])
  const [geneErr, setGeneErr] = useState(false)

  // Debounced gene autocomplete. A failure just stops suggesting.
  useEffect(() => {
    const term = q.gene.trim()
    if (term.length < 2) { setGeneOpts([]); return }
    const t = setTimeout(async () => {
      try {
        const gs = await api.genes(term)
        setGeneOpts(gs.map((g) => g.symbol))
        setGeneErr(false)
      } catch {
        setGeneOpts([])
        setGeneErr(true)
      }
    }, 250)
    return () => clearTimeout(t)
  }, [q.gene])

  const submitSearch = async () => {
    const t = text.trim()
    if (!t) return
    setNlError(null)
    if (looksLikeId(t)) {
      onResolve({ variant: t, build: 'GRCh38' })
      return
    }
    if (!health?.nl_enabled) {
      setNlError(
        'Free-text parsing is unavailable. Enter an rsID or HGVS string, or use Manual input.',
      )
      return
    }
    setParsing(true)
    try {
      const r = await api.nl(t)
      onResolve(r.query, r)
    } catch (e) {
      setNlError(e instanceof ApiError ? e.message : 'Could not parse that request.')
    } finally {
      setParsing(false)
    }
  }

  const submitManual = () => {
    if (!q.variant.trim()) return
    onResolve({
      variant: q.variant.trim(),
      gene: q.gene.trim() || null,
      window_bp: q.window_bp,
      build: q.build,
      ancestry: q.ancestry,
      common_maf: q.common_maf,
      cross_check: q.cross_check,
    })
  }

  const inputRef = useRef<HTMLInputElement>(null)
  useEffect(() => { inputRef.current?.focus() }, [mode])

  // Ghost examples cycle through the accepted input forms, the way a search box does. Only
  // while the box is empty: a placeholder moving under text someone is composing is noise,
  // and the element is aria-hidden from the input's own label anyway.
  const forms = PLACEHOLDERS.filter((p) => !p.nl || health?.nl_enabled)
  const [slot, setSlot] = useState(0)
  const [fading, setFading] = useState(false)
  const idle = mode === 'search' && !text && !busy
  const showGhost = idle

  // Copy the input's own text metrics onto the ghost rather than restate them in CSS: the
  // input's size comes from Mantine's theme and a hero override, so any constant here is a
  // guess that drifts. The ghost must sit exactly where the real text will appear, or the
  // swap to typed text jumps.
  const ghostRef = useRef<HTMLSpanElement>(null)
  useLayoutEffect(() => {
    const el = inputRef.current
    const g = ghostRef.current
    if (!el || !g) return
    const cs = getComputedStyle(el)
    g.style.fontSize = cs.fontSize
    g.style.fontFamily = cs.fontFamily
    g.style.letterSpacing = cs.letterSpacing
    g.style.paddingLeft = cs.paddingLeft
    g.style.paddingRight = cs.paddingRight
  }, [showGhost, hero, slot])

  // Not gated on prefers-reduced-motion. That setting exists for vestibular triggers,
  // which are changes to perceived size, shape or position; a cross-fade changes only
  // opacity and is what the guidance recommends motion be REPLACED with.
  useEffect(() => {
    if (!idle || forms.length < 2) return
    let swap: ReturnType<typeof setTimeout>
    const tick = setInterval(() => {
      setFading(true)
      swap = setTimeout(() => {
        setSlot((i) => (i + 1) % forms.length)
        setFading(false)
      }, FADE_MS)
    }, ROTATE_MS)
    return () => {
      clearInterval(tick)
      clearTimeout(swap)     // or a half-finished swap lands after the box is gone
    }
  }, [idle, forms.length])

  // Wrap inline rather than via a component defined in this body. A component declared
  // here gets a fresh function identity on every render, and React compares element types
  // by reference: a new identity is a new type, so the whole subtree unmounts and remounts
  // on each keystroke, destroying the input and its focus.
  const body = (
    <>
      {hero ? (
        <Group justify="flex-end" gap={4} mb={4}>
          <Select
            size="xs"
            variant="unstyled"
            aria-label="Query entry mode"
            value={mode}
            onChange={(v) => v && setMode(v as Mode)}
            data={[
              { label: 'Search', value: 'search' },
              { label: 'Manual input', value: 'manual' },
            ]}
            allowDeselect={false}
            styles={{ input: { fontSize: 12, color: 'var(--om-text-dim)', textAlign: 'right' } }}
            w={110}
          />
        </Group>
      ) : (
        <Group justify="space-between" className="om-section-title" wrap="nowrap">
          <span>Query</span>
          <SegmentedControl
            size="xs"
            value={mode}
            onChange={(v) => setMode(v as Mode)}
            data={[
              { label: 'Search', value: 'search' },
              { label: 'Manual input', value: 'manual' },
            ]}
            aria-label="Query entry mode"
          />
        </Group>
      )}

      <div style={{ padding: hero ? 0 : 8 }}>
        {mode === 'search' ? (
          <>
            <Group gap={6} wrap="nowrap" align="flex-start">
              {/* The ghost is a real element, not the native placeholder. ::placeholder
                  does not reliably animate opacity in Chrome, so transitioning it swapped
                  the text with no tween: a flash. A span tweens like anything else. It is
                  aria-hidden and pointer-events:none, so the input keeps its own label and
                  stays clickable through it. */}
              <div style={{ flex: 1, position: 'relative' }}>
                {showGhost && (
                  <span
                    ref={ghostRef}
                    aria-hidden
                    className={`om-ghost-text om-mono${fading ? ' om-ghost-out' : ''}`}
                  >
                    {forms[slot]?.text ?? EXAMPLE}
                  </span>
                )}
              <TextInput
                ref={inputRef}
                style={{ width: '100%' }}
                className={hero ? 'om-hero-input' : undefined}
                classNames={{ input: 'om-mono' }}
                aria-label="Variant: HGVS, rsID, or free text"
                value={text}
                disabled={busy}
                onChange={(e) => setText(e.currentTarget.value)}
                onKeyDown={(e) => e.key === 'Enter' && submitSearch()}
              />
              </div>
              <Button
                onClick={submitSearch}
                loading={parsing}
                disabled={busy || !text.trim()}
                size={hero ? 'md' : 'xs'}
                style={hero ? { height: 44 } : undefined}
              >
                Resolve
              </Button>
              <Tooltip label={`Use the ABCC8 example: ${EXAMPLE}`} withArrow>
                <ActionIcon
                  variant="default"
                  size={hero ? 44 : 'lg'}
                  aria-label="Fill in the ABCC8 example variant"
                  onClick={() => setText(EXAMPLE)}
                >
                  <Text size="xs">ex</Text>
                </ActionIcon>
              </Tooltip>
            </Group>

            <Text size="xs" c="dimmed" mt={4}>
              HGVS with a transcript (NM_…:c.…), an rsID, or the variant described in words.
            </Text>

            {nlError && (
              <Alert color="red" mt={6} title="Could not parse" role="alert">
                <Text size="xs">{nlError}</Text>
              </Alert>
            )}

          </>

        ) : (
          <>
            <Group gap={6} align="flex-end" wrap="wrap">
              <TextInput
                style={{ flex: '2 1 260px' }}
                ref={inputRef}
                classNames={{ input: 'om-mono' }}
                label="variant"
                description="HGVS or rsID, verbatim"
                placeholder={EXAMPLE}
                value={q.variant}
                onChange={(e) => setQ({ ...q, variant: e.currentTarget.value })}
                onKeyDown={(e) => e.key === 'Enter' && submitManual()}
              />
              <Autocomplete
                style={{ flex: '1 1 130px' }}
                label="gene"
                description={geneErr ? 'suggestions offline' : 'hint only'}
                placeholder="ABCC8"
                data={geneOpts}
                value={q.gene}
                onChange={(v) => setQ({ ...q, gene: v })}
              />
              <NumberInput
                style={{ flex: '1 1 120px' }}
                label="window_bp"
                description="± around variant"
                min={1000}
                max={2_000_000}
                step={25_000}
                thousandSeparator=","
                value={q.window_bp}
                onChange={(v) => setQ({ ...q, window_bp: typeof v === 'number' ? v : 250_000 })}
              />
              <Select
                style={{ flex: '0 0 110px' }}
                label="build"
                data={['GRCh38', 'GRCh37']}
                value={q.build}
                allowDeselect={false}
                onChange={(v) => setQ({ ...q, build: v ?? 'GRCh38' })}
              />
              <Select
                style={{ flex: '0 0 110px' }}
                label="ancestry"
                placeholder="Global"
                clearable
                data={[...ANCESTRIES]}
                value={q.ancestry}
                onChange={(v) => setQ({ ...q, ancestry: v })}
              />
              <NumberInput
                style={{ flex: '0 0 110px' }}
                label="common_maf"
                min={0}
                max={0.5}
                step={0.01}
                decimalScale={3}
                value={q.common_maf}
                onChange={(v) => setQ({ ...q, common_maf: typeof v === 'number' ? v : 0.05 })}
              />
              <Checkbox
                mb={6}
                size="xs"
                label="cross_check"
                checked={q.cross_check}
                onChange={(e) => setQ({ ...q, cross_check: e.currentTarget.checked })}
              />
              <Button mb={2} onClick={submitManual} disabled={busy || !q.variant.trim()}>
                Resolve
              </Button>
            </Group>
            {q.build === 'GRCh37' && (
              <Text size="xs" c="dimmed" mt={6}>
                Panels are computed on GRCh38. A GRCh37 position is converted via ClinVar's assembly
                mapping and the conversion is labelled on the resolved variant.
              </Text>
            )}
          </>
        )}
      </div>
    </>
  )

  return hero ? <div>{body}</div> : <Paper mb="sm">{body}</Paper>
}

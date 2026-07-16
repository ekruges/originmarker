import { useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react'
import {
  ActionIcon, Alert, Autocomplete, Button, Checkbox, Group, NumberInput,
  Paper, SegmentedControl, Select, Text, TextInput, Tooltip,
} from '@mantine/core'
import { ANCESTRIES, api, ApiError, type Health, type NLResponse, type StructuredQuery } from './api'

const EXAMPLE = 'NM_000352.6(ABCC8):c.3989-9G>A'

/** The rotation: every accepted input form, over a spread of real loci.
 *
 *  Every identifier here was resolved against the live APIs before being listed. An
 *  example that does not resolve is worse than no example, and a genomic HGVS
 *  (NC_000011.10:g.17397055C>T) was listed here until it turned out ClinVar's search
 *  cannot find one, so the box was advertising a form that always failed.
 *
 *  `nl` marks the ones that need the model, and they are dropped when it is switched off.
 */
const PLACEHOLDERS: { text: string; nl?: boolean }[] = [
  { text: EXAMPLE },                                    // HGVS with gene, ABCC8, chr11
  { text: 'rs334' },                                    // sickle cell, HBB, chr11
  { text: 'rs113993960' },                              // CF F508del, CFTR, chr7
  { text: 'rs6025 in Europeans' },                      // factor V Leiden, F5, chr1
  { text: 'rs1800562' },                                // haemochromatosis, HFE, chr6
  { text: 'rs80338939 with a 500kb window' },           // GJB2 deafness, chr13
  { text: 'rs28941770' },                               // Tay-Sachs, HEXA, chr15
  { text: 'rs61750240' },                               // Rett, MECP2, chrX
  { text: 'NM_000518.5(HBB):c.20A>T, MAF at least 0.1' },
  { text: 'VCV000009088' },                             // ClinVar accession
  { text: 'the sickle cell mutation, in Africans', nl: true },
]

/** Fisher-Yates. The first example a visitor sees should not always be the same one. */
function shuffled<T>(xs: T[]): T[] {
  const a = [...xs]
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1))
    ;[a[i], a[j]] = [a[j], a[i]]
  }
  return a
}

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

  /** Everything typed here goes to /api/nl. There is one parser and it is the server's.
   *
   *  This used to shortcut anything "identifier-shaped" straight to /api/resolve, deciding
   *  with a client-side regex. It was a substring test, so "NM_...:c.20A>T, MAF at least
   *  0.1" counted as an identifier and the whole line, modifiers included, was sent as the
   *  variant. It also refused "rs334 in Europeans" when free text was off, though reading
   *  a modifier needs no model.
   *
   *  /api/nl costs nothing for text carrying an identifier: it reads it by regex, and only
   *  reaches the model when nothing else can. The rate limiter meters the model, not this.
   */
  const submitSearch = async () => {
    const t = text.trim()
    if (!t) return
    setNlError(null)
    setParsing(true)
    try {
      const r = await api.nl(t)
      // Provenance only when a model actually chose the variant, so a typed identifier
      // does not carry a caveat about a model that never ran.
      onResolve(r.query, r.used_llm ? r : undefined)
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

  /** Fill the box with one of the advertised examples, chosen at random.
   *
   *  Drawn from `forms`, not from PLACEHOLDERS: that is the same set the ghost rotates
   *  through, already filtered, so the button can never offer a free-text example on an
   *  instance where free text is switched off. Never picks what is already in the box,
   *  or a click would sometimes look broken.
   */
  const fillExample = () => {
    const pool = forms.filter((f) => f.text !== text)
    const pick = pool[Math.floor(Math.random() * pool.length)] ?? forms[0]
    if (!pick) return
    setText(pick.text)
    setNlError(null)
    inputRef.current?.focus()
  }

  // Ghost examples cycle through the accepted input forms, the way a search box does. Only
  // while the box is empty: a placeholder moving under text someone is composing is noise,
  // and the element is aria-hidden from the input's own label anyway.
  // Shuffled once per mount, not per render: a new order on every keystroke would reshuffle
  // the box under the reader. Recomputed when nl_enabled arrives, since that changes the set.
  const forms = useMemo(
    () => shuffled(PLACEHOLDERS.filter((p) => !p.nl || health?.nl_enabled)),
    [health?.nl_enabled],
  )
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
              <Tooltip label="Fill in a random example" withArrow>
                <ActionIcon
                  variant="default"
                  size={hero ? 44 : 'lg'}
                  disabled={busy}
                  aria-label="Fill in a random example query"
                  onClick={fillExample}
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

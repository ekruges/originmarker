import { useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react'
import {
  ActionIcon, Alert, Autocomplete, Button, Checkbox, CloseButton, Combobox, Group,
  NumberInput, Paper, SegmentedControl, Select, Text, TextInput, Tooltip, useCombobox,
} from '@mantine/core'
import {
  ANCESTRIES, api, ApiError, type Health, type NLResponse, type PrimerBuild,
  type StructuredQuery,
} from './api'
import { clearHistory, forgetQuery, readHistory, recordQuery, type Entry } from './history'
import { PrimerFields } from './PrimerOptions'

const EXAMPLE = 'NM_000352.6(ABCC8):c.3989-9G>A'

/** The scope a build gets when nobody picks one. Must stay panelbuilder's own default for
 *  StructuredQuery.primer_scope: it is the value this form shows before anyone touches it,
 *  and the form sends whatever it shows, so a drift here changes what gets built. */
const DEFAULT_SCOPE = 'starred'

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

/** A stored query can be longer than the box. It truncates rather than wrapping the row to
 *  two lines or pushing the count and the ex off the edge. */
const ROW_TEXT = {
  flex: 1, minWidth: 0, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap',
} as const

/** History rows worth offering for what is in the box: everything, until there is something
 *  to narrow on. Matched case-insensitively, which is a search over a list, not the verbatim
 *  identity `history.ts` stores entries under. */
const matching = (es: Entry[], text: string) => {
  const t = text.trim().toLowerCase()
  return t ? es.filter((e) => e.text.toLowerCase().includes(t)) : es
}

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

  // Read once, at mount: localStorage is synchronous, and re-reading it per render would
  // put a disk hit on every keystroke for a list this component already holds.
  const [entries, setEntries] = useState<Entry[]>(readHistory)
  const combobox = useCombobox({ onDropdownClose: () => combobox.resetSelectedOption() })
  const shown = matching(entries, text)

  // The selected index is a ref into the LIVE option list, so any change to `shown` leaves
  // it aimed at whichever row now sits at that index: Enter would then submit a query the
  // user never highlighted. Reset it wherever the rows change, not per writer of text or
  // entries. Left open over an emptied list, the input would also claim aria-expanded over
  // a listbox that is gone.
  useEffect(() => {
    combobox.resetSelectedOption()
    if (!shown.length) combobox.closeDropdown()
    // `shown` is a pure function of these two. useCombobox returns a fresh object every
    // render, so listing it here would reset the selection mid-keystroke and kill the arrows.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [text, entries])

  // Manual mode: every StructuredQuery knob, no parsing in the way.
  const [q, setQ] = useState<Required<Omit<StructuredQuery,
    'gene' | 'ancestry' | 'primer_scope' | 'primer_settings'>> & {
    gene: string
    ancestry: string | null
  }>({
    variant: '', gene: '', window_bp: 250_000, build: 'GRCh38',
    ancestry: null, common_maf: 0.05, cross_check: true,
  })

  // Untouched until the fields are edited. The numbers are never held here: they come off
  // health, from the engine that will use them, so this side has no defaults to drift.
  const [primer, setPrimer] = useState<PrimerBuild | null>(null)

  // What the form shows, and therefore exactly what a build is asked for: sent whenever the
  // form is drawn, touched or not, so the screen cannot say one scope while the server picks
  // another. Null where the server states no defaults, and then it decides everything.
  const primerShown: PrimerBuild | null = health?.primer_defaults
    ? primer ?? { scope: DEFAULT_SCOPE, settings: health.primer_defaults }
    : null

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
  const submitSearch = async (pick?: string) => {
    const t = (pick ?? text).trim()
    if (!t) return
    combobox.closeDropdown()
    setNlError(null)
    setParsing(true)
    try {
      const r = await api.nl(t)
      // Recorded once the parse stands, so the half-typed identifier that failed on the way
      // to this one does not push a real query out of a capped list.
      setEntries(recordQuery(t, r.query.variant))
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
      primer_scope: primerShown?.scope,
      primer_settings: primerShown?.settings,
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
              {/* Not in a portal: the dropdown's own buttons are then reachable by Tab from
                  the input, which is the only keyboard path to Clear all. */}
              <Combobox
                store={combobox}
                withinPortal={false}
                position="bottom-start"
                width="target"
                onOptionSubmit={(v) => { setText(v); submitSearch(v) }}
              >
                <Combobox.Target withExpandedAttribute withKeyboardNavigation={false}>
              <TextInput
                ref={inputRef}
                style={{ width: '100%' }}
                className={hero ? 'om-hero-input' : undefined}
                classNames={{ input: 'om-mono' }}
                aria-label="Variant: HGVS, rsID, or free text"
                value={text}
                disabled={busy}
                onChange={(e) => {
                  // Reset here, not in an effect keyed on `text`: the store's index is a ref
                  // into the option list, an effect runs after paint, and an Enter arriving
                  // in between reads an index aimed at whichever row now sits there. That
                  // submits a variant the user never highlighted.
                  combobox.resetSelectedOption()
                  setText(e.currentTarget.value)
                }}
                // Opened by click, never by focus: the box focuses itself on mount and the
                // ex button focuses it again, and neither is someone asking for a list.
                onClick={() => shown.length && combobox.openDropdown()}
                // Every key handled here, with Combobox.Target's own navigation switched
                // off. Sharing the input between two handlers means relying on which one
                // cloneElement keeps, and the answer was neither: the arrows silently did
                // nothing while both sides looked correct in isolation.
                onKeyDown={(e) => {
                  const open = combobox.dropdownOpened && shown.length > 0
                  const active = open ? combobox.getSelectedOptionIndex() : -1
                  if (e.key === 'ArrowDown') {
                    e.preventDefault()
                    if (!open && shown.length) combobox.openDropdown()
                    combobox.selectNextOption()
                    return
                  }
                  if (e.key === 'ArrowUp') {
                    e.preventDefault()
                    if (open) combobox.selectPreviousOption()
                    return
                  }
                  if (e.key === 'Escape') {
                    if (open) { e.preventDefault(); combobox.closeDropdown() }
                    return
                  }
                  if (e.key === 'Enter') {
                    e.preventDefault()
                    // A highlighted row wins; a plain Enter submits what was typed. Reading
                    // the row out of `shown` rather than the index alone: the index is a ref
                    // into a list that re-filters as the user types.
                    // Trustworthy because onChange resets the index in the same tick it
                    // changes the list: nothing can land in between aiming it at a row that
                    // moved. An effect could not promise that, and did not.
                    const row = active >= 0 ? shown[active] : undefined
                    if (row) { setText(row.text); submitSearch(row.text) } else submitSearch()
                    return
                  }
                  // The row's own ex is a pointer affordance. This is its keyboard twin.
                  if (e.key === 'Delete' && shown[active]) {
                    e.preventDefault()
                    setEntries(forgetQuery(shown[active].text))
                  }
                }}
              />
                </Combobox.Target>
                {shown.length > 0 && (
                  <Combobox.Dropdown>
                    <Combobox.Options>
                      {shown.map((e) => (
                        <Combobox.Option
                          key={e.text}
                          value={e.text}
                          // An option's children are presentational, so the count and the ex
                          // are not read from the row: the row has to say it itself.
                          aria-label={e.count === undefined ? e.text : `${e.text}, ${e.count} candidates`}
                        >
                          <Group gap={8} wrap="nowrap">
                            <Text size="xs" className="om-mono" style={ROW_TEXT}>{e.text}</Text>
                            {e.count !== undefined && (
                              <Text
                                size="xs"
                                className="om-mono"
                                c="dimmed"
                                title={`${e.count} candidates in the panel built from this query`}
                              >
                                {e.count}
                              </Text>
                            )}
                            <CloseButton
                              size="xs"
                              c="black"
                              tabIndex={-1}
                              aria-hidden
                              title={`Forget ${e.text}`}
                              // Or the row underneath reads the click as "run this query".
                              onClick={(ev) => {
                                ev.stopPropagation()
                                setEntries(forgetQuery(e.text))
                              }}
                            />
                          </Group>
                        </Combobox.Option>
                      ))}
                    </Combobox.Options>
                    <Combobox.Footer>
                      <Group justify="space-between" wrap="nowrap" gap={8}>
                        <Text size="xs" c="dimmed">Kept in this browser only</Text>
                        <Button
                          variant="subtle"
                          size="compact-xs"
                          onClick={() => setEntries(clearHistory())}
                        >
                          Clear all
                        </Button>
                      </Group>
                    </Combobox.Footer>
                  </Combobox.Dropdown>
                )}
              </Combobox>
              </div>
              <Button
                onClick={() => submitSearch()}
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
            {/* Only against a server that states its defaults. Without them there is no
                primer feature to configure, and a form here would collect settings that
                nothing downstream reads. */}
            {primerShown && (
              <div
                style={{
                  marginTop: 10,
                  paddingTop: 8,
                  borderTop: '1px solid var(--om-border)',
                }}
              >
                <Text size="xs" fw={600} mb={5}>Primers</Text>
                <PrimerFields value={primerShown} onChange={setPrimer} />
              </div>
            )}
          </>
        )}
      </div>
    </>
  )

  return hero ? <div>{body}</div> : <Paper mb="sm">{body}</Paper>
}

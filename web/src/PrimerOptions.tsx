import { useEffect, useState } from 'react'
import { Alert, Anchor, Button, Group, NumberInput, Popover, Select, Text } from '@mantine/core'
import {
  PRIMER_DOCS, PRIMER_SCOPES, PRIMER_SETTING_KEYS, PRIMER_SIZE_CAP,
  type PrimerBuild, type PrimerScope,
} from './api'

/** What each scope designs for, as panelbuilder.PRIMER_SCOPES defines the set. */
const SCOPE_LABEL: Record<PrimerScope, string> = {
  starred: 'Markers meeting the flanking criteria',
  recommended: 'Every shortlisted marker',
  none: 'No primers',
}

export type GroupKey = 'tm' | 'size' | 'gc' | 'product' | 'salt' | 'mask'

interface Spec {
  group: GroupKey
  label: string
  min: number
  max: number
  step: number
  decimals: number
}

/**
 * The knobs this form draws, with their bounds and their words.
 *
 * A deliberate subset of PRIMER_SETTING_KEYS. The rest are structural rather than tuning:
 * max_ns_accepted is 0 because a primer over an unknown base is a fabricated primer (R1),
 * and the Tm formula and salt correction choose the model the stated Tm means anything
 * under. They ride a rebuild untouched and stay settable through the API.
 */
export const PRIMER_FIELDS = {
  min_tm: { group: 'tm', label: 'min Tm', min: 40, max: 85, step: 0.5, decimals: 1 },
  opt_tm: { group: 'tm', label: 'opt Tm', min: 40, max: 85, step: 0.5, decimals: 1 },
  max_tm: { group: 'tm', label: 'max Tm', min: 40, max: 85, step: 0.5, decimals: 1 },
  max_pair_diff_tm: { group: 'tm', label: 'max pair ΔTm', min: 0, max: 20, step: 0.5, decimals: 1 },
  min_size: { group: 'size', label: 'min length', min: 15, max: PRIMER_SIZE_CAP, step: 1, decimals: 0 },
  opt_size: { group: 'size', label: 'opt length', min: 15, max: PRIMER_SIZE_CAP, step: 1, decimals: 0 },
  max_size: { group: 'size', label: 'max length', min: 15, max: PRIMER_SIZE_CAP, step: 1, decimals: 0 },
  min_gc: { group: 'gc', label: 'min GC %', min: 0, max: 100, step: 1, decimals: 1 },
  max_gc: { group: 'gc', label: 'max GC %', min: 0, max: 100, step: 1, decimals: 1 },
  gc_clamp: { group: 'gc', label: 'GC clamp', min: 0, max: 5, step: 1, decimals: 0 },
  max_poly_x: { group: 'gc', label: 'max poly-X', min: 1, max: 10, step: 1, decimals: 0 },
  min_product: { group: 'product', label: 'min product', min: 50, max: 3000, step: 25, decimals: 0 },
  max_product: { group: 'product', label: 'max product', min: 50, max: 3000, step: 25, decimals: 0 },
  salt_monovalent: { group: 'salt', label: 'monovalent mM', min: 0, max: 500, step: 5, decimals: 1 },
  salt_divalent: { group: 'salt', label: 'divalent mM', min: 0, max: 50, step: 0.1, decimals: 2 },
  dntp_conc: { group: 'salt', label: 'dNTP mM', min: 0, max: 20, step: 0.1, decimals: 2 },
  dna_conc: { group: 'salt', label: 'DNA nM', min: 0, max: 2000, step: 10, decimals: 1 },
  mask_maf: { group: 'mask', label: 'mask MAF floor', min: 0, max: 0.499, step: 0.005, decimals: 4 },
  target_pad: { group: 'mask', label: 'target pad bp', min: 0, max: 300, step: 10, decimals: 0 },
} satisfies Partial<Record<(typeof PRIMER_SETTING_KEYS)[number], Spec>>

export type FieldKey = keyof typeof PRIMER_FIELDS
export const FIELD_KEYS = Object.keys(PRIMER_FIELDS) as FieldKey[]

/**
 * Section headings, and nothing else.
 *
 * The constraint each group carries is real and is written down, but not here: anyone who
 * opens this form knows what a Tm is, and six paragraphs above six rows of number boxes is a
 * form nobody can see. They live in the docs, under one link at the top of this form. The
 * docs table is built from PRIMER_FIELDS itself against a Record<FieldKey, string>, so a
 * knob added here and not documented there fails the typecheck rather than shipping blank.
 */
export const PRIMER_GROUPS: { key: GroupKey; title: string }[] = [
  { key: 'tm', title: 'Melting temperature (°C)' },
  { key: 'size', title: 'Length (bases)' },
  { key: 'gc', title: 'Composition' },
  { key: 'product', title: 'Product (bp)' },
  { key: 'salt', title: 'Reaction conditions' },
  { key: 'mask', title: 'Mask' },
]

export const samePrimerBuild = (a: PrimerBuild, b: PrimerBuild) =>
  a.scope === b.scope && PRIMER_SETTING_KEYS.every((k) => a.settings[k] === b.settings[k])

interface FieldsProps {
  value: PrimerBuild
  onChange: (b: PrimerBuild) => void
  /** Read-only, not disabled: these numbers are also the panel's provenance, and a greyed
   *  form states them less legibly than a live one. */
  readOnly?: boolean
}

/** Every primer setting this form draws, in one grid. Mounted by the results chip and by
 *  manual input, so the two paths cannot offer different knobs. */
export function PrimerFields({ value, onChange, readOnly }: FieldsProps) {
  const set = (k: FieldKey, v: string | number) =>
    typeof v === 'number' && onChange({ ...value, settings: { ...value.settings, [k]: v } })
  return (
    <div>
      <Select
        size="xs"
        w={280}
        label="Design primers for"
        readOnly={readOnly}
        allowDeselect={false}
        value={value.scope}
        onChange={(v) => v && onChange({ ...value, scope: v as PrimerScope })}
        data={PRIMER_SCOPES.map((s) => ({ value: s, label: SCOPE_LABEL[s] }))}
      />
      {value.scope !== 'none' && (
        <Text size="xs" c="dimmed" mt={8}>
          Every field, its bounds and what it constrains:{' '}
          <Anchor href={PRIMER_DOCS} size="xs">primer reference</Anchor>
        </Text>
      )}
      {value.scope !== 'none' && PRIMER_GROUPS.map((s) => (
        <div key={s.key} style={{ marginTop: 12 }}>
          <Text size="xs" fw={600} mb={5}>{s.title}</Text>
          <Group gap={6} wrap="wrap">
            {FIELD_KEYS.filter((k) => PRIMER_FIELDS[k].group === s.key).map((k) => (
              <NumberInput
                key={k}
                size="xs"
                w={104}
                classNames={{ input: 'om-mono' }}
                label={PRIMER_FIELDS[k].label}
                readOnly={readOnly}
                min={PRIMER_FIELDS[k].min}
                max={PRIMER_FIELDS[k].max}
                step={PRIMER_FIELDS[k].step}
                // Clamped on blur as well as on the arrows: a typed number outside the
                // bounds is otherwise submitted, and the cap is the whole point of one.
                clampBehavior="blur"
                decimalScale={PRIMER_FIELDS[k].decimals}
                value={value.settings[k]}
                onChange={(v) => set(k, v)}
              />
            ))}
          </Group>
        </div>
      ))}
    </div>
  )
}

/**
 * The advanced affordance: what this panel was built under, and a rebuild.
 *
 * Without `onRebuild` the page it sits on has no way to build, so the same settings render
 * read-only and point at the path that can. They are provenance either way.
 */
export function PrimerChip({
  build, onRebuild, onVerify, verify,
}: {
  build: PrimerBuild
  onRebuild?: (b: PrimerBuild) => void
  /** Undefined where the server has no UCSC key: the box then says the check is off rather
   *  than offering a button that cannot run. */
  onVerify?: () => void
  verify?: { running: boolean; stage: string; error?: string }
}) {
  const [open, setOpen] = useState(false)
  const [draft, setDraft] = useState(build)
  // A rebuild lands a new panel with new params, and the draft must follow it rather than
  // keep editing settings the page no longer shows.
  useEffect(() => setDraft(build), [build])
  const unchanged = samePrimerBuild(draft, build)

  return (
    <Popover
      opened={open}
      onChange={setOpen}
      width={560}
      position="bottom-start"
      withArrow
      shadow="md"
      trapFocus
    >
      <Popover.Target>
        <Button
          variant="default"
          size="compact-xs"
          radius="xl"
          fw={400}
          aria-expanded={open}
          aria-label="Primer settings"
          onClick={() => setOpen((o) => !o)}
        >
          Primers · {build.scope} · Tm {build.settings.opt_tm} °C {open ? '▴' : '▾'}
        </Button>
      </Popover.Target>
      <Popover.Dropdown mah={520} style={{ overflowY: 'auto' }}>
        <PrimerFields value={draft} onChange={setDraft} readOnly={!onRebuild} />
        {/* Verification is its own act, after the panel and never inside it: UCSC allows one
            request every 15 seconds, so this is minutes of waiting for an answer the panel
            does not need in order to be usable. The button says the cost before it is paid. */}
        <div style={{ borderTop: '1px solid var(--om-border)', marginTop: 12, paddingTop: 10 }}>
          <Text size="xs" fw={600} mb={4}>Check against the genome</Text>
          {onVerify ? (
            <>
              <Text size="xs" c="dimmed" mb={6}>
                Checks each pair against the whole of hg38. About 15 s per pair.{' '}
                <Anchor href={PRIMER_DOCS} size="xs">What it does</Anchor>
              </Text>
              {verify?.error && (
                <Alert color="red" p={6} mb={6} role="alert">
                  <Text size="xs">{verify.error}</Text>
                </Alert>
              )}
              <Group gap={8} align="center">
                <Button size="xs" variant="default" loading={verify?.running} onClick={onVerify}>
                  {verify?.running ? 'Checking' : 'Check pairs'}
                </Button>
                {verify?.running && (
                  <Text size="xs" c="dimmed" className="om-mono">{verify.stage}</Text>
                )}
              </Group>
            </>
          ) : (
            <Text size="xs" c="dimmed">
              Not available on this instance: it needs a UCSC API key. Pairs stay marked as
              not checked, which is what they are.{' '}
              <Anchor href={PRIMER_DOCS} size="xs">Setting one up</Anchor>
            </Text>
          )}
        </div>

        {onRebuild ? (
          <Group justify="flex-end" gap={6} mt={12}>
            <Button size="xs" variant="default" disabled={unchanged} onClick={() => setDraft(build)}>
              Reset
            </Button>
            <Button
              size="xs"
              disabled={unchanged}
              onClick={() => { setOpen(false); onRebuild(draft) }}
            >
              Rebuild
            </Button>
          </Group>
        ) : (
          <Text size="xs" c="dimmed" mt={12}>
            The settings this panel was built under. Change them under Manual input and build
            again.
          </Text>
        )}
      </Popover.Dropdown>
    </Popover>
  )
}

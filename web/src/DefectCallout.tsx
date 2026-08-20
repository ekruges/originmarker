/**
 * Every chromosomal change in one sample, as a grid of small chips that open.
 *
 * WHY IT IS NOT ONE BLOCK PER CHANGE ANY MORE. A real run produced 64 regions in a single sample,
 * and the previous layout gave each of them a headline, a row of evidence chips and a paragraph.
 * Fifty-eight of those paragraphs were the same sentence, because the parental origin on a
 * uniparental genome is inherited from one genome-level call and is therefore identical for every
 * region. The result was several screens of text carrying about three facts, with the confidence
 * number buried in the middle of it.
 *
 * So: the headline counts first, then one chip per region carrying the four things a reader scans
 * for, which are WHERE, WHAT, WHOSE and HOW SURE. Everything else is behind the chip and opens on
 * a click. Reasoning shared by more than one region is stated ONCE underneath rather than repeated
 * per region, which is where most of the words went.
 *
 * The confidence number is on the face of the chip, in the band's own colour, because it is the
 * thing that decides whether a row may be used and it was previously the hardest thing to find.
 */
import { useState } from 'react'
import { Text } from '@mantine/core'
import type { Defect } from './defects.ts'
import { headline, bandColour, BAND_WORD, originBlockedByClass } from './defects.ts'

const mb = (n: number): string => (n / 1e6).toFixed(2)

/** Short label for the parent, so the chip face stays one line. */
const parentMark = (d: Defect): string => (d.origin === 'unclear' || !d.origin ? '?'
  : d.origin === 'maternal' ? 'M' : d.origin === 'paternal' ? 'P' : '?')

const KIND_SHORT: Record<string, string> = {
  'cnn-loh': 'CNN-LOH',
  'copy-loss': 'LOSS',
  'copy-gain': 'GAIN',
  loss: 'LOSS',
  gain: 'GAIN',
  'segmental-duplication': 'DUP',
  'segmental-upd': 'UPD',
  isodisomy: 'ISO',
  triploidy: '3N',
  complex: 'CPLX',
}
const kindShort = (k: string): string => KIND_SHORT[k] ?? k.replace(/[-_]/g, ' ').toUpperCase()

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

/** Sentences, kept with their terminator so rejoining them reads normally. */
const sentences = (s: string): string[] =>
  (s.match(/[^.]+\.?/g) ?? []).map((x) => x.trim()).filter(Boolean)

/**
 * Split every reason into the part unique to one region and the part many regions share.
 *
 * The shared part is real content and is not dropped, it is stated once. On a uniparental genome
 * that is the whole explanation of the parental call, which is identical for every region by
 * construction, and repeating it per region is what made the panel unreadable.
 */
function splitReasons(defects: readonly Defect[]): { unique: string[][]; shared: string[] } {
  const count = new Map<string, number>()
  const per = defects.map((d) => sentences(d.why ?? ''))
  for (const ss of per) for (const s of new Set(ss)) count.set(s, (count.get(s) ?? 0) + 1)
  const isShared = (s: string) => (count.get(s) ?? 0) > 1 && defects.length > 1
  const shared: string[] = []
  for (const ss of per) for (const s of ss) if (isShared(s) && !shared.includes(s)) shared.push(s)
  return { unique: per.map((ss) => ss.filter((s) => !isShared(s))), shared }
}

function OneDefect({ d, why }: { d: Defect; why: string[] }) {
  const [open, setOpen] = useState(false)
  const colour = bandColour(d)
  const named = d.origin && d.origin !== 'unclear'
  /**
   * DIGITS ONLY WHERE THEY MEAN SOMETHING.
   *
   * Bands A and B carry their number. C and D do not, and that is the lesson of the PGT-A
   * mosaicism episode rather than a style choice: a low-specificity intermediate category gets
   * quoted downstream regardless of the label attached to it, and the field's corrective was to
   * withdraw the category rather than improve the caveat. Band D sits around 0.63 measured on its
   * own range, which is thirteen points over guessing, and printing it to three decimals asserts a
   * precision the quantity does not have. The measurement is still in the expanded detail and in
   * the export; it is the row face that stops carrying it.
   */
  const showsDigits = d.band === 'A' || d.band === 'B'
  const conf = showsDigits && d.confidence !== undefined && Number.isFinite(d.confidence)
    ? d.confidence.toFixed(2) : null
  // What stands in place of the digits, per band. F names NO parent: an injection series measured
  // it recovering the right one 0.51 to 0.56 of the time, which is chance, because every call reaching
  // it has an unresolved class and a gain inverts the sign that loss and copy-neutral share. The
  // row is still graded, which is what keeps a refusal the same kind of output as an answer.
  const weakWord = d.band === 'C' ? 'weak'
    : d.band === 'D' ? 'not evaluable'
      : d.band === 'F' ? 'no parent' : null
  return (
    <div style={{ border: `1px solid ${open ? colour : 'var(--om-defect)'}`, background: 'var(--om-defect-chip)' }}>
      <button
        type="button"
        onClick={() => setOpen((v) => !v)}
        aria-expanded={open}
        title={headline(d)}
        style={{
          display: 'flex', alignItems: 'baseline', gap: 7, width: '100%',
          border: 0, background: 'none', cursor: 'pointer', padding: '3px 7px',
          font: 'inherit', textAlign: 'left', lineHeight: 1.5,
        }}
      >
        <span style={{
          fontFamily: 'var(--om-mono)', fontSize: 12, fontWeight: 600, whiteSpace: 'nowrap',
        }}
        >
          chr{d.chrom}
        </span>
        <span style={{ fontFamily: 'var(--om-mono)', fontSize: 12, fontVariantNumeric: 'tabular-nums' }}>
          {mb(d.endBp - d.startBp)}Mb
        </span>
        <span style={{ fontSize: 11, fontWeight: 700, letterSpacing: '0.03em', color: 'var(--om-defect)' }}>
          {kindShort(d.kind)}
        </span>
        <span style={{ flex: 1 }} />
        {/* WHOSE and HOW SURE, together, in the band's colour. The two things a reader scans for. */}
        <span style={{ fontSize: 12, fontWeight: 700, color: named ? colour : 'var(--om-text-dim)' }}>
          {parentMark(d)}
        </span>
        {d.band === 'inherited' ? (
          // A verdict, not a number. The margin says how strong the genome-level call was; there
          // is no per-event accuracy to print and printing one invites it being quoted as one.
          <span style={{ fontSize: 11, fontWeight: 700, color: colour }}>
            inherited{d.inheritedMargin !== undefined && Number.isFinite(d.inheritedMargin)
              ? ` ${d.inheritedMargin.toFixed(1)}x` : ''}
          </span>
        ) : conf ? (
          <span style={{
            fontFamily: 'var(--om-mono)', fontSize: 12, fontWeight: 700, color: colour,
            fontVariantNumeric: 'tabular-nums',
          }}
          >
            {conf}
          </span>
        ) : (
          // A CLASS THAT CANNOT HAVE AN ORIGIN IS NOT A MEASUREMENT THAT FAILED. A triploidy or a
          // genome-level complex call has no parent to name at any quality, which is a different
          // statement from "we could not measure it here", and the two must not read alike.
          <span style={{ fontSize: 11, color: 'var(--om-text-dim)' }}>
            {weakWord ?? (d.originBlocked || originBlockedByClass(d.kind) ? 'n/a' : 'no call')}
          </span>
        )}
        {d.band && d.band !== 'inherited' ? (
          <span style={{
            fontSize: 10, fontWeight: 700, color: colour, border: `1px solid ${colour}`,
            padding: '0 3px', lineHeight: 1.4,
          }}
          >
            {d.band}
          </span>
        ) : null}
      </button>

      {open ? (
        <div style={{ padding: '2px 7px 7px', borderTop: '1px solid var(--om-defect)' }}>
          <Text style={{ fontSize: 13, fontWeight: 700, lineHeight: 1.3, color: colour, margin: '6px 0 6px' }}>
            {headline(d)}
          </Text>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 5 }}>
            <Chip label="at" value={d.locus} mono />
            {d.interval ? <Chip label="edge" value={d.interval} mono /> : null}
            {d.informative !== undefined
              ? <Chip label="markers" value={String(d.informative)} mono /> : null}
            {d.posterior !== undefined && Number.isFinite(d.posterior)
              ? <Chip label="posterior" value={d.posterior.toFixed(3)} mono /> : null}
            {d.exclusive !== undefined
              ? <Chip label="exclusive" value={String(d.exclusive)} mono /> : null}
            {d.z !== undefined ? <Chip label="z" value={d.z.toFixed(2)} mono /> : null}
            {d.impliedF !== undefined
              ? <Chip label="fraction" value={d.impliedF.toFixed(2)} mono /> : null}
            {d.limitedBy && d.limitedBy !== 'none'
              ? <Chip label="limited by" value={d.limitedBy.replace(/-/g, ' ')} /> : null}
            {d.basis ? <Chip label="basis" value={d.basis.replace('-', ' ')} /> : null}
            {d.channel ? <Chip label="channel" value={d.channel} /> : null}
            {d.stage ? <Chip label="material" value={d.stage} /> : null}
            {d.dropout !== undefined
              ? <Chip label="dropout" value={d.dropout.toFixed(3)} mono /> : null}
          </div>
          {why.length ? (
            <Text size="sm" mt={6} style={{ maxWidth: 780, lineHeight: 1.5 }}>{why.join(' ')}</Text>
          ) : null}
        </div>
      ) : null}
    </div>
  )
}

export function DefectCallout({ defects }: { defects: Defect[] }) {
  if (!defects.length) return null
  const { unique, shared } = splitReasons(defects)

  // The counts a reader wants before any coordinate: how many, whose, and how many are strong
  // enough to use. Bands A and B are the reportable ones, so they are counted separately from the
  // rest rather than left for the reader to total.
  const byParent = new Map<string, number>()
  for (const d of defects) {
    const k = d.origin && d.origin !== 'unclear' ? d.origin : 'undetermined'
    byParent.set(k, (byParent.get(k) ?? 0) + 1)
  }
  const reportable = defects.filter((d) => d.band === 'A' || d.band === 'B').length
  // Rows whose parent was INHERITED from one genome-level call rather than measured here. A reader
  // counts rows, so the count has to be stated before the rows are shown.
  const inherited = defects.filter((d) => d.band === 'inherited')
  const named = defects.length - (byParent.get('undetermined') ?? 0)
  const parts = [...byParent.entries()]
    .sort((a, b) => b[1] - a[1])
    .map(([k, n]) => `${n} ${k}`)

  return (
    <div style={{
      border: '1px solid var(--om-defect)', borderLeft: '4px solid var(--om-defect)',
      background: 'var(--om-defect-bg)', padding: '11px 14px', margin: '10px 0 4px',
    }}
    >
      <Text style={{ fontSize: 22, fontWeight: 700, color: 'var(--om-defect)', lineHeight: 1.15 }}>
        {defects.length} chromosomal change{defects.length === 1 ? '' : 's'}
      </Text>
      <Text style={{ fontSize: 15, fontWeight: 600, lineHeight: 1.35, marginTop: 2 }}>
        {parts.join(' · ')}
      </Text>
      <Text size="sm" c="dimmed" mt={2}>
        {named
          ? `${reportable} of ${defects.length} measured at band A or B, the bands meant for `
            + `reporting${inherited.length ? `; ${inherited.length} inherited from one `
              + 'genome-level call' : ''}`
          : 'none carries a parental origin'}
        . Click a change for its evidence.
      </Text>

      {inherited.length > 1 ? (
        <div style={{
          marginTop: 9, border: '1px solid var(--om-defect)', borderLeft: '3px solid var(--om-defect)',
          background: 'var(--om-defect-chip)', padding: '7px 10px',
        }}
        >
          <Text style={{ fontSize: 13, fontWeight: 700, color: 'var(--om-defect)', lineHeight: 1.3 }}>
            ONE parental determination, not {inherited.length}
          </Text>
          <Text size="sm" mt={2} style={{ maxWidth: 780, lineHeight: 1.5 }}>
            {inherited.length} of these {defects.length} changes take their parent from a single
            genome-level call, not from anything measured on the change itself. They are one
            determination inherited {inherited.length} times, and counting them as {inherited.length}
            independent confirmations would overstate the evidence by that factor. If the
            genome-level call is wrong, all {inherited.length} are wrong together.
          </Text>
        </div>
      ) : null}

      <div style={{
        marginTop: 9, display: 'grid', gap: 4,
        gridTemplateColumns: 'repeat(auto-fill, minmax(265px, 1fr))', alignItems: 'start',
      }}
      >
        {defects.map((d, i) => (
          <OneDefect key={`${d.chrom}:${d.startBp}:${d.kind}`} d={d} why={unique[i] ?? []} />
        ))}
      </div>

      {shared.length ? (
        <div style={{ marginTop: 10, borderTop: '1px solid var(--om-defect)', paddingTop: 7 }}>
          <Text style={{ fontSize: 12, fontWeight: 700, color: 'var(--om-defect)' }}>
            Applies to more than one change above
          </Text>
          <Text size="sm" mt={3} style={{ maxWidth: 780, lineHeight: 1.5 }}>
            {shared.join(' ')}
          </Text>
        </div>
      ) : null}

      <Text size="xs" c="dimmed" mt={9} style={{ maxWidth: 780, lineHeight: 1.5 }}>
        Bands: A {BAND_WORD.A}, B {BAND_WORD.B}, C {BAND_WORD.C}, D {BAND_WORD.D}. Position comes
        from the intensity channel and needs no parental genotype, so it is reported for every
        region. Naming the parent needs a genotyped parent or a genome already known to carry one
        parent; where the run has neither the origin reads not determined and the coordinates stand
        on their own rather than being withheld.
      </Text>
    </div>
  )
}

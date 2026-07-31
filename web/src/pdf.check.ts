// Self-check for the PDF writer. Run: node src/pdf.check.ts
//
// A malformed PDF fails in a reader, not here, so what is checked is the structure a reader
// walks: the xref offsets have to land on their own objects, every stream length has to be the
// byte count, and no line of text may exceed the width it was wrapped to.
import assert from 'node:assert/strict'
import { LETTER, Pdf, stringWidth, wrap } from './pdf.ts'
import { syngamyLogText } from './logfile.ts'

// --- 1. metrics are Adobe's, not invented ------------------------------------------------------
assert.equal(stringWidth('W', 'Helvetica', 1000), 944)
assert.equal(stringWidth('i', 'Helvetica', 1000), 222)
assert.equal(stringWidth('W', 'Helvetica-Bold', 1000), 944)
assert.equal(stringWidth('i', 'Helvetica-Bold', 1000), 278)
assert.equal(stringWidth('MMMM', 'Courier', 10), 24, 'Courier is fixed pitch')
assert.ok(stringWidth('Syngamy', 'Helvetica-Bold', 14) > stringWidth('Syngamy', 'Helvetica', 14))

// --- 2. wrapping never overflows, including on a word wider than the column --------------------
for (const [s, w] of [
  ['Parent of origin was determined from SNP array genotypes using OriginMarker.', 120],
  ['a'.repeat(400), 60],
  ['e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855', 90],
  ['', 100],
] as const) {
  for (const line of wrap(s, 'Helvetica', 8, w)) {
    assert.ok(stringWidth(line, 'Helvetica', 8) <= w + 1e-9, `overflow: ${line}`)
  }
}
assert.deepEqual(wrap('one\ntwo', 'Helvetica', 8, 500), ['one', 'two'], 'newlines are kept')

// --- 3. the file a reader walks ----------------------------------------------------------------
const pdf = new Pdf(LETTER)
pdf.setTitle('check')
pdf.setFont('Helvetica-Bold', 14)
pdf.drawString(40, 700, 'Syngamy (test) - a (nested) paren \\ and a slash')
pdf.setFillColor('#337ab7')
pdf.rect(40, 600, 200, 20, true)
pdf.circle(300, 610, 12)
pdf.polys([[[40, 500], [60, 500], [50, 520]]])
pdf.drawRightString(572, 480, 'right edge')
pdf.showPage()
pdf.setFont('Courier', 6)
pdf.drawString(40, 700, 'page two, with a wide character: — and ≥')
const bytes = new Uint8Array(await pdf.save().arrayBuffer())
const text = Buffer.from(bytes).toString('latin1')

assert.ok(text.startsWith('%PDF-1.4\n'))
assert.ok(text.endsWith('%%EOF\n'))
assert.equal(pdf.pageCount, 2)
assert.match(text, /\/Count 2\b/)

// Every xref offset must land on the object it claims, or a reader gives up on the file.
const startxref = Number(/startxref\n(\d+)/.exec(text)![1])
assert.equal(text.slice(startxref, startxref + 4), 'xref')
const entries = [...text.slice(startxref).matchAll(/^(\d{10}) 00000 n $/gm)].map((m) => m[1])
assert.ok(entries.length >= 8, `too few xref entries: ${entries.length}`)
entries.forEach((off, i) => {
  assert.ok(text.startsWith(`${i + 1} 0 obj`, Number(off)),
    `xref entry ${i + 1} points at ${text.slice(Number(off), Number(off) + 12)!}`)
})

// A /Length that disagrees with the stream is the other way a reader gives up.
for (const m of text.matchAll(/<< \/Length (\d+) >>\nstream\n/g)) {
  const start = m.index! + m[0].length
  assert.equal(text.slice(start, start + Number(m[1])).length, Number(m[1]))
  assert.equal(text.slice(start + Number(m[1]), start + Number(m[1]) + 10), '\nendstream')
}

// Parens and backslashes are escaped, so the string operators stay balanced.
assert.ok(text.includes('\\(nested\\)') && text.includes('\\\\ and'))
// Nothing outside Latin-1 survives into the byte stream.
assert.ok(bytes.every((b) => b <= 0xff))
assert.ok(text.includes('- and >='), 'an em dash and a >= fold onto ASCII')

// --- 4. the log export carries its footer ------------------------------------------------------
const log = syngamyLogText(
  [{ tag: 'READ', text: 'a.csv 1.4 MB' }, { tag: 'DONE', text: 'biparental' }],
  { tool: 'OriginMarker 1.0.0', started: '2026-01-01T00:00:00Z' },
  { now: '2026-01-01T00:01:00Z', url: 'http://x/#/syngamy', agent: 'node' },
)
assert.ok(log.startsWith('[READ] a.csv 1.4 MB\n[DONE] biparental'))
assert.match(log, /^# syngamy build\/debug \| tool=OriginMarker 1\.0\.0 \|.*\| lines=2 \|/m)

console.log('pdf.check.ts OK')

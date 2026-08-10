// Self-check for the inferred-array export. Run: node src/inferredArray.check.ts
//
// The point of the export is that another tool can read it, so this reads it back with the real
// ingest rather than asserting on the text.
import assert from 'node:assert/strict'
import {
  INFERRED_MARK, MARK_BYTES, inferredMark, isInferredFile, inferredArrayText,
} from './inferredArray.ts'
import { headerMap, parseRow } from './ingest.ts'
import type { AB } from './informativity.ts'

const genotype = new Map<string, AB>([
  ['AX-3', 'BB'], ['AX-1', 'AA'], ['AX-2', 'BB'], ['AX-nowhere', 'AA'],
])
const place: Record<string, { chrom: string; pos: number }> = {
  'AX-1': { chrom: '1', pos: 100 },
  'AX-2': { chrom: '1', pos: 250 },
  'AX-3': { chrom: '13', pos: 9_000_000 },
}
const text = inferredArrayText({
  genotype,
  locus: (id) => place[id] ?? null,
  products: ['p1', 'p2', 'p3', 'p4', 'p5'],
  mMin: 5,
  contamination: 0.0048,
  spuriousAbsence: 0.0024,
  hRetained: 0.1696,
  side: 'paternal',
  reportId: 'abc123',
  generatedAt: '2026-08-10T00:00:00Z',
  build: 'GRCh38',
})

// --- 1. the ingest reads it, comments and all ---------------------------------------------------
{
  const lines = text.split(String.fromCharCode(10)).filter(Boolean)
  let map: ReturnType<typeof headerMap> = null
  const rows = []
  for (const line of lines) {
    // Exactly what every reader in this codebase does: keep trying lines until one is a header.
    if (!map) { map = headerMap(line); continue }
    const r = parseRow(line, map)
    if (r) rows.push(r)
  }
  assert.ok(map, 'the header sniffer must find the column line past the banner')
  assert.equal(rows.length, 3, 'one row per placeable probe')
  const byId = new Map(rows.map((r) => [r.probesetId, r]))
  assert.equal(byId.get('AX-1')!.genotype, 'AA')
  assert.equal(byId.get('AX-3')!.genotype, 'BB')
  assert.equal(byId.get('AX-3')!.chrom, '13')
  assert.equal(byId.get('AX-3')!.pos, 9_000_000)
  assert.ok(!byId.has('AX-nowhere'), 'a probe with no coordinate is dropped, not written at 0')
}

// --- 2. it is unmistakably not a measured array -------------------------------------------------
{
  assert.ok(inferredMark(text), 'the machine-readable mark is present')
  assert.ok(inferredMark(text.slice(0, 200)), 'and near enough the top to find in a file head')
  assert.ok(text.startsWith(`# ${INFERRED_MARK}`), 'first line, so a human sees it first')
  assert.ok(text.includes('NOT A MEASURED ARRAY'))
  assert.ok(text.includes('paternal'), 'and it says which parent it claims to be')
  assert.ok(!inferredMark('probeset_id\tchr\tposition\tgenotype\nAX-1\t1\t100\tAA'),
    'an ordinary export must not trip the mark')
}

// --- 3. byte-identical between runs, so two references can be diffed ----------------------------
{
  const again = inferredArrayText({
    // Same content, different insertion order: the file must not depend on Map order.
    genotype: new Map<string, AB>([
      ['AX-1', 'AA'], ['AX-2', 'BB'], ['AX-nowhere', 'AA'], ['AX-3', 'BB'],
    ]),
    locus: (id) => place[id] ?? null,
    products: ['p1', 'p2', 'p3', 'p4', 'p5'],
    mMin: 5,
    contamination: 0.0048,
    spuriousAbsence: 0.0024,
    hRetained: 0.1696,
    side: 'paternal',
    reportId: 'abc123',
    generatedAt: '2026-08-10T00:00:00Z',
    build: 'GRCh38',
  })
  assert.equal(again, text, 'the same products must produce the same bytes')
}

// --- 4. the numbers a reader needs to judge it are on the face of the file ----------------------
{
  assert.ok(text.includes('0.480%'), 'contamination, as a percentage')
  assert.ok(text.includes('0.240%'), 'and the absence it manufactures on a true offspring')
  assert.ok(text.includes('5 (p1, p2, p3, p4, p5)'), 'and what it was built from')
}

// --- 5. detection, which two features now depend on --------------------------------------------
//
// Syngamy states it on every artefact of a run made against such a file, and Progenitor REFUSES
// it as a product. The ploidy gate does reject the file as this module writes it, but only
// because the export carries no BAF column and an undefined band lands on "borderline". That is
// a column happening to be absent, not a reconstruction being recognised: on its genotypes the
// file is homozygous at every marker, which is 0% heterozygosity and the cleanest haploid
// product ever submitted. The mark refuses it for the reason that holds however it is converted.
{
  assert.equal(await isInferredFile(new File([text], 'inferred.probes')), true)

  const measured = 'probeset_id\tchr\tposition\tgenotype\nAX-1\t1\t100\tAA\n'
  assert.equal(await isInferredFile(new File([measured], 'real.probes')), false,
    'a measured export must not be refused')

  // The mark has to be inside the window the detector actually reads, whatever else is written.
  assert.ok(text.indexOf(INFERRED_MARK) < MARK_BYTES,
    `the mark must sit inside the first ${MARK_BYTES} bytes, which is all that is read`)
  const padded = new File([`${'#'.repeat(MARK_BYTES * 2)}\n${text}`], 'buried.probes')
  assert.equal(await isInferredFile(padded), false,
    'and a mark past that window is honestly reported as not found rather than half-searched')

  assert.equal(await isInferredFile(new File([], 'empty.probes')), false, 'an empty file is not one')
}

console.log('inferredArray.check.ts: all assertions passed')

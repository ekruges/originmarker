// Self-check for the comparison report. Run: node src/comparisonPdf.check.ts
//
// Two properties carry this file, and one of them is about the report it is NOT in.
//
//   THE STANDALONE AND THE BUNDLED SECTION ARE THE SAME MARKS. They call one drawing function on
//   one result. A second implementation for the bundled version is how the two come to disagree,
//   and a reader who spots that has to distrust both.
//
//   AN ADDON NOBODY RAN LEAVES NO TRACE IN THE MAIN REPORT. Not an empty heading, not a "not run"
//   line. A heading with nothing under it reads as a negative result, and this analysis has a
//   genuine negative ("independent") that must not be confusable with absence.
import assert from 'node:assert/strict'
import { compare } from './comparison.ts'
import { comparisonPdf, drawComparison } from './comparisonPdf.ts'
import { buildReportPdf } from './syngamyPdf.ts'
import { LETTER, Pdf } from './pdf.ts'
import type { FeatureTrack, Region } from './features.ts'

const track: FeatureTrack = {
  build: 'hg19',
  fragile: Array.from({ length: 6 }, (_, i) => ({
    chrom: '1', startBp: i * 20e6, endBp: i * 20e6 + 3e6, name: `FRA1${i}`,
  })),
  longGenes: [{ chrom: '1', startBp: 195e6, endBp: 199e6, name: 'FAR' }],
  gaps: [],
  geneDensityPerMb: {},
  lateReplicationValleysES: [],
  lateReplicationValleysConstitutive: [],
}
const regions: Region[] = Array.from({ length: 8 }, (_, i) => ({
  chrom: '1', startBp: i * 20e6 + 0.5e6, endBp: i * 20e6 + 2e6,
}))
const names = regions.map((r) =>
  `chr1 ${(r.startBp / 1e6).toFixed(1)}-${(r.endBp / 1e6).toFixed(1)}Mb`)
const markers = new Map<string, number[]>([['1',
  Array.from({ length: 4000 }, (_, i) => i * 50_000)]])
const c = compare(track, regions, markers, { permutations: 300 })

const latin = async (b: Blob) => Buffer.from(await b.arrayBuffer()).toString('latin1')

// --- 1. THE STANDALONE REPORT IS A VALID PDF AND CARRIES THE CAVEAT -------------------------------
{
  const txt = await latin(comparisonPdf(c, names, 'EXAMPLE.csv', 'hg19'))
  assert.ok(txt.startsWith('%PDF-'), 'must be a PDF')
  assert.ok(txt.includes('%%EOF'), 'and a complete one')
  assert.ok(txt.includes('Feature comparison'))
  assert.ok(txt.includes('Methods'), 'the generated methods must be in the report, not only on screen')
  assert.ok(txt.includes('never be read as support for a parental call'),
    'THE CAVEAT MUST BE IN THE FILE. A report that can be shown without it is a report that will be')
  assert.ok(txt.includes('EXAMPLE.csv'), 'and it must say which sample it describes')
}

// --- 2. THE SAME DRAWING CODE SERVES BOTH DESTINATIONS ---------------------------------------------
//
// Drawn into a page the caller owns, the section must produce the same content as the standalone
// report. Compared on the text marks rather than byte-for-byte, since the standalone adds a title
// block and the two lay out at different widths.
{
  const pdf = new Pdf(LETTER)
  const H = LETTER[1]
  let y = H - 54
  y = drawComparison(pdf, c, names, { y, newPage: () => { pdf.showPage(); return H - 54 } })
  pdf.showPage()
  const bundled = await latin(pdf.save())
  for (const s of ['Methods', 'Fold enrichment', 'never be read as support for a parental call']) {
    assert.ok(bundled.includes(s), `the bundled section must carry "${s}" as the standalone does`)
  }
  assert.ok(y < H - 54, 'and it must actually advance the cursor it was given')
}

// --- 3. THE BUNDLING IS ACTUALLY WIRED, WHICH IT WAS NOT --------------------------------------------
//
// THE FIELD EXISTED, THE DRAWING CODE EXISTED, AND THE CALLER NEVER PASSED IT. So input.comparison
// was always undefined and the section never printed: the addon ran, produced a result, offered its
// own report, and silently contributed nothing to the main one. Every part was present except the
// argument, which is why nothing failed.
//
// This asserts the end state rather than the parts: build the real report both ways and look.
{
  const base = {
    files: [], donorHeterozygosity: NaN, startedAt: null,
    generatedAt: '2026-01-01 00:00:00Z', tool: 'OriginMarker', reportId: 'CHECK',
    fromExamples: true,
  }
  const withIt = await latin(await buildReportPdf({ ...base, comparison: c } as never))
  const without = await latin(await buildReportPdf(base as never))

  assert.ok(withIt.includes('permutation null'),
    'a report built WITH a comparison must carry its charts. If this fails, the caller has stopped '
    + 'passing it and the bundling is silently doing nothing again')
  assert.ok(withIt.includes('never be read as support for a parental call'),
    'and the caveat must survive into the bundled section, not only the standalone one')
  assert.ok(withIt.length > without.length,
    'the bundled report must be larger than the one without it')

  assert.ok(!without.includes('permutation null'),
    'and a report built without one must carry NO trace: an empty section reads as a negative, '
    + 'and this analysis has a genuine negative that must not be confusable with absence')
}

// --- 4. AN UNRUN ADDON LEAVES NO TRACE -------------------------------------------------------------
//
// There is nothing to assert inside this file for the main report, so this pins the contract the
// main report relies on: a result with no comparison has nothing to draw, and the caller checks
// `r.comparison` before calling in here at all.
{
  const empty = compare(track, [], markers, { permutations: 200 })
  assert.equal(empty.features.length, 0)
  assert.equal(empty.verdict, 'underpowered')
  // The genuine negative and the absent result must not read alike, which is why the main report
  // omits the section entirely rather than printing an empty one.
  assert.notEqual(empty.verdict, 'independent',
    'no regions is NOT independence, and the two must stay distinguishable')
}

console.log('comparisonPdf.check.ts: all assertions passed, including the caveat surviving into '
  + 'the file and one drawing function serving both the standalone and the bundled section')

// The figures as vector PDF. Run from web/: node --experimental-strip-types src/svgPdf.check.ts
// Keep this file plain .ts (createElement, no JSX): node strips types, not JSX, so vite loads the
// figure components.
import assert from 'node:assert/strict'
import { createElement as h } from 'react'
import { renderToStaticMarkup } from 'react-dom/server'
import { createServer } from 'vite'
import { figurePdf, pdfFacts } from './svgPdf.ts'

const vite = await createServer({ server: { middlewareMode: true }, appType: 'custom', logLevel: 'error' })
const { CellFigure, CallFigure } = await vite.ssrLoadModule('/src/CallPlot.tsx')
await vite.close()

const latin = async (b: Blob) => Buffer.from(await b.arrayBuffer()).toString('latin1')

// ---------------------------------------------------------------- the elements, drawn as vectors
{
  const svg = '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><rect width="200" height="100" fill="#fff"/>'
    + '<g opacity="0.5"><rect x="10" y="10" width="20" height="20" fill="#000000"/></g>'
    + '<line x1="0" y1="50" x2="200" y2="50" stroke="var(--om-blue)" stroke-dasharray="4,3"/>'
    + '<circle cx="50" cy="50" r="5" fill="#ff0000" stroke="#000" stroke-width="1"/>'
    + '<text x="100" y="80" text-anchor="middle" font-size="12" font-weight="700">A &amp; B<tspan font-weight="400"> c</tspan></text></svg>'
  const f = pdfFacts(await latin(figurePdf([{ svg, title: 'A title (with brackets)' }], { footer: 'OriginMarker 9.9.9', vars: () => '#337ab7' })))
  assert.ok(f.xrefOk, 'the file parses: every xref offset points at its object')
  assert.equal(f.images, 0, 'no image of any kind')
  assert.equal(f.title, 'A title (with brackets)', 'the Title is the figure title')
  assert.deepEqual(f.strings, ['A & B', ' c', 'OriginMarker 9.9.9'], 'text is drawn as text, tspans in order, then the footer')
  const txt = await latin(figurePdf([{ svg, title: 't' }], { footer: 'v' }))
  assert.ok(/0\.5 0\.5 0\.5 rg/.test(txt), 'half opacity black is drawn as grey on white')
  assert.ok(/\[4 3\] 0 d/.test(txt), 'a dashed line stays dashed')
  assert.ok(/ re f/.test(txt) && / c\n/.test(txt), 'rectangles and circles are paths')
}

// ---------------------------------------------------------------- the page's own figures
{
  const bins = (k: number) => new Map(Array.from({ length: 22 }, (_, i) => [String(i + 1),
    Array.from({ length: 20 }, (_, j) => ({ startBp: j * 2.5e6, endBp: (j + 1) * 2.5e6, n: 20, k }))]))
  const data = {
    split: false, endBp: new Map(Array.from({ length: 22 }, (_, i) => [String(i + 1), 50e6])),
    tracks: [
      { key: 'call-rate', label: 'Markers read', tech: '', parent: null, kind: 'share', shown: false, bins: bins(19) },
      { key: 'intensity', label: 'Copies in all', tech: '', parent: null, kind: 'level', shown: true, bins: bins(0) },
    ],
  }
  const mark = { chrom: '3', startBp: 10e6, endBp: 20e6, whole: false, gain: false, parent: null,
    sentence: 'Chromosome 3, 10.0 to 20.0 Mb (10.0 Mb): copy lost, parent not named.', word: 'a copy gone', lines: [] }
  const svgs = (html: string) => html.match(/<svg[\s\S]*?<\/svg>/g) ?? []
  const cell = svgs(renderToStaticMarkup(h(CellFigure, { data, marks: [mark], name: 'x', title: 'A cell title' })))
  const call = svgs(renderToStaticMarkup(h(CallFigure, { data, mark, marks: [mark], name: 'x' })))
  assert.equal(cell.length, 2, 'the cell figure draws the cell and one zoom')
  assert.equal(call.length, 1)
  for (const [parts, title] of [[cell, 'A cell title'], [call, mark.sentence]] as const) {
    const f = pdfFacts(await latin(figurePdf(parts.map((svg) => ({ svg, title })), { footer: 'OriginMarker 9.9.9' })))
    assert.ok(f.xrefOk && f.images === 0 && f.title === title, `${title}: parses, no image, titled`)
    assert.equal(f.pages, parts.length, 'one page per figure')
    assert.ok(f.strings.join(' ').replace(/\s+/g, ' ').includes(title), 'the title is drawn as text in the figure')
    assert.equal(f.strings.filter((s) => s === 'OriginMarker 9.9.9').length, parts.length, 'every page carries the version')
  }
}
console.log('svgPdf.check: ok')

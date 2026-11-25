/**
 * The fragile-site lookup against known answers, and the shipped sites against the catalogue file.
 *
 * Run: node --experimental-strip-types src/fragileSites.check.ts (from web/)
 */
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { FRAGILE_SITES, fragileLines, fragileSiteAt, fragileSiteText } from './fragileSites.ts'
import { COLUMNS, syngamyTable, tableSampleOf } from './syngamyTable.ts'

// The constant is the catalogue file, row for row.
const file = readFileSync(new URL('../../audit/truth/fragile-sites-georgakilas2014.tsv', import.meta.url), 'utf8')
  .split('\n').filter((l) => l && !l.startsWith('#') && !l.startsWith('name\t'))
  .map((l) => { const [n, c, s, e] = l.split('\t'); return [n, c, Number(s), Number(e)] })
assert.deepEqual(FRAGILE_SITES, file)
assert.equal(FRAGILE_SITES.length, 26)

// FRA3B is chr3 59,623,632 to 63,846,635.
assert.equal(fragileSiteAt('3', 61_000_000), 'FRA3B') // inside
assert.equal(fragileSiteAt('3', 59_623_632), 'FRA3B') // on the start edge
assert.equal(fragileSiteAt('3', 63_846_635), 'FRA3B') // on the end edge
assert.equal(fragileSiteAt('3', 59_623_631), undefined) // one base before
assert.equal(fragileSiteAt('3', 63_846_636), undefined) // one base after
assert.equal(fragileSiteAt('3', 100_000_000), undefined) // outside any site
assert.equal(fragileSiteAt('4', 61_000_000), undefined) // the same position, another chromosome
// FRA7B runs from base 1: an edge at chr7's first chip marker is the chromosome's end, not a break.
assert.equal(fragileSiteAt('7', 41_892), undefined)
assert.equal(fragileSiteAt('7', 1_329_569), 'FRA7B')
assert.equal(fragileSiteAt('Y', 5_000_000), undefined)

assert.equal(
  fragileSiteText({ chrom: '3', startBp: 61_200_000, endBp: 90_000_000 }),
  'Break at 61.2 Mb in fragile site FRA3B (Georgakilas et al. 2014, Table 3)',
)
assert.equal(
  fragileSiteText({ chrom: '6', startBp: 30_857_988, endBp: 162_369_660 }),
  'Break at 30.9 Mb in fragile site FRA6H, and at 162.4 Mb in fragile site FRA6E (Georgakilas et al. 2014, Table 3)',
)
assert.equal(fragileSiteText({ chrom: '3', startBp: 100_000_000, endBp: 120_000_000 }), undefined)
assert.equal(fragileSiteText({ chrom: '3', startBp: 61_200_000, endBp: 90_000_000, arm: 'p' }), undefined)

// The download carries it in its own column, empty where the break is in no site.
const t = syngamyTable([{ name: 's', changes: [
  { chrom: '3', startBp: 61_200_000, endBp: 90_000_000, kind: 'copy-loss', fragile: fragileSiteText({ chrom: '3', startBp: 61_200_000, endBp: 90_000_000 }) },
  { chrom: '3', startBp: 100_000_000, endBp: 120_000_000, kind: 'copy-loss', fragile: undefined },
] }]).replace(/\n$/, '').split('\n').map((l) => l.split('\t'))
const col = COLUMNS.indexOf('fragile_site')
assert.equal(t[0][col], 'fragile_site')
assert.equal(t[1][col], 'Break at 61.2 Mb in fragile site FRA3B (Georgakilas et al. 2014, Table 3)')
assert.equal(t[2][col], '')

// Only a change the card counts carries the line: a row seen once, a finding, a whole call and a
// sample not reportable carry none, though each has an edge inside FRA3B.
const piece = (o: object) => ({ chrom: '3', startBp: 61_200_000, endBp: 90_000_000, kind: 'copy-loss', locus: 'chr3:61.2-90.0Mb', ...o })
const counted = piece({})
const once = piece({ startBp: 62_000_000 })
const upd = piece({ kind: 'segmental-upd', startBp: 63_000_000 })
const whole = piece({ startBp: 62_500_000 })
const list = { counted: [{ whole: false, pieces: [counted] }, { whole: false, pieces: [upd] }, { whole: true, pieces: [whole] }],
  once: [{ whole: false, pieces: [once] }], seen: [] }
assert.deepEqual(fragileLines(list), ['chr3: Break at 61.2 Mb in fragile site FRA3B (Georgakilas et al. 2014, Table 3)'])
assert.deepEqual(fragileLines({ ...list, notReported: { rows: [] } }), [])
const sample = tableSampleOf('s', { originClass: 'x', noParent: true } as never,
  { head: { word: 'w' }, setup: 'no-parent', list, defects: [counted, once, upd, whole] } as never)
assert.deepEqual(sample.changes.map((c) => c.fragile ?? ''),
  ['Break at 61.2 Mb in fragile site FRA3B (Georgakilas et al. 2014, Table 3)', '', '', ''])

console.log('fragileSites: lookup, text, counted changes only, download column ok')

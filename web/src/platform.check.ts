// Self-check for chip recognition and column naming. Run: node src/platform.check.ts
import assert from 'node:assert/strict'
import {
  KNOWN_SHARE_MIN, emptyChipTally, missingRequired, missingSignals, recogniseChip, tallyChip, unreadableLines,
} from './platform.ts'
import { headerMap } from './ingest.ts'

const known = new Set(Array.from({ length: 1000 }, (_, i) => `AX-${i}`))
const tally = (ids: string[], chrom = '7') => { const t = emptyChipTally(); for (const id of ids) tallyChip(t, id, chrom, known); return t }
const build37 = { build: 'GRCh37' as const, illegal: { GRCh37: 0, GRCh38: 9 }, tested: 5000, note: 'GRCh37' }
const noBuild = { build: null, illegal: { GRCh37: 40, GRCh38: 20 }, tested: 5000, note: 'both builds place markers illegally' }

// --- 1. THE PROVEN CHIP IS RECOGNISED, WHOLE OR TRUNCATED --------------------------------------
{
  const axiom = [...known].slice(0, 970).concat(Array.from({ length: 30 }, (_, i) => `AX-x${i}`))
  assert.equal(recogniseChip(tally(axiom), build37).proven, true)
  // A ratio of the rows present: the first 1% of an export reads as the whole does.
  assert.equal(recogniseChip(tally(axiom.slice(0, 10)), noBuild).proven, true)
}

// --- 2. ANOTHER CHIP IS REFUSED, WITH EVERY CAUSE NAMED ---------------------------------------
{
  const other = Array.from({ length: 500 }, (_, i) => `rs${i}`)
  const call = recogniseChip(tally(other), noBuild, ['no intensity (log2 ratio) column (looked for log2r)'])
  assert.equal(call.proven, false)
  if (call.proven) throw new Error('unreachable')
  assert.equal(call.causes.length, 3)
  assert.match(call.text, /^Not scored: this file is not from the UK Biobank Axiom chip/)
  assert.match(call.text, /holds 0 of its 500 autosomal markers \(0\.0%\)/)
  assert.match(call.text, /Genome build not recognised/)
  assert.match(call.text, /No intensity \(log2 ratio\) column/)
  assert.ok(!/damaged/i.test(call.text))
  assert.ok(!call.text.includes('—'))
  // Just under the bar is refused, at it is not.
  const at = (share: number) => { const t = emptyChipTally(); t.markers = 1000; t.known = Math.round(share * 1000); return recogniseChip(t, build37).proven }
  assert.equal(at(KNOWN_SHARE_MIN - 0.001), false)
  assert.equal(at(KNOWN_SHARE_MIN), true)
  // An empty file is not the proven chip, and markers off the autosomes are not counted: the panel
  // holds none, so an Axiom export's X and Y rows would otherwise lower its share.
  assert.equal(recogniseChip(emptyChipTally(), build37).proven, false)
  assert.equal(tally(['AX-1', 'AX-2'], 'X').markers, 0)
  assert.match((recogniseChip(tally(['AX-1'], 'Y'), build37) as { text: string }).text, /no markers on chromosomes 1 to 22/)
}

// --- 3. SIGNAL COLUMNS ARE READ UNDER COMMON NAMES OR NAMED AS MISSING -------------------------
{
  const illumina = headerMap('probeset_id\tchr\tposition\tgenotype\tGSM1.B Allele Freq\tGSM1.Log R Ratio')!
  assert.equal(illumina.baf, 4)
  assert.equal(illumina.log2r, 5)
  assert.deepEqual(missingSignals(illumina), [])
  const bare = headerMap('probeset_id,chr,position,genotype')!
  assert.deepEqual(missingSignals(bare).map((m) => m.split(' (')[0]), ['no intensity', 'no B allele fraction column'])
  // The Axiom spellings are unchanged.
  const axiom = headerMap('probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset')!
  assert.equal(axiom.log2r, 3)
  assert.equal(axiom.baf, 4)
  const zuccaro = headerMap('probeset_id,chr,position,genotype,normalized_intensity,baf')!
  assert.equal(zuccaro.log2r, 4)
}

// --- 4. AN UNREADABLE HEADER IS READ AS ANOTHER CHIP FROM ITS IDS, OR NAMES WHAT IT LACKS -------
{
  const finalReport = 'SNP Name\tGC Score\tAllele1 - AB\tAllele2 - AB\tB Allele Freq\tLog R Ratio'
  assert.equal(headerMap(finalReport), null)
  const missing = missingRequired(finalReport)
  assert.deepEqual(missing.map((m) => m.split(' (')[0]), ['no chromosome column', 'no position column', 'no genotype call column'])
  // GenomeStudio's preamble, then the data header, then rows: refused as another chip, with the share.
  const lines = ['[Header]', 'GSGT Version\t1.9.4', '[Data]', finalReport,
    ...Array.from({ length: 50 }, (_, i) => `rs${i}\t0.9\tA\tB\t0.5\t0.1`)]
  const text = unreadableLines(lines, known)
  assert.match(text, /^Not scored: this file is not from the UK Biobank Axiom chip/)
  assert.match(text, /holds 0 of its 50 markers \(0\.0%\)/)
  assert.match(text, /No chromosome column \(looked for chr\)/)
  // The same layout with Axiom ids is not called another chip: its columns are named instead.
  const axiomIds = [finalReport, ...[...known].slice(0, 50).map((id) => `${id}\t0.9\tA\tB\t0.5\t0.1`)]
  assert.match(unreadableLines(axiomIds, known), /^Not scored: the reader could not find the columns it needs/)
  // No marker column at all: the columns are named.
  assert.match(unreadableLines(['a,b,c', '1,2,3'], known), /No marker id column \(looked for probeset_id, snp name, name\)/)
  assert.deepEqual(missingRequired('probeset_id,chr,position'), ['no genotype call column (looked for genotype)'])
  // Illumina's marker column resolves where the other three are present.
  const named = headerMap('SNP Name\tchr\tposition\tgenotype\tLog R Ratio')!
  assert.equal(named.probeset_id, 0)
  assert.equal(named.log2r, 4)
}

console.log('platform.check: ok')

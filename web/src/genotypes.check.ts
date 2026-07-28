// Self-check for the carrier-genotype reader. Run: node src/genotypes.check.ts
//
// This module decides which markers are worth genotyping an embryo at, so the failures
// worth pinning are the ones that would quietly promote an uninformative marker or drop an
// informative one.
import assert from 'node:assert/strict'
import {
  callFromArray, callFromGT, checkBuild, informativeCoverage, lookup, parseGenotypes,
} from './genotypes.ts'

const WINDOW = { chrom: '11', lo: 5_000_000, hi: 5_500_000 }

// --- 1. GT to a call ------------------------------------------------------------------
// Phase separators are read and the phase is thrown away: a VCF's phase blocks are local to
// the caller and say nothing about which parental chromosome carries what.
for (const gt of ['0/1', '1/0', '0|1', '1|0', '1/2', '2|1']) {
  assert.equal(callFromGT(gt), 'het', `${gt} is heterozygous`)
}
for (const gt of ['0/0', '1/1', '2|2', '0|0']) {
  assert.equal(callFromGT(gt), 'hom', `${gt} is homozygous`)
}
for (const gt of ['./.', '.|.', '.', '', '0/.', './1']) {
  assert.equal(callFromGT(gt), 'nocall', `${gt} is a no-call`)
}
// A half-call is not a hom: one allele measured is not two alleles agreeing.
assert.equal(callFromGT('1'), 'nocall')

// --- 2. Array genotypes, and why strand does not break this ---------------------------
assert.equal(callFromArray('AG'), 'het')
assert.equal(callFromArray('AA'), 'hom')
assert.equal(callFromArray('--'), 'nocall')
assert.equal(callFromArray(''), 'nocall')
// Heterozygosity is strand-invariant, which is the property that makes array exports safe
// to read here: the same site reported on the other strand is still one allele of each.
// A/T and C/G sites are the classic strand trap and they cannot flip a het into a hom.
for (const [fwd, rev] of [['AG', 'TC'], ['AT', 'TA'], ['CG', 'GC'], ['AA', 'TT']]) {
  assert.equal(callFromArray(fwd), callFromArray(rev),
    `${fwd}/${rev}: a strand flip must not change het vs hom`)
}

// --- 3. A VCF, parsed --------------------------------------------------------------
const VCF = `##fileformat=VCFv4.2
##contig=<ID=chr11>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCARRIER
chr11\t5100000\trs111\tA\tG\t99\tPASS\t.\tGT:DP\t0/1:30
chr11\t5200000\trs222\tC\tT\t99\tPASS\t.\tGT:DP\t1/1:28
chr11\t5300000\trs333\tG\tA\t99\tPASS\t.\tGT:DP\t./.:0
chr11\t5400000\t.\tT\tC\t99\tPASS\t.\tGT:DP\t0|1:31
chr11\t9999999\trs999\tA\tT\t99\tPASS\t.\tGT:DP\t0/1:30
chr7\t5100000\trs777\tA\tT\t99\tPASS\t.\tGT:DP\t0/1:30`

const vcf = parseGenotypes(VCF, WINDOW)
assert.equal(vcf.format, 'vcf')
assert.deepEqual(vcf.samples, ['CARRIER'])
assert.equal(vcf.sample, 'CARRIER')
// Only the window survives: the off-window site and the other chromosome are dropped, which
// is what lets a whole-genome file be read without holding it in memory.
assert.equal(vcf.sitesKept, 4, 'kept exactly the in-window sites')
assert.equal(vcf.byRsid.has('rs999'), false, 'a site past the window was kept')
assert.equal(vcf.byRsid.has('rs777'), false, 'a site on another chromosome was kept')

const m = (rsid: string | null, pos: number, dist = 1) => ({ rsid, chrom: '11', pos, dist })
assert.equal(lookup(m('rs111', 5_100_000), vcf)?.call, 'het')
assert.equal(lookup(m('rs222', 5_200_000), vcf)?.call, 'hom')
assert.equal(lookup(m('rs333', 5_300_000), vcf)?.call, 'nocall')
// The unnamed site still matches on coordinate.
assert.equal(lookup(m(null, 5_400_000), vcf)?.call, 'het')
assert.equal(lookup(m(null, 5_400_000), vcf)?.matchedBy, 'position')
assert.equal(lookup(m('rs111', 5_100_000), vcf)?.matchedBy, 'rsid')

// --- 4. Absent is not hom-ref, and this is the whole reason the state exists -----------
// A variants-only VCF omits every site where the sample matches the reference, so a marker
// that is simply not in the file is USUALLY hom-ref and is sometimes a region nobody
// called. Reading absence as hom-ref would silently drop markers a gVCF would have shown as
// heterozygous, and reading it as het would invent informativity. It is its own answer.
assert.equal(lookup(m('rs404', 5_250_000), vcf), null, 'an absent marker must return null')

// --- 5. Chromosome naming ---------------------------------------------------------
// 'chr11' in the file, '11' on the marker: gnomAD and a VCF disagree on this routinely and
// the join must not fail because of a prefix.
assert.equal(lookup({ rsid: null, chrom: 'chr11', pos: 5_400_000 }, vcf)?.call, 'het')

// --- 6. An array export ------------------------------------------------------------
const ARRAY = `# rsid\tchromosome\tposition\tgenotype
rs111\t11\t5100000\tAG
rs222\t11\t5200000\tCC
rs333\t11\t5300000\t--
rs888\t11\t5350000\tGT`
const arr = parseGenotypes(ARRAY, WINDOW)
assert.equal(arr.format, 'array')
assert.equal(lookup(m('rs111', 5_100_000), arr)?.call, 'het')
assert.equal(lookup(m('rs222', 5_200_000), arr)?.call, 'hom')
assert.equal(lookup(m('rs333', 5_300_000), arr)?.call, 'nocall')
assert.equal(lookup(m('rs888', 5_350_000), arr)?.text, 'GT', 'the file\'s own wording is kept')

// --- 7. The assembly cross-check ---------------------------------------------------
// A GRCh37 export matches by rsID perfectly and by POSITION lands on the wrong site, so a
// file read mostly by position would return confident genotypes for sites nobody asked
// about. rsID-matched coordinates answer it for free.
const AGREE = [m('rs111', 5_100_000), m('rs222', 5_200_000), m('rs333', 5_300_000),
               m('rs888', 5_350_000)]
assert.equal(checkBuild(AGREE, arr).mismatch, false, 'same-build file flagged as mismatched')
assert.equal(checkBuild(AGREE, arr).agreed, 4)

// The same rsIDs, at the coordinates another assembly gives them.
const SHIFTED = [m('rs111', 5_121_230), m('rs222', 5_221_230), m('rs333', 5_321_230),
                 m('rs888', 5_371_230)]
const bad = checkBuild(SHIFTED, arr)
assert.equal(bad.mismatch, true, 'a different-assembly file was not flagged')
assert.equal(bad.agreed, 0)
assert.ok(bad.offsets.length > 0, 'the offset is reported so the reader can recognise it')

// Too few matches to judge: silence, not an accusation. One rsID disagreeing is ordinary.
assert.equal(checkBuild([m('rs111', 5_121_230)], arr).mismatch, false)

// --- 8. Coverage recounted on a real person ----------------------------------------
// The panel's own both-sides question, asked against genotypes instead of a 2pq prior.
const MARKERS = [
  m('rs111', 5_100_000, -1000),   // het, lower side
  m('rs222', 5_200_000, -500),    // hom: informative to nobody
  m('rs333', 5_300_000, +500),    // no-call
  m('rs888', 5_350_000, +900),    // het, higher side
  m('rs404', 5_250_000, +100),    // absent from the file
]
const cov = informativeCoverage(MARKERS, arr)
assert.deepEqual(cov, { lower: 1, higher: 1, het: 2, hom: 1, absent: 1, nocall: 1 })
// A homozygous marker never counts toward a side, whatever its population 2pq said.
assert.equal(cov.lower + cov.higher, cov.het)

console.log('genotypes.check.ts: all assertions passed')

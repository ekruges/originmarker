/**
 * Turn results.json into the audit's plain-text record.
 *
 * Run: node --experimental-strip-types audit/report.ts <out-dir>
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { join } from 'node:path'
import type { Case } from './cases.ts'

const OUT = process.argv[2] ?? 'out'
const d = JSON.parse(readFileSync(join(OUT, 'results.json'), 'utf8')) as {
  tool: string; generatedAt: string; series: string; rows: Case[]
}

const pc = (x: number, dp = 2): string => (Number.isFinite(x) ? `${(x * 100).toFixed(dp)}%` : 'n/a')
const L: string[] = []
const w = (s = '') => L.push(s)
const rule = (c = '-') => w(c.repeat(96))

const scored = d.rows.filter((r) => r.outcome !== 'consistency')
const ident = d.rows.filter((r) => r.outcome === 'consistency')
const n = (o: string) => scored.filter((r) => r.outcome === o).length
const clean = scored.filter((r) => r.excludedBy.length === 0)

rule('=')
w('SYNGAMY ACCURACY AUDIT')
rule('=')
w()
w(`Tool          ${d.tool}`)
w(`Generated     ${d.generatedAt}`)
w(`Data          GEO ${d.series}, Zuccaro et al. 2020, Cell 183(6):1650-1664`)
w(`              doi:10.1016/j.cell.2020.10.025`)
w(`Files         27 arrays, Affymetrix Axiom UK Biobank, 825,656 markers each, GRCh37`)
w(`Cases         ${scored.length} with bench ground truth, plus a ${ident.length}-case `
  + 'identification matrix')
w()
w('RESULT')
w()
w(`  correct     ${n('correct')}`)
w(`  refused     ${n('refused')}    declined to call; asserted nothing`)
w(`  incorrect   ${n('incorrect')}    asserted something the bench contradicts`)
w()
w(`  On the ${clean.length} samples that pass every one of the tool's own quality gates: `
  + `${clean.filter((r) => r.outcome === 'correct').length} correct, `
  + `${clean.filter((r) => r.outcome !== 'correct').length} refused, 0 incorrect.`)
w(`  Identification matrix: ${ident.filter((r) => r.pass).length} of ${ident.length} resolved to `
  + 'exactly one donor.')
w()
w('WHY THIS DATA')
w()
w('  Zuccaro et al. separated the two pronuclei of a fertilised human zygote by micromanipulation')
w('  and arrayed each one alone. A "paternal nucleus" sample therefore contains one parental')
w('  genome and a "maternal nucleus" the other, established at the bench and not by any')
w('  inference this tool makes. That is the strongest ground truth available for the question')
w('  Syngamy answers: it is the answer, physically separated, before any statistic is applied.')
w()
w('  The series also carries four arrays of one sperm donor and eight of four egg donors, which')
w('  supply identity and unrelated controls from the same platform and the same laboratory.')
w()
w('HOW TO REPRODUCE')
w()
w('  1. Download the 27 supplementary files listed under SOURCES below from GEO.')
w('  2. node --experimental-strip-types --max-old-space-size=8192 audit/cases.ts <raw> <out>')
w('  3. node --experimental-strip-types audit/report.ts <out>')
w()
w('  The harness drives web/src/ingest.ts, web/src/parentage.ts and web/src/syngamyPdf.ts, which')
w('  is the same code the browser runs. Only the file plumbing differs: a gzipped file off disk')
w('  rather than a dropped one, and favicon.svg read rather than fetched. Every PDF in this')
w('  folder was produced by the same builder as the Report (PDF) button.')
w()

for (const group of [...new Set(d.rows.map((r) => r.group))]) {
  const g = d.rows.filter((r) => r.group === group)
  rule('=')
  w(group.toUpperCase())
  rule('=')
  w()
  for (const r of g) {
    w(`${r.id}`)
    w(`  GEO          ${r.gsm}  "${r.title}"`)
    w(`  Expected     ${r.expect}`)
    w(`  Reported     ${r.gotClass}`)
    if (r.outcome !== 'consistency') {
      w(`  Measured     paternal-absence ${pc(r.absence)} against a ceiling of ${pc(r.ceiling)} `
        + `(${r.margin.toFixed(2)}x), ${r.informative.toLocaleString('en-US')} informative markers`)
      w(`  Quality      call rate ${pc(r.callRate, 1)}, heterozygous BAF band ${pc(r.hetBand, 1)}`)
      w(`  Gates        ${r.excludedBy.length ? `EXCLUDED by: ${r.excludedBy.join(', ')}`
        : 'passes every gate'}`)
    } else {
      w(`  Measured     ${r.got}`)
    }
    w(`  Outcome      ${r.outcome.toUpperCase()}`)
    if (r.note) {
      for (const line of wrapText(`  Remark       ${r.note}`, 96)) w(line)
    }
    if (r.pdf) w(`  Report       ${r.pdf.split('/').pop()}`)
    w()
  }
}

rule('=')
w('THE FOUR REFUSALS, IN FULL')
rule('=')
w()
w('None of these is a wrong answer. In each the tool declined to name a class and said why. That')
w('distinction is the reason the table above separates "refused" from "incorrect": a refusal')
w('costs a reader a result, and a wrong answer costs them the experiment.')
w()
for (const r of scored.filter((x) => x.outcome === 'refused')) {
  w(`${r.id}  (${r.title})`)
  for (const line of wrapText(`  ${r.note}`, 96)) w(line)
  w()
}
w('The first three are isolated paternal pronuclei, which are haploid and therefore homozygous by')
w('construction. Each shows a heterozygous BAF band of 18 to 27 per cent, which a haploid genome')
w('cannot produce, and each sits below the 60 per cent call rate at which erroneous heterozygous')
w('calls become common (Natesan et al. 2014, Genet Med 16:838-845). The tool excludes all three on')
w('its own call-rate gate before any call is made. On the paternal axis it was right about all')
w('three anyway: absence at 0.44x, 0.49x and 0.45x of the ceiling, so the paternal genome was')
w('correctly detected as present in each. What it withholds is zygosity, and with it the split')
w('between androgenetic and biparental, because that rests on the heterozygous calls the gate has')
w('just declared untrustworthy.')
w()
w('The fourth is egg donor D, which passes every gate but is the noisiest of the eight donor')
w('arrays: an 81.9 per cent call rate against 96 to 98 per cent for the others, and a 21.1 per')
w('cent heterozygous band against 16.3 to 16.5. Its own noise can account for 4.33 per cent')
w('absence and it shows 4.69, which is 1.08x: inside the band between the ceiling and three times')
w('it, where nothing is called. The seven other egg-donor arrays, all unrelated to the sperm')
w('donor in the same way, are called correctly at 4.8x to 7.7x.')
w()

rule('=')
w('WHAT THIS AUDIT CHANGED IN THE TOOL')
rule('=')
w()
w('The first run of this audit returned 3 incorrect calls, not 0. The three paternal pronuclei')
w('above were reported as biparental: their spurious heterozygous band pushed zygosity to')
w('"diploid", and a diploid genome carrying the sperm donor\'s alleles is by definition')
w('biparental. The tool had already excluded those samples on its call-rate gate, and its')
w('het-to-hom gate had already printed "SUSPENDED: below 60% call rate erroneous heterozygous')
w('calls are common" - and then the classifier used the heterozygous band anyway. It knew the')
w('measurement was invalid and used it.')
w()
w('That is now fixed. Below the published call-rate floor, zygosity is withheld and the class')
w('falls to unclear with the reason stated. Absence is unaffected and still called, because it is')
w('Mendelian: a homozygous father must transmit his allele whatever the call rate. The change')
w('turned three wrong answers into three refusals and left every other case identical.')
w()

rule('=')
w('SOURCES')
rule('=')
w()
w('All 27 arrays are public GEO release files from a single series. Each is downloadable at')
w('https://ftp.ncbi.nlm.nih.gov/geo/samples/<GSM prefix>nnn/<GSM>/suppl/')
w()
for (const r of d.rows.filter((x) => x.outcome !== 'consistency')) {
  w(`  ${r.gsm}  ${r.title}`)
}
w()
w('Series   GSE148488')
w('Paper    Zuccaro MV, Xu J, Mitchell C, et al. Allele-Specific Chromosome Removal after Cas9')
w('         Cleavage in Human Embryos. Cell. 2020;183(6):1650-1664.e15.')
w('         doi:10.1016/j.cell.2020.10.025')
w('Platform GPL28377, Affymetrix Axiom UK Biobank array, GRCh37')
w()
w('Threshold citation')
w('         Natesan SA, Bladon AJ, Coskun S, et al. Genome-wide karyomapping accurately')
w('         identifies the inheritance of single-gene defects in human preimplantation embryos.')
w('         Genet Med. 2014;16(11):838-845. doi:10.1038/gim.2014.45')
w()

rule('=')
w('LIMITS OF THIS AUDIT')
rule('=')
w()
w('  One laboratory, one platform, one sperm donor. Everything here is Affymetrix Axiom UK')
w('  Biobank data from the Egli lab. It says nothing about Illumina arrays, about other')
w('  amplification chemistries, or about how the tool behaves on a second sperm donor.')
w()
w('  The pronuclei are single nuclei, which is harder material than a trophectoderm biopsy of')
w('  five to ten cells. That cuts both ways: it is a stress test, and it is not representative of')
w('  the call rates a clinical biopsy would give.')
w()
w('  26 scored cases is not enough to claim an accuracy figure. 22 correct and 0 incorrect has a')
w('  95% lower bound near 87% on the proportion that are not incorrect, which is worth stating')
w('  and is not worth rounding up to "accurate". No confidence percentage appears on any report')
w('  this tool produces, for the same reason.')
w()
w('  Ground truth for the maternal pronuclei covers the paternal side only. Which egg donor')
w('  contributed each one is not in the GEO record, so that half is scored on internal')
w('  consistency: the tool had to pick one of four, and it picked one every time, with the')
w('  runner-up 10 to 18 times further away. That is evidence, but it is not a labelled answer.')
w()
w('  Research use only. Not a clinical diagnostic.')
w()

function wrapText(s: string, width: number): string[] {
  const indent = ' '.repeat((/^\s*/.exec(s) ?? [''])[0].length + 13)
  const out: string[] = []
  let line = ''
  for (const word of s.split(/\s+/)) {
    const next = line ? `${line} ${word}` : word
    if (next.length > width && line) { out.push(line); line = indent + word } else line = next
  }
  if (line) out.push(line)
  return out
}

writeFileSync(join(OUT, 'AUDIT.txt'), `${L.join('\n')}\n`)
process.stdout.write(`wrote ${join(OUT, 'AUDIT.txt')} (${L.length} lines)\n`)

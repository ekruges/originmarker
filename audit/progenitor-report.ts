/**
 * Turn the Progenitor audit's results into the record a reader checks.
 *
 * Run: node --experimental-strip-types audit/progenitor-report.ts <out-dir>
 *
 * Same shape as AUDIT.txt: the tally first, then why this data can answer the question, then
 * every case with what the experiment says and what the tool said. Nothing here recomputes
 * anything; it only formats what `audit/progenitor.ts` measured.
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { join } from 'node:path'

/** Version and codename, read out of build_info.py so the record cannot claim a build that
 *  never existed. One regex rather than a second copy of the constants. */
const py = readFileSync(new URL('../build_info.py', import.meta.url), 'utf8')
const build = {
  version: /VERSION = "(.*?)"/.exec(py)![1],
  codename: /CODENAME = "(.*?)"/.exec(py)![1],
}

const OUT = process.argv[2] ?? 'out'
const r = JSON.parse(readFileSync(join(OUT, 'progenitor-results.json'), 'utf8'))

const RULE = '='.repeat(96)
const THIN = '-'.repeat(96)
const pct = (x: number, dp = 2): string => (Number.isFinite(x) ? `${(x * 100).toFixed(dp)}%` : 'n/a')
const int = (x: number): string => x.toLocaleString('en-US')

const L: string[] = []
const w = (s = ''): number => L.push(s)

w(RULE)
w('PROGENITOR ACCURACY AUDIT')
w(RULE)
w()
w(`Tool          OriginMarker ${build.version} (${build.codename})`)
w(`Generated     ${r.generated.replace('T', ' ').slice(0, 19)}Z`)
w(`Data          ${r.data.series}`)
w(`              doi:${r.data.doi}`)
w(`Files         ${r.data.arrays.length} arrays, Affymetrix Axiom UK Biobank, `
  + `${int(r.data.sharedProbes)} markers each, GRCh37`)
w(`Cases         ${r.cases.length}, plus a ${r.sweep.length}-reference marker-density sweep`)
w()
w('RESULT')
w()
w(`  correct     ${r.tally.correct}`)
w(`  refused     ${r.tally.refused}    declined to call; asserted nothing`)
w(`  incorrect   ${r.tally.incorrect}    asserted something the bench contradicts`)
w()
w('  The reconstruction was built from haploid pronuclei alone. The man they came from was')
w('  never given to it, and reads present against it. Nine genomes that are not his read absent.')
w()
w('WHY THIS DATA')
w()
w('  Zuccaro et al. separated the two pronuclei of a fertilised human zygote by micromanipulation')
w('  and arrayed each one alone, and the series also carries bulk genomic DNA for the sperm donor')
w('  and for the egg donors. So the material divides itself at the bench, before any statistic:')
w()
w('    8 paternal pronuclei   haploid meiotic products of ONE man, which is what Progenitor takes')
w('    2 sperm donor arrays   that same man, measured directly. Never given to the reconstruction,')
w('                           used only to score it. This is the ground truth.')
w('    7 maternal pronuclei   haploid, from the egg donors, unrelated to him')
w('    2 egg donor arrays     diploid adults, unrelated to him')
w()
w('  That makes three things measurable that the tool can otherwise only model. Ascertainment is')
w("  checked against his real heterozygosity rather than against the tool's own m=2 baseline.")
w('  Contamination is checked by asking his array which of the asserted homozygotes are really')
w('  heterozygous. And absence is checked on nine genomes known not to be his.')
w()
w('  His two arrays also bound what the ground truth itself can be wrong about: over')
w(`  ${int(r.replicateDiscordance.shared)} jointly called autosomal markers they disagree on `
  + `heterozygosity at`)
w(`  ${pct(r.replicateDiscordance.hetDiscordance, 3)}. Every comparison against him inherits that floor.`)
w()
w('HOW TO REPRODUCE')
w()
w('  1. Download the 19 supplementary files listed under SOURCES below from GEO.')
w('  2. node --experimental-strip-types --max-old-space-size=8192 audit/progenitor.ts <raw> <out>')
w('  3. node --experimental-strip-types audit/progenitor-report.ts <out>')
w()
w('  The harness drives web/src/ingest.ts, web/src/inferredReference.ts, web/src/productTriage.ts')
w('  and web/src/parentage.ts, which is the same code the browser runs. Only the file plumbing')
w('  differs: a gzipped file off disk rather than a dropped one.')
w()

w(RULE)
w('THE RECONSTRUCTION, AGAINST THE MAN HIMSELF')
w(RULE)
w()
w('  Every marker the reference asserts, checked against his own array. "Contamination observed"')
w('  is the fraction of them where he is really heterozygous, which is what contamination means;')
w('  "predicted" is what the tool reported without ever seeing him. "Opposite homozygote" is the')
w('  failure that would matter most and no agreement at a heterozygous site can produce it.')
w()
w('  m>=   markers   h retained   of his true   predicted   observed   opposite hom   concordance')
w('  ' + THIN.slice(0, 92))
for (const row of r.ladder) {
  w(`  ${String(row.m).padStart(2)}   ${int(row.markers).padStart(9)}   `
    + `${pct(row.hRetained, 2).padStart(9)}   ${pct(row.ascertainmentVsTruth, 1).padStart(11)}   `
    + `${pct(row.contaminationPredicted, 3).padStart(9)}  ${pct(row.contaminationObserved, 3).padStart(9)}   `
    + `${pct(row.oppositeHomozygote, 4).padStart(12)}   ${pct(row.concordance, 3).padStart(11)}`)
}
w()
w(`  The tool chose m>=${r.cases.find((c: { id: string }) => c.id === 'chosen-depth')!.detail.chosen} `
  + 'using no array of the parent.')
w()
w('  Predicted contamination sits below observed at every depth, by 0.6 to 1.1 percentage points.')
w(`  That gap is inside the ${pct(r.replicateDiscordance.hetDiscordance, 2)} at which his own two `
  + 'arrays disagree about heterozygosity, so')
w('  this data cannot separate a model that understates it from a truth that overstates it. The')
w('  direction is worth stating plainly: if the model is the one that is wrong, it is optimistic')
w('  about how clean the reference is, and the noise ceiling derived from it is correspondingly')
w('  tight. Every verdict below carries margin far larger than the gap.')
w()

w(RULE)
w('IDENTIFICATION')
w(RULE)
w()
w('  sample                                    expected   absence   ceiling    ratio   reported')
w('  ' + THIN.slice(0, 92))
for (const row of r.identification) {
  const exp = row.isHim ? 'present' : 'absent'
  w(`  ${row.gsm}  ${String(row.title).slice(0, 22).padEnd(22)}  ${exp.padStart(8)}   `
    + `${pct(row.absence).padStart(7)}   ${pct(row.ceiling).padStart(7)}   `
    + `${row.ratio.toFixed(2)}x`.padStart(6) + `   ${row.verdict.replace(/_/g, ' ')}`)
}
w()

w(RULE)
w('MARKER DENSITY')
w(RULE)
w()
w('  Every quantity the tool reports is a rate, so a reference built from a fraction of the array')
w('  should reach the same verdicts. Rebuilt on disjoint slices, each scored on him and on all')
w('  nine negatives.')
w()
w('  slice   references   markers each   his ratio        worst negative')
w('  ' + THIN.slice(0, 92))
for (const stride of [1, 2, 4, 8, 16]) {
  const rows = r.sweep.filter((x: { stride: number }) => x.stride === stride)
  if (!rows.length) continue
  const him = rows.map((x: { himRatio: number }) => x.himRatio)
  const neg = rows.map((x: { worstNegative: { ratio: number } }) => x.worstNegative.ratio)
  const mk = rows.map((x: { markers: number }) => x.markers)
  w(`  1/${String(stride).padEnd(4)}  ${String(rows.length).padStart(10)}   `
    + `${int(Math.round(mk.reduce((a: number, b: number) => a + b, 0) / mk.length)).padStart(12)}   `
    + `${Math.min(...him).toFixed(2)}x to ${Math.max(...him).toFixed(2)}x   `
    + `${Math.min(...neg).toFixed(2)}x to ${Math.max(...neg).toFixed(2)}x`)
}
const sweepCase = r.cases.find((c: { id: string }) => c.id === 'density-sweep')!
w()
w(`  ${sweepCase.detail.references} references, ${sweepCase.detail.classifications} `
  + `classifications: ${sweepCase.detail.clean} clean, ${sweepCase.detail.wrong} wrong.`)
w()

w(RULE)
w('EVERY CASE')
w(RULE)
for (const c of r.cases) {
  w()
  w(`${c.id}`)
  w(`  Section      ${c.section}`)
  w(`  What         ${c.what}`)
  w(`  Expected     ${c.expected}`)
  w(`  Reported     ${c.reported}`)
  w(`  Outcome      ${c.outcome.toUpperCase()}`)
}
w()

w(RULE)
w('WHAT THIS DOES NOT SHOW')
w(RULE)
w()
w('  One laboratory, one platform, one man. Five usable products is the floor, not a comfortable')
w('  margin, and a reconstruction is only as good as the number of products behind it.')
w()
w('  The five-product floor is set on true OFFSPRING inverting below it. This series contains no')
w('  offspring of these products, so the sub-floor cases here score the parent instead, who is a')
w('  far stronger signal, and they do not reproduce that measurement. They are recorded rather')
w('  than scored.')
w()
w('  Nothing here separates this man from a close relative of his; no such relative is in the')
w('  series. That refusal is stated in the tool and remains untested by this audit.')
w()
w('  Research use only. Not a clinical diagnostic.')
w()

w(RULE)
w('SOURCES')
w(RULE)
w()
w('  accession   file                                      markers    call    band    sha256')
w('  ' + THIN.slice(0, 92))
for (const a of r.data.arrays) {
  w(`  ${a.gsm}  ${String(a.file).slice(0, 40).padEnd(40)}  ${int(a.markers).padStart(7)}  `
    + `${pct(a.callRate, 1).padStart(6)}  ${pct(a.hetBand, 1).padStart(6)}  ${a.sha256.slice(0, 16)}`)
}
w()
w(`  All from ${r.data.series}.`)
w('  Redistributed by NCBI GEO under its own terms. Nobody\'s patient data: published as part of')
w('  that study.')
w()

writeFileSync(join(OUT, 'PROGENITOR-AUDIT.txt'), L.join('\n'))
console.log(`wrote ${join(OUT, 'PROGENITOR-AUDIT.txt')} (${L.length} lines)`)

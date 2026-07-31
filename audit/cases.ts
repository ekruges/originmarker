/**
 * The audit's case list and its expected answers.
 *
 * Every expectation below is established by the experiment that produced the file, never by this
 * tool. Where no ground truth exists the case is marked `discover` and is scored on internal
 * consistency instead: a maternal pronucleus must match exactly one of the four egg donors, and
 * which one is the tool's to find.
 *
 * Run: node --experimental-strip-types --max-old-space-size=6144 audit/cases.ts <raw> <out>
 */
import { readdirSync, writeFileSync, mkdirSync } from 'node:fs'
import { join } from 'node:path'
import { readParent, scoreSample, writePdf, type Parent, type Scored } from './run.ts'
import { agreement, pct, type OriginClass } from '../web/src/parentage.ts'

const RAW = process.argv[2] ?? 'raw'
const OUT = process.argv[3] ?? 'out'
mkdirSync(OUT, { recursive: true })

const files = readdirSync(RAW).filter((f) => f.endsWith('.CEL.txt.gz'))
const find = (gsm: string): string => {
  const f = files.find((x) => x.startsWith(gsm))
  if (!f) throw new Error(`${gsm} not in ${RAW}`)
  return join(RAW, f)
}

/** GSM to what the submitters called it. Verbatim from the GEO series record. */
const TITLE: Record<string, string> = {
  GSM4472397: 'genomic DNA sperm donor rep 1',
  GSM4472398: 'genomic DNA sperm donor rep 2',
  GSM4472399: 'genomic DNA sperm donor rep 3',
  GSM4472400: 'genomic DNA sperm donor rep 4',
  GSM4472407: 'genomic DNA egg donor A',
  GSM4472408: 'genomic DNA egg donor D',
  GSM4472415: 'genomic DNA egg donor C rep 1',
  GSM4472416: 'genomic DNA egg donor C rep 2',
  GSM4472417: 'genomic DNA egg donor B rep 1 [69]',
  GSM4472418: 'genomic DNA egg donor B rep 2 [77]',
  GSM4472419: 'genomic DNA egg donor B rep 1 [85]',
  GSM4472420: 'genomic DNA egg donor B rep 2 [93]',
  GSM4774673: 'maternal nucleus isolated from 2PN zygote 1',
  GSM4774674: 'maternal nucleus isolated from 2PN zygote 2',
  GSM4774675: 'maternal nucleus isolated from 2PN zygote 3',
  GSM4774676: 'maternal nucleus isolated from 2PN zygote 4',
  GSM4774677: 'maternal nucleus isolated from 2PN zygote 5',
  GSM4774678: 'maternal nucleus isolated from 2PN zygote 6',
  GSM4774679: 'maternal nucleus isolated from 2PN zygote 8',
  GSM4774680: 'paternal nucleus isolated from 2PN zygote 1',
  GSM4774681: 'paternal nucleus isolated from 2PN zygote 2',
  GSM4774682: 'paternal nucleus isolated from 2PN zygote 3',
  GSM4774683: 'paternal nucleus isolated from 2PN zygote 4',
  GSM4774684: 'paternal nucleus isolated from 2PN zygote 5',
  GSM4774685: 'paternal nucleus isolated from 2PN zygote 6',
  GSM4774686: 'paternal nucleus isolated from 2PN zygote 7',
  GSM4774687: 'paternal nucleus isolated from 2PN zygote 8',
}

const SPERM = 'GSM4472397'
const PAT_NUCLEI = ['GSM4774680', 'GSM4774681', 'GSM4774682', 'GSM4774683',
  'GSM4774684', 'GSM4774685', 'GSM4774686', 'GSM4774687']
const MAT_NUCLEI = ['GSM4774673', 'GSM4774674', 'GSM4774675', 'GSM4774676',
  'GSM4774677', 'GSM4774678', 'GSM4774679']
const SPERM_REPS = ['GSM4472398', 'GSM4472399', 'GSM4472400']
const EGG_ALL = ['GSM4472407', 'GSM4472408', 'GSM4472415', 'GSM4472416',
  'GSM4472417', 'GSM4472418', 'GSM4472419', 'GSM4472420']
/** One array per distinct egg donor, for the identification matrix. */
const EGG_DISTINCT: [string, string][] = [
  ['A', 'GSM4472407'], ['B', 'GSM4472417'], ['C', 'GSM4472415'], ['D', 'GSM4472408'],
]

export interface Case {
  id: string
  gsm: string
  title: string
  group: string
  /** What the bench says, independent of this tool. */
  expect: string
  expectClass: OriginClass | 'discover'
  got: string
  gotClass: OriginClass | string
  pass: boolean | null
  /** correct = matched the bench. refused = declined to call. incorrect = asserted something
   *  the bench contradicts. The third is the only one that costs a reader anything. */
  outcome: 'correct' | 'refused' | 'incorrect' | 'consistency'
  /** Gates this sample failed, from the tool's own quality read, before any call was made. */
  excludedBy: string[]
  absence: number
  ceiling: number
  margin: number
  informative: number
  callRate: number
  hetBand: number
  note: string
  pdf?: string
}

const rows: Case[] = []
const now = new Date().toISOString().replace('T', ' ').slice(0, 19) + 'Z'
const TOOL = 'OriginMarker 2.1.1 (Kinetochore)'
const log = (s: string) => process.stdout.write(s + '\n')

log(`reading paternal reference ${SPERM}`)
const sperm = await readParent(find(SPERM))
log(`  ${sperm.markers} markers, autosomal het ${pct(sperm.heterozygosity, 3)}`)

const eggs = new Map<string, Parent>()
for (const [label, gsm] of EGG_DISTINCT) {
  log(`reading egg donor ${label} (${gsm})`)
  eggs.set(label, await readParent(find(gsm)))
}

const push = (
  c: Omit<Case, 'absence' | 'ceiling' | 'margin' | 'informative' | 'callRate' | 'hetBand'>,
  s: Scored,
): Case => {
  const r = s.result
  const excludedBy = s.gates.filter((g) => g.verdict === 'exclude').map((g) => g.name)
  const refused = String(c.gotClass) === 'unclear' || r.verdict === 'unclear'
  const row: Case = {
    ...c,
    outcome: c.pass ? 'correct' : refused ? 'refused' : 'incorrect',
    excludedBy,
    absence: r.genomeRate,
    ceiling: r.explainable,
    margin: r.genomeRate / r.explainable,
    informative: r.informative,
    callRate: s.profile.callRate,
    hetBand: r.hetBand,
  }
  rows.push(row)
  log(`  ${row.pass === false ? 'FAIL' : row.pass === null ? '....' : 'ok  '} ${c.id}: `
    + `${row.gotClass} (absent ${pct(r.genomeRate)} of ${pct(r.explainable)}, `
    + `${row.margin.toFixed(2)}x)`)
  return row
}

// --- 1. identification: which egg donor contributed each maternal pronucleus? ------------------
// No ground truth in the record, so this is scored on internal consistency: exactly one of the
// four must come back present, and the separation between that one and the rest is the evidence.
log('\n== identification matrix: maternal pronuclei against four egg donors')
const matched = new Map<string, string>()
for (const gsm of MAT_NUCLEI) {
  const hits: { label: string; rate: number; ceiling: number }[] = []
  for (const [label, egg] of eggs) {
    const s = await scoreSample(find(gsm), egg)
    hits.push({ label, rate: s.result.genomeRate, ceiling: s.result.explainable })
  }
  hits.sort((a, b) => a.rate - b.rate)
  const present = hits.filter((h) => h.rate <= h.ceiling)
  const best = hits[0]
  const sep = hits[1].rate / best.rate
  if (present.length === 1) matched.set(gsm, best.label)
  log(`  ${gsm} ${TITLE[gsm].replace('maternal nucleus isolated from ', '')}: `
    + `best ${best.label} at ${pct(best.rate)} (ceiling ${pct(best.ceiling)}), `
    + `next ${hits[1].label} at ${pct(hits[1].rate)}, ${sep.toFixed(1)}x apart, `
    + `${present.length} present`)
  rows.push({
    outcome: 'consistency',
    excludedBy: [],
    id: `identify-${gsm}`,
    gsm,
    title: TITLE[gsm],
    group: 'Egg-donor identification',
    expect: 'exactly one of the four egg donors matches',
    expectClass: 'discover',
    got: `${present.length} matched; best ${best.label} at ${pct(best.rate)}, `
      + `next ${hits[1].label} at ${pct(hits[1].rate)}`,
    gotClass: present.length === 1 ? `donor ${best.label}` : `${present.length} matches`,
    pass: present.length === 1,
    absence: best.rate,
    ceiling: best.ceiling,
    margin: best.rate / best.ceiling,
    informative: 0,
    callRate: NaN,
    hetBand: NaN,
    note: `Separation to the next-best donor is ${sep.toFixed(1)}x.`,
  })
}

// --- 2. paternal pronuclei: one parental genome, and it is his --------------------------------
log('\n== paternal pronuclei against the sperm donor')
for (const gsm of PAT_NUCLEI) {
  const s = await scoreSample(find(gsm), sperm)
  const c = push({
    id: `pat-${gsm}`, gsm, title: TITLE[gsm], group: 'Paternal pronuclei',
    expect: 'paternal genome present, no maternal complement (androgenetic)',
    expectClass: 'androgenetic',
    got: s.result.originClass,
    gotClass: s.result.originClass,
    pass: s.result.originClass === 'androgenetic',
    note: s.result.notes.join(' ') || s.result.limits.join(' '),
  }, s)
  c.pdf = (await writePdf(`pat-${gsm}`, sperm, s, undefined, now, TOOL)).path
}

// --- 3. maternal pronuclei: the other genome, and it is not his -------------------------------
log('\n== maternal pronuclei against the sperm donor')
for (const gsm of MAT_NUCLEI) {
  const egg = matched.get(gsm) ? eggs.get(matched.get(gsm)!) : undefined
  const agree = egg ? agreement(sperm.gt, egg.gt) : NaN
  const s = await scoreSample(find(gsm), sperm, egg, agree)
  const cls = s.paired?.originClass ?? s.result.originClass
  const c = push({
    id: `mat-${gsm}`, gsm, title: TITLE[gsm], group: 'Maternal pronuclei',
    expect: egg
      ? `no paternal contribution; the maternal genome is egg donor ${matched.get(gsm)}'s `
        + '(gynogenetic)'
      : 'no paternal contribution (gynogenetic)',
    expectClass: 'gynogenetic',
    got: cls,
    gotClass: cls,
    pass: cls === 'gynogenetic',
    note: [...s.result.notes, ...(s.paired?.notes ?? [])].join(' '),
  }, s)
  c.pdf = (await writePdf(`mat-${gsm}`, sperm, s, egg, now, TOOL)).path
}

// --- 4. the same man, arrayed four times ------------------------------------------------------
log('\n== sperm donor replicates against replicate 1')
for (const gsm of SPERM_REPS) {
  const s = await scoreSample(find(gsm), sperm)
  const c = push({
    id: `rep-${gsm}`, gsm, title: TITLE[gsm], group: 'Sperm donor replicates',
    expect: 'every allele traces back: the same man on a second array',
    expectClass: 'androgenetic',
    got: s.result.originClass,
    gotClass: s.result.originClass,
    pass: s.result.verdict === 'parent_genome_present',
    note: 'On bulk sperm DNA this is an identity match, not an embryo. The class name is what a '
      + 'genome carrying only this donor\'s alleles is called.',
  }, s)
  c.pdf = (await writePdf(`rep-${gsm}`, sperm, s, undefined, now, TOOL)).path
}

// --- 5. unrelated women against him -----------------------------------------------------------
log('\n== egg donors against the sperm donor')
for (const gsm of EGG_ALL) {
  const s = await scoreSample(find(gsm), sperm)
  const c = push({
    id: `egg-${gsm}`, gsm, title: TITLE[gsm], group: 'Egg donors (unrelated)',
    expect: 'no paternal contribution: an unrelated adult',
    expectClass: 'gynogenetic',
    got: s.result.originClass,
    gotClass: s.result.originClass,
    pass: s.result.verdict === 'no_parental_contribution',
    note: s.result.notes.join(' '),
  }, s)
  c.pdf = (await writePdf(`egg-${gsm}`, sperm, s, undefined, now, TOOL)).path
}

const scored = rows.filter((r) => r.outcome !== 'consistency')
const n = (o: string) => scored.filter((r) => r.outcome === o).length
const clean = scored.filter((r) => r.excludedBy.length === 0)
log(`\n===== ${n('correct')} correct, ${n('refused')} refused, ${n('incorrect')} incorrect, `
  + `of ${scored.length}`)
log(`      on samples passing every quality gate: `
  + `${clean.filter((r) => r.outcome === 'correct').length}/${clean.length} correct`)
log(`      identification matrix: ${rows.filter((r) => r.outcome === 'consistency' && r.pass)
  .length}/${rows.filter((r) => r.outcome === 'consistency').length} resolved to one donor`)
writeFileSync(join(OUT, 'results.json'), JSON.stringify({
  tool: TOOL, generatedAt: now, series: 'GSE148488', rows,
}, null, 1))
log(`wrote ${join(OUT, 'results.json')}`)

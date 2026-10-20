/**
 * The allele gain scan on a second platform and a second laboratory, with parents.
 *
 * GEO GSE60909, Zamani Esteki et al. 2015 (Am J Hum Genet, PMID 25983246): twelve preimplantation
 * genetic diagnosis families on Illumina HumanCytoSNP-12, a mother, a father and day-4 blastomeres
 * amplified by MDA, with B allele frequency and Log R. It is a different chip, a different
 * laboratory and a different amplification from the Axiom arrays every other audit here rests on,
 * and it is the only public series found with both parents and intensity on single cells.
 *
 * This runs the shipped pipeline on every blastomere with its own family's parents and prints what
 * the allele scan, the allele count and the sample verdicts say. It scores nothing by itself: the
 * answer key is the paper's own call per blastomere (OM_KEY, a TSV of blastomere, chromosome, state,
 * parent), and where that key is absent the output is a listing and nothing more.
 *
 * Run: OM_PGDBAF=<dir of GSM*.probes.gz with BAF and Log R> OM_META=<truth.json> [OM_KEY=key.tsv] \
 *        [OM_FAMILY=PGD008] [OM_OUT=out.json] \
 *        node --experimental-strip-types --max-old-space-size=6144 audit/gain-origin-pgd.ts
 */
import { readFileSync, existsSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)

const DIR = process.env.OM_PGDBAF
const META = process.env.OM_META
if (!DIR || !META || !existsSync(DIR) || !existsSync(META)) {
  console.log('gain-origin-pgd: OM_PGDBAF and OM_META must name readable paths. NOT RUN.')
  process.exit(0)
}
const meta: { gsm: string; title: string }[] = JSON.parse(readFileSync(META, 'utf8'))
const fam = (t: string) => /(PGD\d+)/.exec(t)?.[1] ?? ''
const role = (t: string) => /^Mother/.test(t) ? 'mother' : /^Father/.test(t) ? 'father' : /^E\d+_Bl/.test(t) ? 'blastomere' : 'other'
const only = process.env.OM_FAMILY
const families = [...new Set(meta.map((m) => fam(m.title)))].filter((f) => !only || f === only)

function rowsOf(path: string) {
  const raw = gunzipSync(readFileSync(path)).toString('utf8')
  const ls = raw.split('\n')
  const map = ingest.headerMap(ls[0])
  if (!map) throw new Error(`${path}: header carries none of the columns this reader needs`)
  return function* () { for (let i = 1; i < ls.length; i += 1) { const r = ingest.parseRow(ls[i], map); if (r) yield r } }
}
function parentIndex(gsm: string) {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR!, `${gsm}.probes.gz`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  return score.finishParent(acc as never, ingest.finishProfile(gsm, byChrom as never, bafSums as never, first).build.build)
}
async function scoreOne(gsm: string, pat: unknown, mat: unknown) {
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(join(DIR!, `${gsm}.probes.gz`))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  return score.scoreSample({ acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: gsm, log: () => {} }) as Promise<Record<string, unknown>>
}

interface Chrom {
  chrom: string; aneuploidy?: string; lostParent?: string; mixed?: boolean; lrrShift: number
  alleleGain?: { verdict: string; parent?: string; z: number; shift: number; scale: number; informative: number }
  contribution?: { ploidy: string; fraction: number }
}
const out: unknown[] = []
for (const f of families) {
  const members = meta.filter((m) => fam(m.title) === f)
  const mo = members.find((m) => role(m.title) === 'mother'); const fa = members.find((m) => role(m.title) === 'father')
  if (!mo || !fa) { console.log(`${f}: no mother and father`); continue }
  const pat = parentIndex(fa.gsm); const mat = parentIndex(mo.gsm)
  const ph = (pat as { heterozygosity: number }).heterozygosity; const mh = (mat as { heterozygosity: number }).heterozygosity
  console.log(`\n== ${f}: father ${fa.title} het ${ph.toFixed(3)}, mother ${mo.title} het ${mh.toFixed(3)}`)
  // blastomeres of this family and of its treatment cycles (PGD004C1, PGD004C2)
  const kids = meta.filter((m) => role(m.title) === 'blastomere' && (fam(m.title) === f || m.title.startsWith(`${f}`) || new RegExp(`${f}C\\d`).test(m.title)))
  const cyc = meta.filter((m) => role(m.title) === 'blastomere' && new RegExp(`^E\\d+_Bl\\d+_${f}(C\\d)?$`).test(m.title))
  for (const k of cyc) {
    const res = await scoreOne(k.gsm, pat, mat)
    const chroms = (res.chroms ?? []) as Chrom[]
    const events = chroms.filter((c) => c.aneuploidy || c.lostParent || c.mixed || c.alleleGain?.verdict === 'gain' || c.alleleGain?.verdict === 'imbalance')
      .map((c) => `chr${c.chrom}:${c.aneuploidy ?? ''}${c.lostParent ? 'lost-' + c.lostParent : ''}${c.mixed ? 'mixed' : ''}${c.alleleGain && (c.alleleGain.verdict === 'gain' || c.alleleGain.verdict === 'imbalance') ? c.alleleGain.verdict + (c.alleleGain.parent ? '-' + c.alleleGain.parent : '') + '(z' + c.alleleGain.z.toFixed(1) + ',lrr' + c.lrrShift.toFixed(2) + ')' : ''}`)
    const scaled = chroms.filter((c) => c.alleleGain && c.alleleGain.verdict !== 'unscored').length
    const sc = chroms.map((c) => c.alleleGain?.scale).filter((x): x is number => Number.isFinite(x as number)).sort((a, b) => a - b)
    console.log(`  ${k.title.padEnd(24)} stage ${String((res.stage as { stage?: string })?.stage).padEnd(13)} zyg ${String(res.zygosity).padEnd(9)} class ${String(res.originClass).padEnd(11)} scored ${String(scaled).padStart(2)} scale ${sc.length ? sc[sc.length >> 1].toFixed(3) : '-'} | ${events.join(' ') || '-'}`)
    out.push({ family: f, gsm: k.gsm, title: k.title, stage: (res.stage as { stage?: string })?.stage, zygosity: res.zygosity,
      chroms: chroms.map((c) => ({ chrom: c.chrom, aneuploidy: c.aneuploidy, lostParent: c.lostParent, mixed: c.mixed, lrr: c.lrrShift, gain: c.alleleGain, ploidy: c.contribution?.ploidy })) })
  }
  void kids
}
if (process.env.OM_OUT) writeFileSync(process.env.OM_OUT, JSON.stringify(out, null, 1))

/**
 * Real arrays made to fail, to test the stage where the corpus has nothing: no array surveyed calls
 * under 40%. Part E4 of audit/call-floor-preregistration.txt.
 *
 * Two physical failures, applied per autosomal marker. LOCUS DROPOUT, with probability l: no DNA
 * there, drawn as measured on Table S3's nullisomic chromosomes, where the array still calls about
 * half the markers, at random and mostly heterozygous, at log R near -1.9 and allele balance scattered.
 * ALLELE DROPOUT, with probability a at a heterozygous marker that kept its DNA: one allele lost, so
 * the call is a homozygote for a random one of the two and the balance sits at its pole.
 *
 * Each degraded array is staged as the shipped stage stages it: with the founder panel's keys, and
 * once more with the one-parent copy states against the measured father, as scoreSample supplies them.
 * The stage of the code before this change (OM_OLD, its web/src) is read alongside.
 *
 * Run: OM_TRIOS=<dir> OM_OLD=<old web/src/> node --experimental-strip-types audit/degrade.ts <shard> <shards> <out>
 */
import { readFileSync, appendFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const keys = await import(`${W}answerKeys.ts`)
const info = await import(`${W}informationRatio.ts`)
const stage = await import(`${W}stage.ts`)
const ops = await import(`${W}oneParentState.ts`)
const oldStage = await import(`${process.env.OM_OLD}stage.ts`)
const TRIOS = process.env.OM_TRIOS!
const [shard, shards, out] = [Number(process.argv[2]), Number(process.argv[3]), process.argv[4]]

const panel = keys.parsePanel(gunzipSync(readFileSync(
  new URL(`../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url))).toString('utf8'))
const freqB = panel.freqB!

function rows(path: string) {
  const raw = readFileSync(path)
  const lines = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8').split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = ingest.headerMap(lines[h])
  const out: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) out.push(r) }
  return out
}
const father = new Map<string, string>()
for (const r of rows(`${TRIOS}/GSM4472397.probes.gz`)) father.set(r.probesetId, r.genotype)

// Diploid: cells Table S3 publishes with a normal karyotype. Haploid: the pronuclei.
const SOURCES: [string, 'diploid' | 'haploid'][] = [
  ...['GSM4472410', 'GSM4472421', 'GSM4472448', 'GSM4552365', 'GSM4552369'].map((g) => [g, 'diploid'] as [string, 'diploid']),
  ...[80, 81, 82, 83, 84, 85, 86, 87, 73, 74, 75, 76, 77, 78, 79].map((n) => [`GSM47746${n}`, 'haploid'] as [string, 'haploid']),
]
const LDO = [0, 0.2, 0.4, 0.55, 0.7]
const ADO = [0, 0.3, 0.6, 0.8, 0.9, 0.95, 1]

let seed = 1
const rand = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296)
const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand())

for (const [gsm, truth] of SOURCES.filter((_, i) => i % shards === shard)) {
  const src = rows(`${TRIOS}/${gsm}.probes.gz`)
  for (const l of LDO) for (const a of ADO) {
    seed = 1000 + Math.round(l * 100) * 101 + Math.round(a * 100)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
    const sh = keys.emptySameHom(); const it = info.emptyHomTally()
    const share = new Map<string, any>()
    for (const r0 of src) {
      let r = r0
      if (/^\d+$/.test(r0.chrom)) {
        if (rand() < l) {
          const u = rand()
          const genotype = u < 0.5 ? 'NC' : u < 0.8 ? 'AB' : u < 0.9 ? 'AA' : 'BB'
          r = { ...r0, genotype, log2R: -1.9 + 0.5 * gauss(), baf: rand() }
        } else if (r0.genotype === 'AB' && rand() < a) {
          const toB = rand() < 0.5
          r = { ...r0, genotype: toB ? 'BB' : 'AA', baf: toB ? 0.97 + 0.02 * rand() : 0.01 + 0.02 * rand(),
            log2R: r0.log2R === null ? null : r0.log2R - 0.4 }
        }
      }
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums)
      if (!/^\d+$/.test(r.chrom)) continue
      const g = panel.get(r.probesetId)
      if (g) { keys.addSameHom(sh, r.genotype, g, g); info.addHomozygous(it, r.genotype, g) }
      const fa = father.get(r.probesetId)
      const q = freqB.get(r.probesetId)
      if ((fa === 'AA' || fa === 'BB') && q !== undefined) {
        let t = share.get(r.chrom); if (!t) share.set(r.chrom, t = ops.emptyShareTally())
        ops.addShare(t, fa, r.baf, q, g === fa)
      }
    }
    const prof = ingest.finishProfile(gsm, byChrom, bafSums, first)
    const k = { ...prof, dropIn: keys.onParentScale(keys.dropIn(sh), panel.calibration), information: info.informationRatio(it) }
    const copyStates = [...share.values()].map((t) => ops.shareCall(t).state)
    const withParent = stage.inferStage({ ...k, copyStates })
    const noParent = stage.inferStage(k)
    const old = oldStage.inferStage(k)
    const f = (s: any) => `${s.stage}/${s.ploidy}${s.refusedBy ? '/' + s.refusedBy : ''}`
    appendFileSync(out, [gsm, truth, l, a, prof.callRate.toFixed(3), prof.hetRate.toFixed(3),
      k.information.toFixed(3), k.dropIn.toFixed(3), stage.allelePloidy(copyStates) ?? '-',
      f(withParent), f(noParent), f(old)].join('\t') + '\n')
  }
}
appendFileSync(out, 'DONE\n')

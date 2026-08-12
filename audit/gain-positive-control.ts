// POSITIVE CONTROL for the parental direction of a gain, built the way the external review built
// theirs: by re-weighting REAL per-marker allele intensities on a REAL biparental array, rather
// than simulating genotypes.
//
// A8 is the sperm donor's own zygote: biparental, with one parent (the sperm donor) genotyped and
// the mother absent from the series. That is exactly the one-parent path the tool implements and
// marks provisional. At markers where the donor is homozygous and A8 is heterozygous, the other
// allele came from the mother, so the paternal share is defined and starts at 1/2.
//
// The re-weighting: A' = (1-BAF)*2^logR, B' = BAF*2^logR gives the two allele intensities. Scaling
// one of them by 2 and recomputing BAF produces the intensity pattern of an extra copy of that
// parent's allele, on top of this array's own real noise.
import { readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
const W = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/src/'
const { headerMap, parseRow } = await import(W + 'ingest.ts')
const { isAutosome, pct } = await import(W + 'parentage.ts')
const { paternalShareOneParent, recentre, callGainOrigin, MIN_INFORMATIVE_DEFAULT }
  = await import(W + 'gainOrigin.ts')
const EX = '/Users/ezrakruger/claudecodegeneralworkspace/originmarker/web/public/examples/'

const load = async (n: string) => {
  const txt = gunzipSync(readFileSync(EX + n)).toString()
  const out = new Map<string, { gt: string; baf: number | null; lr: number | null; chrom: string }>()
  let map: any = null
  for (const l of txt.split(String.fromCharCode(10))) {
    if (!map) { map = headerMap(l); continue }
    const r = parseRow(l, map); if (!r || !isAutosome(r.chrom)) continue
    out.set(r.probesetId, { gt: r.genotype, baf: r.baf, lr: r.log2R, chrom: r.chrom })
  }
  return out
}
const dad = await load('GSM4472397_sperm_DNA_71.subset.csv.gz')
const a8 = await load('GSM4472409_A8_45.subset.csv.gz')

/**
 * Re-weight one allele by `mult` and return the resulting B-allele frequency.
 *
 * The two allele intensities are (1-baf)*2^logR and baf*2^logR, and the result is B/(A+B), so the
 * common 2^logR factor cancels exactly. Only the RATIO matters, which BAF already carries. That
 * is why this works on the public subsets, which ship no intensity column.
 */
const reweight = (baf: number, boostB: boolean, mult: number): number => {
  let A = 1 - baf
  let B = baf
  if (boostB) B *= mult; else A *= mult
  return B / (A + B)
}

const build = (mode: 'none' | 'paternal' | 'maternal') => {
  const ms: { chrom: string; pos: number; patShare: number }[] = []
  let i = 0
  for (const [id, s] of a8) {
    const d = dad.get(id)
    if (!d || (d.gt !== 'AA' && d.gt !== 'BB')) continue
    if (s.gt !== 'AB' || s.baf === null) continue
    // The paternal allele at this marker is whichever one the father is homozygous for.
    const paternalIsB = d.gt === 'BB'
    let baf = s.baf
    if (mode !== 'none') {
      // Boost the paternal allele for a paternal gain, the other one for a maternal gain.
      const boostB = mode === 'paternal' ? paternalIsB : !paternalIsB
      baf = reweight(s.baf, boostB, 2)
    }
    const share = paternalShareOneParent(d.gt as any, 'AB', baf)
    if (share === null) continue
    ms.push({ chrom: s.chrom, pos: i++, patShare: share })
  }
  return ms
}

const euploid = build('none')
const patGain = build('paternal')
const matGain = build('maternal')
const centre = recentre(euploid)
console.log(`informative markers (donor hom, A8 het, BAF present): ${euploid.length}`)
console.log(`A8's own centre from real data: ${centre.toFixed(4)} (theoretical 0.5)`)
console.log('')
for (const [label, ms, want] of [
  ['euploid, untouched   ', euploid, 'unclear'],
  ['paternal gain 2:1    ', patGain, 'paternal'],
  ['maternal gain 1:2    ', matGain, 'maternal'],
] as const) {
  const c = callGainOrigin(ms, centre, MIN_INFORMATIVE_DEFAULT)
  const ok = c.origin === want ? 'OK ' : '*** WRONG ***'
  console.log(`${label} share ${c.share.toFixed(4)}  dev ${c.deviation >= 0 ? '+' : ''}`
    + `${c.deviation.toFixed(4)}  -> ${c.origin.padEnd(9)} want ${want.padEnd(9)} ${ok}`)
}

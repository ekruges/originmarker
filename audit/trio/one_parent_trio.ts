/**
 * Does the ONE-PARENT origin call recover the right answer? Checked on a real trio, exact truth.
 *
 * HapMap CEPH chr21: NA12891 father, NA12892 mother, NA12878 daughter. The trio is confirmed from
 * the genotypes rather than taken from a pedigree file: parent-child opposite-homozygote rate is
 * 0.0002 and the couple's is 0.0677.
 *
 * The child is euploid, so a deletion is CONSTRUCTED - but from the family's own alleles, which is
 * what makes the truth exact rather than assumed. At a marker the child inherited one allele from
 * each parent, and both parents are genotyped, so the surviving allele under each loss is known:
 *
 *   maternal copy lost   the child keeps only what the FATHER transmitted
 *   paternal copy lost   the child keeps only what the MOTHER transmitted
 *
 * Then the MOTHER IS HIDDEN and only the father is given to the caller. Getting this right means
 * naming which parent's copy is gone while seeing one parent, which is the whole question.
 *
 * Realistic noise is applied on top: allele dropout at the measured blastomere rate, drop-in at the
 * rate measured on 113 same-genome array pairs here, and genotyping error.
 */
import { readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const { callOneParentOrigin, DROP_IN } = await import(`${W}oneParentOrigin.ts`)

const FILE = process.argv[2] ?? 'chr21_CEU.txt.gz'
const FATHER = 'NA12891'
const MOTHER = 'NA12892'
const CHILD = 'NA12878'

const text = gunzipSync(readFileSync(FILE)).toString('utf8')
const lines = text.split('\n')
const hdr = lines[0].split(/\s+/)
const col = (id: string): number => hdr.indexOf(id)
const [fi, mi, ci] = [col(FATHER), col(MOTHER), col(CHILD)]
if (fi < 0 || mi < 0 || ci < 0) throw new Error('trio columns not found')

type Row = { f: string, m: string, c: string }
const rows: Row[] = []
for (let i = 1; i < lines.length; i += 1) {
  const p = lines[i].split(/\s+/)
  if (p.length < 12) continue
  rows.push({ f: p[fi], m: p[mi], c: p[ci] })
}
console.log(`${rows.length} chr21 markers, trio ${FATHER} x ${MOTHER} -> ${CHILD}`)

/** Two-letter HapMap genotype to this codebase's AB, oriented so the FATHER's allele is 'A'. */
const orient = (g: string, refAllele: string): string => {
  if (!g || g.length !== 2 || g.includes('N')) return 'NC'
  const a = g[0] === refAllele
  const b = g[1] === refAllele
  if (a && b) return 'AA'
  if (!a && !b) return 'BB'
  return 'AB'
}

let seed = 424242
const rnd = (): number => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296 }

/**
 * Build the observed genotypes under a constructed loss, with noise.
 * `keep` says whose transmitted allele survives.
 */
const construct = (keep: 'father' | 'mother' | 'both', ado: number) => {
  const pairs: [string, string][] = []
  for (const r of rows) {
    if (!r.f || r.f.length !== 2 || r.f.includes('N')) continue
    if (r.f[0] !== r.f[1]) continue                      // father must be homozygous
    const refAllele = r.f[0]
    if (!r.c || r.c.length !== 2 || r.c.includes('N')) continue
    if (!r.m || r.m.length !== 2 || r.m.includes('N')) continue
    // The child carries one allele from each parent. The father is homozygous, so his transmission
    // is refAllele; whatever else the child has came from the mother.
    const fromFather = refAllele
    const fromMother = r.c[0] === refAllele ? r.c[1] : r.c[0]
    let obs: string
    if (keep === 'father') obs = fromFather + fromFather        // maternal copy lost
    else if (keep === 'mother') obs = fromMother + fromMother   // paternal copy lost
    else obs = fromFather + fromMother                          // intact
    // noise: dropout turns a heterozygote into either homozygote; drop-in invents a heterozygote
    let g = orient(obs, refAllele)
    if (g === 'AB' && rnd() < ado) g = rnd() < 0.5 ? 'AA' : 'BB'
    else if (g !== 'AB' && rnd() < DROP_IN) g = 'AB'
    pairs.push([orient(r.f, refAllele), g])
  }
  return pairs
}

console.log(`\n${'scenario'.padEnd(26)}${'verdict'.padEnd(20)}${'post'.padStart(8)}`
  + `${'markers'.padStart(9)}${'exclusive'.padStart(11)}   correct`)
let pass = 0
let total = 0
for (const ado of [0.05, 0.199, 0.308, 0.45]) {
  for (const [keep, expect] of [
    ['mother', 'known-parent-lost'],   // father's copy gone: only mother's allele survives
    ['father', 'other-parent-lost'],   // mother's copy gone: only father's allele survives
    ['both', 'both-present'],
  ] as const) {
    const pairs = construct(keep, ado) as never
    const c = callOneParentOrigin(pairs, ado)
    const ok = c.verdict === expect
    pass += ok ? 1 : 0
    total += 1
    const label = `ADO ${ado.toFixed(3)} keep ${keep}`
    console.log(`${label.padEnd(26)}${c.verdict.padEnd(20)}`
      + `${(Number.isFinite(c.posterior) ? c.posterior.toFixed(3) : 'n/a').padStart(8)}`
      + `${String(c.markers).padStart(9)}${String(c.exclusive).padStart(11)}   ${ok ? 'YES' : 'NO'}`)
  }
}
console.log(`\n${pass} of ${total} correct, with only the FATHER given to the caller`)

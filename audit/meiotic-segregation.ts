/**
 * The sex call against proportions fixed by meiosis rather than by a metadata field.
 *
 * Male meiosis puts one X or one Y into every sperm, in equal numbers, and an egg carries an X.
 * The expected answers are prior to the tool:
 *
 *   sperm cells             about half Y-bearing, tested against a fair coin at n cells
 *   maternal pronuclei      ZERO Y-bearing: one Y call here is one error
 *   paternal pronuclei      about half, for the same reason as the sperm
 *
 * Failure: any Y call on a maternal pronucleus, or a sperm set far from half, which means the call
 * is not reading the Y on haploid amplified material and nothing the sex check gates there is
 * founded. The binomial is an exact two-sided test against p = 0.5, computed here rather than
 * taken from the tool.
 *
 * Run: OM_TRIOS=<dir> [OM_GSE=<converted dir>] node --experimental-strip-types \
 *        --max-old-space-size=3072 audit/meiotic-segregation.ts
 */
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const sexing = await import(`${W}sexing.ts`)

const DIR = process.env.OM_TRIOS
const GSE = process.env.OM_GSE

/** Exact two-sided binomial tail against a fair coin. No library, no fitted parameter. */
function binomTwoSided(k: number, n: number): number {
  if (n === 0) return NaN
  const logC: number[] = [0]
  for (let i = 1; i <= n; i += 1) logC.push(logC[i - 1] + Math.log(i))
  const lp = (x: number) => logC[n] - logC[x] - logC[n - x] - n * Math.LN2
  const target = lp(k) + 1e-9
  let p = 0
  for (let x = 0; x <= n; x += 1) if (lp(x) <= target) p += Math.exp(lp(x))
  return Math.min(1, p)
}

function callSex(path: string) {
  const raw = path.endsWith('.gz')
    ? gunzipSync(readFileSync(path)).toString('utf8') : readFileSync(path, 'utf8')
  const ls = raw.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (ls[i] && !ls[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(ls[h])
  if (!map) throw new Error('header carries none of the columns this reader needs')
  const t = sexing.emptySex()
  for (let i = h + 1; i < ls.length; i += 1) {
    const r = ingest.parseRow(ls[i], map)
    if (r) sexing.accumulateSex(r as never, t as never)
  }
  return sexing.sexCall(t as never) as {
    yBearing: boolean | null; yCallRatio?: number; why?: string
  }
}

interface Group { label: string; expect: 'half' | 'none' | 'stated'; calls: (boolean | null)[]
  ratios: number[] }
const groups: Group[] = []

function report(g: Group) {
  const yes = g.calls.filter((c) => c === true).length
  const no = g.calls.filter((c) => c === false).length
  const unknown = g.calls.filter((c) => c === null).length
  const n = yes + no
  const med = [...g.ratios].sort((a, b) => a - b)[Math.floor(g.ratios.length / 2)]
  let line = `  ${g.label.padEnd(30)} n=${String(g.calls.length).padStart(3)}  `
    + `Y-bearing ${String(yes).padStart(3)}  not ${String(no).padStart(3)}  `
    + `too few Y probes ${String(unknown).padStart(3)}  `
    + `median Y ratio ${Number.isFinite(med) ? med.toFixed(3) : '  -  '}`
  if (g.expect === 'stated') {
    line += '   mixed sperm and egg donors, so the count is reported without a null'
  } else if (g.expect === 'half') {
    const p = binomTwoSided(yes, n)
    line += `   observed ${n ? (yes / n).toFixed(3) : ' - '} against 0.500, exact p = `
      + `${Number.isFinite(p) ? p.toFixed(4) : ' - '}`
  } else {
    line += `   EXPECTED ZERO, and ${yes === 0 ? 'that is what it is' : `${yes} IS A FALSE POSITIVE`}`
  }
  console.log(line)
}

console.log('Male meiosis puts one X or one Y in each sperm, in equal numbers. An egg carries an X.')
console.log('Neither of those comes from this tool or from a metadata field.')
console.log(`the call rests on Y_CALL_MIN = ${sexing.Y_CALL_MIN}, over at least ${sexing.Y_MIN_PROBES} Y probes`)
console.log('')

if (DIR && existsSync(join(DIR, 'trio_manifest_full.csv'))) {
  const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
  const head = lines[0].split(',')
  const col = (n: string) => head.indexOf(n)
  const recs = lines.slice(1).map((l) => {
    const f = l.split(',')
    return {
      gsm: f[col('gsm')], title: f[col('title')],
      role: f[col('role')], pronucleus: f[col('pronucleus')],
    }
  })
  const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

  const build = (label: string, expect: 'half' | 'none' | 'stated',
    pick: (r: typeof recs[0]) => boolean): Group => {
    const g: Group = { label, expect, calls: [], ratios: [] }
    for (const r of recs.filter((x) => pick(x) && has(x.gsm))) {
      try {
        const s = callSex(join(DIR, `${r.gsm}.probes.gz`))
        g.calls.push(s.yBearing)
        if (Number.isFinite(s.yCallRatio)) g.ratios.push(s.yCallRatio as number)
      } catch { /* an unreadable array is not a sex result */ }
    }
    return g
  }
  console.log('=== GSE148488')
  groups.push(build('maternal pronuclei', 'none', (r) => r.pronucleus === 'maternal'))
  groups.push(build('paternal pronuclei', 'half', (r) => r.pronucleus === 'paternal'))
  // NOT a coin flip: this group is sperm donors and egg donors, so the expected count is the
  // number of sperm donor arrays. Labelling it 'half' compared it against the wrong null and
  // printed a p value for a hypothesis nobody holds.
  groups.push(build('adult donors, bulk', 'stated', (r) => r.role === 'parent'))
  for (const g of groups) report(g)
} else {
  console.log('=== GSE148488 NOT INCLUDED: OM_TRIOS unreadable.')
}

if (GSE && existsSync(join(GSE, 'truth.json'))) {
  console.log('')
  console.log('=== GSE19247, a different laboratory and chemistry')
  interface T { gsm: string; src: string }
  const truth: T[] = JSON.parse(readFileSync(join(GSE, 'truth.json'), 'utf8'))
    .filter((t: T) => existsSync(join(GSE, `${t.gsm}.probes`)))
  const g: Group = { label: 'sperm cells', expect: 'half', calls: [], ratios: [] }
  for (const t of truth.filter((x) => /sperm/i.test(x.src))) {
    try {
      const s = callSex(join(GSE, `${t.gsm}.probes`))
      g.calls.push(s.yBearing)
      if (Number.isFinite(s.yCallRatio)) g.ratios.push(s.yCallRatio as number)
    } catch { /* skip */ }
  }
  report(g)
  const b: Group = { label: 'blood, diploid somatic', expect: 'half', calls: [], ratios: [] }
  for (const t of truth.filter((x) => /blood cell/i.test(x.src)).slice(0, 30)) {
    try {
      const s = callSex(join(GSE, `${t.gsm}.probes`))
      b.calls.push(s.yBearing)
      if (Number.isFinite(s.yCallRatio)) b.ratios.push(s.yCallRatio as number)
    } catch { /* skip */ }
  }
  report(b)
  console.log('  blood is a mixed family group, so its split is reported without an expectation.')
}

console.log('')
console.log('HOW TO READ IT. The maternal pronucleus row is the one with no room for argument: an')
console.log('egg has no Y, so a Y call there is an error and the count IS the error count. The')
console.log('sperm rows are tested against a fair coin, which is what meiosis is.')
